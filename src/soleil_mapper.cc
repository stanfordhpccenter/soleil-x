#include <iostream>
#include <fstream>
#include <regex>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "mappers/default_mapper.h"
#include "realm/logging.h"

#include "config_schema.h"
#include "soleil_mapper.h"

using namespace Legion;
using namespace Legion::Mapping;

//=============================================================================

static Realm::Logger LOG("soleil_mapper");

#define CHECK(cond, ...)                        \
  do {                                          \
    if (!(cond)) {                              \
      LOG.error(__VA_ARGS__);                   \
      exit(1);                                  \
    }                                           \
  } while(0)

#define STARTS_WITH(str, prefix)                \
  (strncmp((str), (prefix), sizeof(prefix) - 1) == 0)

static const Config* find_config(const Task* task) {
  for(; task != NULL; task = task->parent_task) {
    if (strcmp(task->get_task_name(), "work") == 0) {
      const char* ptr = static_cast<const char*>(task->args);
      return reinterpret_cast<const Config*>(ptr + sizeof(uint64_t));
    }
  }
  return NULL;
}

//=============================================================================

// Maps super-tiles to ranks, in row-major order. Within a super-tile, tiles
// are assigned unique processor IDs, in row-major order. The mapper will have
// to match those IDs to real processors (not necessarily one per ID).
class SampleMapping {
public:
  SampleMapping(const Config& config, AddressSpace first_rank)
    : tiles_per_rank_{config.Mapping.tilesPerRank[0],
                      config.Mapping.tilesPerRank[1],
                      config.Mapping.tilesPerRank[2]},
      ranks_per_dim_{config.Mapping.tiles[0] / config.Mapping.tilesPerRank[0],
                     config.Mapping.tiles[1] / config.Mapping.tilesPerRank[1],
                     config.Mapping.tiles[2] / config.Mapping.tilesPerRank[2]},
      first_rank_(first_rank) {}
  AddressSpace get_rank(unsigned x, unsigned y, unsigned z) const {
    return first_rank_ +
      (x / tiles_per_rank_[0]) * ranks_per_dim_[1] * ranks_per_dim_[2] +
      (y / tiles_per_rank_[1]) * ranks_per_dim_[2] +
      (z / tiles_per_rank_[2]);
  }
  unsigned get_proc_id(unsigned x, unsigned y, unsigned z) const {
    return
      (x % tiles_per_rank_[0]) * tiles_per_rank_[1] * tiles_per_rank_[2] +
      (y % tiles_per_rank_[1]) * tiles_per_rank_[2] +
      (z % tiles_per_rank_[2]);
  }
  unsigned num_ranks() const {
    return ranks_per_dim_[0] * ranks_per_dim_[1] * ranks_per_dim_[2];
  }
  unsigned x_tiles() const {
    return tiles_per_rank_[0] * ranks_per_dim_[0];
  }
  unsigned y_tiles() const {
    return tiles_per_rank_[1] * ranks_per_dim_[1];
  }
  unsigned z_tiles() const {
    return tiles_per_rank_[2] * ranks_per_dim_[2];
  }
  unsigned num_tiles() const {
    return x_tiles() * y_tiles() * z_tiles();
  }
private:
  unsigned tiles_per_rank_[3];
  unsigned ranks_per_dim_[3];
  AddressSpace first_rank_;
};

//=============================================================================

class SoleilMapper : public DefaultMapper {
public:
  SoleilMapper(MapperRuntime* rt, Machine machine, Processor local)
    : DefaultMapper(rt, machine, local, "soleil_mapper"),
      all_procs_(remote_cpus.size()) {
    // Set the umask of the process to clear S_IWGRP and S_IWOTH.
    umask(022);
    // Assign ranks sequentially to samples, each sample getting one rank for
    // each super-tile.
    AddressSpace reqd_ranks = 0;
    unsigned num_samples = 0;
    auto process_config = [&](char* config_file) {
      Config config = parse_config(config_file);
      CHECK(config.Mapping.tiles[0] > 0 &&
            config.Mapping.tiles[1] > 0 &&
            config.Mapping.tiles[2] > 0 &&
            config.Mapping.tilesPerRank[0] > 0 &&
            config.Mapping.tilesPerRank[1] > 0 &&
            config.Mapping.tilesPerRank[2] > 0 &&
            config.Mapping.tiles[0] % config.Mapping.tilesPerRank[0] == 0 &&
            config.Mapping.tiles[1] % config.Mapping.tilesPerRank[1] == 0 &&
            config.Mapping.tiles[2] % config.Mapping.tilesPerRank[2] == 0,
            "Invalid tiling for sample %d", num_samples);
      sample_mappings_.emplace_back(config, reqd_ranks);
      unsigned num_ranks = sample_mappings_.back().num_ranks()
      if (reqd_ranks <= node_id && node_id < reqd_ranks + num_ranks) {
	local_sample_id_ = num_samples;
      }
      reqd_ranks += num_ranks;
      num_samples++;
    };
    // Locate all config files specified on the command-line arguments.
    InputArgs args = Runtime::get_input_args();
    for (int i = 0; i < args.argc; ++i) {
      if (strcmp(args.argv[i], "-i") == 0 && i < args.argc-1) {
        process_config(args.argv[i+1]);
      } else if (strcmp(args.argv[i], "-I") == 0 && i < args.argc-1) {
        std::ifstream csv_file(args.argv[i+1]);
        std::string json_filename;
        while (std::getline(csv_file, json_filename)) {
          process_config(&json_filename[0]);
        }
      }
    }
    // Verify that we have enough ranks.
    unsigned supplied_ranks = remote_cpus.size();
    CHECK(reqd_ranks <= supplied_ranks,
          "%d rank(s) required, but %d rank(s) supplied to Legion",
          reqd_ranks, supplied_ranks);
    if (reqd_ranks < supplied_ranks) {
      LOG.warning() << supplied_ranks << " rank(s) supplied to Legion,"
                    << " but only " << reqd_ranks << " required"
    }
    // Cache processor information.
    Machine::ProcessorQuery query(machine);
    for (auto it = query.begin(); it != query.end(); it++) {
      AddressSpace rank = it->address_space();
      Processor::Kind kind = it->kind();
      get_procs(rank, kind).push_back(*it);
    }
  }

public:
  virtual Processor default_policy_select_initial_processor(
                              MapperContext ctx,
                              const Task& task) {
    // DOM sweep & boundary tasks are individually launched; find the tile on
    // which they're centered and send them to the rank responsible for that.
    // TODO: Cache the decision.
    bool is_sweep = STARTS_WITH(task.get_task_name(), "sweep_");
    bool is_bound = STARTS_WITH(task.get_task_name(), "bound_");
    if (is_sweep || is_bound) {
      CHECK(task.parent_task != NULL &&
            strcmp(task.parent_task->get_task_name(), "work") == 0,
            "DOM tasks should only be launched from the work task directly.");
      // Retrieve sample information from parent work task.
      const Config* config = find_config(&task);
      assert(config != NULL);
      unsigned sample_id = config->Mapping.sampleId;
      assert(sample_id < sample_mappings_.size());
      const SampleMapping& mapping = sample_mappings_[sample_id];
      // Find the tile this task launch is centered on.
      DomainPoint tile = DomainPoint::nil();
      for (const RegionRequirement& req : task.regions) {
        assert(req.region.exists());
        const char* name = NULL;
        runtime->retrieve_name(ctx, req.region.get_field_space(), name);
        if ((is_sweep && strcmp(name, "Radiation_columns") == 0) ||
            (is_bound && strcmp(name, "face") == 0)) {
          tile = runtime->get_logical_region_color_point(ctx, req.region);
          break;
        }
      }
      CHECK(!tile.is_null(),
            "Cannot retrieve tile from DOM task launch -- did you change the"
            " names of the field spaces?");
      CHECK(tile.get_dim() == 3,
            "DOM task launches should only use the top-level tiling.");
      // A task that updates the far boundary on some direction is called with
      // a face tile one-over on that direction, so we have to subtract 1 to
      // find the actual tile the launch is centered on.
      if (is_bound) {
        std::regex regex("bound_([xyz])_(lo|hi)");
        std::cmatch match;
        CHECK(std::regex_match(task.get_task_name(), match, regex),
              "Cannot parse name of DOM boundary task -- did you change it?");
        if (match[2].str().compare("hi") == 0) {
          if (match[1].str().compare("x") == 0) { tile[0]--; }
          if (match[1].str().compare("y") == 0) { tile[1]--; }
          if (match[1].str().compare("z") == 0) { tile[2]--; }
        }
      }
      // Assign rank according to the precomputed mapping, then round-robin
      // over all the processors of the preffered kind within that rank.
      CHECK(0 <= tile[0] && tile[0] < config->Mapping.tiles[0] &&
            0 <= tile[1] && tile[1] < config->Mapping.tiles[1] &&
            0 <= tile[2] && tile[2] < config->Mapping.tiles[2],
            "DOM task launches should only use the top-level tiling.");
      AddressSpace target_rank = mapping.get_rank(tile[0], tile[1], tile[2]);
      VariantInfo info =
        default_find_preferred_variant(task, ctx, false/*needs tight*/);
      const std::vector<Processor>& procs =
        get_procs(target_rank, info.proc_kind);
      unsigned target_proc_id = mapping.get_proc_id(tile[0], tile[1], tile[2]);
      Processor target_proc = procs[target_proc_id % procs.size()];
      LOG.debug() << "Sample " << sample_id << ":"
                  << " Sequential launch:"
                  << " Task '" << task.get_task_name() << "'"
		  << " on tile " << tile
		  << " mapped to rank " << target_rank
                  << " processor " << target_proc;
      return target_proc;
    }
    // Send each work task to the first in the set of ranks allocated to the
    // corresponding sample.
    if (strcmp(task.get_task_name(), "work") == 0) {
      const Config* config = find_config(&task);
      assert(config != NULL);
      unsigned sample_id = config->Mapping.sampleId;
      assert(sample_id < sample_mappings_.size());
      AddressSpace target_rank = sample_mappings_[sample_id].get_rank(0,0,0);
      Processor target_proc = remote_cpus[target_rank];
      LOG.debug() << "Sample " << sample_id << ":"
                  << " Sequential launch:"
                  << " Task 'work'"
		  << " mapped to rank " << target_rank
                  << " processor " << target_proc;
      return target_proc;
    }
    // For other tasks, defer to the default mapping policy.
    return DefaultMapper::default_policy_select_initial_processor(ctx, task);
  }

  // Assign priorities to sweep tasks such that we prioritize the tile that has
  // more dependencies downstream.
  virtual TaskPriority default_policy_select_task_priority(
                              MapperContext ctx,
                              const Task& task) {
    if (!STARTS_WITH(task.get_task_name(), "sweep_")) {
      return DefaultMapper::default_policy_select_task_priority(ctx, task);
    }
    // Compute direction of sweep.
    int sweep_id = atoi(task.get_task_name() + sizeof("sweep_") - 1) - 1;
    CHECK(0 <= sweep_id && sweep_id <= 7,
          "Task %s: invalid sweep id", task.get_task_name());
    bool x_rev = (sweep_id >> 0) & 1;
    bool y_rev = (sweep_id >> 1) & 1;
    bool z_rev = (sweep_id >> 2) & 1;
    // Find the tile this task launch is centered on.
    DomainPoint tile = DomainPoint::nil();
    for (const RegionRequirement& req : task.regions) {
      assert(req.region.exists());
      const char* name = NULL;
      runtime->retrieve_name(ctx, req.region.get_field_space(), name);
      if (strcmp(name, "Radiation_columns") == 0) {
        tile = runtime->get_logical_region_color_point(ctx, req.region);
        break;
      }
    }
    CHECK(!tile.is_null(),
          "Cannot retrieve tile from DOM task launch -- did you change the"
          " names of the field spaces?");
    const SampleMapping& local_mapping = sample_mappings_[local_sample_id_];
    CHECK(tile.get_dim() == 3 &&
          0 <= tile[0] && tile[0] < local_mapping.x_tiles() &&
          0 <= tile[1] && tile[1] < local_mapping.y_tiles() &&
          0 <= tile[2] && tile[2] < local_mapping.z_tiles(),
          "DOM task launches should only use the top-level tiling.");
    // Assign priority according to the number of diagonals between this launch
    // and the end of the domain.
    int priority =
      (x_rev ? tile[0] : local_mapping.x_tiles() - tile[0] - 1) +
      (y_rev ? tile[1] : local_mapping.y_tiles() - tile[1] - 1) +
      (z_rev ? tile[2] : local_mapping.z_tiles() - tile[2] - 1);
    LOG.debug() << "Sample " << local_sample_id_ << ":"
                << " Task '" << task.get_task_name() << "'"
                << " on tile " << tile
                << " given priority " << priority;
    return priority;
  }

  // TODO: Select appropriate memories for instances that will be communicated,
  // (e.g. parallelizer-created ghost partitions), such as RDMA memory,
  // zero-copy memory.
  virtual Memory default_policy_select_target_memory(
                              MapperContext ctx,
                              Processor target_proc,
                              const RegionRequirement& req) {
    return DefaultMapper::default_policy_select_target_memory
      (ctx, target_proc, req);
  }

  // Disable an optimization done by the default mapper (attempts to reuse an
  // instance that covers a superset of the requested index space, by searching
  // higher up the partition tree).
  virtual LogicalRegion default_policy_select_instance_region(
                              MapperContext ctx,
                              Memory target_memory,
                              const RegionRequirement& req,
                              const LayoutConstraintSet& constraints,
                              bool force_new_instances,
                              bool meets_constraints) {
    return req.region;
  }

  // Farm index-space launches made by work tasks across all the ranks
  // allocated to the corresponding sample.
  // TODO: Cache the decision.
  virtual void slice_task(const MapperContext ctx,
                          const Task& task,
                          const SliceTaskInput& input,
                          SliceTaskOutput& output) {
    CHECK(task.parent_task != NULL &&
          strcmp(task.parent_task->get_task_name(), "work") == 0,
          "Index-space launches only allowed in work task");
    output.verify_correctness = false;
    // Retrieve sample information from parent task.
    const Config* config = find_config(&task);
    assert(config != NULL);
    unsigned sample_id = config->Mapping.sampleId;
    assert(sample_id < sample_mappings_.size());
    const SampleMapping& mapping = sample_mappings_[sample_id];
    CHECK(input.domain.get_dim() == 3 &&
          input.domain.lo()[0] == 0 &&
          input.domain.lo()[1] == 0 &&
          input.domain.lo()[2] == 0 &&
          input.domain.hi()[0] == config->Mapping.tiles[0] - 1 &&
          input.domain.hi()[1] == config->Mapping.tiles[1] - 1 &&
          input.domain.hi()[2] == config->Mapping.tiles[2] - 1,
          "Index-space launches in the work task should only use the"
          " top-level tiling.");
    // Allocate tasks among all the processors of the same kind as the original
    // target, on each rank allocated to this sample.
    for (unsigned x = 0; x < config->Mapping.tiles[0]; ++x) {
      for (unsigned y = 0; y < config->Mapping.tiles[1]; ++y) {
        for (unsigned z = 0; z < config->Mapping.tiles[2]; ++z) {
          AddressSpace target_rank = mapping.get_rank(x, y, z);
          const std::vector<Processor>& procs =
            get_procs(target_rank, task.target_proc.kind());
          unsigned target_proc_id = mapping.get_proc_id(x, y, z);
          Processor target_proc = procs[target_proc_id % procs.size()];
          output.slices.emplace_back(Rect<3>(Point<3>(x,y,z), Point<3>(x,y,z)),
                                     target_proc,
                                     false /*recurse*/,
                                     false /*stealable*/);
          LOG.debug() << "Sample " << sample_id << ":"
                      << " Index-space launch:"
                      << " Task '" << task.get_task_name() << "'"
		      << " on tile (" << x << "," << y << "," << z << ")"
		      << " mapped to rank " << target_rank
		      << " processor " << target_proc;
        }
      }
    }
  }

private:
  std::vector<Processor>& get_procs(AddressSpace rank, Processor::Kind kind) {
    assert(rank < all_procs_.size());
    auto& rank_procs = all_procs_[rank];
    if (kind >= rank_procs.size()) {
      rank_procs.resize(kind + 1);
    }
    return rank_procs[kind];
  }

private:
  std::vector<SampleMapping> sample_mappings_;
  std::vector<std::vector<std::vector<Processor> > > all_procs_;
  unsigned local_sample_id_;
};

//=============================================================================

static void create_mappers(Machine machine,
                           HighLevelRuntime* runtime,
                           const std::set<Processor>& local_procs) {
  for (Processor proc : local_procs) {
    SoleilMapper* mapper =
      new SoleilMapper(runtime->get_mapper_runtime(), machine, proc);
    runtime->replace_default_mapper(mapper, proc);
  }
}

void register_mappers() {
  Runtime::add_registration_callback(create_mappers);
}
