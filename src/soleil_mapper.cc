#include <array>
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

#define EQUALS(s1, s2) (strcmp((s1), (s2)) == 0)

#define STARTS_WITH(str, prefix)                \
  (strncmp((str), (prefix), sizeof(prefix) - 1) == 0)

static const void* first_arg(const Task& task) {
  const char* ptr = static_cast<const char*>(task.args);
  // Skip over Regent-added arguments.
  // XXX: This assumes Regent's calling convention won't change.
  return static_cast<const void*>(ptr + sizeof(uint64_t));
}

//=============================================================================

// Maps super-tiles to ranks, in row-major order. Within a super-tile, tiles
// are assigned unique processor IDs, in row-major order. The mapper will have
// to match those IDs to real processors (not necessarily one per ID).
class SampleMapping {
public:
  SampleMapping(const Config& config, AddressSpace first_rank)
    : tiles_per_rank_{static_cast<unsigned>(config.Mapping.tilesPerRank[0]),
                      static_cast<unsigned>(config.Mapping.tilesPerRank[1]),
                      static_cast<unsigned>(config.Mapping.tilesPerRank[2])},
      ranks_per_dim_{static_cast<unsigned>(config.Mapping.tiles[0]
                                           / config.Mapping.tilesPerRank[0]),
                     static_cast<unsigned>(config.Mapping.tiles[1]
                                           / config.Mapping.tilesPerRank[1]),
                     static_cast<unsigned>(config.Mapping.tiles[2]
                                           / config.Mapping.tilesPerRank[2])},
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
    auto process_config = [&](const Config& config) {
      CHECK(config.Mapping.tiles[0] > 0 &&
            config.Mapping.tiles[1] > 0 &&
            config.Mapping.tiles[2] > 0 &&
            config.Mapping.tilesPerRank[0] > 0 &&
            config.Mapping.tilesPerRank[1] > 0 &&
            config.Mapping.tilesPerRank[2] > 0 &&
            config.Mapping.tiles[0] % config.Mapping.tilesPerRank[0] == 0 &&
            config.Mapping.tiles[1] % config.Mapping.tilesPerRank[1] == 0 &&
            config.Mapping.tiles[2] % config.Mapping.tilesPerRank[2] == 0,
            "Invalid tiling for sample %lu", sample_mappings_.size() + 1);
      sample_mappings_.emplace_back(config, reqd_ranks);
      reqd_ranks += sample_mappings_.back().num_ranks();
    };
    // Locate all config files specified on the command-line arguments.
    InputArgs args = Runtime::get_input_args();
    for (int i = 0; i < args.argc; ++i) {
      if (EQUALS(args.argv[i], "-i") && i < args.argc-1) {
        Config config;
        parse_Config(&config, args.argv[i+1]);
        process_config(config);
      } else if (EQUALS(args.argv[i], "-m") && i < args.argc-1) {
        MultiConfig mc;
        parse_MultiConfig(&mc, args.argv[i+1]);
        process_config(mc.configs[0]);
        process_config(mc.configs[1]);
      }
    }
    // Verify that we have enough ranks.
    unsigned supplied_ranks = remote_cpus.size();
    CHECK(reqd_ranks <= supplied_ranks,
          "%d rank(s) required, but %d rank(s) supplied to Legion",
          reqd_ranks, supplied_ranks);
    if (reqd_ranks < supplied_ranks) {
      LOG.warning() << supplied_ranks << " rank(s) supplied to Legion,"
                    << " but only " << reqd_ranks << " required";
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
    // For tasks that are individually launched, find the tile on which they're
    // centered and send them to the rank responsible for that.
    // TODO: Cache the decision.
    if (STARTS_WITH(task.get_task_name(), "sweep_") ||
        STARTS_WITH(task.get_task_name(), "TradeQueue_pull")) {
      unsigned sample_id = find_sample_id(ctx, task);
      const SampleMapping& mapping = sample_mappings_[sample_id];
      DomainPoint tile = find_tile(ctx, task, mapping);
      VariantInfo info =
        default_find_preferred_variant(task, ctx, false/*needs tight*/);
      Processor target_proc = select_proc(tile, info.proc_kind, mapping);
      LOG.debug() << "Sample " << sample_id << ":"
                  << " Sequential launch:"
                  << " Task " << task.get_task_name()
                  << " on tile " << tile
                  << " mapped to processor " << target_proc;
      return target_proc;
    }
    // Send each work task to the first in the set of ranks allocated to the
    // corresponding sample.
    else if (STARTS_WITH(task.get_task_name(), "work")) {
      unsigned sample_id = static_cast<unsigned>(-1);
      if (EQUALS(task.get_task_name(), "workSingle")) {
        const Config* config = static_cast<const Config*>(first_arg(task));
        sample_id = static_cast<unsigned>(config->Mapping.sampleId);
        assert(sample_id < sample_mappings_.size());
      } else if (EQUALS(task.get_task_name(), "workDual")) {
        const MultiConfig* mc =
          static_cast<const MultiConfig*>(first_arg(task));
        sample_id = static_cast<unsigned>(mc->configs[0].Mapping.sampleId);
        assert(sample_id < sample_mappings_.size());
      } else {
        CHECK(false, "Unexpected work task name: %s", task.get_task_name());
      }
      const SampleMapping& mapping = sample_mappings_[sample_id];
      Processor target_proc =
        select_proc(Point<3>(0,0,0), Processor::LOC_PROC, mapping);
      LOG.debug() << "Sample " << sample_id << ":"
                  << " Sequential launch:"
                  << " Task " << task.get_task_name()
                  << " mapped to processor " << target_proc;
      return target_proc;
    }
    // For index space tasks, defer to the default mapping policy, and
    // slice_task will decide the final mapping.
    else if (task.is_index_space) {
      return DefaultMapper::default_policy_select_initial_processor(ctx, task);
    }
    // For certain whitelisted tasks, defer to the default mapping policy.
    else if (EQUALS(task.get_task_name(), "main") ||
             STARTS_WITH(task.get_task_name(), "Console_Write") ||
             STARTS_WITH(task.get_task_name(), "Probe_Write") ||
             EQUALS(task.get_task_name(), "createDir") ||
             EQUALS(task.get_task_name(), "cache_grid_translation") ||
             EQUALS(task.get_task_name(), "initialize_angles") ||
             EQUALS(task.get_task_name(), "__dummy") ||
             STARTS_WITH(task.get_task_name(), "__binary_")) {
      return DefaultMapper::default_policy_select_initial_processor(ctx, task);
    }
    // For other tasks, fail & notify the user.
    else {
      CHECK(false, "Unhandled non-index space task %s", task.get_task_name());
      return Processor::NO_PROC;
    }
  }

  virtual TaskPriority default_policy_select_task_priority(
                              MapperContext ctx,
                              const Task& task) {
    // Unless handled specially below, all tasks have the same priority.
    int priority = 0;
    // Assign priorities to sweep tasks such that we prioritize the tile that
    // has more dependencies downstream (count the number of diagonals between
    // the launch tile and the end of the domain).
    if (STARTS_WITH(task.get_task_name(), "sweep_")) {
      unsigned sample_id = find_sample_id(ctx, task);
      const SampleMapping& mapping = sample_mappings_[sample_id];
      std::array<bool,3> dir = parse_direction(task);
      DomainPoint tile = find_tile(ctx, task, mapping);
      priority =
	(dir[0] ? mapping.x_tiles() - tile[0] - 1 : tile[0]) +
	(dir[1] ? mapping.y_tiles() - tile[1] - 1 : tile[1]) +
	(dir[2] ? mapping.z_tiles() - tile[2] - 1 : tile[2]) ;
      LOG.debug() << "Sample " << sample_id << ":"
                  << " Task " << task.get_task_name()
                  << " on tile " << tile
                  << " given priority " << priority;
    }
    // Increase priority of tasks on the critical path of the fluid solve.
    if (STARTS_WITH(task.get_task_name(), "Flow_ComputeVelocityGradient") ||
        STARTS_WITH(task.get_task_name(), "Flow_UpdateGhostVelocityGradient") ||
	STARTS_WITH(task.get_task_name(), "Flow_GetFlux") ||
	STARTS_WITH(task.get_task_name(), "Flow_UpdateUsingFlux")) {
      unsigned sample_id = find_sample_id(ctx, task);
      priority = 1;
      LOG.debug() << "Sample " << sample_id << ":"
                  << " Task " << task.get_task_name()
                  << " given priority " << priority;
    }
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

  // Farm index space launches made by work tasks across all the ranks
  // allocated to the corresponding sample.
  // TODO: Cache the decision.
  virtual void slice_task(const MapperContext ctx,
                          const Task& task,
                          const SliceTaskInput& input,
                          SliceTaskOutput& output) {
    output.verify_correctness = false;
    unsigned sample_id = find_sample_id(ctx, task);
    const SampleMapping& mapping = sample_mappings_[sample_id];
    Domain domain = input.domain;
    // Certain tasks are launched on a 2D domain. Extend each domain point to a
    // 3D tile, by filling in the missing dimension, to decide how to map them.
    unsigned dim = unsigned(-1);
    bool dir = false;
    if (domain.get_dim() == 2) {
      dim = parse_dimension(task);
      dir = parse_direction(task)[dim];
      if (STARTS_WITH(task.get_task_name(), "initialize_faces_") ||
          STARTS_WITH(task.get_task_name(), "bound_")) {
        // Do nothing
      } else if (STARTS_WITH(task.get_task_name(), "cache_intensity_")) {
        // We want to run these tasks on the opposite end of the domain implied
        // by their name.
        dir = !dir;
      } else {
        CHECK(false, "Unexpected 2D domain on index space launch of task %s",
              task.get_task_name());
      }
    } else {
      CHECK(domain.get_dim() == 3 &&
            0 <= domain.lo()[0] && domain.hi()[0] < mapping.x_tiles() &&
            0 <= domain.lo()[1] && domain.hi()[1] < mapping.y_tiles() &&
            0 <= domain.lo()[2] && domain.hi()[2] < mapping.z_tiles(),
            "Unexpected 3D domain on index space launch of task %s",
            task.get_task_name());
    }
    // Allocate tasks among all the processors of the same kind as the original
    // target, on each rank allocated to this sample.
    for (Domain::DomainPointIterator it(domain); it; it++) {
      DomainPoint tile = it.p;
      if (domain.get_dim() == 2) {
        unsigned coord =
          (dim == 0) ? (dir ? 0 : mapping.x_tiles()-1) :
          (dim == 1) ? (dir ? 0 : mapping.y_tiles()-1) :
         /*dim == 2*/  (dir ? 0 : mapping.z_tiles()-1) ;
        tile =
          (dim == 0) ? Point<3>(coord, it.p[0], it.p[1]) :
          (dim == 1) ? Point<3>(it.p[0], coord, it.p[1]) :
         /*dim == 2*/  Point<3>(it.p[0], it.p[1], coord) ;
      }
      Processor target_proc =
        select_proc(tile, task.target_proc.kind(), mapping);
      output.slices.emplace_back(Domain(it.p, it.p), target_proc,
                                 false/*recurse*/, false/*stealable*/);
      LOG.debug() << "Sample " << sample_id << ":"
                  << " Index space launch:"
                  << " Task " << task.get_task_name()
                  << " on domain point " << it.p
                  << " tile " << tile
                  << " mapped to processor " << target_proc;
    }
  }

  virtual void map_copy(const MapperContext ctx,
                        const Copy& copy,
                        const MapCopyInput& input,
                        MapCopyOutput& output) {
    // For HDF copies, defer to the default mapping policy.
    if (EQUALS(copy.parent_task->get_task_name(), "dumpTile") ||
        EQUALS(copy.parent_task->get_task_name(), "loadTile")) {
      DefaultMapper::map_copy(ctx, copy, input, output);
      return;
    }
    // Sanity checks
    // TODO: Check that this is on the fluid grid.
    CHECK(copy.src_indirect_requirements.empty() &&
          copy.dst_indirect_requirements.empty() &&
          !copy.is_index_space &&
          copy.src_requirements.size() == 1 &&
          copy.dst_requirements.size() == 1 &&
          copy.src_requirements[0].region.exists() &&
          copy.dst_requirements[0].region.exists() &&
          !copy.dst_requirements[0].is_restricted() &&
          copy.src_requirements[0].privilege_fields.size() == 1 &&
          copy.dst_requirements[0].privilege_fields.size() == 1 &&
          input.src_instances[0].empty() &&
          // NOTE: The runtime should be passing the existing fluid instances
          // on the destination nodes as usable destinations, but doesn't, so
          // we have to perform an explicit runtime call. If this behavior ever
          // changes, this check will make sure we find out.
          input.dst_instances[0].empty(),
          "Unexpected arguments on explicit copy");
    // Retrieve copy details.
    // We map according to the destination of the copy. We expand the
    // destination domain to the full tile, to make sure we reuse the existing
    // instances.
    const RegionRequirement& src_req = copy.src_requirements[0];
    const RegionRequirement& dst_req = copy.dst_requirements[0];
    unsigned sample_id = find_sample_id(ctx, dst_req);
    const SampleMapping& mapping = sample_mappings_[sample_id];
    LogicalRegion src_region = src_req.region;
    LogicalRegion dst_region = dst_req.region;
    CHECK(runtime->get_index_space_depth
            (ctx, src_region.get_index_space()) == 2 &&
          runtime->get_index_space_depth
            (ctx, dst_region.get_index_space()) == 4,
          "Unexpected bounds on explicit copy");
    dst_region =
      runtime->get_parent_logical_region(ctx,
        runtime->get_parent_logical_partition(ctx, dst_region));
    DomainPoint src_tile =
      runtime->get_logical_region_color_point(ctx, src_region);
    DomainPoint dst_tile =
      runtime->get_logical_region_color_point(ctx, dst_region);
    CHECK(src_tile.get_dim() == 3 &&
          dst_tile.get_dim() == 3 &&
          src_tile[0] == dst_tile[0] &&
          src_tile[1] == dst_tile[1] &&
          src_tile[2] == dst_tile[2] &&
          0 <= dst_tile[0] && dst_tile[0] < mapping.x_tiles() &&
          0 <= dst_tile[1] && dst_tile[1] < mapping.y_tiles() &&
          0 <= dst_tile[2] && dst_tile[2] < mapping.z_tiles(),
          "Unexpected bounds on explicit copy");
    // Always use a virtual instance for the source.
    output.src_instances[0].clear();
    output.src_instances[0].push_back
      (PhysicalInstance::get_virtual_instance());
    // Write the data directly on the best memory for the task that will be
    // using it.
    // XXX: We assume that if we have GPUs, then the GPU variants will be used.
    output.dst_instances[0].clear();
    output.dst_instances[0].emplace_back();
    Processor::Kind proc_kind =
      (local_gpus.size() > 0) ? Processor::TOC_PROC :
      (local_omps.size() > 0) ? Processor::OMP_PROC :
                                Processor::LOC_PROC ;
    Processor target_proc = select_proc(dst_tile, proc_kind, mapping);
    Memory target_memory =
      default_policy_select_target_memory(ctx, target_proc, dst_req);
    LayoutConstraintSet dst_constraints;
    dst_constraints.add_constraint
      (FieldConstraint(dst_req.privilege_fields,
                       false/*contiguous*/, false/*inorder*/));
    CHECK(runtime->find_physical_instance
            (ctx, target_memory, dst_constraints,
             std::vector<LogicalRegion>{dst_region},
             output.dst_instances[0][0],
             true/*acquire*/, false/*tight_region_bounds*/),
          "Could not locate destination instance for explicit copy");
  }

private:
  unsigned find_sample_id(const MapperContext ctx,
                          const RegionRequirement& req) const {
    LogicalRegion region = req.region.exists() ? req.region
      : runtime->get_parent_logical_region(ctx, req.partition);
    region = get_root(ctx, region);
    const void* info = NULL;
    size_t info_size = 0;
    bool success = runtime->retrieve_semantic_information
      (ctx, region, SAMPLE_ID_TAG, info, info_size,
       false/*can_fail*/, true/*wait_until_ready*/);
    CHECK(success, "Missing SAMPLE_ID_TAG semantic information on region");
    assert(info_size == sizeof(unsigned));
    unsigned sample_id = *static_cast<const unsigned*>(info);
    assert(sample_id < sample_mappings_.size());
    return sample_id;
  }

  unsigned find_sample_id(const MapperContext ctx,
                          const Task& task) const {
    CHECK(!task.regions.empty(),
          "No region argument on launch of task %s", task.get_task_name());
    return find_sample_id(ctx, task.regions[0]);
  }

  // XXX: Always using the first region argument to figure out the tile.
  DomainPoint find_tile(const MapperContext ctx,
                        const Task& task,
                        const SampleMapping& mapping) const {
    CHECK(!task.regions.empty(),
          "No region argument on launch of task %s", task.get_task_name());
    const RegionRequirement& req = task.regions[0];
    assert(req.region.exists());
    DomainPoint tile =
      runtime->get_logical_region_color_point(ctx, req.region);
    CHECK(tile.get_dim() == 3 &&
          0 <= tile[0] && tile[0] < mapping.x_tiles() &&
          0 <= tile[1] && tile[1] < mapping.y_tiles() &&
          0 <= tile[2] && tile[2] < mapping.z_tiles(),
          "Launch of task %s using incorrect tiling", task.get_task_name());
    return tile;
  }

  std::array<bool,3> parse_direction(const Task& task) const {
    std::array<bool,3> dir = {{true, true, true}};
    std::regex regex("\\w*_([1-8]|lo|hi)");
    std::cmatch match;
    CHECK(std::regex_match(task.get_task_name(), match, regex),
          "Cannot parse quadrant info from task name: %s",
          task.get_task_name());
    if (match[1].str().compare("lo") == 0) {
      // Do nothing
    } else if (match[1].str().compare("hi") == 0) {
      dir[parse_dimension(task)] = false;
    } else {
      unsigned quadrant = std::stoul(match[1].str()) - 1;
      dir[0] = 1 - ((quadrant >> 0) & 1);
      dir[1] = 1 - ((quadrant >> 1) & 1);
      dir[2] = 1 - ((quadrant >> 2) & 1);
    }
    return dir;
  }

  unsigned parse_dimension(const Task& task) const {
    std::regex regex("\\w*_([xyz])_([1-8]|lo|hi)");
    std::cmatch match;
    CHECK(std::regex_match(task.get_task_name(), match, regex),
          "Cannot parse dimension from task name: %s",
          task.get_task_name());
    return
      (match[1].str().compare("x") == 0) ? 0 :
      (match[1].str().compare("y") == 0) ? 1 :
     /*match[1].str().compare("z") == 0)*/ 2 ;
  }

  // Assign rank according to the precomputed mapping, then round-robin over
  // all the processors of the desired kind within that rank.
  // NOTE: This function doesn't sanity check its input.
  Processor select_proc(const DomainPoint& tile,
                        Processor::Kind kind,
                        const SampleMapping& mapping) {
    AddressSpace rank = mapping.get_rank(tile[0], tile[1], tile[2]);
    const std::vector<Processor>& procs = get_procs(rank, kind);
    unsigned proc_id = mapping.get_proc_id(tile[0], tile[1], tile[2]);
    return procs[proc_id % procs.size()];
  }

  std::vector<Processor>& get_procs(AddressSpace rank, Processor::Kind kind) {
    assert(rank < all_procs_.size());
    auto& rank_procs = all_procs_[rank];
    if (kind >= rank_procs.size()) {
      rank_procs.resize(kind + 1);
    }
    return rank_procs[kind];
  }

  LogicalRegion get_root(const MapperContext ctx, LogicalRegion region) const {
    while (runtime->has_parent_logical_partition(ctx, region)) {
      region =
        runtime->get_parent_logical_region(ctx,
          runtime->get_parent_logical_partition(ctx, region));
    }
    return region;
  }

private:
  std::vector<SampleMapping> sample_mappings_;
  std::vector<std::vector<std::vector<Processor> > > all_procs_;
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
