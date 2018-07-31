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
      if (strcmp(args.argv[i], "-i") == 0 && i < args.argc-1) {
        Config config;
        parse_Config(&config, args.argv[i+1]);
        process_config(config);
      } else if (strcmp(args.argv[i], "-m") == 0 && i < args.argc-1) {
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
    // DOM sweep & boundary tasks are individually launched; find the tile on
    // which they're centered and send them to the rank responsible for that.
    // TODO: Cache the decision.
    bool is_sweep = STARTS_WITH(task.get_task_name(), "sweep_");
    bool is_bound = STARTS_WITH(task.get_task_name(), "bound_");
    if (is_sweep || is_bound) {
      CHECK(task.parent_task != NULL &&
            STARTS_WITH(task.parent_task->get_task_name(), "work"),
            "DOM tasks should only be launched from the work task directly.");
      // Retrieve sample information.
      unsigned sample_id = find_sample_id(ctx, task);
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
      CHECK(0 <= tile[0] && tile[0] < mapping.x_tiles() &&
            0 <= tile[1] && tile[1] < mapping.y_tiles() &&
            0 <= tile[2] && tile[2] < mapping.z_tiles(),
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
                  << " Task " << task.get_task_name()
                  << " on tile " << tile
                  << " mapped to rank " << target_rank
                  << " processor " << target_proc;
      return target_proc;
    }
    // Send each work task to the first in the set of ranks allocated to the
    // corresponding sample.
    if (STARTS_WITH(task.get_task_name(), "work")) {
      unsigned sample_id = static_cast<unsigned>(-1);
      if (strcmp(task.get_task_name(), "workSingle") == 0) {
        const Config* config = static_cast<const Config*>(first_arg(task));
        sample_id = static_cast<unsigned>(config->Mapping.sampleId);
        assert(sample_id < sample_mappings_.size());
      } else if (strcmp(task.get_task_name(), "workDual") == 0) {
        const MultiConfig* mc =
          static_cast<const MultiConfig*>(first_arg(task));
        sample_id = static_cast<unsigned>(mc->configs[0].Mapping.sampleId);
        assert(sample_id < sample_mappings_.size());
      } else {
        CHECK(false, "Unexpected work task name: %s", task.get_task_name());
      }

      AddressSpace target_rank = sample_mappings_[sample_id].get_rank(0,0,0);
      Processor target_proc = remote_cpus[target_rank];
      LOG.debug() << "Sample " << sample_id << ":"
                  << " Sequential launch:"
                  << " Task work"
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
    // Retrieve sample information.
    unsigned sample_id = find_sample_id(ctx, task);
    const SampleMapping& mapping = sample_mappings_[sample_id];
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
    CHECK(tile.get_dim() == 3 &&
          0 <= tile[0] && tile[0] < mapping.x_tiles() &&
          0 <= tile[1] && tile[1] < mapping.y_tiles() &&
          0 <= tile[2] && tile[2] < mapping.z_tiles(),
          "DOM task launches should only use the top-level tiling.");
    // Assign priority according to the number of diagonals between this launch
    // and the end of the domain.
    int priority =
      (x_rev ? tile[0] : mapping.x_tiles() - tile[0] - 1) +
      (y_rev ? tile[1] : mapping.y_tiles() - tile[1] - 1) +
      (z_rev ? tile[2] : mapping.z_tiles() - tile[2] - 1);
    LOG.debug() << "Sample " << sample_id << ":"
                << " Task " << task.get_task_name()
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
          STARTS_WITH(task.parent_task->get_task_name(), "work"),
          "Index-space launches only allowed in work task");
    output.verify_correctness = false;
    // Retrieve sample information.
    unsigned sample_id = find_sample_id(ctx, task);
    const SampleMapping& mapping = sample_mappings_[sample_id];
    CHECK(input.domain.get_dim() == 3 &&
          input.domain.lo()[0] == 0 &&
          input.domain.lo()[1] == 0 &&
          input.domain.lo()[2] == 0 &&
          input.domain.hi()[0] == mapping.x_tiles() - 1 &&
          input.domain.hi()[1] == mapping.y_tiles() - 1 &&
          input.domain.hi()[2] == mapping.z_tiles() - 1,
          "Index-space launches in the work task should only use the"
          " top-level tiling.");
    // Allocate tasks among all the processors of the same kind as the original
    // target, on each rank allocated to this sample.
    for (unsigned x = 0; x < mapping.x_tiles(); ++x) {
      for (unsigned y = 0; y < mapping.y_tiles(); ++y) {
        for (unsigned z = 0; z < mapping.z_tiles(); ++z) {
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
                      << " Task " << task.get_task_name()
                      << " on tile (" << x << "," << y << "," << z << ")"
                      << " mapped to rank " << target_rank
                      << " processor " << target_proc;
        }
      }
    }
  }

  void print_instance_info(const MapperContext ctx,
                           LogicalRegion region,
                           const PhysicalInstance& inst) {
    std::cout << "[" << node_id << "] ";
    std::cout << "    instance:";
    std::cout << " exists " << inst.exists();
    std::cout << " normal " << inst.is_normal_instance();
    std::cout << " virtual " << inst.is_virtual_instance();
    std::cout << " reduction " << inst.is_reduction_instance();
    std::cout << " external " << inst.is_external_instance();
    if (inst.is_normal_instance()) {
      std::cout << " memory " << inst.get_location();
      std::cout << " domain " << runtime->get_index_space_domain(ctx, region.get_index_space());
    }
    std::set<FieldID> fields;
    inst.get_fields(fields);
    std::cout << " fields ";
    for (FieldID fid : fields) {
      const char* field_name = NULL;
      runtime->retrieve_name(ctx, region.get_field_space(), fid, field_name);
      std::cout << " " << field_name;
    }
    std::cout << std::endl;
  }

  void print_region_req_info(const MapperContext ctx,
                             const RegionRequirement& req) {
    assert(req.region.exists());
    std::cout << "[" << node_id << "] ";
    std::cout << "  region req:";
    std::cout << " region " << req.region;
    const char* fs_name = NULL;
    runtime->retrieve_name(ctx, req.region.get_field_space(), fs_name);
    std::cout << " fieldspace " << fs_name;
    std::cout << " tile " << runtime->get_logical_region_color_point(ctx, req.region);
    std::cout << " domain " << runtime->get_index_space_domain(ctx, req.region.get_index_space());
    std::cout << " fields";
      for (FieldID fid : req.privilege_fields) {
        const char* field_name = NULL;
        runtime->retrieve_name(ctx, req.region.get_field_space(), fid, field_name);
        std::cout << " " << field_name;
      }
    std::cout << std::endl;
  }

  void print_op_info(const MapperContext ctx,
                     const std::vector<RegionRequirement>& requirements,
                     const std::vector<std::vector<PhysicalInstance> >& instances) {
    for (size_t idx = 0; idx < requirements.size(); ++idx) {
      const RegionRequirement& req = requirements[idx];
      assert(req.region.exists());
      print_region_req_info(ctx, req);
      std::cout << "[" << node_id << "] " << "  usable instances:" << std::endl;
      for (const PhysicalInstance& inst : instances[idx]) {
        print_instance_info(ctx, req.region, inst);
      }
      std::cout << "[" << node_id << "] " << "  all instances:" << std::endl;
      for (Memory mem : Machine::MemoryQuery(machine)) {
        LayoutConstraintSet constraints;
        constraints.add_constraint
          (FieldConstraint(req.privilege_fields,
                           false/*contiguous*/,
                           false/*inorder*/));
        std::vector<LogicalRegion> regions = {req.region};
        PhysicalInstance inst;
        if (runtime->find_physical_instance(ctx, mem, constraints, regions, inst,
                                            false /*acquire*/,
                                            false /*tight_region_bounds*/)) {
          print_instance_info(ctx, req.region, inst);
        }
      }
    }
  }

  virtual void map_copy(const MapperContext ctx,
                        const Copy& copy,
                        const MapCopyInput& input,
                        MapCopyOutput& output) {
    // // Sanity checks
    // // TODO: Check that this is on the fluid grid.
    // CHECK(STARTS_WITH(copy.parent_task->get_task_name(), "work"),
    //       "Explicit copies only allowed in work task");
    // CHECK(copy.src_indirect_requirements.empty() &&
    //       copy.dst_indirect_requirements.empty() &&
    //       !copy.is_index_space &&
    //       copy.src_requirements.size() == 1 &&
    //       copy.dst_requirements.size() == 1 &&
    //       copy.src_requirements[0].privilege == READ_PRIV &&
    //       copy.dst_requirements[0].privilege == WRITE_PRIV &&
    //       copy.src_requirements[0].region.exists() &&
    //       copy.dst_requirements[0].region.exists() &&
    //       !copy.dst_requirements[0].is_restricted() &&
    //       copy.src_requirements[0].privilege_fields.size() == 1 &&
    //       copy.dst_requirements[0].privilege_fields.size() == 1 &&
    //       input.src_instances[0].empty() &&
    //       input.dst_instances[0].size() <= 1,
    //       "Unexpected arguments on explicit copy");
    // // Retrieve copy details
    // const RegionRequirement& dst_req = copy.dst_requirements[0];
    // FieldID dst_fld = *(dst_req.privilege_fields.begin());
    // DomainPoint dst_tile =
    //   runtime->get_logical_region_color_point(ctx, dst_req.region);
    // // Always use a virtual instance for the source.
    // output.src_instances[0].clear();
    // output.src_instances[0].push_back
    //   (PhysicalInstance::get_virtual_instance());
    // // Place the destination instance directly on the best memory for the task
    // // that will be using this data, on the remote node.
    // output.dst_instances[0].clear();

    // AddressSpace target_rank =
    //   mapping.get_rank(dst_tile[0], dst_tile[1], dst_tile[2]);
    // unsigned target_proc_id =
    //   mapping.get_proc_id(dst_tile[0], dst_tile[1], dst_tile[2]);
    // // We assume that, if we have GPUs, then the GPU variant will be selected.
    // Processor::Kind proc_kind =
    //   (local_gpus.size() > 0) ? Processor::TOC_PROC :
    //   (local_omps.size() > 0) ? Processor::OMP_PROC :
    //   Processor::LOC_PROC);

    // const std::vector<Processor>& procs = get_procs(target_rank, proc_kind);
    // Processor target_proc = procs[target_proc_id % procs.size()];
    // Memory target_memory =
    //   default_policy_select_target_memory(ctx, target_proc, dst_req);
    // // If we have mapped this copy in the past, we should be able to reuse the
    // // previously created instance.
    // if (!input.dst_instances[0].empty()) {
    //   bool acquired = runtime->acquire_instances(ctx, input.dst_instances[0]);
    //   assert(acquired);
    //   const PhysicalInstance& inst = *(input.dst_instances[0].begin());
    //   assert(inst.get_location() == target_memory &&
    //          inst.has_field(dst_fld) &&
    //          inst.is_normal_instance());
    //   output.dst_instances[0].push_back(inst);
    //   return;
    // }
    // // This is the first time we're mapping this task; create a new instance to
    // // house the data on the destination node.
    // LayoutConstraintSet constraints;
    // default_policy_select_constraints
    //   (ctx, constraints, target_memory, dst_req);
    // constraints.add_constraint
    //   (FieldConstraint(std::vector<FieldID>{dst_fld},
    //                    false/*contiguous*/, false/*inorder*/));
    // output.dst_instances[0].emplace_back();
    // CHECK(default_make_instance(ctx, target_memory, constraints,
    //                             output.dst_instances[0].back(), COPY_MAPPING,
    //                             false/*force_new*/, true/*meets*/, dst_req),
    //       "Failed to create remote instance for copy");

    std::cout << "[" << node_id << "] " << "COPY IN TASK " << copy.parent_task->get_task_name() << std::endl;
    std::cout << "[" << node_id << "] " << "INPUT SOURCES:" << std::endl;
    print_op_info(ctx, copy.src_requirements, input.src_instances);
    std::cout << "[" << node_id << "] " << "INPUT DESTINATIONS:" << std::endl;
    print_op_info(ctx, copy.dst_requirements, input.dst_instances);
    DefaultMapper::map_copy(ctx, copy, input, output);
    std::cout << "[" << node_id << "] " << "OUTPUT SOURCES:" << std::endl;
    print_op_info(ctx, copy.src_requirements, output.src_instances);
    std::cout << "[" << node_id << "] " << "OUTPUT DESTINATIONS:" << std::endl;
    print_op_info(ctx, copy.dst_requirements, output.dst_instances);
  }

  virtual void map_task(const MapperContext ctx,
                        const Task& task,
                        const MapTaskInput& input,
                        MapTaskOutput& output) {
    if (STARTS_WITH(task.get_task_name(), "Flow_InitializeCell")) {
      std::cout << "[" << node_id << "] " << "MAPPING " << task.get_task_name() << std::endl;
      std::cout << "[" << node_id << "] " << "INPUT:" << std::endl;
      print_op_info(ctx, task.regions, input.valid_instances);
    }
    DefaultMapper::map_task(ctx, task, input, output);
    if (STARTS_WITH(task.get_task_name(), "Flow_InitializeCell")) {
      std::cout << "[" << node_id << "] " << "OUTPUT:" << std::endl;
      print_op_info(ctx, task.regions, output.chosen_instances);
      std::cout << "[" << node_id << "] " << "CHOSEN PROC: " << output.target_procs[0] << std::endl;
    }
  }

private:
  unsigned find_sample_id(const MapperContext ctx, const Task& task) const {
    for (const RegionRequirement& req : task.regions) {
      LogicalRegion region = req.region.exists() ? req.region
        : runtime->get_parent_logical_region(ctx, req.partition);
      region = get_root(ctx, region);
      const void* info = NULL;
      size_t info_size = 0;
      bool success = runtime->retrieve_semantic_information
        (ctx, region, SAMPLE_ID_TAG, info, info_size,
         false /*can_fail*/, true /*wait_until_ready*/);
      CHECK(success, "Missing SAMPLE_ID_TAG semantic information on region");
      assert(info_size == sizeof(unsigned));
      unsigned sample_id = *static_cast<const unsigned*>(info);
      assert(sample_id < sample_mappings_.size());
      return sample_id;
    }
    CHECK(false, "No region argument on task launch");
    return static_cast<unsigned>(-1);
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
