#include <array>
#include <deque>
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
// DOCUMENTATION
//=============================================================================

// Assume we're running 2 CPU-only samples, A and B, tiled as follows:
//   tiles(A) = [1,2,1], tilesPerRank(A) = [1,1,1]
//   tiles(B) = [6,1,1], tilesPerRank(B) = [3,1,1]
// Based on this configuration, we calculate the following:
//   #shards(A) = 2, #splintersPerShard(A) = 1
//   #shards(B) = 2, #splintersPerShard(B) = 3
// Each shard is placed on a separate rank in row-major order, so we will need
// 4 ranks in total. The splinters within each shard are allocated in
// round-robin, row-major order to the processors on the corresponding rank
// (some processors may receive more than 1 splinter). Assume each rank has 2
// CPU processors. Then the mapping will be as follows:
//   Sample    Tile -> Shard Splinter -> Rank CPU
//   --------------------------------------------
//        A [0,0,0] ->     0        0 ->   0    0
//   --------------------------------------------
//        A [0,1,0] ->     1        0 ->   1    0
//   --------------------------------------------
//        B [0,0,0] ->     0        0 ->   2    0
//        B [1,0,0] ->     0        1 ->   2    1
//        B [2,0,0] ->     0        2 ->   2    0
//   --------------------------------------------
//        B [3,0,0] ->     1        0 ->   3    0
//        B [4,0,0] ->     1        1 ->   3    1
//        B [5,0,0] ->     1        2 ->   3    0

//=============================================================================
// HELPER CODE
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
// INTRA-SAMPLE MAPPING
//=============================================================================

typedef unsigned SplinterID;

class SampleMapping;

class SplinteringFunctor : public ShardingFunctor {
private:
  static ShardingID NEXT_ID;
public:
  SplinteringFunctor(Runtime* rt, SampleMapping& parent)
    : id(NEXT_ID++), parent_(parent) {
    rt->register_sharding_functor(id, this, true);
  }
public:
  AddressSpace get_rank(const DomainPoint &point);
  virtual SplinterID splinter(const DomainPoint &point) = 0;
public:
  const ShardingID id;
protected:
  SampleMapping& parent_;
};

ShardingID SplinteringFunctor::NEXT_ID = 12345;

class SampleMapping {
public:
  class Tiling3DFunctor;
  class Tiling2DFunctor;
  class TileZeroFunctor;

public:
  SampleMapping(Runtime* rt, const Config& config, AddressSpace first_rank)
    : tiles_per_rank_{static_cast<unsigned>(config.Mapping.tilesPerRank[0]),
                      static_cast<unsigned>(config.Mapping.tilesPerRank[1]),
                      static_cast<unsigned>(config.Mapping.tilesPerRank[2])},
      ranks_per_dim_{static_cast<unsigned>(config.Mapping.tiles[0]
                                           / config.Mapping.tilesPerRank[0]),
                     static_cast<unsigned>(config.Mapping.tiles[1]
                                           / config.Mapping.tilesPerRank[1]),
                     static_cast<unsigned>(config.Mapping.tiles[2]
                                           / config.Mapping.tilesPerRank[2])},
      first_rank_(first_rank),
      tiling_3d_functor_(new Tiling3DFunctor(rt, *this)),
      tiling_2d_functors_{{new Tiling2DFunctor(rt, *this, 0, false),
                           new Tiling2DFunctor(rt, *this, 0, true )},
                          {new Tiling2DFunctor(rt, *this, 1, false),
                           new Tiling2DFunctor(rt, *this, 1, true )},
                          {new Tiling2DFunctor(rt, *this, 2, false),
                           new Tiling2DFunctor(rt, *this, 2, true )}},
      tile_zero_functor_(new TileZeroFunctor(rt, *this)) {}
  SampleMapping(const SampleMapping& rhs) = delete;
  SampleMapping& operator=(const SampleMapping& rhs) = delete;

public:
  AddressSpace get_rank(ShardID shard_id) const {
    return first_rank_ + shard_id;
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
  Tiling3DFunctor* tiling_3d_functor() {
    return tiling_3d_functor_;
  }
  Tiling2DFunctor* tiling_2d_functor(int dim, bool dir) {
    assert(0 <= dim && dim < 3);
    return tiling_2d_functors_[dim][dir];
  }
  TileZeroFunctor* tile_zero_functor() {
    return tile_zero_functor_;
  }

public:
  // Maps tasks in a 3D index space launch according to the default tiling
  // logic (see description above).
  class Tiling3DFunctor : public SplinteringFunctor {
  public:
    Tiling3DFunctor(Runtime* rt, SampleMapping& parent)
      : SplinteringFunctor(rt, parent) {}
  public:
    virtual ShardID shard(const DomainPoint& point,
                          const Domain& full_space,
                          const size_t total_shards) {
      assert(point.get_dim() == 3);
      CHECK(0 <= point[0] && point[0] < parent_.x_tiles() &&
            0 <= point[1] && point[1] < parent_.y_tiles() &&
            0 <= point[2] && point[2] < parent_.z_tiles(),
            "Unexpected point on index space launch");
      return
        (point[0] / parent_.tiles_per_rank_[0]) * parent_.ranks_per_dim_[1]
                                                * parent_.ranks_per_dim_[2] +
        (point[1] / parent_.tiles_per_rank_[1]) * parent_.ranks_per_dim_[2] +
        (point[2] / parent_.tiles_per_rank_[2]);
    }
    virtual SplinterID splinter(const DomainPoint &point) {
      assert(point.get_dim() == 3);
      CHECK(0 <= point[0] && point[0] < parent_.x_tiles() &&
            0 <= point[1] && point[1] < parent_.y_tiles() &&
            0 <= point[2] && point[2] < parent_.z_tiles(),
            "Unexpected point on index space launch");
      return
        (point[0] % parent_.tiles_per_rank_[0]) * parent_.tiles_per_rank_[1]
                                                * parent_.tiles_per_rank_[2] +
        (point[1] % parent_.tiles_per_rank_[1]) * parent_.tiles_per_rank_[2] +
        (point[2] % parent_.tiles_per_rank_[2]);
    }
  };

  // Maps tasks in a 2D index space launch, by extending each domain point to a
  // 3D tile and deferring to the default strategy.
  // Parameter `dim` controls which dimension to add.
  // Parameter `dir` controls which extreme of that dimension to set.
  class Tiling2DFunctor : public SplinteringFunctor {
  public:
    Tiling2DFunctor(Runtime* rt, SampleMapping& parent,
                    unsigned dim, bool dir)
      : SplinteringFunctor(rt, parent), dim_(dim), dir_(dir) {}
  public:
    virtual ShardID shard(const DomainPoint& point,
                          const Domain& full_space,
                          const size_t total_shards) {
      return parent_.tiling_3d_functor_->shard
        (to_point_3d(point), full_space, total_shards);
    }
    virtual SplinterID splinter(const DomainPoint &point) {
      return parent_.tiling_3d_functor_->splinter(to_point_3d(point));
    }
  private:
    DomainPoint to_point_3d(const DomainPoint& point) const {
      assert(point.get_dim() == 2);
      unsigned coord =
        (dim_ == 0) ? (dir_ ? 0 : parent_.x_tiles()-1) :
        (dim_ == 1) ? (dir_ ? 0 : parent_.y_tiles()-1) :
       /*dim_ == 2*/  (dir_ ? 0 : parent_.z_tiles()-1) ;
      return
        (dim_ == 0) ? Point<3>(coord, point[0], point[1]) :
        (dim_ == 1) ? Point<3>(point[0], coord, point[1]) :
       /*dim_ == 2*/  Point<3>(point[0], point[1], coord) ;
    }
  private:
    unsigned dim_;
    bool dir_;
  };

  // Maps every task to the shard & splinter corresponding to the first tile.
  class TileZeroFunctor : public SplinteringFunctor {
  public:
    TileZeroFunctor(Runtime* rt, SampleMapping& parent)
      : SplinteringFunctor(rt, parent) {}
  public:
    virtual ShardID shard(const DomainPoint& point,
                          const Domain& full_space,
                          const size_t total_shards) {
      return parent_.tiling_3d_functor_->shard
        (Point<3>(0,0,0), full_space, total_shards);
    }
    virtual SplinterID splinter(const DomainPoint &point) {
      return parent_.tiling_3d_functor_->splinter(Point<3>(0,0,0));
    }
  };

private:
  unsigned tiles_per_rank_[3];
  unsigned ranks_per_dim_[3];
  AddressSpace first_rank_;
  Tiling3DFunctor* tiling_3d_functor_;
  Tiling2DFunctor* tiling_2d_functors_[3][2];
  TileZeroFunctor* tile_zero_functor_;
};

AddressSpace SplinteringFunctor::get_rank(const DomainPoint &point) {
  return parent_.get_rank(shard(point, Domain(), 0));
}

//=============================================================================
// MAPPER CLASS: CONSTRUCTOR
//=============================================================================

class SoleilMapper : public DefaultMapper {
public:
  SoleilMapper(Runtime* rt, Machine machine, Processor local)
    : DefaultMapper(rt->get_mapper_runtime(), machine, local, "soleil_mapper"),
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
      sample_mappings_.emplace_back(rt, config, reqd_ranks);
    };
    // Locate all config files specified on the command-line arguments.
    InputArgs args = Runtime::get_input_args();
    for (int i = 0; i < args.argc; ++i) {
      if (EQUALS(args.argv[i], "-i") && i < args.argc-1) {
        Config config;
        parse_Config(&config, args.argv[i+1]);
        process_config(config);
        reqd_ranks += sample_mappings_.back().num_ranks();
      } else if (EQUALS(args.argv[i], "-m") && i < args.argc-1) {
        MultiConfig mc;
        parse_MultiConfig(&mc, args.argv[i+1]);
        process_config(mc.configs[0]);
        unsigned num_ranks_0 = sample_mappings_.back().num_ranks();
        if (!mc.collocateSections) {
          reqd_ranks += num_ranks_0;
        }
        process_config(mc.configs[1]);
        unsigned num_ranks_1 = sample_mappings_.back().num_ranks();
        if (mc.collocateSections) {
          reqd_ranks += std::max(num_ranks_0, num_ranks_1);
        } else {
          reqd_ranks += num_ranks_1;
        }
      }
    }
    // Verify that we have enough ranks.
    unsigned supplied_ranks = remote_cpus.size();
    CHECK(reqd_ranks <= supplied_ranks,
          "%u rank(s) required, but %u rank(s) supplied to Legion",
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
    // Verify machine configuration.
    for (AddressSpace rank = 0; rank < remote_cpus.size(); ++rank) {
      CHECK(get_procs(rank, Processor::IO_PROC).size() > 0,
            "No IO processor on rank %u", rank);
    }
  }

//=============================================================================
// MAPPER CLASS: MAPPING LOGIC
//=============================================================================

private:
  std::vector<unsigned> find_sample_ids(const MapperContext ctx,
                                        const Task& task) const {
    CHECK(task.has_parent_task(),
          "Unexpected task in find_sample_ids: %s", task.get_task_name());
    std::vector<unsigned> sample_ids;
    // Tasks with Config as 1st argument: read config.Mapping.sampleId
    if (EQUALS(task.get_task_name(), "workSingle")) {
      const Config* config = static_cast<const Config*>(first_arg(task));
      sample_ids.push_back(static_cast<unsigned>(config->Mapping.sampleId));
    }
    // Tasks with MultiConfig as 1st argument: read configs[*].Mapping.sampleId
    else if (EQUALS(task.get_task_name(), "workDual")) {
      const MultiConfig* mc = static_cast<const MultiConfig*>(first_arg(task));
      sample_ids.push_back
        (static_cast<unsigned>(mc->configs[0].Mapping.sampleId));
      sample_ids.push_back
        (static_cast<unsigned>(mc->configs[1].Mapping.sampleId));
    }
    // Other tasks: go up one level to the work task
    else {
      sample_ids = find_sample_ids(ctx, *(task.get_parent_task()));
    }
    // Sanity checks
    assert(!sample_ids.empty());
    for (unsigned sample_id : sample_ids) {
      assert(sample_id < sample_mappings_.size());
    }
    return sample_ids;
  }

  unsigned find_sample_id(const MapperContext ctx, const Task& task) const {
    return find_sample_ids(ctx, task)[0];
  }

  SplinteringFunctor* pick_functor(const MapperContext ctx,
                                   const Task& task) {
    unsigned sample_id = find_sample_id(ctx, task);
    SampleMapping& mapping = sample_mappings_[sample_id];
    // 3D index space tasks
    if ((task.is_index_space && task.index_domain.get_dim() == 3) ||
        task.index_point.get_dim() == 3) {
      return mapping.tiling_3d_functor();
    }
    // 2D index space tasks
    else if ((task.is_index_space && task.index_domain.get_dim() == 2) ||
             task.index_point.get_dim() == 2) {
      unsigned dim = parse_dimension(task);
      bool dir = parse_direction(task)[dim];
      if (STARTS_WITH(task.get_task_name(), "initialize_faces_") ||
          STARTS_WITH(task.get_task_name(), "bound_")) {
        // Do nothing.
      } else if (STARTS_WITH(task.get_task_name(), "cache_intensity_")) {
        // We want to run these tasks on the opposite end of the domain implied
        // by their name.
        dir = !dir;
      } else {
        CHECK(false, "Unexpected 2D domain on index space launch of task %s",
              task.get_task_name());
      }
      return mapping.tiling_2d_functor(dim, dir);
    }
    // Sample-specific tasks that are launched individually
    else {
      CHECK(task.index_point.get_dim() <= 1,
            "Unexpected index information on task %s", task.get_task_name());
      return mapping.tile_zero_functor();
    }
  }

//=============================================================================
// MAPPER CLASS: MAJOR OVERRIDES
//=============================================================================

public:
  // Control-replicate work tasks.
  virtual void select_task_options(const MapperContext ctx,
                                   const Task& task,
                                   TaskOptions& output) {
    DefaultMapper::select_task_options(ctx, task, output);
    output.replicate =
      EQUALS(task.get_task_name(), "workSingle") ||
      EQUALS(task.get_task_name(), "workDual");
  }

  // Enable tracing.
  virtual void memoize_operation(const MapperContext ctx,
                                 const Mappable& mappable,
                                 const MemoizeInput& input,
                                 MemoizeOutput& output) {
    output.memoize = true;
  }

  virtual void default_policy_rank_processor_kinds(
                              MapperContext ctx,
                              const Task& task,
                              std::vector<Processor::Kind>& ranking) {
    // Work tasks: map to IO processors, so they don't get blocked by tiny
    // CPU tasks.
    if (EQUALS(task.get_task_name(), "workSingle") ||
        EQUALS(task.get_task_name(), "workDual")) {
      ranking.push_back(Processor::IO_PROC);
    }
    // Other tasks: defer to the default mapping policy
    else {
      DefaultMapper::default_policy_rank_processor_kinds(ctx, task, ranking);
    }
  }

#ifndef NO_LEGION_CONTROL_REPLICATION
  // Replicate each work task over all ranks assigned to the corresponding
  // sample(s).
  virtual void map_replicate_task(const MapperContext ctx,
                                  const Task& task,
                                  const MapTaskInput& input,
                                  const MapTaskOutput& default_output,
                                  MapReplicateTaskOutput& output) {
    // Read configuration.
    assert(!runtime->is_MPI_interop_configured(ctx));
    assert(EQUALS(task.get_task_name(), "workSingle") ||
           EQUALS(task.get_task_name(), "workDual"));
    VariantInfo info =
      default_find_preferred_variant(task, ctx, false/*needs_tight_bound*/);
    CHECK(task.regions.empty() && info.is_replicable,
          "Unexpected features on work task");
    std::vector<unsigned> sample_ids = find_sample_ids(ctx, task);
    // Create a replicant on the first CPU processor of each sample's ranks.
    for (unsigned sample_id : sample_ids) {
      const SampleMapping& mapping = sample_mappings_[sample_id];
      for (ShardID shard_id = 0; shard_id < mapping.num_ranks(); ++shard_id) {
        AddressSpace rank = mapping.get_rank(shard_id);
        Processor target_proc = get_procs(rank, info.proc_kind)[0];
        output.task_mappings.push_back(default_output);
        output.task_mappings.back().chosen_variant = info.variant;
        output.task_mappings.back().target_procs.push_back(target_proc);
        output.control_replication_map.push_back(target_proc);
      }
    }
  }
#endif

  // NOTE: Will only run if Legion is compiled with dynamic control replication.
  virtual void select_sharding_functor(const MapperContext ctx,
                                       const Task& task,
                                       const SelectShardingFunctorInput& input,
                                       SelectShardingFunctorOutput& output) {
    output.chosen_functor = pick_functor(ctx, task)->id;
  }

  virtual Processor default_policy_select_initial_processor(
                              MapperContext ctx,
                              const Task& task) {
    // Index space tasks: defer to the default mapping policy; slice_task will
    // eventually be called to do the mapping properly
    if (task.is_index_space) {
      return DefaultMapper::default_policy_select_initial_processor(ctx, task);
    }
    // Main task: defer to the default mapping policy
    else if (EQUALS(task.get_task_name(), "main")) {
      return DefaultMapper::default_policy_select_initial_processor(ctx, task);
    }
    // Other tasks: consult functor
    else {
      VariantInfo info =
        default_find_preferred_variant(task, ctx, false/*needs_tight_bound*/);
      SplinteringFunctor* functor = pick_functor(ctx, task);
      Processor target_proc =
        select_proc(task.index_point, info.proc_kind, functor);
      return target_proc;
    }
  }

  virtual void slice_task(const MapperContext ctx,
                          const Task& task,
                          const SliceTaskInput& input,
                          SliceTaskOutput& output) {
    output.verify_correctness = false;
    VariantInfo info =
      default_find_preferred_variant(task, ctx, false/*needs_tight_bound*/);
    SplinteringFunctor* functor = pick_functor(ctx, task);
    for (Domain::DomainPointIterator it(input.domain); it; it++) {
      Processor target_proc = select_proc(it.p, info.proc_kind, functor);
      output.slices.emplace_back(Domain(it.p, it.p), target_proc,
                                 false/*recurse*/, false/*stealable*/);
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
      assert(!task.regions.empty() && task.regions[0].region.exists());
      unsigned sample_id = find_sample_id(ctx, task);
      const SampleMapping& mapping = sample_mappings_[sample_id];
      std::array<bool,3> dir = parse_direction(task);
      DomainPoint tile =
        runtime->get_logical_region_color_point(ctx, task.regions[0].region);
      priority =
        (dir[0] ? mapping.x_tiles() - tile[0] - 1 : tile[0]) +
        (dir[1] ? mapping.y_tiles() - tile[1] - 1 : tile[1]) +
        (dir[2] ? mapping.z_tiles() - tile[2] - 1 : tile[2]) ;
    }
    // Increase priority of tasks on the critical path of the fluid solve.
    if (STARTS_WITH(task.get_task_name(), "Flow_ComputeVelocityGradient") ||
        STARTS_WITH(task.get_task_name(), "Flow_UpdateGhostVelocityGradient") ||
        STARTS_WITH(task.get_task_name(), "Flow_GetFlux") ||
        STARTS_WITH(task.get_task_name(), "Flow_UpdateUsingFlux")) {
      priority = 1;
    }
    return priority;
  }

  // Send each fill to the first rank of the first section.
  // NOTE: Will only run if Legion is compiled with dynamic control replication.
  virtual void select_sharding_functor(const MapperContext ctx,
                                       const Fill& fill,
                                       const SelectShardingFunctorInput& input,
                                       SelectShardingFunctorOutput& output) {
    CHECK(fill.parent_task != NULL &&
          (EQUALS(fill.parent_task->get_task_name(), "workDual") ||
           EQUALS(fill.parent_task->get_task_name(), "workSingle")) &&
          !fill.is_index_space &&
          fill.requirement.region.exists() &&
          runtime->get_index_space_depth
            (ctx, fill.requirement.region.get_index_space()) == 0,
          "Unexpected argument on fill");
    unsigned sample_id = find_sample_id(ctx, *(fill.parent_task));
    SampleMapping& mapping = sample_mappings_[sample_id];
    output.chosen_functor = mapping.tile_zero_functor()->id;
  }

  // Send each dependent partitioning operation to the first rank of the first
  // section.
  // NOTE: Will only run if Legion is compiled with dynamic control replication.
  virtual void select_sharding_functor(const MapperContext ctx,
                                       const Partition& partition,
                                       const SelectShardingFunctorInput& input,
                                       SelectShardingFunctorOutput& output) {
    CHECK(partition.parent_task != NULL &&
          (EQUALS(partition.parent_task->get_task_name(), "workDual") ||
           EQUALS(partition.parent_task->get_task_name(), "workSingle")) &&
          !partition.is_index_space &&
          partition.requirement.region.exists() &&
          runtime->get_index_space_depth
            (ctx, partition.requirement.region.get_index_space()) == 0,
          "Unexpected argument on partition");
    unsigned sample_id = find_sample_id(ctx, *(partition.parent_task));
    SampleMapping& mapping = sample_mappings_[sample_id];
    output.chosen_functor = mapping.tile_zero_functor()->id;
  }

//=============================================================================
// MAPPER CLASS: MINOR OVERRIDES
//=============================================================================

public:
  // TODO: Select appropriate memories for instances that will be communicated,
  // (e.g. parallelizer-created ghost partitions), such as RDMA memory,
  // zero-copy memory.
  // virtual Memory default_policy_select_target_memory(...) { ... }

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

  // Disable an optimization done by the default mapper (extends the set of
  // eligible processors to include all the processors of the same type on the
  // target node).
  virtual void default_policy_select_target_processors(
                              MapperContext ctx,
                              const Task &task,
                              std::vector<Processor> &target_procs) {
    target_procs.push_back(task.target_proc);
  }

  // Shouldn't have to shard any of the following operations.
  virtual void select_sharding_functor(const MapperContext ctx,
                                       const Copy& copy,
                                       const SelectShardingFunctorInput& input,
                                       SelectShardingFunctorOutput& output) {
    CHECK(false, "Unsupported: Sharded Copy");
  }
  virtual void select_sharding_functor(const MapperContext ctx,
                                       const Close& close,
                                       const SelectShardingFunctorInput& input,
                                       SelectShardingFunctorOutput& output) {
    CHECK(false, "Unsupported: Sharded Close");
  }
  virtual void select_sharding_functor(const MapperContext ctx,
                                       const Acquire& acquire,
                                       const SelectShardingFunctorInput& input,
                                       SelectShardingFunctorOutput& output) {
    CHECK(false, "Unsupported: Sharded Acquire");
  }
  virtual void select_sharding_functor(const MapperContext ctx,
                                       const Release& release,
                                       const SelectShardingFunctorInput& input,
                                       SelectShardingFunctorOutput& output) {
    CHECK(false, "Unsupported: Sharded Release");
  }
  virtual void select_sharding_functor(const MapperContext ctx,
                                       const MustEpoch& epoch,
                                       const SelectShardingFunctorInput& input,
                                       MustEpochShardingFunctorOutput& output) {
    CHECK(false, "Unsupported: Sharded MustEpoch");
  }

//=============================================================================
// MAPPER CLASS: HELPER METHODS
//=============================================================================

private:
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
      dir[0] = 1 - ((quadrant >> 2) & 1);
      dir[1] = 1 - ((quadrant >> 1) & 1);
      dir[2] = 1 - ((quadrant >> 0) & 1);
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

  // NOTE: This function doesn't sanity check its input.
  Processor select_proc(const DomainPoint& tile,
                        Processor::Kind kind,
                        SplinteringFunctor* functor) {
    AddressSpace rank = functor->get_rank(tile);
    const std::vector<Processor>& procs = get_procs(rank, kind);
    SplinterID splinter_id = functor->splinter(tile);
    return procs[splinter_id % procs.size()];
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

//=============================================================================
// MAPPER CLASS: MEMBER VARIABLES
//=============================================================================

private:
  std::deque<SampleMapping> sample_mappings_;
  std::vector<std::vector<std::vector<Processor> > > all_procs_;
};

//=============================================================================
// MAPPER REGISTRATION
//=============================================================================

static void create_mappers(Machine machine,
                           Runtime* rt,
                           const std::set<Processor>& local_procs) {
  rt->replace_default_mapper(new SoleilMapper(rt, machine, *(local_procs.begin())));
}

void register_mappers() {
  Runtime::add_registration_callback(create_mappers);
}
