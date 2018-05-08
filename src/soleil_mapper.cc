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

//=============================================================================

struct SampleMapping {
  AddressSpace first_rank;
};

// Maps every sample to a disjoint set of ranks. Each sample gets one rank for
// every tile, as specified in its config file.
class SoleilMapper : public DefaultMapper {
public:
  SoleilMapper(MapperRuntime* rt,
               Machine machine,
               Processor local)
    : DefaultMapper(rt, machine, local, "soleil_mapper") {
    // Set the umask of the process to clear S_IWGRP and S_IWOTH.
    umask(022);
    // Assign ranks sequentially to samples, each sample getting one rank for
    // each tile.
    InputArgs args = Runtime::get_input_args();
    unsigned allocated_ranks = 0;
    auto process_config = [&](char* config_file) {
      Config config = parse_config(config_file);
      CHECK(config.Mapping.xTiles > 0 &&
            config.Mapping.yTiles > 0 &&
            config.Mapping.zTiles > 0,
            "Invalid tiling");
      sample_mappings_.push_back(SampleMapping{
        .first_rank = allocated_ranks,
      });
      allocated_ranks += config.Mapping.xTiles
                       * config.Mapping.yTiles
                       * config.Mapping.zTiles;
    };
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
    unsigned total_ranks = remote_cpus.size();
    CHECK(allocated_ranks <= total_ranks,
          "%d rank(s) required, but %d rank(s) supplied to Legion",
          allocated_ranks, total_ranks);
    if (allocated_ranks < total_ranks) {
      LOG.warning() << total_ranks - allocated_ranks << " rank(s) are unused";
    }
    // TODO: Verify we're running with 1 OpenMP/CPU processor per rank.
  }

public:
  virtual Processor default_policy_select_initial_processor(
                              MapperContext ctx,
                              const Task& task) {
    // Send each work task to the first in the set of ranks allocated to the
    // corresponding sample.
    if (strcmp(task.get_task_name(), "work") == 0) {
      const Config* config = find_config(&task);
      unsigned sample_id = config->Mapping.sampleId;
      assert(sample_id < sample_mappings_.size());
      AddressSpace target_rank = sample_mappings_[sample_id].first_rank;
      Processor target_proc = remote_cpus[target_rank];
      LOG.debug() << "Sample " << sample_id << ":"
                  << " Sequential launch: Work task"
		  << " for sample " << sample_id
		  << " mapped to rank " << target_rank;
      return target_proc;
    }
    // For other tasks, defer to the default mapping policy.
    return DefaultMapper::default_policy_select_initial_processor(ctx, task);
  }

  // TODO: Assign priorities to sweep tasks such that we prioritize the tile
  // that has more dependencies downstream.
  virtual TaskPriority default_policy_select_task_priority(
                              MapperContext ctx,
                              const Task& task) {
    return DefaultMapper::default_policy_select_task_priority(ctx, task);
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
    CHECK(input.domain.get_dim() == 3,
          "Index-space launches should only use the top-level tiling.");
    output.verify_correctness = false;
    // Retrieve sample information from parent task.
    const Config* config = find_config(&task);
    unsigned sample_id = config->Mapping.sampleId;
    assert(sample_id < sample_mappings_.size());
    const SampleMapping& mapping = sample_mappings_[sample_id];
    // Select the 1st (and only) processor of the same kind as the original
    // target, on each rank.
    const std::vector<Processor>& procs =
      remote_procs(task.target_proc.kind());
    // Distribute tiles to the ranks in row-major order.
    for (Domain::DomainPointIterator itr(input.domain); itr; itr++) {
      CHECK(itr.p[0] >= 0 && itr.p[0] < config->Mapping.xTiles &&
	    itr.p[1] >= 0 && itr.p[1] < config->Mapping.yTiles &&
	    itr.p[2] >= 0 && itr.p[2] < config->Mapping.zTiles,
            "Out-of-bounds point on index-space launch");
      int rank = mapping.first_rank
               + itr.p[0] * config->Mapping.yTiles * config->Mapping.zTiles
               + itr.p[1] * config->Mapping.zTiles
               + itr.p[2];
      output.slices.emplace_back(Domain(itr.p, itr.p),
                                 procs[rank],
                                 false /*recurse*/,
                                 false /*stealable*/);
      LOG.debug() << "Sample " << sample_id << ":"
                  << " Index-space launch: Task " << task.get_task_name()
                  << " on tile " << itr.p
                  << " mapped to rank " << rank;
    }
  }

private:
  // TODO: This interface only returns the first processor of the desired kind
  // on each machine; also handle the case where there's more (e.g. >1 GPUs).
  const std::vector<Processor>& remote_procs(Processor::Kind kind) {
    switch (kind) {
    case TOC_PROC: return remote_gpus;
    case LOC_PROC: return remote_cpus;
    case IO_PROC:  return remote_ios;
    case PROC_SET: return remote_procsets;
    case OMP_PROC: return remote_omps;
    case PY_PROC:  return remote_pys;
    default:       assert(false);
    }
    // Should never reach here
    return remote_cpus;
  }

  const Config* find_config(const Task* task) {
    for(; task != NULL; task = task->parent_task) {
      if (strcmp(task->get_task_name(), "work") == 0) {
        const char* ptr = static_cast<const char*>(task->args);
        return reinterpret_cast<const Config*>(ptr + sizeof(uint64_t));
      }
    }
    return NULL;
  }

private:
  std::vector<SampleMapping> sample_mappings_;
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
