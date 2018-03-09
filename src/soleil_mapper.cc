#include <string.h>

#include "mappers/default_mapper.h"
#include "realm/logging.h"

#include "config_schema.h"
#include "soleil_mapper.h"

using namespace Legion;
using namespace Legion::Mapping;

//=============================================================================

static Realm::Logger LOG("soleil_mapper");

//=============================================================================

struct SampleMapping {
  AddressSpace first_node;
  unsigned num_tiles;
};

// Maps every sample to a disjoint set of nodes. Each sample gets one node for
// every tile, as specified in its config file.
class SoleilMapper : public DefaultMapper {
public:
  SoleilMapper(MapperRuntime* rt,
               Machine machine,
               Processor local)
    : DefaultMapper(rt, machine, local, "soleil_mapper") {
    InputArgs args = Runtime::get_input_args();
    unsigned allocated_nodes = 0;
    for (int i = 0; i < args.argc; ++i) {
      if (strcmp(args.argv[i], "-i") == 0 && i < args.argc-1) {
        Config config = parse_config(args.argv[i+1]);
        int num_tiles = config.Mapping.xTiles *
                        config.Mapping.yTiles *
                        config.Mapping.zTiles;
        assert(num_tiles > 0);
        sample_mappings_.push_back(SampleMapping{
          .first_node = allocated_nodes,
          .num_tiles = static_cast<unsigned>(num_tiles)
        });
        allocated_nodes += num_tiles;
      }
    }
    unsigned total_nodes = remote_cpus.size();
    if (allocated_nodes != total_nodes) {
      LOG.error("%d node(s) required, but %d node(s) supplied to Legion.",
                allocated_nodes, total_nodes);
      assert(false);
    }
    // TODO: Verify we're running with 1 OpenMP/CPU processor per node.
  }

public:
  // Send each work task to the first in the set of nodes allocated to the
  // corresponding sample.
  virtual Processor default_policy_select_initial_processor(
                              MapperContext ctx,
                              const Task& task) {
    if (strcmp(task.get_task_name(), "work") != 0) {
      return DefaultMapper::default_policy_select_initial_processor(ctx, task);
    }
    const char* ptr = static_cast<const char*>(task.args);
    const Config* config =
      reinterpret_cast<const Config*>(ptr + sizeof(uint64_t));
    int sample_id = config->Mapping.sampleId;
    assert(sample_id >= 0);
    AddressSpace target_node = sample_mappings_[sample_id].first_node;
    Processor target_proc = remote_cpus[target_node];
    return target_proc;
  }

  // TODO: Select RDMA memory for parallelizer-created ghost partitions.
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

  // Farm index-space launches made by work tasks across all the nodes
  // allocated to the corresponding sample.
  // TODO: Cache the decision.
  virtual void slice_task(const MapperContext ctx,
                          const Task& task,
                          const SliceTaskInput& input,
                          SliceTaskOutput& output) {
    if (task.parent_task == NULL ||
        strcmp(task.parent_task->get_task_name(), "work") != 0) {
      DefaultMapper::slice_task(ctx, task, input, output);
      return;
    }
    output.verify_correctness = false;
    const char* ptr = static_cast<const char*>(task.parent_task->args);
    const Config* config =
      reinterpret_cast<const Config*>(ptr + sizeof(uint64_t));
    int sample_id = config->Mapping.sampleId;
    assert(sample_id >= 0);
    const SampleMapping& mapping = sample_mappings_[sample_id];
    if (mapping.num_tiles != input.domain.get_volume()) {
      LOG.error("Index-space launches in the work task should only use the "
                "top-level tiling.");
      assert(false);
    }
    std::vector<Processor> procs;
    for (unsigned i = 0; i < mapping.num_tiles; ++i) {
      procs.push_back(remote_cpus[mapping.first_node + i]);
    }
    DomainT<3,coord_t> point_space = input.domain;
    Point<3,coord_t> blocking =
      default_select_num_blocks<3>(mapping.num_tiles, point_space.bounds);
    default_decompose_points<3>(input.domain, procs, blocking,
                                false/*recurse*/, false/*stealable*/,
                                output.slices);
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
