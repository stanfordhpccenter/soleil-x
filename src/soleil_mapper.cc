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

struct TaskMapping {
  Processor::Kind proc_kind;
  unsigned first_proc_idx;
  unsigned num_procs;
};

struct ShardMapping {
  Processor::Kind proc_kind;
  unsigned first_node_idx;
  unsigned num_nodes;
#ifdef CONTROL_REPLICATION
  ShardingID sharding_functor_id;
#endif
};

struct SampleMapping {
  ShardMapping shard_mapping;
  TaskMapping task_mapping;
};

struct MachineModel {
  MachineModel(Machine machine);

  ShardMapping generate_shard_mapping(const TaskMapping &task_mapping) const;

  Processor pick_io_proc(unsigned node_id) const;
  Processor pick_cpu(unsigned node_id) const;

  std::vector<Processor> cpus_list;
  std::vector<Processor> ocpus_list;
  std::vector<Processor> gpus_list;
  std::vector<Processor> io_procs_list;
  std::vector<Memory> sysmem_list;
  unsigned num_nodes;
};

class SoleilMapper : public DefaultMapper {
public:
  SoleilMapper(MapperRuntime* rt,
               Machine machine,
               Processor local,
               bool use_dcr,
               MachineModel *_machine_model,
               std::vector<SampleMapping> *_mappings)
    : DefaultMapper(rt, machine, local, "soleil_mapper"), use_dcr_(use_dcr),
      machine_model_(*_machine_model), sample_mappings_(*_mappings) {
  }

  virtual void select_task_options(const MapperContext    ctx,
                                   const Task&            task,
                                   TaskOptions&     output)
  {
    DefaultMapper::select_task_options(ctx, task, output);
#ifdef CONTROL_REPLICATION
    output.replicate = false;
    if (use_dcr_ && strcmp(task.get_task_name(), "work") == 0)
      output.replicate = true;
#endif
  }

public:
  virtual Processor default_policy_select_initial_processor(
                              MapperContext ctx,
                              const Task& task) {
    if (strcmp(task.get_task_name(), "work") == 0)
    {
      const Config* config = find_config(&task);
      unsigned sample_id = config->Mapping.sampleId;
      assert(sample_id < sample_mappings_.size());
      const ShardMapping &shard_mapping =
        sample_mappings_[sample_id].shard_mapping;
      if (shard_mapping.proc_kind == Processor::IO_PROC)
        return machine_model_.pick_io_proc(shard_mapping.first_node_idx);
      else
      {
        assert(shard_mapping.proc_kind == Processor::LOC_PROC);
        return machine_model_.pick_cpu(shard_mapping.first_node_idx);
      }
    }
    else if (task.parent_task != NULL)
    {
      return task.parent_task->target_proc;
    }
    // For other tasks, defer to the default mapping policy.
    return DefaultMapper::default_policy_select_initial_processor(ctx, task);
  }

  virtual void default_policy_select_target_processors(
                                MapperContext ctx,
                                const Task &task,
                                std::vector<Processor> &target_procs)
  {
    target_procs.push_back(task.target_proc);
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

  typedef std::pair<unsigned,TaskID> UniqueTaskId;
  typedef std::pair<Domain,UniqueTaskId> SliceCacheKey;

  // Farm index-space launches made by work tasks across all the ranks
  // allocated to the corresponding sample.
  // TODO: Cache the decision.
  virtual void slice_task(const MapperContext ctx,
                          const Task& task,
                          const SliceTaskInput& input,
                          SliceTaskOutput& output)
  {
    if (input.domain.get_dim() != 3)
    {
      DefaultMapper::slice_task(ctx, task, input, output);
      return;
    }

    CHECK(task.parent_task != NULL &&
          strcmp(task.parent_task->get_task_name(), "work") == 0,
          "Index-space launches only allowed in work task");

    // Retrieve sample information from parent task.
    const Config* config = find_config(&task);
    unsigned sample_id = config->Mapping.sampleId;
    assert(sample_id < sample_mappings_.size());

    output.verify_correctness = false;

    SliceCacheKey key(input.domain, UniqueTaskId(sample_id, task.task_id));
    std::map<SliceCacheKey,std::vector<TaskSlice> >::const_iterator
      finder = cached_slices_.find(key);
    if (finder != cached_slices_.end()) {
      output.slices = finder->second;
      return;
    }

    const TaskMapping& mapping = sample_mappings_[sample_id].task_mapping;
    bool has_variant =
      have_proc_kind_variant(ctx, task.task_id, mapping.proc_kind);
    unsigned num_tiles = mapping.num_procs;
    const std::vector<Processor> *procs = mapping.proc_kind == Processor::OMP_PROC
                                          ? &machine_model_.ocpus_list
                                          : &machine_model_.gpus_list;
    if (!has_variant && mapping.proc_kind == Processor::TOC_PROC)
    {
      has_variant =
        have_proc_kind_variant(ctx, task.task_id, Processor::OMP_PROC);
      if (has_variant)
        procs = &machine_model_.ocpus_list;
    }
    if (!has_variant)
      procs = &machine_model_.cpus_list;
    assert(procs != NULL);
    for (Domain::DomainPointIterator itr(input.domain); itr; itr++)
    {
      DomainPoint p = itr.p;
      p[0] = p[0] < config->Mapping.xTiles ? p[0] : config->Mapping.xTiles - 1;
      p[1] = p[1] < config->Mapping.yTiles ? p[1] : config->Mapping.yTiles - 1;
      p[2] = p[2] < config->Mapping.zTiles ? p[2] : config->Mapping.zTiles - 1;
      size_t linearized = p[0] +
                          p[1] * config->Mapping.xTiles +
                          p[2] * config->Mapping.xTiles * config->Mapping.yTiles;

      size_t proc_idx = mapping.first_proc_idx + linearized;

      if (procs->size() >= num_tiles)
      {
        assert(procs->size() % num_tiles == 0);
        proc_idx *= procs->size() / num_tiles;
      }
      else // procs->size() < num_tiles
      {
        assert(num_tiles % procs->size() == 0);
        proc_idx /= num_tiles / procs->size();
      }
      Processor target_proc = (*procs)[proc_idx % procs->size()];
      output.slices.emplace_back(Domain(itr.p, itr.p), target_proc, false, false);
    }
    cached_slices_[key] = output.slices;
  }

#ifdef CONTROL_REPLICATION
  virtual void map_replicate_task(const MapperContext      ctx,
                                  const Task&              task,
                                  const MapTaskInput&      input,
                                  const MapTaskOutput&     def_output,
                                  MapReplicateTaskOutput&  output)
  {
    assert(strcmp(task.get_task_name(), "work") == 0);
    assert(!runtime->is_MPI_interop_configured(ctx));

    const Config* config = find_config(&task);
    unsigned sample_id = config->Mapping.sampleId;
    assert(sample_id < sample_mappings_.size());
    const ShardMapping& mapping = sample_mappings_[sample_id].shard_mapping;

    Processor::Kind target_kind = mapping.proc_kind;
    // Get the variant that we are going to use to map this task
    VariantInfo chosen = default_find_preferred_variant(task, ctx,
                      true/*needs tight bound*/, true/*cache*/, target_kind);
    assert(chosen.is_replicable);

    unsigned first_node = mapping.first_node_idx;
    unsigned num_nodes = mapping.num_nodes;
    output.task_mappings.resize(num_nodes, def_output);
    output.control_replication_map.resize(num_nodes);

    for (size_t idx = 0; idx < num_nodes; ++idx)
    {
      unsigned node_id = (first_node + idx) % machine_model_.num_nodes;
      Processor target_proc = target_kind == Processor::IO_PROC ?
                              machine_model_.pick_io_proc(node_id) :
                              machine_model_.pick_cpu(node_id);
      output.control_replication_map[idx] = target_proc;
      output.task_mappings[idx].chosen_variant = chosen.variant;
      output.task_mappings[idx].target_procs.push_back(target_proc);
    }
  }
#endif

#ifdef CONTROL_REPLICATION
  virtual void select_sharding_functor(
                                 const MapperContext                ctx,
                                 const Task&                        task,
                                 const SelectShardingFunctorInput&  input,
                                       SelectShardingFunctorOutput& output)
  {
    const Config* config = find_config(&task);
    if (config != NULL && task.is_index_space && task.index_domain.dim == 3)
    {
      const SampleMapping &sample_mapping =
        sample_mappings_[config->Mapping.sampleId];
      output.chosen_functor =
        sample_mapping.shard_mapping.sharding_functor_id;
    }
    else
      DefaultMapper::select_sharding_functor(ctx, task, input, output);
  }
#endif

private:

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
  bool use_dcr_;
  MachineModel &machine_model_;
  std::vector<SampleMapping> &sample_mappings_;
  std::map<SliceCacheKey,std::vector<TaskSlice> > cached_slices_;
};

//=============================================================================

// Cache part of the machine model that Realm provides for mapping purposes
MachineModel::MachineModel(Machine machine)
{
  static const unsigned NUM_PROC_KINDS = 4;
  std::vector<Processor>* p_procs[NUM_PROC_KINDS] =
    { &cpus_list, &ocpus_list, &gpus_list, &io_procs_list };
  static const Processor::Kind query_kinds[] =
    { Processor::LOC_PROC, Processor::OMP_PROC,
      Processor::TOC_PROC, Processor::IO_PROC };

  for (unsigned idx = 0; idx < NUM_PROC_KINDS; ++idx)
  {
    Machine::ProcessorQuery procs_query(machine);
    procs_query.only_kind(query_kinds[idx]);
    for (Machine::ProcessorQuery::iterator it = procs_query.begin();
          it != procs_query.end(); it++)
      p_procs[idx]->push_back(*it);
  }

  {
    Machine::MemoryQuery mems_query(machine);
    mems_query.only_kind(Memory::SYSTEM_MEM);
    for (Machine::MemoryQuery::iterator it = mems_query.begin();
          it != mems_query.end(); it++)
      sysmem_list.push_back(*it);
  }
  num_nodes = sysmem_list.size();
}

ShardMapping
MachineModel::generate_shard_mapping(const TaskMapping &task_mapping) const
{
  std::vector<unsigned> next_shard_proc(0, num_nodes);
  next_shard_proc.resize(num_nodes);
  ShardMapping shard_mapping;

  if (io_procs_list.size() > 0)
    shard_mapping.proc_kind = Processor::IO_PROC;
  else
    shard_mapping.proc_kind = Processor::LOC_PROC;

  unsigned num_procs_per_node = task_mapping.proc_kind == Processor::OMP_PROC ?
                                ocpus_list.size() / num_nodes :
                                gpus_list.size() / num_nodes;
  shard_mapping.first_node_idx =
    task_mapping.first_proc_idx / num_procs_per_node;
  shard_mapping.num_nodes =
    (task_mapping.num_procs + num_procs_per_node - 1) / num_procs_per_node;

  return shard_mapping;
}

Processor
MachineModel::pick_io_proc(unsigned node_id) const
{
  unsigned num_io_procs_per_node = io_procs_list.size() / num_nodes;
  return io_procs_list[num_io_procs_per_node * node_id];
}

Processor
MachineModel::pick_cpu(unsigned node_id) const
{
  unsigned num_cpus_per_node = cpus_list.size() / num_nodes;
  return cpus_list[num_cpus_per_node * node_id];
}

#ifdef CONTROL_REPLICATION
class TilingFunctor : public ShardingFunctor {
  public:
    TilingFunctor(const Config &config, size_t num_shards);
    virtual ~TilingFunctor(void);
  public:
    virtual ShardID shard(const DomainPoint &point,
                          const Domain &full_space,
                          const size_t total_shards);
  protected:
    Config config;
    size_t num_shards;
    size_t num_points_per_shard;
};

TilingFunctor::TilingFunctor(const Config &_config, size_t _num_shards)
  : ShardingFunctor(), config(_config), num_shards(_num_shards)
{
  num_points_per_shard = config.Mapping.xTiles *
                         config.Mapping.yTiles *
                         config.Mapping.zTiles;
  num_points_per_shard =
    (num_points_per_shard + num_shards - 1) / num_shards;
}

TilingFunctor::~TilingFunctor()
{
}

ShardID TilingFunctor::shard(const DomainPoint &point,
                             const Domain &full_space,
                             const size_t total_shards)
{
  assert(point.dim == 3);
  DomainPoint p = point;
  p[0] = p[0] < config.Mapping.xTiles ? p[0] : config.Mapping.xTiles - 1;
  p[1] = p[1] < config.Mapping.yTiles ? p[1] : config.Mapping.yTiles - 1;
  p[2] = p[2] < config.Mapping.zTiles ? p[2] : config.Mapping.zTiles - 1;
  coord_t linearized = p[0] +
                       p[1] * config.Mapping.xTiles +
                       p[2] * config.Mapping.xTiles * config.Mapping.yTiles;
  ShardID shard_id = (ShardID)(linearized / num_points_per_shard);
  assert(shard_id < total_shards);
  return shard_id;
}
#endif

static bool parse_configurations(std::vector<Config> &configs)
{
  // Set the umask of the process to clear S_IWGRP and S_IWOTH.
  umask(022);

  InputArgs args = Runtime::get_input_args();
  unsigned sampleId = 0;
  bool use_dcr = false;
  for (int i = 0; i < args.argc; ++i) {
    if (strcmp(args.argv[i], "-i") == 0 && i < args.argc-1) {
      configs.push_back(parse_config(args.argv[i+1]));
      configs.back().Mapping.sampleId = sampleId++;
    } else if (strcmp(args.argv[i], "-I") == 0 && i < args.argc-1) {
      std::ifstream csv_file(args.argv[i+1]);
      std::string json_filename;
      while (std::getline(csv_file, json_filename)) {
        configs.push_back(parse_config(json_filename.c_str()));
        configs.back().Mapping.sampleId = sampleId++;
      }
    } else if (strcmp(args.argv[i], "-use-dcr") == 0) {
      use_dcr = true;
    }
  }
  return use_dcr;
}

static size_t config_cost(const Config &config)
{
  size_t grid_size =
    (size_t)config.Grid.xNum * config.Grid.yNum * config.Grid.zNum;
  size_t num_tiles =
    (size_t)config.Mapping.xTiles * config.Mapping.yTiles * config.Mapping.zTiles;
  return grid_size / num_tiles;
}

struct CmpConfig {
  bool operator() (const Config *lhs, const Config *rhs) {
    return config_cost(*lhs) > config_cost(*rhs);
  }
};

static void generate_static_mapping(const MachineModel &machine_model,
                                    const std::vector<Config> &configs,
                                    std::vector<SampleMapping> &mappings)
{
  std::vector<const Config*> p_configs;
  for (unsigned idx = 0; idx < configs.size(); ++idx)
    p_configs.push_back(&configs[idx]);
  std::sort(p_configs.begin(), p_configs.end(), CmpConfig());

  mappings.resize(p_configs.size());
  unsigned last_cpu_idx = 0;
  unsigned last_gpu_idx = 0;
  for (unsigned idx = 0; idx < p_configs.size(); ++idx)
  {
    const Config &config = *p_configs[idx];
    TaskMapping &mapping = mappings[config.Mapping.sampleId].task_mapping;
    mapping.num_procs = config.Mapping.xTiles *
                        config.Mapping.yTiles *
                        config.Mapping.zTiles;
    const char *target = (const char*)config.Mapping.target;
    bool target_cpu = strcmp(target, "CPU") == 0;
    bool target_gpu = !target_cpu && strcmp(target, "GPU") == 0;
    bool use_gpu = target_gpu && machine_model.gpus_list.size() > 0;

    if (!(target_cpu || target_gpu))
    {
      fprintf(stderr, "Unknown mapping target: %s\n", config.Mapping.target);
      assert(false);
    }
    else if (target_gpu && !use_gpu)
    {
      fprintf(stderr,
          "Sample %u requested GPU mapping, but no GPU is found."
          "Falling back to CPU mapping.\n", config.Mapping.sampleId);
    }
    if (!use_gpu && machine_model.ocpus_list.size() == 0)
    {
      fprintf(stderr, "OpenMP processors must exist. "
          "Please run with flag -ll:ocpu.\n");
      assert(false);
    }

    if (!use_gpu)
    {
      mapping.first_proc_idx = last_cpu_idx;
      last_cpu_idx =
        (last_cpu_idx + mapping.num_procs) % machine_model.ocpus_list.size();
      mapping.proc_kind = Processor::OMP_PROC;
    }
    else
    {
      mapping.first_proc_idx = last_gpu_idx;
      last_gpu_idx =
        (last_gpu_idx + mapping.num_procs) % machine_model.gpus_list.size();
      mapping.proc_kind = Processor::TOC_PROC;
    }

    mappings[config.Mapping.sampleId].shard_mapping =
      machine_model.generate_shard_mapping(mapping);
    fprintf(stderr, "Sample %u task first_proc %u task num_procs %u shard proc_kind %s "
        "shard fist_node %u shard num_nodes %u shard proc_kind %s\n",
        config.Mapping.sampleId,
        mapping.first_proc_idx, mapping.num_procs,
        mapping.proc_kind == Processor::OMP_PROC ? "cpu" : "gpu",
        mappings[config.Mapping.sampleId].shard_mapping.first_node_idx,
        mappings[config.Mapping.sampleId].shard_mapping.num_nodes,
        mappings[config.Mapping.sampleId].shard_mapping.proc_kind
        == Processor::IO_PROC ? "io" : "cpu");
  }
}

#ifdef CONTROL_REPLICATION
static void register_sharding_functors(HighLevelRuntime* runtime,
                                       const std::vector<Config> &configs,
                                       std::vector<SampleMapping> &mappings)
{
  const unsigned START_FUNCTOR_ID = 123;
  for (unsigned idx = 0; idx < configs.size(); ++idx)
  {
    unsigned functor_id = START_FUNCTOR_ID + idx;
    const Config &config = configs[idx];
    ShardMapping &shard_mapping = mappings[idx].shard_mapping;
    runtime->register_sharding_functor(functor_id,
        new TilingFunctor(config, shard_mapping.num_nodes));
    shard_mapping.sharding_functor_id = functor_id;
  }
}
#endif

static void create_mappers(Machine machine,
                           HighLevelRuntime* runtime,
                           const std::set<Processor>& local_procs) {

  MachineModel *machine_model = new MachineModel(machine);
  std::vector<Config> *configs = new std::vector<Config>();
  std::vector<SampleMapping> *mappings = new std::vector<SampleMapping>();

  bool use_dcr = parse_configurations(*configs);
  generate_static_mapping(*machine_model, *configs, *mappings);
#ifdef CONTROL_REPLICATION
  if (use_dcr) register_sharding_functors(runtime, *configs, *mappings);
#endif

  for (Processor proc : local_procs) {
    SoleilMapper* mapper =
      new SoleilMapper(runtime->get_mapper_runtime(), machine, proc,
                       use_dcr, machine_model, mappings);
    runtime->replace_default_mapper(mapper, proc);
  }
}

void register_mappers() {
  Runtime::add_registration_callback(create_mappers);
}
