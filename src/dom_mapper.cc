#include <iostream>
#include <regex>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "mappers/default_mapper.h"
#include "realm/logging.h"

#include "config_schema.h"
#include "dom_mapper.h"

using namespace Legion;
using namespace Legion::Mapping;

//=============================================================================

static Realm::Logger LOG("dom_mapper");

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

enum
{
  DOM_SHARDING_FUNCTOR_ID = 1,
};

class DOMMapper : public DefaultMapper {
public:
  DOMMapper(MapperRuntime* rt,
            Machine machine,
            Processor local,
            std::vector<Processor>* pl,
            const Config &c)
    : DefaultMapper(rt, machine, local, "dom_mapper"), procs_list(*pl), config(c) {
    num_nodes = procs_list.size() / local_cpus.size();
    assert(procs_list.size() % local_cpus.size() == 0);
    num_points_per_shard = config.Mapping.xTiles *
                           config.Mapping.yTiles *
                           config.Mapping.zTiles;
    num_points_per_shard =
      (num_points_per_shard + num_nodes - 1) / num_nodes;
  }

  virtual void slice_task(const MapperContext ctx,
                          const Task& task,
                          const SliceTaskInput& input,
                          SliceTaskOutput& output) {
    if (input.domain.get_dim() == 3)
    {
      for (Domain::DomainPointIterator itr(input.domain); itr; itr++)
      {
        DomainPoint p = itr.p;
        p[0] = p[0] < config.Mapping.xTiles ? p[0] : config.Mapping.xTiles - 1;
        p[1] = p[1] < config.Mapping.yTiles ? p[1] : config.Mapping.yTiles - 1;
        p[2] = p[2] < config.Mapping.zTiles ? p[2] : config.Mapping.zTiles - 1;
        size_t linearized = p[0] +
                            p[1] * config.Mapping.xTiles +
                            p[2] * config.Mapping.xTiles * config.Mapping.yTiles;
        linearized %= num_points_per_shard;
        assert(0 <= linearized && linearized < local_cpus.size());
        output.slices.emplace_back(Domain(itr.p, itr.p), local_cpus[linearized],
                                   false, false);
      }
    }
    else
      DefaultMapper::slice_task(ctx, task, input, output);
  }

  virtual LogicalRegion default_policy_select_instance_region(
                                MapperContext ctx, Memory target_memory,
                                const RegionRequirement &req,
                                const LayoutConstraintSet &constraints,
                                bool force_new_instances, 
                                bool meets_constraints)
  {
    return req.region;
  }

  virtual void select_sharding_functor(
                                 const MapperContext                ctx,
                                 const Task&                        task,
                                 const SelectShardingFunctorInput&  input,
                                       SelectShardingFunctorOutput& output)
  {
    if (task.is_index_space && task.index_domain.dim == 3)
      output.chosen_functor = DOM_SHARDING_FUNCTOR_ID;
    else
      DefaultMapper::select_sharding_functor(ctx, task, input, output);
  }

  virtual void default_policy_select_target_processors(
                                    MapperContext ctx,
                                    const Task &task,
                                    std::vector<Processor> &target_procs)
  {
    target_procs.push_back(task.target_proc);
  }
private:
  const std::vector<Processor> &procs_list;
  size_t num_nodes;
  size_t num_points_per_shard;
  Config config;
};

//=============================================================================

class GlobalTilingFunctor : public ShardingFunctor {
  public:
    GlobalTilingFunctor(const Config &config);
    virtual ~GlobalTilingFunctor(void);
  public:
    virtual ShardID shard(const DomainPoint &point,
                          const Domain &full_space,
                          const size_t total_shards);
  protected:
    Config config;
};

GlobalTilingFunctor::GlobalTilingFunctor(const Config &_config)
  : ShardingFunctor(), config(_config)
{
}

GlobalTilingFunctor::~GlobalTilingFunctor()
{
}

ShardID GlobalTilingFunctor::shard(const DomainPoint &point,
                                   const Domain &full_space,
                                   const size_t total_shards)
{
  assert(point.dim == 3);
  DomainPoint p = point;
  p[0] = p[0] < config.Mapping.xTiles ? p[0] : config.Mapping.xTiles - 1;
  p[1] = p[1] < config.Mapping.yTiles ? p[1] : config.Mapping.yTiles - 1;
  p[2] = p[2] < config.Mapping.zTiles ? p[2] : config.Mapping.zTiles - 1;
  coord_t num_points_per_shard = config.Mapping.xTiles *
                                 config.Mapping.yTiles *
                                 config.Mapping.zTiles;
  num_points_per_shard =
    (num_points_per_shard + total_shards - 1) / total_shards;
  coord_t linearized = p[0] +
                       p[1] * config.Mapping.xTiles +
                       p[2] * config.Mapping.xTiles * config.Mapping.yTiles;
  ShardID shard_id = (ShardID)(linearized / num_points_per_shard);
  assert(shard_id < total_shards);
  return shard_id;
}

//=============================================================================

static void create_mappers(Machine machine,
                           HighLevelRuntime* runtime,
                           const std::set<Processor>& local_procs) {
  InputArgs args = Runtime::get_input_args();
  assert(args.argc >= 2);
  Config config = parse_config(args.argv[1]);

  std::vector<Processor>* procs_list = new std::vector<Processor>();

  Machine::ProcessorQuery procs_query(machine);
  procs_query.only_kind(Processor::LOC_PROC);
  for (Machine::ProcessorQuery::iterator it = procs_query.begin();
        it != procs_query.end(); it++)
    procs_list->push_back(*it);

  for (Processor proc : local_procs) {
    DOMMapper* mapper =
      new DOMMapper(runtime->get_mapper_runtime(), machine, proc, procs_list, config);
    runtime->replace_default_mapper(mapper, proc);
  }

  runtime->register_sharding_functor(DOM_SHARDING_FUNCTOR_ID,
      new GlobalTilingFunctor(config));
}

void register_mappers() {
  Runtime::add_registration_callback(create_mappers);
}
