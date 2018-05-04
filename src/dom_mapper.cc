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
            const Config &c)
    : DefaultMapper(rt, machine, local, "dom_mapper"), config(c) {
  }

  virtual void slice_task(const MapperContext ctx,
                          const Task& task,
                          const SliceTaskInput& input,
                          SliceTaskOutput& output) {
    if (input.domain.get_dim() == 3)
    {
      DomainT<3,coord_t> point_space = input.domain;
      Point<3,coord_t> num_blocks =
        default_select_num_blocks<3>(local_cpus.size(), point_space.bounds);
      default_decompose_points<3>(point_space, local_cpus,
          num_blocks, false/*recurse*/,
          false, output.slices);
      for (unsigned idx = 0; idx < output.slices.size(); ++idx)
      {
        TaskSlice &slice = output.slices[idx];
        DomainPoint p = slice.domain.lo();
        coord_t linearized = p[0] +
                             p[1] * config.Mapping.xTiles +
                             p[2] * config.Mapping.xTiles * config.Mapping.yTiles;
        slice.proc = local_cpus[linearized % local_cpus.size()];
      }
    }
    else
      DefaultMapper::slice_task(ctx, task, input, output);
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

  for (Processor proc : local_procs) {
    DOMMapper* mapper =
      new DOMMapper(runtime->get_mapper_runtime(), machine, proc, config);
    runtime->replace_default_mapper(mapper, proc);
  }

  runtime->register_sharding_functor(DOM_SHARDING_FUNCTOR_ID,
      new GlobalTilingFunctor(config));
}

void register_mappers() {
  Runtime::add_registration_callback(create_mappers);
}
