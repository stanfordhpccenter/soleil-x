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

class DOMMapper : public DefaultMapper {
public:
  DOMMapper(MapperRuntime* rt,
            Machine machine,
            Processor local,
            std::vector<Processor>* pl,
            const Config &c)
    : DefaultMapper(rt, machine, local, "dom_mapper"),
      procs_list(*pl), config(c) {
  }

  virtual Processor default_policy_select_initial_processor(
                              MapperContext ctx,
                              const Task& task) {
    if (task.regions.size() == 0)
      return task.orig_proc;
    else
      return DefaultMapper::default_policy_select_initial_processor(ctx, task);
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
        coord_t linearized = p[0] +
                             p[1] * config.Mapping.xTiles +
                             p[2] * config.Mapping.xTiles * config.Mapping.yTiles;
        output.slices.emplace_back(Domain(itr.p, itr.p),
                                   procs_list[linearized % procs_list.size()],
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

  virtual void default_policy_select_target_processors(
                                    MapperContext ctx,
                                    const Task &task,
                                    std::vector<Processor> &target_procs)
  {
    target_procs.push_back(task.target_proc);
  }

private:
  const std::vector<Processor> &procs_list;
  Config config;
};

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
      new DOMMapper(runtime->get_mapper_runtime(), machine, proc,
                    procs_list, config);
    runtime->replace_default_mapper(mapper, proc);
  }
}

void register_mappers() {
  Runtime::add_registration_callback(create_mappers);
}
