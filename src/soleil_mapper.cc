#include <string.h>

#include "mappers/default_mapper.h"
#include "realm/logging.h"

#include "config_schema.h"
#include "soleil_mapper.h"

using namespace Legion;
using namespace Legion::Mapping;

Realm::Logger mapper_log("soleil_mapper");

class SoleilMapper : public DefaultMapper {
public:
  SoleilMapper(MapperRuntime* rt,
               Machine machine,
               Processor local)
    : DefaultMapper(rt, machine, local, "soleil_mapper") {
    mapper_log.spew("Soleil mapper created");
  }
};

static void create_mappers(Machine machine,
                           HighLevelRuntime* runtime,
                           const std::set<Processor>& local_procs) {
  for (Processor proc : local_procs) {
    SoleilMapper* mapper = new SoleilMapper(runtime->get_mapper_runtime(),
                                            machine,
                                            proc);
    runtime->replace_default_mapper(mapper, proc);
  }
}

void register_mappers() {
  HighLevelRuntime::add_registration_callback(create_mappers);
}
