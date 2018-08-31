#include "ProblemDescDB.hpp"
#include "LibraryEnvironment.hpp"
#include "DakotaModel.hpp"
#include "DakotaInterface.hpp"
#include "PluginSerialDirectApplicInterface.hpp"

#include "legion.h"

#include "soleil.h"




using namespace Legion;



void run_dakota_parse(const char* dakota_input_file) {
  Dakota::ProgramOptions opts;
  opts.input_file(dakota_input_file);
  Dakota::LibraryEnvironment env(opts);
  serial_interface_plugin(env);
  env.execute();
}

/** client code doesn't require access to detailed Dakota data (such as Model-based parallel
    configuration information) to construct the DirectApplicInterface.
    This example plugs-in a derived serial direct application
    interface instance ("plugin_rosenbrock"). */
void serial_interface_plugin(Dakota::LibraryEnvironment& env) {
  std::string model_type(""); // demo: empty string will match any model type
  std::string interf_type("direct");
  std::string an_driver("soleil");



  Dakota::ProblemDescDB& problem_db = env.problem_description_db();
  Dakota::Interface* serial_iface =
    new SIM::SerialDirectApplicInterface(problem_db);
  // pass addditional stuff to cache here



  bool plugged_in =
    env.plugin_interface(model_type, interf_type, an_driver, serial_iface);

  if (!plugged_in) {
    Cerr << "Error: no serial interface plugin performed.  Check "
	 << "compatibility between parallel\n       configuration and "
	 << "selected analysis_driver." << std::endl;
    Dakota::abort_handler(-1);
  }
}

















task my_regent_task(x : int)
  regentlib.c.printf("Hello from Regent! (value %d)\n", x)
end


// do these in main

terra initSample(config : &Config, num : int, outDirBase : &int8)
  config.Mapping.sampleId = num
  C.snprintf(config.Mapping.outDir, 256, "%s/sample%d", outDirBase, num)
  UTIL.createDir(config.Mapping.outDir)
end

__demand(__inner)
task main()
  var args = regentlib.c.legion_runtime_get_input_args()
  var outDirBase = '.'
  for i = 1, args.argc do
    if C.strcmp(args.argv[i], '-o') == 0 and i < args.argc-1 then
      outDirBase = args.argv[i+1]
    end
  end
  var launched = 0
  for i = 1, args.argc do
    if C.strcmp(args.argv[i], '-i') == 0 and i < args.argc-1 then
      var config : Config[1]
      SCHEMA.parse_Config([&Config](config), args.argv[i+1])
      initSample([&Config](config), launched, outDirBase)
      launched += 1
      workSingle(config[0])
    elseif C.strcmp(args.argv[i], '-m') == 0 and i < args.argc-1 then
      var mc : MultiConfig[1]
      SCHEMA.parse_MultiConfig([&MultiConfig](mc), args.argv[i+1])
      initSample([&Config](mc[0].configs), launched, outDirBase)
      initSample([&Config](mc[0].configs) + 1, launched + 1, outDirBase)
      launched += 2
      workDual(mc[0])
    end
  end
  if launched < 1 then
    var stderr = C.fdopen(2, 'w')
    C.fprintf(stderr, "No testcases supplied.\n")
    C.fflush(stderr)
    C.exit(1)
  end
end









enum {
  TOP_LEVEL_TASK_ID,
};


void top_level_task(const Task *task,
                    const std::vector<PhysicalRegion> &regions,
                    Context ctx, Runtime *runtime) {
  my_regent_task_launcher launcher;
  launcher.add_argument_r(region, region, fields);
  launcher.add_argument_x(12345);
  launcher.execute(runtime, ctx);
}




int main(int argc, char **argv) {
  Runtime::set_top_level_task_id(TOP_LEVEL_TASK_ID);
  TaskVariantRegistrar registrar(TOP_LEVEL_TASK_ID, "top_level");
  registrar.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
  Runtime::preregister_task_variant<top_level_task>(registrar, "top_level");


  Runtime::add_registration_callback(update_mappers);
  soleil_h_register();
  return Runtime::start(argc, argv);
}
