import "regent"

local SCHEMA = terralib.includec("schema.h")

local Config = SCHEMA.Config

task Probe_WriteHeader()
  regentlib.c.printf('hello\n')
end

task workSingle(config : Config)
  for i = 0,config.IO.probes.length do
    Probe_WriteHeader()
  end
end

task main()
  var config : Config[1]
  config[0].IO.probes.length = 0
  workSingle(config[0])
end

regentlib.saveobj(main, "soleil.o", "object")
