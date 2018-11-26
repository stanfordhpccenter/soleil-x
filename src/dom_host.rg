-- Runs dom.rg standalone.
-- Reads configuration options in the same format as main simulation.
-- Uses default values for Ib and sigma.

-------------------------------------------------------------------------------
-- Imports
-------------------------------------------------------------------------------

import 'regent'

local C = regentlib.c
local SCHEMA = terralib.includec("config_schema.h")
local UTIL = require 'util'

local pow = regentlib.pow(double)

-------------------------------------------------------------------------------
-- Compile-time constants
-------------------------------------------------------------------------------

local PI = 3.1415926535898
local SB = 5.67e-8

local MAX_ANGLES_PER_QUAD = 44

-------------------------------------------------------------------------------
-- Proxy radiation grid
-------------------------------------------------------------------------------

struct Point_columns {
  G : double;
  S : double;
  Ib : double;
  sigma : double;
}

-------------------------------------------------------------------------------
-- Import DOM module
-------------------------------------------------------------------------------

local DOM = (require 'dom-desugared')(MAX_ANGLES_PER_QUAD, Point_columns, SCHEMA)
local DOM_INST = DOM.mkInstance()

-------------------------------------------------------------------------------
-- Proxy tasks
-------------------------------------------------------------------------------

local __demand(__leaf) -- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task writeIntensity(points : region(ispace(int3d), Point_columns))
where
  reads(points.G)
do
  var limits = points.bounds
  var f = UTIL.openFile("intensity.dat", "w")
  for i = limits.lo.x, limits.hi.x+1 do
    for j = limits.lo.y, limits.hi.y+1 do
      for k = limits.lo.z, limits.hi.z+1 do
        C.fprintf(f,' %.15e \n', points[{i,j,k}].G)
      end
      C.fprintf(f,'\n')
    end
    C.fprintf(f,'\n')
  end
  C.fclose(f)
end

-------------------------------------------------------------------------------
-- Proxy main
-------------------------------------------------------------------------------

local __forbid(__optimize) __demand(__inner, __replicable)
task work(config : SCHEMA.Config)
  -- Declare externally-managed regions
  var is_points = ispace(int3d, {config.Radiation.u.DOM.xNum,
                                 config.Radiation.u.DOM.yNum,
                                 config.Radiation.u.DOM.zNum})
  var points = region(is_points, Point_columns)
  var tiles = ispace(int3d, {config.Mapping.tiles[0],
                             config.Mapping.tiles[1],
                             config.Mapping.tiles[2]})
  var p_points =
    [UTIL.mkPartitionEqually(int3d, int3d, Point_columns)]
    (points, tiles, int3d{0,0,0}, int3d{0,0,0});
  -- Declare DOM-managed regions
  [DOM_INST.DeclSymbols(config, tiles)];
  [DOM_INST.InitRegions(config, tiles, p_points)];
  -- Prepare fake inputs
  fill(points.Ib, (SB/PI) * pow(1000.0,4.0))
  fill(points.sigma, 5.0);
  -- Invoke DOM solver
  [DOM_INST.ComputeRadiationField(config, tiles, p_points)];
  -- Output results
  writeIntensity(points)
end

local __demand(__inner)
task main()
  var args = C.legion_runtime_get_input_args()
  var stderr = C.fdopen(2, 'w')
  if args.argc < 2 then
    C.fprintf(stderr, "Usage: %s config.json\n", args.argv[0])
    C.fflush(stderr)
    C.exit(1)
  end
  var config : SCHEMA.Config[1]
  SCHEMA.parse_Config([&SCHEMA.Config](config), args.argv[1])
  regentlib.assert(config[0].Radiation.type == SCHEMA.RadiationModel_DOM,
                   'Configuration file must use DOM radiation model')
  work(config[0])
end

regentlib.saveobj(main, 'dom_host.o', 'object')
