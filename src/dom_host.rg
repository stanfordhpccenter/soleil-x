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
  p_to_s3d_1 : int3d;
  p_to_s3d_2 : int3d;
  p_to_s3d_3 : int3d;
  p_to_s3d_4 : int3d;
  p_to_s3d_5 : int3d;
  p_to_s3d_6 : int3d;
  p_to_s3d_7 : int3d;
  p_to_s3d_8 : int3d;
  s3d_to_p_1 : int3d;
  s3d_to_p_2 : int3d;
  s3d_to_p_3 : int3d;
  s3d_to_p_4 : int3d;
  s3d_to_p_5 : int3d;
  s3d_to_p_6 : int3d;
  s3d_to_p_7 : int3d;
  s3d_to_p_8 : int3d;
}

-------------------------------------------------------------------------------
-- Import DOM module
-------------------------------------------------------------------------------

local DOM = (require 'dom-desugared')(MAX_ANGLES_PER_QUAD, Point_columns, SCHEMA.Config)
local DOM_INST = DOM.mkInstance()

-------------------------------------------------------------------------------
-- Proxy tasks
-------------------------------------------------------------------------------

local task writeIntensity(points : region(ispace(int3d), Point_columns))
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

local __forbid(__optimize)
task main()
  -- Read configuration
  var args = C.legion_runtime_get_input_args()
  var stderr = C.fdopen(2, 'w')
  if args.argc < 2 then
    C.fprintf(stderr, "Usage: %s config.json\n", args.argv[0])
    C.fflush(stderr)
    C.exit(1)
  end
  var config : SCHEMA.Config[1]
  SCHEMA.parse_Config([&SCHEMA.Config](config), args.argv[1])
  -- Declare externally-managed regions
  var is_points = ispace(int3d, {config[0].Radiation.xNum,
                                 config[0].Radiation.yNum,
                                 config[0].Radiation.zNum})
  var points = region(is_points, Point_columns)
  var tiles = ispace(int3d, {config[0].Mapping.tiles[0],
                             config[0].Mapping.tiles[1],
                             config[0].Mapping.tiles[2]})
  var p_points =
    [UTIL.mkPartitionEqually(int3d, int3d, Point_columns)](points, tiles);
  -- Declare DOM-managed regions
  [DOM_INST.DeclSymbols(rexpr config[0] end, tiles)];
  [DOM_INST.InitRegions(rexpr config[0] end, tiles, p_points)];
  -- Prepare fake inputs
  fill(points.Ib, (SB/PI) * pow(1000.0,4.0))
  fill(points.sigma, 5.0);
  -- Invoke DOM solver
  [DOM_INST.ComputeRadiationField(rexpr config[0] end, tiles, p_points)];
  -- Output results
  writeIntensity(points)
end

regentlib.saveobj(main, 'dom_host.o', 'object')
