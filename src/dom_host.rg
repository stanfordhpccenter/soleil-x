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

local task writeOutput(points : region(ispace(int3d), Point_columns),
                       Grid_xWidth : double,
                       Grid_yWidth : double,
                       Grid_zWidth : double,
                       Grid_xNum   : double,
                       Grid_yNum   : double,
                       Grid_zNum   : double)
where
  reads(points.{G})
do

  var dx = Grid_xWidth/Grid_xNum
  var dy = Grid_yWidth/Grid_yNum
  var dz = Grid_zWidth/Grid_zNum

  var fp = C.fopen ("volume_solution.txt", "w+");
  var limits = points.bounds
  C.fprintf(fp, "%d \t %d \t %d\n", limits.hi.x+1, limits.hi.y+1, limits.hi.z+1)
  C.fprintf(fp, "i \t j \t k \t x \t\t y \t\t z \t\t G\n")

  for ind in points.ispace do
    C.fprintf(fp, "%d \t %d \t %d \t %f \t %f \t %f \t %f\n",
                   ind.x, ind.y, ind.z,
                   (ind.x+0.5)*dx, (ind.y+0.5)*dy, (ind.z+0.5)*dz,
                   points[ind].G)
  end

  --for i = limits.lo.x, limits.hi.x+1 do
  --  for j = limits.lo.y, limits.hi.y+1 do
  --    for k = limits.lo.z, limits.hi.z+1 do
  --      C.fprintf(fp, "%d \t %d \t %d \t %f \t %f \t %f \t %f\n",
  --                     i, j, k,
  --                     (i+0.5)*dx, (j+0.5)*dy, (k+0.5)*dz,
  --                     points[{i,j,k}].G)
  --    end
  --  end
  --end

  C.fclose(fp);

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
    (points, tiles, 0, 0, 0);
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
  writeOutput(points, 
              config.Grid.xWidth, config.Grid.yWidth, config.Grid.zWidth,
              config.Radiation.u.DOM.xNum, config.Radiation.u.DOM.yNum, config.Radiation.u.DOM.zNum)
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
