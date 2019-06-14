-- Runs dom.rg standalone.
-- Reads configuration options in the same format as main simulation.
-- Uses default values for Ib and sigma.

-------------------------------------------------------------------------------
-- Imports
-------------------------------------------------------------------------------

import 'regent'

local C = regentlib.c
local SCHEMA = terralib.includec("config_schema.h")
local UTIL = require 'util-desugared'

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
  centerCoordinates : double[3];
  cellWidth : double[3];
  G : double;
  S : double;
  Ib : double;
  sigma : double;
}

-------------------------------------------------------------------------------
-- Mesh routines
-------------------------------------------------------------------------------

__demand(__inline)
task uniform_cell_width(x_min : double,
                        x_max : double,
                        Nx    : uint64) : double
  return (x_max-x_min)/Nx
end

__demand(__inline)
task uniform_cell_center(x_min : double,
                         x_max : double,
                         Nx    : uint64,
                         i     : uint64) : double
  var dx = uniform_cell_width(x_min, x_max, Nx)
  return x_min + i*dx + dx/2.0
end

-------------------------------------------------------------------------------
-- Import DOM module
-------------------------------------------------------------------------------

local DOM = (require 'dom-desugared')(MAX_ANGLES_PER_QUAD, Point_columns, SCHEMA)
local DOM_INST = DOM.mkInstance()

-------------------------------------------------------------------------------
-- Proxy tasks
-------------------------------------------------------------------------------

__demand(__leaf) -- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task Radiation_InitializeGeometry(Radiation : region(ispace(int3d), Point_columns),
                                  Grid_xNum : int32, Grid_xOrigin : double, Grid_xWidth : double,
                                  Grid_yNum : int32, Grid_yOrigin : double, Grid_yWidth : double,
                                  Grid_zNum : int32, Grid_zOrigin : double, Grid_zWidth : double)
where
  reads writes(Radiation.{centerCoordinates, cellWidth})
do
  for cell in Radiation do
    var x_neg_boundary = Grid_xOrigin
    var x_pos_boundary = Grid_xOrigin + Grid_xWidth
    var x_idx_interior = cell.x
    cell.centerCoordinates[0] = uniform_cell_center(x_neg_boundary, x_pos_boundary, Grid_xNum, x_idx_interior)
    cell.cellWidth[0]         = uniform_cell_width( x_neg_boundary, x_pos_boundary, Grid_xNum)

    var y_neg_boundary = Grid_yOrigin
    var y_pos_boundary = Grid_yOrigin + Grid_yWidth
    var y_idx_interior = cell.y
    cell.centerCoordinates[1] = uniform_cell_center(y_neg_boundary, y_pos_boundary, Grid_yNum, y_idx_interior)
    cell.cellWidth[1]         = uniform_cell_width( y_neg_boundary, y_pos_boundary, Grid_yNum)

    var z_neg_boundary = Grid_zOrigin
    var z_pos_boundary = Grid_zOrigin + Grid_zWidth
    var z_idx_interior = cell.z
    cell.centerCoordinates[2] = uniform_cell_center(z_neg_boundary, z_pos_boundary, Grid_zNum, z_idx_interior)
    cell.cellWidth[2]         = uniform_cell_width( z_neg_boundary, z_pos_boundary, Grid_zNum)
  end
end

__demand(__leaf) -- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
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

__forbid(__optimize) __demand(__inner, __replicable)
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
    [UTIL.mkPartitionByTile(int3d, int3d, Point_columns)]
    (points, tiles, int3d{0,0,0}, int3d{0,0,0});
  -- Declare DOM-managed regions
  [DOM_INST.DeclSymbols(config, tiles)];
  [DOM_INST.InitRegions(config, tiles, p_points)];
  for c in tiles do
    Radiation_InitializeGeometry(p_points[c],
                                 config.Radiation.u.DOM.xNum, config.Grid.origin[0], config.Grid.xWidth,
                                 config.Radiation.u.DOM.yNum, config.Grid.origin[1], config.Grid.yWidth,
                                 config.Radiation.u.DOM.zNum, config.Grid.origin[2], config.Grid.zWidth)
  end
  -- Prepare fake inputs
  fill(points.Ib, (SB/PI) * pow(1000.0,4.0))
  fill(points.sigma, config.Radiation.u.DOM.qa + config.Radiation.u.DOM.qs);
  -- Invoke DOM solver
  [DOM_INST.ComputeRadiationField(config, tiles, p_points)];
  -- Output results
  writeIntensity(points)
end

__demand(__inner)
task main()
  var args = C.legion_runtime_get_input_args()
  var stderr = C.fdopen(2, 'w')
  if args.argc < 2 then
    C.fprintf(stderr, "Usage: %s config.json\n", args.argv[0])
    C.fflush(stderr)
    C.exit(1)
  end
  var config : SCHEMA.Config
  SCHEMA.parse_Config(&config, args.argv[1])
  regentlib.assert(config.Radiation.type == SCHEMA.RadiationModel_DOM,
                   'Configuration file must use DOM radiation model')
  work(config)
end

regentlib.saveobj(main, 'dom_host.o', 'object')
