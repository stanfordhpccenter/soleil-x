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
local cos = regentlib.cos(double)

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
-- MESH ROUTINES
-------------------------------------------------------------------------------
-- Description:
--     Linear interpolation, given the line defined by the points 
--     (x=alpha, y=a) and (x=beta, y=b) find the y location of the
--     point on the line (x=xi, y=?) 
-- Input:
--     xi = location on x axis
--     alpha = lower point on x axis
--     beta =  upper point on x axis 
--     a = lower point on y axis 
--     b = upper point on y axis 
-- Output:
--     y location on line at x=xi
local __demand(__inline)
task linear_interpolation(xi : double,
                          alpha : double,
                          beta  : double,
                          a     : double,
                          b     : double) : double
  return (b-a)/(beta-alpha)*(xi-alpha) + a
end


-- Description:
--     Generate the cell width of a nonuniform mesh
-- Input:
--     x_min = domain minimum
--     x_max = domain maximum 
--     Nx = number of cells between x_min and x_max
-- Output:
--     width of the non-uniform mesh cell
local __demand(__inline)
task uniform_cell_width(x_min : double,
                        x_max : double,
                        Nx    : uint64) : double
  return (x_max-x_min)/Nx
end


-- Description:
--     Generate the cell center on a uniform mesh
-- Input:
--     x_min = domain minimum
--     x_max = domain maximum 
--     Nx = number of cells between x_min and x_max
--     i  = cell index between x_min and x_max
--         Note: i = 0 has x_min as left face and 
--               i = Nx-1 has x_max as right face
--               so no ghost cells accounted for here
-- Output:
--     location of cell center
local __demand(__inline)
task uniform_cell_center(x_min : double,
                         x_max : double,
                         Nx    : uint64,
                         i     : uint64) : double
  var dx = uniform_cell_width(x_min, x_max, Nx)
  return x_min + i*dx + dx/2.0
end

-- Description:
--     Generate the location of the face in the negative direction
-- Input:
--     x_min = domain minimum
--     x_max = domain maximum 
--     Nx = number of cells between x_min and x_max
--     i  = cell index between x_min and x_max
--         Note: i = 0 has x_min as negative direction (left) face and 
--               i = Nx-1 has x_max as positive direction (right) face
--               so no ghost cells accounted for here
-- Output:
--     location of face in the negative direction
local __demand(__inline)
task uniform_cell_neg_face(x_min : double,
                           x_max : double,
                           Nx    : uint64,
                           i     : uint64) : double
  var dx = uniform_cell_width(x_min, x_max, Nx)
  var x_center = uniform_cell_center(x_min, x_max, Nx, i)
  return x_center - 0.5*dx
end

-- Description:
--     Generate the location of the face in the postive direction
-- Input:
--     x_min = domain minimum
--     x_max = domain maximum 
--     Nx = number of cells between x_min and x_max
--     i  = cell index between x_min and x_max
--         Note: i = 0 has x_min as negative direction (left) face and 
--               i = Nx-1 has x_max as positive direction (right) face
--               so no ghost cells accounted for here
-- Output:
--     location of face in the postive direction
local __demand(__inline)
task uniform_cell_pos_face(x_min : double,
                           x_max : double,
                           Nx    : uint64,
                           i     : uint64) : double
  var dx = uniform_cell_width(x_min, x_max, Nx)
  var x_center = uniform_cell_center(x_min, x_max, Nx, i)
  return x_center + 0.5*dx
end


-- Description:
--     non-linear map point (x) on the interval (x_min, x_max) using
--     a cosine  
-- Input:
--     x = location on uniform mesh
--     x_min = domain minimum
--     x_max = domain maximum 
-- Output:
--     x location on a non-uniform mesh
local __demand(__inline)
task transform_uniform_to_nonuniform(x : double,
                                     x_min : double,
                                     x_max : double) : double
  -- map x onto the interval -1 to 1
  var x_scaled_minus1_to_plus1 = linear_interpolation(x, x_min, x_max, -1.0, 1.0)

  -- map non-uniformly onto the interval -1 to 1
  var x_non_uniform_minus1_to_plus1 = -1.0*cos(PI*(x_scaled_minus1_to_plus1+1.0)/2.0)

  -- map non-uniform sample back to origional interval x_min to x_max
  return  linear_interpolation(x_non_uniform_minus1_to_plus1, -1.0, 1.0, x_min, x_max)
end

-- Description:
--     Generate the location of the face in the negative direction
--     on a non uniform mesh
-- Input:
--     x_min = domain minimum
--     x_max = domain maximum 
--     Nx = number of cells between x_min and x_max
--     i  = cell index between x_min and x_max
--         Note: i = 0 has x_min as negative direction (left) face and 
--               i = Nx-1 has x_max as positive direction (right) face
--               so no ghost cells accounted for here
-- Output:
--     location of face in the negative direction
local __demand(__inline)
task nonuniform_cell_neg_face(x_min : double,
                              x_max : double,
                              Nx    : uint64,
                              i     : uint64) : double
  var x_uniform_neg_face = uniform_cell_neg_face(x_min, x_max, Nx, i)
  return transform_uniform_to_nonuniform(x_uniform_neg_face, x_min, x_max)
end

-- Description:
--     Generate the location of the face in the postive direction
--     on a non uniform mesh
-- Input:
--     x_min = domain minimum
--     x_max = domain maximum 
--     Nx = number of cells between x_min and x_max
--     i  = cell index between x_min and x_max
--         Note: i = 0 has x_min as negative direction (left) face and 
--               i = Nx-1 has x_max as positive direction (right) face
--               so no ghost cells accounted for here
-- Output:
--     location of face in the postive direction
local __demand(__inline)
task nonuniform_cell_pos_face(x_min : double,
                              x_max : double,
                              Nx    : uint64,
                              i     : uint64) : double
  var x_uniform_pos_face = uniform_cell_pos_face(x_min, x_max, Nx, i)
  return transform_uniform_to_nonuniform(x_uniform_pos_face, x_min, x_max)
end


-- Description:
--     Generate the cell center of a nonuniform mesh
-- Input:
--     x_min = domain minimum
--     x_max = domain maximum 
--     Nx = number of cells between x_min and x_max
--     i  = cell index between x_min and x_max
--         Note: i = 0 has x_min as left face and 
--               i = Nx-1 has x_max as right face
--               so no ghost cells accounted for here
-- Output:
--     x location on a non-uniform mesh
local __demand(__inline)
task nonuniform_cell_center(x_min : double,
                            x_max : double,
                            Nx    : uint64,
                            i     : uint64) : double
  var x_non_uniform_neg_face = nonuniform_cell_neg_face(x_min, x_max, Nx, i)
  var x_non_uniform_pos_face = nonuniform_cell_pos_face(x_min, x_max, Nx, i)

  var x_non_uniform_center = 0.5*(x_non_uniform_neg_face + x_non_uniform_pos_face)

  return x_non_uniform_center 
end


-- Description:
--     Generate the cell width of a nonuniform mesh
-- Input:
--     x_min = domain minimum
--     x_max = domain maximum 
--     Nx = number of cells between x_min and x_max
--     i  = cell index between x_min and x_max
--         Note: i = 0 has x_min as left face and 
--               i = Nx-1 has x_max as right face
--               so no ghost cells accounted for here
-- Output:
--     width of the non-uniform mesh cell
local __demand(__inline)
task nonuniform_cell_width(x_min : double,
                           x_max : double,
                           Nx    : uint64,
                           i     : uint64) : double
  var x_non_uniform_neg_face = nonuniform_cell_neg_face(x_min, x_max, Nx, i)
  var x_non_uniform_pos_face = nonuniform_cell_pos_face(x_min, x_max, Nx, i)

  var x_non_uniform_dx = x_non_uniform_pos_face - x_non_uniform_neg_face

  return x_non_uniform_dx 
end
-------------------------------------------------------------------------------
-- Import DOM module
-------------------------------------------------------------------------------

local DOM = (require 'dom-desugared')(MAX_ANGLES_PER_QUAD, Point_columns, SCHEMA)
local DOM_INST = DOM.mkInstance()

-------------------------------------------------------------------------------
-- Proxy tasks
-------------------------------------------------------------------------------


local task Radiation_InitializeGeometry(Radiation : region(ispace(int3d), Point_columns),
                                        Grid_xType : SCHEMA.GridType, Grid_yType : SCHEMA.GridType, Grid_zType : SCHEMA.GridType,
                                        Grid_xNum : int32, Grid_xOrigin : double, Grid_xWidth : double,
                                        Grid_yNum : int32, Grid_yOrigin : double, Grid_yWidth : double,
                                        Grid_zNum : int32, Grid_zOrigin : double, Grid_zWidth : double)
where
  reads writes(Radiation.{centerCoordinates, cellWidth})
do
  -- Find cell center coordinates and cell width of interior cells
  for cell in Radiation do

      var x_neg_boundary = Grid_xOrigin
      var x_pos_boundary = Grid_xOrigin + Grid_xWidth
      var x_idx_interior = cell.x 
      if (Grid_xType == SCHEMA.GridType_Stretched) then
        cell.centerCoordinates[0] = nonuniform_cell_center(x_neg_boundary, x_pos_boundary, Grid_xNum, x_idx_interior)
        cell.cellWidth[0]         = nonuniform_cell_width( x_neg_boundary, x_pos_boundary, Grid_xNum, x_idx_interior)
      elseif (Grid_xType == SCHEMA.GridType_Uniform) then
        cell.centerCoordinates[0] = uniform_cell_center(x_neg_boundary, x_pos_boundary, Grid_xNum, x_idx_interior)
        cell.cellWidth[0]         = uniform_cell_width( x_neg_boundary, x_pos_boundary, Grid_xNum)
      end

      var y_neg_boundary = Grid_yOrigin
      var y_pos_boundary = Grid_yOrigin + Grid_yWidth
      var y_idx_interior = cell.y
      if (Grid_yType == SCHEMA.GridType_Stretched) then
        cell.centerCoordinates[1] = nonuniform_cell_center(y_neg_boundary, y_pos_boundary, Grid_yNum, y_idx_interior)
        cell.cellWidth[1]         = nonuniform_cell_width( y_neg_boundary, y_pos_boundary, Grid_yNum, y_idx_interior)
      elseif (Grid_yType == SCHEMA.GridType_Uniform) then
        cell.centerCoordinates[1] = uniform_cell_center(y_neg_boundary, y_pos_boundary, Grid_yNum, y_idx_interior)
        cell.cellWidth[1]         = uniform_cell_width( y_neg_boundary, y_pos_boundary, Grid_yNum)
      end

      var z_neg_boundary = Grid_zOrigin
      var z_pos_boundary = Grid_zOrigin + Grid_zWidth
      var z_idx_interior = cell.z
      if (Grid_zType == SCHEMA.GridType_Stretched) then
        cell.centerCoordinates[2] = nonuniform_cell_center(z_neg_boundary, z_pos_boundary, Grid_zNum, z_idx_interior)
        cell.cellWidth[2]         = nonuniform_cell_width( z_neg_boundary, z_pos_boundary, Grid_zNum, z_idx_interior)
      elseif (Grid_zType == SCHEMA.GridType_Uniform) then
        cell.centerCoordinates[2] = uniform_cell_center(z_neg_boundary, z_pos_boundary, Grid_zNum, z_idx_interior)
        cell.cellWidth[2]         = uniform_cell_width( z_neg_boundary, z_pos_boundary, Grid_zNum)
      end

  end

end

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

local task writeOutput(points : region(ispace(int3d), Point_columns))
where
  reads(points.{centerCoordinates, cellWidth, G})
do

  var fp = C.fopen ("volume_solution.txt", "w+");
  var limits = points.bounds
  C.fprintf(fp, "%d \t %d \t %d\n", limits.hi.x+1, limits.hi.y+1, limits.hi.z+1)
  C.fprintf(fp, "i \t j \t k \t x \t\t y \t\t z  \t\t dx \t\t dy \t\t dz \t\t G\n")

  for ind in points.ispace do
    C.fprintf(fp, "%d \t %d \t %d \t %f \t %f \t %f \t %f \t %f \t %f \t %f\n",
                   ind.x, ind.y, ind.z,
                   points[ind].centerCoordinates[0], points[ind].centerCoordinates[1], points[ind].centerCoordinates[2],
                   points[ind].cellWidth[0], points[ind].cellWidth[1], points[ind].cellWidth[2],
                   points[ind].G)
  end

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

  --for c in tiles do
  --  Radiation_InitializeGeometry(p_points[c],
  --                               config.Grid.xType, config.Grid.yType, config.Grid.zType,
  --                               config.Radiation.u.DOM.xNum, config.Grid.origin[0], config.Grid.xWidth,
  --                               config.Radiation.u.DOM.yNum, config.Grid.origin[1], config.Grid.yWidth,
  --                               config.Radiation.u.DOM.zNum, config.Grid.origin[2], config.Grid.zWidth)
  --end
  Radiation_InitializeGeometry(points,
                               config.Grid.xType, config.Grid.yType, config.Grid.zType,
                               config.Radiation.u.DOM.xNum, config.Grid.origin[0], config.Grid.xWidth,
                               config.Radiation.u.DOM.yNum, config.Grid.origin[1], config.Grid.yWidth,
                               config.Radiation.u.DOM.zNum, config.Grid.origin[2], config.Grid.zWidth)

  -- Prepare fake inputs
  --fill(points.Ib, (SB/PI) * pow(1000.0,4.0))
  --fill(points.sigma, 5.0);
  fill(points.Ib, 0.0)
  fill(points.sigma, config.Radiation.u.DOM.qa + config.Radiation.u.DOM.qs);

  -- Invoke DOM solver
  [DOM_INST.ComputeRadiationField(config, tiles, p_points)];
  -- Output results
  --writeIntensity(points)
  writeOutput(points)

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
