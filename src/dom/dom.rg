import 'regent'

local A = require 'admiral'

-------------------------------------------------------------------------------
-- MODULE PARAMETERS
-------------------------------------------------------------------------------

return function(radiationRel)

local points = radiationRel:regionSymbol()
local p_points = radiationRel:primPartSymbol()
local pointsType = radiationRel:regionType()
local tiles = A.primColors()

-------------------------------------------------------------------------------
-- COMPILE-TIME COMPUTATION
-------------------------------------------------------------------------------

-- Read configuration

local config
local i = 1
while i <= #arg do
  if arg[i] == '-i' and i < #arg then
    config = loadfile(arg[i+1])()
    break
  end
  i = i + 1
end
if not config then
  print('config file required (-i <file> option)')
  os.exit(1)
end

-- C imports

local c     = regentlib.c
local cmath = terralib.includec('math.h')

-- Some math definitions

local min = regentlib.fmin
local max = regentlib.fmax
local pow = regentlib.pow(double)
local pi  = 2.0*cmath.acos(0.0)

-- Quadrature file name
local quad_file = 'LMquads/'..config.numAngles..'.txt'

-- Grid size

local Nx = radiationRel:Dims()[1]
local Ny = radiationRel:Dims()[2]
local Nz = radiationRel:Dims()[3]

local N_angles = config.numAngles

-- Domain size

local Lx = config.xWidth
local Ly = config.yWidth
local Lz = config.zWidth

-- Grid spacing

local dx = Lx/Nx
local dy = Ly/Ny
local dz = Lz/Nz

local dAx = dy*dz;
local dAy = dx*dz;
local dAz = dx*dy;
local dV = dx*dy*dz;

-- Grid tiling

local ntx,nty,ntz = A.primPartDims()

-- Albedo

local qa = config.qa
local qs = config.qs
local omega = qs/(qa+qs)

-- Wall emissivity

local emiss_east  = config.emiss_east  or 1.0
local emiss_west  = config.emiss_west  or 1.0
local emiss_south = config.emiss_south or 1.0
local emiss_north = config.emiss_north or 1.0
local emiss_up    = config.emiss_up    or 1.0
local emiss_down  = config.emiss_down  or 1.0

-- Wall temperatures

local SB = 5.67e-8

local T_east  = config.T_east  or 300.0
local T_west  = config.T_west  or 300.0
local T_south = config.T_south or 300.0
local T_north = config.T_north or 300.0
local T_up    = config.T_up    or 300.0
local T_down  = config.T_down  or 300.0

-- TODO: Read location of irradiated window

-- Procedure parameters

local tol   = 1e-6   -- solution tolerance
local gamma = 0.5    -- 1 for step differencing, 0.5 for diamond differencing

local terra read_val(f : &c.FILE, val : &double)
  return c.fscanf(f, '%lf\n', &val[0])
end

-------------------------------------------------------------------------------
-- MODULE-LOCAL FIELD SPACES
-------------------------------------------------------------------------------

-- Internal cell values are essentially private,
-- face values are what need to be passed to downstream neighbor
-- Update cell value, then update downstream face values

-- quadrature information
local fspace angle {
  xi  : double,
  eta : double,
  mu  : double,
  w   : double,
}

local fspace face {
  I : double[N_angles],    -- face intensity per angle
}

-------------------------------------------------------------------------------
-- MODULE-LOCAL TASKS
-------------------------------------------------------------------------------

-- Initialize angle quads
local task initialize_angles(angles : region(ispace(int1d), angle))
where
  reads writes (angles.{xi, eta, mu, w})
do

  -- Read angle_value information from file.

  var val : double[1]

  var f = c.fopen([quad_file], 'rb')

  read_val(f, val) -- gets rid of num angles

  for a in angles do
    read_val(f, val)
    a.xi = val[0]
  end

  for a in angles do
    read_val(f, val)
    a.eta = val[0]
  end

  for a in angles do
    read_val(f, val)
    a.mu = val[0]
  end

  for a in angles do
    read_val(f, val)
    a.w = val[0]
  end

  c.fclose(f)

end

----------------------------------------------------------------------------


local task make_private_partition_x(faces : region(ispace(int3d), face),
                                        x_tiles : ispace(int3d))

  var coloring = c.legion_domain_point_coloring_create()
  for tile in x_tiles do

    var lo = int3d { x = tile.x     * Nx / ntx + 1,
                     y = tile.y     * Ny / nty,
                     z = tile.z     * Nz / ntz}
    var hi = int3d { x = (tile.x+1) * Nx / ntx - 1,
                     y = (tile.y+1) * Ny / nty - 1,
                     z = (tile.z+1) * Nz / ntz - 1}

    if hi.x < Nx+1 then               

      var rect = rect3d {lo = lo, hi = hi}
      c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
    end                
    
  end
  var p = partition(disjoint, faces, coloring, x_tiles)
  c.legion_domain_point_coloring_destroy(coloring)
  return p

end

-- 1 - 2 - 3
local task make_shared_partition_x(faces : region(ispace(int3d), face),
                                        x_tiles : ispace(int3d))

  var coloring = c.legion_domain_point_coloring_create()
  for tile in x_tiles do


    var lo = int3d { x = tile.x     * Nx / ntx,
                     y = tile.y     * Ny / nty,
                     z = tile.z     * Nz / ntz}
    var hi = int3d { x = lo.x,
                     y = (tile.y+1) * Ny / nty - 1,
                     z = (tile.z+1) * Nz / ntz - 1}             

    var rect = rect3d {lo = lo, hi = hi}
    c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
  end
  var p = partition(disjoint, faces, coloring, x_tiles)
  c.legion_domain_point_coloring_destroy(coloring)
  return p

end

----------------------------------------------------------------------------


local task make_private_partition_y(faces : region(ispace(int3d), face),
                                        y_tiles : ispace(int3d))

  var coloring = c.legion_domain_point_coloring_create()
  for tile in y_tiles do

    var lo = int3d { x = tile.x     * Nx / ntx,
                     y = tile.y     * Ny / nty + 1,
                     z = tile.z     * Nz / ntz}
    var hi = int3d { x = (tile.x+1) * Nx / ntx - 1,
                     y = (tile.y+1) * Ny / nty - 1,
                     z = (tile.z+1) * Nz / ntz - 1}

    if hi.y < Ny+1 then             

      var rect = rect3d {lo = lo, hi = hi}
      c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
    end                
    
  end
  var p = partition(disjoint, faces, coloring, y_tiles)
  c.legion_domain_point_coloring_destroy(coloring)
  return p

end

-- 1 - 2 - 3
local task make_shared_partition_y(faces : region(ispace(int3d), face),
                                        y_tiles : ispace(int3d))

  var coloring = c.legion_domain_point_coloring_create()
  for tile in y_tiles do


    var lo = int3d { x = tile.x     * Nx / ntx,
                     y = tile.y     * Ny / nty,
                     z = tile.z     * Nz / ntz}
    var hi = int3d { x = (tile.x+1) * Nx / ntx - 1,
                     y = lo.y,
                     z = (tile.z+1) * Nz / ntz - 1}     

    var rect = rect3d {lo = lo, hi = hi}
    c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
  end
  var p = partition(disjoint, faces, coloring, y_tiles)
  c.legion_domain_point_coloring_destroy(coloring)
  return p

end

----------------------------------------------------------------------------


local task make_private_partition_z(faces : region(ispace(int3d), face),
                                        z_tiles : ispace(int3d))

  var coloring = c.legion_domain_point_coloring_create()
  for tile in z_tiles do

    var lo = int3d { x = tile.x     * Nx / ntx,
                     y = tile.y     * Ny / nty,
                     z = tile.z     * Nz / ntz + 1}
    var hi = int3d { x = (tile.x+1) * Nx / ntx - 1,
                     y = (tile.y+1) * Ny / nty - 1,
                     z = (tile.z+1) * Nz / ntz - 1}

    if hi.z < Nz+1 then                 

      var rect = rect3d {lo = lo, hi = hi}
      c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
    end                
    
  end
  var p = partition(disjoint, faces, coloring, z_tiles)
  c.legion_domain_point_coloring_destroy(coloring)
  return p

end

-- 1 - 2 - 3
local task make_shared_partition_z(faces : region(ispace(int3d), face),
                                        z_tiles : ispace(int3d))

  var coloring = c.legion_domain_point_coloring_create()
  for tile in z_tiles do


    var lo = int3d { x = tile.x     * Nx / ntx,
                     y = tile.y     * Ny / nty,
                     z = tile.z     * Nz / ntz}
    var hi = int3d { x = (tile.x+1) * Nx / ntx - 1,
                     y = (tile.y+1) * Ny / nty - 1,
                     z = lo.z}              

    var rect = rect3d {lo = lo, hi = hi}
    c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
  end
  var p = partition(disjoint, faces, coloring, z_tiles)
  c.legion_domain_point_coloring_destroy(coloring)
  return p

end

----------------------------------------------------------------------------


local task make_interior_partition_x_hi(faces : region(ispace(int3d), face),
                                        x_tiles : ispace(int3d))

  var coloring = c.legion_domain_point_coloring_create()
  for tile in x_tiles do

    -- include extra face in last partition
    var val : int = -1
    if tile.x == x_tiles.bounds.hi.x-1 then val = 0 end

    var lo = int3d { x = tile.x     * Nx / ntx,
                     y = tile.y     * Ny / nty,
                     z = tile.z     * Nz / ntz}
    var hi = int3d { x = (tile.x+1) * Nx / ntx + val,
                     y = (tile.y+1) * Ny / nty - 1,
                     z = (tile.z+1) * Nz / ntz - 1}

    -- Create an empty partition
    if hi.x >= Nx+1 then
      lo.x = 1
      hi.x = 0
    end

    var rect = rect3d {lo = lo, hi = hi}
    c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
  end
  var p = partition(disjoint, faces, coloring, x_tiles)
  c.legion_domain_point_coloring_destroy(coloring)
  return p
end

local task make_interior_partition_x_lo(faces : region(ispace(int3d), face),
                                        x_tiles : ispace(int3d))

  var coloring = c.legion_domain_point_coloring_create()
  for tile in x_tiles do

    -- include extra face in first partition
    var val : int = 1
    if tile.x == 1 then val = 0 end

    var lo = int3d { x = (tile.x-1) * Nx / ntx + val,
                     y = tile.y     * Ny / nty,
                     z = tile.z     * Nz / ntz}
    var hi = int3d { x = tile.x     * Nx / ntx,
                     y = (tile.y+1) * Ny / nty - 1,
                     z = (tile.z+1) * Nz / ntz - 1}

    -- Create an empty partition
    if lo.x < 0 then
      lo.x = 1
      hi.x = 0
    end

    var rect = rect3d {lo = lo, hi = hi}
    c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
  end
  var p = partition(disjoint, faces, coloring, x_tiles)
  c.legion_domain_point_coloring_destroy(coloring)
  return p
end

local task make_interior_partition_y_hi(faces : region(ispace(int3d), face),
                                        y_tiles : ispace(int3d))

  var coloring = c.legion_domain_point_coloring_create()
  for tile in y_tiles do

    -- include extra face in last partition
    var val : int = -1
    if tile.y == y_tiles.bounds.hi.y-1 then val = 0 end

    var lo = int3d { x = tile.x     * Nx / ntx,
                     y = tile.y     * Ny / nty,
                     z = tile.z     * Nz / ntz}
    var hi = int3d { x = (tile.x+1) * Nx / ntx - 1,
                     y = (tile.y+1) * Ny / nty + val,
                     z = (tile.z+1) * Nz / ntz - 1}

    -- Create an empty partition
    if hi.y >= Ny+1 then
      lo.y = 1
      hi.y = 0
    end

    var rect = rect3d {lo = lo, hi = hi}
    c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
  end
  var p = partition(disjoint, faces, coloring, y_tiles)
  c.legion_domain_point_coloring_destroy(coloring)
  return p
end

local task make_interior_partition_y_lo(faces : region(ispace(int3d), face),
                                        y_tiles : ispace(int3d))

  var coloring = c.legion_domain_point_coloring_create()
  for tile in y_tiles do

    -- include extra face in first partition
    var val : int = 1
    if tile.y == 1 then val = 0 end

    var lo = int3d { x = tile.x     * Nx / ntx,
                     y = (tile.y-1) * Ny / nty + val,
                     z = tile.z     * Nz / ntz}
    var hi = int3d { x = (tile.x+1) * Nx / ntx - 1,
                     y = tile.y     * Ny / nty,
                     z = (tile.z+1) * Nz / ntz - 1}

    -- Create an empty partition
    if lo.y < 0 then
      lo.y = 1
      hi.y = 0
    end

    var rect = rect3d {lo = lo, hi = hi}
    c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
  end
  var p = partition(disjoint, faces, coloring, y_tiles)
  c.legion_domain_point_coloring_destroy(coloring)
  return p
end

local task make_interior_partition_z_hi(faces : region(ispace(int3d), face),
                                        z_tiles : ispace(int3d))

  var coloring = c.legion_domain_point_coloring_create()
  for tile in z_tiles do

    -- include extra face in last partition
    var val : int = -1
    if tile.z == z_tiles.bounds.hi.z-1 then val = 0 end

    var lo = int3d { x = tile.x     * Nx / ntx,
                     y = tile.y     * Ny / nty,
                     z = tile.z     * Nz / ntz}
    var hi = int3d { x = (tile.x+1) * Nx / ntx - 1,
                     y = (tile.y+1) * Ny / nty - 1,
                     z = (tile.z+1) * Nz / ntz + val}

    -- Create an empty partition
    if hi.z >= Nz+1 then
      lo.z = 1
      hi.z = 0
    end

    var rect = rect3d {lo = lo, hi = hi}
    c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
  end
  var p = partition(disjoint, faces, coloring, z_tiles)
  c.legion_domain_point_coloring_destroy(coloring)
  return p
end

local task make_interior_partition_z_lo(faces : region(ispace(int3d), face),
                                        z_tiles : ispace(int3d))

  var coloring = c.legion_domain_point_coloring_create()
  for tile in z_tiles do

    -- include extra face in first partition
    var val : int = 1
    if tile.z == 1 then val = 0 end

    var lo = int3d { x = tile.x     * Nx / ntx,
                     y = tile.y     * Ny / nty,
                     z = (tile.z-1) * Nz / ntz + val}
    var hi = int3d { x = (tile.x+1) * Nx / ntx - 1,
                     y = (tile.y+1) * Ny / nty - 1,
                     z = tile.z     * Nz / ntz}

    -- Create an empty partition
    if lo.z < 0 then
      lo.z = 1
      hi.z = 0
    end

    var rect = rect3d {lo = lo, hi = hi}
    c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
  end
  var p = partition(disjoint, faces, coloring, z_tiles)
  c.legion_domain_point_coloring_destroy(coloring)
  return p
end

-- Loop over all angles and grid cells to compute the source term
-- for the current iteration.
local task source_term(points : pointsType,
                       angles : region(ispace(int1d), angle))
where
  reads (points.{Iiter_1, Iiter_2, Iiter_3, Iiter_4,
                 Iiter_5, Iiter_6, Iiter_7, Iiter_8,
                 Ib, sigma},
         angles.w),
  reads writes (points.S)
do
  for p in points do
    p.S = (1.0-omega) * p.sigma * p.Ib
    for m = 0, N_angles do
      p.S += omega * p.sigma/(4.0*pi) * angles[m].w * p.Iiter_1[m]
           + omega * p.sigma/(4.0*pi) * angles[m].w * p.Iiter_2[m]
           + omega * p.sigma/(4.0*pi) * angles[m].w * p.Iiter_3[m]
           + omega * p.sigma/(4.0*pi) * angles[m].w * p.Iiter_4[m]
           + omega * p.sigma/(4.0*pi) * angles[m].w * p.Iiter_5[m]
           + omega * p.sigma/(4.0*pi) * angles[m].w * p.Iiter_6[m]
           + omega * p.sigma/(4.0*pi) * angles[m].w * p.Iiter_7[m]
           + omega * p.sigma/(4.0*pi) * angles[m].w * p.Iiter_8[m]
    end
  end
end

local task west_bound(faces_1 : region(ispace(int3d), face),
                      faces_2 : region(ispace(int3d), face),
                      faces_3 : region(ispace(int3d), face),
                      faces_4 : region(ispace(int3d), face),
                      faces_5 : region(ispace(int3d), face),
                      faces_6 : region(ispace(int3d), face),
                      faces_7 : region(ispace(int3d), face),
                      faces_8 : region(ispace(int3d), face),
                      angles : region(ispace(int1d), angle))
where
  reads (angles.{w, xi, eta, mu}),
  reads writes (faces_1.I, faces_2.I, faces_3.I, faces_4.I,
                faces_5.I, faces_6.I, faces_7.I, faces_8.I)
do

  -- Get array bounds

  var limits = faces_1.bounds

  -- Temporary variables

  var reflect : double = 0.0
  var epsw    : double = emiss_west
  var Tw      : double = T_west

  var value : double = epsw*SB*pow(Tw,4.0)/pi

  for j = limits.lo.y, limits.hi.y + 1 do
    for k = limits.lo.z, limits.hi.z + 1 do

      -- Calculate reflect

      reflect = 0
      for m = 0, N_angles do
        if angles[m].xi < 0 then
          var face_value : double = 0.0
          if angles[m].eta > 0 and angles[m].mu > 0 then
            face_value = faces_5[{limits.lo.x,j,k}].I[m]
          elseif angles[m].eta > 0 and angles[m].mu <= 0 then
            face_value = faces_6[{limits.lo.x,j,k}].I[m]
          elseif angles[m].eta <= 0 and angles[m].mu > 0 then
            face_value = faces_7[{limits.lo.x,j,k}].I[m]
          else
            face_value = faces_8[{limits.lo.x,j,k}].I[m]
          end
          reflect += (1.0-epsw)/pi*angles[m].w*cmath.fabs(angles[m].xi)*face_value
        end
      end

      -- Set Ifx values using reflect

      var value : double = epsw*SB*pow(Tw,4.0)/pi + reflect
      for m = 0, N_angles do
        if angles[m].xi > 0 then
          if angles[m].eta > 0 and angles[m].mu > 0 then
            faces_1[{limits.lo.x,j,k}].I[m] = value
          elseif angles[m].eta > 0 and angles[m].mu <= 0 then
            faces_2[{limits.lo.x,j,k}].I[m] = value
          elseif angles[m].eta <= 0 and angles[m].mu > 0 then
            faces_3[{limits.lo.x,j,k}].I[m] = value
          else
            faces_4[{limits.lo.x,j,k}].I[m] = value
          end
        end
      end

    end
  end

end

local task east_bound(faces_1 : region(ispace(int3d), face),
                      faces_2 : region(ispace(int3d), face),
                      faces_3 : region(ispace(int3d), face),
                      faces_4 : region(ispace(int3d), face),
                      faces_5 : region(ispace(int3d), face),
                      faces_6 : region(ispace(int3d), face),
                      faces_7 : region(ispace(int3d), face),
                      faces_8 : region(ispace(int3d), face),
                      angles : region(ispace(int1d), angle))
where
  reads (angles.{w, xi, eta, mu}),
  reads writes (faces_1.I, faces_2.I, faces_3.I, faces_4.I,
                faces_5.I, faces_6.I, faces_7.I, faces_8.I)
do

  -- Get array bounds

  var limits = faces_1.bounds

  -- Temporary variables

  var reflect : double = 0.0
  var epsw    : double = emiss_east
  var Tw      : double = T_east

  var value : double = epsw*SB*pow(Tw,4.0)/pi

  for j = limits.lo.y, limits.hi.y + 1 do
    for k = limits.lo.z, limits.hi.z + 1 do

      -- Calculate reflect

      reflect = 0
      for m = 0, N_angles do
        if angles[m].xi > 0 then
          var face_value : double = 0.0
          if angles[m].eta > 0 and angles[m].mu > 0 then
            face_value = faces_1[{limits.hi.x,j,k}].I[m]
          elseif angles[m].eta > 0 and angles[m].mu <= 0 then
            face_value = faces_2[{limits.hi.x,j,k}].I[m]
          elseif angles[m].eta <= 0 and angles[m].mu > 0 then
            face_value = faces_3[{limits.hi.x,j,k}].I[m]
          else
            face_value = faces_4[{limits.hi.x,j,k}].I[m]
          end
          reflect += (1.0-epsw)/pi*angles[m].w*angles[m].xi*face_value
        end
      end

      -- Set Ifx values using reflect

      var value : double = epsw*SB*pow(Tw,4.0)/pi + reflect
      for m = 0, N_angles do
        if angles[m].xi < 0 then
          if angles[m].eta > 0 and angles[m].mu > 0 then
            faces_5[{limits.hi.x,j,k}].I[m] = value
          elseif angles[m].eta > 0 and angles[m].mu <= 0 then
            faces_6[{limits.hi.x,j,k}].I[m] = value
          elseif angles[m].eta <= 0 and angles[m].mu > 0 then
            faces_7[{limits.hi.x,j,k}].I[m] = value
          else
            faces_8[{limits.hi.x,j,k}].I[m] = value
          end
        end
      end

    end
  end

end

local task north_bound(faces_1 : region(ispace(int3d), face),
                      faces_2 : region(ispace(int3d), face),
                      faces_3 : region(ispace(int3d), face),
                      faces_4 : region(ispace(int3d), face),
                      faces_5 : region(ispace(int3d), face),
                      faces_6 : region(ispace(int3d), face),
                      faces_7 : region(ispace(int3d), face),
                      faces_8 : region(ispace(int3d), face),
                      angles : region(ispace(int1d), angle))
where
  reads (angles.{w, xi, eta, mu}),
  reads writes (faces_1.I, faces_2.I, faces_3.I, faces_4.I,
                faces_5.I, faces_6.I, faces_7.I, faces_8.I)
do

  -- Get array bounds

  var limits = faces_1.bounds

  -- Temporary variables

  var reflect : double = 0.0
  var epsw    : double = emiss_north
  var Tw      : double = T_north

  var value : double = epsw*SB*pow(Tw,4.0)/pi

  for i = limits.lo.x, limits.hi.x + 1 do
    for k = limits.lo.z, limits.hi.z + 1 do

      -- Calculate reflect

      reflect = 0
      for m = 0, N_angles do
        if angles[m].eta > 0 then
          var face_value : double = 0.0
          if angles[m].xi > 0 and angles[m].mu > 0 then
            face_value = faces_1[{i,limits.hi.y,k}].I[m]
          elseif angles[m].xi > 0 and angles[m].mu <= 0 then
            face_value = faces_2[{i,limits.hi.y,k}].I[m]
          elseif angles[m].xi <= 0 and angles[m].mu > 0 then
            face_value = faces_5[{i,limits.hi.y,k}].I[m]
          else
            face_value = faces_6[{i,limits.hi.y,k}].I[m]
          end
          reflect += (1.0-epsw)/pi*angles[m].w*angles[m].eta*face_value
        end
      end

      -- Set Ify values using reflect

      var value : double = epsw*SB*pow(Tw,4.0)/pi + reflect
      for m = 0, N_angles do
        if angles[m].eta < 0 then

          if angles[m].xi > 0 and angles[m].mu > 0 then
            faces_3[{i,limits.hi.y,k}].I[m] = value
          elseif angles[m].xi > 0 and angles[m].mu <= 0 then
            faces_4[{i,limits.hi.y,k}].I[m] = value
          elseif angles[m].xi <= 0 and angles[m].mu > 0 then
            faces_7[{i,limits.hi.y,k}].I[m] = value
          else
            faces_8[{i,limits.hi.y,k}].I[m] = value
          end
        end
      end

    end
  end

end

local task south_bound(faces_1 : region(ispace(int3d), face),
                      faces_2 : region(ispace(int3d), face),
                      faces_3 : region(ispace(int3d), face),
                      faces_4 : region(ispace(int3d), face),
                      faces_5 : region(ispace(int3d), face),
                      faces_6 : region(ispace(int3d), face),
                      faces_7 : region(ispace(int3d), face),
                      faces_8 : region(ispace(int3d), face),
                      angles : region(ispace(int1d), angle))
where
  reads (angles.{w, xi, eta, mu}),
  reads writes (faces_1.I, faces_2.I, faces_3.I, faces_4.I,
                faces_5.I, faces_6.I, faces_7.I, faces_8.I)
do

  -- Get array bounds

  var limits = faces_1.bounds

  -- Temporary variables

  var reflect : double = 0.0
  var epsw    : double = emiss_south
  var Tw      : double = T_south

  var value : double = epsw*SB*pow(Tw,4.0)/pi

  for i = limits.lo.x, limits.hi.x + 1 do
    for k = limits.lo.z, limits.hi.z + 1 do

      -- Calculate reflect

      reflect = 0
      for m = 0, N_angles do
        if angles[m].eta < 0 then
          var face_value : double = 0.0
          if angles[m].xi > 0 and angles[m].mu > 0 then
            face_value = faces_3[{i,limits.lo.y,k}].I[m]
          elseif angles[m].xi > 0 and angles[m].mu <= 0 then
            face_value = faces_4[{i,limits.lo.y,k}].I[m]
          elseif angles[m].xi <= 0 and angles[m].mu > 0 then
            face_value = faces_7[{i,limits.lo.y,k}].I[m]
          else
            face_value = faces_8[{i,limits.lo.y,k}].I[m]
          end
          reflect += (1.0-epsw)/pi*angles[m].w*cmath.fabs(angles[m].eta)*face_value
        end
      end

      -- Set Ify values using reflect

      var value : double = epsw*SB*pow(Tw,4.0)/pi + reflect
      for m = 0, N_angles do
        if angles[m].eta > 0 then

          if angles[m].xi > 0 and angles[m].mu > 0 then
            faces_1[{i,limits.lo.y,k}].I[m] = value
          elseif angles[m].xi > 0 and angles[m].mu <= 0 then
            faces_2[{i,limits.lo.y,k}].I[m] = value
          elseif angles[m].xi <= 0 and angles[m].mu > 0 then
            faces_5[{i,limits.lo.y,k}].I[m] = value
          else
            faces_6[{i,limits.lo.y,k}].I[m] = value
          end
        end
      end

    end
  end

end

local task up_bound(faces_1 : region(ispace(int3d), face),
                      faces_2 : region(ispace(int3d), face),
                      faces_3 : region(ispace(int3d), face),
                      faces_4 : region(ispace(int3d), face),
                      faces_5 : region(ispace(int3d), face),
                      faces_6 : region(ispace(int3d), face),
                      faces_7 : region(ispace(int3d), face),
                      faces_8 : region(ispace(int3d), face),
                      angles : region(ispace(int1d), angle))
where
  reads (angles.{w, xi, eta, mu}),
  reads writes (faces_1.I, faces_2.I, faces_3.I, faces_4.I,
                faces_5.I, faces_6.I, faces_7.I, faces_8.I)
do

  -- Get array bounds

  var limits = faces_1.bounds

  -- Temporary variables

  var reflect : double = 0.0
  var epsw    : double = emiss_up
  var Tw      : double = T_up

  var value : double = epsw*SB*pow(Tw,4.0)/pi

  for i = limits.lo.x, limits.hi.x + 1 do
    for j = limits.lo.y, limits.hi.y + 1 do

      -- Calculate reflect

      reflect = 0
      for m = 0, N_angles do
        if angles[m].mu < 0 then
          var face_value : double = 0.0
          if angles[m].xi > 0 and angles[m].eta > 0 then
            face_value = faces_2[{i,j,limits.lo.z}].I[m]
          elseif angles[m].xi > 0 and angles[m].eta <= 0 then
            face_value = faces_4[{i,j,limits.lo.z}].I[m]
          elseif angles[m].xi <= 0 and angles[m].eta > 0 then
            face_value = faces_6[{i,j,limits.lo.z}].I[m]
          else
            face_value = faces_8[{i,j,limits.lo.z}].I[m]
          end
          reflect += (1.0-epsw)/pi*angles[m].w*cmath.fabs(angles[m].mu)*face_value
        end
      end

      -- Set Ifz values using reflect

      var value : double = epsw*SB*pow(Tw,4.0)/pi + reflect
      for m = 0, N_angles do
        if angles[m].mu > 0 then

          if angles[m].xi > 0 and angles[m].eta > 0 then
            faces_1[{i,j,limits.lo.z}].I[m] = value
          elseif angles[m].xi > 0 and angles[m].eta <= 0 then
            faces_3[{i,j,limits.lo.z}].I[m] = value
          elseif angles[m].xi <= 0 and angles[m].eta > 0 then
            faces_5[{i,j,limits.lo.z}].I[m] = value
          else
            faces_7[{i,j,limits.lo.z}].I[m] = value
          end
        end
      end

    end
  end

end

local task down_bound(faces_1 : region(ispace(int3d), face),
                      faces_2 : region(ispace(int3d), face),
                      faces_3 : region(ispace(int3d), face),
                      faces_4 : region(ispace(int3d), face),
                      faces_5 : region(ispace(int3d), face),
                      faces_6 : region(ispace(int3d), face),
                      faces_7 : region(ispace(int3d), face),
                      faces_8 : region(ispace(int3d), face),
                      angles : region(ispace(int1d), angle))
where
  reads (angles.{w, xi, eta, mu}),
  reads writes (faces_1.I, faces_2.I, faces_3.I, faces_4.I,
                faces_5.I, faces_6.I, faces_7.I, faces_8.I)
do

  -- Get array bounds

  var limits = faces_1.bounds

  -- Temporary variables

  var reflect : double = 0.0
  var epsw    : double = emiss_down
  var Tw      : double = T_down

  var value : double = epsw*SB*pow(Tw,4.0)/pi

  for i = limits.lo.x, limits.hi.x + 1 do
    for j = limits.lo.y, limits.hi.y + 1 do

      -- Calculate reflect

      reflect = 0
      for m = 0, N_angles do
        if angles[m].mu > 0 then
          var face_value : double = 0.0
          if angles[m].xi > 0 and angles[m].eta > 0 then
            face_value = faces_1[{i,j,limits.hi.z}].I[m]
          elseif angles[m].xi > 0 and angles[m].eta <= 0 then
            face_value = faces_3[{i,j,limits.hi.z}].I[m]
          elseif angles[m].xi <= 0 and angles[m].eta > 0 then
            face_value = faces_5[{i,j,limits.hi.z}].I[m]
          else
            face_value = faces_7[{i,j,limits.hi.z}].I[m]
          end
          reflect += (1.0-epsw)/pi*angles[m].w*angles[m].mu*face_value
        end
      end

      -- Set Ifz values using reflect

      var value : double = epsw*SB*pow(Tw,4.0)/pi + reflect
      for m = 0, N_angles do
        if angles[m].mu < 0 then

          if angles[m].xi > 0 and angles[m].eta > 0 then
            faces_2[{i,j,limits.hi.z}].I[m] = value
          elseif angles[m].xi > 0 and angles[m].eta <= 0 then
            faces_4[{i,j,limits.hi.z}].I[m] = value
          elseif angles[m].xi <= 0 and angles[m].eta > 0 then
            faces_6[{i,j,limits.hi.z}].I[m] = value
          else
            faces_8[{i,j,limits.hi.z}].I[m] = value
          end
        end
      end
    end
  end
end

-- Sweeps are all exactly the same except for array read/written to in points

local task sweep_1(points : pointsType,
                   x_faces : region(ispace(int3d), face),
                   y_faces : region(ispace(int3d), face),
                   z_faces : region(ispace(int3d), face),
                   shared_x_faces_upwind : region(ispace(int3d), face),
                   shared_x_faces_downwind : region(ispace(int3d), face),
                   shared_y_faces_upwind : region(ispace(int3d), face),
                   shared_y_faces_downwind : region(ispace(int3d), face),
                   shared_z_faces_upwind : region(ispace(int3d), face),
                   shared_z_faces_downwind : region(ispace(int3d), face),
                   angles : region(ispace(int1d), angle),
                   xi : int64, eta : int64, mu : int64)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         shared_x_faces_upwind.I, shared_y_faces_upwind.I, shared_z_faces_upwind.I),
  reads writes(points.I_1, x_faces.I, y_faces.I, z_faces.I, 
    shared_x_faces_downwind.I, shared_y_faces_downwind.I, shared_z_faces_downwind.I)
do

  -- Determine sweep direction and bounds

  var limits = points.bounds

  var dindx  : int64 = 1
  var startx : int64 = limits.lo.x
  var endx   : int64 = limits.hi.x + 1

  var dindy  : int64 = 1
  var starty : int64 = limits.lo.y
  var endy   : int64 = limits.hi.y + 1

  var dindz  : int64 = 1
  var startz : int64 = limits.lo.z
  var endz   : int64 = limits.hi.z + 1

  if xi < 0 then
    dindx = -1
    startx = limits.hi.x
    endx = limits.lo.x - 1
  end

  if eta < 0 then
    dindy = -1
    starty = limits.hi.y
    endy = limits.lo.y - 1
  end

  if mu < 0 then
    dindz = -1
    startz = limits.hi.z
    endz = limits.lo.z - 1
  end


  -- Outer loop over all angles.
  for m = 0, N_angles do

    if (angles[m].xi * xi > 0 or (angles[m].xi == 0 and xi < 0)) and 
      (angles[m].eta * eta > 0 or (angles[m].eta == 0 and eta < 0)) and 
      (angles[m].mu * mu > 0 or (angles[m].mu == 0 and mu < 0)) then

      -- Use our direction and increments for the sweep.

      for k = startz,endz,dindz do
        for j = starty,endy,dindy do
          for i = startx,endx,dindx do

            -- indx and indy are the upwind indices
            var indx : int64 = i - min(dindx,0)
            var indy : int64 = j - min(dindy,0)
            var indz : int64 = k - min(dindz,0)

            -- Determine if necessary to use ghost partition

            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x or indx > x_faces.bounds.hi.x then
              upwind_x_value = shared_x_faces_upwind[{0,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            ---

            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y or indy > y_faces.bounds.hi.y then
              upwind_y_value = shared_y_faces_upwind[{i,0,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z or indz > z_faces.bounds.hi.z then
              upwind_z_value = shared_z_faces_upwind[{i,j,0}].I[m]
            else
              upwind_z_value = z_faces[{i,j,indz}].I[m]
            end

            -- Integrate to compute cell-centered value of I.

            points[{i,j,k}].I_1[m] = (points[{i,j,k}].S * dV
                                        + cmath.fabs(angles[m].xi) * dAx * upwind_x_value/gamma
                                        + cmath.fabs(angles[m].eta) * dAy * upwind_y_value/gamma
                                        + cmath.fabs(angles[m].mu) * dAz * upwind_z_value/gamma)
              /(points[{i,j,k}].sigma * dV
                  + cmath.fabs(angles[m].xi) * dAx/gamma
                  + cmath.fabs(angles[m].eta) * dAy/gamma
                  + cmath.fabs(angles[m].mu) * dAz/gamma)

            -- Compute intensities on downwind faces

            var x_face_val = (points[{i,j,k}].I_1[m] - (1-gamma)*upwind_x_value)/gamma
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{0, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_1[m] - (1-gamma)*upwind_y_value)/gamma
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, 0, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_1[m] - (1-gamma)*upwind_z_value)/gamma
            if (indz + dindz) > z_faces.bounds.hi.z or (indz + dindz) < z_faces.bounds.lo.z then
              shared_z_faces_downwind[{i, j, 0}].I[m] = z_face_val
            else
              z_faces[{i, j, indz+dindz}].I[m] = z_face_val
            end
            
          end
        end
      end
    end
  end
end

local task sweep_2(points : pointsType,
                   x_faces : region(ispace(int3d), face),
                   y_faces : region(ispace(int3d), face),
                   z_faces : region(ispace(int3d), face),
                   shared_x_faces_upwind : region(ispace(int3d), face),
                   shared_x_faces_downwind : region(ispace(int3d), face),
                   shared_y_faces_upwind : region(ispace(int3d), face),
                   shared_y_faces_downwind : region(ispace(int3d), face),
                   shared_z_faces_upwind : region(ispace(int3d), face),
                   shared_z_faces_downwind : region(ispace(int3d), face),
                   angles : region(ispace(int1d), angle),
                   xi : int64, eta : int64, mu : int64)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         shared_x_faces_upwind.I, shared_y_faces_upwind.I, shared_z_faces_upwind.I),
  reads writes(points.I_2, x_faces.I, y_faces.I, z_faces.I, 
    shared_x_faces_downwind.I, shared_y_faces_downwind.I, shared_z_faces_downwind.I)
do

  -- Determine sweep direction and bounds

  var limits = points.bounds

  var dindx  : int64 = 1
  var startx : int64 = limits.lo.x
  var endx   : int64 = limits.hi.x + 1

  var dindy  : int64 = 1
  var starty : int64 = limits.lo.y
  var endy   : int64 = limits.hi.y + 1

  var dindz  : int64 = 1
  var startz : int64 = limits.lo.z
  var endz   : int64 = limits.hi.z + 1

  if xi < 0 then
    dindx = -1
    startx = limits.hi.x
    endx = limits.lo.x - 1
  end

  if eta < 0 then
    dindy = -1
    starty = limits.hi.y
    endy = limits.lo.y - 1
  end

  if mu < 0 then
    dindz = -1
    startz = limits.hi.z
    endz = limits.lo.z - 1
  end


  -- Outer loop over all angles.
  for m = 0, N_angles do

    if (angles[m].xi * xi > 0 or (angles[m].xi == 0 and xi < 0)) and 
      (angles[m].eta * eta > 0 or (angles[m].eta == 0 and eta < 0)) and 
      (angles[m].mu * mu > 0 or (angles[m].mu == 0 and mu < 0)) then

      -- Use our direction and increments for the sweep.

      for k = startz,endz,dindz do
        for j = starty,endy,dindy do
          for i = startx,endx,dindx do

            -- indx and indy are the upwind indices
            var indx : int64 = i - min(dindx,0)
            var indy : int64 = j - min(dindy,0)
            var indz : int64 = k - min(dindz,0)

            -- Determine if necessary to use ghost partition

            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x or indx > x_faces.bounds.hi.x then
              upwind_x_value = shared_x_faces_upwind[{0,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            ---

            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y or indy > y_faces.bounds.hi.y then
              upwind_y_value = shared_y_faces_upwind[{i,0,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z or indz > z_faces.bounds.hi.z then
              upwind_z_value = shared_z_faces_upwind[{i,j,0}].I[m]
            else
              upwind_z_value = z_faces[{i,j,indz}].I[m]
            end

            -- Integrate to compute cell-centered value of I.

            points[{i,j,k}].I_2[m] = (points[{i,j,k}].S * dV
                                        + cmath.fabs(angles[m].xi) * dAx * upwind_x_value/gamma
                                        + cmath.fabs(angles[m].eta) * dAy * upwind_y_value/gamma
                                        + cmath.fabs(angles[m].mu) * dAz * upwind_z_value/gamma)
              /(points[{i,j,k}].sigma * dV
                  + cmath.fabs(angles[m].xi) * dAx/gamma
                  + cmath.fabs(angles[m].eta) * dAy/gamma
                  + cmath.fabs(angles[m].mu) * dAz/gamma)

            -- Compute intensities on downwind faces

            var x_face_val = (points[{i,j,k}].I_2[m] - (1-gamma)*upwind_x_value)/gamma
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{0, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_2[m] - (1-gamma)*upwind_y_value)/gamma
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, 0, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_2[m] - (1-gamma)*upwind_z_value)/gamma
            if (indz + dindz) > z_faces.bounds.hi.z or (indz + dindz) < z_faces.bounds.lo.z then
              shared_z_faces_downwind[{i, j, 0}].I[m] = z_face_val
            else
              z_faces[{i, j, indz+dindz}].I[m] = z_face_val
            end
            
          end
        end
      end
    end
  end
end

local task sweep_3(points : pointsType,
                   x_faces : region(ispace(int3d), face),
                   y_faces : region(ispace(int3d), face),
                   z_faces : region(ispace(int3d), face),
                   shared_x_faces_upwind : region(ispace(int3d), face),
                   shared_x_faces_downwind : region(ispace(int3d), face),
                   shared_y_faces_upwind : region(ispace(int3d), face),
                   shared_y_faces_downwind : region(ispace(int3d), face),
                   shared_z_faces_upwind : region(ispace(int3d), face),
                   shared_z_faces_downwind : region(ispace(int3d), face),
                   angles : region(ispace(int1d), angle),
                   xi : int64, eta : int64, mu : int64)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         shared_x_faces_upwind.I, shared_y_faces_upwind.I, shared_z_faces_upwind.I),
  reads writes(points.I_3, x_faces.I, y_faces.I, z_faces.I, 
    shared_x_faces_downwind.I, shared_y_faces_downwind.I, shared_z_faces_downwind.I)
do

  -- Determine sweep direction and bounds

  var limits = points.bounds

  var dindx  : int64 = 1
  var startx : int64 = limits.lo.x
  var endx   : int64 = limits.hi.x + 1

  var dindy  : int64 = 1
  var starty : int64 = limits.lo.y
  var endy   : int64 = limits.hi.y + 1

  var dindz  : int64 = 1
  var startz : int64 = limits.lo.z
  var endz   : int64 = limits.hi.z + 1

  if xi < 0 then
    dindx = -1
    startx = limits.hi.x
    endx = limits.lo.x - 1
  end

  if eta < 0 then
    dindy = -1
    starty = limits.hi.y
    endy = limits.lo.y - 1
  end

  if mu < 0 then
    dindz = -1
    startz = limits.hi.z
    endz = limits.lo.z - 1
  end


  -- Outer loop over all angles.
  for m = 0, N_angles do

    if (angles[m].xi * xi > 0 or (angles[m].xi == 0 and xi < 0)) and 
      (angles[m].eta * eta > 0 or (angles[m].eta == 0 and eta < 0)) and 
      (angles[m].mu * mu > 0 or (angles[m].mu == 0 and mu < 0)) then

      -- Use our direction and increments for the sweep.

      for k = startz,endz,dindz do
        for j = starty,endy,dindy do
          for i = startx,endx,dindx do

            -- indx and indy are the upwind indices
            var indx : int64 = i - min(dindx,0)
            var indy : int64 = j - min(dindy,0)
            var indz : int64 = k - min(dindz,0)

            -- Determine if necessary to use ghost partition

            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x or indx > x_faces.bounds.hi.x then
              upwind_x_value = shared_x_faces_upwind[{0,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            ---

            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y or indy > y_faces.bounds.hi.y then
              upwind_y_value = shared_y_faces_upwind[{i,0,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z or indz > z_faces.bounds.hi.z then
              upwind_z_value = shared_z_faces_upwind[{i,j,0}].I[m]
            else
              upwind_z_value = z_faces[{i,j,indz}].I[m]
            end

            -- Integrate to compute cell-centered value of I.

            points[{i,j,k}].I_3[m] = (points[{i,j,k}].S * dV
                                        + cmath.fabs(angles[m].xi) * dAx * upwind_x_value/gamma
                                        + cmath.fabs(angles[m].eta) * dAy * upwind_y_value/gamma
                                        + cmath.fabs(angles[m].mu) * dAz * upwind_z_value/gamma)
              /(points[{i,j,k}].sigma * dV
                  + cmath.fabs(angles[m].xi) * dAx/gamma
                  + cmath.fabs(angles[m].eta) * dAy/gamma
                  + cmath.fabs(angles[m].mu) * dAz/gamma)

            -- Compute intensities on downwind faces

            var x_face_val = (points[{i,j,k}].I_3[m] - (1-gamma)*upwind_x_value)/gamma
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{0, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_3[m] - (1-gamma)*upwind_y_value)/gamma
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, 0, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_3[m] - (1-gamma)*upwind_z_value)/gamma
            if (indz + dindz) > z_faces.bounds.hi.z or (indz + dindz) < z_faces.bounds.lo.z then
              shared_z_faces_downwind[{i, j, 0}].I[m] = z_face_val
            else
              z_faces[{i, j, indz+dindz}].I[m] = z_face_val
            end
            
          end
        end
      end
    end
  end
end

local task sweep_4(points : pointsType,
                   x_faces : region(ispace(int3d), face),
                   y_faces : region(ispace(int3d), face),
                   z_faces : region(ispace(int3d), face),
                   shared_x_faces_upwind : region(ispace(int3d), face),
                   shared_x_faces_downwind : region(ispace(int3d), face),
                   shared_y_faces_upwind : region(ispace(int3d), face),
                   shared_y_faces_downwind : region(ispace(int3d), face),
                   shared_z_faces_upwind : region(ispace(int3d), face),
                   shared_z_faces_downwind : region(ispace(int3d), face),
                   angles : region(ispace(int1d), angle),
                   xi : int64, eta : int64, mu : int64)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         shared_x_faces_upwind.I, shared_y_faces_upwind.I, shared_z_faces_upwind.I),
  reads writes(points.I_4, x_faces.I, y_faces.I, z_faces.I, 
    shared_x_faces_downwind.I, shared_y_faces_downwind.I, shared_z_faces_downwind.I)
do

  -- Determine sweep direction and bounds

  var limits = points.bounds

  var dindx  : int64 = 1
  var startx : int64 = limits.lo.x
  var endx   : int64 = limits.hi.x + 1

  var dindy  : int64 = 1
  var starty : int64 = limits.lo.y
  var endy   : int64 = limits.hi.y + 1

  var dindz  : int64 = 1
  var startz : int64 = limits.lo.z
  var endz   : int64 = limits.hi.z + 1

  if xi < 0 then
    dindx = -1
    startx = limits.hi.x
    endx = limits.lo.x - 1
  end

  if eta < 0 then
    dindy = -1
    starty = limits.hi.y
    endy = limits.lo.y - 1
  end

  if mu < 0 then
    dindz = -1
    startz = limits.hi.z
    endz = limits.lo.z - 1
  end


  -- Outer loop over all angles.
  for m = 0, N_angles do

    if (angles[m].xi * xi > 0 or (angles[m].xi == 0 and xi < 0)) and 
      (angles[m].eta * eta > 0 or (angles[m].eta == 0 and eta < 0)) and 
      (angles[m].mu * mu > 0 or (angles[m].mu == 0 and mu < 0)) then

      -- Use our direction and increments for the sweep.

      for k = startz,endz,dindz do
        for j = starty,endy,dindy do
          for i = startx,endx,dindx do

            -- indx and indy are the upwind indices
            var indx : int64 = i - min(dindx,0)
            var indy : int64 = j - min(dindy,0)
            var indz : int64 = k - min(dindz,0)

            -- Determine if necessary to use ghost partition

            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x or indx > x_faces.bounds.hi.x then
              upwind_x_value = shared_x_faces_upwind[{0,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            ---

            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y or indy > y_faces.bounds.hi.y then
              upwind_y_value = shared_y_faces_upwind[{i,0,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z or indz > z_faces.bounds.hi.z then
              upwind_z_value = shared_z_faces_upwind[{i,j,0}].I[m]
            else
              upwind_z_value = z_faces[{i,j,indz}].I[m]
            end

            -- Integrate to compute cell-centered value of I.

            points[{i,j,k}].I_4[m] = (points[{i,j,k}].S * dV
                                        + cmath.fabs(angles[m].xi) * dAx * upwind_x_value/gamma
                                        + cmath.fabs(angles[m].eta) * dAy * upwind_y_value/gamma
                                        + cmath.fabs(angles[m].mu) * dAz * upwind_z_value/gamma)
              /(points[{i,j,k}].sigma * dV
                  + cmath.fabs(angles[m].xi) * dAx/gamma
                  + cmath.fabs(angles[m].eta) * dAy/gamma
                  + cmath.fabs(angles[m].mu) * dAz/gamma)

            -- Compute intensities on downwind faces

            var x_face_val = (points[{i,j,k}].I_4[m] - (1-gamma)*upwind_x_value)/gamma
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{0, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_4[m] - (1-gamma)*upwind_y_value)/gamma
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, 0, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_4[m] - (1-gamma)*upwind_z_value)/gamma
            if (indz + dindz) > z_faces.bounds.hi.z or (indz + dindz) < z_faces.bounds.lo.z then
              shared_z_faces_downwind[{i, j, 0}].I[m] = z_face_val
            else
              z_faces[{i, j, indz+dindz}].I[m] = z_face_val
            end
            
          end
        end
      end
    end
  end
end

local task sweep_5(points : pointsType,
                   x_faces : region(ispace(int3d), face),
                   y_faces : region(ispace(int3d), face),
                   z_faces : region(ispace(int3d), face),
                   shared_x_faces_upwind : region(ispace(int3d), face),
                   shared_x_faces_downwind : region(ispace(int3d), face),
                   shared_y_faces_upwind : region(ispace(int3d), face),
                   shared_y_faces_downwind : region(ispace(int3d), face),
                   shared_z_faces_upwind : region(ispace(int3d), face),
                   shared_z_faces_downwind : region(ispace(int3d), face),
                   angles : region(ispace(int1d), angle),
                   xi : int64, eta : int64, mu : int64)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         shared_x_faces_upwind.I, shared_y_faces_upwind.I, shared_z_faces_upwind.I),
  reads writes(points.I_5, x_faces.I, y_faces.I, z_faces.I, 
    shared_x_faces_downwind.I, shared_y_faces_downwind.I, shared_z_faces_downwind.I)
do

  -- Determine sweep direction and bounds

  var limits = points.bounds

  var dindx  : int64 = 1
  var startx : int64 = limits.lo.x
  var endx   : int64 = limits.hi.x + 1

  var dindy  : int64 = 1
  var starty : int64 = limits.lo.y
  var endy   : int64 = limits.hi.y + 1

  var dindz  : int64 = 1
  var startz : int64 = limits.lo.z
  var endz   : int64 = limits.hi.z + 1

  if xi < 0 then
    dindx = -1
    startx = limits.hi.x
    endx = limits.lo.x - 1
  end

  if eta < 0 then
    dindy = -1
    starty = limits.hi.y
    endy = limits.lo.y - 1
  end

  if mu < 0 then
    dindz = -1
    startz = limits.hi.z
    endz = limits.lo.z - 1
  end


  -- Outer loop over all angles.
  for m = 0, N_angles do

    if (angles[m].xi * xi > 0 or (angles[m].xi == 0 and xi < 0)) and 
      (angles[m].eta * eta > 0 or (angles[m].eta == 0 and eta < 0)) and 
      (angles[m].mu * mu > 0 or (angles[m].mu == 0 and mu < 0)) then

      -- Use our direction and increments for the sweep.

      for k = startz,endz,dindz do
        for j = starty,endy,dindy do
          for i = startx,endx,dindx do

            -- indx and indy are the upwind indices
            var indx : int64 = i - min(dindx,0)
            var indy : int64 = j - min(dindy,0)
            var indz : int64 = k - min(dindz,0)

            -- Determine if necessary to use ghost partition

            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x or indx > x_faces.bounds.hi.x then
              upwind_x_value = shared_x_faces_upwind[{0,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            ---

            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y or indy > y_faces.bounds.hi.y then
              upwind_y_value = shared_y_faces_upwind[{i,0,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z or indz > z_faces.bounds.hi.z then
              upwind_z_value = shared_z_faces_upwind[{i,j,0}].I[m]
            else
              upwind_z_value = z_faces[{i,j,indz}].I[m]
            end

            -- Integrate to compute cell-centered value of I.

            points[{i,j,k}].I_5[m] = (points[{i,j,k}].S * dV
                                        + cmath.fabs(angles[m].xi) * dAx * upwind_x_value/gamma
                                        + cmath.fabs(angles[m].eta) * dAy * upwind_y_value/gamma
                                        + cmath.fabs(angles[m].mu) * dAz * upwind_z_value/gamma)
              /(points[{i,j,k}].sigma * dV
                  + cmath.fabs(angles[m].xi) * dAx/gamma
                  + cmath.fabs(angles[m].eta) * dAy/gamma
                  + cmath.fabs(angles[m].mu) * dAz/gamma)

            -- Compute intensities on downwind faces

            var x_face_val = (points[{i,j,k}].I_5[m] - (1-gamma)*upwind_x_value)/gamma
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{0, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_5[m] - (1-gamma)*upwind_y_value)/gamma
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, 0, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_5[m] - (1-gamma)*upwind_z_value)/gamma
            if (indz + dindz) > z_faces.bounds.hi.z or (indz + dindz) < z_faces.bounds.lo.z then
              shared_z_faces_downwind[{i, j, 0}].I[m] = z_face_val
            else
              z_faces[{i, j, indz+dindz}].I[m] = z_face_val
            end
            
          end
        end
      end
    end
  end
end

local task sweep_6(points : pointsType,
                   x_faces : region(ispace(int3d), face),
                   y_faces : region(ispace(int3d), face),
                   z_faces : region(ispace(int3d), face),
                   shared_x_faces_upwind : region(ispace(int3d), face),
                   shared_x_faces_downwind : region(ispace(int3d), face),
                   shared_y_faces_upwind : region(ispace(int3d), face),
                   shared_y_faces_downwind : region(ispace(int3d), face),
                   shared_z_faces_upwind : region(ispace(int3d), face),
                   shared_z_faces_downwind : region(ispace(int3d), face),
                   angles : region(ispace(int1d), angle),
                   xi : int64, eta : int64, mu : int64)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         shared_x_faces_upwind.I, shared_y_faces_upwind.I, shared_z_faces_upwind.I),
  reads writes(points.I_6, x_faces.I, y_faces.I, z_faces.I, 
    shared_x_faces_downwind.I, shared_y_faces_downwind.I, shared_z_faces_downwind.I)
do

  -- Determine sweep direction and bounds

  var limits = points.bounds

  var dindx  : int64 = 1
  var startx : int64 = limits.lo.x
  var endx   : int64 = limits.hi.x + 1

  var dindy  : int64 = 1
  var starty : int64 = limits.lo.y
  var endy   : int64 = limits.hi.y + 1

  var dindz  : int64 = 1
  var startz : int64 = limits.lo.z
  var endz   : int64 = limits.hi.z + 1

  if xi < 0 then
    dindx = -1
    startx = limits.hi.x
    endx = limits.lo.x - 1
  end

  if eta < 0 then
    dindy = -1
    starty = limits.hi.y
    endy = limits.lo.y - 1
  end

  if mu < 0 then
    dindz = -1
    startz = limits.hi.z
    endz = limits.lo.z - 1
  end


  -- Outer loop over all angles.
  for m = 0, N_angles do

    if (angles[m].xi * xi > 0 or (angles[m].xi == 0 and xi < 0)) and 
      (angles[m].eta * eta > 0 or (angles[m].eta == 0 and eta < 0)) and 
      (angles[m].mu * mu > 0 or (angles[m].mu == 0 and mu < 0)) then

      -- Use our direction and increments for the sweep.

      for k = startz,endz,dindz do
        for j = starty,endy,dindy do
          for i = startx,endx,dindx do

            -- indx and indy are the upwind indices
            var indx : int64 = i - min(dindx,0)
            var indy : int64 = j - min(dindy,0)
            var indz : int64 = k - min(dindz,0)

            -- Determine if necessary to use ghost partition

            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x or indx > x_faces.bounds.hi.x then
              upwind_x_value = shared_x_faces_upwind[{0,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            ---

            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y or indy > y_faces.bounds.hi.y then
              upwind_y_value = shared_y_faces_upwind[{i,0,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z or indz > z_faces.bounds.hi.z then
              upwind_z_value = shared_z_faces_upwind[{i,j,0}].I[m]
            else
              upwind_z_value = z_faces[{i,j,indz}].I[m]
            end

            -- Integrate to compute cell-centered value of I.

            points[{i,j,k}].I_6[m] = (points[{i,j,k}].S * dV
                                        + cmath.fabs(angles[m].xi) * dAx * upwind_x_value/gamma
                                        + cmath.fabs(angles[m].eta) * dAy * upwind_y_value/gamma
                                        + cmath.fabs(angles[m].mu) * dAz * upwind_z_value/gamma)
              /(points[{i,j,k}].sigma * dV
                  + cmath.fabs(angles[m].xi) * dAx/gamma
                  + cmath.fabs(angles[m].eta) * dAy/gamma
                  + cmath.fabs(angles[m].mu) * dAz/gamma)

            -- Compute intensities on downwind faces

            var x_face_val = (points[{i,j,k}].I_6[m] - (1-gamma)*upwind_x_value)/gamma
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{0, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_6[m] - (1-gamma)*upwind_y_value)/gamma
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, 0, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_6[m] - (1-gamma)*upwind_z_value)/gamma
            if (indz + dindz) > z_faces.bounds.hi.z or (indz + dindz) < z_faces.bounds.lo.z then
              shared_z_faces_downwind[{i, j, 0}].I[m] = z_face_val
            else
              z_faces[{i, j, indz+dindz}].I[m] = z_face_val
            end
            
          end
        end
      end
    end
  end
end

local task sweep_7(points : pointsType,
                   x_faces : region(ispace(int3d), face),
                   y_faces : region(ispace(int3d), face),
                   z_faces : region(ispace(int3d), face),
                   shared_x_faces_upwind : region(ispace(int3d), face),
                   shared_x_faces_downwind : region(ispace(int3d), face),
                   shared_y_faces_upwind : region(ispace(int3d), face),
                   shared_y_faces_downwind : region(ispace(int3d), face),
                   shared_z_faces_upwind : region(ispace(int3d), face),
                   shared_z_faces_downwind : region(ispace(int3d), face),
                   angles : region(ispace(int1d), angle),
                   xi : int64, eta : int64, mu : int64)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         shared_x_faces_upwind.I, shared_y_faces_upwind.I, shared_z_faces_upwind.I),
  reads writes(points.I_7, x_faces.I, y_faces.I, z_faces.I, 
    shared_x_faces_downwind.I, shared_y_faces_downwind.I, shared_z_faces_downwind.I)
do

  -- Determine sweep direction and bounds

  var limits = points.bounds

  var dindx  : int64 = 1
  var startx : int64 = limits.lo.x
  var endx   : int64 = limits.hi.x + 1

  var dindy  : int64 = 1
  var starty : int64 = limits.lo.y
  var endy   : int64 = limits.hi.y + 1

  var dindz  : int64 = 1
  var startz : int64 = limits.lo.z
  var endz   : int64 = limits.hi.z + 1

  if xi < 0 then
    dindx = -1
    startx = limits.hi.x
    endx = limits.lo.x - 1
  end

  if eta < 0 then
    dindy = -1
    starty = limits.hi.y
    endy = limits.lo.y - 1
  end

  if mu < 0 then
    dindz = -1
    startz = limits.hi.z
    endz = limits.lo.z - 1
  end


  -- Outer loop over all angles.
  for m = 0, N_angles do

    if (angles[m].xi * xi > 0 or (angles[m].xi == 0 and xi < 0)) and 
      (angles[m].eta * eta > 0 or (angles[m].eta == 0 and eta < 0)) and 
      (angles[m].mu * mu > 0 or (angles[m].mu == 0 and mu < 0)) then

      -- Use our direction and increments for the sweep.

      for k = startz,endz,dindz do
        for j = starty,endy,dindy do
          for i = startx,endx,dindx do

            -- indx and indy are the upwind indices
            var indx : int64 = i - min(dindx,0)
            var indy : int64 = j - min(dindy,0)
            var indz : int64 = k - min(dindz,0)

            -- Determine if necessary to use ghost partition

            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x or indx > x_faces.bounds.hi.x then
              upwind_x_value = shared_x_faces_upwind[{0,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            ---

            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y or indy > y_faces.bounds.hi.y then
              upwind_y_value = shared_y_faces_upwind[{i,0,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z or indz > z_faces.bounds.hi.z then
              upwind_z_value = shared_z_faces_upwind[{i,j,0}].I[m]
            else
              upwind_z_value = z_faces[{i,j,indz}].I[m]
            end

            -- Integrate to compute cell-centered value of I.

            points[{i,j,k}].I_7[m] = (points[{i,j,k}].S * dV
                                        + cmath.fabs(angles[m].xi) * dAx * upwind_x_value/gamma
                                        + cmath.fabs(angles[m].eta) * dAy * upwind_y_value/gamma
                                        + cmath.fabs(angles[m].mu) * dAz * upwind_z_value/gamma)
              /(points[{i,j,k}].sigma * dV
                  + cmath.fabs(angles[m].xi) * dAx/gamma
                  + cmath.fabs(angles[m].eta) * dAy/gamma
                  + cmath.fabs(angles[m].mu) * dAz/gamma)

            -- Compute intensities on downwind faces

            var x_face_val = (points[{i,j,k}].I_7[m] - (1-gamma)*upwind_x_value)/gamma
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{0, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_7[m] - (1-gamma)*upwind_y_value)/gamma
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, 0, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_7[m] - (1-gamma)*upwind_z_value)/gamma
            if (indz + dindz) > z_faces.bounds.hi.z or (indz + dindz) < z_faces.bounds.lo.z then
              shared_z_faces_downwind[{i, j, 0}].I[m] = z_face_val
            else
              z_faces[{i, j, indz+dindz}].I[m] = z_face_val
            end
            
          end
        end
      end
    end
  end
end

local task sweep_8(points : pointsType,
                   x_faces : region(ispace(int3d), face),
                   y_faces : region(ispace(int3d), face),
                   z_faces : region(ispace(int3d), face),
                   shared_x_faces_upwind : region(ispace(int3d), face),
                   shared_x_faces_downwind : region(ispace(int3d), face),
                   shared_y_faces_upwind : region(ispace(int3d), face),
                   shared_y_faces_downwind : region(ispace(int3d), face),
                   shared_z_faces_upwind : region(ispace(int3d), face),
                   shared_z_faces_downwind : region(ispace(int3d), face),
                   angles : region(ispace(int1d), angle),
                   xi : int64, eta : int64, mu : int64)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         shared_x_faces_upwind.I, shared_y_faces_upwind.I, shared_z_faces_upwind.I),
  reads writes(points.I_8, x_faces.I, y_faces.I, z_faces.I, 
    shared_x_faces_downwind.I, shared_y_faces_downwind.I, shared_z_faces_downwind.I)
do

  -- Determine sweep direction and bounds

  var limits = points.bounds

  var dindx  : int64 = 1
  var startx : int64 = limits.lo.x
  var endx   : int64 = limits.hi.x + 1

  var dindy  : int64 = 1
  var starty : int64 = limits.lo.y
  var endy   : int64 = limits.hi.y + 1

  var dindz  : int64 = 1
  var startz : int64 = limits.lo.z
  var endz   : int64 = limits.hi.z + 1

  if xi < 0 then
    dindx = -1
    startx = limits.hi.x
    endx = limits.lo.x - 1
  end

  if eta < 0 then
    dindy = -1
    starty = limits.hi.y
    endy = limits.lo.y - 1
  end

  if mu < 0 then
    dindz = -1
    startz = limits.hi.z
    endz = limits.lo.z - 1
  end


  -- Outer loop over all angles.
  for m = 0, N_angles do

    if (angles[m].xi * xi > 0 or (angles[m].xi == 0 and xi < 0)) and 
      (angles[m].eta * eta > 0 or (angles[m].eta == 0 and eta < 0)) and 
      (angles[m].mu * mu > 0 or (angles[m].mu == 0 and mu < 0)) then

      -- Use our direction and increments for the sweep.

      for k = startz,endz,dindz do
        for j = starty,endy,dindy do
          for i = startx,endx,dindx do

            -- indx and indy are the upwind indices
            var indx : int64 = i - min(dindx,0)
            var indy : int64 = j - min(dindy,0)
            var indz : int64 = k - min(dindz,0)

            -- Determine if necessary to use ghost partition

            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x or indx > x_faces.bounds.hi.x then
              upwind_x_value = shared_x_faces_upwind[{0,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            ---

            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y or indy > y_faces.bounds.hi.y then
              upwind_y_value = shared_y_faces_upwind[{i,0,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z or indz > z_faces.bounds.hi.z then
              upwind_z_value = shared_z_faces_upwind[{i,j,0}].I[m]
            else
              upwind_z_value = z_faces[{i,j,indz}].I[m]
            end

            -- Integrate to compute cell-centered value of I.

            points[{i,j,k}].I_8[m] = (points[{i,j,k}].S * dV
                                        + cmath.fabs(angles[m].xi) * dAx * upwind_x_value/gamma
                                        + cmath.fabs(angles[m].eta) * dAy * upwind_y_value/gamma
                                        + cmath.fabs(angles[m].mu) * dAz * upwind_z_value/gamma)
              /(points[{i,j,k}].sigma * dV
                  + cmath.fabs(angles[m].xi) * dAx/gamma
                  + cmath.fabs(angles[m].eta) * dAy/gamma
                  + cmath.fabs(angles[m].mu) * dAz/gamma)

            -- Compute intensities on downwind faces

            var x_face_val = (points[{i,j,k}].I_8[m] - (1-gamma)*upwind_x_value)/gamma
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{0, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_8[m] - (1-gamma)*upwind_y_value)/gamma
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, 0, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_8[m] - (1-gamma)*upwind_z_value)/gamma
            if (indz + dindz) > z_faces.bounds.hi.z or (indz + dindz) < z_faces.bounds.lo.z then
              shared_z_faces_downwind[{i, j, 0}].I[m] = z_face_val
            else
              z_faces[{i, j, indz+dindz}].I[m] = z_face_val
            end
            
          end
        end
      end
    end
  end
end

-- Compute the residual after each iteration and return the value.
local task residual(points : pointsType)
where
  reads (points.{I_1, I_2, I_3, I_4, I_5, I_6, I_7, I_8,
                 Iiter_1, Iiter_2, Iiter_3, Iiter_4,
                 Iiter_5, Iiter_6, Iiter_7, Iiter_8})
do
  var res : double = 0.0

  for p in points do
    for m = 0, N_angles do

      if p.I_1[m] > 0 then
        res += (1.0/(Nx*Ny*Nz*(N_angles)))
          * pow((p.I_1[m]-p.Iiter_1[m]),2.0)
          / pow((p.I_1[m]),2.0)
      end

      if p.I_2[m] > 0 then
        res += (1.0/(Nx*Ny*Nz*(N_angles)))
          * pow((p.I_2[m]-p.Iiter_2[m]),2.0)
          / pow((p.I_2[m]),2.0)
      end

      if p.I_3[m] > 0 then
        res += (1.0/(Nx*Ny*Nz*(N_angles)))
          * pow((p.I_3[m]-p.Iiter_3[m]),2.0)
          / pow((p.I_3[m]),2.0)
      end

      if p.I_4[m] > 0 then
        res += (1.0/(Nx*Ny*Nz*(N_angles)))
          * pow((p.I_4[m]-p.Iiter_4[m]),2.0)
          / pow((p.I_4[m]),2.0)
      end

      if p.I_5[m] > 0 then
        res += (1.0/(Nx*Ny*Nz*(N_angles)))
          * pow((p.I_5[m]-p.Iiter_5[m]),2.0)
          / pow((p.I_5[m]),2.0)
      end

      if p.I_6[m] > 0 then
        res += (1.0/(Nx*Ny*Nz*(N_angles)))
          * pow((p.I_6[m]-p.Iiter_6[m]),2.0)
          / pow((p.I_6[m]),2.0)
      end

      if p.I_7[m] > 0 then
        res += (1.0/(Nx*Ny*Nz*(N_angles)))
          * pow((p.I_7[m]-p.Iiter_7[m]),2.0)
          / pow((p.I_7[m]),2.0)
      end

      if p.I_8[m] > 0 then
        res += (1.0/(Nx*Ny*Nz*(N_angles)))
          * pow((p.I_8[m]-p.Iiter_8[m]),2.0)
          / pow((p.I_8[m]),2.0)
      end

    end
  end

  return res
end

-- Update the intensity before moving to the next iteration.
local task update(points : pointsType)
where
  reads (points.{I_1, I_2, I_3, I_4, I_5, I_6, I_7, I_8}),
  reads writes (points.{Iiter_1, Iiter_2, Iiter_3, Iiter_4,
                        Iiter_5, Iiter_6, Iiter_7, Iiter_8})
do
  for p in points do
    for m = 0, N_angles do
      p.Iiter_1[m] = p.I_1[m]
      p.Iiter_2[m] = p.I_2[m]
      p.Iiter_3[m] = p.I_3[m]
      p.Iiter_4[m] = p.I_4[m]
      p.Iiter_5[m] = p.I_5[m]
      p.Iiter_6[m] = p.I_6[m]
      p.Iiter_7[m] = p.I_7[m]
      p.Iiter_8[m] = p.I_8[m]
    end
  end
end

-- Reduce the intensity to summation over all angles
local task reduce_intensity(points : pointsType,
                            angles : region(ispace(int1d), angle))
where
  reads (points.{I_1, I_2, I_3, I_4, I_5, I_6, I_7, I_8}, angles.w),
  reads writes (points.G)
do
  for p in points do
    for m = 0, N_angles do
      p.G += angles[m].w * p.I_1[m]
           + angles[m].w * p.I_2[m]
           + angles[m].w * p.I_3[m]
           + angles[m].w * p.I_4[m]
           + angles[m].w * p.I_5[m]
           + angles[m].w * p.I_6[m]
           + angles[m].w * p.I_7[m]
           + angles[m].w * p.I_8[m]
    end
  end
end

local task write_intensity(points : pointsType)
where
  reads (points.G)
do
  var limits = points.bounds
  var f = c.fopen("intensity.dat", "w")
  for i = limits.lo.x, limits.hi.x+1 do
    for j = limits.lo.y, limits.hi.y+1 do
      for k = limits.lo.z, limits.hi.z+1 do
        c.fprintf(f,' %.6e ', points[{i,j,k}].G)
      end
      c.fprintf(f,'\n')
    end
    c.fprintf(f,'\n')
  end
  c.fclose(f)
end

-- for debugging
local task print_final_intensities(points : pointsType)
where
  reads (points.{I_1, I_2, I_3, I_4, I_5, I_6, I_7, I_8})

do

    var limits = points.bounds
    for i = limits.lo.x, limits.hi.x+1 do
      for j = limits.lo.y, limits.hi.y+1 do
        for k = limits.lo.z, limits.hi.z+1 do
            for m = 0, N_angles do

              if points[{i,j,k}].I_1[m] > 0 then
                c.printf("1 x=%d,y=%d,z=%d,angle=%d I = %lf", i, j, k, m, points[{i,j,k}].I_1[m])
              end 

              if points[{i,j,k}].I_2[m] > 0 then
                c.printf(" 2 x=%d,y=%d,z=%d,angle=%d I = %lf", i, j, k, m, points[{i,j,k}].I_2[m])
              end 

              if points[{i,j,k}].I_3[m] > 0 then
                c.printf(" 3 x=%d,y=%d,z=%d,angle=%d I = %lf", i, j, k, m, points[{i,j,k}].I_3[m])
              end 

              if points[{i,j,k}].I_4[m] > 0 then
                c.printf(" 4 x=%d,y=%d,z=%d,angle=%d I = %lf", i, j, k, m, points[{i,j,k}].I_4[m])
              end 

              if points[{i,j,k}].I_5[m] > 0 then
                c.printf(" 5 x=%d,y=%d,z=%d,angle=%d I = %lf", i, j, k, m, points[{i,j,k}].I_5[m])
              end 

              if points[{i,j,k}].I_6[m] > 0 then
                c.printf(" 6 x=%d,y=%d,z=%d,angle=%d I = %lf", i, j, k, m, points[{i,j,k}].I_6[m])
              end 

              if points[{i,j,k}].I_7[m] > 0 then
                c.printf(" 7 x=%d,y=%d,z=%d,angle=%d I = %lf", i, j, k, m, points[{i,j,k}].I_7[m])
              end 

              if points[{i,j,k}].I_8[m] > 0 then
                c.printf(" 8 x=%d,y=%d,z=%d,angle=%d I = %lf", i, j, k, m, points[{i,j,k}].I_8[m])
              end 

            end
            c.printf("\n")
          end
          c.printf("\n")
      end
      c.printf("\n")
    end

end

-------------------------------------------------------------------------------
-- EXPORTED QUOTES
-------------------------------------------------------------------------------

local exports = {}

-- Symbols shared between quotes
local angles = regentlib.newsymbol('angles')
local x_tiles = regentlib.newsymbol('x_tiles')
local y_tiles = regentlib.newsymbol('y_tiles')
local z_tiles = regentlib.newsymbol('z_tiles')

local s_x_faces = {}
local s_y_faces = {}
local s_z_faces = {}

local p_x_faces_1 = regentlib.newsymbol('p_x_faces_1')
local p_y_faces_1 = regentlib.newsymbol('p_y_faces_1')
local p_z_faces_1 = regentlib.newsymbol('p_z_faces_1')

s_x_faces[1] = regentlib.newsymbol('s_x_faces_1')
s_y_faces[1] = regentlib.newsymbol('s_y_faces_1')
s_z_faces[1] = regentlib.newsymbol('s_z_faces_1')

local p_x_faces_2 = regentlib.newsymbol('p_x_faces_2')
local p_y_faces_2 = regentlib.newsymbol('p_y_faces_2')
local p_z_faces_2 = regentlib.newsymbol('p_z_faces_2')

s_x_faces[2] = regentlib.newsymbol('s_x_faces_2')
s_y_faces[2] = regentlib.newsymbol('s_y_faces_2')
s_z_faces[2] = regentlib.newsymbol('s_z_faces_2')

local p_x_faces_3 = regentlib.newsymbol('p_x_faces_3')
local p_y_faces_3 = regentlib.newsymbol('p_y_faces_3')
local p_z_faces_3 = regentlib.newsymbol('p_z_faces_3')

s_x_faces[3] = regentlib.newsymbol('s_x_faces_3')
s_y_faces[3] = regentlib.newsymbol('s_y_faces_3')
s_z_faces[3] = regentlib.newsymbol('s_z_faces_3')

local p_x_faces_4 = regentlib.newsymbol('p_x_faces_4')
local p_y_faces_4 = regentlib.newsymbol('p_y_faces_4')
local p_z_faces_4 = regentlib.newsymbol('p_z_faces_4')

s_x_faces[4] = regentlib.newsymbol('s_x_faces_4')
s_y_faces[4] = regentlib.newsymbol('s_y_faces_4')
s_z_faces[4] = regentlib.newsymbol('s_z_faces_4')

local p_x_faces_5 = regentlib.newsymbol('p_x_faces_5')
local p_y_faces_5 = regentlib.newsymbol('p_y_faces_5')
local p_z_faces_5 = regentlib.newsymbol('p_z_faces_5')

s_x_faces[5] = regentlib.newsymbol('s_x_faces_5')
s_y_faces[5] = regentlib.newsymbol('s_y_faces_5')
s_z_faces[5] = regentlib.newsymbol('s_z_faces_5')

local p_x_faces_6 = regentlib.newsymbol('p_x_faces_6')
local p_y_faces_6 = regentlib.newsymbol('p_y_faces_6')
local p_z_faces_6 = regentlib.newsymbol('p_z_faces_6')

s_x_faces[6] = regentlib.newsymbol('s_x_faces_6')
s_y_faces[6] = regentlib.newsymbol('s_y_faces_6')
s_z_faces[6] = regentlib.newsymbol('s_z_faces_6')

local p_x_faces_7 = regentlib.newsymbol('p_x_faces_7')
local p_y_faces_7 = regentlib.newsymbol('p_y_faces_7')
local p_z_faces_7 = regentlib.newsymbol('p_z_faces_7')

s_x_faces[7] = regentlib.newsymbol('s_x_faces_7')
s_y_faces[7] = regentlib.newsymbol('s_y_faces_7')
s_z_faces[7] = regentlib.newsymbol('s_z_faces_7')

local p_x_faces_8 = regentlib.newsymbol('p_x_faces_8')
local p_y_faces_8 = regentlib.newsymbol('p_y_faces_8')
local p_z_faces_8 = regentlib.newsymbol('p_z_faces_8')

s_x_faces[8] = regentlib.newsymbol('s_x_faces_8')
s_y_faces[8] = regentlib.newsymbol('s_y_faces_8')
s_z_faces[8] = regentlib.newsymbol('s_z_faces_8')

local res = regentlib.newsymbol(double, 'res')

exports.InitModule = rquote

  -- Regions for faces (+1 in one direction)
  var grid_x = ispace(int3d, {x = Nx+1, y = Ny,   z = Nz})
  var grid_y = ispace(int3d, {x = Nx,   y = Ny+1, z = Nz})
  var grid_z = ispace(int3d, {x = Nx,   y = Ny,   z = Nz+1})

  var x_faces_1 = region(grid_x, face)
  var x_faces_2 = region(grid_x, face)
  var x_faces_3 = region(grid_x, face)
  var x_faces_4 = region(grid_x, face)
  var x_faces_5 = region(grid_x, face)
  var x_faces_6 = region(grid_x, face)
  var x_faces_7 = region(grid_x, face)
  var x_faces_8 = region(grid_x, face)

  var y_faces_1 = region(grid_y, face)
  var y_faces_2 = region(grid_y, face)
  var y_faces_3 = region(grid_y, face)
  var y_faces_4 = region(grid_y, face)
  var y_faces_5 = region(grid_y, face)
  var y_faces_6 = region(grid_y, face)
  var y_faces_7 = region(grid_y, face)
  var y_faces_8 = region(grid_y, face)

  var z_faces_1 = region(grid_z, face)
  var z_faces_2 = region(grid_z, face)
  var z_faces_3 = region(grid_z, face)
  var z_faces_4 = region(grid_z, face)
  var z_faces_5 = region(grid_z, face)
  var z_faces_6 = region(grid_z, face)
  var z_faces_7 = region(grid_z, face)
  var z_faces_8 = region(grid_z, face)


  -- 1D Region for angle values
  var angle_indices = ispace(int1d, N_angles)
  var [angles] = region(angle_indices, angle)

  -- Partition faces
  -- extra tile required for shared edge
  var [x_tiles] = ispace(int3d, {x = ntx+1, y = nty,   z = ntz  })
  var [y_tiles] = ispace(int3d, {x = ntx,   y = nty+1, z = ntz  })
  var [z_tiles] = ispace(int3d, {x = ntx,   y = nty,   z = ntz+1})

  c.printf("ntx = %d, nty = %d, ntz = %d \n", ntx, nty, ntz)
  c.printf("Nx = %d, Ny = %d, Nz = %d \n", Nx, Ny, Nz)

  -- Partition x_faces_1 private/shared
  -- Partition private and shared into tiles
  
  var [p_x_faces_1] = make_private_partition_x(x_faces_1, x_tiles)
  var [p_y_faces_1] = make_private_partition_y(y_faces_1, y_tiles)
  var [p_z_faces_1] = make_private_partition_z(z_faces_1, z_tiles)

  var [s_x_faces[1]] = make_shared_partition_x(x_faces_1, x_tiles)
  var [s_y_faces[1]] = make_shared_partition_y(y_faces_1, y_tiles)
  var [s_z_faces[1]] = make_shared_partition_z(z_faces_1, z_tiles)

  var [p_x_faces_2] = make_private_partition_x(x_faces_2, x_tiles)
  var [p_y_faces_2] = make_private_partition_y(y_faces_2, y_tiles)
  var [p_z_faces_2] = make_private_partition_z(z_faces_2, z_tiles)

  var [s_x_faces[2]] = make_shared_partition_x(x_faces_2, x_tiles)
  var [s_y_faces[2]] = make_shared_partition_y(y_faces_2, y_tiles)
  var [s_z_faces[2]] = make_shared_partition_z(z_faces_2, z_tiles)

  var [p_x_faces_3] = make_private_partition_x(x_faces_3, x_tiles)
  var [p_y_faces_3] = make_private_partition_y(y_faces_3, y_tiles)
  var [p_z_faces_3] = make_private_partition_z(z_faces_3, z_tiles)

  var [s_x_faces[3]] = make_shared_partition_x(x_faces_3, x_tiles)
  var [s_y_faces[3]] = make_shared_partition_y(y_faces_3, y_tiles)
  var [s_z_faces[3]] = make_shared_partition_z(z_faces_3, z_tiles)

  var [p_x_faces_4] = make_private_partition_x(x_faces_4, x_tiles)
  var [p_y_faces_4] = make_private_partition_y(y_faces_4, y_tiles)
  var [p_z_faces_4] = make_private_partition_z(z_faces_4, z_tiles)

  var [s_x_faces[4]] = make_shared_partition_x(x_faces_4, x_tiles)
  var [s_y_faces[4]] = make_shared_partition_y(y_faces_4, y_tiles)
  var [s_z_faces[4]] = make_shared_partition_z(z_faces_4, z_tiles)

  var [p_x_faces_5] = make_private_partition_x(x_faces_5, x_tiles)
  var [p_y_faces_5] = make_private_partition_y(y_faces_5, y_tiles)
  var [p_z_faces_5] = make_private_partition_z(z_faces_5, z_tiles)

  var [s_x_faces[5]] = make_shared_partition_x(x_faces_5, x_tiles)
  var [s_y_faces[5]] = make_shared_partition_y(y_faces_5, y_tiles)
  var [s_z_faces[5]] = make_shared_partition_z(z_faces_5, z_tiles)

  var [p_x_faces_6] = make_private_partition_x(x_faces_6, x_tiles)
  var [p_y_faces_6] = make_private_partition_y(y_faces_6, y_tiles)
  var [p_z_faces_6] = make_private_partition_z(z_faces_6, z_tiles)

  var [s_x_faces[6]] = make_shared_partition_x(x_faces_6, x_tiles)
  var [s_y_faces[6]] = make_shared_partition_y(y_faces_6, y_tiles)
  var [s_z_faces[6]] = make_shared_partition_z(z_faces_6, z_tiles)

  var [p_x_faces_7] = make_private_partition_x(x_faces_7, x_tiles)
  var [p_y_faces_7] = make_private_partition_y(y_faces_7, y_tiles)
  var [p_z_faces_7] = make_private_partition_z(z_faces_7, z_tiles)

  var [s_x_faces[7]] = make_shared_partition_x(x_faces_7, x_tiles)
  var [s_y_faces[7]] = make_shared_partition_y(y_faces_7, y_tiles)
  var [s_z_faces[7]] = make_shared_partition_z(z_faces_7, z_tiles)

  var [p_x_faces_8] = make_private_partition_x(x_faces_8, x_tiles)
  var [p_y_faces_8] = make_private_partition_y(y_faces_8, y_tiles)
  var [p_z_faces_8] = make_private_partition_z(z_faces_8, z_tiles)

  var [s_x_faces[8]] = make_shared_partition_x(x_faces_8, x_tiles)
  var [s_y_faces[8]] = make_shared_partition_y(y_faces_8, y_tiles)
  var [s_z_faces[8]] = make_shared_partition_z(z_faces_8, z_tiles)


  -- Initialize constant values
  initialize_angles(angles)

  -- Declare variables that would go in the main loop, but for static SPMD
  var [res] = 1.0

end

exports.ComputeRadiationField = rquote


  var t   : int64  = 1

  -- Compute until convergence
  while (res > tol) do

    -- Update the source term (in this problem, isotropic)
    for color in tiles do
      source_term(p_points[color], angles)
    end

    -- Update the grid boundary intensities
    -- Update x faces
    for j = 0, nty do
      for k = 0, ntz do
        -- Avoid empty partitions (index 0 for lo, index ntx for hi)
        west_bound([s_x_faces[1]][{0,j,k}],
                   [s_x_faces[2]][{0,j,k}],
                   [s_x_faces[3]][{0,j,k}],
                   [s_x_faces[4]][{0,j,k}],
                   [s_x_faces[5]][{0,j,k}],
                   [s_x_faces[6]][{0,j,k}],
                   [s_x_faces[7]][{0,j,k}],
                   [s_x_faces[8]][{0,j,k}],
                   angles)

        east_bound([s_x_faces[1]][{ntx,j,k}],
                   [s_x_faces[2]][{ntx,j,k}],
                   [s_x_faces[3]][{ntx,j,k}],
                   [s_x_faces[4]][{ntx,j,k}],
                   [s_x_faces[5]][{ntx,j,k}],
                   [s_x_faces[6]][{ntx,j,k}],
                   [s_x_faces[7]][{ntx,j,k}],
                   [s_x_faces[8]][{ntx,j,k}],
                   angles)
      end
    end

    -- Update y faces
    for i = 0, ntx do
      for k = 0, ntz do
        south_bound([s_y_faces[1]][{i,0,k}],
                    [s_y_faces[2]][{i,0,k}],
                    [s_y_faces[3]][{i,0,k}],
                    [s_y_faces[4]][{i,0,k}],
                    [s_y_faces[5]][{i,0,k}],
                    [s_y_faces[6]][{i,0,k}],
                    [s_y_faces[7]][{i,0,k}],
                    [s_y_faces[8]][{i,0,k}],
                    angles)

        north_bound([s_y_faces[1]][{i,nty,k}],
                    [s_y_faces[2]][{i,nty,k}],
                    [s_y_faces[3]][{i,nty,k}],
                    [s_y_faces[4]][{i,nty,k}],
                    [s_y_faces[5]][{i,nty,k}],
                    [s_y_faces[6]][{i,nty,k}],
                    [s_y_faces[7]][{i,nty,k}],
                    [s_y_faces[8]][{i,nty,k}],
                    angles)
      end
    end

    -- Update z faces
    for i = 0, ntx do
      for j = 0, nty do
        up_bound  ([s_z_faces[1]][{i,j,0}],
                   [s_z_faces[2]][{i,j,0}],
                   [s_z_faces[3]][{i,j,0}], 
                   [s_z_faces[4]][{i,j,0}], 
                   [s_z_faces[5]][{i,j,0}], 
                   [s_z_faces[6]][{i,j,0}], 
                   [s_z_faces[7]][{i,j,0}], 
                   [s_z_faces[8]][{i,j,0}], 
                   angles)

        down_bound([s_z_faces[1]][{i,j,ntz}],
                   [s_z_faces[2]][{i,j,ntz}],
                   [s_z_faces[3]][{i,j,ntz}], 
                   [s_z_faces[4]][{i,j,ntz}], 
                   [s_z_faces[5]][{i,j,ntz}], 
                   [s_z_faces[6]][{i,j,ntz}], 
                   [s_z_faces[7]][{i,j,ntz}], 
                   [s_z_faces[8]][{i,j,ntz}], 
                   angles)
      end
    end

    --Perform the sweep for computing new intensities
    --Quadrant 1 - +x, +y, +z
    for i = 0, ntx do
      for j = 0, nty do
        for k = 0, ntz do
          sweep_1(p_points[{i,j,k}],
                  p_x_faces_1[{i,j,k}], p_y_faces_1[{i,j,k}], p_z_faces_1[{i,j,k}],
                  [s_x_faces[1]][{i,j,k}], [s_x_faces[1]][{i+1,j,k}],
                  [s_y_faces[1]][{i,j,k}], [s_y_faces[1]][{i,j+1,k}],
                  [s_z_faces[1]][{i,j,k}], [s_z_faces[1]][{i,j,k+1}],
                  angles, 1, 1, 1)
        end
      end
    end

    -- Quadrant 2 - +x, +y, -z
    for i = 0, ntx do
      for j = 0, nty do
        for k = ntz-1, -1, -1 do
          sweep_2(p_points[{i,j,k}],
                  p_x_faces_2[{i,j,k}], p_y_faces_2[{i,j,k}], p_z_faces_2[{i,j,k}],
                  [s_x_faces[2]][{i,j,k}], [s_x_faces[2]][{i+1,j,k}],
                  [s_y_faces[2]][{i,j,k}], [s_y_faces[2]][{i,j+1,k}],
                  [s_z_faces[2]][{i,j,k+1}], [s_z_faces[2]][{i,j,k}],
                  angles, 1, 1, -1)
        end
      end
    end

    -- Quadrant 3 - +x, -y, +z
    for i = 0, ntx do
      for j = nty-1, -1, -1 do
        for k = 0, ntz do
          sweep_3(p_points[{i,j,k}],
                  p_x_faces_3[{i,j,k}], p_y_faces_3[{i,j,k}], p_z_faces_3[{i,j,k}],
                  [s_x_faces[3]][{i,j,k}], [s_x_faces[3]][{i+1,j,k}],
                  [s_y_faces[3]][{i,j+1,k}], [s_y_faces[3]][{i,j,k}],
                  [s_z_faces[3]][{i,j,k}], [s_z_faces[3]][{i,j,k+1}],
                  angles, 1, -1, 1)
        end
      end
    end

    -- Quadrant 4 - +x, -y, -z
    for i = 0, ntx do
      for j = nty-1, -1, -1 do
        for k = ntz-1, -1, -1 do
          sweep_4(p_points[{i,j,k}],
                  p_x_faces_4[{i,j,k}], p_y_faces_4[{i,j,k}], p_z_faces_4[{i,j,k}],
                  [s_x_faces[4]][{i,j,k}], [s_x_faces[4]][{i+1,j,k}],
                  [s_y_faces[4]][{i,j+1,k}], [s_y_faces[4]][{i,j,k}],
                  [s_z_faces[4]][{i,j,k+1}], [s_z_faces[4]][{i,j,k}],
                  angles, 1, -1, -1)
        end
      end
    end

    -- Quadrant 5 - -x, +y, +z
    for i = ntx-1, -1, -1 do
      for j = 0, nty do
        for k = 0, ntz do
          sweep_5(p_points[{i,j,k}],
                  p_x_faces_5[{i,j,k}], p_y_faces_5[{i,j,k}], p_z_faces_5[{i,j,k}],
                  [s_x_faces[5]][{i+1,j,k}], [s_x_faces[5]][{i,j,k}],
                  [s_y_faces[5]][{i,j,k}], [s_y_faces[5]][{i,j+1,k}],
                  [s_z_faces[5]][{i,j,k}], [s_z_faces[5]][{i,j,k+1}],
                  angles, -1, 1, 1)
        end
      end
    end

    -- Quadrant 6 - -x, +y, -z
    for i = ntx-1, -1, -1 do
      for j = 0, nty do
        for k = ntz-1, -1, -1 do
          sweep_6(p_points[{i,j,k}],
                  p_x_faces_6[{i,j,k}], p_y_faces_6[{i,j,k}], p_z_faces_6[{i,j,k}],
                  [s_x_faces[6]][{i+1,j,k}], [s_x_faces[6]][{i,j,k}],
                  [s_y_faces[6]][{i,j,k}], [s_y_faces[6]][{i,j+1,k}],
                  [s_z_faces[6]][{i,j,k+1}], [s_z_faces[6]][{i,j,k}],
                  angles, -1, 1, -1)
        end
      end
    end

    -- Quadrant 7 - -x, -y, +z
    for i = ntx-1, -1, -1 do
      for j = nty-1, -1, -1 do
        for k = 0, ntz do
          sweep_7(p_points[{i,j,k}],
                  p_x_faces_7[{i,j,k}], p_y_faces_7[{i,j,k}], p_z_faces_7[{i,j,k}],
                  [s_x_faces[7]][{i+1,j,k}], [s_x_faces[7]][{i,j,k}],
                  [s_y_faces[7]][{i,j+1,k}], [s_y_faces[7]][{i,j,k}],
                  [s_z_faces[7]][{i,j,k}], [s_z_faces[7]][{i,j,k+1}],
                  angles, -1, -1, 1)
        end
      end
    end

    -- Quadrant 8 - -x, -y, -z
    for i = ntx-1, -1, -1 do
      for j = nty-1, -1, -1 do
        for k = ntz-1, -1, -1 do
          sweep_8(p_points[{i,j,k}],
                  p_x_faces_8[{i,j,k}], p_y_faces_8[{i,j,k}], p_z_faces_8[{i,j,k}],
                  [s_x_faces[8]][{i+1,j,k}], [s_x_faces[8]][{i,j,k}],
                  [s_y_faces[8]][{i,j+1,k}], [s_y_faces[8]][{i,j,k}],
                  [s_z_faces[8]][{i,j,k+1}], [s_z_faces[8]][{i,j,k}],
                  angles, -1, -1, -1)
        end
      end
    end

    -- Compute the residual
    res = 0.0
    for color in tiles do
      res += residual(p_points[color])
    end
    res = cmath.sqrt(res)

    -- Update the intensities and the iteration number
    for color in tiles do
      update(p_points[color])
    end

    if (t == 1) then
      c.printf("\n")
      c.printf(" Iteration     Residual         \n")
      c.printf(" ------------------------------ \n")
    end
    c.printf( "   %3d    %.15e \n", t, res)

    t = t + 1

  end

  -- Reduce intensity
  for color in tiles do
    reduce_intensity(p_points[color], angles)
  end
  
  -- Debugging
  write_intensity(points)
  -- print_final_intensities(points)

end

-------------------------------------------------------------------------------
-- MODULE EXPORTS
-------------------------------------------------------------------------------

return exports
end
