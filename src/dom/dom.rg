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

local task make_interior_partition_x_hi(faces : region(ispace(int3d), x_face),
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

local task make_interior_partition_x_lo(faces : region(ispace(int3d), x_face),
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

local task make_interior_partition_y_hi(faces : region(ispace(int3d), y_face),
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

local task make_interior_partition_y_lo(faces : region(ispace(int3d), y_face),
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

local task make_interior_partition_z_hi(faces : region(ispace(int3d), z_face),
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

local task make_interior_partition_z_lo(faces : region(ispace(int3d), z_face),
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

  var limits = faces.bounds

  -- Temporary variables

  var reflect : double = 0.0
  var epsw    : double = emiss_west
  var Tw      : double = T_west

  for j = limits.lo.y, limits.hi.y + 1 do
    for k = limits.lo.z, limits.hi.z + 1 do

      -- Calculate reflect

      reflect = 0
      for m = 0, N_angles do
        if angles[m].xi < 0 then
          var face_value : double = 0.0
          if angles[m].eta > 0 and angles[m].mu > 0 then
            face_value = faces_5[{limits.lo.x,j,k}].I[m]
          elseif angles[m].eta > 0 and angles[m].mu < 0 then
            face_value = faces_6[{limits.lo.x,j,k}].I[m]
          elseif angles[m].eta < 0 and angles[m].mu > 0 then
            face_value = faces_7[{limits.lo.x,j,k}].I[m]
          else
            face_value = faces_8[{limits.lo.x,j,k}].I[m]
          end
          reflect += (1.0-epsw)/pi*angles[m].w*cmath.fabs(angles[m].xi)*face_value
        end
      end

      -- Set Ifx values using reflect

      for m = 0, N_angles do
        if angles[m].xi > 0 then
          var value : double = epsw*SB*pow(Tw,4.0)/pi + reflect

          if angles[m].eta > 0 and angles[m].mu > 0 then
            faces_1[{limits.lo.x,j,k}].I[m] = value
          elseif angles[m].eta > 0 and angles[m].mu < 0 then
            faces_2[{limits.lo.x,j,k}].I[m] = value
          elseif angles[m].eta < 0 and angles[m].mu > 0 then
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

  var limits = faces.bounds

  -- Temporary variables

  var reflect : double = 0.0
  var epsw    : double = emiss_east
  var Tw      : double = T_east

  for j = limits.lo.y, limits.hi.y + 1 do
    for k = limits.lo.z, limits.hi.z + 1 do

      -- Calculate reflect

      reflect = 0
      for m = 0, N_angles do
        if angles[m].xi > 0 then
          var face_value : double = 0.0
          if angles[m].eta > 0 and angles[m].mu > 0 then
            face_value = faces_1[{limits.hi.x,j,k}].I[m]
          elseif angles[m].eta > 0 and angles[m].mu < 0 then
            face_value = faces_2[{limits.hi.x,j,k}].I[m]
          elseif angles[m].eta < 0 and angles[m].mu > 0 then
            face_value = faces_3[{limits.hi.x,j,k}].I[m]
          else
            face_value = faces_4[{limits.hi.x,j,k}].I[m]
          end
          reflect += (1.0-epsw)/pi*angles[m].w*angles[m].xi*face_value
        end
      end

      -- Set Ifx values using reflect

      for m = 0, N_angles do
        if angles[m].xi < 0 then
          var value : double = epsw*SB*pow(Tw,4.0)/pi + reflect

          if angles[m].eta > 0 and angles[m].mu > 0 then
            faces_5[{limits.hi.x,j,k}].I[m] = value
          elseif angles[m].eta > 0 and angles[m].mu < 0 then
            faces_6[{limits.hi.x,j,k}].I[m] = value
          elseif angles[m].eta < 0 and angles[m].mu > 0 then
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

  var limits = faces.bounds

  -- Temporary variables

  var reflect : double = 0.0
  var epsw    : double = emiss_north
  var Tw      : double = T_north

  for i = limits.lo.x, limits.hi.x + 1 do
    for k = limits.lo.z, limits.hi.z + 1 do

      -- Calculate reflect

      reflect = 0
      for m = 0, N_angles do
        if angles[m].eta > 0 then
          var face_value : double = 0.0
          if angles[m].xi > 0 and angles[m].mu > 0 then
            face_value = faces_1[{i,limits.hi.y,k}].I[m]
          elseif angles[m].xi > 0 and angles[m].mu < 0 then
            face_value = faces_2[{i,limits.hi.y,k}].I[m]
          elseif angles[m].xi < 0 and angles[m].mu > 0 then
            face_value = faces_5[{i,limits.hi.y,k}].I[m]
          else
            face_value = faces_6[{i,limits.hi.y,k}].I[m]
          end
          reflect += (1.0-epsw)/pi*angles[m].w*angles[m].eta*face_value
        end
      end

      -- Set Ify values using reflect

      for m = 0, N_angles do
        if angles[m].eta < 0 then
          var value : double = epsw*SB*pow(Tw,4.0)/pi + reflect

          if angles[m].xi > 0 and angles[m].mu > 0 then
            faces_3[{i,limits.hi.y,k}].I[m] = value
          elseif angles[m].xi > 0 and angles[m].mu < 0 then
            faces_4[{i,limits.hi.y,k}].I[m] = value
          elseif angles[m].xi < 0 and angles[m].mu > 0 then
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

  var limits = faces.bounds

  -- Temporary variables

  var reflect : double = 0.0
  var epsw    : double = emiss_south
  var Tw      : double = T_south

  for i = limits.lo.x, limits.hi.x + 1 do
    for k = limits.lo.z, limits.hi.z + 1 do

      -- Calculate reflect

      reflect = 0
      for m = 0, N_angles do
        if angles[m].eta < 0 then
          var face_value : double = 0.0
          if angles[m].xi > 0 and angles[m].mu > 0 then
            face_value = faces_3[{i,limits.lo.y,k}].I[m]
          elseif angles[m].xi > 0 and angles[m].mu < 0 then
            face_value = faces_4[{i,limits.lo.y,k}].I[m]
          elseif angles[m].xi < 0 and angles[m].mu > 0 then
            face_value = faces_7[{i,limits.lo.y,k}].I[m]
          else
            face_value = faces_8[{i,limits.lo.y,k}].I[m]
          end
          reflect += (1.0-epsw)/pi*angles[m].w*cmath.fabs(angles[m].eta)*face_value
        end
      end

      -- Set Ify values using reflect

      for m = 0, N_angles do
        if angles[m].eta > 0 then
          var value : double = epsw*SB*pow(Tw,4.0)/pi + reflect

          if angles[m].xi > 0 and angles[m].mu > 0 then
            faces_1[{i,limits.lo.y,k}].I[m] = value
          elseif angles[m].xi > 0 and angles[m].mu < 0 then
            faces_2[{i,limits.lo.y,k}].I[m] = value
          elseif angles[m].xi < 0 and angles[m].mu > 0 then
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

  var limits = faces.bounds

  -- Temporary variables

  var reflect : double = 0.0
  var epsw    : double = emiss_up
  var Tw      : double = T_up

  for i = limits.lo.x, limits.hi.x + 1 do
    for j = limits.lo.y, limits.hi.y + 1 do

      -- Calculate reflect

      reflect = 0
      for m = 0, N_angles do
        if angles[m].mu < 0 then
          var face_value : double = 0.0
          if angles[m].xi > 0 and angles[m].eta > 0 then
            face_value = faces_2[{i,j,limits.lo.z}].I[m]
          elseif angles[m].xi > 0 and angles[m].eta < 0 then
            face_value = faces_4[{i,j,limits.lo.z}].I[m]
          elseif angles[m].xi < 0 and angles[m].eta > 0 then
            face_value = faces_6[{i,j,limits.lo.z}].I[m]
          else
            face_value = faces_8[{i,j,limits.lo.z}].I[m]
          end
          reflect += (1.0-epsw)/pi*angles[m].w*cmath.fabs(angles[m].mu)*face_value
        end
      end

      -- Set Ifz values using reflect

      for m = 0, N_angles do
        if angles[m].mu > 0 then
          var value : double = epsw*SB*pow(Tw,4.0)/pi + reflect

          if angles[m].xi > 0 and angles[m].eta > 0 then
            faces_1[{i,j,limits.lo.z}].I[m] = value
          elseif angles[m].xi > 0 and angles[m].eta < 0 then
            faces_3[{i,j,limits.lo.z}].I[m] = value
          elseif angles[m].xi < 0 and angles[m].eta > 0 then
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

  var limits = faces.bounds

  -- Temporary variables

  var reflect : double = 0.0
  var epsw    : double = emiss_down
  var Tw      : double = T_down

  for i = limits.lo.x, limits.hi.x + 1 do
    for j = limits.lo.y, limits.hi.y + 1 do

      -- Calculate reflect

      reflect = 0
      for m = 0, N_angles do
        if angles[m].mu > 0 then
          var face_value : double = 0.0
          if angles[m].xi > 0 and angles[m].eta > 0 then
            face_value = faces_1[{i,j,limits.hi.z}].I[m]
          elseif angles[m].xi > 0 and angles[m].eta < 0 then
            face_value = faces_3[{i,j,limits.hi.z}].I[m]
          elseif angles[m].xi < 0 and angles[m].eta > 0 then
            face_value = faces_5[{i,j,limits.hi.z}].I[m]
          else
            face_value = faces_7[{i,j,limits.hi.z}].I[m]
          end
          reflect += (1.0-epsw)/pi*angles[m].w*angles[m].mu*face_value
        end
      end

      -- Set Ifz values using reflect

      for m = 0, N_angles do
        if angles[m].mu < 0 then
          var value : double = epsw*SB*pow(Tw,4.0)/pi + reflect

          if angles[m].xi > 0 and angles[m].eta > 0 then
            faces_2[{i,j,limits.hi.z}].I[m] = value
          elseif angles[m].xi > 0 and angles[m].eta < 0 then
            faces_4[{i,j,limits.hi.z}].I[m] = value
          elseif angles[m].xi < 0 and angles[m].eta > 0 then
            faces_6[{i,j,limits.hi.z}].I[m] = value
          else
            faces_8[{i,j,limits.hi.z}].I[m] = value
          end
        end
      end

    end
  end

end

-- single sweep since only 1 array with 4 regions
local task sweep(points : pointsType,
                   x_faces : region(ispace(int3d), face),
                   y_faces : region(ispace(int3d), face),
                   z_faces : region(ispace(int3d), face),
                   ghost_x_faces : region(ispace(int3d), face),
                   ghost_y_faces : region(ispace(int3d), face),
                   ghost_z_faces : region(ispace(int3d), face),
                   angles : region(ispace(int1d), angle),
                   startx : int64, endx : int64, dindx : int64,
                   starty : int64, endy : int64, dindy : int64,
                   startz : int64, endz : int64, dindz : int64,
                   xi : int64, eta : int64, mu : int64)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         ghost_x_faces.I, ghost_y_faces.I, ghost_z_faces.I),
  reads writes(points.I, x_faces.I, y_faces.I, z_faces.I)
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

  if xi <= 0 then
    dindx = -1
    startx = limits.hi.x
    endx = limits.lo.x - 1
  end

  if eta <= 0 then
    dindx = -1
    startx = limits.hi.y
    endx = limits.lo.y - 1
  end

  if mu <= 0 then
    dindx = -1
    startx = limits.hi.z
    endx = limits.lo.z - 1
  end


  -- Outer loop over all angles.
  for m = 0, N_angles do

    -- todo: if angle_values[m].xi > 0 and angle_values[m].eta > 0 and angle_values[m].mu > 0 then

    -- Use our direction and increments for the sweep.

    for k = startz,endz,dindz do
      for j = starty,endy,dindy do
        for i = startx,endx,dindx do

          -- indx and indy are the upwind indices
          var indx : int64 = i - min(dindx,0)
          var indy : int64 = j - min(dindy,0)
          var indz : int64 = k - min(dindz,0)

          -- Determine if necessary to use ghost partition

          var ghost_x_limits = ghost_x_faces.bounds
          var upwind_x_value : double = 0.0
          if indx < x_faces.bounds.lo.x then
            upwind_x_value = ghost_x_faces[{ghost_x_limits.hi.x,j,k}].I[m]
          else if indx > x_faces.bounds.hi.x then
            upwind_x_value = ghost_x_faces[{i,j,ghost_x_limits.lo.x}].I[m]
          else
            upwind_x_value = x_faces[{indx,j,k}].I[m]
          end

          var ghost_y_limits = ghost_y_faces.bounds
          var upwind_y_value : double = 0.0
          if indy < y_faces.bounds.lo.y then
            upwind_y_value = ghost_y_faces[{i,ghost_y_limits.hi.y,k}].I[m]
          else if indy > y_faces.bounds.hi.y then
            upwind_y_value = ghost_y_faces[{i,j,ghost_y_limits.lo.y}].I[m]
          else
            upwind_y_value = y_faces[{i,indy,k}].I[m]
          end

          var ghost_z_limits = ghost_z_faces.bounds
          var upwind_z_value : double = 0.0
          if indz < z_faces.bounds.lo.z then
            upwind_z_value = ghost_z_faces[{i,j,ghost_z_limits.hi.z}].I[m]
          else if indz > z_faces.bounds.hi.z then
            upwind_z_value = ghost_z_faces[{i,j,ghost_z_limits.lo.z}].I[m]
          else
            upwind_z_value = z_faces[{i,j,indz}].I[m]
          end

          -- Integrate to compute cell-centered value of I.

          points[{i,j,k}].I[m] = (points[{i,j,k}].S * dV
                                      + cmath.fabs(angles[m].xi) * dAx * upwind_x_value/gamma
                                      + cmath.fabs(angles[m].eta) * dAy * upwind_y_value/gamma
                                      + cmath.fabs(angles[m].mu) * dAz * upwind_z_value/gamma)
            /(points[{i,j,k}].sigma * dV
                + cmath.fabs(angles[m].xi) * dAx/gamma
                + cmath.fabs(angles[m].eta) * dAy/gamma
                + cmath.fabs(angles[m].mu) * dAz/gamma)

          -- Compute intensities on downwind faces

          x_faces[{indx+dindx, j, k}].I[m] = (points[{i,j,k}].I[m] - (1-gamma)*upwind_x_value)/gamma
          y_faces[{i, indy+dindy, k}].I[m] = (points[{i,j,k}].I[m] - (1-gamma)*upwind_y_value)/gamma
          z_faces[{i, j, indz+dindz}].I[m] = (points[{i,j,k}].I[m] - (1-gamma)*upwind_z_value)/gamma
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

-------------------------------------------------------------------------------
-- EXPORTED QUOTES
-------------------------------------------------------------------------------

local exports = {}

-- Symbols shared between quotes
local angles = regentlib.newsymbol('angles')
local x_tiles = regentlib.newsymbol('x_tiles')
local y_tiles = regentlib.newsymbol('y_tiles')
local z_tiles = regentlib.newsymbol('z_tiles')
local p_x_faces_hi = regentlib.newsymbol('p_x_faces_hi')
local p_y_faces_hi = regentlib.newsymbol('p_y_faces_hi')
local p_z_faces_hi = regentlib.newsymbol('p_z_faces_hi')
local p_x_faces_lo = regentlib.newsymbol('p_x_faces_lo')
local p_y_faces_lo = regentlib.newsymbol('p_y_faces_lo')
local p_z_faces_lo = regentlib.newsymbol('p_z_faces_lo')
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
  -- extra tile required for ghost
  var [x_tiles] = ispace(int3d, {x = ntx+1, y = nty,   z = ntz  })
  var [y_tiles] = ispace(int3d, {x = ntx,   y = nty+1, z = ntz  })
  var [z_tiles] = ispace(int3d, {x = ntx,   y = nty,   z = ntz+1})

  var [p_x_faces_1] = make_interior_partition_x_lo(x_faces, x_tiles)
  var [p_y_faces_1] = make_interior_partition_y_lo(y_faces, y_tiles)
  var [p_z_faces_1] = make_interior_partition_z_lo(z_faces, z_tiles)

  var [p_x_faces_2] = make_interior_partition_x_lo(x_faces, x_tiles)
  var [p_y_faces_2] = make_interior_partition_y_lo(y_faces, y_tiles)
  var [p_z_faces_2] = make_interior_partition_z_hi(z_faces, z_tiles)

  var [p_x_faces_3] = make_interior_partition_x_lo(x_faces, x_tiles)
  var [p_y_faces_3] = make_interior_partition_y_hi(y_faces, y_tiles)
  var [p_z_faces_3] = make_interior_partition_z_lo(z_faces, z_tiles)

  var [p_x_faces_4] = make_interior_partition_x_lo(x_faces, x_tiles)
  var [p_y_faces_4] = make_interior_partition_y_hi(y_faces, y_tiles)
  var [p_z_faces_4] = make_interior_partition_z_hi(z_faces, z_tiles)

  var [p_x_faces_5] = make_interior_partition_x_hi(x_faces, x_tiles)
  var [p_y_faces_5] = make_interior_partition_y_lo(y_faces, y_tiles)
  var [p_z_faces_5] = make_interior_partition_z_lo(z_faces, z_tiles)

  var [p_x_faces_6] = make_interior_partition_x_hi(x_faces, x_tiles)
  var [p_y_faces_6] = make_interior_partition_y_lo(y_faces, y_tiles)
  var [p_z_faces_6] = make_interior_partition_z_hi(z_faces, z_tiles)

  var [p_x_faces_7] = make_interior_partition_x_hi(x_faces, x_tiles)
  var [p_y_faces_7] = make_interior_partition_y_hi(y_faces, y_tiles)
  var [p_z_faces_7] = make_interior_partition_z_lo(z_faces, z_tiles)

  var [p_x_faces_8] = make_interior_partition_x_hi(x_faces, x_tiles)
  var [p_y_faces_8] = make_interior_partition_y_hi(y_faces, y_tiles)
  var [p_z_faces_8] = make_interior_partition_z_hi(z_faces, z_tiles)

  

  -- Initialize constant values
  initialize_angles(angles)

  -- Declare variables that would go in the main loop, but for static SPMD
  var [res] = 1.0

end

exports.ComputeRadiationField = rquote

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
        west_bound(p_x_faces_hi[{0  ,j,k}], angles)
        east_bound(p_x_faces_lo[{ntx,j,k}], angles)
      end
    end
    -- Update y faces
    for i = 0, ntx do
      for k = 0, ntz do
        south_bound(p_y_faces_hi[{i,0,  k}], angles)
        north_bound(p_y_faces_lo[{i,nty,k}], angles)
      end
    end
    -- Update z faces
    for i = 0, ntx do
      for j = 0, nty do
        up_bound  (p_z_faces_hi[{i,j,0  }], angles)
        down_bound(p_z_faces_lo[{i,j,ntz}], angles)
      end
    end

    -- Perform the sweep for computing new intensities
    -- Quadrant 1 - +x, +y, +z
    for i = 0, ntx do
      for j = 0, nty do
        for k = 0, ntz do
          sweep(p_points[{i,j,k}],
                  p_x_faces_1[{i+1,j,k}], p_y_faces_1[{i,j+1,k}], p_z_faces_1[{i,j,k+1}],
                  p_x_faces_1[{i,  j,k}], p_y_faces_1[{i,j,  k}], p_z_faces_1[{i,j,k  }],
                  angles)
        end
      end
    end
    -- Quadrant 2 - +x, +y, -z
    for i = 0, ntx do
      for j = 0, nty do
        for k = ntz-1, -1, -1 do
          sweep(p_points[{i,j,k}],
                  p_x_faces_2[{i+1,j,k}], p_y_faces_2[{i,j+1,k}], p_z_faces_2[{i,j,k  }],
                  p_x_faces_2[{i,  j,k}], p_y_faces_2[{i,j,  k}], p_z_faces_2[{i,j,k+1}],
                  angles)
        end
      end
    end
    -- Quadrant 3 - +x, -y, +z
    for i = 0, ntx do
      for j = nty-1, -1, -1 do
        for k = 0, ntz do
          sweep(p_points[{i,j,k}],
                  p_x_faces_3[{i+1,j,k}], p_y_faces_3[{i,j,  k}], p_z_faces_3[{i,j,k+1}],
                  p_x_faces_3[{i,  j,k}], p_y_faces_3[{i,j+1,k}], p_z_faces_3[{i,j,k  }],
                  angles)
        end
      end
    end
    -- Quadrant 4 - +x, -y, -z
    for i = 0, ntx do
      for j = nty-1, -1, -1 do
        for k = ntz-1, -1, -1 do
          sweep(p_points[{i,j,k}],
                  p_x_faces_4[{i+1,j,k}], p_y_faces_4[{i,j,  k}], p_z_faces_4[{i,j,k  }],
                  p_x_faces_4[{i,  j,k}], p_y_faces_4[{i,j+1,k}], p_z_faces_4[{i,j,k+1}],
                  angles)
        end
      end
    end
    -- Quadrant 5 - -x, +y, +z
    for i = ntx-1, -1, -1 do
      for j = 0, nty do
        for k = 0, ntz do
          sweep(p_points[{i,j,k}],
                  p_x_faces_5[{i,  j,k}], p_y_faces_5[{i,j+1,k}], p_z_faces_5[{i,j,k+1}],
                  p_x_faces_5[{i+1,j,k}], p_y_faces_5[{i,j,  k}], p_z_faces_5[{i,j,k  }],
                  angles)
        end
      end
    end
    -- Quadrant 6 - -x, +y, -z
    for i = ntx-1, -1, -1 do
      for j = 0, nty do
        for k = ntz-1, -1, -1 do
          sweep(p_points[{i,j,k}],
                  p_x_faces_6[{i,  j,k}], p_y_faces_6[{i,j+1,k}], p_z_faces_6[{i,j,k  }],
                  p_x_faces_6[{i+1,j,k}], p_y_faces_6[{i,j,  k}], p_z_faces_6[{i,j,k+1}],
                  angles)
        end
      end
    end
    -- Quadrant 7 - -x, -y, +z
    for i = ntx-1, -1, -1 do
      for j = nty-1, -1, -1 do
        for k = 0, ntz do
          sweep(p_points[{i,j,k}],
                  p_x_faces_6[{i,  j,k}], p_y_faces_hi[{i,j,  k}], p_z_faces_lo[{i,j,k+1}],
                  p_x_faces_hi[{i+1,j,k}], p_y_faces_hi[{i,j+1,k}], p_z_faces_lo[{i,j,k  }],
                  angles)
        end
      end
    end
    -- Quadrant 8 - -x, -y, -z
    for i = ntx-1, -1, -1 do
      for j = nty-1, -1, -1 do
        for k = ntz-1, -1, -1 do
          sweep_8(p_points[{i,j,k}],
                  p_x_faces_hi[{i,  j,k}], p_y_faces_hi[{i,j,  k}], p_z_faces_hi[{i,j,k  }],
                  p_x_faces_hi[{i+1,j,k}], p_y_faces_hi[{i,j+1,k}], p_z_faces_hi[{i,j,k+1}],
                  angles)
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

  end

  -- Reduce intensity
  for color in tiles do
    reduce_intensity(p_points[color], angles)
  end

end

-------------------------------------------------------------------------------
-- MODULE EXPORTS
-------------------------------------------------------------------------------

return exports
end
