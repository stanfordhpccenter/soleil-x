import 'regent'

local A = require 'admiral'

-------------------------------------------------------------------------------
-- MODULE PARAMETERS
-------------------------------------------------------------------------------

return function(radiationRel, NUM_ANGLES,
                NxGlobal, NyGlobal, NzGlobal,
                dxGlobal, dyGlobal, dzGlobal)

local points = radiationRel:regionSymbol()
local p_points = radiationRel:primPartSymbol()
local pointsType = radiationRel:regionType()
local tiles = A.primColors()

local Nx = NxGlobal:varSymbol()
local Ny = NyGlobal:varSymbol()
local Nz = NzGlobal:varSymbol()

local dx = dxGlobal:varSymbol()
local dy = dyGlobal:varSymbol()
local dz = dzGlobal:varSymbol()

local config = A.configSymbol()

-------------------------------------------------------------------------------
-- COMPILE-TIME COMPUTATION
-------------------------------------------------------------------------------

-- C imports

local c     = regentlib.c
local cmath = terralib.includec('math.h')

-- Some math definitions

local min = regentlib.fmin
local max = regentlib.fmax
local pow = regentlib.pow(double)
local pi  = 2.0*cmath.acos(0.0)

-- Quadrature file name
local quad_file = 'LMquads/'..NUM_ANGLES..'.txt'

-- Wall temperatures

local SB = 5.67e-8

-- TODO: Read location of irradiated window

-- Procedure parameters

local tol   = 1e-6   -- solution tolerance
local gamma = 0.5    -- 1 for step differencing, 0.5 for diamond differencing

local terra read_val(f : &c.FILE, val : &double)
  return c.fscanf(f, '%lf\n', &val[0])
end
A.registerFun(read_val, 'read_val')

-------------------------------------------------------------------------------
-- MODULE-LOCAL FIELD SPACES
-------------------------------------------------------------------------------

-- Internal cell values are essentially private,
-- face values are what need to be passed to downstream neighbor
-- Update cell value, then update downstream face values

-- quadrature information
local struct angle {
  xi  : double,
  eta : double,
  mu  : double,
  w   : double,
}
A.registerStruct(angle)

local struct face {
  I : double[NUM_ANGLES],
}
A.registerStruct(face)

-------------------------------------------------------------------------------
-- MODULE-LOCAL TASKS
-------------------------------------------------------------------------------

-- Initialize face values
local task initialize_faces(faces : region(ispace(int3d), face))
where
  reads writes (faces.I)
do
  for m = 0, NUM_ANGLES do
    faces.I[m] = 0.0
  end
end

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
A.registerTask(initialize_angles, 'initialize_angles')

local task make_interior_partition_x_hi(faces : region(ispace(int3d), face),
                                        x_tiles : ispace(int3d),
                                        Nx : int, Ny : int, Nz : int,
                                        ntx : int, nty : int, ntz : int)

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
A.registerTask(make_interior_partition_x_hi, 'make_interior_partition_x_hi')

local task make_interior_partition_x_lo(faces : region(ispace(int3d), face),
                                        x_tiles : ispace(int3d),
                                        Nx : int, Ny : int, Nz : int,
                                        ntx : int, nty : int, ntz : int)

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
A.registerTask(make_interior_partition_x_lo, 'make_interior_partition_x_lo')

local task make_interior_partition_y_hi(faces : region(ispace(int3d), face),
                                        y_tiles : ispace(int3d),
                                        Nx : int, Ny : int, Nz : int,
                                        ntx : int, nty : int, ntz : int)

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
A.registerTask(make_interior_partition_y_hi, 'make_interior_partition_y_hi')

local task make_interior_partition_y_lo(faces : region(ispace(int3d), face),
                                        y_tiles : ispace(int3d),
                                        Nx : int, Ny : int, Nz : int,
                                        ntx : int, nty : int, ntz : int)

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
A.registerTask(make_interior_partition_y_lo, 'make_interior_partition_y_lo')

local task make_interior_partition_z_hi(faces : region(ispace(int3d), face),
                                        z_tiles : ispace(int3d),
                                        Nx : int, Ny : int, Nz : int,
                                        ntx : int, nty : int, ntz : int)

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
A.registerTask(make_interior_partition_z_hi, 'make_interior_partition_z_hi')

local task make_interior_partition_z_lo(faces : region(ispace(int3d), face),
                                        z_tiles : ispace(int3d),
                                        Nx : int, Ny : int, Nz : int,
                                        ntx : int, nty : int, ntz : int)

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
A.registerTask(make_interior_partition_z_lo, 'make_interior_partition_z_lo')

-- Loop over all angles and grid cells to compute the source term
-- for the current iteration.
local task source_term(points : pointsType,
                       angles : region(ispace(int1d), angle),
                       omega : double)
where
  reads (points.{Iiter_1, Iiter_2, Iiter_3, Iiter_4,
                 Iiter_5, Iiter_6, Iiter_7, Iiter_8,
                 Ib, sigma},
         angles.w),
  reads writes (points.S)
do
  for p in points do
    p.S = (1.0-omega) * p.sigma * p.Ib
    for m = 0, NUM_ANGLES do
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
A.registerTask(source_term, 'source_term')

local task west_bound(faces_1 : region(ispace(int3d), face),
                      faces_2 : region(ispace(int3d), face),
                      faces_3 : region(ispace(int3d), face),
                      faces_4 : region(ispace(int3d), face),
                      faces_5 : region(ispace(int3d), face),
                      faces_6 : region(ispace(int3d), face),
                      faces_7 : region(ispace(int3d), face),
                      faces_8 : region(ispace(int3d), face),
                      angles : region(ispace(int1d), angle),
                      emissWest : double,
                      tempWest : double)
where
  reads (angles.{w, xi, eta, mu}),
  reads writes (faces_1.I, faces_2.I, faces_3.I, faces_4.I,
                faces_5.I, faces_6.I, faces_7.I, faces_8.I)
do

  -- Get array bounds

  var limits = faces_1.bounds

  -- Temporary variables

  var reflect : double = 0.0
  var epsw    : double = emissWest
  var Tw      : double = tempWest

  for j = limits.lo.y, limits.hi.y + 1 do
    for k = limits.lo.z, limits.hi.z + 1 do

      -- Calculate reflect

      reflect = 0
      for m = 0, NUM_ANGLES do
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

      for m = 0, NUM_ANGLES do
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
A.registerTask(west_bound, 'west_bound')

local task east_bound(faces_1 : region(ispace(int3d), face),
                      faces_2 : region(ispace(int3d), face),
                      faces_3 : region(ispace(int3d), face),
                      faces_4 : region(ispace(int3d), face),
                      faces_5 : region(ispace(int3d), face),
                      faces_6 : region(ispace(int3d), face),
                      faces_7 : region(ispace(int3d), face),
                      faces_8 : region(ispace(int3d), face),
                      angles : region(ispace(int1d), angle),
                      emissEast : double,
                      tempEast : double)
where
  reads (angles.{w, xi, eta, mu}),
  reads writes (faces_1.I, faces_2.I, faces_3.I, faces_4.I,
                faces_5.I, faces_6.I, faces_7.I, faces_8.I)
do

  -- Get array bounds

  var limits = faces_1.bounds

  -- Temporary variables

  var reflect : double = 0.0
  var epsw    : double = emissEast
  var Tw      : double = tempEast

  for j = limits.lo.y, limits.hi.y + 1 do
    for k = limits.lo.z, limits.hi.z + 1 do

      -- Calculate reflect

      reflect = 0
      for m = 0, NUM_ANGLES do
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

      for m = 0, NUM_ANGLES do
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
A.registerTask(east_bound, 'east_bound')

local task north_bound(faces_1 : region(ispace(int3d), face),
                       faces_2 : region(ispace(int3d), face),
                       faces_3 : region(ispace(int3d), face),
                       faces_4 : region(ispace(int3d), face),
                       faces_5 : region(ispace(int3d), face),
                       faces_6 : region(ispace(int3d), face),
                       faces_7 : region(ispace(int3d), face),
                       faces_8 : region(ispace(int3d), face),
                       angles : region(ispace(int1d), angle),
                       emissNorth : double,
                       tempNorth : double)
where
  reads (angles.{w, xi, eta, mu}),
  reads writes (faces_1.I, faces_2.I, faces_3.I, faces_4.I,
                faces_5.I, faces_6.I, faces_7.I, faces_8.I)
do

  -- Get array bounds

  var limits = faces_1.bounds

  -- Temporary variables

  var reflect : double = 0.0
  var epsw    : double = emissNorth
  var Tw      : double = tempNorth

  for i = limits.lo.x, limits.hi.x + 1 do
    for k = limits.lo.z, limits.hi.z + 1 do

      -- Calculate reflect

      reflect = 0
      for m = 0, NUM_ANGLES do
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

      for m = 0, NUM_ANGLES do
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
A.registerTask(north_bound, 'north_bound')

local task south_bound(faces_1 : region(ispace(int3d), face),
                       faces_2 : region(ispace(int3d), face),
                       faces_3 : region(ispace(int3d), face),
                       faces_4 : region(ispace(int3d), face),
                       faces_5 : region(ispace(int3d), face),
                       faces_6 : region(ispace(int3d), face),
                       faces_7 : region(ispace(int3d), face),
                       faces_8 : region(ispace(int3d), face),
                       angles : region(ispace(int1d), angle),
                       emissSouth : double,
                       tempSouth : double)
where
  reads (angles.{w, xi, eta, mu}),
  reads writes (faces_1.I, faces_2.I, faces_3.I, faces_4.I,
                faces_5.I, faces_6.I, faces_7.I, faces_8.I)
do

  -- Get array bounds

  var limits = faces_1.bounds

  -- Temporary variables

  var reflect : double = 0.0
  var epsw    : double = emissSouth
  var Tw      : double = tempSouth

  for i = limits.lo.x, limits.hi.x + 1 do
    for k = limits.lo.z, limits.hi.z + 1 do

      -- Calculate reflect

      reflect = 0
      for m = 0, NUM_ANGLES do
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

      for m = 0, NUM_ANGLES do
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
A.registerTask(south_bound, 'south_bound')

local task up_bound(faces_1 : region(ispace(int3d), face),
                    faces_2 : region(ispace(int3d), face),
                    faces_3 : region(ispace(int3d), face),
                    faces_4 : region(ispace(int3d), face),
                    faces_5 : region(ispace(int3d), face),
                    faces_6 : region(ispace(int3d), face),
                    faces_7 : region(ispace(int3d), face),
                    faces_8 : region(ispace(int3d), face),
                    angles : region(ispace(int1d), angle),
                    emissUp : double,
                    tempUp : double)
where
  reads (angles.{w, xi, eta, mu}),
  reads writes (faces_1.I, faces_2.I, faces_3.I, faces_4.I,
                faces_5.I, faces_6.I, faces_7.I, faces_8.I)
do

  -- Get array bounds

  var limits = faces_1.bounds

  -- Temporary variables

  var reflect : double = 0.0
  var epsw    : double = emissUp
  var Tw      : double = tempUp

  for i = limits.lo.x, limits.hi.x + 1 do
    for j = limits.lo.y, limits.hi.y + 1 do

      -- Calculate reflect

      reflect = 0
      for m = 0, NUM_ANGLES do
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

      for m = 0, NUM_ANGLES do
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
A.registerTask(up_bound, 'up_bound')

local task down_bound(faces_1 : region(ispace(int3d), face),
                      faces_2 : region(ispace(int3d), face),
                      faces_3 : region(ispace(int3d), face),
                      faces_4 : region(ispace(int3d), face),
                      faces_5 : region(ispace(int3d), face),
                      faces_6 : region(ispace(int3d), face),
                      faces_7 : region(ispace(int3d), face),
                      faces_8 : region(ispace(int3d), face),
                      angles : region(ispace(int1d), angle),
                      emissDown : double,
                      tempDown : double)
where
  reads (angles.{w, xi, eta, mu}),
  reads writes (faces_1.I, faces_2.I, faces_3.I, faces_4.I,
                faces_5.I, faces_6.I, faces_7.I, faces_8.I)
do

  -- Get array bounds

  var limits = faces_1.bounds

  -- Temporary variables

  var reflect : double = 0.0
  var epsw    : double = emissDown
  var Tw      : double = tempDown

  for i = limits.lo.x, limits.hi.x + 1 do
    for j = limits.lo.y, limits.hi.y + 1 do

      -- Calculate reflect

      reflect = 0
      for m = 0, NUM_ANGLES do
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

      for m = 0, NUM_ANGLES do
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
A.registerTask(down_bound, 'down_bound')

local task sweep_1(points : pointsType,
                   x_faces : region(ispace(int3d), face),
                   y_faces : region(ispace(int3d), face),
                   z_faces : region(ispace(int3d), face),
                   ghost_x_faces : region(ispace(int3d), face),
                   ghost_y_faces : region(ispace(int3d), face),
                   ghost_z_faces : region(ispace(int3d), face),
                   angles : region(ispace(int1d), angle),
                   xi : int64, eta : int64, mu : int64,
                   dx : double, dy : double, dz : double)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         ghost_x_faces.I, ghost_y_faces.I, ghost_z_faces.I),
  reads writes(points.I_1, x_faces.I, y_faces.I, z_faces.I)
do
  var dAx = dy*dz;
  var dAy = dx*dz;
  var dAz = dx*dy;
  var dV = dx*dy*dz;

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
  for m = 0, NUM_ANGLES do

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

            var ghost_x_limits = ghost_x_faces.bounds
            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x then
              upwind_x_value = ghost_x_faces[{ghost_x_limits.hi.x,j,k}].I[m]
            elseif indx > x_faces.bounds.hi.x then
              upwind_x_value = ghost_x_faces[{ghost_x_limits.lo.x,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            var ghost_y_limits = ghost_y_faces.bounds
            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y then
              upwind_y_value = ghost_y_faces[{i,ghost_y_limits.hi.y,k}].I[m]
            elseif indy > y_faces.bounds.hi.y then
              upwind_y_value = ghost_y_faces[{i,ghost_y_limits.lo.y,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var ghost_z_limits = ghost_z_faces.bounds
            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z then
              upwind_z_value = ghost_z_faces[{i,j,ghost_z_limits.hi.z}].I[m]
            elseif indz > z_faces.bounds.hi.z then
              upwind_z_value = ghost_z_faces[{i,j,ghost_z_limits.lo.z}].I[m]
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

            x_faces[{indx+dindx, j, k}].I[m] = (points[{i,j,k}].I_1[m] - (1-gamma)*upwind_x_value)/gamma
            y_faces[{i, indy+dindy, k}].I[m] = (points[{i,j,k}].I_1[m] - (1-gamma)*upwind_y_value)/gamma
            z_faces[{i, j, indz+dindz}].I[m] = (points[{i,j,k}].I_1[m] - (1-gamma)*upwind_z_value)/gamma
          end
        end
      end
    end
  end
end
A.registerTask(sweep_1, 'sweep_1')

local task sweep_2(points : pointsType,
                   x_faces : region(ispace(int3d), face),
                   y_faces : region(ispace(int3d), face),
                   z_faces : region(ispace(int3d), face),
                   ghost_x_faces : region(ispace(int3d), face),
                   ghost_y_faces : region(ispace(int3d), face),
                   ghost_z_faces : region(ispace(int3d), face),
                   angles : region(ispace(int1d), angle),
                   xi : int64, eta : int64, mu : int64,
                   dx : double, dy : double, dz : double)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         ghost_x_faces.I, ghost_y_faces.I, ghost_z_faces.I),
  reads writes(points.I_2, x_faces.I, y_faces.I, z_faces.I)
do
  var dAx = dy*dz;
  var dAy = dx*dz;
  var dAz = dx*dy;
  var dV = dx*dy*dz;

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
  for m = 0, NUM_ANGLES do

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

            var ghost_x_limits = ghost_x_faces.bounds
            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x then
              upwind_x_value = ghost_x_faces[{ghost_x_limits.hi.x,j,k}].I[m]
            elseif indx > x_faces.bounds.hi.x then
              upwind_x_value = ghost_x_faces[{ghost_x_limits.lo.x,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            var ghost_y_limits = ghost_y_faces.bounds
            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y then
              upwind_y_value = ghost_y_faces[{i,ghost_y_limits.hi.y,k}].I[m]
            elseif indy > y_faces.bounds.hi.y then
              upwind_y_value = ghost_y_faces[{i,ghost_y_limits.lo.y,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var ghost_z_limits = ghost_z_faces.bounds
            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z then
              upwind_z_value = ghost_z_faces[{i,j,ghost_z_limits.hi.z}].I[m]
            elseif indz > z_faces.bounds.hi.z then
              upwind_z_value = ghost_z_faces[{i,j,ghost_z_limits.lo.z}].I[m]
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

            x_faces[{indx+dindx, j, k}].I[m] = (points[{i,j,k}].I_2[m] - (1-gamma)*upwind_x_value)/gamma
            y_faces[{i, indy+dindy, k}].I[m] = (points[{i,j,k}].I_2[m] - (1-gamma)*upwind_y_value)/gamma
            z_faces[{i, j, indz+dindz}].I[m] = (points[{i,j,k}].I_2[m] - (1-gamma)*upwind_z_value)/gamma
          end
        end
      end
    end
  end
end
A.registerTask(sweep_2, 'sweep_2')

local task sweep_3(points : pointsType,
                   x_faces : region(ispace(int3d), face),
                   y_faces : region(ispace(int3d), face),
                   z_faces : region(ispace(int3d), face),
                   ghost_x_faces : region(ispace(int3d), face),
                   ghost_y_faces : region(ispace(int3d), face),
                   ghost_z_faces : region(ispace(int3d), face),
                   angles : region(ispace(int1d), angle),
                   xi : int64, eta : int64, mu : int64,
                   dx : double, dy : double, dz : double)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         ghost_x_faces.I, ghost_y_faces.I, ghost_z_faces.I),
  reads writes(points.I_3, x_faces.I, y_faces.I, z_faces.I)
do
  var dAx = dy*dz;
  var dAy = dx*dz;
  var dAz = dx*dy;
  var dV = dx*dy*dz;

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
  for m = 0, NUM_ANGLES do

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

            var ghost_x_limits = ghost_x_faces.bounds
            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x then
              upwind_x_value = ghost_x_faces[{ghost_x_limits.hi.x,j,k}].I[m]
            elseif indx > x_faces.bounds.hi.x then
              upwind_x_value = ghost_x_faces[{ghost_x_limits.lo.x,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            var ghost_y_limits = ghost_y_faces.bounds
            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y then
              upwind_y_value = ghost_y_faces[{i,ghost_y_limits.hi.y,k}].I[m]
            elseif indy > y_faces.bounds.hi.y then
              upwind_y_value = ghost_y_faces[{i,ghost_y_limits.lo.y,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var ghost_z_limits = ghost_z_faces.bounds
            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z then
              upwind_z_value = ghost_z_faces[{i,j,ghost_z_limits.hi.z}].I[m]
            elseif indz > z_faces.bounds.hi.z then
              upwind_z_value = ghost_z_faces[{i,j,ghost_z_limits.lo.z}].I[m]
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

            x_faces[{indx+dindx, j, k}].I[m] = (points[{i,j,k}].I_3[m] - (1-gamma)*upwind_x_value)/gamma
            y_faces[{i, indy+dindy, k}].I[m] = (points[{i,j,k}].I_3[m] - (1-gamma)*upwind_y_value)/gamma
            z_faces[{i, j, indz+dindz}].I[m] = (points[{i,j,k}].I_3[m] - (1-gamma)*upwind_z_value)/gamma
          end
        end
      end
    end
  end
end
A.registerTask(sweep_3, 'sweep_3')

local task sweep_4(points : pointsType,
                   x_faces : region(ispace(int3d), face),
                   y_faces : region(ispace(int3d), face),
                   z_faces : region(ispace(int3d), face),
                   ghost_x_faces : region(ispace(int3d), face),
                   ghost_y_faces : region(ispace(int3d), face),
                   ghost_z_faces : region(ispace(int3d), face),
                   angles : region(ispace(int1d), angle),
                   xi : int64, eta : int64, mu : int64,
                   dx : double, dy : double, dz : double)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         ghost_x_faces.I, ghost_y_faces.I, ghost_z_faces.I),
  reads writes(points.I_4, x_faces.I, y_faces.I, z_faces.I)
do
  var dAx = dy*dz;
  var dAy = dx*dz;
  var dAz = dx*dy;
  var dV = dx*dy*dz;

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
  for m = 0, NUM_ANGLES do

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

            var ghost_x_limits = ghost_x_faces.bounds
            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x then
              upwind_x_value = ghost_x_faces[{ghost_x_limits.hi.x,j,k}].I[m]
            elseif indx > x_faces.bounds.hi.x then
              upwind_x_value = ghost_x_faces[{ghost_x_limits.lo.x,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            var ghost_y_limits = ghost_y_faces.bounds
            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y then
              upwind_y_value = ghost_y_faces[{i,ghost_y_limits.hi.y,k}].I[m]
            elseif indy > y_faces.bounds.hi.y then
              upwind_y_value = ghost_y_faces[{i,ghost_y_limits.lo.y,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var ghost_z_limits = ghost_z_faces.bounds
            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z then
              upwind_z_value = ghost_z_faces[{i,j,ghost_z_limits.hi.z}].I[m]
            elseif indz > z_faces.bounds.hi.z then
              upwind_z_value = ghost_z_faces[{i,j,ghost_z_limits.lo.z}].I[m]
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

            x_faces[{indx+dindx, j, k}].I[m] = (points[{i,j,k}].I_4[m] - (1-gamma)*upwind_x_value)/gamma
            y_faces[{i, indy+dindy, k}].I[m] = (points[{i,j,k}].I_4[m] - (1-gamma)*upwind_y_value)/gamma
            z_faces[{i, j, indz+dindz}].I[m] = (points[{i,j,k}].I_4[m] - (1-gamma)*upwind_z_value)/gamma
          end
        end
      end
    end
  end
end
A.registerTask(sweep_4, 'sweep_4')

local task sweep_5(points : pointsType,
                   x_faces : region(ispace(int3d), face),
                   y_faces : region(ispace(int3d), face),
                   z_faces : region(ispace(int3d), face),
                   ghost_x_faces : region(ispace(int3d), face),
                   ghost_y_faces : region(ispace(int3d), face),
                   ghost_z_faces : region(ispace(int3d), face),
                   angles : region(ispace(int1d), angle),
                   xi : int64, eta : int64, mu : int64,
                   dx : double, dy : double, dz : double)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         ghost_x_faces.I, ghost_y_faces.I, ghost_z_faces.I),
  reads writes(points.I_5, x_faces.I, y_faces.I, z_faces.I)
do
  var dAx = dy*dz;
  var dAy = dx*dz;
  var dAz = dx*dy;
  var dV = dx*dy*dz;

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
  for m = 0, NUM_ANGLES do

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

            var ghost_x_limits = ghost_x_faces.bounds
            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x then
              upwind_x_value = ghost_x_faces[{ghost_x_limits.hi.x,j,k}].I[m]
            elseif indx > x_faces.bounds.hi.x then
              upwind_x_value = ghost_x_faces[{ghost_x_limits.lo.x,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            var ghost_y_limits = ghost_y_faces.bounds
            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y then
              upwind_y_value = ghost_y_faces[{i,ghost_y_limits.hi.y,k}].I[m]
            elseif indy > y_faces.bounds.hi.y then
              upwind_y_value = ghost_y_faces[{i,ghost_y_limits.lo.y,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var ghost_z_limits = ghost_z_faces.bounds
            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z then
              upwind_z_value = ghost_z_faces[{i,j,ghost_z_limits.hi.z}].I[m]
            elseif indz > z_faces.bounds.hi.z then
              upwind_z_value = ghost_z_faces[{i,j,ghost_z_limits.lo.z}].I[m]
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

            x_faces[{indx+dindx, j, k}].I[m] = (points[{i,j,k}].I_5[m] - (1-gamma)*upwind_x_value)/gamma
            y_faces[{i, indy+dindy, k}].I[m] = (points[{i,j,k}].I_5[m] - (1-gamma)*upwind_y_value)/gamma
            z_faces[{i, j, indz+dindz}].I[m] = (points[{i,j,k}].I_5[m] - (1-gamma)*upwind_z_value)/gamma
          end
        end
      end
    end
  end
end
A.registerTask(sweep_5, 'sweep_5')

local task sweep_6(points : pointsType,
                   x_faces : region(ispace(int3d), face),
                   y_faces : region(ispace(int3d), face),
                   z_faces : region(ispace(int3d), face),
                   ghost_x_faces : region(ispace(int3d), face),
                   ghost_y_faces : region(ispace(int3d), face),
                   ghost_z_faces : region(ispace(int3d), face),
                   angles : region(ispace(int1d), angle),
                   xi : int64, eta : int64, mu : int64,
                   dx : double, dy : double, dz : double)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         ghost_x_faces.I, ghost_y_faces.I, ghost_z_faces.I),
  reads writes(points.I_6, x_faces.I, y_faces.I, z_faces.I)
do
  var dAx = dy*dz;
  var dAy = dx*dz;
  var dAz = dx*dy;
  var dV = dx*dy*dz;

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
  for m = 0, NUM_ANGLES do

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

            var ghost_x_limits = ghost_x_faces.bounds
            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x then
              upwind_x_value = ghost_x_faces[{ghost_x_limits.hi.x,j,k}].I[m]
            elseif indx > x_faces.bounds.hi.x then
              upwind_x_value = ghost_x_faces[{ghost_x_limits.lo.x,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            var ghost_y_limits = ghost_y_faces.bounds
            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y then
              upwind_y_value = ghost_y_faces[{i,ghost_y_limits.hi.y,k}].I[m]
            elseif indy > y_faces.bounds.hi.y then
              upwind_y_value = ghost_y_faces[{i,ghost_y_limits.lo.y,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var ghost_z_limits = ghost_z_faces.bounds
            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z then
              upwind_z_value = ghost_z_faces[{i,j,ghost_z_limits.hi.z}].I[m]
            elseif indz > z_faces.bounds.hi.z then
              upwind_z_value = ghost_z_faces[{i,j,ghost_z_limits.lo.z}].I[m]
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

            x_faces[{indx+dindx, j, k}].I[m] = (points[{i,j,k}].I_6[m] - (1-gamma)*upwind_x_value)/gamma
            y_faces[{i, indy+dindy, k}].I[m] = (points[{i,j,k}].I_6[m] - (1-gamma)*upwind_y_value)/gamma
            z_faces[{i, j, indz+dindz}].I[m] = (points[{i,j,k}].I_6[m] - (1-gamma)*upwind_z_value)/gamma
          end
        end
      end
    end
  end
end
A.registerTask(sweep_6, 'sweep_6')

local task sweep_7(points : pointsType,
                   x_faces : region(ispace(int3d), face),
                   y_faces : region(ispace(int3d), face),
                   z_faces : region(ispace(int3d), face),
                   ghost_x_faces : region(ispace(int3d), face),
                   ghost_y_faces : region(ispace(int3d), face),
                   ghost_z_faces : region(ispace(int3d), face),
                   angles : region(ispace(int1d), angle),
                   xi : int64, eta : int64, mu : int64,
                   dx : double, dy : double, dz : double)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         ghost_x_faces.I, ghost_y_faces.I, ghost_z_faces.I),
  reads writes(points.I_7, x_faces.I, y_faces.I, z_faces.I)
do
  var dAx = dy*dz;
  var dAy = dx*dz;
  var dAz = dx*dy;
  var dV = dx*dy*dz;

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
  for m = 0, NUM_ANGLES do

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

            var ghost_x_limits = ghost_x_faces.bounds
            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x then
              upwind_x_value = ghost_x_faces[{ghost_x_limits.hi.x,j,k}].I[m]
            elseif indx > x_faces.bounds.hi.x then
              upwind_x_value = ghost_x_faces[{ghost_x_limits.lo.x,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            var ghost_y_limits = ghost_y_faces.bounds
            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y then
              upwind_y_value = ghost_y_faces[{i,ghost_y_limits.hi.y,k}].I[m]
            elseif indy > y_faces.bounds.hi.y then
              upwind_y_value = ghost_y_faces[{i,ghost_y_limits.lo.y,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var ghost_z_limits = ghost_z_faces.bounds
            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z then
              upwind_z_value = ghost_z_faces[{i,j,ghost_z_limits.hi.z}].I[m]
            elseif indz > z_faces.bounds.hi.z then
              upwind_z_value = ghost_z_faces[{i,j,ghost_z_limits.lo.z}].I[m]
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

            x_faces[{indx+dindx, j, k}].I[m] = (points[{i,j,k}].I_7[m] - (1-gamma)*upwind_x_value)/gamma
            y_faces[{i, indy+dindy, k}].I[m] = (points[{i,j,k}].I_7[m] - (1-gamma)*upwind_y_value)/gamma
            z_faces[{i, j, indz+dindz}].I[m] = (points[{i,j,k}].I_7[m] - (1-gamma)*upwind_z_value)/gamma
          end
        end
      end
    end
  end
end
A.registerTask(sweep_7, 'sweep_7')

local task sweep_8(points : pointsType,
                   x_faces : region(ispace(int3d), face),
                   y_faces : region(ispace(int3d), face),
                   z_faces : region(ispace(int3d), face),
                   ghost_x_faces : region(ispace(int3d), face),
                   ghost_y_faces : region(ispace(int3d), face),
                   ghost_z_faces : region(ispace(int3d), face),
                   angles : region(ispace(int1d), angle),
                   xi : int64, eta : int64, mu : int64,
                   dx : double, dy : double, dz : double)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         ghost_x_faces.I, ghost_y_faces.I, ghost_z_faces.I),
  reads writes(points.I_8, x_faces.I, y_faces.I, z_faces.I)
do
  var dAx = dy*dz;
  var dAy = dx*dz;
  var dAz = dx*dy;
  var dV = dx*dy*dz;

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
  for m = 0, NUM_ANGLES do

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

            var ghost_x_limits = ghost_x_faces.bounds
            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x then
              upwind_x_value = ghost_x_faces[{ghost_x_limits.hi.x,j,k}].I[m]
            elseif indx > x_faces.bounds.hi.x then
              upwind_x_value = ghost_x_faces[{ghost_x_limits.lo.x,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            var ghost_y_limits = ghost_y_faces.bounds
            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y then
              upwind_y_value = ghost_y_faces[{i,ghost_y_limits.hi.y,k}].I[m]
            elseif indy > y_faces.bounds.hi.y then
              upwind_y_value = ghost_y_faces[{i,ghost_y_limits.lo.y,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var ghost_z_limits = ghost_z_faces.bounds
            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z then
              upwind_z_value = ghost_z_faces[{i,j,ghost_z_limits.hi.z}].I[m]
            elseif indz > z_faces.bounds.hi.z then
              upwind_z_value = ghost_z_faces[{i,j,ghost_z_limits.lo.z}].I[m]
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

            x_faces[{indx+dindx, j, k}].I[m] = (points[{i,j,k}].I_8[m] - (1-gamma)*upwind_x_value)/gamma
            y_faces[{i, indy+dindy, k}].I[m] = (points[{i,j,k}].I_8[m] - (1-gamma)*upwind_y_value)/gamma
            z_faces[{i, j, indz+dindz}].I[m] = (points[{i,j,k}].I_8[m] - (1-gamma)*upwind_z_value)/gamma
          end
        end
      end
    end
  end
end
A.registerTask(sweep_8, 'sweep_8')

-- Compute the residual after each iteration and return the value.
local task residual(points : pointsType, Nx : int, Ny : int, Nz : int)
where
  reads (points.{I_1, I_2, I_3, I_4, I_5, I_6, I_7, I_8,
                 Iiter_1, Iiter_2, Iiter_3, Iiter_4,
                 Iiter_5, Iiter_6, Iiter_7, Iiter_8})
do
  var res : double = 0.0

  for p in points do
    for m = 0, NUM_ANGLES do

      if p.I_1[m] > 0 then
        res += (1.0/(Nx*Ny*Nz*(NUM_ANGLES)))
          * pow((p.I_1[m]-p.Iiter_1[m]),2.0)
          / pow((p.I_1[m]),2.0)
      end

      if p.I_2[m] > 0 then
        res += (1.0/(Nx*Ny*Nz*(NUM_ANGLES)))
          * pow((p.I_2[m]-p.Iiter_2[m]),2.0)
          / pow((p.I_2[m]),2.0)
      end

      if p.I_3[m] > 0 then
        res += (1.0/(Nx*Ny*Nz*(NUM_ANGLES)))
          * pow((p.I_3[m]-p.Iiter_3[m]),2.0)
          / pow((p.I_3[m]),2.0)
      end

      if p.I_4[m] > 0 then
        res += (1.0/(Nx*Ny*Nz*(NUM_ANGLES)))
          * pow((p.I_4[m]-p.Iiter_4[m]),2.0)
          / pow((p.I_4[m]),2.0)
      end

      if p.I_5[m] > 0 then
        res += (1.0/(Nx*Ny*Nz*(NUM_ANGLES)))
          * pow((p.I_5[m]-p.Iiter_5[m]),2.0)
          / pow((p.I_5[m]),2.0)
      end

      if p.I_6[m] > 0 then
        res += (1.0/(Nx*Ny*Nz*(NUM_ANGLES)))
          * pow((p.I_6[m]-p.Iiter_6[m]),2.0)
          / pow((p.I_6[m]),2.0)
      end

      if p.I_7[m] > 0 then
        res += (1.0/(Nx*Ny*Nz*(NUM_ANGLES)))
          * pow((p.I_7[m]-p.Iiter_7[m]),2.0)
          / pow((p.I_7[m]),2.0)
      end

      if p.I_8[m] > 0 then
        res += (1.0/(Nx*Ny*Nz*(NUM_ANGLES)))
          * pow((p.I_8[m]-p.Iiter_8[m]),2.0)
          / pow((p.I_8[m]),2.0)
      end

    end
  end

  return res
end
A.registerTask(residual, 'residual')

-- Update the intensity before moving to the next iteration.
local task update(points : pointsType)
where
  reads (points.{I_1, I_2, I_3, I_4, I_5, I_6, I_7, I_8}),
  reads writes (points.{Iiter_1, Iiter_2, Iiter_3, Iiter_4,
                        Iiter_5, Iiter_6, Iiter_7, Iiter_8})
do
  for p in points do
    for m = 0, NUM_ANGLES do
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
A.registerTask(update, 'update')

-- Reduce the intensity to summation over all angles
local task reduce_intensity(points : pointsType,
                            angles : region(ispace(int1d), angle))
where
  reads (points.{I_1, I_2, I_3, I_4, I_5, I_6, I_7, I_8}, angles.w),
  reads writes (points.G)
do
  for p in points do
    for m = 0, NUM_ANGLES do
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
A.registerTask(reduce_intensity, 'reduce_intensity')

-------------------------------------------------------------------------------
-- EXPORTED QUOTES
-------------------------------------------------------------------------------

local exports = {}

-- Symbols shared between quotes
local ntx = regentlib.newsymbol('ntx')
local nty = regentlib.newsymbol('nty')
local ntz = regentlib.newsymbol('ntz')

local angles = regentlib.newsymbol('angles')
local x_tiles = regentlib.newsymbol('x_tiles')
local y_tiles = regentlib.newsymbol('y_tiles')
local z_tiles = regentlib.newsymbol('z_tiles')

local p_x_faces_1 = regentlib.newsymbol('p_x_faces_1')
local p_y_faces_1 = regentlib.newsymbol('p_y_faces_1')
local p_z_faces_1 = regentlib.newsymbol('p_z_faces_1')

local p_x_faces_2 = regentlib.newsymbol('p_x_faces_2')
local p_y_faces_2 = regentlib.newsymbol('p_y_faces_2')
local p_z_faces_2 = regentlib.newsymbol('p_z_faces_2')

local p_x_faces_3 = regentlib.newsymbol('p_x_faces_3')
local p_y_faces_3 = regentlib.newsymbol('p_y_faces_3')
local p_z_faces_3 = regentlib.newsymbol('p_z_faces_3')

local p_x_faces_4 = regentlib.newsymbol('p_x_faces_4')
local p_y_faces_4 = regentlib.newsymbol('p_y_faces_4')
local p_z_faces_4 = regentlib.newsymbol('p_z_faces_4')

local p_x_faces_5 = regentlib.newsymbol('p_x_faces_5')
local p_y_faces_5 = regentlib.newsymbol('p_y_faces_5')
local p_z_faces_5 = regentlib.newsymbol('p_z_faces_5')

local p_x_faces_6 = regentlib.newsymbol('p_x_faces_6')
local p_y_faces_6 = regentlib.newsymbol('p_y_faces_6')
local p_z_faces_6 = regentlib.newsymbol('p_z_faces_6')

local p_x_faces_7 = regentlib.newsymbol('p_x_faces_7')
local p_y_faces_7 = regentlib.newsymbol('p_y_faces_7')
local p_z_faces_7 = regentlib.newsymbol('p_z_faces_7')

local p_x_faces_8 = regentlib.newsymbol('p_x_faces_8')
local p_y_faces_8 = regentlib.newsymbol('p_y_faces_8')
local p_z_faces_8 = regentlib.newsymbol('p_z_faces_8')

local res = regentlib.newsymbol(double, 'res')

exports.InitModule = rquote

  var [ntx] = config.Grid.xTiles
  var [nty] = config.Grid.yTiles
  var [ntz] = config.Grid.zTiles

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
  var angle_indices = ispace(int1d, NUM_ANGLES)
  var [angles] = region(angle_indices, angle)

  -- Partition faces
  -- extra tile required for ghost
  var [x_tiles] = ispace(int3d, {x = ntx+1, y = nty,   z = ntz  })
  var [y_tiles] = ispace(int3d, {x = ntx,   y = nty+1, z = ntz  })
  var [z_tiles] = ispace(int3d, {x = ntx,   y = nty,   z = ntz+1})

  var [p_x_faces_1] = make_interior_partition_x_lo(x_faces_1, x_tiles, Nx, Ny, Nz, ntx, nty, ntz)
  var [p_y_faces_1] = make_interior_partition_y_lo(y_faces_1, y_tiles, Nx, Ny, Nz, ntx, nty, ntz)
  var [p_z_faces_1] = make_interior_partition_z_lo(z_faces_1, z_tiles, Nx, Ny, Nz, ntx, nty, ntz)

  var [p_x_faces_2] = make_interior_partition_x_lo(x_faces_2, x_tiles, Nx, Ny, Nz, ntx, nty, ntz)
  var [p_y_faces_2] = make_interior_partition_y_lo(y_faces_2, y_tiles, Nx, Ny, Nz, ntx, nty, ntz)
  var [p_z_faces_2] = make_interior_partition_z_hi(z_faces_2, z_tiles, Nx, Ny, Nz, ntx, nty, ntz)

  var [p_x_faces_3] = make_interior_partition_x_lo(x_faces_3, x_tiles, Nx, Ny, Nz, ntx, nty, ntz)
  var [p_y_faces_3] = make_interior_partition_y_hi(y_faces_3, y_tiles, Nx, Ny, Nz, ntx, nty, ntz)
  var [p_z_faces_3] = make_interior_partition_z_lo(z_faces_3, z_tiles, Nx, Ny, Nz, ntx, nty, ntz)

  var [p_x_faces_4] = make_interior_partition_x_lo(x_faces_4, x_tiles, Nx, Ny, Nz, ntx, nty, ntz)
  var [p_y_faces_4] = make_interior_partition_y_hi(y_faces_4, y_tiles, Nx, Ny, Nz, ntx, nty, ntz)
  var [p_z_faces_4] = make_interior_partition_z_hi(z_faces_4, z_tiles, Nx, Ny, Nz, ntx, nty, ntz)

  var [p_x_faces_5] = make_interior_partition_x_hi(x_faces_5, x_tiles, Nx, Ny, Nz, ntx, nty, ntz)
  var [p_y_faces_5] = make_interior_partition_y_lo(y_faces_5, y_tiles, Nx, Ny, Nz, ntx, nty, ntz)
  var [p_z_faces_5] = make_interior_partition_z_lo(z_faces_5, z_tiles, Nx, Ny, Nz, ntx, nty, ntz)

  var [p_x_faces_6] = make_interior_partition_x_hi(x_faces_6, x_tiles, Nx, Ny, Nz, ntx, nty, ntz)
  var [p_y_faces_6] = make_interior_partition_y_lo(y_faces_6, y_tiles, Nx, Ny, Nz, ntx, nty, ntz)
  var [p_z_faces_6] = make_interior_partition_z_hi(z_faces_6, z_tiles, Nx, Ny, Nz, ntx, nty, ntz)

  var [p_x_faces_7] = make_interior_partition_x_hi(x_faces_7, x_tiles, Nx, Ny, Nz, ntx, nty, ntz)
  var [p_y_faces_7] = make_interior_partition_y_hi(y_faces_7, y_tiles, Nx, Ny, Nz, ntx, nty, ntz)
  var [p_z_faces_7] = make_interior_partition_z_lo(z_faces_7, z_tiles, Nx, Ny, Nz, ntx, nty, ntz)

  var [p_x_faces_8] = make_interior_partition_x_hi(x_faces_8, x_tiles, Nx, Ny, Nz, ntx, nty, ntz)
  var [p_y_faces_8] = make_interior_partition_y_hi(y_faces_8, y_tiles, Nx, Ny, Nz, ntx, nty, ntz)
  var [p_z_faces_8] = make_interior_partition_z_hi(z_faces_8, z_tiles, Nx, Ny, Nz, ntx, nty, ntz)

  -- Initialize face values
  initialize_faces(x_faces_1)
  initialize_faces(x_faces_2)
  initialize_faces(x_faces_3)
  initialize_faces(x_faces_4)
  initialize_faces(x_faces_5)
  initialize_faces(x_faces_6)
  initialize_faces(x_faces_7)
  initialize_faces(x_faces_8)

  initialize_faces(y_faces_1)
  initialize_faces(y_faces_2)
  initialize_faces(y_faces_3)
  initialize_faces(y_faces_4)
  initialize_faces(y_faces_5)
  initialize_faces(y_faces_6)
  initialize_faces(y_faces_7)
  initialize_faces(y_faces_8)

  initialize_faces(z_faces_1)
  initialize_faces(z_faces_2)
  initialize_faces(z_faces_3)
  initialize_faces(z_faces_4)
  initialize_faces(z_faces_5)
  initialize_faces(z_faces_6)
  initialize_faces(z_faces_7)
  initialize_faces(z_faces_8)

  -- Initialize constant values
  initialize_angles(angles)

  -- Declare variables that would go in the main loop, but for static SPMD
  var [res] = 1.0

end

exports.ComputeRadiationField = rquote

  var t : int64  = 1
  var omega = config.Radiation.qs/(config.Radiation.qa+config.Radiation.qs)

  -- Compute until convergence
  while (res > tol) do

    -- Update the source term (in this problem, isotropic)
    for color in tiles do
      source_term(p_points[color], angles, omega)
    end

    -- Update the grid boundary intensities
    -- Update x faces
    for j = 0, nty do
      for k = 0, ntz do
        -- Avoid empty partitions (index 0 for lo, index ntx for hi)
        west_bound(p_x_faces_1[{1  ,j,k}],
                   p_x_faces_2[{1  ,j,k}],
                   p_x_faces_3[{1  ,j,k}],
                   p_x_faces_4[{1  ,j,k}],
                   p_x_faces_5[{0  ,j,k}],
                   p_x_faces_6[{0  ,j,k}],
                   p_x_faces_7[{0  ,j,k}],
                   p_x_faces_8[{0  ,j,k}],
                   angles,
                   config.Radiation.emissWest,
                   config.Radiation.tempWest)

        east_bound(p_x_faces_1[{ntx  ,j,k}],
                   p_x_faces_2[{ntx  ,j,k}],
                   p_x_faces_3[{ntx  ,j,k}],
                   p_x_faces_4[{ntx  ,j,k}],
                   p_x_faces_5[{ntx-1,j,k}],
                   p_x_faces_6[{ntx-1,j,k}],
                   p_x_faces_7[{ntx-1,j,k}],
                   p_x_faces_8[{ntx-1,j,k}],
                   angles,
                   config.Radiation.emissEast,
                   config.Radiation.tempEast)
      end
    end

    -- Update y faces
    for i = 0, ntx do
      for k = 0, ntz do
        south_bound(p_y_faces_1[{i,1,  k}],
                    p_y_faces_2[{i,1,  k}],
                    p_y_faces_3[{i,0,  k}],
                    p_y_faces_4[{i,0,  k}],
                    p_y_faces_5[{i,1,  k}],
                    p_y_faces_6[{i,1,  k}],
                    p_y_faces_7[{i,0,  k}],
                    p_y_faces_8[{i,0,  k}],
                    angles,
                    config.Radiation.emissSouth,
                    config.Radiation.tempSouth)

        north_bound(p_y_faces_1[{i,nty,k}],
                    p_y_faces_2[{i,nty,k}],
                    p_y_faces_3[{i,nty-1,k}],
                    p_y_faces_4[{i,nty-1,k}],
                    p_y_faces_5[{i,nty,k}],
                    p_y_faces_6[{i,nty,k}],
                    p_y_faces_7[{i,nty-1,k}],
                    p_y_faces_8[{i,nty-1,k}],
                    angles,
                    config.Radiation.emissNorth,
                    config.Radiation.tempNorth)
      end
    end

    -- Update z faces
    for i = 0, ntx do
      for j = 0, nty do
        up_bound  (p_z_faces_1[{i,j,1  }],
                   p_z_faces_2[{i,j,0  }],
                   p_z_faces_3[{i,j,1  }],
                   p_z_faces_4[{i,j,0  }],
                   p_z_faces_5[{i,j,1  }],
                   p_z_faces_6[{i,j,0  }],
                   p_z_faces_7[{i,j,1  }],
                   p_z_faces_8[{i,j,0  }],
                   angles,
                   config.Radiation.emissUp,
                   config.Radiation.tempUp)

        down_bound(p_z_faces_1[{i,j,ntz}],
                   p_z_faces_2[{i,j,ntz-1}],
                   p_z_faces_3[{i,j,ntz}],
                   p_z_faces_4[{i,j,ntz-1}],
                   p_z_faces_5[{i,j,ntz}],
                   p_z_faces_6[{i,j,ntz-1}],
                   p_z_faces_7[{i,j,ntz}],
                   p_z_faces_8[{i,j,ntz-1}],
                   angles,
                   config.Radiation.emissDown,
                   config.Radiation.tempDown)
      end
    end

    -- Perform the sweep for computing new intensities
    -- Quadrant 1 - +x, +y, +z
    for i = 0, ntx do
      for j = 0, nty do
        for k = 0, ntz do
          sweep_1(p_points[{i,j,k}],
                  p_x_faces_1[{i+1,j,k}], p_y_faces_1[{i,j+1,k}], p_z_faces_1[{i,j,k+1}],
                  p_x_faces_1[{i,  j,k}], p_y_faces_1[{i,j,  k}], p_z_faces_1[{i,j,k  }],
                  angles, 1, 1, 1, dx, dy, dz)
        end
      end
    end

    -- Quadrant 2 - +x, +y, -z
    for i = 0, ntx do
      for j = 0, nty do
        for k = ntz-1, -1, -1 do
          sweep_2(p_points[{i,j,k}],
                  p_x_faces_2[{i+1,j,k}], p_y_faces_2[{i,j+1,k}], p_z_faces_2[{i,j,k  }],
                  p_x_faces_2[{i,  j,k}], p_y_faces_2[{i,j,  k}], p_z_faces_2[{i,j,k+1}],
                  angles, 1, 1, -1, dx, dy, dz)
        end
      end
    end

    -- Quadrant 3 - +x, -y, +z
    for i = 0, ntx do
      for j = nty-1, -1, -1 do
        for k = 0, ntz do
          sweep_3(p_points[{i,j,k}],
                  p_x_faces_3[{i+1,j,k}], p_y_faces_3[{i,j,  k}], p_z_faces_3[{i,j,k+1}],
                  p_x_faces_3[{i,  j,k}], p_y_faces_3[{i,j+1,k}], p_z_faces_3[{i,j,k  }],
                  angles, 1, -1, 1, dx, dy, dz)
        end
      end
    end

    -- Quadrant 4 - +x, -y, -z
    for i = 0, ntx do
      for j = nty-1, -1, -1 do
        for k = ntz-1, -1, -1 do
          sweep_4(p_points[{i,j,k}],
                  p_x_faces_4[{i+1,j,k}], p_y_faces_4[{i,j,  k}], p_z_faces_4[{i,j,k  }],
                  p_x_faces_4[{i,  j,k}], p_y_faces_4[{i,j+1,k}], p_z_faces_4[{i,j,k+1}],
                  angles, 1, -1, -1, dx, dy, dz)
        end
      end
    end

    -- Quadrant 5 - -x, +y, +z
    for i = ntx-1, -1, -1 do
      for j = 0, nty do
        for k = 0, ntz do
          sweep_5(p_points[{i,j,k}],
                  p_x_faces_5[{i,  j,k}], p_y_faces_5[{i,j+1,k}], p_z_faces_5[{i,j,k+1}],
                  p_x_faces_5[{i+1,j,k}], p_y_faces_5[{i,j,  k}], p_z_faces_5[{i,j,k  }],
                  angles, -1, 1, 1, dx, dy, dz)
        end
      end
    end

    -- Quadrant 6 - -x, +y, -z
    for i = ntx-1, -1, -1 do
      for j = 0, nty do
        for k = ntz-1, -1, -1 do
          sweep_6(p_points[{i,j,k}],
                  p_x_faces_6[{i,  j,k}], p_y_faces_6[{i,j+1,k}], p_z_faces_6[{i,j,k  }],
                  p_x_faces_6[{i+1,j,k}], p_y_faces_6[{i,j,  k}], p_z_faces_6[{i,j,k+1}],
                  angles, -1, 1, -1, dx, dy, dz)
        end
      end
    end

    -- Quadrant 7 - -x, -y, +z
    for i = ntx-1, -1, -1 do
      for j = nty-1, -1, -1 do
        for k = 0, ntz do
          sweep_7(p_points[{i,j,k}],
                  p_x_faces_7[{i,  j,k}], p_y_faces_7[{i,j,  k}], p_z_faces_7[{i,j,k+1}],
                  p_x_faces_7[{i+1,j,k}], p_y_faces_7[{i,j+1,k}], p_z_faces_7[{i,j,k  }],
                  angles, -1, -1, 1, dx, dy, dz)
        end
      end
    end

    -- Quadrant 8 - -x, -y, -z
    for i = ntx-1, -1, -1 do
      for j = nty-1, -1, -1 do
        for k = ntz-1, -1, -1 do
          sweep_8(p_points[{i,j,k}],
                  p_x_faces_8[{i,  j,k}], p_y_faces_8[{i,j,  k}], p_z_faces_8[{i,j,k  }],
                  p_x_faces_8[{i+1,j,k}], p_y_faces_8[{i,j+1,k}], p_z_faces_8[{i,j,k+1}],
                  angles, -1, -1, -1, dx, dy, dz)
        end
      end
    end

    -- Compute the residual
    res = 0.0
    for color in tiles do
      res += residual(p_points[color], Nx, Ny, Nz)
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

end

-------------------------------------------------------------------------------
-- MODULE EXPORTS
-------------------------------------------------------------------------------

return exports
end
