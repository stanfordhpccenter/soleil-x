import 'regent'

-------------------------------------------------------------------------------
-- MODULE PARAMETERS
-------------------------------------------------------------------------------

return function(NUM_ANGLES, pointsFSpace, Config) local MODULE = {}

-------------------------------------------------------------------------------
-- COMPILE-TIME COMPUTATION
-------------------------------------------------------------------------------

-- C imports

local C = regentlib.c

local MAPPER = terralib.includec("soleil_mapper.h")
local UTIL = require 'util'

-- Math imports

local fabs = regentlib.fabs(double)
local max = regentlib.fmax
local min = regentlib.fmin
local pow = regentlib.pow(double)
local sqrt = regentlib.sqrt(double)

-- Math Constants

local pi = 3.1415926535898

-- Wall temperatures

local SB = 5.67e-8

-- Procedure parameters

local tol   = 1e-6   -- solution tolerance
local gamma = 0.5    -- 1 for step differencing, 0.5 for diamond differencing

local terra open_quad_file() : &C.FILE
  var fname = [&int8](C.malloc(256))
  C.snprintf(fname, 256, ['%s/src/LMquads/'..NUM_ANGLES..'.txt'],
             C.getenv('SOLEIL_DIR'))
  var f = C.fopen(fname, 'rb')
  if f == nil then
    var stderr = C.fdopen(2, 'w')
    C.fprintf(stderr, 'Error opening angle file %s\n', fname)
    C.fflush(stderr)
    C.exit(1)
  end
  C.free(fname)
  return f
end

local terra read_double(f : &C.FILE) : double
  var val : double
  if C.fscanf(f, '%lf\n', &val) < 1 then
    C.printf('Error while reading angle file\n')
    C.exit(1)
  end
  return val
end

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

local struct face {
  I : double[NUM_ANGLES],
  is_private : int1d, -- 1 = private, 0 = shared
  tile : int3d,
  diagonal : int1d,
}

local struct tile_info {
  diagonal : int1d,
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

  var f = open_quad_file()

  read_double(f) -- gets rid of num angles

  for a in angles do
    a.xi = read_double(f)
  end

  for a in angles do
    a.eta = read_double(f)
  end

  for a in angles do
    a.mu = read_double(f)
  end

  for a in angles do
    a.w = read_double(f)
  end

  C.fclose(f)

end

local __demand(__inline)
task ite(b : bool, x : int, y : int)
  var res = y
  if b then res = x end
  return res
end

-- Coloring example: 6x6 square, 2x2 tiling, +1-1 direction
--
-- points: 6x6 square
--
-- tile:
--   (0,1)(0,1)(0,1)(1,1)(1,1)(1,1)
--   (0,1)(0,1)(0,1)(1,1)(1,1)(1,1)
-- ^ (0,1)(0,1)(0,1)(1,1)(1,1)(1,1)
-- | (0,0)(0,0)(0,0)(1,0)(1,0)(1,0)
-- y (0,0)(0,0)(0,0)(1,0)(1,0)(1,0)
-- | (0,0)(0,0)(0,0)(1,0)(1,0)(1,0)
--   --x-->
--
-- x-faces: 7x6 square
--
-- is_private:
-- 0110110
-- 0110110
-- 0110110
-- 0110110
-- 0110110
-- 0110110
--
-- private x-faces:
-- diagonal: tile:
-- _00_11_   _____(0,1)(0,1)_____(1,1)(1,1)_____
-- _00_11_   _____(0,1)(0,1)_____(1,1)(1,1)_____
-- _00_11_   _____(0,1)(0,1)_____(1,1)(1,1)_____
-- _11_22_   _____(0,0)(0,0)_____(1,0)(1,0)_____
-- _11_22_   _____(0,0)(0,0)_____(1,0)(1,0)_____
-- _11_22_   _____(0,0)(0,0)_____(1,0)(1,0)_____
--
-- shared x-faces:
-- diagonal: tile:
-- 0__1__2   (0,1)__________(1,1)__________(2,1)
-- 0__1__2   (0,1)__________(1,1)__________(2,1)
-- 0__1__2   (0,1)__________(1,1)__________(2,1)
-- 1__2__3   (0,0)__________(1,0)__________(2,0)
-- 1__2__3   (0,0)__________(1,0)__________(2,0)
-- 1__2__3   (0,0)__________(1,0)__________(2,0)
--
-- y-faces: 6x7 square
--
-- is_private:
-- 111111
-- 000000
-- 000000
-- 111111
-- 000000
-- 000000
-- 111111
--
-- private y-faces:
-- diagonal: tile:
-- ______    ______________________________
-- 000111    (0,1)(0,1)(0,1)(1,1)(1,1)(1,1)
-- 000111    (0,1)(0,1)(0,1)(1,1)(1,1)(1,1)
-- ______    ______________________________
-- 111222    (0,0)(0,0)(0,0)(1,0)(1,0)(1,0)
-- 111222    (0,0)(0,0)(0,0)(1,0)(1,0)(1,0)
-- ______    ______________________________
--
-- shared y-faces:
-- diagonal: tile:
-- 000111    (0,2)(0,2)(0,2)(1,2)(1,2)(1,2)
-- ______    ______________________________
-- ______    ______________________________
-- 111222    (0,1)(0,1)(0,1)(1,1)(1,1)(1,1)
-- ______    ______________________________
-- ______    ______________________________
-- 222333    (0,0)(0,0)(0,0)(1,0)(1,0)(1,0)

local task color_faces(faces : region(ispace(int3d), face),
                       Nx : int, Ny : int, Nz : int,
                       ntx : int, nty : int, ntz : int,
                       dimension : int, sweepDir : bool[3])
where
  writes (faces.{is_private, diagonal, tile})
do
  for idx in faces do
    faces[idx].is_private = 1
    var x_tile = idx.x / (Nx/ntx)
    var y_tile = idx.y / (Ny/nty)
    var z_tile = idx.z / (Nz/ntz)
    faces[idx].tile = int3d{x_tile, y_tile, z_tile}
    if dimension == 0 then
      if idx.x % (Nx/ntx) == 0 then
        faces[idx].is_private = 0
        if not sweepDir[0] then x_tile -= 1 end
      end
    elseif dimension == 1 then
      if idx.y % (Ny/nty) == 0 then
        faces[idx].is_private = 0
        if not sweepDir[1] then y_tile -= 1 end
      end
    elseif dimension == 2 then
      if idx.z % (Nz/ntz) == 0 then
        faces[idx].is_private = 0
        if not sweepDir[2] then z_tile -= 1 end
      end
    else regentlib.assert(false, '') end
    faces[idx].diagonal =
      ite(sweepDir[0], x_tile, ntx-1-x_tile) +
      ite(sweepDir[1], y_tile, nty-1-y_tile) +
      ite(sweepDir[2], z_tile, ntz-1-z_tile)
  end
end

local task fill_tile_info(r_tiles : region(ispace(int3d), tile_info),
                          sweepDir : bool[3])
where
  writes (r_tiles.diagonal)
do
  var ntx = r_tiles.bounds.hi.x + 1
  var nty = r_tiles.bounds.hi.y + 1
  var ntz = r_tiles.bounds.hi.z + 1
  for t in r_tiles do
    r_tiles[t].diagonal =
      ite(sweepDir[0], t.x, ntx-1-t.x) +
      ite(sweepDir[1], t.y, nty-1-t.y) +
      ite(sweepDir[2], t.z, ntz-1-t.z)
  end
end

-- Loop over all angles and grid cells to compute the source term
-- for the current iteration.
local task source_term(points : region(ispace(int3d), pointsFSpace),
                       angles : region(ispace(int1d), angle),
                       omega : double)
where
  reads (points.{Iiter_1, Iiter_2, Iiter_3, Iiter_4,
                 Iiter_5, Iiter_6, Iiter_7, Iiter_8,
                 Ib, sigma},
         angles.w),
  reads writes (points.S)
do
  __demand(__openmp)
  for p in points do
    p.S = (1.0-omega) * p.sigma * p.Ib
    for m = 0, NUM_ANGLES do
      p.S += omega * p.sigma/(4.0*pi) * angles[m].w
          * (p.Iiter_1[m]
           + p.Iiter_2[m]
           + p.Iiter_3[m]
           + p.Iiter_4[m]
           + p.Iiter_5[m]
           + p.Iiter_6[m]
           + p.Iiter_7[m]
           + p.Iiter_8[m])
    end
  end
end

local task bound_x_lo(faces_1 : region(ispace(int3d), face),
                      faces_2 : region(ispace(int3d), face),
                      faces_3 : region(ispace(int3d), face),
                      faces_4 : region(ispace(int3d), face),
                      faces_5 : region(ispace(int3d), face),
                      faces_6 : region(ispace(int3d), face),
                      faces_7 : region(ispace(int3d), face),
                      faces_8 : region(ispace(int3d), face),
                      angles : region(ispace(int1d), angle),
                      config : Config)
where
  reads (angles.{w, xi, eta, mu}),
  reads writes (faces_1.I, faces_2.I, faces_3.I, faces_4.I,
                faces_5.I, faces_6.I, faces_7.I, faces_8.I)
do

  -- Get array bounds

  var limits = faces_1.bounds

  -- Temporary variables

  var epsw = config.Radiation.xLoEmiss
  var Tw = config.Radiation.xLoTemp

  for k = limits.lo.z, limits.hi.z + 1 do
    for j = limits.lo.y, limits.hi.y + 1 do
      var value = 0.0

      -- Calculate reflect

      if epsw < 1 then
        for m = 0, NUM_ANGLES do
          if angles[m].xi < 0 then
            var face_value = 0.0
            if angles[m].eta > 0 and angles[m].mu > 0 then
              face_value = faces_5[{limits.lo.x,j,k}].I[m]
            elseif angles[m].eta > 0 and angles[m].mu <= 0 then
              face_value = faces_6[{limits.lo.x,j,k}].I[m]
            elseif angles[m].eta <= 0 and angles[m].mu > 0 then
              face_value = faces_7[{limits.lo.x,j,k}].I[m]
            else
              face_value = faces_8[{limits.lo.x,j,k}].I[m]
            end
            value += (1.0-epsw)/pi*angles[m].w*fabs(angles[m].xi)*face_value
          end
        end
      end

      -- Add blackbody radiation

      if j >= config.Radiation.xLoWindow.fromCell[0] and
         k >= config.Radiation.xLoWindow.fromCell[1] and
         j <= config.Radiation.xLoWindow.uptoCell[0] and
         k <= config.Radiation.xLoWindow.uptoCell[1] then
        value += epsw*SB*pow(Tw,4.0)/pi
      end

      -- Set Ifx values

      for m = 0, NUM_ANGLES do
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

local task bound_x_hi(faces_1 : region(ispace(int3d), face),
                      faces_2 : region(ispace(int3d), face),
                      faces_3 : region(ispace(int3d), face),
                      faces_4 : region(ispace(int3d), face),
                      faces_5 : region(ispace(int3d), face),
                      faces_6 : region(ispace(int3d), face),
                      faces_7 : region(ispace(int3d), face),
                      faces_8 : region(ispace(int3d), face),
                      angles : region(ispace(int1d), angle),
                      config : Config)
where
  reads (angles.{w, xi, eta, mu}),
  reads writes (faces_1.I, faces_2.I, faces_3.I, faces_4.I,
                faces_5.I, faces_6.I, faces_7.I, faces_8.I)
do

  -- Get array bounds

  var limits = faces_1.bounds

  -- Temporary variables

  var epsw = config.Radiation.xHiEmiss
  var Tw = config.Radiation.xHiTemp

  for k = limits.lo.z, limits.hi.z + 1 do
    for j = limits.lo.y, limits.hi.y + 1 do
      var value = 0.0

      -- Calculate reflect

      if epsw < 1 then
        for m = 0, NUM_ANGLES do
          if angles[m].xi > 0 then
            var face_value = 0.0
            if angles[m].eta > 0 and angles[m].mu > 0 then
              face_value = faces_1[{limits.hi.x,j,k}].I[m]
            elseif angles[m].eta > 0 and angles[m].mu <= 0 then
              face_value = faces_2[{limits.hi.x,j,k}].I[m]
            elseif angles[m].eta <= 0 and angles[m].mu > 0 then
              face_value = faces_3[{limits.hi.x,j,k}].I[m]
            else
              face_value = faces_4[{limits.hi.x,j,k}].I[m]
            end
            value += (1.0-epsw)/pi*angles[m].w*angles[m].xi*face_value
          end
        end
      end

      -- Add blackbody radiation

      if j >= config.Radiation.xHiWindow.fromCell[0] and
         k >= config.Radiation.xHiWindow.fromCell[1] and
         j <= config.Radiation.xHiWindow.uptoCell[0] and
         k <= config.Radiation.xHiWindow.uptoCell[1] then
        value += epsw*SB*pow(Tw,4.0)/pi
      end

      -- Set Ifx values

      for m = 0, NUM_ANGLES do
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

local task bound_y_hi(faces_1 : region(ispace(int3d), face),
                      faces_2 : region(ispace(int3d), face),
                      faces_3 : region(ispace(int3d), face),
                      faces_4 : region(ispace(int3d), face),
                      faces_5 : region(ispace(int3d), face),
                      faces_6 : region(ispace(int3d), face),
                      faces_7 : region(ispace(int3d), face),
                      faces_8 : region(ispace(int3d), face),
                      angles : region(ispace(int1d), angle),
                      config : Config)
where
  reads (angles.{w, xi, eta, mu}),
  reads writes (faces_1.I, faces_2.I, faces_3.I, faces_4.I,
                faces_5.I, faces_6.I, faces_7.I, faces_8.I)
do

  -- Get array bounds

  var limits = faces_1.bounds

  -- Temporary variables

  var epsw = config.Radiation.yHiEmiss
  var Tw = config.Radiation.yHiTemp

  for k = limits.lo.z, limits.hi.z + 1 do
    for i = limits.lo.x, limits.hi.x + 1 do
      var value = 0.0

      -- Calculate reflect

      if epsw < 1 then
        for m = 0, NUM_ANGLES do
          if angles[m].eta > 0 then
            var face_value = 0.0
            if angles[m].xi > 0 and angles[m].mu > 0 then
              face_value = faces_1[{i,limits.hi.y,k}].I[m]
            elseif angles[m].xi > 0 and angles[m].mu <= 0 then
              face_value = faces_2[{i,limits.hi.y,k}].I[m]
            elseif angles[m].xi <= 0 and angles[m].mu > 0 then
              face_value = faces_5[{i,limits.hi.y,k}].I[m]
            else
              face_value = faces_6[{i,limits.hi.y,k}].I[m]
            end
            value += (1.0-epsw)/pi*angles[m].w*angles[m].eta*face_value
          end
        end
      end

      -- Add blackbody radiation

      if i >= config.Radiation.yHiWindow.fromCell[0] and
         k >= config.Radiation.yHiWindow.fromCell[1] and
         i <= config.Radiation.yHiWindow.uptoCell[0] and
         k <= config.Radiation.yHiWindow.uptoCell[1] then
        value += epsw*SB*pow(Tw,4.0)/pi
      end

      -- Set Ify values

      for m = 0, NUM_ANGLES do
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

local task bound_y_lo(faces_1 : region(ispace(int3d), face),
                      faces_2 : region(ispace(int3d), face),
                      faces_3 : region(ispace(int3d), face),
                      faces_4 : region(ispace(int3d), face),
                      faces_5 : region(ispace(int3d), face),
                      faces_6 : region(ispace(int3d), face),
                      faces_7 : region(ispace(int3d), face),
                      faces_8 : region(ispace(int3d), face),
                      angles : region(ispace(int1d), angle),
                      config : Config)
where
  reads (angles.{w, xi, eta, mu}),
  reads writes (faces_1.I, faces_2.I, faces_3.I, faces_4.I,
                faces_5.I, faces_6.I, faces_7.I, faces_8.I)
do

  -- Get array bounds

  var limits = faces_1.bounds

  -- Temporary variables

  var epsw = config.Radiation.yLoEmiss
  var Tw = config.Radiation.yLoTemp

  for k = limits.lo.z, limits.hi.z + 1 do
    for i = limits.lo.x, limits.hi.x + 1 do
      var value = 0.0

      -- Calculate reflect

      if epsw < 1 then
        for m = 0, NUM_ANGLES do
          if angles[m].eta < 0 then
            var face_value = 0.0
            if angles[m].xi > 0 and angles[m].mu > 0 then
              face_value = faces_3[{i,limits.lo.y,k}].I[m]
            elseif angles[m].xi > 0 and angles[m].mu <= 0 then
              face_value = faces_4[{i,limits.lo.y,k}].I[m]
            elseif angles[m].xi <= 0 and angles[m].mu > 0 then
              face_value = faces_7[{i,limits.lo.y,k}].I[m]
            else
              face_value = faces_8[{i,limits.lo.y,k}].I[m]
            end
            value += (1.0-epsw)/pi*angles[m].w*fabs(angles[m].eta)*face_value
          end
        end
      end

      -- Add blackbody radiation

      if i >= config.Radiation.yLoWindow.fromCell[0] and
         k >= config.Radiation.yLoWindow.fromCell[1] and
         i <= config.Radiation.yLoWindow.uptoCell[0] and
         k <= config.Radiation.yLoWindow.uptoCell[1] then
        value += epsw*SB*pow(Tw,4.0)/pi
      end

      -- Set Ify values

      for m = 0, NUM_ANGLES do
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

local task bound_z_lo(faces_1 : region(ispace(int3d), face),
                      faces_2 : region(ispace(int3d), face),
                      faces_3 : region(ispace(int3d), face),
                      faces_4 : region(ispace(int3d), face),
                      faces_5 : region(ispace(int3d), face),
                      faces_6 : region(ispace(int3d), face),
                      faces_7 : region(ispace(int3d), face),
                      faces_8 : region(ispace(int3d), face),
                      angles : region(ispace(int1d), angle),
                      config : Config)
where
  reads (angles.{w, xi, eta, mu}),
  reads writes (faces_1.I, faces_2.I, faces_3.I, faces_4.I,
                faces_5.I, faces_6.I, faces_7.I, faces_8.I)
do

  -- Get array bounds

  var limits = faces_1.bounds

  -- Temporary variables

  var epsw = config.Radiation.zLoEmiss
  var Tw = config.Radiation.zLoTemp

  for j = limits.lo.y, limits.hi.y + 1 do
    for i = limits.lo.x, limits.hi.x + 1 do
      var value = 0.0

      -- Calculate reflect

      if epsw < 1 then
        for m = 0, NUM_ANGLES do
          if angles[m].mu < 0 then
            var face_value = 0.0
            if angles[m].xi > 0 and angles[m].eta > 0 then
              face_value = faces_2[{i,j,limits.lo.z}].I[m]
            elseif angles[m].xi > 0 and angles[m].eta <= 0 then
              face_value = faces_4[{i,j,limits.lo.z}].I[m]
            elseif angles[m].xi <= 0 and angles[m].eta > 0 then
              face_value = faces_6[{i,j,limits.lo.z}].I[m]
            else
              face_value = faces_8[{i,j,limits.lo.z}].I[m]
            end
            value += (1.0-epsw)/pi*angles[m].w*fabs(angles[m].mu)*face_value
          end
        end
      end

      -- Add blackbody radiation

      if i >= config.Radiation.zLoWindow.fromCell[0] and
         j >= config.Radiation.zLoWindow.fromCell[1] and
         i <= config.Radiation.zLoWindow.uptoCell[0] and
         j <= config.Radiation.zLoWindow.uptoCell[1] then
        value += epsw*SB*pow(Tw,4.0)/pi
      end

      -- Set Ifz values

      for m = 0, NUM_ANGLES do
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

local task bound_z_hi(faces_1 : region(ispace(int3d), face),
                      faces_2 : region(ispace(int3d), face),
                      faces_3 : region(ispace(int3d), face),
                      faces_4 : region(ispace(int3d), face),
                      faces_5 : region(ispace(int3d), face),
                      faces_6 : region(ispace(int3d), face),
                      faces_7 : region(ispace(int3d), face),
                      faces_8 : region(ispace(int3d), face),
                      angles : region(ispace(int1d), angle),
                      config : Config)
where
  reads (angles.{w, xi, eta, mu}),
  reads writes (faces_1.I, faces_2.I, faces_3.I, faces_4.I,
                faces_5.I, faces_6.I, faces_7.I, faces_8.I)
do

  -- Get array bounds

  var limits = faces_1.bounds

  -- Temporary variables

  var epsw = config.Radiation.zHiEmiss
  var Tw = config.Radiation.zHiTemp

  for j = limits.lo.y, limits.hi.y + 1 do
    for i = limits.lo.x, limits.hi.x + 1 do
      var value = 0.0

      -- Calculate reflect

      if epsw < 1 then
        for m = 0, NUM_ANGLES do
          if angles[m].mu > 0 then
            var face_value = 0.0
            if angles[m].xi > 0 and angles[m].eta > 0 then
              face_value = faces_1[{i,j,limits.hi.z}].I[m]
            elseif angles[m].xi > 0 and angles[m].eta <= 0 then
              face_value = faces_3[{i,j,limits.hi.z}].I[m]
            elseif angles[m].xi <= 0 and angles[m].eta > 0 then
              face_value = faces_5[{i,j,limits.hi.z}].I[m]
            else
              face_value = faces_7[{i,j,limits.hi.z}].I[m]
            end
            value += (1.0-epsw)/pi*angles[m].w*angles[m].mu*face_value
          end
        end
      end

      -- Add blackbody radiation

      if i >= config.Radiation.zHiWindow.fromCell[0] and
         j >= config.Radiation.zHiWindow.fromCell[1] and
         i <= config.Radiation.zHiWindow.uptoCell[0] and
         j <= config.Radiation.zHiWindow.uptoCell[1] then
        value += epsw*SB*pow(Tw,4.0)/pi
      end

      -- Set Ifz values

      for m = 0, NUM_ANGLES do
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

local task sweep_1(points : region(ispace(int3d), pointsFSpace),
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
                   xi : int64, eta : int64, mu : int64,
                   dx : double, dy : double, dz : double)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         shared_x_faces_upwind.I, shared_y_faces_upwind.I, shared_z_faces_upwind.I),
  reads writes(points.I_1, x_faces.I, y_faces.I, z_faces.I,
    shared_x_faces_downwind.I, shared_y_faces_downwind.I, shared_z_faces_downwind.I)
do

  var dAx = dy*dz
  var dAy = dx*dz
  var dAz = dx*dy
  var dV = dx*dy*dz

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

  -- xi,eta,mu can only be 1 or -1 since they are directional indicators
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


  -- Use our direction and increments for the sweep.

  for k = startz,endz,dindz do
    for j = starty,endy,dindy do
      for i = startx,endx,dindx do

      -- Loop over all angles.
      for m = 0, NUM_ANGLES do

        if (angles[m].xi * xi > 0 or (angles[m].xi == 0 and xi < 0)) and
          (angles[m].eta * eta > 0 or (angles[m].eta == 0 and eta < 0)) and
          (angles[m].mu * mu > 0 or (angles[m].mu == 0 and mu < 0)) then

            -- indx and indy are the upwind indices
            var indx : int64 = i - min(dindx,0)
            var indy : int64 = j - min(dindy,0)
            var indz : int64 = k - min(dindz,0)

            -- Determine if necessary to use ghost partition

            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x or indx > x_faces.bounds.hi.x then
              upwind_x_value = shared_x_faces_upwind[{indx,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y or indy > y_faces.bounds.hi.y then
              upwind_y_value = shared_y_faces_upwind[{i,indy,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z or indz > z_faces.bounds.hi.z then
              upwind_z_value = shared_z_faces_upwind[{i,j,indz}].I[m]
            else
              upwind_z_value = z_faces[{i,j,indz}].I[m]
            end

            -- Integrate to compute cell-centered value of I.
            points[{i,j,k}].I_1[m] = (points[{i,j,k}].S * dV
                                        + fabs(angles[m].xi) * dAx * upwind_x_value/gamma
                                        + fabs(angles[m].eta) * dAy * upwind_y_value/gamma
                                        + fabs(angles[m].mu) * dAz * upwind_z_value/gamma)
                                    /(points[{i,j,k}].sigma * dV
                                        + fabs(angles[m].xi) * dAx/gamma
                                        + fabs(angles[m].eta) * dAy/gamma
                                        + fabs(angles[m].mu) * dAz/gamma)



            -- Compute intensities on downwind faces

            var x_face_val = (points[{i,j,k}].I_1[m] - (1-gamma)*upwind_x_value)/gamma
            if (x_face_val < 0) then x_face_val = 0 end
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{indx+dindx, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_1[m] - (1-gamma)*upwind_y_value)/gamma
            if (y_face_val < 0) then y_face_val = 0 end
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, indy + dindy, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_1[m] - (1-gamma)*upwind_z_value)/gamma
            if (z_face_val < 0) then z_face_val = 0 end
            if (indz + dindz) > z_faces.bounds.hi.z or (indz + dindz) < z_faces.bounds.lo.z then
              shared_z_faces_downwind[{i, j, indz + dindz}].I[m] = z_face_val
            else
              z_faces[{i, j, indz+dindz}].I[m] = z_face_val
            end

          end
        end
      end
    end
  end
end

local task sweep_2(points : region(ispace(int3d), pointsFSpace),
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
                   xi : int64, eta : int64, mu : int64,
                   dx : double, dy : double, dz : double)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         shared_x_faces_upwind.I, shared_y_faces_upwind.I, shared_z_faces_upwind.I),
  reads writes(points.I_2, x_faces.I, y_faces.I, z_faces.I,
    shared_x_faces_downwind.I, shared_y_faces_downwind.I, shared_z_faces_downwind.I)
do
  var dAx = dy*dz
  var dAy = dx*dz
  var dAz = dx*dy
  var dV = dx*dy*dz

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

  -- Use our direction and increments for the sweep.

  for k = startz,endz,dindz do
    for j = starty,endy,dindy do
      for i = startx,endx,dindx do

      -- Loop over all angles.
      for m = 0, NUM_ANGLES do

        if (angles[m].xi * xi > 0 or (angles[m].xi == 0 and xi < 0)) and
          (angles[m].eta * eta > 0 or (angles[m].eta == 0 and eta < 0)) and
          (angles[m].mu * mu > 0 or (angles[m].mu == 0 and mu < 0)) then

            -- indx and indy are the upwind indices
            var indx : int64 = i - min(dindx,0)
            var indy : int64 = j - min(dindy,0)
            var indz : int64 = k - min(dindz,0)

            -- Determine if necessary to use ghost partition

            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x or indx > x_faces.bounds.hi.x then
              upwind_x_value = shared_x_faces_upwind[{indx,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            ---

            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y or indy > y_faces.bounds.hi.y then
              upwind_y_value = shared_y_faces_upwind[{i,indy,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z or indz > z_faces.bounds.hi.z then
              upwind_z_value = shared_z_faces_upwind[{i,j,indz}].I[m]
            else
              upwind_z_value = z_faces[{i,j,indz}].I[m]
            end

            -- Integrate to compute cell-centered value of I.

            points[{i,j,k}].I_2[m] = (points[{i,j,k}].S * dV
                                        + fabs(angles[m].xi) * dAx * upwind_x_value/gamma
                                        + fabs(angles[m].eta) * dAy * upwind_y_value/gamma
                                        + fabs(angles[m].mu) * dAz * upwind_z_value/gamma)
                                    /(points[{i,j,k}].sigma * dV
                                        + fabs(angles[m].xi) * dAx/gamma
                                        + fabs(angles[m].eta) * dAy/gamma
                                        + fabs(angles[m].mu) * dAz/gamma)

            -- Compute intensities on downwind faces

            var x_face_val = (points[{i,j,k}].I_2[m] - (1-gamma)*upwind_x_value)/gamma
            if (x_face_val < 0) then x_face_val = 0 end
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{indx + dindx, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_2[m] - (1-gamma)*upwind_y_value)/gamma
            if (y_face_val < 0) then y_face_val = 0 end
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, indy + dindy, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_2[m] - (1-gamma)*upwind_z_value)/gamma
            if (z_face_val < 0) then z_face_val = 0 end
            if (indz + dindz) > z_faces.bounds.hi.z or (indz + dindz) < z_faces.bounds.lo.z then
              shared_z_faces_downwind[{i, j, indz + dindz}].I[m] = z_face_val
            else
              z_faces[{i, j, indz+dindz}].I[m] = z_face_val
            end

          end
        end
      end
    end
  end
end

local task sweep_3(points : region(ispace(int3d), pointsFSpace),
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
                   xi : int64, eta : int64, mu : int64,
                   dx : double, dy : double, dz : double)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         shared_x_faces_upwind.I, shared_y_faces_upwind.I, shared_z_faces_upwind.I),
  reads writes(points.I_3, x_faces.I, y_faces.I, z_faces.I,
    shared_x_faces_downwind.I, shared_y_faces_downwind.I, shared_z_faces_downwind.I)
do
  var dAx = dy*dz
  var dAy = dx*dz
  var dAz = dx*dy
  var dV = dx*dy*dz

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


  -- Use our direction and increments for the sweep.

  for k = startz,endz,dindz do
    for j = starty,endy,dindy do
      for i = startx,endx,dindx do

      -- Loop over all angles.
      for m = 0, NUM_ANGLES do

        if (angles[m].xi * xi > 0 or (angles[m].xi == 0 and xi < 0)) and
          (angles[m].eta * eta > 0 or (angles[m].eta == 0 and eta < 0)) and
          (angles[m].mu * mu > 0 or (angles[m].mu == 0 and mu < 0)) then

            -- indx and indy are the upwind indices
            var indx : int64 = i - min(dindx,0)
            var indy : int64 = j - min(dindy,0)
            var indz : int64 = k - min(dindz,0)

            -- Determine if necessary to use ghost partition

            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x or indx > x_faces.bounds.hi.x then
              upwind_x_value = shared_x_faces_upwind[{indx,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            ---

            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y or indy > y_faces.bounds.hi.y then
              upwind_y_value = shared_y_faces_upwind[{i,indy,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z or indz > z_faces.bounds.hi.z then
              upwind_z_value = shared_z_faces_upwind[{i,j,indz}].I[m]
            else
              upwind_z_value = z_faces[{i,j,indz}].I[m]
            end

            -- Integrate to compute cell-centered value of I.

            points[{i,j,k}].I_3[m] = (points[{i,j,k}].S * dV
                                        + fabs(angles[m].xi) * dAx * upwind_x_value/gamma
                                        + fabs(angles[m].eta) * dAy * upwind_y_value/gamma
                                        + fabs(angles[m].mu) * dAz * upwind_z_value/gamma)
                                    /(points[{i,j,k}].sigma * dV
                                        + fabs(angles[m].xi) * dAx/gamma
                                        + fabs(angles[m].eta) * dAy/gamma
                                        + fabs(angles[m].mu) * dAz/gamma)

            -- Compute intensities on downwind faces

            var x_face_val = (points[{i,j,k}].I_3[m] - (1-gamma)*upwind_x_value)/gamma
            if (x_face_val < 0) then x_face_val = 0 end
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{indx + dindx, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_3[m] - (1-gamma)*upwind_y_value)/gamma
            if (y_face_val < 0) then y_face_val = 0 end
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, indy + dindy, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_3[m] - (1-gamma)*upwind_z_value)/gamma
            if (z_face_val < 0) then z_face_val = 0 end
            if (indz + dindz) > z_faces.bounds.hi.z or (indz + dindz) < z_faces.bounds.lo.z then
              shared_z_faces_downwind[{i, j, indz + dindz}].I[m] = z_face_val
            else
              z_faces[{i, j, indz+dindz}].I[m] = z_face_val
            end

          end
        end
      end
    end
  end
end

local task sweep_4(points : region(ispace(int3d), pointsFSpace),
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
                   xi : int64, eta : int64, mu : int64,
                   dx : double, dy : double, dz : double)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         shared_x_faces_upwind.I, shared_y_faces_upwind.I, shared_z_faces_upwind.I),
  reads writes(points.I_4, x_faces.I, y_faces.I, z_faces.I,
    shared_x_faces_downwind.I, shared_y_faces_downwind.I, shared_z_faces_downwind.I)
do
  var dAx = dy*dz
  var dAy = dx*dz
  var dAz = dx*dy
  var dV = dx*dy*dz

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


  -- Use our direction and increments for the sweep.

  for k = startz,endz,dindz do
    for j = starty,endy,dindy do
      for i = startx,endx,dindx do

      -- Loop over all angles.
      for m = 0, NUM_ANGLES do

        if (angles[m].xi * xi > 0 or (angles[m].xi == 0 and xi < 0)) and
          (angles[m].eta * eta > 0 or (angles[m].eta == 0 and eta < 0)) and
          (angles[m].mu * mu > 0 or (angles[m].mu == 0 and mu < 0)) then

            -- indx and indy are the upwind indices
            var indx : int64 = i - min(dindx,0)
            var indy : int64 = j - min(dindy,0)
            var indz : int64 = k - min(dindz,0)

            -- Determine if necessary to use ghost partition

            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x or indx > x_faces.bounds.hi.x then
              upwind_x_value = shared_x_faces_upwind[{indx,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            ---

            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y or indy > y_faces.bounds.hi.y then
              upwind_y_value = shared_y_faces_upwind[{i,indy,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z or indz > z_faces.bounds.hi.z then
              upwind_z_value = shared_z_faces_upwind[{i,j,indz}].I[m]
            else
              upwind_z_value = z_faces[{i,j,indz}].I[m]
            end

            -- Integrate to compute cell-centered value of I.

            points[{i,j,k}].I_4[m] = (points[{i,j,k}].S * dV
                                        + fabs(angles[m].xi) * dAx * upwind_x_value/gamma
                                        + fabs(angles[m].eta) * dAy * upwind_y_value/gamma
                                        + fabs(angles[m].mu) * dAz * upwind_z_value/gamma)
                                    /(points[{i,j,k}].sigma * dV
                                        + fabs(angles[m].xi) * dAx/gamma
                                        + fabs(angles[m].eta) * dAy/gamma
                                        + fabs(angles[m].mu) * dAz/gamma)

            -- Compute intensities on downwind faces

            var x_face_val = (points[{i,j,k}].I_4[m] - (1-gamma)*upwind_x_value)/gamma
            if (x_face_val < 0) then x_face_val = 0 end
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{indx + dindx, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_4[m] - (1-gamma)*upwind_y_value)/gamma
            if (y_face_val < 0) then y_face_val = 0 end
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, indy + dindy, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_4[m] - (1-gamma)*upwind_z_value)/gamma
            if (z_face_val < 0) then z_face_val = 0 end
            if (indz + dindz) > z_faces.bounds.hi.z or (indz + dindz) < z_faces.bounds.lo.z then
              shared_z_faces_downwind[{i, j, indz + dindz}].I[m] = z_face_val
            else
              z_faces[{i, j, indz+dindz}].I[m] = z_face_val
            end

          end
        end
      end
    end
  end
end

local task sweep_5(points : region(ispace(int3d), pointsFSpace),
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
                   xi : int64, eta : int64, mu : int64,
                   dx : double, dy : double, dz : double)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         shared_x_faces_upwind.I, shared_y_faces_upwind.I, shared_z_faces_upwind.I),
  reads writes(points.I_5, x_faces.I, y_faces.I, z_faces.I,
    shared_x_faces_downwind.I, shared_y_faces_downwind.I, shared_z_faces_downwind.I)
do
  var dAx = dy*dz
  var dAy = dx*dz
  var dAz = dx*dy
  var dV = dx*dy*dz

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


  -- Use our direction and increments for the sweep.

  for k = startz,endz,dindz do
    for j = starty,endy,dindy do
      for i = startx,endx,dindx do

      -- Loop over all angles.
      for m = 0, NUM_ANGLES do

        if (angles[m].xi * xi > 0 or (angles[m].xi == 0 and xi < 0)) and
          (angles[m].eta * eta > 0 or (angles[m].eta == 0 and eta < 0)) and
          (angles[m].mu * mu > 0 or (angles[m].mu == 0 and mu < 0)) then

            -- indx and indy are the upwind indices
            var indx : int64 = i - min(dindx,0)
            var indy : int64 = j - min(dindy,0)
            var indz : int64 = k - min(dindz,0)

            -- Determine if necessary to use ghost partition

            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x or indx > x_faces.bounds.hi.x then
              upwind_x_value = shared_x_faces_upwind[{indx,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            ---

            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y or indy > y_faces.bounds.hi.y then
              upwind_y_value = shared_y_faces_upwind[{i,indy,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z or indz > z_faces.bounds.hi.z then
              upwind_z_value = shared_z_faces_upwind[{i,j,indz}].I[m]
            else
              upwind_z_value = z_faces[{i,j,indz}].I[m]
            end

            -- Integrate to compute cell-centered value of I.

            points[{i,j,k}].I_5[m] = (points[{i,j,k}].S * dV
                                        + fabs(angles[m].xi) * dAx * upwind_x_value/gamma
                                        + fabs(angles[m].eta) * dAy * upwind_y_value/gamma
                                        + fabs(angles[m].mu) * dAz * upwind_z_value/gamma)
                                    /(points[{i,j,k}].sigma * dV
                                        + fabs(angles[m].xi) * dAx/gamma
                                        + fabs(angles[m].eta) * dAy/gamma
                                        + fabs(angles[m].mu) * dAz/gamma)

            -- Compute intensities on downwind faces

            var x_face_val = (points[{i,j,k}].I_5[m] - (1-gamma)*upwind_x_value)/gamma
            if (x_face_val < 0) then x_face_val = 0 end
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{indx + dindx, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_5[m] - (1-gamma)*upwind_y_value)/gamma
            if (y_face_val < 0) then y_face_val = 0 end
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, indy + dindy, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_5[m] - (1-gamma)*upwind_z_value)/gamma
            if (z_face_val < 0) then z_face_val = 0 end
            if (indz + dindz) > z_faces.bounds.hi.z or (indz + dindz) < z_faces.bounds.lo.z then
              shared_z_faces_downwind[{i, j, indz + dindz}].I[m] = z_face_val
            else
              z_faces[{i, j, indz+dindz}].I[m] = z_face_val
            end

          end
        end
      end
    end
  end
end

local task sweep_6(points : region(ispace(int3d), pointsFSpace),
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
                   xi : int64, eta : int64, mu : int64,
                   dx : double, dy : double, dz : double)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         shared_x_faces_upwind.I, shared_y_faces_upwind.I, shared_z_faces_upwind.I),
  reads writes(points.I_6, x_faces.I, y_faces.I, z_faces.I,
    shared_x_faces_downwind.I, shared_y_faces_downwind.I, shared_z_faces_downwind.I)
do
  var dAx = dy*dz
  var dAy = dx*dz
  var dAz = dx*dy
  var dV = dx*dy*dz

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


  -- Use our direction and increments for the sweep.

  for k = startz,endz,dindz do
    for j = starty,endy,dindy do
      for i = startx,endx,dindx do

      -- Loop over all angles.
      for m = 0, NUM_ANGLES do

        if (angles[m].xi * xi > 0 or (angles[m].xi == 0 and xi < 0)) and
          (angles[m].eta * eta > 0 or (angles[m].eta == 0 and eta < 0)) and
          (angles[m].mu * mu > 0 or (angles[m].mu == 0 and mu < 0)) then

            -- indx and indy are the upwind indices
            var indx : int64 = i - min(dindx,0)
            var indy : int64 = j - min(dindy,0)
            var indz : int64 = k - min(dindz,0)

            -- Determine if necessary to use ghost partition

            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x or indx > x_faces.bounds.hi.x then
              upwind_x_value = shared_x_faces_upwind[{indx,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            ---

            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y or indy > y_faces.bounds.hi.y then
              upwind_y_value = shared_y_faces_upwind[{i,indy,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z or indz > z_faces.bounds.hi.z then
              upwind_z_value = shared_z_faces_upwind[{i,j,indz}].I[m]
            else
              upwind_z_value = z_faces[{i,j,indz}].I[m]
            end

            -- Integrate to compute cell-centered value of I.

            points[{i,j,k}].I_6[m] = (points[{i,j,k}].S * dV
                                        + fabs(angles[m].xi) * dAx * upwind_x_value/gamma
                                        + fabs(angles[m].eta) * dAy * upwind_y_value/gamma
                                        + fabs(angles[m].mu) * dAz * upwind_z_value/gamma)
                                    /(points[{i,j,k}].sigma * dV
                                        + fabs(angles[m].xi) * dAx/gamma
                                        + fabs(angles[m].eta) * dAy/gamma
                                        + fabs(angles[m].mu) * dAz/gamma)

            -- Compute intensities on downwind faces

            var x_face_val = (points[{i,j,k}].I_6[m] - (1-gamma)*upwind_x_value)/gamma
            if (x_face_val < 0) then x_face_val = 0 end
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{indx + dindx, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_6[m] - (1-gamma)*upwind_y_value)/gamma
            if (y_face_val < 0) then y_face_val = 0 end
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, indy + dindy, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_6[m] - (1-gamma)*upwind_z_value)/gamma
            if (z_face_val < 0) then z_face_val = 0 end
            if (indz + dindz) > z_faces.bounds.hi.z or (indz + dindz) < z_faces.bounds.lo.z then
              shared_z_faces_downwind[{i, j, indz + dindz}].I[m] = z_face_val
            else
              z_faces[{i, j, indz+dindz}].I[m] = z_face_val
            end

          end
        end
      end
    end
  end
end

local task sweep_7(points : region(ispace(int3d), pointsFSpace),
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
                   xi : int64, eta : int64, mu : int64,
                   dx : double, dy : double, dz : double)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         shared_x_faces_upwind.I, shared_y_faces_upwind.I, shared_z_faces_upwind.I),
  reads writes(points.I_7, x_faces.I, y_faces.I, z_faces.I,
    shared_x_faces_downwind.I, shared_y_faces_downwind.I, shared_z_faces_downwind.I)
do
  var dAx = dy*dz
  var dAy = dx*dz
  var dAz = dx*dy
  var dV = dx*dy*dz

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


  -- Use our direction and increments for the sweep.

  for k = startz,endz,dindz do
    for j = starty,endy,dindy do
      for i = startx,endx,dindx do

      -- Loop over all angles.
      for m = 0, NUM_ANGLES do

        if (angles[m].xi * xi > 0 or (angles[m].xi == 0 and xi < 0)) and
          (angles[m].eta * eta > 0 or (angles[m].eta == 0 and eta < 0)) and
          (angles[m].mu * mu > 0 or (angles[m].mu == 0 and mu < 0)) then

            -- indx and indy are the upwind indices
            var indx : int64 = i - min(dindx,0)
            var indy : int64 = j - min(dindy,0)
            var indz : int64 = k - min(dindz,0)

            -- Determine if necessary to use ghost partition

            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x or indx > x_faces.bounds.hi.x then
              upwind_x_value = shared_x_faces_upwind[{indx,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            ---

            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y or indy > y_faces.bounds.hi.y then
              upwind_y_value = shared_y_faces_upwind[{i,indy,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z or indz > z_faces.bounds.hi.z then
              upwind_z_value = shared_z_faces_upwind[{i,j,indz}].I[m]
            else
              upwind_z_value = z_faces[{i,j,indz}].I[m]
            end

            -- Integrate to compute cell-centered value of I.

            points[{i,j,k}].I_7[m] = (points[{i,j,k}].S * dV
                                        + fabs(angles[m].xi) * dAx * upwind_x_value/gamma
                                        + fabs(angles[m].eta) * dAy * upwind_y_value/gamma
                                        + fabs(angles[m].mu) * dAz * upwind_z_value/gamma)
                                    /(points[{i,j,k}].sigma * dV
                                        + fabs(angles[m].xi) * dAx/gamma
                                        + fabs(angles[m].eta) * dAy/gamma
                                        + fabs(angles[m].mu) * dAz/gamma)

            -- Compute intensities on downwind faces

            var x_face_val = (points[{i,j,k}].I_7[m] - (1-gamma)*upwind_x_value)/gamma
            if (x_face_val < 0) then x_face_val = 0 end
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{indx + dindx, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_7[m] - (1-gamma)*upwind_y_value)/gamma
            if (y_face_val < 0) then y_face_val = 0 end
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, indy + dindy, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_7[m] - (1-gamma)*upwind_z_value)/gamma
            if (z_face_val < 0) then z_face_val = 0 end
            if (indz + dindz) > z_faces.bounds.hi.z or (indz + dindz) < z_faces.bounds.lo.z then
              shared_z_faces_downwind[{i, j, indz + dindz}].I[m] = z_face_val
            else
              z_faces[{i, j, indz+dindz}].I[m] = z_face_val
            end

          end
        end
      end
    end
  end
end

local task sweep_8(points : region(ispace(int3d), pointsFSpace),
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
                   xi : int64, eta : int64, mu : int64,
                   dx : double, dy : double, dz : double)
where
  reads (angles.{xi, eta, mu}, points.{S, sigma},
         shared_x_faces_upwind.I, shared_y_faces_upwind.I, shared_z_faces_upwind.I),
  reads writes(points.I_8, x_faces.I, y_faces.I, z_faces.I,
    shared_x_faces_downwind.I, shared_y_faces_downwind.I, shared_z_faces_downwind.I)
do
  var dAx = dy*dz
  var dAy = dx*dz
  var dAz = dx*dy
  var dV = dx*dy*dz

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


  -- Use our direction and increments for the sweep.

  for k = startz,endz,dindz do
    for j = starty,endy,dindy do
      for i = startx,endx,dindx do

      -- Loop over all angles.
      for m = 0, NUM_ANGLES do

        if (angles[m].xi * xi > 0 or (angles[m].xi == 0 and xi < 0)) and
          (angles[m].eta * eta > 0 or (angles[m].eta == 0 and eta < 0)) and
          (angles[m].mu * mu > 0 or (angles[m].mu == 0 and mu < 0)) then

            -- indx and indy are the upwind indices
            var indx : int64 = i - min(dindx,0)
            var indy : int64 = j - min(dindy,0)
            var indz : int64 = k - min(dindz,0)

            -- Determine if necessary to use ghost partition

            var upwind_x_value : double = 0.0
            if indx < x_faces.bounds.lo.x or indx > x_faces.bounds.hi.x then
              upwind_x_value = shared_x_faces_upwind[{indx,j,k}].I[m]
            else
              upwind_x_value = x_faces[{indx,j,k}].I[m]
            end

            ---

            var upwind_y_value : double = 0.0
            if indy < y_faces.bounds.lo.y or indy > y_faces.bounds.hi.y then
              upwind_y_value = shared_y_faces_upwind[{i,indy,k}].I[m]
            else
              upwind_y_value = y_faces[{i,indy,k}].I[m]
            end

            var upwind_z_value : double = 0.0
            if indz < z_faces.bounds.lo.z or indz > z_faces.bounds.hi.z then
              upwind_z_value = shared_z_faces_upwind[{i,j,indz}].I[m]
            else
              upwind_z_value = z_faces[{i,j,indz}].I[m]
            end

            -- Integrate to compute cell-centered value of I.

            points[{i,j,k}].I_8[m] = (points[{i,j,k}].S * dV
                                        + fabs(angles[m].xi) * dAx * upwind_x_value/gamma
                                        + fabs(angles[m].eta) * dAy * upwind_y_value/gamma
                                        + fabs(angles[m].mu) * dAz * upwind_z_value/gamma)
                                    /(points[{i,j,k}].sigma * dV
                                        + fabs(angles[m].xi) * dAx/gamma
                                        + fabs(angles[m].eta) * dAy/gamma
                                        + fabs(angles[m].mu) * dAz/gamma)

            -- Compute intensities on downwind faces

            var x_face_val = (points[{i,j,k}].I_8[m] - (1-gamma)*upwind_x_value)/gamma
            if (x_face_val < 0) then x_face_val = 0 end
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{indx + dindx, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_8[m] - (1-gamma)*upwind_y_value)/gamma
            if (y_face_val < 0) then y_face_val = 0 end
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, indy + dindy, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_8[m] - (1-gamma)*upwind_z_value)/gamma
            if (z_face_val < 0) then z_face_val = 0 end
            if (indz + dindz) > z_faces.bounds.hi.z or (indz + dindz) < z_faces.bounds.lo.z then
              shared_z_faces_downwind[{i, j, indz + dindz}].I[m] = z_face_val
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
local task residual(points : region(ispace(int3d), pointsFSpace),
                    Nx : int, Ny : int, Nz : int)
where
  reads (points.{I_1, I_2, I_3, I_4, I_5, I_6, I_7, I_8,
                 Iiter_1, Iiter_2, Iiter_3, Iiter_4,
                 Iiter_5, Iiter_6, Iiter_7, Iiter_8})
do
  var res : double = 0.0

  __demand(__openmp)
  for p in points do
    for m = 0, NUM_ANGLES do
      var v1 = p.I_1[m]
      if v1 > 0 then
        var v2 = v1 - p.Iiter_1[m]
        res += (v2 * v2) / (v1 * v1)
      end
    end
  end
  __demand(__openmp)
  for p in points do
    for m = 0, NUM_ANGLES do
      var v1 = p.I_2[m]
      if v1 > 0 then
        var v2 = v1 - p.Iiter_2[m]
        res += (v2 * v2) / (v1 * v1)
      end
    end
  end
  __demand(__openmp)
  for p in points do
    for m = 0, NUM_ANGLES do
      var v1 = p.I_3[m]
      if v1 > 0 then
        var v2 = v1 - p.Iiter_3[m]
        res += (v2 * v2) / (v1 * v1)
      end
    end
  end
  __demand(__openmp)
  for p in points do
    for m = 0, NUM_ANGLES do
      var v1 = p.I_4[m]
      if v1 > 0 then
        var v2 = v1 - p.Iiter_4[m]
        res += (v2 * v2) / (v1 * v1)
      end
    end
  end
  __demand(__openmp)
  for p in points do
    for m = 0, NUM_ANGLES do
      var v1 = p.I_5[m]
      if v1 > 0 then
        var v2 = v1 - p.Iiter_5[m]
        res += (v2 * v2) / (v1 * v1)
      end
    end
  end
  __demand(__openmp)
  for p in points do
    for m = 0, NUM_ANGLES do
      var v1 = p.I_6[m]
      if v1 > 0 then
        var v2 = v1 - p.Iiter_6[m]
        res += (v2 * v2) / (v1 * v1)
      end
    end
  end
  __demand(__openmp)
  for p in points do
    for m = 0, NUM_ANGLES do
      var v1 = p.I_7[m]
      if v1 > 0 then
        var v2 = v1 - p.Iiter_7[m]
        res += (v2 * v2) / (v1 * v1)
      end
    end
  end
  __demand(__openmp)
  for p in points do
    for m = 0, NUM_ANGLES do
      var v1 = p.I_8[m]
      if v1 > 0 then
        var v2 = v1 - p.Iiter_8[m]
        res += (v2 * v2) / (v1 * v1)
      end
    end
  end

  return res
end

-- Update the intensity before moving to the next iteration.
local task update(points : region(ispace(int3d), pointsFSpace))
where
  reads (points.{I_1, I_2, I_3, I_4, I_5, I_6, I_7, I_8}),
  reads writes (points.{Iiter_1, Iiter_2, Iiter_3, Iiter_4,
                        Iiter_5, Iiter_6, Iiter_7, Iiter_8})
do
  -- Request that outer loop is not vectorized, to ensure the inner one gets
  -- vectorized instead.
  __forbid(__vectorize) __demand(__openmp)
  for p in points do
    for m = 0, NUM_ANGLES do
      p.Iiter_1[m] = p.I_1[m]
    end
  end
  __forbid(__vectorize) __demand(__openmp)
  for p in points do
    for m = 0, NUM_ANGLES do
      p.Iiter_2[m] = p.I_2[m]
    end
  end
  __forbid(__vectorize) __demand(__openmp)
  for p in points do
    for m = 0, NUM_ANGLES do
      p.Iiter_3[m] = p.I_3[m]
    end
  end
  __forbid(__vectorize) __demand(__openmp)
  for p in points do
    for m = 0, NUM_ANGLES do
      p.Iiter_4[m] = p.I_4[m]
    end
  end
  __forbid(__vectorize) __demand(__openmp)
  for p in points do
    for m = 0, NUM_ANGLES do
      p.Iiter_5[m] = p.I_5[m]
    end
  end
  __forbid(__vectorize) __demand(__openmp)
  for p in points do
    for m = 0, NUM_ANGLES do
      p.Iiter_6[m] = p.I_6[m]
    end
  end
  __forbid(__vectorize) __demand(__openmp)
  for p in points do
    for m = 0, NUM_ANGLES do
      p.Iiter_7[m] = p.I_7[m]
    end
  end
  __forbid(__vectorize) __demand(__openmp)
  for p in points do
    for m = 0, NUM_ANGLES do
      p.Iiter_8[m] = p.I_8[m]
    end
  end
end

-- Reduce the intensity to summation over all angles
local task reduce_intensity(points : region(ispace(int3d), pointsFSpace),
                            angles : region(ispace(int1d), angle))
where
  reads (points.{I_1, I_2, I_3, I_4, I_5, I_6, I_7, I_8}, angles.w),
  reads writes (points.G)
do
  __demand(__openmp)
  for p in points do
    for m = 0, NUM_ANGLES do
      p.G += angles[m].w * (p.I_1[m])
       + angles[m].w * (p.I_2[m])
       + angles[m].w * (p.I_3[m])
       + angles[m].w * (p.I_4[m])
       + angles[m].w * (p.I_5[m])
       + angles[m].w * (p.I_6[m])
       + angles[m].w * (p.I_7[m])
       + angles[m].w * (p.I_8[m])
    end
  end
end

-------------------------------------------------------------------------------
-- FULL SIMULATION QUOTES
-------------------------------------------------------------------------------

function MODULE.mkInstance() local INSTANCE = {}

  -- Symbols shared between quotes

  local Nx = regentlib.newsymbol('Nx')
  local Ny = regentlib.newsymbol('Ny')
  local Nz = regentlib.newsymbol('Nz')

  local ntx = regentlib.newsymbol('ntx')
  local nty = regentlib.newsymbol('nty')
  local ntz = regentlib.newsymbol('ntz')

  local grid_x = regentlib.newsymbol('grid_x')
  local grid_y = regentlib.newsymbol('grid_y')
  local grid_z = regentlib.newsymbol('grid_z')

  local x_faces = {
    regentlib.newsymbol('x_faces_1'),
    regentlib.newsymbol('x_faces_2'),
    regentlib.newsymbol('x_faces_3'),
    regentlib.newsymbol('x_faces_4'),
    regentlib.newsymbol('x_faces_5'),
    regentlib.newsymbol('x_faces_6'),
    regentlib.newsymbol('x_faces_7'),
    regentlib.newsymbol('x_faces_8'),
  }

  local y_faces = {
    regentlib.newsymbol('y_faces_1'),
    regentlib.newsymbol('y_faces_2'),
    regentlib.newsymbol('y_faces_3'),
    regentlib.newsymbol('y_faces_4'),
    regentlib.newsymbol('y_faces_5'),
    regentlib.newsymbol('y_faces_6'),
    regentlib.newsymbol('y_faces_7'),
    regentlib.newsymbol('y_faces_8'),
  }

  local z_faces = {
    regentlib.newsymbol('z_faces_1'),
    regentlib.newsymbol('z_faces_2'),
    regentlib.newsymbol('z_faces_3'),
    regentlib.newsymbol('z_faces_4'),
    regentlib.newsymbol('z_faces_5'),
    regentlib.newsymbol('z_faces_6'),
    regentlib.newsymbol('z_faces_7'),
    regentlib.newsymbol('z_faces_8'),
  }

  local angles = regentlib.newsymbol('angles')
  local tiles_private = regentlib.newsymbol('tiles_private')
  local x_tiles_shared = regentlib.newsymbol('x_tiles_shared')
  local y_tiles_shared = regentlib.newsymbol('y_tiles_shared')
  local z_tiles_shared = regentlib.newsymbol('z_tiles_shared')

  local s_x_faces = {
    regentlib.newsymbol('s_x_faces_1'),
    regentlib.newsymbol('s_x_faces_2'),
    regentlib.newsymbol('s_x_faces_3'),
    regentlib.newsymbol('s_x_faces_4'),
    regentlib.newsymbol('s_x_faces_5'),
    regentlib.newsymbol('s_x_faces_6'),
    regentlib.newsymbol('s_x_faces_7'),
    regentlib.newsymbol('s_x_faces_8'),
  }
  local s_y_faces = {
    regentlib.newsymbol('s_y_faces_1'),
    regentlib.newsymbol('s_y_faces_2'),
    regentlib.newsymbol('s_y_faces_3'),
    regentlib.newsymbol('s_y_faces_4'),
    regentlib.newsymbol('s_y_faces_5'),
    regentlib.newsymbol('s_y_faces_6'),
    regentlib.newsymbol('s_y_faces_7'),
    regentlib.newsymbol('s_y_faces_8'),
  }
  local s_z_faces = {
    regentlib.newsymbol('s_z_faces_1'),
    regentlib.newsymbol('s_z_faces_2'),
    regentlib.newsymbol('s_z_faces_3'),
    regentlib.newsymbol('s_z_faces_4'),
    regentlib.newsymbol('s_z_faces_5'),
    regentlib.newsymbol('s_z_faces_6'),
    regentlib.newsymbol('s_z_faces_7'),
    regentlib.newsymbol('s_z_faces_8'),
  }

  local p_x_faces_equal = {
    regentlib.newsymbol('p_x_faces_equal_1'),
    regentlib.newsymbol('p_x_faces_equal_2'),
    regentlib.newsymbol('p_x_faces_equal_3'),
    regentlib.newsymbol('p_x_faces_equal_4'),
    regentlib.newsymbol('p_x_faces_equal_5'),
    regentlib.newsymbol('p_x_faces_equal_6'),
    regentlib.newsymbol('p_x_faces_equal_7'),
    regentlib.newsymbol('p_x_faces_equal_8'),
  }
  local p_y_faces_equal = {
    regentlib.newsymbol('p_y_faces_equal_1'),
    regentlib.newsymbol('p_y_faces_equal_2'),
    regentlib.newsymbol('p_y_faces_equal_3'),
    regentlib.newsymbol('p_y_faces_equal_4'),
    regentlib.newsymbol('p_y_faces_equal_5'),
    regentlib.newsymbol('p_y_faces_equal_6'),
    regentlib.newsymbol('p_y_faces_equal_7'),
    regentlib.newsymbol('p_y_faces_equal_8'),
  }
  local p_z_faces_equal = {
    regentlib.newsymbol('p_z_faces_equal_1'),
    regentlib.newsymbol('p_z_faces_equal_2'),
    regentlib.newsymbol('p_z_faces_equal_3'),
    regentlib.newsymbol('p_z_faces_equal_4'),
    regentlib.newsymbol('p_z_faces_equal_5'),
    regentlib.newsymbol('p_z_faces_equal_6'),
    regentlib.newsymbol('p_z_faces_equal_7'),
    regentlib.newsymbol('p_z_faces_equal_8'),
  }

  local p_x_faces = {
    regentlib.newsymbol('p_x_faces_1'),
    regentlib.newsymbol('p_x_faces_2'),
    regentlib.newsymbol('p_x_faces_3'),
    regentlib.newsymbol('p_x_faces_4'),
    regentlib.newsymbol('p_x_faces_5'),
    regentlib.newsymbol('p_x_faces_6'),
    regentlib.newsymbol('p_x_faces_7'),
    regentlib.newsymbol('p_x_faces_8'),
  }
  local p_y_faces = {
    regentlib.newsymbol('p_y_faces_1'),
    regentlib.newsymbol('p_y_faces_2'),
    regentlib.newsymbol('p_y_faces_3'),
    regentlib.newsymbol('p_y_faces_4'),
    regentlib.newsymbol('p_y_faces_5'),
    regentlib.newsymbol('p_y_faces_6'),
    regentlib.newsymbol('p_y_faces_7'),
    regentlib.newsymbol('p_y_faces_8'),
  }
  local p_z_faces = {
    regentlib.newsymbol('p_z_faces_1'),
    regentlib.newsymbol('p_z_faces_2'),
    regentlib.newsymbol('p_z_faces_3'),
    regentlib.newsymbol('p_z_faces_4'),
    regentlib.newsymbol('p_z_faces_5'),
    regentlib.newsymbol('p_z_faces_6'),
    regentlib.newsymbol('p_z_faces_7'),
    regentlib.newsymbol('p_z_faces_8'),
  }

  function INSTANCE.DeclSymbols(config)

    local decl_symbols = rquote
      var sampleId = config.Mapping.sampleId

      -- Number of points in each dimension
      var [Nx] = config.Radiation.xNum
      var [Ny] = config.Radiation.yNum
      var [Nz] = config.Radiation.zNum

      -- Number of tiles in each dimension
      var [ntx] = config.Mapping.tiles[0]
      var [nty] = config.Mapping.tiles[1]
      var [ntz] = config.Mapping.tiles[2]

      -- Regions for faces (+1 in one direction since one more face than points)
      var [grid_x] = ispace(int3d, {x = Nx+1, y = Ny,   z = Nz  })
      var [grid_y] = ispace(int3d, {x = Nx,   y = Ny+1, z = Nz  })
      var [grid_z] = ispace(int3d, {x = Nx,   y = Ny,   z = Nz+1})

      -- 1D Region for angle values
      var angle_indices = ispace(int1d, NUM_ANGLES)
      var [angles] = region(angle_indices, angle);
      [UTIL.mkRegionTagAttach(angles, MAPPER.SAMPLE_ID_TAG, sampleId, int)];
    end

    for i = 1, 8 do
      decl_symbols = rquote
        [decl_symbols];
        var sampleId = config.Mapping.sampleId
        var [x_faces[i]] = region([grid_x], face);
        [UTIL.mkRegionTagAttach(x_faces[i], MAPPER.SAMPLE_ID_TAG, sampleId, int)];
        var [y_faces[i]] = region([grid_y], face);
        [UTIL.mkRegionTagAttach(y_faces[i], MAPPER.SAMPLE_ID_TAG, sampleId, int)];
        var [z_faces[i]] = region([grid_z], face);
        [UTIL.mkRegionTagAttach(z_faces[i], MAPPER.SAMPLE_ID_TAG, sampleId, int)];
      end
    end

    decl_symbols = rquote
      [decl_symbols];
      -- extra tile required for shared edge
      var [tiles_private]  = ispace(int3d, {x = ntx,   y = nty,   z = ntz  })
      var [x_tiles_shared] = ispace(int3d, {x = ntx+1, y = nty,   z = ntz  })
      var [y_tiles_shared] = ispace(int3d, {x = ntx,   y = nty+1, z = ntz  })
      var [z_tiles_shared] = ispace(int3d, {x = ntx,   y = nty,   z = ntz+1})
    end

    -- Partition faces

    local directions = {
      rexpr array(true,true,true) end,
      rexpr array(true,true,false) end,
      rexpr array(true,false,true) end,
      rexpr array(true,false,false) end,
      rexpr array(false,true,true) end,
      rexpr array(false,true,false) end,
      rexpr array(false,false,true) end,
      rexpr array(false,false,false) end,
    }

    for i = 1, 8 do

      local x_by_privacy = regentlib.newsymbol('x_by_privacy_' .. tostring(i))
      local p_x = regentlib.newsymbol('p_x_' .. tostring(i))
      local s_x = regentlib.newsymbol('s_x_' .. tostring(i))

      local y_by_privacy = regentlib.newsymbol('y_by_privacy_' .. tostring(i))
      local p_y = regentlib.newsymbol('p_y_' .. tostring(i))
      local s_y = regentlib.newsymbol('s_y_' .. tostring(i))

      local z_by_privacy = regentlib.newsymbol('z_by_privacy_' .. tostring(i))
      local p_z = regentlib.newsymbol('p_z_' .. tostring(i))
      local s_z = regentlib.newsymbol('s_z_' .. tostring(i))

      decl_symbols = rquote
        [decl_symbols];
        -- x
        var [p_x_faces_equal[i]] = partition(equal, [x_faces[i]], [tiles_private])
        for c in [tiles_private] do
          color_faces([p_x_faces_equal[i]][c], Nx, Ny, Nz, ntx, nty, ntz, 0, [directions[i]])
        end
        var [x_by_privacy] = partition([x_faces[i]].is_private, ispace(int1d,2))
        var [p_x] = [x_by_privacy][1]
        var [p_x_faces[i]] = partition([p_x].tile, [tiles_private])
        var [s_x] = [x_by_privacy][0]
        var [s_x_faces[i]] = partition([s_x].tile, [x_tiles_shared])
        -- y
        var [p_y_faces_equal[i]] = partition(equal, [y_faces[i]], [tiles_private])
        for c in [tiles_private] do
          color_faces([p_y_faces_equal[i]][c], Nx, Ny, Nz, ntx, nty, ntz, 1, [directions[i]])
        end
        var [y_by_privacy] = partition([y_faces[i]].is_private, ispace(int1d,2))
        var [p_y] = [y_by_privacy][1]
        var [p_y_faces[i]] = partition([p_y].tile, [tiles_private])
        var [s_y] = [y_by_privacy][0]
        var [s_y_faces[i]] = partition([s_y].tile, [y_tiles_shared])
        -- z
        var [p_z_faces_equal[i]] = partition(equal, [z_faces[i]], [tiles_private])
        for c in [tiles_private] do
          color_faces([p_z_faces_equal[i]][c], Nx, Ny, Nz, ntx, nty, ntz, 2, [directions[i]])
        end
        var [z_by_privacy] = partition([z_faces[i]].is_private, ispace(int1d,2))
        var [p_z] = [z_by_privacy][1]
        var [p_z_faces[i]] = partition([p_z].tile, [tiles_private])
        var [s_z] = [z_by_privacy][0]
        var [s_z_faces[i]] = partition([s_z].tile, [z_tiles_shared])
      end
    end

    return decl_symbols

  end -- DeclSymbols

  function INSTANCE.InitRegions() return rquote

    -- Initialize angle values
    initialize_angles(angles)

  end end -- InitRegions

  function INSTANCE.ComputeRadiationField(config, tiles, p_points) return rquote

    var dx = config.Grid.xWidth / config.Radiation.xNum
    var dy = config.Grid.yWidth / config.Radiation.yNum
    var dz = config.Grid.zWidth / config.Radiation.zNum

    var omega = config.Radiation.qs/(config.Radiation.qa+config.Radiation.qs)

    -- Compute until convergence
    var res : double = 1.0
    while (res > tol) do

      -- Update the source term (in this problem, isotropic)
      for t in tiles do
        source_term(p_points[t], angles, omega)
      end

      -- Update the grid boundary intensities
      -- TODO: Should launch these on just the boundaries
      for j = 0, nty do
        for k = 0, ntz do
          bound_x_lo([s_x_faces[1]][{0,j,k}],
                     [s_x_faces[2]][{0,j,k}],
                     [s_x_faces[3]][{0,j,k}],
                     [s_x_faces[4]][{0,j,k}],
                     [s_x_faces[5]][{0,j,k}],
                     [s_x_faces[6]][{0,j,k}],
                     [s_x_faces[7]][{0,j,k}],
                     [s_x_faces[8]][{0,j,k}],
                     angles,
                     config)
        end
      end

      for j = 0, nty do
        for k = 0, ntz do
          bound_x_hi([s_x_faces[1]][{ntx,j,k}],
                     [s_x_faces[2]][{ntx,j,k}],
                     [s_x_faces[3]][{ntx,j,k}],
                     [s_x_faces[4]][{ntx,j,k}],
                     [s_x_faces[5]][{ntx,j,k}],
                     [s_x_faces[6]][{ntx,j,k}],
                     [s_x_faces[7]][{ntx,j,k}],
                     [s_x_faces[8]][{ntx,j,k}],
                     angles,
                     config)
        end
      end

      -- Update y faces
      for i = 0, ntx do
        for k = 0, ntz do
          bound_y_lo([s_y_faces[1]][{i,0,k}],
                     [s_y_faces[2]][{i,0,k}],
                     [s_y_faces[3]][{i,0,k}],
                     [s_y_faces[4]][{i,0,k}],
                     [s_y_faces[5]][{i,0,k}],
                     [s_y_faces[6]][{i,0,k}],
                     [s_y_faces[7]][{i,0,k}],
                     [s_y_faces[8]][{i,0,k}],
                     angles,
                     config)
        end
      end

      for i = 0, ntx do
        for k = 0, ntz do
          bound_y_hi([s_y_faces[1]][{i,nty,k}],
                     [s_y_faces[2]][{i,nty,k}],
                     [s_y_faces[3]][{i,nty,k}],
                     [s_y_faces[4]][{i,nty,k}],
                     [s_y_faces[5]][{i,nty,k}],
                     [s_y_faces[6]][{i,nty,k}],
                     [s_y_faces[7]][{i,nty,k}],
                     [s_y_faces[8]][{i,nty,k}],
                     angles,
                     config)
        end
      end

      -- Update z faces
      for i = 0, ntx do
        for j = 0, nty do
          bound_z_lo([s_z_faces[1]][{i,j,0}],
                     [s_z_faces[2]][{i,j,0}],
                     [s_z_faces[3]][{i,j,0}],
                     [s_z_faces[4]][{i,j,0}],
                     [s_z_faces[5]][{i,j,0}],
                     [s_z_faces[6]][{i,j,0}],
                     [s_z_faces[7]][{i,j,0}],
                     [s_z_faces[8]][{i,j,0}],
                     angles,
                     config)
        end
      end

      for i = 0, ntx do
        for j = 0, nty do
          bound_z_hi([s_z_faces[1]][{i,j,ntz}],
                     [s_z_faces[2]][{i,j,ntz}],
                     [s_z_faces[3]][{i,j,ntz}],
                     [s_z_faces[4]][{i,j,ntz}],
                     [s_z_faces[5]][{i,j,ntz}],
                     [s_z_faces[6]][{i,j,ntz}],
                     [s_z_faces[7]][{i,j,ntz}],
                     [s_z_faces[8]][{i,j,ntz}],
                     angles,
                     config)
        end
      end

      --Perform the sweep for computing new intensities
      --Quadrant 1 - +x, +y, +z
      for i = 0, ntx do
        for j = 0, nty do
          for k = 0, ntz do
            sweep_1(p_points[{i,j,k}],
                    [p_x_faces[1]][{i,j,k}], [p_y_faces[1]][{i,j,k}], [p_z_faces[1]][{i,j,k}],
                    [s_x_faces[1]][{i,j,k}], [s_x_faces[1]][{i+1,j,k}],
                    [s_y_faces[1]][{i,j,k}], [s_y_faces[1]][{i,j+1,k}],
                    [s_z_faces[1]][{i,j,k}], [s_z_faces[1]][{i,j,k+1}],
                    angles, 1, 1, 1, dx, dy, dz)
          end
        end
      end

      -- Quadrant 2 - +x, +y, -z
      for i = 0, ntx do
        for j = 0, nty do
          for k = ntz-1, -1, -1 do
            sweep_2(p_points[{i,j,k}],
                    [p_x_faces[2]][{i,j,k}], [p_y_faces[2]][{i,j,k}], [p_z_faces[2]][{i,j,k}],
                    [s_x_faces[2]][{i,j,k}], [s_x_faces[2]][{i+1,j,k}],
                    [s_y_faces[2]][{i,j,k}], [s_y_faces[2]][{i,j+1,k}],
                    [s_z_faces[2]][{i,j,k+1}], [s_z_faces[2]][{i,j,k}],
                    angles, 1, 1, -1, dx, dy, dz)
          end
        end
      end

      -- Quadrant 3 - +x, -y, +z
      for i = 0, ntx do
        for j = nty-1, -1, -1 do
          for k = 0, ntz do
            sweep_3(p_points[{i,j,k}],
                    [p_x_faces[3]][{i,j,k}], [p_y_faces[3]][{i,j,k}], [p_z_faces[3]][{i,j,k}],
                    [s_x_faces[3]][{i,j,k}], [s_x_faces[3]][{i+1,j,k}],
                    [s_y_faces[3]][{i,j+1,k}], [s_y_faces[3]][{i,j,k}],
                    [s_z_faces[3]][{i,j,k}], [s_z_faces[3]][{i,j,k+1}],
                    angles, 1, -1, 1, dx, dy, dz)
          end
        end
      end

      -- Quadrant 4 - +x, -y, -z
      for i = 0, ntx do
        for j = nty-1, -1, -1 do
          for k = ntz-1, -1, -1 do
            sweep_4(p_points[{i,j,k}],
                    [p_x_faces[4]][{i,j,k}], [p_y_faces[4]][{i,j,k}], [p_z_faces[4]][{i,j,k}],
                    [s_x_faces[4]][{i,j,k}], [s_x_faces[4]][{i+1,j,k}],
                    [s_y_faces[4]][{i,j+1,k}], [s_y_faces[4]][{i,j,k}],
                    [s_z_faces[4]][{i,j,k+1}], [s_z_faces[4]][{i,j,k}],
                    angles, 1, -1, -1, dx, dy, dz)
          end
        end
      end

      -- Quadrant 5 - -x, +y, +z
      for i = ntx-1, -1, -1 do
        for j = 0, nty do
          for k = 0, ntz do
            sweep_5(p_points[{i,j,k}],
                    [p_x_faces[5]][{i,j,k}], [p_y_faces[5]][{i,j,k}], [p_z_faces[5]][{i,j,k}],
                    [s_x_faces[5]][{i+1,j,k}], [s_x_faces[5]][{i,j,k}],
                    [s_y_faces[5]][{i,j,k}], [s_y_faces[5]][{i,j+1,k}],
                    [s_z_faces[5]][{i,j,k}], [s_z_faces[5]][{i,j,k+1}],
                    angles, -1, 1, 1, dx, dy, dz)
          end
        end
      end

      -- Quadrant 6 - -x, +y, -z
      for i = ntx-1, -1, -1 do
        for j = 0, nty do
          for k = ntz-1, -1, -1 do
            sweep_6(p_points[{i,j,k}],
                    [p_x_faces[6]][{i,j,k}], [p_y_faces[6]][{i,j,k}], [p_z_faces[6]][{i,j,k}],
                    [s_x_faces[6]][{i+1,j,k}], [s_x_faces[6]][{i,j,k}],
                    [s_y_faces[6]][{i,j,k}], [s_y_faces[6]][{i,j+1,k}],
                    [s_z_faces[6]][{i,j,k+1}], [s_z_faces[6]][{i,j,k}],
                    angles, -1, 1, -1, dx, dy, dz)
          end
        end
      end

      -- Quadrant 7 - -x, -y, +z
      for i = ntx-1, -1, -1 do
        for j = nty-1, -1, -1 do
          for k = 0, ntz do
            sweep_7(p_points[{i,j,k}],
                    [p_x_faces[7]][{i,j,k}], [p_y_faces[7]][{i,j,k}], [p_z_faces[7]][{i,j,k}],
                    [s_x_faces[7]][{i+1,j,k}], [s_x_faces[7]][{i,j,k}],
                    [s_y_faces[7]][{i,j+1,k}], [s_y_faces[7]][{i,j,k}],
                    [s_z_faces[7]][{i,j,k}], [s_z_faces[7]][{i,j,k+1}],
                    angles, -1, -1, 1, dx, dy, dz)
          end
        end
      end

      -- Quadrant 8 - -x, -y, -z
      for i = ntx-1, -1, -1 do
        for j = nty-1, -1, -1 do
          for k = ntz-1, -1, -1 do
            sweep_8(p_points[{i,j,k}],
                    [p_x_faces[8]][{i,j,k}], [p_y_faces[8]][{i,j,k}], [p_z_faces[8]][{i,j,k}],
                    [s_x_faces[8]][{i+1,j,k}], [s_x_faces[8]][{i,j,k}],
                    [s_y_faces[8]][{i,j+1,k}], [s_y_faces[8]][{i,j,k}],
                    [s_z_faces[8]][{i,j,k+1}], [s_z_faces[8]][{i,j,k}],
                    angles, -1, -1, -1, dx, dy, dz)
          end
        end
      end

      -- Compute the residual
      res = 0.0
      for t in tiles do
        res += residual(p_points[t], Nx, Ny, Nz)
      end
      res = sqrt(res/(Nx*Ny*Nz*(NUM_ANGLES)))

      -- Update the intensities
      for t in tiles do
        update(p_points[t])
      end

    end

    -- Reduce intensity
    for t in tiles do
      reduce_intensity(p_points[t], angles)
    end

  end end -- ComputeRadiationField

return INSTANCE end -- mkInstance

-------------------------------------------------------------------------------
-- MODULE END
-------------------------------------------------------------------------------

return MODULE end
