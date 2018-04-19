import 'regent'

-------------------------------------------------------------------------------
-- MODULE PARAMETERS
-------------------------------------------------------------------------------

return function(NUM_ANGLES, pointsFSpace)

-------------------------------------------------------------------------------
-- COMPILE-TIME COMPUTATION
-------------------------------------------------------------------------------

-- C imports

local C = regentlib.c

-- Math imports

local fabs = regentlib.fabs(double)
local max = regentlib.fmax
local min = regentlib.fmin
local pow = regentlib.pow(double)
local sqrt = regentlib.sqrt(double)

-- Math Constants

local pi = 3.1415926535898

-- Quadrature file name

local quad_file = 'LMquads/'..NUM_ANGLES..'.txt'

-- Wall temperatures

local SB = 5.67e-8

-- Procedure parameters

local tol   = 1e-6   -- solution tolerance
local gamma = 0.5    -- 1 for step differencing, 0.5 for diamond differencing

local terra open_quad_file() : &C.FILE
  var f = C.fopen([quad_file], 'rb')
  if f == nil then
    C.printf('Error opening angle file\n')
    C.exit(1)
  end
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

-------------------------------------------------------------------------------
-- MODULE-LOCAL TASKS
-------------------------------------------------------------------------------

-- Initialize face values
local task initialize_faces(faces : region(ispace(int3d), face))
where
  reads writes (faces.I)
do
  for f in faces do
    for m = 0, NUM_ANGLES do
      f.I[m] = 0.0
    end
  end
end

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

-- Nx = 8 (num x) but contains 9 since Nx + 1
-- ntx = 2 (tiles)
-- x_tiles contains ntx+1 in x direction for extra shared (3 tiles)
-- s -- p -- p -- p -- s -- p -- p -- p -- s
-- 0                   4                   8

-- Nx = 6
-- ntx = 2 (tiles)
-- s -- p -- p -- s -- p -- p -- s
-- 0              3              6

-- shared = i % (Nx/ntx) == 0
-- private

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

local task bound_x_lo(faces_1 : region(ispace(int3d), face),
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
          elseif angles[m].eta > 0 and angles[m].mu <= 0 then
            face_value = faces_6[{limits.lo.x,j,k}].I[m]
          elseif angles[m].eta <= 0 and angles[m].mu > 0 then
            face_value = faces_7[{limits.lo.x,j,k}].I[m]
          else
            face_value = faces_8[{limits.lo.x,j,k}].I[m]
          end
          reflect += (1.0-epsw)/pi*angles[m].w*fabs(angles[m].xi)*face_value
        end
      end

      -- Set Ifx values using reflect

      var value : double = epsw*SB*pow(Tw,4.0)/pi + reflect
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
          elseif angles[m].xi > 0 and angles[m].mu <= 0 then
            face_value = faces_4[{i,limits.lo.y,k}].I[m]
          elseif angles[m].xi <= 0 and angles[m].mu > 0 then
            face_value = faces_7[{i,limits.lo.y,k}].I[m]
          else
            face_value = faces_8[{i,limits.lo.y,k}].I[m]
          end
          reflect += (1.0-epsw)/pi*angles[m].w*fabs(angles[m].eta)*face_value
        end
      end

      -- Set Ify values using reflect

      var value : double = epsw*SB*pow(Tw,4.0)/pi + reflect
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
          reflect += (1.0-epsw)/pi*angles[m].w*fabs(angles[m].mu)*face_value
        end
      end

      -- Set Ifz values using reflect

      var value : double = epsw*SB*pow(Tw,4.0)/pi + reflect
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
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{indx+dindx, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_1[m] - (1-gamma)*upwind_y_value)/gamma
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, indy + dindy, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_1[m] - (1-gamma)*upwind_z_value)/gamma
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
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{indx + dindx, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_2[m] - (1-gamma)*upwind_y_value)/gamma
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, indy + dindy, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_2[m] - (1-gamma)*upwind_z_value)/gamma
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
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{indx + dindx, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_3[m] - (1-gamma)*upwind_y_value)/gamma
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, indy + dindy, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_3[m] - (1-gamma)*upwind_z_value)/gamma
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
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{indx + dindx, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_4[m] - (1-gamma)*upwind_y_value)/gamma
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, indy + dindy, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_4[m] - (1-gamma)*upwind_z_value)/gamma
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
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{indx + dindx, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_5[m] - (1-gamma)*upwind_y_value)/gamma
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, indy + dindy, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_5[m] - (1-gamma)*upwind_z_value)/gamma
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
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{indx + dindx, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_6[m] - (1-gamma)*upwind_y_value)/gamma
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, indy + dindy, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_6[m] - (1-gamma)*upwind_z_value)/gamma
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
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{indx + dindx, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_7[m] - (1-gamma)*upwind_y_value)/gamma
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, indy + dindy, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_7[m] - (1-gamma)*upwind_z_value)/gamma
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
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{indx + dindx, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end

            var y_face_val = (points[{i,j,k}].I_8[m] - (1-gamma)*upwind_y_value)/gamma
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, indy + dindy, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end

            var z_face_val = (points[{i,j,k}].I_8[m] - (1-gamma)*upwind_z_value)/gamma
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

-- Update the intensity before moving to the next iteration.
local task update(points : region(ispace(int3d), pointsFSpace))
where
  reads (points.{I_1, I_2, I_3, I_4, I_5, I_6, I_7, I_8}),
  reads writes (points.{Iiter_1, Iiter_2, Iiter_3, Iiter_4,
                        Iiter_5, Iiter_6, Iiter_7, Iiter_8})
do
  __demand(__openmp)
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

local Exports = {}

-- Symbols shared between quotes

local Nx = regentlib.newsymbol('Nx')
local Ny = regentlib.newsymbol('Ny')
local Nz = regentlib.newsymbol('Nz')

local ntx = regentlib.newsymbol('ntx')
local nty = regentlib.newsymbol('nty')
local ntz = regentlib.newsymbol('ntz')

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
local diagonals_private = regentlib.newsymbol('diagonals_private')
local diagonals_shared = regentlib.newsymbol('diagonals_shared')

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

local s_x_faces_by_tile = {
  regentlib.newsymbol('s_x_faces_by_tile_1'),
  regentlib.newsymbol('s_x_faces_by_tile_2'),
  regentlib.newsymbol('s_x_faces_by_tile_3'),
  regentlib.newsymbol('s_x_faces_by_tile_4'),
  regentlib.newsymbol('s_x_faces_by_tile_5'),
  regentlib.newsymbol('s_x_faces_by_tile_6'),
  regentlib.newsymbol('s_x_faces_by_tile_7'),
  regentlib.newsymbol('s_x_faces_by_tile_8'),
}
local s_y_faces_by_tile = {
  regentlib.newsymbol('s_y_faces_by_tile_1'),
  regentlib.newsymbol('s_y_faces_by_tile_2'),
  regentlib.newsymbol('s_y_faces_by_tile_3'),
  regentlib.newsymbol('s_y_faces_by_tile_4'),
  regentlib.newsymbol('s_y_faces_by_tile_5'),
  regentlib.newsymbol('s_y_faces_by_tile_6'),
  regentlib.newsymbol('s_y_faces_by_tile_7'),
  regentlib.newsymbol('s_y_faces_by_tile_8'),
}
local s_z_faces_by_tile = {
  regentlib.newsymbol('s_z_faces_by_tile_1'),
  regentlib.newsymbol('s_z_faces_by_tile_2'),
  regentlib.newsymbol('s_z_faces_by_tile_3'),
  regentlib.newsymbol('s_z_faces_by_tile_4'),
  regentlib.newsymbol('s_z_faces_by_tile_5'),
  regentlib.newsymbol('s_z_faces_by_tile_6'),
  regentlib.newsymbol('s_z_faces_by_tile_7'),
  regentlib.newsymbol('s_z_faces_by_tile_8'),
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

function Exports.DeclSymbols(config) return rquote

  -- Number of points in each dimension
  var [Nx] = config.Radiation.xNum
  var [Ny] = config.Radiation.yNum
  var [Nz] = config.Radiation.zNum

  -- Number of tiles in each dimension
  var [ntx] = config.Mapping.xTiles
  var [nty] = config.Mapping.yTiles
  var [ntz] = config.Mapping.zTiles

  -- Regions for faces (+1 in one direction since one more face than points)
  var grid_x = ispace(int3d, {x = Nx+1, y = Ny,   z = Nz  })
  var grid_y = ispace(int3d, {x = Nx,   y = Ny+1, z = Nz  })
  var grid_z = ispace(int3d, {x = Nx,   y = Ny,   z = Nz+1})

  var [x_faces[1]] = region(grid_x, face)
  var [x_faces[2]] = region(grid_x, face)
  var [x_faces[3]] = region(grid_x, face)
  var [x_faces[4]] = region(grid_x, face)
  var [x_faces[5]] = region(grid_x, face)
  var [x_faces[6]] = region(grid_x, face)
  var [x_faces[7]] = region(grid_x, face)
  var [x_faces[8]] = region(grid_x, face)

  var [y_faces[1]] = region(grid_y, face)
  var [y_faces[2]] = region(grid_y, face)
  var [y_faces[3]] = region(grid_y, face)
  var [y_faces[4]] = region(grid_y, face)
  var [y_faces[5]] = region(grid_y, face)
  var [y_faces[6]] = region(grid_y, face)
  var [y_faces[7]] = region(grid_y, face)
  var [y_faces[8]] = region(grid_y, face)

  var [z_faces[1]] = region(grid_z, face)
  var [z_faces[2]] = region(grid_z, face)
  var [z_faces[3]] = region(grid_z, face)
  var [z_faces[4]] = region(grid_z, face)
  var [z_faces[5]] = region(grid_z, face)
  var [z_faces[6]] = region(grid_z, face)
  var [z_faces[7]] = region(grid_z, face)
  var [z_faces[8]] = region(grid_z, face)

  -- 1D Region for angle values
  var angle_indices = ispace(int1d, NUM_ANGLES)
  var [angles] = region(angle_indices, angle)

  -- Partition faces

  -- extra tile required for shared edge
  var [tiles_private]  = ispace(int3d, {x = ntx,   y = nty,   z = ntz  })
  var [x_tiles_shared] = ispace(int3d, {x = ntx+1, y = nty,   z = ntz  })
  var [y_tiles_shared] = ispace(int3d, {x = ntx,   y = nty+1, z = ntz  })
  var [z_tiles_shared] = ispace(int3d, {x = ntx,   y = nty,   z = ntz+1})
  var [diagonals_private] = ispace(int1d, ntx-1 + nty-1 + ntz-1 + 1)
  var [diagonals_shared]  = ispace(int1d, ntx-1 + nty-1 + ntz-1 + 2)

  -- x

  color_faces([x_faces[1]], Nx, Ny, Nz, ntx, nty, ntz, 0, array(true,true,true))
  var x_1_by_privacy = partition([x_faces[1]].is_private, ispace(int1d,2))
  var p_x_1 = x_1_by_privacy[1]
  var p_x_1_by_diagonal = partition(p_x_1.diagonal, diagonals_private)
  var p_x_1_by_tile = partition(p_x_1.tile, tiles_private)
  var [p_x_faces[1]] = cross_product(p_x_1_by_diagonal, p_x_1_by_tile)
  var s_x_1 = x_1_by_privacy[0]
  var s_x_1_by_diagonal = partition(s_x_1.diagonal, diagonals_shared)
  var [s_x_faces_by_tile[1]] = partition(s_x_1.tile, x_tiles_shared)
  var [s_x_faces[1]] = cross_product(s_x_1_by_diagonal, [s_x_faces_by_tile[1]])

  color_faces([x_faces[2]], Nx, Ny, Nz, ntx, nty, ntz, 0, array(true,true,false))
  var x_2_by_privacy = partition([x_faces[2]].is_private, ispace(int1d,2))
  var p_x_2 = x_2_by_privacy[1]
  var p_x_2_by_diagonal = partition(p_x_2.diagonal, diagonals_private)
  var p_x_2_by_tile = partition(p_x_2.tile, tiles_private)
  var [p_x_faces[2]] = cross_product(p_x_2_by_diagonal, p_x_2_by_tile)
  var s_x_2 = x_2_by_privacy[0]
  var s_x_2_by_diagonal = partition(s_x_2.diagonal, diagonals_shared)
  var [s_x_faces_by_tile[2]] = partition(s_x_2.tile, x_tiles_shared)
  var [s_x_faces[2]] = cross_product(s_x_2_by_diagonal, [s_x_faces_by_tile[2]])

  color_faces([x_faces[3]], Nx, Ny, Nz, ntx, nty, ntz, 0, array(true,false,true))
  var x_3_by_privacy = partition([x_faces[3]].is_private, ispace(int1d,2))
  var p_x_3 = x_3_by_privacy[1]
  var p_x_3_by_diagonal = partition(p_x_3.diagonal, diagonals_private)
  var p_x_3_by_tile = partition(p_x_3.tile, tiles_private)
  var [p_x_faces[3]] = cross_product(p_x_3_by_diagonal, p_x_3_by_tile)
  var s_x_3 = x_3_by_privacy[0]
  var s_x_3_by_diagonal = partition(s_x_3.diagonal, diagonals_shared)
  var [s_x_faces_by_tile[3]] = partition(s_x_3.tile, x_tiles_shared)
  var [s_x_faces[3]] = cross_product(s_x_3_by_diagonal, [s_x_faces_by_tile[3]])

  color_faces([x_faces[4]], Nx, Ny, Nz, ntx, nty, ntz, 0, array(true,false,false))
  var x_4_by_privacy = partition([x_faces[4]].is_private, ispace(int1d,2))
  var p_x_4 = x_4_by_privacy[1]
  var p_x_4_by_diagonal = partition(p_x_4.diagonal, diagonals_private)
  var p_x_4_by_tile = partition(p_x_4.tile, tiles_private)
  var [p_x_faces[4]] = cross_product(p_x_4_by_diagonal, p_x_4_by_tile)
  var s_x_4 = x_4_by_privacy[0]
  var s_x_4_by_diagonal = partition(s_x_4.diagonal, diagonals_shared)
  var [s_x_faces_by_tile[4]] = partition(s_x_4.tile, x_tiles_shared)
  var [s_x_faces[4]] = cross_product(s_x_4_by_diagonal, [s_x_faces_by_tile[4]])

  color_faces([x_faces[5]], Nx, Ny, Nz, ntx, nty, ntz, 0, array(false,true,true))
  var x_5_by_privacy = partition([x_faces[5]].is_private, ispace(int1d,2))
  var p_x_5 = x_5_by_privacy[1]
  var p_x_5_by_diagonal = partition(p_x_5.diagonal, diagonals_private)
  var p_x_5_by_tile = partition(p_x_5.tile, tiles_private)
  var [p_x_faces[5]] = cross_product(p_x_5_by_diagonal, p_x_5_by_tile)
  var s_x_5 = x_5_by_privacy[0]
  var s_x_5_by_diagonal = partition(s_x_5.diagonal, diagonals_shared)
  var [s_x_faces_by_tile[5]] = partition(s_x_5.tile, x_tiles_shared)
  var [s_x_faces[5]] = cross_product(s_x_5_by_diagonal, [s_x_faces_by_tile[5]])

  color_faces([x_faces[6]], Nx, Ny, Nz, ntx, nty, ntz, 0, array(false,true,false))
  var x_6_by_privacy = partition([x_faces[6]].is_private, ispace(int1d,2))
  var p_x_6 = x_6_by_privacy[1]
  var p_x_6_by_diagonal = partition(p_x_6.diagonal, diagonals_private)
  var p_x_6_by_tile = partition(p_x_6.tile, tiles_private)
  var [p_x_faces[6]] = cross_product(p_x_6_by_diagonal, p_x_6_by_tile)
  var s_x_6 = x_6_by_privacy[0]
  var s_x_6_by_diagonal = partition(s_x_6.diagonal, diagonals_shared)
  var [s_x_faces_by_tile[6]] = partition(s_x_6.tile, x_tiles_shared)
  var [s_x_faces[6]] = cross_product(s_x_6_by_diagonal, [s_x_faces_by_tile[6]])

  color_faces([x_faces[7]], Nx, Ny, Nz, ntx, nty, ntz, 0, array(false,false,true))
  var x_7_by_privacy = partition([x_faces[7]].is_private, ispace(int1d,2))
  var p_x_7 = x_7_by_privacy[1]
  var p_x_7_by_diagonal = partition(p_x_7.diagonal, diagonals_private)
  var p_x_7_by_tile = partition(p_x_7.tile, tiles_private)
  var [p_x_faces[7]] = cross_product(p_x_7_by_diagonal, p_x_7_by_tile)
  var s_x_7 = x_7_by_privacy[0]
  var s_x_7_by_diagonal = partition(s_x_7.diagonal, diagonals_shared)
  var [s_x_faces_by_tile[7]] = partition(s_x_7.tile, x_tiles_shared)
  var [s_x_faces[7]] = cross_product(s_x_7_by_diagonal, [s_x_faces_by_tile[7]])

  color_faces([x_faces[8]], Nx, Ny, Nz, ntx, nty, ntz, 0, array(false,false,false))
  var x_8_by_privacy = partition([x_faces[8]].is_private, ispace(int1d,2))
  var p_x_8 = x_8_by_privacy[1]
  var p_x_8_by_diagonal = partition(p_x_8.diagonal, diagonals_private)
  var p_x_8_by_tile = partition(p_x_8.tile, tiles_private)
  var [p_x_faces[8]] = cross_product(p_x_8_by_diagonal, p_x_8_by_tile)
  var s_x_8 = x_8_by_privacy[0]
  var s_x_8_by_diagonal = partition(s_x_8.diagonal, diagonals_shared)
  var [s_x_faces_by_tile[8]] = partition(s_x_8.tile, x_tiles_shared)
  var [s_x_faces[8]] = cross_product(s_x_8_by_diagonal, [s_x_faces_by_tile[8]])

  -- y

  color_faces([y_faces[1]], Nx, Ny, Nz, ntx, nty, ntz, 1, array(true,true,true))
  var y_1_by_privacy = partition([y_faces[1]].is_private, ispace(int1d,2))
  var p_y_1 = y_1_by_privacy[1]
  var p_y_1_by_diagonal = partition(p_y_1.diagonal, diagonals_private)
  var p_y_1_by_tile = partition(p_y_1.tile, tiles_private)
  var [p_y_faces[1]] = cross_product(p_y_1_by_diagonal, p_y_1_by_tile)
  var s_y_1 = y_1_by_privacy[0]
  var s_y_1_by_diagonal = partition(s_y_1.diagonal, diagonals_shared)
  var [s_y_faces_by_tile[1]] = partition(s_y_1.tile, y_tiles_shared)
  var [s_y_faces[1]] = cross_product(s_y_1_by_diagonal, [s_y_faces_by_tile[1]])

  color_faces([y_faces[2]], Nx, Ny, Nz, ntx, nty, ntz, 1, array(true,true,false))
  var y_2_by_privacy = partition([y_faces[2]].is_private, ispace(int1d,2))
  var p_y_2 = y_2_by_privacy[1]
  var p_y_2_by_diagonal = partition(p_y_2.diagonal, diagonals_private)
  var p_y_2_by_tile = partition(p_y_2.tile, tiles_private)
  var [p_y_faces[2]] = cross_product(p_y_2_by_diagonal, p_y_2_by_tile)
  var s_y_2 = y_2_by_privacy[0]
  var s_y_2_by_diagonal = partition(s_y_2.diagonal, diagonals_shared)
  var [s_y_faces_by_tile[2]] = partition(s_y_2.tile, y_tiles_shared)
  var [s_y_faces[2]] = cross_product(s_y_2_by_diagonal, [s_y_faces_by_tile[2]])

  color_faces([y_faces[3]], Nx, Ny, Nz, ntx, nty, ntz, 1, array(true,false,true))
  var y_3_by_privacy = partition([y_faces[3]].is_private, ispace(int1d,2))
  var p_y_3 = y_3_by_privacy[1]
  var p_y_3_by_diagonal = partition(p_y_3.diagonal, diagonals_private)
  var p_y_3_by_tile = partition(p_y_3.tile, tiles_private)
  var [p_y_faces[3]] = cross_product(p_y_3_by_diagonal, p_y_3_by_tile)
  var s_y_3 = y_3_by_privacy[0]
  var s_y_3_by_diagonal = partition(s_y_3.diagonal, diagonals_shared)
  var [s_y_faces_by_tile[3]] = partition(s_y_3.tile, y_tiles_shared)
  var [s_y_faces[3]] = cross_product(s_y_3_by_diagonal, [s_y_faces_by_tile[3]])

  color_faces([y_faces[4]], Nx, Ny, Nz, ntx, nty, ntz, 1, array(true,false,false))
  var y_4_by_privacy = partition([y_faces[4]].is_private, ispace(int1d,2))
  var p_y_4 = y_4_by_privacy[1]
  var p_y_4_by_diagonal = partition(p_y_4.diagonal, diagonals_private)
  var p_y_4_by_tile = partition(p_y_4.tile, tiles_private)
  var [p_y_faces[4]] = cross_product(p_y_4_by_diagonal, p_y_4_by_tile)
  var s_y_4 = y_4_by_privacy[0]
  var s_y_4_by_diagonal = partition(s_y_4.diagonal, diagonals_shared)
  var [s_y_faces_by_tile[4]] = partition(s_y_4.tile, y_tiles_shared)
  var [s_y_faces[4]] = cross_product(s_y_4_by_diagonal, [s_y_faces_by_tile[4]])

  color_faces([y_faces[5]], Nx, Ny, Nz, ntx, nty, ntz, 1, array(false,true,true))
  var y_5_by_privacy = partition([y_faces[5]].is_private, ispace(int1d,2))
  var p_y_5 = y_5_by_privacy[1]
  var p_y_5_by_diagonal = partition(p_y_5.diagonal, diagonals_private)
  var p_y_5_by_tile = partition(p_y_5.tile, tiles_private)
  var [p_y_faces[5]] = cross_product(p_y_5_by_diagonal, p_y_5_by_tile)
  var s_y_5 = y_5_by_privacy[0]
  var s_y_5_by_diagonal = partition(s_y_5.diagonal, diagonals_shared)
  var [s_y_faces_by_tile[5]] = partition(s_y_5.tile, y_tiles_shared)
  var [s_y_faces[5]] = cross_product(s_y_5_by_diagonal, [s_y_faces_by_tile[5]])

  color_faces([y_faces[6]], Nx, Ny, Nz, ntx, nty, ntz, 1, array(false,true,false))
  var y_6_by_privacy = partition([y_faces[6]].is_private, ispace(int1d,2))
  var p_y_6 = y_6_by_privacy[1]
  var p_y_6_by_diagonal = partition(p_y_6.diagonal, diagonals_private)
  var p_y_6_by_tile = partition(p_y_6.tile, tiles_private)
  var [p_y_faces[6]] = cross_product(p_y_6_by_diagonal, p_y_6_by_tile)
  var s_y_6 = y_6_by_privacy[0]
  var s_y_6_by_diagonal = partition(s_y_6.diagonal, diagonals_shared)
  var [s_y_faces_by_tile[6]] = partition(s_y_6.tile, y_tiles_shared)
  var [s_y_faces[6]] = cross_product(s_y_6_by_diagonal, [s_y_faces_by_tile[6]])

  color_faces([y_faces[7]], Nx, Ny, Nz, ntx, nty, ntz, 1, array(false,false,true))
  var y_7_by_privacy = partition([y_faces[7]].is_private, ispace(int1d,2))
  var p_y_7 = y_7_by_privacy[1]
  var p_y_7_by_diagonal = partition(p_y_7.diagonal, diagonals_private)
  var p_y_7_by_tile = partition(p_y_7.tile, tiles_private)
  var [p_y_faces[7]] = cross_product(p_y_7_by_diagonal, p_y_7_by_tile)
  var s_y_7 = y_7_by_privacy[0]
  var s_y_7_by_diagonal = partition(s_y_7.diagonal, diagonals_shared)
  var [s_y_faces_by_tile[7]] = partition(s_y_7.tile, y_tiles_shared)
  var [s_y_faces[7]] = cross_product(s_y_7_by_diagonal, [s_y_faces_by_tile[7]])

  color_faces([y_faces[8]], Nx, Ny, Nz, ntx, nty, ntz, 1, array(false,false,false))
  var y_8_by_privacy = partition([y_faces[8]].is_private, ispace(int1d,2))
  var p_y_8 = y_8_by_privacy[1]
  var p_y_8_by_diagonal = partition(p_y_8.diagonal, diagonals_private)
  var p_y_8_by_tile = partition(p_y_8.tile, tiles_private)
  var [p_y_faces[8]] = cross_product(p_y_8_by_diagonal, p_y_8_by_tile)
  var s_y_8 = y_8_by_privacy[0]
  var s_y_8_by_diagonal = partition(s_y_8.diagonal, diagonals_shared)
  var [s_y_faces_by_tile[8]] = partition(s_y_8.tile, y_tiles_shared)
  var [s_y_faces[8]] = cross_product(s_y_8_by_diagonal, [s_y_faces_by_tile[8]])

  -- z

  color_faces([z_faces[1]], Nx, Ny, Nz, ntx, nty, ntz, 2, array(true,true,true))
  var z_1_by_privacy = partition([z_faces[1]].is_private, ispace(int1d,2))
  var p_z_1 = z_1_by_privacy[1]
  var p_z_1_by_diagonal = partition(p_z_1.diagonal, diagonals_private)
  var p_z_1_by_tile = partition(p_z_1.tile, tiles_private)
  var [p_z_faces[1]] = cross_product(p_z_1_by_diagonal, p_z_1_by_tile)
  var s_z_1 = z_1_by_privacy[0]
  var s_z_1_by_diagonal = partition(s_z_1.diagonal, diagonals_shared)
  var [s_z_faces_by_tile[1]] = partition(s_z_1.tile, z_tiles_shared)
  var [s_z_faces[1]] = cross_product(s_z_1_by_diagonal, [s_z_faces_by_tile[1]])

  color_faces([z_faces[2]], Nx, Ny, Nz, ntx, nty, ntz, 2, array(true,true,false))
  var z_2_by_privacy = partition([z_faces[2]].is_private, ispace(int1d,2))
  var p_z_2 = z_2_by_privacy[1]
  var p_z_2_by_diagonal = partition(p_z_2.diagonal, diagonals_private)
  var p_z_2_by_tile = partition(p_z_2.tile, tiles_private)
  var [p_z_faces[2]] = cross_product(p_z_2_by_diagonal, p_z_2_by_tile)
  var s_z_2 = z_2_by_privacy[0]
  var s_z_2_by_diagonal = partition(s_z_2.diagonal, diagonals_shared)
  var [s_z_faces_by_tile[2]] = partition(s_z_2.tile, z_tiles_shared)
  var [s_z_faces[2]] = cross_product(s_z_2_by_diagonal, [s_z_faces_by_tile[2]])

  color_faces([z_faces[3]], Nx, Ny, Nz, ntx, nty, ntz, 2, array(true,false,true))
  var z_3_by_privacy = partition([z_faces[3]].is_private, ispace(int1d,2))
  var p_z_3 = z_3_by_privacy[1]
  var p_z_3_by_diagonal = partition(p_z_3.diagonal, diagonals_private)
  var p_z_3_by_tile = partition(p_z_3.tile, tiles_private)
  var [p_z_faces[3]] = cross_product(p_z_3_by_diagonal, p_z_3_by_tile)
  var s_z_3 = z_3_by_privacy[0]
  var s_z_3_by_diagonal = partition(s_z_3.diagonal, diagonals_shared)
  var [s_z_faces_by_tile[3]] = partition(s_z_3.tile, z_tiles_shared)
  var [s_z_faces[3]] = cross_product(s_z_3_by_diagonal, [s_z_faces_by_tile[3]])

  color_faces([z_faces[4]], Nx, Ny, Nz, ntx, nty, ntz, 2, array(true,false,false))
  var z_4_by_privacy = partition([z_faces[4]].is_private, ispace(int1d,2))
  var p_z_4 = z_4_by_privacy[1]
  var p_z_4_by_diagonal = partition(p_z_4.diagonal, diagonals_private)
  var p_z_4_by_tile = partition(p_z_4.tile, tiles_private)
  var [p_z_faces[4]] = cross_product(p_z_4_by_diagonal, p_z_4_by_tile)
  var s_z_4 = z_4_by_privacy[0]
  var s_z_4_by_diagonal = partition(s_z_4.diagonal, diagonals_shared)
  var [s_z_faces_by_tile[4]] = partition(s_z_4.tile, z_tiles_shared)
  var [s_z_faces[4]] = cross_product(s_z_4_by_diagonal, [s_z_faces_by_tile[4]])

  color_faces([z_faces[5]], Nx, Ny, Nz, ntx, nty, ntz, 2, array(false,true,true))
  var z_5_by_privacy = partition([z_faces[5]].is_private, ispace(int1d,2))
  var p_z_5 = z_5_by_privacy[1]
  var p_z_5_by_diagonal = partition(p_z_5.diagonal, diagonals_private)
  var p_z_5_by_tile = partition(p_z_5.tile, tiles_private)
  var [p_z_faces[5]] = cross_product(p_z_5_by_diagonal, p_z_5_by_tile)
  var s_z_5 = z_5_by_privacy[0]
  var s_z_5_by_diagonal = partition(s_z_5.diagonal, diagonals_shared)
  var [s_z_faces_by_tile[5]] = partition(s_z_5.tile, z_tiles_shared)
  var [s_z_faces[5]] = cross_product(s_z_5_by_diagonal, [s_z_faces_by_tile[5]])

  color_faces([z_faces[6]], Nx, Ny, Nz, ntx, nty, ntz, 2, array(false,true,false))
  var z_6_by_privacy = partition([z_faces[6]].is_private, ispace(int1d,2))
  var p_z_6 = z_6_by_privacy[1]
  var p_z_6_by_diagonal = partition(p_z_6.diagonal, diagonals_private)
  var p_z_6_by_tile = partition(p_z_6.tile, tiles_private)
  var [p_z_faces[6]] = cross_product(p_z_6_by_diagonal, p_z_6_by_tile)
  var s_z_6 = z_6_by_privacy[0]
  var s_z_6_by_diagonal = partition(s_z_6.diagonal, diagonals_shared)
  var [s_z_faces_by_tile[6]] = partition(s_z_6.tile, z_tiles_shared)
  var [s_z_faces[6]] = cross_product(s_z_6_by_diagonal, [s_z_faces_by_tile[6]])

  color_faces([z_faces[7]], Nx, Ny, Nz, ntx, nty, ntz, 2, array(false,false,true))
  var z_7_by_privacy = partition([z_faces[7]].is_private, ispace(int1d,2))
  var p_z_7 = z_7_by_privacy[1]
  var p_z_7_by_diagonal = partition(p_z_7.diagonal, diagonals_private)
  var p_z_7_by_tile = partition(p_z_7.tile, tiles_private)
  var [p_z_faces[7]] = cross_product(p_z_7_by_diagonal, p_z_7_by_tile)
  var s_z_7 = z_7_by_privacy[0]
  var s_z_7_by_diagonal = partition(s_z_7.diagonal, diagonals_shared)
  var [s_z_faces_by_tile[7]] = partition(s_z_7.tile, z_tiles_shared)
  var [s_z_faces[7]] = cross_product(s_z_7_by_diagonal, [s_z_faces_by_tile[7]])

  color_faces([z_faces[8]], Nx, Ny, Nz, ntx, nty, ntz, 2, array(false,false,false))
  var z_8_by_privacy = partition([z_faces[8]].is_private, ispace(int1d,2))
  var p_z_8 = z_8_by_privacy[1]
  var p_z_8_by_diagonal = partition(p_z_8.diagonal, diagonals_private)
  var p_z_8_by_tile = partition(p_z_8.tile, tiles_private)
  var [p_z_faces[8]] = cross_product(p_z_8_by_diagonal, p_z_8_by_tile)
  var s_z_8 = z_8_by_privacy[0]
  var s_z_8_by_diagonal = partition(s_z_8.diagonal, diagonals_shared)
  var [s_z_faces_by_tile[8]] = partition(s_z_8.tile, z_tiles_shared)
  var [s_z_faces[8]] = cross_product(s_z_8_by_diagonal, [s_z_faces_by_tile[8]])

end end

function Exports.InitRegions() return rquote

  -- Initialize face values
  initialize_faces([x_faces[1]])
  initialize_faces([x_faces[2]])
  initialize_faces([x_faces[3]])
  initialize_faces([x_faces[4]])
  initialize_faces([x_faces[5]])
  initialize_faces([x_faces[6]])
  initialize_faces([x_faces[7]])
  initialize_faces([x_faces[8]])

  initialize_faces([y_faces[1]])
  initialize_faces([y_faces[2]])
  initialize_faces([y_faces[3]])
  initialize_faces([y_faces[4]])
  initialize_faces([y_faces[5]])
  initialize_faces([y_faces[6]])
  initialize_faces([y_faces[7]])
  initialize_faces([y_faces[8]])

  initialize_faces([z_faces[1]])
  initialize_faces([z_faces[2]])
  initialize_faces([z_faces[3]])
  initialize_faces([z_faces[4]])
  initialize_faces([z_faces[5]])
  initialize_faces([z_faces[6]])
  initialize_faces([z_faces[7]])
  initialize_faces([z_faces[8]])

  -- Initialize constant values
  initialize_angles(angles)

end end

function Exports.ComputeRadiationField(config, tiles, p_points) return rquote

  var dx = config.Grid.xWidth / config.Radiation.xNum
  var dy = config.Grid.yWidth / config.Radiation.yNum
  var dz = config.Grid.zWidth / config.Radiation.zNum

  var iter = 1
  var omega = config.Radiation.qs/(config.Radiation.qa+config.Radiation.qs)

  -- Compute until convergence
  var res = 1.0
  while (res > tol) do

    -- Update the source term (in this problem, isotropic)
    for t in tiles do
      source_term(p_points[t], angles, omega)
    end

    -- Update the grid boundary intensities
    -- TODO: Should launch these on just the boundaries
    for j = 0, nty do
      for k = 0, ntz do
        bound_x_lo([s_x_faces_by_tile[1]][{0,j,k}],
                   [s_x_faces_by_tile[2]][{0,j,k}],
                   [s_x_faces_by_tile[3]][{0,j,k}],
                   [s_x_faces_by_tile[4]][{0,j,k}],
                   [s_x_faces_by_tile[5]][{0,j,k}],
                   [s_x_faces_by_tile[6]][{0,j,k}],
                   [s_x_faces_by_tile[7]][{0,j,k}],
                   [s_x_faces_by_tile[8]][{0,j,k}],
                   angles,
                   config.Radiation.emissWest,
                   config.Radiation.tempWest)

        bound_x_hi([s_x_faces_by_tile[1]][{ntx,j,k}],
                   [s_x_faces_by_tile[2]][{ntx,j,k}],
                   [s_x_faces_by_tile[3]][{ntx,j,k}],
                   [s_x_faces_by_tile[4]][{ntx,j,k}],
                   [s_x_faces_by_tile[5]][{ntx,j,k}],
                   [s_x_faces_by_tile[6]][{ntx,j,k}],
                   [s_x_faces_by_tile[7]][{ntx,j,k}],
                   [s_x_faces_by_tile[8]][{ntx,j,k}],
                   angles,
                   config.Radiation.emissEast,
                   config.Radiation.tempEast)
      end
    end

    -- Update y faces
    for i = 0, ntx do
      for k = 0, ntz do
        bound_y_lo([s_y_faces_by_tile[1]][{i,0,k}],
                   [s_y_faces_by_tile[2]][{i,0,k}],
                   [s_y_faces_by_tile[3]][{i,0,k}],
                   [s_y_faces_by_tile[4]][{i,0,k}],
                   [s_y_faces_by_tile[5]][{i,0,k}],
                   [s_y_faces_by_tile[6]][{i,0,k}],
                   [s_y_faces_by_tile[7]][{i,0,k}],
                   [s_y_faces_by_tile[8]][{i,0,k}],
                   angles,
                   config.Radiation.emissSouth,
                   config.Radiation.tempSouth)

        bound_y_hi([s_y_faces_by_tile[1]][{i,nty,k}],
                   [s_y_faces_by_tile[2]][{i,nty,k}],
                   [s_y_faces_by_tile[3]][{i,nty,k}],
                   [s_y_faces_by_tile[4]][{i,nty,k}],
                   [s_y_faces_by_tile[5]][{i,nty,k}],
                   [s_y_faces_by_tile[6]][{i,nty,k}],
                   [s_y_faces_by_tile[7]][{i,nty,k}],
                   [s_y_faces_by_tile[8]][{i,nty,k}],
                   angles,
                   config.Radiation.emissNorth,
                   config.Radiation.tempNorth)
      end
    end

    -- Update z faces
    for i = 0, ntx do
      for j = 0, nty do
        bound_z_lo([s_z_faces_by_tile[1]][{i,j,0}],
                   [s_z_faces_by_tile[2]][{i,j,0}],
                   [s_z_faces_by_tile[3]][{i,j,0}],
                   [s_z_faces_by_tile[4]][{i,j,0}],
                   [s_z_faces_by_tile[5]][{i,j,0}],
                   [s_z_faces_by_tile[6]][{i,j,0}],
                   [s_z_faces_by_tile[7]][{i,j,0}],
                   [s_z_faces_by_tile[8]][{i,j,0}],
                   angles,
                   config.Radiation.emissDown,
                   config.Radiation.tempDown)

        bound_z_hi([s_z_faces_by_tile[1]][{i,j,ntz}],
                   [s_z_faces_by_tile[2]][{i,j,ntz}],
                   [s_z_faces_by_tile[3]][{i,j,ntz}],
                   [s_z_faces_by_tile[4]][{i,j,ntz}],
                   [s_z_faces_by_tile[5]][{i,j,ntz}],
                   [s_z_faces_by_tile[6]][{i,j,ntz}],
                   [s_z_faces_by_tile[7]][{i,j,ntz}],
                   [s_z_faces_by_tile[8]][{i,j,ntz}],
                   angles,
                   config.Radiation.emissUp,
                   config.Radiation.tempUp)
      end
    end

    --Perform the sweep for computing new intensities
    --Quadrant 1 - +x, +y, +z
    for d in diagonals_private do
      __demand(__parallel)
      for t in [p_x_faces[1]][d].colors do
        sweep_1(p_points[t],
                [p_x_faces[1]][d][t], [p_y_faces[1]][d][t], [p_z_faces[1]][d][t],
                [s_x_faces[1]][d][t], [s_x_faces[1]][d+1][t+int3d{1,0,0}],
                [s_y_faces[1]][d][t], [s_y_faces[1]][d+1][t+int3d{0,1,0}],
                [s_z_faces[1]][d][t], [s_z_faces[1]][d+1][t+int3d{0,0,1}],
                angles, 1, 1, 1, dx, dy, dz)
      end
    end

    -- Quadrant 2 - +x, +y, -z
    for d in diagonals_private do
      __demand(__parallel)
      for t in [p_x_faces[2]][d].colors do
        sweep_2(p_points[t],
                [p_x_faces[2]][d][t], [p_y_faces[2]][d][t], [p_z_faces[2]][d][t],
                [s_x_faces[2]][d][t], [s_x_faces[2]][d+1][t+int3d{1,0,0}],
                [s_y_faces[2]][d][t], [s_y_faces[2]][d+1][t+int3d{0,1,0}],
                [s_z_faces[2]][d][t+int3d{0,0,1}], [s_z_faces[2]][d+1][t],
                angles, 1, 1, -1, dx, dy, dz)
      end
    end

    -- Quadrant 3 - +x, -y, +z
    for d in diagonals_private do
      __demand(__parallel)
      for t in [p_x_faces[3]][d].colors do
        sweep_3(p_points[t],
                [p_x_faces[3]][d][t], [p_y_faces[3]][d][t], [p_z_faces[3]][d][t],
                [s_x_faces[3]][d][t], [s_x_faces[3]][d+1][t+int3d{1,0,0}],
                [s_y_faces[3]][d][t+int3d{0,1,0}], [s_y_faces[3]][d+1][t],
                [s_z_faces[3]][d][t], [s_z_faces[3]][d+1][t+int3d{0,0,1}],
                angles, 1, -1, 1, dx, dy, dz)
      end
    end

    -- Quadrant 4 - +x, -y, -z
    for d in diagonals_private do
      __demand(__parallel)
      for t in [p_x_faces[4]][d].colors do
        sweep_4(p_points[t],
                [p_x_faces[4]][d][t], [p_y_faces[4]][d][t], [p_z_faces[4]][d][t],
                [s_x_faces[4]][d][t], [s_x_faces[4]][d+1][t+int3d{1,0,0}],
                [s_y_faces[4]][d][t+int3d{0,1,0}], [s_y_faces[4]][d+1][t],
                [s_z_faces[4]][d][t+int3d{0,0,1}], [s_z_faces[4]][d+1][t],
                angles, 1, -1, -1, dx, dy, dz)
      end
    end

    -- Quadrant 5 - -x, +y, +z
    for d in diagonals_private do
      __demand(__parallel)
      for t in [p_x_faces[5]][d].colors do
        sweep_5(p_points[t],
                [p_x_faces[5]][d][t], [p_y_faces[5]][d][t], [p_z_faces[5]][d][t],
                [s_x_faces[5]][d][t+int3d{1,0,0}], [s_x_faces[5]][d+1][t],
                [s_y_faces[5]][d][t], [s_y_faces[5]][d+1][t+int3d{0,1,0}],
                [s_z_faces[5]][d][t], [s_z_faces[5]][d+1][t+int3d{0,0,1}],
                angles, -1, 1, 1, dx, dy, dz)
      end
    end

    -- Quadrant 6 - -x, +y, -z
    for d in diagonals_private do
      __demand(__parallel)
      for t in [p_x_faces[6]][d].colors do
        sweep_6(p_points[t],
                [p_x_faces[6]][d][t], [p_y_faces[6]][d][t], [p_z_faces[6]][d][t],
                [s_x_faces[6]][d][t+int3d{1,0,0}], [s_x_faces[6]][d+1][t],
                [s_y_faces[6]][d][t], [s_y_faces[6]][d+1][t+int3d{0,1,0}],
                [s_z_faces[6]][d][t+int3d{0,0,1}], [s_z_faces[6]][d+1][t],
                angles, -1, 1, -1, dx, dy, dz)
      end
    end

    -- Quadrant 7 - -x, -y, +z
    for d in diagonals_private do
      __demand(__parallel)
      for t in [p_x_faces[7]][d].colors do
        sweep_7(p_points[t],
                [p_x_faces[7]][d][t], [p_y_faces[7]][d][t], [p_z_faces[7]][d][t],
                [s_x_faces[7]][d][t+int3d{1,0,0}], [s_x_faces[7]][d+1][t],
                [s_y_faces[7]][d][t+int3d{0,1,0}], [s_y_faces[7]][d+1][t],
                [s_z_faces[7]][d][t], [s_z_faces[7]][d+1][t+int3d{0,0,1}],
                angles, -1, -1, 1, dx, dy, dz)
      end
    end

    -- Quadrant 8 - -x, -y, -z
    for d in diagonals_private do
      __demand(__parallel)
      for t in [p_x_faces[8]][d].colors do
        sweep_8(p_points[t],
                [p_x_faces[8]][d][t], [p_y_faces[8]][d][t], [p_z_faces[8]][d][t],
                [s_x_faces[8]][d][t+int3d{1,0,0}], [s_x_faces[8]][d+1][t],
                [s_y_faces[8]][d][t+int3d{0,1,0}], [s_y_faces[8]][d+1][t],
                [s_z_faces[8]][d][t+int3d{0,0,1}], [s_z_faces[8]][d+1][t],
                angles, -1, -1, -1, dx, dy, dz)
      end
    end

    -- Compute the residual
    res = 0.0
    for t in tiles do
      res += residual(p_points[t], Nx, Ny, Nz)
    end
    res = sqrt(res)

    -- Update the intensities and the iteration number
    for t in tiles do
      update(p_points[t])
    end

    if (iter == 1) then
      C.printf("\n")
      C.printf(" Iteration     Residual         \n")
      C.printf(" ------------------------------ \n")
    end
    C.printf( "   %3d    %.15e \n", iter, res)

    iter += 1

  end

  -- Reduce intensity
  for t in tiles do
    reduce_intensity(p_points[t], angles)
  end

end end

-------------------------------------------------------------------------------
-- MODULE EXPORTS
-------------------------------------------------------------------------------

return Exports
end
