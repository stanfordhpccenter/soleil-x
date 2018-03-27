import 'regent'

-------------------------------------------------------------------------------
-- MODULE PARAMETERS
-------------------------------------------------------------------------------

return function(NUM_ANGLES, pointsFSpace)

-------------------------------------------------------------------------------
-- COMPILE-TIME COMPUTATION
-------------------------------------------------------------------------------

-- C imports

local c = regentlib.c

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

local terra open_quad_file() : &c.FILE
  var f = c.fopen([quad_file], 'rb')
  if f == nil then
    c.printf('Error opening angle file\n')
    c.exit(1)
  end
  return f
end

local terra read_double(f : &c.FILE) : double
  var val : double
  if c.fscanf(f, '%lf\n', &val) < 1 then
    c.printf('Error while reading angle file\n')
    c.exit(1)
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
  private_color : int3d,      -- Used for partition_by_field
  shared_color : int3d,       -- Used for partition_by_field
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

  c.fclose(f)

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

local task color_faces_x(faces : region(ispace(int3d), face),
                                        Nx : int, Ny : int, Nz : int,
                                        ntx : int, nty : int, ntz : int)
  where
    reads writes (faces)
  do

  var limits = faces.bounds

  for i = limits.lo.x, limits.hi.x + 1 do
    var x_tile = i / (Nx/ntx)
    for j = limits.lo.y, limits.hi.y + 1 do
      var y_tile = j / (Ny/nty)
      for k = limits.lo.z, limits.hi.z + 1 do
        var z_tile = k / (Nz/ntz)

        var color = {x = x_tile, y = y_tile, z = z_tile}
        if i % (Nx/ntx) == 0 then
          faces[{i,j,k}].shared_color = color
          faces[{i,j,k}].private_color = {x=-1, y=-1, z=-1}
        else
          faces[{i,j,k}].shared_color = {x=-1, y=-1, z=-1}
          faces[{i,j,k}].private_color = color
        end
      end
    end
  end
end

local task color_faces_y(faces : region(ispace(int3d), face),
                                        Nx : int, Ny : int, Nz : int,
                                        ntx : int, nty : int, ntz : int)
  where
    reads writes (faces)
  do

  var limits = faces.bounds

  for i = limits.lo.x, limits.hi.x + 1 do
    var x_tile = i / (Nx/ntx)
    for j = limits.lo.y, limits.hi.y + 1 do
      var y_tile = j / (Ny/nty)
      for k = limits.lo.z, limits.hi.z + 1 do
        var z_tile = k / (Nz/ntz)

        var color = {x = x_tile, y = y_tile, z = z_tile}
        if j % (Ny/nty) == 0 then
          faces[{i,j,k}].shared_color = color
          faces[{i,j,k}].private_color = {x=-1, y=-1, z=-1}
        else
          faces[{i,j,k}].shared_color = {x=-1, y=-1, z=-1}
          faces[{i,j,k}].private_color = color
        end
      end
    end
  end
end

local task color_faces_z(faces : region(ispace(int3d), face),
                                        Nx : int, Ny : int, Nz : int,
                                        ntx : int, nty : int, ntz : int)
  where
    reads writes (faces)
  do

  var limits = faces.bounds

  for i = limits.lo.x, limits.hi.x + 1 do
    var x_tile = i / (Nx/ntx)
    for j = limits.lo.y, limits.hi.y + 1 do
      var y_tile = j / (Ny/nty)
      for k = limits.lo.z, limits.hi.z + 1 do
        var z_tile = k / (Nz/ntz)

        var color = {x = x_tile, y = y_tile, z = z_tile}
        if k % (Nz/ntz) == 0 then
          faces[{i,j,k}].shared_color = color
          faces[{i,j,k}].private_color = {x=-1, y=-1, z=-1}
        else
          faces[{i,j,k}].shared_color = {x=-1, y=-1, z=-1}
          faces[{i,j,k}].private_color = color
        end
      end
    end
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

  var value : double = epsw*SB*pow(Tw,4.0)/pi

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

  var value : double = epsw*SB*pow(Tw,4.0)/pi

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

  var value : double = epsw*SB*pow(Tw,4.0)/pi

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

  var value : double = epsw*SB*pow(Tw,4.0)/pi

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

  var value : double = epsw*SB*pow(Tw,4.0)/pi

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

  var value : double = epsw*SB*pow(Tw,4.0)/pi

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

local task write_intensity(points : region(ispace(int3d), pointsFSpace))
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
local task print_final_intensities(points : region(ispace(int3d), pointsFSpace))
where
  reads (points.{I_1, I_2, I_3, I_4, I_5, I_6, I_7, I_8})

do

    var limits = points.bounds
    for i = limits.lo.x, limits.hi.x+1 do
      for j = limits.lo.y, limits.hi.y+1 do
        for k = limits.lo.z, limits.hi.z+1 do
            for m = 0, NUM_ANGLES do

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

local s_x_faces = {}
local s_y_faces = {}
local s_z_faces = {}

local p_x_faces = {}
local p_y_faces = {}
local p_z_faces = {}

p_x_faces[1] = regentlib.newsymbol('p_x_faces_1')
p_y_faces[1] = regentlib.newsymbol('p_y_faces_1')
p_z_faces[1] = regentlib.newsymbol('p_z_faces_1')

s_x_faces[1] = regentlib.newsymbol('s_x_faces_1')
s_y_faces[1] = regentlib.newsymbol('s_y_faces_1')
s_z_faces[1] = regentlib.newsymbol('s_z_faces_1')

p_x_faces[2] = regentlib.newsymbol('p_x_faces_2')
p_y_faces[2] = regentlib.newsymbol('p_y_faces_2')
p_z_faces[2] = regentlib.newsymbol('p_z_faces_2')

s_x_faces[2] = regentlib.newsymbol('s_x_faces_2')
s_y_faces[2] = regentlib.newsymbol('s_y_faces_2')
s_z_faces[2] = regentlib.newsymbol('s_z_faces_2')

p_x_faces[3] = regentlib.newsymbol('p_x_faces_3')
p_y_faces[3] = regentlib.newsymbol('p_y_faces_3')
p_z_faces[3] = regentlib.newsymbol('p_z_faces_3')

s_x_faces[3] = regentlib.newsymbol('s_x_faces_3')
s_y_faces[3] = regentlib.newsymbol('s_y_faces_3')
s_z_faces[3] = regentlib.newsymbol('s_z_faces_3')

p_x_faces[4] = regentlib.newsymbol('p_x_faces_4')
p_y_faces[4] = regentlib.newsymbol('p_y_faces_4')
p_z_faces[4] = regentlib.newsymbol('p_z_faces_4')

s_x_faces[4] = regentlib.newsymbol('s_x_faces_4')
s_y_faces[4] = regentlib.newsymbol('s_y_faces_4')
s_z_faces[4] = regentlib.newsymbol('s_z_faces_4')

p_x_faces[5] = regentlib.newsymbol('p_x_faces_5')
p_y_faces[5] = regentlib.newsymbol('p_y_faces_5')
p_z_faces[5] = regentlib.newsymbol('p_z_faces_5')

s_x_faces[5] = regentlib.newsymbol('s_x_faces_5')
s_y_faces[5] = regentlib.newsymbol('s_y_faces_5')
s_z_faces[5] = regentlib.newsymbol('s_z_faces_5')

p_x_faces[6] = regentlib.newsymbol('p_x_faces_6')
p_y_faces[6] = regentlib.newsymbol('p_y_faces_6')
p_z_faces[6] = regentlib.newsymbol('p_z_faces_6')

s_x_faces[6] = regentlib.newsymbol('s_x_faces_6')
s_y_faces[6] = regentlib.newsymbol('s_y_faces_6')
s_z_faces[6] = regentlib.newsymbol('s_z_faces_6')

p_x_faces[7] = regentlib.newsymbol('p_x_faces_7')
p_y_faces[7] = regentlib.newsymbol('p_y_faces_7')
p_z_faces[7] = regentlib.newsymbol('p_z_faces_7')

s_x_faces[7] = regentlib.newsymbol('s_x_faces_7')
s_y_faces[7] = regentlib.newsymbol('s_y_faces_7')
s_z_faces[7] = regentlib.newsymbol('s_z_faces_7')

p_x_faces[8] = regentlib.newsymbol('p_x_faces_8')
p_y_faces[8] = regentlib.newsymbol('p_y_faces_8')
p_z_faces[8] = regentlib.newsymbol('p_z_faces_8')

s_x_faces[8] = regentlib.newsymbol('s_x_faces_8')
s_y_faces[8] = regentlib.newsymbol('s_y_faces_8')
s_z_faces[8] = regentlib.newsymbol('s_z_faces_8')

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
  var grid_x = ispace(int3d, {x = Nx+1, y = Ny,   z = Nz})
  var grid_y = ispace(int3d, {x = Nx,   y = Ny+1, z = Nz})
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
  var [tiles_private] = ispace(int3d, {x = ntx, y = nty,   z = ntz  })
  var [x_tiles_shared] = ispace(int3d, {x = ntx+1, y = nty,   z = ntz  })
  var [y_tiles_shared] = ispace(int3d, {x = ntx,   y = nty+1, z = ntz  })
  var [z_tiles_shared] = ispace(int3d, {x = ntx,   y = nty,   z = ntz+1})

  c.printf("ntx = %d, nty = %d, ntz = %d \n", ntx, nty, ntz)
  c.printf("Nx = %d, Ny = %d, Nz = %d \n", Nx, Ny, Nz)

  -- x
  
  color_faces_x([x_faces[1]], Nx, Ny, Nz, ntx, nty, ntz)         
  var [p_x_faces[1]] = partition([x_faces[1]].private_color, [tiles_private])
  var [s_x_faces[1]] = partition([x_faces[1]].shared_color, [x_tiles_shared])

  color_faces_x([x_faces[2]], Nx, Ny, Nz, ntx, nty, ntz)
  var [p_x_faces[2]] = partition([x_faces[2]].private_color, [tiles_private])
  var [s_x_faces[2]] = partition([x_faces[2]].shared_color, [x_tiles_shared])

  color_faces_x([x_faces[3]], Nx, Ny, Nz, ntx, nty, ntz)
  var [p_x_faces[3]] = partition([x_faces[3]].private_color, [tiles_private])
  var [s_x_faces[3]] = partition([x_faces[3]].shared_color, [x_tiles_shared])

  color_faces_x([x_faces[4]], Nx, Ny, Nz, ntx, nty, ntz)
  var [p_x_faces[4]] = partition([x_faces[4]].private_color, [tiles_private])
  var [s_x_faces[4]] = partition([x_faces[4]].shared_color, [x_tiles_shared])

  color_faces_x([x_faces[5]], Nx, Ny, Nz, ntx, nty, ntz)
  var [p_x_faces[5]] = partition([x_faces[5]].private_color, [tiles_private])
  var [s_x_faces[5]] = partition([x_faces[5]].shared_color, [x_tiles_shared])

  color_faces_x([x_faces[6]], Nx, Ny, Nz, ntx, nty, ntz)
  var [p_x_faces[6]] = partition([x_faces[6]].private_color, [tiles_private])
  var [s_x_faces[6]] = partition([x_faces[6]].shared_color, [x_tiles_shared])

  color_faces_x([x_faces[7]], Nx, Ny, Nz, ntx, nty, ntz)
  var [p_x_faces[7]] = partition([x_faces[7]].private_color, [tiles_private])
  var [s_x_faces[7]] = partition([x_faces[7]].shared_color, [x_tiles_shared])

  color_faces_x([x_faces[8]], Nx, Ny, Nz, ntx, nty, ntz)
  var [p_x_faces[8]] = partition([x_faces[8]].private_color, [tiles_private])
  var [s_x_faces[8]] = partition([x_faces[8]].shared_color, [x_tiles_shared])

  -- y

  color_faces_y([y_faces[1]], Nx, Ny, Nz, ntx, nty, ntz)
  var [p_y_faces[1]] = partition([y_faces[1]].private_color, [tiles_private])
  var [s_y_faces[1]] = partition([y_faces[1]].shared_color, [y_tiles_shared])

  color_faces_y([y_faces[2]], Nx, Ny, Nz, ntx, nty, ntz)
  var [p_y_faces[2]] = partition([y_faces[2]].private_color, [tiles_private])
  var [s_y_faces[2]] = partition([y_faces[2]].shared_color, [y_tiles_shared])

  color_faces_y([y_faces[3]], Nx, Ny, Nz, ntx, nty, ntz)
  var [p_y_faces[3]] = partition([y_faces[3]].private_color, [tiles_private])
  var [s_y_faces[3]] = partition([y_faces[3]].shared_color, [y_tiles_shared])

  color_faces_y([y_faces[4]], Nx, Ny, Nz, ntx, nty, ntz)
  var [p_y_faces[4]] = partition([y_faces[4]].private_color, [tiles_private])
  var [s_y_faces[4]] = partition([y_faces[4]].shared_color, [y_tiles_shared])

  color_faces_y([y_faces[5]], Nx, Ny, Nz, ntx, nty, ntz)
  var [p_y_faces[5]] = partition([y_faces[5]].private_color, [tiles_private])
  var [s_y_faces[5]] = partition([y_faces[5]].shared_color, [y_tiles_shared])

  color_faces_y([y_faces[6]], Nx, Ny, Nz, ntx, nty, ntz)
  var [p_y_faces[6]] = partition([y_faces[6]].private_color, [tiles_private])
  var [s_y_faces[6]] = partition([y_faces[6]].shared_color, [y_tiles_shared])

  color_faces_y([y_faces[7]], Nx, Ny, Nz, ntx, nty, ntz)
  var [p_y_faces[7]] = partition([y_faces[7]].private_color, [tiles_private])
  var [s_y_faces[7]] = partition([y_faces[7]].shared_color, [y_tiles_shared])

  color_faces_y([y_faces[8]], Nx, Ny, Nz, ntx, nty, ntz)
  var [p_y_faces[8]] = partition([y_faces[8]].private_color, [tiles_private])
  var [s_y_faces[8]] = partition([y_faces[8]].shared_color, [y_tiles_shared])

  -- z

  color_faces_z([z_faces[1]], Nx, Ny, Nz, ntx, nty, ntz)
  var [p_z_faces[1]] = partition([z_faces[1]].private_color, [tiles_private])
  var [s_z_faces[1]] = partition([z_faces[1]].shared_color, [z_tiles_shared])

  color_faces_z([z_faces[2]], Nx, Ny, Nz, ntx, nty, ntz)
  var [p_z_faces[2]] = partition([z_faces[2]].private_color, [tiles_private])
  var [s_z_faces[2]] = partition([z_faces[2]].shared_color, [z_tiles_shared])

  color_faces_z([z_faces[3]], Nx, Ny, Nz, ntx, nty, ntz)
  var [p_z_faces[3]] = partition([z_faces[3]].private_color, [tiles_private])
  var [s_z_faces[3]] = partition([z_faces[3]].shared_color, [z_tiles_shared])

  color_faces_z([z_faces[4]], Nx, Ny, Nz, ntx, nty, ntz)
  var [p_z_faces[4]] = partition([z_faces[4]].private_color, [tiles_private])
  var [s_z_faces[4]] = partition([z_faces[4]].shared_color, [z_tiles_shared])

  color_faces_z([z_faces[5]], Nx, Ny, Nz, ntx, nty, ntz)
  var [p_z_faces[5]] = partition([z_faces[5]].private_color, [tiles_private])
  var [s_z_faces[5]] = partition([z_faces[5]].shared_color, [z_tiles_shared])

  color_faces_z([z_faces[6]], Nx, Ny, Nz, ntx, nty, ntz)
  var [p_z_faces[6]] = partition([z_faces[6]].private_color, [tiles_private])
  var [s_z_faces[6]] = partition([z_faces[6]].shared_color, [z_tiles_shared])

  color_faces_z([z_faces[7]], Nx, Ny, Nz, ntx, nty, ntz)
  var [p_z_faces[7]] = partition([z_faces[7]].private_color, [tiles_private])
  var [s_z_faces[7]] = partition([z_faces[7]].shared_color, [z_tiles_shared])

  color_faces_z([z_faces[8]], Nx, Ny, Nz, ntx, nty, ntz)
  var [p_z_faces[8]] = partition([z_faces[8]].private_color, [tiles_private])
  var [s_z_faces[8]] = partition([z_faces[8]].shared_color, [z_tiles_shared])

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

  var t : int64  = 1
  var omega = config.Radiation.qs/(config.Radiation.qa+config.Radiation.qs)

  -- Compute until convergence
  var res = 1.0
  while (res > tol) do

    -- Update the source term (in this problem, isotropic)
    for color in tiles do
      source_term(p_points[color], angles, omega)
    end

    -- Update the grid boundary intensities
    for j = 0, nty do
      for k = 0, ntz do
        west_bound([s_x_faces[1]][{0,j,k}],
                   [s_x_faces[2]][{0,j,k}],
                   [s_x_faces[3]][{0,j,k}],
                   [s_x_faces[4]][{0,j,k}],
                   [s_x_faces[5]][{0,j,k}],
                   [s_x_faces[6]][{0,j,k}],
                   [s_x_faces[7]][{0,j,k}],
                   [s_x_faces[8]][{0,j,k}],
                   angles,
                   config.Radiation.emissWest,
                   config.Radiation.tempWest)

        east_bound([s_x_faces[1]][{ntx,j,k}],
                   [s_x_faces[2]][{ntx,j,k}],
                   [s_x_faces[3]][{ntx,j,k}],
                   [s_x_faces[4]][{ntx,j,k}],
                   [s_x_faces[5]][{ntx,j,k}],
                   [s_x_faces[6]][{ntx,j,k}],
                   [s_x_faces[7]][{ntx,j,k}],
                   [s_x_faces[8]][{ntx,j,k}],
                   angles,
                   config.Radiation.emissEast,
                   config.Radiation.tempEast)
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
                    angles,
                    config.Radiation.emissSouth,
                    config.Radiation.tempSouth)

        north_bound([s_y_faces[1]][{i,nty,k}],
                    [s_y_faces[2]][{i,nty,k}],
                    [s_y_faces[3]][{i,nty,k}],
                    [s_y_faces[4]][{i,nty,k}],
                    [s_y_faces[5]][{i,nty,k}],
                    [s_y_faces[6]][{i,nty,k}],
                    [s_y_faces[7]][{i,nty,k}],
                    [s_y_faces[8]][{i,nty,k}],
                    angles,
                    config.Radiation.emissNorth,
                    config.Radiation.tempNorth)
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
                   angles,
                   config.Radiation.emissUp,
                   config.Radiation.tempUp)

        down_bound([s_z_faces[1]][{i,j,ntz}],
                   [s_z_faces[2]][{i,j,ntz}],
                   [s_z_faces[3]][{i,j,ntz}], 
                   [s_z_faces[4]][{i,j,ntz}], 
                   [s_z_faces[5]][{i,j,ntz}], 
                   [s_z_faces[6]][{i,j,ntz}], 
                   [s_z_faces[7]][{i,j,ntz}], 
                   [s_z_faces[8]][{i,j,ntz}], 
                   angles,
                   config.Radiation.emissDown,
                   config.Radiation.tempDown)
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
    for color in tiles do
      res += residual(p_points[color], Nx, Ny, Nz)
    end
    res = sqrt(res)

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
  -- write_intensity(p_points)
  -- print_final_intensities(p_points)

end end

-------------------------------------------------------------------------------
-- MODULE EXPORTS
-------------------------------------------------------------------------------

return Exports
end
