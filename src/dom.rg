import 'regent'

-- Coloring example: 6x6 square, 2x2 tiling, +1-1 direction
--
-- Internal cell values are essentially private,
-- face values are what need to be passed to downstream neighbor.
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

-------------------------------------------------------------------------------
-- MODULE PARAMETERS
-------------------------------------------------------------------------------

return function(MAX_ANGLES_PER_QUAD, pointsFSpace, Config) local MODULE = {}

-------------------------------------------------------------------------------
-- IMPORTS
-------------------------------------------------------------------------------

local C = regentlib.c
local UTIL = require 'util'

local fabs = regentlib.fabs(double)
local max = regentlib.fmax
local min = regentlib.fmin
local pow = regentlib.pow(double)
local sqrt = regentlib.sqrt(double)

-------------------------------------------------------------------------------
-- CONSTANTS
-------------------------------------------------------------------------------

local PI = 3.1415926535898
local SB = 5.67e-8
local TOLERANCE = 1e-6 -- solution tolerance
local GAMMA = 0.5 -- 1 for step differencing, 0.5 for diamond differencing

-------------------------------------------------------------------------------
-- HELPER FUNCTIONS
-------------------------------------------------------------------------------

local terra open_quad_file(num_angles : int) : &C.FILE
  var fname = [&int8](C.malloc(256))
  C.snprintf(fname, 256, '%s/src/LMquads/%d.txt',
             C.getenv('SOLEIL_DIR'), num_angles)
  var f = UTIL.openFile(fname, 'rb')
  C.free(fname)
  return f
end

local terra read_double(f : &C.FILE) : double
  var val : double
  if C.fscanf(f, '%lf\n', &val) < 1 then
    var stderr = C.fdopen(2, 'w')
    C.fprintf(stderr, 'Error while reading angle file\n')
    C.fflush(stderr)
    C.exit(1)
  end
  return val
end

-------------------------------------------------------------------------------
-- MODULE-LOCAL FIELD SPACES
-------------------------------------------------------------------------------

local struct angle {
  xi  : double,
  eta : double,
  mu  : double,
  w   : double,
}

local struct face {
  I : double[MAX_ANGLES_PER_QUAD],
  is_private : int1d, -- 1 = private, 0 = shared
  tile : int3d,
}

-------------------------------------------------------------------------------
-- QUADRANT MACROS
-------------------------------------------------------------------------------

local directions = terralib.newlist{
  terralib.newlist{ true,  true,  true},
  terralib.newlist{ true,  true, false},
  terralib.newlist{ true, false,  true},
  terralib.newlist{ true, false, false},
  terralib.newlist{false,  true,  true},
  terralib.newlist{false,  true, false},
  terralib.newlist{false, false,  true},
  terralib.newlist{false, false, false},
}

local intensityFields = terralib.newlist{
  'I_1', 'I_2', 'I_3', 'I_4', 'I_5', 'I_6', 'I_7', 'I_8'
}

local iterIntFields = terralib.newlist{
  'Iiter_1', 'Iiter_2', 'Iiter_3', 'Iiter_4', 'Iiter_5', 'Iiter_6', 'Iiter_7', 'Iiter_8'
}

-- 1..8, regentlib.rexpr -> regentlib.rexpr
local function angleInQuadrant(q, angle)
  return rexpr
    [directions[q][1] and rexpr angle.xi  >= 0 end or rexpr angle.xi  <= 0 end] and
    [directions[q][2] and rexpr angle.eta >= 0 end or rexpr angle.eta <= 0 end] and
    [directions[q][3] and rexpr angle.mu  >= 0 end or rexpr angle.mu  <= 0 end]
  end
end

local __demand(__inline)
task quadrantSize(q : int, num_angles : int)
  return num_angles/8 + max(0, min(1, num_angles%8 - q + 1))
end

-------------------------------------------------------------------------------
-- MODULE-LOCAL TASKS
-------------------------------------------------------------------------------

local angles = UTIL.generate(8, function()
  return regentlib.newsymbol(region(ispace(int1d), angle))
end)

local -- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task initialize_angles([angles],
                       config : Config)
where
  [angles:map(function(a) return terralib.newlist{
     regentlib.privilege(regentlib.reads, a, 'xi'),
     regentlib.privilege(regentlib.reads, a, 'eta'),
     regentlib.privilege(regentlib.reads, a, 'mu'),
     regentlib.privilege(regentlib.reads, a, 'w'),
     regentlib.privilege(regentlib.writes, a, 'xi'),
     regentlib.privilege(regentlib.writes, a, 'eta'),
     regentlib.privilege(regentlib.writes, a, 'mu'),
     regentlib.privilege(regentlib.writes, a, 'w'),
   } end):flatten()]
do
  -- Open angles file
  var num_angles = config.Radiation.angles
  regentlib.assert(
    MAX_ANGLES_PER_QUAD * 8 >= num_angles,
    'Too many angles; recompile with larger MAX_ANGLES_PER_QUAD')
  var f = open_quad_file(num_angles)
  -- Throw away num angles header
  read_double(f)
  -- Read fields round-robin into angle quadrants
  for m = 0, MAX_ANGLES_PER_QUAD do
    @ESCAPE for q = 1, 8 do @EMIT
      if m*8 + q - 1 == num_angles then break end
      [angles[q]][m].xi = read_double(f)
    @TIME end @EPACSE
  end
  for m = 0, MAX_ANGLES_PER_QUAD do
    @ESCAPE for q = 1, 8 do @EMIT
      if m*8 + q - 1 == num_angles then break end
      [angles[q]][m].eta = read_double(f)
    @TIME end @EPACSE
  end
  for m = 0, MAX_ANGLES_PER_QUAD do
    @ESCAPE for q = 1, 8 do @EMIT
      if m*8 + q - 1 == num_angles then break end
      [angles[q]][m].mu = read_double(f)
    @TIME end @EPACSE
  end
  for m = 0, MAX_ANGLES_PER_QUAD do
    @ESCAPE for q = 1, 8 do @EMIT
      if m*8 + q - 1 == num_angles then break end
      [angles[q]][m].w = read_double(f)
    @TIME end @EPACSE
  end
  -- Check that angles are partitioned correctly into quadrants.
  for m = 0, MAX_ANGLES_PER_QUAD do
    @ESCAPE for q = 1, 8 do @EMIT
      if m*8 + q - 1 == num_angles then break end
      regentlib.assert([angleInQuadrant(q, rexpr [angles[q]][m] end)],
                       'Angle in wrong quadrant')
    @TIME end @EPACSE
  end
  -- Close angles file.
  C.fclose(f)
end

local -- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task color_faces(faces : region(ispace(int3d), face),
                 Nx : int, Ny : int, Nz : int,
                 ntx : int, nty : int, ntz : int,
                 dimension : int, sweepDir : bool[3])
where
  writes(faces.{is_private, tile})
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
  end
end

local -- MANUALLY PARALLELIZED, NO CUDA
task source_term(points : region(ispace(int3d), pointsFSpace),
                 [angles],
                 config : Config,
                 omega : double)
where
  reads(points.[iterIntFields], points.{Ib, sigma}),
  [angles:map(function(a)
     return regentlib.privilege(regentlib.reads, a, 'w')
   end)],
  writes(points.S)
do
  var num_angles = config.Radiation.angles
  __demand(__openmp)
  for p in points do
    var S = (1.0-omega) * p.sigma * p.Ib;
    @ESCAPE for q = 1, 8 do @EMIT
      for m = 0, quadrantSize(q, num_angles) do
        S += omega
           * p.sigma/(4.0*PI)
           * [angles[q]][m].w
           * p.[iterIntFields[q]][m]
      end
    @TIME end @EPACSE
    p.S = S
  end
end

-- 1..6 -> regentlib.task
local function mkBound(wall)

  local emissField = terralib.newlist{
    'xLoEmiss', 'xHiEmiss', 'yLoEmiss', 'yHiEmiss', 'zLoEmiss', 'zHiEmiss'
  }[wall]
  local tempField = terralib.newlist{
    'xLoTemp', 'xHiTemp', 'yLoTemp', 'yHiTemp', 'zLoTemp', 'zHiTemp'
  }[wall]
  local windowField = terralib.newlist{
    'xLoWindow', 'xHiWindow', 'yLoWindow', 'yHiWindow', 'zLoWindow', 'zHiWindow'
  }[wall]
  local incomingQuadrants = terralib.newlist{
    terralib.newlist{5, 6, 7, 8}, -- xi < 0
    terralib.newlist{1, 2, 3, 4}, -- xi > 0
    terralib.newlist{3, 4, 7, 8}, -- eta < 0
    terralib.newlist{1, 2, 5, 6}, -- eta > 0
    terralib.newlist{2, 4, 6, 8}, -- mu < 0
    terralib.newlist{1, 3, 5, 7}, -- mu > 0
  }[wall]
  local outgoingQuadrants = terralib.newlist{
    terralib.newlist{1, 2, 3, 4}, -- xi > 0
    terralib.newlist{5, 6, 7, 8}, -- xi < 0
    terralib.newlist{1, 2, 5, 6}, -- eta > 0
    terralib.newlist{3, 4, 7, 8}, -- eta < 0
    terralib.newlist{1, 3, 5, 7}, -- mu > 0
    terralib.newlist{2, 4, 6, 8}, -- mu < 0
  }[wall]

  local faces = UTIL.generate(8, function()
    return regentlib.newsymbol(region(ispace(int3d), face))
  end)

  local function inCells(dim0, dim1, fromCell, uptoCell)
    return rexpr
      fromCell[0] <= dim0 and dim0 <= uptoCell[0] and
      fromCell[1] <= dim1 and dim1 <= uptoCell[1]
    end
  end

  local -- MANUALLY PARALLELIZED, NO CUDA
  task bound([faces],
             [angles],
             config : Config)
  where
    [faces:map(function(f) return terralib.newlist{
       regentlib.privilege(regentlib.reads, f, 'I'),
       regentlib.privilege(regentlib.writes, f, 'I'),
     } end):flatten()],
    [angles:map(function(a) return terralib.newlist{
       regentlib.privilege(regentlib.reads, a, 'xi'),
       regentlib.privilege(regentlib.reads, a, 'eta'),
       regentlib.privilege(regentlib.reads, a, 'mu'),
       regentlib.privilege(regentlib.reads, a, 'w'),
     } end):flatten()]
  do
    var Nx = config.Radiation.xNum
    var Ny = config.Radiation.yNum
    var Nz = config.Radiation.zNum
    var epsw = config.Radiation.[emissField]
    var Tw = config.Radiation.[tempField]
    var fromCell = config.Radiation.[windowField].fromCell
    var uptoCell = config.Radiation.[windowField].uptoCell
    var num_angles = config.Radiation.angles
    __demand(__openmp)
    for idx in [faces[1]] do
      -- Only update cells on the boundary
      if [terralib.newlist{
            rexpr idx.x == 0  end,
            rexpr idx.x == Nx end,
            rexpr idx.y == 0  end,
            rexpr idx.y == Ny end,
            rexpr idx.z == 0  end,
            rexpr idx.z == Nz end,
          }[wall]] then
        var value = 0.0
        -- Calculate reflected intensity
        if epsw < 1.0 then
          @ESCAPE for _,q in ipairs(incomingQuadrants) do @EMIT
            for m = 0, quadrantSize(q, num_angles) do
              value +=
                (1.0-epsw)/PI * [angles[q]][m].w * [faces[q]][idx].I[m]
                * fabs([terralib.newlist{
                          rexpr [angles[q]][m].xi  end,
                          rexpr [angles[q]][m].xi  end,
                          rexpr [angles[q]][m].eta end,
                          rexpr [angles[q]][m].eta end,
                          rexpr [angles[q]][m].mu  end,
                          rexpr [angles[q]][m].mu  end,
                        }[wall]])
            end
          @TIME end @EPACSE
        end
        -- Add blackbody radiation
        if [terralib.newlist{
              inCells(rexpr idx.y end, rexpr idx.z end, fromCell, uptoCell),
              inCells(rexpr idx.y end, rexpr idx.z end, fromCell, uptoCell),
              inCells(rexpr idx.x end, rexpr idx.z end, fromCell, uptoCell),
              inCells(rexpr idx.x end, rexpr idx.z end, fromCell, uptoCell),
              inCells(rexpr idx.x end, rexpr idx.y end, fromCell, uptoCell),
              inCells(rexpr idx.x end, rexpr idx.y end, fromCell, uptoCell),
            }[wall]] then
          value += epsw*SB*pow(Tw,4.0)/PI
        end
        -- Set outgoing intensity values
        @ESCAPE for _,q in ipairs(outgoingQuadrants) do @EMIT
          for m = 0, quadrantSize(q, num_angles) do
            if [terralib.newlist{
                  rexpr [angles[q]][m].xi  > 0 end,
                  rexpr [angles[q]][m].xi  < 0 end,
                  rexpr [angles[q]][m].eta > 0 end,
                  rexpr [angles[q]][m].eta < 0 end,
                  rexpr [angles[q]][m].mu  > 0 end,
                  rexpr [angles[q]][m].mu  < 0 end,
                }[wall]] then
              [faces[q]][idx].I[m] = value
            end
          end
        @TIME end @EPACSE
      end
    end
  end

  local name = terralib.newlist{
    'bound_x_lo', 'bound_x_hi',
    'bound_y_lo', 'bound_y_hi',
    'bound_z_lo', 'bound_z_hi',
  }[wall]
  bound:set_name(name)
  bound:get_primary_variant():get_ast().name[1] = name -- XXX: Dangerous
  return bound

end -- mkBound

local bound_x_lo = mkBound(1)
local bound_x_hi = mkBound(2)
local bound_y_lo = mkBound(3)
local bound_y_hi = mkBound(4)
local bound_z_lo = mkBound(5)
local bound_z_hi = mkBound(6)

-- 1..8 -> regentlib.task
local function mkSweep(q)

  local fld = intensityFields[q]
  local bnd = regentlib.newsymbol()

  local startx = directions[q][1] and rexpr     bnd.lo.x end or rexpr     bnd.hi.x end
  local dindx  = directions[q][1] and rexpr            1 end or rexpr           -1 end
  local endx   = directions[q][1] and rexpr bnd.hi.x + 1 end or rexpr bnd.lo.x - 1 end

  local starty = directions[q][2] and rexpr     bnd.lo.y end or rexpr     bnd.hi.y end
  local dindy  = directions[q][2] and rexpr            1 end or rexpr           -1 end
  local endy   = directions[q][2] and rexpr bnd.hi.y + 1 end or rexpr bnd.lo.y - 1 end

  local startz = directions[q][3] and rexpr     bnd.lo.z end or rexpr     bnd.hi.z end
  local dindz  = directions[q][3] and rexpr            1 end or rexpr           -1 end
  local endz   = directions[q][3] and rexpr bnd.hi.z + 1 end or rexpr bnd.lo.z - 1 end

  local -- MANUALLY PARALLELIZED, NO OPENMP, NO CUDA
  task sweep(points : region(ispace(int3d), pointsFSpace),
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
             config : Config,
             dx : double, dy : double, dz : double)
  where
    reads(angles.{xi, eta, mu}, points.{S, sigma},
          shared_x_faces_upwind.I, shared_y_faces_upwind.I, shared_z_faces_upwind.I),
    reads writes(points.[fld], x_faces.I, y_faces.I, z_faces.I,
                 shared_x_faces_downwind.I, shared_y_faces_downwind.I, shared_z_faces_downwind.I)
  do
    var num_angles = config.Radiation.angles
    var dAx = dy*dz
    var dAy = dx*dz
    var dAz = dx*dy
    var dV = dx*dy*dz
    var [bnd] = points.bounds
    -- Use our direction and increments for the sweep
    for k = startz,endz,dindz do
      for j = starty,endy,dindy do
        for i = startx,endx,dindx do
          -- Loop over this quadrant's angles.
          for m = 0, quadrantSize(q, num_angles) do
            -- Derive upwind indices
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
            -- Integrate to compute cell-centered value of I
            points[{i,j,k}].[fld][m] = (points[{i,j,k}].S * dV
                                        + fabs(angles[m].xi) * dAx * upwind_x_value/GAMMA
                                        + fabs(angles[m].eta) * dAy * upwind_y_value/GAMMA
                                        + fabs(angles[m].mu) * dAz * upwind_z_value/GAMMA)
                                     / (points[{i,j,k}].sigma * dV
                                        + fabs(angles[m].xi) * dAx/GAMMA
                                        + fabs(angles[m].eta) * dAy/GAMMA
                                        + fabs(angles[m].mu) * dAz/GAMMA)
            -- Compute intensities on downwind faces
            var x_face_val = (points[{i,j,k}].[fld][m] - (1-GAMMA)*upwind_x_value)/GAMMA
            if (x_face_val < 0) then x_face_val = 0 end
            if (indx + dindx) > x_faces.bounds.hi.x or (indx + dindx) < x_faces.bounds.lo.x then
              shared_x_faces_downwind[{indx + dindx, j, k}].I[m] = x_face_val
            else
              x_faces[{indx+dindx, j, k}].I[m] = x_face_val
            end
            var y_face_val = (points[{i,j,k}].[fld][m] - (1-GAMMA)*upwind_y_value)/GAMMA
            if (y_face_val < 0) then y_face_val = 0 end
            if (indy + dindy) > y_faces.bounds.hi.y or (indy + dindy) < y_faces.bounds.lo.y then
              shared_y_faces_downwind[{i, indy + dindy, k}].I[m] = y_face_val
            else
              y_faces[{i, indy+dindy, k}].I[m] = y_face_val
            end
            var z_face_val = (points[{i,j,k}].[fld][m] - (1-GAMMA)*upwind_z_value)/GAMMA
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

  local name = 'sweep_'..tostring(q)
  sweep:set_name(name)
  sweep:get_primary_variant():get_ast().name[1] = name -- XXX: Dangerous
  return sweep

end -- mkSweep

local sweeps = terralib.newlist{
  mkSweep(1),
  mkSweep(2),
  mkSweep(3),
  mkSweep(4),
  mkSweep(5),
  mkSweep(6),
  mkSweep(7),
  mkSweep(8),
}

local -- MANUALLY PARALLELIZED, NO CUDA
task residual(points : region(ispace(int3d), pointsFSpace),
              config : Config,
              Nx : int, Ny : int, Nz : int)
where
  reads(points.[intensityFields], points.[iterIntFields])
do
  var num_angles = config.Radiation.angles
  var res : double = 0.0;
  @ESCAPE for q = 1, 8 do @EMIT
    __demand(__openmp)
    for p in points do
      for m = 0, quadrantSize(q, num_angles) do
        var v1 = p.[intensityFields[q]][m]
        if v1 > 0 then
          var v2 = v1 - p.[iterIntFields[q]][m]
          res += (v2 * v2) / (v1 * v1)
        end
      end
    end
  @TIME end @EPACSE
  return res
end

local -- MANUALLY PARALLELIZED, NO CUDA
task update(points : region(ispace(int3d), pointsFSpace),
            config : Config)
where
  reads(points.[intensityFields]),
  reads writes(points.[iterIntFields])
do
  var num_angles = config.Radiation.angles;
  @ESCAPE for q = 1, 8 do @EMIT
    -- Request that outer loop is not vectorized, to ensure the inner one gets
    -- vectorized instead.
    __forbid(__vectorize) __demand(__openmp)
    for p in points do
      for m = 0, quadrantSize(q, num_angles) do
        p.[iterIntFields[q]][m] = p.[intensityFields[q]][m]
      end
    end
  @TIME end @EPACSE
end

local -- MANUALLY PARALLELIZED, NO CUDA
task reduce_intensity(points : region(ispace(int3d), pointsFSpace),
                      [angles],
                      config : Config)
where
  reads(points.[intensityFields]),
  [angles:map(function(a) return
     regentlib.privilege(regentlib.reads, a, 'w')
   end)],
  reads writes(points.G)
do
  var num_angles = config.Radiation.angles
  __demand(__openmp)
  for p in points do
    @ESCAPE for q = 1, 8 do @EMIT
      for m = 0, quadrantSize(q, num_angles) do
        p.G += [angles[q]][m].w * p.[intensityFields[q]][m]
      end
    @TIME end @EPACSE
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

  local x_faces = UTIL.generate(8, regentlib.newsymbol)
  local y_faces = UTIL.generate(8, regentlib.newsymbol)
  local z_faces = UTIL.generate(8, regentlib.newsymbol)
  local angles = UTIL.generate(8, regentlib.newsymbol)

  local tiles_private = regentlib.newsymbol('tiles_private')
  local x_tiles_shared = regentlib.newsymbol('x_tiles_shared')
  local y_tiles_shared = regentlib.newsymbol('y_tiles_shared')
  local z_tiles_shared = regentlib.newsymbol('z_tiles_shared')

  local s_x_faces = UTIL.generate(8, regentlib.newsymbol)
  local s_y_faces = UTIL.generate(8, regentlib.newsymbol)
  local s_z_faces = UTIL.generate(8, regentlib.newsymbol)

  local p_x_faces = UTIL.generate(8, regentlib.newsymbol)
  local p_y_faces = UTIL.generate(8, regentlib.newsymbol)
  local p_z_faces = UTIL.generate(8, regentlib.newsymbol)

  function INSTANCE.DeclSymbols(config) return rquote

    -- Number of points in each dimension
    var [Nx] = config.Radiation.xNum
    var [Ny] = config.Radiation.yNum
    var [Nz] = config.Radiation.zNum
    -- Number of tiles in each dimension
    var [ntx] = config.Mapping.tiles[0]
    var [nty] = config.Mapping.tiles[1]
    var [ntz] = config.Mapping.tiles[2]

    -- Regions for faces (+1 in one direction since one more face than points)
    var grid_x = ispace(int3d, {Nx+1, Ny,   Nz  })
    var grid_y = ispace(int3d, {Nx,   Ny+1, Nz  })
    var grid_z = ispace(int3d, {Nx,   Ny,   Nz+1});
    @ESCAPE for q = 1, 8 do @EMIT
      var [x_faces[q]] = region(grid_x, face)
      var [y_faces[q]] = region(grid_y, face)
      var [z_faces[q]] = region(grid_z, face)
    @TIME end @EPACSE

    -- Regions for angle values
    var angle_indices = ispace(int1d, MAX_ANGLES_PER_QUAD);
    @ESCAPE for q = 1, 8 do @EMIT
      var [angles[q]] = region(angle_indices, angle)
    @TIME end @EPACSE

    -- Partition faces
    var [tiles_private]  = ispace(int3d, {ntx,   nty,   ntz  })
    -- extra tile required for shared edge
    var [x_tiles_shared] = ispace(int3d, {ntx+1, nty,   ntz  })
    var [y_tiles_shared] = ispace(int3d, {ntx,   nty+1, ntz  })
    var [z_tiles_shared] = ispace(int3d, {ntx,   nty,   ntz+1});
    @ESCAPE for q = 1, 8 do @EMIT
      -- x
      var p_x_faces_equal = partition(equal, [x_faces[q]], [tiles_private])
      for c in [tiles_private] do
        color_faces(p_x_faces_equal[c], Nx, Ny, Nz, ntx, nty, ntz, 0, array([directions[q]]))
      end
      var x_by_privacy = partition([x_faces[q]].is_private, ispace(int1d,2))
      var p_x = x_by_privacy[1]
      var [p_x_faces[q]] = partition(p_x.tile, [tiles_private])
      var s_x = x_by_privacy[0]
      var [s_x_faces[q]] = partition(s_x.tile, [x_tiles_shared])
      -- y
      var p_y_faces_equal = partition(equal, [y_faces[q]], [tiles_private])
      for c in [tiles_private] do
        color_faces(p_y_faces_equal[c], Nx, Ny, Nz, ntx, nty, ntz, 1, array([directions[q]]))
      end
      var y_by_privacy = partition([y_faces[q]].is_private, ispace(int1d,2))
      var p_y = y_by_privacy[1]
      var [p_y_faces[q]] = partition(p_y.tile, [tiles_private])
      var s_y = y_by_privacy[0]
      var [s_y_faces[q]] = partition(s_y.tile, [y_tiles_shared])
      -- z
      var p_z_faces_equal = partition(equal, [z_faces[q]], [tiles_private])
      for c in [tiles_private] do
        color_faces(p_z_faces_equal[c], Nx, Ny, Nz, ntx, nty, ntz, 2, array([directions[q]]))
      end
      var z_by_privacy = partition([z_faces[q]].is_private, ispace(int1d,2))
      var p_z = z_by_privacy[1]
      var [p_z_faces[q]] = partition(p_z.tile, [tiles_private])
      var s_z = z_by_privacy[0]
      var [s_z_faces[q]] = partition(s_z.tile, [z_tiles_shared])
    @TIME end @EPACSE

  end end -- DeclSymbols

  function INSTANCE.InitRegions(config) return rquote

    -- Initialize angle values
    initialize_angles([angles], config)

  end end -- InitRegions

  function INSTANCE.ComputeRadiationField(config, tiles, p_points) return rquote

    var dx = config.Grid.xWidth / config.Radiation.xNum
    var dy = config.Grid.yWidth / config.Radiation.yNum
    var dz = config.Grid.zWidth / config.Radiation.zNum

    var omega = config.Radiation.qs/(config.Radiation.qa+config.Radiation.qs)

    -- Compute until convergence
    var res : double = 1.0
    while (res > TOLERANCE) do

      -- Update the source term (in this problem, isotropic)
      for t in tiles do
        source_term(p_points[t], [angles], config, omega)
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
                     [angles],
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
                     [angles],
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
                     [angles],
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
                     [angles],
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
                     [angles],
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
                     [angles],
                     config)
        end
      end

      --Perform the sweep for computing new intensities
      --Quadrant 1 - +x, +y, +z
      for i = 0, ntx do
        for j = 0, nty do
          for k = 0, ntz do
            [sweeps[1]](p_points[{i,j,k}],
                        [p_x_faces[1]][{i,j,k}], [p_y_faces[1]][{i,j,k}], [p_z_faces[1]][{i,j,k}],
                        [s_x_faces[1]][{i,j,k}], [s_x_faces[1]][{i+1,j,k}],
                        [s_y_faces[1]][{i,j,k}], [s_y_faces[1]][{i,j+1,k}],
                        [s_z_faces[1]][{i,j,k}], [s_z_faces[1]][{i,j,k+1}],
                        [angles[1]], config, dx, dy, dz)
          end
        end
      end

      -- Quadrant 2 - +x, +y, -z
      for i = 0, ntx do
        for j = 0, nty do
          for k = ntz-1, -1, -1 do
            [sweeps[2]](p_points[{i,j,k}],
                        [p_x_faces[2]][{i,j,k}], [p_y_faces[2]][{i,j,k}], [p_z_faces[2]][{i,j,k}],
                        [s_x_faces[2]][{i,j,k}], [s_x_faces[2]][{i+1,j,k}],
                        [s_y_faces[2]][{i,j,k}], [s_y_faces[2]][{i,j+1,k}],
                        [s_z_faces[2]][{i,j,k+1}], [s_z_faces[2]][{i,j,k}],
                        [angles[2]], config, dx, dy, dz)
          end
        end
      end

      -- Quadrant 3 - +x, -y, +z
      for i = 0, ntx do
        for j = nty-1, -1, -1 do
          for k = 0, ntz do
            [sweeps[3]](p_points[{i,j,k}],
                        [p_x_faces[3]][{i,j,k}], [p_y_faces[3]][{i,j,k}], [p_z_faces[3]][{i,j,k}],
                        [s_x_faces[3]][{i,j,k}], [s_x_faces[3]][{i+1,j,k}],
                        [s_y_faces[3]][{i,j+1,k}], [s_y_faces[3]][{i,j,k}],
                        [s_z_faces[3]][{i,j,k}], [s_z_faces[3]][{i,j,k+1}],
                        [angles[3]], config, dx, dy, dz)
          end
        end
      end

      -- Quadrant 4 - +x, -y, -z
      for i = 0, ntx do
        for j = nty-1, -1, -1 do
          for k = ntz-1, -1, -1 do
            [sweeps[4]](p_points[{i,j,k}],
                        [p_x_faces[4]][{i,j,k}], [p_y_faces[4]][{i,j,k}], [p_z_faces[4]][{i,j,k}],
                        [s_x_faces[4]][{i,j,k}], [s_x_faces[4]][{i+1,j,k}],
                        [s_y_faces[4]][{i,j+1,k}], [s_y_faces[4]][{i,j,k}],
                        [s_z_faces[4]][{i,j,k+1}], [s_z_faces[4]][{i,j,k}],
                        [angles[4]], config, dx, dy, dz)
          end
        end
      end

      -- Quadrant 5 - -x, +y, +z
      for i = ntx-1, -1, -1 do
        for j = 0, nty do
          for k = 0, ntz do
            [sweeps[5]](p_points[{i,j,k}],
                        [p_x_faces[5]][{i,j,k}], [p_y_faces[5]][{i,j,k}], [p_z_faces[5]][{i,j,k}],
                        [s_x_faces[5]][{i+1,j,k}], [s_x_faces[5]][{i,j,k}],
                        [s_y_faces[5]][{i,j,k}], [s_y_faces[5]][{i,j+1,k}],
                        [s_z_faces[5]][{i,j,k}], [s_z_faces[5]][{i,j,k+1}],
                        [angles[5]], config, dx, dy, dz)
          end
        end
      end

      -- Quadrant 6 - -x, +y, -z
      for i = ntx-1, -1, -1 do
        for j = 0, nty do
          for k = ntz-1, -1, -1 do
            [sweeps[6]](p_points[{i,j,k}],
                        [p_x_faces[6]][{i,j,k}], [p_y_faces[6]][{i,j,k}], [p_z_faces[6]][{i,j,k}],
                        [s_x_faces[6]][{i+1,j,k}], [s_x_faces[6]][{i,j,k}],
                        [s_y_faces[6]][{i,j,k}], [s_y_faces[6]][{i,j+1,k}],
                        [s_z_faces[6]][{i,j,k+1}], [s_z_faces[6]][{i,j,k}],
                        [angles[6]], config, dx, dy, dz)
          end
        end
      end

      -- Quadrant 7 - -x, -y, +z
      for i = ntx-1, -1, -1 do
        for j = nty-1, -1, -1 do
          for k = 0, ntz do
            [sweeps[7]](p_points[{i,j,k}],
                        [p_x_faces[7]][{i,j,k}], [p_y_faces[7]][{i,j,k}], [p_z_faces[7]][{i,j,k}],
                        [s_x_faces[7]][{i+1,j,k}], [s_x_faces[7]][{i,j,k}],
                        [s_y_faces[7]][{i,j+1,k}], [s_y_faces[7]][{i,j,k}],
                        [s_z_faces[7]][{i,j,k}], [s_z_faces[7]][{i,j,k+1}],
                        [angles[7]], config, dx, dy, dz)
          end
        end
      end

      -- Quadrant 8 - -x, -y, -z
      for i = ntx-1, -1, -1 do
        for j = nty-1, -1, -1 do
          for k = ntz-1, -1, -1 do
            [sweeps[8]](p_points[{i,j,k}],
                        [p_x_faces[8]][{i,j,k}], [p_y_faces[8]][{i,j,k}], [p_z_faces[8]][{i,j,k}],
                        [s_x_faces[8]][{i+1,j,k}], [s_x_faces[8]][{i,j,k}],
                        [s_y_faces[8]][{i,j+1,k}], [s_y_faces[8]][{i,j,k}],
                        [s_z_faces[8]][{i,j,k+1}], [s_z_faces[8]][{i,j,k}],
                        [angles[8]], config, dx, dy, dz)
          end
        end
      end

      -- Compute the residual
      res = 0.0
      for t in tiles do
        res += residual(p_points[t], config, Nx, Ny, Nz)
      end
      res = sqrt(res/(Nx*Ny*Nz*config.Radiation.angles))

      -- Update the intensities
      for t in tiles do
        update(p_points[t], config)
      end

    end

    -- Reduce intensity
    for t in tiles do
      reduce_intensity(p_points[t], [angles], config)
    end

  end end -- ComputeRadiationField

return INSTANCE end -- mkInstance

-------------------------------------------------------------------------------
-- MODULE END
-------------------------------------------------------------------------------

return MODULE end
