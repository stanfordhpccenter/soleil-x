import 'regent'

-------------------------------------------------------------------------------
-- MODULE PARAMETERS
-------------------------------------------------------------------------------

return function(MAX_ANGLES_PER_QUAD, Point_columns, Config) local MODULE = {}

-------------------------------------------------------------------------------
-- IMPORTS
-------------------------------------------------------------------------------

local C = regentlib.c
local MAPPER = terralib.includec("soleil_mapper.h")
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

local struct Angle_columns {
  xi  : double;
  eta : double;
  mu  : double;
  w   : double;
}

local struct Face_columns {
  I      : regentlib.array(double, MAX_ANGLES_PER_QUAD);
  I_prev : regentlib.array(double, MAX_ANGLES_PER_QUAD);
}

-- A sub-point holds information specific to a cell center and angle.
local struct SubPoint_columns {
  I : double;
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
  return regentlib.newsymbol(region(ispace(int1d), Angle_columns))
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

local __demand(__cuda) -- MANUALLY PARALLELIZED
task initialize_points(points : region(ispace(int3d), Point_columns))
where
  writes(points.{G, S})
do
  __demand(__openmp)
  for p in points do
    p.G = 0.0
    p.S = 0.0
  end
end

-- Sub-points within a tile are laid out in the order that the sweep code will
-- process them: (x,y,z) point coordinates are grouped by diagonal, and angle
-- values are contiguous per point. E.g. on a 2x2 grid with 2 angles, subpoints
-- would be laid out as follows (order of elements is mxy):
-- 000 100 200 010 110 210 001 101 201 011 111 211
-- |<diag. 0>| |<    diagonal 1     >| |<diag. 2>|
-- This task fills in the mappings that allow us to move from one ordering to
-- the other. Given a 1d subpoint index s, we would proceed as follows to find
-- the 3d point p it corresponds to:
-- * Split the 1d index point s into its 4 coordinates m,x,y,z.
-- * Follow the s3d_to_p_Q field for the appropriate quadrant Q on (x,y,z).

local p_to_s3d = terralib.newlist{
  'p_to_s3d_1', 'p_to_s3d_2', 'p_to_s3d_3', 'p_to_s3d_4',
  'p_to_s3d_5', 'p_to_s3d_6', 'p_to_s3d_7', 'p_to_s3d_8'
}
local s3d_to_p = terralib.newlist{
  's3d_to_p_1', 's3d_to_p_2', 's3d_to_p_3', 's3d_to_p_4',
  's3d_to_p_5', 's3d_to_p_6', 's3d_to_p_7', 's3d_to_p_8'
}

local -- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task cache_grid_translation(points : region(ispace(int3d), Point_columns))
where
  writes(points.[p_to_s3d], points.[s3d_to_p])
do
  var Tx = points.bounds.hi.x - points.bounds.lo.x + 1
  var Ty = points.bounds.hi.y - points.bounds.lo.y + 1
  var Tz = points.bounds.hi.z - points.bounds.lo.z + 1;
  @ESCAPE for q = 1, 8 do @EMIT
    -- Start one-before the first grid-order index
    var grid = int3d{-1,0,0}
    for d = 0, (Tx-1)+(Ty-1)+(Tz-1)+1 do
      -- Set diagonal-order index to the smallest for this diagonal
      var diag : int3d
      diag.x = min(d, Tx-1)
      diag.y = min(d-diag.x, Ty-1)
      diag.z = d-diag.x-diag.y
      while true do
        -- Advance grid-order index
        if grid.x < Tx-1 then
          grid.x += 1
        elseif grid.y < Ty-1 then
          grid.x = 0
          grid.y += 1
        elseif grid.z < Tz-1 then
          grid.x = 0
          grid.y = 0
          grid.z += 1
        else regentlib.assert(false, 'Internal error') end
        -- Store mapping for this pair of indices
        var real = int3d{
          [directions[q][1] and rexpr diag.x end or rexpr Tx-diag.x-1 end],
          [directions[q][2] and rexpr diag.y end or rexpr Ty-diag.y-1 end],
          [directions[q][3] and rexpr diag.z end or rexpr Tz-diag.z-1 end]}
        points[real + points.bounds.lo].[p_to_s3d[q]] = grid + points.bounds.lo
        points[grid + points.bounds.lo].[s3d_to_p[q]] = real + points.bounds.lo
        -- Advance diagonal-order index
        if diag.x > 0 and diag.y < Ty-1 then
          diag.x -= 1
          diag.y += 1
        elseif diag.z < min(d, Tz-1) then
          diag.z += 1
          diag.x = min(d-diag.z, Tx-1)
          diag.y = d-diag.z-diag.x
        else
          -- We've run out of indices on this diagonal, continue to next one
          break
        end
      end
    end
    regentlib.assert(grid.x == Tx-1 and grid.y == Ty-1 and grid.z == Tz-1,
                     'Internal error')
  @TIME end @EPACSE
end

local __demand(__cuda) -- MANUALLY PARALLELIZED
task initialize_sub_points(sub_points : region(ispace(int1d), SubPoint_columns))
where
  writes(sub_points.I)
do
  __demand(__openmp)
  for s1d in sub_points do
    s1d.I = 0.0
  end
end

-- 'x'|'y'|'z', 1..8 -> regentlib.task
local function mkInitializeFaces(dim, q)

  local __demand(__cuda) -- MANUALLY PARALLELIZED
  task initialize_faces(faces : region(ispace(int2d), Face_columns),
                        config : Config)
  where
    reads writes(faces.I)
  do
    var num_angles = config.Radiation.angles
    for m = 0, quadrantSize(q, num_angles) do
      __demand(__openmp)
      for f in faces do
        f.I[m] = 0.0
      end
    end
  end

  local name = 'initialize_faces_'..dim..'_'..tostring(q)
  initialize_faces:set_name(name)
  initialize_faces:get_primary_variant():get_ast().name[1] = name
  return initialize_faces

end -- mkInitializeFaces

local initialize_faces = {
  x = UTIL.range(1,8):map(function(q) return mkInitializeFaces('x', q) end),
  y = UTIL.range(1,8):map(function(q) return mkInitializeFaces('y', q) end),
  z = UTIL.range(1,8):map(function(q) return mkInitializeFaces('z', q) end),
}

local __demand(__cuda) -- MANUALLY PARALLELIZED
task source_term(points : region(ispace(int3d), Point_columns),
                 config : Config)
where
  reads(points.{Ib, sigma, G}),
  writes(points.S)
do
  var omega = config.Radiation.qs/(config.Radiation.qa+config.Radiation.qs)
  __demand(__openmp)
  for p in points do
    p.S = (1.0-omega) * p.sigma * p.Ib + omega * p.sigma/(4.0*PI) * p.G
  end
end

-- 'x'|'y'|'z', 1..8 -> regentlib.task
local function mkCacheIntensity(dim, q)

  local __demand(__cuda) -- MANUALLY PARALLELIZED
  task cache_intensity(faces : region(ispace(int2d), Face_columns),
                       config : Config)
  where
    reads(faces.I),
    reads writes(faces.I_prev)
  do
    var num_angles = config.Radiation.angles
    for m = 0, quadrantSize(q, num_angles) do
      __demand(__openmp)
      for f in faces do
        f.I_prev[m] = f.I[m]
      end
    end
  end

  local name = 'cache_intensity_'..dim..'_'..tostring(q)
  cache_intensity:set_name(name)
  cache_intensity:get_primary_variant():get_ast().name[1] = name
  return cache_intensity

end -- mkCacheIntensity

local cache_intensity = {
  x = UTIL.range(1,8):map(function(q) return mkCacheIntensity('x', q) end),
  y = UTIL.range(1,8):map(function(q) return mkCacheIntensity('y', q) end),
  z = UTIL.range(1,8):map(function(q) return mkCacheIntensity('z', q) end),
}

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
    return regentlib.newsymbol(region(ispace(int2d), Face_columns))
  end)

  local __demand(__cuda) -- MANUALLY PARALLELIZED
  task bound([faces],
             [angles],
             config : Config)
  where
    [incomingQuadrants:map(function(q)
       return regentlib.privilege(regentlib.reads, faces[q], 'I_prev')
     end)],
    [outgoingQuadrants:map(function(q) return terralib.newlist{
       regentlib.privilege(regentlib.reads, faces[q], 'I'),
       regentlib.privilege(regentlib.writes, faces[q], 'I')
     } end):flatten()],
    [angles:map(function(a) return terralib.newlist{
       regentlib.privilege(regentlib.reads, a, 'xi'),
       regentlib.privilege(regentlib.reads, a, 'eta'),
       regentlib.privilege(regentlib.reads, a, 'mu'),
       regentlib.privilege(regentlib.reads, a, 'w'),
     } end):flatten()]
  do
    var epsw = config.Radiation.[emissField]
    var Tw = config.Radiation.[tempField]
    var fromCell = config.Radiation.[windowField].fromCell
    var uptoCell = config.Radiation.[windowField].uptoCell
    var num_angles = config.Radiation.angles
    __demand(__openmp)
    for idx in [faces[1]] do
      var a = idx.x
      var b = idx.y
      var value = 0.0
      -- Calculate reflected intensity
      if epsw < 1.0 then
        @ESCAPE for _,q in ipairs(incomingQuadrants) do @EMIT
          for m = 0, quadrantSize(q, num_angles) do
            value +=
              (1.0-epsw)/PI * [angles[q]][m].w * [faces[q]][idx].I_prev[m]
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
      if fromCell[0] <= a and a <= uptoCell[0] and
         fromCell[1] <= b and b <= uptoCell[1] then
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

  local __demand(__cuda) -- MANUALLY PARALLELIZED
  task sweep(points : region(ispace(int3d), Point_columns),
             sub_points : region(ispace(int1d), SubPoint_columns),
             x_faces : region(ispace(int2d), Face_columns),
             y_faces : region(ispace(int2d), Face_columns),
             z_faces : region(ispace(int2d), Face_columns),
             angles : region(ispace(int1d), Angle_columns),
             config : Config)
  where
    reads(angles.{xi, eta, mu}, points.{S, sigma, [p_to_s3d[q]]}),
    reads writes(sub_points.I, x_faces.I, y_faces.I, z_faces.I)
  do
    var Tx = points.bounds.hi.x - points.bounds.lo.x + 1
    var Ty = points.bounds.hi.y - points.bounds.lo.y + 1
    var Tz = points.bounds.hi.z - points.bounds.lo.z + 1
    var dx = config.Grid.xWidth / config.Radiation.xNum
    var dy = config.Grid.yWidth / config.Radiation.yNum
    var dz = config.Grid.zWidth / config.Radiation.zNum
    var dAx = dy*dz
    var dAy = dx*dz
    var dAz = dx*dy
    var dV = dx*dy*dz
    var num_angles = config.Radiation.angles
    var [bnd] = points.bounds
    var res = 0.0
    for z = startz,endz,dindz do
      for y = starty,endy,dindy do
        for x = startx,endx,dindx do
          var s3d = points[{x,y,z}].[p_to_s3d[q]]
          var s1d = MAX_ANGLES_PER_QUAD * s3d.x
                  + MAX_ANGLES_PER_QUAD * Tx    * s3d.y
                  + MAX_ANGLES_PER_QUAD * Tx    * Ty    * s3d.z
          for m = 0, quadrantSize(q, num_angles) do
            -- Read upwind face values
            var x_value = x_faces[{  y,z}].I[m]
            var y_value = y_faces[{x,  z}].I[m]
            var z_value = z_faces[{x,y  }].I[m]
            -- Integrate to compute cell-centered value of I
            var oldI = sub_points[s1d + m].I
            var newI = (points[{x,y,z}].S * dV
                        + fabs(angles[m].xi)  * dAx * x_value/GAMMA
                        + fabs(angles[m].eta) * dAy * y_value/GAMMA
                        + fabs(angles[m].mu)  * dAz * z_value/GAMMA)
                     / (points[{x,y,z}].sigma * dV
                        + fabs(angles[m].xi)  * dAx/GAMMA
                        + fabs(angles[m].eta) * dAy/GAMMA
                        + fabs(angles[m].mu)  * dAz/GAMMA)
            if newI > 0.0 then
              res += pow(newI-oldI,2) / pow(newI,2)
            end
            sub_points[s1d + m].I = newI
            -- Compute intensities on downwind faces
            x_faces[{  y,z}].I[m] = max(0.0, (newI-(1-GAMMA)*x_value)/GAMMA)
            y_faces[{x,  z}].I[m] = max(0.0, (newI-(1-GAMMA)*y_value)/GAMMA)
            z_faces[{x,y  }].I[m] = max(0.0, (newI-(1-GAMMA)*z_value)/GAMMA)
          end
        end
      end
    end
    return res
  end

  local name = 'sweep_'..tostring(q)
  sweep:set_name(name)
  sweep:get_primary_variant():get_ast().name[1] = name -- XXX: Dangerous
  return sweep

end -- mkSweep

local sweep = UTIL.range(1,8):map(function(q) return mkSweep(q) end)

local sub_points = UTIL.generate(8, function()
  return regentlib.newsymbol(region(ispace(int1d), SubPoint_columns))
end)

local __demand(__cuda) -- MANUALLY PARALLELIZED
task reduce_intensity(points : region(ispace(int3d), Point_columns),
                      [sub_points],
                      [angles],
                      config : Config)
where
  reads(points.[p_to_s3d]),
  [sub_points:map(function(s)
     return regentlib.privilege(regentlib.reads, s, 'I')
   end)],
  [angles:map(function(a)
     return regentlib.privilege(regentlib.reads, a, 'w')
   end)],
  reads writes(points.G)
do
  var Tx = points.bounds.hi.x - points.bounds.lo.x + 1
  var Ty = points.bounds.hi.y - points.bounds.lo.y + 1
  var Tz = points.bounds.hi.z - points.bounds.lo.z + 1
  var num_angles = config.Radiation.angles
  __demand(__openmp)
  for p in points do
    p.G = 0.0
  end
  @ESCAPE for q = 1, 8 do @EMIT
    __demand(__openmp)
    for p in points do
      var G = 0.0
      var s3d = p.[p_to_s3d[q]]
      var s1d = MAX_ANGLES_PER_QUAD * s3d.x
              + MAX_ANGLES_PER_QUAD * Tx    * s3d.y
              + MAX_ANGLES_PER_QUAD * Tx    * Ty    * s3d.z
      for m = 0, quadrantSize(q, num_angles) do
        G += [angles[q]][m].w * [sub_points[q]][s1d + m].I
      end
      p.G += G
    end
  @TIME end @EPACSE
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
  local Tx = regentlib.newsymbol('Tx')
  local Ty = regentlib.newsymbol('Ty')
  local Tz = regentlib.newsymbol('Tz')

  local sub_points = UTIL.generate(8, regentlib.newsymbol)
  local x_faces = UTIL.generate(8, regentlib.newsymbol)
  local y_faces = UTIL.generate(8, regentlib.newsymbol)
  local z_faces = UTIL.generate(8, regentlib.newsymbol)
  local angles = UTIL.generate(8, regentlib.newsymbol)

  local x_tiles = regentlib.newsymbol('x_tiles')
  local y_tiles = regentlib.newsymbol('y_tiles')
  local z_tiles = regentlib.newsymbol('z_tiles')

  local p_sub_points = UTIL.generate(8, regentlib.newsymbol)
  local p_x_faces = UTIL.generate(8, regentlib.newsymbol)
  local p_y_faces = UTIL.generate(8, regentlib.newsymbol)
  local p_z_faces = UTIL.generate(8, regentlib.newsymbol)

  function INSTANCE.DeclSymbols(config, tiles) return rquote

    var sampleId = config.Mapping.sampleId

    -- Number of points in each dimension
    var [Nx] = config.Radiation.xNum
    var [Ny] = config.Radiation.yNum
    var [Nz] = config.Radiation.zNum
    -- Number of tiles in each dimension
    var [ntx] = config.Mapping.tiles[0]
    var [nty] = config.Mapping.tiles[1]
    var [ntz] = config.Mapping.tiles[2]
    -- Sanity-check partitioning
    regentlib.assert(Nx % ntx == 0, "Uneven partitioning of radiation grid on x")
    regentlib.assert(Ny % nty == 0, "Uneven partitioning of radiation grid on y")
    regentlib.assert(Nz % ntz == 0, "Uneven partitioning of radiation grid on z")
    -- Number of points in each tile
    var [Tx] = Nx / ntx
    var [Ty] = Ny / nty
    var [Tz] = Nz / ntz

    -- Region for points
    -- (managed by the host code)

    -- Regions for sub-points
    -- Conceptually int4d, but rolled into 1 dimension to make CUDA code
    -- generation easier. The effective storage order is Z > Y > X > M.
    var is_sub_points = ispace(int1d, MAX_ANGLES_PER_QUAD*Nx*Ny*Nz);
    @ESCAPE for q = 1, 8 do @EMIT
      var [sub_points[q]] = region(is_sub_points, SubPoint_columns);
      [UTIL.emitRegionTagAttach(sub_points[q], MAPPER.SAMPLE_ID_TAG, sampleId, int)];
    @TIME end @EPACSE

    -- Regions for faces
    var grid_x = ispace(int2d, {   Ny,Nz})
    var grid_y = ispace(int2d, {Nx,   Nz})
    var grid_z = ispace(int2d, {Nx,Ny   });
    @ESCAPE for q = 1, 8 do @EMIT
      var [x_faces[q]] = region(grid_x, Face_columns);
      [UTIL.emitRegionTagAttach(x_faces[q], MAPPER.SAMPLE_ID_TAG, sampleId, int)];
      var [y_faces[q]] = region(grid_y, Face_columns);
      [UTIL.emitRegionTagAttach(y_faces[q], MAPPER.SAMPLE_ID_TAG, sampleId, int)];
      var [z_faces[q]] = region(grid_z, Face_columns);
      [UTIL.emitRegionTagAttach(z_faces[q], MAPPER.SAMPLE_ID_TAG, sampleId, int)];
    @TIME end @EPACSE

    -- Regions for angles
    @ESCAPE for q = 1, 8 do @EMIT
      var is = ispace(int1d, quadrantSize(q, config.Radiation.angles))
      var [angles[q]] = region(is, Angle_columns);
      [UTIL.emitRegionTagAttach(angles[q], MAPPER.SAMPLE_ID_TAG, sampleId, int)];
    @TIME end @EPACSE

    -- Partition points
    -- (done by the host code)

    -- Partition sub-points
    @ESCAPE for q = 1, 8 do @EMIT
      var [p_sub_points[q]] =
        [UTIL.mkPartitionEqually(int1d, int3d, SubPoint_columns)]
        ([sub_points[q]], tiles)
    @TIME end @EPACSE

    -- Partition faces
    var [x_tiles] = ispace(int2d, {    nty,ntz})
    var [y_tiles] = ispace(int2d, {ntx,    ntz})
    var [z_tiles] = ispace(int2d, {ntx,nty    });
    @ESCAPE for q = 1, 8 do @EMIT
      -- x
      var x_coloring = regentlib.c.legion_domain_point_coloring_create()
      for c in x_tiles do
        var a = c.x
        var b = c.y
        var rect = rect2d{lo = int2d{Ty*a,       Tz*b      },
                          hi = int2d{Ty*(a+1)-1, Tz*(b+1)-1}}
        regentlib.c.legion_domain_point_coloring_color_domain(x_coloring, c, rect)
      end
      var [p_x_faces[q]] = partition(disjoint, [x_faces[q]], x_coloring, x_tiles)
      regentlib.c.legion_domain_point_coloring_destroy(x_coloring);
      -- y
      var y_coloring = regentlib.c.legion_domain_point_coloring_create()
      for c in y_tiles do
        var a = c.x
        var b = c.y
        var rect = rect2d{lo = int2d{Tx*a,       Tz*b      },
                          hi = int2d{Tx*(a+1)-1, Tz*(b+1)-1}}
        regentlib.c.legion_domain_point_coloring_color_domain(y_coloring, c, rect)
      end
      var [p_y_faces[q]] = partition(disjoint, [y_faces[q]], y_coloring, y_tiles)
      regentlib.c.legion_domain_point_coloring_destroy(y_coloring);
      -- z
      var z_coloring = regentlib.c.legion_domain_point_coloring_create()
      for c in z_tiles do
        var a = c.x
        var b = c.y
        var rect = rect2d{lo = int2d{Tx*a,       Ty*b      },
                          hi = int2d{Tx*(a+1)-1, Ty*(b+1)-1}}
        regentlib.c.legion_domain_point_coloring_color_domain(z_coloring, c, rect)
      end
      var [p_z_faces[q]] = partition(disjoint, [z_faces[q]], z_coloring, z_tiles)
      regentlib.c.legion_domain_point_coloring_destroy(z_coloring);
    @TIME end @EPACSE

  end end -- DeclSymbols

  function INSTANCE.InitRegions(config, tiles, p_points) return rquote

    -- Initialize points
    for c in tiles do
      initialize_points(p_points[c])
    end
    for c in tiles do
      cache_grid_translation(p_points[c])
    end

    -- Initialize sub-points
    @ESCAPE for q = 1, 8 do @EMIT
      for c in tiles do
        initialize_sub_points([p_sub_points[q]][c])
      end
    @TIME end @EPACSE

    -- Initialize faces
    @ESCAPE for q = 1, 8 do @EMIT
      for c in x_tiles do
        [initialize_faces['x'][q]]([p_x_faces[q]][c], config)
      end
      for c in y_tiles do
        [initialize_faces['y'][q]]([p_y_faces[q]][c], config)
      end
      for c in z_tiles do
        [initialize_faces['z'][q]]([p_z_faces[q]][c], config)
      end
    @TIME end @EPACSE

    -- Initialize angles
    initialize_angles([angles], config)

  end end -- InitRegions

  function INSTANCE.ComputeRadiationField(config, tiles, p_points) return rquote

    -- Initialize intensity.
    for c in tiles do
      reduce_intensity(p_points[c],
                       [p_sub_points:map(function(s) return rexpr s[c] end end)],
                       [angles],
                       config)
    end

    -- Compute until convergence.
    var res : double = 1.0
    while res > TOLERANCE do

      -- Update the source term.
      for c in tiles do
        source_term(p_points[c], config)
      end

      -- Cache the face intensity values from the previous iteration (those
      -- values represent the final downwind values).
      @ESCAPE for q = 1, 8 do @EMIT
        for c in x_tiles do
          [cache_intensity['x'][q]]([p_x_faces[q]][c], config)
        end
        for c in y_tiles do
          [cache_intensity['y'][q]]([p_y_faces[q]][c], config)
        end
        for c in z_tiles do
          [cache_intensity['z'][q]]([p_z_faces[q]][c], config)
        end
      @TIME end @EPACSE

      -- Update face intensity values, to represent initial upwind values for
      -- this iteration.
      for c in x_tiles do
        bound_x_lo([p_x_faces:map(function(f) return rexpr f[c] end end)],
                   [angles],
                   config)
      end
      for c in x_tiles do
        bound_x_hi([p_x_faces:map(function(f) return rexpr f[c] end end)],
                   [angles],
                   config)
      end
      for c in y_tiles do
        bound_y_lo([p_y_faces:map(function(f) return rexpr f[c] end end)],
                   [angles],
                   config)
      end
      for c in y_tiles do
        bound_y_hi([p_y_faces:map(function(f) return rexpr f[c] end end)],
                   [angles],
                   config)
      end
      for c in z_tiles do
        bound_z_lo([p_z_faces:map(function(f) return rexpr f[c] end end)],
                   [angles],
                   config)
      end
      for c in z_tiles do
        bound_z_hi([p_z_faces:map(function(f) return rexpr f[c] end end)],
                   [angles],
                   config)
      end

      -- Perform the sweep for computing new intensities.
      res = 0.0
      -- Quadrant 1 - +x, +y, +z
      for i = 0, ntx do
        for j = 0, nty do
          for k = 0, ntz do
            res +=
              [sweep[1]](p_points[{i,j,k}],
                         [p_sub_points[1]][{i,j,k}],
                         [p_x_faces[1]][{  j,k}],
                         [p_y_faces[1]][{i,  k}],
                         [p_z_faces[1]][{i,j  }],
                         [angles[1]],
                         config)
          end
        end
      end
      -- Quadrant 2 - +x, +y, -z
      for i = 0, ntx do
        for j = 0, nty do
          for k = ntz-1, -1, -1 do
            res +=
              [sweep[2]](p_points[{i,j,k}],
                         [p_sub_points[2]][{i,j,k}],
                         [p_x_faces[2]][{  j,k}],
                         [p_y_faces[2]][{i,  k}],
                         [p_z_faces[2]][{i,j  }],
                         [angles[2]],
                         config)
          end
        end
      end
      -- Quadrant 3 - +x, -y, +z
      for i = 0, ntx do
        for j = nty-1, -1, -1 do
          for k = 0, ntz do
            res +=
              [sweep[3]](p_points[{i,j,k}],
                         [p_sub_points[3]][{i,j,k}],
                         [p_x_faces[3]][{  j,k}],
                         [p_y_faces[3]][{i,  k}],
                         [p_z_faces[3]][{i,j  }],
                         [angles[3]],
                         config)
          end
        end
      end
      -- Quadrant 4 - +x, -y, -z
      for i = 0, ntx do
        for j = nty-1, -1, -1 do
          for k = ntz-1, -1, -1 do
            res +=
              [sweep[4]](p_points[{i,j,k}],
                         [p_sub_points[4]][{i,j,k}],
                         [p_x_faces[4]][{  j,k}],
                         [p_y_faces[4]][{i,  k}],
                         [p_z_faces[4]][{i,j  }],
                         [angles[4]],
                         config)
          end
        end
      end
      -- Quadrant 5 - -x, +y, +z
      for i = ntx-1, -1, -1 do
        for j = 0, nty do
          for k = 0, ntz do
            res +=
              [sweep[5]](p_points[{i,j,k}],
                         [p_sub_points[5]][{i,j,k}],
                         [p_x_faces[5]][{  j,k}],
                         [p_y_faces[5]][{i,  k}],
                         [p_z_faces[5]][{i,j  }],
                         [angles[5]],
                         config)
          end
        end
      end
      -- Quadrant 6 - -x, +y, -z
      for i = ntx-1, -1, -1 do
        for j = 0, nty do
          for k = ntz-1, -1, -1 do
            res +=
              [sweep[6]](p_points[{i,j,k}],
                         [p_sub_points[6]][{i,j,k}],
                         [p_x_faces[6]][{  j,k}],
                         [p_y_faces[6]][{i,  k}],
                         [p_z_faces[6]][{i,j  }],
                         [angles[6]],
                         config)
          end
        end
      end
      -- Quadrant 7 - -x, -y, +z
      for i = ntx-1, -1, -1 do
        for j = nty-1, -1, -1 do
          for k = 0, ntz do
            res +=
              [sweep[7]](p_points[{i,j,k}],
                         [p_sub_points[7]][{i,j,k}],
                         [p_x_faces[7]][{  j,k}],
                         [p_y_faces[7]][{i,  k}],
                         [p_z_faces[7]][{i,j  }],
                         [angles[7]],
                         config)
          end
        end
      end
      -- Quadrant 8 - -x, -y, -z
      for i = ntx-1, -1, -1 do
        for j = nty-1, -1, -1 do
          for k = ntz-1, -1, -1 do
            res +=
              [sweep[8]](p_points[{i,j,k}],
                         [p_sub_points[8]][{i,j,k}],
                         [p_x_faces[8]][{  j,k}],
                         [p_y_faces[8]][{i,  k}],
                         [p_z_faces[8]][{i,j  }],
                         [angles[8]],
                         config)
          end
        end
      end

      -- Compute the residual.
      res = sqrt(res/(Nx*Ny*Nz*config.Radiation.angles))

      -- Update intensity.
      for c in tiles do
        reduce_intensity(p_points[c],
                         [p_sub_points:map(function(s) return rexpr s[c] end end)],
                         [angles],
                         config)
      end

    end -- while res > TOLERANCE

  end end -- ComputeRadiationField

return INSTANCE end -- mkInstance

-------------------------------------------------------------------------------
-- MODULE END
-------------------------------------------------------------------------------

return MODULE end
