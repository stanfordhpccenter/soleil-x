import 'regent'

-------------------------------------------------------------------------------
-- MODULE PARAMETERS
-------------------------------------------------------------------------------

return function(MAX_ANGLES_PER_QUAD, Point, Config) local MODULE = {}

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

-- 'x'|'y'|'z', 1..8 -> regentlib.task
local function mkInitializeFaces(dim, q)

  local -- MANUALLY PARALLELIZED, NO CUDA
  task initialize_faces(faces : region(ispace(int2d), Face_columns),
                        config : Config)
  where
    reads writes(faces.I)
  do
    var num_angles = config.Radiation.angles
    __demand(__openmp)
    for f in faces do
      for m = 0, quadrantSize(q, num_angles) do
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

local -- MANUALLY PARALLELIZED, NO CUDA
task source_term(points : region(ispace(int3d), Point),
                 [angles],
                 config : Config,
                 omega : double)
where
  reads(points.[intensityFields], points.{Ib, sigma}),
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
           * p.[intensityFields[q]][m]
      end
    @TIME end @EPACSE
    p.S = S
  end
end

-- 'x'|'y'|'z', 1..8 -> regentlib.task
local function mkCacheIntensity(dim, q)

  local -- MANUALLY PARALLELIZED, NO CUDA
  task cache_intensity(faces : region(ispace(int2d), Face_columns),
                       config : Config)
  where
    reads(faces.I),
    reads writes(faces.I_prev)
  do
    var num_angles = config.Radiation.angles
    __demand(__openmp)
    for f in faces do
      for m = 0, quadrantSize(q, num_angles) do
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

  local -- MANUALLY PARALLELIZED, NO CUDA
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

  local -- MANUALLY PARALLELIZED, NO CUDA
  task sweep(points : region(ispace(int3d), Point),
             x_faces : region(ispace(int2d), Face_columns),
             y_faces : region(ispace(int2d), Face_columns),
             z_faces : region(ispace(int2d), Face_columns),
             angles : region(ispace(int1d), Angle_columns),
             config : Config)
  where
    reads(angles.{xi, eta, mu}, points.{S, sigma}),
    reads writes(points.[fld], x_faces.I, y_faces.I, z_faces.I)
  do
    var dx = config.Grid.xWidth / config.Radiation.xNum
    var dy = config.Grid.yWidth / config.Radiation.yNum
    var dz = config.Grid.zWidth / config.Radiation.zNum
    var dAx = dy*dz
    var dAy = dx*dz
    var dAz = dx*dy
    var dV = dx*dy*dz
    var [bnd] = points.bounds
    var res = 0.0
    -- Loop over this quadrant's angles
    __demand(__openmp)
    for m in angles do
      -- Use our direction and increments for the sweep
      for k = startz,endz,dindz do
        for j = starty,endy,dindy do
          for i = startx,endx,dindx do
            -- Read upwind face values
            var x_value = x_faces[{  j,k}].I[int(m)]
            var y_value = y_faces[{i,  k}].I[int(m)]
            var z_value = z_faces[{i,j  }].I[int(m)]
            -- Integrate to compute cell-centered value of I
            var oldI = points[{i,j,k}].[fld][int(m)]
            var newI = (points[{i,j,k}].S * dV
                        + fabs(m.xi)  * dAx * x_value/GAMMA
                        + fabs(m.eta) * dAy * y_value/GAMMA
                        + fabs(m.mu)  * dAz * z_value/GAMMA)
                     / (points[{i,j,k}].sigma * dV
                        + fabs(m.xi)  * dAx/GAMMA
                        + fabs(m.eta) * dAy/GAMMA
                        + fabs(m.mu)  * dAz/GAMMA)
            if newI > 0.0 then
              res += pow(newI-oldI,2) / pow(newI,2)
            end
            points[{i,j,k}].[fld][int(m)] = newI
            -- Compute intensities on downwind faces
            x_faces[{  j,k}].I[int(m)] = max(0.0, (newI-(1-GAMMA)*x_value)/GAMMA)
            y_faces[{i,  k}].I[int(m)] = max(0.0, (newI-(1-GAMMA)*y_value)/GAMMA)
            z_faces[{i,j  }].I[int(m)] = max(0.0, (newI-(1-GAMMA)*z_value)/GAMMA)
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

local -- MANUALLY PARALLELIZED, NO CUDA
task reduce_intensity(points : region(ispace(int3d), Point),
                      [angles],
                      config : Config)
where
  reads(points.[intensityFields]),
  [angles:map(function(a)
     return regentlib.privilege(regentlib.reads, a, 'w')
   end)],
  writes(points.G)
do
  var num_angles = config.Radiation.angles
  __demand(__openmp)
  for p in points do
    var G = 0.0;
    @ESCAPE for q = 1, 8 do @EMIT
      for m = 0, quadrantSize(q, num_angles) do
        G += [angles[q]][m].w * p.[intensityFields[q]][m]
      end
    @TIME end @EPACSE
    p.G = G
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

  local angles = UTIL.generate(8, regentlib.newsymbol)

  local x_faces = UTIL.generate(8, regentlib.newsymbol)
  local y_faces = UTIL.generate(8, regentlib.newsymbol)
  local z_faces = UTIL.generate(8, regentlib.newsymbol)

  local x_tiles = regentlib.newsymbol('x_tiles')
  local y_tiles = regentlib.newsymbol('y_tiles')
  local z_tiles = regentlib.newsymbol('z_tiles')

  local p_x_faces = UTIL.generate(8, regentlib.newsymbol)
  local p_y_faces = UTIL.generate(8, regentlib.newsymbol)
  local p_z_faces = UTIL.generate(8, regentlib.newsymbol)

  function INSTANCE.DeclSymbols(config) return rquote

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

    -- Regions for faces
    var grid_x = ispace(int2d, {   Ny,Nz})
    var grid_y = ispace(int2d, {Nx,   Nz})
    var grid_z = ispace(int2d, {Nx,Ny   });
    @ESCAPE for q = 1, 8 do @EMIT
      var [x_faces[q]] = region(grid_x, Face_columns);
      [UTIL.mkRegionTagAttach(x_faces[q], MAPPER.SAMPLE_ID_TAG, sampleId, int)];
      var [y_faces[q]] = region(grid_y, Face_columns);
      [UTIL.mkRegionTagAttach(y_faces[q], MAPPER.SAMPLE_ID_TAG, sampleId, int)];
      var [z_faces[q]] = region(grid_z, Face_columns);
      [UTIL.mkRegionTagAttach(z_faces[q], MAPPER.SAMPLE_ID_TAG, sampleId, int)];
    @TIME end @EPACSE

    -- Regions for angle values
    @ESCAPE for q = 1, 8 do @EMIT
      var is = ispace(int1d, quadrantSize(q, config.Radiation.angles))
      var [angles[q]] = region(is, Angle_columns);
      [UTIL.mkRegionTagAttach(angles[q], MAPPER.SAMPLE_ID_TAG, sampleId, int)];
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
        var rect = rect2d{lo = int2d{(Ny/nty)*a,       (Nz/ntz)*b      },
                          hi = int2d{(Ny/nty)*(a+1)-1, (Nz/ntz)*(b+1)-1}}
        regentlib.c.legion_domain_point_coloring_color_domain(x_coloring, c, rect)
      end
      var [p_x_faces[q]] = partition(disjoint, [x_faces[q]], x_coloring, x_tiles)
      regentlib.c.legion_domain_point_coloring_destroy(x_coloring);
      -- y
      var y_coloring = regentlib.c.legion_domain_point_coloring_create()
      for c in y_tiles do
        var a = c.x
        var b = c.y
        var rect = rect2d{lo = int2d{(Nx/ntx)*a,       (Nz/ntz)*b      },
                          hi = int2d{(Nx/ntx)*(a+1)-1, (Nz/ntz)*(b+1)-1}}
        regentlib.c.legion_domain_point_coloring_color_domain(y_coloring, c, rect)
      end
      var [p_y_faces[q]] = partition(disjoint, [y_faces[q]], y_coloring, y_tiles)
      regentlib.c.legion_domain_point_coloring_destroy(y_coloring);
      -- z
      var z_coloring = regentlib.c.legion_domain_point_coloring_create()
      for c in z_tiles do
        var a = c.x
        var b = c.y
        var rect = rect2d{lo = int2d{(Nx/ntx)*a,       (Ny/nty)*b      },
                          hi = int2d{(Nx/ntx)*(a+1)-1, (Ny/nty)*(b+1)-1}}
        regentlib.c.legion_domain_point_coloring_color_domain(z_coloring, c, rect)
      end
      var [p_z_faces[q]] = partition(disjoint, [z_faces[q]], z_coloring, z_tiles)
      regentlib.c.legion_domain_point_coloring_destroy(z_coloring);
    @TIME end @EPACSE

  end end -- DeclSymbols

  function INSTANCE.InitRegions(config) return rquote

    initialize_angles([angles], config);
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

  end end -- InitRegions

  function INSTANCE.ComputeRadiationField(config, tiles, p_points) return rquote

    var omega = config.Radiation.qs/(config.Radiation.qa+config.Radiation.qs)

    -- Compute until convergence
    var res : double = 1.0
    while res > TOLERANCE do

      -- Update the source term (in this problem, isotropic)
      for c in tiles do
        source_term(p_points[c], [angles], config, omega)
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

      --Perform the sweep for computing new intensities
      res = 0.0
      -- Quadrant 1 - +x, +y, +z
      for i = 0, ntx do
        for j = 0, nty do
          for k = 0, ntz do
            res +=
              [sweep[1]](p_points[{i,j,k}],
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
                         [p_x_faces[8]][{  j,k}],
                         [p_y_faces[8]][{i,  k}],
                         [p_z_faces[8]][{i,j  }],
                         [angles[8]],
                         config)
          end
        end
      end

      -- Compute the residual
      res = sqrt(res/(Nx*Ny*Nz*config.Radiation.angles))

    end -- while res > TOLERANCE

    -- Reduce intensity
    for c in tiles do
      reduce_intensity(p_points[c], [angles], config)
    end

  end end -- ComputeRadiationField

return INSTANCE end -- mkInstance

-------------------------------------------------------------------------------
-- MODULE END
-------------------------------------------------------------------------------

return MODULE end
