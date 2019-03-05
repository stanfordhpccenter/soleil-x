import 'regent'

-------------------------------------------------------------------------------
-- MODULE PARAMETERS
-------------------------------------------------------------------------------

return function(MAX_ANGLES_PER_QUAD, Point_columns, SCHEMA) local MODULE = {}

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
  I      : double[MAX_ANGLES_PER_QUAD];
  I_prev : double[MAX_ANGLES_PER_QUAD];
}

local struct GridMap_columns {
  p_to_s3d : int3d;
  s3d_to_p : int3d;
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

-- 1..6, regentlib.rexpr -> regentlib.rexpr
local function isWallNormal(wall, angle)
  return terralib.newlist{
    rexpr angle.xi ==  1.0 and angle.eta ==  0.0 and angle.mu ==  0.0 end,
    rexpr angle.xi == -1.0 and angle.eta ==  0.0 and angle.mu ==  0.0 end,
    rexpr angle.xi ==  0.0 and angle.eta ==  1.0 and angle.mu ==  0.0 end,
    rexpr angle.xi ==  0.0 and angle.eta == -1.0 and angle.mu ==  0.0 end,
    rexpr angle.xi ==  0.0 and angle.eta ==  0.0 and angle.mu ==  1.0 end,
    rexpr angle.xi ==  0.0 and angle.eta ==  0.0 and angle.mu == -1.0 end,
  }[wall]
end

-------------------------------------------------------------------------------
-- MODULE-LOCAL TASKS
-------------------------------------------------------------------------------

local angles = UTIL.generate(8, function()
  return regentlib.newsymbol(region(ispace(int1d), Angle_columns))
end)

local -- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task initialize_angles([angles],
                       config : SCHEMA.Config)
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
  var num_angles = config.Radiation.u.DOM.angles
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
  -- Close angles file.
  C.fclose(f);
  -- Check that angles are partitioned correctly into quadrants.
  @ESCAPE for q = 1, 8 do @EMIT
    for m = 0, quadrantSize(q, num_angles) do
      regentlib.assert([angleInQuadrant(q, rexpr [angles[q]][m] end)],
                       'Angle in wrong quadrant')
    end
  @TIME end @EPACSE
  -- Check that normals exist for all walls.
  var normalExists = array(false, false, false, false, false, false);
  @ESCAPE for q = 1, 8 do @EMIT
    for m = 0, quadrantSize(q, num_angles) do
      @ESCAPE for wall = 1, 6 do @EMIT
        if [isWallNormal(wall, rexpr [angles[q]][m] end)] then
          normalExists[wall-1] = true
        end
      @TIME end @EPACSE
    end
  @TIME end @EPACSE
  regentlib.assert(normalExists[0], 'Normal missing for wall xLo')
  regentlib.assert(normalExists[1], 'Normal missing for wall xHi')
  regentlib.assert(normalExists[2], 'Normal missing for wall yLo')
  regentlib.assert(normalExists[3], 'Normal missing for wall yHi')
  regentlib.assert(normalExists[4], 'Normal missing for wall zLo')
  regentlib.assert(normalExists[5], 'Normal missing for wall zHi')
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

local __demand(__cuda) -- MANUALLY PARALLELIZED
task initialize_geometry(points : region(ispace(int3d), Point_columns),
                         config : SCHEMA.Config)
where
  reads writes(points.{cellWidth})
do
  __demand(__openmp)
  for p in points do
    p.cellWidth[0] = config.Grid.xWidth / config.Radiation.u.DOM.xNum
    p.cellWidth[1] = config.Grid.yWidth / config.Radiation.u.DOM.yNum
    p.cellWidth[2] = config.Grid.zWidth / config.Radiation.u.DOM.zNum
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
-- * Follow the s3d_to_p field.

local -- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task cache_grid_translation(grid_map : region(ispace(int3d), GridMap_columns),
                            sub_point_offsets : region(ispace(int1d), bool),
                            diagonals : ispace(int1d))
where
  writes(grid_map.{p_to_s3d, s3d_to_p})
do
  var Tx = grid_map.bounds.hi.x + 1
  var Ty = grid_map.bounds.hi.y + 1
  var Tz = grid_map.bounds.hi.z + 1
  regentlib.assert(
    grid_map.bounds.lo.x == 0 and
    grid_map.bounds.lo.y == 0 and
    grid_map.bounds.lo.z == 0 and
    int(sub_point_offsets.bounds.lo) == 0 and
    int(sub_point_offsets.bounds.hi + 1) == MAX_ANGLES_PER_QUAD*Tx*Ty*Tz and
    int(diagonals.bounds.lo) == 0 and
    int(diagonals.bounds.hi) == (Tx-1)+(Ty-1)+(Tz-1),
    'Internal error')
  var coloring = regentlib.c.legion_domain_point_coloring_create()
  var rect_start = 0
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
      grid_map[diag].p_to_s3d = grid
      grid_map[grid].s3d_to_p = diag
      -- Advance diagonal-order index
      if diag.x > 0 and diag.y < Ty-1 then
        diag.x -= 1
        diag.y += 1
      elseif diag.z < min(d, Tz-1) then
        diag.z += 1
        diag.x = min(d-diag.z, Tx-1)
        diag.y = d-diag.z-diag.x
      else
        -- We've run out of indices on this diagonal, color it on the sub-point
        -- offsets and continue to the next one
        var rect_end = MAX_ANGLES_PER_QUAD
                     + MAX_ANGLES_PER_QUAD * grid.x
                     + MAX_ANGLES_PER_QUAD * Tx     * grid.y
                     + MAX_ANGLES_PER_QUAD * Tx     * Ty     * grid.z
        regentlib.c.legion_domain_point_coloring_color_domain(
          coloring, int1d(d), rect1d{ lo = rect_start, hi = rect_end - 1 })
        rect_start = rect_end
        break
      end
    end
  end
  regentlib.assert(grid.x == Tx-1 and
                   grid.y == Ty-1 and
                   grid.z == Tz-1 and
                   rect_start == int(sub_point_offsets.bounds.hi + 1),
                   'Internal error')
  -- Construct & return partition of sub-point offsets
  var p = partition(disjoint, sub_point_offsets, coloring, diagonals)
  regentlib.c.legion_domain_point_coloring_destroy(coloring)
  return p
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
                        config : SCHEMA.Config)
  where
    reads writes(faces.I)
  do
    var num_angles = config.Radiation.u.DOM.angles
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

local __demand(__cuda) -- MANUALLY PARALLELIZED
task source_term(points : region(ispace(int3d), Point_columns),
                 config : SCHEMA.Config)
where
  reads(points.{Ib, sigma, G}),
  writes(points.S)
do
  var omega =
    config.Radiation.u.DOM.qs /
    (config.Radiation.u.DOM.qa + config.Radiation.u.DOM.qs)
  __demand(__openmp)
  for p in points do
    p.S = (1.0-omega) * p.sigma * p.Ib + omega * p.sigma/(4.0*PI) * p.G
  end
end

-- 'x'|'y'|'z', 1..8 -> regentlib.task
local function mkCacheIntensity(dim, q)

  local __demand(__cuda) -- MANUALLY PARALLELIZED
  task cache_intensity(faces : region(ispace(int2d), Face_columns),
                       config : SCHEMA.Config)
  where
    reads(faces.I),
    reads writes(faces.I_prev)
  do
    var num_angles = config.Radiation.u.DOM.angles
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
  local intensityField = terralib.newlist{
    'xLoIntensity', 'xHiIntensity',
    'yLoIntensity', 'yHiIntensity',
    'zLoIntensity', 'zHiIntensity'
  }[wall]
  local windowField = terralib.newlist{
    'xLoWindow', 'xHiWindow', 'yLoWindow', 'yHiWindow', 'zLoWindow', 'zHiWindow'
  }[wall]
  local incomingQuadrants = terralib.newlist{
    terralib.newlist{5, 6, 7, 8}, -- xi <= 0
    terralib.newlist{1, 2, 3, 4}, -- xi >= 0
    terralib.newlist{3, 4, 7, 8}, -- eta <= 0
    terralib.newlist{1, 2, 5, 6}, -- eta >= 0
    terralib.newlist{2, 4, 6, 8}, -- mu <= 0
    terralib.newlist{1, 3, 5, 7}, -- mu >= 0
  }[wall]
  local outgoingQuadrants = terralib.newlist{
    terralib.newlist{1, 2, 3, 4}, -- xi >= 0
    terralib.newlist{5, 6, 7, 8}, -- xi <= 0
    terralib.newlist{1, 2, 5, 6}, -- eta >= 0
    terralib.newlist{3, 4, 7, 8}, -- eta <= 0
    terralib.newlist{1, 3, 5, 7}, -- mu >= 0
    terralib.newlist{2, 4, 6, 8}, -- mu <= 0
  }[wall]

  local faces = UTIL.generate(8, function()
    return regentlib.newsymbol(region(ispace(int2d), Face_columns))
  end)

  local __demand(__cuda) -- MANUALLY PARALLELIZED
  task bound([faces],
             [angles],
             config : SCHEMA.Config)
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
    var epsw = config.Radiation.u.DOM.[emissField]
    var Tw = config.Radiation.u.DOM.[tempField]
    var incidentI = config.Radiation.u.DOM.[intensityField]
    var fromCell = config.Radiation.u.DOM.[windowField].fromCell
    var uptoCell = config.Radiation.u.DOM.[windowField].uptoCell
    var num_angles = config.Radiation.u.DOM.angles
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
      value += epsw*SB*pow(Tw,4.0)/PI;
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
            var I = value
            -- Add incident radiation on the wall normal
            if fromCell[0] <= a and a <= uptoCell[0] and
               fromCell[1] <= b and b <= uptoCell[1] and
               [isWallNormal(wall, rexpr [angles[q]][m] end)] then
              I += incidentI
            end
            [faces[q]][idx].I[m] = I
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

  local __demand(__cuda) -- MANUALLY PARALLELIZED
  task sweep(points : region(ispace(int3d), Point_columns),
             sub_points : region(ispace(int1d), SubPoint_columns),
             grid_map : region(ispace(int3d), GridMap_columns),
             sub_point_offsets : region(ispace(int1d), bool),
             diagonals : ispace(int1d),
             p_sub_point_offsets : partition(disjoint, sub_point_offsets, diagonals),
             x_faces : region(ispace(int2d), Face_columns),
             y_faces : region(ispace(int2d), Face_columns),
             z_faces : region(ispace(int2d), Face_columns),
             angles : region(ispace(int1d), Angle_columns),
             config : SCHEMA.Config)
  where
    reads(angles.{xi, eta, mu}, points.{cellWidth, S, sigma}, grid_map.s3d_to_p),
    reads writes(sub_points.I, x_faces.I, y_faces.I, z_faces.I)
  do
    var Tx = points.bounds.hi.x - points.bounds.lo.x + 1
    var Ty = points.bounds.hi.y - points.bounds.lo.y + 1
    var Tz = points.bounds.hi.z - points.bounds.lo.z + 1
    regentlib.assert(
      int(sub_points.bounds.hi - sub_points.bounds.lo + 1)
      == MAX_ANGLES_PER_QUAD*Tx*Ty*Tz and
      grid_map.bounds.lo.x == 0 and grid_map.bounds.hi.x + 1 == Tx and
      grid_map.bounds.lo.y == 0 and grid_map.bounds.hi.y + 1 == Ty and
      grid_map.bounds.lo.z == 0 and grid_map.bounds.hi.z + 1 == Tz and
      int(sub_point_offsets.bounds.lo) == 0 and
      int(sub_point_offsets.bounds.hi + 1) == MAX_ANGLES_PER_QUAD*Tx*Ty*Tz and
      x_faces.bounds.hi.x - x_faces.bounds.lo.x + 1 == Ty and
      x_faces.bounds.hi.y - x_faces.bounds.lo.y + 1 == Tz and
      y_faces.bounds.hi.x - y_faces.bounds.lo.x + 1 == Tx and
      y_faces.bounds.hi.y - y_faces.bounds.lo.y + 1 == Tz and
      z_faces.bounds.hi.x - z_faces.bounds.lo.x + 1 == Tx and
      z_faces.bounds.hi.y - z_faces.bounds.lo.y + 1 == Ty,
      'Internal error')
    var num_angles = config.Radiation.u.DOM.angles
    var res = 0.0
    -- Launch in order of intra-tile diagonals
    for d = int(diagonals.bounds.lo), int(diagonals.bounds.hi+1) do
      __demand(__openmp)
      for s1d_off in p_sub_point_offsets[d] do
        -- Compute sub-point index, translate to point index
        var m = int(s1d_off) % MAX_ANGLES_PER_QUAD
        if m < quadrantSize(q, num_angles) then
          var s3d_off = int3d{s1d_off / MAX_ANGLES_PER_QUAD % Tx,
                              s1d_off / MAX_ANGLES_PER_QUAD / Tx % Ty,
                              s1d_off / MAX_ANGLES_PER_QUAD / Tx / Ty}
          var p_off = grid_map[s3d_off].s3d_to_p
          p_off = int3d{
            [directions[q][1] and rexpr p_off.x end or rexpr Tx-p_off.x-1 end],
            [directions[q][2] and rexpr p_off.y end or rexpr Ty-p_off.y-1 end],
            [directions[q][3] and rexpr p_off.z end or rexpr Tz-p_off.z-1 end]}
          var s1d = sub_points.bounds.lo + s1d_off
          var p = points.bounds.lo + p_off
          -- Read upwind face values
          var x_value = x_faces[{    p.y,p.z}].I[m]
          var y_value = y_faces[{p.x,    p.z}].I[m]
          var z_value = z_faces[{p.x,p.y    }].I[m]
          -- Integrate to compute cell-centered value of I
          var oldI = sub_points[s1d].I
          -- TODO update for non uniform mesh
          var dx = points[p].cellWidth[0]
          var dy = points[p].cellWidth[1]
          var dz = points[p].cellWidth[2]
          var dAx = dy*dz
          var dAy = dx*dz
          var dAz = dx*dy
          var dV  = dx*dy*dz
          var newI = (points[p].S * dV
                      + fabs(angles[m].xi)  * dAx * x_value/GAMMA
                      + fabs(angles[m].eta) * dAy * y_value/GAMMA
                      + fabs(angles[m].mu)  * dAz * z_value/GAMMA)
                   / (points[p].sigma * dV
                      + fabs(angles[m].xi)  * dAx/GAMMA
                      + fabs(angles[m].eta) * dAy/GAMMA
                      + fabs(angles[m].mu)  * dAz/GAMMA)
          if newI > 0.0 then
            -- TODO update for non uniform mesh
            res += ( pow(newI-oldI,2) / pow(newI,2) ) * dV
          end
          sub_points[s1d].I = newI
          -- Compute intensities on downwind faces
          x_faces[{    p.y,p.z}].I[m] = max(0.0, (newI-(1-GAMMA)*x_value)/GAMMA)
          y_faces[{p.x,    p.z}].I[m] = max(0.0, (newI-(1-GAMMA)*y_value)/GAMMA)
          z_faces[{p.x,p.y    }].I[m] = max(0.0, (newI-(1-GAMMA)*z_value)/GAMMA)
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
                      grid_map : region(ispace(int3d), GridMap_columns),
                      [angles],
                      config : SCHEMA.Config)
where
  reads(grid_map.p_to_s3d),
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
  regentlib.assert(
    grid_map.bounds.lo.x == 0 and grid_map.bounds.hi.x + 1 == Tx and
    grid_map.bounds.lo.y == 0 and grid_map.bounds.hi.y + 1 == Ty and
    grid_map.bounds.lo.z == 0 and grid_map.bounds.hi.z + 1 == Tz,
    'Internal error');
  @ESCAPE for q = 1, 8 do @EMIT
    regentlib.assert(
      int([sub_points[q]].bounds.hi - [sub_points[q]].bounds.lo + 1)
      == MAX_ANGLES_PER_QUAD*Tx*Ty*Tz,
      'Internal error')
  @TIME end @EPACSE
  var num_angles = config.Radiation.u.DOM.angles
  __demand(__openmp)
  for p in points do
    p.G = 0.0
  end
  @ESCAPE for q = 1, 8 do @EMIT
    __demand(__openmp)
    for p in points do
      var G = 0.0
      var p_off = p - points.bounds.lo
      p_off = int3d{
        [directions[q][1] and rexpr p_off.x end or rexpr Tx-p_off.x-1 end],
        [directions[q][2] and rexpr p_off.y end or rexpr Ty-p_off.y-1 end],
        [directions[q][3] and rexpr p_off.z end or rexpr Tz-p_off.z-1 end]}
      var s3d_off = grid_map[p_off].p_to_s3d
      var s1d_off = MAX_ANGLES_PER_QUAD * s3d_off.x
                  + MAX_ANGLES_PER_QUAD * Tx        * s3d_off.y
                  + MAX_ANGLES_PER_QUAD * Tx        * Ty        * s3d_off.z
      var s1d = [sub_points[q]].bounds.lo + s1d_off
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
  local p_sub_points = UTIL.generate(8, regentlib.newsymbol)

  local x_faces = UTIL.generate(8, regentlib.newsymbol)
  local y_faces = UTIL.generate(8, regentlib.newsymbol)
  local z_faces = UTIL.generate(8, regentlib.newsymbol)
  local x_tiles = regentlib.newsymbol('x_tiles')
  local y_tiles = regentlib.newsymbol('y_tiles')
  local z_tiles = regentlib.newsymbol('z_tiles')
  local p_x_faces = UTIL.generate(8, regentlib.newsymbol)
  local p_y_faces = UTIL.generate(8, regentlib.newsymbol)
  local p_z_faces = UTIL.generate(8, regentlib.newsymbol)

  local angles = UTIL.generate(8, regentlib.newsymbol)

  local grid_map = regentlib.newsymbol('grid_map')
  local sub_point_offsets = regentlib.newsymbol('sub_point_offsets')
  local diagonals = regentlib.newsymbol('diagonals')
  local p_sub_point_offsets = regentlib.newsymbol('p_sub_point_offsets')

  -- NOTE: This quote is included into the main simulation whether or not
  -- we're using DOM, so the values will be garbage if type ~= DOM.
  function INSTANCE.DeclSymbols(config, tiles) return rquote

    var sampleId = config.Mapping.sampleId

    -- Number of tiles in each dimension
    var [ntx] = config.Mapping.tiles[0]
    var [nty] = config.Mapping.tiles[1]
    var [ntz] = config.Mapping.tiles[2]
    -- Number of points in each dimension
    var [Nx] = ntx
    var [Ny] = nty
    var [Nz] = ntz
    if config.Radiation.type == SCHEMA.RadiationModel_DOM then
      Nx = config.Radiation.u.DOM.xNum
      Ny = config.Radiation.u.DOM.yNum
      Nz = config.Radiation.u.DOM.zNum
    end
    -- Sanity-check partitioning
    regentlib.assert(Nx % ntx == 0, "Uneven partitioning of radiation grid on x")
    regentlib.assert(Ny % nty == 0, "Uneven partitioning of radiation grid on y")
    regentlib.assert(Nz % ntz == 0, "Uneven partitioning of radiation grid on z")
    -- Number of points in each tile
    var [Tx] = Nx / ntx
    var [Ty] = Ny / nty
    var [Tz] = Nz / ntz

    -- Regions for points
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
    var num_angles = 8
    if config.Radiation.type == SCHEMA.RadiationModel_DOM then
      num_angles = config.Radiation.u.DOM.angles
    end
    @ESCAPE for q = 1, 8 do @EMIT
      var is_angles = ispace(int1d, quadrantSize(q, num_angles))
      var [angles[q]] = region(is_angles, Angle_columns);
      [UTIL.emitRegionTagAttach(angles[q], MAPPER.SAMPLE_ID_TAG, sampleId, int)];
    @TIME end @EPACSE

    -- Regions for intra-tile information
    var is_sub_point_offsets = ispace(int1d, MAX_ANGLES_PER_QUAD*Tx*Ty*Tz)
    var [sub_point_offsets] = region(is_sub_point_offsets, bool);
    [UTIL.emitRegionTagAttach(sub_point_offsets, MAPPER.SAMPLE_ID_TAG, sampleId, int)];
    var is_grid_map = ispace(int3d, {Tx,Ty,Tz})
    var [grid_map] = region(is_grid_map, GridMap_columns);
    [UTIL.emitRegionTagAttach(grid_map, MAPPER.SAMPLE_ID_TAG, sampleId, int)];

    -- Partition points
    -- (done by the host code)

    -- Partition sub-points
    @ESCAPE for q = 1, 8 do @EMIT
      var [p_sub_points[q]] =
        [UTIL.mkPartitionEqually(int1d, int3d, SubPoint_columns)]
        ([sub_points[q]], tiles, 0)
    @TIME end @EPACSE

    -- Partition faces
    var [x_tiles] = ispace(int2d, {    nty,ntz})
    var [y_tiles] = ispace(int2d, {ntx,    ntz})
    var [z_tiles] = ispace(int2d, {ntx,nty    });
    @ESCAPE for q = 1, 8 do @EMIT
      var [p_x_faces[q]] =
        [UTIL.mkPartitionEqually(int2d, int2d, Face_columns)]
        ([x_faces[q]], x_tiles, 0, 0)
      var [p_y_faces[q]] =
        [UTIL.mkPartitionEqually(int2d, int2d, Face_columns)]
        ([y_faces[q]], y_tiles, 0, 0)
      var [p_z_faces[q]] =
        [UTIL.mkPartitionEqually(int2d, int2d, Face_columns)]
        ([z_faces[q]], z_tiles, 0, 0)
    @TIME end @EPACSE

    -- Cache intra-tile information
    var [diagonals] = ispace(int1d, (Tx-1)+(Ty-1)+(Tz-1)+1)
    var [p_sub_point_offsets] =
      cache_grid_translation(grid_map, sub_point_offsets, diagonals)

  end end -- DeclSymbols

  function INSTANCE.InitRegions(config, tiles, p_points) return rquote

    -- Initialize points
    for c in tiles do
      initialize_points(p_points[c])
    end

    -- TEST Initialize geometry
    for c in tiles do
      initialize_geometry(p_points[c], config)
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
    initialize_angles([angles], config);

  end end -- InitRegions

  function INSTANCE.ComputeRadiationField(config, tiles, p_points) return rquote

    -- Initialize intensity.
    for c in tiles do
      reduce_intensity(p_points[c],
                       [p_sub_points:map(function(s) return rexpr s[c] end end)],
                       grid_map,
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
      res = 0.0;
      @ESCAPE for q = 1, 8 do @EMIT
        for i = [directions[q][1] and rexpr   0 end or rexpr ntx-1 end],
                [directions[q][1] and rexpr ntx end or rexpr    -1 end],
                [directions[q][1] and rexpr   1 end or rexpr    -1 end] do
          for j = [directions[q][2] and rexpr   0 end or rexpr nty-1 end],
                  [directions[q][2] and rexpr nty end or rexpr    -1 end],
                  [directions[q][2] and rexpr   1 end or rexpr    -1 end] do
            for k = [directions[q][3] and rexpr   0 end or rexpr ntz-1 end],
                    [directions[q][3] and rexpr ntz end or rexpr    -1 end],
                    [directions[q][3] and rexpr   1 end or rexpr    -1 end] do
              res +=
                [sweep[q]](p_points[{i,j,k}],
                           [p_sub_points[q]][{i,j,k}],
                           grid_map,
                           sub_point_offsets,
                           diagonals,
                           p_sub_point_offsets,
                           [p_x_faces[q]][{  j,k}],
                           [p_y_faces[q]][{i,  k}],
                           [p_z_faces[q]][{i,j  }],
                           [angles[q]],
                           config)
            end
          end
        end
      @TIME end @EPACSE

      -- Compute the residual.
      res = sqrt(res/(config.Grid.xWidth*config.Grid.yWidth*config.Grid.zWidth*config.Radiation.u.DOM.angles))

      -- Update intensity.
      for c in tiles do
        reduce_intensity(p_points[c],
                         [p_sub_points:map(function(s) return rexpr s[c] end end)],
                         grid_map,
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
