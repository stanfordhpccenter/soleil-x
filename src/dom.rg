import 'regent'

-------------------------------------------------------------------------------
-- MODULE PARAMETERS
-------------------------------------------------------------------------------

return function(MAX_ANGLES_PER_QUAD, Point_columns, SCHEMA) local MODULE = {}

-------------------------------------------------------------------------------
-- IMPORTS
-------------------------------------------------------------------------------

local UTIL = require 'util'

-------------------------------------------------------------------------------
-- MODULE-LOCAL FIELD SPACES
-------------------------------------------------------------------------------

local struct Angle_columns {
  xi  : double;
  eta : double;
  mu  : double;
  w   : double;
}

-------------------------------------------------------------------------------
-- QUADRANT MACROS
-------------------------------------------------------------------------------

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
  var num_angles = config.Radiation.u.DOM.angles;
  @ESCAPE for q = 1, 8 do @EMIT
    for m = 0, quadrantSize(q, num_angles) do
      @ESCAPE for wall = 1, 6 do @EMIT
        if [isWallNormal(wall, rexpr [angles[q]][m] end)] then
        end
      @TIME end @EPACSE
    end
  @TIME end @EPACSE
end

-------------------------------------------------------------------------------
-- MODULE END
-------------------------------------------------------------------------------

return MODULE end
