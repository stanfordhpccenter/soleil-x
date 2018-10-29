import 'regent'

-------------------------------------------------------------------------------
-- MODULE PARAMETERS
-------------------------------------------------------------------------------

return function(MAX_ANGLES_PER_QUAD, Point_columns, SCHEMA) local MODULE = {}

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

local angles = regentlib.newsymbol(region(ispace(int1d), Angle_columns))

local -- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task initialize_angles([angles],
                       config : SCHEMA.Config)
where
  reads writes(angles.{xi, eta, mu, w})
do
  -- Open angles file
  var num_angles = config.Radiation.u.DOM.angles;
  for m = 0, 12345 do
    @ESCAPE for wall = 1, 6 do @EMIT
      if [isWallNormal(wall, rexpr angles[m] end)] then
      end
    @TIME end @EPACSE
  end
end

-------------------------------------------------------------------------------
-- MODULE END
-------------------------------------------------------------------------------

return MODULE end
