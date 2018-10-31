import 'regent'

-- It seems like we need to generate a certain amount of allocations for this to trigger.
local SCHEMA = {}
for i = 0, 40000 do SCHEMA[i] = {"asdf", "qwer" .. tostring(i)} end

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
  }[wall]
end

-------------------------------------------------------------------------------
-- MODULE-LOCAL TASKS
-------------------------------------------------------------------------------

task initialize_angles(angles : region(ispace(int1d), int))
  for m = 0, 12345 do
    [(function() local __quotes = terralib.newlist() for wall = 1, 3 do __quotes:insert(rquote
      if [isWallNormal(wall, rexpr angles[m] end)] then
      end
    end) end return __quotes end)()]
  end
end
