import 'regent'

local A = require 'admiral'

-------------------------------------------------------------------------------
-- Module parameters
-------------------------------------------------------------------------------

return function(particlesRel, absorptivity, heatCapacity, radiationIntensity)

-------------------------------------------------------------------------------
-- Local tasks
-------------------------------------------------------------------------------

local acos = regentlib.acos(double)
local pow = regentlib.pow(double)
local pi = rexpr 2.0 * acos(0.0) end

local __demand(__parallel) task AddRadiation
  (particles : particlesRel:regionType(),
   absorptivity : double,
   heatCapacity : double,
   radiationIntensity : double)
where
  reads(particles.{density, diameter}),
  reads writes(particles.temperature_t)
do
  for p in particles do
    var crossSectionArea = pi * pow(p.diameter, 2.0) / 4.0
    var volume = pi * pow(p.diameter, 3.0) / 6.0
    var mass = volume * p.density
    var absorbedRadiationIntensity =
      absorptivity * radiationIntensity * crossSectionArea
    p.temperature_t += absorbedRadiationIntensity / (mass * heatCapacity)
  end
end
A.registerTask(AddRadiation, 'AddRadiation')

-------------------------------------------------------------------------------
-- Exported quotes
-------------------------------------------------------------------------------

local exports = {}

exports.AddRadiation = rquote
  AddRadiation([particlesRel:regionSymbol()],
               [absorptivity:varSymbol()],
               [heatCapacity:varSymbol()],
               [radiationIntensity:varSymbol()])
end

-------------------------------------------------------------------------------
-- Module exports
-------------------------------------------------------------------------------

return exports
end
