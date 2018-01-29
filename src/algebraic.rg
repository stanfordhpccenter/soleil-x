import 'regent'

local A = require 'admiral'

-------------------------------------------------------------------------------
-- Module parameters
-------------------------------------------------------------------------------

return function(particlesRel)

-------------------------------------------------------------------------------
-- Local tasks
-------------------------------------------------------------------------------

local acos = regentlib.acos(double)
local pow = regentlib.pow(double)
local pi = rexpr 2.0 * acos(0.0) end

local __demand(__parallel, __cuda)
task AddRadiation(particles : particlesRel:regionType(),
                  config : A.configStruct())
where
  reads(particles.{density, diameter}),
  reads writes(particles.temperature_t)
do
  __demand(__openmp)
  for p in particles do
    var crossSectionArea = pi * pow(p.diameter, 2.0) / 4.0
    var volume = pi * pow(p.diameter, 3.0) / 6.0
    var mass = volume * p.density
    var absorbedRadiationIntensity =
      config.Particles.absorptivity
      * config.Radiation.intensity
      * crossSectionArea
    p.temperature_t +=
      absorbedRadiationIntensity / (mass * config.Particles.heatCapacity)
  end
end
A.registerTask(AddRadiation, 'AddRadiation')

-------------------------------------------------------------------------------
-- Exported quotes
-------------------------------------------------------------------------------

local exports = {}

exports.AddRadiation = rquote
  AddRadiation([particlesRel:regionSymbol()],
               [A.configSymbol()])
end

-------------------------------------------------------------------------------
-- Module exports
-------------------------------------------------------------------------------

return exports
end
