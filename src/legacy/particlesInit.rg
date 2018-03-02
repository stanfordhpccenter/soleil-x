import 'regent'

local A = require 'admiral'

-------------------------------------------------------------------------------
-- Module parameters
-------------------------------------------------------------------------------

return function(particlesRel, cellsRel,
                xBnumGlobal, yBnumGlobal, zBnumGlobal)

-------------------------------------------------------------------------------
-- Local tasks
-------------------------------------------------------------------------------

local __demand(__parallel) task InitParticlesUniform
  (particles : particlesRel:regionType(),
   cells : cellsRel:regionType(),
   config : A.configStruct(),
   xBnum : int,
   yBnum : int,
   zBnum : int)
where
  reads writes(particles),
  reads(cells.{velocity, centerCoordinates})
do
  var pBase = 0
  for p in particles do
    pBase = [int32](p)
    break
  end
  var lo : int3d = cells.bounds.lo
  lo.x = max(lo.x, xBnum)
  lo.y = max(lo.y, yBnum)
  lo.z = max(lo.z, zBnum)
  var hi : int3d = cells.bounds.hi
  hi.x = min(hi.x, config.Grid.xNum + xBnum - 1)
  hi.y = min(hi.y, config.Grid.yNum + yBnum - 1)
  hi.z = min(hi.z, config.Grid.zNum + zBnum - 1)
  var xSize = hi.x - lo.x + 1
  var ySize = hi.y - lo.y + 1
  var particlesPerTask =
    config.Particles.initNum
    / (config.Grid.xTiles * config.Grid.yTiles * config.Grid.zTiles)
  __demand(__openmp)
  for p in particles do
    if [int32](p) - pBase < particlesPerTask then
      p.__valid = true
      var relIdx = [int32](p) - pBase
      var c : int3d = { lo.x + relIdx % xSize,
                        lo.y + relIdx / xSize % ySize,
                        lo.z + relIdx / xSize / ySize }
      p.cell = c
      p.position = cells[p.cell].centerCoordinates
      p.velocity = cells[p.cell].velocity
      p.density = config.Particles.density
      p.temperature = config.Particles.initTemperature
      p.diameter = config.Particles.diameterMean
    end
  end
end
A.registerTask(InitParticlesUniform, 'InitParticlesUniform')

-------------------------------------------------------------------------------
-- Exported quotes
-------------------------------------------------------------------------------

local exports = {}

exports.InitParticlesUniform = rquote
  InitParticlesUniform([particlesRel:regionSymbol()],
                       [cellsRel:regionSymbol()],
                       [A.configSymbol()],
                       [xBnumGlobal:varSymbol()],
                       [yBnumGlobal:varSymbol()],
                       [zBnumGlobal:varSymbol()])
end

-------------------------------------------------------------------------------
-- Module exports
-------------------------------------------------------------------------------

return exports
end
