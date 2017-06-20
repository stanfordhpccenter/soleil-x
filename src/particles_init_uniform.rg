import 'regent'

-------------------------------------------------------------------------------
-- Module parameters
-------------------------------------------------------------------------------

return function(particlesRel, cellsRel)

-------------------------------------------------------------------------------
-- Compile-time computation
-------------------------------------------------------------------------------

local config
local i = 1
while i <= #arg do
  if arg[i] == '-i' and i < #arg then
    config = loadfile(arg[i+1])()
    break
  end
  i = i + 1
end
if not config then
  print("config file required (-i <file> option)")
  os.exit(1)
end

local NUM_PRIM_PARTS = 1

if regentlib.config["parallelize"] then
  local dop = regentlib.config['parallelize-dop']
  local NX,NY,NZ = dop:match('^(%d+),(%d+),(%d+)$')
  if NX then
    NUM_PRIM_PARTS = NX * NY * NZ
  else
    NUM_PRIM_PARTS = assert(tonumber(dop))
  end
end

local particlesPerTask = math.ceil(config.num / NUM_PRIM_PARTS)
local density = config.density
local initialTemperature = config.initialTemperature
local diameter_mean = config.diameter_mean

local xBoundaryWidth = (config.xBCLeft == 'periodic') and 0 or 1
local yBoundaryWidth = (config.yBCLeft == 'periodic') and 0 or 1
local zBoundaryWidth = (config.zBCLeft == 'periodic') and 0 or 1
local xInteriorWidth = config.xnum
local yInteriorWidth = config.ynum
local zInteriorWidth = config.znum

-------------------------------------------------------------------------------
-- Local tasks
-------------------------------------------------------------------------------

local __demand(__parallel) task InitParticlesUniform
  (particles : particlesRel:regionType(), cells : cellsRel:regionType())
where
  reads writes(particles),
  reads(cells.{velocity, centerCoordinates})
do
  var pBase = 0
  for p in particles do
    pBase = __raw(p).value
    break
  end
  var lo : int3d = cells.bounds.lo
  lo.x = max(lo.x, xBoundaryWidth)
  lo.y = max(lo.y, yBoundaryWidth)
  lo.z = max(lo.z, zBoundaryWidth)
  var hi : int3d = cells.bounds.hi
  hi.x = min(hi.x, xInteriorWidth + xBoundaryWidth - 1)
  hi.y = min(hi.y, yInteriorWidth + yBoundaryWidth - 1)
  hi.z = min(hi.z, zInteriorWidth + zBoundaryWidth - 1)
  var xSize = hi.x - lo.x + 1
  var ySize = hi.y - lo.y + 1
  --__demand(__openmp)
  for p in particles do
    if __raw(p).value - pBase < particlesPerTask then
      p.__valid = true
      var relIdx = __raw(p).value - pBase
      var c : int3d = { lo.x + relIdx % xSize,
                        lo.y + relIdx / xSize % ySize,
                        lo.z + relIdx / xSize / ySize }
      p.cell = c
      p.position = cells[p.cell].centerCoordinates
      p.particle_velocity = cells[p.cell].velocity
      p.density = density
      p.particle_temperature = initialTemperature
      p.diameter = diameter_mean
    end
  end
end

-------------------------------------------------------------------------------
-- Exported quotes
-------------------------------------------------------------------------------

local exports = {}

exports.InitParticlesUniform = rquote
  InitParticlesUniform([particlesRel:regionSymbol()],
                       [cellsRel:regionSymbol()])
end

-------------------------------------------------------------------------------
-- Module exports
-------------------------------------------------------------------------------

return exports
end
