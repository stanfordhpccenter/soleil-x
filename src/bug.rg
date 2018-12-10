import "regent"

-------------------------------------------------------------------------------
-- IMPORTS
-------------------------------------------------------------------------------

local C = regentlib.c
local MAPPER = terralib.includec("soleil_mapper.h")
local SCHEMA = terralib.includec("config_schema.h")
local UTIL = require 'util'

local MAX_ANGLES_PER_QUAD = 44


-------------------------------------------------------------------------------
-- DATA STRUCTURES
-------------------------------------------------------------------------------

local Config = SCHEMA.Config

local struct Fluid_columns {
  rho : double;
  pressure : double;
  velocity : double[3];
  centerCoordinates : double[3];
  velocityGradientX : double[3];
  velocityGradientY : double[3];
  velocityGradientZ : double[3];
  temperature : double;
  rhoVelocity : double[3];
  rhoEnergy : double;
  rho_old : double;
  rhoVelocity_old : double[3];
  rhoEnergy_old : double;
  rho_new : double;
  rhoVelocity_new : double[3];
  rhoEnergy_new : double;
  rho_t : double;
  rhoVelocity_t : double[3];
  rhoEnergy_t : double;
  rhoFluxX : double;
  rhoVelocityFluxX : double[3];
  rhoEnergyFluxX : double;
  rhoFluxY : double;
  rhoVelocityFluxY : double[3];
  rhoEnergyFluxY : double;
  rhoFluxZ : double;
  rhoVelocityFluxZ : double[3];
  rhoEnergyFluxZ : double;
  dissipation : double;
  dissipationFlux : double;
  to_Radiation : int3d;
  dudtBoundary : double;
  dTdtBoundary : double;
  velocity_old_NSCBC : double[3];
  temperature_old_NSCBC : double;
  velocity_inc : double[3];
  temperature_inc : double;
}

local Fluid_primitives = terralib.newlist({
  'rho',
  'pressure',
  'velocity',
  'temperature',
})



-------------------------------------------------------------------------------
-- EXTERNAL MODULE IMPORTS
-------------------------------------------------------------------------------


-------------------------------------------------------------------------------
-- CONSTANTS
-------------------------------------------------------------------------------

local PI = 3.1415926535898
local SB = 5.67e-08


-------------------------------------------------------------------------------
-- MACROS
-------------------------------------------------------------------------------



-------------------------------------------------------------------------------
-- OTHER ROUTINES
-------------------------------------------------------------------------------


__demand(__parallel, __cuda)
task Flow_InitializeCell(Fluid : region(ispace(int3d), Fluid_columns))
where
  writes(Fluid.centerCoordinates),
  writes(Fluid.dissipation),
  writes(Fluid.dissipationFlux),
  writes(Fluid.pressure),
  writes(Fluid.rho),
  writes(Fluid.rhoEnergy),
  writes(Fluid.{rhoEnergyFluxX, rhoEnergyFluxY, rhoEnergyFluxZ}),
  writes(Fluid.rhoEnergy_new),
  writes(Fluid.rhoEnergy_old),
  writes(Fluid.rhoEnergy_t),
  writes(Fluid.rhoFluxX),
  writes(Fluid.rhoFluxY),
  writes(Fluid.rhoFluxZ),
  writes(Fluid.rhoVelocity),
  writes(Fluid.{rhoVelocityFluxX, rhoVelocityFluxY, rhoVelocityFluxZ}),
  writes(Fluid.rhoVelocity_new),
  writes(Fluid.rhoVelocity_old),
  writes(Fluid.rhoVelocity_t),
  writes(Fluid.rho_new),
  writes(Fluid.rho_old),
  writes(Fluid.rho_t),
  writes(Fluid.temperature),
  writes(Fluid.velocity),
  writes(Fluid.{velocityGradientX, velocityGradientY, velocityGradientZ}),
  writes(Fluid.dudtBoundary),
  writes(Fluid.dTdtBoundary),
  writes(Fluid.velocity_old_NSCBC),
  writes(Fluid.temperature_old_NSCBC),
  writes(Fluid.velocity_inc),
  writes(Fluid.temperature_inc)
do
  __demand(__openmp)
  for c in Fluid do
    Fluid[c].rho = 0.0
    Fluid[c].pressure = 0.0
    Fluid[c].velocity = array(0.0, 0.0, 0.0)
    Fluid[c].centerCoordinates = array(0.0, 0.0, 0.0)
    Fluid[c].velocityGradientX = array(0.0, 0.0, 0.0)
    Fluid[c].velocityGradientY = array(0.0, 0.0, 0.0)
    Fluid[c].velocityGradientZ = array(0.0, 0.0, 0.0)
    Fluid[c].temperature = 0.0
    Fluid[c].rhoVelocity = array(0.0, 0.0, 0.0)
    Fluid[c].rhoEnergy = 0.0
    Fluid[c].rho_old = 0.0
    Fluid[c].rhoVelocity_old = array(0.0, 0.0, 0.0)
    Fluid[c].rhoEnergy_old = 0.0
    Fluid[c].rho_new = 0.0
    Fluid[c].rhoVelocity_new = array(0.0, 0.0, 0.0)
    Fluid[c].rhoEnergy_new = 0.0
    Fluid[c].rho_t = 0.0
    Fluid[c].rhoVelocity_t = array(0.0, 0.0, 0.0)
    Fluid[c].rhoEnergy_t = 0.0
    Fluid[c].rhoFluxX = 0.0
    Fluid[c].rhoVelocityFluxX = array(0.0, 0.0, 0.0)
    Fluid[c].rhoEnergyFluxX = 0.0
    Fluid[c].rhoFluxY = 0.0
    Fluid[c].rhoVelocityFluxY = array(0.0, 0.0, 0.0)
    Fluid[c].rhoEnergyFluxY = 0.0
    Fluid[c].rhoFluxZ = 0.0
    Fluid[c].rhoVelocityFluxZ = array(0.0, 0.0, 0.0)
    Fluid[c].rhoEnergyFluxZ = 0.0
    Fluid[c].dissipation = 0.0
    Fluid[c].dissipationFlux = 0.0
    Fluid[c].dudtBoundary = 0.0
    Fluid[c].dTdtBoundary = 0.0
    Fluid[c].velocity_old_NSCBC = array(0.0, 0.0, 0.0)
    Fluid[c].temperature_old_NSCBC = 0.0
    Fluid[c].velocity_inc = array(0.0, 0.0, 0.0)
    Fluid[c].temperature_inc = 0.0
  end
end

__demand(__parallel, __cuda)
task Flow_InitializeCenterCoordinates(Fluid : region(ispace(int3d), Fluid_columns),
                                      Grid_xBnum : int32, Grid_xNum : int32, Grid_xOrigin : double, Grid_xWidth : double,
                                      Grid_yBnum : int32, Grid_yNum : int32, Grid_yOrigin : double, Grid_yWidth : double,
                                      Grid_zBnum : int32, Grid_zNum : int32, Grid_zOrigin : double, Grid_zWidth : double)
where
  writes(Fluid.centerCoordinates)
do
C.printf("Grid_xBnum %d\n", Grid_xBnum);
-- __demand(__openmp)
for c in Fluid do
var centerCoordinates : double[3];
centerCoordinates[0] = Grid_xOrigin + (Grid_xWidth/Grid_xNum) * (c.x-Grid_xBnum+0.5)
centerCoordinates[1] = Grid_yOrigin + (Grid_yWidth/Grid_yNum) * (c.y-Grid_yBnum+0.5)
centerCoordinates[2] = Grid_zOrigin + (Grid_zWidth/Grid_zNum) * (c.z-Grid_zBnum+0.5)
Fluid[c].centerCoordinates = centerCoordinates
C.printf("Grid_xOrigin %lg + (Grid_xWidth %lg / Grid_xNum %d) * (c.x %lg - Grid_xBnum %d + 0.5) = %lg\n",
Grid_xOrigin, Grid_xWidth, Grid_xNum, c.x, Grid_xBnum, centerCoordinates[0])
C.printf("Grid_xBnum %d c.x %lg Grid_xNum %d Grid_xWidth %lg Grid_xOrigin %lg\n",
Grid_xBnum, c.x, Grid_xNum, Grid_xWidth, Grid_xOrigin)
C.printf("--- Grid_xBnum %d result %lg\n", Grid_xBnum, centerCoordinates[0]);
C.printf("--- cxyz %lg %lg %lg\n", c.x, c.y, c.z)
  end
end




-------------------------------------------------------------------------------
-- TOP-LEVEL INTERFACE
-------------------------------------------------------------------------------
local function mkInstance() local INSTANCE = {}

-----------------------------------------------------------------------------
-- Symbols shared between quotes
-----------------------------------------------------------------------------

local startTime = regentlib.newsymbol()
local Grid = {
xCellWidth = regentlib.newsymbol(),
yCellWidth = regentlib.newsymbol(),
zCellWidth = regentlib.newsymbol(),
cellVolume = regentlib.newsymbol(),
xBnum = regentlib.newsymbol(),
yBnum = regentlib.newsymbol(),
zBnum = regentlib.newsymbol(),
xRealOrigin = regentlib.newsymbol(),
yRealOrigin = regentlib.newsymbol(),
zRealOrigin = regentlib.newsymbol(),
}
local BC = {
xPosSign = regentlib.newsymbol(double[3]),
xNegSign = regentlib.newsymbol(double[3]),
xPosVelocity = regentlib.newsymbol(double[3]),
xNegVelocity = regentlib.newsymbol(double[3]),
xPosTemperature = regentlib.newsymbol(double),
xNegTemperature = regentlib.newsymbol(double),
yPosSign = regentlib.newsymbol(double[3]),
yNegSign = regentlib.newsymbol(double[3]),
yPosVelocity = regentlib.newsymbol(double[3]),
yNegVelocity = regentlib.newsymbol(double[3]),
yPosTemperature = regentlib.newsymbol(double),
yNegTemperature = regentlib.newsymbol(double),
zPosSign = regentlib.newsymbol(double[3]),
zNegSign = regentlib.newsymbol(double[3]),
zPosVelocity = regentlib.newsymbol(double[3]),
zNegVelocity = regentlib.newsymbol(double[3]),
zPosTemperature = regentlib.newsymbol(double),
zNegTemperature = regentlib.newsymbol(double),
xBCParticles = regentlib.newsymbol(SCHEMA.ParticlesBC),
yBCParticles = regentlib.newsymbol(SCHEMA.ParticlesBC),
zBCParticles = regentlib.newsymbol(SCHEMA.ParticlesBC),
}
local NX = regentlib.newsymbol()
local NY = regentlib.newsymbol()
local NZ = regentlib.newsymbol()
local numTiles = regentlib.newsymbol()

local Integrator_deltaTime = regentlib.newsymbol()
local Integrator_simTime = regentlib.newsymbol()
local Integrator_timeStep = regentlib.newsymbol()
local Integrator_exitCond = regentlib.newsymbol()
local Particles_number = regentlib.newsymbol()

local Fluid = regentlib.newsymbol()
local Fluid_copy = regentlib.newsymbol()
local Particles = regentlib.newsymbol()
local Particles_copy = regentlib.newsymbol()
local TradeQueue = UTIL.generate(26, regentlib.newsymbol)
local Radiation = regentlib.newsymbol()
local tiles = regentlib.newsymbol()
local p_Fluid = regentlib.newsymbol()
local p_Fluid_copy = regentlib.newsymbol()
local p_Particles = regentlib.newsymbol()
local p_Particles_copy = regentlib.newsymbol()
local p_TradeQueue = UTIL.generate(26, regentlib.newsymbol)
local p_Radiation = regentlib.newsymbol()


-----------------------------------------------------------------------------
-- Exported symbols
-----------------------------------------------------------------------------

INSTANCE.Grid = Grid
INSTANCE.Fluid = Fluid
INSTANCE.tiles = tiles
INSTANCE.p_Fluid = p_Fluid

-----------------------------------------------------------------------------
-- Symbol declaration & initialization
-----------------------------------------------------------------------------

function INSTANCE.DeclSymbols(config) return rquote

---------------------------------------------------------------------------
-- Declare & initialize state variables
---------------------------------------------------------------------------

-- Cell step size (TODO: Change when we go to non-uniform meshes)
var [Grid.xCellWidth] = config.Grid.xWidth / config.Grid.xNum
var [Grid.yCellWidth] = config.Grid.yWidth / config.Grid.yNum
var [Grid.zCellWidth] = config.Grid.zWidth / config.Grid.zNum
var [Grid.cellVolume] = Grid.xCellWidth * Grid.yCellWidth * Grid.zCellWidth
var [Grid.xBnum] = 1
var [Grid.yBnum] = 1
var [Grid.zBnum] = 1
var [NX] = config.Mapping.tiles[0]
var [NY] = config.Mapping.tiles[1]
var [NZ] = config.Mapping.tiles[2]
var [numTiles] = NX * NY * NZ

---------------------------------------------------------------------------
-- Create Regions and Partitions
---------------------------------------------------------------------------

-- Create Fluid Regions
var is_Fluid = ispace(int3d, {x = config.Grid.xNum + 2*Grid.xBnum,
y = config.Grid.yNum + 2*Grid.yBnum,
z = config.Grid.zNum + 2*Grid.zBnum})
var [Fluid] = region(is_Fluid, Fluid_columns);
-- Partitioning domain
var [tiles] = ispace(int3d, {NX,NY,NZ})

-- Fluid Partitioning
var [p_Fluid] =
[UTIL.mkPartitionEqually(int3d, int3d, Fluid_columns)]
(Fluid, tiles, Grid.xBnum, Grid.yBnum, Grid.zBnum)


end end -- DeclSymbols

-----------------------------------------------------------------------------
-- Region initialization
-----------------------------------------------------------------------------

function INSTANCE.InitRegions(config) return rquote

Flow_InitializeCell(Fluid)
Flow_InitializeCenterCoordinates(Fluid,
Grid.xBnum, config.Grid.xNum, config.Grid.origin[0], config.Grid.xWidth,
Grid.yBnum, config.Grid.yNum, config.Grid.origin[1], config.Grid.yWidth,
Grid.zBnum, config.Grid.zNum, config.Grid.origin[2], config.Grid.zWidth)


end end -- InitRegions

return INSTANCE end -- mkInstance






local SIM = mkInstance()

__forbid(__optimize) 
task workSingle(config : Config)
  [SIM.DeclSymbols(config)];
  [SIM.InitRegions(config)];
end

__demand(__inner)
task main()
var args = regentlib.c.legion_runtime_get_input_args()
var outDirBase = '.'
for i = 1, args.argc do
if C.strcmp(args.argv[i], '-o') == 0 and i < args.argc-1 then
outDirBase = args.argv[i+1]
end
end
var launched = 0
for i = 1, args.argc do
if C.strcmp(args.argv[i], '-i') == 0 and i < args.argc-1 then
var config : Config[1]
SCHEMA.parse_Config([&Config](config), args.argv[i+1])
workSingle(config[0])
end
end
end

-------------------------------------------------------------------------------
-- COMPILATION CALL
-------------------------------------------------------------------------------

regentlib.saveobj(main, "bug.o", "object")
