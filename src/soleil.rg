import "regent"

-------------------------------------------------------------------------------
-- IMPORTS
-------------------------------------------------------------------------------

local C = regentlib.c
local SCHEMA = terralib.includec("config_schema.h")
local UTIL = require 'util'

local ceil = regentlib.ceil(double)

-------------------------------------------------------------------------------
-- DATA STRUCTURES
-------------------------------------------------------------------------------

local Config = SCHEMA.Config
local MultiConfig = SCHEMA.MultiConfig

local struct Particles_columns {
  cell : int3d;
  position : double[3];
  velocity : double[3];
  temperature : double;
  diameter : double;
  density : double;
  deltaVelocityOverRelaxationTime : double[3];
  deltaTemperatureTerm : double;
  position_old : double[3];
  velocity_old : double[3];
  temperature_old : double;
  position_new : double[3];
  velocity_new : double[3];
  temperature_new : double;
  velocity_t : double[3];
  temperature_t : double;
  __valid : bool;
  __xfer_dir : int8;
  __xfer_slot : int;
}

local struct Fluid_columns {
  rho : double;
  pressure : double;
  velocity : double[3];
  centerCoordinates : double[3];
  velocityGradientX : double[3];
  velocityGradientY : double[3];
  velocityGradientZ : double[3];
  temperature : double;
  rhoEnthalpy : double;
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

-------------------------------------------------------------------------------
-- I/O ROUTINES
-------------------------------------------------------------------------------

-- regentlib.rexpr, regentlib.rexpr, regentlib.rexpr* -> regentlib.rquote
local function emitConsoleWrite(config, format, ...)
  local args = terralib.newlist{...}
  return rquote
    var consoleFile = [&int8](C.malloc(256))
    C.snprintf(consoleFile, 256, '%s/console.txt', config.Mapping.outDir)
    var console = UTIL.openFile(consoleFile, 'a')
    C.free(consoleFile)
    C.fprintf(console, format, [args])
    C.fflush(console)
    C.fclose(console)
  end
end

-- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task Console_WriteHeader(config : Config)
  [emitConsoleWrite(config, 'Iter\t'..
                            'Sim Time\t'..
                            'Wall t\t'..
                            'Delta Time\t'..
                            'Avg Press\t'..
                            'Avg Temp\t'..
                            'Average KE\t'..
                            '#Part\t'..
                            'Particle T\n')];
end

-- regentlib.rexpr, regentlib.rexpr, regentlib.rexpr, regentlib.rexpr*
--   -> regentlib.rquote
local function emitProbeWrite(config, probeId, format, ...)
  local args = terralib.newlist{...}
  return rquote
    var filename = [&int8](C.malloc(256))
    C.snprintf(filename, 256, '%s/probe%d.csv', config.Mapping.outDir, probeId)
    var file = UTIL.openFile(filename, 'a')
    C.free(filename)
    C.fprintf(file, format, [args])
    C.fflush(file)
    C.fclose(file)
  end
end

-- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task Probe_WriteHeader(config : Config,
                       probeId : int)
  [emitProbeWrite(config, probeId, 'Iter\t'..
                                   'AvgFluidT\t'..
                                   'AvgParticleT\t'..
                                   'AvgCellOfParticleT\n')];
end

-- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task Probe_Write(config : Config,
                 probeId : int,
                 Integrator_timeStep : int,
                 avgFluidT : double,
                 avgParticleT : double,
                 avgCellOfParticleT : double)
  [emitProbeWrite(config, probeId, '%d\t'..
                                   '%e\t'..
                                   '%e\t'..
                                   '%e\n',
                  Integrator_timeStep,
                  avgFluidT,
                  avgParticleT,
                  avgCellOfParticleT)];
end

-------------------------------------------------------------------------------
-- MAIN SIMULATION
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
  local Particles = regentlib.newsymbol()
  local tiles = regentlib.newsymbol()
  local p_Fluid = regentlib.newsymbol()
  local p_Particles = regentlib.newsymbol()

  -----------------------------------------------------------------------------
  -- Exported symbols
  -----------------------------------------------------------------------------

  INSTANCE.Grid = Grid
  INSTANCE.Integrator_deltaTime = Integrator_deltaTime
  INSTANCE.Integrator_exitCond = Integrator_exitCond
  INSTANCE.Fluid = Fluid
  INSTANCE.Particles = Particles
  INSTANCE.tiles = tiles
  INSTANCE.p_Fluid = p_Fluid
  INSTANCE.p_Particles = p_Particles

  -----------------------------------------------------------------------------
  -- Symbol declaration & initialization
  -----------------------------------------------------------------------------

  function INSTANCE.DeclSymbols(config) return rquote

    ---------------------------------------------------------------------------
    -- Preparation
    ---------------------------------------------------------------------------

    -- Start timer
    var [startTime] = C.legion_get_current_time_in_micros() / 1000;

    -- Write console header
    Console_WriteHeader(config)

    -- Write probe file headers
    for i = 0,config.IO.probes.length do
      Probe_WriteHeader(config, i)
    end

    ---------------------------------------------------------------------------
    -- Declare & initialize state variables
    ---------------------------------------------------------------------------

    -- Cell step size (TODO: Change when we go to non-uniform meshes)
    var [Grid.xCellWidth] = config.Grid.xWidth / config.Grid.xNum
    var [Grid.yCellWidth] = config.Grid.yWidth / config.Grid.yNum
    var [Grid.zCellWidth] = config.Grid.zWidth / config.Grid.zNum
    var [Grid.cellVolume] = Grid.xCellWidth * Grid.yCellWidth * Grid.zCellWidth

    -- Determine number of ghost cells in each direction
    -- 0 ghost cells if periodic and 1 otherwise
    var [Grid.xBnum] = 1
    var [Grid.yBnum] = 1
    var [Grid.zBnum] = 1
    if config.BC.xBCLeft == SCHEMA.FlowBC_Periodic then Grid.xBnum = 0 end
    if config.BC.yBCLeft == SCHEMA.FlowBC_Periodic then Grid.yBnum = 0 end
    if config.BC.zBCLeft == SCHEMA.FlowBC_Periodic then Grid.zBnum = 0 end

    -- Compute real origin, accounting for ghost cells
    var [Grid.xRealOrigin] = (config.Grid.origin[0]-(Grid.xCellWidth*Grid.xBnum))
    var [Grid.yRealOrigin] = (config.Grid.origin[1]-(Grid.yCellWidth*Grid.yBnum))
    var [Grid.zRealOrigin] = (config.Grid.origin[2]-(Grid.zCellWidth*Grid.zBnum))

    var [NX] = config.Mapping.tiles[0]
    var [NY] = config.Mapping.tiles[1]
    var [NZ] = config.Mapping.tiles[2]
    var [numTiles] = NX * NY * NZ

    var [Integrator_exitCond] = true
    var [Integrator_simTime] = 0.0
    var [Integrator_timeStep] = 0
    var [Integrator_deltaTime] = 0.0

    var [Particles_number] = int64(0)

    ---------------------------------------------------------------------------
    -- Create Regions and Partitions
    ---------------------------------------------------------------------------

    var sampleId = config.Mapping.sampleId

    -- Create Fluid Regions
    var is_Fluid = ispace(int3d, {x = config.Grid.xNum + 2*Grid.xBnum,
                                  y = config.Grid.yNum + 2*Grid.yBnum,
                                  z = config.Grid.zNum + 2*Grid.zBnum})
    var [Fluid] = region(is_Fluid, Fluid_columns);

    -- Create Particles Regions
    regentlib.assert(config.Particles.maxNum % numTiles == 0,
                     'Uneven partitioning of particles')
    var maxParticlesPerTile = ceil((config.Particles.maxNum / numTiles) * config.Particles.maxSkew)
    var is_Particles = ispace(int1d, maxParticlesPerTile * numTiles)
    var [Particles] = region(is_Particles, Particles_columns);

    -- Partitioning domain
    var [tiles] = ispace(int3d, {NX,NY,NZ})

    -- Fluid Partitioning
    var [p_Fluid] =
      [UTIL.mkPartitionEqually(int3d, int3d, Fluid_columns)]
      (Fluid, tiles, Grid.xBnum, Grid.yBnum, Grid.zBnum)

    -- Particles Partitioning
    var [p_Particles] =
      [UTIL.mkPartitionEqually(int1d, int3d, Particles_columns)]
      (Particles, tiles, 0)

  end end -- DeclSymbols

return INSTANCE end -- mkInstance

-------------------------------------------------------------------------------
-- TOP-LEVEL INTERFACE
-------------------------------------------------------------------------------

local SIM = mkInstance()

__forbid(__optimize) __demand(__inner)
task workSingle(config : Config)
  [SIM.DeclSymbols(config)];
end

terra initSample(config : &Config, num : int, outDirBase : &int8)
  config.Mapping.sampleId = num
  C.snprintf(config.Mapping.outDir, 256, "%s/sample%d", outDirBase, num)
  UTIL.createDir(config.Mapping.outDir)
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
      initSample([&Config](config), launched, outDirBase)
      launched += 1
      workSingle(config[0])
    end
  end
  if launched < 1 then
    var stderr = C.fdopen(2, 'w')
    C.fprintf(stderr, "No testcases supplied.\n")
    C.fflush(stderr)
    C.exit(1)
  end
end

-------------------------------------------------------------------------------
-- COMPILATION CALL
-------------------------------------------------------------------------------

regentlib.saveobj(main, "soleil.o", "object")
