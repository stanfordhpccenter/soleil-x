local Exports = {}

-- Unions & enumeration constants
Exports.FlowBC = Enum('Periodic','Symmetry','AdiabaticWall','IsothermalWall','NSCBC_SubsonicInflow','NSCBC_SubsonicOutflow','NonUniformTemperatureWall')
Exports.ParticlesBC = Enum('Periodic','Bounce','Disappear')
Exports.ViscosityModel = Enum('Constant','PowerLaw','Sutherland')
Exports.FlowInitCase = Enum('Uniform','Random','Restart','Perturbed','TaylorGreen2DVortex','TaylorGreen3DVortex')
Exports.ParticlesInitCase = Enum('Random','Restart','Uniform')
Exports.RadiationType = Enum('OFF','Algebraic','DOM')
Exports.WallHeatModel = Union{
  Constant = {
    temperature = double,
  },
  Parabola = {
    T_left = double,
    T_right = double,
    T_mid = double,
  },
}
Exports.InflowProfile = Union{
  Constant = {
    velocity = double,
  },
  DuctProfile = {
    meanVelocity = double,
  },
}
Exports.PertubationModel = Union{
  OFF = {},
  Random = {
    fromCell = Array(3,int),
    uptoCell = Array(3,int),
  },
}
local Window = {
  fromCell = Array(2,int),
  uptoCell = Array(2,int),
}

-- Main config struct
Exports.Config = {
  Mapping = {
    -- number of tiles in which to split the domain, on the x,y,z dimensions
    tiles = Array(3,int),
    -- number of tiles to allocate to each rank, on the x,y,z dimensions
    tilesPerRank = Array(3,int),
    -- unique id assigned to each sample, according to its order in the command
    -- line (first sample is 0, second is 1 etc.); the initial value of this
    -- option is irrelevant, it will be overriden by the code
    sampleId = int,
    -- output directory for each sample; the initial value of this option is
    -- irrelevant, it will be overriden by the code
    outDir = String(256),
    -- expected wall-clock execution time, in minutes
    wallTime = int,
  },
  Grid = {
    -- number of cells in the fluid grid, on the x,y,z dimensions
    xNum = int,
    yNum = int,
    zNum = int,
    -- coordinates of the fluid grid's origin, in meters
    origin = Array(3,double),
    -- width of the fluid grid, on the x,y,z dimensions, in meters
    xWidth = double,
    yWidth = double,
    zWidth = double,
  },
  BC = {
    xBCLeft = Exports.FlowBC,
    xBCLeftVel = Array(3,double),
    xBCLeftHeat = Exports.WallHeatModel,
    xBCLeftInflowProfile = Exports.InflowProfile,
    xBCRight = Exports.FlowBC,
    xBCRightVel = Array(3,double),
    xBCRightHeat = Exports.WallHeatModel,
    -- Pressure that the sub-sonic outlet relaxes
    xBCRightP_inf = double,
    yBCLeft = Exports.FlowBC,
    yBCLeftVel = Array(3,double),
    yBCLeftHeat = Exports.WallHeatModel,
    yBCRight = Exports.FlowBC,
    yBCRightVel = Array(3,double),
    yBCRightHeat = Exports.WallHeatModel,
    zBCLeft = Exports.FlowBC,
    zBCLeftVel = Array(3,double),
    zBCLeftHeat = Exports.WallHeatModel,
    zBCRight = Exports.FlowBC,
    zBCRightVel = Array(3,double),
    zBCRightHeat = Exports.WallHeatModel,
  },
  Integrator = {
    finalTime = double,
    restartIter = int,
    restartTime = double,
    maxIter = int,
    cfl = double,
    fixedDeltaTime = double,
  },
  Flow = {
    gasConstant = double,
    gamma = double,
    prandtl = double,
    viscosityModel = Exports.ViscosityModel,
    constantVisc = double,
    powerlawViscRef = double,
    powerlawTempRef = double,
    sutherlandViscRef = double,
    sutherlandTempRef = double,
    sutherlandSRef = double,
    initCase = Exports.FlowInitCase,
    restartDir = String(256),
    initParams = Array(5,double),
    bodyForce = Array(3,double),
    turbForcing = bool,
    pertubation = Exports.PertubationModel,
  },
  Particles = {
    initCase = Exports.ParticlesInitCase,
    restartDir = String(256),
    initNum = int,
    maxNum = int,
    restitutionCoeff = double,
    convectiveCoeff = double,
    absorptivity = double,
    heatCapacity = double,
    initTemperature = double,
    density = double,
    diameterMean = double,
    bodyForce = Array(3,double),
    maxSkew = double,
    maxXferNum = int,
    collisions = bool,
  },
  Radiation = {
    type = Exports.RadiationType,
    intensity = double,
    qa = double,
    qs = double,
    -- number of cells in the radiation grid, on the x,y,z dimensions
    xNum = int,
    yNum = int,
    zNum = int,
    xHiEmiss = double,
    xLoEmiss = double,
    yHiEmiss = double,
    yLoEmiss = double,
    zHiEmiss = double,
    zLoEmiss = double,
    xHiTemp = double,
    xLoTemp = double,
    yHiTemp = double,
    yLoTemp = double,
    zHiTemp = double,
    zLoTemp = double,
    xHiWindow = Window,
    xLoWindow = Window,
    yHiWindow = Window,
    yLoWindow = Window,
    zHiWindow = Window,
    zLoWindow = Window,
    angles = int,
  },
  IO = {
    -- whether to write restart files (requires compiling with HDF support)
    wrtRestart = bool,
    -- how often to write restart files
    restartEveryTimeSteps = int,
    -- Temperature probes
    probes = UpTo(5,{
      coords = Array(3,int),
      frequency = int,
    }),
  },
}

return Exports
