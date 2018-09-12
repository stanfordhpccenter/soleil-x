local Exports = {}

-- Helper definitions
Exports.Volume = {
  fromCell = Array(3,int),
  uptoCell = Array(3,int),
}
Exports.Window = {
  fromCell = Array(2,int),
  uptoCell = Array(2,int),
}

-- Unions & enumeration constants
Exports.FlowBC = Enum('Periodic','Symmetry','AdiabaticWall','IsothermalWall','NSCBC_SubsonicInflow','NSCBC_SubsonicOutflow','NonUniformTemperatureWall')
Exports.ParticlesBC = Enum('Periodic','Bounce','Disappear')
Exports.ViscosityModel = Enum('Constant','PowerLaw','Sutherland')
Exports.FlowInitCase = Enum('Uniform','Random','Restart','Perturbed','TaylorGreen2DVortex','TaylorGreen3DVortex')
Exports.ParticlesInitCase = Enum('Random','Restart','Uniform')
Exports.TempProfile = Union{
  Constant = {
    temperature = double,
  },
  Parabola = {
    T_left = double,
    T_right = double,
    T_mid = double,
  },
  Incoming = {},
}
Exports.InflowProfile = Union{
  Constant = {
    velocity = double,
  },
  Duct = {
    meanVelocity = double,
  },
  Incoming = {
    addedVelocity = double,
  },
}
Exports.TurbForcingModel = Union{
  OFF = {},
  HIT = {
    G = double,
    t_o = double,
    K_o = double,
  },
}
Exports.FeedModel = Union{
  OFF = {},
  Incoming = {
    addedVelocity = Array(3,double),
  },
}
Exports.RadiationModel = Union{
  OFF = {},
  Algebraic = {
    intensity = double,
    absorptivity = double,
  },
  DOM = {
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
    xHiWindow = Exports.Window,
    xLoWindow = Exports.Window,
    yHiWindow = Exports.Window,
    yLoWindow = Exports.Window,
    zHiWindow = Exports.Window,
    zLoWindow = Exports.Window,
    angles = int,
  },
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
    xBCLeftHeat = Exports.TempProfile,
    xBCLeftInflowProfile = Exports.InflowProfile,
    xBCRight = Exports.FlowBC,
    xBCRightVel = Array(3,double),
    xBCRightHeat = Exports.TempProfile,
    -- Pressure that the sub-sonic outlet relaxes
    xBCRightP_inf = double,
    yBCLeft = Exports.FlowBC,
    yBCLeftVel = Array(3,double),
    yBCLeftHeat = Exports.TempProfile,
    yBCRight = Exports.FlowBC,
    yBCRightVel = Array(3,double),
    yBCRightHeat = Exports.TempProfile,
    zBCLeft = Exports.FlowBC,
    zBCLeftVel = Array(3,double),
    zBCLeftHeat = Exports.TempProfile,
    zBCRight = Exports.FlowBC,
    zBCRightVel = Array(3,double),
    zBCRightHeat = Exports.TempProfile,
  },
  Integrator = {
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
    turbForcing = Exports.TurbForcingModel,
  },
  Particles = {
    initCase = Exports.ParticlesInitCase,
    restartDir = String(256),
    initNum = int,
    maxNum = int,
    restitutionCoeff = double,
    convectiveCoeff = double,
    heatCapacity = double,
    initTemperature = double,
    density = double,
    diameterMean = double,
    bodyForce = Array(3,double),
    maxSkew = double,
    maxXferNum = int,
    collisions = bool,
    feeding = Exports.FeedModel,
  },
  Radiation = Exports.RadiationModel,
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

-- Dual-section simulation config
Exports.MultiConfig = {
  -- Case configurations for the two sections
  configs = Array(2,Exports.Config),
  -- Volume to copy from every timestep (in the 1st section)
  copySrc = Exports.Volume,
  -- Volume to copy into every timestep (in the 2nd section)
  copyTgt = Exports.Volume,
}

return Exports
