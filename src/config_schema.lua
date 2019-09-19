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
Exports.GridType = Enum('Uniform','Stretched')
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
    meanVelocity = Array(3,double),
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
    -- number of cells in the radiation grid
    xNum = int,
    yNum = int,
    zNum = int,
    -- number of quadrature points
    angles = int,
    -- wall emissivity [0.0-1.0]
    xHiEmiss = double,
    xLoEmiss = double,
    yHiEmiss = double,
    yLoEmiss = double,
    zHiEmiss = double,
    zLoEmiss = double,
    -- wall blackbody temperature [K]
    xHiTemp = double,
    xLoTemp = double,
    yHiTemp = double,
    yLoTemp = double,
    zHiTemp = double,
    zLoTemp = double,
    -- incoming wall intensity [W/m^2]
    -- power per unit of particle area (as projected on the wall)
    -- assumed monochromatic and collimated
    -- only applied over the quadrature point that is normal to the wall
    xHiIntensity = double,
    xLoIntensity = double,
    yHiIntensity = double,
    yLoIntensity = double,
    zHiIntensity = double,
    zLoIntensity = double,
    -- illuminated window on each wall
    xHiWindow = Exports.Window,
    xLoWindow = Exports.Window,
    yHiWindow = Exports.Window,
    yLoWindow = Exports.Window,
    zHiWindow = Exports.Window,
    zLoWindow = Exports.Window,
  },
}

-- Main config struct
Exports.Config = {
  Mapping = {
    -- number of tiles in which to split the domain
    tiles = Array(3,int),
    -- number of tiles to allocate to each rank
    tilesPerRank = Array(3,int),
    -- unique id assigned to each sample, according to its order in the command
    -- line (first sample is 0, second is 1 etc.); the initial value of this
    -- option is irrelevant, it will be overriden by the code
    sampleId = int,
    -- output directory for each sample; the initial value of this option is
    -- irrelevant, it will be overriden by the code
    outDir = String(256),
    -- expected wall-clock execution time [minutes]
    wallTime = int,
  },
  Grid = {
    -- number of cells in the fluid grid
    xNum = int,
    yNum = int,
    zNum = int,
    -- coordinates of the fluid grid's origin [m]
    origin = Array(3,double),
    -- width of the fluid grid [m]
    xWidth = double,
    yWidth = double,
    zWidth = double,
    -- grid type in each direction
    xType = Exports.GridType,
    yType = Exports.GridType,
    zType = Exports.GridType,
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
    startIter = int,
    startTime = double,
    maxIter = int,
    cfl = double,
    fixedDeltaTime = double,
    -- what order RK method to use [2-4]
    rkOrder = int,
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
    initParams = Array(6,double),
    bodyForce = Array(3,double),
    turbForcing = Exports.TurbForcingModel,
  },
  Particles = {
    initCase = Exports.ParticlesInitCase,
    restartDir = String(256),
    initNum = int64,
    maxNum = int64,
    restitutionCoeff = double,
    convectiveCoeff = double,
    heatCapacity = double,
    initTemperature = double,
    density = double,
    diameterMean = double,
    bodyForce = Array(3,double),
    maxSkew = double,
    escapeRatioPerDir = double,
    collisions = bool,
    twoWayCoupled = bool,
    feeding = Exports.FeedModel,
    -- how many timesteps to advance the fluid before every particle solve
    staggerFactor = int,
    parcelSize = int,
  },
  Radiation = Exports.RadiationModel,
  IO = {
    -- whether to write restart files (requires compiling with HDF support)
    wrtRestart = bool,
    -- how often to write restart files
    restartEveryTimeSteps = int,
    -- temperature probes
    probes = UpTo(5, Exports.Volume),
  },
}

-- Dual-section simulation config
Exports.MultiConfig = {
  -- case configurations for the two sections
  configs = Array(2,Exports.Config),
  -- volume to copy from every timestep (in the 1st section)
  copySrc = Exports.Volume,
  -- volume to copy into every timestep (in the 2nd section)
  copyTgt = Exports.Volume,
  -- whether to place the tiles of the two sections on the same set of ranks
  collocateSections = bool,
  -- How often to copy values from one section to the other
  copyEveryTimeSteps = int,
}

return Exports
