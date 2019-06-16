local Exports = {}

-- Helper definitions
Exports.Volume = Struct{
  fromCell = Array(3,Int(4)),
  uptoCell = Array(3,Int(4)),
}
Exports.Window = Struct{
  fromCell = Array(2,Int(4)),
  uptoCell = Array(2,Int(4)),
}

-- Unions & enumeration constants
Exports.FlowBC = Enum('Periodic','Symmetry','AdiabaticWall','IsothermalWall','NSCBC_SubsonicInflow','NSCBC_SubsonicOutflow','NonUniformTemperatureWall')
Exports.ParticlesBC = Enum('Periodic','Bounce','Disappear')
Exports.ViscosityModel = Enum('Constant','PowerLaw','Sutherland')
Exports.FlowInitCase = Enum('Uniform','Random','Restart','Perturbed','TaylorGreen2DVortex','TaylorGreen3DVortex')
Exports.ParticlesInitCase = Enum('Random','Restart','Uniform')
Exports.GridType = Enum('Uniform','Stretched')
Exports.TempProfile = Union{
  Constant = Struct{
    temperature = Float(8),
  },
  Parabola = Struct{
    T_left = Float(8),
    T_right = Float(8),
    T_mid = Float(8),
  },
  Incoming = Struct{},
}
Exports.InflowProfile = Union{
  Constant = Struct{
    velocity = Float(8),
  },
  Duct = Struct{
    meanVelocity = Float(8),
  },
  Incoming = Struct{
    addedVelocity = Float(8),
  },
}
Exports.TurbForcingModel = Union{
  OFF = Struct{},
  HIT = Struct{
    meanVelocity = Array(3,Float(8)),
    G = Float(8),
    t_o = Float(8),
    K_o = Float(8),
  },
}
Exports.FeedModel = Union{
  OFF = Struct{},
  Incoming = Struct{
    addedVelocity = Array(3,Float(8)),
  },
}
Exports.RadiationModel = Union{
  OFF = Struct{},
  Algebraic = Struct{
    intensity = Float(8),
    absorptivity = Float(8),
  },
  DOM = Struct{
    qa = Float(8),
    qs = Float(8),
    -- number of cells in the radiation grid
    xNum = Int(4),
    yNum = Int(4),
    zNum = Int(4),
    -- number of quadrature points
    angles = Int(4),
    -- wall emissivity [0.0-1.0]
    xHiEmiss = Float(8),
    xLoEmiss = Float(8),
    yHiEmiss = Float(8),
    yLoEmiss = Float(8),
    zHiEmiss = Float(8),
    zLoEmiss = Float(8),
    -- wall blackbody temperature [K]
    xHiTemp = Float(8),
    xLoTemp = Float(8),
    yHiTemp = Float(8),
    yLoTemp = Float(8),
    zHiTemp = Float(8),
    zLoTemp = Float(8),
    -- incoming wall intensity [W/m^2]
    -- power per unit of particle area (as projected on the wall)
    -- assumed monochromatic and collimated
    -- only applied over the quadrature point that is normal to the wall
    xHiIntensity = Float(8),
    xLoIntensity = Float(8),
    yHiIntensity = Float(8),
    yLoIntensity = Float(8),
    zHiIntensity = Float(8),
    zLoIntensity = Float(8),
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
Exports.Config = Struct{
  Mapping = Struct{
    -- number of tiles in which to split the domain
    tiles = Array(3,Int(4)),
    -- number of tiles to allocate to each rank
    tilesPerRank = Array(3,Int(4)),
    -- unique id assigned to each sample, according to its order in the command
    -- line (first sample is 0, second is 1 etc.); the initial value of this
    -- option is irrelevant, it will be overriden by the code
    sampleId = Int(4),
    -- output directory for each sample; the initial value of this option is
    -- irrelevant, it will be overriden by the code
    outDir = String(256),
    -- expected wall-clock execution time [minutes]
    wallTime = Int(4),
  },
  Grid = Struct{
    -- number of cells in the fluid grid
    xNum = Int(4),
    yNum = Int(4),
    zNum = Int(4),
    -- coordinates of the fluid grid's origin [m]
    origin = Array(3,Float(8)),
    -- width of the fluid grid [m]
    xWidth = Float(8),
    yWidth = Float(8),
    zWidth = Float(8),
    -- grid type in each direction
    xType = Exports.GridType,
    yType = Exports.GridType,
    zType = Exports.GridType,
  },
  BC = Struct{
    xBCLeft = Exports.FlowBC,
    xBCLeftVel = Array(3,Float(8)),
    xBCLeftHeat = Exports.TempProfile,
    xBCLeftInflowProfile = Exports.InflowProfile,
    xBCRight = Exports.FlowBC,
    xBCRightVel = Array(3,Float(8)),
    xBCRightHeat = Exports.TempProfile,
    -- Pressure that the sub-sonic outlet relaxes
    xBCRightP_inf = Float(8),
    yBCLeft = Exports.FlowBC,
    yBCLeftVel = Array(3,Float(8)),
    yBCLeftHeat = Exports.TempProfile,
    yBCRight = Exports.FlowBC,
    yBCRightVel = Array(3,Float(8)),
    yBCRightHeat = Exports.TempProfile,
    zBCLeft = Exports.FlowBC,
    zBCLeftVel = Array(3,Float(8)),
    zBCLeftHeat = Exports.TempProfile,
    zBCRight = Exports.FlowBC,
    zBCRightVel = Array(3,Float(8)),
    zBCRightHeat = Exports.TempProfile,
  },
  Integrator = Struct{
    startIter = Int(4),
    startTime = Float(8),
    maxIter = Int(4),
    cfl = Float(8),
    fixedDeltaTime = Float(8),
    -- what order RK method to use [2-4]
    rkOrder = Int(4),
  },
  Flow = Struct{
    gasConstant = Float(8),
    gamma = Float(8),
    prandtl = Float(8),
    viscosityModel = Exports.ViscosityModel,
    constantVisc = Float(8),
    powerlawViscRef = Float(8),
    powerlawTempRef = Float(8),
    sutherlandViscRef = Float(8),
    sutherlandTempRef = Float(8),
    sutherlandSRef = Float(8),
    initCase = Exports.FlowInitCase,
    restartDir = String(256),
    initParams = Array(6,Float(8)),
    bodyForce = Array(3,Float(8)),
    turbForcing = Exports.TurbForcingModel,
  },
  Particles = Struct{
    initCase = Exports.ParticlesInitCase,
    restartDir = String(256),
    initNum = Int(8),
    maxNum = Int(8),
    restitutionCoeff = Float(8),
    convectiveCoeff = Float(8),
    heatCapacity = Float(8),
    initTemperature = Float(8),
    density = Float(8),
    diameterMean = Float(8),
    bodyForce = Array(3,Float(8)),
    maxSkew = Float(8),
    escapeRatioPerDir = Float(8),
    collisions = Bool(),
    feeding = Exports.FeedModel,
    -- how many timesteps to advance the fluid before every particle solve
    staggerFactor = Int(4),
    parcelSize = Int(4),
  },
  Radiation = Exports.RadiationModel,
  IO = Struct{
    -- whether to write restart files (requires compiling with HDF support)
    wrtRestart = Bool(),
    -- how often to write restart files
    restartEveryTimeSteps = Int(4),
    -- temperature probes
    probes = UpTo(5, Exports.Volume),
  },
}

-- Dual-section simulation config
Exports.MultiConfig = Struct{
  -- case configurations for the two sections
  configs = Array(2,Exports.Config),
  -- volume to copy from every timestep (in the 1st section)
  copySrc = Exports.Volume,
  -- volume to copy into every timestep (in the 2nd section)
  copyTgt = Exports.Volume,
  -- whether to place the tiles of the two sections on the same set of ranks
  collocateSections = Bool(),
  -- How often to copy values from one section to the other
  copyEveryTimeSteps = Int(4),
}

return Exports
