local Exports = {}

-- Unions & enumeration constants
Exports.FlowBC = Union{
  Periodic = {},
  Symmetry = {},
  AdiabaticWall = {
    velocity = Array(3,double),
  },
  IsothermalWall = {
    velocity = Array(3,double),
    temperature = double,
  },
}
Exports.ViscosityModel = Union{
  Constant = {
    viscosity = double,
  },
  PowerLaw = {
    viscRef = double,
    tempRef = double,
  },
  Sutherland = {
    viscRef = double,
    tempRef = double,
    sRef = double,
  },
}
Exports.FlowInitCase = Union{
  Uniform = {
    rho = double,
    pressure = double,
    velocity = Array(3,double),
  },
  Restart = {
    dir = String(256),
    iter = int,
    time = double,
  },
  Perturbed = {
    rho = double,
    pressure = double,
    velocity = Array(3,double),
  },
  TaylorGreen2DVortex = {
    density = double,
    pressure = double,
    velocity = double,
  },
  TaylorGreen3DVortex = {
    density = double,
    pressure = double,
    velocity = double,
  },
}
Exports.ParticlesInitCase = Union{
  Random = {},
  Restart = {
    dir = String(256),
    density = double,
  },
  Uniform = {
    num = int,
    density = double,
    temperature = double,
    diameter = double,
  },
}
Exports.RadiationModel = Union{
  OFF = {},
  Algebraic = {
    intensity = double,
  },
  DOM = {
    qa = double,
    qs = double,
    -- number of cells in the radiation grid, on the x,y,z dimensions
    xNum = int,
    yNum = int,
    zNum = int,
    emissEast = double,
    emissWest = double,
    emissSouth = double,
    emissNorth = double,
    emissUp = double,
    emissDown = double,
    tempEast = double,
    tempWest = double,
    tempSouth = double,
    tempNorth = double,
    tempUp = double,
    tempDown = double,
  },
}

-- Main config struct
Exports.Config = {
  Mapping = {
    -- number of tiles in which to split the domain, on the x,y,z dimensions;
    -- each tile will occupy a different node
    xTiles = int,
    yTiles = int,
    zTiles = int,
    -- unique id assigned to each sample, according to its order in the command
    -- line (first sample is 0, second is 1 etc.); the initial value of this
    -- option is irrelevant, it will be overriden by the code
    sampleId = int,
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
    xBCRight = Exports.FlowBC,
    yBCLeft = Exports.FlowBC,
    yBCRight = Exports.FlowBC,
    zBCLeft = Exports.FlowBC,
    zBCRight = Exports.FlowBC,
  },
  Integrator = {
    finalTime = double,
    maxIter = int,
    cfl = double,
    fixedDeltaTime = double,
  },
  Flow = {
    gasConstant = double,
    gamma = double,
    prandtl = double,
    viscosityModel = Exports.ViscosityModel,
    init = Exports.FlowInitCase,
    bodyForce = Array(3,double),
    turbForcing = bool,
  },
  Particles = {
    init = Exports.ParticlesInitCase,
    maxNum = int,
    restitutionCoeff = double,
    convectiveCoeff = double,
    absorptivity = double,
    heatCapacity = double,
    bodyForce = Array(3,double),
    maxSkew = double,
    maxXferNum = int,
  },
  Radiation = Exports.RadiationModel,
  IO = {
    -- whether to write restart files (requires compiling with HDF support)
    wrtRestart = bool,
    -- how often to write restart files
    restartEveryTimeSteps = int,
    -- how often to write intermediate statistics to the console
    consoleFrequency = int,
    headerFrequency = int,
    probes = UpTo(5,{
      coords = Array(3,int),
      frequency = int,
    }),
  },
}

return Exports
