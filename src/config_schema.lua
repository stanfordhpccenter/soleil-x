local Exports = {}

-- Enumeration constants
Exports.FlowBC = Enum('Periodic','Symmetry','AdiabaticWall','IsothermalWall','NSCBC_SubsonicInflow','NSCBC_SubsonicOutflow')
Exports.ParticleBC = Enum('Permeable','Solid')
Exports.ViscosityModel = Enum('Constant','PowerLaw','Sutherland')
Exports.FlowInitCase = Enum('Uniform','Restart','Perturbed','TaylorGreen2DVortex','TaylorGreen3DVortex')
Exports.ParticlesInitCase = Enum('Random','Restart','Uniform')
Exports.RadiationType = Enum('OFF','Algebraic','DOM')
Exports.PertubationModel = Union{
  OFF = {},
  Random = {
    fromCell = Array(3,int),
    toCell = Array(3,int),
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
    xBCLeftVel = Array(3,double),
    xBCLeftTemp = double,
    xBCRight = Exports.FlowBC,
    xBCRightVel = Array(3,double),
    xBCRightTemp = double,
    -- Pressure that the sub-sonic outlet relaxes
    xBCRightP_inf = double,
    yBCLeft = Exports.FlowBC,
    yBCLeftVel = Array(3,double),
    yBCLeftTemp = double,
    yBCRight = Exports.FlowBC,
    yBCRightVel = Array(3,double),
    yBCRightTemp = double,
    zBCLeft = Exports.FlowBC,
    zBCLeftVel = Array(3,double),
    zBCLeftTemp = double,
    zBCRight = Exports.FlowBC,
    zBCRightVel = Array(3,double),
    zBCRightTemp = double,
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
