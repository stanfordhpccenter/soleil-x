local Exports = {}

-- Enumeration constants
Exports.FlowBC = {
  Periodic = 0,
  Symmetry = 1,
  AdiabaticWall = 2,
  IsothermalWall = 3,
}
Exports.ParticleBC = {
  Permeable = 0,
  Solid = 1,
}
Exports.ViscosityModel = {
  Constant = 0,
  PowerLaw = 1,
  Sutherland = 2,
}
Exports.FlowInitCase = {
  Uniform = 0,
  Restart = 1,
  Perturbed = 2,
  TaylorGreen2DVortex = 3,
  TaylorGreen3DVortex = 4,
}
Exports.OnOrOff = {
  OFF = 0,
  ON = 1,
}
Exports.ParticlesInitCase = {
  Random = 0,
  Restart = 1,
  Uniform = 2,
}
Exports.RadiationType = {
  OFF = 0,
  Algebraic = 1,
  DOM = 2,
}

-- Main config struct
Exports.Config = {
  Mapping = {
    xTiles = int,
    yTiles = int,
    zTiles = int,
    sampleId = int,
  },
  Grid = {
    xNum = int,
    yNum = int,
    zNum = int,
    origin = double[3],
    xWidth = double,
    yWidth = double,
    zWidth = double,
  },
  BC = {
    xBCLeft = Exports.FlowBC,
    xBCLeftVel = double[3],
    xBCLeftTemp = double,
    xBCRight = Exports.FlowBC,
    xBCRightVel = double[3],
    xBCRightTemp = double,
    yBCLeft = Exports.FlowBC,
    yBCLeftVel = double[3],
    yBCLeftTemp = double,
    yBCRight = Exports.FlowBC,
    yBCRightVel = double[3],
    yBCRightTemp = double,
    zBCLeft = Exports.FlowBC,
    zBCLeftVel = double[3],
    zBCLeftTemp = double,
    zBCRight = Exports.FlowBC,
    zBCRightVel = double[3],
    zBCRightTemp = double,
  },
  Integrator = {
    finalTime = double,
    restartIter = int,
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
    initParams = double[5],
    bodyForce = double[3],
    turbForcing = Exports.OnOrOff,
  },
  Particles = {
    initCase = Exports.ParticlesInitCase,
    initNum = int,
    maxNum = int,
    restitutionCoeff = double,
    convectiveCoeff = double,
    absorptivity = double,
    heatCapacity = double,
    initTemperature = double,
    density = double,
    diameterMean = double,
    bodyForce = double[3],
    maxSkew = double,
    maxXferNum = int,
  },
  Radiation = {
    type = Exports.RadiationType,
    intensity = double,
    qa = double,
    qs = double,
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
    wrtRestart = Exports.OnOrOff,
    restartEveryTimeSteps = int,
    consoleFrequency = int,
    headerFrequency = int,
  },
}

return Exports
