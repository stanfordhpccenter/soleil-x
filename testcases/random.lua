return {

  -- Flow Initialization  Options
  initCase    = 'Uniform', -- Uniform, Restart, TaylorGreen2DVortex, TaylorGreen3DVortex
  initParams  = {0.4472728, 0.4001904, 0.3713178, 0.6807681, 0.4900506},
  bodyForce   = {0.7751856, 0.7816449, 0.1264908},
  turbForceCoeff = 0.6875105,
  turbForcing = 'OFF',

  restartIter = 111,

  -- Grid Options
  xnum = 124320, -- 148 * 840
  ynum = 320364, -- 396 * 809
  znum = 201584, -- 344 * 586
  -- dop = 148 * 396 * 344
  origin = {0.6987052, 0.9980223, 0.7057003},
  xWidth = 0.5988125,
  yWidth = 0.8466831,
  zWidth = 0.8310359,
  -- BCs: 'periodic,' 'symmetry,' 'adiabatic_wall,' or 'isothermal_wall'
  xBCLeft  = 'periodic',
  xBCLeftVel = {0.2099739, 0.4410183, 0.4906633},
  xBCLeftTemp = 0.4044863,
  xBCRight = 'periodic',
  xBCRightVel = {0.1156091, 0.6859607, 0.4446281},
  xBCRightTemp = 0.4196384,
  yBCLeft  = 'adiabatic_wall',
  yBCLeftVel = {0.9918297, 0.8624091, 0.3291094},
  yBCLeftTemp = 0.8230053,
  yBCRight = 'adiabatic_wall',
  yBCRightVel = {0.4206601, 0.6641698, 0.5674205},
  yBCRightTemp = 0.7810168,
  zBCLeft  = 'symmetry',
  zBCLeftVel = {0.1792673, 0.8027272, 0.3218365},
  zBCLeftTemp = 0.3692598,
  zBCRight = 'symmetry',
  zBCRightVel = {0.9480248, 0.9376553, 0.4133212},
  zBCRightTemp = 0.3741956,

  -- Time Integration Options
  final_time            = 0.8029553,
  max_iter              = 222,
  cfl                   = 0.4106583,

  -- File Output Options
  wrtRestart = 'ON',
  consoleFrequency = 333,
  restartEveryTimeSteps = 444,
  headerFrequency       = 555,

  -- Fluid Options
  gasConstant = 0.7272462,
  gamma = 0.8198141,
  prandtl = 0.1565651,
  viscosity_model = 'Sutherland', -- 'Constant', 'PowerLaw', or 'Sutherland'
  constant_visc = 0.2454194,
  powerlaw_visc_ref = 0.9683056,
  powerlaw_temp_ref = 0.9913615,
  suth_visc_ref = 0.5815515,
  suth_temp_ref = 0.6935851,
  suth_s_ref = 0.2772065,

  -- Particle Options
  modeParticles = 'ON',
  initParticles = 'Uniform', -- 'Random' or 'Restart'
  particleType = 'Free', -- Fixed or Free
  twoWayCoupling = 'ON',
  num = 120,
  maximum_num = 240,
  restitutionCoefficient = 0.1487899,
  convectiveCoefficient = 0.8720309,
  heatCapacity = 0.3767944,
  initialTemperature = 0.7490512,
  density = 0.4769553,
  diameter_mean = 0.5250806,
  diameter_maxDeviation = 0.9746212,
  bodyForceParticles = {0.5899043, 0.3012763, 0.6355317},
  absorptivity = 0.6796266,

  -- Radiation Options
  radiationType = 'Algebraic', -- 'Algebraic', 'DOM', 'MCRT', 'OFF'
  radiationIntensity = 0.258018,
  zeroAvgHeatSource = 'OFF',
}
