-- This is a Lua config file for the Soleil code.

return {

  xnum = 1024, -- number of cells in the x-direction
  ynum = 512, -- number of cells in the y-direction
  znum = 512, -- number of cells in the z-direction

  -- if you increase the cell number and the calculation diverges
  -- right away, decrease the time step on the next line

  delta_time = 1.0e-4,

  -- Control max iterations here. Set to a very high number if
  -- you want to run the full calculation (it will stop once
  -- it hits the 20.0 second max time set below)

  max_iter = 20,

  -- Output options. All output is off by default, but we
  -- will need to turn it on to check results/make visualizations

  consoleFrequency = 5,  -- Iterations between console output of statistics
  wrtRestart = 'OFF',
  wrtVolumeSolution = 'OFF',
  wrt1DSlice = 'OFF',
  wrtParticleEvolution = 'OFF',
  particleEvolutionIndex = 0,
  outputEveryTimeSteps  = 50,
  restartEveryTimeSteps = 50,

  -------------------------------------------
  --[ SHOULD NOT NEED TO MODIFY BELOW HERE]--
  -------------------------------------------

  -- Flow Initialization  Options --
  initCase     = 'TaylorGreen3DVortex', -- Uniform, Restart, TaylorGreen2DVortex, TaylorGreen3DVortex
  initParams = {1.0,100.0,2.0,0.0,0.0}, -- for TGV: first three are density, pressure, velocity
  bodyForce = {0,0.0,0}, -- body force in x, y, z
  turbForcing = 'OFF',          -- Turn turbulent forcing on or off
  turbForceCoeff = 0.0,         -- Turbulent linear forcing coefficient (f = A*rho*u)
  restartIter = 10000,

  -- Grid Options -- PERIODICITY
  origin = {0.0, 0.0, 0.0}, -- spatial origin of the computational domain
  xWidth = 6.283185307179586,
  yWidth = 6.283185307179586,
  zWidth = 6.283185307179586,
  -- BCs on each boundary: 'periodic,' 'symmetry,' or 'wall'
  xBCLeft  = 'periodic',
  xBCLeftVel = {0.0, 0.0, 0.0},
  xBCLeftTemp = 0.0,
  xBCRight = 'periodic',
  xBCRightVel = {0.0, 0.0, 0.0},
  xBCRightTemp = 0.0,
  yBCLeft  = 'periodic',
  yBCLeftVel = {0.0, 0.0, 0.0},
  yBCLeftTemp = 0.0,
  yBCRight = 'periodic',
  yBCRightVel = {0.0, 0.0, 0.0},
  yBCRightTemp = 0.0,
  zBCLeft  = 'periodic',
  zBCLeftVel = {0.0, 0.0, 0.0},
  zBCLeftTemp = 0.0,
  zBCRight = 'periodic',
  zBCRightVel = {0.0, 0.0, 0.0},
  zBCRightTemp = 0.0,

  --Time Integration Options --
  cfl                   = 1.0, -- Negative CFL implies that we will used fixed delta T
  final_time            = 20.00001,

  --- File Output Options --
  headerFrequency       = 200000,
  outputFormat = 'Tecplot', --Tecplot or Python

  -- Fluid Options --
  gasConstant = 20.4128,
  gamma = 1.4,
  viscosity_model = 'PowerLaw', -- Constant, PowerLaw, Sutherland
  constant_visc = 0.004491,          -- Value for a constant viscosity [kg/m/s]
  powerlaw_visc_ref = 0.00044,
  powerlaw_temp_ref = 1.0,
  prandtl = 0.7,
  suth_visc_ref = 1.716E-5,     -- Sutherland's Law reference viscosity [kg/m/s]
  suth_temp_ref = 273.15,       -- Sutherland's Law reference temperature [K]
  suth_s_ref = 110.4,           -- Sutherland's Law S constant [K]

  -- Particle Options --
  -- completely disable particles, including all data
  modeParticles = 'ON',
  initParticles = 'Uniform', -- 'Random' or 'Restart'
  restartParticleIter = 0,
  particleType = 'Free', -- Fixed or Free
  twoWayCoupling = 'OFF', -- ON or OFF
  num = 16000000,
  restitutionCoefficient = 1.0,
  convectiveCoefficient = 0.7, -- W m^-2 K^-1
  heatCapacity = 0.7, -- J Kg^-1 K^-1
  initialTemperature = 20, -- K
  density = 8900, -- kg/m^3
  diameter_mean = 5e-3, -- m
  diameter_maxDeviation = 1e-3, -- m, for statistical distribution
  bodyForceParticles = {0.0,0.0,0.0},
  absorptivity = 1.0, -- Equal to emissivity in thermal equilibrium
  maximum_num = 16000000, -- upper bound on particles with insertion
  insertion_rate = 0, -- per face and per time step
  insertion_mode = {0,0,0,0,0,0}, --bool, MinX MaxX MinY MaxY MinZ MaxZ
  deletion_mode = {0,0,0,0,0,0}, --bool, MinX MaxX MinY MaxY MinZ MaxZ
  -- (Kirchhoff law of thermal radiation)

  -- Radiation Options --
  radiationType = 'Algebraic', -- ON or OFF
  radiationIntensity = 1e3,
  zeroAvgHeatSource = 'OFF'

}
