import "regent"

-------------------------------------------------------------------------------
-- IMPORTS
-------------------------------------------------------------------------------

local C = terralib.includecstring[[
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
]]
local MAPPER = terralib.includec("soleil_mapper.h")
local SCHEMA = terralib.includec("config_schema.h")
local UTIL = require 'util'

local acos = regentlib.acos(double)
local ceil = regentlib.ceil(double)
local cos = regentlib.cos(double)
local exp = regentlib.exp(double)
local fabs = regentlib.fabs(double)
local floor = regentlib.floor(double)
local fmod = regentlib.fmod(double)
local pow = regentlib.pow(double)
local sin = regentlib.sin(double)
local sqrt = regentlib.sqrt(double)
local log = regentlib.log(double)

-------------------------------------------------------------------------------
-- COMPILE-TIME CONFIGURATION
-------------------------------------------------------------------------------

local USE_HDF = assert(os.getenv('USE_HDF')) ~= '0'

local MAX_ANGLES_PER_QUAD = 44

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
  position_ghost : double[3];
  velocity_ghost : double[3];
  velocity_t_ghost : double[3];
  position_t : double[3];
  velocity_t : double[3];
  temperature_t : double;
  __valid : bool;
}

local Particles_primitives = terralib.newlist({
  'position',
  'velocity',
  'temperature',
  'diameter',
  'density',
  '__valid',
})

local Particles_subStepTemp = terralib.newlist({
  'deltaVelocityOverRelaxationTime',
  'deltaTemperatureTerm',
  'position_ghost',
  'velocity_ghost',
  'velocity_t_ghost',
  'position_t',
  'velocity_t',
  'temperature_t',
})

local Particles_subStepConserved =
  UTIL.setToList(
    UTIL.setSubList(
      UTIL.listToSet(UTIL.fieldNames(Particles_columns)),
      Particles_subStepTemp))

local TradeQueue_columns =
  UTIL.deriveStruct('TradeQueue_columns',
                    Particles_columns,
                    Particles_subStepConserved,
                    {__source=int1d, __target=int1d})

local CopyQueue_columns =
  UTIL.deriveStruct('CopyQueue_columns',
                    Particles_columns,
                    Particles_primitives)

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
  convectiveSpectralRadius : double;
  viscousSpectralRadius : double;
  heatConductionSpectralRadius : double;
  rhoVelocity : double[3];
  rhoEnergy : double;
  rhoBoundary : double;
  rhoVelocityBoundary : double[3];
  rhoEnergyBoundary : double;
  velocityBoundary : double[3];
  pressureBoundary : double;
  temperatureBoundary : double;
  velocityGradientXBoundary : double[3];
  velocityGradientYBoundary : double[3];
  velocityGradientZBoundary : double[3];
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
  PD : double;
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

local Fluid_primitives = terralib.newlist({
  'rho',
  'pressure',
  'velocity',
  'temperature',
})

struct Radiation_columns {
  G : double;
  S : double;
  Ib : double;
  sigma : double;
  acc_d2 : double;
  acc_d2t4 : double;
}

-------------------------------------------------------------------------------
-- EXTERNAL MODULE IMPORTS
-------------------------------------------------------------------------------

local DOM = (require 'dom-desugared')(MAX_ANGLES_PER_QUAD, Radiation_columns, SCHEMA)

-------------------------------------------------------------------------------
-- CONSTANTS
-------------------------------------------------------------------------------

local PI = 3.1415926535898
local SB = 5.67e-08

local SIZEOF_PARTICLE = sizeof(Particles_columns)

local terra validFieldOffset()
  var x : Particles_columns
  return [int64]([&int8](&(x.__valid)) - [&int8](&x))
end
local VALID_FIELD_OFFSET = validFieldOffset()

-------------------------------------------------------------------------------
-- MACROS
-------------------------------------------------------------------------------

-- TODO: Define macros for in_boundary, neg_depth etc.

local __demand(__inline)
task drand48_r(rngState : &C.drand48_data)
  var res : double[1]
  C.drand48_r(rngState, [&double](res))
  return res[0]
end

__demand(__inline)
task vs_mul(a : double[3], b : double)
  return array(a[0] * b, a[1] * b, a[2] * b)
end

__demand(__inline)
task vs_div(a : double[3], b : double)
  return array(a[0] / b, a[1] / b, a[2] / b)
end

__demand(__inline)
task dot(a : double[3], b : double[3])
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
end

__demand(__inline)
task vv_add(a : double[3], b : double[3])
  return array(a[0] + b[0], a[1] + b[1], a[2] + b[2])
end

__demand(__inline)
task vv_sub(a : double[3], b : double[3])
  return array(a[0] - b[0], a[1] - b[1], a[2] - b[2])
end

__demand(__inline)
task vv_mul(a : double[3], b : double[3])
  return array(a[0] * b[0], a[1] * b[1], a[2] * b[2])
end

__demand(__inline)
task vv_div(a : double[3], b : double[3])
  return array(a[0] / b[0], a[1] / b[1], a[2] / b[2])
end

-------------------------------------------------------------------------------
-- I/O ROUTINES
-------------------------------------------------------------------------------

local __demand(__inline)
task Fluid_dump(colors : ispace(int3d),
                dirname : &int8,
                Fluid : region(ispace(int3d),Fluid_columns),
                Fluid_copy : region(ispace(int3d),Fluid_columns),
                p_Fluid : partition(disjoint, Fluid, colors),
                p_Fluid_copy : partition(disjoint, Fluid_copy, colors))
where
  reads(Fluid.[Fluid_primitives]),
  reads writes(Fluid_copy.[Fluid_primitives]),
  Fluid * Fluid_copy
do
  regentlib.assert(false, 'Recompile with USE_HDF=1')
end

local __demand(__inline)
task Fluid_load(colors : ispace(int3d),
                dirname : &int8,
                Fluid : region(ispace(int3d),Fluid_columns),
                Fluid_copy : region(ispace(int3d),Fluid_columns),
                p_Fluid : partition(disjoint, Fluid, colors),
                p_Fluid_copy : partition(disjoint, Fluid_copy, colors))
where
  reads writes(Fluid.[Fluid_primitives]),
  reads writes(Fluid_copy.[Fluid_primitives]),
  Fluid * Fluid_copy
do
  regentlib.assert(false, 'Recompile with USE_HDF=1')
end

local __demand(__inline)
task Particles_dump(colors : ispace(int3d),
                    dirname : &int8,
                    Particles : region(ispace(int1d),Particles_columns),
                    Particles_copy : region(ispace(int1d),Particles_columns),
                    p_Particles : partition(disjoint, Particles, colors),
                    p_Particles_copy : partition(disjoint, Particles_copy, colors))
where
  reads(Particles.[Particles_primitives]),
  reads writes(Particles_copy.[Particles_primitives]),
  Particles * Particles_copy
do
  regentlib.assert(false, 'Recompile with USE_HDF=1')
end

local __demand(__inline)
task Particles_load(colors : ispace(int3d),
                    dirname : &int8,
                    Particles : region(ispace(int1d),Particles_columns),
                    Particles_copy : region(ispace(int1d),Particles_columns),
                    p_Particles : partition(disjoint, Particles, colors),
                    p_Particles_copy : partition(disjoint, Particles_copy, colors))
where
  reads writes(Particles.[Particles_primitives]),
  reads writes(Particles_copy.[Particles_primitives]),
  Particles * Particles_copy
do
  regentlib.assert(false, 'Recompile with USE_HDF=1')
end

if USE_HDF then
  local HDF = require "hdf_helper"
  Fluid_dump, Fluid_load = HDF.mkHDFTasks(
    int3d, int3d, Fluid_columns, Fluid_primitives)
  Particles_dump, Particles_load = HDF.mkHDFTasks(
    int1d, int3d, Particles_columns, Particles_primitives)
end

-- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task Probes_write(Fluid : region(ispace(int3d), Fluid_columns),
                  exitCond : bool,
                  Integrator_timeStep : int,
                  config : Config)
where
  reads(Fluid.temperature)
do
  var bounds = Fluid.bounds
  for i = 0,config.IO.probes.length do
    var probe = config.IO.probes.values[i]
    var coords = probe.coords
    if (exitCond or Integrator_timeStep % probe.frequency == 0)
    and bounds.lo.x <= coords[0] and coords[0] <= bounds.hi.x
    and bounds.lo.y <= coords[1] and coords[1] <= bounds.hi.y
    and bounds.lo.z <= coords[2] and coords[2] <= bounds.hi.z then
      var filename = [&int8](C.malloc(256))
      C.snprintf(filename, 256, '%s/probe%d.csv', config.Mapping.outDir, i)
      var file = UTIL.openFile(filename, 'a')
      C.free(filename)
      var temp = Fluid[int3d{coords[0],coords[1],coords[2]}].temperature
      C.fprintf(file, '%d\t%lf\n', Integrator_timeStep, temp)
      C.fclose(file)
    end
  end
end

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

task Console_write(config : Config,
                   Integrator_timeStep : int,
                   Integrator_simTime : double,
                   startTime : uint64,
                   Integrator_deltaTime : double,
                   Flow_averagePressure : double,
                   Flow_averageTemperature : double,
                   Flow_averageKineticEnergy : double,
                   Particles_number : int64,
                   Particles_averageTemperature : double)
  var currTime = regentlib.c.legion_get_current_time_in_micros() / 1000;
  [emitConsoleWrite(config, '%d\t'..
                            '%e\t'..
                            '%llu.%03llu\t'..
                            '%e\t'..
                            '%e\t'..
                            '%e\t'..
                            '%e\t'..
                            '%lld\t'..
                            '%e\n',
                    Integrator_timeStep,
                    Integrator_simTime,
                    rexpr (currTime - startTime) / 1000 end,
                    rexpr (currTime - startTime) % 1000 end,
                    Integrator_deltaTime,
                    Flow_averagePressure,
                    Flow_averageTemperature,
                    Flow_averageKineticEnergy,
                    Particles_number,
                    Particles_averageTemperature)];
end

-------------------------------------------------------------------------------
-- OTHER ROUTINES
-------------------------------------------------------------------------------

__demand(__inline)
task GetDynamicViscosity(temperature : double,
                         Flow_constantVisc : double,
                         Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                         Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                         Flow_viscosityModel : SCHEMA.ViscosityModel)
  var viscosity = 0.0
  if (Flow_viscosityModel == SCHEMA.ViscosityModel_Constant) then
    viscosity = Flow_constantVisc
  else
    if (Flow_viscosityModel == SCHEMA.ViscosityModel_PowerLaw) then
      viscosity = (Flow_powerlawViscRef*pow((temperature/Flow_powerlawTempRef), double(0.75)))
    else
      viscosity = ((Flow_sutherlandViscRef*pow((temperature/Flow_sutherlandTempRef), (3.0/2.0)))*((Flow_sutherlandTempRef+Flow_sutherlandSRef)/(temperature+Flow_sutherlandSRef)))
    end
  end
  return viscosity
end

__demand(__parallel) -- NO CUDA
task InitParticlesUniform(Particles : region(ispace(int1d), Particles_columns),
                          Fluid : region(ispace(int3d), Fluid_columns),
                          config : Config,
                          Grid_xBnum : int32, Grid_yBnum : int32, Grid_zBnum : int32)
where
  reads(Fluid.{centerCoordinates, velocity}),
  reads(Particles.cell),
  writes(Particles.{__valid, cell, position, velocity, density, temperature, diameter})
do
  var pBase = 0
  for p in Particles do -- this loop trips up the CUDA codegen
    pBase = int32(p)
    break
  end
  var lo = Fluid.bounds.lo
  lo.x = max(lo.x, Grid_xBnum)
  lo.y = max(lo.y, Grid_yBnum)
  lo.z = max(lo.z, Grid_zBnum)
  var hi = Fluid.bounds.hi
  hi.x = min(hi.x, ((config.Grid.xNum+Grid_xBnum)-1))
  hi.y = min(hi.y, ((config.Grid.yNum+Grid_yBnum)-1))
  hi.z = min(hi.z, ((config.Grid.zNum+Grid_zBnum)-1))
  var xSize = ((hi.x-lo.x)+1)
  var ySize = ((hi.y-lo.y)+1)
  var zSize = ((hi.z-lo.z)+1)
  var particlesPerTask = config.Particles.initNum / (config.Mapping.tiles[0]*config.Mapping.tiles[1]*config.Mapping.tiles[2])
  var Particles_density = config.Particles.density
  var Particles_initTemperature = config.Particles.initTemperature
  var Particles_diameterMean = config.Particles.diameterMean
  __demand(__openmp)
  for p in Particles do
    if ((int32(p)-pBase)<particlesPerTask) then
      p.__valid = true
      var relIdx = (int32(p)-pBase)
      var c = int3d({(lo.x+(relIdx%xSize)), (lo.y+((relIdx/xSize)%ySize)), (lo.z+((relIdx/xSize)/ySize%zSize))})
      p.cell = c
      p.position = Fluid[p.cell].centerCoordinates
      p.velocity = Fluid[p.cell].velocity
      p.density = Particles_density
      p.temperature = Particles_initTemperature
      p.diameter = Particles_diameterMean
    end
  end
end

__demand(__parallel, __cuda)
task AddRadiation(Particles : region(ispace(int1d), Particles_columns),
                  config : Config)
where
  reads(Particles.{density, diameter, __valid}),
  reads writes(Particles.temperature_t)
do
  var absorptivity = config.Radiation.u.Algebraic.absorptivity
  var intensity = config.Radiation.u.Algebraic.intensity
  var heatCapacity = config.Particles.heatCapacity
  __demand(__openmp)
  for p in Particles do
    if p.__valid then
      var crossSectionArea = PI*pow(p.diameter,2.0)/4.0
      var volume = PI*pow(p.diameter,3.0)/6.0
      var mass = volume*p.density
      var absorbedRadiationIntensity = absorptivity*intensity*crossSectionArea
      p.temperature_t += absorbedRadiationIntensity/(mass*heatCapacity)
    end
  end
end

__demand(__parallel, __cuda)
task Particles_initValidField(Particles : region(ispace(int1d), Particles_columns))
where
  writes(Particles.__valid)
do
  __demand(__openmp)
  for p in Particles do
    p.__valid = false
  end
end

__demand(__parallel, __cuda)
task SetCoarseningField(Fluid : region(ispace(int3d), Fluid_columns),
                        Grid_xBnum : int32, Grid_xNum : int32,
                        Grid_yBnum : int32, Grid_yNum : int32,
                        Grid_zBnum : int32, Grid_zNum : int32,
                        Radiation_xNum : int32,
                        Radiation_yNum : int32,
                        Radiation_zNum : int32)
where
  writes(Fluid.to_Radiation)
do
  var xFactor = (Grid_xNum/Radiation_xNum)
  var yFactor = (Grid_yNum/Radiation_yNum)
  var zFactor = (Grid_zNum/Radiation_zNum)
  __demand(__openmp)
  for f in Fluid do
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(f).x)), 0)>0) or (max(int32((int3d(f).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(f).y)), 0)>0)) or (max(int32((int3d(f).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(f).z)), 0)>0)) or (max(int32((int3d(f).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      Fluid[f].to_Radiation = int3d{(int3d(f).x-Grid_xBnum)/xFactor,
                                    (int3d(f).y-Grid_yBnum)/yFactor,
                                    (int3d(f).z-Grid_zBnum)/zFactor}
    else
      Fluid[f].to_Radiation = int3d({uint64(0ULL), uint64(0ULL), uint64(0ULL)})
    end
  end
end

__demand(__parallel, __cuda)
task Flow_InitializeCell(Fluid : region(ispace(int3d), Fluid_columns))
where
  writes(Fluid.PD),
  writes(Fluid.centerCoordinates),
  writes(Fluid.convectiveSpectralRadius),
  writes(Fluid.dissipation),
  writes(Fluid.dissipationFlux),
  writes(Fluid.heatConductionSpectralRadius),
  writes(Fluid.pressure),
  writes(Fluid.pressureBoundary),
  writes(Fluid.rho),
  writes(Fluid.rhoBoundary),
  writes(Fluid.rhoEnergy),
  writes(Fluid.rhoEnergyBoundary),
  writes(Fluid.{rhoEnergyFluxX, rhoEnergyFluxY, rhoEnergyFluxZ}),
  writes(Fluid.rhoEnergy_new),
  writes(Fluid.rhoEnergy_old),
  writes(Fluid.rhoEnergy_t),
  writes(Fluid.rhoEnthalpy),
  writes(Fluid.rhoFluxX),
  writes(Fluid.rhoFluxY),
  writes(Fluid.rhoFluxZ),
  writes(Fluid.rhoVelocity),
  writes(Fluid.rhoVelocityBoundary),
  writes(Fluid.{rhoVelocityFluxX, rhoVelocityFluxY, rhoVelocityFluxZ}),
  writes(Fluid.rhoVelocity_new),
  writes(Fluid.rhoVelocity_old),
  writes(Fluid.rhoVelocity_t),
  writes(Fluid.rho_new),
  writes(Fluid.rho_old),
  writes(Fluid.rho_t),
  writes(Fluid.temperature),
  writes(Fluid.temperatureBoundary),
  writes(Fluid.velocity),
  writes(Fluid.velocityBoundary),
  writes(Fluid.{velocityGradientX, velocityGradientY, velocityGradientZ}),
  writes(Fluid.{velocityGradientXBoundary, velocityGradientYBoundary, velocityGradientZBoundary}),
  writes(Fluid.viscousSpectralRadius),
  writes(Fluid.dudtBoundary),
  writes(Fluid.dTdtBoundary),
  writes(Fluid.velocity_old_NSCBC),
  writes(Fluid.temperature_old_NSCBC),
  writes(Fluid.velocity_inc),
  writes(Fluid.temperature_inc)
do
  __demand(__openmp)
  for c in Fluid do
    Fluid[c].rho = 0.0
    Fluid[c].pressure = 0.0
    Fluid[c].velocity = array(0.0, 0.0, 0.0)
    Fluid[c].centerCoordinates = array(0.0, 0.0, 0.0)
    Fluid[c].velocityGradientX = array(0.0, 0.0, 0.0)
    Fluid[c].velocityGradientY = array(0.0, 0.0, 0.0)
    Fluid[c].velocityGradientZ = array(0.0, 0.0, 0.0)
    Fluid[c].temperature = 0.0
    Fluid[c].rhoEnthalpy = 0.0
    Fluid[c].convectiveSpectralRadius = 0.0
    Fluid[c].viscousSpectralRadius = 0.0
    Fluid[c].heatConductionSpectralRadius = 0.0
    Fluid[c].rhoVelocity = array(0.0, 0.0, 0.0)
    Fluid[c].rhoEnergy = 0.0
    Fluid[c].rhoBoundary = 0.0
    Fluid[c].rhoVelocityBoundary = array(0.0, 0.0, 0.0)
    Fluid[c].rhoEnergyBoundary = 0.0
    Fluid[c].velocityBoundary = array(0.0, 0.0, 0.0)
    Fluid[c].pressureBoundary = 0.0
    Fluid[c].temperatureBoundary = 0.0
    Fluid[c].velocityGradientXBoundary = array(0.0, 0.0, 0.0)
    Fluid[c].velocityGradientYBoundary = array(0.0, 0.0, 0.0)
    Fluid[c].velocityGradientZBoundary = array(0.0, 0.0, 0.0)
    Fluid[c].rho_old = 0.0
    Fluid[c].rhoVelocity_old = array(0.0, 0.0, 0.0)
    Fluid[c].rhoEnergy_old = 0.0
    Fluid[c].rho_new = 0.0
    Fluid[c].rhoVelocity_new = array(0.0, 0.0, 0.0)
    Fluid[c].rhoEnergy_new = 0.0
    Fluid[c].rho_t = 0.0
    Fluid[c].rhoVelocity_t = array(0.0, 0.0, 0.0)
    Fluid[c].rhoEnergy_t = 0.0
    Fluid[c].rhoFluxX = 0.0
    Fluid[c].rhoVelocityFluxX = array(0.0, 0.0, 0.0)
    Fluid[c].rhoEnergyFluxX = 0.0
    Fluid[c].rhoFluxY = 0.0
    Fluid[c].rhoVelocityFluxY = array(0.0, 0.0, 0.0)
    Fluid[c].rhoEnergyFluxY = 0.0
    Fluid[c].rhoFluxZ = 0.0
    Fluid[c].rhoVelocityFluxZ = array(0.0, 0.0, 0.0)
    Fluid[c].rhoEnergyFluxZ = 0.0
    Fluid[c].PD = 0.0
    Fluid[c].dissipation = 0.0
    Fluid[c].dissipationFlux = 0.0
    Fluid[c].dudtBoundary = 0.0
    Fluid[c].dTdtBoundary = 0.0
    Fluid[c].velocity_old_NSCBC = array(0.0, 0.0, 0.0)
    Fluid[c].temperature_old_NSCBC = 0.0
    Fluid[c].velocity_inc = array(0.0, 0.0, 0.0)
    Fluid[c].temperature_inc = 0.0
  end
end

__demand(__parallel, __cuda)
task Flow_InitializeCenterCoordinates(Fluid : region(ispace(int3d), Fluid_columns),
                                      Grid_xBnum : int32, Grid_xNum : int32, Grid_xOrigin : double, Grid_xWidth : double,
                                      Grid_yBnum : int32, Grid_yNum : int32, Grid_yOrigin : double, Grid_yWidth : double,
                                      Grid_zBnum : int32, Grid_zNum : int32, Grid_zOrigin : double, Grid_zWidth : double)
where
  writes(Fluid.centerCoordinates)
do
  __demand(__openmp)
  for c in Fluid do
    Fluid[c].centerCoordinates = array(Grid_xOrigin + (Grid_xWidth/Grid_xNum) * (c.x-Grid_xBnum+0.5),
                                       Grid_yOrigin + (Grid_yWidth/Grid_yNum) * (c.y-Grid_yBnum+0.5),
                                       Grid_zOrigin + (Grid_zWidth/Grid_zNum) * (c.z-Grid_zBnum+0.5))
  end
end

__demand(__parallel, __cuda)
task Flow_InitializeUniform(Fluid : region(ispace(int3d), Fluid_columns), Flow_initParams : double[5])
where
  writes(Fluid.{rho, pressure, velocity})
do
  __demand(__openmp)
  for c in Fluid do
    Fluid[c].rho = Flow_initParams[0]
    Fluid[c].pressure = Flow_initParams[1]
    Fluid[c].velocity = array(Flow_initParams[2], Flow_initParams[3], Flow_initParams[4])
  end
end

__demand(__parallel) -- NO CUDA, NO OPENMP
task Flow_InitializeRandom(Fluid : region(ispace(int3d), Fluid_columns),
                           Flow_initParams : double[5])
where
  writes(Fluid.{rho, pressure, velocity})
do
  var magnitude = Flow_initParams[2]
  var rngState : C.drand48_data[1]
  var rngStatePtr = [&C.drand48_data](rngState)
  C.srand48_r(regentlib.c.legion_get_current_time_in_nanos(), rngStatePtr)
  for c in Fluid do
    Fluid[c].rho = Flow_initParams[0]
    Fluid[c].pressure = Flow_initParams[1]
    Fluid[c].velocity = array(2 * (drand48_r(rngStatePtr) - 0.5) * magnitude,
                              2 * (drand48_r(rngStatePtr) - 0.5) * magnitude,
                              2 * (drand48_r(rngStatePtr) - 0.5) * magnitude)
  end
end

-- CHANGE do not compute xy instead just pass in cell center since it is computed before this task will be called
__demand(__parallel, __cuda)
task Flow_InitializeTaylorGreen2D(Fluid : region(ispace(int3d), Fluid_columns),
                                  Flow_initParams : double[5],
                                  Grid_xBnum : int32, Grid_xNum : int32, Grid_xOrigin : double, Grid_xWidth : double,
                                  Grid_yBnum : int32, Grid_yNum : int32, Grid_yOrigin : double, Grid_yWidth : double,
                                  Grid_zBnum : int32, Grid_zNum : int32, Grid_zOrigin : double, Grid_zWidth : double)
where
  writes(Fluid.{rho, velocity, pressure})
do
  __demand(__openmp)
  for c in Fluid do
    var taylorGreenDensity = Flow_initParams[0]
    var taylorGreenPressure = Flow_initParams[1]
    var taylorGreenVelocity = Flow_initParams[2]
    var xy = [double[3]](array((Grid_xOrigin+((Grid_xWidth/double(Grid_xNum))*(double((int3d(c).x-uint64(Grid_xBnum)))+0.5))), (Grid_yOrigin+((Grid_yWidth/double(Grid_yNum))*(double((int3d(c).y-uint64(Grid_yBnum)))+0.5))), (Grid_zOrigin+((Grid_zWidth/double(Grid_zNum))*(double((int3d(c).z-uint64(Grid_zBnum)))+0.5)))))
    var coorZ = 0
    Fluid[c].rho = taylorGreenDensity
    Fluid[c].velocity = vs_mul([double[3]](array(((sin(xy[0])*cos(xy[1]))*cos(coorZ)), (((-cos(xy[0]))*sin(xy[1]))*cos(coorZ)), 0.0)), taylorGreenVelocity)
    var factorA = (cos((2.0*double(coorZ)))+2.0)
    var factorB = (cos((2.0*xy[0]))+cos((2.0*xy[1])))
    Fluid[c].pressure = (taylorGreenPressure+((((taylorGreenDensity*pow(taylorGreenVelocity, 2.0))/16.0)*factorA)*factorB))
  end
end

-- CHANGE do not compute xy instead just pass in cell center since it is computed before this task will be called
__demand(__parallel, __cuda)
task Flow_InitializeTaylorGreen3D(Fluid : region(ispace(int3d), Fluid_columns),
                                  Flow_initParams : double[5],
                                  Grid_xBnum : int32, Grid_xNum : int32, Grid_xOrigin : double, Grid_xWidth : double,
                                  Grid_yBnum : int32, Grid_yNum : int32, Grid_yOrigin : double, Grid_yWidth : double,
                                  Grid_zBnum : int32, Grid_zNum : int32, Grid_zOrigin : double, Grid_zWidth : double)
where
  writes(Fluid.{rho, velocity, pressure})
do
  __demand(__openmp)
  for c in Fluid do
    var taylorGreenDensity = Flow_initParams[0]
    var taylorGreenPressure = Flow_initParams[1]
    var taylorGreenVelocity = Flow_initParams[2]
    var xy = [double[3]](array((Grid_xOrigin+((Grid_xWidth/double(Grid_xNum))*(double((int3d(c).x-uint64(Grid_xBnum)))+0.5))), (Grid_yOrigin+((Grid_yWidth/double(Grid_yNum))*(double((int3d(c).y-uint64(Grid_yBnum)))+0.5))), (Grid_zOrigin+((Grid_zWidth/double(Grid_zNum))*(double((int3d(c).z-uint64(Grid_zBnum)))+0.5)))))
    Fluid[c].rho = taylorGreenDensity
    Fluid[c].velocity = vs_mul([double[3]](array(((sin(xy[0])*cos(xy[1]))*cos(xy[2])), (((-cos(xy[0]))*sin(xy[1]))*cos(xy[2])), 0.0)), taylorGreenVelocity)
    var factorA = (cos((2.0*xy[2]))+2.0)
    var factorB = (cos((2.0*xy[0]))+cos((2.0*xy[1])))
    Fluid[c].pressure = (taylorGreenPressure+((((taylorGreenDensity*pow(taylorGreenVelocity, 2.0))/16.0)*factorA)*factorB))
  end
end

__demand(__parallel) -- NO CUDA, NO OPENMP
task Flow_InitializePerturbed(Fluid : region(ispace(int3d), Fluid_columns),
                              Flow_initParams : double[5])
where
  writes(Fluid.{rho, pressure, velocity})
do
  var rngState : C.drand48_data[1]
  var rngStatePtr = [&C.drand48_data](rngState)
  C.srand48_r(regentlib.c.legion_get_current_time_in_nanos(), rngStatePtr)
  for c in Fluid do
    Fluid[c].rho = Flow_initParams[0]
    Fluid[c].pressure = Flow_initParams[1]
    Fluid[c].velocity = array(Flow_initParams[2] + (drand48_r(rngStatePtr)-0.5)*10.0,
                              Flow_initParams[3] + (drand48_r(rngStatePtr)-0.5)*10.0,
                              Flow_initParams[4] + (drand48_r(rngStatePtr)-0.5)*10.0)
  end
end

-- Modify initial rho, pressure, velocity to be consistent with NSCBC inflow outflow condition, also set up stuff for RHS of inflow
__demand(__parallel, __cuda)
task Flow_InitializeGhostNSCBC(Fluid : region(ispace(int3d), Fluid_columns),
                               config : Config,
                               Flow_gasConstant : double,
                               Flow_constantVisc : double,
                               Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                               Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                               Flow_viscosityModel : SCHEMA.ViscosityModel,
                               Grid_xBnum : int32, Grid_xNum : int32,
                               Grid_yBnum : int32, Grid_yNum : int32,
                               Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{rho, velocity, pressure, temperature, centerCoordinates}),
  reads writes(Fluid.{rho, velocity, pressure}),
  writes(Fluid.{velocity_old_NSCBC, temperature_old_NSCBC, dudtBoundary, dTdtBoundary, velocity_inc, temperature_inc})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  -- Domain origin
  var Grid_xOrigin = config.Grid.origin[0]
  var Grid_yOrigin = config.Grid.origin[1]
  var Grid_zOrigin = config.Grid.origin[2]
  -- Domain Size
  var Grid_xWidth = config.Grid.xWidth
  var Grid_yWidth = config.Grid.yWidth
  var Grid_zWidth = config.Grid.zWidth
  -- Cell step size
  var Grid_xCellWidth = (Grid_xWidth/Grid_xNum)
  var Grid_yCellWidth = (Grid_yWidth/Grid_yNum)
  var Grid_zCellWidth = (Grid_zWidth/Grid_zNum)
  -- Compute real origin and width accounting for ghost cells
  var Grid_xRealOrigin = (Grid_xOrigin-(Grid_xCellWidth*Grid_xBnum))
  var Grid_yRealOrigin = (Grid_yOrigin-(Grid_yCellWidth*Grid_yBnum))
  var Grid_zRealOrigin = (Grid_zOrigin-(Grid_zCellWidth*Grid_zBnum))
  -- Inflow values
  var BC_xBCLeftHeat_type = config.BC.xBCLeftHeat.type
  var BC_xBCLeftHeat_Constant_temperature = config.BC.xBCLeftHeat.u.Constant.temperature
  var BC_xBCLeftInflowProfile_type = config.BC.xBCLeftInflowProfile.type
  var BC_xBCLeftInflowProfile_Constant_velocity = config.BC.xBCLeftInflowProfile.u.Constant.velocity
  var BC_xBCLeftInflowProfile_Duct_meanVelocity = config.BC.xBCLeftInflowProfile.u.Duct.meanVelocity
  var BC_xBCLeftInflowProfile_Incoming_addedVelocity = config.BC.xBCLeftInflowProfile.u.Incoming.addedVelocity
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0)
    var xPosGhost  = (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)
    var yNegGhost = (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)
    var yPosGhost  = (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)
    var zNegGhost = (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)
    var zPosGhost  = (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)

    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if  NSCBC_inflow_cell then -- x- boundary cell
      -- Assume subsonic inflow
      var c_bnd = int3d(c)
      var c_int = ((c+{1, 0, 0})%Fluid.bounds)

      var velocity = array(0.0, 0.0, 0.0)
      if BC_xBCLeftInflowProfile_type == SCHEMA.InflowProfile_Constant then
        velocity[0] = BC_xBCLeftInflowProfile_Constant_velocity
      elseif BC_xBCLeftInflowProfile_type == SCHEMA.InflowProfile_Duct then
        var y = Fluid[c].centerCoordinates[1]
        var z = Fluid[c].centerCoordinates[2]
        var y_dist_to_wall = 0.0
        var y_local = 0.0
        if y < (Grid_yWidth/ 2.0) then
          y_dist_to_wall = y
          y_local = (Grid_yWidth/ 2.0) - y
        else
          y_dist_to_wall = Grid_yWidth - y
          y_local = y - (Grid_yWidth/ 2.0)
        end
        var z_dist_to_wall = 0.0
        var z_local = 0.0
        if z < (Grid_zWidth/ 2.0) then
          z_dist_to_wall = z
          z_local = (Grid_zWidth/ 2.0) - z
        else
          z_dist_to_wall = Grid_zWidth - z
          z_local = z - (Grid_zWidth/ 2.0)
        end
        var d = 0.0
        var d_max = 0.0
        if y_dist_to_wall < z_dist_to_wall then
          d = y_dist_to_wall
          d_max = (Grid_yWidth/ 2.0)
        else
          d = z_dist_to_wall
          d_max = (Grid_zWidth/ 2.0)
        end
        var meanVelocity = BC_xBCLeftInflowProfile_Duct_meanVelocity
        var mu = GetDynamicViscosity(Fluid[c].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)
        var Re = Fluid[c].rho*meanVelocity*Grid_yWidth / mu
        var n = -1.7 + 1.8*log(Re)
        velocity[0] = meanVelocity*pow((d/d_max), (1.0/n))
      else -- BC_xBCLeftInflowProfile_type == SCHEMA.InflowProfile_Incoming
        -- This value will be overwritten by the incoming fluid, so just set
        -- it to something reasonable.
        velocity = Fluid[c_int].velocity
        -- HACK: The inflow boundary code later overrides the velocity field
        -- based on the incoming values, so we set a fake incoming value, that
        -- would have produced the fake velocity setting above.
        var velocity_inc = velocity
        velocity_inc[0] -= BC_xBCLeftInflowProfile_Incoming_addedVelocity
        Fluid[c_bnd].velocity_inc = velocity_inc
      end
      Fluid[c_bnd].velocity = velocity

      -- Just copy over the density from the interior
      Fluid[c_bnd].rho = Fluid[c_int].rho

      var temperature : double
      if BC_xBCLeftHeat_type == SCHEMA.TempProfile_Constant then
        temperature = BC_xBCLeftHeat_Constant_temperature

        -- Use the specified temperature to find the correct pressure for current density from EOS
        Fluid[c_bnd].pressure = temperature*Flow_gasConstant*Fluid[c_bnd].rho

      -- elseif BC_xBCLeftHeat_type == SCHEMA.TempProfile_Parabola then
      --   regentlib.assert(false, 'Parabola heat model not supported')
      else -- BC_xBCLeftHeat_type == SCHEMA.TempProfile_Incoming
        -- This value will be overwritten by the incoming fluid, so just set
        -- it to something reasonable.
        Fluid[c_bnd].pressure = Fluid[c_int].pressure
        -- Use equation of state to find temperature of cell
        temperature = (Fluid[c_bnd].pressure/(Flow_gasConstant*Fluid[c_bnd].rho))
        -- HACK: The inflow boundary code later overrides the temperature field
        -- based on the incoming values, so we set a fake incoming value, that
        -- would have produced the fake pressure setting above.
        Fluid[c_bnd].temperature_inc = temperature
      end

      -- for time stepping RHS of INFLOW
      Fluid[c_bnd].velocity_old_NSCBC = velocity
      Fluid[c_bnd].temperature_old_NSCBC = temperature
      Fluid[c_bnd].dudtBoundary = 0.0
      Fluid[c_bnd].dTdtBoundary= 0.0
    end
    if NSCBC_outflow_cell then -- x+ boundary
      -- Assume subsonic outflow
      var c_bnd = int3d(c)
      var c_int = ((c+{-1, 0, 0})%Fluid.bounds)

      -- just copy over the interior variable values
      Fluid[c_bnd].rho       = Fluid[c_int].rho
      Fluid[c_bnd].velocity  = Fluid[c_int].velocity
      Fluid[c_bnd].pressure  = Fluid[c_int].pressure

    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdateConservedFromPrimitive(Fluid : region(ispace(int3d), Fluid_columns),
                                       Flow_gamma : double,
                                       Flow_gasConstant : double,
                                       Grid_xBnum : int32, Grid_xNum : int32,
                                       Grid_yBnum : int32, Grid_yNum : int32,
                                       Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{rho, velocity, pressure}),
  writes(Fluid.{rhoVelocity, rhoEnergy})
do
  __demand(__openmp)
  for c in Fluid do
    -- if interior cell
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      var tmpTemperature = (Fluid[c].pressure/(Flow_gasConstant*Fluid[c].rho))
      var velocity = Fluid[c].velocity
      Fluid[c].rhoVelocity = vs_mul(Fluid[c].velocity, Fluid[c].rho)
      var cv = (Flow_gasConstant/(Flow_gamma-1.0))
      Fluid[c].rhoEnergy = Fluid[c].rho*((cv*tmpTemperature)+(0.5*dot(velocity, velocity)))
    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdateConservedFromPrimitiveGhostNSCBC(Fluid : region(ispace(int3d), Fluid_columns),
                                                 config : Config,
                                                 Flow_gamma : double,
                                                 Flow_gasConstant : double,
                                                 Grid_xBnum : int32, Grid_xNum : int32,
                                                 Grid_yBnum : int32, Grid_yNum : int32,
                                                 Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{rho, velocity, pressure}),
  writes(Fluid.{rhoVelocity, rhoEnergy})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0)
    var xPosGhost  = (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)
    var yNegGhost = (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)
    var yPosGhost  = (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)
    var zNegGhost = (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)
    var zPosGhost  = (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)

    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if (NSCBC_inflow_cell or NSCBC_outflow_cell) then
      var tmpTemperature = (Fluid[c].pressure/(Flow_gasConstant*Fluid[c].rho))
      var velocity = Fluid[c].velocity
      Fluid[c].rhoVelocity = vs_mul(Fluid[c].velocity, Fluid[c].rho)
      var cv = (Flow_gasConstant/(Flow_gamma-1.0))
      Fluid[c].rhoEnergy = Fluid[c].rho*((cv*tmpTemperature)+(0.5*dot(velocity, velocity)))
    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdateAuxiliaryVelocity(Fluid : region(ispace(int3d), Fluid_columns),
                                  Grid_xBnum : int32, Grid_xNum : int32,
                                  Grid_yBnum : int32, Grid_yNum : int32,
                                  Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{rho, rhoVelocity}),
  writes(Fluid.{velocity})
do
  __demand(__openmp)
  for c in Fluid do
    -- If interior cells
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      var velocity = vs_div(Fluid[c].rhoVelocity, Fluid[c].rho)
      Fluid[c].velocity = velocity
    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdateAuxiliaryVelocityGhostNSCBC(Fluid : region(ispace(int3d), Fluid_columns),
                                            config : Config,
                                            Flow_constantVisc : double,
                                            Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                                            Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                                            Flow_viscosityModel : SCHEMA.ViscosityModel,
                                            Grid_xBnum : int32, Grid_xNum : int32,
                                            Grid_yBnum : int32, Grid_yNum : int32,
                                            Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{rho, rhoVelocity, temperature, centerCoordinates, velocity_inc}),
  writes(Fluid.{velocity})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  -- Domain origin
  var Grid_xOrigin = config.Grid.origin[0]
  var Grid_yOrigin = config.Grid.origin[1]
  var Grid_zOrigin = config.Grid.origin[2]
  -- Domain Size
  var Grid_xWidth = config.Grid.xWidth
  var Grid_yWidth = config.Grid.yWidth
  var Grid_zWidth = config.Grid.zWidth
  -- Cell step size
  var Grid_xCellWidth = (Grid_xWidth/Grid_xNum)
  var Grid_yCellWidth = (Grid_yWidth/Grid_yNum)
  var Grid_zCellWidth = (Grid_zWidth/Grid_zNum)
  -- Compute real origin and width accounting for ghost cells
  var Grid_xRealOrigin = (Grid_xOrigin-(Grid_xCellWidth*Grid_xBnum))
  var Grid_yRealOrigin = (Grid_yOrigin-(Grid_yCellWidth*Grid_yBnum))
  var Grid_zRealOrigin = (Grid_zOrigin-(Grid_zCellWidth*Grid_zBnum))
  -- Inflow values
  var BC_xBCLeftInflowProfile_type = config.BC.xBCLeftInflowProfile.type
  var BC_xBCLeftInflowProfile_Constant_velocity = config.BC.xBCLeftInflowProfile.u.Constant.velocity
  var BC_xBCLeftInflowProfile_Duct_meanVelocity = config.BC.xBCLeftInflowProfile.u.Duct.meanVelocity
  var BC_xBCLeftInflowProfile_Incoming_addedVelocity = config.BC.xBCLeftInflowProfile.u.Incoming.addedVelocity
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0)
    var xPosGhost  = (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)
    var yNegGhost = (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)
    var yPosGhost  = (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)
    var zNegGhost = (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)
    var zPosGhost  = (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)

    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if (NSCBC_inflow_cell) then
      var velocity = array(0.0, 0.0, 0.0)
      if BC_xBCLeftInflowProfile_type == SCHEMA.InflowProfile_Constant then
        velocity[0] = BC_xBCLeftInflowProfile_Constant_velocity
      elseif BC_xBCLeftInflowProfile_type == SCHEMA.InflowProfile_Duct then
        var y = Fluid[c].centerCoordinates[1]
        var z = Fluid[c].centerCoordinates[2]
        var y_dist_to_wall = 0.0
        var y_local = 0.0
        if y < (Grid_yWidth/ 2.0) then
          y_dist_to_wall = y
          y_local = (Grid_yWidth/ 2.0) - y
        else
          y_dist_to_wall = Grid_yWidth - y
          y_local = y - (Grid_yWidth/ 2.0)
        end
        var z_dist_to_wall = 0.0
        var z_local = 0.0
        if z < (Grid_zWidth/ 2.0) then
          z_dist_to_wall = z
          z_local = (Grid_zWidth/ 2.0) - z
        else
          z_dist_to_wall = Grid_zWidth - z
          z_local = z - (Grid_zWidth/ 2.0)
        end
        var d = 0.0
        var d_max = 0.0
        if y_dist_to_wall < z_dist_to_wall then
          d = y_dist_to_wall
          d_max = (Grid_yWidth/ 2.0)
        else
          d = z_dist_to_wall
          d_max = (Grid_zWidth/ 2.0)
        end
        var meanVelocity = BC_xBCLeftInflowProfile_Duct_meanVelocity
        var mu = GetDynamicViscosity(Fluid[c].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)
        var Re = Fluid[c].rho*meanVelocity*Grid_yWidth / mu
        var n = -1.7 + 1.8*log(Re)
        velocity[0] = meanVelocity*pow((d/d_max), (1.0/n))
      else -- BC_xBCLeftInflowProfile_type == SCHEMA.InflowProfile_Incoming
        velocity = Fluid[c].velocity_inc
        velocity[0] += BC_xBCLeftInflowProfile_Incoming_addedVelocity
      end
      Fluid[c].velocity = velocity
    end

    if (NSCBC_outflow_cell) then
      var velocity = vs_div(Fluid[c].rhoVelocity, Fluid[c].rho)
      Fluid[c].velocity = velocity
    end
  end
end

-- given a valid rho, v, and T in the interior cells + the specified BC's-> compute the conserved variabs in the ghost cells
__demand(__parallel, __cuda)
task Flow_UpdateGhostConservedStep1(Fluid : region(ispace(int3d), Fluid_columns),
                                    config : Config,
                                    BC_xNegTemperature : double, BC_xNegVelocity : double[3],
                                    BC_xPosTemperature : double, BC_xPosVelocity : double[3],
                                    BC_xNegSign : double[3], BC_xPosSign : double[3],
                                    BC_yNegTemperature : double, BC_yNegVelocity : double[3],
                                    BC_yPosTemperature : double, BC_yPosVelocity : double[3],
                                    BC_yNegSign : double[3], BC_yPosSign : double[3],
                                    BC_zNegTemperature : double, BC_zNegVelocity : double[3],
                                    BC_zPosTemperature : double, BC_zPosVelocity : double[3],
                                    BC_zNegSign : double[3], BC_zPosSign : double[3],
                                    Flow_gamma : double,
                                    Flow_gasConstant : double,
                                    Flow_constantVisc : double,
                                    Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                                    Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                                    Flow_viscosityModel : SCHEMA.ViscosityModel,
                                    Grid_xBnum : int32, Grid_xNum : int32,
                                    Grid_yBnum : int32, Grid_yNum : int32,
                                    Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{rho, pressure, temperature, rhoVelocity, rhoEnergy, centerCoordinates, velocity_inc, temperature_inc}),
  writes(Fluid.{rhoBoundary, rhoEnergyBoundary, rhoVelocityBoundary})
do
  var BC_xBCLeft  = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  var BC_yBCLeft  = config.BC.yBCLeft
  var BC_yBCRight = config.BC.yBCRight
  var BC_zBCLeft  = config.BC.zBCLeft
  var BC_zBCRight = config.BC.zBCRight
  var BC_yBCLeftHeat_T_left  = config.BC.yBCLeftHeat.u.Parabola.T_left
  var BC_yBCLeftHeat_T_mid   = config.BC.yBCLeftHeat.u.Parabola.T_mid
  var BC_yBCLeftHeat_T_right = config.BC.yBCLeftHeat.u.Parabola.T_right
  var BC_yBCRightHeat_T_left  = config.BC.yBCRightHeat.u.Parabola.T_left
  var BC_yBCRightHeat_T_mid   = config.BC.yBCRightHeat.u.Parabola.T_mid
  var BC_yBCRightHeat_T_right = config.BC.yBCRightHeat.u.Parabola.T_right
  var BC_zBCLeftHeat_T_left  = config.BC.zBCLeftHeat.u.Parabola.T_left
  var BC_zBCLeftHeat_T_mid   = config.BC.zBCLeftHeat.u.Parabola.T_mid
  var BC_zBCLeftHeat_T_right = config.BC.zBCLeftHeat.u.Parabola.T_right
  var BC_zBCRightHeat_T_left  = config.BC.zBCRightHeat.u.Parabola.T_left
  var BC_zBCRightHeat_T_mid   = config.BC.zBCRightHeat.u.Parabola.T_mid
  var BC_zBCRightHeat_T_right = config.BC.zBCRightHeat.u.Parabola.T_right
  -- Domain origin
  var Grid_xOrigin = config.Grid.origin[0]
  var Grid_yOrigin = config.Grid.origin[1]
  var Grid_zOrigin = config.Grid.origin[2]
  -- Domain Size
  var Grid_xWidth = config.Grid.xWidth
  var Grid_yWidth = config.Grid.yWidth
  var Grid_zWidth = config.Grid.zWidth
  -- Cell step size
  var Grid_xCellWidth = (Grid_xWidth/Grid_xNum)
  var Grid_yCellWidth = (Grid_yWidth/Grid_yNum)
  var Grid_zCellWidth = (Grid_zWidth/Grid_zNum)
  -- Compute real origin and width accounting for ghost cells
  var Grid_xRealOrigin = (Grid_xOrigin-(Grid_xCellWidth*Grid_xBnum))
  var Grid_yRealOrigin = (Grid_yOrigin-(Grid_yCellWidth*Grid_yBnum))
  var Grid_zRealOrigin = (Grid_zOrigin-(Grid_zCellWidth*Grid_zBnum))
  -- Inflow values
  var BC_xBCLeftHeat_type = config.BC.xBCLeftHeat.type
  var BC_xBCLeftHeat_Constant_temperature = config.BC.xBCLeftHeat.u.Constant.temperature
  var BC_xBCLeftInflowProfile_type = config.BC.xBCLeftInflowProfile.type
  var BC_xBCLeftInflowProfile_Constant_velocity = config.BC.xBCLeftInflowProfile.u.Constant.velocity
  var BC_xBCLeftInflowProfile_Duct_meanVelocity = config.BC.xBCLeftInflowProfile.u.Duct.meanVelocity
  var BC_xBCLeftInflowProfile_Incoming_addedVelocity = config.BC.xBCLeftInflowProfile.u.Incoming.addedVelocity
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0)
    var xPosGhost = (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)
    var yNegGhost = (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)
    var yPosGhost = (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)
    var zNegGhost = (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)
    var zPosGhost = (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)
    var ghost_cell = (xNegGhost or xPosGhost or
                      yNegGhost or yPosGhost or
                      zNegGhost or zPosGhost )
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if xNegGhost then
      var c_bnd = int3d(c)
      var c_int = ((c+{1, 0, 0})%Fluid.bounds)
      var cv = (Flow_gasConstant/(Flow_gamma-1.0))

      if NSCBC_inflow_cell then
        var rho = Fluid[c_bnd].rho
        var velocity = array(0.0, 0.0, 0.0)
        if BC_xBCLeftInflowProfile_type == SCHEMA.InflowProfile_Constant then
          velocity[0] = BC_xBCLeftInflowProfile_Constant_velocity
        elseif BC_xBCLeftInflowProfile_type == SCHEMA.InflowProfile_Duct then
          var y = Fluid[c].centerCoordinates[1]
          var z = Fluid[c].centerCoordinates[2]
          var y_dist_to_wall = 0.0
          var y_local = 0.0
          if y < (Grid_yWidth/ 2.0) then
            y_dist_to_wall = y
            y_local = (Grid_yWidth/ 2.0) - y
          else
            y_dist_to_wall = Grid_yWidth - y
            y_local = y - (Grid_yWidth/ 2.0)
          end
          var z_dist_to_wall = 0.0
          var z_local = 0.0
          if z < (Grid_zWidth/ 2.0) then
            z_dist_to_wall = z
            z_local = (Grid_zWidth/ 2.0) - z
          else
            z_dist_to_wall = Grid_zWidth - z
            z_local = z - (Grid_zWidth/ 2.0)
          end
          var d = 0.0
          var d_max = 0.0
          if y_dist_to_wall < z_dist_to_wall then
            d = y_dist_to_wall
            d_max = (Grid_yWidth/ 2.0)
          else
            d = z_dist_to_wall
            d_max = (Grid_zWidth/ 2.0)
          end
          var meanVelocity = BC_xBCLeftInflowProfile_Duct_meanVelocity
          var mu = GetDynamicViscosity(Fluid[c].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)
          var Re = Fluid[c].rho*meanVelocity*Grid_yWidth / mu
          var n = -1.7 + 1.8*log(Re)
          velocity[0] = meanVelocity*pow((d/d_max), (1.0/n))
        else -- BC_xBCLeftInflowProfile_type == SCHEMA.InflowProfile_Incoming
          velocity = Fluid[c].velocity_inc
          velocity[0] += BC_xBCLeftInflowProfile_Incoming_addedVelocity
        end
        var temperature : double
        if BC_xBCLeftHeat_type == SCHEMA.TempProfile_Constant then
          temperature = BC_xBCLeftHeat_Constant_temperature
        -- elseif BC_xBCLeftHeat_type == SCHEMA.TempProfile_Parabola then
        --   regentlib.assert(false, 'Parabola heat model not supported')
        else -- BC_xBCLeftHeat_type == SCHEMA.TempProfile_Incoming
          temperature = Fluid[c].temperature_inc
        end
        Fluid[c_bnd].rhoBoundary = rho
        Fluid[c_bnd].rhoVelocityBoundary = vs_mul(velocity, rho)
        Fluid[c_bnd].rhoEnergyBoundary = rho*((cv*temperature)+(0.5*dot(velocity, velocity)))
      else
        var sign = BC_xNegSign
        var bnd_velocity = BC_xNegVelocity
        var bnd_temperature = BC_xNegTemperature

        var rho = 0.0
        var temperature = 0.0
        var velocity = array(0.0, 0.0, 0.0)
        var wall_temperature = 0.0

        velocity = vv_add(vv_mul(vs_div(Fluid[c_int].rhoVelocity, Fluid[c_int].rho), sign), bnd_velocity)

        wall_temperature = Fluid[c_int].temperature -- adibatic wall
        if (bnd_temperature>0.0) then -- isothermal wall
          wall_temperature = bnd_temperature
        end
        temperature = (2.0*wall_temperature)-Fluid[c_int].temperature

        rho = Fluid[c_int].pressure/(Flow_gasConstant*temperature)

        Fluid[c_bnd].rhoBoundary = rho
        Fluid[c_bnd].rhoVelocityBoundary = vs_mul(velocity, rho)
        Fluid[c_bnd].rhoEnergyBoundary = rho*((cv*temperature)+(0.5*dot(velocity, velocity)))
      end

    end
    if xPosGhost then
      var c_bnd = int3d(c)
      var c_int = ((c+{-1, 0, 0})%Fluid.bounds)
      var cv = (Flow_gasConstant/(Flow_gamma-1.0))

      if NSCBC_outflow_cell then
        Fluid[c_bnd].rhoBoundary = Fluid[c_bnd].rho
        Fluid[c_bnd].rhoVelocityBoundary = Fluid[c_bnd].rhoVelocity
        Fluid[c_bnd].rhoEnergyBoundary = Fluid[c_bnd].rhoEnergy
      else
        var sign = BC_xPosSign
        var bnd_velocity = BC_xPosVelocity
        var bnd_temperature = BC_xPosTemperature

        var rho = 0.0
        var temperature = 0.0
        var velocity = array(0.0, 0.0, 0.0)
        var wall_temperature = 0.0

        velocity = vv_add(vv_mul(vs_div(Fluid[c_int].rhoVelocity, Fluid[c_int].rho), sign), bnd_velocity)

        wall_temperature = Fluid[c_int].temperature
        if (bnd_temperature>0.0) then
          wall_temperature = bnd_temperature
        end
        temperature = (2.0*wall_temperature)-Fluid[c_int].temperature

        rho = (Fluid[c_int].pressure/(Flow_gasConstant*temperature))

        Fluid[c_bnd].rhoBoundary = rho
        Fluid[c_bnd].rhoVelocityBoundary = vs_mul(velocity, rho)
        Fluid[c_bnd].rhoEnergyBoundary = rho*((cv*temperature)+(0.5*dot(velocity, velocity)))
      end
    end
    if yNegGhost then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, 1, 0})%Fluid.bounds)
      var cv = (Flow_gasConstant/(Flow_gamma-1.0))

      if (BC_yBCLeft == SCHEMA.FlowBC_NonUniformTemperatureWall) then
        var sign = BC_yNegSign
        var bnd_velocity = BC_yNegVelocity -- velocity at face/boundary

        var rho = 0.0
        var temperature = 0.0
        var velocity = array(0.0, 0.0, 0.0)
        var wall_temperature = 0.0

        velocity = vv_add(vv_mul(vs_div(Fluid[c_int].rhoVelocity, Fluid[c_int].rho), sign), bnd_velocity)

        var c_1 = 2.0/(Grid_xWidth*Grid_xWidth)*( (BC_yBCLeftHeat_T_right - BC_yBCLeftHeat_T_left) - 2.0*(BC_yBCLeftHeat_T_mid - BC_yBCLeftHeat_T_left))
        var c_2 = 4.0/(Grid_xWidth)*((BC_yBCLeftHeat_T_mid - BC_yBCLeftHeat_T_left) - 1.0/4.0*(BC_yBCLeftHeat_T_right - BC_yBCLeftHeat_T_left))
        var c_3 = BC_yBCLeftHeat_T_left
        wall_temperature = c_1*Fluid[c_bnd].centerCoordinates[0]*Fluid[c_bnd].centerCoordinates[0] + c_2*Fluid[c_bnd].centerCoordinates[0] + c_3
        if wall_temperature < 0.0 then --unphysical.... set wall themperature to zero
          wall_temperature = 0.0
        end
        temperature = (2.0*wall_temperature)-Fluid[c_int].temperature

        rho = Fluid[c_int].pressure/(Flow_gasConstant*temperature)

        Fluid[c_bnd].rhoBoundary = rho
        Fluid[c_bnd].rhoVelocityBoundary = vs_mul(velocity, rho)
        Fluid[c_bnd].rhoEnergyBoundary = rho*((cv*temperature)+(0.5*dot(velocity, velocity)))
      else
        var sign = BC_yNegSign
        var bnd_velocity = BC_yNegVelocity
        var bnd_temperature = BC_yNegTemperature
        var rho = 0.0
        var temp_wall = 0.0
        var temperature = 0.0
        var velocity = array(0.0, 0.0, 0.0)
        var velocity__3527 = array(0.0, 0.0, 0.0)
        velocity__3527 = vv_add(vv_mul(vs_div(Fluid[c_int].rhoVelocity, Fluid[c_int].rho), sign), bnd_velocity)
        temp_wall = Fluid[c_int].temperature
        if (bnd_temperature>0.0) then
          temp_wall = bnd_temperature
        end
        temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
        rho = (Fluid[c_int].pressure/(Flow_gasConstant*temperature))
        Fluid[c_bnd].rhoBoundary = rho
        Fluid[c_bnd].rhoVelocityBoundary = vs_mul(velocity__3527, rho)
        Fluid[c_bnd].rhoEnergyBoundary = (rho*((cv*temperature)+(0.5*dot(velocity__3527, velocity__3527))))
      end

    end
    if yPosGhost then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, -1, 0})%Fluid.bounds)
      var cv = (Flow_gasConstant/(Flow_gamma-1.0))

      if (BC_yBCRight == SCHEMA.FlowBC_NonUniformTemperatureWall) then
        var sign = BC_yPosSign
        var bnd_velocity = BC_yPosVelocity

        var rho = 0.0
        var temperature = 0.0
        var velocity = array(0.0, 0.0, 0.0)
        var wall_temperature = 0.0

        velocity = vv_add(vv_mul(vs_div(Fluid[c_int].rhoVelocity, Fluid[c_int].rho), sign), bnd_velocity)

        var c_1 = 2.0/(Grid_xWidth*Grid_xWidth)*( (BC_yBCRightHeat_T_right - BC_yBCRightHeat_T_left) - 2.0*(BC_yBCRightHeat_T_mid - BC_yBCRightHeat_T_left))
        var c_2 = 4.0/(Grid_xWidth)*((BC_yBCRightHeat_T_mid - BC_yBCRightHeat_T_left) - 1.0/4.0*(BC_yBCRightHeat_T_right - BC_yBCRightHeat_T_left))
        var c_3 = BC_yBCRightHeat_T_left
        wall_temperature = c_1*Fluid[c_bnd].centerCoordinates[0]*Fluid[c_bnd].centerCoordinates[0] + c_2*Fluid[c_bnd].centerCoordinates[0] + c_3
        if wall_temperature < 0.0 then --unphysical.... set wall themperature to zero
          wall_temperature = 0.0
        end
        temperature = ((2.0*wall_temperature)-Fluid[c_int].temperature)

        rho = (Fluid[c_int].pressure/(Flow_gasConstant*temperature))

        Fluid[c_bnd].rhoBoundary = rho
        Fluid[c_bnd].rhoVelocityBoundary = vs_mul(velocity, rho)
        Fluid[c_bnd].rhoEnergyBoundary = (rho*((cv*temperature)+(0.5*dot(velocity, velocity))))
      else
        var sign = BC_yPosSign
        var bnd_velocity = BC_yPosVelocity
        var bnd_temperature = BC_yPosTemperature
        var rho = 0.0
        var temp_wall = 0.0
        var temperature = 0.0
        var velocity = array(0.0, 0.0, 0.0)
        var velocity__3538 = array(0.0, 0.0, 0.0)
        velocity__3538 = vv_add(vv_mul(vs_div(Fluid[c_int].rhoVelocity, Fluid[c_int].rho), sign), bnd_velocity)
        temp_wall = Fluid[c_int].temperature
        if (bnd_temperature>0.0) then
          temp_wall = bnd_temperature
        end
        temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
        rho = (Fluid[c_int].pressure/(Flow_gasConstant*temperature))
        Fluid[c_bnd].rhoBoundary = rho
        Fluid[c_bnd].rhoVelocityBoundary = vs_mul(velocity__3538, rho)
        Fluid[c_bnd].rhoEnergyBoundary = (rho*((cv*temperature)+(0.5*dot(velocity__3538, velocity__3538))))
      end
    end
    if zNegGhost then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, 0, 1})%Fluid.bounds)
      var cv = (Flow_gasConstant/(Flow_gamma-1.0))

      if (BC_zBCLeft == SCHEMA.FlowBC_NonUniformTemperatureWall) then
        var sign = BC_zNegSign
        var bnd_velocity = BC_zNegVelocity

        var rho = 0.0
        var temperature = 0.0
        var velocity = array(0.0, 0.0, 0.0)
        var wall_temperature = 0.0

        velocity = vv_add(vv_mul(vs_div(Fluid[c_int].rhoVelocity, Fluid[c_int].rho), sign), bnd_velocity)

        var c_1 = 2.0/(Grid_xWidth*Grid_xWidth)*( (BC_zBCLeftHeat_T_right - BC_zBCLeftHeat_T_left) - 2.0*(BC_zBCLeftHeat_T_mid - BC_zBCLeftHeat_T_left))
        var c_2 = 4.0/(Grid_xWidth)*((BC_zBCLeftHeat_T_mid - BC_zBCLeftHeat_T_left) - 1.0/4.0*(BC_zBCLeftHeat_T_right - BC_zBCLeftHeat_T_left))
        var c_3 = BC_zBCLeftHeat_T_left
        wall_temperature = c_1*Fluid[c_bnd].centerCoordinates[0]*Fluid[c_bnd].centerCoordinates[0] + c_2*Fluid[c_bnd].centerCoordinates[0] + c_3
        if wall_temperature < 0.0 then --unphysical.... set wall themperature to zero
          wall_temperature = 0.0
        end
        temperature = ((2.0*wall_temperature)-Fluid[c_int].temperature)

        rho = (Fluid[c_int].pressure/(Flow_gasConstant*temperature))

        Fluid[c_bnd].rhoBoundary = rho
        Fluid[c_bnd].rhoVelocityBoundary = vs_mul(velocity, rho)
        Fluid[c_bnd].rhoEnergyBoundary = (rho*((cv*temperature)+(0.5*dot(velocity, velocity))))
      else
        var sign = BC_zNegSign
        var bnd_velocity = BC_zNegVelocity
        var bnd_temperature = BC_zNegTemperature
        var rho = 0.0
        var temp_wall = 0.0
        var temperature = 0.0
        var velocity = array(0.0, 0.0, 0.0)
        var velocity__3549 = array(0.0, 0.0, 0.0)
        velocity__3549 = vv_add(vv_mul(vs_div(Fluid[c_int].rhoVelocity, Fluid[c_int].rho), sign), bnd_velocity)
        temp_wall = Fluid[c_int].temperature
        if (bnd_temperature>0.0) then
          temp_wall = bnd_temperature
        end
        temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
        rho = (Fluid[c_int].pressure/(Flow_gasConstant*temperature))
        Fluid[c_bnd].rhoBoundary = rho
        Fluid[c_bnd].rhoVelocityBoundary = vs_mul(velocity__3549, rho)
        Fluid[c_bnd].rhoEnergyBoundary = (rho*((cv*temperature)+(0.5*dot(velocity__3549, velocity__3549))))
      end
    end
    if zPosGhost then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, 0, -1})%Fluid.bounds)
      var cv = (Flow_gasConstant/(Flow_gamma-1.0))

      if (BC_zBCRight == SCHEMA.FlowBC_NonUniformTemperatureWall)then
        var sign = BC_zPosSign
        var bnd_velocity = BC_zPosVelocity

        var rho = 0.0
        var temperature = 0.0
        var velocity = array(0.0, 0.0, 0.0)
        var wall_temperature = 0.0

        velocity = vv_add(vv_mul(vs_div(Fluid[c_int].rhoVelocity, Fluid[c_int].rho), sign), bnd_velocity)

        var c_1 = 2.0/(Grid_xWidth*Grid_xWidth)*( (BC_zBCRightHeat_T_right - BC_zBCRightHeat_T_left) - 2.0*(BC_zBCRightHeat_T_mid - BC_zBCRightHeat_T_left))
        var c_2 = 4.0/(Grid_xWidth)*((BC_zBCRightHeat_T_mid - BC_zBCRightHeat_T_left) - 1.0/4.0*(BC_zBCRightHeat_T_right - BC_zBCRightHeat_T_left))
        var c_3 = BC_zBCRightHeat_T_left
        wall_temperature = c_1*Fluid[c_bnd].centerCoordinates[0]*Fluid[c_bnd].centerCoordinates[0] + c_2*Fluid[c_bnd].centerCoordinates[0] + c_3
        if wall_temperature < 0.0 then --unphysical.... set wall themperature to zero
          wall_temperature = 0.0
        end
        temperature = ((2.0*wall_temperature)-Fluid[c_int].temperature)

        rho = (Fluid[c_int].pressure/(Flow_gasConstant*temperature))

        Fluid[c_bnd].rhoBoundary = rho
        Fluid[c_bnd].rhoVelocityBoundary = vs_mul(velocity, rho)
        Fluid[c_bnd].rhoEnergyBoundary = (rho*((cv*temperature)+(0.5*dot(velocity, velocity))))
      else
        var sign = BC_zPosSign
        var bnd_velocity = BC_zPosVelocity
        var bnd_temperature = BC_zPosTemperature
        var rho = 0.0
        var temp_wall = 0.0
        var temperature = 0.0
        var velocity = array(0.0, 0.0, 0.0)
        var velocity__3560 = array(0.0, 0.0, 0.0)
        velocity__3560 = vv_add(vv_mul(vs_div(Fluid[c_int].rhoVelocity, Fluid[c_int].rho), sign), bnd_velocity)
        temp_wall = Fluid[c_int].temperature
        if (bnd_temperature>0.0) then
          temp_wall = bnd_temperature
        end
        temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
        rho = (Fluid[c_int].pressure/(Flow_gasConstant*temperature))
        Fluid[c_bnd].rhoBoundary = rho
        Fluid[c_bnd].rhoVelocityBoundary = vs_mul(velocity__3560, rho)
        Fluid[c_bnd].rhoEnergyBoundary = (rho*((cv*temperature)+(0.5*dot(velocity__3560, velocity__3560))))
      end
    end
  end
end

-- given a valid rho*v in the interior cells + the specified BC's -> compute v in the ghost cells
__demand(__parallel, __cuda)
task Flow_UpdateGhostConservedStep2(Fluid : region(ispace(int3d), Fluid_columns),
                                    config : Config,
                                    Grid_xBnum : int32, Grid_xNum : int32,
                                    Grid_yBnum : int32, Grid_yNum : int32,
                                    Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{rhoBoundary, rhoVelocityBoundary, rhoEnergyBoundary}),
  writes(Fluid.{rho, rhoVelocity, rhoEnergy})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0)
    var xPosGhost  = (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)
    var yNegGhost = (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)
    var yPosGhost  = (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)
    var zNegGhost = (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)
    var zPosGhost  = (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)
    var ghost_cell = (xNegGhost or xPosGhost or
                      yNegGhost or yPosGhost or
                      zNegGhost or zPosGhost )
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if ghost_cell then
      Fluid[c].rho = Fluid[c].rhoBoundary
      Fluid[c].rhoVelocity = Fluid[c].rhoVelocityBoundary
      Fluid[c].rhoEnergy = Fluid[c].rhoEnergyBoundary
    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdateGhostVelocityStep1(Fluid : region(ispace(int3d), Fluid_columns),
                                   config : Config,
                                   BC_xNegVelocity : double[3], BC_xPosVelocity : double[3], BC_xNegSign : double[3], BC_xPosSign : double[3],
                                   BC_yNegVelocity : double[3], BC_yPosVelocity : double[3], BC_yNegSign : double[3], BC_yPosSign : double[3],
                                   BC_zNegVelocity : double[3], BC_zPosVelocity : double[3], BC_zNegSign : double[3], BC_zPosSign : double[3],
                                   Grid_xBnum : int32, Grid_xNum : int32,
                                   Grid_yBnum : int32, Grid_yNum : int32,
                                   Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.velocity),
  writes(Fluid.velocityBoundary)
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  __demand(__openmp)
  for c in Fluid do

    var xNegGhost = (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0)
    var xPosGhost = (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)
    var yNegGhost = (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)
    var yPosGhost = (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)
    var zNegGhost = (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)
    var zPosGhost = (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)
    var ghost_cell = (xNegGhost or xPosGhost or
                      yNegGhost or yPosGhost or
                      zNegGhost or zPosGhost )
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if (ghost_cell and not (NSCBC_inflow_cell or NSCBC_outflow_cell)) then
      if xNegGhost then
        var c_bnd = int3d(c)
        var c_int = ((c+{1, 0, 0})%Fluid.bounds)
        var sign = BC_xNegSign
        var bnd_velocity = BC_xNegVelocity
        Fluid[c_bnd].velocityBoundary = vv_add(vv_mul(Fluid[c_int].velocity, sign), bnd_velocity)
      end
      if xPosGhost then
        var c_bnd = int3d(c)
        var c_int = ((c+{-1, 0, 0})%Fluid.bounds)
        var sign = BC_xPosSign
        var bnd_velocity = BC_xPosVelocity
        Fluid[c_bnd].velocityBoundary = vv_add(vv_mul(Fluid[c_int].velocity, sign), bnd_velocity)
      end
      if (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0) then
        var c_bnd = int3d(c)
        var c_int = ((c+{0, 1, 0})%Fluid.bounds)
        var sign = BC_yNegSign
        var bnd_velocity = BC_yNegVelocity
        Fluid[c_bnd].velocityBoundary = vv_add(vv_mul(Fluid[c_int].velocity, sign), bnd_velocity)
      end
      if (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0) then
        var c_bnd = int3d(c)
        var c_int = ((c+{0, -1, 0})%Fluid.bounds)
        var sign = BC_yPosSign
        var bnd_velocity = BC_yPosVelocity
        Fluid[c_bnd].velocityBoundary = vv_add(vv_mul(Fluid[c_int].velocity, sign), bnd_velocity)
      end
      if (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0) then
        var c_bnd = int3d(c)
        var c_int = ((c+{0, 0, 1})%Fluid.bounds)
        var sign = BC_zNegSign
        var bnd_velocity = BC_zNegVelocity
        Fluid[c_bnd].velocityBoundary = vv_add(vv_mul(Fluid[c_int].velocity, sign), bnd_velocity)
      end
      if (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0) then
        var c_bnd = int3d(c)
        var c_int = ((c+{0, 0, -1})%Fluid.bounds)
        var sign = BC_zPosSign
        var bnd_velocity = BC_zPosVelocity
        Fluid[c_bnd].velocityBoundary = vv_add(vv_mul(Fluid[c_int].velocity, sign), bnd_velocity)
      end
    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdateGhostVelocityStep2(Fluid : region(ispace(int3d), Fluid_columns),
                                   config : Config,
                                   Grid_xBnum : int32, Grid_xNum : int32,
                                   Grid_yBnum : int32, Grid_yNum : int32,
                                   Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.velocityBoundary),
  writes(Fluid.velocity)
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0)
    var xPosGhost  = (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)
    var yNegGhost = (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)
    var yPosGhost  = (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)
    var zNegGhost = (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)
    var zPosGhost  = (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)
    var ghost_cell = (xNegGhost or xPosGhost or
                      yNegGhost or yPosGhost or
                      zNegGhost or zPosGhost )
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if (ghost_cell and not (NSCBC_inflow_cell or NSCBC_outflow_cell)) then
      Fluid[c].velocity = Fluid[c].velocityBoundary
    end
  end
end

__demand(__parallel, __cuda)
task Flow_ComputeVelocityGradientAll(Fluid : region(ispace(int3d), Fluid_columns),
                                     Grid_xBnum : int32, Grid_xCellWidth : double, Grid_xNum : int32,
                                     Grid_yBnum : int32, Grid_yCellWidth : double, Grid_yNum : int32,
                                     Grid_zBnum : int32, Grid_zCellWidth : double, Grid_zNum : int32)
where
  reads(Fluid.velocity),
  writes(Fluid.{velocityGradientX, velocityGradientY, velocityGradientZ})
do
  __demand(__openmp)
  for c in Fluid do
    -- if interior cell
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      Fluid[c].velocityGradientX = vs_div(vs_mul(vv_sub(Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocity, Fluid[((c+{-1, 0, 0})%Fluid.bounds)].velocity), 0.5), Grid_xCellWidth)
      Fluid[c].velocityGradientY = vs_div(vs_mul(vv_sub(Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocity, Fluid[((c+{0, -1, 0})%Fluid.bounds)].velocity), 0.5), Grid_yCellWidth)
      Fluid[c].velocityGradientZ = vs_div(vs_mul(vv_sub(Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocity, Fluid[((c+{0, 0, -1})%Fluid.bounds)].velocity), 0.5), Grid_zCellWidth)
    end
  end
end

__demand(__parallel, __cuda)
task Flow_ComputeVelocityGradientGhostNSCBC(Fluid : region(ispace(int3d), Fluid_columns),
                                            config : Config,
                                            Grid_xBnum : int32, Grid_xCellWidth : double, Grid_xNum : int32,
                                            Grid_yBnum : int32, Grid_yCellWidth : double, Grid_yNum : int32,
                                            Grid_zBnum : int32, Grid_zCellWidth : double, Grid_zNum : int32)
where
  reads(Fluid.velocity),
  writes(Fluid.{velocityGradientX, velocityGradientY, velocityGradientZ})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0)
    var xPosGhost  = (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)
    var yNegGhost = (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)
    var yPosGhost  = (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)
    var zNegGhost = (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)
    var zPosGhost  = (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)
    var ghost_cell = (xNegGhost or xPosGhost or
                      yNegGhost or yPosGhost or
                      zNegGhost or zPosGhost )
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if NSCBC_inflow_cell  then
      -- forward one sided difference
      Fluid[c].velocityGradientX = vs_div(vv_sub(Fluid[(c+{1, 0, 0})].velocity, Fluid[c].velocity), Grid_xCellWidth)

      -- centeral difference
      Fluid[c].velocityGradientY = vs_div(vs_mul(vv_sub(Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocity, Fluid[((c+{0, -1, 0})%Fluid.bounds)].velocity), 0.5), Grid_yCellWidth)
      Fluid[c].velocityGradientZ = vs_div(vs_mul(vv_sub(Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocity, Fluid[((c+{0, 0, -1})%Fluid.bounds)].velocity), 0.5), Grid_zCellWidth)
    end

    if NSCBC_outflow_cell  then
      -- backward one sided difference
      Fluid[c].velocityGradientX = vs_div(vv_sub(Fluid[c].velocity, Fluid[(c+{-1, 0, 0})].velocity), Grid_xCellWidth)

      -- centeral difference
      Fluid[c].velocityGradientY = vs_div(vs_mul(vv_sub(Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocity, Fluid[((c+{0, -1, 0})%Fluid.bounds)].velocity), 0.5), Grid_yCellWidth)
      Fluid[c].velocityGradientZ = vs_div(vs_mul(vv_sub(Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocity, Fluid[((c+{0, 0, -1})%Fluid.bounds)].velocity), 0.5), Grid_zCellWidth)

    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdateAuxiliaryThermodynamics(Fluid : region(ispace(int3d), Fluid_columns),
                                        Flow_gamma : double,
                                        Flow_gasConstant : double,
                                        Grid_xBnum : int32, Grid_xNum : int32,
                                        Grid_yBnum : int32, Grid_yNum : int32,
                                        Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{rho, velocity, rhoEnergy}),
  writes(Fluid.{pressure, temperature})
do
  __demand(__openmp)
  for c in Fluid do
    -- if interior cells
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      var kineticEnergy = ((0.5*Fluid[c].rho)*dot(Fluid[c].velocity, Fluid[c].velocity))
      var pressure = ((Flow_gamma-1.0)*(Fluid[c].rhoEnergy-kineticEnergy))
      Fluid[c].pressure = pressure
      Fluid[c].temperature = (pressure/(Flow_gasConstant*Fluid[c].rho))
    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdateAuxiliaryThermodynamicsGhostNSCBC(Fluid : region(ispace(int3d), Fluid_columns),
                                                  config : Config,
                                                  Flow_gamma : double,
                                                  Flow_gasConstant : double,
                                                  Grid_xBnum : int32, Grid_xNum : int32,
                                                  Grid_yBnum : int32, Grid_yNum : int32,
                                                  Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{rho, velocity, rhoEnergy, temperature_inc}),
  writes(Fluid.{pressure, temperature})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCLeftHeat_type = config.BC.xBCLeftHeat.type
  var BC_xBCLeftHeat_Constant_temperature = config.BC.xBCLeftHeat.u.Constant.temperature
  var BC_xBCRight = config.BC.xBCRight
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0)
    var xPosGhost  = (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)
    var yNegGhost = (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)
    var yPosGhost  = (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)
    var zNegGhost = (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)
    var zPosGhost  = (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)
    var ghost_cell = (xNegGhost or xPosGhost or
                      yNegGhost or yPosGhost or
                      zNegGhost or zPosGhost )
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if (ghost_cell) then
      if (NSCBC_inflow_cell)  then
        var kineticEnergy = (0.5*Fluid[c].rho) * dot(Fluid[c].velocity,Fluid[c].velocity)
        Fluid[c].pressure = (Flow_gamma-1.0) * (Fluid[c].rhoEnergy-kineticEnergy)
        var temperature : double
        if BC_xBCLeftHeat_type == SCHEMA.TempProfile_Constant then
          temperature = BC_xBCLeftHeat_Constant_temperature
        -- elseif BC_xBCLeftHeat_type == SCHEMA.TempProfile_Parabola then
        --   regentlib.assert(false, 'Parabola heat model not supported')
        else -- BC_xBCLeftHeat_type == SCHEMA.TempProfile_Incoming
          temperature = Fluid[c].temperature_inc
        end
        Fluid[c].temperature = temperature
      end

      if (NSCBC_outflow_cell)  then
        var kineticEnergy = ((0.5*Fluid[c].rho)*dot(Fluid[c].velocity, Fluid[c].velocity))
        var pressure = ((Flow_gamma-1.0)*(Fluid[c].rhoEnergy-kineticEnergy))
        Fluid[c].pressure = pressure
        Fluid[c].temperature = (pressure/(Flow_gasConstant*Fluid[c].rho))
      end
    end

  end
end

__demand(__parallel, __cuda)
task Flow_UpdateGhostThermodynamicsStep1(Fluid : region(ispace(int3d), Fluid_columns),
                                         config : Config,
                                         Flow_gamma : double,
                                         Flow_gasConstant : double,
                                         BC_xNegTemperature : double, BC_xPosTemperature : double,
                                         BC_yNegTemperature : double, BC_yPosTemperature : double,
                                         BC_zNegTemperature : double, BC_zPosTemperature : double,
                                         Grid_xBnum : int32, Grid_xNum : int32,
                                         Grid_yBnum : int32, Grid_yNum : int32,
                                         Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{pressure, temperature, centerCoordinates}),
  writes(Fluid.{pressureBoundary, temperatureBoundary})
do
  var BC_xBCLeft  = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  var BC_yBCLeft  = config.BC.yBCLeft
  var BC_yBCRight = config.BC.yBCRight
  var BC_zBCLeft  = config.BC.zBCLeft
  var BC_zBCRight = config.BC.zBCRight
  var Grid_xWidth = config.Grid.xWidth
  var BC_yBCLeftHeat_T_left  = config.BC.yBCLeftHeat.u.Parabola.T_left
  var BC_yBCLeftHeat_T_mid   = config.BC.yBCLeftHeat.u.Parabola.T_mid
  var BC_yBCLeftHeat_T_right = config.BC.yBCLeftHeat.u.Parabola.T_right
  var BC_yBCRightHeat_T_left  = config.BC.yBCRightHeat.u.Parabola.T_left
  var BC_yBCRightHeat_T_mid   = config.BC.yBCRightHeat.u.Parabola.T_mid
  var BC_yBCRightHeat_T_right = config.BC.yBCRightHeat.u.Parabola.T_right
  var BC_zBCLeftHeat_T_left  = config.BC.zBCLeftHeat.u.Parabola.T_left
  var BC_zBCLeftHeat_T_mid   = config.BC.zBCLeftHeat.u.Parabola.T_mid
  var BC_zBCLeftHeat_T_right = config.BC.zBCLeftHeat.u.Parabola.T_right
  var BC_zBCRightHeat_T_left  = config.BC.zBCRightHeat.u.Parabola.T_left
  var BC_zBCRightHeat_T_mid   = config.BC.zBCRightHeat.u.Parabola.T_mid
  var BC_zBCRightHeat_T_right = config.BC.zBCRightHeat.u.Parabola.T_right
  __demand(__openmp)
  for c in Fluid do

    var xNegGhost = (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0)
    var xPosGhost = (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)
    var yNegGhost = (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)
    var yPosGhost = (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)
    var zNegGhost = (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)
    var zPosGhost = (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)
    var ghost_cell = (xNegGhost or xPosGhost or
                      yNegGhost or yPosGhost or
                      zNegGhost or zPosGhost )
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if ghost_cell and not (NSCBC_inflow_cell or NSCBC_outflow_cell) then
      if xNegGhost then
        var c_bnd = int3d(c)
        var c_int = ((c+{1, 0, 0})%Fluid.bounds)
        var temperature = 0.0
        var pressure = 0.0

        pressure = Fluid[c_int].pressure

        var bnd_temperature  = BC_xNegTemperature
        var wall_temperature = Fluid[c_int].temperature
        if (bnd_temperature>0.0) then
          wall_temperature = bnd_temperature
        end
        temperature = (2.0*wall_temperature)-Fluid[c_int].temperature

        Fluid[c_bnd].pressureBoundary = pressure
        Fluid[c_bnd].temperatureBoundary = temperature

      end
      if xPosGhost then
        var c_bnd = int3d(c)
        var c_int = ((c+{-1, 0, 0})%Fluid.bounds)
        var temperature = 0.0
        var pressure = 0.0

        pressure = Fluid[c_int].pressure
        var bnd_temperature = BC_xPosTemperature
        var wall_temperature = Fluid[c_int].temperature
        if (bnd_temperature>0.0) then
          wall_temperature = bnd_temperature
        end
        temperature = ((2.0*wall_temperature)-Fluid[c_int].temperature)

        Fluid[c_bnd].pressureBoundary = pressure
        Fluid[c_bnd].temperatureBoundary = temperature
      end
      if yNegGhost then
        var c_bnd = int3d(c)
        var c_int = ((c+{0, 1, 0})%Fluid.bounds)

        if (BC_yBCLeft == SCHEMA.FlowBC_NonUniformTemperatureWall) then
          var c_1 = 2.0/(Grid_xWidth*Grid_xWidth)*( (BC_yBCLeftHeat_T_right - BC_yBCLeftHeat_T_left) - 2.0*(BC_yBCLeftHeat_T_mid - BC_yBCLeftHeat_T_left))
          var c_2 = 4.0/(Grid_xWidth)*((BC_yBCLeftHeat_T_mid - BC_yBCLeftHeat_T_left) - 1.0/4.0*(BC_yBCLeftHeat_T_right - BC_yBCLeftHeat_T_left))
          var c_3 = BC_yBCLeftHeat_T_left
          var wall_temperature = c_1*Fluid[c_bnd].centerCoordinates[0]*Fluid[c_bnd].centerCoordinates[0] + c_2*Fluid[c_bnd].centerCoordinates[0] + c_3
          if wall_temperature < 0.0 then --unphysical.... set wall themperature to zero
            wall_temperature = 0.0
          end

          var temperature = ((2.0*wall_temperature)-Fluid[c_int].temperature)
          Fluid[c_bnd].pressureBoundary = Fluid[c_int].pressure
          Fluid[c_bnd].temperatureBoundary = temperature
        else
          var bnd_temperature = BC_yNegTemperature
          var temp_wall = 0.0
          var temperature = 0.0
          temp_wall = Fluid[c_int].temperature
          if (bnd_temperature>0.0) then
            temp_wall = bnd_temperature
          end
          temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
          Fluid[c_bnd].pressureBoundary = Fluid[c_int].pressure
          Fluid[c_bnd].temperatureBoundary = temperature
        end
      end
      if yPosGhost then
        var c_bnd = int3d(c)
        var c_int = ((c+{0, -1, 0})%Fluid.bounds)
        if (BC_yBCRight == SCHEMA.FlowBC_NonUniformTemperatureWall) then
          var c_1 = 2.0/(Grid_xWidth*Grid_xWidth)*( (BC_yBCLeftHeat_T_right - BC_yBCLeftHeat_T_left) - 2.0*(BC_yBCLeftHeat_T_mid - BC_yBCLeftHeat_T_left))
          var c_2 = 4.0/(Grid_xWidth)*((BC_yBCLeftHeat_T_mid - BC_yBCLeftHeat_T_left) - 1.0/4.0*(BC_yBCLeftHeat_T_right - BC_yBCLeftHeat_T_left))
          var c_3 = BC_yBCLeftHeat_T_left
          var wall_temperature = c_1*Fluid[c_bnd].centerCoordinates[0]*Fluid[c_bnd].centerCoordinates[0] + c_2*Fluid[c_bnd].centerCoordinates[0] + c_3
          if wall_temperature < 0.0 then --unphysical.... set wall themperature to zero
            wall_temperature = 0.0
          end

          var temperature = ((2.0*wall_temperature)-Fluid[c_int].temperature)
          Fluid[c_bnd].pressureBoundary = Fluid[c_int].pressure
          Fluid[c_bnd].temperatureBoundary = temperature
        else
          var bnd_temperature = BC_yPosTemperature
          var temp_wall = 0.0
          var temperature = 0.0
          temp_wall = Fluid[c_int].temperature
          if (bnd_temperature>0.0) then
            temp_wall = bnd_temperature
          end
          temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
          Fluid[c_bnd].pressureBoundary = Fluid[c_int].pressure
          Fluid[c_bnd].temperatureBoundary = temperature
        end
      end
      if zNegGhost then
        var c_bnd = int3d(c)
        var c_int = ((c+{0, 0, 1})%Fluid.bounds)
        if (BC_zBCLeft == SCHEMA.FlowBC_NonUniformTemperatureWall) then
          var c_1 = 2.0/(Grid_xWidth*Grid_xWidth)*( (BC_yBCLeftHeat_T_right - BC_yBCLeftHeat_T_left) - 2.0*(BC_yBCLeftHeat_T_mid - BC_yBCLeftHeat_T_left))
          var c_2 = 4.0/(Grid_xWidth)*((BC_yBCLeftHeat_T_mid - BC_yBCLeftHeat_T_left) - 1.0/4.0*(BC_yBCLeftHeat_T_right - BC_yBCLeftHeat_T_left))
          var c_3 = BC_yBCLeftHeat_T_left
          var wall_temperature = c_1*Fluid[c_bnd].centerCoordinates[0]*Fluid[c_bnd].centerCoordinates[0] + c_2*Fluid[c_bnd].centerCoordinates[0] + c_3
          if wall_temperature < 0.0 then --unphysical.... set wall themperature to zero
            wall_temperature = 0.0
          end

          var temperature = ((2.0*wall_temperature)-Fluid[c_int].temperature)
          Fluid[c_bnd].pressureBoundary = Fluid[c_int].pressure
          Fluid[c_bnd].temperatureBoundary = temperature
        else
          var bnd_temperature = BC_zNegTemperature
          var temp_wall = 0.0
          var temperature = 0.0
          temp_wall = Fluid[c_int].temperature
          if (bnd_temperature>0.0) then
            temp_wall = bnd_temperature
          end
          temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
          Fluid[c_bnd].pressureBoundary = Fluid[c_int].pressure
          Fluid[c_bnd].temperatureBoundary = temperature
        end
      end
      if zPosGhost then
        var c_bnd = int3d(c)
        var c_int = ((c+{0, 0, -1})%Fluid.bounds)
        if (BC_zBCRight == SCHEMA.FlowBC_NonUniformTemperatureWall) then
          var c_1 = 2.0/(Grid_xWidth*Grid_xWidth)*( (BC_yBCLeftHeat_T_right - BC_yBCLeftHeat_T_left) - 2.0*(BC_yBCLeftHeat_T_mid - BC_yBCLeftHeat_T_left))
          var c_2 = 4.0/(Grid_xWidth)*((BC_yBCLeftHeat_T_mid - BC_yBCLeftHeat_T_left) - 1.0/4.0*(BC_yBCLeftHeat_T_right - BC_yBCLeftHeat_T_left))
          var c_3 = BC_yBCLeftHeat_T_left
          var wall_temperature = c_1*Fluid[c_bnd].centerCoordinates[0]*Fluid[c_bnd].centerCoordinates[0] + c_2*Fluid[c_bnd].centerCoordinates[0] + c_3
          if wall_temperature < 0.0 then --unphysical.... set wall themperature to zero
            wall_temperature = 0.0
          end

          var temperature = ((2.0*wall_temperature)-Fluid[c_int].temperature)
          Fluid[c_bnd].pressureBoundary = Fluid[c_int].pressure
          Fluid[c_bnd].temperatureBoundary = temperature

        else
          var bnd_temperature = BC_zPosTemperature
          var temp_wall = 0.0
          var temperature = 0.0
          temp_wall = Fluid[c_int].temperature
          if (bnd_temperature>0.0) then
            temp_wall = bnd_temperature
          end
          temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
          Fluid[c_bnd].pressureBoundary = Fluid[c_int].pressure
          Fluid[c_bnd].temperatureBoundary = temperature
        end
      end
    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdateGhostThermodynamicsStep2(Fluid : region(ispace(int3d), Fluid_columns),
                                         config : Config,
                                         Grid_xBnum : int32, Grid_xNum : int32,
                                         Grid_yBnum : int32, Grid_yNum : int32,
                                         Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{pressureBoundary, temperatureBoundary}),
  writes(Fluid.{pressure, temperature})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0)
    var xPosGhost  = (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)
    var yNegGhost = (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)
    var yPosGhost  = (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)
    var zNegGhost = (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)
    var zPosGhost  = (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)

    var ghost_cell = (xNegGhost or xPosGhost or
                      yNegGhost or yPosGhost or
                      zNegGhost or zPosGhost )

    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if ghost_cell and not (NSCBC_inflow_cell or NSCBC_outflow_cell) then
      Fluid[c].pressure = Fluid[c].pressureBoundary
      Fluid[c].temperature = Fluid[c].temperatureBoundary
    end
  end
end

__demand(__parallel, __cuda)
task Particles_CalculateNumber(Particles : region(ispace(int1d), Particles_columns))
where
  reads(Particles.__valid)
do
  var acc = int64(0)
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      acc += 1
    end
  end
  return acc
end

__demand(__parallel, __cuda)
task CalculateAveragePressure(Fluid : region(ispace(int3d), Fluid_columns),
                              Grid_cellVolume : double,
                              Grid_xBnum : int32, Grid_xNum : int32,
                              Grid_yBnum : int32, Grid_yNum : int32,
                              Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.pressure)
do
  var acc = 0.0
  __demand(__openmp)
  for c in Fluid do
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      acc += (Fluid[c].pressure*Grid_cellVolume)
    end
  end
  return acc
end

__demand(__parallel, __cuda)
task CalculateAverageTemperature(Fluid : region(ispace(int3d), Fluid_columns),
                                 Grid_cellVolume : double,
                                 Grid_xBnum : int32, Grid_xNum : int32,
                                 Grid_yBnum : int32, Grid_yNum : int32,
                                 Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.temperature)
do
  var acc = 0.0
  __demand(__openmp)
  for c in Fluid do
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      acc += (Fluid[c].temperature*Grid_cellVolume)
    end
  end
  return acc
end

__demand(__parallel, __cuda)
task CalculateAverageKineticEnergy(Fluid : region(ispace(int3d), Fluid_columns),
                                   Grid_cellVolume : double,
                                   Grid_xBnum : int32, Grid_xNum : int32,
                                   Grid_yBnum : int32, Grid_yNum : int32,
                                   Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{rho, velocity})
do
  var acc = 0.0
  __demand(__openmp)
  for c in Fluid do
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      var kineticEnergy = ((0.5*Fluid[c].rho)*dot(Fluid[c].velocity, Fluid[c].velocity))
      acc += (kineticEnergy*Grid_cellVolume)
    end
  end
  return acc
end

__demand(__parallel, __cuda)
task CalculateMinTemperature(Fluid : region(ispace(int3d), Fluid_columns),
                             Grid_xBnum : int32, Grid_xNum : int32,
                             Grid_yBnum : int32, Grid_yNum : int32,
                             Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.temperature)
do
  var acc = math.huge
  __demand(__openmp)
  for c in Fluid do
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      acc min= Fluid[c].temperature
    end
  end
  return acc
end

__demand(__parallel, __cuda)
task CalculateMaxTemperature(Fluid : region(ispace(int3d), Fluid_columns),
                             Grid_xBnum : int32, Grid_xNum : int32,
                             Grid_yBnum : int32, Grid_yNum : int32,
                             Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.temperature)
do
  var acc = -math.huge
  __demand(__openmp)
  for c in Fluid do
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      acc max= Fluid[c].temperature
    end
  end
  return acc
end

__demand(__parallel, __cuda)
task Particles_IntegrateQuantities(Particles : region(ispace(int1d), Particles_columns))
where
  reads(Particles.{temperature, __valid})
do
  var acc = 0.0
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      acc += Particles[p].temperature
    end
  end
  return acc
end

__demand(__inline)
task GetSoundSpeed(temperature : double, Flow_gamma : double, Flow_gasConstant : double)
  return sqrt(((Flow_gamma*Flow_gasConstant)*temperature))
end

__demand(__parallel, __cuda)
task CalculateMaxMachNumber(Fluid : region(ispace(int3d), Fluid_columns),
                            config : Config,
                            Flow_gamma : double,
                            Flow_gasConstant : double,
                            Grid_xBnum : int32, Grid_xNum : int32,
                            Grid_yBnum : int32, Grid_yNum : int32,
                            Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{velocity, temperature})
do
  var acc = -math.huge
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  __demand(__openmp)
  for c in Fluid do

    var xNegGhost = (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0)
    var xPosGhost = (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)
    var yNegGhost = (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)
    var yPosGhost = (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)
    var zNegGhost = (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)
    var zPosGhost = (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)
    var ghost_cell = (xNegGhost or xPosGhost or
                      yNegGhost or yPosGhost or
                      zNegGhost or zPosGhost )
    var interior_cell = not (ghost_cell)
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    if interior_cell or NSCBC_inflow_cell or NSCBC_outflow_cell then
      var c_sound = GetSoundSpeed(Fluid[c].temperature, Flow_gamma, Flow_gasConstant)
      acc max= (Fluid[c].velocity[0]*Fluid[c].velocity[0] + Fluid[c].velocity[1]*Fluid[c].velocity[1] + Fluid[c].velocity[2]*Fluid[c].velocity[2]) / c_sound
    end
  end
  return acc
end

__demand(__parallel, __cuda)
task CalculateConvectiveSpectralRadius(Fluid : region(ispace(int3d), Fluid_columns),
                                       Flow_gamma : double,
                                       Flow_gasConstant : double,
                                       Grid_dXYZInverseSquare : double,
                                       Grid_xCellWidth : double, Grid_yCellWidth : double, Grid_zCellWidth : double)
where
  reads(Fluid.{velocity, temperature}),
  writes(Fluid.convectiveSpectralRadius)
do
  var acc = -math.huge
  __demand(__openmp)
  for c in Fluid do
    var tmp =  ((((fabs(Fluid[c].velocity[0])/Grid_xCellWidth)+(fabs(Fluid[c].velocity[1])/Grid_yCellWidth))+(fabs(Fluid[c].velocity[2])/Grid_zCellWidth))+(GetSoundSpeed(Fluid[c].temperature, Flow_gamma, Flow_gasConstant)*sqrt(Grid_dXYZInverseSquare)))
    Fluid[c].convectiveSpectralRadius = tmp
    acc max= tmp
  end
  return acc
end


__demand(__parallel, __cuda)
task CalculateViscousSpectralRadius(Fluid : region(ispace(int3d), Fluid_columns),
                                    Flow_constantVisc : double,
                                    Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                                    Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                                    Flow_viscosityModel : SCHEMA.ViscosityModel,
                                    Grid_dXYZInverseSquare : double)
where
  reads(Fluid.{rho, temperature}),
  writes(Fluid.viscousSpectralRadius)
do
  var acc = -math.huge
  __demand(__openmp)
  for c in Fluid do
    var dynamicViscosity = GetDynamicViscosity(Fluid[c].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)
    var tmp = ((((2.0*dynamicViscosity)/Fluid[c].rho)*Grid_dXYZInverseSquare)*4.0)
    Fluid[c].viscousSpectralRadius = tmp
    acc max= tmp
  end
  return acc
end

__demand(__parallel, __cuda)
task CalculateHeatConductionSpectralRadius(Fluid : region(ispace(int3d), Fluid_columns),
                                           Flow_constantVisc : double,
                                           Flow_gamma : double,
                                           Flow_gasConstant : double,
                                           Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                                           Flow_prandtl : double,
                                           Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                                           Flow_viscosityModel : SCHEMA.ViscosityModel,
                                           Grid_dXYZInverseSquare : double)
where
  reads(Fluid.{rho, temperature}),
  writes(Fluid.heatConductionSpectralRadius)
do
  var acc = -math.huge
  __demand(__openmp)
  for c in Fluid do
    var dynamicViscosity = GetDynamicViscosity(Fluid[c].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)
    var cv = (Flow_gasConstant/(Flow_gamma-1.0))
    var cp = (Flow_gamma*cv)
    var kappa = ((cp/Flow_prandtl)*dynamicViscosity)
    var tmp = (((kappa/(cv*Fluid[c].rho))*Grid_dXYZInverseSquare)*4.0)
    Fluid[c].heatConductionSpectralRadius = tmp
    acc max= tmp
  end
  return acc
end

__demand(__parallel, __cuda)
task Flow_InitializeTemporaries(Fluid : region(ispace(int3d), Fluid_columns))
where
  reads(Fluid.{rho, rhoEnergy, rhoVelocity}),
  writes(Fluid.{rhoEnergy_new, rhoEnergy_old, rhoVelocity_new, rhoVelocity_old, rho_new, rho_old})
do
  __demand(__openmp)
  for c in Fluid do
    Fluid[c].rho_old = Fluid[c].rho
    Fluid[c].rhoVelocity_old = Fluid[c].rhoVelocity
    Fluid[c].rhoEnergy_old = Fluid[c].rhoEnergy
    Fluid[c].rho_new = Fluid[c].rho
    Fluid[c].rhoVelocity_new = Fluid[c].rhoVelocity
    Fluid[c].rhoEnergy_new = Fluid[c].rhoEnergy
  end
end

__demand(__parallel, __cuda)
task Particles_InitializeTemporaries(Particles : region(ispace(int1d), Particles_columns))
where
  reads(Particles.{position, velocity, temperature, __valid}),
  writes(Particles.{position_new, position_old, temperature_new, temperature_old, velocity_new, velocity_old})
do
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      Particles[p].position_old = Particles[p].position
      Particles[p].velocity_old = Particles[p].velocity
      Particles[p].temperature_old = Particles[p].temperature
      Particles[p].position_new = Particles[p].position
      Particles[p].velocity_new = Particles[p].velocity
      Particles[p].temperature_new = Particles[p].temperature
    end
  end
end

__demand(__parallel, __cuda)
task Flow_InitializeTimeDerivatives(Fluid : region(ispace(int3d), Fluid_columns))
where
  reads(Fluid.{pressure, rhoEnergy}),
  writes(Fluid.{rho_t, rhoVelocity_t, rhoEnergy_t, rhoEnthalpy})
do
  __demand(__openmp)
  for c in Fluid do
    Fluid[c].rho_t = 0.0
    Fluid[c].rhoVelocity_t = array(0.0, 0.0, 0.0)
    Fluid[c].rhoEnergy_t = 0.0
    Fluid[c].rhoEnthalpy = (Fluid[c].rhoEnergy+Fluid[c].pressure)
  end
end

__demand(__parallel, __cuda)
task Particles_InitializeTimeDerivatives(Particles : region(ispace(int1d), Particles_columns))
where
  reads(Particles.__valid),
  writes(Particles.{position_t, velocity_t, temperature_t})
do
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      Particles[p].position_t = array(0.0, 0.0, 0.0)
      Particles[p].velocity_t = array(0.0, 0.0, 0.0)
      Particles[p].temperature_t = 0.0
    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdateGhostVelocityGradientStep1(Fluid : region(ispace(int3d), Fluid_columns),
                                           config : Config,
                                           BC_xNegSign : double[3], BC_yNegSign : double[3], BC_zNegSign : double[3],
                                           BC_xPosSign : double[3], BC_yPosSign : double[3], BC_zPosSign : double[3],
                                           Grid_xBnum : int32, Grid_xCellWidth : double, Grid_xNum : int32,
                                           Grid_yBnum : int32, Grid_yCellWidth : double, Grid_yNum : int32,
                                           Grid_zBnum : int32, Grid_zCellWidth : double, Grid_zNum : int32)
where
  reads(Fluid.{velocityGradientX, velocityGradientY, velocityGradientZ}),
  writes(Fluid.{velocityGradientXBoundary, velocityGradientYBoundary, velocityGradientZBoundary})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0)
    var xPosGhost  = (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)
    var yNegGhost = (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)
    var yPosGhost  = (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)
    var zNegGhost = (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)
    var zPosGhost  = (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)
    var ghost_cell = (xNegGhost or xPosGhost or
                      yNegGhost or yPosGhost or
                      zNegGhost or zPosGhost )
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if ghost_cell and not (NSCBC_inflow_cell or NSCBC_outflow_cell) then
      if (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) then -- x- boundary
        var c_bnd = int3d(c)
        var c_int = ((c+{1, 0, 0})%Fluid.bounds)
        var sign = BC_xNegSign
        Fluid[c_bnd].velocityGradientXBoundary = vv_mul(sign, Fluid[c_int].velocityGradientX)
        Fluid[c_bnd].velocityGradientYBoundary = vv_mul(sign, Fluid[c_int].velocityGradientY)
        Fluid[c_bnd].velocityGradientZBoundary = vv_mul(sign, Fluid[c_int].velocityGradientZ)
      end
      if (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0) then -- x+ boundary
        var c_bnd = int3d(c)
        var c_int = ((c+{-1, 0, 0})%Fluid.bounds)
        var sign = BC_xPosSign
        Fluid[c_bnd].velocityGradientXBoundary = vv_mul(sign, Fluid[c_int].velocityGradientX)
        Fluid[c_bnd].velocityGradientYBoundary = vv_mul(sign, Fluid[c_int].velocityGradientY)
        Fluid[c_bnd].velocityGradientZBoundary = vv_mul(sign, Fluid[c_int].velocityGradientZ)
      end
      if (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0) then
        var c_bnd = int3d(c)
        var c_int = ((c+{0, 1, 0})%Fluid.bounds)
        var sign = BC_yNegSign
        Fluid[c_bnd].velocityGradientXBoundary = vv_mul(sign, Fluid[c_int].velocityGradientX)
        Fluid[c_bnd].velocityGradientYBoundary = vv_mul(sign, Fluid[c_int].velocityGradientY)
        Fluid[c_bnd].velocityGradientZBoundary = vv_mul(sign, Fluid[c_int].velocityGradientZ)
      end
      if (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0) then
        var c_bnd = int3d(c)
        var c_int = ((c+{0, -1, 0})%Fluid.bounds)
        var sign = BC_yPosSign
        Fluid[c_bnd].velocityGradientXBoundary = vv_mul(sign, Fluid[c_int].velocityGradientX)
        Fluid[c_bnd].velocityGradientYBoundary = vv_mul(sign, Fluid[c_int].velocityGradientY)
        Fluid[c_bnd].velocityGradientZBoundary = vv_mul(sign, Fluid[c_int].velocityGradientZ)
      end
      if (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0) then
        var c_bnd = int3d(c)
        var c_int = ((c+{0, 0, 1})%Fluid.bounds)
        var sign = BC_zNegSign
        Fluid[c_bnd].velocityGradientXBoundary = vv_mul(sign, Fluid[c_int].velocityGradientX)
        Fluid[c_bnd].velocityGradientYBoundary = vv_mul(sign, Fluid[c_int].velocityGradientY)
        Fluid[c_bnd].velocityGradientZBoundary = vv_mul(sign, Fluid[c_int].velocityGradientZ)
      end
      if (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0) then
        var c_bnd = int3d(c)
        var c_int = ((c+{0, 0, -1})%Fluid.bounds)
        var sign = BC_zPosSign
        Fluid[c_bnd].velocityGradientXBoundary = vv_mul(sign, Fluid[c_int].velocityGradientX)
        Fluid[c_bnd].velocityGradientYBoundary = vv_mul(sign, Fluid[c_int].velocityGradientY)
        Fluid[c_bnd].velocityGradientZBoundary = vv_mul(sign, Fluid[c_int].velocityGradientZ)
      end
    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdateGhostVelocityGradientStep2(Fluid : region(ispace(int3d), Fluid_columns),
                                           config : Config,
                                           Grid_xBnum : int32, Grid_xNum : int32,
                                           Grid_yBnum : int32, Grid_yNum : int32,
                                           Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{velocityGradientXBoundary, velocityGradientYBoundary, velocityGradientZBoundary}),
  writes(Fluid.{velocityGradientX, velocityGradientY, velocityGradientZ})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0)
    var xPosGhost  = (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)
    var yNegGhost = (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)
    var yPosGhost  = (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)
    var zNegGhost = (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)
    var zPosGhost  = (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)
    var ghost_cell = (xNegGhost or xPosGhost or
                      yNegGhost or yPosGhost or
                      zNegGhost or zPosGhost )
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if ghost_cell and not (NSCBC_inflow_cell or NSCBC_outflow_cell) then
        Fluid[c].velocityGradientX = Fluid[c].velocityGradientXBoundary
        Fluid[c].velocityGradientY = Fluid[c].velocityGradientYBoundary
        Fluid[c].velocityGradientZ = Fluid[c].velocityGradientZBoundary
    end
  end
end

__demand(__inline)
task CenteredInviscidFlux(c_l : int3d,
                          c_r : int3d,
                          Fluid : region(ispace(int3d), Fluid_columns))
where
  reads(Fluid.{rho, velocity, pressure, rhoVelocity, rhoEnthalpy})
do
  var rhoFactorDiagonal = 0.0
  var rhoVelocityFactorDiagonal = array(0.0, 0.0, 0.0)
  var rhoEnergyFactorDiagonal = 0.0
  var fpdiag = 0.0
  rhoFactorDiagonal = (0.5*((Fluid[c_l].rho*Fluid[c_l].velocity[0])+(Fluid[c_r].rho*Fluid[c_r].velocity[0])))
  rhoVelocityFactorDiagonal = vs_mul(vv_add(vs_mul(Fluid[c_l].rhoVelocity, Fluid[c_l].velocity[0]), vs_mul(Fluid[c_r].rhoVelocity, Fluid[c_r].velocity[0])), 0.5)
  rhoEnergyFactorDiagonal = (0.5*((Fluid[c_l].rhoEnthalpy*Fluid[c_l].velocity[0])+(Fluid[c_r].rhoEnthalpy*Fluid[c_r].velocity[0])))
  fpdiag += (0.5*(Fluid[c_l].pressure+Fluid[c_r].pressure))
  var rhoFactorSkew = 0.0
  var rhoVelocityFactorSkew = array(0.0, 0.0, 0.0)
  var rhoEnergyFactorSkew = 0.0
  var tmp = 0.0
  tmp = (0.5*Fluid[c_r].velocity[0])
  rhoFactorSkew += (Fluid[c_l].rho*tmp)
  var tmp__6137 = vs_mul(Fluid[c_l].rhoVelocity, tmp)
  var v = rhoVelocityFactorSkew
  v[0] += tmp__6137[0]
  v[1] += tmp__6137[1]
  v[2] += tmp__6137[2]
  rhoVelocityFactorSkew = v
  rhoEnergyFactorSkew += (Fluid[c_l].rhoEnthalpy*tmp)
  tmp = (0.5*Fluid[c_l].velocity[0])
  rhoFactorSkew += (Fluid[c_r].rho*tmp)
  var tmp__6139 = vs_mul(Fluid[c_r].rhoVelocity, tmp)
  var v__6140 = rhoVelocityFactorSkew
  v__6140[0] += tmp__6139[0]
  v__6140[1] += tmp__6139[1]
  v__6140[2] += tmp__6139[2]
  rhoVelocityFactorSkew = v__6140
  rhoEnergyFactorSkew += (Fluid[c_r].rhoEnthalpy*tmp)
  var s = 0.5
  var rhoFlux_temp = ((s*rhoFactorDiagonal)+((1.0-s)*rhoFactorSkew))
  var rhoVelocityFlux_temp = vv_add(vs_mul(rhoVelocityFactorDiagonal, s), vs_mul(rhoVelocityFactorSkew, (1.0-s)))
  var rhoEnergyFlux_temp = ((s*rhoEnergyFactorDiagonal)+((1.0-s)*rhoEnergyFactorSkew))
  rhoVelocityFlux_temp[0] += fpdiag
  return array(rhoFlux_temp, rhoVelocityFlux_temp[0], rhoVelocityFlux_temp[1], rhoVelocityFlux_temp[2], rhoEnergyFlux_temp)
end

__demand(__inline)
task CenteredInviscidFlux_(c_l : int3d,
                           c_r : int3d,
                           Fluid : region(ispace(int3d), Fluid_columns))
where
  reads(Fluid.{rho, pressure, velocity, rhoVelocity, rhoEnthalpy})
do
  var rhoFactorDiagonal = 0.0
  var rhoVelocityFactorDiagonal = array(0.0, 0.0, 0.0)
  var rhoEnergyFactorDiagonal = 0.0
  var fpdiag = 0.0
  rhoFactorDiagonal = (0.5*((Fluid[c_l].rho*Fluid[c_l].velocity[1])+(Fluid[c_r].rho*Fluid[c_r].velocity[1])))
  rhoVelocityFactorDiagonal = vs_mul(vv_add(vs_mul(Fluid[c_l].rhoVelocity, Fluid[c_l].velocity[1]), vs_mul(Fluid[c_r].rhoVelocity, Fluid[c_r].velocity[1])), 0.5)
  rhoEnergyFactorDiagonal = (0.5*((Fluid[c_l].rhoEnthalpy*Fluid[c_l].velocity[1])+(Fluid[c_r].rhoEnthalpy*Fluid[c_r].velocity[1])))
  fpdiag += (0.5*(Fluid[c_l].pressure+Fluid[c_r].pressure))
  var rhoFactorSkew = 0.0
  var rhoVelocityFactorSkew = array(0.0, 0.0, 0.0)
  var rhoEnergyFactorSkew = 0.0
  var tmp = 0.0
  tmp = (0.5*Fluid[c_r].velocity[1])
  rhoFactorSkew += (Fluid[c_l].rho*tmp)
  var tmp__6344 = vs_mul(Fluid[c_l].rhoVelocity, tmp)
  var v = rhoVelocityFactorSkew
  v[0] += tmp__6344[0]
  v[1] += tmp__6344[1]
  v[2] += tmp__6344[2]
  rhoVelocityFactorSkew = v
  rhoEnergyFactorSkew += (Fluid[c_l].rhoEnthalpy*tmp)
  tmp = (0.5*Fluid[c_l].velocity[1])
  rhoFactorSkew += (Fluid[c_r].rho*tmp)
  var tmp__6346 = vs_mul(Fluid[c_r].rhoVelocity, tmp)
  var v__6347 = rhoVelocityFactorSkew
  v__6347[0] += tmp__6346[0]
  v__6347[1] += tmp__6346[1]
  v__6347[2] += tmp__6346[2]
  rhoVelocityFactorSkew = v__6347
  rhoEnergyFactorSkew += (Fluid[c_r].rhoEnthalpy*tmp)
  var s = 0.5
  var rhoFlux_temp = ((s*rhoFactorDiagonal)+((1.0-s)*rhoFactorSkew))
  var rhoVelocityFlux_temp = vv_add(vs_mul(rhoVelocityFactorDiagonal, s), vs_mul(rhoVelocityFactorSkew, (1.0-s)))
  var rhoEnergyFlux_temp = ((s*rhoEnergyFactorDiagonal)+((1.0-s)*rhoEnergyFactorSkew))
  rhoVelocityFlux_temp[1] += fpdiag
  return array(rhoFlux_temp, rhoVelocityFlux_temp[0], rhoVelocityFlux_temp[1], rhoVelocityFlux_temp[2], rhoEnergyFlux_temp)
end

__demand(__inline)
task CenteredInviscidFlux__(c_l : int3d,
                            c_r : int3d,
                            Fluid : region(ispace(int3d), Fluid_columns))
where
  reads(Fluid.{rho, pressure, velocity, rhoVelocity, rhoEnthalpy})
do
  var rhoFactorDiagonal = 0.0
  var rhoVelocityFactorDiagonal = array(0.0, 0.0, 0.0)
  var rhoEnergyFactorDiagonal = 0.0
  var fpdiag = 0.0
  rhoFactorDiagonal = (0.5*((Fluid[c_l].rho*Fluid[c_l].velocity[2])+(Fluid[c_r].rho*Fluid[c_r].velocity[2])))
  rhoVelocityFactorDiagonal = vs_mul(vv_add(vs_mul(Fluid[c_l].rhoVelocity, Fluid[c_l].velocity[2]), vs_mul(Fluid[c_r].rhoVelocity, Fluid[c_r].velocity[2])), 0.5)
  rhoEnergyFactorDiagonal = (0.5*((Fluid[c_l].rhoEnthalpy*Fluid[c_l].velocity[2])+(Fluid[c_r].rhoEnthalpy*Fluid[c_r].velocity[2])))
  fpdiag += (0.5*(Fluid[c_l].pressure+Fluid[c_r].pressure))
  var rhoFactorSkew = 0.0
  var rhoVelocityFactorSkew = array(0.0, 0.0, 0.0)
  var rhoEnergyFactorSkew = 0.0
  var tmp = 0.0
  tmp = (0.5*Fluid[c_r].velocity[2])
  rhoFactorSkew += (Fluid[c_l].rho*tmp)
  var tmp__6551 = vs_mul(Fluid[c_l].rhoVelocity, tmp)
  var v = rhoVelocityFactorSkew
  v[0] += tmp__6551[0]
  v[1] += tmp__6551[1]
  v[2] += tmp__6551[2]
  rhoVelocityFactorSkew = v
  rhoEnergyFactorSkew += (Fluid[c_l].rhoEnthalpy*tmp)
  tmp = (0.5*Fluid[c_l].velocity[2])
  rhoFactorSkew += (Fluid[c_r].rho*tmp)
  var tmp__6553 = vs_mul(Fluid[c_r].rhoVelocity, tmp)
  var v__6554 = rhoVelocityFactorSkew
  v__6554[0] += tmp__6553[0]
  v__6554[1] += tmp__6553[1]
  v__6554[2] += tmp__6553[2]
  rhoVelocityFactorSkew = v__6554
  rhoEnergyFactorSkew += (Fluid[c_r].rhoEnthalpy*tmp)
  var s = 0.5
  var rhoFlux_temp = ((s*rhoFactorDiagonal)+((1.0-s)*rhoFactorSkew))
  var rhoVelocityFlux_temp = vv_add(vs_mul(rhoVelocityFactorDiagonal, s), vs_mul(rhoVelocityFactorSkew, (1.0-s)))
  var rhoEnergyFlux_temp = ((s*rhoEnergyFactorDiagonal)+((1.0-s)*rhoEnergyFactorSkew))
  rhoVelocityFlux_temp[2] += fpdiag
  return array(rhoFlux_temp, rhoVelocityFlux_temp[0], rhoVelocityFlux_temp[1], rhoVelocityFlux_temp[2], rhoEnergyFlux_temp)
end

__demand(__parallel, __cuda)
task Flow_AddGetFlux(Fluid : region(ispace(int3d), Fluid_columns),
                     config : Config,
                     Flow_constantVisc : double,
                     Flow_gamma : double,
                     Flow_gasConstant : double,
                     Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                     Flow_prandtl : double,
                     Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                     Flow_viscosityModel : SCHEMA.ViscosityModel,
                     Grid_xBnum : int32, Grid_xCellWidth : double, Grid_xNum : int32,
                     Grid_yBnum : int32, Grid_yCellWidth : double, Grid_yNum : int32,
                     Grid_zBnum : int32, Grid_zCellWidth : double, Grid_zNum : int32)
where
  reads(Fluid.{rho, pressure, velocity, rhoVelocity, rhoEnthalpy, temperature}),
  reads(Fluid.{velocityGradientX, velocityGradientY, velocityGradientZ}),
  reads writes(Fluid.{rhoEnergyFluxX, rhoEnergyFluxY, rhoEnergyFluxZ}),
  reads writes(Fluid.{rhoFluxX, rhoFluxY, rhoFluxZ}),
  reads writes(Fluid.{rhoVelocityFluxX, rhoVelocityFluxY, rhoVelocityFluxZ})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  var recip_xCellWidth = 1 / (Grid_xCellWidth * 0.5)
  var recip_yCellWidth = 1 / (Grid_yCellWidth * 0.5)
  var recip_zCellWidth = 1 / (Grid_zCellWidth * 0.5)

  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0)
    var xPosGhost  = (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)
    var yNegGhost = (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)
    var yPosGhost  = (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)
    var zNegGhost = (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)
    var zPosGhost  = (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)
    var ghost_cell = (xNegGhost or xPosGhost or
                      yNegGhost or yPosGhost or
                      zNegGhost or zPosGhost )
    var interior_cell = not (ghost_cell)
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if interior_cell or xNegGhost  then
      var stencil = (c + {1, 0, 0}) % Fluid.bounds
      var flux = CenteredInviscidFlux(int3d(c), (c + {1, 0, 0}) % Fluid.bounds, Fluid)
      Fluid[c].rhoFluxX = flux[0]
      Fluid[c].rhoVelocityFluxX = array(flux[1], flux[2], flux[3])
      Fluid[c].rhoEnergyFluxX = flux[4]

      var temperature = Fluid[c].temperature
      var temperature_stencil = Fluid[stencil].temperature
      var muFace = 0.5 * GetDynamicViscosity(temperature,
                                             Flow_constantVisc,
                                             Flow_powerlawTempRef,
                                             Flow_powerlawViscRef,
                                             Flow_sutherlandSRef,
                                             Flow_sutherlandTempRef,
                                             Flow_sutherlandViscRef,
                                             Flow_viscosityModel)
                 + 0.5 * GetDynamicViscosity(temperature_stencil,
                                             Flow_constantVisc,
                                             Flow_powerlawTempRef,
                                             Flow_powerlawViscRef,
                                             Flow_sutherlandSRef,
                                             Flow_sutherlandTempRef,
                                             Flow_sutherlandViscRef,
                                             Flow_viscosityModel)

      var velocity_stencil = Fluid[stencil].velocity
      var velocity = Fluid[c].velocity

      var velocityFace = vs_mul(vv_add(velocity, velocity_stencil), 0.5)
      var velocityX_YFace = 0.5 * (Fluid[c].velocityGradientY[0] + Fluid[stencil].velocityGradientY[0])
      var velocityX_ZFace = 0.5 * (Fluid[c].velocityGradientZ[0] + Fluid[stencil].velocityGradientZ[0])
      var velocityY_YFace = 0.5 * (Fluid[c].velocityGradientY[1] + Fluid[stencil].velocityGradientY[1])
      var velocityZ_ZFace = 0.5 * (Fluid[c].velocityGradientZ[2] + Fluid[stencil].velocityGradientZ[2])

      var velocityX_XFace   = (0.5 * (velocity_stencil[0] - velocity[0])) * recip_xCellWidth
      var velocityY_XFace   = (0.5 * (velocity_stencil[1] - velocity[1])) * recip_xCellWidth
      var velocityZ_XFace   = (0.5 * (velocity_stencil[2] - velocity[2])) * recip_xCellWidth
      var temperature_XFace = (0.5 * (temperature_stencil - temperature)) * recip_xCellWidth

      var sigmaXX = ((muFace*(((4.0*velocityX_XFace)-(2.0*velocityY_YFace))-(2.0*velocityZ_ZFace)))/3.0)
      var sigmaYX = (muFace*(velocityY_XFace+velocityX_YFace))
      var sigmaZX = (muFace*(velocityZ_XFace+velocityX_ZFace))
      var usigma = (((velocityFace[0]*sigmaXX)+(velocityFace[1]*sigmaYX))+(velocityFace[2]*sigmaZX))
      var cp = ((Flow_gamma*Flow_gasConstant)/(Flow_gamma-1.0))
      var heatFlux = ((-((cp*muFace)/Flow_prandtl))*temperature_XFace)
      Fluid[c].rhoVelocityFluxX[0] += (-sigmaXX)
      Fluid[c].rhoVelocityFluxX[1] += (-sigmaYX)
      Fluid[c].rhoVelocityFluxX[2] += (-sigmaZX)
      Fluid[c].rhoEnergyFluxX += (-(usigma-heatFlux))
    end
    if interior_cell or yNegGhost  then
      var stencil = (c + {0, 1, 0}) % Fluid.bounds
      var flux = CenteredInviscidFlux_(int3d(c), (c + {0, 1, 0}) % Fluid.bounds, Fluid)
      Fluid[c].rhoFluxY = flux[0]
      Fluid[c].rhoVelocityFluxY = array(flux[1], flux[2], flux[3])
      Fluid[c].rhoEnergyFluxY = flux[4]

      var temperature = Fluid[c].temperature
      var temperature_stencil = Fluid[stencil].temperature
      var muFace = 0.5 * GetDynamicViscosity(temperature,
                                             Flow_constantVisc,
                                             Flow_powerlawTempRef,
                                             Flow_powerlawViscRef,
                                             Flow_sutherlandSRef,
                                             Flow_sutherlandTempRef,
                                             Flow_sutherlandViscRef,
                                             Flow_viscosityModel)
                 + 0.5 * GetDynamicViscosity(temperature_stencil,
                                             Flow_constantVisc,
                                             Flow_powerlawTempRef,
                                             Flow_powerlawViscRef,
                                             Flow_sutherlandSRef,
                                             Flow_sutherlandTempRef,
                                             Flow_sutherlandViscRef,
                                             Flow_viscosityModel)

      var velocity_stencil = Fluid[stencil].velocity
      var velocity = Fluid[c].velocity

      var velocityFace = vs_mul(vv_add(velocity, velocity_stencil), 0.5)
      var velocityY_XFace = 0.5 * (Fluid[c].velocityGradientX[1] + Fluid[stencil].velocityGradientX[1])
      var velocityY_ZFace = 0.5 * (Fluid[c].velocityGradientZ[1] + Fluid[stencil].velocityGradientZ[1])
      var velocityX_XFace = 0.5 * (Fluid[c].velocityGradientX[0] + Fluid[stencil].velocityGradientX[0])
      var velocityZ_ZFace = 0.5 * (Fluid[c].velocityGradientZ[2] + Fluid[stencil].velocityGradientZ[2])

      var velocityX_YFace   = (0.5 * (velocity_stencil[0] - velocity[0])) * recip_yCellWidth
      var velocityY_YFace   = (0.5 * (velocity_stencil[1] - velocity[1])) * recip_yCellWidth
      var velocityZ_YFace   = (0.5 * (velocity_stencil[2] - velocity[2])) * recip_yCellWidth
      var temperature_YFace = (0.5 * (temperature_stencil - temperature)) * recip_yCellWidth

      var sigmaXY = (muFace*(velocityX_YFace+velocityY_XFace))
      var sigmaYY = ((muFace*(((4.0*velocityY_YFace)-(2.0*velocityX_XFace))-(2.0*velocityZ_ZFace)))/3.0)
      var sigmaZY = (muFace*(velocityZ_YFace+velocityY_ZFace))
      var usigma = (((velocityFace[0]*sigmaXY)+(velocityFace[1]*sigmaYY))+(velocityFace[2]*sigmaZY))
      var cp = ((Flow_gamma*Flow_gasConstant)/(Flow_gamma-1.0))
      var heatFlux = ((-((cp*muFace)/Flow_prandtl))*temperature_YFace)
      Fluid[c].rhoVelocityFluxY[0] += (-sigmaXY)
      Fluid[c].rhoVelocityFluxY[1] += (-sigmaYY)
      Fluid[c].rhoVelocityFluxY[2] += (-sigmaZY)
      Fluid[c].rhoEnergyFluxY += (-(usigma-heatFlux))
    end
    if interior_cell or zNegGhost then
      var stencil = (c + {0, 0, 1}) % Fluid.bounds
      var flux = CenteredInviscidFlux__(int3d(c), (c + {0, 0, 1}) % Fluid.bounds, Fluid)
      Fluid[c].rhoFluxZ = flux[0]
      Fluid[c].rhoVelocityFluxZ = array(flux[1], flux[2], flux[3])
      Fluid[c].rhoEnergyFluxZ = flux[4]

      var temperature = Fluid[c].temperature
      var temperature_stencil = Fluid[stencil].temperature
      var muFace = 0.5 * GetDynamicViscosity(temperature,
                                             Flow_constantVisc,
                                             Flow_powerlawTempRef,
                                             Flow_powerlawViscRef,
                                             Flow_sutherlandSRef,
                                             Flow_sutherlandTempRef,
                                             Flow_sutherlandViscRef,
                                             Flow_viscosityModel)
                 + 0.5 * GetDynamicViscosity(temperature_stencil,
                                             Flow_constantVisc,
                                             Flow_powerlawTempRef,
                                             Flow_powerlawViscRef,
                                             Flow_sutherlandSRef,
                                             Flow_sutherlandTempRef,
                                             Flow_sutherlandViscRef,
                                             Flow_viscosityModel)

      var velocity_stencil = Fluid[stencil].velocity
      var velocity = Fluid[c].velocity

      var velocityFace = vs_mul(vv_add(Fluid[c].velocity, Fluid[stencil].velocity), 0.5)
      var velocityZ_XFace = 0.5 * (Fluid[c].velocityGradientX[2] + Fluid[stencil].velocityGradientX[2])
      var velocityZ_YFace = 0.5 * (Fluid[c].velocityGradientY[2] + Fluid[stencil].velocityGradientY[2])
      var velocityX_XFace = 0.5 * (Fluid[c].velocityGradientX[0] + Fluid[stencil].velocityGradientX[0])
      var velocityY_YFace = 0.5 * (Fluid[c].velocityGradientY[1] + Fluid[stencil].velocityGradientY[1])

      var velocityX_ZFace   = (0.5 * (velocity_stencil[0] - velocity[0])) * recip_zCellWidth
      var velocityY_ZFace   = (0.5 * (velocity_stencil[1] - velocity[1])) * recip_zCellWidth
      var velocityZ_ZFace   = (0.5 * (velocity_stencil[2] - velocity[2])) * recip_zCellWidth
      var temperature_ZFace = (0.5 * (temperature_stencil - temperature)) * recip_zCellWidth

      var sigmaXZ = (muFace*(velocityX_ZFace+velocityZ_XFace))
      var sigmaYZ = (muFace*(velocityY_ZFace+velocityZ_YFace))
      var sigmaZZ = ((muFace*(((4.0*velocityZ_ZFace)-(2.0*velocityX_XFace))-(2.0*velocityY_YFace)))/3.0)
      var usigma = (((velocityFace[0]*sigmaXZ)+(velocityFace[1]*sigmaYZ))+(velocityFace[2]*sigmaZZ))
      var cp = ((Flow_gamma*Flow_gasConstant)/(Flow_gamma-1.0))
      var heatFlux = ((-((cp*muFace)/Flow_prandtl))*temperature_ZFace)
      Fluid[c].rhoVelocityFluxZ[0] += (-sigmaXZ)
      Fluid[c].rhoVelocityFluxZ[1] += (-sigmaYZ)
      Fluid[c].rhoVelocityFluxZ[2] += (-sigmaZZ)
      Fluid[c].rhoEnergyFluxZ += (-(usigma-heatFlux))
    end
  end
end

__demand(__parallel, __cuda)
task Flow_AddUpdateUsingFlux(Fluid : region(ispace(int3d), Fluid_columns),
                             Grid_xBnum : int32, Grid_xCellWidth : double, Grid_xNum : int32,
                             Grid_yBnum : int32, Grid_yCellWidth : double, Grid_yNum : int32,
                             Grid_zBnum : int32, Grid_zCellWidth : double, Grid_zNum : int32)
where
  reads(Fluid.{rhoFluxX, rhoFluxY, rhoFluxZ}),
  reads(Fluid.{rhoVelocityFluxX, rhoVelocityFluxY, rhoVelocityFluxZ}),
  reads(Fluid.{rhoEnergyFluxX, rhoEnergyFluxY, rhoEnergyFluxZ}),
  reads writes(Fluid.{rho_t, rhoVelocity_t, rhoEnergy_t})
do
  __demand(__openmp)
  for c in Fluid do
    -- If interior cell
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then

      var stencil1 = ((c+{-1, 0, 0})%Fluid.bounds)
      var stencil2 = ((c+{0, -1, 0})%Fluid.bounds)
      var stencil3 = ((c+{0, 0, -1})%Fluid.bounds)

      Fluid[c].rho_t += ((-(Fluid[c].rhoFluxX-Fluid[stencil1].rhoFluxX))/Grid_xCellWidth)
      var tmp = vs_div(vs_mul(vv_sub(Fluid[c].rhoVelocityFluxX, Fluid[stencil1].rhoVelocityFluxX), double((-1))), Grid_xCellWidth)
      var v = Fluid[c].rhoVelocity_t
      v[0] += tmp[0]
      v[1] += tmp[1]
      v[2] += tmp[2]
      Fluid[c].rhoVelocity_t = v
      Fluid[c].rhoEnergy_t += ((-(Fluid[c].rhoEnergyFluxX-Fluid[stencil1].rhoEnergyFluxX))/Grid_xCellWidth)

      Fluid[c].rho_t += ((-(Fluid[c].rhoFluxY-Fluid[stencil2].rhoFluxY))/Grid_yCellWidth)
      var tmp__7144 = vs_div(vs_mul(vv_sub(Fluid[c].rhoVelocityFluxY, Fluid[stencil2].rhoVelocityFluxY), double((-1))), Grid_yCellWidth)
      var v__7145 = Fluid[c].rhoVelocity_t
      v__7145[0] += tmp__7144[0]
      v__7145[1] += tmp__7144[1]
      v__7145[2] += tmp__7144[2]
      Fluid[c].rhoVelocity_t = v__7145
      Fluid[c].rhoEnergy_t += ((-(Fluid[c].rhoEnergyFluxY-Fluid[stencil2].rhoEnergyFluxY))/Grid_yCellWidth)

      Fluid[c].rho_t += ((-(Fluid[c].rhoFluxZ-Fluid[stencil3].rhoFluxZ))/Grid_zCellWidth)
      var tmp__7146 = vs_div(vs_mul(vv_sub(Fluid[c].rhoVelocityFluxZ, Fluid[stencil3].rhoVelocityFluxZ), double((-1))), Grid_zCellWidth)
      var v__7147 = Fluid[c].rhoVelocity_t
      v__7147[0] += tmp__7146[0]
      v__7147[1] += tmp__7146[1]
      v__7147[2] += tmp__7146[2]
      Fluid[c].rhoVelocity_t = v__7147
      Fluid[c].rhoEnergy_t += ((-(Fluid[c].rhoEnergyFluxZ-Fluid[stencil3].rhoEnergyFluxZ))/Grid_zCellWidth)

    end
  end
end

__demand(__parallel, __cuda)
task Flow_AddGetFluxGhostNSCBC(Fluid : region(ispace(int3d), Fluid_columns),
                               config : Config,
                               Flow_constantVisc : double,
                               Flow_gamma : double,
                               Flow_gasConstant : double,
                               Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                               Flow_prandtl : double,
                               Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                               Flow_viscosityModel : SCHEMA.ViscosityModel,
                               Grid_xBnum : int32, Grid_xCellWidth : double, Grid_xNum : int32,
                               Grid_yBnum : int32, Grid_yCellWidth : double, Grid_yNum : int32,
                               Grid_zBnum : int32, Grid_zCellWidth : double, Grid_zNum : int32)
where
  reads(Fluid.{rho, pressure, velocity, rhoVelocity, rhoEnthalpy, temperature}),
  reads(Fluid.{velocityGradientX, velocityGradientY, velocityGradientZ}),
  reads writes(Fluid.{rhoEnergyFluxY, rhoEnergyFluxZ}),
  reads writes(Fluid.{rhoFluxY, rhoFluxZ}),
  reads writes(Fluid.{rhoVelocityFluxY, rhoVelocityFluxZ})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0)
    var xPosGhost  = (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)
    var yNegGhost = (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)
    var yPosGhost  = (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)
    var zNegGhost = (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)
    var zPosGhost  = (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)
    var ghost_cell = (xNegGhost or xPosGhost or
                      yNegGhost or yPosGhost or
                      zNegGhost or zPosGhost )
    var interior_cell = not (ghost_cell)
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if NSCBC_inflow_cell or NSCBC_outflow_cell  then
      -- y fluxes
      var flux = CenteredInviscidFlux_(int3d(c), ((c+{0, 1, 0})%Fluid.bounds), Fluid)
      Fluid[c].rhoFluxY = flux[0]
      Fluid[c].rhoVelocityFluxY = array(flux[1], flux[2], flux[3])
      Fluid[c].rhoEnergyFluxY = flux[4]
      var muFace = (0.5*(GetDynamicViscosity(Fluid[c].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)+GetDynamicViscosity(Fluid[((c+{0, 1, 0})%Fluid.bounds)].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)))
      var velocityFace = array(0.0, 0.0, 0.0)
      var velocityY_XFace = 0.0
      var velocityY_ZFace = 0.0
      var velocityX_XFace = 0.0
      var velocityZ_ZFace = 0.0
      velocityFace = vs_mul(vv_add(Fluid[c].velocity, Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocity), 0.5)
      velocityY_XFace = (0.5*(Fluid[c].velocityGradientX[1]+Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocityGradientX[1]))
      velocityY_ZFace = (0.5*(Fluid[c].velocityGradientZ[1]+Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocityGradientZ[1]))
      velocityX_XFace = (0.5*(Fluid[c].velocityGradientX[0]+Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocityGradientX[0]))
      velocityZ_ZFace = (0.5*(Fluid[c].velocityGradientZ[2]+Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocityGradientZ[2]))
      var velocityX_YFace = 0.0
      var velocityY_YFace = 0.0
      var velocityZ_YFace = 0.0
      var temperature_YFace = 0.0
      velocityX_YFace = (0.5*(Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocity[0]-Fluid[c].velocity[0]))
      velocityY_YFace = (0.5*(Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocity[1]-Fluid[c].velocity[1]))
      velocityZ_YFace = (0.5*(Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocity[2]-Fluid[c].velocity[2]))
      temperature_YFace = (0.5*(Fluid[((c+{0, 1, 0})%Fluid.bounds)].temperature-Fluid[c].temperature))
      velocityX_YFace *= (1/(Grid_yCellWidth*0.5))
      velocityY_YFace *= (1/(Grid_yCellWidth*0.5))
      velocityZ_YFace *= (1/(Grid_yCellWidth*0.5))
      temperature_YFace *= (1/(Grid_yCellWidth*0.5))
      var sigmaXY = (muFace*(velocityX_YFace+velocityY_XFace))
      var sigmaYY = ((muFace*(((4.0*velocityY_YFace)-(2.0*velocityX_XFace))-(2.0*velocityZ_ZFace)))/3.0)
      var sigmaZY = (muFace*(velocityZ_YFace+velocityY_ZFace))
      var usigma = (((velocityFace[0]*sigmaXY)+(velocityFace[1]*sigmaYY))+(velocityFace[2]*sigmaZY))
      var cp = ((Flow_gamma*Flow_gasConstant)/(Flow_gamma-1.0))
      var heatFlux = ((-((cp*muFace)/Flow_prandtl))*temperature_YFace)
      Fluid[c].rhoVelocityFluxY[0] += (-sigmaXY)
      Fluid[c].rhoVelocityFluxY[1] += (-sigmaYY)
      Fluid[c].rhoVelocityFluxY[2] += (-sigmaZY)
      Fluid[c].rhoEnergyFluxY += (-(usigma-heatFlux))
    end
    if (NSCBC_inflow_cell or NSCBC_outflow_cell) then
      -- z fluxes
      var flux = CenteredInviscidFlux__(int3d(c), ((c+{0, 0, 1})%Fluid.bounds), Fluid)
      Fluid[c].rhoFluxZ = flux[0]
      Fluid[c].rhoVelocityFluxZ = array(flux[1], flux[2], flux[3])
      Fluid[c].rhoEnergyFluxZ = flux[4]
      var muFace = (0.5*(GetDynamicViscosity(Fluid[c].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)+GetDynamicViscosity(Fluid[((c+{0, 0, 1})%Fluid.bounds)].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)))
      var velocityFace = array(0.0, 0.0, 0.0)
      var velocityZ_XFace = 0.0
      var velocityZ_YFace = 0.0
      var velocityX_XFace = 0.0
      var velocityY_YFace = 0.0
      velocityFace = vs_mul(vv_add(Fluid[c].velocity, Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocity), 0.5)
      velocityZ_XFace = (0.5*(Fluid[c].velocityGradientX[2]+Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocityGradientX[2]))
      velocityZ_YFace = (0.5*(Fluid[c].velocityGradientY[2]+Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocityGradientY[2]))
      velocityX_XFace = (0.5*(Fluid[c].velocityGradientX[0]+Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocityGradientX[0]))
      velocityY_YFace = (0.5*(Fluid[c].velocityGradientY[1]+Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocityGradientY[1]))
      var velocityX_ZFace = 0.0
      var velocityY_ZFace = 0.0
      var velocityZ_ZFace = 0.0
      var temperature_ZFace = 0.0
      velocityX_ZFace = (0.5*(Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocity[0]-Fluid[c].velocity[0]))
      velocityY_ZFace = (0.5*(Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocity[1]-Fluid[c].velocity[1]))
      velocityZ_ZFace = (0.5*(Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocity[2]-Fluid[c].velocity[2]))
      temperature_ZFace = (0.5*(Fluid[((c+{0, 0, 1})%Fluid.bounds)].temperature-Fluid[c].temperature))
      velocityX_ZFace *= (1/(Grid_zCellWidth*0.5))
      velocityY_ZFace *= (1/(Grid_zCellWidth*0.5))
      velocityZ_ZFace *= (1/(Grid_zCellWidth*0.5))
      temperature_ZFace *= (1/(Grid_zCellWidth*0.5))
      var sigmaXZ = (muFace*(velocityX_ZFace+velocityZ_XFace))
      var sigmaYZ = (muFace*(velocityY_ZFace+velocityZ_YFace))
      var sigmaZZ = ((muFace*(((4.0*velocityZ_ZFace)-(2.0*velocityX_XFace))-(2.0*velocityY_YFace)))/3.0)
      var usigma = (((velocityFace[0]*sigmaXZ)+(velocityFace[1]*sigmaYZ))+(velocityFace[2]*sigmaZZ))
      var cp = ((Flow_gamma*Flow_gasConstant)/(Flow_gamma-1.0))
      var heatFlux = ((-((cp*muFace)/Flow_prandtl))*temperature_ZFace)
      Fluid[c].rhoVelocityFluxZ[0] += (-sigmaXZ)
      Fluid[c].rhoVelocityFluxZ[1] += (-sigmaYZ)
      Fluid[c].rhoVelocityFluxZ[2] += (-sigmaZZ)
      Fluid[c].rhoEnergyFluxZ += (-(usigma-heatFlux))
    end
  end
end

__demand(__parallel, __cuda)
task Flow_AddUpdateUsingFluxGhostNSCBC(Fluid : region(ispace(int3d), Fluid_columns),
                                       config : Config,
                                       Flow_gamma : double, Flow_gasConstant : double,
                                       Flow_prandtl : double,
                                       maxMach : double,
                                       Flow_lengthScale : double,
                                       Flow_constantVisc : double,
                                       Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                                       Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                                       Flow_viscosityModel : SCHEMA.ViscosityModel,
                                       BC_xPosP_inf : double,
                                       Grid_xBnum : int32, Grid_xCellWidth : double, Grid_xNum : int32,
                                       Grid_yBnum : int32, Grid_yCellWidth : double, Grid_yNum : int32,
                                       Grid_zBnum : int32, Grid_zCellWidth : double, Grid_zNum : int32)
where
  reads(Fluid.{rho, velocity, pressure, temperature, rhoVelocity, velocityGradientX, velocityGradientY, velocityGradientZ, rhoEnergy, dudtBoundary, dTdtBoundary}),
  reads(Fluid.{rhoFluxX, rhoFluxY, rhoFluxZ}),
  reads(Fluid.{rhoVelocityFluxX, rhoVelocityFluxY, rhoVelocityFluxZ}),
  reads(Fluid.{rhoEnergyFluxX, rhoEnergyFluxY, rhoEnergyFluxZ}),
  reads writes(Fluid.{rho_t, rhoVelocity_t, rhoEnergy_t})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0)
    var xPosGhost  = (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)
    var yNegGhost = (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)
    var yPosGhost  = (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)
    var zNegGhost = (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)
    var zPosGhost  = (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)
    var ghost_cell = (xNegGhost or xPosGhost or
                      yNegGhost or yPosGhost or
                      zNegGhost or zPosGhost )
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if ghost_cell then
      if NSCBC_inflow_cell then
        -- add y and z fluxes for inflow cells
        Fluid[c].rho_t += ((-(Fluid[c].rhoFluxY-Fluid[((c+{0, -1, 0})%Fluid.bounds)].rhoFluxY))/Grid_yCellWidth)
        Fluid[c].rho_t += ((-(Fluid[c].rhoFluxZ-Fluid[((c+{0, 0, -1})%Fluid.bounds)].rhoFluxZ))/Grid_zCellWidth)

        -- Add in the x flux using NSCBC
        var c_bnd = int3d(c)
        var c_int = ((c+{1, 0, 0})%Fluid.bounds)

        -- compute amplitudes of waves leaving the domain
        var c_sound = GetSoundSpeed(Fluid[c_bnd].temperature, Flow_gamma, Flow_gasConstant)
        var lambda_1 = Fluid[c_bnd].velocity[0] - c_sound
        var dP_dx = (Fluid[c_int].pressure    - Fluid[c_bnd].pressure)    /  Grid_xCellWidth
        var du_dx = (Fluid[c_int].velocity[0] - Fluid[c_bnd].velocity[0]) /  Grid_xCellWidth
        var L1 = lambda_1*(dP_dx - Fluid[c_bnd].rho*c_sound*du_dx)

        -- compute amplitudes of waves entering the domain
        var L5 = L1 - 2*Fluid[c_bnd].rho*c_sound*Fluid[c_bnd].dudtBoundary
        var L2 = 0.5*(Flow_gamma - 1.0)*(L5+L1) + (Fluid[c_bnd].rho*c_sound*c_sound/Fluid[c_bnd].temperature)*Fluid[c_bnd].dTdtBoundary

        -- update RHS of transport equation for boundary cell
        var d1 = 1/(c_sound*c_sound)*(L2+0.5*(L1+L5))

        -- Set RHS to update the density in the ghost inflow cells
        Fluid[c_bnd].rho_t += - d1
      end

      if NSCBC_outflow_cell then
        -- update y and z fluxes for outflow cells
        Fluid[c].rho_t += ((-(Fluid[c].rhoFluxY-Fluid[((c+{0, -1, 0})%Fluid.bounds)].rhoFluxY))/Grid_yCellWidth)
        var tmp__7144 = vs_div(vs_mul(vv_sub(Fluid[c].rhoVelocityFluxY, Fluid[((c+{0, -1, 0})%Fluid.bounds)].rhoVelocityFluxY), double((-1))), Grid_yCellWidth)
        var v__7145 = Fluid[c].rhoVelocity_t
        v__7145[0] += tmp__7144[0]
        v__7145[1] += tmp__7144[1]
        v__7145[2] += tmp__7144[2]
        Fluid[c].rhoVelocity_t = v__7145
        Fluid[c].rhoEnergy_t += ((-(Fluid[c].rhoEnergyFluxY-Fluid[((c+{0, -1, 0})%Fluid.bounds)].rhoEnergyFluxY))/Grid_yCellWidth)

        Fluid[c].rho_t += ((-(Fluid[c].rhoFluxZ-Fluid[((c+{0, 0, -1})%Fluid.bounds)].rhoFluxZ))/Grid_zCellWidth)
        var tmp__7146 = vs_div(vs_mul(vv_sub(Fluid[c].rhoVelocityFluxZ, Fluid[((c+{0, 0, -1})%Fluid.bounds)].rhoVelocityFluxZ), double((-1))), Grid_zCellWidth)
        var v__7147 = Fluid[c].rhoVelocity_t
        v__7147[0] += tmp__7146[0]
        v__7147[1] += tmp__7146[1]
        v__7147[2] += tmp__7146[2]
        Fluid[c].rhoVelocity_t = v__7147
        Fluid[c].rhoEnergy_t += ((-(Fluid[c].rhoEnergyFluxZ-Fluid[((c+{0, 0, -1})%Fluid.bounds)].rhoEnergyFluxZ))/Grid_zCellWidth)

        -- Add in the x fluxes using NSCBC for outflow
        var c_bnd = int3d(c)
        var c_int = ((c+{-1, 0, 0})%Fluid.bounds)

        var sigma = 0.25 -- Specified constant
        var c_sound = GetSoundSpeed(Fluid[c_bnd].temperature, Flow_gamma, Flow_gasConstant)
        var K = sigma*(1.0-maxMach*maxMach)*c_sound/Flow_lengthScale

        var L1 = K*(Fluid[c_bnd].pressure - BC_xPosP_inf)

        var lambda_2 = Fluid[c_bnd].velocity[0]
        var lambda_3 = Fluid[c_bnd].velocity[0]
        var lambda_4 = Fluid[c_bnd].velocity[0]
        var lambda_5 = Fluid[c_bnd].velocity[0] + c_sound

        var drho_dx = (Fluid[c_bnd].rho - Fluid[c_int].rho) /  Grid_xCellWidth
        var dp_dx   = (Fluid[c_bnd].pressure    - Fluid[c_int].pressure   ) /  Grid_xCellWidth
        var du_dx   = (Fluid[c_bnd].velocity[0] - Fluid[c_int].velocity[0]) /  Grid_xCellWidth
        var dv_dx   = (Fluid[c_bnd].velocity[1] - Fluid[c_int].velocity[1]) /  Grid_xCellWidth
        var dw_dx   = (Fluid[c_bnd].velocity[2] - Fluid[c_int].velocity[2]) /  Grid_xCellWidth

        var L2 = lambda_2*(c_sound*c_sound*drho_dx - dp_dx)
        var L3 = lambda_3*(dv_dx)
        var L4 = lambda_4*(dw_dx)
        var L5 = lambda_5*(dp_dx + Fluid[c_bnd].rho*c_sound*du_dx)

        var d1 = 1.0/(c_sound*c_sound)*(L2 + 0.5*(L5 + L1))
        var d2 = 0.5*(L5 + L1)
        var d3 = 1.0/(2.0*Fluid[c_bnd].rho*c_sound)*(L5 - L1)
        var d4 = L3
        var d5 = L4

        var mu_pos = GetDynamicViscosity(Fluid[c_bnd].temperature,
                                        Flow_constantVisc,
                                        Flow_powerlawTempRef, Flow_powerlawViscRef,
                                        Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef,
                                        Flow_viscosityModel)
        var tau11_pos = mu_pos*( Fluid[c_bnd].velocityGradientX[0] + Fluid[c_bnd].velocityGradientX[0] - (2.0/3.0)*(Fluid[c_bnd].velocityGradientX[0] + Fluid[c_bnd].velocityGradientY[1] + Fluid[c_bnd].velocityGradientZ[2]) )
        var tau21_pos = mu_pos*( Fluid[c_bnd].velocityGradientX[1] + Fluid[c_bnd].velocityGradientY[0] )
        var tau31_pos = mu_pos*( Fluid[c_bnd].velocityGradientX[2] + Fluid[c_bnd].velocityGradientZ[0] )

        var mu_neg = GetDynamicViscosity(Fluid[c_int].temperature,
                                        Flow_constantVisc,
                                        Flow_powerlawTempRef, Flow_powerlawViscRef,
                                        Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef,
                                        Flow_viscosityModel)
        var tau11_neg = mu_neg*( Fluid[c_int].velocityGradientX[0] + Fluid[c_int].velocityGradientX[0] - (2.0/3.0)*(Fluid[c_int].velocityGradientX[0] + Fluid[c_int].velocityGradientY[1] + Fluid[c_int].velocityGradientZ[2]) )
        var tau21_neg = mu_neg*( Fluid[c_int].velocityGradientX[1] + Fluid[c_int].velocityGradientY[0] )
        var tau31_neg = mu_neg*( Fluid[c_int].velocityGradientX[2] + Fluid[c_int].velocityGradientZ[0] )

        -- Stuff for momentum equations
        var dtau11_dx = (tau11_pos - tau11_neg) / (Grid_xCellWidth)
        var dtau21_dx = (tau21_pos - tau21_neg) / (Grid_xCellWidth)
        var dtau31_dx = (tau31_pos - tau31_neg) / (Grid_xCellWidth)

        -- Stuff for energy equation
        var mu = GetDynamicViscosity(Fluid[c_bnd].temperature,
                                     Flow_constantVisc,
                                     Flow_powerlawTempRef, Flow_powerlawViscRef,
                                     Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef,
                                     Flow_viscosityModel)
        var tau_12 =  mu*( Fluid[c_bnd].velocityGradientY[0] + Fluid[c_bnd].velocityGradientX[1] )
        var tau_13 =  mu*( Fluid[c_bnd].velocityGradientZ[0] + Fluid[c_bnd].velocityGradientX[2] )
        var energy_term_x = (Fluid[c_bnd].velocity[0]*tau11_pos - Fluid[c_int].velocity[0]*tau11_neg) / (Grid_xCellWidth) + c.velocityGradientX[1]*tau_12 + c.velocityGradientX[2]*tau_13

        -- Update the RHS of conservation equaions with x fluxes
        Fluid[c_bnd].rho_t += - d1
        Fluid[c_bnd].rhoVelocity_t[0] += -Fluid[c_bnd].velocity[0]*d1 - Fluid[c_bnd].rho*d3 + dtau11_dx
        Fluid[c_bnd].rhoVelocity_t[1] += -Fluid[c_bnd].velocity[1]*d1 - Fluid[c_bnd].rho*d4 + dtau21_dx
        Fluid[c_bnd].rhoVelocity_t[2] += -Fluid[c_bnd].velocity[2]*d1 - Fluid[c_bnd].rho*d5 + dtau31_dx
        Fluid[c_bnd].rhoEnergy_t += -0.5*(Fluid[c_bnd].velocity[0]*Fluid[c_bnd].velocity[0] + Fluid[c_bnd].velocity[1]*Fluid[c_bnd].velocity[1] + Fluid[c_bnd].velocity[2]*Fluid[c_bnd].velocity[2])*d1 - d2/(Flow_gamma-1.0) - Fluid[c_bnd].rhoVelocity[0]*d3 - Fluid[c_bnd].rhoVelocity[1]*d4 - Fluid[c_bnd].rhoVelocity[2]*d5 + energy_term_x
      end
    end
  end
end

-- Update the time derivative values needed for subsonic inflow
__demand(__parallel, __cuda)
task Flow_UpdateNSCBCGhostCellTimeDerivatives(Fluid : region(ispace(int3d), Fluid_columns),
                                              config : Config,
                                              Grid_xBnum : int32, Grid_xNum : int32,
                                              Grid_yBnum : int32, Grid_yNum : int32,
                                              Grid_zBnum : int32, Grid_zNum : int32,
                                              Integrator_deltaTime : double)
where
  reads(Fluid.{velocity, temperature}),
  writes(Fluid.{dudtBoundary, dTdtBoundary}),
  reads writes(Fluid.{velocity_old_NSCBC, temperature_old_NSCBC})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0)
    var xPosGhost  = (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)
    var yNegGhost = (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)
    var yPosGhost  = (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)
    var zNegGhost = (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)
    var zPosGhost  = (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)
    var ghost_cell = (xNegGhost or xPosGhost or
                      yNegGhost or yPosGhost or
                      zNegGhost or zPosGhost )
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if NSCBC_inflow_cell then
      Fluid[c].dudtBoundary = (Fluid[c].velocity[0] - Fluid[c].velocity_old_NSCBC[0]) / Integrator_deltaTime
      Fluid[c].dTdtBoundary = (Fluid[c].temperature - Fluid[c].temperature_old_NSCBC) / Integrator_deltaTime

      Fluid[c].velocity_old_NSCBC    = Fluid[c].velocity
      Fluid[c].temperature_old_NSCBC = Fluid[c].temperature
    end

  end
end

__demand(__parallel, __cuda)
task Flow_AddBodyForces(Fluid : region(ispace(int3d), Fluid_columns),
                        Flow_bodyForce : double[3],
                        Grid_xBnum : int32, Grid_xNum : int32,
                        Grid_yBnum : int32, Grid_yNum : int32,
                        Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{rho, velocity}),
  reads writes(Fluid.{rhoEnergy_t, rhoVelocity_t})
do
  __demand(__openmp)
  for c in Fluid do
    -- if interior cell
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      var tmp = vs_mul(Flow_bodyForce, Fluid[c].rho)
      var v = Fluid[c].rhoVelocity_t
      v[0] += tmp[0]
      v[1] += tmp[1]
      v[2] += tmp[2]
      Fluid[c].rhoVelocity_t = v
      Fluid[c].rhoEnergy_t += (Fluid[c].rho*dot(Flow_bodyForce, Fluid[c].velocity))
    end
  end
end

__demand(__parallel, __cuda)
task Flow_AddBodyForcesGhostNSCBC(Fluid : region(ispace(int3d), Fluid_columns),
                                  config : Config,
                                  Flow_bodyForce : double[3],
                                  Grid_xBnum : int32, Grid_xNum : int32,
                                  Grid_yBnum : int32, Grid_yNum : int32,
                                  Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{rho, velocity}),
  reads writes(Fluid.{rhoEnergy_t, rhoVelocity_t})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0)
    var xPosGhost  = (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)
    var yNegGhost = (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)
    var yPosGhost  = (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)
    var zNegGhost = (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)
    var zPosGhost  = (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if NSCBC_inflow_cell or NSCBC_outflow_cell then
      var tmp = vs_mul(Flow_bodyForce, Fluid[c].rho)
      var v = Fluid[c].rhoVelocity_t
      v[0] += tmp[0]
      v[1] += tmp[1]
      v[2] += tmp[2]
      Fluid[c].rhoVelocity_t = v
      Fluid[c].rhoEnergy_t += (Fluid[c].rho*dot(Flow_bodyForce, Fluid[c].velocity))
    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdatePD(Fluid : region(ispace(int3d), Fluid_columns),
                   Grid_xBnum : int32, Grid_xNum : int32,
                   Grid_yBnum : int32, Grid_yNum : int32,
                   Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{pressure, velocityGradientX, velocityGradientY, velocityGradientZ}),
  writes(Fluid.PD)
do
  __demand(__openmp)
  for c in Fluid do
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      var divU = 0.0
      divU = ((Fluid[c].velocityGradientX[0]+Fluid[c].velocityGradientY[1])+Fluid[c].velocityGradientZ[2])
      Fluid[c].PD = (divU*Fluid[c].pressure)
    end
  end
end

__demand(__parallel, __cuda)
task Flow_ResetDissipation(Fluid : region(ispace(int3d), Fluid_columns))
where
  writes(Fluid.dissipation)
do
  __demand(__openmp)
  for c in Fluid do
    Fluid[c].dissipation = 0.0
  end
end

__demand(__parallel, __cuda)
task Flow_ComputeDissipationX(Fluid : region(ispace(int3d), Fluid_columns),
                              Flow_constantVisc : double,
                              Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                              Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                              Flow_viscosityModel : SCHEMA.ViscosityModel,
                              Grid_xBnum : int32, Grid_xNum : int32, Grid_xCellWidth : double,
                              Grid_yBnum : int32, Grid_yNum : int32,
                              Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{velocity, temperature, velocityGradientY, velocityGradientZ}),
  writes(Fluid.dissipationFlux)
do
  __demand(__openmp)
  for c in Fluid do
    if ((not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) or (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)==1)) then
      var muFace = (0.5*(GetDynamicViscosity(Fluid[c].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)+GetDynamicViscosity(Fluid[((c+{1, 0, 0})%Fluid.bounds)].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)))
      var velocityFace = array(0.0, 0.0, 0.0)
      var velocityX_YFace = 0.0
      var velocityX_ZFace = 0.0
      var velocityY_YFace = 0.0
      var velocityZ_ZFace = 0.0
      velocityFace = vs_mul(vv_add(Fluid[c].velocity, Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocity), 0.5)
      velocityX_YFace = (0.5*(Fluid[c].velocityGradientY[0]+Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocityGradientY[0]))
      velocityX_ZFace = (0.5*(Fluid[c].velocityGradientZ[0]+Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocityGradientZ[0]))
      velocityY_YFace = (0.5*(Fluid[c].velocityGradientY[1]+Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocityGradientY[1]))
      velocityZ_ZFace = (0.5*(Fluid[c].velocityGradientZ[2]+Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocityGradientZ[2]))
      var velocityX_XFace = 0.0
      var velocityY_XFace = 0.0
      var velocityZ_XFace = 0.0
      var temperature_XFace = 0.0
      velocityX_XFace = (0.5*(Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocity[0]-Fluid[c].velocity[0]))
      velocityY_XFace = (0.5*(Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocity[1]-Fluid[c].velocity[1]))
      velocityZ_XFace = (0.5*(Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocity[2]-Fluid[c].velocity[2]))
      temperature_XFace = (0.5*(Fluid[((c+{1, 0, 0})%Fluid.bounds)].temperature-Fluid[c].temperature))
      velocityX_XFace *= (1/(Grid_xCellWidth*0.5))
      velocityY_XFace *= (1/(Grid_xCellWidth*0.5))
      velocityZ_XFace *= (1/(Grid_xCellWidth*0.5))
      temperature_XFace *= (1/(Grid_xCellWidth*0.5))
      var sigmaXX = ((muFace*(((4.0*velocityX_XFace)-(2.0*velocityY_YFace))-(2.0*velocityZ_ZFace)))/3.0)
      var sigmaYX = (muFace*(velocityY_XFace+velocityX_YFace))
      var sigmaZX = (muFace*(velocityZ_XFace+velocityX_ZFace))
      var usigma = (((velocityFace[0]*sigmaXX)+(velocityFace[1]*sigmaYX))+(velocityFace[2]*sigmaZX))
      Fluid[c].dissipationFlux = usigma
    end
  end
end

-- CHANGE to reduces
__demand(__parallel, __cuda)
task Flow_UpdateDissipationX(Fluid : region(ispace(int3d), Fluid_columns),
                             Grid_xBnum : int32, Grid_xNum : int32, Grid_xCellWidth : double,
                             Grid_yBnum : int32, Grid_yNum : int32,
                             Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.dissipationFlux),
  reads writes(Fluid.dissipation)
do
  __demand(__openmp)
  for c in Fluid do
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      Fluid[c].dissipation += ((Fluid[c].dissipationFlux-Fluid[((c+{-1, 0, 0})%Fluid.bounds)].dissipationFlux)/Grid_xCellWidth)
    end
  end
end

__demand(__parallel, __cuda)
task Flow_ComputeDissipationY(Fluid : region(ispace(int3d), Fluid_columns),
                              Flow_constantVisc : double,
                              Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                              Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                              Flow_viscosityModel : SCHEMA.ViscosityModel,
                              Grid_xBnum : int32, Grid_xNum : int32,
                              Grid_yBnum : int32, Grid_yNum : int32, Grid_yCellWidth : double,
                              Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{velocity, temperature, velocityGradientX, velocityGradientZ}),
  writes(Fluid.dissipationFlux)
do
  __demand(__openmp)
  for c in Fluid do
    if ((not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)==1)) then
      var muFace = (0.5*(GetDynamicViscosity(Fluid[c].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)+GetDynamicViscosity(Fluid[((c+{0, 1, 0})%Fluid.bounds)].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)))
      var velocityFace = array(0.0, 0.0, 0.0)
      var velocityY_XFace = 0.0
      var velocityY_ZFace = 0.0
      var velocityX_XFace = 0.0
      var velocityZ_ZFace = 0.0
      velocityFace = vs_mul(vv_add(Fluid[c].velocity, Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocity), 0.5)
      velocityY_XFace = (0.5*(Fluid[c].velocityGradientX[1]+Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocityGradientX[1]))
      velocityY_ZFace = (0.5*(Fluid[c].velocityGradientZ[1]+Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocityGradientZ[1]))
      velocityX_XFace = (0.5*(Fluid[c].velocityGradientX[0]+Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocityGradientX[0]))
      velocityZ_ZFace = (0.5*(Fluid[c].velocityGradientZ[2]+Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocityGradientZ[2]))
      var velocityX_YFace = 0.0
      var velocityY_YFace = 0.0
      var velocityZ_YFace = 0.0
      var temperature_YFace = 0.0
      velocityX_YFace = (0.5*(Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocity[0]-Fluid[c].velocity[0]))
      velocityY_YFace = (0.5*(Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocity[1]-Fluid[c].velocity[1]))
      velocityZ_YFace = (0.5*(Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocity[2]-Fluid[c].velocity[2]))
      temperature_YFace = (0.5*(Fluid[((c+{0, 1, 0})%Fluid.bounds)].temperature-Fluid[c].temperature))
      velocityX_YFace *= (1/(Grid_yCellWidth*0.5))
      velocityY_YFace *= (1/(Grid_yCellWidth*0.5))
      velocityZ_YFace *= (1/(Grid_yCellWidth*0.5))
      temperature_YFace *= (1/(Grid_yCellWidth*0.5))
      var sigmaXY = (muFace*(velocityX_YFace+velocityY_XFace))
      var sigmaYY = ((muFace*(((4.0*velocityY_YFace)-(2.0*velocityX_XFace))-(2.0*velocityZ_ZFace)))/3.0)
      var sigmaZY = (muFace*(velocityZ_YFace+velocityY_ZFace))
      var usigma = (((velocityFace[0]*sigmaXY)+(velocityFace[1]*sigmaYY))+(velocityFace[2]*sigmaZY))
      Fluid[c].dissipationFlux = usigma
    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdateDissipationY(Fluid : region(ispace(int3d), Fluid_columns),
                             Grid_xBnum : int32, Grid_xNum : int32,
                             Grid_yBnum : int32, Grid_yNum : int32, Grid_yCellWidth : double,
                             Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.dissipationFlux),
  reads writes(Fluid.dissipation)
do
  __demand(__openmp)
  for c in Fluid do
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      Fluid[c].dissipation += ((Fluid[c].dissipationFlux-Fluid[((c+{0, -1, 0})%Fluid.bounds)].dissipationFlux)/Grid_yCellWidth)
    end
  end
end

__demand(__parallel, __cuda)
task Flow_ComputeDissipationZ(Fluid : region(ispace(int3d), Fluid_columns),
                              Flow_constantVisc : double,
                              Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                              Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                              Flow_viscosityModel : SCHEMA.ViscosityModel,
                              Grid_xBnum : int32, Grid_xNum : int32,
                              Grid_yBnum : int32, Grid_yNum : int32,
                              Grid_zBnum : int32, Grid_zNum : int32, Grid_zCellWidth : double)
where
  reads(Fluid.{velocity, temperature, velocityGradientX, velocityGradientY}),
  writes(Fluid.dissipationFlux)
do
  __demand(__openmp)
  for c in Fluid do
    if ((not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)==1)) then
      var muFace = (0.5*(GetDynamicViscosity(Fluid[c].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)+GetDynamicViscosity(Fluid[((c+{0, 0, 1})%Fluid.bounds)].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)))
      var velocityFace = array(0.0, 0.0, 0.0)
      var velocityZ_XFace = 0.0
      var velocityZ_YFace = 0.0
      var velocityX_XFace = 0.0
      var velocityY_YFace = 0.0
      velocityFace = vs_mul(vv_add(Fluid[c].velocity, Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocity), 0.5)
      velocityZ_XFace = (0.5*(Fluid[c].velocityGradientX[2]+Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocityGradientX[2]))
      velocityZ_YFace = (0.5*(Fluid[c].velocityGradientY[2]+Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocityGradientY[2]))
      velocityX_XFace = (0.5*(Fluid[c].velocityGradientX[0]+Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocityGradientX[0]))
      velocityY_YFace = (0.5*(Fluid[c].velocityGradientY[1]+Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocityGradientY[1]))
      var velocityX_ZFace = 0.0
      var velocityY_ZFace = 0.0
      var velocityZ_ZFace = 0.0
      var temperature_ZFace = 0.0
      velocityX_ZFace = (0.5*(Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocity[0]-Fluid[c].velocity[0]))
      velocityY_ZFace = (0.5*(Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocity[1]-Fluid[c].velocity[1]))
      velocityZ_ZFace = (0.5*(Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocity[2]-Fluid[c].velocity[2]))
      temperature_ZFace = (0.5*(Fluid[((c+{0, 0, 1})%Fluid.bounds)].temperature-Fluid[c].temperature))
      velocityX_ZFace *= (1/(Grid_zCellWidth*0.5))
      velocityY_ZFace *= (1/(Grid_zCellWidth*0.5))
      velocityZ_ZFace *= (1/(Grid_zCellWidth*0.5))
      temperature_ZFace *= (1/(Grid_zCellWidth*0.5))
      var sigmaXZ = (muFace*(velocityX_ZFace+velocityZ_XFace))
      var sigmaYZ = (muFace*(velocityY_ZFace+velocityZ_YFace))
      var sigmaZZ = ((muFace*(((4.0*velocityZ_ZFace)-(2.0*velocityX_XFace))-(2.0*velocityY_YFace)))/3.0)
      var usigma = (((velocityFace[0]*sigmaXZ)+(velocityFace[1]*sigmaYZ))+(velocityFace[2]*sigmaZZ))
      Fluid[c].dissipationFlux = usigma
    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdateDissipationZ(Fluid : region(ispace(int3d), Fluid_columns),
                             Grid_xBnum : int32, Grid_xNum : int32,
                             Grid_yBnum : int32, Grid_yNum : int32,
                             Grid_zBnum : int32, Grid_zNum : int32, Grid_zCellWidth : double)
where
  reads(Fluid.dissipationFlux),
  reads writes(Fluid.dissipation)
do
  __demand(__openmp)
  for c in Fluid do
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      Fluid[c].dissipation += ((Fluid[c].dissipationFlux-Fluid[((c+{0, 0, -1})%Fluid.bounds)].dissipationFlux)/Grid_zCellWidth)
    end
  end
end

__demand(__parallel, __cuda)
task CalculateAveragePD(Fluid : region(ispace(int3d), Fluid_columns),
                        Grid_cellVolume : double,
                        Grid_xBnum : int32, Grid_xNum : int32,
                        Grid_yBnum : int32, Grid_yNum : int32,
                        Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.PD)
do
  var acc = 0.0
  __demand(__openmp)
  for c in Fluid do
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      acc += (Fluid[c].PD*Grid_cellVolume)
    end
  end
  return acc
end

__demand(__parallel, __cuda)
task CalculateAverageDissipation(Fluid : region(ispace(int3d), Fluid_columns),
                                 Grid_cellVolume : double,
                                 Grid_xBnum : int32, Grid_xNum : int32,
                                 Grid_yBnum : int32, Grid_yNum : int32,
                                 Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.dissipation)
do
  var acc = 0.0
  __demand(__openmp)
  for c in Fluid do
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      acc += (Fluid[c].dissipation*Grid_cellVolume)
    end
  end
  return acc
end

__demand(__parallel, __cuda)
task CalculateAverageK(Fluid : region(ispace(int3d), Fluid_columns),
                       Grid_cellVolume : double,
                       Grid_xBnum : int32, Grid_xNum : int32,
                       Grid_yBnum : int32, Grid_yNum : int32,
                       Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{rho, velocity})
do
  var acc = 0.0
  __demand(__openmp)
  for c in Fluid do
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      acc += (((0.5*Fluid[c].rho)*dot(Fluid[c].velocity, Fluid[c].velocity))*Grid_cellVolume)
    end
  end
  return acc
end

__demand(__parallel, __cuda)
task Flow_AddTurbulentSource(Fluid : region(ispace(int3d), Fluid_columns),
                             Flow_averageDissipation : double,
                             Flow_averageK : double,
                             Flow_averagePD : double,
                             Grid_cellVolume : double,
                             Grid_xBnum : int32, Grid_xNum : int32,
                             Grid_yBnum : int32, Grid_yNum : int32,
                             Grid_zBnum : int32, Grid_zNum : int32,
                             config : Config)
where
  reads(Fluid.{rho, velocity}),
  reads writes(Fluid.{rhoVelocity_t, rhoEnergy_t})
do
  var W = Flow_averagePD + Flow_averageDissipation
  var G = config.Flow.turbForcing.u.HIT.G
  var t_o = config.Flow.turbForcing.u.HIT.t_o
  var K_o = config.Flow.turbForcing.u.HIT.K_o
  var A = (-W-G*(Flow_averageK-K_o)/t_o) / (2.0*Flow_averageK)
  var acc = 0.0
  __demand(__openmp)
  for c in Fluid do
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      var force = vs_mul(Fluid[c].velocity, Fluid[c].rho*A)
      Fluid[c].rhoVelocity_t = vv_add(Fluid[c].rhoVelocity_t, force)
      Fluid[c].rhoEnergy_t += dot(force, Fluid[c].velocity)
      acc += dot(force, Fluid[c].velocity) * Grid_cellVolume
    end
  end
  return acc
end

-- CHANGE to reduces+?
__demand(__parallel, __cuda)
task Flow_AdjustTurbulentSource(Fluid : region(ispace(int3d), Fluid_columns),
                                Flow_averageFe : double,
                                Grid_xBnum : int32, Grid_xNum : int32,
                                Grid_yBnum : int32, Grid_yNum : int32,
                                Grid_zBnum : int32, Grid_zNum : int32)
where
  reads writes(Fluid.rhoEnergy_t)
do
  __demand(__openmp)
  for c in Fluid do
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      Fluid[c].rhoEnergy_t += (-Flow_averageFe)
    end
  end
end

-------------------------------------------------------------------------------
-- PARTICLE MOVEMENT
-------------------------------------------------------------------------------

__demand(__inline)
task locate(pos : double[3],
            Grid_xBnum : int32, Grid_xNum : int32, Grid_xOrigin : double, Grid_xWidth : double,
            Grid_yBnum : int32, Grid_yNum : int32, Grid_yOrigin : double, Grid_yWidth : double,
            Grid_zBnum : int32, Grid_zNum : int32, Grid_zOrigin : double, Grid_zWidth : double)
  var xcw = Grid_xWidth/Grid_xNum
  var xro = Grid_xOrigin-Grid_xBnum*xcw
  var xpos = floor((pos[0]-xro)/xcw)
  var xrnum = Grid_xNum+2*Grid_xBnum
  var xidx = max(0, min(xrnum-1, xpos))
  var ycw = Grid_yWidth/Grid_yNum
  var yro = Grid_yOrigin-Grid_yBnum*ycw
  var ypos = floor((pos[1]-yro)/ycw)
  var yrnum = Grid_yNum+2*Grid_yBnum
  var yidx = max(0, min(yrnum-1, ypos))
  var zcw = Grid_zWidth/Grid_zNum
  var zro = Grid_zOrigin-Grid_zBnum*zcw
  var zpos = floor((pos[2]-zro)/zcw)
  var zrnum = Grid_zNum+2*Grid_zBnum
  var zidx = max(0, min(zrnum-1, zpos))
  return int3d{xidx, yidx, zidx}
end

__demand(__cuda) -- MANUALLY PARALLELIZED
task Particles_LocateInCells(Particles : region(ispace(int1d), Particles_columns),
                             Grid_xBnum : int32, Grid_xNum : int32, Grid_xOrigin : double, Grid_xWidth : double,
                             Grid_yBnum : int32, Grid_yNum : int32, Grid_yOrigin : double, Grid_yWidth : double,
                             Grid_zBnum : int32, Grid_zNum : int32, Grid_zOrigin : double, Grid_zWidth : double)
where
  reads(Particles.{position, __valid}),
  writes(Particles.cell)
do
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      Particles[p].cell = locate(Particles[p].position,
                                 Grid_xBnum, Grid_xNum, Grid_xOrigin, Grid_xWidth,
                                 Grid_yBnum, Grid_yNum, Grid_yOrigin, Grid_yWidth,
                                 Grid_zBnum, Grid_zNum, Grid_zOrigin, Grid_zWidth)
    end
  end
end

__demand(__inline)
task Fluid_elemColor(idx : int3d,
                     Grid_xBnum : int32, Grid_xNum : int32, NX : int32,
                     Grid_yBnum : int32, Grid_yNum : int32, NY : int32,
                     Grid_zBnum : int32, Grid_zNum : int32, NZ : int32)
  idx.x = min(max(idx.x, Grid_xBnum), Grid_xNum+Grid_xBnum-1)
  idx.y = min(max(idx.y, Grid_yBnum), Grid_yNum+Grid_yBnum-1)
  idx.z = min(max(idx.z, Grid_zBnum), Grid_zNum+Grid_zBnum-1)
  return int3d{(idx.x-Grid_xBnum)/(Grid_xNum/NX),
               (idx.y-Grid_yBnum)/(Grid_yNum/NY),
               (idx.z-Grid_zBnum)/(Grid_zNum/NZ)}
end

-- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task Particles_CheckPartitioning(color : int3d,
                                 Particles : region(ispace(int1d), Particles_columns),
                                 Grid_xBnum : int32, Grid_xNum : int32, NX : int32,
                                 Grid_yBnum : int32, Grid_yNum : int32, NY : int32,
                                 Grid_zBnum : int32, Grid_zNum : int32, NZ : int32)
where
  reads(Particles.{cell, __valid})
do
  for p in Particles do
    if p.__valid then
      regentlib.assert(color == Fluid_elemColor(p.cell,
                                                Grid_xBnum, Grid_xNum, NX,
                                                Grid_yBnum, Grid_yNum, NY,
                                                Grid_zBnum, Grid_zNum, NZ),
                       'Invalid particle partitioning')
    end
  end
end

local colorOffsets = terralib.newlist({
  rexpr int3d({ 0,  0,  1}) end,
  rexpr int3d({ 0,  0, -1}) end,
  rexpr int3d({ 0,  1,  0}) end,
  rexpr int3d({ 0,  1,  1}) end,
  rexpr int3d({ 0,  1, -1}) end,
  rexpr int3d({ 0, -1,  0}) end,
  rexpr int3d({ 0, -1,  1}) end,
  rexpr int3d({ 0, -1, -1}) end,
  rexpr int3d({ 1,  0,  0}) end,
  rexpr int3d({ 1,  0,  1}) end,
  rexpr int3d({ 1,  0, -1}) end,
  rexpr int3d({ 1,  1,  0}) end,
  rexpr int3d({ 1,  1,  1}) end,
  rexpr int3d({ 1,  1, -1}) end,
  rexpr int3d({ 1, -1,  0}) end,
  rexpr int3d({ 1, -1,  1}) end,
  rexpr int3d({ 1, -1, -1}) end,
  rexpr int3d({-1,  0,  0}) end,
  rexpr int3d({-1,  0,  1}) end,
  rexpr int3d({-1,  0, -1}) end,
  rexpr int3d({-1,  1,  0}) end,
  rexpr int3d({-1,  1,  1}) end,
  rexpr int3d({-1,  1, -1}) end,
  rexpr int3d({-1, -1,  0}) end,
  rexpr int3d({-1, -1,  1}) end,
  rexpr int3d({-1, -1, -1}) end,
})

local tradeQueues = UTIL.generate(26, function()
  return regentlib.newsymbol(region(ispace(int1d), TradeQueue_columns))
end)
local tradeQueuePtrs = UTIL.generate(26, regentlib.newsymbol)

__demand(__cuda) -- MANUALLY PARALLELIZED
task TradeQueue_clearSource([tradeQueues])
where
  [tradeQueues:map(function(q)
     return regentlib.privilege(regentlib.writes, q, '__source')
   end)]
do
  @ESCAPE for k = 1,26 do local q = tradeQueues[k] @EMIT
    __demand(__openmp)
    for qPtr in q do
      qPtr.__source = -1
    end
  @TIME end @EPACSE
end

-- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task TradeQueue_fillSource(partColor : int3d,
                           Particles : region(ispace(int1d), Particles_columns),
                           [tradeQueues],
                           Grid_xBnum : int32, Grid_xNum : int32, NX : int32,
                           Grid_yBnum : int32, Grid_yNum : int32, NY : int32,
                           Grid_zBnum : int32, Grid_zNum : int32, NZ : int32)
where
  reads(Particles.{cell, __valid}),
  [tradeQueues:map(function(q)
     return regentlib.privilege(regentlib.writes, q, '__source')
   end)]
do
  @ESCAPE for k = 1,26 do local q = tradeQueues[k] local j = tradeQueuePtrs[k] @EMIT
    var [j] = q.bounds.lo
  @TIME end @EPACSE
  for i in Particles do
    if Particles[i].__valid then
      var elemColor = Fluid_elemColor(Particles[i].cell,
                                      Grid_xBnum, Grid_xNum, NX,
                                      Grid_yBnum, Grid_yNum, NY,
                                      Grid_zBnum, Grid_zNum, NZ)
      if elemColor ~= partColor then
        var transferred = false;
        @ESCAPE for k = 1,26 do local q = tradeQueues[k] local j = tradeQueuePtrs[k] @EMIT
          if not transferred and
             elemColor == (partColor + [colorOffsets[k]] + {NX,NY,NZ}) % {NX,NY,NZ} then
            regentlib.assert(j <= q.bounds.hi,
                             'Ran out of space in particle transfer queue')
            q[j].__source = i
            j += 1
            transferred = true
          end
        @TIME end @EPACSE
        regentlib.assert(transferred, 'Particle moved past expected stencil')
      end
    end
  end
end

__demand(__cuda) -- MANUALLY PARALLELIZED
task TradeQueue_push(Particles : region(ispace(int1d), Particles_columns),
                     [tradeQueues])
where
  reads(Particles.[Particles_subStepConserved]),
  writes(Particles.__valid),
  [tradeQueues:map(function(q)
     return Particles_subStepConserved:map(function(f)
       return regentlib.privilege(regentlib.writes, q, f)
     end)
   end):flatten()],
  [tradeQueues:map(function(q)
     return regentlib.privilege(regentlib.reads, q, '__source')
   end)]
do
  @ESCAPE for k = 1,26 do local q = tradeQueues[k] @EMIT
    __demand(__openmp)
    for j in q do
      var i = q[j].__source
      if [int](i) >= 0 then
        @ESCAPE for _,fld in ipairs(Particles_subStepConserved) do @EMIT
          q[j].[fld] = Particles[i].[fld]
        @TIME end @EPACSE
        Particles[i].__valid = false
      else
        q[j].__valid = false
      end
    end
  @TIME end @EPACSE
end

-- MANUALLY_PARALLELIZED, NO CUDA, NO OPENMP
task TradeQueue_fillTarget(Particles : region(ispace(int1d), Particles_columns),
                           [tradeQueues])
where
  reads(Particles.__valid),
  [tradeQueues:map(function(q)
     return regentlib.privilege(regentlib.reads, q, '__valid')
   end)],
  [tradeQueues:map(function(q)
     return regentlib.privilege(regentlib.writes, q, '__target')
   end)]
do
  var i = Particles.bounds.lo;
  @ESCAPE for k = 1,26 do local q = tradeQueues[k] @EMIT
    for j in q do
      if q[j].__valid then
        while i <= Particles.bounds.hi and Particles[i].__valid do
          i += 1
        end
        regentlib.assert(i <= Particles.bounds.hi,
                         'Ran out of space in particles region')
        q[j].__target = i
        i += 1
      end
    end
  @TIME end @EPACSE
end

__demand(__cuda) -- MANUALLY PARALLELIZED
task TradeQueue_pull(Particles : region(ispace(int1d), Particles_columns),
                     [tradeQueues])
where
  [tradeQueues:map(function(q)
     return Particles_subStepConserved:map(function(f)
       return regentlib.privilege(regentlib.reads, q, f)
     end)
   end):flatten()],
  [tradeQueues:map(function(q)
     return regentlib.privilege(regentlib.reads, q, '__target')
   end)],
  writes(Particles.[Particles_subStepConserved])
do
  @ESCAPE for k = 1,26 do local q = tradeQueues[k] @EMIT
    __demand(__openmp)
    for j in q do
      if q[j].__valid then
        var i = q[j].__target;
        @ESCAPE for _,fld in ipairs(Particles_subStepConserved) do @EMIT
          Particles[i].[fld] = q[j].[fld]
        @TIME end @EPACSE
      end
    end
  @TIME end @EPACSE
end

__demand(__inline)
task intersection(a : rect3d, b : SCHEMA.Volume)
  var res = rect3d{ lo = int3d{0,0,0}, hi = int3d{-1,-1,-1} }
  if  a.hi.x >= b.fromCell[0] and b.uptoCell[0] >= a.lo.x
  and a.hi.y >= b.fromCell[1] and b.uptoCell[1] >= a.lo.y
  and a.hi.z >= b.fromCell[2] and b.uptoCell[2] >= a.lo.z then
    res = rect3d{
      lo = int3d{max(a.lo.x,b.fromCell[0]), max(a.lo.y,b.fromCell[1]), max(a.lo.z,b.fromCell[2])},
      hi = int3d{min(a.hi.x,b.uptoCell[0]), min(a.hi.y,b.uptoCell[1]), min(a.hi.z,b.uptoCell[2])}}
  end
  return res
end

__demand(__inline)
task rectSize(a : rect3d)
  return (a.hi.x - a.lo.x + 1) * (a.hi.y - a.lo.y + 1) * (a.hi.z - a.lo.z + 1)
end

__demand(__inline)
task CopyQueue_partSize(fluidPartBounds : rect3d,
                        config : Config,
                        copySrc : SCHEMA.Volume)
  var totalCells = config.Grid.xNum * config.Grid.yNum * config.Grid.zNum
  var copiedCells = rectSize(intersection(fluidPartBounds, copySrc))
  return ceil(
    copiedCells
    / [double](totalCells)
    * config.Particles.maxNum
    * config.Particles.maxSkew
  )
end

-- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task CopyQueue_push(Particles : region(ispace(int1d), Particles_columns),
                    CopyQueue : region(ispace(int1d), CopyQueue_columns),
                    copySrc : SCHEMA.Volume,
                    copySrcOrigin : double[3], copyTgtOrigin : double[3],
                    Fluid0_cellWidth : double[3], Fluid1_cellWidth : double[3])
where
  reads(Particles.[Particles_primitives], Particles.cell),
  writes(CopyQueue.[Particles_primitives])
do
  var p2 = CopyQueue.bounds.lo
  for p1 in Particles do
    if p1.__valid then
      var cell = p1.cell
      if  copySrc.fromCell[0] <= cell.x and cell.x <= copySrc.uptoCell[0]
      and copySrc.fromCell[1] <= cell.y and cell.y <= copySrc.uptoCell[1]
      and copySrc.fromCell[2] <= cell.z and cell.z <= copySrc.uptoCell[2] then
        regentlib.assert(
          p2 <= CopyQueue.bounds.hi,
          'Ran out of space in cross-section particles copy queue')
        CopyQueue[p2].position =
          vv_add(copyTgtOrigin, vv_mul(Fluid1_cellWidth,
            vv_div(vv_sub(p1.position, copySrcOrigin), Fluid0_cellWidth)))
        CopyQueue[p2].velocity = p1.velocity
        CopyQueue[p2].temperature = p1.temperature
        CopyQueue[p2].diameter = p1.diameter
        CopyQueue[p2].density = p1.density
        CopyQueue[p2].__valid = true
        p2 += 1
      end
    end
  end
end

-- NOTE: It is important that Particles are placed first in the arguments list,
-- to make sure the mapper will map this task according to the sample the
-- Particles belong to (the second in a 2-section simulation). The CopyQueue
-- technically belongs to the first section.
-- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task CopyQueue_pull(partColor : int3d,
                    Particles : region(ispace(int1d), Particles_columns),
                    CopyQueue : region(ispace(int1d), CopyQueue_columns),
                    config : Config,
                    Grid_xBnum : int32, Grid_yBnum : int32, Grid_zBnum : int32)
where
  reads(CopyQueue.[Particles_primitives], Particles.__valid),
  writes(Particles.[Particles_primitives], Particles.cell)
do
  var acc : int64 = 0
  var addedVelocity = config.Particles.feeding.u.Incoming.addedVelocity
  var p1 = Particles.bounds.lo
  for p2 in CopyQueue do
    if p2.__valid then
      var cell = locate(p2.position,
                        Grid_xBnum, config.Grid.xNum, config.Grid.origin[0], config.Grid.xWidth,
                        Grid_yBnum, config.Grid.yNum, config.Grid.origin[1], config.Grid.yWidth,
                        Grid_zBnum, config.Grid.zNum, config.Grid.origin[2], config.Grid.zWidth)
      var elemColor = Fluid_elemColor(cell,
                                      Grid_xBnum, config.Grid.xNum, config.Mapping.tiles[0],
                                      Grid_yBnum, config.Grid.yNum, config.Mapping.tiles[1],
                                      Grid_zBnum, config.Grid.zNum, config.Mapping.tiles[2])
      if elemColor == partColor then
        while p1 <= Particles.bounds.hi and Particles[p1].__valid do
          p1 += 1
        end
        regentlib.assert(
          p1 <= Particles.bounds.hi,
          'Ran out of space while copying over particles from other section')
        Particles[p1].position = p2.position
        Particles[p1].velocity = vv_add(p2.velocity, addedVelocity)
        Particles[p1].temperature = p2.temperature
        Particles[p1].diameter = p2.diameter
        Particles[p1].density = p2.density
        Particles[p1].__valid = true
        acc += 1
      end
    end
  end
  return acc
end

-------------------------------------------------------------------------------
-- OTHER ROUTINES
-------------------------------------------------------------------------------

__demand(__inline)
task TrilinearInterpolateVelocity(xyz : double[3],
                                  c000 : double[3],
                                  c100 : double[3],
                                  c010 : double[3],
                                  c110 : double[3],
                                  c001 : double[3],
                                  c101 : double[3],
                                  c011 : double[3],
                                  c111 : double[3],
                                  Grid_xCellWidth : double, Grid_xRealOrigin : double,
                                  Grid_yCellWidth : double, Grid_yRealOrigin : double,
                                  Grid_zCellWidth : double, Grid_zRealOrigin : double)
  var dX = fmod((((xyz[0]-Grid_xRealOrigin)/Grid_xCellWidth)+0.5), 1.0)
  var dY = fmod((((xyz[1]-Grid_yRealOrigin)/Grid_yCellWidth)+0.5), 1.0)
  var dZ = fmod((((xyz[2]-Grid_zRealOrigin)/Grid_zCellWidth)+0.5), 1.0)
  var oneMinusdX = (1.0-dX)
  var oneMinusdY = (1.0-dY)
  var oneMinusdZ = (1.0-dZ)
  var weight00 = vv_add(vs_mul(c000, oneMinusdX), vs_mul(c100, dX))
  var weight10 = vv_add(vs_mul(c010, oneMinusdX), vs_mul(c110, dX))
  var weight01 = vv_add(vs_mul(c001, oneMinusdX), vs_mul(c101, dX))
  var weight11 = vv_add(vs_mul(c011, oneMinusdX), vs_mul(c111, dX))
  var weight0 = vv_add(vs_mul(weight00, oneMinusdY), vs_mul(weight10, dY))
  var weight1 = vv_add(vs_mul(weight01, oneMinusdY), vs_mul(weight11, dY))
  return vv_add(vs_mul(weight0, oneMinusdZ), vs_mul(weight1, dZ))
end

__demand(__inline)
task InterpolateTriVelocity(c : int3d,
                            xyz : double[3],
                            Fluid : region(ispace(int3d), Fluid_columns),
                            Grid_xCellWidth : double, Grid_xRealOrigin : double,
                            Grid_yCellWidth : double, Grid_yRealOrigin : double,
                            Grid_zCellWidth : double, Grid_zRealOrigin : double)
where
  reads(Fluid.{centerCoordinates, velocity})
do
  var i000 = Fluid[c].velocity
  var i001 = Fluid[((c+{ 0, 0, 1})%Fluid.bounds)].velocity
  var i00_ = Fluid[((c+{ 0, 0,-1})%Fluid.bounds)].velocity
  var i010 = Fluid[((c+{ 0, 1, 0})%Fluid.bounds)].velocity
  var i011 = Fluid[((c+{ 0, 1, 1})%Fluid.bounds)].velocity
  var i01_ = Fluid[((c+{ 0, 1,-1})%Fluid.bounds)].velocity
  var i0_0 = Fluid[((c+{ 0,-1, 0})%Fluid.bounds)].velocity
  var i0_1 = Fluid[((c+{ 0,-1, 1})%Fluid.bounds)].velocity
  var i0__ = Fluid[((c+{ 0,-1,-1})%Fluid.bounds)].velocity

  var i100 = Fluid[((c+{ 1, 0, 0})%Fluid.bounds)].velocity
  var i101 = Fluid[((c+{ 1, 0, 1})%Fluid.bounds)].velocity
  var i10_ = Fluid[((c+{ 1, 0,-1})%Fluid.bounds)].velocity
  var i110 = Fluid[((c+{ 1, 1, 0})%Fluid.bounds)].velocity
  var i111 = Fluid[((c+{ 1, 1, 1})%Fluid.bounds)].velocity
  var i11_ = Fluid[((c+{ 1, 1,-1})%Fluid.bounds)].velocity
  var i1_0 = Fluid[((c+{ 1,-1, 0})%Fluid.bounds)].velocity
  var i1_1 = Fluid[((c+{ 1,-1, 1})%Fluid.bounds)].velocity
  var i1__ = Fluid[((c+{ 1,-1,-1})%Fluid.bounds)].velocity

  var i_00 = Fluid[((c+{-1, 0, 0})%Fluid.bounds)].velocity
  var i_01 = Fluid[((c+{-1, 0, 1})%Fluid.bounds)].velocity
  var i_0_ = Fluid[((c+{-1, 0,-1})%Fluid.bounds)].velocity
  var i_10 = Fluid[((c+{-1, 1, 0})%Fluid.bounds)].velocity
  var i_11 = Fluid[((c+{-1, 1, 1})%Fluid.bounds)].velocity
  var i_1_ = Fluid[((c+{-1, 1,-1})%Fluid.bounds)].velocity
  var i__0 = Fluid[((c+{-1,-1, 0})%Fluid.bounds)].velocity
  var i__1 = Fluid[((c+{-1,-1, 1})%Fluid.bounds)].velocity
  var i___ = Fluid[((c+{-1,-1,-1})%Fluid.bounds)].velocity

  var v000 = array(0.0, 0.0, 0.0)
  var v001 = array(0.0, 0.0, 0.0)
  var v010 = array(0.0, 0.0, 0.0)
  var v011 = array(0.0, 0.0, 0.0)
  var v100 = array(0.0, 0.0, 0.0)
  var v101 = array(0.0, 0.0, 0.0)
  var v110 = array(0.0, 0.0, 0.0)
  var v111 = array(0.0, 0.0, 0.0)

  if (xyz[0]>Fluid[c].centerCoordinates[0]) then
    if (xyz[1]>Fluid[c].centerCoordinates[1]) then
      if (xyz[2]>Fluid[c].centerCoordinates[2]) then
        v000 = i000
        v001 = i001
        v010 = i010
        v011 = i011
        v100 = i100
        v101 = i101
        v110 = i110
        v111 = i111
      else
        v000 = i00_
        v001 = i000
        v010 = i01_
        v011 = i010
        v100 = i10_
        v101 = i100
        v110 = i11_
        v111 = i110
      end
    else
      if (xyz[2]>Fluid[c].centerCoordinates[2]) then
        v000 = i0_0
        v001 = i0_1
        v010 = i000
        v011 = i001
        v100 = i1_0
        v101 = i1_1
        v110 = i100
        v111 = i101
      else
        v000 = i0__
        v001 = i0_0
        v010 = i00_
        v011 = i000
        v100 = i1__
        v101 = i1_0
        v110 = i10_
        v111 = i100
      end
    end
  else
    if (xyz[1]>Fluid[c].centerCoordinates[1]) then
      if (xyz[2]>Fluid[c].centerCoordinates[2]) then
        v000 = i_00
        v001 = i_01
        v010 = i_10
        v011 = i_11
        v100 = i000
        v101 = i001
        v110 = i010
        v111 = i011
      else
        v000 = i_0_
        v001 = i_00
        v010 = i_1_
        v011 = i_10
        v100 = i00_
        v101 = i000
        v110 = i01_
        v111 = i010
      end
    else
      if (xyz[2]>Fluid[c].centerCoordinates[2]) then
        v000 = i__0
        v001 = i__1
        v010 = i_00
        v011 = i_01
        v100 = i0_0
        v101 = i0_1
        v110 = i000
        v111 = i001
      else
        v000 = i___
        v001 = i__0
        v010 = i_0_
        v011 = i_00
        v100 = i0__
        v101 = i0_0
        v110 = i00_
        v111 = i000
      end
    end
  end

  return TrilinearInterpolateVelocity(xyz, v000, v100, v010, v110, v001, v101, v011, v111, Grid_xCellWidth, Grid_xRealOrigin, Grid_yCellWidth, Grid_yRealOrigin, Grid_zCellWidth, Grid_zRealOrigin)
end

__demand(__inline)
task TrilinearInterpolateTemp(xyz : double[3],
                              c000 : double,
                              c100 : double,
                              c010 : double,
                              c110 : double,
                              c001 : double,
                              c101 : double,
                              c011 : double,
                              c111 : double,
                              Grid_xCellWidth : double, Grid_xRealOrigin : double,
                              Grid_yCellWidth : double, Grid_yRealOrigin : double,
                              Grid_zCellWidth : double, Grid_zRealOrigin : double)
  var dX = fmod((((xyz[0]-Grid_xRealOrigin)/Grid_xCellWidth)+0.5), 1.0)
  var dY = fmod((((xyz[1]-Grid_yRealOrigin)/Grid_yCellWidth)+0.5), 1.0)
  var dZ = fmod((((xyz[2]-Grid_zRealOrigin)/Grid_zCellWidth)+0.5), 1.0)
  var oneMinusdX = (1.0-dX)
  var oneMinusdY = (1.0-dY)
  var oneMinusdZ = (1.0-dZ)
  var weight00 = ((c000*oneMinusdX)+(c100*dX))
  var weight10 = ((c010*oneMinusdX)+(c110*dX))
  var weight01 = ((c001*oneMinusdX)+(c101*dX))
  var weight11 = ((c011*oneMinusdX)+(c111*dX))
  var weight0 = ((weight00*oneMinusdY)+(weight10*dY))
  var weight1 = ((weight01*oneMinusdY)+(weight11*dY))
  return ((weight0*oneMinusdZ)+(weight1*dZ))
end

__demand(__inline)
task InterpolateTriTemp(c : int3d,
                        xyz : double[3],
                        Fluid : region(ispace(int3d), Fluid_columns),
                        Grid_xCellWidth : double, Grid_xRealOrigin : double,
                        Grid_yCellWidth : double, Grid_yRealOrigin : double,
                        Grid_zCellWidth : double, Grid_zRealOrigin : double)
where
  reads(Fluid.{centerCoordinates, temperature})
do
  var i000 = Fluid[c].temperature
  var i001 = Fluid[((c+{ 0, 0, 1})%Fluid.bounds)].temperature
  var i00_ = Fluid[((c+{ 0, 0,-1})%Fluid.bounds)].temperature
  var i010 = Fluid[((c+{ 0, 1, 0})%Fluid.bounds)].temperature
  var i011 = Fluid[((c+{ 0, 1, 1})%Fluid.bounds)].temperature
  var i01_ = Fluid[((c+{ 0, 1,-1})%Fluid.bounds)].temperature
  var i0_0 = Fluid[((c+{ 0,-1, 0})%Fluid.bounds)].temperature
  var i0_1 = Fluid[((c+{ 0,-1, 1})%Fluid.bounds)].temperature
  var i0__ = Fluid[((c+{ 0,-1,-1})%Fluid.bounds)].temperature

  var i100 = Fluid[((c+{ 1, 0, 0})%Fluid.bounds)].temperature
  var i101 = Fluid[((c+{ 1, 0, 1})%Fluid.bounds)].temperature
  var i10_ = Fluid[((c+{ 1, 0,-1})%Fluid.bounds)].temperature
  var i110 = Fluid[((c+{ 1, 1, 0})%Fluid.bounds)].temperature
  var i111 = Fluid[((c+{ 1, 1, 1})%Fluid.bounds)].temperature
  var i11_ = Fluid[((c+{ 1, 1,-1})%Fluid.bounds)].temperature
  var i1_0 = Fluid[((c+{ 1,-1, 0})%Fluid.bounds)].temperature
  var i1_1 = Fluid[((c+{ 1,-1, 1})%Fluid.bounds)].temperature
  var i1__ = Fluid[((c+{ 1,-1,-1})%Fluid.bounds)].temperature

  var i_00 = Fluid[((c+{-1, 0, 0})%Fluid.bounds)].temperature
  var i_01 = Fluid[((c+{-1, 0, 1})%Fluid.bounds)].temperature
  var i_0_ = Fluid[((c+{-1, 0,-1})%Fluid.bounds)].temperature
  var i_10 = Fluid[((c+{-1, 1, 0})%Fluid.bounds)].temperature
  var i_11 = Fluid[((c+{-1, 1, 1})%Fluid.bounds)].temperature
  var i_1_ = Fluid[((c+{-1, 1,-1})%Fluid.bounds)].temperature
  var i__0 = Fluid[((c+{-1,-1, 0})%Fluid.bounds)].temperature
  var i__1 = Fluid[((c+{-1,-1, 1})%Fluid.bounds)].temperature
  var i___ = Fluid[((c+{-1,-1,-1})%Fluid.bounds)].temperature

  var v000 = 0.0
  var v001 = 0.0
  var v010 = 0.0
  var v011 = 0.0
  var v100 = 0.0
  var v101 = 0.0
  var v110 = 0.0
  var v111 = 0.0

  if (xyz[0]>Fluid[c].centerCoordinates[0]) then
    if (xyz[1]>Fluid[c].centerCoordinates[1]) then
      if (xyz[2]>Fluid[c].centerCoordinates[2]) then
        v000 = i000
        v001 = i001
        v010 = i010
        v011 = i011
        v100 = i100
        v101 = i101
        v110 = i110
        v111 = i111
      else
        v000 = i00_
        v001 = i000
        v010 = i01_
        v011 = i010
        v100 = i10_
        v101 = i100
        v110 = i11_
        v111 = i110
      end
    else
      if (xyz[2]>Fluid[c].centerCoordinates[2]) then
        v000 = i0_0
        v001 = i0_1
        v010 = i000
        v011 = i001
        v100 = i1_0
        v101 = i1_1
        v110 = i100
        v111 = i101
      else
        v000 = i0__
        v001 = i0_0
        v010 = i00_
        v011 = i000
        v100 = i1__
        v101 = i1_0
        v110 = i10_
        v111 = i100
      end
    end
  else
    if (xyz[1]>Fluid[c].centerCoordinates[1]) then
      if (xyz[2]>Fluid[c].centerCoordinates[2]) then
        v000 = i_00
        v001 = i_01
        v010 = i_10
        v011 = i_11
        v100 = i000
        v101 = i001
        v110 = i010
        v111 = i011
      else
        v000 = i_0_
        v001 = i_00
        v010 = i_1_
        v011 = i_10
        v100 = i00_
        v101 = i000
        v110 = i01_
        v111 = i010
      end
    else
      if (xyz[2]>Fluid[c].centerCoordinates[2]) then
        v000 = i__0
        v001 = i__1
        v010 = i_00
        v011 = i_01
        v100 = i0_0
        v101 = i0_1
        v110 = i000
        v111 = i001
      else
        v000 = i___
        v001 = i__0
        v010 = i_0_
        v011 = i_00
        v100 = i0__
        v101 = i0_0
        v110 = i00_
        v111 = i000
      end
    end
  end

  return TrilinearInterpolateTemp(xyz, v000, v100, v010, v110, v001, v101, v011, v111, Grid_xCellWidth, Grid_xRealOrigin, Grid_yCellWidth, Grid_yRealOrigin, Grid_zCellWidth, Grid_zRealOrigin)
end

__demand(__parallel, __cuda)
task Particles_AddFlowCoupling(Particles : region(ispace(int1d), Particles_columns),
                               Fluid : region(ispace(int3d), Fluid_columns),
                               Flow_constantVisc : double,
                               Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                               Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                               Flow_viscosityModel : SCHEMA.ViscosityModel,
                               Grid_xCellWidth : double, Grid_xRealOrigin : double,
                               Grid_yCellWidth : double, Grid_yRealOrigin : double,
                               Grid_zCellWidth : double, Grid_zRealOrigin : double,
                               Particles_convectiveCoeff : double,
                               Particles_heatCapacity : double)
where
  reads(Fluid.{centerCoordinates, velocity, temperature}),
  reads(Particles.{cell, position, velocity, diameter, density, temperature, __valid}),
  reads writes(Particles.{position_t, velocity_t, temperature_t}),
  writes(Particles.{deltaTemperatureTerm, deltaVelocityOverRelaxationTime})
do
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      var flowVelocity = InterpolateTriVelocity(Particles[p].cell, Particles[p].position, Fluid, Grid_xCellWidth, Grid_xRealOrigin, Grid_yCellWidth, Grid_yRealOrigin, Grid_zCellWidth, Grid_zRealOrigin)
      var flowTemperature = InterpolateTriTemp(Particles[p].cell, Particles[p].position, Fluid, Grid_xCellWidth, Grid_xRealOrigin, Grid_yCellWidth, Grid_yRealOrigin, Grid_zCellWidth, Grid_zRealOrigin)
      var flowDynamicViscosity = GetDynamicViscosity(flowTemperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)
      Particles[p].position_t = vv_add(Particles[p].position_t, Particles[p].velocity)
      var particleReynoldsNumber = 0.0
      var relaxationTime = (((Particles[p].density*pow(Particles[p].diameter, 2.0))/(18.0*flowDynamicViscosity))/(1.0+(double(0.15)*pow(particleReynoldsNumber, double(0.687)))))
      var tmp2 = vs_div(vv_sub(flowVelocity, Particles[p].velocity), relaxationTime)
      Particles[p].deltaVelocityOverRelaxationTime = tmp2
      Particles[p].velocity_t = vv_add(Particles[p].velocity_t, tmp2)
      var tmp3 = (((PI*pow(Particles[p].diameter, 2.0))*Particles_convectiveCoeff)*(flowTemperature-Particles[p].temperature))
      Particles[p].deltaTemperatureTerm = tmp3
      Particles[p].temperature_t += (tmp3/((((PI*pow(Particles[p].diameter, 3.0))/6.0)*Particles[p].density)*Particles_heatCapacity))
    end
  end
end

__demand(__parallel, __cuda)
task Particles_AddBodyForces(Particles : region(ispace(int1d), Particles_columns),
                             Particles_bodyForce : double[3])
where
  reads(Particles.__valid),
  reads writes(Particles.velocity_t)
do
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      var tmp = Particles_bodyForce
      var v = Particles[p].velocity_t
      v[0] += tmp[0]
      v[1] += tmp[1]
      v[2] += tmp[2]
      Particles[p].velocity_t = v
    end
  end
end

__demand(__cuda) -- MANUALLY PARALLELIZED
task Radiation_AccumulateParticleValues(Particles : region(ispace(int1d), Particles_columns),
                                        Fluid : region(ispace(int3d), Fluid_columns),
                                        Radiation : region(ispace(int3d), Radiation_columns))
where
  reads(Fluid.to_Radiation),
  reads(Particles.{cell, diameter, temperature, __valid}),
  reads writes(Radiation.{acc_d2, acc_d2t4})
do
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      Radiation[Fluid[Particles[p].cell].to_Radiation].acc_d2 += pow(Particles[p].diameter, 2.0)
      Radiation[Fluid[Particles[p].cell].to_Radiation].acc_d2t4 += (pow(Particles[p].diameter, 2.0)*pow(Particles[p].temperature, 4.0))
    end
  end
end

__demand(__parallel, __cuda)
task Radiation_UpdateFieldValues(Radiation : region(ispace(int3d), Radiation_columns),
                                 Radiation_cellVolume : double,
                                 Radiation_qa : double,
                                 Radiation_qs : double)
where
  reads(Radiation.{acc_d2, acc_d2t4}),
  writes(Radiation.{Ib, sigma})
do
  __demand(__openmp)
  for c in Radiation do
    c.sigma = c.acc_d2*PI*(Radiation_qa+Radiation_qs)/(4.0*Radiation_cellVolume)
    if c.acc_d2 == 0.0 then
      c.Ib = 0.0
    else
      c.Ib = (SB*c.acc_d2t4)/(PI*c.acc_d2)
    end
  end
end

__demand(__cuda) -- MANUALLY PARALLELIZED
task Particles_AbsorbRadiation(Particles : region(ispace(int1d), Particles_columns),
                               Fluid : region(ispace(int3d), Fluid_columns),
                               Radiation : region(ispace(int3d), Radiation_columns),
                               Particles_heatCapacity : double,
                               Radiation_qa : double)
where
  reads(Fluid.to_Radiation),
  reads(Radiation.G),
  reads(Particles.{cell, density, diameter, temperature, __valid}),
  reads writes(Particles.temperature_t)
do
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      var t4 = pow(Particles[p].temperature, 4.0)
      var alpha = ((((PI*Radiation_qa)*pow(Particles[p].diameter, 2.0))*(Radiation[Fluid[Particles[p].cell].to_Radiation].G-((4.0*SB)*t4)))/4.0)
      Particles[p].temperature_t += (alpha/((((PI*pow(Particles[p].diameter, 3.0))/6.0)*Particles[p].density)*Particles_heatCapacity))
    end
  end
end

__demand(__parallel, __cuda)
task Flow_AddParticlesCoupling(Particles : region(ispace(int1d), Particles_columns),
                               Fluid : region(ispace(int3d), Fluid_columns),
                               Grid_cellVolume : double)
where
  reads(Particles.{cell, diameter, density, deltaTemperatureTerm, deltaVelocityOverRelaxationTime, __valid}),
  reads writes(Fluid.{rhoVelocity_t, rhoEnergy_t})
do
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      var tmp = vs_div(vs_mul(Particles[p].deltaVelocityOverRelaxationTime, (-(((PI*pow(Particles[p].diameter, 3.0))/6.0)*Particles[p].density))), Grid_cellVolume)
      var v = Fluid[Particles[p].cell].rhoVelocity_t
      v[0] += tmp[0]
      v[1] += tmp[1]
      v[2] += tmp[2]
      Fluid[Particles[p].cell].rhoVelocity_t = v
      Fluid[Particles[p].cell].rhoEnergy_t += ((-Particles[p].deltaTemperatureTerm)/Grid_cellVolume)
    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdateVars(Fluid : region(ispace(int3d), Fluid_columns),
                     Integrator_deltaTime : double,
                     Integrator_stage : int32)
where
  reads(Fluid.{rho_old, rhoEnergy_old, rhoVelocity_old}),
  reads(Fluid.{rho_t, rhoEnergy_t, rhoVelocity_t}),
  writes(Fluid.{rho, rhoEnergy, rhoVelocity}),
  reads writes(Fluid.{rho_new, rhoEnergy_new, rhoVelocity_new})
do
  __demand(__openmp)
  for c in Fluid do
    var deltaTime = Integrator_deltaTime
    if (Integrator_stage==1) then
      Fluid[c].rho_new += (((1.0/6.0)*deltaTime)*Fluid[c].rho_t)
      Fluid[c].rho = (Fluid[c].rho_old+((0.5*deltaTime)*Fluid[c].rho_t))

      var tmp = vs_mul(Fluid[c].rhoVelocity_t, ((1.0/6.0)*deltaTime))
      var v = Fluid[c].rhoVelocity_new
      v[0] += tmp[0]
      v[1] += tmp[1]
      v[2] += tmp[2]
      Fluid[c].rhoVelocity_new = v
      Fluid[c].rhoVelocity = vv_add(Fluid[c].rhoVelocity_old, vs_mul(Fluid[c].rhoVelocity_t, (0.5*deltaTime)))

      Fluid[c].rhoEnergy_new += (((1.0/6.0)*deltaTime)*Fluid[c].rhoEnergy_t)
      Fluid[c].rhoEnergy = (Fluid[c].rhoEnergy_old+((0.5*deltaTime)*Fluid[c].rhoEnergy_t))
    else
      if (Integrator_stage==2) then
        Fluid[c].rho_new += (((1.0/3.0)*deltaTime)*Fluid[c].rho_t)
        Fluid[c].rho = (Fluid[c].rho_old+((0.5*deltaTime)*Fluid[c].rho_t))
        var tmp = vs_mul(Fluid[c].rhoVelocity_t, ((1.0/3.0)*deltaTime))
        var v = Fluid[c].rhoVelocity_new
        v[0] += tmp[0]
        v[1] += tmp[1]
        v[2] += tmp[2]
        Fluid[c].rhoVelocity_new = v
        Fluid[c].rhoVelocity = vv_add(Fluid[c].rhoVelocity_old, vs_mul(Fluid[c].rhoVelocity_t, (0.5*deltaTime)))
        Fluid[c].rhoEnergy_new += (((1.0/3.0)*deltaTime)*Fluid[c].rhoEnergy_t)
        Fluid[c].rhoEnergy = (Fluid[c].rhoEnergy_old+((0.5*deltaTime)*Fluid[c].rhoEnergy_t))
      else
        if (Integrator_stage==3) then
          Fluid[c].rho_new += (((1.0/3.0)*deltaTime)*Fluid[c].rho_t)
          Fluid[c].rho = (Fluid[c].rho_old+((1.0*deltaTime)*Fluid[c].rho_t))
          var tmp = vs_mul(Fluid[c].rhoVelocity_t, ((1.0/3.0)*deltaTime))
          var v = Fluid[c].rhoVelocity_new
          v[0] += tmp[0]
          v[1] += tmp[1]
          v[2] += tmp[2]
          Fluid[c].rhoVelocity_new = v
          Fluid[c].rhoVelocity = vv_add(Fluid[c].rhoVelocity_old, vs_mul(Fluid[c].rhoVelocity_t, (1.0*deltaTime)))
          Fluid[c].rhoEnergy_new += (((1.0/3.0)*deltaTime)*Fluid[c].rhoEnergy_t)
          Fluid[c].rhoEnergy = (Fluid[c].rhoEnergy_old+((1.0*deltaTime)*Fluid[c].rhoEnergy_t))
        else
          Fluid[c].rho = (Fluid[c].rho_new+(((1.0/6.0)*deltaTime)*Fluid[c].rho_t))
          Fluid[c].rhoVelocity = vv_add(Fluid[c].rhoVelocity_new, vs_mul(Fluid[c].rhoVelocity_t, ((1.0/6.0)*deltaTime)))
          Fluid[c].rhoEnergy = (Fluid[c].rhoEnergy_new+(((1.0/6.0)*deltaTime)*Fluid[c].rhoEnergy_t))
        end
      end
    end
  end
end

__demand(__parallel, __cuda)
task Particles_UpdateVars(Particles : region(ispace(int1d), Particles_columns),
                          Integrator_deltaTime : double,
                          Integrator_stage : int32)
where
  reads(Particles.{position_old, velocity_old, temperature_old}),
  reads(Particles.{position_t, velocity_t, temperature_t}),
  reads(Particles.__valid),
  writes(Particles.{position, temperature, velocity}),
  reads writes(Particles.{position_new, temperature_new, velocity_new})
do
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      var deltaTime = Integrator_deltaTime
      if Integrator_stage == 1 then
        var tmp = vs_mul(Particles[p].position_t, ((1.0/6.0)*deltaTime))
        var v = Particles[p].position_new
        v[0] += tmp[0]
        v[1] += tmp[1]
        v[2] += tmp[2]
        Particles[p].position_new = v
        Particles[p].position = vv_add(Particles[p].position_old, vs_mul(Particles[p].position_t, (0.5*deltaTime)))
        var tmp__11020 = vs_mul(Particles[p].velocity_t, ((1.0/6.0)*deltaTime))
        var v__11021 = Particles[p].velocity_new
        v__11021[0] += tmp__11020[0]
        v__11021[1] += tmp__11020[1]
        v__11021[2] += tmp__11020[2]
        Particles[p].velocity_new = v__11021
        Particles[p].velocity = vv_add(Particles[p].velocity_old, vs_mul(Particles[p].velocity_t, (0.5*deltaTime)))
        Particles[p].temperature_new += (((1.0/6.0)*deltaTime)*Particles[p].temperature_t)
        Particles[p].temperature = (Particles[p].temperature_old+((0.5*deltaTime)*Particles[p].temperature_t))
      elseif Integrator_stage == 2 then
        var tmp = vs_mul(Particles[p].position_t, ((1.0/3.0)*deltaTime))
        var v = Particles[p].position_new
        v[0] += tmp[0]
        v[1] += tmp[1]
        v[2] += tmp[2]
        Particles[p].position_new = v
        Particles[p].position = vv_add(Particles[p].position_old, vs_mul(Particles[p].position_t, (0.5*deltaTime)))
        var tmp__11024 = vs_mul(Particles[p].velocity_t, ((1.0/3.0)*deltaTime))
        var v__11025 = Particles[p].velocity_new
        v__11025[0] += tmp__11024[0]
        v__11025[1] += tmp__11024[1]
        v__11025[2] += tmp__11024[2]
        Particles[p].velocity_new = v__11025
        Particles[p].velocity = vv_add(Particles[p].velocity_old, vs_mul(Particles[p].velocity_t, (0.5*deltaTime)))
        Particles[p].temperature_new += (((1.0/3.0)*deltaTime)*Particles[p].temperature_t)
        Particles[p].temperature = (Particles[p].temperature_old+((0.5*deltaTime)*Particles[p].temperature_t))
      elseif Integrator_stage == 3 then
        var tmp = vs_mul(Particles[p].position_t, ((1.0/3.0)*deltaTime))
        var v = Particles[p].position_new
        v[0] += tmp[0]
        v[1] += tmp[1]
        v[2] += tmp[2]
        Particles[p].position_new = v
        Particles[p].position = vv_add(Particles[p].position_old, vs_mul(Particles[p].position_t, (1.0*deltaTime)))
        var tmp__11028 = vs_mul(Particles[p].velocity_t, ((1.0/3.0)*deltaTime))
        var v__11029 = Particles[p].velocity_new
        v__11029[0] += tmp__11028[0]
        v__11029[1] += tmp__11028[1]
        v__11029[2] += tmp__11028[2]
        Particles[p].velocity_new = v__11029
        Particles[p].velocity = vv_add(Particles[p].velocity_old, vs_mul(Particles[p].velocity_t, (1.0*deltaTime)))
        Particles[p].temperature_new += (((1.0/3.0)*deltaTime)*Particles[p].temperature_t)
        Particles[p].temperature = (Particles[p].temperature_old+((1.0*deltaTime)*Particles[p].temperature_t))
      else -- Integrator_stage == 4
        Particles[p].position = vv_add(Particles[p].position_new, vs_mul(Particles[p].position_t, ((1.0/6.0)*deltaTime)))
        Particles[p].velocity = vv_add(Particles[p].velocity_new, vs_mul(Particles[p].velocity_t, ((1.0/6.0)*deltaTime)))
        Particles[p].temperature = (Particles[p].temperature_new+(((1.0/6.0)*deltaTime)*Particles[p].temperature_t))
      end
    end
  end
end

-- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task Particles_HandleCollisions(Particles : region(ispace(int1d), Particles_columns),
                                Integrator_deltaTime : double, Particles_restitutionCoeff : double )
-- This is an adaption of collisionPrt routine of the Soleil-MPI version
-- TODO: search box implementation
where
  reads(Particles.{position_old, diameter, density, __valid}),
  reads writes(Particles.{position, velocity})
do
  for p1 in Particles do
    if p1.__valid then
      for p2 in Particles do
        if p2.__valid and p1 < p2 then

          -- Relative position of particles
          var x = p2.position[0] - p1.position[0]
          var y = p2.position[1] - p1.position[1]
          var z = p2.position[2] - p1.position[2]

          -- Old relative position of particles
          var xold = p2.position_old[0] - p1.position_old[0]
          var yold = p2.position_old[1] - p1.position_old[1]
          var zold = p2.position_old[2] - p1.position_old[2]


          -- Relative velocity
          var ux = (x-xold)/Integrator_deltaTime
          var uy = (y-yold)/Integrator_deltaTime
          var uz = (z-zold)/Integrator_deltaTime

          -- Relevant scalar products
          var x_scal_u = xold*ux + yold*uy + zold*uz
          var x_scal_x = xold*xold + yold*yold + zold*zold
          var u_scal_u = ux*ux + uy*uy + uz*uz

          -- Critical distance
          var dcrit = 0.5*( p1.diameter + p2.diameter )

          -- Checking if particles are getting away from each other
          if x_scal_u<0.0 then

            -- Checking if particles are in a collision path
            var det = x_scal_u*x_scal_u - u_scal_u*(x_scal_x - dcrit*dcrit)
            if det>0.0 then

              -- Checking if collision occurs in this time step
              var timecol = ( -x_scal_u - sqrt(det) ) / u_scal_u
              if (timecol>0.0 and timecol<Integrator_deltaTime) then

                -- We do have a collision


                -- Mass ratio of particles
                var mr = (p2.density * p2.diameter * p2.diameter * p2.diameter)
                mr = mr/ (p1.density * p1.diameter * p1.diameter * p1.diameter)

                -- Change of velocity and particle location after impact
                -- Note: for now particle restitution coeff is the same for all particles ?
                var du = ( 1.0 + min( Particles_restitutionCoeff, Particles_restitutionCoeff ) ) / (1.0 + mr)*x_scal_u/x_scal_x
                var dx = du * ( Integrator_deltaTime - timecol )

                -- Update velocities

                p1.velocity[0] = p1.velocity[0] + du*xold*mr
                p1.velocity[1] = p1.velocity[1] + du*yold*mr
                p1.velocity[2] = p1.velocity[2] + du*zold*mr

                p2.velocity[0] = p2.velocity[0] - du*xold
                p2.velocity[1] = p2.velocity[1] - du*yold
                p2.velocity[2] = p2.velocity[2] - du*zold

                -- Update positions
                p1.position[0] = p1.position[0] + dx*xold*mr
                p1.position[1] = p1.position[1] + dx*yold*mr
                p1.position[2] = p1.position[2] + dx*zold*mr

                p2.position[0] = p2.position[0] - dx*xold
                p2.position[1] = p2.position[1] - dx*yold
                p2.position[2] = p2.position[2] - dx*zold

              end

            end
          end

        end
      end
    end
  end
end

__demand(__parallel, __cuda)
task Particles_UpdateAuxiliaryStep1(Particles : region(ispace(int1d), Particles_columns),
                                    BC_xBCParticles : SCHEMA.ParticlesBC,
                                    BC_yBCParticles : SCHEMA.ParticlesBC,
                                    BC_zBCParticles : SCHEMA.ParticlesBC,
                                    Grid_xOrigin : double, Grid_xWidth : double,
                                    Grid_yOrigin : double, Grid_yWidth : double,
                                    Grid_zOrigin : double, Grid_zWidth : double,
                                    Particles_restitutionCoeff : double)
where
  reads(Particles.{position, velocity, velocity_t, __valid}),
  reads writes(Particles.{position_ghost, velocity_ghost, velocity_t_ghost})
do
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      Particles[p].position_ghost[0] = Particles[p].position[0]
      Particles[p].position_ghost[1] = Particles[p].position[1]
      Particles[p].position_ghost[2] = Particles[p].position[2]
      Particles[p].velocity_ghost[0] = Particles[p].velocity[0]
      Particles[p].velocity_ghost[1] = Particles[p].velocity[1]
      Particles[p].velocity_ghost[2] = Particles[p].velocity[2]
      Particles[p].velocity_t_ghost[0] = Particles[p].velocity_t[0]
      Particles[p].velocity_t_ghost[1] = Particles[p].velocity_t[1]
      Particles[p].velocity_t_ghost[2] = Particles[p].velocity_t[2]
      if (Particles[p].position[0]<Grid_xOrigin) then
        if BC_xBCParticles == SCHEMA.ParticlesBC_Periodic then
          Particles[p].position_ghost[0] = (Particles[p].position[0]+Grid_xWidth)
        elseif BC_xBCParticles == SCHEMA.ParticlesBC_Bounce then
          Particles[p].position_ghost[0] = Grid_xOrigin
          var impulse = ((-(1.0+Particles_restitutionCoeff))*Particles[p].velocity[0])
          if (impulse<=0.0) then
            Particles[p].velocity_ghost[0] += impulse
          end
          var contact_force = (double(-1)*Particles[p].velocity_t[0])
          if (contact_force>0.0) then
            Particles[p].velocity_t_ghost[0] += contact_force
          end
        else -- BC_xBCParticles == SCHEMA.ParticlesBC_Disappear
          -- Do nothing, let out-of-bounds particles get deleted
        end
      end
      if (Particles[p].position[0]>(Grid_xOrigin+Grid_xWidth)) then
        if BC_xBCParticles == SCHEMA.ParticlesBC_Periodic then
          Particles[p].position_ghost[0] = (Particles[p].position[0]-Grid_xWidth)
        elseif BC_xBCParticles == SCHEMA.ParticlesBC_Bounce then
          Particles[p].position_ghost[0] = (Grid_xOrigin+Grid_xWidth)
          var impulse = ((-(1.0+Particles_restitutionCoeff))*Particles[p].velocity[0])
          if (impulse>=0.0) then
            Particles[p].velocity_ghost[0] += impulse
          end
          var contact_force = (double(-1)*Particles[p].velocity_t[0])
          if (contact_force<0.0) then
            Particles[p].velocity_t_ghost[0] += contact_force
          end
        else -- BC_xBCParticles == SCHEMA.ParticlesBC_Disappear
          -- Do nothing, let out-of-bounds particles get deleted
        end
      end
      if (Particles[p].position[1]<Grid_yOrigin) then
        if BC_yBCParticles == SCHEMA.ParticlesBC_Periodic then
          Particles[p].position_ghost[1] = (Particles[p].position[1]+Grid_yWidth)
        elseif BC_yBCParticles == SCHEMA.ParticlesBC_Bounce then
          Particles[p].position_ghost[1] = Grid_yOrigin
          var impulse = ((-(1.0+Particles_restitutionCoeff))*Particles[p].velocity[1])
          if (impulse<=0.0) then
            Particles[p].velocity_ghost[1] += impulse
          end
          var contact_force = (double(-1)*Particles[p].velocity_t[1])
          if (contact_force>0.0) then
            Particles[p].velocity_t_ghost[1] += contact_force
          end
        else -- BC_yBCParticles == SCHEMA.ParticlesBC_Disappear
          -- Do nothing, let out-of-bounds particles get deleted
        end
      end
      if (Particles[p].position[1]>(Grid_yOrigin+Grid_yWidth)) then
        if BC_yBCParticles == SCHEMA.ParticlesBC_Periodic then
          Particles[p].position_ghost[1] = (Particles[p].position[1]-Grid_yWidth)
        elseif BC_yBCParticles == SCHEMA.ParticlesBC_Bounce then
          Particles[p].position_ghost[1] = (Grid_yOrigin+Grid_yWidth)
          var impulse = ((-(1.0+Particles_restitutionCoeff))*Particles[p].velocity[1])
          if (impulse>=0.0) then
            Particles[p].velocity_ghost[1] += impulse
          end
          var contact_force = (double(-1)*Particles[p].velocity_t[1])
          if (contact_force<0.0) then
            Particles[p].velocity_t_ghost[1] += contact_force
          end
        else -- BC_yBCParticles == SCHEMA.ParticlesBC_Disappear
          -- Do nothing, let out-of-bounds particles get deleted
        end
      end
      if (Particles[p].position[2]<Grid_zOrigin) then
        if BC_zBCParticles == SCHEMA.ParticlesBC_Periodic then
          Particles[p].position_ghost[2] = (Particles[p].position[2]+Grid_zWidth)
        elseif BC_zBCParticles == SCHEMA.ParticlesBC_Bounce then
          Particles[p].position_ghost[2] = Grid_zOrigin
          var impulse = ((-(1.0+Particles_restitutionCoeff))*Particles[p].velocity[2])
          if (impulse<=0.0) then
            Particles[p].velocity_ghost[2] += impulse
          end
          var contact_force = (double(-1)*Particles[p].velocity_t[2])
          if (contact_force>0.0) then
            Particles[p].velocity_t_ghost[2] += contact_force
          end
        else -- BC_zBCParticles == SCHEMA.ParticlesBC_Disappear
          -- Do nothing, let out-of-bounds particles get deleted
        end
      end
      if (Particles[p].position[2]>(Grid_zOrigin+Grid_zWidth)) then
        if BC_zBCParticles == SCHEMA.ParticlesBC_Periodic then
          Particles[p].position_ghost[2] = (Particles[p].position[2]-Grid_zWidth)
        elseif BC_zBCParticles == SCHEMA.ParticlesBC_Bounce then
          Particles[p].position_ghost[2] = (Grid_zOrigin+Grid_zWidth)
          var impulse = ((-(1.0+Particles_restitutionCoeff))*Particles[p].velocity[2])
          if (impulse>=0.0) then
            Particles[p].velocity_ghost[2] += impulse
          end
          var contact_force = (double(-1)*Particles[p].velocity_t[2])
          if (contact_force<0.0) then
            Particles[p].velocity_t_ghost[2] += contact_force
          end
        else -- BC_zBCParticles == SCHEMA.ParticlesBC_Disappear
          -- Do nothing, let out-of-bounds particles get deleted
        end
      end
    end
  end
end

__demand(__parallel, __cuda)
task Particles_UpdateAuxiliaryStep2(Particles : region(ispace(int1d), Particles_columns))
where
  reads(Particles.{position_ghost, velocity_ghost, velocity_t_ghost, __valid}),
  reads writes(Particles.{position, velocity, velocity_t})
do
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      Particles[p].position = Particles[p].position_ghost
      Particles[p].velocity = Particles[p].velocity_ghost
      Particles[p].velocity_t = Particles[p].velocity_t_ghost
    end
  end
end

__demand(__cuda) -- MANUALLY PARALLELIZED
task Particles_DeleteEscapingParticles(Particles : region(ispace(int1d), Particles_columns),
                                       Grid_xOrigin : double, Grid_xWidth : double,
                                       Grid_yOrigin : double, Grid_yWidth : double,
                                       Grid_zOrigin : double, Grid_zWidth : double)
where
  reads(Particles.position),
  reads writes(Particles.__valid)
do
  var acc = int64(0)
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      var min_x = Grid_xOrigin
      var max_x = Grid_xOrigin+Grid_xWidth
      var min_y = Grid_yOrigin
      var max_y = Grid_yOrigin+Grid_yWidth
      var min_z = Grid_zOrigin
      var max_z = Grid_zOrigin+Grid_zWidth
      var pos = Particles[p].position
      if pos[0]>max_x or pos[0]<min_x or pos[1]>max_y or pos[1]<min_y or pos[2]>max_z or pos[2]<min_z then
        Particles[p].__valid = false
        acc += (-1)
      end
    end
  end
  return acc
end

-------------------------------------------------------------------------------
-- MAIN SIMULATION
-------------------------------------------------------------------------------

local function mkInstance() local INSTANCE = {}

  local DOM_INST = DOM.mkInstance()

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
  local BC = {
    xPosSign = regentlib.newsymbol(double[3]),
    xNegSign = regentlib.newsymbol(double[3]),
    xPosVelocity = regentlib.newsymbol(double[3]),
    xNegVelocity = regentlib.newsymbol(double[3]),
    xPosTemperature = regentlib.newsymbol(double),
    xNegTemperature = regentlib.newsymbol(double),
    yPosSign = regentlib.newsymbol(double[3]),
    yNegSign = regentlib.newsymbol(double[3]),
    yPosVelocity = regentlib.newsymbol(double[3]),
    yNegVelocity = regentlib.newsymbol(double[3]),
    yPosTemperature = regentlib.newsymbol(double),
    yNegTemperature = regentlib.newsymbol(double),
    zPosSign = regentlib.newsymbol(double[3]),
    zNegSign = regentlib.newsymbol(double[3]),
    zPosVelocity = regentlib.newsymbol(double[3]),
    zNegVelocity = regentlib.newsymbol(double[3]),
    zPosTemperature = regentlib.newsymbol(double),
    zNegTemperature = regentlib.newsymbol(double),
    xBCParticles = regentlib.newsymbol(SCHEMA.ParticlesBC),
    yBCParticles = regentlib.newsymbol(SCHEMA.ParticlesBC),
    zBCParticles = regentlib.newsymbol(SCHEMA.ParticlesBC),
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
  local Fluid_copy = regentlib.newsymbol()
  local Particles = regentlib.newsymbol()
  local Particles_copy = regentlib.newsymbol()
  local TradeQueue = UTIL.generate(26, regentlib.newsymbol)
  local Radiation = regentlib.newsymbol()
  local tiles = regentlib.newsymbol()
  local p_Fluid = regentlib.newsymbol()
  local p_Fluid_copy = regentlib.newsymbol()
  local p_Particles = regentlib.newsymbol()
  local p_Particles_copy = regentlib.newsymbol()
  local p_TradeQueue = UTIL.generate(26, regentlib.newsymbol)
  local p_Radiation = regentlib.newsymbol()

  -----------------------------------------------------------------------------
  -- Exported symbols
  -----------------------------------------------------------------------------

  INSTANCE.Grid = Grid
  INSTANCE.Integrator_deltaTime = Integrator_deltaTime
  INSTANCE.Integrator_exitCond = Integrator_exitCond
  INSTANCE.Fluid = Fluid
  INSTANCE.Fluid_copy = Fluid_copy
  INSTANCE.Particles = Particles
  INSTANCE.Particles_copy = Particles_copy
  INSTANCE.Radiation = Radiation
  INSTANCE.tiles = tiles
  INSTANCE.p_Fluid = p_Fluid
  INSTANCE.p_Fluid_copy = p_Fluid_copy
  INSTANCE.p_Particles = p_Particles
  INSTANCE.p_Particles_copy = p_Particles_copy
  INSTANCE.p_Radiation = p_Radiation

  -----------------------------------------------------------------------------
  -- Symbol declaration & initialization
  -----------------------------------------------------------------------------

  function INSTANCE.DeclSymbols(config) return rquote

    ---------------------------------------------------------------------------
    -- Preparation
    ---------------------------------------------------------------------------

    -- Start timer
    var [startTime] = regentlib.c.legion_get_current_time_in_micros() / 1000;

    -- Write console header
    [emitConsoleWrite(config, 'Iter\t'..
                              'Sim Time\t'..
                              'Wall t\t'..
                              'Delta Time\t'..
                              'Avg Press\t'..
                              'Avg Temp\t'..
                              'Average KE\t'..
                              '#Part\t'..
                              'Particle T\n')];

    ---------------------------------------------------------------------------
    -- Declare & initialize state variables
    ---------------------------------------------------------------------------

    -- Cell step size (TODO: Change when we go to non-uniform meshes)
    var [Grid.xCellWidth] = config.Grid.xWidth / config.Grid.xNum
    var [Grid.yCellWidth] = config.Grid.yWidth / config.Grid.yNum
    var [Grid.zCellWidth] = config.Grid.zWidth / config.Grid.zNum
    var [Grid.cellVolume] = Grid.xCellWidth * Grid.yCellWidth * Grid.zCellWidth

    var [BC.xPosSign]
    var [BC.xNegSign]
    var [BC.xPosVelocity]
    var [BC.xNegVelocity]
    var [BC.xPosTemperature]
    var [BC.xNegTemperature]

    var [BC.yPosSign]
    var [BC.yNegSign]
    var [BC.yPosVelocity]
    var [BC.yNegVelocity]
    var [BC.yPosTemperature]
    var [BC.yNegTemperature]

    var [BC.zPosSign]
    var [BC.zNegSign]
    var [BC.zPosVelocity]
    var [BC.zNegVelocity]
    var [BC.zPosTemperature]
    var [BC.zNegTemperature]

    var [BC.xBCParticles]
    var [BC.yBCParticles]
    var [BC.zBCParticles]

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

    if config.Radiation.type == SCHEMA.RadiationModel_DOM then
      regentlib.assert(config.Grid.xNum >= config.Radiation.u.DOM.xNum and
                       config.Grid.yNum >= config.Radiation.u.DOM.yNum and
                       config.Grid.zNum >= config.Radiation.u.DOM.zNum,
                       'Radiation grid cannnot be finer than fluid grid')
      regentlib.assert(config.Grid.xNum % config.Radiation.u.DOM.xNum == 0 and
                       config.Grid.yNum % config.Radiation.u.DOM.yNum == 0 and
                       config.Grid.zNum % config.Radiation.u.DOM.zNum == 0,
                       'Inexact radiation grid coarsening factor')
    end

    -- Set up flow BC's in x direction
    if ((config.BC.xBCLeft == SCHEMA.FlowBC_Periodic) and (config.BC.xBCRight == SCHEMA.FlowBC_Periodic)) then
      BC.xBCParticles = SCHEMA.ParticlesBC_Periodic
    elseif ((config.BC.xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow) and (config.BC.xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow)) then
      if config.BC.xBCLeftHeat.type == SCHEMA.TempProfile_Constant then
        -- Do nothing
      elseif config.BC.xBCLeftHeat.type == SCHEMA.TempProfile_Parabola then
        regentlib.assert(false, 'Parabola heat model not supported')
      elseif config.BC.xBCLeftHeat.type == SCHEMA.TempProfile_Incoming then
        -- Do nothing
      else regentlib.assert(false, 'Unhandled case in switch') end
      BC.xBCParticles = SCHEMA.ParticlesBC_Disappear
    else
      if (config.BC.xBCLeft == SCHEMA.FlowBC_Symmetry) then
        BC.xNegSign = array(-1.0, 1.0, 1.0)
        BC.xNegVelocity = array(0.0, 0.0, 0.0)
        BC.xNegTemperature = -1.0
        BC.xBCParticles = SCHEMA.ParticlesBC_Bounce
      elseif (config.BC.xBCLeft == SCHEMA.FlowBC_AdiabaticWall) then
        BC.xNegSign = array(-1.0, -1.0, -1.0)
        BC.xNegVelocity = vs_mul(config.BC.xBCLeftVel, 2.0)
        BC.xNegTemperature = -1.0
        BC.xBCParticles = SCHEMA.ParticlesBC_Bounce
      elseif (config.BC.xBCLeft == SCHEMA.FlowBC_IsothermalWall) then
        BC.xNegSign = array(-1.0, -1.0, -1.0)
        BC.xNegVelocity = vs_mul(config.BC.xBCLeftVel, 2.0)
        if config.BC.xBCLeftHeat.type == SCHEMA.TempProfile_Constant then
          BC.xNegTemperature = config.BC.xBCLeftHeat.u.Constant.temperature
        else
          regentlib.assert(false, 'Only constant heat model supported')
        end
        BC.xBCParticles = SCHEMA.ParticlesBC_Bounce
      else
        regentlib.assert(false, "Boundary conditions in xBCLeft not implemented")
      end

      if (config.BC.xBCRight == SCHEMA.FlowBC_Symmetry) then
        BC.xPosSign = array(-1.0, 1.0, 1.0)
        BC.xPosVelocity = array(0.0, 0.0, 0.0)
        BC.xPosTemperature = -1.0
        BC.xBCParticles = SCHEMA.ParticlesBC_Bounce
      elseif (config.BC.xBCRight == SCHEMA.FlowBC_AdiabaticWall) then
        BC.xPosSign = array(-1.0, -1.0, -1.0)
        BC.xPosVelocity = vs_mul(config.BC.xBCRightVel, 2.0)
        BC.xPosTemperature = -1.0
        BC.xBCParticles = SCHEMA.ParticlesBC_Bounce
      elseif (config.BC.xBCRight == SCHEMA.FlowBC_IsothermalWall) then
        BC.xPosSign = array(-1.0, -1.0, -1.0)
        BC.xPosVelocity = vs_mul(config.BC.xBCRightVel, 2.0)
        if config.BC.xBCRightHeat.type == SCHEMA.TempProfile_Constant then
          BC.xPosTemperature = config.BC.xBCRightHeat.u.Constant.temperature
        else
          regentlib.assert(false, 'Only constant heat model supported')
        end
        BC.xBCParticles = SCHEMA.ParticlesBC_Bounce
      else
        regentlib.assert(false, "Boundary conditions in xBCRight not implemented")
      end
    end

    -- Set up flow BC's in y direction
    if ((config.BC.yBCLeft == SCHEMA.FlowBC_Periodic) and (config.BC.yBCRight == SCHEMA.FlowBC_Periodic)) then
      BC.yBCParticles = SCHEMA.ParticlesBC_Periodic
    else
      if (config.BC.yBCLeft == SCHEMA.FlowBC_Symmetry) then
        BC.yNegSign = array(1.0, -1.0, 1.0)
        BC.yNegVelocity = array(0.0, 0.0, 0.0)
        BC.yNegTemperature = -1.0
        BC.yBCParticles = SCHEMA.ParticlesBC_Bounce
      elseif (config.BC.yBCLeft == SCHEMA.FlowBC_AdiabaticWall) then
        BC.yNegSign = array(-1.0, -1.0, -1.0)
        BC.yNegVelocity = vs_mul(config.BC.yBCLeftVel, 2.0)
        BC.yNegTemperature = -1.0
        BC.yBCParticles = SCHEMA.ParticlesBC_Bounce
      elseif (config.BC.yBCLeft == SCHEMA.FlowBC_IsothermalWall) then
        BC.yNegSign = array(-1.0, -1.0, -1.0)
        BC.yNegVelocity = vs_mul(config.BC.yBCLeftVel, 2.0)
        if config.BC.yBCLeftHeat.type == SCHEMA.TempProfile_Constant then
          BC.yNegTemperature = config.BC.yBCLeftHeat.u.Constant.temperature
        else
          regentlib.assert(false, 'Only constant heat model supported')
        end
        BC.yBCParticles = SCHEMA.ParticlesBC_Bounce
      elseif (config.BC.yBCLeft == SCHEMA.FlowBC_NonUniformTemperatureWall) then
        BC.yNegSign = array(-1.0, -1.0, -1.0)
        BC.yNegVelocity = vs_mul(config.BC.yBCLeftVel, 2.0)
        if not (config.BC.yBCLeftHeat.type == SCHEMA.TempProfile_Parabola) then
          regentlib.assert(false, 'Only parabola heat model supported')
        end
        BC.yBCParticles = SCHEMA.ParticlesBC_Bounce
      else
        regentlib.assert(false, "Boundary conditions in y not implemented")
      end

      if (config.BC.yBCRight == SCHEMA.FlowBC_Symmetry) then
        BC.yPosSign = array(1.0, -1.0, 1.0)
        BC.yPosVelocity = array(0.0, 0.0, 0.0)
        BC.yPosTemperature = -1.0
        BC.yBCParticles = SCHEMA.ParticlesBC_Bounce
      elseif (config.BC.yBCRight == SCHEMA.FlowBC_AdiabaticWall) then
        BC.yPosSign = array(-1.0, -1.0, -1.0)
        BC.yPosVelocity = vs_mul(config.BC.yBCRightVel, 2.0)
        BC.yPosTemperature = -1.0
        BC.yBCParticles = SCHEMA.ParticlesBC_Bounce
      elseif (config.BC.yBCRight == SCHEMA.FlowBC_IsothermalWall) then
        BC.yPosSign = array(-1.0, -1.0, -1.0)
        BC.yPosVelocity = vs_mul(config.BC.yBCRightVel, 2.0)
        if config.BC.yBCRightHeat.type == SCHEMA.TempProfile_Constant then
          BC.yPosTemperature = config.BC.yBCRightHeat.u.Constant.temperature
        else
          regentlib.assert(false, 'Only constant heat model supported')
        end
        BC.yBCParticles = SCHEMA.ParticlesBC_Bounce
      elseif (config.BC.yBCRight == SCHEMA.FlowBC_NonUniformTemperatureWall) then
        BC.yPosSign = array(-1.0, -1.0, -1.0)
        BC.yPosVelocity = vs_mul(config.BC.yBCRightVel, 2.0)
        if not (config.BC.yBCRightHeat.type == SCHEMA.TempProfile_Parabola) then
          regentlib.assert(false, 'Only parabola heat model supported')
        end
        BC.yBCParticles = SCHEMA.ParticlesBC_Bounce
      else
        regentlib.assert(false, "Boundary conditions in y not implemented")
      end
    end

    -- Set up flow BC's in z direction
    if ((config.BC.zBCLeft == SCHEMA.FlowBC_Periodic) and (config.BC.zBCRight == SCHEMA.FlowBC_Periodic)) then
      BC.zBCParticles = SCHEMA.ParticlesBC_Periodic
    else
      if (config.BC.zBCLeft == SCHEMA.FlowBC_Symmetry) then
        BC.zNegSign = array(1.0, 1.0, -1.0)
        BC.zNegVelocity = array(0.0, 0.0, 0.0)
        BC.zNegTemperature = -1.0
        BC.zBCParticles = SCHEMA.ParticlesBC_Bounce
      elseif (config.BC.zBCLeft == SCHEMA.FlowBC_AdiabaticWall) then
        BC.zNegSign = array(-1.0, -1.0, -1.0)
        BC.zNegVelocity = vs_mul(config.BC.zBCLeftVel, 2.0)
        BC.zNegTemperature = -1.0
        BC.zBCParticles = SCHEMA.ParticlesBC_Bounce
      elseif (config.BC.zBCLeft == SCHEMA.FlowBC_IsothermalWall) then
        BC.zNegSign = array(-1.0, -1.0, -1.0)
        BC.zNegVelocity = vs_mul(config.BC.zBCLeftVel, 2.0)
        if config.BC.zBCLeftHeat.type == SCHEMA.TempProfile_Constant then
          BC.zNegTemperature = config.BC.zBCLeftHeat.u.Constant.temperature
        else
          regentlib.assert(false, 'Only constant heat model supported')
        end
        BC.zBCParticles = SCHEMA.ParticlesBC_Bounce
      elseif (config.BC.zBCLeft == SCHEMA.FlowBC_NonUniformTemperatureWall) then
        BC.zNegSign = array(-1.0, -1.0, -1.0)
        BC.zNegVelocity = vs_mul(config.BC.zBCLeftVel, 2.0)
        if not (config.BC.zBCLeftHeat.type == SCHEMA.TempProfile_Parabola) then
          regentlib.assert(false, 'Only parabola heat model supported')
        end
        BC.zBCParticles = SCHEMA.ParticlesBC_Bounce
      else
        regentlib.assert(false, "Boundary conditions in zBCLeft not implemented")
      end

      if (config.BC.zBCRight == SCHEMA.FlowBC_Symmetry) then
        BC.zPosSign = array(1.0, 1.0, -1.0)
        BC.zPosVelocity = array(0.0, 0.0, 0.0)
        BC.zPosTemperature = -1.0
        BC.zBCParticles = SCHEMA.ParticlesBC_Bounce
      elseif (config.BC.zBCRight == SCHEMA.FlowBC_AdiabaticWall) then
        BC.zPosSign = array(-1.0, -1.0, -1.0)
        BC.zPosVelocity = vs_mul(config.BC.zBCRightVel, 2.0)
        BC.zPosTemperature = -1.0
        BC.zBCParticles = SCHEMA.ParticlesBC_Bounce
      elseif (config.BC.zBCRight == SCHEMA.FlowBC_IsothermalWall) then
        BC.zPosSign = array(-1.0, -1.0, -1.0)
        BC.zPosVelocity = vs_mul(config.BC.zBCRightVel, 2.0)
        if config.BC.zBCRightHeat.type == SCHEMA.TempProfile_Constant then
          BC.zPosTemperature = config.BC.zBCRightHeat.u.Constant.temperature
        else
          regentlib.assert(false, 'Only constant heat model supported')
        end
        BC.zBCParticles = SCHEMA.ParticlesBC_Bounce
      elseif (config.BC.zBCRight == SCHEMA.FlowBC_NonUniformTemperatureWall) then
        BC.zPosSign = array(-1.0, -1.0, -1.0)
        BC.zPosVelocity = vs_mul(config.BC.zBCRightVel, 2.0)
        if not (config.BC.zBCRightHeat.type == SCHEMA.TempProfile_Parabola) then
          regentlib.assert(false, 'Only parabola heat model supported')
        end
        BC.zBCParticles = SCHEMA.ParticlesBC_Bounce
      else
        regentlib.assert(false, "Boundary conditions in zBCRight not implemented")
      end
    end

    -- Check if boundary conditions in each direction are either both periodic or both non-periodic
    if (not (((config.BC.xBCLeft == SCHEMA.FlowBC_Periodic) and (config.BC.xBCRight == SCHEMA.FlowBC_Periodic)) or ((not (config.BC.xBCLeft == SCHEMA.FlowBC_Periodic)) and (not (config.BC.xBCRight == SCHEMA.FlowBC_Periodic))))) then
      regentlib.assert(false, "Boundary conditions in x should match for periodicity")
    end
    if (not (((config.BC.yBCLeft == SCHEMA.FlowBC_Periodic) and (config.BC.yBCRight == SCHEMA.FlowBC_Periodic)) or ((not (config.BC.yBCLeft == SCHEMA.FlowBC_Periodic)) and (not (config.BC.yBCRight == SCHEMA.FlowBC_Periodic))))) then
      regentlib.assert(false, "Boundary conditions in y should match for periodicity")
    end
    if (not (((config.BC.zBCLeft == SCHEMA.FlowBC_Periodic) and (config.BC.zBCRight == SCHEMA.FlowBC_Periodic)) or ((not (config.BC.zBCLeft == SCHEMA.FlowBC_Periodic)) and (not (config.BC.zBCRight == SCHEMA.FlowBC_Periodic))))) then
      regentlib.assert(false, "Boundary conditions in z should match for periodicity")
    end

    ---------------------------------------------------------------------------
    -- Create Regions and Partitions
    ---------------------------------------------------------------------------

    var sampleId = config.Mapping.sampleId

    -- Create Fluid Regions
    var is_Fluid = ispace(int3d, {x = config.Grid.xNum + 2*Grid.xBnum,
                                  y = config.Grid.yNum + 2*Grid.yBnum,
                                  z = config.Grid.zNum + 2*Grid.zBnum})
    var [Fluid] = region(is_Fluid, Fluid_columns);
    [UTIL.emitRegionTagAttach(Fluid, MAPPER.SAMPLE_ID_TAG, sampleId, int)];
    var [Fluid_copy] = region(is_Fluid, Fluid_columns);
    [UTIL.emitRegionTagAttach(Fluid_copy, MAPPER.SAMPLE_ID_TAG, sampleId, int)];

    -- Create Particles Regions
    var maxParticlesPerTile = ceil((config.Particles.maxNum/numTiles)*config.Particles.maxSkew)
    var is_Particles = ispace(int1d, maxParticlesPerTile * numTiles)
    var [Particles] = region(is_Particles, Particles_columns);
    [UTIL.emitRegionTagAttach(Particles, MAPPER.SAMPLE_ID_TAG, sampleId, int)];
    var [Particles_copy] = region(is_Particles, Particles_columns);
    [UTIL.emitRegionTagAttach(Particles_copy, MAPPER.SAMPLE_ID_TAG, sampleId, int)];
    var is_TradeQueue = ispace(int1d, config.Particles.maxXferNum * numTiles);
    @ESCAPE for k = 1,26 do @EMIT
      var [TradeQueue[k]] = region(is_TradeQueue, TradeQueue_columns);
      [UTIL.emitRegionTagAttach(TradeQueue[k], MAPPER.SAMPLE_ID_TAG, sampleId, int)];
    @TIME end @EPACSE

    -- Create Radiation Regions
    var rad_x = NX
    var rad_y = NY
    var rad_z = NZ
    if config.Radiation.type == SCHEMA.RadiationModel_DOM then
      rad_x = config.Radiation.u.DOM.xNum
      rad_y = config.Radiation.u.DOM.yNum
      rad_z = config.Radiation.u.DOM.zNum
    end
    var is_Radiation = ispace(int3d, {x = rad_x, y = rad_y, z = rad_z})
    var [Radiation] = region(is_Radiation, Radiation_columns);
    [UTIL.emitRegionTagAttach(Radiation, MAPPER.SAMPLE_ID_TAG, sampleId, int)];

    -- Partitioning domain
    var [tiles] = ispace(int3d, {NX,NY,NZ})

    -- Fluid Partitioning
    regentlib.assert(config.Grid.xNum % NX == 0, "Uneven partitioning of fluid grid on x")
    regentlib.assert(config.Grid.yNum % NY == 0, "Uneven partitioning of fluid grid on y")
    regentlib.assert(config.Grid.zNum % NZ == 0, "Uneven partitioning of fluid grid on z")
    var coloring_Fluid = regentlib.c.legion_domain_point_coloring_create()
    for c in tiles do
      var rect = rect3d{lo = int3d{x = Grid.xBnum+(config.Grid.xNum/NX)*c.x,       y = Grid.yBnum+(config.Grid.yNum/NY)*c.y,       z = Grid.zBnum+(config.Grid.zNum/NZ)*c.z      },
                        hi = int3d{x = Grid.xBnum+(config.Grid.xNum/NX)*(c.x+1)-1, y = Grid.yBnum+(config.Grid.yNum/NY)*(c.y+1)-1, z = Grid.zBnum+(config.Grid.zNum/NZ)*(c.z+1)-1}}
      if c.x == 0 then
        rect.lo.x -= Grid.xBnum
      end
      if c.x == NX-1 then
        rect.hi.x += Grid.xBnum
      end
      if c.y == 0 then
        rect.lo.y -= Grid.yBnum
      end
      if c.y == NY-1 then
        rect.hi.y += Grid.yBnum
      end
      if c.z == 0 then
        rect.lo.z -= Grid.zBnum
      end
      if c.z == NZ-1 then
        rect.hi.z += Grid.zBnum
      end
      regentlib.c.legion_domain_point_coloring_color_domain(coloring_Fluid, c, rect)
    end
    var [p_Fluid] = partition(disjoint, Fluid, coloring_Fluid, tiles)
    var [p_Fluid_copy] = partition(disjoint, Fluid_copy, coloring_Fluid, tiles)
    regentlib.c.legion_domain_point_coloring_destroy(coloring_Fluid)

    -- Particles Partitioning
    regentlib.assert(config.Particles.maxNum % numTiles == 0, "Uneven partitioning of particles")
    var coloring_Particles = regentlib.c.legion_domain_point_coloring_create()
    for z = 0,NZ do
      for y = 0,NY do
        for x = 0,NX do
          var rBase : int64
          for rStart in Particles do
            rBase = rStart + (z*NX*NY + y*NX + x) * maxParticlesPerTile
            break
          end
          regentlib.c.legion_domain_point_coloring_color_domain(coloring_Particles, int3d{x,y,z}, rect1d{rBase,rBase+maxParticlesPerTile-1})
        end
      end
    end
    var [p_Particles] = partition(disjoint, Particles, coloring_Particles, tiles)
    var [p_Particles_copy] = partition(disjoint, Particles_copy, coloring_Particles, tiles)
    regentlib.c.legion_domain_point_coloring_destroy(coloring_Particles)
    var coloring_TradeQueue = regentlib.c.legion_domain_point_coloring_create()
    for z = 0,NZ do
      for y = 0,NY do
        for x = 0,NX do
          var rBase : int64
          for rStart in [TradeQueue[1]] do
            rBase = rStart + (z*NX*NY + y*NX + x) * config.Particles.maxXferNum
            break
          end
          regentlib.c.legion_domain_point_coloring_color_domain(coloring_TradeQueue, int3d{x,y,z}, rect1d{rBase,rBase+config.Particles.maxXferNum-1})
        end
      end
    end
    @ESCAPE for k = 1,26 do @EMIT
      var [p_TradeQueue[k]] = partition(disjoint, [TradeQueue[k]], coloring_TradeQueue, tiles)
    @TIME end @EPACSE
    regentlib.c.legion_domain_point_coloring_destroy(coloring_TradeQueue)

    -- Radiation Partitioning
    var [p_Radiation] =
      [UTIL.mkPartitionEqually(int3d, int3d, Radiation_columns)]
      (Radiation, tiles);

    ---------------------------------------------------------------------------
    -- DOM code declarations
    ---------------------------------------------------------------------------

    [DOM_INST.DeclSymbols(config, tiles)];

  end end -- DeclSymbols

  -----------------------------------------------------------------------------
  -- Region initialization
  -----------------------------------------------------------------------------

  function INSTANCE.InitRegions(config) return rquote

    Particles_initValidField(Particles)
    if config.Radiation.type == SCHEMA.RadiationModel_DOM then
      SetCoarseningField(Fluid,
                         Grid.xBnum, config.Grid.xNum,
                         Grid.yBnum, config.Grid.yNum,
                         Grid.zBnum, config.Grid.zNum,
                         config.Radiation.u.DOM.xNum,
                         config.Radiation.u.DOM.yNum,
                         config.Radiation.u.DOM.zNum)
    end
    Flow_InitializeCell(Fluid)
    Flow_InitializeCenterCoordinates(Fluid,
                                     Grid.xBnum, config.Grid.xNum, config.Grid.origin[0], config.Grid.xWidth,
                                     Grid.yBnum, config.Grid.yNum, config.Grid.origin[1], config.Grid.yWidth,
                                     Grid.zBnum, config.Grid.zNum, config.Grid.origin[2], config.Grid.zWidth)

    if config.Flow.initCase == SCHEMA.FlowInitCase_Uniform then
      Flow_InitializeUniform(Fluid, config.Flow.initParams)
    elseif config.Flow.initCase == SCHEMA.FlowInitCase_Random then
      Flow_InitializeRandom(Fluid, config.Flow.initParams)
    elseif config.Flow.initCase == SCHEMA.FlowInitCase_TaylorGreen2DVortex then
      Flow_InitializeTaylorGreen2D(Fluid,
                                   config.Flow.initParams,
                                   Grid.xBnum, config.Grid.xNum, config.Grid.origin[0], config.Grid.xWidth,
                                   Grid.yBnum, config.Grid.yNum, config.Grid.origin[1], config.Grid.yWidth,
                                   Grid.zBnum, config.Grid.zNum, config.Grid.origin[2], config.Grid.zWidth)
    elseif config.Flow.initCase == SCHEMA.FlowInitCase_TaylorGreen3DVortex then
      Flow_InitializeTaylorGreen3D(Fluid,
                                   config.Flow.initParams,
                                   Grid.xBnum, config.Grid.xNum, config.Grid.origin[0], config.Grid.xWidth,
                                   Grid.yBnum, config.Grid.yNum, config.Grid.origin[1], config.Grid.yWidth,
                                   Grid.zBnum, config.Grid.zNum, config.Grid.origin[2], config.Grid.zWidth)
    elseif config.Flow.initCase == SCHEMA.FlowInitCase_Perturbed then
      Flow_InitializePerturbed(Fluid, config.Flow.initParams)
    elseif config.Flow.initCase == SCHEMA.FlowInitCase_Restart then
      Fluid_load(tiles, config.Flow.restartDir, Fluid, Fluid_copy, p_Fluid, p_Fluid_copy)
      Integrator_timeStep = config.Integrator.restartIter
      Integrator_simTime = config.Integrator.restartTime
    else regentlib.assert(false, 'Unhandled case in switch') end

    -- initialize ghost cells to their specified values in NSCBC case
    if ((config.BC.xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow) and (config.BC.xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow)) then
        Flow_InitializeGhostNSCBC(Fluid,
                                  config,
                                  config.Flow.gasConstant,
                                  config.Flow.constantVisc,
                                  config.Flow.powerlawTempRef, config.Flow.powerlawViscRef,
                                  config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef,
                                  config.Flow.viscosityModel,
                                  Grid.xBnum, config.Grid.xNum,
                                  Grid.yBnum, config.Grid.yNum,
                                  Grid.zBnum, config.Grid.zNum)
    end

    -- update interior cells from initialized primitive values
    Flow_UpdateConservedFromPrimitive(Fluid,
                                      config.Flow.gamma,
                                      config.Flow.gasConstant,
                                      Grid.xBnum, config.Grid.xNum,
                                      Grid.yBnum, config.Grid.yNum,
                                      Grid.zBnum, config.Grid.zNum)
    if ((config.BC.xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow) and (config.BC.xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow)) then
      Flow_UpdateConservedFromPrimitiveGhostNSCBC(Fluid,
                                                  config,
                                                  config.Flow.gamma, config.Flow.gasConstant,
                                                  Grid.xBnum, config.Grid.xNum,
                                                  Grid.yBnum, config.Grid.yNum,
                                                  Grid.zBnum, config.Grid.zNum)
    end

    Flow_UpdateAuxiliaryVelocity(Fluid, Grid.xBnum, config.Grid.xNum, Grid.yBnum, config.Grid.yNum, Grid.zBnum, config.Grid.zNum)
    if ((config.BC.xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow) and (config.BC.xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow)) then
      Flow_UpdateAuxiliaryVelocityGhostNSCBC(Fluid,
                                             config,
                                             config.Flow.constantVisc,
                                             config.Flow.powerlawTempRef, config.Flow.powerlawViscRef,
                                             config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef,
                                             config.Flow.viscosityModel,
                                             Grid.xBnum, config.Grid.xNum,
                                             Grid.yBnum, config.Grid.yNum,
                                             Grid.zBnum, config.Grid.zNum)
    end
    Flow_UpdateGhostVelocityStep1(Fluid,
                                  config,
                                  BC.xNegVelocity, BC.xPosVelocity, BC.xNegSign, BC.xPosSign,
                                  BC.yNegVelocity, BC.yPosVelocity, BC.yNegSign, BC.yPosSign,
                                  BC.zNegVelocity, BC.zPosVelocity, BC.zNegSign, BC.zPosSign,
                                  Grid.xBnum, config.Grid.xNum,
                                  Grid.yBnum, config.Grid.yNum,
                                  Grid.zBnum, config.Grid.zNum)
    Flow_UpdateGhostVelocityStep2(Fluid,
                                  config,
                                  Grid.xBnum, config.Grid.xNum,
                                  Grid.yBnum, config.Grid.yNum,
                                  Grid.zBnum, config.Grid.zNum)

    Flow_UpdateAuxiliaryThermodynamics(Fluid, config.Flow.gamma, config.Flow.gasConstant, Grid.xBnum, config.Grid.xNum, Grid.yBnum, config.Grid.yNum, Grid.zBnum, config.Grid.zNum)
    if ((config.BC.xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow) and (config.BC.xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow)) then
      Flow_UpdateAuxiliaryThermodynamicsGhostNSCBC(Fluid,
                                                   config,
                                                   config.Flow.gamma,
                                                   config.Flow.gasConstant,
                                                   Grid.xBnum, config.Grid.xNum,
                                                   Grid.yBnum, config.Grid.yNum,
                                                   Grid.zBnum, config.Grid.zNum)
    end
    Flow_UpdateGhostThermodynamicsStep1(Fluid,
                                        config,
                                        config.Flow.gamma,
                                        config.Flow.gasConstant,
                                        BC.xNegTemperature, BC.xPosTemperature,
                                        BC.yNegTemperature, BC.yPosTemperature,
                                        BC.zNegTemperature, BC.zPosTemperature,
                                        Grid.xBnum, config.Grid.xNum,
                                        Grid.yBnum, config.Grid.yNum,
                                        Grid.zBnum, config.Grid.zNum)
    Flow_UpdateGhostThermodynamicsStep2(Fluid,
                                        config,
                                        Grid.xBnum, config.Grid.xNum,
                                        Grid.yBnum, config.Grid.yNum,
                                        Grid.zBnum, config.Grid.zNum)

    Flow_UpdateGhostConservedStep1(Fluid,
                                   config,
                                   BC.xNegTemperature, BC.xNegVelocity, BC.xPosTemperature, BC.xPosVelocity, BC.xNegSign, BC.xPosSign,
                                   BC.yNegTemperature, BC.yNegVelocity, BC.yPosTemperature, BC.yPosVelocity, BC.yNegSign, BC.yPosSign,
                                   BC.zNegTemperature, BC.zNegVelocity, BC.zPosTemperature, BC.zPosVelocity, BC.zNegSign, BC.zPosSign,
                                   config.Flow.gamma, config.Flow.gasConstant,
                                   config.Flow.constantVisc,
                                   config.Flow.powerlawTempRef, config.Flow.powerlawViscRef,
                                   config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef,
                                   config.Flow.viscosityModel,
                                   Grid.xBnum, config.Grid.xNum,
                                   Grid.yBnum, config.Grid.yNum,
                                   Grid.zBnum, config.Grid.zNum)
    Flow_UpdateGhostConservedStep2(Fluid,
                                   config,
                                   Grid.xBnum, config.Grid.xNum,
                                   Grid.yBnum, config.Grid.yNum,
                                   Grid.zBnum, config.Grid.zNum)

    -- Initialize particles
    if config.Particles.initCase == SCHEMA.ParticlesInitCase_Random then
      regentlib.assert(false, "Random particle initialization is disabled")
    elseif config.Particles.initCase == SCHEMA.ParticlesInitCase_Restart then
      Particles_load(tiles, config.Particles.restartDir, Particles, Particles_copy, p_Particles, p_Particles_copy)
      for c in tiles do
        Particles_LocateInCells(p_Particles[c],
                                Grid.xBnum, config.Grid.xNum, config.Grid.origin[0], config.Grid.xWidth,
                                Grid.yBnum, config.Grid.yNum, config.Grid.origin[1], config.Grid.yWidth,
                                Grid.zBnum, config.Grid.zNum, config.Grid.origin[2], config.Grid.zWidth)
      end
      for c in tiles do
        Particles_CheckPartitioning(c,
                                    p_Particles[c],
                                    Grid.xBnum, config.Grid.xNum, NX,
                                    Grid.yBnum, config.Grid.yNum, NY,
                                    Grid.zBnum, config.Grid.zNum, NZ)
      end
      Particles_number += Particles_CalculateNumber(Particles)
    elseif config.Particles.initCase == SCHEMA.ParticlesInitCase_Uniform then
      regentlib.assert(config.Particles.initNum <= config.Particles.maxNum,
                       "Not enough space for initial number of particles")
      InitParticlesUniform(Particles, Fluid, config, Grid.xBnum, Grid.yBnum, Grid.zBnum)
      Particles_number = (config.Particles.initNum / numTiles) * numTiles
    else regentlib.assert(false, 'Unhandled case in switch') end

    -- Initialize radiation
    if config.Radiation.type == SCHEMA.RadiationModel_OFF then
      -- Do nothing
    elseif config.Radiation.type == SCHEMA.RadiationModel_Algebraic then
      -- Do nothing
    elseif config.Radiation.type == SCHEMA.RadiationModel_DOM then
      [DOM_INST.InitRegions(config, tiles, p_Radiation)];
    else regentlib.assert(false, 'Unhandled case in switch') end

  end end -- InitRegions

  -----------------------------------------------------------------------------
  -- Main time-step loop header
  -----------------------------------------------------------------------------

  function INSTANCE.MainLoopHeader(config) return rquote

    -- Calculate exit condition
    Integrator_exitCond =
      Integrator_timeStep >= config.Integrator.maxIter

    -- Determine time step size
    if config.Integrator.cfl < 0.0 then
      Integrator_deltaTime = config.Integrator.fixedDeltaTime
    else
      var Integrator_maxConvectiveSpectralRadius = 0.0
      var Integrator_maxViscousSpectralRadius = 0.0
      var Integrator_maxHeatConductionSpectralRadius = 0.0
      var Grid_dXYZInverseSquare =
        1.0/Grid.xCellWidth/Grid.xCellWidth +
        1.0/Grid.yCellWidth/Grid.yCellWidth +
        1.0/Grid.zCellWidth/Grid.zCellWidth
      Integrator_maxConvectiveSpectralRadius max=
        CalculateConvectiveSpectralRadius(Fluid,
                                          config.Flow.gamma, config.Flow.gasConstant,
                                          Grid_dXYZInverseSquare, Grid.xCellWidth, Grid.yCellWidth, Grid.zCellWidth)
      Integrator_maxViscousSpectralRadius max=
        CalculateViscousSpectralRadius(Fluid,
                                       config.Flow.constantVisc,
                                       config.Flow.powerlawTempRef, config.Flow.powerlawViscRef,
                                       config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef,
                                       config.Flow.viscosityModel,
                                       Grid_dXYZInverseSquare)
      Integrator_maxHeatConductionSpectralRadius max=
        CalculateHeatConductionSpectralRadius(Fluid,
                                              config.Flow.constantVisc,
                                              config.Flow.gamma, config.Flow.gasConstant,
                                              config.Flow.powerlawTempRef, config.Flow.powerlawViscRef,
                                              config.Flow.prandtl,
                                              config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef,
                                              config.Flow.viscosityModel,
                                              Grid_dXYZInverseSquare)
      Integrator_deltaTime = (config.Integrator.cfl/max(Integrator_maxConvectiveSpectralRadius, max(Integrator_maxViscousSpectralRadius, Integrator_maxHeatConductionSpectralRadius)))
    end

  end end -- MainLoopHeader

  -----------------------------------------------------------------------------
  -- Per-time-step I/O
  -----------------------------------------------------------------------------

  function INSTANCE.PerformIO(config) return rquote

    -- Write to console
    var Flow_averagePressure = 0.0
    var Flow_averageTemperature = 0.0
    var Flow_averageKineticEnergy = 0.0
    var Particles_averageTemperature = 0.0
    Flow_averagePressure += CalculateAveragePressure(Fluid, Grid.cellVolume, Grid.xBnum, config.Grid.xNum, Grid.yBnum, config.Grid.yNum, Grid.zBnum, config.Grid.zNum)
    Flow_averageTemperature += CalculateAverageTemperature(Fluid, Grid.cellVolume, Grid.xBnum, config.Grid.xNum, Grid.yBnum, config.Grid.yNum, Grid.zBnum, config.Grid.zNum)
    Flow_averageKineticEnergy += CalculateAverageKineticEnergy(Fluid, Grid.cellVolume, Grid.xBnum, config.Grid.xNum, Grid.yBnum, config.Grid.yNum, Grid.zBnum, config.Grid.zNum)
    Particles_averageTemperature += Particles_IntegrateQuantities(Particles)
    Flow_averagePressure = (Flow_averagePressure/(((config.Grid.xNum*config.Grid.yNum)*config.Grid.zNum)*Grid.cellVolume))
    Flow_averageTemperature = (Flow_averageTemperature/(((config.Grid.xNum*config.Grid.yNum)*config.Grid.zNum)*Grid.cellVolume))
    Flow_averageKineticEnergy = (Flow_averageKineticEnergy/(((config.Grid.xNum*config.Grid.yNum)*config.Grid.zNum)*Grid.cellVolume))
    Particles_averageTemperature = (Particles_averageTemperature/Particles_number)
    Console_write(config,
                  Integrator_timeStep,
                  Integrator_simTime,
                  startTime,
                  Integrator_deltaTime,
                  Flow_averagePressure,
                  Flow_averageTemperature,
                  Flow_averageKineticEnergy,
                  Particles_number,
                  Particles_averageTemperature)

    -- Dump restart files
    if config.IO.wrtRestart then
      if Integrator_exitCond or Integrator_timeStep % config.IO.restartEveryTimeSteps == 0 then
        var dirname = [&int8](C.malloc(256))
        C.snprintf(dirname, 256, '%s/fluid_iter%010d', config.Mapping.outDir, Integrator_timeStep)
        Fluid_dump(tiles, dirname, Fluid, Fluid_copy, p_Fluid, p_Fluid_copy)
        C.snprintf(dirname, 256, '%s/particles_iter%010d', config.Mapping.outDir, Integrator_timeStep)
        Particles_dump(tiles, dirname, Particles, Particles_copy, p_Particles, p_Particles_copy)
        C.free(dirname)
      end
    end

    -- Write probe files
    if config.IO.probes.length > 0 then
      for c in tiles do
        Probes_write(p_Fluid[c], Integrator_exitCond, Integrator_timeStep, config)
      end
    end

  end end -- PerformIO

  -----------------------------------------------------------------------------
  -- Main time-step loop body
  -----------------------------------------------------------------------------

  function INSTANCE.MainLoopBody(config, CopyQueue) return rquote

    -- Feed particles
    if config.Particles.feeding.type == SCHEMA.FeedModel_OFF then
      -- Do nothing
    elseif config.Particles.feeding.type == SCHEMA.FeedModel_Incoming then
      for c in tiles do
        Particles_number +=
          CopyQueue_pull(c,
                         p_Particles[c],
                         CopyQueue,
                         config,
                         Grid.xBnum, Grid.yBnum, Grid.zBnum)
      end
    else regentlib.assert(false, 'Unhandled case in switch') end

    -- Set iteration-specific fields that persist across RK4 sub-steps
    Flow_InitializeTemporaries(Fluid)
    Particles_InitializeTemporaries(Particles)

    -- RK4 sub-time-stepping loop
    var Integrator_time_old = Integrator_simTime
    for Integrator_stage = 1,5 do

      Flow_ComputeVelocityGradientAll(Fluid,
                                      Grid.xBnum, Grid.xCellWidth, config.Grid.xNum,
                                      Grid.yBnum, Grid.yCellWidth, config.Grid.yNum,
                                      Grid.zBnum, Grid.zCellWidth, config.Grid.zNum)
      if ((config.BC.xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow) and (config.BC.xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow)) then
        Flow_ComputeVelocityGradientGhostNSCBC(Fluid,
                                               config,
                                               Grid.xBnum, Grid.xCellWidth, config.Grid.xNum,
                                               Grid.yBnum, Grid.yCellWidth, config.Grid.yNum,
                                               Grid.zBnum, Grid.zCellWidth, config.Grid.zNum)
      end
      Flow_UpdateGhostVelocityGradientStep1(Fluid,
                                            config,
                                            BC.xNegSign, BC.yNegSign, BC.zNegSign,
                                            BC.xPosSign, BC.yPosSign, BC.zPosSign,
                                            Grid.xBnum, Grid.xCellWidth, config.Grid.xNum,
                                            Grid.yBnum, Grid.yCellWidth, config.Grid.yNum,
                                            Grid.zBnum, Grid.zCellWidth, config.Grid.zNum)
      Flow_UpdateGhostVelocityGradientStep2(Fluid,
                                            config,
                                            Grid.xBnum, config.Grid.xNum,
                                            Grid.yBnum, config.Grid.yNum,
                                            Grid.zBnum, config.Grid.zNum)

      Flow_InitializeTimeDerivatives(Fluid)
      Particles_InitializeTimeDerivatives(Particles)

      Flow_AddGetFlux(Fluid,
                      config,
                      config.Flow.constantVisc,
                      config.Flow.gamma, config.Flow.gasConstant,
                      config.Flow.powerlawTempRef, config.Flow.powerlawViscRef,
                      config.Flow.prandtl,
                      config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef,
                      config.Flow.viscosityModel,
                      Grid.xBnum, Grid.xCellWidth, config.Grid.xNum,
                      Grid.yBnum, Grid.yCellWidth, config.Grid.yNum,
                      Grid.zBnum, Grid.zCellWidth, config.Grid.zNum)
      Flow_AddUpdateUsingFlux(Fluid,
                              Grid.xBnum, Grid.xCellWidth, config.Grid.xNum,
                              Grid.yBnum, Grid.yCellWidth, config.Grid.yNum,
                              Grid.zBnum, Grid.zCellWidth, config.Grid.zNum)

      if ((config.BC.xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow) and (config.BC.xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow)) then
        Flow_AddGetFluxGhostNSCBC(Fluid,
                                  config,
                                  config.Flow.constantVisc,
                                  config.Flow.gamma,
                                  config.Flow.gasConstant,
                                  config.Flow.powerlawTempRef, config.Flow.powerlawViscRef,
                                  config.Flow.prandtl,
                                  config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef,
                                  config.Flow.viscosityModel,
                                  Grid.xBnum, Grid.xCellWidth, config.Grid.xNum,
                                  Grid.yBnum, Grid.yCellWidth, config.Grid.yNum,
                                  Grid.zBnum, Grid.zCellWidth, config.Grid.zNum)

        var maxMach = -math.huge
        maxMach max= CalculateMaxMachNumber(Fluid,
                                            config,
                                            config.Flow.gamma, config.Flow.gasConstant,
                                            Grid.xBnum, config.Grid.xNum,
                                            Grid.yBnum, config.Grid.yNum,
                                            Grid.zBnum, config.Grid.zNum)

        var Flow_lengthScale = config.Grid.xWidth

        Flow_AddUpdateUsingFluxGhostNSCBC(Fluid,
                                          config,
                                          config.Flow.gamma, config.Flow.gasConstant,
                                          config.Flow.prandtl,
                                          maxMach,
                                          Flow_lengthScale,
                                          config.Flow.constantVisc,
                                          config.Flow.powerlawTempRef, config.Flow.powerlawViscRef,
                                          config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef,
                                          config.Flow.viscosityModel,
                                          config.BC.xBCRightP_inf,
                                          Grid.xBnum, Grid.xCellWidth, config.Grid.xNum,
                                          Grid.yBnum, Grid.yCellWidth, config.Grid.yNum,
                                          Grid.zBnum, Grid.zCellWidth, config.Grid.zNum)
      end

      Flow_AddBodyForces(Fluid,
                         config.Flow.bodyForce,
                         Grid.xBnum, config.Grid.xNum,
                         Grid.yBnum, config.Grid.yNum,
                         Grid.zBnum, config.Grid.zNum)
      if ((config.BC.xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow) and (config.BC.xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow)) then
        Flow_AddBodyForcesGhostNSCBC(Fluid,
                                     config,
                                     config.Flow.bodyForce,
                                     Grid.xBnum, config.Grid.xNum,
                                     Grid.yBnum, config.Grid.yNum,
                                     Grid.zBnum, config.Grid.zNum)
      end

      -- Add turbulent forcing
      if config.Flow.turbForcing.type == SCHEMA.TurbForcingModel_HIT then
        var Flow_averagePD = 0.0
        var Flow_averageDissipation = 0.0
        var Flow_averageFe = 0.0
        var Flow_averageK = 0.0
        Flow_UpdatePD(Fluid, Grid.xBnum, config.Grid.xNum, Grid.yBnum, config.Grid.yNum, Grid.zBnum, config.Grid.zNum)
        Flow_ResetDissipation(Fluid)
        Flow_ComputeDissipationX(Fluid, config.Flow.constantVisc, config.Flow.powerlawTempRef, config.Flow.powerlawViscRef, config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef, config.Flow.viscosityModel, Grid.xBnum, config.Grid.xNum, Grid.xCellWidth, Grid.yBnum, config.Grid.yNum, Grid.zBnum, config.Grid.zNum)
        Flow_UpdateDissipationX(Fluid, Grid.xBnum, config.Grid.xNum, Grid.xCellWidth, Grid.yBnum, config.Grid.yNum, Grid.zBnum, config.Grid.zNum)
        Flow_ComputeDissipationY(Fluid, config.Flow.constantVisc, config.Flow.powerlawTempRef, config.Flow.powerlawViscRef, config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef, config.Flow.viscosityModel, Grid.xBnum, config.Grid.xNum, Grid.yBnum, config.Grid.yNum, Grid.yCellWidth, Grid.zBnum, config.Grid.zNum)
        Flow_UpdateDissipationY(Fluid, Grid.xBnum, config.Grid.xNum, Grid.yBnum, config.Grid.yNum, Grid.yCellWidth, Grid.zBnum, config.Grid.zNum)
        Flow_ComputeDissipationZ(Fluid, config.Flow.constantVisc, config.Flow.powerlawTempRef, config.Flow.powerlawViscRef, config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef, config.Flow.viscosityModel, Grid.xBnum, config.Grid.xNum, Grid.yBnum, config.Grid.yNum, Grid.zBnum, config.Grid.zNum, Grid.zCellWidth)
        Flow_UpdateDissipationZ(Fluid, Grid.xBnum, config.Grid.xNum, Grid.yBnum, config.Grid.yNum, Grid.zBnum, config.Grid.zNum, Grid.zCellWidth)
        Flow_averagePD += CalculateAveragePD(Fluid, Grid.cellVolume, Grid.xBnum, config.Grid.xNum, Grid.yBnum, config.Grid.yNum, Grid.zBnum, config.Grid.zNum)
        Flow_averagePD = (Flow_averagePD/(((config.Grid.xNum*config.Grid.yNum)*config.Grid.zNum)*Grid.cellVolume))
        Flow_averageDissipation += CalculateAverageDissipation(Fluid, Grid.cellVolume, Grid.xBnum, config.Grid.xNum, Grid.yBnum, config.Grid.yNum, Grid.zBnum, config.Grid.zNum)
        Flow_averageDissipation = (Flow_averageDissipation/(((config.Grid.xNum*config.Grid.yNum)*config.Grid.zNum)*Grid.cellVolume))
        Flow_averageK += CalculateAverageK(Fluid, Grid.cellVolume, Grid.xBnum, config.Grid.xNum, Grid.yBnum, config.Grid.yNum, Grid.zBnum, config.Grid.zNum)
        Flow_averageK = (Flow_averageK/(((config.Grid.xNum*config.Grid.yNum)*config.Grid.zNum)*Grid.cellVolume))
        Flow_averageFe += Flow_AddTurbulentSource(Fluid, Flow_averageDissipation, Flow_averageK, Flow_averagePD, Grid.cellVolume, Grid.xBnum, config.Grid.xNum, Grid.yBnum, config.Grid.yNum, Grid.zBnum, config.Grid.zNum, config)
        Flow_averageFe = (Flow_averageFe/(((config.Grid.xNum*config.Grid.yNum)*config.Grid.zNum)*Grid.cellVolume))
        Flow_AdjustTurbulentSource(Fluid, Flow_averageFe, Grid.xBnum, config.Grid.xNum, Grid.yBnum, config.Grid.yNum, Grid.zBnum, config.Grid.zNum)
      end

      -- Add fluid forces to particles
      Particles_AddFlowCoupling(Particles, Fluid, config.Flow.constantVisc, config.Flow.powerlawTempRef, config.Flow.powerlawViscRef, config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef, config.Flow.viscosityModel, Grid.xCellWidth, Grid.xRealOrigin, Grid.yCellWidth, Grid.yRealOrigin, Grid.zCellWidth, Grid.zRealOrigin, config.Particles.convectiveCoeff, config.Particles.heatCapacity)
      Particles_AddBodyForces(Particles, config.Particles.bodyForce)

      -- Add radiation
      if config.Radiation.type == SCHEMA.RadiationModel_OFF then
        -- Do nothing
      elseif config.Radiation.type == SCHEMA.RadiationModel_Algebraic then
        AddRadiation(Particles, config)
      elseif config.Radiation.type == SCHEMA.RadiationModel_DOM then
        fill(Radiation.acc_d2, 0.0)
        fill(Radiation.acc_d2t4, 0.0)
        for c in tiles do
          Radiation_AccumulateParticleValues(p_Particles[c], p_Fluid[c], p_Radiation[c])
        end
        var Radiation_xCellWidth = (config.Grid.xWidth/config.Radiation.u.DOM.xNum)
        var Radiation_yCellWidth = (config.Grid.yWidth/config.Radiation.u.DOM.yNum)
        var Radiation_zCellWidth = (config.Grid.zWidth/config.Radiation.u.DOM.zNum)
        var Radiation_cellVolume = Radiation_xCellWidth * Radiation_yCellWidth * Radiation_zCellWidth
        Radiation_UpdateFieldValues(Radiation,
                                    Radiation_cellVolume,
                                    config.Radiation.u.DOM.qa,
                                    config.Radiation.u.DOM.qs);
        [DOM_INST.ComputeRadiationField(config, tiles, p_Radiation)];
        for c in tiles do
          Particles_AbsorbRadiation(p_Particles[c],
                                    p_Fluid[c],
                                    p_Radiation[c],
                                    config.Particles.heatCapacity,
                                    config.Radiation.u.DOM.qa)
        end
      else regentlib.assert(false, 'Unhandled case in switch') end

      -- Add particle forces to fluid
      Flow_AddParticlesCoupling(Particles, Fluid, Grid.cellVolume)

      -- Time step
      Flow_UpdateVars(Fluid, Integrator_deltaTime, Integrator_stage)
      Particles_UpdateVars(Particles, Integrator_deltaTime, Integrator_stage)

      -- Now the new conserved variables values are used so update everything else
      Flow_UpdateAuxiliaryVelocity(Fluid, Grid.xBnum, config.Grid.xNum, Grid.yBnum, config.Grid.yNum, Grid.zBnum, config.Grid.zNum)
      if ((config.BC.xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow) and (config.BC.xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow)) then
        Flow_UpdateAuxiliaryVelocityGhostNSCBC(Fluid,
                                               config,
                                               config.Flow.constantVisc,
                                               config.Flow.powerlawTempRef, config.Flow.powerlawViscRef,
                                               config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef,
                                               config.Flow.viscosityModel,
                                               Grid.xBnum, config.Grid.xNum,
                                               Grid.yBnum, config.Grid.yNum,
                                               Grid.zBnum, config.Grid.zNum)
      end
      Flow_UpdateGhostVelocityStep1(Fluid,
                                    config,
                                    BC.xNegVelocity, BC.xPosVelocity, BC.xNegSign, BC.xPosSign,
                                    BC.yNegVelocity, BC.yPosVelocity, BC.yNegSign, BC.yPosSign,
                                    BC.zNegVelocity, BC.zPosVelocity, BC.zNegSign, BC.zPosSign,
                                    Grid.xBnum, config.Grid.xNum,
                                    Grid.yBnum, config.Grid.yNum,
                                    Grid.zBnum, config.Grid.zNum)
      Flow_UpdateGhostVelocityStep2(Fluid,
                                    config,
                                    Grid.xBnum, config.Grid.xNum,
                                    Grid.yBnum, config.Grid.yNum,
                                    Grid.zBnum, config.Grid.zNum)

      Flow_UpdateAuxiliaryThermodynamics(Fluid,
                                         config.Flow.gamma, config.Flow.gasConstant,
                                         Grid.xBnum, config.Grid.xNum,
                                         Grid.yBnum, config.Grid.yNum,
                                         Grid.zBnum, config.Grid.zNum)
      if ((config.BC.xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow) and (config.BC.xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow)) then
        Flow_UpdateAuxiliaryThermodynamicsGhostNSCBC(Fluid,
                                                     config,
                                                     config.Flow.gamma,
                                                     config.Flow.gasConstant,
                                                     Grid.xBnum, config.Grid.xNum,
                                                     Grid.yBnum, config.Grid.yNum,
                                                     Grid.zBnum, config.Grid.zNum)
      end
      Flow_UpdateGhostThermodynamicsStep1(Fluid,
                                          config,
                                          config.Flow.gamma,
                                          config.Flow.gasConstant,
                                          BC.xNegTemperature, BC.xPosTemperature,
                                          BC.yNegTemperature, BC.yPosTemperature,
                                          BC.zNegTemperature, BC.zPosTemperature,
                                          Grid.xBnum, config.Grid.xNum,
                                          Grid.yBnum, config.Grid.yNum,
                                          Grid.zBnum, config.Grid.zNum)
      Flow_UpdateGhostThermodynamicsStep2(Fluid,
                                          config,
                                          Grid.xBnum, config.Grid.xNum,
                                          Grid.yBnum, config.Grid.yNum,
                                          Grid.zBnum, config.Grid.zNum)

      Flow_UpdateGhostConservedStep1(Fluid,
                                     config,
                                     BC.xNegTemperature, BC.xNegVelocity, BC.xPosTemperature, BC.xPosVelocity, BC.xNegSign, BC.xPosSign,
                                     BC.yNegTemperature, BC.yNegVelocity, BC.yPosTemperature, BC.yPosVelocity, BC.yNegSign, BC.yPosSign,
                                     BC.zNegTemperature, BC.zNegVelocity, BC.zPosTemperature, BC.zPosVelocity, BC.zNegSign, BC.zPosSign,
                                     config.Flow.gamma, config.Flow.gasConstant,
                                     config.Flow.constantVisc,
                                     config.Flow.powerlawTempRef, config.Flow.powerlawViscRef,
                                     config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef,
                                     config.Flow.viscosityModel,
                                     Grid.xBnum, config.Grid.xNum,
                                     Grid.yBnum, config.Grid.yNum,
                                     Grid.zBnum, config.Grid.zNum)
      Flow_UpdateGhostConservedStep2(Fluid,
                                     config,
                                     Grid.xBnum, config.Grid.xNum,
                                     Grid.yBnum, config.Grid.yNum,
                                     Grid.zBnum, config.Grid.zNum)

      -- Handle particle collisions
      -- TODO: Collisions across tiles are not handled.
      if config.Particles.collisions and Integrator_stage==4 then
        for c in tiles do
          Particles_HandleCollisions(p_Particles[c], Integrator_deltaTime, config.Particles.restitutionCoeff)
        end
      end

      -- Handle particle boundary conditions
      Particles_UpdateAuxiliaryStep1(Particles,
                                     BC.xBCParticles,
                                     BC.yBCParticles,
                                     BC.zBCParticles,
                                     config.Grid.origin[0], config.Grid.xWidth,
                                     config.Grid.origin[1], config.Grid.yWidth,
                                     config.Grid.origin[2], config.Grid.zWidth,
                                     config.Particles.restitutionCoeff)
      Particles_UpdateAuxiliaryStep2(Particles)
      for c in tiles do
        Particles_number +=
          Particles_DeleteEscapingParticles(p_Particles[c],
                                            config.Grid.origin[0], config.Grid.xWidth,
                                            config.Grid.origin[1], config.Grid.yWidth,
                                            config.Grid.origin[2], config.Grid.zWidth)
      end

      -- Move particles to new partitions
      for c in tiles do
        Particles_LocateInCells(p_Particles[c],
                                Grid.xBnum, config.Grid.xNum, config.Grid.origin[0], config.Grid.xWidth,
                                Grid.yBnum, config.Grid.yNum, config.Grid.origin[1], config.Grid.yWidth,
                                Grid.zBnum, config.Grid.zNum, config.Grid.origin[2], config.Grid.zWidth)
      end
      if numTiles > 1 then
        for c in tiles do
          TradeQueue_clearSource([UTIL.range(1,26):map(function(k) return rexpr [p_TradeQueue[k]][c] end end)])
        end
        for c in tiles do
          TradeQueue_fillSource(c,
                                p_Particles[c],
                                [UTIL.range(1,26):map(function(k) return rexpr [p_TradeQueue[k]][c] end end)],
                                Grid.xBnum, config.Grid.xNum, NX,
                                Grid.yBnum, config.Grid.yNum, NY,
                                Grid.zBnum, config.Grid.zNum, NZ)
        end
        for c in tiles do
          TradeQueue_push(p_Particles[c],
                          [UTIL.range(1,26):map(function(k) return rexpr [p_TradeQueue[k]][c] end end)])
        end
        for c in tiles do
          TradeQueue_fillTarget(p_Particles[c],
                                [UTIL.range(1,26):map(function(k) return rexpr
                                   [p_TradeQueue[k]][ (c-[colorOffsets[k]]+{NX,NY,NZ}) % {NX,NY,NZ} ]
                                 end end)])
        end
        for c in tiles do
          TradeQueue_pull(p_Particles[c],
                          [UTIL.range(1,26):map(function(k) return rexpr
                             [p_TradeQueue[k]][ (c-[colorOffsets[k]]+{NX,NY,NZ}) % {NX,NY,NZ} ]
                           end end)])
        end
      end

      Integrator_simTime = (Integrator_time_old+((0.5*(1+(Integrator_stage/3)))*Integrator_deltaTime))

    end -- RK4 sub-time-stepping

    -- update time derivatives at boundary for NSCBC
    if ((config.BC.xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow) and (config.BC.xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow)) then
      Flow_UpdateNSCBCGhostCellTimeDerivatives(Fluid,
                                               config,
                                               Grid.xBnum, config.Grid.xNum,
                                               Grid.yBnum, config.Grid.yNum,
                                               Grid.zBnum, config.Grid.zNum,
                                               Integrator_deltaTime)
    end

    Integrator_timeStep += 1

  end end -- MainLoopBody

  -----------------------------------------------------------------------------
  -- Cleanup code
  -----------------------------------------------------------------------------

  function INSTANCE.Cleanup(config) return rquote

    -- Wait for everything above to finish
    __fence(__execution, __block)

    -- Report final time
    var endTime = regentlib.c.legion_get_current_time_in_micros() / 1000;
    [emitConsoleWrite(config,
                      'Total time: %llu.%03llu seconds\n',
                      rexpr (endTime - startTime) / 1000 end,
                      rexpr (endTime - startTime) % 1000 end)];

  end end -- Cleanup

return INSTANCE end -- mkInstance

-------------------------------------------------------------------------------
-- TOP-LEVEL INTERFACE
-------------------------------------------------------------------------------

local CopyQueue = regentlib.newsymbol()
local FakeCopyQueue = regentlib.newsymbol()

local function parallelizeFor(sim, stmts)
  return rquote
    __parallelize_with
      sim.p_Fluid, sim.p_Particles, sim.p_Radiation, sim.tiles,
      image(sim.Fluid, sim.p_Particles, [sim.Particles].cell) <= sim.p_Fluid
    do [stmts] end
  end
end

local SIM = mkInstance()

__forbid(__optimize) __demand(__inner)
task workSingle(config : Config)
  [SIM.DeclSymbols(config)];
  var is_FakeCopyQueue = ispace(int1d, 0)
  var [FakeCopyQueue] = region(is_FakeCopyQueue, CopyQueue_columns);
  [UTIL.emitRegionTagAttach(FakeCopyQueue, MAPPER.SAMPLE_ID_TAG, -1, int)];
  [parallelizeFor(SIM, rquote
    [SIM.InitRegions(config)];
    while true do
      [SIM.MainLoopHeader(config)];
      [SIM.PerformIO(config)];
      if SIM.Integrator_exitCond then
        break
      end
      [SIM.MainLoopBody(config, FakeCopyQueue)];
    end
  end)];
  [SIM.Cleanup(config)];
end

local SIM0 = mkInstance()
local SIM1 = mkInstance()

__forbid(__optimize) __demand(__inner)
task workDual(mc : MultiConfig)
  -- Declare symbols
  [SIM0.DeclSymbols(rexpr mc.configs[0] end)];
  [SIM1.DeclSymbols(rexpr mc.configs[1] end)];
  var is_FakeCopyQueue = ispace(int1d, 0)
  var [FakeCopyQueue] = region(is_FakeCopyQueue, CopyQueue_columns);
  [UTIL.emitRegionTagAttach(FakeCopyQueue, MAPPER.SAMPLE_ID_TAG, -1, int)];
  var copySrcOrigin = array(
    SIM0.Grid.xRealOrigin + mc.copySrc.fromCell[0] * SIM0.Grid.xCellWidth,
    SIM0.Grid.yRealOrigin + mc.copySrc.fromCell[1] * SIM0.Grid.yCellWidth,
    SIM0.Grid.zRealOrigin + mc.copySrc.fromCell[2] * SIM0.Grid.zCellWidth)
  var copyTgtOrigin = array(
    SIM1.Grid.xRealOrigin + mc.copyTgt.fromCell[0] * SIM1.Grid.xCellWidth,
    SIM1.Grid.yRealOrigin + mc.copyTgt.fromCell[1] * SIM1.Grid.yCellWidth,
    SIM1.Grid.zRealOrigin + mc.copyTgt.fromCell[2] * SIM1.Grid.zCellWidth)
  var Fluid0_cellWidth = array(SIM0.Grid.xCellWidth, SIM0.Grid.yCellWidth, SIM0.Grid.zCellWidth)
  var Fluid1_cellWidth = array(SIM1.Grid.xCellWidth, SIM1.Grid.yCellWidth, SIM1.Grid.zCellWidth)
  var CopyQueue_ptr : int64 = 0
  var coloring_CopyQueue = regentlib.c.legion_domain_point_coloring_create()
  for c in SIM0.tiles do
    var partSize = CopyQueue_partSize(SIM0.p_Fluid[c].bounds,
                                      mc.configs[0],
                                      mc.copySrc)
    regentlib.c.legion_domain_point_coloring_color_domain(
      coloring_CopyQueue, c, rect1d{CopyQueue_ptr,CopyQueue_ptr+partSize-1})
    CopyQueue_ptr += partSize
  end
  var is_CopyQueue = ispace(int1d, CopyQueue_ptr)
  var [CopyQueue] = region(is_CopyQueue, CopyQueue_columns);
  [UTIL.emitRegionTagAttach(CopyQueue, MAPPER.SAMPLE_ID_TAG, rexpr mc.configs[0].Mapping.sampleId end, int)];
  var p_CopyQueue = partition(disjoint, CopyQueue, coloring_CopyQueue, SIM0.tiles)
  regentlib.c.legion_domain_point_coloring_destroy(coloring_CopyQueue)
  -- Check 2-section configuration
  regentlib.assert(
    -- copySrc is a valid volume
    0 <= mc.copySrc.fromCell[0] and
    0 <= mc.copySrc.fromCell[1] and
    0 <= mc.copySrc.fromCell[2] and
    mc.copySrc.fromCell[0] <= mc.copySrc.uptoCell[0] and
    mc.copySrc.fromCell[1] <= mc.copySrc.uptoCell[1] and
    mc.copySrc.fromCell[2] <= mc.copySrc.uptoCell[2] and
    mc.copySrc.uptoCell[0] < mc.configs[0].Grid.xNum + 2 * SIM0.Grid.xBnum and
    mc.copySrc.uptoCell[1] < mc.configs[0].Grid.yNum + 2 * SIM0.Grid.yBnum and
    mc.copySrc.uptoCell[2] < mc.configs[0].Grid.zNum + 2 * SIM0.Grid.zBnum and
    -- copyTgt is a valid volume
    0 <= mc.copyTgt.fromCell[0] and
    0 <= mc.copyTgt.fromCell[1] and
    0 <= mc.copyTgt.fromCell[2] and
    mc.copyTgt.fromCell[0] <= mc.copyTgt.uptoCell[0] and
    mc.copyTgt.fromCell[1] <= mc.copyTgt.uptoCell[1] and
    mc.copyTgt.fromCell[2] <= mc.copyTgt.uptoCell[2] and
    mc.copyTgt.uptoCell[0] < mc.configs[1].Grid.xNum + 2 * SIM1.Grid.xBnum and
    mc.copyTgt.uptoCell[1] < mc.configs[1].Grid.yNum + 2 * SIM1.Grid.yBnum and
    mc.copySrc.uptoCell[2] < mc.configs[1].Grid.zNum + 2 * SIM1.Grid.zBnum and
    -- volumes have the same size
    mc.copySrc.uptoCell[0] - mc.copySrc.fromCell[0] ==
    mc.copyTgt.uptoCell[0] - mc.copyTgt.fromCell[0] and
    mc.copySrc.uptoCell[1] - mc.copySrc.fromCell[1] ==
    mc.copyTgt.uptoCell[1] - mc.copyTgt.fromCell[1] and
    mc.copySrc.uptoCell[2] - mc.copySrc.fromCell[2] ==
    mc.copyTgt.uptoCell[2] - mc.copyTgt.fromCell[2],
    'Invalid volume copy configuration');
  -- Initialize regions & partitions
  [parallelizeFor(SIM0, SIM0.InitRegions(rexpr mc.configs[0] end))];
  [parallelizeFor(SIM1, SIM1.InitRegions(rexpr mc.configs[1] end))];
  var srcOrigin = int3d{mc.copySrc.fromCell[0], mc.copySrc.fromCell[1], mc.copySrc.fromCell[2]}
  var tgtOrigin = int3d{mc.copyTgt.fromCell[0], mc.copyTgt.fromCell[1], mc.copyTgt.fromCell[2]}
  var srcColoring = regentlib.c.legion_domain_point_coloring_create()
  for c in SIM1.tiles do
    var tgtRect = intersection(SIM1.p_Fluid[c].bounds, mc.copyTgt)
    if rectSize(tgtRect) > 0 then
      var srcRect = rect3d{lo = tgtRect.lo - tgtOrigin + srcOrigin,
                           hi = tgtRect.hi - tgtOrigin + srcOrigin}
      regentlib.c.legion_domain_point_coloring_color_domain(srcColoring, c, srcRect)
    end
  end
  var p_Fluid0_src = partition(disjoint, SIM0.Fluid, srcColoring, SIM1.tiles)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring)
  var tgtColoring = regentlib.c.legion_domain_point_coloring_create()
  regentlib.c.legion_domain_point_coloring_color_domain(tgtColoring, int1d(0),
    rect3d{lo = int3d{mc.copyTgt.fromCell[0], mc.copyTgt.fromCell[1], mc.copyTgt.fromCell[2]},
           hi = int3d{mc.copyTgt.uptoCell[0], mc.copyTgt.uptoCell[1], mc.copyTgt.uptoCell[2]}})
  var p_Fluid1_isCopied = partition(disjoint, SIM1.Fluid, tgtColoring, ispace(int1d,1))
  regentlib.c.legion_domain_point_coloring_destroy(tgtColoring)
  var p_Fluid1_tgt = cross_product(SIM1.p_Fluid, p_Fluid1_isCopied)
  -- Main simulation loop
  while true do
    -- Perform preliminary actions before each timestep
    [parallelizeFor(SIM0, SIM0.MainLoopHeader(rexpr mc.configs[0] end))];
    [parallelizeFor(SIM1, SIM1.MainLoopHeader(rexpr mc.configs[1] end))];
    -- Make sure both simulations are using the same timestep
    SIM0.Integrator_deltaTime = min(SIM0.Integrator_deltaTime, SIM1.Integrator_deltaTime)
    SIM1.Integrator_deltaTime = min(SIM0.Integrator_deltaTime, SIM1.Integrator_deltaTime);
    [parallelizeFor(SIM0, SIM0.PerformIO(rexpr mc.configs[0] end))];
    [parallelizeFor(SIM1, SIM1.PerformIO(rexpr mc.configs[1] end))];
    if SIM0.Integrator_exitCond or SIM1.Integrator_exitCond then
      break
    end
    -- Run 1 iteration of first section
    [parallelizeFor(SIM0, SIM0.MainLoopBody(rexpr mc.configs[0] end, FakeCopyQueue))];
    -- Copy fluid values to second section
    for c in SIM1.tiles do
      var src = p_Fluid0_src[c]
      var tgt = p_Fluid1_tgt[c][0]
      copy(src.temperature, tgt.temperature_inc)
      copy(src.velocity, tgt.velocity_inc)
    end
    -- Copy particles to second section
    if mc.configs[1].Particles.feeding.type == SCHEMA.FeedModel_Incoming then
      fill(CopyQueue.position, array(-1.0, -1.0, -1.0))
      fill(CopyQueue.velocity, array(-1.0, -1.0, -1.0))
      fill(CopyQueue.temperature, -1.0)
      fill(CopyQueue.diameter, -1.0)
      fill(CopyQueue.density, -1.0)
      fill(CopyQueue.__valid, false)
      for c in SIM0.tiles do
        CopyQueue_push(SIM0.p_Particles[c],
                       p_CopyQueue[c],
                       mc.copySrc,
                       copySrcOrigin, copyTgtOrigin,
                       Fluid0_cellWidth, Fluid1_cellWidth)
      end
    end
    -- Run 1 iteration of second section
    [parallelizeFor(SIM1, SIM1.MainLoopBody(rexpr mc.configs[1] end, CopyQueue))];
  end
  -- Cleanups
  [SIM0.Cleanup(rexpr mc.configs[0] end)];
  [SIM1.Cleanup(rexpr mc.configs[1] end)];
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
    elseif C.strcmp(args.argv[i], '-m') == 0 and i < args.argc-1 then
      var mc : MultiConfig[1]
      SCHEMA.parse_MultiConfig([&MultiConfig](mc), args.argv[i+1])
      initSample([&Config](mc[0].configs), launched, outDirBase)
      initSample([&Config](mc[0].configs) + 1, launched + 1, outDirBase)
      launched += 2
      workDual(mc[0])
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

regentlib.saveobj(main, "soleil.o", "object", MAPPER.register_mappers)
