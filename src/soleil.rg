import "regent"

-------------------------------------------------------------------------------
-- IMPORTS
-------------------------------------------------------------------------------

local C = terralib.includecstring[[
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
]]
local MAPPER = terralib.includec("soleil_mapper.h")
local SCHEMA = terralib.includec("config_schema.h")

local acos = regentlib.acos(double)
local ceil = regentlib.ceil(double)
local cos = regentlib.cos(double)
local fabs = regentlib.fabs(double)
local fmod = regentlib.fmod(double)
local pow = regentlib.pow(double)
local sin = regentlib.sin(double)
local sqrt = regentlib.sqrt(double)

-------------------------------------------------------------------------------
-- COMPILE-TIME CONFIGURATION
-------------------------------------------------------------------------------

local USE_HDF = assert(os.getenv('USE_HDF')) ~= '0'

local NUM_ANGLES = 14

-------------------------------------------------------------------------------
-- DATA STRUCTURES
-------------------------------------------------------------------------------

local Config = SCHEMA.Config

local struct particles_columns {
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
  kineticEnergy : double;
  sgsEnergy : double;
  sgsEddyViscosity : double;
  sgsEddyKappa : double;
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
}

struct Radiation_columns {
  I_1 : double[NUM_ANGLES];
  I_2 : double[NUM_ANGLES];
  I_3 : double[NUM_ANGLES];
  I_4 : double[NUM_ANGLES];
  I_5 : double[NUM_ANGLES];
  I_6 : double[NUM_ANGLES];
  I_7 : double[NUM_ANGLES];
  I_8 : double[NUM_ANGLES];
  Iiter_1 : double[NUM_ANGLES];
  Iiter_2 : double[NUM_ANGLES];
  Iiter_3 : double[NUM_ANGLES];
  Iiter_4 : double[NUM_ANGLES];
  Iiter_5 : double[NUM_ANGLES];
  Iiter_6 : double[NUM_ANGLES];
  Iiter_7 : double[NUM_ANGLES];
  Iiter_8 : double[NUM_ANGLES];
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

local DOM = (require 'dom/dom')(NUM_ANGLES, Radiation_columns)

-------------------------------------------------------------------------------
-- CONSTANTS
-------------------------------------------------------------------------------

local PI = 3.1415926535898

local FILE_PTR_SIZE = sizeof(&C.FILE)

-------------------------------------------------------------------------------
-- MACROS
-------------------------------------------------------------------------------

-- TODO: Define macros for in_boundary, neg_depth etc.

-------------------------------------------------------------------------------
-- I/O ROUTINES
-------------------------------------------------------------------------------

local task Fluid_dump(colors : ispace(int3d),
                      filename : rawstring,
                      r : region(ispace(int3d),Fluid_columns),
                      s : region(ispace(int3d),Fluid_columns),
                      p_r : partition(disjoint, r, colors),
                      p_s : partition(disjoint, s, colors))
  regentlib.assert(false, 'Recompile with USE_HDF=1')
end
local task Fluid_load(colors : ispace(int3d),
                      filename : rawstring,
                      r : region(ispace(int3d),Fluid_columns),
                      s : region(ispace(int3d),Fluid_columns),
                      p_r : partition(disjoint, r, colors),
                      p_s : partition(disjoint, s, colors))
  regentlib.assert(false, 'Recompile with USE_HDF=1')
end
local task particles_dump(colors : ispace(int3d),
                          filename : rawstring,
                          r : region(ispace(int1d),particles_columns),
                          s : region(ispace(int1d),particles_columns),
                          p_r : partition(disjoint, r, colors),
                          p_s : partition(disjoint, s, colors))
  regentlib.assert(false, 'Recompile with USE_HDF=1')
end
local task particles_load(colors : ispace(int3d),
                          filename : rawstring,
                          r : region(ispace(int1d),particles_columns),
                          s : region(ispace(int1d),particles_columns),
                          p_r : partition(disjoint, r, colors),
                          p_s : partition(disjoint, s, colors))
  regentlib.assert(false, 'Recompile with USE_HDF=1')
end

if USE_HDF then
  local HDF = require "hdf_helper"
  Fluid_dump, Fluid_load = HDF.mkHDFTasks(
    int3d, int3d, Fluid_columns,
    {"rho","pressure","velocity"})
  particles_dump, particles_load = HDF.mkHDFTasks(
    int1d, int3d, particles_columns,
    {"cell","position","velocity","temperature","diameter","__valid"})
end

-------------------------------------------------------------------------------
-- OTHER ROUTINES
-------------------------------------------------------------------------------

__demand(__parallel)
task InitParticlesUniform(particles : region(ispace(int1d), particles_columns),
                          cells : region(ispace(int3d), Fluid_columns),
                          config : Config,
                          xBnum : int32, yBnum : int32, zBnum : int32)
where
  reads(cells.{centerCoordinates, velocity}),
  reads writes(particles)
do
  var pBase = 0
  for p in particles do
    pBase = int32(p)
    break
  end
  var lo = cells.bounds.lo
  lo.x = max(lo.x, xBnum)
  lo.y = max(lo.y, yBnum)
  lo.z = max(lo.z, zBnum)
  var hi = cells.bounds.hi
  hi.x = min(hi.x, ((config.Grid.xNum+xBnum)-1))
  hi.y = min(hi.y, ((config.Grid.yNum+yBnum)-1))
  hi.z = min(hi.z, ((config.Grid.zNum+zBnum)-1))
  var xSize = ((hi.x-lo.x)+1)
  var ySize = ((hi.y-lo.y)+1)
  var particlesPerTask = (config.Particles.initNum/((config.Mapping.xTiles*config.Mapping.yTiles)*config.Mapping.zTiles))
  __demand(__openmp)
  for p in particles do
    if ((int32(p)-pBase)<particlesPerTask) then
      p.__valid = true
      var relIdx = (int32(p)-pBase)
      var c = int3d({(lo.x+(relIdx%xSize)), (lo.y+((relIdx/xSize)%ySize)), (lo.z+((relIdx/xSize)/ySize))})
      p.cell = c
      p.position = cells[p.cell].centerCoordinates
      p.velocity = cells[p.cell].velocity
      p.density = config.Particles.density
      p.temperature = config.Particles.initTemperature
      p.diameter = config.Particles.diameterMean
    end
  end
end

__demand(__parallel, __cuda)
task AddRadiation(particles : region(ispace(int1d), particles_columns),
                  config : Config)
where
  reads(particles.{density, diameter}),
  reads writes(particles.temperature_t)
do
  var absorptivity = config.Particles.absorptivity
  var intensity = config.Radiation.intensity
  var heatCapacity = config.Particles.heatCapacity
  __demand(__openmp)
  for p in particles do
    var crossSectionArea = (((2*acos(0))*pow(p.diameter, 2))/4)
    var volume = (((2*acos(0))*pow(p.diameter, 3))/6)
    var mass = (volume*p.density)
    var absorbedRadiationIntensity = ((absorptivity*intensity)*crossSectionArea)
    p.temperature_t += (absorbedRadiationIntensity/(mass*heatCapacity))
  end
end

__demand(__parallel, __cuda)
task particles_initValidField(r : region(ispace(int1d), particles_columns))
where
  writes(r.__valid)
do
  __demand(__openmp)
  for e in r do
    e.__valid = false
  end
end

__demand(__parallel, __cuda)
task SetCoarseningField(Fluid : region(ispace(int3d), Fluid_columns),
                        Grid_xBnum : int32, Grid_xNum : int32,
                        Grid_yBnum : int32, Grid_yNum : int32,
                        Grid_zBnum : int32, Grid_zNum : int32,
                        Radiation_xBnum : int32, Radiation_xNum : int32,
                        Radiation_yBnum : int32, Radiation_yNum : int32,
                        Radiation_zBnum : int32, Radiation_zNum : int32)
where
  reads writes(Fluid.to_Radiation)
do
  __demand(__openmp)
  for f in Fluid do
    var xFactor = (Grid_xNum/Radiation_xNum)
    var yFactor = (Grid_yNum/Radiation_yNum)
    var zFactor = (Radiation_zNum/Radiation_zNum)
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(f).x)), 0)>0) or (max(int32((int3d(f).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(f).y)), 0)>0)) or (max(int32((int3d(f).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(f).z)), 0)>0)) or (max(int32((int3d(f).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      Fluid[f].to_Radiation = int3d({(((int3d(f).x-uint64(Grid_xBnum))/uint64(xFactor))+uint64(Radiation_xBnum)), (((int3d(f).y-uint64(Grid_yBnum))/uint64(yFactor))+uint64(Radiation_yBnum)), (((int3d(f).z-uint64(Grid_zBnum))/uint64(zFactor))+uint64(Radiation_zBnum))})
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
  writes(Fluid.kineticEnergy),
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
  writes(Fluid.sgsEddyKappa),
  writes(Fluid.sgsEddyViscosity),
  writes(Fluid.sgsEnergy),
  writes(Fluid.temperature),
  writes(Fluid.temperatureBoundary),
  writes(Fluid.velocity),
  writes(Fluid.velocityBoundary),
  writes(Fluid.{velocityGradientX, velocityGradientY, velocityGradientZ}),
  writes(Fluid.{velocityGradientXBoundary, velocityGradientYBoundary, velocityGradientZBoundary}),
  writes(Fluid.viscousSpectralRadius)
do
  __demand(__openmp)
  for c in Fluid do
    Fluid[c].rho = 0.0
    Fluid[c].pressure = 0.0
    Fluid[c].velocity = [double[3]](array(0.0, 0.0, 0.0))
    Fluid[c].centerCoordinates = [double[3]](array(0.0, 0.0, 0.0))
    Fluid[c].velocityGradientX = [double[3]](array(0.0, 0.0, 0.0))
    Fluid[c].velocityGradientY = [double[3]](array(0.0, 0.0, 0.0))
    Fluid[c].velocityGradientZ = [double[3]](array(0.0, 0.0, 0.0))
    Fluid[c].temperature = 0.0
    Fluid[c].rhoEnthalpy = 0.0
    Fluid[c].kineticEnergy = 0.0
    Fluid[c].sgsEnergy = 0.0
    Fluid[c].sgsEddyViscosity = 0.0
    Fluid[c].sgsEddyKappa = 0.0
    Fluid[c].convectiveSpectralRadius = 0.0
    Fluid[c].viscousSpectralRadius = 0.0
    Fluid[c].heatConductionSpectralRadius = 0.0
    Fluid[c].rhoVelocity = [double[3]](array(0.0, 0.0, 0.0))
    Fluid[c].rhoEnergy = 0.0
    Fluid[c].rhoBoundary = 0.0
    Fluid[c].rhoVelocityBoundary = [double[3]](array(0.0, 0.0, 0.0))
    Fluid[c].rhoEnergyBoundary = 0.0
    Fluid[c].velocityBoundary = [double[3]](array(0.0, 0.0, 0.0))
    Fluid[c].pressureBoundary = 0.0
    Fluid[c].temperatureBoundary = 0.0
    Fluid[c].velocityGradientXBoundary = [double[3]](array(0.0, 0.0, 0.0))
    Fluid[c].velocityGradientYBoundary = [double[3]](array(0.0, 0.0, 0.0))
    Fluid[c].velocityGradientZBoundary = [double[3]](array(0.0, 0.0, 0.0))
    Fluid[c].rho_old = 0.0
    Fluid[c].rhoVelocity_old = [double[3]](array(0.0, 0.0, 0.0))
    Fluid[c].rhoEnergy_old = 0.0
    Fluid[c].rho_new = 0.0
    Fluid[c].rhoVelocity_new = [double[3]](array(0.0, 0.0, 0.0))
    Fluid[c].rhoEnergy_new = 0.0
    Fluid[c].rho_t = 0.0
    Fluid[c].rhoVelocity_t = [double[3]](array(0.0, 0.0, 0.0))
    Fluid[c].rhoEnergy_t = 0.0
    Fluid[c].rhoFluxX = 0.0
    Fluid[c].rhoVelocityFluxX = [double[3]](array(0.0, 0.0, 0.0))
    Fluid[c].rhoEnergyFluxX = 0.0
    Fluid[c].rhoFluxY = 0.0
    Fluid[c].rhoVelocityFluxY = [double[3]](array(0.0, 0.0, 0.0))
    Fluid[c].rhoEnergyFluxY = 0.0
    Fluid[c].rhoFluxZ = 0.0
    Fluid[c].rhoVelocityFluxZ = [double[3]](array(0.0, 0.0, 0.0))
    Fluid[c].rhoEnergyFluxZ = 0.0
    Fluid[c].PD = 0.0
    Fluid[c].dissipation = 0.0
    Fluid[c].dissipationFlux = 0.0
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
    var xy = [double[3]](array((Grid_xOrigin+((Grid_xWidth/double(Grid_xNum))*(double((int3d(c).x-uint64(Grid_xBnum)))+double(0.5)))), (Grid_yOrigin+((Grid_yWidth/double(Grid_yNum))*(double((int3d(c).y-uint64(Grid_yBnum)))+double(0.5)))), (Grid_zOrigin+((Grid_zWidth/double(Grid_zNum))*(double((int3d(c).z-uint64(Grid_zBnum)))+double(0.5))))))
    Fluid[c].centerCoordinates = [double[3]](array(double(xy[0]), double(xy[1]), double(xy[2])))
  end
end

__demand(__parallel, __cuda)
task Flow_InitializeUniform(Fluid : region(ispace(int3d), Fluid_columns), Flow_initParams : double[5])
where
  writes(Fluid.{rho, pressure}),
  reads writes(Fluid.velocity)
do
  __demand(__openmp)
  for c in Fluid do
    Fluid[c].rho = Flow_initParams[0]
    Fluid[c].pressure = Flow_initParams[1]
    Fluid[c].velocity[0] = Flow_initParams[2]
    Fluid[c].velocity[1] = Flow_initParams[3]
    Fluid[c].velocity[2] = Flow_initParams[4]
  end
end

terra vs_mul_double_3(a : double[3],b : double) : double[3]
  return array([&double](a)[0] * b, [&double](a)[1] * b, [&double](a)[2] * b)
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
    var xy = [double[3]](array((Grid_xOrigin+((Grid_xWidth/double(Grid_xNum))*(double((int3d(c).x-uint64(Grid_xBnum)))+double(0.5)))), (Grid_yOrigin+((Grid_yWidth/double(Grid_yNum))*(double((int3d(c).y-uint64(Grid_yBnum)))+double(0.5)))), (Grid_zOrigin+((Grid_zWidth/double(Grid_zNum))*(double((int3d(c).z-uint64(Grid_zBnum)))+double(0.5))))))
    var coorZ = 0
    Fluid[c].rho = taylorGreenDensity
    Fluid[c].velocity = vs_mul_double_3([double[3]](array(((sin(xy[0])*cos(xy[1]))*cos(coorZ)), (((-cos(xy[0]))*sin(xy[1]))*cos(coorZ)), 0.0)), taylorGreenVelocity)
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
    var xy = [double[3]](array((Grid_xOrigin+((Grid_xWidth/double(Grid_xNum))*(double((int3d(c).x-uint64(Grid_xBnum)))+double(0.5)))), (Grid_yOrigin+((Grid_yWidth/double(Grid_yNum))*(double((int3d(c).y-uint64(Grid_yBnum)))+double(0.5)))), (Grid_zOrigin+((Grid_zWidth/double(Grid_zNum))*(double((int3d(c).z-uint64(Grid_zBnum)))+double(0.5))))))
    Fluid[c].rho = taylorGreenDensity
    Fluid[c].velocity = vs_mul_double_3([double[3]](array(((sin(xy[0])*cos(xy[1]))*cos(xy[2])), (((-cos(xy[0]))*sin(xy[1]))*cos(xy[2])), 0.0)), taylorGreenVelocity)
    var factorA = (cos((2.0*xy[2]))+2.0)
    var factorB = (cos((2.0*xy[0]))+cos((2.0*xy[1])))
    Fluid[c].pressure = (taylorGreenPressure+((((taylorGreenDensity*pow(taylorGreenVelocity, 2.0))/16.0)*factorA)*factorB))
  end
end

__demand(__parallel)
task Flow_InitializePerturbed(Fluid : region(ispace(int3d), Fluid_columns),
                              Flow_initParams : double[5])
where
  writes(Fluid.{rho, pressure}),
  reads writes(Fluid.velocity)
do
  for c in Fluid do
    Fluid[c].rho = Flow_initParams[0]
    Fluid[c].pressure = Flow_initParams[1]
    Fluid[c].velocity[0] = (Flow_initParams[2]+(((double(C.rand())/2147483647)-double(0.5))*10.0))
    Fluid[c].velocity[1] = (Flow_initParams[3]+(((double(C.rand())/2147483647)-double(0.5))*10.0))
    Fluid[c].velocity[2] = (Flow_initParams[4]+(((double(C.rand())/2147483647)-double(0.5))*10.0))
  end
end

terra dot_double_3(a : double[3],b : double[3]) : double
  return [&double](a)[0] * [&double](b)[0] + [&double](a)[1] * [&double](b)[1] + [&double](a)[2] * [&double](b)[2]
end

__demand(__parallel, __cuda)
task Flow_UpdateConservedFromPrimitive(Fluid : region(ispace(int3d), Fluid_columns),
                                       Flow_gamma : double,
                                       Flow_gasConstant : double,
                                       Grid_xBnum : int32, Grid_xNum : int32, Grid_yBnum : int32, Grid_yNum : int32, Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{rho, velocity, pressure, sgsEnergy}),
  writes(Fluid.{rhoVelocity, rhoEnergy})
do
  __demand(__openmp)
  for c in Fluid do
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      var tmpTemperature = (Fluid[c].pressure/(Flow_gasConstant*Fluid[c].rho))
      var velocity = Fluid[c].velocity
      Fluid[c].rhoVelocity = vs_mul_double_3(Fluid[c].velocity, Fluid[c].rho)
      var cv = (Flow_gasConstant/(Flow_gamma-1.0))
      Fluid[c].rhoEnergy = ((Fluid[c].rho*((cv*tmpTemperature)+(double(0.5)*dot_double_3(velocity, velocity))))+Fluid[c].sgsEnergy)
    end
  end
end

terra vs_div_double_3(a : double[3],b : double) : double[3]
  return array([&double](a)[0] / b, [&double](a)[1] / b, [&double](a)[2] / b)
end

__demand(__parallel, __cuda)
task Flow_UpdateAuxiliaryVelocity(Fluid : region(ispace(int3d), Fluid_columns),
                                  Grid_xBnum : int32, Grid_xNum : int32,
                                  Grid_yBnum : int32, Grid_yNum : int32,
                                  Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{rho, rhoVelocity}),
  writes(Fluid.{velocity, kineticEnergy})
do
  __demand(__openmp)
  for c in Fluid do
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      var velocity = vs_div_double_3(Fluid[c].rhoVelocity, Fluid[c].rho)
      Fluid[c].velocity = velocity
      Fluid[c].kineticEnergy = ((double(0.5)*Fluid[c].rho)*dot_double_3(velocity, velocity))
    end
  end
end

terra vv_mul_double_3(a : double[3],b : double[3]) : double[3]
  return array([&double](a)[0] * [&double](b)[0], [&double](a)[1] * [&double](b)[1], [&double](a)[2] * [&double](b)[2])
end

terra vv_add_double_3(a : double[3],b : double[3]) : double[3]
  return array([&double](a)[0] + [&double](b)[0], [&double](a)[1] + [&double](b)[1], [&double](a)[2] + [&double](b)[2])
end

__demand(__parallel, __cuda)
task Flow_UpdateGhostConservedStep1(Fluid : region(ispace(int3d), Fluid_columns),
                                    BC_xNegTemperature : double, BC_xNegVelocity : double[3], BC_xPosTemperature : double, BC_xPosVelocity : double[3], BC_xSign : double[3],
                                    BC_yNegTemperature : double, BC_yNegVelocity : double[3], BC_yPosTemperature : double, BC_yPosVelocity : double[3], BC_ySign : double[3],
                                    BC_zNegTemperature : double, BC_zNegVelocity : double[3], BC_zPosTemperature : double, BC_zPosVelocity : double[3], BC_zSign : double[3],
                                    Flow_gamma : double,
                                    Flow_gasConstant : double,
                                    Grid_xBnum : int32, Grid_xNum : int32,
                                    Grid_yBnum : int32, Grid_yNum : int32,
                                    Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{rho, pressure, rhoVelocity, temperature}),
  reads writes(Fluid.{rhoBoundary, rhoEnergyBoundary, rhoVelocityBoundary})
do
  __demand(__openmp)
  for c in Fluid do
    if (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{1, 0, 0})%Fluid.bounds)
      var sign = BC_xSign
      var bnd_velocity = BC_xNegVelocity
      var bnd_temperature = BC_xNegTemperature
      var rho = double(0.0)
      var temp_wall = double(0.0)
      var temperature = double(0.0)
      var velocity = [double[3]](array(0.0, 0.0, 0.0))
      var cv = (Flow_gasConstant/(Flow_gamma-1.0))
      var velocity__3505 = [double[3]](array(0.0, 0.0, 0.0))
      velocity__3505 = vv_add_double_3(vv_mul_double_3(vs_div_double_3(Fluid[c_int].rhoVelocity, Fluid[c_int].rho), sign), bnd_velocity)
      temp_wall = Fluid[c_int].temperature
      if (bnd_temperature>0.0) then
        temp_wall = bnd_temperature
      end
      temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
      rho = (Fluid[c_int].pressure/(Flow_gasConstant*temperature))
      Fluid[c_bnd].rhoBoundary = rho
      Fluid[c_bnd].rhoVelocityBoundary = vs_mul_double_3(velocity__3505, rho)
      Fluid[c_bnd].rhoEnergyBoundary = (rho*((cv*temperature)+(double(0.5)*dot_double_3(velocity__3505, velocity__3505))))
    end
    if (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{-1, 0, 0})%Fluid.bounds)
      var sign = BC_xSign
      var bnd_velocity = BC_xPosVelocity
      var bnd_temperature = BC_xPosTemperature
      var rho = double(0.0)
      var temp_wall = double(0.0)
      var temperature = double(0.0)
      var velocity = [double[3]](array(0.0, 0.0, 0.0))
      var cv = (Flow_gasConstant/(Flow_gamma-1.0))
      var velocity__3516 = [double[3]](array(0.0, 0.0, 0.0))
      velocity__3516 = vv_add_double_3(vv_mul_double_3(vs_div_double_3(Fluid[c_int].rhoVelocity, Fluid[c_int].rho), sign), bnd_velocity)
      temp_wall = Fluid[c_int].temperature
      if (bnd_temperature>0.0) then
        temp_wall = bnd_temperature
      end
      temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
      rho = (Fluid[c_int].pressure/(Flow_gasConstant*temperature))
      Fluid[c_bnd].rhoBoundary = rho
      Fluid[c_bnd].rhoVelocityBoundary = vs_mul_double_3(velocity__3516, rho)
      Fluid[c_bnd].rhoEnergyBoundary = (rho*((cv*temperature)+(double(0.5)*dot_double_3(velocity__3516, velocity__3516))))
    end
    if (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, 1, 0})%Fluid.bounds)
      var sign = BC_ySign
      var bnd_velocity = BC_yNegVelocity
      var bnd_temperature = BC_yNegTemperature
      var rho = double(0.0)
      var temp_wall = double(0.0)
      var temperature = double(0.0)
      var velocity = [double[3]](array(0.0, 0.0, 0.0))
      var cv = (Flow_gasConstant/(Flow_gamma-1.0))
      var velocity__3527 = [double[3]](array(0.0, 0.0, 0.0))
      velocity__3527 = vv_add_double_3(vv_mul_double_3(vs_div_double_3(Fluid[c_int].rhoVelocity, Fluid[c_int].rho), sign), bnd_velocity)
      temp_wall = Fluid[c_int].temperature
      if (bnd_temperature>0.0) then
        temp_wall = bnd_temperature
      end
      temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
      rho = (Fluid[c_int].pressure/(Flow_gasConstant*temperature))
      Fluid[c_bnd].rhoBoundary = rho
      Fluid[c_bnd].rhoVelocityBoundary = vs_mul_double_3(velocity__3527, rho)
      Fluid[c_bnd].rhoEnergyBoundary = (rho*((cv*temperature)+(double(0.5)*dot_double_3(velocity__3527, velocity__3527))))
    end
    if (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, -1, 0})%Fluid.bounds)
      var sign = BC_ySign
      var bnd_velocity = BC_yPosVelocity
      var bnd_temperature = BC_yPosTemperature
      var rho = double(0.0)
      var temp_wall = double(0.0)
      var temperature = double(0.0)
      var velocity = [double[3]](array(0.0, 0.0, 0.0))
      var cv = (Flow_gasConstant/(Flow_gamma-1.0))
      var velocity__3538 = [double[3]](array(0.0, 0.0, 0.0))
      velocity__3538 = vv_add_double_3(vv_mul_double_3(vs_div_double_3(Fluid[c_int].rhoVelocity, Fluid[c_int].rho), sign), bnd_velocity)
      temp_wall = Fluid[c_int].temperature
      if (bnd_temperature>0.0) then
        temp_wall = bnd_temperature
      end
      temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
      rho = (Fluid[c_int].pressure/(Flow_gasConstant*temperature))
      Fluid[c_bnd].rhoBoundary = rho
      Fluid[c_bnd].rhoVelocityBoundary = vs_mul_double_3(velocity__3538, rho)
      Fluid[c_bnd].rhoEnergyBoundary = (rho*((cv*temperature)+(double(0.5)*dot_double_3(velocity__3538, velocity__3538))))
    end
    if (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, 0, 1})%Fluid.bounds)
      var sign = BC_zSign
      var bnd_velocity = BC_zNegVelocity
      var bnd_temperature = BC_zNegTemperature
      var rho = double(0.0)
      var temp_wall = double(0.0)
      var temperature = double(0.0)
      var velocity = [double[3]](array(0.0, 0.0, 0.0))
      var cv = (Flow_gasConstant/(Flow_gamma-1.0))
      var velocity__3549 = [double[3]](array(0.0, 0.0, 0.0))
      velocity__3549 = vv_add_double_3(vv_mul_double_3(vs_div_double_3(Fluid[c_int].rhoVelocity, Fluid[c_int].rho), sign), bnd_velocity)
      temp_wall = Fluid[c_int].temperature
      if (bnd_temperature>0.0) then
        temp_wall = bnd_temperature
      end
      temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
      rho = (Fluid[c_int].pressure/(Flow_gasConstant*temperature))
      Fluid[c_bnd].rhoBoundary = rho
      Fluid[c_bnd].rhoVelocityBoundary = vs_mul_double_3(velocity__3549, rho)
      Fluid[c_bnd].rhoEnergyBoundary = (rho*((cv*temperature)+(double(0.5)*dot_double_3(velocity__3549, velocity__3549))))
    end
    if (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, 0, -1})%Fluid.bounds)
      var sign = BC_zSign
      var bnd_velocity = BC_zPosVelocity
      var bnd_temperature = BC_zPosTemperature
      var rho = double(0.0)
      var temp_wall = double(0.0)
      var temperature = double(0.0)
      var velocity = [double[3]](array(0.0, 0.0, 0.0))
      var cv = (Flow_gasConstant/(Flow_gamma-1.0))
      var velocity__3560 = [double[3]](array(0.0, 0.0, 0.0))
      velocity__3560 = vv_add_double_3(vv_mul_double_3(vs_div_double_3(Fluid[c_int].rhoVelocity, Fluid[c_int].rho), sign), bnd_velocity)
      temp_wall = Fluid[c_int].temperature
      if (bnd_temperature>0.0) then
        temp_wall = bnd_temperature
      end
      temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
      rho = (Fluid[c_int].pressure/(Flow_gasConstant*temperature))
      Fluid[c_bnd].rhoBoundary = rho
      Fluid[c_bnd].rhoVelocityBoundary = vs_mul_double_3(velocity__3560, rho)
      Fluid[c_bnd].rhoEnergyBoundary = (rho*((cv*temperature)+(double(0.5)*dot_double_3(velocity__3560, velocity__3560))))
    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdateGhostConservedStep2(Fluid : region(ispace(int3d), Fluid_columns),
                                    Grid_xBnum : int32, Grid_xNum : int32,
                                    Grid_yBnum : int32, Grid_yNum : int32,
                                    Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{rhoBoundary, rhoVelocityBoundary, rhoEnergyBoundary}),
  writes(Fluid.{rho, rhoVelocity, rhoEnergy})
do
  __demand(__openmp)
  for c in Fluid do
    if ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)) then
      Fluid[c].rho = Fluid[c].rhoBoundary
      Fluid[c].rhoVelocity = Fluid[c].rhoVelocityBoundary
      Fluid[c].rhoEnergy = Fluid[c].rhoEnergyBoundary
    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdateGhostVelocityStep1(Fluid : region(ispace(int3d), Fluid_columns),
                                   BC_xNegVelocity : double[3], BC_xPosVelocity : double[3], BC_xSign : double[3],
                                   BC_yNegVelocity : double[3], BC_yPosVelocity : double[3], BC_ySign : double[3],
                                   BC_zNegVelocity : double[3], BC_zPosVelocity : double[3], BC_zSign : double[3],
                                   Grid_xBnum : int32, Grid_xNum : int32,
                                   Grid_yBnum : int32, Grid_yNum : int32,
                                   Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.velocity),
  writes(Fluid.velocityBoundary)
do
  __demand(__openmp)
  for c in Fluid do
    if (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{1, 0, 0})%Fluid.bounds)
      var sign = BC_xSign
      var bnd_velocity = BC_xNegVelocity
      Fluid[c_bnd].velocityBoundary = vv_add_double_3(vv_mul_double_3(Fluid[c_int].velocity, sign), bnd_velocity)
    end
    if (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{-1, 0, 0})%Fluid.bounds)
      var sign = BC_xSign
      var bnd_velocity = BC_xPosVelocity
      Fluid[c_bnd].velocityBoundary = vv_add_double_3(vv_mul_double_3(Fluid[c_int].velocity, sign), bnd_velocity)
    end
    if (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, 1, 0})%Fluid.bounds)
      var sign = BC_ySign
      var bnd_velocity = BC_yNegVelocity
      Fluid[c_bnd].velocityBoundary = vv_add_double_3(vv_mul_double_3(Fluid[c_int].velocity, sign), bnd_velocity)
    end
    if (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, -1, 0})%Fluid.bounds)
      var sign = BC_ySign
      var bnd_velocity = BC_yPosVelocity
      Fluid[c_bnd].velocityBoundary = vv_add_double_3(vv_mul_double_3(Fluid[c_int].velocity, sign), bnd_velocity)
    end
    if (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, 0, 1})%Fluid.bounds)
      var sign = BC_zSign
      var bnd_velocity = BC_zNegVelocity
      Fluid[c_bnd].velocityBoundary = vv_add_double_3(vv_mul_double_3(Fluid[c_int].velocity, sign), bnd_velocity)
    end
    if (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, 0, -1})%Fluid.bounds)
      var sign = BC_zSign
      var bnd_velocity = BC_zPosVelocity
      Fluid[c_bnd].velocityBoundary = vv_add_double_3(vv_mul_double_3(Fluid[c_int].velocity, sign), bnd_velocity)
    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdateGhostVelocityStep2(Fluid : region(ispace(int3d), Fluid_columns),
                                   Grid_xBnum : int32, Grid_xNum : int32,
                                   Grid_yBnum : int32, Grid_yNum : int32,
                                   Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.velocityBoundary),
  writes(Fluid.velocity)
do
  __demand(__openmp)
  for c in Fluid do
    if ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)) then
      Fluid[c].velocity = Fluid[c].velocityBoundary
    end
  end
end

terra vv_sub_double_3(a : double[3],b : double[3]) : double[3]
  return array([&double](a)[0] - [&double](b)[0], [&double](a)[1] - [&double](b)[1], [&double](a)[2] - [&double](b)[2])
end

__demand(__parallel, __cuda)
task Flow_ComputeVelocityGradientAll(Fluid : region(ispace(int3d), Fluid_columns),
                                     Grid_xBnum : int32, Grid_xCellWidth : double, Grid_xNum : int32,
                                     Grid_yBnum : int32, Grid_yCellWidth : double, Grid_yNum : int32,
                                     Grid_zBnum : int32, Grid_zCellWidth : double, Grid_zNum : int32)
where
  reads(Fluid.velocity),
  writes(Fluid.velocityGradientX),
  writes(Fluid.velocityGradientY),
  writes(Fluid.velocityGradientZ)
do
  __demand(__openmp)
  for c in Fluid do
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      Fluid[c].velocityGradientX = vs_div_double_3(vs_mul_double_3(vv_sub_double_3(Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocity, Fluid[((c+{-1, 0, 0})%Fluid.bounds)].velocity), double(0.5)), Grid_xCellWidth)
      Fluid[c].velocityGradientY = vs_div_double_3(vs_mul_double_3(vv_sub_double_3(Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocity, Fluid[((c+{0, -1, 0})%Fluid.bounds)].velocity), double(0.5)), Grid_yCellWidth)
      Fluid[c].velocityGradientZ = vs_div_double_3(vs_mul_double_3(vv_sub_double_3(Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocity, Fluid[((c+{0, 0, -1})%Fluid.bounds)].velocity), double(0.5)), Grid_zCellWidth)
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
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      var kineticEnergy = ((double(0.5)*Fluid[c].rho)*dot_double_3(Fluid[c].velocity, Fluid[c].velocity))
      var pressure = ((Flow_gamma-1.0)*(Fluid[c].rhoEnergy-kineticEnergy))
      Fluid[c].pressure = pressure
      Fluid[c].temperature = (pressure/(Flow_gasConstant*Fluid[c].rho))
    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdateGhostThermodynamicsStep1(Fluid : region(ispace(int3d), Fluid_columns),
                                         BC_xNegTemperature : double, BC_xPosTemperature : double,
                                         BC_yNegTemperature : double, BC_yPosTemperature : double,
                                         BC_zNegTemperature : double, BC_zPosTemperature : double,
                                         Grid_xBnum : int32, Grid_xNum : int32,
                                         Grid_yBnum : int32, Grid_yNum : int32,
                                         Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{pressure, temperature}),
  writes(Fluid.{pressureBoundary, temperatureBoundary})
do
  __demand(__openmp)
  for c in Fluid do
    if (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{1, 0, 0})%Fluid.bounds)
      var bnd_temperature = BC_xNegTemperature
      var temp_wall = double(0.0)
      var temperature = double(0.0)
      temp_wall = Fluid[c_int].temperature
      if (bnd_temperature>0.0) then
        temp_wall = bnd_temperature
      end
      temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
      Fluid[c_bnd].pressureBoundary = Fluid[c_int].pressure
      Fluid[c_bnd].temperatureBoundary = temperature
    end
    if (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{-1, 0, 0})%Fluid.bounds)
      var bnd_temperature = BC_xPosTemperature
      var temp_wall = double(0.0)
      var temperature = double(0.0)
      temp_wall = Fluid[c_int].temperature
      if (bnd_temperature>0.0) then
        temp_wall = bnd_temperature
      end
      temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
      Fluid[c_bnd].pressureBoundary = Fluid[c_int].pressure
      Fluid[c_bnd].temperatureBoundary = temperature
    end
    if (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, 1, 0})%Fluid.bounds)
      var bnd_temperature = BC_yNegTemperature
      var temp_wall = double(0.0)
      var temperature = double(0.0)
      temp_wall = Fluid[c_int].temperature
      if (bnd_temperature>0.0) then
        temp_wall = bnd_temperature
      end
      temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
      Fluid[c_bnd].pressureBoundary = Fluid[c_int].pressure
      Fluid[c_bnd].temperatureBoundary = temperature
    end
    if (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, -1, 0})%Fluid.bounds)
      var bnd_temperature = BC_yPosTemperature
      var temp_wall = double(0.0)
      var temperature = double(0.0)
      temp_wall = Fluid[c_int].temperature
      if (bnd_temperature>0.0) then
        temp_wall = bnd_temperature
      end
      temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
      Fluid[c_bnd].pressureBoundary = Fluid[c_int].pressure
      Fluid[c_bnd].temperatureBoundary = temperature
    end
    if (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, 0, 1})%Fluid.bounds)
      var bnd_temperature = BC_zNegTemperature
      var temp_wall = double(0.0)
      var temperature = double(0.0)
      temp_wall = Fluid[c_int].temperature
      if (bnd_temperature>0.0) then
        temp_wall = bnd_temperature
      end
      temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
      Fluid[c_bnd].pressureBoundary = Fluid[c_int].pressure
      Fluid[c_bnd].temperatureBoundary = temperature
    end
    if (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, 0, -1})%Fluid.bounds)
      var bnd_temperature = BC_zPosTemperature
      var temp_wall = double(0.0)
      var temperature = double(0.0)
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

__demand(__parallel, __cuda)
task Flow_UpdateGhostThermodynamicsStep2(Fluid : region(ispace(int3d), Fluid_columns),
                                         Grid_xBnum : int32, Grid_xNum : int32,
                                         Grid_yBnum : int32, Grid_yNum : int32,
                                         Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{pressureBoundary, temperatureBoundary}),
  writes(Fluid.{pressure, temperature})
do
  __demand(__openmp)
  for c in Fluid do
    if ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)) then
      Fluid[c].pressure = Fluid[c].pressureBoundary
      Fluid[c].temperature = Fluid[c].temperatureBoundary
    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdateGhostFieldsStep1(Fluid : region(ispace(int3d), Fluid_columns),
                                 BC_xNegTemperature : double, BC_xNegVelocity : double[3], BC_xPosTemperature : double, BC_xPosVelocity : double[3], BC_xSign : double[3],
                                 BC_yNegTemperature : double, BC_yNegVelocity : double[3], BC_yPosTemperature : double, BC_yPosVelocity : double[3], BC_ySign : double[3],
                                 BC_zNegTemperature : double, BC_zNegVelocity : double[3], BC_zPosTemperature : double, BC_zPosVelocity : double[3], BC_zSign : double[3],
                                 Flow_gamma : double,
                                 Flow_gasConstant : double,
                                 Grid_xBnum : int32, Grid_xNum : int32,
                                 Grid_yBnum : int32, Grid_yNum : int32,
                                 Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{rho, rhoVelocity, pressure, temperature}),
  writes(Fluid.{rhoBoundary, velocityBoundary, pressureBoundary, rhoVelocityBoundary, rhoEnergyBoundary, temperatureBoundary})
do
  __demand(__openmp)
  for c in Fluid do
    if (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{1, 0, 0})%Fluid.bounds)
      var sign = BC_xSign
      var bnd_velocity = BC_xNegVelocity
      var bnd_temperature = BC_xNegTemperature
      var rho = double(0.0)
      var temp_wall = double(0.0)
      var temperature = double(0.0)
      var velocity = [double[3]](array(0.0, 0.0, 0.0))
      var cv = (Flow_gasConstant/(Flow_gamma-1.0))
      velocity = vv_add_double_3(vv_mul_double_3(vs_div_double_3(Fluid[c_int].rhoVelocity, Fluid[c_int].rho), sign), bnd_velocity)
      temp_wall = Fluid[c_int].temperature
      if (bnd_temperature>0.0) then
        temp_wall = bnd_temperature
      end
      temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
      rho = (Fluid[c_int].pressure/(Flow_gasConstant*temperature))
      Fluid[c_bnd].rhoBoundary = rho
      Fluid[c_bnd].rhoVelocityBoundary = vs_mul_double_3(velocity, rho)
      Fluid[c_bnd].rhoEnergyBoundary = (rho*((cv*temperature)+(double(0.5)*dot_double_3(velocity, velocity))))
      Fluid[c_bnd].velocityBoundary = velocity
      Fluid[c_bnd].pressureBoundary = Fluid[c_int].pressure
      Fluid[c_bnd].temperatureBoundary = temperature
    end
    if (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{-1, 0, 0})%Fluid.bounds)
      var sign = BC_xSign
      var bnd_velocity = BC_xPosVelocity
      var bnd_temperature = BC_xPosTemperature
      var rho = double(0.0)
      var temp_wall = double(0.0)
      var temperature = double(0.0)
      var velocity = [double[3]](array(0.0, 0.0, 0.0))
      var cv = (Flow_gasConstant/(Flow_gamma-1.0))
      velocity = vv_add_double_3(vv_mul_double_3(vs_div_double_3(Fluid[c_int].rhoVelocity, Fluid[c_int].rho), sign), bnd_velocity)
      temp_wall = Fluid[c_int].temperature
      if (bnd_temperature>0.0) then
        temp_wall = bnd_temperature
      end
      temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
      rho = (Fluid[c_int].pressure/(Flow_gasConstant*temperature))
      Fluid[c_bnd].rhoBoundary = rho
      Fluid[c_bnd].rhoVelocityBoundary = vs_mul_double_3(velocity, rho)
      Fluid[c_bnd].rhoEnergyBoundary = (rho*((cv*temperature)+(double(0.5)*dot_double_3(velocity, velocity))))
      Fluid[c_bnd].velocityBoundary = velocity
      Fluid[c_bnd].pressureBoundary = Fluid[c_int].pressure
      Fluid[c_bnd].temperatureBoundary = temperature
    end
    if (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, 1, 0})%Fluid.bounds)
      var sign = BC_ySign
      var bnd_velocity = BC_yNegVelocity
      var bnd_temperature = BC_yNegTemperature
      var rho = double(0.0)
      var temp_wall = double(0.0)
      var temperature = double(0.0)
      var velocity = [double[3]](array(0.0, 0.0, 0.0))
      var cv = (Flow_gasConstant/(Flow_gamma-1.0))
      velocity = vv_add_double_3(vv_mul_double_3(vs_div_double_3(Fluid[c_int].rhoVelocity, Fluid[c_int].rho), sign), bnd_velocity)
      temp_wall = Fluid[c_int].temperature
      if (bnd_temperature>0.0) then
        temp_wall = bnd_temperature
      end
      temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
      rho = (Fluid[c_int].pressure/(Flow_gasConstant*temperature))
      Fluid[c_bnd].rhoBoundary = rho
      Fluid[c_bnd].rhoVelocityBoundary = vs_mul_double_3(velocity, rho)
      Fluid[c_bnd].rhoEnergyBoundary = (rho*((cv*temperature)+(double(0.5)*dot_double_3(velocity, velocity))))
      Fluid[c_bnd].velocityBoundary = velocity
      Fluid[c_bnd].pressureBoundary = Fluid[c_int].pressure
      Fluid[c_bnd].temperatureBoundary = temperature
    end
    if (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, -1, 0})%Fluid.bounds)
      var sign = BC_ySign
      var bnd_velocity = BC_yPosVelocity
      var bnd_temperature = BC_yPosTemperature
      var rho = double(0.0)
      var temp_wall = double(0.0)
      var temperature = double(0.0)
      var velocity = [double[3]](array(0.0, 0.0, 0.0))
      var cv = (Flow_gasConstant/(Flow_gamma-1.0))
      velocity = vv_add_double_3(vv_mul_double_3(vs_div_double_3(Fluid[c_int].rhoVelocity, Fluid[c_int].rho), sign), bnd_velocity)
      temp_wall = Fluid[c_int].temperature
      if (bnd_temperature>0.0) then
        temp_wall = bnd_temperature
      end
      temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
      rho = (Fluid[c_int].pressure/(Flow_gasConstant*temperature))
      Fluid[c_bnd].rhoBoundary = rho
      Fluid[c_bnd].rhoVelocityBoundary = vs_mul_double_3(velocity, rho)
      Fluid[c_bnd].rhoEnergyBoundary = (rho*((cv*temperature)+(double(0.5)*dot_double_3(velocity, velocity))))
      Fluid[c_bnd].velocityBoundary = velocity
      Fluid[c_bnd].pressureBoundary = Fluid[c_int].pressure
      Fluid[c_bnd].temperatureBoundary = temperature
    end
    if (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, 0, 1})%Fluid.bounds)
      var sign = BC_zSign
      var bnd_velocity = BC_zNegVelocity
      var bnd_temperature = BC_zNegTemperature
      var rho = double(0.0)
      var temp_wall = double(0.0)
      var temperature = double(0.0)
      var velocity = [double[3]](array(0.0, 0.0, 0.0))
      var cv = (Flow_gasConstant/(Flow_gamma-1.0))
      velocity = vv_add_double_3(vv_mul_double_3(vs_div_double_3(Fluid[c_int].rhoVelocity, Fluid[c_int].rho), sign), bnd_velocity)
      temp_wall = Fluid[c_int].temperature
      if (bnd_temperature>0.0) then
        temp_wall = bnd_temperature
      end
      temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
      rho = (Fluid[c_int].pressure/(Flow_gasConstant*temperature))
      Fluid[c_bnd].rhoBoundary = rho
      Fluid[c_bnd].rhoVelocityBoundary = vs_mul_double_3(velocity, rho)
      Fluid[c_bnd].rhoEnergyBoundary = (rho*((cv*temperature)+(double(0.5)*dot_double_3(velocity, velocity))))
      Fluid[c_bnd].velocityBoundary = velocity
      Fluid[c_bnd].pressureBoundary = Fluid[c_int].pressure
      Fluid[c_bnd].temperatureBoundary = temperature
    end
    if (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, 0, -1})%Fluid.bounds)
      var sign = BC_zSign
      var bnd_velocity = BC_zPosVelocity
      var bnd_temperature = BC_zPosTemperature
      var rho = double(0.0)
      var temp_wall = double(0.0)
      var temperature = double(0.0)
      var velocity = [double[3]](array(0.0, 0.0, 0.0))
      var cv = (Flow_gasConstant/(Flow_gamma-1.0))
      velocity = vv_add_double_3(vv_mul_double_3(vs_div_double_3(Fluid[c_int].rhoVelocity, Fluid[c_int].rho), sign), bnd_velocity)
      temp_wall = Fluid[c_int].temperature
      if (bnd_temperature>0.0) then
        temp_wall = bnd_temperature
      end
      temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
      rho = (Fluid[c_int].pressure/(Flow_gasConstant*temperature))
      Fluid[c_bnd].rhoBoundary = rho
      Fluid[c_bnd].rhoVelocityBoundary = vs_mul_double_3(velocity, rho)
      Fluid[c_bnd].rhoEnergyBoundary = (rho*((cv*temperature)+(double(0.5)*dot_double_3(velocity, velocity))))
      Fluid[c_bnd].velocityBoundary = velocity
      Fluid[c_bnd].pressureBoundary = Fluid[c_int].pressure
      Fluid[c_bnd].temperatureBoundary = temperature
    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdateGhostFieldsStep2(Fluid : region(ispace(int3d), Fluid_columns),
                                 Grid_xBnum : int32, Grid_xNum : int32,
                                 Grid_yBnum : int32, Grid_yNum : int32,
                                 Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{rhoBoundary, rhoVelocityBoundary, pressureBoundary, rhoEnergyBoundary, temperatureBoundary}),
  writes(Fluid.{rho, pressure, rhoEnergy, rhoVelocity, temperature})
do
  __demand(__openmp)
  for c in Fluid do
    if ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)) then
      Fluid[c].rho = Fluid[c].rhoBoundary
      Fluid[c].rhoVelocity = Fluid[c].rhoVelocityBoundary
      Fluid[c].rhoEnergy = Fluid[c].rhoEnergyBoundary
      Fluid[c].pressure = Fluid[c].pressureBoundary
      Fluid[c].temperature = Fluid[c].temperatureBoundary
    end
  end
end

__demand(__parallel, __cuda)
task Particles_InitializeDensity(particles : region(ispace(int1d), particles_columns),
                                 Particles_density : double)
where
  reads(particles.__valid),
  writes(particles.density)
do
  __demand(__openmp)
  for p in particles do
    if particles[p].__valid then
      particles[p].density = Particles_density
    end
  end
end

__demand(__parallel)
task Particles_CalculateNumber(particles : region(ispace(int1d), particles_columns)) : int64
where
  reads(particles.__valid)
do
  var acc = int64(0)
  __demand(__openmp)
  for p in particles do
    if particles[p].__valid then
      acc += int64(1)
    end
  end
  return acc
end

__demand(__parallel)
task CalculateAveragePressure(Fluid : region(ispace(int3d), Fluid_columns),
                              Grid_cellVolume : double,
                              Grid_xBnum : int32, Grid_xNum : int32,
                              Grid_yBnum : int32, Grid_yNum : int32,
                              Grid_zBnum : int32, Grid_zNum : int32) : double
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

__demand(__parallel)
task CalculateAverageTemperature(Fluid : region(ispace(int3d), Fluid_columns),
                                 Grid_cellVolume : double,
                                 Grid_xBnum : int32, Grid_xNum : int32,
                                 Grid_yBnum : int32, Grid_yNum : int32,
                                 Grid_zBnum : int32, Grid_zNum : int32) : double
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

__demand(__parallel)
task CalculateAverageKineticEnergy(Fluid : region(ispace(int3d), Fluid_columns),
                                   Grid_cellVolume : double,
                                   Grid_xBnum : int32, Grid_xNum : int32,
                                   Grid_yBnum : int32, Grid_yNum : int32,
                                   Grid_zBnum : int32, Grid_zNum : int32) : double
where
  reads(Fluid.kineticEnergy)
do
  var acc = 0.0
  __demand(__openmp)
  for c in Fluid do
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      acc += (Fluid[c].kineticEnergy*Grid_cellVolume)
    end
  end
  return acc
end

__demand(__parallel)
task CalculateMinTemperature(Fluid : region(ispace(int3d), Fluid_columns),
                             Grid_xBnum : int32, Grid_xNum : int32,
                             Grid_yBnum : int32, Grid_yNum : int32,
                             Grid_zBnum : int32, Grid_zNum : int32) : double
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

__demand(__parallel)
task CalculateMaxTemperature(Fluid : region(ispace(int3d), Fluid_columns),
                             Grid_xBnum : int32, Grid_xNum : int32,
                             Grid_yBnum : int32, Grid_yNum : int32,
                             Grid_zBnum : int32, Grid_zNum : int32) : double
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

__demand(__parallel)
task Particles_IntegrateQuantities(particles : region(ispace(int1d), particles_columns)) : double
where
  reads(particles.{temperature, __valid})
do
  var acc = 0.0
  __demand(__openmp)
  for p in particles do
    if particles[p].__valid then
      acc += particles[p].temperature
    end
  end
  return acc
end

__demand(__parallel, __cuda)
task Radiation_InitializeCell(Radiation : region(ispace(int3d), Radiation_columns))
where
  reads(Radiation.G), writes(Radiation.G),
  reads(Radiation.I_1), writes(Radiation.I_1),
  reads(Radiation.I_2), writes(Radiation.I_2),
  reads(Radiation.I_3), writes(Radiation.I_3),
  reads(Radiation.I_4), writes(Radiation.I_4),
  reads(Radiation.I_5), writes(Radiation.I_5),
  reads(Radiation.I_6), writes(Radiation.I_6),
  reads(Radiation.I_7), writes(Radiation.I_7),
  reads(Radiation.I_8), writes(Radiation.I_8),
  reads(Radiation.Iiter_1), writes(Radiation.Iiter_1),
  reads(Radiation.Iiter_2), writes(Radiation.Iiter_2),
  reads(Radiation.Iiter_3), writes(Radiation.Iiter_3),
  reads(Radiation.Iiter_4), writes(Radiation.Iiter_4),
  reads(Radiation.Iiter_5), writes(Radiation.Iiter_5),
  reads(Radiation.Iiter_6), writes(Radiation.Iiter_6),
  reads(Radiation.Iiter_7), writes(Radiation.Iiter_7),
  reads(Radiation.Iiter_8), writes(Radiation.Iiter_8),
  reads(Radiation.S), writes(Radiation.S)
do
  __demand(__openmp)
  for c in Radiation do
    for m : int32 = 0, NUM_ANGLES do
      Radiation[c].I_1[m] = 0.0
      Radiation[c].I_2[m] = 0.0
      Radiation[c].I_3[m] = 0.0
      Radiation[c].I_4[m] = 0.0
      Radiation[c].I_5[m] = 0.0
      Radiation[c].I_6[m] = 0.0
      Radiation[c].I_7[m] = 0.0
      Radiation[c].I_8[m] = 0.0
      Radiation[c].Iiter_1[m] = 0.0
      Radiation[c].Iiter_2[m] = 0.0
      Radiation[c].Iiter_3[m] = 0.0
      Radiation[c].Iiter_4[m] = 0.0
      Radiation[c].Iiter_5[m] = 0.0
      Radiation[c].Iiter_6[m] = 0.0
      Radiation[c].Iiter_7[m] = 0.0
      Radiation[c].Iiter_8[m] = 0.0
    end
    Radiation[c].G = 0.0
    Radiation[c].S = 0.0
  end
end

__demand(__inline)
task GetSoundSpeed(temperature : double, Flow_gamma : double, Flow_gasConstant : double) : double
  return sqrt(((Flow_gamma*Flow_gasConstant)*temperature))
end

__demand(__parallel)
task CalculateConvectiveSpectralRadius(Fluid : region(ispace(int3d), Fluid_columns),
                                       Flow_gamma : double,
                                       Flow_gasConstant : double,
                                       Grid_dXYZInverseSquare : double,
                                       Grid_xCellWidth : double, Grid_yCellWidth : double, Grid_zCellWidth : double) : double
where
  reads(Fluid.{velocity, temperature}),
  reads writes(Fluid.convectiveSpectralRadius)
do
  var acc = -math.huge
  __demand(__openmp)
  for c in Fluid do
    Fluid[c].convectiveSpectralRadius = ((((fabs(Fluid[c].velocity[0])/Grid_xCellWidth)+(fabs(Fluid[c].velocity[1])/Grid_yCellWidth))+(fabs(Fluid[c].velocity[2])/Grid_zCellWidth))+(GetSoundSpeed(Fluid[c].temperature, Flow_gamma, Flow_gasConstant)*sqrt(Grid_dXYZInverseSquare)))
    acc max= Fluid[c].convectiveSpectralRadius
  end
  return acc
end

__demand(__inline)
task GetDynamicViscosity(temperature : double,
                         Flow_constantVisc : double,
                         Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                         Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                         Flow_viscosityModel : SCHEMA.ViscosityModel) : double
  var viscosity = double(0.0)
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

__demand(__parallel)
task CalculateViscousSpectralRadius(Fluid : region(ispace(int3d), Fluid_columns),
                                    Flow_constantVisc : double,
                                    Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                                    Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                                    Flow_viscosityModel : SCHEMA.ViscosityModel,
                                    Grid_dXYZInverseSquare : double) : double
where
  reads(Fluid.{rho, temperature, sgsEddyViscosity}),
  reads writes(Fluid.viscousSpectralRadius)
do
  var acc = -math.huge
  __demand(__openmp)
  for c in Fluid do
    var dynamicViscosity = GetDynamicViscosity(Fluid[c].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)
    var eddyViscosity = Fluid[c].sgsEddyViscosity
    Fluid[c].viscousSpectralRadius = ((((2.0*(dynamicViscosity+eddyViscosity))/Fluid[c].rho)*Grid_dXYZInverseSquare)*4.0)
    acc max= Fluid[c].viscousSpectralRadius
  end
  return acc
end

__demand(__parallel)
task CalculateHeatConductionSpectralRadius(Fluid : region(ispace(int3d), Fluid_columns),
                                           Flow_constantVisc : double,
                                           Flow_gamma : double,
                                           Flow_gasConstant : double,
                                           Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                                           Flow_prandtl : double,
                                           Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                                           Flow_viscosityModel : SCHEMA.ViscosityModel,
                                           Grid_dXYZInverseSquare : double) : double
where
  reads(Fluid.{rho, temperature, sgsEddyKappa}),
  reads writes(Fluid.heatConductionSpectralRadius)
do
  var acc = -math.huge
  __demand(__openmp)
  for c in Fluid do
    var dynamicViscosity = GetDynamicViscosity(Fluid[c].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)
    var cv = (Flow_gasConstant/(Flow_gamma-1.0))
    var cp = (Flow_gamma*cv)
    var kappa = ((cp/Flow_prandtl)*dynamicViscosity)
    Fluid[c].heatConductionSpectralRadius = ((((kappa+Fluid[c].sgsEddyKappa)/(cv*Fluid[c].rho))*Grid_dXYZInverseSquare)*4.0)
    acc max= Fluid[c].heatConductionSpectralRadius
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
task Particles_InitializeTemporaries(particles : region(ispace(int1d), particles_columns))
where
  reads(particles.{position, velocity, temperature, __valid}),
  writes(particles.{position_new, position_old, temperature_new, temperature_old, velocity_new, velocity_old})
do
  __demand(__openmp)
  for p in particles do
    if particles[p].__valid then
      particles[p].position_old = particles[p].position
      particles[p].velocity_old = particles[p].velocity
      particles[p].temperature_old = particles[p].temperature
      particles[p].position_new = particles[p].position
      particles[p].velocity_new = particles[p].velocity
      particles[p].temperature_new = particles[p].temperature
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
    Fluid[c].rho_t = double(0.0)
    Fluid[c].rhoVelocity_t = [double[3]](array(0.0, 0.0, 0.0))
    Fluid[c].rhoEnergy_t = double(0.0)
    Fluid[c].rhoEnthalpy = (Fluid[c].rhoEnergy+Fluid[c].pressure)
  end
end

__demand(__parallel, __cuda)
task Particles_InitializeTimeDerivatives(particles : region(ispace(int1d), particles_columns))
where
  reads(particles.__valid),
  writes(particles.{position_t, velocity_t, temperature_t})
do
  __demand(__openmp)
  for p in particles do
    if particles[p].__valid then
      particles[p].position_t = [double[3]](array(0.0, 0.0, 0.0))
      particles[p].velocity_t = [double[3]](array(0.0, 0.0, 0.0))
      particles[p].temperature_t = 0.0
    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdateGhostVelocityGradientStep1(Fluid : region(ispace(int3d), Fluid_columns),
                                           BC_xSign : double[3], BC_ySign : double[3], BC_zSign : double[3],
                                           Grid_xBnum : int32, Grid_xNum : int32,
                                           Grid_yBnum : int32, Grid_yNum : int32,
                                           Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{velocityGradientX, velocityGradientY, velocityGradientZ}),
  writes(Fluid.{velocityGradientXBoundary, velocityGradientYBoundary, velocityGradientZBoundary})
do
  __demand(__openmp)
  for c in Fluid do
    if (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{1, 0, 0})%Fluid.bounds)
      var sign = BC_xSign
      Fluid[c_bnd].velocityGradientXBoundary = vv_mul_double_3(sign, Fluid[c_int].velocityGradientX)
      Fluid[c_bnd].velocityGradientYBoundary = vv_mul_double_3(sign, Fluid[c_int].velocityGradientY)
      Fluid[c_bnd].velocityGradientZBoundary = vv_mul_double_3(sign, Fluid[c_int].velocityGradientZ)
    end
    if (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{-1, 0, 0})%Fluid.bounds)
      var sign = BC_xSign
      Fluid[c_bnd].velocityGradientXBoundary = vv_mul_double_3(sign, Fluid[c_int].velocityGradientX)
      Fluid[c_bnd].velocityGradientYBoundary = vv_mul_double_3(sign, Fluid[c_int].velocityGradientY)
      Fluid[c_bnd].velocityGradientZBoundary = vv_mul_double_3(sign, Fluid[c_int].velocityGradientZ)
    end
    if (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, 1, 0})%Fluid.bounds)
      var sign = BC_ySign
      Fluid[c_bnd].velocityGradientXBoundary = vv_mul_double_3(sign, Fluid[c_int].velocityGradientX)
      Fluid[c_bnd].velocityGradientYBoundary = vv_mul_double_3(sign, Fluid[c_int].velocityGradientY)
      Fluid[c_bnd].velocityGradientZBoundary = vv_mul_double_3(sign, Fluid[c_int].velocityGradientZ)
    end
    if (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, -1, 0})%Fluid.bounds)
      var sign = BC_ySign
      Fluid[c_bnd].velocityGradientXBoundary = vv_mul_double_3(sign, Fluid[c_int].velocityGradientX)
      Fluid[c_bnd].velocityGradientYBoundary = vv_mul_double_3(sign, Fluid[c_int].velocityGradientY)
      Fluid[c_bnd].velocityGradientZBoundary = vv_mul_double_3(sign, Fluid[c_int].velocityGradientZ)
    end
    if (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, 0, 1})%Fluid.bounds)
      var sign = BC_zSign
      Fluid[c_bnd].velocityGradientXBoundary = vv_mul_double_3(sign, Fluid[c_int].velocityGradientX)
      Fluid[c_bnd].velocityGradientYBoundary = vv_mul_double_3(sign, Fluid[c_int].velocityGradientY)
      Fluid[c_bnd].velocityGradientZBoundary = vv_mul_double_3(sign, Fluid[c_int].velocityGradientZ)
    end
    if (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0) then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, 0, -1})%Fluid.bounds)
      var sign = BC_zSign
      Fluid[c_bnd].velocityGradientXBoundary = vv_mul_double_3(sign, Fluid[c_int].velocityGradientX)
      Fluid[c_bnd].velocityGradientYBoundary = vv_mul_double_3(sign, Fluid[c_int].velocityGradientY)
      Fluid[c_bnd].velocityGradientZBoundary = vv_mul_double_3(sign, Fluid[c_int].velocityGradientZ)
    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdateGhostVelocityGradientStep2(Fluid : region(ispace(int3d), Fluid_columns),
                                           Grid_xBnum : int32, Grid_xNum : int32,
                                           Grid_yBnum : int32, Grid_yNum : int32,
                                           Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{velocityGradientXBoundary, velocityGradientYBoundary, velocityGradientZBoundary}),
  writes(Fluid.{velocityGradientX, velocityGradientY, velocityGradientZ})
do
  __demand(__openmp)
  for c in Fluid do
    if ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0)) then
      Fluid[c].velocityGradientX = Fluid[c].velocityGradientXBoundary
      Fluid[c].velocityGradientY = Fluid[c].velocityGradientYBoundary
      Fluid[c].velocityGradientZ = Fluid[c].velocityGradientZBoundary
    end
  end
end

__demand(__inline)
task CenteredInviscidFlux(c_l : int3d,
                          c_r : int3d,
                          Fluid : region(ispace(int3d), Fluid_columns)) : double[5]
where
  reads(Fluid.{rho, velocity, pressure, rhoVelocity, rhoEnthalpy})
do
  var rhoFactorDiagonal = double(0.0)
  var rhoVelocityFactorDiagonal = [double[3]](array(0.0, 0.0, 0.0))
  var rhoEnergyFactorDiagonal = double(0.0)
  var fpdiag = double(0.0)
  rhoFactorDiagonal = (double(0.5)*((Fluid[c_l].rho*Fluid[c_l].velocity[0])+(Fluid[c_r].rho*Fluid[c_r].velocity[0])))
  rhoVelocityFactorDiagonal = vs_mul_double_3(vv_add_double_3(vs_mul_double_3(Fluid[c_l].rhoVelocity, Fluid[c_l].velocity[0]), vs_mul_double_3(Fluid[c_r].rhoVelocity, Fluid[c_r].velocity[0])), double(0.5))
  rhoEnergyFactorDiagonal = (double(0.5)*((Fluid[c_l].rhoEnthalpy*Fluid[c_l].velocity[0])+(Fluid[c_r].rhoEnthalpy*Fluid[c_r].velocity[0])))
  fpdiag += (double(0.5)*(Fluid[c_l].pressure+Fluid[c_r].pressure))
  var rhoFactorSkew = double(0.0)
  var rhoVelocityFactorSkew = [double[3]](array(0.0, 0.0, 0.0))
  var rhoEnergyFactorSkew = double(0.0)
  var tmp = double(0.0)
  tmp = (double(0.5)*Fluid[c_r].velocity[0])
  rhoFactorSkew += (Fluid[c_l].rho*tmp)
  var tmp__6137 = vs_mul_double_3(Fluid[c_l].rhoVelocity, tmp)
  var v = rhoVelocityFactorSkew
  v[0] += tmp__6137[0]
  v[1] += tmp__6137[1]
  v[2] += tmp__6137[2]
  rhoVelocityFactorSkew = v
  rhoEnergyFactorSkew += (Fluid[c_l].rhoEnthalpy*tmp)
  tmp = (double(0.5)*Fluid[c_l].velocity[0])
  rhoFactorSkew += (Fluid[c_r].rho*tmp)
  var tmp__6139 = vs_mul_double_3(Fluid[c_r].rhoVelocity, tmp)
  var v__6140 = rhoVelocityFactorSkew
  v__6140[0] += tmp__6139[0]
  v__6140[1] += tmp__6139[1]
  v__6140[2] += tmp__6139[2]
  rhoVelocityFactorSkew = v__6140
  rhoEnergyFactorSkew += (Fluid[c_r].rhoEnthalpy*tmp)
  var s = double(0.5)
  var rhoFlux_temp = ((s*rhoFactorDiagonal)+((1.0-s)*rhoFactorSkew))
  var rhoVelocityFlux_temp = vv_add_double_3(vs_mul_double_3(rhoVelocityFactorDiagonal, s), vs_mul_double_3(rhoVelocityFactorSkew, (1.0-s)))
  var rhoEnergyFlux_temp = ((s*rhoEnergyFactorDiagonal)+((1.0-s)*rhoEnergyFactorSkew))
  rhoVelocityFlux_temp[0] += fpdiag
  return array(rhoFlux_temp, rhoVelocityFlux_temp[0], rhoVelocityFlux_temp[1], rhoVelocityFlux_temp[2], rhoEnergyFlux_temp)
end

__demand(__inline)
task CenteredInviscidFlux_(c_l : int3d,
                           c_r : int3d,
                           Fluid : region(ispace(int3d), Fluid_columns)) : double[5]
where
  reads(Fluid.{rho, pressure, velocity, rhoVelocity, rhoEnthalpy})
do
  var rhoFactorDiagonal = double(0.0)
  var rhoVelocityFactorDiagonal = [double[3]](array(0.0, 0.0, 0.0))
  var rhoEnergyFactorDiagonal = double(0.0)
  var fpdiag = double(0.0)
  rhoFactorDiagonal = (double(0.5)*((Fluid[c_l].rho*Fluid[c_l].velocity[1])+(Fluid[c_r].rho*Fluid[c_r].velocity[1])))
  rhoVelocityFactorDiagonal = vs_mul_double_3(vv_add_double_3(vs_mul_double_3(Fluid[c_l].rhoVelocity, Fluid[c_l].velocity[1]), vs_mul_double_3(Fluid[c_r].rhoVelocity, Fluid[c_r].velocity[1])), double(0.5))
  rhoEnergyFactorDiagonal = (double(0.5)*((Fluid[c_l].rhoEnthalpy*Fluid[c_l].velocity[1])+(Fluid[c_r].rhoEnthalpy*Fluid[c_r].velocity[1])))
  fpdiag += (double(0.5)*(Fluid[c_l].pressure+Fluid[c_r].pressure))
  var rhoFactorSkew = double(0.0)
  var rhoVelocityFactorSkew = [double[3]](array(0.0, 0.0, 0.0))
  var rhoEnergyFactorSkew = double(0.0)
  var tmp = double(0.0)
  tmp = (double(0.5)*Fluid[c_r].velocity[1])
  rhoFactorSkew += (Fluid[c_l].rho*tmp)
  var tmp__6344 = vs_mul_double_3(Fluid[c_l].rhoVelocity, tmp)
  var v = rhoVelocityFactorSkew
  v[0] += tmp__6344[0]
  v[1] += tmp__6344[1]
  v[2] += tmp__6344[2]
  rhoVelocityFactorSkew = v
  rhoEnergyFactorSkew += (Fluid[c_l].rhoEnthalpy*tmp)
  tmp = (double(0.5)*Fluid[c_l].velocity[1])
  rhoFactorSkew += (Fluid[c_r].rho*tmp)
  var tmp__6346 = vs_mul_double_3(Fluid[c_r].rhoVelocity, tmp)
  var v__6347 = rhoVelocityFactorSkew
  v__6347[0] += tmp__6346[0]
  v__6347[1] += tmp__6346[1]
  v__6347[2] += tmp__6346[2]
  rhoVelocityFactorSkew = v__6347
  rhoEnergyFactorSkew += (Fluid[c_r].rhoEnthalpy*tmp)
  var s = double(0.5)
  var rhoFlux_temp = ((s*rhoFactorDiagonal)+((1.0-s)*rhoFactorSkew))
  var rhoVelocityFlux_temp = vv_add_double_3(vs_mul_double_3(rhoVelocityFactorDiagonal, s), vs_mul_double_3(rhoVelocityFactorSkew, (1.0-s)))
  var rhoEnergyFlux_temp = ((s*rhoEnergyFactorDiagonal)+((1.0-s)*rhoEnergyFactorSkew))
  rhoVelocityFlux_temp[1] += fpdiag
  return array(rhoFlux_temp, rhoVelocityFlux_temp[0], rhoVelocityFlux_temp[1], rhoVelocityFlux_temp[2], rhoEnergyFlux_temp)
end

__demand(__inline)
task CenteredInviscidFlux__(c_l : int3d,
                            c_r : int3d,
                            Fluid : region(ispace(int3d), Fluid_columns)) : double[5]
where
  reads(Fluid.{rho, pressure, velocity, rhoVelocity, rhoEnthalpy})
do
  var rhoFactorDiagonal = double(0.0)
  var rhoVelocityFactorDiagonal = [double[3]](array(0.0, 0.0, 0.0))
  var rhoEnergyFactorDiagonal = double(0.0)
  var fpdiag = double(0.0)
  rhoFactorDiagonal = (double(0.5)*((Fluid[c_l].rho*Fluid[c_l].velocity[2])+(Fluid[c_r].rho*Fluid[c_r].velocity[2])))
  rhoVelocityFactorDiagonal = vs_mul_double_3(vv_add_double_3(vs_mul_double_3(Fluid[c_l].rhoVelocity, Fluid[c_l].velocity[2]), vs_mul_double_3(Fluid[c_r].rhoVelocity, Fluid[c_r].velocity[2])), double(0.5))
  rhoEnergyFactorDiagonal = (double(0.5)*((Fluid[c_l].rhoEnthalpy*Fluid[c_l].velocity[2])+(Fluid[c_r].rhoEnthalpy*Fluid[c_r].velocity[2])))
  fpdiag += (double(0.5)*(Fluid[c_l].pressure+Fluid[c_r].pressure))
  var rhoFactorSkew = double(0.0)
  var rhoVelocityFactorSkew = [double[3]](array(0.0, 0.0, 0.0))
  var rhoEnergyFactorSkew = double(0.0)
  var tmp = double(0.0)
  tmp = (double(0.5)*Fluid[c_r].velocity[2])
  rhoFactorSkew += (Fluid[c_l].rho*tmp)
  var tmp__6551 = vs_mul_double_3(Fluid[c_l].rhoVelocity, tmp)
  var v = rhoVelocityFactorSkew
  v[0] += tmp__6551[0]
  v[1] += tmp__6551[1]
  v[2] += tmp__6551[2]
  rhoVelocityFactorSkew = v
  rhoEnergyFactorSkew += (Fluid[c_l].rhoEnthalpy*tmp)
  tmp = (double(0.5)*Fluid[c_l].velocity[2])
  rhoFactorSkew += (Fluid[c_r].rho*tmp)
  var tmp__6553 = vs_mul_double_3(Fluid[c_r].rhoVelocity, tmp)
  var v__6554 = rhoVelocityFactorSkew
  v__6554[0] += tmp__6553[0]
  v__6554[1] += tmp__6553[1]
  v__6554[2] += tmp__6553[2]
  rhoVelocityFactorSkew = v__6554
  rhoEnergyFactorSkew += (Fluid[c_r].rhoEnthalpy*tmp)
  var s = double(0.5)
  var rhoFlux_temp = ((s*rhoFactorDiagonal)+((1.0-s)*rhoFactorSkew))
  var rhoVelocityFlux_temp = vv_add_double_3(vs_mul_double_3(rhoVelocityFactorDiagonal, s), vs_mul_double_3(rhoVelocityFactorSkew, (1.0-s)))
  var rhoEnergyFlux_temp = ((s*rhoEnergyFactorDiagonal)+((1.0-s)*rhoEnergyFactorSkew))
  rhoVelocityFlux_temp[2] += fpdiag
  return array(rhoFlux_temp, rhoVelocityFlux_temp[0], rhoVelocityFlux_temp[1], rhoVelocityFlux_temp[2], rhoEnergyFlux_temp)
end

__demand(__parallel, __cuda)
task Flow_AddGetFlux(Fluid : region(ispace(int3d), Fluid_columns),
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
  __demand(__openmp)
  for c in Fluid do
    if ((not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) or (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)==1)) then
      var flux = CenteredInviscidFlux(int3d(c), ((c+{1, 0, 0})%Fluid.bounds), Fluid)
      Fluid[c].rhoFluxX = flux[0]
      Fluid[c].rhoVelocityFluxX = array(flux[1], flux[2], flux[3])
      Fluid[c].rhoEnergyFluxX = flux[4]
      var muFace = (double(0.5)*(GetDynamicViscosity(Fluid[c].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)+GetDynamicViscosity(Fluid[((c+{1, 0, 0})%Fluid.bounds)].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)))
      var velocityFace = [double[3]](array(0.0, 0.0, 0.0))
      var velocityX_YFace = double(0.0)
      var velocityX_ZFace = double(0.0)
      var velocityY_YFace = double(0.0)
      var velocityZ_ZFace = double(0.0)
      velocityFace = vs_mul_double_3(vv_add_double_3(Fluid[c].velocity, Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocity), double(0.5))
      velocityX_YFace = (double(0.5)*(Fluid[c].velocityGradientY[0]+Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocityGradientY[0]))
      velocityX_ZFace = (double(0.5)*(Fluid[c].velocityGradientZ[0]+Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocityGradientZ[0]))
      velocityY_YFace = (double(0.5)*(Fluid[c].velocityGradientY[1]+Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocityGradientY[1]))
      velocityZ_ZFace = (double(0.5)*(Fluid[c].velocityGradientZ[2]+Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocityGradientZ[2]))
      var velocityX_XFace = double(0.0)
      var velocityY_XFace = double(0.0)
      var velocityZ_XFace = double(0.0)
      var temperature_XFace = double(0.0)
      velocityX_XFace = (double(0.5)*(Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocity[0]-Fluid[c].velocity[0]))
      velocityY_XFace = (double(0.5)*(Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocity[1]-Fluid[c].velocity[1]))
      velocityZ_XFace = (double(0.5)*(Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocity[2]-Fluid[c].velocity[2]))
      temperature_XFace = (double(0.5)*(Fluid[((c+{1, 0, 0})%Fluid.bounds)].temperature-Fluid[c].temperature))
      velocityX_XFace *= (1/(Grid_xCellWidth*double(0.5)))
      velocityY_XFace *= (1/(Grid_xCellWidth*double(0.5)))
      velocityZ_XFace *= (1/(Grid_xCellWidth*double(0.5)))
      temperature_XFace *= (1/(Grid_xCellWidth*double(0.5)))
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
    if ((not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)==1)) then
      var flux = CenteredInviscidFlux_(int3d(c), ((c+{0, 1, 0})%Fluid.bounds), Fluid)
      Fluid[c].rhoFluxY = flux[0]
      Fluid[c].rhoVelocityFluxY = array(flux[1], flux[2], flux[3])
      Fluid[c].rhoEnergyFluxY = flux[4]
      var muFace = (double(0.5)*(GetDynamicViscosity(Fluid[c].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)+GetDynamicViscosity(Fluid[((c+{0, 1, 0})%Fluid.bounds)].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)))
      var velocityFace = [double[3]](array(0.0, 0.0, 0.0))
      var velocityY_XFace = double(0.0)
      var velocityY_ZFace = double(0.0)
      var velocityX_XFace = double(0.0)
      var velocityZ_ZFace = double(0.0)
      velocityFace = vs_mul_double_3(vv_add_double_3(Fluid[c].velocity, Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocity), double(0.5))
      velocityY_XFace = (double(0.5)*(Fluid[c].velocityGradientX[1]+Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocityGradientX[1]))
      velocityY_ZFace = (double(0.5)*(Fluid[c].velocityGradientZ[1]+Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocityGradientZ[1]))
      velocityX_XFace = (double(0.5)*(Fluid[c].velocityGradientX[0]+Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocityGradientX[0]))
      velocityZ_ZFace = (double(0.5)*(Fluid[c].velocityGradientZ[2]+Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocityGradientZ[2]))
      var velocityX_YFace = double(0.0)
      var velocityY_YFace = double(0.0)
      var velocityZ_YFace = double(0.0)
      var temperature_YFace = double(0.0)
      velocityX_YFace = (double(0.5)*(Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocity[0]-Fluid[c].velocity[0]))
      velocityY_YFace = (double(0.5)*(Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocity[1]-Fluid[c].velocity[1]))
      velocityZ_YFace = (double(0.5)*(Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocity[2]-Fluid[c].velocity[2]))
      temperature_YFace = (double(0.5)*(Fluid[((c+{0, 1, 0})%Fluid.bounds)].temperature-Fluid[c].temperature))
      velocityX_YFace *= (1/(Grid_yCellWidth*double(0.5)))
      velocityY_YFace *= (1/(Grid_yCellWidth*double(0.5)))
      velocityZ_YFace *= (1/(Grid_yCellWidth*double(0.5)))
      temperature_YFace *= (1/(Grid_yCellWidth*double(0.5)))
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
    if ((not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)==1)) then
      var flux = CenteredInviscidFlux__(int3d(c), ((c+{0, 0, 1})%Fluid.bounds), Fluid)
      Fluid[c].rhoFluxZ = flux[0]
      Fluid[c].rhoVelocityFluxZ = array(flux[1], flux[2], flux[3])
      Fluid[c].rhoEnergyFluxZ = flux[4]
      var muFace = (double(0.5)*(GetDynamicViscosity(Fluid[c].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)+GetDynamicViscosity(Fluid[((c+{0, 0, 1})%Fluid.bounds)].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)))
      var velocityFace = [double[3]](array(0.0, 0.0, 0.0))
      var velocityZ_XFace = double(0.0)
      var velocityZ_YFace = double(0.0)
      var velocityX_XFace = double(0.0)
      var velocityY_YFace = double(0.0)
      velocityFace = vs_mul_double_3(vv_add_double_3(Fluid[c].velocity, Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocity), double(0.5))
      velocityZ_XFace = (double(0.5)*(Fluid[c].velocityGradientX[2]+Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocityGradientX[2]))
      velocityZ_YFace = (double(0.5)*(Fluid[c].velocityGradientY[2]+Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocityGradientY[2]))
      velocityX_XFace = (double(0.5)*(Fluid[c].velocityGradientX[0]+Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocityGradientX[0]))
      velocityY_YFace = (double(0.5)*(Fluid[c].velocityGradientY[1]+Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocityGradientY[1]))
      var velocityX_ZFace = double(0.0)
      var velocityY_ZFace = double(0.0)
      var velocityZ_ZFace = double(0.0)
      var temperature_ZFace = double(0.0)
      velocityX_ZFace = (double(0.5)*(Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocity[0]-Fluid[c].velocity[0]))
      velocityY_ZFace = (double(0.5)*(Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocity[1]-Fluid[c].velocity[1]))
      velocityZ_ZFace = (double(0.5)*(Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocity[2]-Fluid[c].velocity[2]))
      temperature_ZFace = (double(0.5)*(Fluid[((c+{0, 0, 1})%Fluid.bounds)].temperature-Fluid[c].temperature))
      velocityX_ZFace *= (1/(Grid_zCellWidth*double(0.5)))
      velocityY_ZFace *= (1/(Grid_zCellWidth*double(0.5)))
      velocityZ_ZFace *= (1/(Grid_zCellWidth*double(0.5)))
      temperature_ZFace *= (1/(Grid_zCellWidth*double(0.5)))
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
    var v = 0
    if (v==1) then
      var tmp1 = Fluid[((c+{-1, 0, 0})%Fluid.bounds)].velocity[0]
      var tmp2 = Fluid[((c+{0, -1, 0})%Fluid.bounds)].velocity[1]
      var tmp3 = Fluid[((c+{0, 0, -1})%Fluid.bounds)].velocity[2]
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
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      Fluid[c].rho_t += ((-(Fluid[c].rhoFluxX-Fluid[((c+{-1, 0, 0})%Fluid.bounds)].rhoFluxX))/Grid_xCellWidth)
      var tmp = vs_div_double_3(vs_mul_double_3(vv_sub_double_3(Fluid[c].rhoVelocityFluxX, Fluid[((c+{-1, 0, 0})%Fluid.bounds)].rhoVelocityFluxX), double((-1))), Grid_xCellWidth)
      var v = Fluid[c].rhoVelocity_t
      v[0] += tmp[0]
      v[1] += tmp[1]
      v[2] += tmp[2]
      Fluid[c].rhoVelocity_t = v
      Fluid[c].rhoEnergy_t += ((-(Fluid[c].rhoEnergyFluxX-Fluid[((c+{-1, 0, 0})%Fluid.bounds)].rhoEnergyFluxX))/Grid_xCellWidth)
      Fluid[c].rho_t += ((-(Fluid[c].rhoFluxY-Fluid[((c+{0, -1, 0})%Fluid.bounds)].rhoFluxY))/Grid_yCellWidth)
      var tmp__7144 = vs_div_double_3(vs_mul_double_3(vv_sub_double_3(Fluid[c].rhoVelocityFluxY, Fluid[((c+{0, -1, 0})%Fluid.bounds)].rhoVelocityFluxY), double((-1))), Grid_yCellWidth)
      var v__7145 = Fluid[c].rhoVelocity_t
      v__7145[0] += tmp__7144[0]
      v__7145[1] += tmp__7144[1]
      v__7145[2] += tmp__7144[2]
      Fluid[c].rhoVelocity_t = v__7145
      Fluid[c].rhoEnergy_t += ((-(Fluid[c].rhoEnergyFluxY-Fluid[((c+{0, -1, 0})%Fluid.bounds)].rhoEnergyFluxY))/Grid_yCellWidth)
      Fluid[c].rho_t += ((-(Fluid[c].rhoFluxZ-Fluid[((c+{0, 0, -1})%Fluid.bounds)].rhoFluxZ))/Grid_zCellWidth)
      var tmp__7146 = vs_div_double_3(vs_mul_double_3(vv_sub_double_3(Fluid[c].rhoVelocityFluxZ, Fluid[((c+{0, 0, -1})%Fluid.bounds)].rhoVelocityFluxZ), double((-1))), Grid_zCellWidth)
      var v__7147 = Fluid[c].rhoVelocity_t
      v__7147[0] += tmp__7146[0]
      v__7147[1] += tmp__7146[1]
      v__7147[2] += tmp__7146[2]
      Fluid[c].rhoVelocity_t = v__7147
      Fluid[c].rhoEnergy_t += ((-(Fluid[c].rhoEnergyFluxZ-Fluid[((c+{0, 0, -1})%Fluid.bounds)].rhoEnergyFluxZ))/Grid_zCellWidth)
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
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      var tmp = vs_mul_double_3(Flow_bodyForce, Fluid[c].rho)
      var v = Fluid[c].rhoVelocity_t
      v[0] += tmp[0]
      v[1] += tmp[1]
      v[2] += tmp[2]
      Fluid[c].rhoVelocity_t = v
      Fluid[c].rhoEnergy_t += (Fluid[c].rho*dot_double_3(Flow_bodyForce, Fluid[c].velocity))
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
  reads writes(Fluid.PD)
do
  __demand(__openmp)
  for c in Fluid do
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      var divU = double(0.0)
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
  reads writes(Fluid.dissipationFlux)
do
  __demand(__openmp)
  for c in Fluid do
    if ((not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) or (max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)==1)) then
      var muFace = (double(0.5)*(GetDynamicViscosity(Fluid[c].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)+GetDynamicViscosity(Fluid[((c+{1, 0, 0})%Fluid.bounds)].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)))
      var velocityFace = [double[3]](array(0.0, 0.0, 0.0))
      var velocityX_YFace = double(0.0)
      var velocityX_ZFace = double(0.0)
      var velocityY_YFace = double(0.0)
      var velocityZ_ZFace = double(0.0)
      velocityFace = vs_mul_double_3(vv_add_double_3(Fluid[c].velocity, Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocity), double(0.5))
      velocityX_YFace = (double(0.5)*(Fluid[c].velocityGradientY[0]+Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocityGradientY[0]))
      velocityX_ZFace = (double(0.5)*(Fluid[c].velocityGradientZ[0]+Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocityGradientZ[0]))
      velocityY_YFace = (double(0.5)*(Fluid[c].velocityGradientY[1]+Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocityGradientY[1]))
      velocityZ_ZFace = (double(0.5)*(Fluid[c].velocityGradientZ[2]+Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocityGradientZ[2]))
      var velocityX_XFace = double(0.0)
      var velocityY_XFace = double(0.0)
      var velocityZ_XFace = double(0.0)
      var temperature_XFace = double(0.0)
      velocityX_XFace = (double(0.5)*(Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocity[0]-Fluid[c].velocity[0]))
      velocityY_XFace = (double(0.5)*(Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocity[1]-Fluid[c].velocity[1]))
      velocityZ_XFace = (double(0.5)*(Fluid[((c+{1, 0, 0})%Fluid.bounds)].velocity[2]-Fluid[c].velocity[2]))
      temperature_XFace = (double(0.5)*(Fluid[((c+{1, 0, 0})%Fluid.bounds)].temperature-Fluid[c].temperature))
      velocityX_XFace *= (1/(Grid_xCellWidth*double(0.5)))
      velocityY_XFace *= (1/(Grid_xCellWidth*double(0.5)))
      velocityZ_XFace *= (1/(Grid_xCellWidth*double(0.5)))
      temperature_XFace *= (1/(Grid_xCellWidth*double(0.5)))
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
  reads writes(Fluid.dissipationFlux)
do
  __demand(__openmp)
  for c in Fluid do
    if ((not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)==1)) then
      var muFace = (double(0.5)*(GetDynamicViscosity(Fluid[c].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)+GetDynamicViscosity(Fluid[((c+{0, 1, 0})%Fluid.bounds)].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)))
      var velocityFace = [double[3]](array(0.0, 0.0, 0.0))
      var velocityY_XFace = double(0.0)
      var velocityY_ZFace = double(0.0)
      var velocityX_XFace = double(0.0)
      var velocityZ_ZFace = double(0.0)
      velocityFace = vs_mul_double_3(vv_add_double_3(Fluid[c].velocity, Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocity), double(0.5))
      velocityY_XFace = (double(0.5)*(Fluid[c].velocityGradientX[1]+Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocityGradientX[1]))
      velocityY_ZFace = (double(0.5)*(Fluid[c].velocityGradientZ[1]+Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocityGradientZ[1]))
      velocityX_XFace = (double(0.5)*(Fluid[c].velocityGradientX[0]+Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocityGradientX[0]))
      velocityZ_ZFace = (double(0.5)*(Fluid[c].velocityGradientZ[2]+Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocityGradientZ[2]))
      var velocityX_YFace = double(0.0)
      var velocityY_YFace = double(0.0)
      var velocityZ_YFace = double(0.0)
      var temperature_YFace = double(0.0)
      velocityX_YFace = (double(0.5)*(Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocity[0]-Fluid[c].velocity[0]))
      velocityY_YFace = (double(0.5)*(Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocity[1]-Fluid[c].velocity[1]))
      velocityZ_YFace = (double(0.5)*(Fluid[((c+{0, 1, 0})%Fluid.bounds)].velocity[2]-Fluid[c].velocity[2]))
      temperature_YFace = (double(0.5)*(Fluid[((c+{0, 1, 0})%Fluid.bounds)].temperature-Fluid[c].temperature))
      velocityX_YFace *= (1/(Grid_yCellWidth*double(0.5)))
      velocityY_YFace *= (1/(Grid_yCellWidth*double(0.5)))
      velocityZ_YFace *= (1/(Grid_yCellWidth*double(0.5)))
      temperature_YFace *= (1/(Grid_yCellWidth*double(0.5)))
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
      var muFace = (double(0.5)*(GetDynamicViscosity(Fluid[c].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)+GetDynamicViscosity(Fluid[((c+{0, 0, 1})%Fluid.bounds)].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)))
      var velocityFace = [double[3]](array(0.0, 0.0, 0.0))
      var velocityZ_XFace = double(0.0)
      var velocityZ_YFace = double(0.0)
      var velocityX_XFace = double(0.0)
      var velocityY_YFace = double(0.0)
      velocityFace = vs_mul_double_3(vv_add_double_3(Fluid[c].velocity, Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocity), double(0.5))
      velocityZ_XFace = (double(0.5)*(Fluid[c].velocityGradientX[2]+Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocityGradientX[2]))
      velocityZ_YFace = (double(0.5)*(Fluid[c].velocityGradientY[2]+Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocityGradientY[2]))
      velocityX_XFace = (double(0.5)*(Fluid[c].velocityGradientX[0]+Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocityGradientX[0]))
      velocityY_YFace = (double(0.5)*(Fluid[c].velocityGradientY[1]+Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocityGradientY[1]))
      var velocityX_ZFace = double(0.0)
      var velocityY_ZFace = double(0.0)
      var velocityZ_ZFace = double(0.0)
      var temperature_ZFace = double(0.0)
      velocityX_ZFace = (double(0.5)*(Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocity[0]-Fluid[c].velocity[0]))
      velocityY_ZFace = (double(0.5)*(Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocity[1]-Fluid[c].velocity[1]))
      velocityZ_ZFace = (double(0.5)*(Fluid[((c+{0, 0, 1})%Fluid.bounds)].velocity[2]-Fluid[c].velocity[2]))
      temperature_ZFace = (double(0.5)*(Fluid[((c+{0, 0, 1})%Fluid.bounds)].temperature-Fluid[c].temperature))
      velocityX_ZFace *= (1/(Grid_zCellWidth*double(0.5)))
      velocityY_ZFace *= (1/(Grid_zCellWidth*double(0.5)))
      velocityZ_ZFace *= (1/(Grid_zCellWidth*double(0.5)))
      temperature_ZFace *= (1/(Grid_zCellWidth*double(0.5)))
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

__demand(__parallel)
task CalculateAveragePD(Fluid : region(ispace(int3d), Fluid_columns),
                        Grid_cellVolume : double,
                        Grid_xBnum : int32, Grid_xNum : int32,
                        Grid_yBnum : int32, Grid_yNum : int32,
                        Grid_zBnum : int32, Grid_zNum : int32) : double
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

__demand(__parallel)
task CalculateAverageDissipation(Fluid : region(ispace(int3d), Fluid_columns),
                                 Grid_cellVolume : double,
                                 Grid_xBnum : int32, Grid_xNum : int32,
                                 Grid_yBnum : int32, Grid_yNum : int32,
                                 Grid_zBnum : int32, Grid_zNum : int32) : double
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

__demand(__parallel)
task CalculateAverageK(Fluid : region(ispace(int3d), Fluid_columns),
                       Grid_cellVolume : double,
                       Grid_xBnum : int32, Grid_xNum : int32,
                       Grid_yBnum : int32, Grid_yNum : int32,
                       Grid_zBnum : int32, Grid_zNum : int32) : double
where
  reads(Fluid.{rho, velocity})
do
  var acc = 0.0
  __demand(__openmp)
  for c in Fluid do
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      acc += (((double(0.5)*Fluid[c].rho)*dot_double_3(Fluid[c].velocity, Fluid[c].velocity))*Grid_cellVolume)
    end
  end
  return acc
end

__demand(__parallel)
task Flow_AddTurbulentSource(Fluid : region(ispace(int3d), Fluid_columns),
                             Flow_averageDissipation : double,
                             Flow_averageK : double,
                             Flow_averagePD : double,
                             Grid_cellVolume : double,
                             Grid_xBnum : int32, Grid_xNum : int32,
                             Grid_yBnum : int32, Grid_yNum : int32,
                             Grid_zBnum : int32, Grid_zNum : int32) : double
where
  reads(Fluid.{rho, velocity}),
  reads writes(Fluid.{rhoVelocity_t, rhoEnergy_t})
do
  var acc = 0.0
  __demand(__openmp)
  for c in Fluid do
    if (not ((((((max(int32((uint64(Grid_xBnum)-int3d(c).x)), 0)>0) or (max(int32((int3d(c).x-uint64(((Grid_xNum+Grid_xBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_yBnum)-int3d(c).y)), 0)>0)) or (max(int32((int3d(c).y-uint64(((Grid_yNum+Grid_yBnum)-1)))), 0)>0)) or (max(int32((uint64(Grid_zBnum)-int3d(c).z)), 0)>0)) or (max(int32((int3d(c).z-uint64(((Grid_zNum+Grid_zBnum)-1)))), 0)>0))) then
      var W = double(0.0)
      var A = double(0.0)
      var G = double(0.0)
      var t_o = double(0.0)
      var K_o = double(0.0)
      var force = [double[3]](array(0.0, 0.0, 0.0))
      W = (Flow_averagePD+Flow_averageDissipation)
      G = 300.0
      t_o = double(3.00889e-06)
      K_o = double(66.27348)
      A = (((-W)-((G*(Flow_averageK-K_o))/t_o))/(2.0*Flow_averageK))
      force = vs_mul_double_3(Fluid[c].velocity, (Fluid[c].rho*A))
      var tmp = force
      var v = Fluid[c].rhoVelocity_t
      v[0] += tmp[0]
      v[1] += tmp[1]
      v[2] += tmp[2]
      Fluid[c].rhoVelocity_t = v
      Fluid[c].rhoEnergy_t += dot_double_3(force, Fluid[c].velocity)
      acc += (dot_double_3(force, Fluid[c].velocity)*Grid_cellVolume)
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

__demand(__inline)
task locate(pos : double[3],
            BC_xBCPeriodic : bool, BC_yBCPeriodic : bool, BC_zBCPeriodic : bool,
            Grid_xBnum : int32, Grid_xNum : int32, Grid_xOrigin : double, Grid_xWidth : double,
            Grid_yBnum : int32, Grid_yNum : int32, Grid_yOrigin : double, Grid_yWidth : double,
            Grid_zBnum : int32, Grid_zNum : int32, Grid_zOrigin : double, Grid_zWidth : double) : int3d
  var xcw = (Grid_xWidth/double(Grid_xNum))
  var xro = (Grid_xOrigin-(double(Grid_xBnum)*xcw))
  var xpos = ((pos[0]-xro)/xcw)
  var xrnum = (Grid_xNum+(2*Grid_xBnum))
  var xidx : uint64
  if BC_xBCPeriodic then
    xidx = (uint64((fmod(xpos, double(xrnum))+double(xrnum)))%uint64(xrnum))
  else
    xidx = uint64(max(0.0, min(double((xrnum-1)), xpos)))
  end
  var ycw = (Grid_yWidth/double(Grid_yNum))
  var yro = (Grid_yOrigin-(double(Grid_yBnum)*ycw))
  var ypos = ((pos[1]-yro)/ycw)
  var yrnum = (Grid_yNum+(2*Grid_yBnum))
  var yidx : uint64
  if BC_yBCPeriodic then
    yidx = (uint64((fmod(ypos, double(yrnum))+double(yrnum)))%uint64(yrnum))
  else
    yidx = uint64(max(0.0, min(double((yrnum-1)), ypos)))
  end
  var zcw = (Grid_zWidth/double(Grid_zNum))
  var zro = (Grid_zOrigin-(double(Grid_zBnum)*zcw))
  var zpos = ((pos[2]-zro)/zcw)
  var zrnum = (Grid_zNum+(2*Grid_zBnum))
  var zidx : uint64
  if BC_zBCPeriodic then
    zidx = (uint64((fmod(zpos, double(zrnum))+double(zrnum)))%uint64(zrnum))
  else
    zidx = uint64(max(0.0, min(double((zrnum-1)), zpos)))
  end
  return int3d({xidx, yidx, zidx})
end

__demand(__cuda)
task Particles_LocateInCells(particles : region(ispace(int1d), particles_columns),
                             BC_xBCPeriodic : bool, BC_yBCPeriodic : bool, BC_zBCPeriodic : bool,
                             Grid_xBnum : int32, Grid_xNum : int32, Grid_xOrigin : double, Grid_xWidth : double,
                             Grid_yBnum : int32, Grid_yNum : int32, Grid_yOrigin : double, Grid_yWidth : double,
                             Grid_zBnum : int32, Grid_zNum : int32, Grid_zOrigin : double, Grid_zWidth : double)
where
  reads(particles.{position, __valid}),
  writes(particles.cell)
do
  __demand(__openmp)
  for p in particles do
    if particles[p].__valid then
      particles[p].cell = locate(particles[p].position, BC_xBCPeriodic, BC_yBCPeriodic, BC_zBCPeriodic, Grid_xBnum, Grid_xNum, Grid_xOrigin, Grid_xWidth, Grid_yBnum, Grid_yNum, Grid_yOrigin, Grid_yWidth, Grid_zBnum, Grid_zNum, Grid_zOrigin, Grid_zWidth)
    end
  end
end

__demand(__inline)
task Fluid_elemColor(idx : int3d,
                     xNum  : int32, yNum  : int32, zNum  : int32,
                     xBnum : int32, yBnum : int32, zBnum : int32,
                     NX_   : int32, NY_   : int32, NZ_   : int32) : int3d
  idx.x = min(max(idx.x, xBnum), ((xNum+xBnum)-1))
  idx.y = min(max(idx.y, yBnum), ((yNum+yBnum)-1))
  idx.z = min(max(idx.z, zBnum), ((zNum+zBnum)-1))
  return int3d({((idx.x-xBnum)/(xNum/NX_)), ((idx.y-yBnum)/(yNum/NY_)), ((idx.z-zBnum)/(zNum/NZ_))})
end

terra particles_pushElement(dst : &opaque,idx : int32,src : particles_columns) : {}
  var ptr = [&int8](dst) + idx * 376
  C.memcpy([&opaque](ptr), [&opaque](&src), [uint64](376))
end

terra particles_getBasePointer(pr : regentlib.c.legion_physical_region_t,fid : uint32,runtime : regentlib.c.legion_runtime_t) : &opaque
  var acc = regentlib.c.legion_physical_region_get_field_accessor_array_1d(pr, fid)
  var lr = regentlib.c.legion_physical_region_get_logical_region(pr)
  var domain = regentlib.c.legion_index_space_get_domain(runtime, lr.index_space)
  var rect = regentlib.c.legion_domain_get_rect_1d(domain)
  var subrect : regentlib.c.legion_rect_1d_t
  var offsets : regentlib.c.legion_byte_offset_t[1]
  var p = regentlib.c.legion_accessor_array_1d_raw_rect_ptr(acc, rect, &subrect, &[&regentlib.c.legion_byte_offset_t](offsets)[0])
  regentlib.c.legion_accessor_array_1d_destroy(acc)
  return p
end

terra particles_getOffset() : int64
  var x : particles_columns
  return [&int8](&x.__valid) - [&int8](&x)
end

task particles_pushAll(partColor : int3d,
                       r : region(ispace(int1d), particles_columns),
                       q0 : region(ispace(int1d), int8[376]),
                       q1 : region(ispace(int1d), int8[376]),
                       q2 : region(ispace(int1d), int8[376]),
                       q3 : region(ispace(int1d), int8[376]),
                       q4 : region(ispace(int1d), int8[376]),
                       q5 : region(ispace(int1d), int8[376]),
                       q6 : region(ispace(int1d), int8[376]),
                       q7 : region(ispace(int1d), int8[376]),
                       q8 : region(ispace(int1d), int8[376]),
                       q9 : region(ispace(int1d), int8[376]),
                       q10 : region(ispace(int1d), int8[376]),
                       q11 : region(ispace(int1d), int8[376]),
                       q12 : region(ispace(int1d), int8[376]),
                       q13 : region(ispace(int1d), int8[376]),
                       q14 : region(ispace(int1d), int8[376]),
                       q15 : region(ispace(int1d), int8[376]),
                       q16 : region(ispace(int1d), int8[376]),
                       q17 : region(ispace(int1d), int8[376]),
                       q18 : region(ispace(int1d), int8[376]),
                       q19 : region(ispace(int1d), int8[376]),
                       q20 : region(ispace(int1d), int8[376]),
                       q21 : region(ispace(int1d), int8[376]),
                       q22 : region(ispace(int1d), int8[376]),
                       q23 : region(ispace(int1d), int8[376]),
                       q24 : region(ispace(int1d), int8[376]),
                       q25 : region(ispace(int1d), int8[376]),
                       rngXNum : int32,  rngYNum  : int32, rngZNum  : int32,
                       rngXbnum : int32, rngYbnum : int32, rngZbnum : int32,
                       NX_ : int32, NY_ : int32, NZ_ : int32)
where
  reads(r),
  writes(r.__valid),
  reads writes(q0),
  reads writes(q1),
  reads writes(q2),
  reads writes(q3),
  reads writes(q4),
  reads writes(q5),
  reads writes(q6),
  reads writes(q7),
  reads writes(q8),
  reads writes(q9),
  reads writes(q10),
  reads writes(q11),
  reads writes(q12),
  reads writes(q13),
  reads writes(q14),
  reads writes(q15),
  reads writes(q16),
  reads writes(q17),
  reads writes(q18),
  reads writes(q19),
  reads writes(q20),
  reads writes(q21),
  reads writes(q22),
  reads writes(q23),
  reads writes(q24),
  reads writes(q25)
do
  for qPtr in q0 do
    q0[qPtr][368LL] = int8(false)
  end
  var qBasePtr0 = particles_getBasePointer(__physical(q0)[0], __fields(q0)[0], __runtime())
  for qPtr in q1 do
    q1[qPtr][368LL] = int8(false)
  end
  var qBasePtr1 = particles_getBasePointer(__physical(q1)[0], __fields(q1)[0], __runtime())
  for qPtr in q2 do
    q2[qPtr][368LL] = int8(false)
  end
  var qBasePtr2 = particles_getBasePointer(__physical(q2)[0], __fields(q2)[0], __runtime())
  for qPtr in q3 do
    q3[qPtr][368LL] = int8(false)
  end
  var qBasePtr3 = particles_getBasePointer(__physical(q3)[0], __fields(q3)[0], __runtime())
  for qPtr in q4 do
    q4[qPtr][368LL] = int8(false)
  end
  var qBasePtr4 = particles_getBasePointer(__physical(q4)[0], __fields(q4)[0], __runtime())
  for qPtr in q5 do
    q5[qPtr][368LL] = int8(false)
  end
  var qBasePtr5 = particles_getBasePointer(__physical(q5)[0], __fields(q5)[0], __runtime())
  for qPtr in q6 do
    q6[qPtr][368LL] = int8(false)
  end
  var qBasePtr6 = particles_getBasePointer(__physical(q6)[0], __fields(q6)[0], __runtime())
  for qPtr in q7 do
    q7[qPtr][368LL] = int8(false)
  end
  var qBasePtr7 = particles_getBasePointer(__physical(q7)[0], __fields(q7)[0], __runtime())
  for qPtr in q8 do
    q8[qPtr][368LL] = int8(false)
  end
  var qBasePtr8 = particles_getBasePointer(__physical(q8)[0], __fields(q8)[0], __runtime())
  for qPtr in q9 do
    q9[qPtr][368LL] = int8(false)
  end
  var qBasePtr9 = particles_getBasePointer(__physical(q9)[0], __fields(q9)[0], __runtime())
  for qPtr in q10 do
    q10[qPtr][368LL] = int8(false)
  end
  var qBasePtr10 = particles_getBasePointer(__physical(q10)[0], __fields(q10)[0], __runtime())
  for qPtr in q11 do
    q11[qPtr][368LL] = int8(false)
  end
  var qBasePtr11 = particles_getBasePointer(__physical(q11)[0], __fields(q11)[0], __runtime())
  for qPtr in q12 do
    q12[qPtr][368LL] = int8(false)
  end
  var qBasePtr12 = particles_getBasePointer(__physical(q12)[0], __fields(q12)[0], __runtime())
  for qPtr in q13 do
    q13[qPtr][368LL] = int8(false)
  end
  var qBasePtr13 = particles_getBasePointer(__physical(q13)[0], __fields(q13)[0], __runtime())
  for qPtr in q14 do
    q14[qPtr][368LL] = int8(false)
  end
  var qBasePtr14 = particles_getBasePointer(__physical(q14)[0], __fields(q14)[0], __runtime())
  for qPtr in q15 do
    q15[qPtr][368LL] = int8(false)
  end
  var qBasePtr15 = particles_getBasePointer(__physical(q15)[0], __fields(q15)[0], __runtime())
  for qPtr in q16 do
    q16[qPtr][368LL] = int8(false)
  end
  var qBasePtr16 = particles_getBasePointer(__physical(q16)[0], __fields(q16)[0], __runtime())
  for qPtr in q17 do
    q17[qPtr][368LL] = int8(false)
  end
  var qBasePtr17 = particles_getBasePointer(__physical(q17)[0], __fields(q17)[0], __runtime())
  for qPtr in q18 do
    q18[qPtr][368LL] = int8(false)
  end
  var qBasePtr18 = particles_getBasePointer(__physical(q18)[0], __fields(q18)[0], __runtime())
  for qPtr in q19 do
    q19[qPtr][368LL] = int8(false)
  end
  var qBasePtr19 = particles_getBasePointer(__physical(q19)[0], __fields(q19)[0], __runtime())
  for qPtr in q20 do
    q20[qPtr][368LL] = int8(false)
  end
  var qBasePtr20 = particles_getBasePointer(__physical(q20)[0], __fields(q20)[0], __runtime())
  for qPtr in q21 do
    q21[qPtr][368LL] = int8(false)
  end
  var qBasePtr21 = particles_getBasePointer(__physical(q21)[0], __fields(q21)[0], __runtime())
  for qPtr in q22 do
    q22[qPtr][368LL] = int8(false)
  end
  var qBasePtr22 = particles_getBasePointer(__physical(q22)[0], __fields(q22)[0], __runtime())
  for qPtr in q23 do
    q23[qPtr][368LL] = int8(false)
  end
  var qBasePtr23 = particles_getBasePointer(__physical(q23)[0], __fields(q23)[0], __runtime())
  for qPtr in q24 do
    q24[qPtr][368LL] = int8(false)
  end
  var qBasePtr24 = particles_getBasePointer(__physical(q24)[0], __fields(q24)[0], __runtime())
  for qPtr in q25 do
    q25[qPtr][368LL] = int8(false)
  end
  var qBasePtr25 = particles_getBasePointer(__physical(q25)[0], __fields(q25)[0], __runtime())
  for rPtr in r do
    if rPtr.__valid then
      var elemColor = Fluid_elemColor(rPtr.cell, rngXNum, rngYNum, rngZNum, rngXbnum, rngYbnum, rngZbnum, NX_, NY_, NZ_)
      if (elemColor~=partColor) then
        do
          var colorOff = int3d({0, 0, 1})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q0 do
              if (not bool(q0[qPtr][368LL])) then
                particles_pushElement(qBasePtr0, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q0[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({0, 0, -1})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q1 do
              if (not bool(q1[qPtr][368LL])) then
                particles_pushElement(qBasePtr1, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q1[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({0, 1, 0})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q2 do
              if (not bool(q2[qPtr][368LL])) then
                particles_pushElement(qBasePtr2, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q2[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({0, 1, 1})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q3 do
              if (not bool(q3[qPtr][368LL])) then
                particles_pushElement(qBasePtr3, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q3[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({0, 1, -1})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q4 do
              if (not bool(q4[qPtr][368LL])) then
                particles_pushElement(qBasePtr4, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q4[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({0, -1, 0})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q5 do
              if (not bool(q5[qPtr][368LL])) then
                particles_pushElement(qBasePtr5, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q5[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({0, -1, 1})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q6 do
              if (not bool(q6[qPtr][368LL])) then
                particles_pushElement(qBasePtr6, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q6[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({0, -1, -1})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q7 do
              if (not bool(q7[qPtr][368LL])) then
                particles_pushElement(qBasePtr7, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q7[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({1, 0, 0})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q8 do
              if (not bool(q8[qPtr][368LL])) then
                particles_pushElement(qBasePtr8, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q8[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({1, 0, 1})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q9 do
              if (not bool(q9[qPtr][368LL])) then
                particles_pushElement(qBasePtr9, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q9[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({1, 0, -1})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q10 do
              if (not bool(q10[qPtr][368LL])) then
                particles_pushElement(qBasePtr10, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q10[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({1, 1, 0})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q11 do
              if (not bool(q11[qPtr][368LL])) then
                particles_pushElement(qBasePtr11, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q11[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({1, 1, 1})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q12 do
              if (not bool(q12[qPtr][368LL])) then
                particles_pushElement(qBasePtr12, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q12[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({1, 1, -1})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q13 do
              if (not bool(q13[qPtr][368LL])) then
                particles_pushElement(qBasePtr13, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q13[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({1, -1, 0})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q14 do
              if (not bool(q14[qPtr][368LL])) then
                particles_pushElement(qBasePtr14, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q14[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({1, -1, 1})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q15 do
              if (not bool(q15[qPtr][368LL])) then
                particles_pushElement(qBasePtr15, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q15[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({1, -1, -1})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q16 do
              if (not bool(q16[qPtr][368LL])) then
                particles_pushElement(qBasePtr16, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q16[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({-1, 0, 0})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q17 do
              if (not bool(q17[qPtr][368LL])) then
                particles_pushElement(qBasePtr17, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q17[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({-1, 0, 1})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q18 do
              if (not bool(q18[qPtr][368LL])) then
                particles_pushElement(qBasePtr18, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q18[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({-1, 0, -1})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q19 do
              if (not bool(q19[qPtr][368LL])) then
                particles_pushElement(qBasePtr19, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q19[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({-1, 1, 0})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q20 do
              if (not bool(q20[qPtr][368LL])) then
                particles_pushElement(qBasePtr20, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q20[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({-1, 1, 1})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q21 do
              if (not bool(q21[qPtr][368LL])) then
                particles_pushElement(qBasePtr21, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q21[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({-1, 1, -1})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q22 do
              if (not bool(q22[qPtr][368LL])) then
                particles_pushElement(qBasePtr22, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q22[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({-1, -1, 0})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q23 do
              if (not bool(q23[qPtr][368LL])) then
                particles_pushElement(qBasePtr23, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q23[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({-1, -1, 1})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q24 do
              if (not bool(q24[qPtr][368LL])) then
                particles_pushElement(qBasePtr24, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q24[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        do
          var colorOff = int3d({-1, -1, -1})
          if (rPtr.__valid and (elemColor==(((partColor+colorOff)+{NX_, NY_, NZ_})%{NX_, NY_, NZ_}))) then
            var idx = 0
            for qPtr in q25 do
              if (not bool(q25[qPtr][368LL])) then
                particles_pushElement(qBasePtr25, idx, r[rPtr])
                rPtr.__valid = false
                regentlib.assert(bool(q25[qPtr][368LL]), "Element did not get copied properly")
                break
              end
              idx += 1
            end
            regentlib.assert((not rPtr.__valid), "Transfer queue ran out of space")
          end
        end
        regentlib.assert((not rPtr.__valid), "Element moved past predicted stencil")
      end
    end
  end
end

terra particles_pullElement(src : &int8) : particles_columns
  var dst : particles_columns
  C.memcpy([&opaque](&dst), [&opaque](src), [uint64](376))
  return dst
end

task particles_pullAll(color : int3d,
                       r : region(ispace(int1d), particles_columns),
                       q0 : region(ispace(int1d), int8[376]),
                       q1 : region(ispace(int1d), int8[376]),
                       q2 : region(ispace(int1d), int8[376]),
                       q3 : region(ispace(int1d), int8[376]),
                       q4 : region(ispace(int1d), int8[376]),
                       q5 : region(ispace(int1d), int8[376]),
                       q6 : region(ispace(int1d), int8[376]),
                       q7 : region(ispace(int1d), int8[376]),
                       q8 : region(ispace(int1d), int8[376]),
                       q9 : region(ispace(int1d), int8[376]),
                       q10 : region(ispace(int1d), int8[376]),
                       q11 : region(ispace(int1d), int8[376]),
                       q12 : region(ispace(int1d), int8[376]),
                       q13 : region(ispace(int1d), int8[376]),
                       q14 : region(ispace(int1d), int8[376]),
                       q15 : region(ispace(int1d), int8[376]),
                       q16 : region(ispace(int1d), int8[376]),
                       q17 : region(ispace(int1d), int8[376]),
                       q18 : region(ispace(int1d), int8[376]),
                       q19 : region(ispace(int1d), int8[376]),
                       q20 : region(ispace(int1d), int8[376]),
                       q21 : region(ispace(int1d), int8[376]),
                       q22 : region(ispace(int1d), int8[376]),
                       q23 : region(ispace(int1d), int8[376]),
                       q24 : region(ispace(int1d), int8[376]),
                       q25 : region(ispace(int1d), int8[376]))
where
  reads writes(r),
  reads writes(q0),
  reads writes(q1),
  reads writes(q2),
  reads writes(q3),
  reads writes(q4),
  reads writes(q5),
  reads writes(q6),
  reads writes(q7),
  reads writes(q8),
  reads writes(q9),
  reads writes(q10),
  reads writes(q11),
  reads writes(q12),
  reads writes(q13),
  reads writes(q14),
  reads writes(q15),
  reads writes(q16),
  reads writes(q17),
  reads writes(q18),
  reads writes(q19),
  reads writes(q20),
  reads writes(q21),
  reads writes(q22),
  reads writes(q23),
  reads writes(q24),
  reads writes(q25)
do
  for qPtr in q0 do
    if bool(q0[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q0[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q1 do
    if bool(q1[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q1[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q2 do
    if bool(q2[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q2[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q3 do
    if bool(q3[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q3[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q4 do
    if bool(q4[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q4[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q5 do
    if bool(q5[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q5[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q6 do
    if bool(q6[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q6[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q7 do
    if bool(q7[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q7[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q8 do
    if bool(q8[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q8[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q9 do
    if bool(q9[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q9[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q10 do
    if bool(q10[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q10[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q11 do
    if bool(q11[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q11[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q12 do
    if bool(q12[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q12[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q13 do
    if bool(q13[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q13[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q14 do
    if bool(q14[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q14[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q15 do
    if bool(q15[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q15[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q16 do
    if bool(q16[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q16[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q17 do
    if bool(q17[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q17[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q18 do
    if bool(q18[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q18[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q19 do
    if bool(q19[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q19[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q20 do
    if bool(q20[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q20[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q21 do
    if bool(q21[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q21[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q22 do
    if bool(q22[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q22[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q23 do
    if bool(q23[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q23[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q24 do
    if bool(q24[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q24[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
  for qPtr in q25 do
    if bool(q25[qPtr][368LL]) then
      var copied = false
      for rPtr in r do
        if (not rPtr.__valid) then
          r[rPtr] = particles_pullElement([&int8](q25[qPtr]))
          copied = true
          regentlib.assert(r[rPtr].__valid, "Pulled particle was not copied correctly")
          break
        end
      end
      regentlib.assert(copied, "Ran out of space on sub-partition")
    end
  end
end

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
                                  Grid_zCellWidth : double, Grid_zRealOrigin : double) : double[3]
  var dX = fmod((((xyz[0]-Grid_xRealOrigin)/Grid_xCellWidth)+double(0.5)), 1.0)
  var dY = fmod((((xyz[1]-Grid_yRealOrigin)/Grid_yCellWidth)+double(0.5)), 1.0)
  var dZ = fmod((((xyz[2]-Grid_zRealOrigin)/Grid_zCellWidth)+double(0.5)), 1.0)
  var oneMinusdX = (1.0-dX)
  var oneMinusdY = (1.0-dY)
  var oneMinusdZ = (1.0-dZ)
  var weight00 = vv_add_double_3(vs_mul_double_3(c000, oneMinusdX), vs_mul_double_3(c100, dX))
  var weight10 = vv_add_double_3(vs_mul_double_3(c010, oneMinusdX), vs_mul_double_3(c110, dX))
  var weight01 = vv_add_double_3(vs_mul_double_3(c001, oneMinusdX), vs_mul_double_3(c101, dX))
  var weight11 = vv_add_double_3(vs_mul_double_3(c011, oneMinusdX), vs_mul_double_3(c111, dX))
  var weight0 = vv_add_double_3(vs_mul_double_3(weight00, oneMinusdY), vs_mul_double_3(weight10, dY))
  var weight1 = vv_add_double_3(vs_mul_double_3(weight01, oneMinusdY), vs_mul_double_3(weight11, dY))
  return vv_add_double_3(vs_mul_double_3(weight0, oneMinusdZ), vs_mul_double_3(weight1, dZ))
end

__demand(__inline)
task InterpolateTriVelocity(c : int3d,
                            xyz : double[3],
                            Fluid : region(ispace(int3d), Fluid_columns),
                            Grid_xCellWidth : double, Grid_xRealOrigin : double,
                            Grid_yCellWidth : double, Grid_yRealOrigin : double,
                            Grid_zCellWidth : double, Grid_zRealOrigin : double) : double[3]
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
                              Grid_zCellWidth : double, Grid_zRealOrigin : double) : double
  var dX = fmod((((xyz[0]-Grid_xRealOrigin)/Grid_xCellWidth)+double(0.5)), 1.0)
  var dY = fmod((((xyz[1]-Grid_yRealOrigin)/Grid_yCellWidth)+double(0.5)), 1.0)
  var dZ = fmod((((xyz[2]-Grid_zRealOrigin)/Grid_zCellWidth)+double(0.5)), 1.0)
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
                        Grid_zCellWidth : double, Grid_zRealOrigin : double) : double
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
task Particles_AddFlowCoupling(particles : region(ispace(int1d), particles_columns),
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
  reads(particles.{cell, position, velocity, diameter, density, temperature, __valid}),
  reads writes(particles.{position_t, velocity_t, temperature_t, deltaTemperatureTerm, deltaVelocityOverRelaxationTime})
do
  __demand(__openmp)
  for p in particles do
    if particles[p].__valid then
      var flowVelocity = InterpolateTriVelocity(particles[p].cell, particles[p].position, Fluid, Grid_xCellWidth, Grid_xRealOrigin, Grid_yCellWidth, Grid_yRealOrigin, Grid_zCellWidth, Grid_zRealOrigin)
      var flowTemperature = InterpolateTriTemp(particles[p].cell, particles[p].position, Fluid, Grid_xCellWidth, Grid_xRealOrigin, Grid_yCellWidth, Grid_yRealOrigin, Grid_zCellWidth, Grid_zRealOrigin)
      var flowDynamicViscosity = GetDynamicViscosity(flowTemperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)
      var tmp = particles[p].velocity
      var v = particles[p].position_t
      v[0] += tmp[0]
      v[1] += tmp[1]
      v[2] += tmp[2]
      particles[p].position_t = v
      var particleReynoldsNumber = 0.0
      var relaxationTime = (((particles[p].density*pow(particles[p].diameter, 2.0))/(18.0*flowDynamicViscosity))/(1.0+(double(0.15)*pow(particleReynoldsNumber, double(0.687)))))
      particles[p].deltaVelocityOverRelaxationTime = vs_div_double_3(vv_sub_double_3(flowVelocity, particles[p].velocity), relaxationTime)
      particles[p].deltaTemperatureTerm = (((PI*pow(particles[p].diameter, 2.0))*Particles_convectiveCoeff)*(flowTemperature-particles[p].temperature))
      var tmp__10496 = particles[p].deltaVelocityOverRelaxationTime
      var v__10497 = particles[p].velocity_t
      v__10497[0] += tmp__10496[0]
      v__10497[1] += tmp__10496[1]
      v__10497[2] += tmp__10496[2]
      particles[p].velocity_t = v__10497
      particles[p].temperature_t += (particles[p].deltaTemperatureTerm/((((PI*pow(particles[p].diameter, 3.0))/6.0)*particles[p].density)*Particles_heatCapacity))
    end
  end
end

__demand(__parallel, __cuda)
task Particles_AddBodyForces(particles : region(ispace(int1d), particles_columns),
                             Particles_bodyForce : double[3])
where
  reads(particles.__valid),
  reads writes(particles.velocity_t)
do
  __demand(__openmp)
  for p in particles do
    if particles[p].__valid then
      var tmp = Particles_bodyForce
      var v = particles[p].velocity_t
      v[0] += tmp[0]
      v[1] += tmp[1]
      v[2] += tmp[2]
      particles[p].velocity_t = v
    end
  end
end

__demand(__parallel, __cuda)
task Radiation_ClearAccumulators(Radiation : region(ispace(int3d), Radiation_columns))
where
  reads writes(Radiation.{acc_d2, acc_d2t4})
do
  __demand(__openmp)
  for c in Radiation do
    Radiation[c].acc_d2 = 0.0
    Radiation[c].acc_d2t4 = 0.0
  end
end

__demand(__cuda)
task Radiation_AccumulateParticleValues(particles : region(ispace(int1d), particles_columns),
                                        Fluid : region(ispace(int3d), Fluid_columns),
                                        Radiation : region(ispace(int3d), Radiation_columns))
where
  reads(Fluid.to_Radiation),
  reads(particles.{cell, diameter, temperature, __valid}),
  reads writes(Radiation.{acc_d2, acc_d2t4})
do
  __demand(__openmp)
  for p in particles do
    if particles[p].__valid then
      Radiation[Fluid[particles[p].cell].to_Radiation].acc_d2 += pow(particles[p].diameter, 2.0)
      Radiation[Fluid[particles[p].cell].to_Radiation].acc_d2t4 += (pow(particles[p].diameter, 2.0)*pow(particles[p].temperature, 4.0))
    end
  end
end

__demand(__parallel, __cuda)
task Radiation_UpdateFieldValues(Radiation : region(ispace(int3d), Radiation_columns),
                                 Radiation_cellVolume : double,
                                 Radiation_qa : double,
                                 Radiation_qs : double)
where
  reads writes(Radiation.{Ib, sigma}),
  reads(Radiation.{acc_d2, acc_d2t4})
do
  __demand(__openmp)
  for c in Radiation do
    Radiation[c].sigma = (((Radiation[c].acc_d2*PI)*(Radiation_qa+Radiation_qs))/(4.0*Radiation_cellVolume))
    if (Radiation[c].acc_d2==0.0) then
      Radiation[c].Ib = 0.0
    else
      Radiation[c].Ib = ((double(5.67e-08)*Radiation[c].acc_d2t4)/(PI*Radiation[c].acc_d2))
    end
  end
end

__demand(__cuda)
task Particles_AbsorbRadiation(particles : region(ispace(int1d), particles_columns),
                               Fluid : region(ispace(int3d), Fluid_columns),
                               Radiation : region(ispace(int3d), Radiation_columns),
                               Particles_heatCapacity : double,
                               Radiation_qa : double)
where
  reads(Fluid.to_Radiation),
  reads(Radiation.G),
  reads(particles.{cell, density, diameter, temperature, __valid}),
  reads writes(particles.temperature_t)
do
  __demand(__openmp)
  for p in particles do
    if particles[p].__valid then
      var t4 = pow(particles[p].temperature, 4.0)
      var alpha = ((((PI*Radiation_qa)*pow(particles[p].diameter, 2.0))*(Radiation[Fluid[particles[p].cell].to_Radiation].G-((4.0*double(5.67e-08))*t4)))/4.0)
      particles[p].temperature_t += (alpha/((((PI*pow(particles[p].diameter, 3.0))/6.0)*particles[p].density)*Particles_heatCapacity))
    end
  end
end

__demand(__parallel, __cuda)
task Flow_AddParticlesCoupling(particles : region(ispace(int1d), particles_columns),
                               Fluid : region(ispace(int3d), Fluid_columns),
                               Grid_cellVolume : double)
where
  reads writes(Fluid.{rhoVelocity_t, rhoEnergy_t}),
  reads(particles.{cell, diameter, density, deltaTemperatureTerm, deltaVelocityOverRelaxationTime, __valid})
do
  __demand(__openmp)
  for p in particles do
    if particles[p].__valid then
      var tmp = vs_div_double_3(vs_mul_double_3(particles[p].deltaVelocityOverRelaxationTime, (-(((PI*pow(particles[p].diameter, 3.0))/6.0)*particles[p].density))), Grid_cellVolume)
      var v = Fluid[particles[p].cell].rhoVelocity_t
      v[0] += tmp[0]
      v[1] += tmp[1]
      v[2] += tmp[2]
      Fluid[particles[p].cell].rhoVelocity_t = v
      Fluid[particles[p].cell].rhoEnergy_t += ((-particles[p].deltaTemperatureTerm)/Grid_cellVolume)
    end
  end
end

__demand(__parallel, __cuda)
task Flow_UpdateVars(Fluid : region(ispace(int3d), Fluid_columns),
                     Integrator_deltaTime : double,
                     Integrator_stage : int32)
where
  reads(Fluid.{rhoEnergy_old, rhoEnergy_new}),
  reads(Fluid.{rhoVelocity_old, rhoVelocity_new}),
  reads(Fluid.{rho_old, rho_new}),
  reads(Fluid.{rho_t, rhoVelocity_t, rhoEnergy_t}),
  reads writes(Fluid.{rho, rho_new}),
  reads writes(Fluid.{rhoEnergy, rhoEnergy_new}),
  reads writes(Fluid.{rhoVelocity, rhoVelocity_new})
do
  __demand(__openmp)
  for c in Fluid do
    var deltaTime = Integrator_deltaTime
    if (Integrator_stage==1) then
      Fluid[c].rho_new += (((1.0/6.0)*deltaTime)*Fluid[c].rho_t)
      Fluid[c].rho = (Fluid[c].rho_old+((double(0.5)*deltaTime)*Fluid[c].rho_t))
      var tmp = vs_mul_double_3(Fluid[c].rhoVelocity_t, ((1.0/6.0)*deltaTime))
      var v = Fluid[c].rhoVelocity_new
      v[0] += tmp[0]
      v[1] += tmp[1]
      v[2] += tmp[2]
      Fluid[c].rhoVelocity_new = v
      Fluid[c].rhoVelocity = vv_add_double_3(Fluid[c].rhoVelocity_old, vs_mul_double_3(Fluid[c].rhoVelocity_t, (double(0.5)*deltaTime)))
      Fluid[c].rhoEnergy_new += (((1.0/6.0)*deltaTime)*Fluid[c].rhoEnergy_t)
      Fluid[c].rhoEnergy = (Fluid[c].rhoEnergy_old+((double(0.5)*deltaTime)*Fluid[c].rhoEnergy_t))
    else
      if (Integrator_stage==2) then
        Fluid[c].rho_new += (((1.0/3.0)*deltaTime)*Fluid[c].rho_t)
        Fluid[c].rho = (Fluid[c].rho_old+((double(0.5)*deltaTime)*Fluid[c].rho_t))
        var tmp = vs_mul_double_3(Fluid[c].rhoVelocity_t, ((1.0/3.0)*deltaTime))
        var v = Fluid[c].rhoVelocity_new
        v[0] += tmp[0]
        v[1] += tmp[1]
        v[2] += tmp[2]
        Fluid[c].rhoVelocity_new = v
        Fluid[c].rhoVelocity = vv_add_double_3(Fluid[c].rhoVelocity_old, vs_mul_double_3(Fluid[c].rhoVelocity_t, (double(0.5)*deltaTime)))
        Fluid[c].rhoEnergy_new += (((1.0/3.0)*deltaTime)*Fluid[c].rhoEnergy_t)
        Fluid[c].rhoEnergy = (Fluid[c].rhoEnergy_old+((double(0.5)*deltaTime)*Fluid[c].rhoEnergy_t))
      else
        if (Integrator_stage==3) then
          Fluid[c].rho_new += (((1.0/3.0)*deltaTime)*Fluid[c].rho_t)
          Fluid[c].rho = (Fluid[c].rho_old+((1.0*deltaTime)*Fluid[c].rho_t))
          var tmp = vs_mul_double_3(Fluid[c].rhoVelocity_t, ((1.0/3.0)*deltaTime))
          var v = Fluid[c].rhoVelocity_new
          v[0] += tmp[0]
          v[1] += tmp[1]
          v[2] += tmp[2]
          Fluid[c].rhoVelocity_new = v
          Fluid[c].rhoVelocity = vv_add_double_3(Fluid[c].rhoVelocity_old, vs_mul_double_3(Fluid[c].rhoVelocity_t, (1.0*deltaTime)))
          Fluid[c].rhoEnergy_new += (((1.0/3.0)*deltaTime)*Fluid[c].rhoEnergy_t)
          Fluid[c].rhoEnergy = (Fluid[c].rhoEnergy_old+((1.0*deltaTime)*Fluid[c].rhoEnergy_t))
        else
          Fluid[c].rho = (Fluid[c].rho_new+(((1.0/6.0)*deltaTime)*Fluid[c].rho_t))
          Fluid[c].rhoVelocity = vv_add_double_3(Fluid[c].rhoVelocity_new, vs_mul_double_3(Fluid[c].rhoVelocity_t, ((1.0/6.0)*deltaTime)))
          Fluid[c].rhoEnergy = (Fluid[c].rhoEnergy_new+(((1.0/6.0)*deltaTime)*Fluid[c].rhoEnergy_t))
        end
      end
    end
  end
end

__demand(__parallel, __cuda)
task Particles_UpdateVars(particles : region(ispace(int1d), particles_columns),
                          Integrator_deltaTime : double,
                          Integrator_stage : int32)
where
  reads(particles.{position_old, velocity_old, temperature_old}),
  reads(particles.{position_t, velocity_t, temperature_t, __valid}),
  reads writes(particles.{position, position_new}),
  reads writes(particles.{temperature, temperature_new}),
  reads writes(particles.{velocity, velocity_new})
do
  __demand(__openmp)
  for p in particles do
    if particles[p].__valid then
      var deltaTime = Integrator_deltaTime
      if (Integrator_stage==1) then
        var tmp = vs_mul_double_3(particles[p].position_t, ((1.0/6.0)*deltaTime))
        var v = particles[p].position_new
        v[0] += tmp[0]
        v[1] += tmp[1]
        v[2] += tmp[2]
        particles[p].position_new = v
        particles[p].position = vv_add_double_3(particles[p].position_old, vs_mul_double_3(particles[p].position_t, (double(0.5)*deltaTime)))
        var tmp__11020 = vs_mul_double_3(particles[p].velocity_t, ((1.0/6.0)*deltaTime))
        var v__11021 = particles[p].velocity_new
        v__11021[0] += tmp__11020[0]
        v__11021[1] += tmp__11020[1]
        v__11021[2] += tmp__11020[2]
        particles[p].velocity_new = v__11021
        particles[p].velocity = vv_add_double_3(particles[p].velocity_old, vs_mul_double_3(particles[p].velocity_t, (double(0.5)*deltaTime)))
        particles[p].temperature_new += (((1.0/6.0)*deltaTime)*particles[p].temperature_t)
        particles[p].temperature = (particles[p].temperature_old+((double(0.5)*deltaTime)*particles[p].temperature_t))
      else
        if (Integrator_stage==2) then
          var tmp = vs_mul_double_3(particles[p].position_t, ((1.0/3.0)*deltaTime))
          var v = particles[p].position_new
          v[0] += tmp[0]
          v[1] += tmp[1]
          v[2] += tmp[2]
          particles[p].position_new = v
          particles[p].position = vv_add_double_3(particles[p].position_old, vs_mul_double_3(particles[p].position_t, (double(0.5)*deltaTime)))
          var tmp__11024 = vs_mul_double_3(particles[p].velocity_t, ((1.0/3.0)*deltaTime))
          var v__11025 = particles[p].velocity_new
          v__11025[0] += tmp__11024[0]
          v__11025[1] += tmp__11024[1]
          v__11025[2] += tmp__11024[2]
          particles[p].velocity_new = v__11025
          particles[p].velocity = vv_add_double_3(particles[p].velocity_old, vs_mul_double_3(particles[p].velocity_t, (double(0.5)*deltaTime)))
          particles[p].temperature_new += (((1.0/3.0)*deltaTime)*particles[p].temperature_t)
          particles[p].temperature = (particles[p].temperature_old+((double(0.5)*deltaTime)*particles[p].temperature_t))
        else
          if (Integrator_stage==3) then
            var tmp = vs_mul_double_3(particles[p].position_t, ((1.0/3.0)*deltaTime))
            var v = particles[p].position_new
            v[0] += tmp[0]
            v[1] += tmp[1]
            v[2] += tmp[2]
            particles[p].position_new = v
            particles[p].position = vv_add_double_3(particles[p].position_old, vs_mul_double_3(particles[p].position_t, (1.0*deltaTime)))
            var tmp__11028 = vs_mul_double_3(particles[p].velocity_t, ((1.0/3.0)*deltaTime))
            var v__11029 = particles[p].velocity_new
            v__11029[0] += tmp__11028[0]
            v__11029[1] += tmp__11028[1]
            v__11029[2] += tmp__11028[2]
            particles[p].velocity_new = v__11029
            particles[p].velocity = vv_add_double_3(particles[p].velocity_old, vs_mul_double_3(particles[p].velocity_t, (1.0*deltaTime)))
            particles[p].temperature_new += (((1.0/3.0)*deltaTime)*particles[p].temperature_t)
            particles[p].temperature = (particles[p].temperature_old+((1.0*deltaTime)*particles[p].temperature_t))
          else
            particles[p].position = vv_add_double_3(particles[p].position_new, vs_mul_double_3(particles[p].position_t, ((1.0/6.0)*deltaTime)))
            particles[p].velocity = vv_add_double_3(particles[p].velocity_new, vs_mul_double_3(particles[p].velocity_t, ((1.0/6.0)*deltaTime)))
            particles[p].temperature = (particles[p].temperature_new+(((1.0/6.0)*deltaTime)*particles[p].temperature_t))
          end
        end
      end
    end
  end
end

__demand(__parallel, __cuda)
task Particles_UpdateAuxiliaryStep1(particles : region(ispace(int1d), particles_columns),
                                    BC_xBCLeftParticles : int32, BC_xBCRightParticles : int32,
                                    BC_yBCLeftParticles : int32, BC_yBCRightParticles : int32,
                                    BC_zBCLeftParticles : int32, BC_zBCRightParticles : int32,
                                    Grid_xOrigin : double, Grid_xWidth : double,
                                    Grid_yOrigin : double, Grid_yWidth : double,
                                    Grid_zOrigin : double, Grid_zWidth : double,
                                    Particles_restitutionCoeff : double)
where
  reads(particles.{position, velocity, velocity_t, __valid}),
  reads writes(particles.{position_ghost, velocity_ghost, velocity_t_ghost})
do
  __demand(__openmp)
  for p in particles do
    if particles[p].__valid then
      particles[p].position_ghost[0] = particles[p].position[0]
      particles[p].position_ghost[1] = particles[p].position[1]
      particles[p].position_ghost[2] = particles[p].position[2]
      particles[p].velocity_ghost[0] = particles[p].velocity[0]
      particles[p].velocity_ghost[1] = particles[p].velocity[1]
      particles[p].velocity_ghost[2] = particles[p].velocity[2]
      particles[p].velocity_t_ghost[0] = particles[p].velocity_t[0]
      particles[p].velocity_t_ghost[1] = particles[p].velocity_t[1]
      particles[p].velocity_t_ghost[2] = particles[p].velocity_t[2]
      if (particles[p].position[0]<Grid_xOrigin) then
        if (BC_xBCLeftParticles==0) then
          particles[p].position_ghost[0] = (particles[p].position[0]+Grid_xWidth)
        else
          particles[p].position_ghost[0] = Grid_xOrigin
          var impulse = ((-(1.0+Particles_restitutionCoeff))*particles[p].velocity[0])
          if (impulse<=0.0) then
            particles[p].velocity_ghost[0] += impulse
          end
          var contact_force = (double(-1)*particles[p].velocity_t[0])
          if (contact_force>0.0) then
            particles[p].velocity_t_ghost[0] += contact_force
          end
        end
      end
      if (particles[p].position[0]>(Grid_xOrigin+Grid_xWidth)) then
        if (BC_xBCRightParticles==0) then
          particles[p].position_ghost[0] = (particles[p].position[0]-Grid_xWidth)
        else
          particles[p].position_ghost[0] = (Grid_xOrigin+Grid_xWidth)
          var impulse = ((-(1.0+Particles_restitutionCoeff))*particles[p].velocity[0])
          if (impulse>=0.0) then
            particles[p].velocity_ghost[0] += impulse
          end
          var contact_force = (double(-1)*particles[p].velocity_t[0])
          if (contact_force<0.0) then
            particles[p].velocity_t_ghost[0] += contact_force
          end
        end
      end
      if (particles[p].position[1]<Grid_yOrigin) then
        if (BC_yBCLeftParticles==0) then
          particles[p].position_ghost[1] = (particles[p].position[1]+Grid_yWidth)
        else
          particles[p].position_ghost[1] = Grid_yOrigin
          var impulse = ((-(1.0+Particles_restitutionCoeff))*particles[p].velocity[1])
          if (impulse<=0.0) then
            particles[p].velocity_ghost[1] += impulse
          end
          var contact_force = (double(-1)*particles[p].velocity_t[1])
          if (contact_force>0.0) then
            particles[p].velocity_t_ghost[1] += contact_force
          end
        end
      end
      if (particles[p].position[1]>(Grid_yOrigin+Grid_yWidth)) then
        if (BC_yBCRightParticles==0) then
          particles[p].position_ghost[1] = (particles[p].position[1]-Grid_yWidth)
        else
          particles[p].position_ghost[1] = (Grid_yOrigin+Grid_yWidth)
          var impulse = ((-(1.0+Particles_restitutionCoeff))*particles[p].velocity[1])
          if (impulse>=0.0) then
            particles[p].velocity_ghost[1] += impulse
          end
          var contact_force = (double(-1)*particles[p].velocity_t[1])
          if (contact_force<0.0) then
            particles[p].velocity_t_ghost[1] += contact_force
          end
        end
      end
      if (particles[p].position[2]<Grid_zOrigin) then
        if (BC_zBCLeftParticles==0) then
          particles[p].position_ghost[2] = (particles[p].position[2]+Grid_zWidth)
        else
          particles[p].position_ghost[2] = Grid_zOrigin
          var impulse = ((-(1.0+Particles_restitutionCoeff))*particles[p].velocity[2])
          if (impulse<=0.0) then
            particles[p].velocity_ghost[2] += impulse
          end
          var contact_force = (double(-1)*particles[p].velocity_t[2])
          if (contact_force>0.0) then
            particles[p].velocity_t_ghost[2] += contact_force
          end
        end
      end
      if (particles[p].position[2]>(Grid_zOrigin+Grid_zWidth)) then
        if (BC_zBCRightParticles==0) then
          particles[p].position_ghost[2] = (particles[p].position[2]-Grid_zWidth)
        else
          particles[p].position_ghost[2] = (Grid_zOrigin+Grid_zWidth)
          var impulse = ((-(1.0+Particles_restitutionCoeff))*particles[p].velocity[2])
          if (impulse>=0.0) then
            particles[p].velocity_ghost[2] += impulse
          end
          var contact_force = (double(-1)*particles[p].velocity_t[2])
          if (contact_force<0.0) then
            particles[p].velocity_t_ghost[2] += contact_force
          end
        end
      end
    end
  end
end

__demand(__parallel, __cuda)
task Particles_UpdateAuxiliaryStep2(particles : region(ispace(int1d), particles_columns))
where
  reads(particles.{position_ghost, velocity_ghost, velocity_t_ghost, __valid}),
  reads writes(particles.{position, velocity, velocity_t})
do
  __demand(__openmp)
  for p in particles do
    if particles[p].__valid then
      particles[p].position = particles[p].position_ghost
      particles[p].velocity = particles[p].velocity_ghost
      particles[p].velocity_t = particles[p].velocity_t_ghost
    end
  end
end

task Particles_DeleteEscapingParticles(particles : region(ispace(int1d), particles_columns),
                                       Grid_xRealOrigin : double, Grid_xRealWidth : double,
                                       Grid_yRealOrigin : double, Grid_yRealWidth : double,
                                       Grid_zRealOrigin : double, Grid_zRealWidth : double) : int64
where
  reads(particles.position),
  reads writes(particles.__valid)
do
  var acc = int64(0)
  __demand(__openmp)
  for p in particles do
    if particles[p].__valid then
      var min_x = Grid_xRealOrigin
      var max_x = (Grid_xRealOrigin+Grid_xRealWidth)
      var min_y = Grid_yRealOrigin
      var max_y = (Grid_yRealOrigin+Grid_yRealWidth)
      var min_z = Grid_zRealOrigin
      var max_z = (Grid_zRealOrigin+Grid_zRealWidth)
      var pos = particles[p].position
      if ((((((pos[0]>max_x) or (pos[0]<min_x)) or (pos[1]>max_y)) or (pos[1]<min_y)) or (pos[2]>max_z)) or (pos[2]<min_z)) then
        particles[p].__valid = false
        acc += (-int64(1))
      end
    end
  end
  return acc
end

-------------------------------------------------------------------------------
-- MAIN FUNCTION
-------------------------------------------------------------------------------

__forbid(__optimize)
task work(config : Config)
  var NX = config.Mapping.xTiles
  var NY = config.Mapping.yTiles
  var NZ = config.Mapping.zTiles
  var Grid_xNum = config.Grid.xNum
  var Grid_yNum = config.Grid.yNum
  var Grid_zNum = config.Grid.zNum
  var Grid_xOrigin = config.Grid.origin[0]
  var Grid_yOrigin = config.Grid.origin[1]
  var Grid_zOrigin = config.Grid.origin[2]
  var Grid_xWidth = config.Grid.xWidth
  var Grid_yWidth = config.Grid.yWidth
  var Grid_zWidth = config.Grid.zWidth
  var Grid_xCellWidth = (Grid_xWidth/Grid_xNum)
  var Grid_yCellWidth = (Grid_yWidth/Grid_yNum)
  var Grid_zCellWidth = (Grid_zWidth/Grid_zNum)
  var Grid_cellVolume = ((Grid_xCellWidth*Grid_yCellWidth)*Grid_zCellWidth)
  var Grid_dXYZInverseSquare = (((((1/Grid_xCellWidth)*1)/Grid_xCellWidth)+(((1/Grid_yCellWidth)*1)/Grid_yCellWidth))+(((1/Grid_zCellWidth)*1)/Grid_zCellWidth))
  var BC_xBCPeriodic = (config.BC.xBCLeft == SCHEMA.FlowBC_Periodic)
  var BC_xSign = array(double(0.1), double(0.1), double(0.1))
  var BC_xPosVelocity = array(double(0.1), double(0.1), double(0.1))
  var BC_xNegVelocity = array(double(0.1), double(0.1), double(0.1))
  var BC_xPosTemperature = 0.0
  var BC_xNegTemperature = 0.0
  var BC_yBCPeriodic = (config.BC.yBCLeft == SCHEMA.FlowBC_Periodic)
  var BC_ySign = array(double(0.1), double(0.1), double(0.1))
  var BC_yPosVelocity = array(double(0.1), double(0.1), double(0.1))
  var BC_yNegVelocity = array(double(0.1), double(0.1), double(0.1))
  var BC_yPosTemperature = 0.0
  var BC_yNegTemperature = 0.0
  var BC_zBCPeriodic = (config.BC.zBCLeft == SCHEMA.FlowBC_Periodic)
  var BC_zSign = array(double(0.1), double(0.1), double(0.1))
  var BC_zPosVelocity = array(double(0.1), double(0.1), double(0.1))
  var BC_zNegVelocity = array(double(0.1), double(0.1), double(0.1))
  var BC_zPosTemperature = 0.0
  var BC_zNegTemperature = 0.0
  var BC_xBCLeftParticles = int32(-1)
  var BC_xBCRightParticles = int32(-1)
  var BC_yBCLeftParticles = int32(-1)
  var BC_yBCRightParticles = int32(-1)
  var BC_zBCLeftParticles = int32(-1)
  var BC_zBCRightParticles = int32(-1)
  var Grid_xBnum = 1
  if BC_xBCPeriodic then Grid_xBnum = 0 end
  var Grid_yBnum = 1
  if BC_yBCPeriodic then Grid_yBnum = 0 end
  var Grid_zBnum = 1
  if BC_zBCPeriodic then Grid_zBnum = 0 end
  var Grid_xRealOrigin = (Grid_xOrigin-(Grid_xCellWidth*Grid_xBnum))
  var Grid_yRealOrigin = (Grid_yOrigin-(Grid_yCellWidth*Grid_yBnum))
  var Grid_zRealOrigin = (Grid_zOrigin-(Grid_zCellWidth*Grid_zBnum))
  var Grid_xRealWidth = (Grid_xWidth+(2*(Grid_xCellWidth*Grid_xBnum)))
  var Grid_yRealWidth = (Grid_yWidth+(2*(Grid_yCellWidth*Grid_yBnum)))
  var Grid_zRealWidth = (Grid_zWidth+(2*(Grid_zCellWidth*Grid_zBnum)))
  var Integrator_simTime = 0.0
  var Integrator_time_old = 0.0
  var Integrator_timeStep = 0
  var Integrator_deltaTime = double(0.0001)
  var Integrator_stage = 0
  var Integrator_maxConvectiveSpectralRadius = 0.0
  var Integrator_maxViscousSpectralRadius = 0.0
  var Integrator_maxHeatConductionSpectralRadius = 0.0
  var Flow_gasConstant = config.Flow.gasConstant
  var Flow_gamma = config.Flow.gamma
  var Flow_prandtl = config.Flow.prandtl
  var Flow_viscosityModel = config.Flow.viscosityModel
  var Flow_constantVisc = config.Flow.constantVisc
  var Flow_powerlawViscRef = config.Flow.powerlawViscRef
  var Flow_powerlawTempRef = config.Flow.powerlawTempRef
  var Flow_sutherlandViscRef = config.Flow.sutherlandViscRef
  var Flow_sutherlandTempRef = config.Flow.sutherlandTempRef
  var Flow_sutherlandSRef = config.Flow.sutherlandSRef
  var Flow_initParams = config.Flow.initParams
  var Flow_bodyForce = config.Flow.bodyForce
  var Flow_averagePD = 0.0
  var Flow_averageDissipation = 0.0
  var Flow_averageFe = 0.0
  var Flow_averageK = 0.0
  var Particles_maxNum = config.Particles.maxNum
  var Particles_restitutionCoeff = config.Particles.restitutionCoeff
  var Particles_convectiveCoeff = config.Particles.convectiveCoeff
  var Particles_heatCapacity = config.Particles.heatCapacity
  var Particles_density = config.Particles.density
  var Particles_bodyForce = config.Particles.bodyForce
  var Particles_maxSkew = config.Particles.maxSkew
  var Particles_maxXferNum = config.Particles.maxXferNum
  var Particles_number = int64(0)
  var Radiation_qa = config.Radiation.qa
  var Radiation_qs = config.Radiation.qs
  var Radiation_xNum = config.Radiation.xNum
  var Radiation_yNum = config.Radiation.yNum
  var Radiation_zNum = config.Radiation.zNum
  var Radiation_xBnum = 0
  var Radiation_yBnum = 0
  var Radiation_zBnum = 0
  var Radiation_xPeriodic = false
  var Radiation_yPeriodic = false
  var Radiation_zPeriodic = false
  var Radiation_xCellWidth = (Grid_xWidth/Radiation_xNum)
  var Radiation_yCellWidth = (Grid_yWidth/Radiation_yNum)
  var Radiation_zCellWidth = (Grid_zWidth/Radiation_zNum)
  var Radiation_cellVolume = ((Radiation_xCellWidth*Radiation_yCellWidth)*Radiation_zCellWidth)
  var is = ispace(int3d, int3d({x = (Grid_xNum+(2*Grid_xBnum)), y = (Grid_yNum+(2*Grid_yBnum)), z = (Grid_zNum+(2*Grid_zBnum))}))
  var Fluid = region(is, Fluid_columns)
  var Fluid_copy = region(is, Fluid_columns)
  var is__11726 = ispace(int1d, int1d((ceil(((Particles_maxNum/((NX*NY)*NZ))*Particles_maxSkew))*((NX*NY)*NZ))))
  var particles = region(is__11726, particles_columns)
  var particles_copy = region(is__11726, particles_columns)
  var is__11729 = ispace(int3d, int3d({x = (Radiation_xNum+(2*Radiation_xBnum)), y = (Radiation_yNum+(2*Radiation_yBnum)), z = (Radiation_zNum+(2*Radiation_zBnum))}))
  var Radiation = region(is__11729, Radiation_columns)
  var Radiation_copy = region(is__11729, Radiation_columns)
  var primColors = ispace(int3d, int3d({NX, NY, NZ}))
  regentlib.assert(((Grid_xNum%NX)==0), "Uneven partitioning of fluid grid on x")
  regentlib.assert(((Grid_yNum%NY)==0), "Uneven partitioning of fluid grid on y")
  regentlib.assert(((Grid_zNum%NZ)==0), "Uneven partitioning of fluid grid on z")
  var coloring = regentlib.c.legion_domain_point_coloring_create()
  for c in primColors do
    var rect = rect3d({lo = int3d({x = (Grid_xBnum+((Grid_xNum/NX)*c.x)), y = (Grid_yBnum+((Grid_yNum/NY)*c.y)), z = (Grid_zBnum+((Grid_zNum/NZ)*c.z))}), hi = int3d({x = ((Grid_xBnum+((Grid_xNum/NX)*(c.x+1)))-1), y = ((Grid_yBnum+((Grid_yNum/NY)*(c.y+1)))-1), z = ((Grid_zBnum+((Grid_zNum/NZ)*(c.z+1)))-1)})})
    if (c.x==0) then
      rect.lo.x -= Grid_xBnum
    end
    if (c.x==(NX-1)) then
      rect.hi.x += Grid_xBnum
    end
    if (c.y==0) then
      rect.lo.y -= Grid_yBnum
    end
    if (c.y==(NY-1)) then
      rect.hi.y += Grid_yBnum
    end
    if (c.z==0) then
      rect.lo.z -= Grid_zBnum
    end
    if (c.z==(NZ-1)) then
      rect.hi.z += Grid_zBnum
    end
    regentlib.c.legion_domain_point_coloring_color_domain(coloring, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect))
  end
  var Fluid_primPart = partition(disjoint, Fluid, coloring, primColors)
  var Fluid_copy_primPart = partition(disjoint, Fluid_copy, coloring, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(coloring)
  regentlib.assert(((Particles_maxNum%((NX*NY)*NZ))==0), "Uneven partitioning of particles")
  var coloring__11738 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var rBase : int64
        for rStart in particles do
          rBase = int64((rStart+(((((z*NX)*NY)+(y*NX))+x)*ceil(((Particles_maxNum/((NX*NY)*NZ))*Particles_maxSkew)))))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(coloring__11738, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({rBase, ((rBase+ceil(((Particles_maxNum/((NX*NY)*NZ))*Particles_maxSkew)))-1)})))
      end
    end
  end
  var particles_primPart = partition(disjoint, particles, coloring__11738, primColors)
  var particles_copy_primPart = partition(disjoint, particles_copy, coloring__11738, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(coloring__11738)
  var particles_queue_0 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_0 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_0 = partition(disjoint, particles_queue_0, srcColoring, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring)
  var dstColoring = regentlib.c.legion_domain_point_coloring_create()
  var colorOff = int3d({0, 0, 1})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_0[(((c-colorOff)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_0 = partition(aliased, particles_queue_0, dstColoring, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring)
  var particles_queue_1 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__11761 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_1 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__11761, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_1 = partition(disjoint, particles_queue_1, srcColoring__11761, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__11761)
  var dstColoring__11768 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__11769 = int3d({0, 0, -1})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_1[(((c-colorOff__11769)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__11768, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_1 = partition(aliased, particles_queue_1, dstColoring__11768, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__11768)
  var particles_queue_2 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__11775 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_2 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__11775, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_2 = partition(disjoint, particles_queue_2, srcColoring__11775, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__11775)
  var dstColoring__11782 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__11783 = int3d({0, 1, 0})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_2[(((c-colorOff__11783)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__11782, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_2 = partition(aliased, particles_queue_2, dstColoring__11782, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__11782)
  var particles_queue_3 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__11789 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_3 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__11789, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_3 = partition(disjoint, particles_queue_3, srcColoring__11789, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__11789)
  var dstColoring__11796 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__11797 = int3d({0, 1, 1})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_3[(((c-colorOff__11797)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__11796, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_3 = partition(aliased, particles_queue_3, dstColoring__11796, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__11796)
  var particles_queue_4 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__11803 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_4 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__11803, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_4 = partition(disjoint, particles_queue_4, srcColoring__11803, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__11803)
  var dstColoring__11810 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__11811 = int3d({0, 1, -1})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_4[(((c-colorOff__11811)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__11810, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_4 = partition(aliased, particles_queue_4, dstColoring__11810, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__11810)
  var particles_queue_5 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__11817 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_5 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__11817, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_5 = partition(disjoint, particles_queue_5, srcColoring__11817, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__11817)
  var dstColoring__11824 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__11825 = int3d({0, -1, 0})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_5[(((c-colorOff__11825)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__11824, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_5 = partition(aliased, particles_queue_5, dstColoring__11824, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__11824)
  var particles_queue_6 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__11831 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_6 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__11831, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_6 = partition(disjoint, particles_queue_6, srcColoring__11831, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__11831)
  var dstColoring__11838 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__11839 = int3d({0, -1, 1})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_6[(((c-colorOff__11839)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__11838, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_6 = partition(aliased, particles_queue_6, dstColoring__11838, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__11838)
  var particles_queue_7 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__11845 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_7 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__11845, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_7 = partition(disjoint, particles_queue_7, srcColoring__11845, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__11845)
  var dstColoring__11852 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__11853 = int3d({0, -1, -1})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_7[(((c-colorOff__11853)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__11852, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_7 = partition(aliased, particles_queue_7, dstColoring__11852, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__11852)
  var particles_queue_8 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__11859 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_8 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__11859, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_8 = partition(disjoint, particles_queue_8, srcColoring__11859, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__11859)
  var dstColoring__11866 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__11867 = int3d({1, 0, 0})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_8[(((c-colorOff__11867)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__11866, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_8 = partition(aliased, particles_queue_8, dstColoring__11866, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__11866)
  var particles_queue_9 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__11873 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_9 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__11873, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_9 = partition(disjoint, particles_queue_9, srcColoring__11873, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__11873)
  var dstColoring__11880 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__11881 = int3d({1, 0, 1})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_9[(((c-colorOff__11881)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__11880, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_9 = partition(aliased, particles_queue_9, dstColoring__11880, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__11880)
  var particles_queue_10 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__11887 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_10 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__11887, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_10 = partition(disjoint, particles_queue_10, srcColoring__11887, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__11887)
  var dstColoring__11894 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__11895 = int3d({1, 0, -1})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_10[(((c-colorOff__11895)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__11894, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_10 = partition(aliased, particles_queue_10, dstColoring__11894, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__11894)
  var particles_queue_11 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__11901 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_11 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__11901, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_11 = partition(disjoint, particles_queue_11, srcColoring__11901, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__11901)
  var dstColoring__11908 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__11909 = int3d({1, 1, 0})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_11[(((c-colorOff__11909)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__11908, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_11 = partition(aliased, particles_queue_11, dstColoring__11908, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__11908)
  var particles_queue_12 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__11915 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_12 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__11915, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_12 = partition(disjoint, particles_queue_12, srcColoring__11915, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__11915)
  var dstColoring__11922 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__11923 = int3d({1, 1, 1})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_12[(((c-colorOff__11923)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__11922, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_12 = partition(aliased, particles_queue_12, dstColoring__11922, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__11922)
  var particles_queue_13 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__11929 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_13 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__11929, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_13 = partition(disjoint, particles_queue_13, srcColoring__11929, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__11929)
  var dstColoring__11936 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__11937 = int3d({1, 1, -1})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_13[(((c-colorOff__11937)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__11936, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_13 = partition(aliased, particles_queue_13, dstColoring__11936, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__11936)
  var particles_queue_14 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__11943 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_14 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__11943, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_14 = partition(disjoint, particles_queue_14, srcColoring__11943, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__11943)
  var dstColoring__11950 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__11951 = int3d({1, -1, 0})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_14[(((c-colorOff__11951)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__11950, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_14 = partition(aliased, particles_queue_14, dstColoring__11950, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__11950)
  var particles_queue_15 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__11957 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_15 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__11957, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_15 = partition(disjoint, particles_queue_15, srcColoring__11957, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__11957)
  var dstColoring__11964 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__11965 = int3d({1, -1, 1})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_15[(((c-colorOff__11965)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__11964, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_15 = partition(aliased, particles_queue_15, dstColoring__11964, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__11964)
  var particles_queue_16 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__11971 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_16 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__11971, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_16 = partition(disjoint, particles_queue_16, srcColoring__11971, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__11971)
  var dstColoring__11978 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__11979 = int3d({1, -1, -1})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_16[(((c-colorOff__11979)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__11978, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_16 = partition(aliased, particles_queue_16, dstColoring__11978, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__11978)
  var particles_queue_17 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__11985 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_17 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__11985, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_17 = partition(disjoint, particles_queue_17, srcColoring__11985, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__11985)
  var dstColoring__11992 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__11993 = int3d({-1, 0, 0})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_17[(((c-colorOff__11993)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__11992, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_17 = partition(aliased, particles_queue_17, dstColoring__11992, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__11992)
  var particles_queue_18 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__11999 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_18 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__11999, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_18 = partition(disjoint, particles_queue_18, srcColoring__11999, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__11999)
  var dstColoring__12006 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__12007 = int3d({-1, 0, 1})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_18[(((c-colorOff__12007)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__12006, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_18 = partition(aliased, particles_queue_18, dstColoring__12006, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__12006)
  var particles_queue_19 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__12013 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_19 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__12013, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_19 = partition(disjoint, particles_queue_19, srcColoring__12013, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__12013)
  var dstColoring__12020 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__12021 = int3d({-1, 0, -1})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_19[(((c-colorOff__12021)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__12020, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_19 = partition(aliased, particles_queue_19, dstColoring__12020, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__12020)
  var particles_queue_20 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__12027 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_20 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__12027, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_20 = partition(disjoint, particles_queue_20, srcColoring__12027, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__12027)
  var dstColoring__12034 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__12035 = int3d({-1, 1, 0})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_20[(((c-colorOff__12035)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__12034, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_20 = partition(aliased, particles_queue_20, dstColoring__12034, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__12034)
  var particles_queue_21 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__12041 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_21 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__12041, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_21 = partition(disjoint, particles_queue_21, srcColoring__12041, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__12041)
  var dstColoring__12048 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__12049 = int3d({-1, 1, 1})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_21[(((c-colorOff__12049)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__12048, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_21 = partition(aliased, particles_queue_21, dstColoring__12048, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__12048)
  var particles_queue_22 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__12055 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_22 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__12055, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_22 = partition(disjoint, particles_queue_22, srcColoring__12055, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__12055)
  var dstColoring__12062 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__12063 = int3d({-1, 1, -1})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_22[(((c-colorOff__12063)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__12062, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_22 = partition(aliased, particles_queue_22, dstColoring__12062, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__12062)
  var particles_queue_23 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__12069 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_23 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__12069, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_23 = partition(disjoint, particles_queue_23, srcColoring__12069, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__12069)
  var dstColoring__12076 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__12077 = int3d({-1, -1, 0})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_23[(((c-colorOff__12077)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__12076, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_23 = partition(aliased, particles_queue_23, dstColoring__12076, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__12076)
  var particles_queue_24 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__12083 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_24 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__12083, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_24 = partition(disjoint, particles_queue_24, srcColoring__12083, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__12083)
  var dstColoring__12090 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__12091 = int3d({-1, -1, 1})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_24[(((c-colorOff__12091)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__12090, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_24 = partition(aliased, particles_queue_24, dstColoring__12090, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__12090)
  var particles_queue_25 = region(ispace(int1d, int1d((Particles_maxXferNum*((NX*NY)*NZ)))), int8[376])
  var srcColoring__12097 = regentlib.c.legion_domain_point_coloring_create()
  for z : int32 = 0, NZ do
    for y : int32 = 0, NY do
      for x : int32 = 0, NX do
        var qBase : int64
        for qStart in particles_queue_25 do
          qBase = int64((qStart+(((((z*NX)*NY)+(y*NX))+x)*Particles_maxXferNum)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(srcColoring__12097, regentlib.c.legion_domain_point_t(int3d({x, y, z})), regentlib.c.legion_domain_t(rect1d({qBase, ((qBase+Particles_maxXferNum)-1)})))
      end
    end
  end
  var particles_qSrcPart_25 = partition(disjoint, particles_queue_25, srcColoring__12097, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(srcColoring__12097)
  var dstColoring__12104 = regentlib.c.legion_domain_point_coloring_create()
  var colorOff__12105 = int3d({-1, -1, -1})
  for c in primColors do
    var srcBase : int64
    for qptr in particles_qSrcPart_25[(((c-colorOff__12105)+{NX, NY, NZ})%{NX, NY, NZ})] do
      srcBase = int64(int1d(qptr))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(dstColoring__12104, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect1d({srcBase, ((srcBase+Particles_maxXferNum)-1)})))
  end
  var particles_qDstPart_25 = partition(aliased, particles_queue_25, dstColoring__12104, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(dstColoring__12104)
  regentlib.assert(((Radiation_xNum%NX)==0), "Uneven partitioning of radiation grid on x")
  regentlib.assert(((Radiation_yNum%NY)==0), "Uneven partitioning of radiation grid on y")
  regentlib.assert(((Radiation_zNum%NZ)==0), "Uneven partitioning of radiation grid on z")
  var coloring__12110 = regentlib.c.legion_domain_point_coloring_create()
  for c in primColors do
    var rect = rect3d({lo = int3d({x = (Radiation_xBnum+((Radiation_xNum/NX)*c.x)), y = (Radiation_yBnum+((Radiation_yNum/NY)*c.y)), z = (Radiation_zBnum+((Radiation_zNum/NZ)*c.z))}), hi = int3d({x = ((Radiation_xBnum+((Radiation_xNum/NX)*(c.x+1)))-1), y = ((Radiation_yBnum+((Radiation_yNum/NY)*(c.y+1)))-1), z = ((Radiation_zBnum+((Radiation_zNum/NZ)*(c.z+1)))-1)})})
    if (c.x==0) then
      rect.lo.x -= Radiation_xBnum
    end
    if (c.x==(NX-1)) then
      rect.hi.x += Radiation_xBnum
    end
    if (c.y==0) then
      rect.lo.y -= Radiation_yBnum
    end
    if (c.y==(NY-1)) then
      rect.hi.y += Radiation_yBnum
    end
    if (c.z==0) then
      rect.lo.z -= Radiation_zBnum
    end
    if (c.z==(NZ-1)) then
      rect.hi.z += Radiation_zBnum
    end
    regentlib.c.legion_domain_point_coloring_color_domain(coloring__12110, regentlib.c.legion_domain_point_t(c), regentlib.c.legion_domain_t(rect))
  end
  var Radiation_primPart = partition(disjoint, Radiation, coloring__12110, primColors)
  var Radiation_copy_primPart = partition(disjoint, Radiation_copy, coloring__12110, primColors)
  regentlib.c.legion_domain_point_coloring_destroy(coloring__12110)
  __parallelize_with Fluid_primPart, particles_primPart, Radiation_primPart, primColors, (image(Fluid, particles_primPart, particles.cell)<=Fluid_primPart) do
    particles_initValidField(particles)
    if ((not ((Grid_xNum%Radiation_xNum)==0)) or ((not ((Grid_yNum%Radiation_yNum)==0)) or (not ((Grid_zNum%Radiation_zNum)==0)))) then
      regentlib.assert(false, "Inexact coarsening factor")
    end
    SetCoarseningField(Fluid, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum, Radiation_xBnum, Radiation_xNum, Radiation_yBnum, Radiation_yNum, Radiation_zBnum, Radiation_zNum)
    if ((config.BC.xBCLeft == SCHEMA.FlowBC_Periodic) and (config.BC.xBCRight == SCHEMA.FlowBC_Periodic)) then
      BC_xSign = array(1.0, 1.0, 1.0)
      BC_xPosVelocity = array(0.0, 0.0, 0.0)
      BC_xNegVelocity = array(0.0, 0.0, 0.0)
      BC_xPosTemperature = -1.0
      BC_xNegTemperature = -1.0
      BC_xBCLeftParticles = 0
      BC_xBCRightParticles = 0
    elseif ((config.BC.xBCLeft == SCHEMA.FlowBC_Symmetry) and (config.BC.xBCRight == SCHEMA.FlowBC_Symmetry)) then
      BC_xSign = array(-1.0, 1.0, 1.0)
      BC_xPosVelocity = array(0.0, 0.0, 0.0)
      BC_xNegVelocity = array(0.0, 0.0, 0.0)
      BC_xPosTemperature = -1.0
      BC_xNegTemperature = -1.0
      BC_xBCLeftParticles = 1
      BC_xBCRightParticles = 1
    elseif ((config.BC.xBCLeft == SCHEMA.FlowBC_AdiabaticWall) and (config.BC.xBCRight == SCHEMA.FlowBC_AdiabaticWall)) then
      BC_xSign = array(-1.0, -1.0, -1.0)
      BC_xPosVelocity = array((2*config.BC.xBCRightVel[0]), (2*config.BC.xBCRightVel[1]), (2*config.BC.xBCRightVel[2]))
      BC_xNegVelocity = array((2*config.BC.xBCLeftVel[0]), (2*config.BC.xBCLeftVel[1]), (2*config.BC.xBCLeftVel[2]))
      BC_xPosTemperature = -1.0
      BC_xNegTemperature = -1.0
      BC_xBCLeftParticles = 1
      BC_xBCRightParticles = 1
    elseif ((config.BC.xBCLeft == SCHEMA.FlowBC_IsothermalWall) and (config.BC.xBCRight == SCHEMA.FlowBC_IsothermalWall)) then
      BC_xSign = array(-1.0, -1.0, -1.0)
      BC_xPosVelocity = array((2*config.BC.xBCRightVel[0]), (2*config.BC.xBCRightVel[1]), (2*config.BC.xBCRightVel[2]))
      BC_xNegVelocity = array((2*config.BC.xBCLeftVel[0]), (2*config.BC.xBCLeftVel[1]), (2*config.BC.xBCLeftVel[2]))
      BC_xPosTemperature = config.BC.xBCRightTemp
      BC_xNegTemperature = config.BC.xBCLeftTemp
      BC_xBCLeftParticles = 1
      BC_xBCRightParticles = 1
    else
      regentlib.assert(false, "Boundary conditions in x not implemented")
    end
    if ((config.BC.yBCLeft == SCHEMA.FlowBC_Periodic) and (config.BC.yBCRight == SCHEMA.FlowBC_Periodic)) then
      BC_ySign = array(1.0, 1.0, 1.0)
      BC_yPosVelocity = array(0.0, 0.0, 0.0)
      BC_yNegVelocity = array(0.0, 0.0, 0.0)
      BC_yPosTemperature = -1.0
      BC_yNegTemperature = -1.0
      BC_yBCLeftParticles = 0
      BC_yBCRightParticles = 0
    elseif ((config.BC.yBCLeft == SCHEMA.FlowBC_Symmetry) and (config.BC.yBCRight == SCHEMA.FlowBC_Symmetry)) then
      BC_ySign = array(1.0, -1.0, 1.0)
      BC_yPosVelocity = array(0.0, 0.0, 0.0)
      BC_yNegVelocity = array(0.0, 0.0, 0.0)
      BC_yPosTemperature = -1.0
      BC_yNegTemperature = -1.0
      BC_yBCLeftParticles = 1
      BC_yBCRightParticles = 1
    elseif ((config.BC.yBCLeft == SCHEMA.FlowBC_AdiabaticWall) and (config.BC.yBCRight == SCHEMA.FlowBC_AdiabaticWall)) then
      BC_ySign = array(-1.0, -1.0, -1.0)
      BC_yPosVelocity = array((2*config.BC.yBCRightVel[0]), (2*config.BC.yBCRightVel[1]), (2*config.BC.yBCRightVel[2]))
      BC_yNegVelocity = array((2*config.BC.yBCLeftVel[0]), (2*config.BC.yBCLeftVel[1]), (2*config.BC.yBCLeftVel[2]))
      BC_yPosTemperature = -1.0
      BC_yNegTemperature = -1.0
      BC_yBCLeftParticles = 1
      BC_yBCRightParticles = 1
    elseif ((config.BC.yBCLeft == SCHEMA.FlowBC_IsothermalWall) and (config.BC.yBCRight == SCHEMA.FlowBC_IsothermalWall)) then
      BC_ySign = array(-1.0, -1.0, -1.0)
      BC_yPosVelocity = array((2*config.BC.yBCRightVel[0]), (2*config.BC.yBCRightVel[1]), (2*config.BC.yBCRightVel[2]))
      BC_yNegVelocity = array((2*config.BC.yBCLeftVel[0]), (2*config.BC.yBCLeftVel[1]), (2*config.BC.yBCLeftVel[2]))
      BC_yPosTemperature = config.BC.yBCRightTemp
      BC_yNegTemperature = config.BC.yBCLeftTemp
      BC_yBCLeftParticles = 1
      BC_yBCRightParticles = 1
    else
      regentlib.assert(false, "Boundary conditions in y not implemented")
    end
    if ((config.BC.zBCLeft == SCHEMA.FlowBC_Periodic) and (config.BC.zBCRight == SCHEMA.FlowBC_Periodic)) then
      BC_zSign = array(1.0, 1.0, 1.0)
      BC_zPosVelocity = array(0.0, 0.0, 0.0)
      BC_zNegVelocity = array(0.0, 0.0, 0.0)
      BC_zPosTemperature = -1.0
      BC_zNegTemperature = -1.0
      BC_zBCLeftParticles = 0
      BC_zBCRightParticles = 0
    elseif ((config.BC.zBCLeft == SCHEMA.FlowBC_Symmetry) and (config.BC.zBCRight == SCHEMA.FlowBC_Symmetry)) then
      BC_zSign = array(1.0, 1.0, -1.0)
      BC_zPosVelocity = array(0.0, 0.0, 0.0)
      BC_zNegVelocity = array(0.0, 0.0, 0.0)
      BC_zPosTemperature = -1.0
      BC_zNegTemperature = -1.0
      BC_zBCLeftParticles = 1
      BC_zBCRightParticles = 1
    elseif ((config.BC.zBCLeft == SCHEMA.FlowBC_AdiabaticWall) and (config.BC.zBCRight == SCHEMA.FlowBC_AdiabaticWall)) then
      BC_zSign = array(-1.0, -1.0, -1.0)
      BC_zPosVelocity = array((2*config.BC.zBCRightVel[0]), (2*config.BC.zBCRightVel[1]), (2*config.BC.zBCRightVel[2]))
      BC_zNegVelocity = array((2*config.BC.zBCLeftVel[0]), (2*config.BC.zBCLeftVel[1]), (2*config.BC.zBCLeftVel[2]))
      BC_zPosTemperature = -1.0
      BC_zNegTemperature = -1.0
      BC_zBCLeftParticles = 1
      BC_zBCRightParticles = 1
    elseif ((config.BC.zBCLeft == SCHEMA.FlowBC_IsothermalWall) and (config.BC.zBCRight == SCHEMA.FlowBC_IsothermalWall)) then
      BC_zSign = array(-1.0, -1.0, -1.0)
      BC_zPosVelocity = array((2*config.BC.zBCRightVel[0]), (2*config.BC.zBCRightVel[1]), (2*config.BC.zBCRightVel[2]))
      BC_zNegVelocity = array((2*config.BC.zBCLeftVel[0]), (2*config.BC.zBCLeftVel[1]), (2*config.BC.zBCLeftVel[2]))
      BC_zPosTemperature = config.BC.zBCRightTemp
      BC_zNegTemperature = config.BC.zBCLeftTemp
      BC_zBCLeftParticles = 1
      BC_zBCRightParticles = 1
    else
      regentlib.assert(false, "Boundary conditions in z not implemented")
    end
    if (not (((config.BC.xBCLeft == SCHEMA.FlowBC_Periodic) and (config.BC.xBCRight == SCHEMA.FlowBC_Periodic)) or ((not (config.BC.xBCLeft == SCHEMA.FlowBC_Periodic)) and (not (config.BC.xBCRight == SCHEMA.FlowBC_Periodic))))) then
      regentlib.assert(false, "Boundary conditions in x should match for periodicity")
    end
    if (not (((config.BC.yBCLeft == SCHEMA.FlowBC_Periodic) and (config.BC.yBCRight == SCHEMA.FlowBC_Periodic)) or ((not (config.BC.yBCLeft == SCHEMA.FlowBC_Periodic)) and (not (config.BC.yBCRight == SCHEMA.FlowBC_Periodic))))) then
      regentlib.assert(false, "Boundary conditions in y should match for periodicity")
    end
    if (not (((config.BC.zBCLeft == SCHEMA.FlowBC_Periodic) and (config.BC.zBCRight == SCHEMA.FlowBC_Periodic)) or ((not (config.BC.zBCLeft == SCHEMA.FlowBC_Periodic)) and (not (config.BC.zBCRight == SCHEMA.FlowBC_Periodic))))) then
      regentlib.assert(false, "Boundary conditions in z should match for periodicity")
    end
    if (config.Flow.initCase == SCHEMA.FlowInitCase_Restart) then
      Integrator_timeStep = config.Integrator.restartIter
    end
    Flow_InitializeCell(Fluid)
    Flow_InitializeCenterCoordinates(Fluid, Grid_xBnum, Grid_xNum, Grid_xOrigin, Grid_xWidth, Grid_yBnum, Grid_yNum, Grid_yOrigin, Grid_yWidth, Grid_zBnum, Grid_zNum, Grid_zOrigin, Grid_zWidth)
    if (config.Flow.initCase == SCHEMA.FlowInitCase_Uniform) then
      Flow_InitializeUniform(Fluid, Flow_initParams)
    end
    if (config.Flow.initCase == SCHEMA.FlowInitCase_TaylorGreen2DVortex) then
      Flow_InitializeTaylorGreen2D(Fluid, Flow_initParams, Grid_xBnum, Grid_xNum, Grid_xOrigin, Grid_xWidth, Grid_yBnum, Grid_yNum, Grid_yOrigin, Grid_yWidth, Grid_zBnum, Grid_zNum, Grid_zOrigin, Grid_zWidth)
    end
    if (config.Flow.initCase == SCHEMA.FlowInitCase_TaylorGreen3DVortex) then
      Flow_InitializeTaylorGreen3D(Fluid, Flow_initParams, Grid_xBnum, Grid_xNum, Grid_xOrigin, Grid_xWidth, Grid_yBnum, Grid_yNum, Grid_yOrigin, Grid_yWidth, Grid_zBnum, Grid_zNum, Grid_zOrigin, Grid_zWidth)
    end
    if (config.Flow.initCase == SCHEMA.FlowInitCase_Perturbed) then
      Flow_InitializePerturbed(Fluid, Flow_initParams)
    end
    if (config.Flow.initCase == SCHEMA.FlowInitCase_Restart) then
      var dirname = [&int8](C.malloc(256))
      C.snprintf(dirname, 256, "fluid_sample%d_iter%d", config.Mapping.sampleId, config.Integrator.restartIter)
      Fluid_load(primColors, dirname, Fluid, Fluid_copy, Fluid_primPart, Fluid_copy_primPart)
      C.free([&opaque](dirname))
    end
    Flow_UpdateConservedFromPrimitive(Fluid, Flow_gamma, Flow_gasConstant, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
    Flow_UpdateAuxiliaryVelocity(Fluid, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
    Flow_UpdateGhostConservedStep1(Fluid, BC_xNegTemperature, BC_xNegVelocity, BC_xPosTemperature, BC_xPosVelocity, BC_xSign, BC_yNegTemperature, BC_yNegVelocity, BC_yPosTemperature, BC_yPosVelocity, BC_ySign, BC_zNegTemperature, BC_zNegVelocity, BC_zPosTemperature, BC_zPosVelocity, BC_zSign, Flow_gamma, Flow_gasConstant, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
    Flow_UpdateGhostConservedStep2(Fluid, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
    Flow_UpdateGhostVelocityStep1(Fluid, BC_xNegVelocity, BC_xPosVelocity, BC_xSign, BC_yNegVelocity, BC_yPosVelocity, BC_ySign, BC_zNegVelocity, BC_zPosVelocity, BC_zSign, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
    Flow_UpdateGhostVelocityStep2(Fluid, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
    Flow_ComputeVelocityGradientAll(Fluid, Grid_xBnum, Grid_xCellWidth, Grid_xNum, Grid_yBnum, Grid_yCellWidth, Grid_yNum, Grid_zBnum, Grid_zCellWidth, Grid_zNum)
    Flow_UpdateAuxiliaryThermodynamics(Fluid, Flow_gamma, Flow_gasConstant, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
    Flow_UpdateGhostThermodynamicsStep1(Fluid, BC_xNegTemperature, BC_xPosTemperature, BC_yNegTemperature, BC_yPosTemperature, BC_zNegTemperature, BC_zPosTemperature, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
    Flow_UpdateGhostThermodynamicsStep2(Fluid, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
    Flow_UpdateGhostFieldsStep1(Fluid, BC_xNegTemperature, BC_xNegVelocity, BC_xPosTemperature, BC_xPosVelocity, BC_xSign, BC_yNegTemperature, BC_yNegVelocity, BC_yPosTemperature, BC_yPosVelocity, BC_ySign, BC_zNegTemperature, BC_zNegVelocity, BC_zPosTemperature, BC_zPosVelocity, BC_zSign, Flow_gamma, Flow_gasConstant, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
    Flow_UpdateGhostFieldsStep2(Fluid, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
    if (config.Particles.initCase == SCHEMA.ParticlesInitCase_Random) then
      regentlib.assert(false, "Random particle initialization is disabled")
    end
    if (config.Particles.initCase == SCHEMA.ParticlesInitCase_Restart) then
      var dirname = [&int8](C.malloc(256))
      C.snprintf(dirname, 256, "particles_sample%d_iter%d", config.Mapping.sampleId, config.Integrator.restartIter)
      particles_load(primColors, dirname, particles, particles_copy, particles_primPart, particles_copy_primPart)
      C.free([&opaque](dirname))
      Particles_InitializeDensity(particles, Particles_density)
      Particles_number += Particles_CalculateNumber(particles)
    end
    if (config.Particles.initCase == SCHEMA.ParticlesInitCase_Uniform) then
      InitParticlesUniform(particles, Fluid, config, Grid_xBnum, Grid_yBnum, Grid_zBnum)
      Particles_number = int64(((config.Particles.initNum/((config.Mapping.xTiles*config.Mapping.yTiles)*config.Mapping.zTiles))*((config.Mapping.xTiles*config.Mapping.yTiles)*config.Mapping.zTiles)))
    end
    [DOM.DeclSymbols(config)];
    if (config.Radiation.type == SCHEMA.RadiationType_DOM) then
      Radiation_InitializeCell(Radiation);
      [DOM.InitRegions()];
    end
    -- Open long-running files
    var probeFiles = [&&C.FILE](C.malloc(config.IO.probes.length * FILE_PTR_SIZE))
    for i = 0,config.IO.probes.length do
      var filename = [&int8](C.malloc(32))
      C.snprintf(filename, 32, 'probe%d.csv', i)
      if C.fopen(filename, 'r') ~= [&C.FILE](0) then
        var stderr = C.fdopen(2, 'w')
        C.fprintf(stderr, 'Probe file %s already exists.\n', filename)
        C.fflush(stderr)
        C.exit(1)
      end
      probeFiles[i] = C.fopen(filename, 'w')
      if probeFiles[i] == [&C.FILE](0) then
        var stderr = C.fdopen(2, 'w')
        C.fprintf(stderr, 'Cannot open file %s for writing.\n', filename)
        C.fflush(stderr)
        C.exit(1)
      end
      C.fprintf(probeFiles[i], 'timeStep\ttemperature\n')
      C.fflush(probeFiles[i])
      C.free(filename)
    end

    -- Main time-step loop
    while true do
      -- Calculate exit condition, but don't exit yet
      var exitCond =
        Integrator_simTime >= config.Integrator.finalTime or
        Integrator_timeStep >= config.Integrator.maxIter
      -- Perform IO (on user-requested timesteps, and always on the last one)
      if exitCond or Integrator_timeStep % config.IO.headerFrequency == 0 then
        var Flow_minTemperature = math.huge
        var Flow_maxTemperature = -math.huge
        Flow_minTemperature min= CalculateMinTemperature(Fluid, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
        Flow_maxTemperature max= CalculateMaxTemperature(Fluid, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
        C.printf("\n")
        C.printf(" Current time step: %2.6e s.\n", Integrator_deltaTime)
        C.printf(" Min Flow Temp: %11.6f K. Max Flow Temp: %11.6f K.\n", Flow_minTemperature, Flow_maxTemperature)
        C.printf(" Current number of particles: %d.\n", Particles_number)
        C.printf("\n")
        C.printf("    Iter     Time(s)   Avg Press    Avg Temp      Avg KE  Particle T\n")
      end
      if exitCond or Integrator_timeStep % config.IO.consoleFrequency == 0 then
        var Flow_averagePressure = 0.0
        var Flow_averageTemperature = 0.0
        var Flow_averageKineticEnergy = 0.0
        var Particles_averageTemperature = 0.0
        Flow_averagePressure += CalculateAveragePressure(Fluid, Grid_cellVolume, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
        Flow_averageTemperature += CalculateAverageTemperature(Fluid, Grid_cellVolume, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
        Flow_averageKineticEnergy += CalculateAverageKineticEnergy(Fluid, Grid_cellVolume, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
        Particles_averageTemperature += Particles_IntegrateQuantities(particles)
        Flow_averagePressure = (Flow_averagePressure/(((Grid_xNum*Grid_yNum)*Grid_zNum)*Grid_cellVolume))
        Flow_averageTemperature = (Flow_averageTemperature/(((Grid_xNum*Grid_yNum)*Grid_zNum)*Grid_cellVolume))
        Flow_averageKineticEnergy = (Flow_averageKineticEnergy/(((Grid_xNum*Grid_yNum)*Grid_zNum)*Grid_cellVolume))
        Particles_averageTemperature = (Particles_averageTemperature/Particles_number)
        C.printf("%8d %11.6f %11.6f %11.6f %11.6f %11.6f\n", Integrator_timeStep, Integrator_simTime, Flow_averagePressure, Flow_averageTemperature, Flow_averageKineticEnergy, Particles_averageTemperature)
      end
      if config.IO.wrtRestart then
        if exitCond or Integrator_timeStep % config.IO.restartEveryTimeSteps == 0 then
          var dirname = [&int8](C.malloc(256))
          C.snprintf(dirname, 256, "fluid_sample%d_iter%d", config.Mapping.sampleId, Integrator_timeStep)
          Fluid_dump(primColors, dirname, Fluid, Fluid_copy, Fluid_primPart, Fluid_copy_primPart)
          C.snprintf(dirname, 256, "particles_sample%d_iter%d", config.Mapping.sampleId, Integrator_timeStep)
          particles_dump(primColors, dirname, particles, particles_copy, particles_primPart, particles_copy_primPart)
          C.free([&opaque](dirname))
        end
      end
      for i = 0,config.IO.probes.length do
        var probe = config.IO.probes.values[i]
        if exitCond or Integrator_timeStep % probe.frequency == 0 then
          var temp = Fluid[int3d{probe.coords[0],probe.coords[1],probe.coords[2]}].temperature
          C.fprintf(probeFiles[i], '%d\t%lf\n', Integrator_timeStep, temp)
          C.fflush(probeFiles[i])
        end
      end

      -- Check exit condition after I/O
      if exitCond then
        break
      end
      if (config.Integrator.cfl<0) then
        Integrator_deltaTime = config.Integrator.fixedDeltaTime
      else
        Integrator_maxConvectiveSpectralRadius max= CalculateConvectiveSpectralRadius(Fluid, Flow_gamma, Flow_gasConstant, Grid_dXYZInverseSquare, Grid_xCellWidth, Grid_yCellWidth, Grid_zCellWidth)
        Integrator_maxViscousSpectralRadius max= CalculateViscousSpectralRadius(Fluid, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel, Grid_dXYZInverseSquare)
        Integrator_maxHeatConductionSpectralRadius max= CalculateHeatConductionSpectralRadius(Fluid, Flow_constantVisc, Flow_gamma, Flow_gasConstant, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_prandtl, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel, Grid_dXYZInverseSquare)
        Integrator_deltaTime = (config.Integrator.cfl/max(Integrator_maxConvectiveSpectralRadius, max(Integrator_maxViscousSpectralRadius, Integrator_maxHeatConductionSpectralRadius)))
      end
      Flow_InitializeTemporaries(Fluid)
      Particles_InitializeTemporaries(particles)
      Integrator_time_old = Integrator_simTime
      Integrator_stage = 1
      while (Integrator_stage<5) do
        Flow_InitializeTimeDerivatives(Fluid)
        Particles_InitializeTimeDerivatives(particles)
        Flow_UpdateGhostVelocityGradientStep1(Fluid, BC_xSign, BC_ySign, BC_zSign, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
        Flow_UpdateGhostVelocityGradientStep2(Fluid, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
        Flow_AddGetFlux(Fluid, Flow_constantVisc, Flow_gamma, Flow_gasConstant, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_prandtl, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel, Grid_xBnum, Grid_xCellWidth, Grid_xNum, Grid_yBnum, Grid_yCellWidth, Grid_yNum, Grid_zBnum, Grid_zCellWidth, Grid_zNum)
        Flow_AddUpdateUsingFlux(Fluid, Grid_xBnum, Grid_xCellWidth, Grid_xNum, Grid_yBnum, Grid_yCellWidth, Grid_yNum, Grid_zBnum, Grid_zCellWidth, Grid_zNum)
        Flow_AddBodyForces(Fluid, Flow_bodyForce, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
        if config.Flow.turbForcing then
          Flow_averagePD = 0.0
          Flow_averageDissipation = 0.0
          Flow_averageFe = 0.0
          Flow_averageK = 0.0
          Flow_UpdatePD(Fluid, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
          Flow_ResetDissipation(Fluid)
          Flow_ComputeDissipationX(Fluid, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel, Grid_xBnum, Grid_xCellWidth, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
          Flow_UpdateDissipationX(Fluid, Grid_xBnum, Grid_xCellWidth, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
          Flow_ComputeDissipationY(Fluid, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yCellWidth, Grid_yNum, Grid_zBnum, Grid_zNum)
          Flow_UpdateDissipationY(Fluid, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yCellWidth, Grid_yNum, Grid_zBnum, Grid_zNum)
          Flow_ComputeDissipationZ(Fluid, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zCellWidth, Grid_zNum)
          Flow_UpdateDissipationZ(Fluid, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zCellWidth, Grid_zNum)
          Flow_averagePD += CalculateAveragePD(Fluid, Grid_cellVolume, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
          Flow_averagePD = (Flow_averagePD/(((Grid_xNum*Grid_yNum)*Grid_zNum)*Grid_cellVolume))
          Flow_averageDissipation += CalculateAverageDissipation(Fluid, Grid_cellVolume, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
          Flow_averageDissipation = (Flow_averageDissipation/(((Grid_xNum*Grid_yNum)*Grid_zNum)*Grid_cellVolume))
          Flow_averageK += CalculateAverageK(Fluid, Grid_cellVolume, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
          Flow_averageK = (Flow_averageK/(((Grid_xNum*Grid_yNum)*Grid_zNum)*Grid_cellVolume))
          Flow_averageFe += Flow_AddTurbulentSource(Fluid, Flow_averageDissipation, Flow_averageK, Flow_averagePD, Grid_cellVolume, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
          Flow_averageFe = (Flow_averageFe/(((Grid_xNum*Grid_yNum)*Grid_zNum)*Grid_cellVolume))
          Flow_AdjustTurbulentSource(Fluid, Flow_averageFe, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
        end
        for c in primColors do
          Particles_LocateInCells(particles_primPart[c], BC_xBCPeriodic, BC_yBCPeriodic, BC_zBCPeriodic, Grid_xBnum, Grid_xNum, Grid_xOrigin, Grid_xWidth, Grid_yBnum, Grid_yNum, Grid_yOrigin, Grid_yWidth, Grid_zBnum, Grid_zNum, Grid_zOrigin, Grid_zWidth)
        end
        for c in primColors do
          particles_pushAll(c, particles_primPart[c], particles_qSrcPart_0[c], particles_qSrcPart_1[c], particles_qSrcPart_2[c], particles_qSrcPart_3[c], particles_qSrcPart_4[c], particles_qSrcPart_5[c], particles_qSrcPart_6[c], particles_qSrcPart_7[c], particles_qSrcPart_8[c], particles_qSrcPart_9[c], particles_qSrcPart_10[c], particles_qSrcPart_11[c], particles_qSrcPart_12[c], particles_qSrcPart_13[c], particles_qSrcPart_14[c], particles_qSrcPart_15[c], particles_qSrcPart_16[c], particles_qSrcPart_17[c], particles_qSrcPart_18[c], particles_qSrcPart_19[c], particles_qSrcPart_20[c], particles_qSrcPart_21[c], particles_qSrcPart_22[c], particles_qSrcPart_23[c], particles_qSrcPart_24[c], particles_qSrcPart_25[c], Grid_xNum, Grid_yNum, Grid_zNum, Grid_xBnum, Grid_yBnum, Grid_zBnum, NX, NY, NZ)
        end
        for c in primColors do
          particles_pullAll(c, particles_primPart[c], particles_qDstPart_0[c], particles_qDstPart_1[c], particles_qDstPart_2[c], particles_qDstPart_3[c], particles_qDstPart_4[c], particles_qDstPart_5[c], particles_qDstPart_6[c], particles_qDstPart_7[c], particles_qDstPart_8[c], particles_qDstPart_9[c], particles_qDstPart_10[c], particles_qDstPart_11[c], particles_qDstPart_12[c], particles_qDstPart_13[c], particles_qDstPart_14[c], particles_qDstPart_15[c], particles_qDstPart_16[c], particles_qDstPart_17[c], particles_qDstPart_18[c], particles_qDstPart_19[c], particles_qDstPart_20[c], particles_qDstPart_21[c], particles_qDstPart_22[c], particles_qDstPart_23[c], particles_qDstPart_24[c], particles_qDstPart_25[c])
        end
        Particles_AddFlowCoupling(particles, Fluid, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel, Grid_xCellWidth, Grid_xRealOrigin, Grid_yCellWidth, Grid_yRealOrigin, Grid_zCellWidth, Grid_zRealOrigin, Particles_convectiveCoeff, Particles_heatCapacity)
        Particles_AddBodyForces(particles, Particles_bodyForce)
        if (config.Radiation.type == SCHEMA.RadiationType_Algebraic) then
          AddRadiation(particles, config)
        end
        if (config.Radiation.type == SCHEMA.RadiationType_DOM) then
          Radiation_ClearAccumulators(Radiation)
          for c in primColors do
            Radiation_AccumulateParticleValues(particles_primPart[c], Fluid_primPart[c], Radiation_primPart[c])
          end
          Radiation_UpdateFieldValues(Radiation, Radiation_cellVolume, Radiation_qa, Radiation_qs);
          [DOM.ComputeRadiationField(config, primColors, Radiation_primPart)];
          for c in primColors do
            Particles_AbsorbRadiation(particles_primPart[c], Fluid_primPart[c], Radiation_primPart[c], Particles_heatCapacity, Radiation_qa)
          end
        end
        Flow_AddParticlesCoupling(particles, Fluid, Grid_cellVolume)
        Flow_UpdateVars(Fluid, Integrator_deltaTime, Integrator_stage)
        Particles_UpdateVars(particles, Integrator_deltaTime, Integrator_stage)
        Flow_UpdateAuxiliaryVelocity(Fluid, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
        Flow_UpdateGhostConservedStep1(Fluid, BC_xNegTemperature, BC_xNegVelocity, BC_xPosTemperature, BC_xPosVelocity, BC_xSign, BC_yNegTemperature, BC_yNegVelocity, BC_yPosTemperature, BC_yPosVelocity, BC_ySign, BC_zNegTemperature, BC_zNegVelocity, BC_zPosTemperature, BC_zPosVelocity, BC_zSign, Flow_gamma, Flow_gasConstant, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
        Flow_UpdateGhostConservedStep2(Fluid, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
        Flow_UpdateGhostVelocityStep1(Fluid, BC_xNegVelocity, BC_xPosVelocity, BC_xSign, BC_yNegVelocity, BC_yPosVelocity, BC_ySign, BC_zNegVelocity, BC_zPosVelocity, BC_zSign, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
        Flow_UpdateGhostVelocityStep2(Fluid, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
        Flow_ComputeVelocityGradientAll(Fluid, Grid_xBnum, Grid_xCellWidth, Grid_xNum, Grid_yBnum, Grid_yCellWidth, Grid_yNum, Grid_zBnum, Grid_zCellWidth, Grid_zNum)
        Flow_UpdateAuxiliaryThermodynamics(Fluid, Flow_gamma, Flow_gasConstant, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
        Flow_UpdateGhostThermodynamicsStep1(Fluid, BC_xNegTemperature, BC_xPosTemperature, BC_yNegTemperature, BC_yPosTemperature, BC_zNegTemperature, BC_zPosTemperature, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
        Flow_UpdateGhostThermodynamicsStep2(Fluid, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
        Particles_UpdateAuxiliaryStep1(particles, BC_xBCLeftParticles, BC_xBCRightParticles, BC_yBCLeftParticles, BC_yBCRightParticles, BC_zBCLeftParticles, BC_zBCRightParticles, Grid_xOrigin, Grid_xWidth, Grid_yOrigin, Grid_yWidth, Grid_zOrigin, Grid_zWidth, Particles_restitutionCoeff)
        Particles_UpdateAuxiliaryStep2(particles)
        Integrator_simTime = (Integrator_time_old+((double(0.5)*(1+(Integrator_stage/3)))*Integrator_deltaTime))
        Integrator_stage = (Integrator_stage+1)
        for c in primColors do
          Particles_number += Particles_DeleteEscapingParticles(particles_primPart[c], Grid_xRealOrigin, Grid_xRealWidth, Grid_yRealOrigin, Grid_yRealWidth, Grid_zRealOrigin, Grid_zRealWidth)
        end
      end
      Integrator_timeStep = (Integrator_timeStep+1)
    end

    -- Close long-running files
    for i = 0,config.IO.probes.length do
      C.fclose(probeFiles[i])
    end
    C.free(probeFiles)

  end
end

task main()
  var args = regentlib.c.legion_runtime_get_input_args()
  var launched = 0
  for i = 1, args.argc do
    if C.strcmp(args.argv[i],"-i") == 0 and i < args.argc-1 then
      var config = SCHEMA.parse_config(args.argv[i+1])
      config.Mapping.sampleId = launched
      work(config)
      launched += 1
    end
  end
  if launched < 1 then
    var stderr = C.fdopen(2, 'w')
    C.fprintf(stderr, "No testcases supplied.\n")
    C.fprintf(stderr, "Usage: %s -i config1.json [-i config2.json ...]\n", args.argv[0])
    C.fflush(stderr)
    C.exit(1)
  end
end

-------------------------------------------------------------------------------
-- COMPILATION CALL
-------------------------------------------------------------------------------

regentlib.saveobj(main, "soleil.o", "object", MAPPER.register_mappers)
