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

local ceil = regentlib.ceil(double)
local floor = regentlib.floor(double)
local fmod = regentlib.fmod(double)
local pow = regentlib.pow(double)

-------------------------------------------------------------------------------
-- COMPILE-TIME CONFIGURATION
-------------------------------------------------------------------------------

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
  velocity_t : double[3];
  temperature_t : double;
  __valid : bool;
  __xfer_dir : int8;
  __xfer_slot : int;
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
  'velocity_t',
  'temperature_t',
})

local Particles_subStepConserved =
  UTIL.setToList(
    UTIL.setSubList(
      UTIL.listToSet(UTIL.fieldNames(Particles_columns)),
      Particles_subStepTemp))

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
  rhoVelocity : double[3];
  rhoEnergy : double;
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
  dissipation : double;
  dissipationFlux : double;
  dudtBoundary : double;
  dTdtBoundary : double;
  velocity_old_NSCBC : double[3];
  temperature_old_NSCBC : double;
  velocity_inc : double[3];
  temperature_inc : double;
}

-------------------------------------------------------------------------------
-- CONSTANTS
-------------------------------------------------------------------------------

local PI = 3.1415926535898

-------------------------------------------------------------------------------
-- MACROS
-------------------------------------------------------------------------------

__demand(__inline)
task vs_mul(a : double[3], b : double)
  return array(a[0] * b, a[1] * b, a[2] * b)
end

__demand(__inline)
task vs_div(a : double[3], b : double)
  return array(a[0] / b, a[1] / b, a[2] / b)
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

__demand(__cuda) -- MANUALLY PARALLELIZED
task Particles_InitializeUniform(Particles : region(ispace(int1d), Particles_columns),
                                 Fluid : region(ispace(int3d), Fluid_columns),
                                 config : Config,
                                 Grid_xBnum : int32, Grid_yBnum : int32, Grid_zBnum : int32)
where
  reads(Fluid.{centerCoordinates, velocity}),
  writes(Particles.{__valid, cell, position, velocity, density, temperature, diameter})
do
  var pBase = Particles.bounds.lo
  var lo = Fluid.bounds.lo
  lo.x = max(lo.x, Grid_xBnum)
  lo.y = max(lo.y, Grid_yBnum)
  lo.z = max(lo.z, Grid_zBnum)
  var hi = Fluid.bounds.hi
  hi.x = min(hi.x, config.Grid.xNum+Grid_xBnum-1)
  hi.y = min(hi.y, config.Grid.yNum+Grid_yBnum-1)
  hi.z = min(hi.z, config.Grid.zNum+Grid_zBnum-1)
  var xSize = hi.x-lo.x+1
  var ySize = hi.y-lo.y+1
  var zSize = hi.z-lo.z+1
  var particlesPerTask = config.Particles.initNum / (config.Mapping.tiles[0]*config.Mapping.tiles[1]*config.Mapping.tiles[2])
  var Particles_density = config.Particles.density
  var Particles_initTemperature = config.Particles.initTemperature
  var Particles_diameterMean = config.Particles.diameterMean
  __demand(__openmp)
  for p in Particles do
    var relIdx = int(p - pBase)
    if relIdx < particlesPerTask then
      p.__valid = true
      var c = lo + int3d{relIdx%xSize, relIdx/xSize%ySize, relIdx/xSize/ySize%zSize}
      p.cell = c
      p.position = Fluid[c].centerCoordinates
      p.velocity = Fluid[c].velocity
      p.density = Particles_density
      p.temperature = Particles_initTemperature
      p.diameter = Particles_diameterMean
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
task Flow_InitializeCell(Fluid : region(ispace(int3d), Fluid_columns))
where
  writes(Fluid.centerCoordinates),
  writes(Fluid.dissipation),
  writes(Fluid.dissipationFlux),
  writes(Fluid.pressure),
  writes(Fluid.rho),
  writes(Fluid.rhoEnergy),
  writes(Fluid.{rhoEnergyFluxX, rhoEnergyFluxY, rhoEnergyFluxZ}),
  writes(Fluid.rhoEnergy_new),
  writes(Fluid.rhoEnergy_old),
  writes(Fluid.rhoEnergy_t),
  writes(Fluid.rhoEnthalpy),
  writes(Fluid.rhoFluxX),
  writes(Fluid.rhoFluxY),
  writes(Fluid.rhoFluxZ),
  writes(Fluid.rhoVelocity),
  writes(Fluid.{rhoVelocityFluxX, rhoVelocityFluxY, rhoVelocityFluxZ}),
  writes(Fluid.rhoVelocity_new),
  writes(Fluid.rhoVelocity_old),
  writes(Fluid.rhoVelocity_t),
  writes(Fluid.rho_new),
  writes(Fluid.rho_old),
  writes(Fluid.rho_t),
  writes(Fluid.temperature),
  writes(Fluid.velocity),
  writes(Fluid.{velocityGradientX, velocityGradientY, velocityGradientZ}),
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
    Fluid[c].rhoVelocity = array(0.0, 0.0, 0.0)
    Fluid[c].rhoEnergy = 0.0
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
        Particles[p1].cell = cell
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
  reads writes(Particles.{velocity_t, temperature_t}),
  writes(Particles.{deltaTemperatureTerm, deltaVelocityOverRelaxationTime})
do
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      var flowVelocity = InterpolateTriVelocity(Particles[p].cell,
                                                Particles[p].position,
                                                Fluid,
                                                Grid_xCellWidth, Grid_xRealOrigin,
                                                Grid_yCellWidth, Grid_yRealOrigin,
                                                Grid_zCellWidth, Grid_zRealOrigin)
      var flowTemperature = InterpolateTriTemp(Particles[p].cell,
                                               Particles[p].position,
                                               Fluid,
                                               Grid_xCellWidth, Grid_xRealOrigin,
                                               Grid_yCellWidth, Grid_yRealOrigin,
                                               Grid_zCellWidth, Grid_zRealOrigin)
      var flowDynamicViscosity = GetDynamicViscosity(flowTemperature,
                                                     Flow_constantVisc,
                                                     Flow_powerlawTempRef, Flow_powerlawViscRef,
                                                     Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef,
                                                     Flow_viscosityModel)
      var particleReynoldsNumber = 0.0
      var relaxationTime =
        Particles[p].density
        * pow(Particles[p].diameter,2.0)
        / (18.0*flowDynamicViscosity)
        / (1.0 + (0.15*pow(particleReynoldsNumber,0.687)))
      var tmp2 = vs_div(vv_sub(flowVelocity, Particles[p].velocity), relaxationTime)
      Particles[p].deltaVelocityOverRelaxationTime = tmp2
      Particles[p].velocity_t = vv_add(Particles[p].velocity_t, tmp2)
      var tmp3 = PI * pow(Particles[p].diameter,2.0) * Particles_convectiveCoeff * (flowTemperature-Particles[p].temperature)
      Particles[p].deltaTemperatureTerm = tmp3
      Particles[p].temperature_t += tmp3/(PI*pow(Particles[p].diameter,3.0)/6.0*Particles[p].density*Particles_heatCapacity)
    end
  end
end

-------------------------------------------------------------------------------
-- MAIN SIMULATION
-------------------------------------------------------------------------------

local function mkInstance() local INSTANCE = {}

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
  local Particles = regentlib.newsymbol()
  local tiles = regentlib.newsymbol()
  local p_Fluid = regentlib.newsymbol()
  local p_Particles = regentlib.newsymbol()

  -----------------------------------------------------------------------------
  -- Exported symbols
  -----------------------------------------------------------------------------

  INSTANCE.Grid = Grid
  INSTANCE.Integrator_deltaTime = Integrator_deltaTime
  INSTANCE.Integrator_exitCond = Integrator_exitCond
  INSTANCE.Fluid = Fluid
  INSTANCE.Particles = Particles
  INSTANCE.tiles = tiles
  INSTANCE.p_Fluid = p_Fluid
  INSTANCE.p_Particles = p_Particles

  -----------------------------------------------------------------------------
  -- Symbol declaration & initialization
  -----------------------------------------------------------------------------

  function INSTANCE.DeclSymbols(config) return rquote

    ---------------------------------------------------------------------------
    -- Declare & initialize state variables
    ---------------------------------------------------------------------------

    -- Cell step size (TODO: Change when we go to non-uniform meshes)
    var [Grid.xCellWidth] = config.Grid.xWidth / config.Grid.xNum
    var [Grid.yCellWidth] = config.Grid.yWidth / config.Grid.yNum
    var [Grid.zCellWidth] = config.Grid.zWidth / config.Grid.zNum
    var [Grid.cellVolume] = Grid.xCellWidth * Grid.yCellWidth * Grid.zCellWidth

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

    -- Create Particles Regions
    regentlib.assert(config.Particles.maxNum % numTiles == 0,
                     'Uneven partitioning of particles')
    var maxParticlesPerTile = ceil((config.Particles.maxNum / numTiles) * config.Particles.maxSkew)
    var is_Particles = ispace(int1d, maxParticlesPerTile * numTiles)
    var [Particles] = region(is_Particles, Particles_columns);
    [UTIL.emitRegionTagAttach(Particles, MAPPER.SAMPLE_ID_TAG, sampleId, int)];

    -- Partitioning domain
    var [tiles] = ispace(int3d, {NX,NY,NZ})

    -- Fluid Partitioning
    var [p_Fluid] =
      [UTIL.mkPartitionEqually(int3d, int3d, Fluid_columns)]
      (Fluid, tiles, Grid.xBnum, Grid.yBnum, Grid.zBnum)

    -- Particles Partitioning
    var [p_Particles] =
      [UTIL.mkPartitionEqually(int1d, int3d, Particles_columns)]
      (Particles, tiles, 0)

    ---------------------------------------------------------------------------
    -- DOM code declarations
    ---------------------------------------------------------------------------

  end end -- DeclSymbols

  -----------------------------------------------------------------------------
  -- Region initialization
  -----------------------------------------------------------------------------

  function INSTANCE.InitRegions(config) return rquote

    Particles_initValidField(Particles)
    Flow_InitializeCell(Fluid)
    Flow_InitializeCenterCoordinates(Fluid,
                                     Grid.xBnum, config.Grid.xNum, config.Grid.origin[0], config.Grid.xWidth,
                                     Grid.yBnum, config.Grid.yNum, config.Grid.origin[1], config.Grid.yWidth,
                                     Grid.zBnum, config.Grid.zNum, config.Grid.origin[2], config.Grid.zWidth)

    -- Initialize particles
    if config.Particles.initCase == SCHEMA.ParticlesInitCase_Random then
      regentlib.assert(false, "Random particle initialization is disabled")
    elseif config.Particles.initCase == SCHEMA.ParticlesInitCase_Uniform then
      regentlib.assert(config.Particles.initNum <= config.Particles.maxNum,
                       "Not enough space for initial number of particles")
      for c in tiles do
        Particles_InitializeUniform(p_Particles[c],
                                    p_Fluid[c],
                                    config,
                                    Grid.xBnum, Grid.yBnum, Grid.zBnum)
      end
      Particles_number = (config.Particles.initNum / numTiles) * numTiles
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
    end

  end end -- MainLoopHeader

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

      -- Add fluid forces to particles
      Particles_AddFlowCoupling(Particles, Fluid, config.Flow.constantVisc, config.Flow.powerlawTempRef, config.Flow.powerlawViscRef, config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef, config.Flow.viscosityModel, Grid.xCellWidth, Grid.xRealOrigin, Grid.yCellWidth, Grid.yRealOrigin, Grid.zCellWidth, Grid.zRealOrigin, config.Particles.convectiveCoeff, config.Particles.heatCapacity)

    end -- RK4 sub-time-stepping

    Integrator_timeStep += 1

  end end -- MainLoopBody

return INSTANCE end -- mkInstance

-------------------------------------------------------------------------------
-- TOP-LEVEL INTERFACE
-------------------------------------------------------------------------------

local CopyQueue = regentlib.newsymbol()
local FakeCopyQueue = regentlib.newsymbol()

local function parallelizeFor(sim, stmts)
  return rquote
    __parallelize_with
      sim.p_Fluid, sim.p_Particles, sim.tiles,
      image(sim.Fluid, sim.p_Particles, [sim.Particles].cell) <= sim.p_Fluid
    do [stmts] end
  end
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
  var timestep = 0
  while true do
    -- Perform preliminary actions before each timestep
    [parallelizeFor(SIM0, SIM0.MainLoopHeader(rexpr mc.configs[0] end))];
    [parallelizeFor(SIM1, SIM1.MainLoopHeader(rexpr mc.configs[1] end))];
    -- Make sure both simulations are using the same timestep
    SIM0.Integrator_deltaTime = min(SIM0.Integrator_deltaTime, SIM1.Integrator_deltaTime)
    SIM1.Integrator_deltaTime = min(SIM0.Integrator_deltaTime, SIM1.Integrator_deltaTime);
    if SIM0.Integrator_exitCond or SIM1.Integrator_exitCond then
      break
    end
    -- Run 1 iteration of first section
    [parallelizeFor(SIM0, SIM0.MainLoopBody(rexpr mc.configs[0] end, FakeCopyQueue))];
    -- Copy fluid & particles to second section
    fill(CopyQueue.__valid, false) -- clear the copyqueue from the previous iteration
    if timestep % mc.copyEveryTimeSteps == 0 then
      fill(CopyQueue.position, array(-1.0, -1.0, -1.0))
      fill(CopyQueue.velocity, array(-1.0, -1.0, -1.0))
      fill(CopyQueue.temperature, -1.0)
      fill(CopyQueue.diameter, -1.0)
      fill(CopyQueue.density, -1.0)
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
    timestep += 1
  end
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
