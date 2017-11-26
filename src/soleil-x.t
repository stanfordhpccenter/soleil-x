-----------------------------------------------------------------------------
--[[
-----------------------------------------------------------------------------

Soleil-X Version 0.0.1
Copyright (C) 2013-2015, Dr. Thomas D. Economon,
                         Dr. Ivan Bermejo-Moreno

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public
License along with this program; if not, write to the Free
Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
nBoston, MA 02110-1301 USA.

-----------------------------------------------------------------------------
]]--
-----------------------------------------------------------------------------

import 'ebb'

local A    = require 'admiral'
local GRID = require 'ebb.domains.grid'
local JSON = require 'json'
local L    = require 'ebblib'
local M    = require 'ebb.src.main'

-----------------------------------------------------------------------------
--[[                          MATH IMPORTS                               ]]--
-----------------------------------------------------------------------------

local C = terralib.includecstring [[
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
]]

----------------------------------------------------------------------------
--[[                      COMPILE-TIME OPTIONS                           ]]--
-----------------------------------------------------------------------------

local function PrintUsageAndExit ()
  print("Usage : $LISZT_PATH/liszt-legion.sh $SOLEIL_PATH/src/soleil-x.t")
  print("          -i <config.json> (** required **)")
  os.exit(1)
end

do
  local configFileName
  local i = 1
  while i <= #arg do
    if arg[i] == '-i' then
      if i == #arg then PrintUsageAndExit() end
      configFileName = arg[i+1]
      break
    end
    i = i + 1
  end
  if not configFileName then PrintUsageAndExit() end
  local f = io.open(configFileName, 'r')
  if not f then PrintUsageAndExit() end
  local content = f:read('*all')
  if not content then PrintUsageAndExit() end
  f:close()
  config = JSON.decode(content)
  if not config then PrintUsageAndExit() end
end

-----------------------------------------------------------------------------
--[[                            CONSTANTS                                ]]--
-----------------------------------------------------------------------------

local pi = 2.0*L.acos(0)
local twoPi = 2.0*pi
local SB = 5.67e-8

-----------------------------------------------------------------------------
--[[                            NAMESPACES                               ]]--
-----------------------------------------------------------------------------

local Grid = {}
local BC = {}
local Integrator = {}
local Flow = {}
local Particles = {}
local Radiation = {}
local Statistics = {}
local IO = {}

-----------------------------------------------------------------------------
--[[                         OPTIONS & GLOBALS                           ]]--
-----------------------------------------------------------------------------

---------------------------------
--[[ Grid (without boundary) ]]--
---------------------------------

-- Number of cells in the x, y, & z directions
Grid.xNum = A.globalFromConfig('Grid.xNum', int)
Grid.yNum = A.globalFromConfig('Grid.yNum', int)
Grid.zNum = A.globalFromConfig('Grid.zNum', int)
-- Number of tiles in each direction
Grid.xTiles = A.readConfig('Grid.xTiles', int)
Grid.yTiles = A.readConfig('Grid.yTiles', int)
Grid.zTiles = A.readConfig('Grid.zTiles', int)
Grid.numTiles = Grid.xTiles * Grid.yTiles * Grid.zTiles
-- Origin of the computational domain (meters)
Grid.origin = A.readConfig('Grid.origin', double[3])
Grid.xOrigin = L.Global('Grid.xOrigin', L.double, Grid.origin[0])
Grid.yOrigin = L.Global('Grid.yOrigin', L.double, Grid.origin[1])
Grid.zOrigin = L.Global('Grid.zOrigin', L.double, Grid.origin[2])
-- Width of the computational domain in the x, y, & z directions (meters)
Grid.xWidth = A.globalFromConfig('Grid.xWidth', double)
Grid.yWidth = A.globalFromConfig('Grid.yWidth', double)
Grid.zWidth = A.globalFromConfig('Grid.zWidth', double)
-- Width of each cell in the x, y, & z directions (meters)
Grid.xCellWidth = L.Global('Grid.xCellWidth', L.double, Grid.xWidth:get() / Grid.xNum:get())
Grid.yCellWidth = L.Global('Grid.yCellWidth', L.double, Grid.yWidth:get() / Grid.yNum:get())
Grid.zCellWidth = L.Global('Grid.zCellWidth', L.double, Grid.zWidth:get() / Grid.zNum:get())
-- Cell volume
Grid.cellVolume =
  L.Global('Grid.cellVolume', L.double,
           Grid.xCellWidth:get() * Grid.yCellWidth:get() * Grid.zCellWidth:get())
Grid.areaInterior =
  Grid.xNum:get() * Grid.yNum:get() * Grid.zNum:get() * Grid.cellVolume:get()
Grid.dXYZInverseSquare =
  L.Global('Grid.dXYZInverseSquare', L.double,
           1.0/Grid.xCellWidth:get() * 1.0/Grid.xCellWidth:get() +
           1.0/Grid.yCellWidth:get() * 1.0/Grid.yCellWidth:get() +
           1.0/Grid.zCellWidth:get() * 1.0/Grid.zCellWidth:get())

-----------------------------
--[[ Boundary Conditions ]]--
-----------------------------

-- Boundary conditions for each face of the block and possible
-- wall velocity, if no-slip.

local FlowBC = A.Enum('Periodic','Symmetry','AdiabaticWall','IsothermalWall')

-- X
BC.xBCLeft         = A.readConfig('BC.xBCLeft', FlowBC)
BC.xBCLeftVel      = A.readConfig('BC.xBCLeftVel', double[3])
BC.xBCLeftTemp     = A.readConfig('BC.xBCLeftTemp', double)
BC.xBCRight        = A.readConfig('BC.xBCRight', FlowBC)
BC.xBCRightVel     = A.readConfig('BC.xBCRightVel', double[3])
BC.xBCRightTemp    = A.readConfig('BC.xBCRightTemp', double)
BC.xBCPeriodic     = L.Global('BC.xBCPeriodic', L.bool,
                              M.COND2EXPR(M.EQ(BC.xBCLeft, FlowBC.Periodic)))
BC.xSign           = L.Global('BC.xSign', L.vec3d, {0.0,0.0,0.0})
BC.xPosVelocity    = L.Global('BC.xPosVelocity', L.vec3d, {0.0,0.0,0.0})
BC.xNegVelocity    = L.Global('BC.xNegVelocity', L.vec3d, {0.0,0.0,0.0})
BC.xPosTemperature = L.Global('BC.xPosTemperature', L.double, 0.0)
BC.xNegTemperature = L.Global('BC.xNegTemperature', L.double, 0.0)

-- Y
BC.yBCLeft         = A.readConfig('BC.yBCLeft', FlowBC)
BC.yBCLeftVel      = A.readConfig('BC.yBCLeftVel', double[3])
BC.yBCLeftTemp     = A.readConfig('BC.yBCLeftTemp', double)
BC.yBCRight        = A.readConfig('BC.yBCRight', FlowBC)
BC.yBCRightVel     = A.readConfig('BC.yBCRightVel', double[3])
BC.yBCRightTemp    = A.readConfig('BC.yBCRightTemp', double)
BC.yBCPeriodic     = L.Global('BC.yBCPeriodic', L.bool,
                              M.COND2EXPR(M.EQ(BC.yBCLeft, FlowBC.Periodic)))
BC.ySign           = L.Global('BC.ySign', L.vec3d, {0.0,0.0,0.0})
BC.yPosVelocity    = L.Global('BC.yPosVelocity', L.vec3d, {0.0,0.0,0.0})
BC.yNegVelocity    = L.Global('BC.yNegVelocity', L.vec3d, {0.0,0.0,0.0})
BC.yPosTemperature = L.Global('BC.yPosTemperature', L.double, 0.0)
BC.yNegTemperature = L.Global('BC.yNegTemperature', L.double, 0.0)

-- Z
BC.zBCLeft         = A.readConfig('BC.zBCLeft', FlowBC)
BC.zBCLeftVel      = A.readConfig('BC.zBCLeftVel', double[3])
BC.zBCLeftTemp     = A.readConfig('BC.zBCLeftTemp', double)
BC.zBCRight        = A.readConfig('BC.zBCRight', FlowBC)
BC.zBCRightVel     = A.readConfig('BC.zBCRightVel', double[3])
BC.zBCRightTemp    = A.readConfig('BC.zBCRightTemp', double)
BC.zBCPeriodic     = L.Global('BC.zBCPeriodic', L.bool,
                              M.COND2EXPR(M.EQ(BC.zBCLeft, FlowBC.Periodic)))
BC.zSign           = L.Global('BC.zSign', L.vec3d, {0.0,0.0,0.0})
BC.zPosVelocity    = L.Global('BC.zPosVelocity', L.vec3d, {0.0,0.0,0.0})
BC.zNegVelocity    = L.Global('BC.zNegVelocity', L.vec3d, {0.0,0.0,0.0})
BC.zPosTemperature = L.Global('BC.zPosTemperature', L.double, 0.0)
BC.zNegTemperature = L.Global('BC.zNegTemperature', L.double, 0.0)

-- Particle boundary conditions

local ParticleBC = A.Enum('Permeable','Solid')

BC.xBCLeftParticles  = L.Global('BC.xBCLeftParticles', L.int, -1)
BC.xBCRightParticles = L.Global('BC.xBCRightParticles', L.int, -1)
BC.yBCLeftParticles  = L.Global('BC.yBCLeftParticles', L.int, -1)
BC.yBCRightParticles = L.Global('BC.yBCRightParticles', L.int, -1)
BC.zBCLeftParticles  = L.Global('BC.zBCLeftParticles', L.int, -1)
BC.zBCRightParticles = L.Global('BC.zBCRightParticles', L.int, -1)

------------------------------
--[[ Grid (with boundary) ]]--
------------------------------

-- Number of boundary cells
Grid.xBnum = L.Global('Grid.xBnum', L.int,
                      M.MIN(1, (BC.xBCLeft - FlowBC.Periodic) *
                               (BC.xBCLeft - FlowBC.Periodic)))
Grid.yBnum = L.Global('Grid.yBnum', L.int,
                      M.MIN(1, (BC.yBCLeft - FlowBC.Periodic) *
                               (BC.yBCLeft - FlowBC.Periodic)))
Grid.zBnum = L.Global('Grid.zBnum', L.int,
                      M.MIN(1, (BC.zBCLeft - FlowBC.Periodic) *
                               (BC.zBCLeft - FlowBC.Periodic)))
-- Boundary width
Grid.xBwidth = Grid.xCellWidth:get() * Grid.xBnum:get()
Grid.yBwidth = Grid.yCellWidth:get() * Grid.yBnum:get()
Grid.zBwidth = Grid.zCellWidth:get() * Grid.zBnum:get()
-- Number of cells
Grid.xRealNum = Grid.xNum:get() + 2*Grid.xBnum:get()
Grid.yRealNum = Grid.yNum:get() + 2*Grid.yBnum:get()
Grid.zRealNum = Grid.zNum:get() + 2*Grid.zBnum:get()
-- Origin
Grid.xRealOrigin = L.Global('Grid.xRealOrigin', L.double, Grid.xOrigin:get() - Grid.xBwidth)
Grid.yRealOrigin = L.Global('Grid.yRealOrigin', L.double, Grid.yOrigin:get() - Grid.yBwidth)
Grid.zRealOrigin = L.Global('Grid.zRealOrigin', L.double, Grid.zOrigin:get() - Grid.zBwidth)
-- Width
Grid.xRealWidth = L.Global('Grid.xRealWidth', L.double, Grid.xWidth:get() + 2*Grid.xBwidth)
Grid.yRealWidth = L.Global('Grid.yRealWidth', L.double, Grid.yWidth:get() + 2*Grid.yBwidth)
Grid.zRealWidth = L.Global('Grid.zRealWidth', L.double, Grid.zWidth:get() + 2*Grid.zBwidth)

--------------------
--[[ Integrator ]]--
--------------------

-- Spatial integration
Integrator.SPLIT = 0.5 -- Splitting parameter

-- Time integration
Integrator.finalTime      = A.readConfig('Integrator.finalTime', double)
Integrator.restartIter    = A.readConfig('Integrator.restartIter', int)
Integrator.maxIter        = A.readConfig('Integrator.maxIter', int)
Integrator.cfl            = A.readConfig('Integrator.cfl', double)
Integrator.fixedDeltaTime = A.readConfig('Integrator.fixedDeltaTime', double)
Integrator.simTime        = L.Global('Integrator.simTime', L.double, 0.0)
Integrator.time_old       = L.Global('Integrator.time_old', L.double, 0.0)
Integrator.timeStep       = L.Global('Integrator.timeStep', L.int, 0)
Integrator.deltaTime      = L.Global('Integrator.deltaTime', L.double, 0.0001)
Integrator.stage          = L.Global('Integrator.stage', L.int, 0)

-- Spectral radii for cfl-based delta time
Integrator.maxConvectiveSpectralRadius =
  L.Global('Integrator.maxConvectiveSpectralRadius', L.double, 0.0)
Integrator.maxViscousSpectralRadius =
  L.Global('Integrator.maxViscousSpectralRadius', L.double, 0.0)
Integrator.maxHeatConductionSpectralRadius =
  L.Global('Integrator.maxHeatConductionSpectralRadius', L.double, 0.0)

--------------
--[[ Flow ]]--
--------------

Flow.gasConstant = A.globalFromConfig('Flow.gasConstant', double)
Flow.gamma       = A.globalFromConfig('Flow.gamma', double)
Flow.prandtl     = A.globalFromConfig('Flow.prandtl', double)

local ViscosityModel = A.Enum('Constant','PowerLaw','Sutherland')

Flow.viscosityModel    = A.globalFromConfig('Flow.viscosityModel', ViscosityModel)
Flow.constantVisc      = A.globalFromConfig('Flow.constantVisc', double)
Flow.powerlawViscRef   = A.globalFromConfig('Flow.powerlawViscRef', double)
Flow.powerlawTempRef   = A.globalFromConfig('Flow.powerlawTempRef', double)
Flow.sutherlandViscRef = A.globalFromConfig('Flow.sutherlandViscRef', double)
Flow.sutherlandTempRef = A.globalFromConfig('Flow.sutherlandTempRef', double)
Flow.sutherlandSRef    = A.globalFromConfig('Flow.sutherlandSRef', double)

local FlowInitCase = A.Enum('Uniform','Restart','Perturbed',
                            'TaylorGreen2DVortex','TaylorGreen3DVortex')
local OnOrOff = A.Enum('OFF','ON')

Flow.initCase    = A.readConfig('Flow.initCase', FlowInitCase)
Flow.initParams  = A.globalFromConfig('Flow.initParams', double[5])
Flow.bodyForce   = A.globalFromConfig('Flow.bodyForce', double[3])
Flow.turbForcing = A.readConfig('turbForcing', OnOrOff)

Flow.averagePressure      = L.Global('Flow.averagePressure', L.double, 0.0)
Flow.averageTemperature   = L.Global('Flow.averageTemperature', L.double, 0.0)
Flow.averageKineticEnergy = L.Global('Flow.averageKineticEnergy', L.double, 0.0)
Flow.minTemperature       = L.Global('Flow.minTemperature', L.double, 0.0)
Flow.maxTemperature       = L.Global('Flow.maxTemperature', L.double, 0.0)
Flow.averagePD            = L.Global('Flow.averagePD', L.double, 0.0)
Flow.averageDissipation   = L.Global('Flow.averageDissipation', L.double, 0.0)
Flow.averageFe            = L.Global('Flow.averageFe', L.double, 0.0)
Flow.averageK             = L.Global('Flow.averageK', L.double, 0.0)

-------------------
--[[ Particles ]]--
-------------------

local ParticlesInitCase = A.Enum('Random','Restart','Uniform')

-- Define the initial number of particles and insertion/deletion
Particles.initCase           = A.readConfig('Particles.initCase', ParticlesInitCase)
Particles.initNum            = A.readConfig('Particles.initNum', int)
Particles.maxNum             = A.globalFromConfig('Particles.maxNum', int)
Particles.restitutionCoeff   = A.globalFromConfig('Particles.restitutionCoeff', double)
Particles.convectiveCoeff    = A.globalFromConfig('Particles.convectiveCoeff', double)
Particles.absorptivity       = A.readConfig('Particles.absorptivity', double)
Particles.heatCapacity       = A.globalFromConfig('Particles.heatCapacity', double)
Particles.initTemperature    = A.readConfig('Particles.initTemperature', double)
Particles.density            = A.globalFromConfig('Particles.density', double)
Particles.diameterMean       = A.readConfig('Particles.diameterMean', double)
Particles.bodyForce          = A.globalFromConfig('Particles.bodyForce', double[3])
Particles.maxSkew            = A.globalFromConfig('Particles.maxSkew', double)
Particles.maxXferNum         = A.globalFromConfig('Particles.maxXferNum', int)
Particles.averageTemperature = L.Global('Particles.averageTemperature', L.double, 0.0)
Particles.number             = L.Global('Particles.number', L.int64, 0)

-------------------
--[[ Radiation ]]--
-------------------

Radiation.TYPE = config.Radiation.TYPE
assert(Radiation.TYPE == 'OFF' or
       Radiation.TYPE == 'Algebraic' or
       Radiation.TYPE == 'DOM')

if Radiation.TYPE == 'Algebraic' then
  Radiation.intensity = A.readConfig('Radiation.intensity', double)
end

if Radiation.TYPE == 'DOM' then
  Radiation.NUM_ANGLES = config.Radiation.NUM_ANGLES
  Radiation.qa         = A.globalFromConfig('Radiation.qa', double)
  Radiation.qs         = A.globalFromConfig('Radiation.qs', double)
  Radiation.xNum       = A.globalFromConfig('Radiation.xNum', int)
  Radiation.yNum       = A.globalFromConfig('Radiation.yNum', int)
  Radiation.zNum       = A.globalFromConfig('Radiation.zNum', int)
  Radiation.xBnum      = L.Global('Radiation.xBnum', L.int, 0)
  Radiation.yBnum      = L.Global('Radiation.yBnum', L.int, 0)
  Radiation.zBnum      = L.Global('Radiation.zBnum', L.int, 0)
  Radiation.xPeriodic  = L.Global('Radiation.xPeriodic', L.bool, false)
  Radiation.yPeriodic  = L.Global('Radiation.yPeriodic', L.bool, false)
  Radiation.zPeriodic  = L.Global('Radiation.zPeriodic', L.bool, false)
  Radiation.xCellWidth = L.Global('Radiation.xCellWidth', L.double,
                                  Grid.xWidth:get() / Radiation.xNum:get())
  Radiation.yCellWidth = L.Global('Radiation.yCellWidth', L.double,
                                  Grid.yWidth:get() / Radiation.yNum:get())
  Radiation.zCellWidth = L.Global('Radiation.zCellWidth', L.double,
                                  Grid.zWidth:get() / Radiation.zNum:get())
  Radiation.cellVolume = L.Global('Radiation.cellVolume', L.double,
                                  Radiation.xCellWidth:get() *
                                  Radiation.yCellWidth:get() *
                                  Radiation.zCellWidth:get())
  Radiation.emissEast  = A.readConfig('Radiation.emissEast', double)
  Radiation.emissWest  = A.readConfig('Radiation.emissWest', double)
  Radiation.emissSouth = A.readConfig('Radiation.emissSouth', double)
  Radiation.emissNorth = A.readConfig('Radiation.emissNorth', double)
  Radiation.emissUp    = A.readConfig('Radiation.emissUp', double)
  Radiation.emissDown  = A.readConfig('Radiation.emissDown', double)
  Radiation.tempEast  = A.readConfig('Radiation.tempEast', double)
  Radiation.tempWest  = A.readConfig('Radiation.tempWest', double)
  Radiation.tempSouth = A.readConfig('Radiation.tempSouth', double)
  Radiation.tempNorth = A.readConfig('Radiation.tempNorth', double)
  Radiation.tempUp    = A.readConfig('Radiation.tempUp', double)
  Radiation.tempDown  = A.readConfig('Radiation.tempDown', double)
end

------------
--[[ IO ]]--
------------

IO.wrtRestart            = A.readConfig('wrtRestart', OnOrOff)
IO.restartEveryTimeSteps = A.readConfig('restartEveryTimeSteps', int)
IO.headerFrequency       = A.readConfig('headerFrequency', int)
IO.consoleFrequency      = A.readConfig('consoleFrequency', int)

-----------------------------------------------------------------------------
--[[                           FLUID GRID                                ]]--
-----------------------------------------------------------------------------

local fluidGrid = GRID.NewGrid{
  name = 'Fluid',
  xNum = Grid.xNum,
  yNum = Grid.yNum,
  zNum = Grid.zNum,
  xOrigin = Grid.xOrigin,
  yOrigin = Grid.yOrigin,
  zOrigin = Grid.zOrigin,
  xWidth = Grid.xWidth,
  yWidth = Grid.yWidth,
  zWidth = Grid.zWidth,
  xBnum = Grid.xBnum,
  yBnum = Grid.yBnum,
  zBnum = Grid.zBnum,
  xPeriodic = BC.xBCPeriodic,
  yPeriodic = BC.yBCPeriodic,
  zPeriodic = BC.zBCPeriodic,
}

-- Primitive variables
fluidGrid:NewField('rho', L.double)
fluidGrid:NewField('pressure', L.double)
fluidGrid:NewField('velocity', L.vec3d)

-- Remaining base variables
fluidGrid:NewField('centerCoordinates', L.vec3d)
fluidGrid:NewField('velocityGradientX', L.vec3d)
fluidGrid:NewField('velocityGradientY', L.vec3d)
fluidGrid:NewField('velocityGradientZ', L.vec3d)
fluidGrid:NewField('temperature', L.double)
fluidGrid:NewField('rhoEnthalpy', L.double)
fluidGrid:NewField('kineticEnergy', L.double)
fluidGrid:NewField('sgsEnergy', L.double)
fluidGrid:NewField('sgsEddyViscosity', L.double)
fluidGrid:NewField('sgsEddyKappa', L.double)
fluidGrid:NewField('convectiveSpectralRadius', L.double)
fluidGrid:NewField('viscousSpectralRadius', L.double)
fluidGrid:NewField('heatConductionSpectralRadius', L.double)

-- Conserved variables
fluidGrid:NewField('rhoVelocity', L.vec3d)
fluidGrid:NewField('rhoEnergy', L.double)

-- Fields for boundary treatment
fluidGrid:NewField('rhoBoundary', L.double)
fluidGrid:NewField('rhoVelocityBoundary', L.vec3d)
fluidGrid:NewField('rhoEnergyBoundary', L.double)
fluidGrid:NewField('velocityBoundary', L.vec3d)
fluidGrid:NewField('pressureBoundary', L.double)
fluidGrid:NewField('temperatureBoundary', L.double)
fluidGrid:NewField('velocityGradientXBoundary', L.vec3d)
fluidGrid:NewField('velocityGradientYBoundary', L.vec3d)
fluidGrid:NewField('velocityGradientZBoundary', L.vec3d)

-- scratch (temporary) fields
-- intermediate value and copies
fluidGrid:NewField('rho_old', L.double)
fluidGrid:NewField('rhoVelocity_old', L.vec3d)
fluidGrid:NewField('rhoEnergy_old', L.double)
fluidGrid:NewField('rho_new', L.double)
fluidGrid:NewField('rhoVelocity_new', L.vec3d)
fluidGrid:NewField('rhoEnergy_new', L.double)
-- time derivatives
fluidGrid:NewField('rho_t', L.double)
fluidGrid:NewField('rhoVelocity_t', L.vec3d)
fluidGrid:NewField('rhoEnergy_t', L.double)
-- fluxes
fluidGrid:NewField('rhoFluxX', L.double)
fluidGrid:NewField('rhoVelocityFluxX', L.vec3d)
fluidGrid:NewField('rhoEnergyFluxX', L.double)
fluidGrid:NewField('rhoFluxY', L.double)
fluidGrid:NewField('rhoVelocityFluxY', L.vec3d)
fluidGrid:NewField('rhoEnergyFluxY', L.double)
fluidGrid:NewField('rhoFluxZ', L.double)
fluidGrid:NewField('rhoVelocityFluxZ', L.vec3d)
fluidGrid:NewField('rhoEnergyFluxZ', L.double)

-- Right hand side of the kinetic energy equation
fluidGrid:NewField('PD', L.double)
fluidGrid:NewField('dissipation', L.double)
fluidGrid:NewField('dissipationFlux', L.double)

-----------------------------------------------------------------------------
--[[                         PARTICLES TABLE                             ]]--
-----------------------------------------------------------------------------

local particles = L.NewRelation {
    name = 'particles',
    mode = 'COUPLED',
    coupled_with = fluidGrid,
    coupling_field = 'cell',
    size = Particles.maxNum,
    max_skew = Particles.maxSkew,
    max_xfer_num = Particles.maxXferNum,
    xfer_stencil = {            { 0, 0, 1}, { 0, 0,-1},
                    { 0, 1, 0}, { 0, 1, 1}, { 0, 1,-1},
                    { 0,-1, 0}, { 0,-1, 1}, { 0,-1,-1},
                    { 1, 0, 0}, { 1, 0, 1}, { 1, 0,-1},
                    { 1, 1, 0}, { 1, 1, 1}, { 1, 1,-1},
                    { 1,-1, 0}, { 1,-1, 1}, { 1,-1,-1},
                    {-1, 0, 0}, {-1, 0, 1}, {-1, 0,-1},
                    {-1, 1, 0}, {-1, 1, 1}, {-1, 1,-1},
                    {-1,-1, 0}, {-1,-1, 1}, {-1,-1,-1}}
}

-- Primitive variables
particles:NewField('position', L.vec3d)
particles:NewField('velocity', L.vec3d)
particles:NewField('temperature', L.double)
particles:NewField('diameter', L.double)
particles:NewField('density', L.double)

-- Remaining base variables
particles:NewField('deltaVelocityOverRelaxationTime', L.vec3d)
particles:NewField('deltaTemperatureTerm', L.double)

-- scratch (temporary) fields
-- intermediate values and copies
particles:NewField('position_old', L.vec3d)
particles:NewField('velocity_old', L.vec3d)
particles:NewField('temperature_old', L.double)
particles:NewField('position_new', L.vec3d)
particles:NewField('velocity_new', L.vec3d)
particles:NewField('temperature_new', L.double)
particles:NewField('position_ghost', L.vec3d)
particles:NewField('velocity_ghost', L.vec3d)
particles:NewField('velocity_t_ghost', L.vec3d)

-- derivatives
particles:NewField('position_t', L.vec3d)
particles:NewField('velocity_t', L.vec3d)
particles:NewField('temperature_t', L.double)

-----------------------------------------------------------------------------
--[[                        RADIATION PREPROCESSING                      ]]--
-----------------------------------------------------------------------------

if Radiation.TYPE == 'DOM' then

  local domGrid = GRID.NewGrid{
    name = 'Radiation',
    xNum = Radiation.xNum,
    yNum = Radiation.yNum,
    zNum = Radiation.zNum,
    xOrigin = Grid.xOrigin,
    yOrigin = Grid.yOrigin,
    zOrigin = Grid.zOrigin,
    xWidth = Grid.xWidth,
    yWidth = Grid.yWidth,
    zWidth = Grid.zWidth,
    xBnum = Radiation.xBnum,
    yBnum = Radiation.yBnum,
    zBnum = Radiation.zBnum,
    xPeriodic = Radiation.xPeriodic,
    yPeriodic = Radiation.yPeriodic,
    zPeriodic = Radiation.zPeriodic,
  }
  fluidGrid:LinkWithCoarse(domGrid, 'to_Radiation')

  -- cell center intensity per angle
  domGrid:NewField('I_1', L.vector(L.double, Radiation.NUM_ANGLES))
  domGrid:NewField('I_2', L.vector(L.double, Radiation.NUM_ANGLES))
  domGrid:NewField('I_3', L.vector(L.double, Radiation.NUM_ANGLES))
  domGrid:NewField('I_4', L.vector(L.double, Radiation.NUM_ANGLES))
  domGrid:NewField('I_5', L.vector(L.double, Radiation.NUM_ANGLES))
  domGrid:NewField('I_6', L.vector(L.double, Radiation.NUM_ANGLES))
  domGrid:NewField('I_7', L.vector(L.double, Radiation.NUM_ANGLES))
  domGrid:NewField('I_8', L.vector(L.double, Radiation.NUM_ANGLES))

  -- iterative intensity per angle
  domGrid:NewField('Iiter_1', L.vector(L.double, Radiation.NUM_ANGLES))
  domGrid:NewField('Iiter_2', L.vector(L.double, Radiation.NUM_ANGLES))
  domGrid:NewField('Iiter_3', L.vector(L.double, Radiation.NUM_ANGLES))
  domGrid:NewField('Iiter_4', L.vector(L.double, Radiation.NUM_ANGLES))
  domGrid:NewField('Iiter_5', L.vector(L.double, Radiation.NUM_ANGLES))
  domGrid:NewField('Iiter_6', L.vector(L.double, Radiation.NUM_ANGLES))
  domGrid:NewField('Iiter_7', L.vector(L.double, Radiation.NUM_ANGLES))
  domGrid:NewField('Iiter_8', L.vector(L.double, Radiation.NUM_ANGLES))

  domGrid:NewField('G',     L.double) -- intensity summation over all angles
  domGrid:NewField('S',     L.double) -- source term
  domGrid:NewField('Ib',    L.double) -- blackbody intensity
  domGrid:NewField('sigma', L.double) -- extinction coefficient

  -- partial sums, over particles inside volume
  domGrid:NewField('acc_d2',   L.double) -- p.diameter^2
  domGrid:NewField('acc_d2t4', L.double) -- p.diameter^2 * p.temperature^4

end

-----------------------------------------------------------------------------
--[[                       BOUNDARY CONFIG CHECK                         ]]--
-----------------------------------------------------------------------------

-- Define offsets for boundary conditions in flow solver
-- The sign variables define the necessary reflections for the
-- different types of BCs. The wall velocity is specified above,
-- and then the velocity adjustment is calculated here and applied
-- to the boundaries below.

-- Define offsets, signs, and velocities for the x BCs

M.IF(M.AND(M.EQ(BC.xBCLeft,  FlowBC.Periodic),
           M.EQ(BC.xBCRight, FlowBC.Periodic)))
  BC.xSign:set({1.0,1.0,1.0})
  BC.xPosVelocity:set({0.0,0.0,0.0})
  BC.xNegVelocity:set({0.0,0.0,0.0})
  BC.xPosTemperature:set(-1.0)
  BC.xNegTemperature:set(-1.0)
  BC.xBCLeftParticles:set(ParticleBC.Permeable)
  BC.xBCRightParticles:set(ParticleBC.Permeable)
M.ELSE() M.IF(M.AND(M.EQ(BC.xBCLeft,  FlowBC.Symmetry),
                    M.EQ(BC.xBCRight, FlowBC.Symmetry)))
  BC.xSign:set({-1.0,1.0,1.0})
  BC.xPosVelocity:set({0.0,0.0,0.0})
  BC.xNegVelocity:set({0.0,0.0,0.0})
  BC.xPosTemperature:set(-1.0)
  BC.xNegTemperature:set(-1.0)
  BC.xBCLeftParticles:set(ParticleBC.Solid)
  BC.xBCRightParticles:set(ParticleBC.Solid)
M.ELSE() M.IF(M.AND(M.EQ(BC.xBCLeft,  FlowBC.AdiabaticWall),
                    M.EQ(BC.xBCRight, FlowBC.AdiabaticWall)))
  BC.xSign:set({-1.0,-1.0,-1.0})
  BC.xPosVelocity:set(M.ARRAY{2.0*BC.xBCRightVel[0],
                              2.0*BC.xBCRightVel[1],
                              2.0*BC.xBCRightVel[2]})
  BC.xNegVelocity:set(M.ARRAY{2.0*BC.xBCLeftVel[0],
                              2.0*BC.xBCLeftVel[1],
                              2.0*BC.xBCLeftVel[2]})
  BC.xPosTemperature:set(-1.0)
  BC.xNegTemperature:set(-1.0)
  BC.xBCLeftParticles:set(ParticleBC.Solid)
  BC.xBCRightParticles:set(ParticleBC.Solid)
M.ELSE() M.IF(M.AND(M.EQ(BC.xBCLeft,  FlowBC.IsothermalWall),
                    M.EQ(BC.xBCRight, FlowBC.IsothermalWall)))
  BC.xSign:set({-1.0,-1.0,-1.0})
  BC.xPosVelocity:set(M.ARRAY{2.0*BC.xBCRightVel[0],
                              2.0*BC.xBCRightVel[1],
                              2.0*BC.xBCRightVel[2]})
  BC.xNegVelocity:set(M.ARRAY{2.0*BC.xBCLeftVel[0],
                              2.0*BC.xBCLeftVel[1],
                              2.0*BC.xBCLeftVel[2]})
  BC.xPosTemperature:set(BC.xBCRightTemp)
  BC.xNegTemperature:set(BC.xBCLeftTemp)
  BC.xBCLeftParticles:set(ParticleBC.Solid)
  BC.xBCRightParticles:set(ParticleBC.Solid)
M.ELSE()
  M.ERROR('Boundary conditions in x not implemented')
M.END() M.END() M.END() M.END()

-- Define offsets, signs, and velocities for the y BCs

M.IF(M.AND(M.EQ(BC.yBCLeft,  FlowBC.Periodic),
           M.EQ(BC.yBCRight, FlowBC.Periodic)))
  BC.ySign:set({1.0,1.0,1.0})
  BC.yPosVelocity:set({0.0,0.0,0.0})
  BC.yNegVelocity:set({0.0,0.0,0.0})
  BC.yPosTemperature:set(-1.0)
  BC.yNegTemperature:set(-1.0)
  BC.yBCLeftParticles:set(ParticleBC.Permeable)
  BC.yBCRightParticles:set(ParticleBC.Permeable)
M.ELSE() M.IF(M.AND(M.EQ(BC.yBCLeft,  FlowBC.Symmetry),
                    M.EQ(BC.yBCRight, FlowBC.Symmetry)))
  BC.ySign:set({1.0,-1.0,1.0})
  BC.yPosVelocity:set({0.0,0.0,0.0})
  BC.yNegVelocity:set({0.0,0.0,0.0})
  BC.yPosTemperature:set(-1.0)
  BC.yNegTemperature:set(-1.0)
  BC.yBCLeftParticles:set(ParticleBC.Solid)
  BC.yBCRightParticles:set(ParticleBC.Solid)
M.ELSE() M.IF(M.AND(M.EQ(BC.yBCLeft,  FlowBC.AdiabaticWall),
                    M.EQ(BC.yBCRight, FlowBC.AdiabaticWall)))
  BC.ySign:set({-1.0,-1.0,-1.0})
  BC.yPosVelocity:set(M.ARRAY{2.0*BC.yBCRightVel[0],
                              2.0*BC.yBCRightVel[1],
                              2.0*BC.yBCRightVel[2]})
  BC.yNegVelocity:set(M.ARRAY{2.0*BC.yBCLeftVel[0],
                              2.0*BC.yBCLeftVel[1],
                              2.0*BC.yBCLeftVel[2]})
  BC.yPosTemperature:set(-1.0)
  BC.yNegTemperature:set(-1.0)
  BC.yBCLeftParticles:set(ParticleBC.Solid)
  BC.yBCRightParticles:set(ParticleBC.Solid)
M.ELSE() M.IF(M.AND(M.EQ(BC.yBCLeft,  FlowBC.IsothermalWall),
                    M.EQ(BC.yBCRight, FlowBC.IsothermalWall)))
  BC.ySign:set({-1.0,-1.0,-1.0})
  BC.yPosVelocity:set(M.ARRAY{2.0*BC.yBCRightVel[0],
                              2.0*BC.yBCRightVel[1],
                              2.0*BC.yBCRightVel[2]})
  BC.yNegVelocity:set(M.ARRAY{2.0*BC.yBCLeftVel[0],
                              2.0*BC.yBCLeftVel[1],
                              2.0*BC.yBCLeftVel[2]})
  BC.yPosTemperature:set(BC.yBCRightTemp)
  BC.yNegTemperature:set(BC.yBCLeftTemp)
  BC.yBCLeftParticles:set(ParticleBC.Solid)
  BC.yBCRightParticles:set(ParticleBC.Solid)
M.ELSE()
  M.ERROR('Boundary conditions in y not implemented')
M.END() M.END() M.END() M.END()

-- Define offsets, signs, and velocities for the z BCs

M.IF(M.AND(M.EQ(BC.zBCLeft,  FlowBC.Periodic),
           M.EQ(BC.zBCRight, FlowBC.Periodic)))
  BC.zSign:set({1.0,1.0,1.0})
  BC.zPosVelocity:set({0.0,0.0,0.0})
  BC.zNegVelocity:set({0.0,0.0,0.0})
  BC.zPosTemperature:set(-1.0)
  BC.zNegTemperature:set(-1.0)
  BC.zBCLeftParticles:set(ParticleBC.Permeable)
  BC.zBCRightParticles:set(ParticleBC.Permeable)
M.ELSE() M.IF(M.AND(M.EQ(BC.zBCLeft,  FlowBC.Symmetry),
                    M.EQ(BC.zBCRight, FlowBC.Symmetry)))
  BC.zSign:set({1.0,1.0,-1.0})
  BC.zPosVelocity:set({0.0,0.0,0.0})
  BC.zNegVelocity:set({0.0,0.0,0.0})
  BC.zPosTemperature:set(-1.0)
  BC.zNegTemperature:set(-1.0)
  BC.zBCLeftParticles:set(ParticleBC.Solid)
  BC.zBCRightParticles:set(ParticleBC.Solid)
M.ELSE() M.IF(M.AND(M.EQ(BC.zBCLeft,  FlowBC.AdiabaticWall),
                    M.EQ(BC.zBCRight, FlowBC.AdiabaticWall)))
  BC.zSign:set({-1.0,-1.0,-1.0})
  BC.zPosVelocity:set(M.ARRAY{2.0*BC.zBCRightVel[0],
                              2.0*BC.zBCRightVel[1],
                              2.0*BC.zBCRightVel[2]})
  BC.zNegVelocity:set(M.ARRAY{2.0*BC.zBCLeftVel[0],
                              2.0*BC.zBCLeftVel[1],
                              2.0*BC.zBCLeftVel[2]})
  BC.zPosTemperature:set(-1.0)
  BC.zNegTemperature:set(-1.0)
  BC.zBCLeftParticles:set(ParticleBC.Solid)
  BC.zBCRightParticles:set(ParticleBC.Solid)
M.ELSE() M.IF(M.AND(M.EQ(BC.zBCLeft,  FlowBC.IsothermalWall),
                    M.EQ(BC.zBCRight, FlowBC.IsothermalWall)))
  BC.zSign:set({-1.0,-1.0,-1.0})
  BC.zPosVelocity:set(M.ARRAY{2.0*BC.zBCRightVel[0],
                              2.0*BC.zBCRightVel[1],
                              2.0*BC.zBCRightVel[2]})
  BC.zNegVelocity:set(M.ARRAY{2.0*BC.zBCLeftVel[0],
                              2.0*BC.zBCLeftVel[1],
                              2.0*BC.zBCLeftVel[2]})
  BC.zPosTemperature:set(BC.zBCRightTemp)
  BC.zNegTemperature:set(BC.zBCLeftTemp)
  BC.zBCLeftParticles:set(ParticleBC.Solid)
  BC.zBCRightParticles:set(ParticleBC.Solid)
M.ELSE()
  M.ERROR('Boundary conditions in z not implemented')
M.END() M.END() M.END() M.END()

-- Check boundary type consistency for the periodic BCs

M.IF(M.NOT(M.OR(M.AND(M.EQ(BC.xBCLeft,  FlowBC.Periodic),
                      M.EQ(BC.xBCRight, FlowBC.Periodic)),
                M.AND(M.NOT(M.EQ(BC.xBCLeft,  FlowBC.Periodic)),
                      M.NOT(M.EQ(BC.xBCRight, FlowBC.Periodic))))))
  M.ERROR("Boundary conditions in x should match for periodicity")
M.END()

M.IF(M.NOT(M.OR(M.AND(M.EQ(BC.yBCLeft,  FlowBC.Periodic),
                      M.EQ(BC.yBCRight, FlowBC.Periodic)),
                M.AND(M.NOT(M.EQ(BC.yBCLeft,  FlowBC.Periodic)),
                      M.NOT(M.EQ(BC.yBCRight, FlowBC.Periodic))))))
  M.ERROR("Boundary conditions in y should match for periodicity")
M.END()

M.IF(M.NOT(M.OR(M.AND(M.EQ(BC.zBCLeft,  FlowBC.Periodic),
                      M.EQ(BC.zBCRight, FlowBC.Periodic)),
                M.AND(M.NOT(M.EQ(BC.zBCLeft,  FlowBC.Periodic)),
                      M.NOT(M.EQ(BC.zBCRight, FlowBC.Periodic))))))
  M.ERROR("Boundary conditions in z should match for periodicity")
M.END()

-----------------------------------------------------------------------------
--[[                       EXTERNAL REGENT MODULES                       ]]--
-----------------------------------------------------------------------------

local PARTICLES_INIT =
  (require 'particlesInit')(particles, fluidGrid,
                            Grid.xBnum, Grid.yBnum, Grid.zBnum)

if Radiation.TYPE == 'Algebraic' then
  ALGEBRAIC = (require 'algebraic')(particles)
end

if Radiation.TYPE == 'DOM' then
  DOM = (require 'dom/dom')(
    domGrid, Radiation.NUM_ANGLES,
    Radiation.xNum, Radiation.yNum, Radiation.zNum,
    Radiation.xCellWidth, Radiation.yCellWidth, Radiation.zCellWidth)
end

-----------------------------------------------------------------------------
--[[                       LOAD DATA FOR RESTART                         ]]--
-----------------------------------------------------------------------------

M.IF(M.EQ(Flow.initCase, FlowInitCase.Restart))
  -- Increment the time step and physical time so the simulation doesn't
  -- repeat from 0. Also, increase the max number of iterations so the solve
  -- doesn't immediately exit.
  Integrator.timeStep:set(Integrator.restartIter)
  -- TODO: No way to pass Integrator.simTime for the restart
M.END()

-----------------------------------------------------------------------------
--[[                       USER-DEFINED FUNCTIONS                        ]]--
-----------------------------------------------------------------------------

-- Norm of a vector
local ebb Norm (v)
  return L.sqrt(L.dot(v, v))
end

-- Compute fluid dynamic viscosity from fluid temperature
local ebb GetDynamicViscosity (temperature)
  var viscosity = L.double(0.0)
  if Flow.viscosityModel == ViscosityModel.Constant then
    viscosity = Flow.constantVisc
  elseif Flow.viscosityModel == ViscosityModel.PowerLaw then
    viscosity = Flow.powerlawViscRef *
      L.pow(temperature/Flow.powerlawTempRef, 0.75)
  elseif Flow.viscosityModel == ViscosityModel.Sutherland then
    viscosity = Flow.sutherlandViscRef *
      L.pow((temperature/Flow.sutherlandTempRef),(3.0/2.0))*
      ((Flow.sutherlandTempRef + Flow.sutherlandSRef)/
         (temperature + Flow.sutherlandSRef))
  else L.assert(false) end
  return viscosity
end

-- Compute fluid flow sound speed based on temperature (a = sqrt(gamma*R*T))
local ebb GetSoundSpeed (temperature)
  return L.sqrt(Flow.gamma * Flow.gasConstant * temperature)
end

-- Functions to retrieve particle area, volume and mass
particles:NewFieldReadFunction('cross_section_area', ebb(p)
  return pi * L.pow(p.diameter, 2) / 4.0
end)
particles:NewFieldReadFunction('volume', ebb(p)
  return pi * L.pow(p.diameter, 3) / 6.0
end)
particles:NewFieldReadFunction('mass', ebb(p)
  return p.volume * p.density
end)

-----------------------------------------------------------------------------
--[[                              EBB MACROS                             ]]--
-----------------------------------------------------------------------------

local ebb TrilinearInterpolateRho (xyz, c000, c100, c010, c110, c001, c101, c011, c111)
  var dX   = L.fmod((xyz[0] - Grid.xRealOrigin)/Grid.xCellWidth + 0.5, 1.0)
  var dY   = L.fmod((xyz[1] - Grid.yRealOrigin)/Grid.yCellWidth + 0.5, 1.0)
  var dZ   = L.fmod((xyz[2] - Grid.zRealOrigin)/Grid.zCellWidth + 0.5, 1.0)

  var oneMinusdX = 1.0 - dX
  var oneMinusdY = 1.0 - dY
  var oneMinusdZ = 1.0 - dZ
  var weight00 = c000 * oneMinusdX + c100 * dX
  var weight10 = c010 * oneMinusdX + c110 * dX
  var weight01 = c001 * oneMinusdX + c101 * dX
  var weight11 = c011 * oneMinusdX + c111 * dX
  var weight0  = weight00 * oneMinusdY + weight10 * dY
  var weight1  = weight01 * oneMinusdY + weight11 * dY

  return weight0 * oneMinusdZ + weight1 * dZ
end

local ebb TrilinearInterpolateVelocity (xyz, c000, c100, c010, c110, c001, c101, c011, c111)
  var dX   = L.fmod((xyz[0] - Grid.xRealOrigin)/Grid.xCellWidth + 0.5, 1.0)
  var dY   = L.fmod((xyz[1] - Grid.yRealOrigin)/Grid.yCellWidth + 0.5, 1.0)
  var dZ   = L.fmod((xyz[2] - Grid.zRealOrigin)/Grid.zCellWidth + 0.5, 1.0)

  var oneMinusdX = 1.0 - dX
  var oneMinusdY = 1.0 - dY
  var oneMinusdZ = 1.0 - dZ
  var weight00 = c000 * oneMinusdX + c100 * dX
  var weight10 = c010 * oneMinusdX + c110 * dX
  var weight01 = c001 * oneMinusdX + c101 * dX
  var weight11 = c011 * oneMinusdX + c111 * dX
  var weight0  = weight00 * oneMinusdY + weight10 * dY
  var weight1  = weight01 * oneMinusdY + weight11 * dY

  return weight0 * oneMinusdZ + weight1 * dZ
end

local ebb TrilinearInterpolateTemp (xyz, c000, c100, c010, c110, c001, c101, c011, c111)
  var dX   = L.fmod((xyz[0] - Grid.xRealOrigin)/Grid.xCellWidth + 0.5, 1.0)
  var dY   = L.fmod((xyz[1] - Grid.yRealOrigin)/Grid.yCellWidth + 0.5, 1.0)
  var dZ   = L.fmod((xyz[2] - Grid.zRealOrigin)/Grid.zCellWidth + 0.5, 1.0)

  var oneMinusdX = 1.0 - dX
  var oneMinusdY = 1.0 - dY
  var oneMinusdZ = 1.0 - dZ
  var weight00 = c000 * oneMinusdX + c100 * dX
  var weight10 = c010 * oneMinusdX + c110 * dX
  var weight01 = c001 * oneMinusdX + c101 * dX
  var weight11 = c011 * oneMinusdX + c111 * dX
  var weight0  = weight00 * oneMinusdY + weight10 * dY
  var weight1  = weight01 * oneMinusdY + weight11 * dY

  return weight0 * oneMinusdZ + weight1 * dZ
end

local ebb InterpolateTriVelocity (c, xyz)
  var velocity000 = L.vec3d({0.0,0.0,0.0})
  var velocity100 = L.vec3d({0.0,0.0,0.0})
  var velocity010 = L.vec3d({0.0,0.0,0.0})
  var velocity110 = L.vec3d({0.0,0.0,0.0})
  var velocity001 = L.vec3d({0.0,0.0,0.0})
  var velocity101 = L.vec3d({0.0,0.0,0.0})
  var velocity011 = L.vec3d({0.0,0.0,0.0})
  var velocity111 = L.vec3d({0.0,0.0,0.0})

  var velocity0 = c.velocity

  if xyz[0] > c.centerCoordinates[0] then
    var velocityb = c( 1,  0,  0).velocity
    if xyz[1] > c.centerCoordinates[1] then
      var velocityaa = c( 0,  1,  0).velocity
      var velocityab = c( 1,  1,  0).velocity
      if xyz[2] > c.centerCoordinates[2] then
        velocity000 = velocity0
        velocity100 = velocityb
        velocity010 = velocityaa
        velocity110 = velocityab
        velocity001 = c( 0,  0,  1).velocity
        velocity101 = c( 1,  0,  1).velocity
        velocity011 = c( 0,  1,  1).velocity
        velocity111 = c( 1,  1,  1).velocity

      else -- xyz[2] <= c.centerCoordinates[2]
        velocity000 = c( 0,  0, -1).velocity
        velocity100 = c( 1,  0, -1).velocity
        velocity010 = c( 0,  1, -1).velocity
        velocity110 = c( 1,  1, -1).velocity
        velocity001 = velocity0
        velocity101 = velocityb
        velocity011 = velocityaa
        velocity111 = velocityab

      end

    else -- xyz[1] <= c.centerCoordinates[1]
      var velocityaa = c( 0, -1,  0).velocity
      var velocityab = c( 1, -1,  0).velocity

      if xyz[2] > c.centerCoordinates[2] then
        velocity000 = velocityaa
        velocity100 = velocityab
        velocity010 = velocity0
        velocity110 = velocityb
        velocity001 = c( 0, -1,  1).velocity
        velocity101 = c( 1, -1,  1).velocity
        velocity011 = c( 0,  0,  1).velocity
        velocity111 = c( 1,  0,  1).velocity

      else -- xyz[2] <= c.centerCoordinates[2]
        velocity000 = c( 0, -1, -1).velocity
        velocity100 = c( 1, -1, -1).velocity
        velocity010 = c( 0,  0, -1).velocity
        velocity110 = c( 1,  0, -1).velocity
        velocity001 = velocityaa
        velocity101 = velocityab
        velocity011 = velocity0
        velocity111 = velocityb

      end
    end

  else -- xyz[0] <= c.centerCoordinates[0]
    var velocitya = c(-1,  0,  0).velocity
    if xyz[1] > c.centerCoordinates[1] then
      var velocityaa = c(-1,  1,  0).velocity
      var velocityab = c( 0,  1,  0).velocity

      if xyz[2] > c.centerCoordinates[2] then
        velocity000 = velocitya
        velocity100 = velocity0
        velocity010 = velocityaa
        velocity110 = velocityab
        velocity001 = c(-1,  0,  1).velocity
        velocity101 = c( 0,  0,  1).velocity
        velocity011 = c(-1,  1,  1).velocity
        velocity111 = c( 0,  1,  1).velocity

      else -- xyz[2] <= c.centerCoordinates[2]
        velocity000 = c(-1,  0, -1).velocity
        velocity100 = c( 0,  0, -1).velocity
        velocity010 = c(-1,  1, -1).velocity
        velocity110 = c( 0,  1, -1).velocity
        velocity001 = velocitya
        velocity101 = velocity0
        velocity011 = velocityaa
        velocity111 = velocityab

      end

    else -- xyz[1] <= c.centerCoordinates[1]
      var velocityaa = c(-1, -1,  0).velocity
      var velocityab = c( 0, -1,  0).velocity

      if xyz[2] > c.centerCoordinates[2] then
        velocity000 = velocityaa
        velocity100 = velocityab
        velocity010 = velocitya
        velocity110 = velocity0
        velocity001 = c(-1, -1,  1).velocity
        velocity101 = c( 0, -1,  1).velocity
        velocity011 = c(-1,  0,  1).velocity
        velocity111 = c( 0,  0,  1).velocity

      else -- xyz[2] <= c.centerCoordinates[2]
        velocity000 = c(-1, -1, -1).velocity
        velocity100 = c( 0, -1, -1).velocity
        velocity010 = c(-1,  0, -1).velocity
        velocity110 = c( 0,  0, -1).velocity
        velocity001 = velocityaa
        velocity101 = velocityab
        velocity011 = velocitya
        velocity111 = velocity0

      end
    end
  end

  return TrilinearInterpolateVelocity (xyz, velocity000, velocity100,
                                            velocity010, velocity110,
                                            velocity001, velocity101,
                                            velocity011, velocity111)
end

local ebb InterpolateTriTemp (c, xyz)
  var temp000 = L.double(0.0)
  var temp100 = L.double(0.0)
  var temp010 = L.double(0.0)
  var temp110 = L.double(0.0)
  var temp001 = L.double(0.0)
  var temp101 = L.double(0.0)
  var temp011 = L.double(0.0)
  var temp111 = L.double(0.0)

  var temp0     = c.temperature

  if xyz[0] > c.centerCoordinates[0] then
    var tempb     = c( 1,  0,  0).temperature
    if xyz[1] > c.centerCoordinates[1] then
      var tempaa = c( 0,  1,  0).temperature
      var tempab = c( 1,  1,  0).temperature
      if xyz[2] > c.centerCoordinates[2] then
        temp000 = temp0
        temp100 = tempb
        temp010 = tempaa
        temp110 = tempab
        temp001 = c( 0,  0,  1).temperature
        temp101 = c( 1,  0,  1).temperature
        temp011 = c( 0,  1,  1).temperature
        temp111 = c( 1,  1,  1).temperature

      else -- xyz[2] <= c.centerCoordinates[2]
        temp000 = c( 0,  0, -1).temperature
        temp100 = c( 1,  0, -1).temperature
        temp010 = c( 0,  1, -1).temperature
        temp110 = c( 1,  1, -1).temperature
        temp001 = temp0
        temp101 = tempb
        temp011 = tempaa
        temp111 = tempab

      end

    else -- xyz[1] <= c.centerCoordinates[1]
      var tempaa = c( 0, -1,  0).temperature
      var tempab = c( 1, -1,  0).temperature

      if xyz[2] > c.centerCoordinates[2] then
        temp000 = tempaa
        temp100 = tempab
        temp010 = temp0
        temp110 = tempb
        temp001 = c( 0, -1,  1).temperature
        temp101 = c( 1, -1,  1).temperature
        temp011 = c( 0,  0,  1).temperature
        temp111 = c( 1,  0,  1).temperature

      else -- xyz[2] <= c.centerCoordinates[2]
        temp000 = c( 0, -1, -1).temperature
        temp100 = c( 1, -1, -1).temperature
        temp010 = c( 0,  0, -1).temperature
        temp110 = c( 1,  0, -1).temperature
        temp001 = tempaa
        temp101 = tempab
        temp011 = temp0
        temp111 = tempb

      end
    end

  else -- xyz[0] <= c.centerCoordinates[0]
    var tempa = c(-1,  0,  0).temperature
    if xyz[1] > c.centerCoordinates[1] then
      var tempaa = c(-1,  1,  0).temperature
      var tempab = c( 0,  1,  0).temperature

      if xyz[2] > c.centerCoordinates[2] then
        temp000 = tempa
        temp100 = temp0
        temp010 = tempaa
        temp110 = tempab
        temp001 = c(-1,  0,  1).temperature
        temp101 = c( 0,  0,  1).temperature
        temp011 = c(-1,  1,  1).temperature
        temp111 = c( 0,  1,  1).temperature

      else -- xyz[2] <= c.centerCoordinates[2]
        temp000 = c(-1,  0, -1).temperature
        temp100 = c( 0,  0, -1).temperature
        temp010 = c(-1,  1, -1).temperature
        temp110 = c( 0,  1, -1).temperature
        temp001 = tempa
        temp101 = temp0
        temp011 = tempaa
        temp111 = tempab

      end

    else -- xyz[1] <= c.centerCoordinates[1]
      var tempaa = c(-1, -1,  0).temperature
      var tempab = c( 0, -1,  0).temperature

      if xyz[2] > c.centerCoordinates[2] then
        temp000 = tempaa
        temp100 = tempab
        temp010 = tempa
        temp110 = temp0
        temp001 = c(-1, -1,  1).temperature
        temp101 = c( 0, -1,  1).temperature
        temp011 = c(-1,  0,  1).temperature
        temp111 = c( 0,  0,  1).temperature

      else -- xyz[2] <= c.centerCoordinates[2]
        temp000 = c(-1, -1, -1).temperature
        temp100 = c( 0, -1, -1).temperature
        temp010 = c(-1,  0, -1).temperature
        temp110 = c( 0,  0, -1).temperature
        temp001 = tempaa
        temp101 = tempab
        temp011 = tempa
        temp111 = temp0

      end
    end
  end

  return TrilinearInterpolateTemp (xyz, temp000, temp100,
                                        temp010, temp110,
                                        temp001, temp101,
                                        temp011, temp111)
end

-----------------------------------------------------------------------------
--[[                            EBB FUNCTIONS                            ]]--
-----------------------------------------------------------------------------

-------
-- FLOW
-------

ebb Flow.InitializeCell (c : fluidGrid)
  c.rho = 0.0
  c.pressure = 0.0
  c.velocity = L.vec3d({0.0, 0.0, 0.0})
  c.centerCoordinates = L.vec3d({0.0, 0.0, 0.0})
  c.velocityGradientX = L.vec3d({0.0, 0.0, 0.0})
  c.velocityGradientY = L.vec3d({0.0, 0.0, 0.0})
  c.velocityGradientZ = L.vec3d({0.0, 0.0, 0.0})
  c.temperature = 0.0
  c.rhoEnthalpy = 0.0
  c.kineticEnergy = 0.0
  c.sgsEnergy = 0.0
  c.sgsEddyViscosity = 0.0
  c.sgsEddyKappa = 0.0
  c.convectiveSpectralRadius = 0.0
  c.viscousSpectralRadius = 0.0
  c.heatConductionSpectralRadius = 0.0
  c.rhoVelocity = L.vec3d({0.0, 0.0, 0.0})
  c.rhoEnergy = 0.0
  c.rhoBoundary = 0.0
  c.rhoVelocityBoundary = L.vec3d({0.0, 0.0, 0.0})
  c.rhoEnergyBoundary = 0.0
  c.velocityBoundary = L.vec3d({0.0, 0.0, 0.0})
  c.pressureBoundary = 0.0
  c.temperatureBoundary = 0.0
  c.velocityGradientXBoundary = L.vec3d({0.0, 0.0, 0.0})
  c.velocityGradientYBoundary = L.vec3d({0.0, 0.0, 0.0})
  c.velocityGradientZBoundary = L.vec3d({0.0, 0.0, 0.0})
  c.rho_old = 0.0
  c.rhoVelocity_old = L.vec3d({0.0, 0.0, 0.0})
  c.rhoEnergy_old = 0.0
  c.rho_new = 0.0
  c.rhoVelocity_new = L.vec3d({0.0, 0.0, 0.0})
  c.rhoEnergy_new = 0.0
  c.rho_t = 0.0
  c.rhoVelocity_t = L.vec3d({0.0, 0.0, 0.0})
  c.rhoEnergy_t = 0.0
  c.rhoFluxX = 0.0
  c.rhoVelocityFluxX = L.vec3d({0.0, 0.0, 0.0})
  c.rhoEnergyFluxX = 0.0
  c.rhoFluxY = 0.0
  c.rhoVelocityFluxY = L.vec3d({0.0, 0.0, 0.0})
  c.rhoEnergyFluxY = 0.0
  c.rhoFluxZ = 0.0
  c.rhoVelocityFluxZ = L.vec3d({0.0, 0.0, 0.0})
  c.rhoEnergyFluxZ = 0.0
  c.PD = 0.0
  c.dissipation = 0.0
  c.dissipationFlux = 0.0
end

-- Initialize flow variables
-- Cell center coordinates are stored in the grid field macro 'center'.
-- Here, we use a field for convenience when outputting to file, but this is
-- to be removed after grid outputing is well defined from within the grid.t
-- module. Similar story with the vertex coordinates (output only).
ebb Flow.InitializeCenterCoordinates (c : fluidGrid)
  var xy = c.center
  c.centerCoordinates = L.vec3d({L.double(xy[0]), L.double(xy[1]), L.double((xy[2]))})
end

ebb Flow.InitializeUniform (c : fluidGrid)
  c.rho         = Flow.initParams[0]
  c.pressure    = Flow.initParams[1]
  c.velocity[0] = Flow.initParams[2]
  c.velocity[1] = Flow.initParams[3]
  c.velocity[2] = Flow.initParams[4]
end

ebb Flow.InitializeTaylorGreen2D (c : fluidGrid)
  -- Define Taylor Green Vortex
  var taylorGreenDensity  = Flow.initParams[0]
  var taylorGreenPressure = Flow.initParams[1]
  var taylorGreenVelocity = Flow.initParams[2]
  -- Initialize
  var xy = c.center
  var coorZ = 0
  c.rho = taylorGreenDensity
  c.velocity =
    taylorGreenVelocity *
    L.vec3d({L.sin(xy[0]) * L.cos(xy[1]) * L.cos(coorZ),
             - L.cos(xy[0]) * L.sin(xy[1]) * L.cos(coorZ),
             0})
  var factorA = L.cos(2.0*coorZ) + 2.0
  var factorB = L.cos(2.0*xy[0]) + L.cos(2.0*xy[1])
  c.pressure =
    taylorGreenPressure +
    taylorGreenDensity * L.pow(taylorGreenVelocity,2) / 16 *
    factorA * factorB
end

ebb Flow.InitializeTaylorGreen3D (c : fluidGrid)
  -- Define Taylor Green Vortex
  var taylorGreenDensity  = Flow.initParams[0]
  var taylorGreenPressure = Flow.initParams[1]
  var taylorGreenVelocity = Flow.initParams[2]
  -- Initialize
  var xy = c.center
  c.rho = taylorGreenDensity
  c.velocity =
    taylorGreenVelocity *
    L.vec3d({L.sin(xy[0]) * L.cos(xy[1]) * L.cos(xy[2]),
             - L.cos(xy[0]) * L.sin(xy[1]) * L.cos(xy[2]),
             0})
  var factorA = L.cos(2.0*xy[2]) + 2.0
  var factorB = L.cos(2.0*xy[0]) + L.cos(2.0*xy[1])
  c.pressure =
    taylorGreenPressure +
    taylorGreenDensity * L.pow(taylorGreenVelocity,2) / 16 *
    factorA * factorB
end

ebb Flow.InitializePerturbed (c : fluidGrid)
  -- This initialization imposes a small random perturbation in
  -- the velocity field used to start up forced turbulence cases
  c.rho         = Flow.initParams[0]
  c.pressure    = Flow.initParams[1]
  c.velocity[0] = Flow.initParams[2] + ((L.rand()-0.5)*10.0)
  c.velocity[1] = Flow.initParams[3] + ((L.rand()-0.5)*10.0)
  c.velocity[2] = Flow.initParams[4] + ((L.rand()-0.5)*10.0)
end

ebb Flow.UpdateConservedFromPrimitive (c : fluidGrid)
  if c.in_interior then
    -- Equation of state: T = p / ( R * rho )
    var tmpTemperature = c.pressure / (Flow.gasConstant * c.rho)
    var velocity = c.velocity
    c.rhoVelocity = c.rho * c.velocity

    -- rhoE = rhoe (= rho * cv * T) + kineticEnergy + sgsEnergy
    var cv = Flow.gasConstant / (Flow.gamma - 1.0)
    c.rhoEnergy =
      c.rho * ( cv * tmpTemperature + 0.5 * L.dot(velocity,velocity) )
      + c.sgsEnergy
  end
end

-- Initialize temporaries
ebb Flow.InitializeTemporaries (c : fluidGrid)
  c.rho_old         = c.rho
  c.rhoVelocity_old = c.rhoVelocity
  c.rhoEnergy_old   = c.rhoEnergy
  c.rho_new         = c.rho
  c.rhoVelocity_new = c.rhoVelocity
  c.rhoEnergy_new   = c.rhoEnergy
end

-- Initialize derivatives
ebb Flow.InitializeTimeDerivatives (c : fluidGrid)
  c.rho_t         = L.double(0.0)
  c.rhoVelocity_t = L.vec3d({0.0, 0.0, 0.0})
  c.rhoEnergy_t   = L.double(0.0)
  -- Initialize enthalpy
  c.rhoEnthalpy = c.rhoEnergy + c.pressure
end

---------------------
-- Particles coupling
---------------------

ebb Flow.AddParticlesCoupling (p : particles)

  -- WARNING: Assumes that deltaVelocityOverRelaxationTime and
  -- deltaTemperatureTerm have been computed previously, and that
  -- we have called the cell_locate kernel for the particles.
  -- (for example, when adding the flow coupling to the particles,
  -- which should be called before in the time stepper)

  -- Add contribution to momentum and energy equations from the previously
  -- computed deltaVelocityOverRelaxationTime and deltaTemperatureTerm
  p.cell.rhoVelocity_t += -p.mass * p.deltaVelocityOverRelaxationTime / Grid.cellVolume
  p.cell.rhoEnergy_t   += -p.deltaTemperatureTerm / Grid.cellVolume
end

--------------
-- Body Forces
--------------

ebb Flow.AddBodyForces (c : fluidGrid)
  if c.in_interior then
    -- Add body forces (accelerations) to the momentum
    c.rhoVelocity_t += c.rho * Flow.bodyForce

    -- Body force contribution to energy equation
    c.rhoEnergy_t += c.rho * L.dot(Flow.bodyForce,c.velocity)
  end
end

ebb Flow.UpdatePD (c : fluidGrid)
  if c.in_interior then
    var divU = L.double(0.0)

    -- compute the divergence of the velocity (trace of the velocity gradient)
    divU = c.velocityGradientX[0] + c.velocityGradientY[1] + c.velocityGradientZ[2]

    -- Compute pressure dilation by multiplying by pressure (assumes homogeneity)
    -- PD = - <u_i P,j> = <Ui,i P >
    c.PD = divU * c.pressure
  end
end

-- Compute viscous fluxes in X direction
ebb Flow.ComputeDissipationX (c : fluidGrid)
  if c.in_interior or c.xneg_depth == 1 then
    -- Consider first boundary element (c.xneg_depth == 1) to define left flux
    -- on first interior cell
    var muFace = 0.5 * (GetDynamicViscosity(c.temperature) +
                        GetDynamicViscosity(c(1,0,0).temperature))
    var velocityFace    = L.vec3d({0.0, 0.0, 0.0})
    var velocityX_YFace = L.double(0.0)
    var velocityX_ZFace = L.double(0.0)
    var velocityY_YFace = L.double(0.0)
    var velocityZ_ZFace = L.double(0.0)

    -- Interpolate velocity and derivatives to face
    velocityFace = 0.5 * ( c.velocity + c(1,0,0).velocity )
    velocityX_YFace = 0.5 * ( c.velocityGradientY[0] +
                              c(1,0,0).velocityGradientY[0] )
    velocityX_ZFace = 0.5 * ( c.velocityGradientZ[0] +
                              c(1,0,0).velocityGradientZ[0] )
    velocityY_YFace = 0.5 * ( c.velocityGradientY[1] +
                              c(1,0,0).velocityGradientY[1] )
    velocityZ_ZFace = 0.5 * ( c.velocityGradientZ[2] +
                              c(1,0,0).velocityGradientZ[2] )

    -- Differentiate at face
    var velocityX_XFace   = L.double(0.0)
    var velocityY_XFace   = L.double(0.0)
    var velocityZ_XFace   = L.double(0.0)
    var temperature_XFace = L.double(0.0)

    velocityX_XFace   = 0.5*( c(1,0,0).velocity[0] - c.velocity[0] )
    velocityY_XFace   = 0.5*( c(1,0,0).velocity[1] - c.velocity[1] )
    velocityZ_XFace   = 0.5*( c(1,0,0).velocity[2] - c.velocity[2] )
    temperature_XFace = 0.5*( c(1,0,0).temperature - c.temperature )

    -- Half cell size due to the 0.5 above
    velocityX_XFace   /= (Grid.xCellWidth*0.5)
    velocityY_XFace   /= (Grid.xCellWidth*0.5)
    velocityZ_XFace   /= (Grid.xCellWidth*0.5)
    temperature_XFace /= (Grid.xCellWidth*0.5)

    -- Tensor components (at face)
    var sigmaXX = muFace * ( 4.0 * velocityX_XFace -
                             2.0 * velocityY_YFace -
                             2.0 * velocityZ_ZFace ) / 3.0
    var sigmaYX = muFace * ( velocityY_XFace + velocityX_YFace )
    var sigmaZX = muFace * ( velocityZ_XFace + velocityX_ZFace )
    var usigma  = velocityFace[0] * sigmaXX +
                  velocityFace[1] * sigmaYX +
                  velocityFace[2] * sigmaZX

    -- Fluxes
    c.dissipationFlux = usigma -- possible just x component?
  end
end

-- Compute viscous fluxes in Y direction
ebb Flow.ComputeDissipationY (c : fluidGrid)
  if c.in_interior or c.yneg_depth == 1 then
    -- Consider first boundary element (c.yneg_depth == 1) to define down flux
    -- on first interior cell
    var muFace = 0.5 * (GetDynamicViscosity(c.temperature) +
                        GetDynamicViscosity(c(0,1,0).temperature))
    var velocityFace    = L.vec3d({0.0, 0.0, 0.0})
    var velocityY_XFace = L.double(0.0)
    var velocityY_ZFace = L.double(0.0)
    var velocityX_XFace = L.double(0.0)
    var velocityZ_ZFace = L.double(0.0)

    -- Interpolate velocity and derivatives to face
    velocityFace = 0.5 * ( c.velocity + c(0,1,0).velocity )
    velocityY_XFace = 0.5 * ( c.velocityGradientX[1] +
                              c(0,1,0).velocityGradientX[1] )
    velocityY_ZFace = 0.5 * ( c.velocityGradientZ[1] +
                              c(0,1,0).velocityGradientZ[1] )
    velocityX_XFace = 0.5 * ( c.velocityGradientX[0] +
                              c(0,1,0).velocityGradientX[0] )
    velocityZ_ZFace = 0.5 * ( c.velocityGradientZ[2] +
                              c(0,1,0).velocityGradientZ[2] )

    -- Differentiate at face
    var velocityX_YFace   = L.double(0.0)
    var velocityY_YFace   = L.double(0.0)
    var velocityZ_YFace   = L.double(0.0)
    var temperature_YFace = L.double(0.0)

    velocityX_YFace   = 0.5*( c(0,1,0).velocity[0] - c.velocity[0] )
    velocityY_YFace   = 0.5*( c(0,1,0).velocity[1] - c.velocity[1] )
    velocityZ_YFace   = 0.5*( c(0,1,0).velocity[2] - c.velocity[2] )
    temperature_YFace = 0.5*( c(0,1,0).temperature - c.temperature )

    -- Half cell size due to the 0.5 above
    velocityX_YFace   /= (Grid.yCellWidth*0.5)
    velocityY_YFace   /= (Grid.yCellWidth*0.5)
    velocityZ_YFace   /= (Grid.yCellWidth*0.5)
    temperature_YFace /= (Grid.yCellWidth*0.5)

    -- Tensor components (at face)
    var sigmaXY = muFace * ( velocityX_YFace + velocityY_XFace )
    var sigmaYY = muFace * ( 4.0 * velocityY_YFace -
                             2.0 * velocityX_XFace -
                             2.0 * velocityZ_ZFace ) / 3.0
    var sigmaZY = muFace * ( velocityZ_YFace + velocityY_ZFace )
    var usigma  = velocityFace[0] * sigmaXY +
                  velocityFace[1] * sigmaYY +
                  velocityFace[2] * sigmaZY

    -- Fluxes
    c.dissipationFlux = usigma
  end
end

-- Compute viscous fluxes in Z direction
ebb Flow.ComputeDissipationZ (c : fluidGrid)
  if c.in_interior or c.zneg_depth == 1 then
    -- Consider first boundary element (c.zneg_depth == 1) to define down flux
    -- on first interior cell
    var muFace = 0.5 * (GetDynamicViscosity(c.temperature) +
                          GetDynamicViscosity(c(0,0,1).temperature))
    var velocityFace    = L.vec3d({0.0, 0.0, 0.0})
    var velocityZ_XFace = L.double(0.0)
    var velocityZ_YFace = L.double(0.0)
    var velocityX_XFace = L.double(0.0)
    var velocityY_YFace = L.double(0.0)

    -- Interpolate velocity and derivatives to face
    velocityFace = 0.5 * ( c.velocity + c(0,0,1).velocity )
    velocityZ_XFace = 0.5 * ( c.velocityGradientX[2] +
                              c(0,0,1).velocityGradientX[2] )
    velocityZ_YFace = 0.5 * ( c.velocityGradientY[2] +
                              c(0,0,1).velocityGradientY[2] )
    velocityX_XFace = 0.5 * ( c.velocityGradientX[0] +
                              c(0,0,1).velocityGradientX[0] )
    velocityY_YFace = 0.5 * ( c.velocityGradientY[1] +
                              c(0,0,1).velocityGradientY[1] )

    -- Differentiate at face
    var velocityX_ZFace   = L.double(0.0)
    var velocityY_ZFace   = L.double(0.0)
    var velocityZ_ZFace   = L.double(0.0)
    var temperature_ZFace = L.double(0.0)

    velocityX_ZFace   = 0.5*( c(0,0,1).velocity[0] - c.velocity[0] )
    velocityY_ZFace   = 0.5*( c(0,0,1).velocity[1] - c.velocity[1] )
    velocityZ_ZFace   = 0.5*( c(0,0,1).velocity[2] - c.velocity[2] )
    temperature_ZFace = 0.5*( c(0,0,1).temperature - c.temperature )

    -- Half cell size due to the 0.5 above
    velocityX_ZFace   /= (Grid.zCellWidth*0.5)
    velocityY_ZFace   /= (Grid.zCellWidth*0.5)
    velocityZ_ZFace   /= (Grid.zCellWidth*0.5)
    temperature_ZFace /= (Grid.zCellWidth*0.5)

    -- Tensor components (at face)
    var sigmaXZ = muFace * ( velocityX_ZFace + velocityZ_XFace )
    var sigmaYZ = muFace * ( velocityY_ZFace + velocityZ_YFace )
    var sigmaZZ = muFace * ( 4.0 * velocityZ_ZFace -
                             2.0 * velocityX_XFace -
                             2.0 * velocityY_YFace ) / 3.0
    var usigma  = velocityFace[0] * sigmaXZ +
                  velocityFace[1] * sigmaYZ +
                  velocityFace[2] * sigmaZZ

    -- Fluxes
    c.dissipationFlux = usigma
  end
end

ebb Flow.UpdateDissipationX (c : fluidGrid)
  if c.in_interior then
    c.dissipation += (c( 0,0,0).dissipationFlux -
                      c(-1,0,0).dissipationFlux)/Grid.xCellWidth
  end
end

ebb Flow.UpdateDissipationY (c : fluidGrid)
  if c.in_interior then
    c.dissipation += (c(0, 0,0).dissipationFlux -
                      c(0,-1,0).dissipationFlux)/Grid.yCellWidth
  end
end

ebb Flow.UpdateDissipationZ (c : fluidGrid)
  if c.in_interior then
    c.dissipation += (c(0,0, 0).dissipationFlux -
                      c(0,0,-1).dissipationFlux)/Grid.zCellWidth
  end
end

ebb Flow.ResetDissipation (c : fluidGrid)
  c.dissipation = 0.0
end

function Flow.UpdateDissipation ()
  fluidGrid:foreach(Flow.ResetDissipation)
  fluidGrid:foreach(Flow.ComputeDissipationX)
  fluidGrid:foreach(Flow.UpdateDissipationX)
  fluidGrid:foreach(Flow.ComputeDissipationY)
  fluidGrid:foreach(Flow.UpdateDissipationY)
  fluidGrid:foreach(Flow.ComputeDissipationZ)
  fluidGrid:foreach(Flow.UpdateDissipationZ)
end

-- WARNING: uniform grid assumption
local ebb CalculateAveragePD (c : fluidGrid)
  if c.in_interior then
    Flow.averagePD += c.PD * Grid.cellVolume
  end
end

local ebb CalculateAverageDissipation (c : fluidGrid)
  if c.in_interior then
    Flow.averageDissipation += c.dissipation * Grid.cellVolume
  end
end

local ebb CalculateAverageK (c : fluidGrid)
  if c.in_interior then
    Flow.averageK += 0.5 * c.rho * L.dot(c.velocity,c.velocity) * Grid.cellVolume
  end
end

function Flow.UpdateTurbulentAverages ()
  fluidGrid:foreach(CalculateAveragePD)
  Flow.averagePD:set(Flow.averagePD:get() / Grid.areaInterior)
  fluidGrid:foreach(CalculateAverageDissipation)
  Flow.averageDissipation:set(Flow.averageDissipation:get()/ Grid.areaInterior)
  fluidGrid:foreach(CalculateAverageK)
  Flow.averageK:set(Flow.averageK:get() / Grid.areaInterior)
end

ebb Flow.AddTurbulentSource (c : fluidGrid)
  if c.in_interior then

    var W   = L.double(0.0)
    var A   = L.double(0.0)
    var G   = L.double(0.0)
    var t_o = L.double(0.0)
    var K_o = L.double(0.0)
    var force = L.vec3d({0.0,0.0,0.0})

    -- Compute W (pressure dilatation term and dissipation)
    W = Flow.averagePD + Flow.averageDissipation

    -- Compute forcing coefficient using gain controller
    -- Inputs: G, t_o, Ko, where G ~ 300.0, t_o ~ L_o / u_o, L_o is domain length,
    -- u_o ~ from Re relationship or sqrt(K_o/rho_o)
    G   = 300.0
    t_o = 3.00889E-06
    K_o = 66.27348

    A =  ( - W - G * ( Flow.averageK - K_o ) / t_o  ) / (2.0 * Flow.averageK)

    -- Compute the turbulent force vector
    force = c.rho * A * c.velocity

    -- Add the forcing terms to the momentum and energy equations
    c.rhoVelocity_t += force
    c.rhoEnergy_t   += L.dot(force,c.velocity)

    -- Store the increment in the average energy source (to be subtracted later)
    Flow.averageFe += L.dot(force,c.velocity) * Grid.cellVolume
  end
end

ebb Flow.AdjustTurbulentSource (c : fluidGrid)
  if c.in_interior then
    -- Remove the average of the forcing term that has been added to the energy
    -- equation so that the flow can reach a statistical steady state.
    -- Note that this has been pre-computed before reaching this kernel (above).
    c.rhoEnergy_t -= Flow.averageFe
  end
end

-- One high level routine that runs all steps
function Flow.AddTurbulentForcing ()
  -- Need to reset these averages somewhere
  Flow.averagePD:set(0.0)
  Flow.averageDissipation:set(0.0)
  Flow.averageFe:set(0.0)
  Flow.averageK:set(0.0)
  fluidGrid:foreach(Flow.UpdatePD)
  Flow.UpdateDissipation()
  -- average PD and EPS
  Flow.UpdateTurbulentAverages()
  -- Compute A & force, f_i
  -- Add rho * A * u_i to momentum, f_i*u_i to energy, accumulate f_i*u_i for average
  fluidGrid:foreach(Flow.AddTurbulentSource)
  -- Update average of the energy source
  Flow.averageFe:set(Flow.averageFe:get()/Grid.areaInterior)
  -- Subtract <f_e> from energy
  fluidGrid:foreach(Flow.AdjustTurbulentSource)
end

-------------------
-- Update functions
-------------------

-- Update flow variables using derivatives
-- Assumes 4th-order Runge-Kutta
ebb Flow.UpdateVars (c : fluidGrid)
  var deltaTime = Integrator.deltaTime
  if Integrator.stage == 1 then
    c.rho_new += (1.0/6.0) * deltaTime * c.rho_t
    c.rho = c.rho_old + 0.5 * deltaTime * c.rho_t
    c.rhoVelocity_new += (1.0/6.0) * deltaTime * c.rhoVelocity_t
    c.rhoVelocity = c.rhoVelocity_old + 0.5 * deltaTime * c.rhoVelocity_t
    c.rhoEnergy_new += (1.0/6.0) * deltaTime * c.rhoEnergy_t
    c.rhoEnergy = c.rhoEnergy_old + 0.5 * deltaTime * c.rhoEnergy_t
  elseif Integrator.stage == 2 then
    c.rho_new += (1.0/3.0) * deltaTime * c.rho_t
    c.rho = c.rho_old + 0.5 * deltaTime * c.rho_t
    c.rhoVelocity_new += (1.0/3.0) * deltaTime * c.rhoVelocity_t
    c.rhoVelocity = c.rhoVelocity_old + 0.5 * deltaTime * c.rhoVelocity_t
    c.rhoEnergy_new += (1.0/3.0) * deltaTime * c.rhoEnergy_t
    c.rhoEnergy = c.rhoEnergy_old + 0.5 * deltaTime * c.rhoEnergy_t
  elseif Integrator.stage == 3 then
    c.rho_new += (1.0/3.0) * deltaTime * c.rho_t
    c.rho = c.rho_old + 1.0 * deltaTime * c.rho_t
    c.rhoVelocity_new += (1.0/3.0) * deltaTime * c.rhoVelocity_t
    c.rhoVelocity = c.rhoVelocity_old + 1.0 * deltaTime * c.rhoVelocity_t
    c.rhoEnergy_new += (1.0/3.0) * deltaTime * c.rhoEnergy_t
    c.rhoEnergy = c.rhoEnergy_old + 1.0 * deltaTime * c.rhoEnergy_t
  else -- Integrator.stage == 4
    c.rho = c.rho_new + (1.0/6.0) * deltaTime * c.rho_t
    c.rhoVelocity = c.rhoVelocity_new + (1.0/6.0) * deltaTime * c.rhoVelocity_t
    c.rhoEnergy = c.rhoEnergy_new + (1.0/6.0) * deltaTime * c.rhoEnergy_t
  end
end

ebb Flow.UpdateAuxiliaryVelocity (c : fluidGrid)
  if c.in_interior then
    var velocity = c.rhoVelocity / c.rho
    c.velocity = velocity
    c.kineticEnergy = 0.5 *  c.rho * L.dot(velocity,velocity)
  end
end

-- Helper function for updating the ghost fields to minimize repeated code
local ebb UpdateGhostFieldsHelper (c_bnd, c_int, sign, bnd_velocity, bnd_temperature)
  -- Temporary variables for computing new halo state
  var rho         = L.double(0.0)
  var temp_wall   = L.double(0.0)
  var temperature = L.double(0.0)
  var velocity    = L.vec3d({0.0, 0.0, 0.0})
  -- Compute the Cv for updating the Energy equation
  var cv = Flow.gasConstant / (Flow.gamma - 1.0)
  -- Compute the new velocity (including any wall conditions)
  velocity = L.times(c_int.rhoVelocity/c_int.rho, sign) + bnd_velocity
  -- Compute the temperature for the halo cell (possibly adiabatic/isothermal)
  temp_wall = c_int.temperature
  if bnd_temperature > 0.0 then
    temp_wall = bnd_temperature
  end
  temperature = 2.0*temp_wall - c_int.temperature
  -- Recompute the density in the halo in case of a temperature change
  -- Pressure is a zero-order extrapolation
  rho = c_int.pressure / ( Flow.gasConstant * temperature )
  -- Update the boundary cell based on the values in the matching interior cell
  c_bnd.rhoBoundary         =  rho
  c_bnd.rhoVelocityBoundary =  rho*velocity
  c_bnd.rhoEnergyBoundary   =  rho * (cv * temperature + 0.5*L.dot(velocity,velocity))
  c_bnd.velocityBoundary    =  velocity
  c_bnd.pressureBoundary    =  c_int.pressure
  c_bnd.temperatureBoundary =  temperature
end

ebb Flow.UpdateGhostFieldsStep1 (c : fluidGrid)
  if c.xneg_depth > 0 then
    UpdateGhostFieldsHelper(c, c( 1,0,0), BC.xSign, BC.xNegVelocity, BC.xNegTemperature)
  end
  if c.xpos_depth > 0 then
    UpdateGhostFieldsHelper(c, c(-1,0,0), BC.xSign, BC.xPosVelocity, BC.xPosTemperature)
  end
  if c.yneg_depth > 0 then
    UpdateGhostFieldsHelper(c, c(0, 1,0), BC.ySign, BC.yNegVelocity, BC.yNegTemperature)
  end
  if c.ypos_depth > 0 then
    UpdateGhostFieldsHelper(c, c(0,-1,0), BC.ySign, BC.yPosVelocity, BC.yPosTemperature)
  end
  if c.zneg_depth > 0 then
    UpdateGhostFieldsHelper(c, c(0,0, 1), BC.zSign, BC.zNegVelocity, BC.zNegTemperature)
  end
  if c.zpos_depth > 0 then
    UpdateGhostFieldsHelper(c, c(0,0,-1), BC.zSign, BC.zPosVelocity, BC.zPosTemperature)
  end
end

ebb Flow.UpdateGhostFieldsStep2 (c : fluidGrid)
  if c.in_boundary then
    c.rho         = c.rhoBoundary
    c.rhoVelocity = c.rhoVelocityBoundary
    c.rhoEnergy   = c.rhoEnergyBoundary
    c.pressure    = c.pressureBoundary
    c.temperature = c.temperatureBoundary
  end
end

function Flow.UpdateGhost ()
  fluidGrid:foreach(Flow.UpdateGhostFieldsStep1)
  fluidGrid:foreach(Flow.UpdateGhostFieldsStep2)
end

-- Helper function for updating the ghost fields to minimize repeated code
local ebb UpdateGhostThermodynamicsHelper (c_bnd, c_int, bnd_temperature)
  -- Temporary variables for computing new halo state
  var temp_wall   = L.double(0.0)
  var temperature = L.double(0.0)
  -- Compute the temperature for the halo cell (possibly adiabatic/isothermal)
  temp_wall = c_int.temperature
  if bnd_temperature > 0.0 then
    temp_wall = bnd_temperature
  end
  temperature = 2.0*temp_wall - c_int.temperature
  -- Update the boundary cell based on the values in the matching interior cell
  c_bnd.pressureBoundary    = c_int.pressure
  c_bnd.temperatureBoundary = temperature
end

ebb Flow.UpdateGhostThermodynamicsStep1 (c : fluidGrid)
  if c.xneg_depth > 0 then
    UpdateGhostThermodynamicsHelper(c, c( 1,0,0), BC.xNegTemperature)
  end
  if c.xpos_depth > 0 then
    UpdateGhostThermodynamicsHelper(c, c(-1,0,0), BC.xPosTemperature)
  end
  if c.yneg_depth > 0 then
    UpdateGhostThermodynamicsHelper(c, c(0, 1,0), BC.yNegTemperature)
  end
  if c.ypos_depth > 0 then
    UpdateGhostThermodynamicsHelper(c, c(0,-1,0), BC.yPosTemperature)
  end
  if c.zneg_depth > 0 then
    UpdateGhostThermodynamicsHelper(c, c(0,0, 1), BC.zNegTemperature)
  end
  if c.zpos_depth > 0 then
    UpdateGhostThermodynamicsHelper(c, c(0,0,-1), BC.zPosTemperature)
  end
end

ebb Flow.UpdateGhostThermodynamicsStep2 (c : fluidGrid)
  if c.in_boundary then
    c.pressure    = c.pressureBoundary
    c.temperature = c.temperatureBoundary
  end
end

function Flow.UpdateGhostThermodynamics ()
  fluidGrid:foreach(Flow.UpdateGhostThermodynamicsStep1)
  fluidGrid:foreach(Flow.UpdateGhostThermodynamicsStep2)
end

-- Helper function for updating the ghost fields to minimize repeated code
local ebb UpdateGhostVelocityHelper (c_bnd, c_int, sign, bnd_velocity)
  -- Update the boundary cell based on the values in the matching interior cell
  c_bnd.velocityBoundary = L.times(c_int.velocity, sign) + bnd_velocity
end

ebb Flow.UpdateGhostVelocityStep1 (c : fluidGrid)
  if c.xneg_depth > 0 then
    UpdateGhostVelocityHelper(c, c( 1,0,0), BC.xSign, BC.xNegVelocity)
  end
  if c.xpos_depth > 0 then
    UpdateGhostVelocityHelper(c, c(-1,0,0), BC.xSign, BC.xPosVelocity)
  end
  if c.yneg_depth > 0 then
    UpdateGhostVelocityHelper(c, c(0, 1,0), BC.ySign, BC.yNegVelocity)
  end
  if c.ypos_depth > 0 then
    UpdateGhostVelocityHelper(c, c(0,-1,0), BC.ySign, BC.yPosVelocity)
  end
  if c.zneg_depth > 0 then
    UpdateGhostVelocityHelper(c, c(0,0, 1), BC.zSign, BC.zNegVelocity)
  end
  if c.zpos_depth > 0 then
    UpdateGhostVelocityHelper(c, c(0,0,-1), BC.zSign, BC.zPosVelocity)
  end
end

ebb Flow.UpdateGhostVelocityStep2 (c : fluidGrid)
  if c.in_boundary then
    c.velocity = c.velocityBoundary
  end
end

function Flow.UpdateGhostVelocity ()
  fluidGrid:foreach(Flow.UpdateGhostVelocityStep1)
  fluidGrid:foreach(Flow.UpdateGhostVelocityStep2)
end

-- Helper function for updating the conservatives to minimize repeated code
local ebb UpdateGhostConservedHelper (c_bnd, c_int, sign, bnd_velocity, bnd_temperature)
  -- Temporary variables for computing new halo state
  var rho         = L.double(0.0)
  var temp_wall   = L.double(0.0)
  var temperature = L.double(0.0)
  var velocity    = L.vec3d({0.0, 0.0, 0.0})
  -- Compute the Cv for updating the Energy equation
  var cv = Flow.gasConstant / (Flow.gamma - 1.0)
  -- Compute the new velocity (including any wall conditions)
  var velocity = L.vec3d({0.0, 0.0, 0.0})
  velocity = L.times(c_int.rhoVelocity/c_int.rho, sign) + bnd_velocity
  -- Compute the temperature for the halo cell (possibly adiabatic/isothermal)
  temp_wall = c_int.temperature
  if bnd_temperature > 0.0 then
    temp_wall = bnd_temperature
  end
  temperature = 2.0*temp_wall - c_int.temperature
  -- Recompute the density in the halo in case of a temperature change
  -- Pressure is a zero-order extrapolation
  rho = c_int.pressure / ( Flow.gasConstant * temperature )
  -- Update the boundary cell based on the values in the matching interior cell
  c_bnd.rhoBoundary         = rho
  c_bnd.rhoVelocityBoundary = rho*velocity
  c_bnd.rhoEnergyBoundary   = rho * (cv * temperature + 0.5*L.dot(velocity,velocity))
end

ebb Flow.UpdateGhostConservedStep1 (c : fluidGrid)
  if c.xneg_depth > 0 then
    UpdateGhostConservedHelper(c, c( 1,0,0), BC.xSign, BC.xNegVelocity, BC.xNegTemperature)
  end
  if c.xpos_depth > 0 then
    UpdateGhostConservedHelper(c, c(-1,0,0), BC.xSign, BC.xPosVelocity, BC.xPosTemperature)
  end
  if c.yneg_depth > 0 then
    UpdateGhostConservedHelper(c, c(0, 1,0), BC.ySign, BC.yNegVelocity, BC.yNegTemperature)
  end
  if c.ypos_depth > 0 then
    UpdateGhostConservedHelper(c, c(0,-1,0), BC.ySign, BC.yPosVelocity, BC.yPosTemperature)
  end
  if c.zneg_depth > 0 then
    UpdateGhostConservedHelper(c, c(0,0, 1), BC.zSign, BC.zNegVelocity, BC.zNegTemperature)
  end
  if c.zpos_depth > 0 then
    UpdateGhostConservedHelper(c, c(0,0,-1), BC.zSign, BC.zPosVelocity, BC.zPosTemperature)
  end
end

ebb Flow.UpdateGhostConservedStep2 (c : fluidGrid)
  if c.in_boundary then
    c.rho         = c.rhoBoundary
    c.rhoVelocity = c.rhoVelocityBoundary
    c.rhoEnergy   = c.rhoEnergyBoundary
  end
end

function Flow.UpdateGhostConserved ()
  fluidGrid:foreach(Flow.UpdateGhostConservedStep1)
  fluidGrid:foreach(Flow.UpdateGhostConservedStep2)
end

ebb Flow.UpdateAuxiliaryThermodynamics (c : fluidGrid)
  if c.in_interior then
    var kineticEnergy = 0.5 * c.rho * L.dot(c.velocity,c.velocity)
    var pressure  = (Flow.gamma - 1.0) *( c.rhoEnergy - kineticEnergy )
    c.pressure    = pressure
    c.temperature = pressure / ( Flow.gasConstant * c.rho )
  end
end

---------------------
-- Velocity gradients
---------------------

-- WARNING: non-uniform grid assumption
ebb Flow.ComputeVelocityGradientAll (c : fluidGrid)
  if c.in_interior then
    c.velocityGradientX = 0.5*(c(1,0,0).velocity - c(-1,0,0).velocity)/Grid.xCellWidth
    c.velocityGradientY = 0.5*(c(0,1,0).velocity - c(0,-1,0).velocity)/Grid.yCellWidth
    c.velocityGradientZ = 0.5*(c(0,0,1).velocity - c(0,0,-1).velocity)/Grid.zCellWidth
  end
end

-- Helper function for updating the boundary gradients to minimize repeated code
local ebb UpdateGhostVelocityGradientHelper (c_bnd, c_int, sign)
  -- Apply sign change and copy gradients from interior to boundary
  c_bnd.velocityGradientXBoundary = L.times(sign, c_int.velocityGradientX)
  c_bnd.velocityGradientYBoundary = L.times(sign, c_int.velocityGradientY)
  c_bnd.velocityGradientZBoundary = L.times(sign, c_int.velocityGradientZ)
end

ebb Flow.UpdateGhostVelocityGradientStep1 (c : fluidGrid)
  if c.xneg_depth > 0 then
    UpdateGhostVelocityGradientHelper(c, c( 1,0,0), BC.xSign)
  end
  if c.xpos_depth > 0 then
    UpdateGhostVelocityGradientHelper(c, c(-1,0,0), BC.xSign)
  end
  if c.yneg_depth > 0 then
    UpdateGhostVelocityGradientHelper(c, c(0, 1,0), BC.ySign)
  end
  if c.ypos_depth > 0 then
    UpdateGhostVelocityGradientHelper(c, c(0,-1,0), BC.ySign)
  end
  if c.zneg_depth > 0 then
    UpdateGhostVelocityGradientHelper(c, c(0,0, 1), BC.zSign)
  end
  if c.zpos_depth > 0 then
    UpdateGhostVelocityGradientHelper(c, c(0,0,-1), BC.zSign)
  end
end

ebb Flow.UpdateGhostVelocityGradientStep2 (c : fluidGrid)
  if c.in_boundary then
    c.velocityGradientX = c.velocityGradientXBoundary
    c.velocityGradientY = c.velocityGradientYBoundary
    c.velocityGradientZ = c.velocityGradientZBoundary
  end
end

local ebb CalculateConvectiveSpectralRadius (c : fluidGrid)
  -- Convective spectral radii
  -- WARNING: uniform grid assumption
  c.convectiveSpectralRadius =
   (L.fabs(c.velocity[0])/Grid.xCellWidth  +
    L.fabs(c.velocity[1])/Grid.yCellWidth  +
    L.fabs(c.velocity[2])/Grid.zCellWidth  +
    GetSoundSpeed(c.temperature) * L.sqrt(Grid.dXYZInverseSquare))
  Integrator.maxConvectiveSpectralRadius max= c.convectiveSpectralRadius
end

local ebb CalculateViscousSpectralRadius (c : fluidGrid)
  -- Viscous spectral radii (including sgs model component)
  var dynamicViscosity = GetDynamicViscosity(c.temperature)
  var eddyViscosity = c.sgsEddyViscosity
  c.viscousSpectralRadius =
   (2.0 * ( dynamicViscosity + eddyViscosity ) /
    c.rho * Grid.dXYZInverseSquare) * 4.0
  Integrator.maxViscousSpectralRadius max= c.viscousSpectralRadius
end

local ebb CalculateHeatConductionSpectralRadius (c : fluidGrid)
  var dynamicViscosity  = GetDynamicViscosity(c.temperature)
  -- Heat conduction spectral radii (including sgs model component)
  var cv    = Flow.gasConstant / (Flow.gamma - 1.0)
  var cp    = Flow.gamma * cv
  var kappa = cp / Flow.prandtl *  dynamicViscosity
  c.heatConductionSpectralRadius =
     ((kappa + c.sgsEddyKappa) / (cv * c.rho) * Grid.dXYZInverseSquare) * 4.0
  Integrator.maxHeatConductionSpectralRadius max= c.heatConductionSpectralRadius
end

function Flow.CalculateSpectralRadii ()
  fluidGrid:foreach(CalculateConvectiveSpectralRadius)
  fluidGrid:foreach(CalculateViscousSpectralRadius)
  fluidGrid:foreach(CalculateHeatConductionSpectralRadius)
end

-------------
-- Statistics
-------------

local ebb CalculateAveragePressure (c : fluidGrid)
  if c.in_interior then
    Flow.averagePressure += c.pressure * Grid.cellVolume
  end
end

local ebb CalculateAverageTemperature (c : fluidGrid)
  if c.in_interior then
    Flow.averageTemperature += c.temperature * Grid.cellVolume
  end
end

local ebb CalculateAverageKineticEnergy (c : fluidGrid)
  if c.in_interior then
    Flow.averageKineticEnergy += c.kineticEnergy * Grid.cellVolume
  end
end

local ebb CalculateMinTemperature (c : fluidGrid)
  if c.in_interior then
    Flow.minTemperature min= c.temperature
  end
end

local ebb CalculateMaxTemperature (c : fluidGrid)
  if c.in_interior then
    Flow.maxTemperature max= c.temperature
  end
end

function Flow.IntegrateQuantities ()
  fluidGrid:foreach(CalculateAveragePressure)
  fluidGrid:foreach(CalculateAverageTemperature)
  fluidGrid:foreach(CalculateAverageKineticEnergy)
  fluidGrid:foreach(CalculateMinTemperature)
  fluidGrid:foreach(CalculateMaxTemperature)
end

------------
-- PARTICLES
------------

ebb Particles.LocateInCells (p : particles)
  p.cell = fluidGrid.locate(p.position)
end

-- Locate particles in cells
function Particles.Locate ()
  particles:foreach(Particles.LocateInCells)
end

-- Initialize temporaries for time stepper
ebb Particles.InitializeTemporaries (p : particles)
  p.position_old    = p.position
  p.velocity_old    = p.velocity
  p.temperature_old = p.temperature
  p.position_new    = p.position
  p.velocity_new    = p.velocity
  p.temperature_new = p.temperature
end

----------------
-- Flow Coupling
----------------

-- Initialize time derivative for each stage of time stepper
ebb Particles.InitializeTimeDerivatives (p : particles)
  p.position_t = L.vec3d({0.0, 0.0, 0.0})
  p.velocity_t = L.vec3d({0.0, 0.0, 0.0})
  p.temperature_t = L.double(0)
end

-- Update particle fields based on flow fields
ebb Particles.AddFlowCoupling (p: particles)

  -- WARNING: assumes we have already located particles

  -- Trilinear interpolation for the flow quantities
  var flowVelocity = InterpolateTriVelocity(p.cell, p.position)
  var flowTemperature = InterpolateTriTemp(p.cell, p.position)
  var flowDynamicViscosity = GetDynamicViscosity(flowTemperature)

  -- Update the particle position using the current velocity
  p.position_t += p.velocity

  -- Relaxation time for small particles
  -- - particles Reynolds number (set to zero for Stokesian)
  var particleReynoldsNumber = 0.0
  --(p.density * Norm(flowVelocity - p.velocity) * p.diameter) / flowDynamicViscosity
  var relaxationTime =
    ( p.density * L.pow(p.diameter,2) / (18.0 * flowDynamicViscosity) ) /
    ( 1.0 + 0.15 * L.pow(particleReynoldsNumber,0.687) )
  p.deltaVelocityOverRelaxationTime =
    (flowVelocity - p.velocity) / relaxationTime
  p.deltaTemperatureTerm =
    pi * L.pow(p.diameter, 2) * Particles.convectiveCoeff *
    (flowTemperature - p.temperature)

  -- Update the particle velocity and temperature
  p.velocity_t += p.deltaVelocityOverRelaxationTime
  p.temperature_t += p.deltaTemperatureTerm / (p.mass * Particles.heatCapacity)

end

--------------
-- Body forces
--------------

ebb Particles.AddBodyForces (p : particles)
  p.velocity_t += Particles.bodyForce
end

------------
-- Radiation
------------

-- Update particle variables using derivatives
ebb Particles.UpdateVars (p : particles)
  var deltaTime = Integrator.deltaTime
  if Integrator.stage == 1 then
    p.position_new += (1.0/6.0) * deltaTime * p.position_t
    p.position = p.position_old + 0.5 * deltaTime * p.position_t
    p.velocity_new += (1.0/6.0) * deltaTime * p.velocity_t
    p.velocity = p.velocity_old + 0.5 * deltaTime * p.velocity_t
    p.temperature_new += (1.0/6.0) * deltaTime * p.temperature_t
    p.temperature = p.temperature_old + 0.5 * deltaTime * p.temperature_t
  elseif Integrator.stage == 2 then
    p.position_new += (1.0/3.0) * deltaTime * p.position_t
    p.position = p.position_old + 0.5 * deltaTime * p.position_t
    p.velocity_new += (1.0/3.0) * deltaTime * p.velocity_t
    p.velocity = p.velocity_old + 0.5 * deltaTime * p.velocity_t
    p.temperature_new += (1.0/3.0) * deltaTime * p.temperature_t
    p.temperature = p.temperature_old + 0.5 * deltaTime * p.temperature_t
  elseif Integrator.stage == 3 then
    p.position_new += (1.0/3.0) * deltaTime * p.position_t
    p.position = p.position_old + 1.0 * deltaTime * p.position_t
    p.velocity_new += (1.0/3.0) * deltaTime * p.velocity_t
    p.velocity = p.velocity_old + 1.0 * deltaTime * p.velocity_t
    p.temperature_new += (1.0/3.0) * deltaTime * p.temperature_t
    p.temperature = p.temperature_old + 1.0 * deltaTime * p.temperature_t
  else -- Integrator.stage == 4
    p.position = p.position_new + (1.0/6.0) * deltaTime * p.position_t
    p.velocity = p.velocity_new + (1.0/6.0) * deltaTime * p.velocity_t
    p.temperature = p.temperature_new + (1.0/6.0) * deltaTime * p.temperature_t
  end
end

ebb Particles.UpdateAuxiliaryStep1 (p : particles)

  -- Initialize position and velocity before we check for wall collisions

  p.position_ghost[0]   = p.position[0]
  p.position_ghost[1]   = p.position[1]
  p.position_ghost[2]   = p.position[2]
  p.velocity_ghost[0]   = p.velocity[0]
  p.velocity_ghost[1]   = p.velocity[1]
  p.velocity_ghost[2]   = p.velocity[2]
  p.velocity_t_ghost[0] = p.velocity_t[0]
  p.velocity_t_ghost[1] = p.velocity_t[1]
  p.velocity_t_ghost[2] = p.velocity_t[2]

  -- Check here for particles exiting the domain. For periodic
  -- boundaries, the particle is transported to the matching periodic
  -- face. For symmetry or wall boundaries, an elastic collision is
  -- assumed. To start, the collision is perfectly elastic.

  -- Left X boundary
  if p.position[0] < Grid.xOrigin then
    if BC.xBCLeftParticles == ParticleBC.Permeable then
      p.position_ghost[0] = p.position[0] + Grid.xWidth
    elseif BC.xBCLeftParticles == ParticleBC.Solid then

      -- Set the position to be on the wall
      p.position_ghost[0] = Grid.xOrigin

      -- Apply an impulse to kick particle away from the wall
      var impulse = - (1.0+Particles.restitutionCoeff) * p.velocity[0]
      if impulse <= 0 then
        p.velocity_ghost[0] += impulse
      end

      -- Add a contact force in case particle rests on the wall
      var contact_force = -1.0*p.velocity_t[0]

      -- To prevent sticky walls, only add contact force if current
      -- force would push the particle through the wall
      if contact_force > 0 then
        p.velocity_t_ghost[0] += contact_force
      end

    else L.assert(false) end
  end

  -- Right X boundary
  if p.position[0] > Grid.xOrigin + Grid.xWidth then
    if BC.xBCRightParticles == ParticleBC.Permeable then
      p.position_ghost[0] = p.position[0] - Grid.xWidth
    elseif BC.xBCRightParticles == ParticleBC.Solid then

      -- Set the position to be on the wall
      p.position_ghost[0] = Grid.xOrigin + Grid.xWidth

      -- Apply an impulse to kick particle away from the wall
      var impulse = - (1.0+Particles.restitutionCoeff) * p.velocity[0]
      if impulse >= 0 then
        p.velocity_ghost[0] += impulse
      end

      -- Add a contact force in case particle rests on the wall
      var contact_force = -1.0*p.velocity_t[0]

      -- To prevent sticky walls, only add contact force if current
      -- force would push the particle through the wall
      if contact_force < 0 then
        p.velocity_t_ghost[0] += contact_force
      end

    else L.assert(false) end
  end

  -- Left Y boundary
  if p.position[1] < Grid.yOrigin then
    if BC.yBCLeftParticles == ParticleBC.Permeable then
      p.position_ghost[1] = p.position[1] + Grid.yWidth
    elseif BC.yBCLeftParticles == ParticleBC.Solid then

      -- Set the position to be on the wall
      p.position_ghost[1] = Grid.yOrigin

      -- Apply an impulse to kick particle away from the wall
      var impulse = - (1.0+Particles.restitutionCoeff) * p.velocity[1]
      if impulse <= 0 then
        p.velocity_ghost[1] += impulse
      end

      -- Add a contact force in case particle rests on the wall
      var contact_force = -1.0*p.velocity_t[1]

      -- To prevent sticky walls, only add contact force if current
      -- force would push the particle through the wall
      if contact_force > 0 then
        p.velocity_t_ghost[1] += contact_force
      end

    else L.assert(false) end
  end

  -- Right Y boundary
  if p.position[1] > Grid.yOrigin + Grid.yWidth then
    if BC.yBCRightParticles == ParticleBC.Permeable then
      p.position_ghost[1] = p.position[1] - Grid.yWidth
    elseif BC.yBCRightParticles == ParticleBC.Solid then

      -- Set the position to be on the wall
      p.position_ghost[1] = Grid.yOrigin + Grid.yWidth

      -- Apply an impulse to kick particle away from the wall
      var impulse = - (1.0+Particles.restitutionCoeff) * p.velocity[1]
      if impulse >= 0 then
        p.velocity_ghost[1] += impulse
      end

      -- Add a contact force in case particle rests on the wall
      var contact_force = -1.0*p.velocity_t[1]

      -- To prevent sticky walls, only add contact force if current
      -- force would push the particle through the wall
      if contact_force < 0 then
        p.velocity_t_ghost[1] += contact_force
      end

    else L.assert(false) end
  end

  -- Left Z boundary
  if p.position[2] < Grid.zOrigin then
    if BC.zBCLeftParticles == ParticleBC.Permeable then
      p.position_ghost[2] = p.position[2] + Grid.zWidth
    elseif BC.zBCLeftParticles == ParticleBC.Solid then

      -- Set the position to be on the wall
      p.position_ghost[2] = Grid.zOrigin

      -- Apply an impulse to kick particle away from the wall
      var impulse = - (1.0+Particles.restitutionCoeff) * p.velocity[2]
      if impulse <= 0 then
        p.velocity_ghost[2] += impulse
      end

      -- Add a contact force in case particle rests on the wall
      var contact_force = -1.0*p.velocity_t[2]

      -- To prevent sticky walls, only add contact force if current
      -- force would push the particle through the wall
      if contact_force > 0 then
        p.velocity_t_ghost[2] += contact_force
      end

    else L.assert(false) end
  end

  -- Right Z boundary
  if p.position[2] > Grid.zOrigin + Grid.zWidth then
    if BC.zBCRightParticles == ParticleBC.Permeable then
      p.position_ghost[2] = p.position[2] - Grid.zWidth
    elseif BC.zBCRightParticles == ParticleBC.Solid then

      -- Set the position to be on the wall
      p.position_ghost[2] = Grid.zOrigin + Grid.zWidth

      -- Apply an impulse to kick particle away from the wall
      var impulse = - (1.0+Particles.restitutionCoeff) * p.velocity[2]
      if impulse >= 0 then
        p.velocity_ghost[2] += impulse
      end

      -- Add a contact force in case particle rests on the wall
      var contact_force = -1.0*p.velocity_t[2]

      -- To prevent sticky walls, only add contact force if current
      -- force would push the particle through the wall
      if contact_force < 0 then
        p.velocity_t_ghost[2] += contact_force
      end

    else L.assert(false) end
  end

end

ebb Particles.UpdateAuxiliaryStep2 (p : particles)
  p.position   = p.position_ghost
  p.velocity   = p.velocity_ghost
  p.velocity_t = p.velocity_t_ghost
end

if Radiation.TYPE == 'DOM' then

  ebb Radiation.InitializeCell (c : domGrid)
    for m = 0,Radiation.NUM_ANGLES do
      c.I_1[m]     = 0.0
      c.I_2[m]     = 0.0
      c.I_3[m]     = 0.0
      c.I_4[m]     = 0.0
      c.I_5[m]     = 0.0
      c.I_6[m]     = 0.0
      c.I_7[m]     = 0.0
      c.I_8[m]     = 0.0
      c.Iiter_1[m] = 0.0
      c.Iiter_2[m] = 0.0
      c.Iiter_3[m] = 0.0
      c.Iiter_4[m] = 0.0
      c.Iiter_5[m] = 0.0
      c.Iiter_6[m] = 0.0
      c.Iiter_7[m] = 0.0
      c.Iiter_8[m] = 0.0
    end
    c.G = 0.0
    c.S = 0.0
  end

  ebb Radiation.ClearAccumulators (c : domGrid)
    c.acc_d2 = 0.0
    c.acc_d2t4 = 0.0
  end

  ebb Radiation.AccumulateParticleValues (p : particles)
    p.cell.to_Radiation.acc_d2 +=
      L.pow(p.diameter,2.0)
    p.cell.to_Radiation.acc_d2t4 +=
      L.pow(p.diameter,2.0) * L.pow(p.temperature,4.0)
  end
  Radiation.AccumulateParticleValues._MANUAL_PARAL = true

  ebb Radiation.UpdateFieldValues (c : domGrid)
    c.sigma = c.acc_d2 * pi
      * (Radiation.qa + Radiation.qs)
      / (4.0 * Radiation.cellVolume)
    if c.acc_d2 == 0.0 then
      c.Ib = 0.0
    else
      c.Ib = SB * c.acc_d2t4 / (pi * c.acc_d2)
    end
  end

  ebb Particles.AbsorbRadiation (p : particles)
    var t4 = L.pow(p.temperature,4.0)
    var alpha = pi * Radiation.qa * L.pow(p.diameter,2.0)
      * (p.cell.to_Radiation.G - 4.0 * SB * t4) / 4.0
    p.temperature_t += alpha / (p.mass * Particles.heatCapacity)
  end
  Particles.AbsorbRadiation._MANUAL_PARAL = true

end

------------
-- Collector
------------

ebb Particles.DeleteEscapingParticles (p: particles)
  var min_x = Grid.xRealOrigin
  var max_x = Grid.xRealOrigin + Grid.xRealWidth
  var min_y = Grid.yRealOrigin
  var max_y = Grid.yRealOrigin + Grid.yRealWidth
  var min_z = Grid.zRealOrigin
  var max_z = Grid.zRealOrigin + Grid.zRealWidth
  var pos = p.position
  if (pos[0] > max_x or pos[0] < min_x  or
      pos[1] > max_y or pos[1] < min_y  or
      pos[2] > max_z or pos[2] < min_z) then
    delete(p)
    Particles.number -= L.int64(1)
  end
end

-- Particle collector
function Particles.Collect ()
  -- For now, delete anything that leaves the domain.
  particles:foreach(Particles.DeleteEscapingParticles)
end

-------------
-- Statistics
-------------

ebb Particles.IntegrateQuantities (p : particles)
  Particles.averageTemperature += p.temperature
end

ebb Particles.CalculateNumber (p : particles)
  Particles.number += L.int64(1)
end

-----------------------------------------------------------------------------
--[[                                MAIN FUNCTIONS                       ]]--
-----------------------------------------------------------------------------

-------
-- FLOW
-------

function Flow.InitializePrimitives ()
  M.IF(M.EQ(Flow.initCase, FlowInitCase.Uniform))
    fluidGrid:foreach(Flow.InitializeUniform)
  M.END()
  M.IF(M.EQ(Flow.initCase, FlowInitCase.TaylorGreen2DVortex))
    fluidGrid:foreach(Flow.InitializeTaylorGreen2D)
  M.END()
  M.IF(M.EQ(Flow.initCase, FlowInitCase.TaylorGreen3DVortex))
    fluidGrid:foreach(Flow.InitializeTaylorGreen3D)
  M.END()
  M.IF(M.EQ(Flow.initCase, FlowInitCase.Perturbed))
    fluidGrid:foreach(Flow.InitializePerturbed)
  M.END()
  M.IF(M.EQ(Flow.initCase, FlowInitCase.Restart))
    fluidGrid:Load({'rho','pressure','velocity'},
                   'restart_fluid_%d.hdf',
                   Integrator.restartIter)
  M.END()
end

function Flow.UpdateGhostVelocityGradient ()
  fluidGrid:foreach(Flow.UpdateGhostVelocityGradientStep1)
  fluidGrid:foreach(Flow.UpdateGhostVelocityGradientStep2)
end

-- Routine that computes the inviscid flux through the face of
-- any two adjacent cells with a centered scheme. The left cell (c_l),
-- right cell (c_r), and coordinate direction (x = 0, y = 1, or z = 2)
-- are the inputs.
local function mkCenteredInviscidFlux (direction)
  local ebb CenteredInviscidFlux (c_l, c_r)

    -- Diagonal terms of inviscid flux
    var rhoFactorDiagonal         = L.double(0.0)
    var rhoVelocityFactorDiagonal = L.vec3d({0.0, 0.0, 0.0})
    var rhoEnergyFactorDiagonal   = L.double(0.0)
    var fpdiag                    = L.double(0.0)

    rhoFactorDiagonal = 0.5 * ( c_l.rho * c_l.velocity[direction] +
                                c_r.rho * c_r.velocity[direction] )
    rhoVelocityFactorDiagonal = 0.5 *
                              ( c_l.rhoVelocity *
                                c_l.velocity[direction] +
                                c_r.rhoVelocity *
                                c_r.velocity[direction] )
    rhoEnergyFactorDiagonal = 0.5 *
                            ( c_l.rhoEnthalpy *
                              c_l.velocity[direction] +
                              c_r.rhoEnthalpy *
                              c_r.velocity[direction] )
    fpdiag += 0.5 * ( c_l.pressure + c_r.pressure )

    -- Skewed terms
    var rhoFactorSkew         = L.double(0.0)
    var rhoVelocityFactorSkew = L.vec3d({0.0, 0.0, 0.0})
    var rhoEnergyFactorSkew   = L.double(0.0)
    var tmp                   = L.double(0.0)

    tmp = 0.5 * c_r.velocity[direction]

    rhoFactorSkew         += c_l.rho * tmp
    rhoVelocityFactorSkew += c_l.rhoVelocity * tmp
    rhoEnergyFactorSkew   += c_l.rhoEnthalpy * tmp

    tmp = 0.5 * c_l.velocity[direction]

    rhoFactorSkew         += c_r.rho * tmp
    rhoVelocityFactorSkew += c_r.rhoVelocity * tmp
    rhoEnergyFactorSkew   += c_r.rhoEnthalpy * tmp

    -- Compute fluxes with prescribed splitting
    var s = Integrator.SPLIT
    var rhoFlux_temp         = s * rhoFactorDiagonal +
                              (1-s) * rhoFactorSkew
    var rhoVelocityFlux_temp = s * rhoVelocityFactorDiagonal +
                              (1-s) * rhoVelocityFactorSkew
    var rhoEnergyFlux_temp   = s * rhoEnergyFactorDiagonal +
                              (1-s) * rhoEnergyFactorSkew
    rhoVelocityFlux_temp[direction] += fpdiag

    -- Return the fluxes in a 5D array
    return {rhoFlux_temp,
            rhoVelocityFlux_temp[0],
            rhoVelocityFlux_temp[1],
            rhoVelocityFlux_temp[2],
            rhoEnergyFlux_temp}
  end
  return CenteredInviscidFlux
end
Flow.CenteredInviscidFluxX = mkCenteredInviscidFlux(0)
Flow.CenteredInviscidFluxY = mkCenteredInviscidFlux(1)
Flow.CenteredInviscidFluxZ = mkCenteredInviscidFlux(2)

-- Compute inviscid and viscous fluxes in all directions.
ebb Flow.AddGetFlux (c : fluidGrid)
  if c.in_interior or c.xneg_depth == 1 then
    ------------
    --- Inviscid
    ------------
    -- Compute the inviscid flux with a centered scheme.
    -- Input the left and right cell states for this face and
    -- the direction index for the flux (x = 0, y = 1, or z = 2).
    var flux = Flow.CenteredInviscidFluxX(c, c(1,0,0))

    -- Store this flux in the cell to the left of the face.
    c.rhoFluxX         =  flux[0]
    c.rhoVelocityFluxX = {flux[1],flux[2],flux[3]}
    c.rhoEnergyFluxX   =  flux[4]

    ------------
    --- Viscous
    ------------
    -- Consider first boundary element (c.xneg_depth == 1) to define left flux
    -- on first interior cell
    var muFace = 0.5 * (GetDynamicViscosity(c.temperature) +
                        GetDynamicViscosity(c(1,0,0).temperature))
    var velocityFace    = L.vec3d({0.0, 0.0, 0.0})
    var velocityX_YFace = L.double(0.0)
    var velocityX_ZFace = L.double(0.0)
    var velocityY_YFace = L.double(0.0)
    var velocityZ_ZFace = L.double(0.0)

    -- Interpolate velocity and derivatives to face
    velocityFace = 0.5 * ( c.velocity + c(1,0,0).velocity )
    velocityX_YFace = 0.5 * ( c.velocityGradientY[0] +
                              c(1,0,0).velocityGradientY[0] )
    velocityX_ZFace = 0.5 * ( c.velocityGradientZ[0] +
                              c(1,0,0).velocityGradientZ[0] )
    velocityY_YFace = 0.5 * ( c.velocityGradientY[1] +
                              c(1,0,0).velocityGradientY[1] )
    velocityZ_ZFace = 0.5 * ( c.velocityGradientZ[2] +
                              c(1,0,0).velocityGradientZ[2] )

    -- Differentiate at face
    var velocityX_XFace   = L.double(0.0)
    var velocityY_XFace   = L.double(0.0)
    var velocityZ_XFace   = L.double(0.0)
    var temperature_XFace = L.double(0.0)

    velocityX_XFace   = 0.5*( c(1,0,0).velocity[0] - c.velocity[0] )
    velocityY_XFace   = 0.5*( c(1,0,0).velocity[1] - c.velocity[1] )
    velocityZ_XFace   = 0.5*( c(1,0,0).velocity[2] - c.velocity[2] )
    temperature_XFace = 0.5*( c(1,0,0).temperature - c.temperature )

    -- Half cell size due to the 0.5 above
    velocityX_XFace   /= (Grid.xCellWidth*0.5)
    velocityY_XFace   /= (Grid.xCellWidth*0.5)
    velocityZ_XFace   /= (Grid.xCellWidth*0.5)
    temperature_XFace /= (Grid.xCellWidth*0.5)

    -- Tensor components (at face)
    var sigmaXX = muFace * ( 4.0 * velocityX_XFace -
                             2.0 * velocityY_YFace -
                             2.0 * velocityZ_ZFace ) / 3.0
    var sigmaYX = muFace * ( velocityY_XFace + velocityX_YFace )
    var sigmaZX = muFace * ( velocityZ_XFace + velocityX_ZFace )
    var usigma  = velocityFace[0] * sigmaXX +
                  velocityFace[1] * sigmaYX +
                  velocityFace[2] * sigmaZX
    var cp = Flow.gamma * Flow.gasConstant / (Flow.gamma - 1.0)
    var heatFlux = - (cp*muFace/Flow.prandtl)*temperature_XFace

    -- Fluxes
    c.rhoVelocityFluxX[0] -= sigmaXX
    c.rhoVelocityFluxX[1] -= sigmaYX
    c.rhoVelocityFluxX[2] -= sigmaZX
    c.rhoEnergyFluxX -= usigma - heatFlux
    -- WARNING: Add SGS terms for LES
  end

  if c.in_interior or c.yneg_depth == 1 then
    ------------
    --- Inviscid
    ------------
    -- Compute the inviscid flux with a centered scheme.
    -- Input the left and right cell states for this face and
    -- the direction index for the flux (x = 0, y = 1, or z = 2).
    var flux = Flow.CenteredInviscidFluxY(c, c(0,1,0))

    -- Store this flux in the cell to the left of the face.
    c.rhoFluxY         =  flux[0]
    c.rhoVelocityFluxY = {flux[1],flux[2],flux[3]}
    c.rhoEnergyFluxY   =  flux[4]

    ------------
    --- Viscous
    ------------
    -- Consider first boundary element (c.yneg_depth == 1) to define down flux
    -- on first interior cell
    var muFace = 0.5 * (GetDynamicViscosity(c.temperature) +
                        GetDynamicViscosity(c(0,1,0).temperature))
    var velocityFace    = L.vec3d({0.0, 0.0, 0.0})
    var velocityY_XFace = L.double(0.0)
    var velocityY_ZFace = L.double(0.0)
    var velocityX_XFace = L.double(0.0)
    var velocityZ_ZFace = L.double(0.0)

    -- Interpolate velocity and derivatives to face
    velocityFace = 0.5 * ( c.velocity + c(0,1,0).velocity )
    velocityY_XFace = 0.5 * ( c.velocityGradientX[1] +
                              c(0,1,0).velocityGradientX[1] )
    velocityY_ZFace = 0.5 * ( c.velocityGradientZ[1] +
                              c(0,1,0).velocityGradientZ[1] )
    velocityX_XFace = 0.5 * ( c.velocityGradientX[0] +
                              c(0,1,0).velocityGradientX[0] )
    velocityZ_ZFace = 0.5 * ( c.velocityGradientZ[2] +
                              c(0,1,0).velocityGradientZ[2] )

    -- Differentiate at face
    var velocityX_YFace   = L.double(0.0)
    var velocityY_YFace   = L.double(0.0)
    var velocityZ_YFace   = L.double(0.0)
    var temperature_YFace = L.double(0.0)

    velocityX_YFace   = 0.5*( c(0,1,0).velocity[0] - c.velocity[0] )
    velocityY_YFace   = 0.5*( c(0,1,0).velocity[1] - c.velocity[1] )
    velocityZ_YFace   = 0.5*( c(0,1,0).velocity[2] - c.velocity[2] )
    temperature_YFace = 0.5*( c(0,1,0).temperature - c.temperature )

    -- Half cell size due to the 0.5 above
    velocityX_YFace   /= (Grid.yCellWidth*0.5)
    velocityY_YFace   /= (Grid.yCellWidth*0.5)
    velocityZ_YFace   /= (Grid.yCellWidth*0.5)
    temperature_YFace /= (Grid.yCellWidth*0.5)

    -- Tensor components (at face)
    var sigmaXY = muFace * ( velocityX_YFace + velocityY_XFace )
    var sigmaYY = muFace * ( 4.0 * velocityY_YFace -
                             2.0 * velocityX_XFace -
                             2.0 * velocityZ_ZFace ) / 3.0
    var sigmaZY = muFace * ( velocityZ_YFace + velocityY_ZFace )
    var usigma  = velocityFace[0] * sigmaXY +
                  velocityFace[1] * sigmaYY +
                  velocityFace[2] * sigmaZY
    var cp = Flow.gamma * Flow.gasConstant / (Flow.gamma - 1.0)
    var heatFlux = - (cp*muFace/Flow.prandtl)*temperature_YFace

    -- Fluxes
    c.rhoVelocityFluxY[0] -= sigmaXY
    c.rhoVelocityFluxY[1] -= sigmaYY
    c.rhoVelocityFluxY[2] -= sigmaZY
    c.rhoEnergyFluxY -= usigma - heatFlux
    -- WARNING: Add SGS terms for LES
  end

  if c.in_interior or c.zneg_depth == 1 then
    ------------
    --- Inviscid
    ------------
    -- Compute the inviscid flux with a centered scheme.
    -- Input the left and right cell states for this face and
    -- the direction index for the flux (x = 0, y = 1, or z = 2).
    var flux = Flow.CenteredInviscidFluxZ(c, c(0,0,1))

    -- Store this flux in the cell to the left of the face.
    c.rhoFluxZ         =  flux[0]
    c.rhoVelocityFluxZ = {flux[1],flux[2],flux[3]}
    c.rhoEnergyFluxZ   =  flux[4]

    ------------
    --- Viscous
    ------------
    -- Consider first boundary element (c.zneg_depth == 1) to define down flux
    -- on first interior cell
    var muFace = 0.5 * (GetDynamicViscosity(c.temperature) +
                        GetDynamicViscosity(c(0,0,1).temperature))
    var velocityFace    = L.vec3d({0.0, 0.0, 0.0})
    var velocityZ_XFace = L.double(0.0)
    var velocityZ_YFace = L.double(0.0)
    var velocityX_XFace = L.double(0.0)
    var velocityY_YFace = L.double(0.0)

    -- Interpolate velocity and derivatives to face
    velocityFace = 0.5 * ( c.velocity + c(0,0,1).velocity )
    velocityZ_XFace = 0.5 * ( c.velocityGradientX[2] +
                              c(0,0,1).velocityGradientX[2] )
    velocityZ_YFace = 0.5 * ( c.velocityGradientY[2] +
                              c(0,0,1).velocityGradientY[2] )
    velocityX_XFace = 0.5 * ( c.velocityGradientX[0] +
                              c(0,0,1).velocityGradientX[0] )
    velocityY_YFace = 0.5 * ( c.velocityGradientY[1] +
                              c(0,0,1).velocityGradientY[1] )

    -- Differentiate at face
    var velocityX_ZFace   = L.double(0.0)
    var velocityY_ZFace   = L.double(0.0)
    var velocityZ_ZFace   = L.double(0.0)
    var temperature_ZFace = L.double(0.0)

    velocityX_ZFace   = 0.5*( c(0,0,1).velocity[0] - c.velocity[0] )
    velocityY_ZFace   = 0.5*( c(0,0,1).velocity[1] - c.velocity[1] )
    velocityZ_ZFace   = 0.5*( c(0,0,1).velocity[2] - c.velocity[2] )
    temperature_ZFace = 0.5*( c(0,0,1).temperature - c.temperature )

    -- Half cell size due to the 0.5 above
    velocityX_ZFace   /= (Grid.zCellWidth*0.5)
    velocityY_ZFace   /= (Grid.zCellWidth*0.5)
    velocityZ_ZFace   /= (Grid.zCellWidth*0.5)
    temperature_ZFace /= (Grid.zCellWidth*0.5)

    -- Tensor components (at face)
    var sigmaXZ = muFace * ( velocityX_ZFace + velocityZ_XFace )
    var sigmaYZ = muFace * ( velocityY_ZFace + velocityZ_YFace )
    var sigmaZZ = muFace * ( 4.0 * velocityZ_ZFace -
                             2.0 * velocityX_XFace -
                             2.0 * velocityY_YFace ) / 3.0
    var usigma  = velocityFace[0] * sigmaXZ +
                  velocityFace[1] * sigmaYZ +
                  velocityFace[2] * sigmaZZ
    var cp = Flow.gamma * Flow.gasConstant / (Flow.gamma - 1.0)
    var heatFlux = - (cp*muFace/Flow.prandtl)*temperature_ZFace

    -- Fluxes
    c.rhoVelocityFluxZ[0] -= sigmaXZ
    c.rhoVelocityFluxZ[1] -= sigmaYZ
    c.rhoVelocityFluxZ[2] -= sigmaZZ
    c.rhoEnergyFluxZ -= usigma - heatFlux
    -- WARNING: Add SGS terms for LES
  end

  -- HACK: obnoxious workaround to establish strict domninance between tasks
  var v = 0
  if v == 1 then
    var tmp1 = c(-1,0,0).velocity[0]
    var tmp2 = c(0,-1,0).velocity[1]
    var tmp3 = c(0,0,-1).velocity[2]
  end
end

ebb Flow.AddUpdateUsingFlux (c : fluidGrid)
  --if c.in_interior or c.xneg_depth == 1 then
  if c.in_interior then --or c.xneg_depth == 1 then
    c.rho_t += -(c( 0,0,0).rhoFluxX -
                 c(-1,0,0).rhoFluxX)/Grid.xCellWidth
    c.rhoVelocity_t += -(c( 0,0,0).rhoVelocityFluxX -
                         c(-1,0,0).rhoVelocityFluxX)/Grid.xCellWidth
    c.rhoEnergy_t += -(c( 0,0,0).rhoEnergyFluxX -
                       c(-1,0,0).rhoEnergyFluxX)/Grid.xCellWidth
  --end
  --if c.in_interior or c.zneg_depth == 1 then
    c.rho_t += -(c(0, 0,0).rhoFluxY -
                 c(0,-1,0).rhoFluxY)/Grid.yCellWidth
    c.rhoVelocity_t += -(c(0, 0,0).rhoVelocityFluxY -
                         c(0,-1,0).rhoVelocityFluxY)/Grid.yCellWidth
    c.rhoEnergy_t += -(c(0, 0,0).rhoEnergyFluxY -
                       c(0,-1,0).rhoEnergyFluxY)/Grid.yCellWidth
  --end
  --if c.in_interior or c.zneg_depth == 1 then
    c.rho_t += -(c(0,0, 0).rhoFluxZ -
                 c(0,0,-1).rhoFluxZ)/Grid.zCellWidth
    c.rhoVelocity_t += -(c(0,0, 0).rhoVelocityFluxZ -
                         c(0,0,-1).rhoVelocityFluxZ)/Grid.zCellWidth
    c.rhoEnergy_t += -(c(0,0, 0).rhoEnergyFluxZ -
                       c(0,0,-1).rhoEnergyFluxZ)/Grid.zCellWidth
  end
end

function Flow.AddFluxes ()
  fluidGrid:foreach(Flow.AddGetFlux)
  fluidGrid:foreach(Flow.AddUpdateUsingFlux)
end

function Flow.ComputeVelocityGradients ()
  fluidGrid:foreach(Flow.ComputeVelocityGradientAll)
end

function Flow.UpdateAuxiliaryVelocityConservedAndGradients ()
  fluidGrid:foreach(Flow.UpdateAuxiliaryVelocity)
  Flow.UpdateGhostConserved()
  Flow.UpdateGhostVelocity()
  Flow.ComputeVelocityGradients()
end

function Flow.UpdateAuxiliary ()
  Flow.UpdateAuxiliaryVelocityConservedAndGradients()
  fluidGrid:foreach(Flow.UpdateAuxiliaryThermodynamics)
  Flow.UpdateGhostThermodynamics()
end


------------
-- PARTICLES
------------

ebb Particles.InitializeDensity (p : particles)
  p.density = Particles.density
end

function Particles.InitializePrimitives ()
  M.IF(M.EQ(Particles.initCase, ParticlesInitCase.Random))
    M.ERROR('Random particle initialization is disabled')
  M.END()
  M.IF(M.EQ(Particles.initCase, ParticlesInitCase.Restart))
    particles:Load(
      {'cell','position','velocity','temperature','diameter'},
      'restart_particles_%d.hdf',
      Integrator.restartIter)
    particles:foreach(Particles.InitializeDensity)
    particles:foreach(Particles.CalculateNumber)
  M.END()
  M.IF(M.EQ(Particles.initCase, ParticlesInitCase.Uniform))
    M.INLINE(PARTICLES_INIT.InitParticlesUniform)
    Particles.number:set((Particles.initNum / Grid.numTiles) * Grid.numTiles)
  M.END()
end

function Particles.UpdateAuxiliary ()
  particles:foreach(Particles.UpdateAuxiliaryStep1)
  particles:foreach(Particles.UpdateAuxiliaryStep2)
end

------------------
-- TIME INTEGRATOR
------------------

function Integrator.SetupTimeStep ()
  fluidGrid:foreach(Flow.InitializeTemporaries)
  particles:foreach(Particles.InitializeTemporaries)
end

function Integrator.ConcludeTimeStep ()
  Particles.Collect()
end

function Integrator.InitializeTimeDerivatives ()
  fluidGrid:foreach(Flow.InitializeTimeDerivatives)
  particles:foreach(Particles.InitializeTimeDerivatives)
end

function Integrator.UpdateAuxiliary ()
  Flow.UpdateAuxiliary()
  Particles.UpdateAuxiliary()
end

function Integrator.UpdateTime ()
  -- HACK
  Integrator.simTime:set(Integrator.time_old:get() +
                         -- stage = 1 => 0.5
                         -- stage = 2 => 0.5
                         -- stage = 3 => 1.0
                         -- stage = 4 => 1.0
                         0.5 * (1 + Integrator.stage:get() / 3) *
                         Integrator.deltaTime:get())
end

function Integrator.InitializeVariables ()

  fluidGrid:foreach(Flow.InitializeCell)

  -- Initialize several grid related entitities
  fluidGrid:foreach(Flow.InitializeCenterCoordinates)

  -- Set initial condition for the flow and all auxiliary flow variables
  Flow.InitializePrimitives()
  fluidGrid:foreach(Flow.UpdateConservedFromPrimitive)
  Flow.UpdateAuxiliary()
  Flow.UpdateGhost()

  -- Initialize the particles (position, velocity, temp, diameter, locate)
  Particles.InitializePrimitives()

end

function Integrator.ComputeDFunctionDt ()

  -- Compute flow convective, viscous, and body force residuals
  Flow.UpdateGhostVelocityGradient()
  Flow.AddFluxes()
  fluidGrid:foreach(Flow.AddBodyForces)

  -- FIXME: turbulent-related tasks should be revised
  M.IF(M.EQ(Flow.turbForcing, OnOrOff.ON))
    Flow.AddTurbulentForcing()
  M.END()

  -- Compute residuals for the particles (locate all particles first)
  Particles.Locate()
  particles:foreach(Particles.AddFlowCoupling)
  particles:foreach(Particles.AddBodyForces)

  if Radiation.TYPE == 'Algebraic' then
    M.INLINE(ALGEBRAIC.AddRadiation)
  elseif Radiation.TYPE == 'DOM' then
    -- Compute radiation field values from particles
    domGrid:foreach(Radiation.ClearAccumulators)
    particles:foreach(Radiation.AccumulateParticleValues)
    domGrid:foreach(Radiation.UpdateFieldValues)
    M.INLINE(DOM.ComputeRadiationField)
    -- Absorb radiation into each particle
    particles:foreach(Particles.AbsorbRadiation)
  end

  -- Compute two-way coupling in momentum and energy
  particles:foreach(Flow.AddParticlesCoupling)

end

function Integrator.UpdateSolution ()
  fluidGrid:foreach(Flow.UpdateVars)
  particles:foreach(Particles.UpdateVars)
end

function Integrator.AdvanceTimeStep ()

  Integrator.SetupTimeStep()
  Integrator.time_old:set(Integrator.simTime:get())

  Integrator.stage:set(1)
  M.WHILE(M.LT(Integrator.stage:get(), 5))
    Integrator.InitializeTimeDerivatives()
    Integrator.ComputeDFunctionDt()
    Integrator.UpdateSolution()
    Integrator.UpdateAuxiliary()
    Integrator.UpdateTime()
    Integrator.stage:set(Integrator.stage:get() + 1)
    -- HACK: Move escaping particle deletion here, to appease the SPMD
    -- transformation. It should be fine to do this multiple times.
    Integrator.ConcludeTimeStep()
  M.END()

  Integrator.timeStep:set(Integrator.timeStep:get() + 1)

end

function Integrator.CalculateDeltaTime ()
  -- Check whether we are imposing a delta time or basing it on the CFL,
  -- i.e. a negative CFL was imposed in the config
  M.IF(M.LT(Integrator.cfl, 0.0))
    -- Impose a fixed time step from the config
    Integrator.deltaTime:set(Integrator.fixedDeltaTime)
  M.ELSE()
    -- Calculate the convective, viscous, and heat spectral radii
    Flow.CalculateSpectralRadii()
    -- Calculate diffusive spectral radius as the maximum between
    -- heat conduction and convective spectral radii
    -- Calculate global spectral radius as the maximum between the convective
    -- and diffusive spectral radii
    -- Delta time using the CFL and max spectral radius for stability
    Integrator.deltaTime:set(
      Integrator.cfl /
        M.MAX(Integrator.maxConvectiveSpectralRadius:get(),
        M.MAX(Integrator.maxViscousSpectralRadius:get(),
              Integrator.maxHeatConductionSpectralRadius:get())))
  M.END()
end

-------------
-- STATISTICS
-------------

function Statistics.ResetSpatialAverages ()
  Flow.averagePressure:set(0.0)
  Flow.averageTemperature:set(0.0)
  Flow.averageKineticEnergy:set(0.0)
  Flow.minTemperature:set(math.huge)
  Flow.maxTemperature:set(-math.huge)
  Flow.averagePD:set(0.0)
  Flow.averageDissipation:set(0.0)
  Particles.averageTemperature:set(0.0)
end

function Statistics.UpdateSpatialAverages ()
  -- Flow
  Flow.averagePressure:set(
    Flow.averagePressure:get() / Grid.areaInterior)
  Flow.averageTemperature:set(
    Flow.averageTemperature:get() / Grid.areaInterior)
  Flow.averageKineticEnergy:set(
    Flow.averageKineticEnergy:get() / Grid.areaInterior)
  -- Particles
  Particles.averageTemperature:set(
    Particles.averageTemperature:get() / Particles.number:get())
end

function Statistics.ComputeSpatialAverages ()
  Statistics.ResetSpatialAverages()
  Flow.IntegrateQuantities()
  particles:foreach(Particles.IntegrateQuantities)
  Statistics.UpdateSpatialAverages()
end

-----
-- IO
-----

function IO.WriteConsoleOutput ()
  M.IF(M.EQ(Integrator.timeStep:get() % IO.consoleFrequency, 0))
    -- Output log headers at a specified frequency
    M.IF(M.EQ(Integrator.timeStep:get() % IO.headerFrequency, 0))
      M.PRINT("\n Current time step: %2.6e s.\n",
              Integrator.deltaTime)
      M.PRINT(" Min Flow Temp: %11.6f K. Max Flow Temp: %11.6f K.\n",
              Flow.minTemperature, Flow.maxTemperature)
      M.PRINT(" Current number of particles: %d.\n", Particles.number)
      M.PRINT("\n")
      M.PRINT("    Iter     Time(s)   Avg Press    Avg Temp      Avg KE  Particle T\n")
    M.END()
    -- Output the current stats to the console for this iteration
    M.PRINT("%8d %11.6f %11.6f %11.6f %11.6f %11.6f\n",
            Integrator.timeStep,
            Integrator.simTime,
            Flow.averagePressure,
            Flow.averageTemperature,
            Flow.averageKineticEnergy,
            Particles.averageTemperature)
  M.END()
end

function IO.WriteFlowRestart ()
  -- Check if it is time to output a flow restart file
  M.IF(M.EQ(Integrator.timeStep:get() % IO.restartEveryTimeSteps, 0))
    -- Write the restart files for density, pressure, and velocity
    fluidGrid:Dump({'rho','pressure','velocity'},
                   'restart_fluid_%d.hdf',
                   Integrator.timeStep:get())
  M.END()
end

function IO.WriteParticleRestart ()
  -- Check if it is time to output a particle restart file
  M.IF(M.EQ(Integrator.timeStep:get() % IO.restartEveryTimeSteps, 0))
    -- Write the restart files for position, velocity, temperature and diameter
    particles:Dump({'cell','position','velocity','temperature','diameter'},
                   'restart_particles_%d.hdf',
                   Integrator.timeStep:get())
  M.END()
end

function IO.WriteOutput ()
  -- Write the console output to the screen
  IO.WriteConsoleOutput()
  -- Write the restart files
  M.IF(M.EQ(IO.wrtRestart, OnOrOff.ON))
    -- Write the flow restart files
    IO.WriteFlowRestart()
    -- Write the particle restart files
    IO.WriteParticleRestart()
  M.END()
end

-----------------------------------------------------------------------------
--[[                            MAIN EXECUTION                           ]]--
-----------------------------------------------------------------------------

-- Initialize all variables

Integrator.InitializeVariables()
Statistics.ComputeSpatialAverages()
if Radiation.TYPE == 'DOM' then
  domGrid:foreach(Radiation.InitializeCell)
  M.INLINE(DOM.InitModule)
end
IO.WriteOutput()

-- Main iteration loop

M.WHILE(M.AND(M.LT(Integrator.simTime:get(), Integrator.finalTime),
              M.LT(Integrator.timeStep:get(), Integrator.maxIter)),
        true)
  Integrator.CalculateDeltaTime()
  Integrator.AdvanceTimeStep()
  if not regentlib.config['flow-spmd'] then
    M.IF(M.EQ(Integrator.timeStep:get() % IO.consoleFrequency, 0))
      Statistics.ComputeSpatialAverages()
      IO.WriteOutput()
    M.END()
  end
M.END()

-- Final stats printing

if regentlib.config['flow-spmd'] then
  Statistics.ComputeSpatialAverages()
end
IO.WriteConsoleOutput()

A.translate(Grid.xTiles, Grid.yTiles, Grid.zTiles)
