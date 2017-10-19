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
Boston, MA 02110-1301 USA.

-----------------------------------------------------------------------------
]]--
-----------------------------------------------------------------------------

import 'ebb'

local A    = require 'admiral'
local Grid = require 'ebb.domains.grid'
local L    = require 'ebblib'
local M    = require 'ebb.src.main'
local PN   = require 'ebb.lib.pathname'

-----------------------------------------------------------------------------
--[[                          MATH IMPORTS                               ]]--
-----------------------------------------------------------------------------

local C = terralib.includecstring [[
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>

double rand_double() {
      double r = (double)rand();
      return r;
}

double rand_unity() {
    double r = (double)rand()/(double)RAND_MAX;
    return r;
}

double rand_gauss() {
	static double V1, V2, S;
	static int phase = 0;
	double X;

	if(phase == 0) {
		do {
			double U1 = (double)rand() / (double)RAND_MAX;
			double U2 = (double)rand() / (double)RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
			} while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;

	return X;
}
]]

-- Use the built in rand() function from Liszt
local rand_float = L.rand

-----------------------------------------------------------------------------
--[[                       COMMAND LINE OPTIONS                          ]]--
-----------------------------------------------------------------------------

local function printUsageAndExit()
  print("Usage : liszt-legion.sh ~/path/to/soleil-x.t <options>")
  print("          -i <parameter file with Soleil-X options> (** required **)")
  os.exit(1)
end

local configFileName = nil

do
  local i = 1
  while i <= #arg do
    if arg[i] == '-i' then
      if i == #arg then printUsageAndExit() end
      configFileName = arg[i+1]
      i = i + 2
    else
      i = i + 1
    end
  end
end
if not configFileName then
  print("Config file name required")
  printUsageAndExit()
end

-- Load up the configuration file.

local config,errorMsg = loadfile(configFileName)
if not config then
  print('Error when reading configuration file:')
  print(errorMsg)
  os.exit(1)
end
local config = config()

-- Set the output directory to the current working directory

local outputdir = PN.pwd_str()

-----------------------------------------------------------------------------
--[[                            CONSTANTS                                ]]--
-----------------------------------------------------------------------------

local pi = 2.0*L.acos(0)
local twoPi = 2.0*pi
local SB = 5.67e-8

-----------------------------------------------------------------------------
--[[                            NAMESPACES                               ]]--
-----------------------------------------------------------------------------

local Flow = {}
local Particles = {}
local Radiation = {}
local TimeIntegrator = {}
local Statistics = {}
local IO = {}

-----------------------------------------------------------------------------
--[[                   INITIALIZE OPTIONS FROM CONFIG                    ]]--
-----------------------------------------------------------------------------

local function joinList(list, sep)
  sep = sep or ' '
  local res = ''
  for i,elem in ipairs(list) do
    if i > 1 then
      res = res..sep
    end
    res = res..tostring(elem)
  end
  return res
end

local function Enum(base, ...)
  local enum = {}
  enum.__values = {...}
  for i,val in ipairs({...}) do
    enum[val] = base + i
  end
  return enum
end

local function parseEnum(name, enum)
  local val = enum[config[name]]
  if not val then
    valuesStr = joinList(enum.__values, ', ')
    error('Configuration value "'..name..'" not defined ('..valuesStr..')')
  end
  return val
end

local function parseBool(name)
  if config[name] == 'ON' then
    return true
  elseif config[name] == 'OFF' then
    return false
  end
  error('Configuration value "'..name..'" not defined (ON,OFF)')
end

local FlowBC = Enum(1000, 'periodic','symmetry','adiabatic_wall','isothermal_wall')
local ParticleBC = Enum(2000, 'Permeable','Solid')

local grid_options = {
  -- Number of cells in the x, y, & z directions
  xnum        = config.xnum,
  ynum        = config.ynum,
  znum        = config.znum,
  -- Origin of the computational domain (meters)
  origin      = config.origin,
  -- Width of the computational domain in the x, y, & z directions (meters)
  xWidth      = config.xWidth,
  yWidth      = config.yWidth,
  zWidth      = config.zWidth,
  -- Width of each cell in the x, y, & z directions (meters)
  xCellWidth  = config.xWidth / config.xnum,
  yCellWidth  = config.yWidth / config.ynum,
  zCellWidth  = config.zWidth / config.znum,
  -- Boundary condition type for each face of the block and possible
  -- wall velocity, if no-slip.
  xBCLeft      = parseEnum('xBCLeft', FlowBC),
  xBCLeftVel   = config.xBCLeftVel,
  xBCLeftTemp  = config.xBCLeftTemp,
  xBCRight     = parseEnum('xBCRight', FlowBC),
  xBCRightVel  = config.xBCRightVel,
  xBCRightTemp = config.xBCRightTemp,
  yBCLeft      = parseEnum('yBCLeft', FlowBC),
  yBCLeftVel   = config.yBCLeftVel,
  yBCLeftTemp  = config.yBCLeftTemp,
  yBCRight     = parseEnum('yBCRight', FlowBC),
  yBCRightVel  = config.yBCRightVel,
  yBCRightTemp = config.yBCRightTemp,
  zBCLeft      = parseEnum('zBCLeft', FlowBC),
  zBCLeftVel   = config.zBCLeftVel,
  zBCLeftTemp  = config.zBCLeftTemp,
  zBCRight     = parseEnum('zBCRight', FlowBC),
  zBCRightVel  = config.zBCRightVel,
  zBCRightTemp = config.zBCRightTemp,
}

-- Define offsets for boundary conditions in flow solver
-- The sign variables define the necessary reflections for the
-- different types of BCs. The wall velocity is specified above,
-- and then the velocity adjustment is calculated here and applied
-- to the boundaries below.

-- Define offsets, signs, and velocities for the x BCs

local x_sign
local xpos_velocity
local xneg_velocity
local xpos_temperature
local xneg_temperature

if grid_options.xBCLeft  == FlowBC.periodic and
   grid_options.xBCRight == FlowBC.periodic then
  x_sign = L.Constant(L.vec3d, {1.0,1.0,1.0})
  xpos_velocity = L.Constant(L.vec3d, {0.0,0.0,0.0})
  xneg_velocity = L.Constant(L.vec3d, {0.0,0.0,0.0})
  xpos_temperature = L.Constant(L.double, -1.0)
  xneg_temperature = L.Constant(L.double, -1.0)
  grid_options.xBCLeftParticles  = ParticleBC.Permeable
  grid_options.xBCRightParticles = ParticleBC.Permeable
elseif grid_options.xBCLeft == FlowBC.symmetry and
       grid_options.xBCRight == FlowBC.symmetry then
  x_sign = L.Constant(L.vec3d, {-1.0,1.0,1.0})
  xpos_velocity = L.Constant(L.vec3d, {0.0,0.0,0.0})
  xneg_velocity = L.Constant(L.vec3d, {0.0,0.0,0.0})
  xpos_temperature = L.Constant(L.double, -1.0)
  xneg_temperature = L.Constant(L.double, -1.0)
  grid_options.xBCLeftParticles  = ParticleBC.Solid
  grid_options.xBCRightParticles = ParticleBC.Solid
elseif grid_options.xBCLeft  == FlowBC.adiabatic_wall and
       grid_options.xBCRight == FlowBC.adiabatic_wall then
  x_sign = L.Constant(L.vec3d, {-1.0,-1.0,-1.0})
  xpos_velocity = L.Constant(L.vec3d, {2.0*grid_options.xBCRightVel[1],
                                       2.0*grid_options.xBCRightVel[2],
                                       2.0*grid_options.xBCRightVel[3]})
  xneg_velocity = L.Constant(L.vec3d, {2.0*grid_options.xBCLeftVel[1],
                                       2.0*grid_options.xBCLeftVel[2],
                                       2.0*grid_options.xBCLeftVel[3]})
  xpos_temperature = L.Constant(L.double, -1.0)
  xneg_temperature = L.Constant(L.double, -1.0)
  grid_options.xBCLeftParticles  = ParticleBC.Solid
  grid_options.xBCRightParticles = ParticleBC.Solid
elseif grid_options.xBCLeft  == FlowBC.isothermal_wall and
       grid_options.xBCRight == FlowBC.isothermal_wall then
  x_sign = L.Constant(L.vec3d, {-1.0,-1.0,-1.0})
  xpos_velocity = L.Constant(L.vec3d, {2.0*grid_options.xBCRightVel[1],
                                       2.0*grid_options.xBCRightVel[2],
                                       2.0*grid_options.xBCRightVel[3]})
  xneg_velocity = L.Constant(L.vec3d, {2.0*grid_options.xBCLeftVel[1],
                                       2.0*grid_options.xBCLeftVel[2],
                                       2.0*grid_options.xBCLeftVel[3]})
  xpos_temperature = L.Constant(L.double, grid_options.xBCRightTemp)
  xneg_temperature = L.Constant(L.double, grid_options.xBCLeftTemp)
  grid_options.xBCLeftParticles  = ParticleBC.Solid
  grid_options.xBCRightParticles = ParticleBC.Solid
else
  error("Boundary conditions in x not implemented")
end

-- Define offsets, signs, and velocities for the y BCs

local y_sign
local ypos_velocity
local yneg_velocity
local ypos_temperature
local yneg_temperature

if grid_options.yBCLeft  == FlowBC.periodic and
   grid_options.yBCRight == FlowBC.periodic then
  y_sign = L.Constant(L.vec3d, {1.0,1.0,1.0})
  ypos_velocity = L.Constant(L.vec3d, {0.0,0.0,0.0})
  yneg_velocity = L.Constant(L.vec3d, {0.0,0.0,0.0})
  ypos_temperature = L.Constant(L.double, -1.0)
  yneg_temperature = L.Constant(L.double, -1.0)
  grid_options.yBCLeftParticles  = ParticleBC.Permeable
  grid_options.yBCRightParticles = ParticleBC.Permeable
elseif grid_options.yBCLeft  == FlowBC.symmetry and
       grid_options.yBCRight == FlowBC.symmetry then
  y_sign = L.Constant(L.vec3d, {1.0,-1.0,1.0})
  ypos_velocity = L.Constant(L.vec3d, {0.0,0.0,0.0})
  yneg_velocity = L.Constant(L.vec3d, {0.0,0.0,0.0})
  ypos_temperature = L.Constant(L.double, -1.0)
  yneg_temperature = L.Constant(L.double, -1.0)
  grid_options.yBCLeftParticles  = ParticleBC.Solid
  grid_options.yBCRightParticles = ParticleBC.Solid
elseif grid_options.yBCLeft  == FlowBC.adiabatic_wall and
       grid_options.yBCRight == FlowBC.adiabatic_wall then
  y_sign = L.Constant(L.vec3d, {-1.0,-1.0,-1.0})
  ypos_velocity = L.Constant(L.vec3d, {2.0*grid_options.yBCRightVel[1],
                                       2.0*grid_options.yBCRightVel[2],
                                       2.0*grid_options.yBCRightVel[3]})
  yneg_velocity = L.Constant(L.vec3d, {2.0*grid_options.yBCLeftVel[1],
                                       2.0*grid_options.yBCLeftVel[2],
                                       2.0*grid_options.yBCLeftVel[3]})
  ypos_temperature = L.Constant(L.double, -1.0)
  yneg_temperature = L.Constant(L.double, -1.0)
  grid_options.yBCLeftParticles  = ParticleBC.Solid
  grid_options.yBCRightParticles = ParticleBC.Solid
elseif grid_options.yBCLeft  == FlowBC.isothermal_wall and
       grid_options.yBCRight == FlowBC.isothermal_wall then
  y_sign = L.Constant(L.vec3d, {-1.0,-1.0,-1.0})
  ypos_velocity = L.Constant(L.vec3d, {2.0*grid_options.yBCRightVel[1],
                                       2.0*grid_options.yBCRightVel[2],
                                       2.0*grid_options.yBCRightVel[3]})
  yneg_velocity = L.Constant(L.vec3d, {2.0*grid_options.yBCLeftVel[1],
                                       2.0*grid_options.yBCLeftVel[2],
                                       2.0*grid_options.yBCLeftVel[3]})
  ypos_temperature = L.Constant(L.double, grid_options.yBCRightTemp)
  yneg_temperature = L.Constant(L.double, grid_options.yBCLeftTemp)
  grid_options.yBCLeftParticles  = ParticleBC.Solid
  grid_options.yBCRightParticles = ParticleBC.Solid
else
  error("Boundary conditions in y not implemented")
end

-- Define offsets, signs, and velocities for the z BCs

local z_sign
local zpos_velocity
local zneg_velocity
local zpos_temperature
local zneg_temperature

if grid_options.zBCLeft  == FlowBC.periodic and
   grid_options.zBCRight == FlowBC.periodic then
  z_sign = L.Constant(L.vec3d, {1.0,1.0,1.0})
  zpos_velocity = L.Constant(L.vec3d, {0.0,0.0,0.0})
  zneg_velocity = L.Constant(L.vec3d, {0.0,0.0,0.0})
  zpos_temperature = L.Constant(L.double, -1.0)
  zneg_temperature = L.Constant(L.double, -1.0)
  grid_options.zBCLeftParticles  = ParticleBC.Permeable
  grid_options.zBCRightParticles = ParticleBC.Permeable
elseif grid_options.zBCLeft == FlowBC.symmetry and
       grid_options.zBCRight == FlowBC.symmetry then
  z_sign = L.Constant(L.vec3d, {1.0,1.0,-1.0})
  zpos_velocity = L.Constant(L.vec3d, {0.0,0.0,0.0})
  zneg_velocity = L.Constant(L.vec3d, {0.0,0.0,0.0})
  zpos_temperature = L.Constant(L.double, -1.0)
  zneg_temperature = L.Constant(L.double, -1.0)
  grid_options.zBCLeftParticles  = ParticleBC.Solid
  grid_options.zBCRightParticles = ParticleBC.Solid
elseif grid_options.zBCLeft  == FlowBC.adiabatic_wall and
       grid_options.zBCRight == FlowBC.adiabatic_wall then
  z_sign = L.Constant(L.vec3d, {-1.0,-1.0,-1.0})
  zpos_velocity = L.Constant(L.vec3d, {2.0*grid_options.zBCRightVel[1],
                                       2.0*grid_options.zBCRightVel[2],
                                       2.0*grid_options.zBCRightVel[3]})
  zneg_velocity = L.Constant(L.vec3d, {2.0*grid_options.zBCLeftVel[1],
                                       2.0*grid_options.zBCLeftVel[2],
                                       2.0*grid_options.zBCLeftVel[3]})
  zpos_temperature = L.Constant(L.double, -1.0)
  zneg_temperature = L.Constant(L.double, -1.0)
  grid_options.zBCLeftParticles  = ParticleBC.Solid
  grid_options.zBCRightParticles = ParticleBC.Solid
elseif grid_options.zBCLeft  == FlowBC.isothermal_wall and
       grid_options.zBCRight == FlowBC.isothermal_wall then
  z_sign = L.Constant(L.vec3d, {-1.0,-1.0,-1.0})
  zpos_velocity = L.Constant(L.vec3d, {2.0*grid_options.zBCRightVel[1],
                                       2.0*grid_options.zBCRightVel[2],
                                       2.0*grid_options.zBCRightVel[3]})
  zneg_velocity = L.Constant(L.vec3d, {2.0*grid_options.zBCLeftVel[1],
                                       2.0*grid_options.zBCLeftVel[2],
                                       2.0*grid_options.zBCLeftVel[3]})
  zpos_temperature = L.Constant(L.double, grid_options.zBCRightTemp)
  zneg_temperature = L.Constant(L.double, grid_options.zBCLeftTemp)
  grid_options.zBCLeftParticles  = ParticleBC.Solid
  grid_options.zBCRightParticles = ParticleBC.Solid
else
  error("Boundary conditions in z not implemented")
end

-- Spatial integration options
local spatial_options = {
  split = 0.5  --  Splitting parameter
}

-- Time integration options
local time_options = {
  final_time            = config.final_time,
  restartIter           = config.restartIter,
  max_iter              = config.max_iter,
  cfl                   = config.cfl,
  delta_time            = config.delta_time,
  restartEveryTimeSteps = config.restartEveryTimeSteps,
  headerFrequency       = config.headerFrequency,
  consoleFrequency      = config.consoleFrequency,
}

local ViscosityModel = Enum(3000, 'Constant','PowerLaw','Sutherland')
local fluid_options = {
  viscosity_model    = parseEnum('viscosity_model', ViscosityModel),
  gasConstant        = config.gasConstant,
  gamma              = config.gamma,
  prandtl            = config.prandtl,
  constant_visc      = config.constant_visc,
  powerlaw_visc_ref  = config.powerlaw_visc_ref,
  powerlaw_temp_ref  = config.powerlaw_temp_ref,
  suth_visc_ref      = config.suth_visc_ref,
  suth_temp_ref      = config.suth_temp_ref,
  suth_s_ref         = config.suth_s_ref,
}

local InitCase = Enum(4000, 'Uniform','Restart','Perturbed','TaylorGreen2DVortex','TaylorGreen3DVortex')
local flow_options = {
  initCase       = parseEnum('initCase', InitCase),
  initParams     = L.Constant(L.vector(L.double,5), config.initParams),
  bodyForce      = L.Constant(L.vec3d, config.bodyForce),
  turbForceCoeff = L.Constant(L.double, config.turbForceCoeff),
  turbForcing    = parseBool('turbForcing'),
}

local InitParticles = Enum(5000, 'Random','Restart','Uniform')
local ParticleType = Enum(6000, 'Fixed','Free')
local particles_options = {
  -- Define the initial number of particles and insertion/deletion
  num            = config.num,
  maximum_num    = config.maximum_num,

  -- Particle characteristics
  restitution_coefficient = L.Constant(L.double, config.restitutionCoefficient),
  convective_coefficient  = L.Constant(L.double, config.convectiveCoefficient),
  heat_capacity           = L.Constant(L.double, config.heatCapacity),
  initialTemperature      = config.initialTemperature,
  density                 = config.density,
  diameter_mean           = config.diameter_mean,
  diameter_maxDeviation   = config.diameter_maxDeviation,
  bodyForce               = L.Constant(L.vec3d, config.bodyForceParticles),
  absorptivity            = config.absorptivity,
  restartParticleIter     = config.restartParticleIter,

  -- Particles mode
  modeParticles  = parseBool('modeParticles'),
  initParticles  = parseEnum('initParticles', InitParticles),
  particleType   = parseEnum('particleType', ParticleType),
  twoWayCoupling = parseBool('twoWayCoupling'),
}

local RadiationType = Enum(7000, 'Algebraic','DOM','MCRT','OFF')
local radiation_options = {
  radiationType      = parseEnum('radiationType', RadiationType),
  zeroAvgHeatSource  = parseBool('zeroAvgHeatSource'),
}
if radiation_options.radiationType ~= RadiationType.OFF
  and not particles_options.modeParticles then
  error('Radiation support requires particles to be enabled')
end
if radiation_options.radiationType == RadiationType.DOM then
  radiation_options.qa = config.qa
  radiation_options.qs = config.qs
  radiation_options.numAngles = config.numAngles
  radiation_options.coarsenFactor = config.coarsenFactor
  radiation_options.xCellWidth = grid_options.xCellWidth * config.coarsenFactor[1]
  radiation_options.yCellWidth = grid_options.yCellWidth * config.coarsenFactor[2]
  radiation_options.zCellWidth = grid_options.zCellWidth * config.coarsenFactor[3]
  radiation_options.cellVolume = radiation_options.xCellWidth *
                                 radiation_options.yCellWidth *
                                 radiation_options.zCellWidth
end

local io_options = {
  wrtRestart           = parseBool('wrtRestart'),
  outputFileNamePrefix = outputdir .. '/',
}

-----------------------------------------------------------------------------
--[[                      BOUNDARY CONFIG CHECK                          ]]--
-----------------------------------------------------------------------------

-- Check boundary type consistency for the periodic BCs
if ( grid_options.xBCLeft  == FlowBC.periodic and
     grid_options.xBCRight ~= FlowBC.periodic ) or
   ( grid_options.xBCLeft  ~= FlowBC.periodic and
     grid_options.xBCRight == FlowBC.periodic ) then
  error("Boundary conditions in x should match for periodicity")
end
if ( grid_options.yBCLeft  == FlowBC.periodic and
     grid_options.yBCRight ~= FlowBC.periodic ) or
   ( grid_options.yBCLeft  ~= FlowBC.periodic and
     grid_options.yBCRight == FlowBC.periodic ) then
  error("Boundary conditions in y should match for periodicity")
end
if ( grid_options.zBCLeft  == FlowBC.periodic and
     grid_options.zBCRight ~= FlowBC.periodic ) or
   ( grid_options.zBCLeft  ~= FlowBC.periodic and
     grid_options.zBCRight == FlowBC.periodic ) then
  error("Boundary conditions in z should match for periodicity")
end
if ( grid_options.xBCLeft  == FlowBC.periodic and
     grid_options.xBCRight == FlowBC.periodic ) then
  xBCPeriodic = true
else
  xBCPeriodic = false
end
if ( grid_options.yBCLeft  == FlowBC.periodic and
     grid_options.yBCRight == FlowBC.periodic ) then
  yBCPeriodic = true
else
  yBCPeriodic = false
end
if ( grid_options.zBCLeft  == FlowBC.periodic and
     grid_options.zBCRight == FlowBC.periodic ) then
  zBCPeriodic = true
else
  zBCPeriodic = false
end

-----------------------------------------------------------------------------
--[[                         GRID PREPROCESSING                          ]]--
-----------------------------------------------------------------------------

-- Declare and initialize grid and related fields
-- As we are second-order, we will initialize the grid
-- with a single layer of halo cells (unless running a
-- periodic case, which is natively handled w/out halos).
local bnum = 1
if xBCPeriodic then xBnum = 0 else xBnum = bnum end
if yBCPeriodic then yBnum = 0 else yBnum = bnum end
if zBCPeriodic then zBnum = 0 else zBnum = bnum end
local xBw = grid_options.xWidth/grid_options.xnum * xBnum
local yBw = grid_options.yWidth/grid_options.ynum * yBnum
local zBw = grid_options.zWidth/grid_options.znum * zBnum
local gridOriginInteriorX = grid_options.origin[1]
local gridOriginInteriorY = grid_options.origin[2]
local gridOriginInteriorZ = grid_options.origin[3]

local fluidGrid = Grid.NewGrid{
  name              = 'Fluid',
  dims              = {grid_options.xnum + 2*xBnum,
                       grid_options.ynum + 2*yBnum,
                       grid_options.znum + 2*zBnum},
  origin            = {grid_options.origin[1] - xBw,
                       grid_options.origin[2] - yBw,
                       grid_options.origin[3] - zBw},
  width             = {grid_options.xWidth + 2*xBw,
                       grid_options.yWidth + 2*yBw,
                       grid_options.zWidth + 2*zBw},
  boundary_depth    = {xBnum, yBnum, zBnum},
  periodic          = {xBCPeriodic, yBCPeriodic, zBCPeriodic}
}

-- Define uniform grid spacing
-- WARNING: These are used for uniform grids and should be replaced by different
-- metrics for non-uniform ones (see other WARNINGS throughout the code)
local grid_originX = L.Constant(L.double, fluidGrid:Origin()[1])
local grid_originY = L.Constant(L.double, fluidGrid:Origin()[2])
local grid_originZ = L.Constant(L.double, fluidGrid:Origin()[3])
local grid_widthX  = L.Constant(L.double, fluidGrid:Width()[1])
local grid_widthY  = L.Constant(L.double, fluidGrid:Width()[2])
local grid_widthZ  = L.Constant(L.double, fluidGrid:Width()[3])
local grid_dx      = L.Constant(L.double, fluidGrid:CellWidth()[1])
local grid_dy      = L.Constant(L.double, fluidGrid:CellWidth()[2])
local grid_dz      = L.Constant(L.double, fluidGrid:CellWidth()[3])

-- Primitive variables
fluidGrid:NewField('rho', L.double)
fluidGrid:NewField('pressure', L.double)
fluidGrid:NewField('velocity', L.vec3d)

-- Remaining primitive variables
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
--[[                       PARTICLE PREPROCESSING                        ]]--
-----------------------------------------------------------------------------

-- Check whether particles are even active in order to avoid allocating
-- any data for the particles.

local particles

if particles_options.modeParticles then

  -- Declare particle relation and fields over the particles.
  -- This is a flexible relation, and thus starts as empty, so we don't
  -- initialize the fields here.

  particles = L.NewRelation {
    name = 'particles',
    mode = 'COUPLED',
    coupled_with = fluidGrid,
    coupling_field = 'cell',
    size = particles_options.maximum_num,
    max_skew = 1.5,
    max_xfer_num = 1000,
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

  particles:NewField('position', L.vec3d)
  particles:NewField('particle_velocity', L.vec3d)
  particles:NewField('density', L.double)
  particles:NewField('particle_temperature', L.double)
  particles:NewField('diameter', L.double)
  particles:NewField('position_ghost', L.vec3d)
  particles:NewField('velocity_ghost', L.vec3d)
  particles:NewField('velocity_t_ghost', L.vec3d)
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
  -- derivatives
  particles:NewField('position_t', L.vec3d)
  particles:NewField('velocity_t', L.vec3d)
  particles:NewField('temperature_t', L.double)

end

-----------------------------------------------------------------------------
--[[                        RADIATION PREPROCESSING                      ]]--
-----------------------------------------------------------------------------

local radiationGrid

if radiation_options.radiationType == RadiationType.DOM then

  radiationGrid = fluidGrid:Coarsen('Radiation', radiation_options.coarsenFactor, {0,0,0})

  -- cell center intensity per angle
  radiationGrid:NewField('I_1', L.vector(L.double,radiation_options.numAngles))
  radiationGrid:NewField('I_2', L.vector(L.double,radiation_options.numAngles))
  radiationGrid:NewField('I_3', L.vector(L.double,radiation_options.numAngles))
  radiationGrid:NewField('I_4', L.vector(L.double,radiation_options.numAngles))
  radiationGrid:NewField('I_5', L.vector(L.double,radiation_options.numAngles))
  radiationGrid:NewField('I_6', L.vector(L.double,radiation_options.numAngles))
  radiationGrid:NewField('I_7', L.vector(L.double,radiation_options.numAngles))
  radiationGrid:NewField('I_8', L.vector(L.double,radiation_options.numAngles))

  -- iterative intensity per angle
  radiationGrid:NewField('Iiter_1', L.vector(L.double,radiation_options.numAngles))
  radiationGrid:NewField('Iiter_2', L.vector(L.double,radiation_options.numAngles))
  radiationGrid:NewField('Iiter_3', L.vector(L.double,radiation_options.numAngles))
  radiationGrid:NewField('Iiter_4', L.vector(L.double,radiation_options.numAngles))
  radiationGrid:NewField('Iiter_5', L.vector(L.double,radiation_options.numAngles))
  radiationGrid:NewField('Iiter_6', L.vector(L.double,radiation_options.numAngles))
  radiationGrid:NewField('Iiter_7', L.vector(L.double,radiation_options.numAngles))
  radiationGrid:NewField('Iiter_8', L.vector(L.double,radiation_options.numAngles))

  radiationGrid:NewField('G',     L.double) -- intensity summation over all angles
  radiationGrid:NewField('S',     L.double) -- source term
  radiationGrid:NewField('Ib',    L.double) -- blackbody intensity
  radiationGrid:NewField('sigma', L.double) -- extinction coefficient

  -- partial sums, over particles inside volume
  radiationGrid:NewField('acc_d2',   L.double) -- p.diameter^2
  radiationGrid:NewField('acc_d2t4', L.double) -- p.diameter^2 * p.temperature^4

end

-----------------------------------------------------------------------------
--[[                           GLOBAL VARIABLES                          ]]--
-----------------------------------------------------------------------------

-- Integration quantities

TimeIntegrator.simTime   = L.Global('TimeIntegrator.simTime', L.double, 0)
TimeIntegrator.timeOld   = L.Global('TimeIntegrator.timeOld', L.double, 0)
TimeIntegrator.timeStep  = L.Global('TimeIntegrator.timeStep', L.int, 0)
TimeIntegrator.deltaTime = L.Global('TimeIntegrator.deltaTime', L.double, 0.0001)
TimeIntegrator.stage     = L.Global('TimeIntegrator.stage', L.int, 0)

-- Statistics quantities

-- Note: - numberOfInteriorCells and areaInterior could be defined as variables
-- from grid instead of Flow. Here Flow is used to avoid adding things to grid
-- externally
Flow.numberOfInteriorCells   = L.Global('Flow.numberOfInteriorCells', L.int64, 0)
Flow.areaInterior            = L.Global('Flow.areaInterior', L.double, 0.0)
Flow.averagePressure         = L.Global('Flow.averagePressure', L.double, 0.0)
Flow.averageTemperature      = L.Global('Flow.averageTemperature', L.double, 0.0)
Flow.averageHeatSource       = L.Global('Flow.averageHeatSource', L.double, 0.0)
Flow.averageKineticEnergy    = L.Global('Flow.averageKineticEnergy', L.double, 0.0)
Flow.minTemperature          = L.Global('Flow.minTemperature', L.double, 0.0)
Flow.maxTemperature          = L.Global('Flow.maxTemperature', L.double, 0.0)
Particles.averageTemperature = L.Global('Particles.averageTemperature', L.double, 0.0)
Particles.number             = L.Global('Particles.number', L.int64, 0)
Particles.limit              = L.Global('Particles.limit', L.int64, 0)

Flow.averagePD          = L.Global('Flow.averagePD', L.double, 0.0)
Flow.averageDissipation = L.Global('Flow.averageDissipation', L.double, 0.0)
Flow.averageFe          = L.Global('Flow.averageFe', L.double, 0.0)
Flow.averageK           = L.Global('Flow.averageK', L.double, 0.0)

-----------------------------------------------------------------------------
--[[                       EXTERNAL REGENT MODULES                       ]]--
-----------------------------------------------------------------------------

local particles_init_uniform
if particles_options.modeParticles then
  particles_init_uniform =
    (require 'particles_init_uniform')(particles, fluidGrid)
end

local radiation
if radiation_options.radiationType == RadiationType.Algebraic then
  radiation = (require 'algebraic')(particles)
elseif radiation_options.radiationType == RadiationType.DOM then
  radiation = (require 'dom/dom')(radiationGrid)
elseif radiation_options.radiationType == RadiationType.MCRT then
  error('MCRT not supported yet')
elseif radiation_options.radiationType == RadiationType.OFF then
  -- do nothing
else assert(false) end

-----------------------------------------------------------------------------
--[[                       LOAD DATA FOR RESTART                         ]]--
-----------------------------------------------------------------------------

if flow_options.initCase == InitCase.Restart then
  -- Increment the time step and physical time so the simulation doesn't
  -- repeat from 0. Also, increase the max number of iterations so the solve
  -- doesn't immediately exit.
  TimeIntegrator.timeStep:set(time_options.restartIter)
  -- TODO: No way to pass TimeIntegrator.simTime for the restart
  time_options.max_iter = time_options.max_iter + time_options.restartIter
end

-----------------------------------------------------------------------------
--[[                       USER-DEFINED FUNCTIONS                        ]]--
-----------------------------------------------------------------------------

-- Norm of a vector
local ebb norm (v)
  return L.sqrt(L.dot(v, v))
end

-- Compute fluid dynamic viscosity from fluid temperature
local ebb GetDynamicViscosity (temperature)
  var viscosity = L.double(0.0)
  if fluid_options.viscosity_model == ViscosityModel.Constant then
    viscosity = fluid_options.constant_visc
  elseif fluid_options.viscosity_model == ViscosityModel.PowerLaw then
    viscosity = fluid_options.powerlaw_visc_ref *
      L.pow(temperature/fluid_options.powerlaw_temp_ref, 0.75)
  elseif fluid_options.viscosity_model == ViscosityModel.Sutherland then
    viscosity = fluid_options.suth_visc_ref *
      L.pow((temperature/fluid_options.suth_temp_ref),(3.0/2.0))*
      ((fluid_options.suth_temp_ref + fluid_options.suth_s_ref)/
         (temperature + fluid_options.suth_s_ref))
  else L.assert(false) end
  return viscosity
end

-- Compute fluid flow sound speed based on temperature (a = sqrt(gamma*R*T))
local ebb GetSoundSpeed (temperature)
  return L.sqrt(fluid_options.gamma * fluid_options.gasConstant * temperature)
end

-- Function to retrieve particle area, volume and mass
-- These are Ebb user-defined functions that behave like a field
if particles_options.modeParticles then
  particles:NewFieldReadFunction('cross_section_area', ebb(p)
    return pi * L.pow(p.diameter, 2) / 4.0
  end)
  particles:NewFieldReadFunction('volume', ebb(p)
    return pi * L.pow(p.diameter, 3) / 6.0
  end)
  particles:NewFieldReadFunction('mass', ebb(p)
    return p.volume * p.density
  end)
end

-- Function for returning a Gaussian random variable
local ebb rand_gauss()
  var x = L.double(0.0)
  for i = 0,25 do
    x += rand_float()
  end
  x -= 25.0 / 2.0
  x /= L.sqrt(25.0 / 12.0)
  return x
end

-- WARNING: update cellVolume computation for non-uniform grids
local cellVolume = L.Constant(L.double,
                              grid_dx:get() * grid_dy:get() * grid_dz:get())
local ebb numberOfInteriorCells ( c : fluidGrid )
  if c.in_interior then
    Flow.numberOfInteriorCells += L.int64(1)
  end
end
local ebb areaInterior ( c : fluidGrid )
  if c.in_interior then
    Flow.areaInterior += cellVolume
  end
end
function Flow.IntegrateGeometricQuantities(cells)
  Flow.numberOfInteriorCells:set(0)
  Flow.areaInterior:set(0)
  cells:foreach(numberOfInteriorCells)
  cells:foreach(areaInterior         )
end

-----------------------------------------------------------------------------
--[[                              EBB MACROS                             ]]--
-----------------------------------------------------------------------------

ebb TrilinearInterpolateRho (xyz, c000, c100, c010, c110, c001, c101, c011, c111 )
  var dX   = L.fmod((xyz[0] - grid_originX)/grid_dx + 0.5, 1.0)
  var dY   = L.fmod((xyz[1] - grid_originY)/grid_dy + 0.5, 1.0)
  var dZ   = L.fmod((xyz[2] - grid_originZ)/grid_dz + 0.5, 1.0)

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

ebb TrilinearInterpolateVelocity (xyz, c000, c100, c010, c110, c001, c101, c011, c111 )
  var dX   = L.fmod((xyz[0] - grid_originX)/grid_dx + 0.5, 1.0)
  var dY   = L.fmod((xyz[1] - grid_originY)/grid_dy + 0.5, 1.0)
  var dZ   = L.fmod((xyz[2] - grid_originZ)/grid_dz + 0.5, 1.0)

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

ebb TrilinearInterpolateTemp (xyz, c000, c100, c010, c110, c001, c101, c011, c111 )
  var dX   = L.fmod((xyz[0] - grid_originX)/grid_dx + 0.5, 1.0)
  var dY   = L.fmod((xyz[1] - grid_originY)/grid_dy + 0.5, 1.0)
  var dZ   = L.fmod((xyz[2] - grid_originZ)/grid_dz + 0.5, 1.0)

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

ebb InterpolateTriVelocity (c, xyz)
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
                                            velocity011, velocity111 )
end

ebb InterpolateTriTemp (c, xyz)
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
                                        temp011, temp111 )
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
  c.rho         = flow_options.initParams[0]
  c.pressure    = flow_options.initParams[1]
  c.velocity[0] = flow_options.initParams[2]
  c.velocity[1] = flow_options.initParams[3]
  c.velocity[2] = flow_options.initParams[4]
end

ebb Flow.InitializeTaylorGreen2D (c : fluidGrid)
  -- Define Taylor Green Vortex
  var taylorGreenDensity  = flow_options.initParams[0]
  var taylorGreenPressure = flow_options.initParams[1]
  var taylorGreenVelocity = flow_options.initParams[2]
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
  var taylorGreenDensity  = flow_options.initParams[0]
  var taylorGreenPressure = flow_options.initParams[1]
  var taylorGreenVelocity = flow_options.initParams[2]
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
  c.rho         = flow_options.initParams[0]
  c.pressure    = flow_options.initParams[1]
  c.velocity[0] = flow_options.initParams[2] + ((rand_float()-0.5)*10.0)
  c.velocity[1] = flow_options.initParams[3] + ((rand_float()-0.5)*10.0)
  c.velocity[2] = flow_options.initParams[4] + ((rand_float()-0.5)*10.0)
end

ebb Flow.UpdateConservedFromPrimitive (c : fluidGrid)
  if c.in_interior then
    -- Equation of state: T = p / ( R * rho )
    var tmpTemperature = c.pressure / (fluid_options.gasConstant * c.rho)
    var velocity = c.velocity
    c.rhoVelocity = c.rho * c.velocity

    -- rhoE = rhoe (= rho * cv * T) + kineticEnergy + sgsEnergy
    var cv = fluid_options.gasConstant / (fluid_options.gamma - 1.0)
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

if particles_options.modeParticles then
  ebb Flow.AddParticlesCoupling (p : particles)

    -- WARNING: Assumes that deltaVelocityOverRelaxationTime and
    -- deltaTemperatureTerm have been computed previously, and that
    -- we have called the cell_locate kernel for the particles.
    -- (for example, when adding the flow coupling to the particles,
    -- which should be called before in the time stepper)

    -- WARNING: Uniform grid assumption
    var cellVolume = grid_dx * grid_dy * grid_dz

    -- Add contribution to momentum and energy equations from the previously
    -- computed deltaVelocityOverRelaxationTime and deltaTemperatureTerm
    p.cell.rhoVelocity_t += -p.mass * p.deltaVelocityOverRelaxationTime / cellVolume
    p.cell.rhoEnergy_t   += -p.deltaTemperatureTerm / cellVolume

    -- In case we want to hold a fixed temperature by subtracting
    -- a constant heat flux from the fluid, compute the avg.
    -- deltaTemperatureTerm to be adjusted later (note change in sign)
    if radiation_options.zeroAvgHeatSource then
      Flow.averageHeatSource += p.deltaTemperatureTerm / cellVolume
    end

  end
end

--------------------------------------------------------------
-- Holding avg. temperature fixed in the presence of radiation
--------------------------------------------------------------

if radiation_options.zeroAvgHeatSource then
  ebb Flow.AdjustHeatSource (c : fluidGrid)
    if c.in_interior then
    -- Remove a constant heat flux in all cells to balance with radiation.
    -- Note that this has been pre-computed before reaching this kernel (above).
    c.rhoEnergy_t += Flow.averageHeatSource
    end
  end
end

--------------
-- Body Forces
--------------

ebb Flow.AddBodyForces (c : fluidGrid)
  if c.in_interior then
    -- Add body forces (accelerations) to the momentum
    c.rhoVelocity_t += c.rho * flow_options.bodyForce

    -- Body force contribution to energy equation
    c.rhoEnergy_t += c.rho * L.dot(flow_options.bodyForce,c.velocity)

    -- Compute average heat source contribution in case we would
    -- like to subtract this later to recover a steady solution with radiation.
    --Flow.averageHeatSource += -c.rho * L.dot(flow_options.bodyForce,c.velocity)
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
    velocityX_XFace   /= (grid_dx*0.5)
    velocityY_XFace   /= (grid_dx*0.5)
    velocityZ_XFace   /= (grid_dx*0.5)
    temperature_XFace /= (grid_dx*0.5)

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
    velocityX_YFace   /= (grid_dy*0.5)
    velocityY_YFace   /= (grid_dy*0.5)
    velocityZ_YFace   /= (grid_dy*0.5)
    temperature_YFace /= (grid_dy*0.5)

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
    velocityX_ZFace   /= (grid_dz*0.5)
    velocityY_ZFace   /= (grid_dz*0.5)
    velocityZ_ZFace   /= (grid_dz*0.5)
    temperature_ZFace /= (grid_dz*0.5)

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
                      c(-1,0,0).dissipationFlux)/grid_dx
  end
end

ebb Flow.UpdateDissipationY (c : fluidGrid)
  if c.in_interior then
    c.dissipation += (c(0, 0,0).dissipationFlux -
                      c(0,-1,0).dissipationFlux)/grid_dy
  end
end

ebb Flow.UpdateDissipationZ (c : fluidGrid)
  if c.in_interior then
    c.dissipation += (c(0,0, 0).dissipationFlux -
                      c(0,0,-1).dissipationFlux)/grid_dz
  end
end

ebb Flow.ResetDissipation (c : fluidGrid)
  c.dissipation = 0.0
end

function Flow.UpdateDissipation (cells)
  fluidGrid:foreach(Flow.ResetDissipation)
  fluidGrid:foreach(Flow.ComputeDissipationX)
  fluidGrid:foreach(Flow.UpdateDissipationX)
  fluidGrid:foreach(Flow.ComputeDissipationY)
  fluidGrid:foreach(Flow.UpdateDissipationY)
  fluidGrid:foreach(Flow.ComputeDissipationZ)
  fluidGrid:foreach(Flow.UpdateDissipationZ)
end


-- WARNING: uniform grid assumption
local ebb averagePD ( c : fluidGrid )
  Flow.averagePD += c.PD * cellVolume
end
local ebb averageDissipation ( c : fluidGrid )
  Flow.averageDissipation += c.dissipation * cellVolume
end
local ebb averageK ( c : fluidGrid )
  Flow.averageK += 0.5 * c.rho * L.dot(c.velocity,c.velocity) * cellVolume
end
function Flow.UpdateTurbulentAverages(cells)

  cells:foreach(averagePD)
  Flow.averagePD:set(
      Flow.averagePD:get()/
      Flow.areaInterior:get())

  cells:foreach(averageDissipation)
  Flow.averageDissipation:set(
      Flow.averageDissipation:get()/
      Flow.areaInterior:get())

  cells:foreach(averageK)
  Flow.averageK:set(
      Flow.averageK:get()/
      Flow.areaInterior:get())
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

  --L.print(Flow.averagePD, Flow.averageDissipation, Flow.averageK, A)

  -- Add the forcing terms to the momentum and energy equations
  c.rhoVelocity_t += force
  c.rhoEnergy_t   += L.dot(force,c.velocity)

  -- Store the increment in the average energy source (to be subtracted later)
  -- WARNING: Uniform grid assumption
  var cellVolume = grid_dx * grid_dy * grid_dz
  Flow.averageFe += L.dot(force,c.velocity) * cellVolume
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
function Flow.AddTurbulentForcing (cells)

  -- Need to reset these averages somewhere

  Flow.averagePD:set(0.0)
  Flow.averageDissipation:set(0.0)
  Flow.averageFe:set(0.0)
  Flow.averageK:set(0.0)

  fluidGrid:foreach(Flow.UpdatePD)
  Flow.UpdateDissipation(cells)

  -- average PD and EPS
  Flow.UpdateTurbulentAverages(cells)

  -- Compute A & force, f_i
  -- Add rho * A * u_i to momentum, f_i*u_i to energy, accumulate f_i*u_i for average
  fluidGrid:foreach(Flow.AddTurbulentSource)

  -- Update average of the energy source
  Flow.averageFe:set(Flow.averageFe:get()/Flow.areaInterior:get())

  -- Subtract <f_e> from energy
  fluidGrid:foreach(Flow.AdjustTurbulentSource)

end

-------------------
-- Update functions
-------------------

-- Update flow variables using derivatives
-- Assumes 4th-order Runge-Kutta
ebb Flow.UpdateVars(c : fluidGrid)
  var deltaTime = TimeIntegrator.deltaTime
  if TimeIntegrator.stage == 1 then
    c.rho_new += (1.0/6.0) * deltaTime * c.rho_t
    c.rho = c.rho_old + 0.5 * deltaTime * c.rho_t
    c.rhoVelocity_new += (1.0/6.0) * deltaTime * c.rhoVelocity_t
    c.rhoVelocity = c.rhoVelocity_old + 0.5 * deltaTime * c.rhoVelocity_t
    c.rhoEnergy_new += (1.0/6.0) * deltaTime * c.rhoEnergy_t
    c.rhoEnergy = c.rhoEnergy_old + 0.5 * deltaTime * c.rhoEnergy_t
  elseif TimeIntegrator.stage == 2 then
    c.rho_new += (1.0/3.0) * deltaTime * c.rho_t
    c.rho = c.rho_old + 0.5 * deltaTime * c.rho_t
    c.rhoVelocity_new += (1.0/3.0) * deltaTime * c.rhoVelocity_t
    c.rhoVelocity = c.rhoVelocity_old + 0.5 * deltaTime * c.rhoVelocity_t
    c.rhoEnergy_new += (1.0/3.0) * deltaTime * c.rhoEnergy_t
    c.rhoEnergy = c.rhoEnergy_old + 0.5 * deltaTime * c.rhoEnergy_t
  elseif TimeIntegrator.stage == 3 then
    c.rho_new += (1.0/3.0) * deltaTime * c.rho_t
    c.rho = c.rho_old + 1.0 * deltaTime * c.rho_t
    c.rhoVelocity_new += (1.0/3.0) * deltaTime * c.rhoVelocity_t
    c.rhoVelocity = c.rhoVelocity_old + 1.0 * deltaTime * c.rhoVelocity_t
    c.rhoEnergy_new += (1.0/3.0) * deltaTime * c.rhoEnergy_t
    c.rhoEnergy = c.rhoEnergy_old + 1.0 * deltaTime * c.rhoEnergy_t
  else -- TimeIntegrator.stage == 4
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
local ebb UpdateGhostFieldsHelper(c_bnd, c_int, sign, bnd_velocity, bnd_temperature)
  -- Temporary variables for computing new halo state
  var rho         = L.double(0.0)
  var temp_wall   = L.double(0.0)
  var temperature = L.double(0.0)
  var velocity    = L.vec3d({0.0, 0.0, 0.0})

  -- Compute the Cv for updating the Energy equation
  var cv = fluid_options.gasConstant / (fluid_options.gamma - 1.0)

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
  rho = c_int.pressure / ( fluid_options.gasConstant * temperature )

  -- Update the boundary cell based on the values in the matching interior cell
  c_bnd.rhoBoundary         =  rho
  c_bnd.rhoVelocityBoundary =  rho*velocity
  c_bnd.rhoEnergyBoundary   =  rho * (cv * temperature +
                                      0.5*L.dot(velocity,velocity))
  c_bnd.velocityBoundary    =  velocity
  c_bnd.pressureBoundary    =  c_int.pressure
  c_bnd.temperatureBoundary =  temperature
end

ebb Flow.UpdateGhostFieldsStep1 (c : fluidGrid)
  if c.xneg_depth > 0 then
    UpdateGhostFieldsHelper(c, c( 1,0,0), x_sign, xneg_velocity, xneg_temperature)
  end
  if c.xpos_depth > 0 then
    UpdateGhostFieldsHelper(c, c(-1,0,0), x_sign, xpos_velocity, xpos_temperature)
  end
  if c.yneg_depth > 0 then
    UpdateGhostFieldsHelper(c, c(0, 1,0), y_sign, yneg_velocity, yneg_temperature)
  end
  if c.ypos_depth > 0 then
    UpdateGhostFieldsHelper(c, c(0,-1,0), y_sign, ypos_velocity, ypos_temperature)
  end
  if c.zneg_depth > 0 then
    UpdateGhostFieldsHelper(c, c(0,0, 1), z_sign, zneg_velocity, zneg_temperature)
  end
  if c.zpos_depth > 0 then
    UpdateGhostFieldsHelper(c, c(0,0,-1), z_sign, zpos_velocity, zpos_temperature)
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
function Flow.UpdateGhost()
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
    UpdateGhostThermodynamicsHelper(c, c( 1,0,0), xneg_temperature)
  end
  if c.xpos_depth > 0 then
    UpdateGhostThermodynamicsHelper(c, c(-1,0,0), xpos_temperature)
  end
  if c.yneg_depth > 0 then
    UpdateGhostThermodynamicsHelper(c, c(0, 1,0), yneg_temperature)
  end
  if c.ypos_depth > 0 then
    UpdateGhostThermodynamicsHelper(c, c(0,-1,0), ypos_temperature)
  end
  if c.zneg_depth > 0 then
    UpdateGhostThermodynamicsHelper(c, c(0,0, 1), zneg_temperature)
  end
  if c.zpos_depth > 0 then
    UpdateGhostThermodynamicsHelper(c, c(0,0,-1), zpos_temperature)
  end
end

ebb Flow.UpdateGhostThermodynamicsStep2 (c : fluidGrid)
  if c.in_boundary then
    c.pressure    = c.pressureBoundary
    c.temperature = c.temperatureBoundary
  end
end

function Flow.UpdateGhostThermodynamics()
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
    UpdateGhostVelocityHelper(c, c( 1,0,0), x_sign, xneg_velocity)
  end
  if c.xpos_depth > 0 then
    UpdateGhostVelocityHelper(c, c(-1,0,0), x_sign, xpos_velocity)
  end
  if c.yneg_depth > 0 then
    UpdateGhostVelocityHelper(c, c(0, 1,0), y_sign, yneg_velocity)
  end
  if c.ypos_depth > 0 then
    UpdateGhostVelocityHelper(c, c(0,-1,0), y_sign, ypos_velocity)
  end
  if c.zneg_depth > 0 then
    UpdateGhostVelocityHelper(c, c(0,0, 1), z_sign, zneg_velocity)
  end
  if c.zpos_depth > 0 then
    UpdateGhostVelocityHelper(c, c(0,0,-1), z_sign, zpos_velocity)
  end
end

ebb Flow.UpdateGhostVelocityStep2 (c : fluidGrid)
  if c.in_boundary then
    c.velocity = c.velocityBoundary
  end
end

function Flow.UpdateGhostVelocity()
  fluidGrid:foreach(Flow.UpdateGhostVelocityStep1)
  fluidGrid:foreach(Flow.UpdateGhostVelocityStep2)
end

-- Helper function for updating the conservatives to minimize repeated code
local ebb UpdateGhostConservedHelper (c_bnd, c_int, sign, bnd_velocity,
                                        bnd_temperature)

  -- Temporary variables for computing new halo state
  var rho         = L.double(0.0)
  var temp_wall   = L.double(0.0)
  var temperature = L.double(0.0)
  var velocity    = L.vec3d({0.0, 0.0, 0.0})

  -- Compute the Cv for updating the Energy equation
  var cv = fluid_options.gasConstant / (fluid_options.gamma - 1.0)

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
  rho = c_int.pressure / ( fluid_options.gasConstant * temperature )

  -- Update the boundary cell based on the values in the matching interior cell
  c_bnd.rhoBoundary         = rho
  c_bnd.rhoVelocityBoundary = rho*velocity
  c_bnd.rhoEnergyBoundary   = rho * (cv * temperature +
                                     0.5*L.dot(velocity,velocity))

end
ebb Flow.UpdateGhostConservedStep1 (c : fluidGrid)
  if c.xneg_depth > 0 then
    UpdateGhostConservedHelper(c, c( 1,0,0), x_sign, xneg_velocity, xneg_temperature)
  end
  if c.xpos_depth > 0 then
    UpdateGhostConservedHelper(c, c(-1,0,0), x_sign, xpos_velocity, xpos_temperature)
  end
  if c.yneg_depth > 0 then
    UpdateGhostConservedHelper(c, c(0, 1,0), y_sign, yneg_velocity, yneg_temperature)
  end
  if c.ypos_depth > 0 then
    UpdateGhostConservedHelper(c, c(0,-1,0), y_sign, ypos_velocity, ypos_temperature)
  end
  if c.zneg_depth > 0 then
    UpdateGhostConservedHelper(c, c(0,0, 1), z_sign, zneg_velocity, zneg_temperature)
  end
  if c.zpos_depth > 0 then
    UpdateGhostConservedHelper(c, c(0,0,-1), z_sign, zpos_velocity, zpos_temperature)
  end
end

ebb Flow.UpdateGhostConservedStep2 (c : fluidGrid)
  if c.in_boundary then
    c.rho         = c.rhoBoundary
    c.rhoVelocity = c.rhoVelocityBoundary
    c.rhoEnergy   = c.rhoEnergyBoundary
  end
end

function Flow.UpdateGhostConserved()
  fluidGrid:foreach(Flow.UpdateGhostConservedStep1)
  fluidGrid:foreach(Flow.UpdateGhostConservedStep2)
end

ebb Flow.UpdateAuxiliaryThermodynamics (c : fluidGrid)
  if c.in_interior then
    var kineticEnergy = 0.5 * c.rho * L.dot(c.velocity,c.velocity)
    var pressure  = (fluid_options.gamma - 1.0) *( c.rhoEnergy - kineticEnergy )
    c.pressure    = pressure
    c.temperature = pressure / ( fluid_options.gasConstant * c.rho )
  end
end

---------------------
-- Velocity gradients
---------------------

-- WARNING: non-uniform grid assumption
ebb Flow.ComputeVelocityGradientAll (c : fluidGrid)
  if c.in_interior then
    c.velocityGradientX = 0.5*(c(1,0,0).velocity - c(-1,0,0).velocity)/grid_dx
    c.velocityGradientY = 0.5*(c(0,1,0).velocity - c(0,-1,0).velocity)/grid_dy
    c.velocityGradientZ = 0.5*(c(0,0,1).velocity - c(0,0,-1).velocity)/grid_dz
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
    UpdateGhostVelocityGradientHelper(c, c( 1,0,0), x_sign)
  end
  if c.xpos_depth > 0 then
    UpdateGhostVelocityGradientHelper(c, c(-1,0,0), x_sign)
  end
  if c.yneg_depth > 0 then
    UpdateGhostVelocityGradientHelper(c, c(0, 1,0), y_sign)
  end
  if c.ypos_depth > 0 then
    UpdateGhostVelocityGradientHelper(c, c(0,-1,0), y_sign)
  end
  if c.zneg_depth > 0 then
    UpdateGhostVelocityGradientHelper(c, c(0,0, 1), z_sign)
  end
  if c.zpos_depth > 0 then
    UpdateGhostVelocityGradientHelper(c, c(0,0,-1), z_sign)
  end
end

ebb Flow.UpdateGhostVelocityGradientStep2 (c : fluidGrid)
  if c.in_boundary then
    c.velocityGradientX = c.velocityGradientXBoundary
    c.velocityGradientY = c.velocityGradientYBoundary
    c.velocityGradientZ = c.velocityGradientZBoundary
  end
end

-- Calculation of spectral radii for clf-based delta time
local maxConvectiveSpectralRadius     = L.Global('maxC', L.double, 0.0)
local maxViscousSpectralRadius        = L.Global('maxV', L.double, 0.0)
local maxHeatConductionSpectralRadius = L.Global('maxH', L.double, 0.0)

-- WARNING: update cellVolume computation for non-uniform grids
local dXYZInverseSquare = L.Constant(L.double,
                                     1.0/grid_dx:get() * 1.0/grid_dx:get() +
                                     1.0/grid_dy:get() * 1.0/grid_dy:get() +
                                     1.0/grid_dz:get() * 1.0/grid_dz:get())
local ebb calculateConvectiveSpectralRadius     ( c : fluidGrid )
  -- Convective spectral radii
  -- WARNING: uniform grid assumption
  c.convectiveSpectralRadius =
   (L.fabs(c.velocity[0])/grid_dx  +
    L.fabs(c.velocity[1])/grid_dy  +
    L.fabs(c.velocity[2])/grid_dz  +
    GetSoundSpeed(c.temperature) * L.sqrt(dXYZInverseSquare))

  maxConvectiveSpectralRadius max= c.convectiveSpectralRadius
end
local ebb calculateViscousSpectralRadius        ( c : fluidGrid )
  -- Viscous spectral radii (including sgs model component)
  var dynamicViscosity = GetDynamicViscosity(c.temperature)
  var eddyViscosity = c.sgsEddyViscosity
  c.viscousSpectralRadius =
   (2.0 * ( dynamicViscosity + eddyViscosity ) /
    c.rho * dXYZInverseSquare) * 4.0

  maxViscousSpectralRadius max= c.viscousSpectralRadius
end
local ebb calculateHeatConductionSpectralRadius ( c : fluidGrid )
  var dynamicViscosity  = GetDynamicViscosity(c.temperature)

  -- Heat conduction spectral radii (including sgs model component)
  var cv = fluid_options.gasConstant / (fluid_options.gamma - 1.0)
  var cp = fluid_options.gamma * cv

  var kappa = cp / fluid_options.prandtl *  dynamicViscosity

  c.heatConductionSpectralRadius =
     ((kappa + c.sgsEddyKappa) / (cv * c.rho) * dXYZInverseSquare) * 4.0
  maxHeatConductionSpectralRadius max= c.heatConductionSpectralRadius
end
function Flow.CalculateSpectralRadii(cells)
  cells:foreach(calculateConvectiveSpectralRadius)
  cells:foreach(calculateViscousSpectralRadius)
  cells:foreach(calculateHeatConductionSpectralRadius)
end

-------------
-- Statistics
-------------

local ebb averagePressure       ( c : fluidGrid )
  if c.in_interior then
    Flow.averagePressure          += c.pressure * cellVolume
  end
end
local ebb averageTemperature    ( c : fluidGrid )
  if c.in_interior then
    Flow.averageTemperature       += c.temperature * cellVolume
  end
end
local ebb averageKineticEnergy  ( c : fluidGrid )
  if c.in_interior then
    Flow.averageKineticEnergy     += c.kineticEnergy * cellVolume
  end
end
local ebb minTemperature        ( c : fluidGrid )
  if c.in_interior then
    Flow.minTemperature         min= c.temperature
  end
end
local ebb maxTemperature        ( c : fluidGrid )
  if c.in_interior then
    Flow.maxTemperature         max= c.temperature
  end
end
function Flow.IntegrateQuantities(cells)
  cells:foreach(averagePressure      )
  cells:foreach(averageTemperature   )
  cells:foreach(averageKineticEnergy )
  cells:foreach(minTemperature       )
  cells:foreach(maxTemperature       )
end


------------
-- PARTICLES
------------

-- Put a guard around the entire particles section so that we don't invoke
-- any of these kernels when the particles are turned off.
if particles_options.modeParticles then

  ebb Particles.LocateInCells( p : particles )
    p.cell = fluidGrid.locate(p.position)
  end

  -- Locate particles in cells
  function Particles.Locate()
    particles:foreach(Particles.LocateInCells)
  end

  -- Initialize temporaries for time stepper
  ebb Particles.InitializeTemporaries (p : particles)
    p.position_old    = p.position
    p.velocity_old    = p.particle_velocity
    p.temperature_old = p.particle_temperature
    p.position_new    = p.position
    p.velocity_new    = p.particle_velocity
    p.temperature_new = p.particle_temperature
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
    if particles_options.particleType == ParticleType.Fixed then
      -- Don't move the particle
    elseif particles_options.particleType == ParticleType.Free then
      p.position_t += p.particle_velocity
    else L.assert(false) end

    -- Relaxation time for small particles
    -- - particles Reynolds number (set to zero for Stokesian)
    var particleReynoldsNumber = 0.0
    --(p.density * norm(flowVelocity - p.velocity) * p.diameter) / flowDynamicViscosity
    var relaxationTime =
      ( p.density * L.pow(p.diameter,2) / (18.0 * flowDynamicViscosity) ) /
      ( 1.0 + 0.15 * L.pow(particleReynoldsNumber,0.687) )
    p.deltaVelocityOverRelaxationTime =
      (flowVelocity - p.particle_velocity) / relaxationTime
    p.deltaTemperatureTerm =
      pi * L.pow(p.diameter, 2) * particles_options.convective_coefficient *
      (flowTemperature - p.particle_temperature)

    -- Update the particle velocity and temperature
    if particles_options.particleType == ParticleType.Fixed then
      p.velocity_t  = {0.0,0.0,0.0} -- Don't move the particle
    elseif particles_options.particleType == ParticleType.Free then
      p.velocity_t += p.deltaVelocityOverRelaxationTime
    else L.assert(false) end
    p.temperature_t += p.deltaTemperatureTerm /
      (p.mass * particles_options.heat_capacity)

  end

  --------------
  -- Body forces
  --------------

  ebb Particles.AddBodyForces (p : particles)
    p.velocity_t += particles_options.bodyForce
  end

  ------------
  -- Radiation
  ------------

  -- Update particle variables using derivatives
  ebb Particles.UpdateVars(p : particles)
    var deltaTime = TimeIntegrator.deltaTime
    if TimeIntegrator.stage == 1 then
      p.position_new += (1.0/6.0) * deltaTime * p.position_t
      p.position = p.position_old + 0.5 * deltaTime * p.position_t
      p.velocity_new += (1.0/6.0) * deltaTime * p.velocity_t
      p.particle_velocity = p.velocity_old + 0.5 * deltaTime * p.velocity_t
      p.temperature_new += (1.0/6.0) * deltaTime * p.temperature_t
      p.particle_temperature = p.temperature_old +
        0.5 * deltaTime * p.temperature_t
    elseif TimeIntegrator.stage == 2 then
      p.position_new += (1.0/3.0) * deltaTime * p.position_t
      p.position = p.position_old + 0.5 * deltaTime * p.position_t
      p.velocity_new += (1.0/3.0) * deltaTime * p.velocity_t
      p.particle_velocity = p.velocity_old + 0.5 * deltaTime * p.velocity_t
      p.temperature_new += (1.0/3.0) * deltaTime * p.temperature_t
      p.particle_temperature = p.temperature_old +
        0.5 * deltaTime * p.temperature_t
    elseif TimeIntegrator.stage == 3 then
      p.position_new += (1.0/3.0) * deltaTime * p.position_t
      p.position = p.position_old + 1.0 * deltaTime * p.position_t
      p.velocity_new += (1.0/3.0) * deltaTime * p.velocity_t
      p.particle_velocity = p.velocity_old + 1.0 * deltaTime * p.velocity_t
      p.temperature_new += (1.0/3.0) * deltaTime * p.temperature_t
      p.particle_temperature = p.temperature_old +
        1.0 * deltaTime * p.temperature_t
    else -- TimeIntegrator.stage == 4
      p.position = p.position_new + (1.0/6.0) * deltaTime * p.position_t
      p.particle_velocity = p.velocity_new + (1.0/6.0) * deltaTime * p.velocity_t
      p.particle_temperature = p.temperature_new +
        (1.0/6.0) * deltaTime * p.temperature_t
    end
  end

  ebb Particles.UpdateAuxiliaryStep1 (p : particles)

    -- Initialize position and velocity before we check for wall collisions

    p.position_ghost[0]   = p.position[0]
    p.position_ghost[1]   = p.position[1]
    p.position_ghost[2]   = p.position[2]
    p.velocity_ghost[0]   = p.particle_velocity[0]
    p.velocity_ghost[1]   = p.particle_velocity[1]
    p.velocity_ghost[2]   = p.particle_velocity[2]
    p.velocity_t_ghost[0] = p.velocity_t[0]
    p.velocity_t_ghost[1] = p.velocity_t[1]
    p.velocity_t_ghost[2] = p.velocity_t[2]

    -- Check here for particles exiting the domain. For periodic
    -- boundaries, the particle is transported to the matching periodic
    -- face. For symmetry or wall boundaries, an elastic collision is
    -- assumed. To start, the collision is perfectly elastic.

    -- Left X boundary
    if p.position[0] < gridOriginInteriorX then
      if grid_options.xBCLeftParticles == ParticleBC.Permeable then
        p.position_ghost[0] = p.position[0] + grid_options.xWidth
      elseif grid_options.xBCLeftParticles == ParticleBC.Solid then

        -- Set the position to be on the wall
        p.position_ghost[0] = gridOriginInteriorX

        -- Apply an impulse to kick particle away from the wall
        var impulse = - (1.0+particles_options.restitution_coefficient) *
          p.particle_velocity[0]
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
    if p.position[0] > gridOriginInteriorX + grid_options.xWidth then
      if grid_options.xBCRightParticles == ParticleBC.Permeable then
        p.position_ghost[0] = p.position[0] - grid_options.xWidth
      elseif grid_options.xBCRightParticles == ParticleBC.Solid then

        -- Set the position to be on the wall
        p.position_ghost[0] = gridOriginInteriorX + grid_options.xWidth

        -- Apply an impulse to kick particle away from the wall
        var impulse = - (1.0+particles_options.restitution_coefficient) *
          p.particle_velocity[0]
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
    if p.position[1] < gridOriginInteriorY then
      if grid_options.yBCLeftParticles == ParticleBC.Permeable then
        p.position_ghost[1] = p.position[1] + grid_options.yWidth
      elseif grid_options.yBCLeftParticles == ParticleBC.Solid then

        -- Set the position to be on the wall
        p.position_ghost[1] = gridOriginInteriorY

        -- Apply an impulse to kick particle away from the wall
        var impulse = - (1.0+particles_options.restitution_coefficient) *
          p.particle_velocity[1]
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
    if p.position[1] > gridOriginInteriorY + grid_options.yWidth then
      if grid_options.yBCRightParticles == ParticleBC.Permeable then
        p.position_ghost[1] = p.position[1] - grid_options.yWidth
      elseif grid_options.yBCRightParticles == ParticleBC.Solid then

        -- Set the position to be on the wall
        p.position_ghost[1] = gridOriginInteriorY + grid_options.yWidth

        -- Apply an impulse to kick particle away from the wall
        var impulse = - (1.0+particles_options.restitution_coefficient) *
          p.particle_velocity[1]
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
    if p.position[2] < gridOriginInteriorZ then
      if grid_options.zBCLeftParticles == ParticleBC.Permeable then
        p.position_ghost[2] = p.position[2] + grid_options.zWidth
      elseif grid_options.zBCLeftParticles == ParticleBC.Solid then

        -- Set the position to be on the wall
        p.position_ghost[2] = gridOriginInteriorZ

        -- Apply an impulse to kick particle away from the wall
        var impulse = - (1.0+particles_options.restitution_coefficient) *
          p.particle_velocity[2]
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
    if p.position[2] > gridOriginInteriorZ + grid_options.zWidth then
      if grid_options.zBCRightParticles == ParticleBC.Permeable then
        p.position_ghost[2] = p.position[2] - grid_options.zWidth
      elseif grid_options.zBCRightParticles == ParticleBC.Solid then

        -- Set the position to be on the wall
        p.position_ghost[2] = gridOriginInteriorZ + grid_options.zWidth

        -- Apply an impulse to kick particle away from the wall
        var impulse = - (1.0+particles_options.restitution_coefficient) *
          p.particle_velocity[2]
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
    p.position          = p.position_ghost
    p.particle_velocity = p.velocity_ghost
    p.velocity_t        = p.velocity_t_ghost
  end

  if radiation_options.radiationType == RadiationType.DOM then

    ebb Radiation.InitializeCell(c : radiationGrid)
      for m = 0,radiation_options.numAngles do
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

    ebb Radiation.ClearAccumulators(c : radiationGrid)
      c.acc_d2 = 0.0
      c.acc_d2t4 = 0.0
    end

    ebb Radiation.AccumulateParticleValues(p : particles)
      p.cell.to_Radiation.acc_d2 +=
        L.pow(p.diameter,2.0)
      p.cell.to_Radiation.acc_d2t4 +=
        L.pow(p.diameter,2.0) * L.pow(p.particle_temperature,4.0)
    end
    Radiation.AccumulateParticleValues._MANUAL_PARAL = true

    ebb Radiation.UpdateFieldValues(c : radiationGrid)
      c.sigma = c.acc_d2 * pi
        * (radiation_options.qa + radiation_options.qs)
        / (4.0 * radiation_options.cellVolume)
      if c.acc_d2 == 0.0 then
        c.Ib = 0.0
      else
        c.Ib = SB * c.acc_d2t4 / (pi * c.acc_d2)
      end
    end

    ebb Particles.AbsorbRadiation (p : particles)
      var t4 = L.pow(p.particle_temperature,4.0)
      var alpha = pi * radiation_options.qa * L.pow(p.diameter,2.0)
        * (p.cell.to_Radiation.G - 4.0 * SB * t4) / 4.0
      p.temperature_t += alpha / (p.mass * particles_options.heat_capacity)
    end
    Particles.AbsorbRadiation._MANUAL_PARAL = true

  end

  ---------
  -- Feeder
  ---------

  -- Convert the cell coordinates to a number in 0..size(interior)-1
  ebb Flow.InteriorCellNumber(c)
    var xid = L.int64(L.xid(c))
    var yid = L.int64(L.yid(c))
    var zid = L.int64(L.zid(c))
    if not xBCPeriodic then xid = xid-L.int64(1) end
    if not yBCPeriodic then yid = yid-L.int64(1) end
    if not zBCPeriodic then zid = zid-L.int64(1) end
    return (zid * L.int64(grid_options.xnum) * L.int64(grid_options.ynum) +
            yid * L.int64(grid_options.xnum) +
            xid)
  end

  -- Pick a diameter value according to a random distribution,
  -- with given mean value and maximum deviation.
  ebb Particles.RandomDiameter()
    return (rand_float() - 0.5) * particles_options.diameter_maxDeviation +
      particles_options.diameter_mean
  end

  -- Calculate particle velocity from underlying flow velocity.
  if particles_options.particleType == ParticleType.Fixed then
    ebb Particles.VelocityFromFlow(cell, position)
      return {0.0, 0.0, 0.0} -- Don't move the particle
    end
  elseif particles_options.particleType == ParticleType.Free then
    ebb Particles.VelocityFromFlow(cell, position)
      return InterpolateTriVelocity(cell, position)
    end
  else assert(false) end

  -- Insert one particle on each cell, with a small probability.
  -- TODO: Inserting exactly at the center, to avoid the need for stencils.
  ebb Flow.InsertParticlesAtRandom(c : fluidGrid)
    if c.in_interior and
       Flow.InteriorCellNumber(c) < Particles.limit and
       rand_float() < 0.01 then
      insert {
        cell = c,
        position = c.center,
        particle_velocity = c.velocity,
        density = particles_options.density,
        particle_temperature = particles_options.initialTemperature,
        diameter = Particles.RandomDiameter()
      } into particles
      Particles.number += 1
    end
  end

  -- Particle feeder
  function Particles.Feed()
    -- For now, insert at random just for testing.
    Particles.limit:set(particles_options.maximum_num - Particles.number:get())
    fluidGrid:foreach(Flow.InsertParticlesAtRandom)
  end

  ebb Particles.DeleteEscapingParticles(p: particles)
    var min_x = grid_originX
    var max_x = grid_originX + grid_widthX
    var min_y = grid_originY
    var max_y = grid_originY + grid_widthY
    var min_z = grid_originZ
    var max_z = grid_originZ + grid_widthZ
    var pos = p.position
    if (pos[0] > max_x or pos[0] < min_x  or
        pos[1] > max_y or pos[1] < min_y  or
        pos[2] > max_z or pos[2] < min_z) then
      delete(p)
      Particles.number -= L.int64(1)
    end
  end

  -- Particle collector
  function Particles.Collect()
    -- For now, delete anything that leaves the domain.
    particles:foreach(Particles.DeleteEscapingParticles)
  end

  -------------
  -- Statistics
  -------------

  ebb Particles.IntegrateQuantities (p : particles)
    Particles.averageTemperature += p.particle_temperature
  end

  ebb Particles.numberOfParticles (p : particles)
    Particles.number += L.int64(1)
  end

end

-----------------------------------------------------------------------------
--[[                                MAIN FUNCTIONS                       ]]--
-----------------------------------------------------------------------------

-------
-- FLOW
-------

function Flow.InitializePrimitives()
  if flow_options.initCase == InitCase.Uniform then
    fluidGrid:foreach(Flow.InitializeUniform)
  elseif flow_options.initCase == InitCase.TaylorGreen2DVortex then
    fluidGrid:foreach(Flow.InitializeTaylorGreen2D)
  elseif flow_options.initCase == InitCase.TaylorGreen3DVortex then
    fluidGrid:foreach(Flow.InitializeTaylorGreen3D)
  elseif flow_options.initCase == InitCase.Perturbed then
    fluidGrid:foreach(Flow.InitializePerturbed)
  elseif flow_options.initCase == InitCase.Restart then
    fluidGrid:Load({'rho','pressure','velocity'},
                    io_options.outputFileNamePrefix .. 'restart_fluid_' ..
                      time_options.restartIter .. '.hdf')
  else assert(false) end
end

function Flow.UpdateGhostVelocityGradient()
  fluidGrid:foreach(Flow.UpdateGhostVelocityGradientStep1)
  fluidGrid:foreach(Flow.UpdateGhostVelocityGradientStep2)
end

-- Routine that computes the inviscid flux through the face of
-- any two adjacent cells with a centered scheme. The left cell (c_l),
-- right cell (c_r), and coordinate direction (x = 0, y = 1, or z = 2)
-- are the inputs.
local function mkCenteredInviscidFlux(direction)
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
    var s = spatial_options.split
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
    velocityX_XFace   /= (grid_dx*0.5)
    velocityY_XFace   /= (grid_dx*0.5)
    velocityZ_XFace   /= (grid_dx*0.5)
    temperature_XFace /= (grid_dx*0.5)

    -- Tensor components (at face)
    var sigmaXX = muFace * ( 4.0 * velocityX_XFace -
                             2.0 * velocityY_YFace -
                             2.0 * velocityZ_ZFace ) / 3.0
    var sigmaYX = muFace * ( velocityY_XFace + velocityX_YFace )
    var sigmaZX = muFace * ( velocityZ_XFace + velocityX_ZFace )
    var usigma  = velocityFace[0] * sigmaXX +
                  velocityFace[1] * sigmaYX +
                  velocityFace[2] * sigmaZX
    var cp = fluid_options.gamma * fluid_options.gasConstant /
             (fluid_options.gamma - 1.0)
    var heatFlux = - (cp*muFace/fluid_options.prandtl)*temperature_XFace

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
    velocityX_YFace   /= (grid_dy*0.5)
    velocityY_YFace   /= (grid_dy*0.5)
    velocityZ_YFace   /= (grid_dy*0.5)
    temperature_YFace /= (grid_dy*0.5)

    -- Tensor components (at face)
    var sigmaXY = muFace * ( velocityX_YFace + velocityY_XFace )
    var sigmaYY = muFace * ( 4.0 * velocityY_YFace -
                             2.0 * velocityX_XFace -
                             2.0 * velocityZ_ZFace ) / 3.0
    var sigmaZY = muFace * ( velocityZ_YFace + velocityY_ZFace )
    var usigma  = velocityFace[0] * sigmaXY +
                  velocityFace[1] * sigmaYY +
                  velocityFace[2] * sigmaZY
    var cp = fluid_options.gamma * fluid_options.gasConstant /
             (fluid_options.gamma - 1.0)
    var heatFlux = - (cp*muFace/fluid_options.prandtl)*temperature_YFace

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
    velocityX_ZFace   /= (grid_dz*0.5)
    velocityY_ZFace   /= (grid_dz*0.5)
    velocityZ_ZFace   /= (grid_dz*0.5)
    temperature_ZFace /= (grid_dz*0.5)

    -- Tensor components (at face)
    var sigmaXZ = muFace * ( velocityX_ZFace + velocityZ_XFace )
    var sigmaYZ = muFace * ( velocityY_ZFace + velocityZ_YFace )
    var sigmaZZ = muFace * ( 4.0 * velocityZ_ZFace -
                             2.0 * velocityX_XFace -
                             2.0 * velocityY_YFace ) / 3.0
    var usigma  = velocityFace[0] * sigmaXZ +
                  velocityFace[1] * sigmaYZ +
                  velocityFace[2] * sigmaZZ
    var cp = fluid_options.gamma * fluid_options.gasConstant /
             (fluid_options.gamma - 1.0)
    var heatFlux = - (cp*muFace/fluid_options.prandtl)*temperature_ZFace

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

ebb Flow.AddUpdateUsingFlux(c : fluidGrid)
  --if c.in_interior or c.xneg_depth == 1 then
  if c.in_interior then --or c.xneg_depth == 1 then
    c.rho_t += -(c( 0,0,0).rhoFluxX -
                 c(-1,0,0).rhoFluxX)/grid_dx
    c.rhoVelocity_t += -(c( 0,0,0).rhoVelocityFluxX -
                         c(-1,0,0).rhoVelocityFluxX)/grid_dx
    c.rhoEnergy_t += -(c( 0,0,0).rhoEnergyFluxX -
                       c(-1,0,0).rhoEnergyFluxX)/grid_dx
  --end
  --if c.in_interior or c.zneg_depth == 1 then
    c.rho_t += -(c(0, 0,0).rhoFluxY -
                 c(0,-1,0).rhoFluxY)/grid_dy
    c.rhoVelocity_t += -(c(0, 0,0).rhoVelocityFluxY -
                         c(0,-1,0).rhoVelocityFluxY)/grid_dy
    c.rhoEnergy_t += -(c(0, 0,0).rhoEnergyFluxY -
                       c(0,-1,0).rhoEnergyFluxY)/grid_dy
  --end
  --if c.in_interior or c.zneg_depth == 1 then
    c.rho_t += -(c(0,0, 0).rhoFluxZ -
                 c(0,0,-1).rhoFluxZ)/grid_dz
    c.rhoVelocity_t += -(c(0,0, 0).rhoVelocityFluxZ -
                         c(0,0,-1).rhoVelocityFluxZ)/grid_dz
    c.rhoEnergy_t += -(c(0,0, 0).rhoEnergyFluxZ -
                       c(0,0,-1).rhoEnergyFluxZ)/grid_dz
  end
end

function Flow.AddFluxes()
  fluidGrid:foreach(Flow.AddGetFlux)
  fluidGrid:foreach(Flow.AddUpdateUsingFlux)
end

function Flow.ComputeVelocityGradients()
  fluidGrid:foreach(Flow.ComputeVelocityGradientAll)
end

function Flow.UpdateAuxiliaryVelocityConservedAndGradients()
  fluidGrid:foreach(Flow.UpdateAuxiliaryVelocity)
  Flow.UpdateGhostConserved()
  Flow.UpdateGhostVelocity()
  Flow.ComputeVelocityGradients()
end

function Flow.UpdateAuxiliary()
  Flow.UpdateAuxiliaryVelocityConservedAndGradients()
  fluidGrid:foreach(Flow.UpdateAuxiliaryThermodynamics)
  Flow.UpdateGhostThermodynamics()
end


------------
-- PARTICLES
------------

-- Put a guard around all particle kernels in case they're inactive.

if particles_options.modeParticles then

  -- Insert one particle at the center of each cell (plus a tiny offset to help
  -- the interpolation verification checks).
  ebb Flow.InsertParticlesUniform(c : fluidGrid)
    if c.in_interior then
      var cellId = Flow.InteriorCellNumber(c)
      var numCells = L.int64(grid_options.xnum) *
                     L.int64(grid_options.ynum) *
                     L.int64(grid_options.znum)
      if cellId == L.int64(0) or
         (cellId-L.int64(1)) * (Particles.limit-L.int64(1)) / (numCells-L.int64(1)) <
           cellId * (Particles.limit-L.int64(1)) / (numCells-L.int64(1)) then
        insert {
          cell = c,
          position = c.center,
          particle_velocity = c.velocity,
          density = particles_options.density,
          particle_temperature = particles_options.initialTemperature,
          diameter = particles_options.diameter_mean
        } into particles
        Particles.number += L.int64(1)
      end
    end
  end

  function Particles.InitializePrimitives()
    if particles_options.initParticles == InitParticles.Random then
      error("Random particle initialization is disabled")
    elseif particles_options.initParticles == InitParticles.Restart then
      particles:Load(
        {'cell','position','particle_velocity','particle_temperature','diameter'},
        io_options.outputFileNamePrefix .. 'restart_particles_' ..
          tostring(particles_options.restartParticleIter) .. '.hdf')
      particles.density:Fill(particles_options.density)
    elseif particles_options.initParticles == InitParticles.Uniform then
      Particles.number:set(particles_options.num)
      M.INLINE(particles_init_uniform.InitParticlesUniform)
    else assert(false) end
  end

  function Particles.UpdateAuxiliary()
    particles:foreach(Particles.UpdateAuxiliaryStep1)
    particles:foreach(Particles.UpdateAuxiliaryStep2)
  end

end -- particles_options.modeParticles


------------------
-- TIME INTEGRATOR
------------------

function TimeIntegrator.SetupTimeStep()
  fluidGrid:foreach(Flow.InitializeTemporaries)
  if particles_options.modeParticles then
    -- Particles.Feed()
    particles:foreach(Particles.InitializeTemporaries)
  end
end

function TimeIntegrator.ConcludeTimeStep()
  if particles_options.modeParticles then
    Particles.Collect()
  end
end

function TimeIntegrator.InitializeTimeDerivatives()
  fluidGrid:foreach(Flow.InitializeTimeDerivatives)
  if particles_options.modeParticles then
    particles:foreach(Particles.InitializeTimeDerivatives)
  end
end

function TimeIntegrator.UpdateAuxiliary()
  Flow.UpdateAuxiliary()
  if particles_options.modeParticles then
    Particles.UpdateAuxiliary()
  end
end

function TimeIntegrator.UpdateTime()
  -- HACK
  TimeIntegrator.simTime:set(TimeIntegrator.timeOld:get() +
                             -- stage = 1 => 0.5
                             -- stage = 2 => 0.5
                             -- stage = 3 => 1.0
                             -- stage = 4 => 1.0
                             0.5 * (1 + TimeIntegrator.stage:get() / 3) *
                             TimeIntegrator.deltaTime:get())
end

function TimeIntegrator.InitializeVariables()

  fluidGrid:foreach(Flow.InitializeCell)

  -- Initialize several grid related entitities
  fluidGrid:foreach(Flow.InitializeCenterCoordinates)

  -- Set initial condition for the flow and all auxiliary flow variables
  Flow.InitializePrimitives()
  fluidGrid:foreach(Flow.UpdateConservedFromPrimitive)
  Flow.UpdateAuxiliary()
  Flow.UpdateGhost()

  -- Initialize the particles (position, velocity, temp, diameter, locate)
  if particles_options.modeParticles then
    Particles.InitializePrimitives()
  end

end

function TimeIntegrator.ComputeDFunctionDt()

  -- Compute flow convective, viscous, and body force residuals
  Flow.UpdateGhostVelocityGradient()
  Flow.AddFluxes()
  if radiation_options.zeroAvgHeatSource then
    Flow.averageHeatSource:set(0.0)
  end
  fluidGrid:foreach(Flow.AddBodyForces)

  -- FIXME: turbulent-related tasks should be revised
  if flow_options.turbForcing then
    Flow.AddTurbulentForcing(fluidGrid.interior)
  end

  -- Compute residuals for the particles (locate all particles first)
  if particles_options.modeParticles then

    Particles.Locate()
    particles:foreach(Particles.AddFlowCoupling)

    if particles_options.particleType == ParticleType.Free then
      particles:foreach(Particles.AddBodyForces)
    end

    if radiation_options.radiationType == RadiationType.Algebraic then
      M.INLINE(radiation.AddRadiation)
    elseif radiation_options.radiationType == RadiationType.DOM then
      -- Compute radiation field values from particles
      radiationGrid:foreach(Radiation.ClearAccumulators)
      particles:foreach(Radiation.AccumulateParticleValues)
      radiationGrid:foreach(Radiation.UpdateFieldValues)
      M.INLINE(radiation.ComputeRadiationField)
      -- Absorb radiation into each particle
      particles:foreach(Particles.AbsorbRadiation)
   end

    -- Compute two-way coupling in momentum and energy
    if particles_options.twoWayCoupling then
      particles:foreach(Flow.AddParticlesCoupling)
    end

  end

  -- In case we want to hold flow temp fixed with radiation active
  --print(Flow.averageHeatSource:get())

  if radiation_options.zeroAvgHeatSource then
    Flow.averageHeatSource:set(Flow.averageHeatSource:get()/
                                 Flow.numberOfInteriorCells:get())
    fluidGrid:foreach(Flow.AdjustHeatSource)
  end

end

function TimeIntegrator.UpdateSolution()
  fluidGrid:foreach(Flow.UpdateVars)
  if particles_options.modeParticles then
    particles:foreach(Particles.UpdateVars)
  end
end

function TimeIntegrator.AdvanceTimeStep()

  TimeIntegrator.SetupTimeStep()
  TimeIntegrator.timeOld:set(TimeIntegrator.simTime:get())

  TimeIntegrator.stage:set(1)
  M.WHILE(M.LT(TimeIntegrator.stage:get(), 5))
    TimeIntegrator.InitializeTimeDerivatives()
    TimeIntegrator.ComputeDFunctionDt()
    TimeIntegrator.UpdateSolution()
    TimeIntegrator.UpdateAuxiliary()
    TimeIntegrator.UpdateTime()
    TimeIntegrator.stage:set(TimeIntegrator.stage:get() + 1)
    -- HACK: Move escaping particle deletion here, to appease the SPMD
    -- transformation. It should be fine to do this multiple times.
    TimeIntegrator.ConcludeTimeStep()
  M.END()

  TimeIntegrator.timeStep:set(TimeIntegrator.timeStep:get() + 1)

end

function TimeIntegrator.CalculateDeltaTime()

  -- Check whether we are imposing a delta time or basing it on the CFL,
  -- i.e. a negative CFL was imposed in the config
  if time_options.cfl < 0 then

    -- Impose a fixed time step from the config
    TimeIntegrator.deltaTime:set(time_options.delta_time)

  else

    -- Calculate the convective, viscous, and heat spectral radii
    Flow.CalculateSpectralRadii(fluidGrid)

    local maxV = maxViscousSpectralRadius:get()
    local maxH = maxHeatConductionSpectralRadius:get()
    local maxC = maxConvectiveSpectralRadius:get()

    -- Calculate diffusive spectral radius as the maximum between
    -- heat conduction and convective spectral radii
    -- Calculate global spectral radius as the maximum between the convective
    -- and diffusive spectral radii
    -- Delta time using the CFL and max spectral radius for stability
    TimeIntegrator.deltaTime:set(time_options.cfl /
                                 M.MAX(maxC,M.MAX(maxV,maxH)))
  end

end


-------------
-- STATISTICS
-------------

function Statistics.ResetSpatialAverages()
  Flow.averagePressure:set(0.0)
  Flow.averageTemperature:set(0.0)
  Flow.averageKineticEnergy:set(0.0)
  Flow.minTemperature:set(math.huge)
  Flow.maxTemperature:set(-math.huge)
  Flow.averagePD:set(0.0)
  Flow.averageDissipation:set(0.0)
  Particles.averageTemperature:set(0.0)
end

function Statistics.UpdateSpatialAverages(grid, particles)
  -- Flow
  Flow.averagePressure:set(
    Flow.averagePressure:get() / Flow.areaInterior:get()
  )
  Flow.averageTemperature:set(
    Flow.averageTemperature:get() / Flow.areaInterior:get()
  )
  Flow.averageKineticEnergy:set(
    Flow.averageKineticEnergy:get() / Flow.areaInterior:get()
  )

  -- Particles
  if particles_options.modeParticles then
    Particles.averageTemperature:set(
      Particles.averageTemperature:get() / Particles.number:get()
    )
  end
end

function Statistics.ComputeSpatialAverages()
  Statistics.ResetSpatialAverages()
  Flow.IntegrateQuantities(fluidGrid)
  if particles_options.modeParticles then
    particles:foreach(Particles.IntegrateQuantities)
  end
  Statistics.UpdateSpatialAverages(fluidGrid, particles)
end


-----
-- IO
-----

function IO.WriteConsoleOutput()
  M.IF(M.EQ(TimeIntegrator.timeStep:get() % time_options.consoleFrequency, 0))
    -- Output log headers at a specified frequency
    M.IF(M.EQ(TimeIntegrator.timeStep:get() % time_options.headerFrequency, 0))
      M.PRINT("\n Current time step: %2.6e s.\n",
              TimeIntegrator.deltaTime)
      M.PRINT(" Min Flow Temp: %11.6f K. Max Flow Temp: %11.6f K.\n",
              Flow.minTemperature, Flow.maxTemperature)
      if particles_options.modeParticles then
        M.PRINT(" Current number of particles: %d.\n", Particles.number)
      end
      M.PRINT("\n")
      M.PRINT("    Iter     Time(s)   Avg Press    Avg Temp      Avg KE  Particle T\n")
    M.END()
    -- Output the current stats to the console for this iteration
    M.PRINT("%8d %11.6f %11.6f %11.6f %11.6f %11.6f\n",
            TimeIntegrator.timeStep,
            TimeIntegrator.simTime,
            Flow.averagePressure,
            Flow.averageTemperature,
            Flow.averageKineticEnergy,
            Particles.averageTemperature)
  M.END()
end

function IO.WriteFlowRestart()
  -- Check if it is time to output a flow restart file
  M.IF(M.EQ(TimeIntegrator.timeStep:get() % time_options.restartEveryTimeSteps, 0))
    -- Write the restart files for density, pressure, and velocity
    fluidGrid:Dump({'rho','pressure','velocity'},
                   io_options.outputFileNamePrefix .. "restart_fluid_%d.hdf",
                   TimeIntegrator.timeStep:get())
  M.END()
end

-- put guards around the particle kernels in case inactive
if particles_options.modeParticles then

  function IO.WriteParticleRestart()
    -- Check if it is time to output a particle restart file
    M.IF(M.EQ(TimeIntegrator.timeStep:get() % time_options.restartEveryTimeSteps, 0))
      -- Write the restart files for position, velocity, temperature and diameter
      particles:Dump({'cell','position','particle_velocity','particle_temperature','diameter'},
                     io_options.outputFileNamePrefix .. "restart_particles_%d.hdf",
                     TimeIntegrator.timeStep:get())
    M.END()
  end

end

function IO.WriteOutput()
  -- Write the console output to the screen
  IO.WriteConsoleOutput()
  -- Write the restart files
  if io_options.wrtRestart then
    -- Write the flow restart files
    IO.WriteFlowRestart()
    -- Write the particle restart files
    if particles_options.modeParticles then
      IO.WriteParticleRestart()
    end
  end
end

-----------------------------------------------------------------------------
--[[                            MAIN EXECUTION                           ]]--
-----------------------------------------------------------------------------

-- Initialize all variables

TimeIntegrator.InitializeVariables()
Flow.IntegrateGeometricQuantities(fluidGrid)
Statistics.ComputeSpatialAverages()
if radiation_options.radiationType == RadiationType.DOM then
  radiationGrid:foreach(Radiation.InitializeCell)
  M.INLINE(radiation.InitModule)
end
IO.WriteOutput()

-- Main iteration loop

M.WHILE(M.AND(M.LT(TimeIntegrator.simTime:get(), time_options.final_time),
              M.LT(TimeIntegrator.timeStep:get(), time_options.max_iter)),
        true)
  TimeIntegrator.CalculateDeltaTime()
  TimeIntegrator.AdvanceTimeStep()
  if not regentlib.config['flow-spmd'] then
    M.IF(M.EQ(TimeIntegrator.timeStep:get() % time_options.consoleFrequency, 0))
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

A.translateAndRun()
