import "regent"

-------------------------------------------------------------------------------
-- IMPORTS
-------------------------------------------------------------------------------

local C = regentlib.c
local MAPPER = terralib.includec("soleil_mapper.h")
local SCHEMA = terralib.includec("config_schema.h")
local UTIL = require 'util-desugared'

local acos = regentlib.acos(double)
local tanh = regentlib.tanh(double)
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
  __xfer_slot : int64;
}

local Particles_primitives = terralib.newlist({
  'position',
  'velocity',
  'temperature',
  'diameter',
  'density',
  '__valid',
})
local Particles_derived = terralib.newlist({
  'cell',
  'deltaVelocityOverRelaxationTime',
  'deltaTemperatureTerm',
})
local Particles_iterTemp = terralib.newlist({
  'position_old',
  'velocity_old',
  'temperature_old',
  'position_new',
  'velocity_new',
  'temperature_new',
})
local Particles_subStepTemp = terralib.newlist({
  'velocity_t',
  'temperature_t',
  '__xfer_dir',
  '__xfer_slot',
})
for _,e in ipairs(Particles_columns.entries) do
  local fld,_ = UTIL.parseStructEntry(e)
  assert(Particles_primitives:find(fld) or
         Particles_derived:find(fld) or
         Particles_iterTemp:find(fld) or
         Particles_subStepTemp:find(fld))
end

local Particles_subStepConserved =
  UTIL.setToList(
    UTIL.setSubList(
      UTIL.listToSet(UTIL.fieldNames(Particles_columns)),
      Particles_subStepTemp))

local TradeQueue_columns =
  UTIL.deriveStruct('TradeQueue_columns',
                    Particles_columns,
                    Particles_subStepConserved)

local CopyQueue_columns =
  UTIL.deriveStruct('CopyQueue_columns',
                    Particles_columns,
                    Particles_primitives)

local struct Fluid_columns {
  rho : double;
  pressure : double;
  velocity : double[3];
  centerCoordinates : double[3];
  cellWidth : double[3];
  velocityGradientX : double[3];
  velocityGradientY : double[3];
  velocityGradientZ : double[3];
  temperature : double;
  temperatureGradient : double[3];
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
  to_Radiation : int3d;
  dudtBoundary : double;
  dTdtBoundary : double;
  velocity_old_NSCBC : double[3];
  temperature_old_NSCBC : double;
  velocity_inc : double[3];
  temperature_inc : double;
}

local Fluid_primitives = terralib.newlist({
  'centerCoordinates',
  'cellWidth',
  'rho',
  'pressure',
  'velocity',
  'temperature',
})

struct Radiation_columns {
  centerCoordinates : double[3];
  cellWidth : double[3];
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

local HDF_FLUID = (require 'hdf_helper')(int3d, int3d, Fluid_columns,
                                         Fluid_primitives,
                                         {timeStep=int,simTime=double})

local HDF_PARTICLES = (require 'hdf_helper')(int1d, int3d, Particles_columns,
                                             Particles_primitives,
                                             {timeStep=int,simTime=double})

-------------------------------------------------------------------------------
-- CONSTANTS
-------------------------------------------------------------------------------

local PI = 3.1415926535898
local SB = 5.67e-08
local DBL_DECIMAL_DIG = 17 -- HACK: normally defined in float.h
local DBL_FORMAT = '%.'..tostring(DBL_DECIMAL_DIG)..'e'

local RK_MIN_ORDER = 2
local RK_MAX_ORDER = 4
-- We only support methods with C[i+1] = A[i+1,i] and A[i,j] = 0 for i != j+1
local RK_B = { -- B[1] B[2] ... B[s]
  [2] = {    0.0,     1.0},
  [3] = {1.0/4.0,     0.0, 3.0/4.0},
  [4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0},
}
local RK_C = { -- C[2] ... C[s]
  [2] = {1.0/2.0},
  [3] = {1.0/3.0, 2.0/3.0},
  [4] = {1.0/2.0, 1.0/2.0,     1.0},
}

-------------------------------------------------------------------------------
-- MACROS
-------------------------------------------------------------------------------

__demand(__inline)
task is_xNegGhost(c : int3d, Grid_xBnum : int)
  return c.x < Grid_xBnum
end

__demand(__inline)
task is_yNegGhost(c : int3d, Grid_yBnum : int)
  return c.y < Grid_yBnum
end

__demand(__inline)
task is_zNegGhost(c : int3d, Grid_zBnum : int)
  return c.z < Grid_zBnum
end

__demand(__inline)
task is_xPosGhost(c : int3d, Grid_xBnum : int, Grid_xNum : int)
  return c.x >= Grid_xNum + Grid_xBnum
end

__demand(__inline)
task is_yPosGhost(c : int3d, Grid_yBnum : int, Grid_yNum : int)
  return c.y >= Grid_yNum + Grid_yBnum
end

__demand(__inline)
task is_zPosGhost(c : int3d, Grid_zBnum : int, Grid_zNum : int)
  return c.z >= Grid_zNum + Grid_zBnum
end

__demand(__inline)
task in_interior(c : int3d,
                 Grid_xBnum : int, Grid_xNum : int,
                 Grid_yBnum : int, Grid_yNum : int,
                 Grid_zBnum : int, Grid_zNum : int)
  return
    Grid_xBnum <= c.x and c.x < Grid_xNum + Grid_xBnum and
    Grid_yBnum <= c.y and c.y < Grid_yNum + Grid_yBnum and
    Grid_zBnum <= c.z and c.z < Grid_zNum + Grid_zBnum
end

__demand(__inline)
task drand48_r(rngState : &C.drand48_data)
  var res : double
  C.drand48_r(rngState, &res)
  return res
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

__demand(__leaf) -- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task Console_WriteHeader(_ : int,
                         config : Config)
  [emitConsoleWrite(config, 'Iteration\t'..
                            'Sim Time\t'..
                            'Wall Time\t'..
                            'Delta Time\t'..
                            'Avg Press\t'..
                            'Avg Temp\t'..
                            'Avg KE\t'..
                            'Particle Num\t'..
                            'Avg Particle T\n')];
  return _
end

__demand(__leaf) -- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task Console_Write(config : Config,
                   Integrator_timeStep : int,
                   Integrator_simTime : double,
                   startTime : uint64,
                   Integrator_deltaTime : double,
                   Flow_averagePressure : double,
                   Flow_averageTemperature : double,
                   Flow_averageKineticEnergy : double,
                   Particles_number : int64,
                   Particles_averageTemperature : double)
  var currTime = C.legion_get_current_time_in_micros() / 1000;
  [emitConsoleWrite(config, '%d\t'..
                            DBL_FORMAT..'\t'..
                            '%llu.%03llu\t'..
                            DBL_FORMAT..'\t'..
                            DBL_FORMAT..'\t'..
                            DBL_FORMAT..'\t'..
                            DBL_FORMAT..'\t'..
                            '%lld\t'..
                            DBL_FORMAT..'\n',
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

-- regentlib.rexpr, regentlib.rexpr, regentlib.rexpr, regentlib.rexpr*
--   -> regentlib.rquote
local function emitProbeWrite(config, probeId, format, ...)
  local args = terralib.newlist{...}
  return rquote
    var filename = [&int8](C.malloc(256))
    C.snprintf(filename, 256, '%s/probe%d.csv', config.Mapping.outDir, probeId)
    var file = UTIL.openFile(filename, 'a')
    C.free(filename)
    C.fprintf(file, format, [args])
    C.fflush(file)
    C.fclose(file)
  end
end

__demand(__leaf, __parallel, __cuda)
task Probe_AvgFluidT(Fluid : region(ispace(int3d), Fluid_columns),
                     probe : SCHEMA.Volume,
                     totalCells : int)
where
  reads(Fluid.temperature)
do
  var fromCell = probe.fromCell
  var uptoCell = probe.uptoCell
  var acc = 0.0
  __demand(__openmp)
  for c in Fluid do
    if fromCell[0] <= c.x and c.x <= uptoCell[0] and
       fromCell[1] <= c.y and c.y <= uptoCell[1] and
       fromCell[2] <= c.z and c.z <= uptoCell[2] then
      acc += Fluid[c].temperature / totalCells
    end
  end
  return acc
end

__demand(__leaf, __parallel, __cuda)
task Probe_CountParticles(Particles : region(ispace(int1d), Particles_columns),
                          probe : SCHEMA.Volume)
where
  reads(Particles.{cell, __valid})
do
  var fromCell = probe.fromCell
  var uptoCell = probe.uptoCell
  var acc = int64(0)
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      var cell = Particles[p].cell
      if fromCell[0] <= cell.x and cell.x <= uptoCell[0] and
         fromCell[1] <= cell.y and cell.y <= uptoCell[1] and
         fromCell[2] <= cell.z and cell.z <= uptoCell[2] then
        acc += 1
      end
    end
  end
  return acc
end

__demand(__leaf, __parallel, __cuda)
task Probe_AvgParticleT(Particles : region(ispace(int1d), Particles_columns),
                        probe : SCHEMA.Volume,
                        totalParticles : int64)
where
  reads(Particles.{cell, temperature, __valid})
do
  var fromCell = probe.fromCell
  var uptoCell = probe.uptoCell
  var acc = 0.0
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      var cell = Particles[p].cell
      if fromCell[0] <= cell.x and cell.x <= uptoCell[0] and
         fromCell[1] <= cell.y and cell.y <= uptoCell[1] and
         fromCell[2] <= cell.z and cell.z <= uptoCell[2] then
        acc += Particles[p].temperature / totalParticles
      end
    end
  end
  return acc
end

__demand(__leaf, __parallel, __cuda)
task Probe_AvgCellOfParticleT(Fluid : region(ispace(int3d), Fluid_columns),
                              Particles : region(ispace(int1d), Particles_columns),
                              probe : SCHEMA.Volume,
                              totalParticles : int64)
where
  reads(Particles.{cell, temperature, __valid}),
  reads(Fluid.temperature)
do
  var fromCell = probe.fromCell
  var uptoCell = probe.uptoCell
  var acc = 0.0
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      var cell = Particles[p].cell
      if fromCell[0] <= cell.x and cell.x <= uptoCell[0] and
         fromCell[1] <= cell.y and cell.y <= uptoCell[1] and
         fromCell[2] <= cell.z and cell.z <= uptoCell[2] then
        acc += Fluid[cell].temperature / totalParticles
      end
    end
  end
  return acc
end

__demand(__leaf) -- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task Probe_WriteHeader(_ : int,
                       config : Config,
                       probeId : int)
  [emitProbeWrite(config, probeId, 'Iter\t'..
                                   'AvgFluidT\t'..
                                   'AvgParticleT\t'..
                                   'AvgCellOfParticleT\n')];
  return _
end

__demand(__leaf) -- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task Probe_Write(_ : int,
                 config : Config,
                 probeId : int,
                 Integrator_timeStep : int,
                 avgFluidT : double,
                 avgParticleT : double,
                 avgCellOfParticleT : double)
  [emitProbeWrite(config, probeId, '%d\t'..
                                   DBL_FORMAT..'\t'..
                                   DBL_FORMAT..'\t'..
                                   DBL_FORMAT..'\n',
                  Integrator_timeStep,
                  avgFluidT,
                  avgParticleT,
                  avgCellOfParticleT)];
  return _
end

__demand(__leaf) -- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task IO_CreateDir(_ : int,
                  dirname : regentlib.string)
  UTIL.createDir(dirname)
  return _
end

-------------------------------------------------------------------------------
-- MESH ROUTINES
-------------------------------------------------------------------------------
-- Description:
--     Linear interpolation, given the line defined by the points
--     (x=alpha, y=a) and (x=beta, y=b) find the y location of the
--     point on the line (x=xi, y=?)
-- Input:
--     xi = location on x axis
--     alpha = lower point on x axis
--     beta =  upper point on x axis
--     a = lower point on y axis
--     b = upper point on y axis
-- Output:
--     y location on line at x=xi
__demand(__inline)
task linear_interpolation(x     : double,
                          alpha : double,
                          beta  : double,
                          a     : double,
                          b     : double) : double
  return (b-a)/(beta-alpha)*(x-alpha) + a
end

-- Description:
--     Generate the cell width of a nonuniform mesh
-- Input:
--     x_min = domain minimum
--     x_max = domain maximum
--     Nx = number of cells between x_min and x_max
-- Output:
--     width of the non-uniform mesh cell
__demand(__inline)
task uniform_cell_width(x_min : double,
                        x_max : double,
                        Nx    : uint64) : double
  return (x_max-x_min)/Nx
end

-- Description:
--     Generate the cell center on a uniform mesh
-- Input:
--     x_min = domain minimum
--     x_max = domain maximum
--     Nx = number of cells between x_min and x_max
--     i  = cell index between x_min and x_max
--         Note: i = 0 has x_min as left face and
--               i = Nx-1 has x_max as right face
-- Output:
--     location of cell center
__demand(__inline)
task uniform_cell_center(x_min : double,
                         x_max : double,
                         Nx    : uint64,
                         i     : uint64) : double
  var dx = uniform_cell_width(x_min, x_max, Nx)
  return x_min + i*dx + dx/2.0
end

-- Description:
--     Generate the location of the face in the negative direction
-- Input:
--     x_min = domain minimum
--     x_max = domain maximum
--     Nx = number of cells between x_min and x_max
--     i  = cell index between x_min and x_max
--         Note: i = 0 has x_min as negative direction (left) face and
--               i = Nx-1 has x_max as positive direction (right) face
-- Output:
--     location of face in the negative direction
__demand(__inline)
task uniform_cell_neg_face(x_min : double,
                           x_max : double,
                           Nx    : uint64,
                           i     : uint64) : double
  var dx = uniform_cell_width(x_min, x_max, Nx)
  var x_center = uniform_cell_center(x_min, x_max, Nx, i)
  return x_center - 0.5*dx
end

-- Description:
--     Generate the location of the face in the postive direction
-- Input:
--     x_min = domain minimum
--     x_max = domain maximum
--     Nx = number of cells between x_min and x_max
--     i  = cell index between x_min and x_max
--         Note: i = 0 has x_min as negative direction (left) face and
--               i = Nx-1 has x_max as positive direction (right) face
-- Output:
--     location of face in the postive direction
__demand(__inline)
task uniform_cell_pos_face(x_min : double,
                           x_max : double,
                           Nx    : uint64,
                           i     : uint64) : double
  var dx = uniform_cell_width(x_min, x_max, Nx)
  var x_center = uniform_cell_center(x_min, x_max, Nx, i)
  return x_center + 0.5*dx
end

-- Description:
--     non-linear map point (x) on the interval (x_min, x_max) using
--     a cosine
-- Input:
--     x = location on uniform mesh
--     x_min = domain minimum
--     x_max = domain maximum
-- Output:
--     x location on a non-uniform mesh
__demand(__inline)
task transform_uniform_to_nonuniform(x     : double,
                                     x_min : double,
                                     x_max : double) : double

  -------------------------
  ---- Cosine Stretching --
  -------------------------
  ---- map x onto the interval -1 to 1
  --var x_scaled_minus1_to_plus1 = linear_interpolation(x, x_min, x_max, -1.0, 1.0)
  ---- map non-uniformly onto the interval -1 to 1
  --var x_non_uniform_minus1_to_plus1 = -1.0*cos(PI*(x_scaled_minus1_to_plus1+1.0)/2.0)
  ---- map non-uniform sample back to origional interval x_min to x_max
  --var x_non_uniform = linear_interpolation(x_non_uniform_minus1_to_plus1, -1.0, 1.0, x_min, x_max)

  ----------------------------------
  -- Hyperbolic Tanget Stretching --
  ----------------------------------
  var stretching_factor = 1.7
  -- map x onto the interval 0 to 1
  var x_scaled_zero_to_plus1 = linear_interpolation(x, x_min, x_max, 0.0, 1.0)

  -- map non-uniformly onto the interval 0 to 1
  var x_transformed     = tanh(stretching_factor*(2.0*x_scaled_zero_to_plus1-1.0))/tanh(stretching_factor)
  var x_transformed_min = tanh(stretching_factor*(2.0*0.0-1.0))/tanh(stretching_factor)
  var x_transformed_max = tanh(stretching_factor*(2.0*1.0-1.0))/tanh(stretching_factor)

  -- map non-uniform sample back to origional interval x_min to x_max
  var x_non_uniform = linear_interpolation(x_transformed, x_transformed_min, x_transformed_max, x_min, x_max)


  return x_non_uniform

end

-- Description:
--     Generate the location of the face in the negative direction
--     on a non uniform mesh
-- Input:
--     x_min = domain minimum
--     x_max = domain maximum
--     Nx = number of cells between x_min and x_max
--     i  = cell index between x_min and x_max
--         Note: i = 0 has x_min as negative direction (left) face and
--               i = Nx-1 has x_max as positive direction (right) face
--               so no ghost cells accounted for here
-- Output:
--     location of face in the negative direction
__demand(__inline)
task nonuniform_cell_neg_face(x_min : double,
                              x_max : double,
                              Nx    : uint64,
                              i     : uint64) : double
  var x_uniform_neg_face = uniform_cell_neg_face(x_min, x_max, Nx, i)
  return transform_uniform_to_nonuniform(x_uniform_neg_face, x_min, x_max)
end

-- Description:
--     Generate the location of the face in the postive direction
--     on a non uniform mesh
-- Input:
--     x_min = domain minimum
--     x_max = domain maximum
--     Nx = number of cells between x_min and x_max
--     i  = cell index between x_min and x_max
--         Note: i = 0 has x_min as negative direction (left) face and
--               i = Nx-1 has x_max as positive direction (right) face
--               so no ghost cells accounted for here
-- Output:
--     location of face in the postive direction
__demand(__inline)
task nonuniform_cell_pos_face(x_min : double,
                              x_max : double,
                              Nx    : uint64,
                              i     : uint64) : double
  var x_uniform_pos_face = uniform_cell_pos_face(x_min, x_max, Nx, i)
  return transform_uniform_to_nonuniform(x_uniform_pos_face, x_min, x_max)
end

-- Description:
--     Generate the cell center of a nonuniform mesh
-- Input:
--     x_min = domain minimum
--     x_max = domain maximum
--     Nx = number of cells between x_min and x_max
--     i  = cell index between x_min and x_max
--         Note: i = 0 has x_min as left face and
--               i = Nx-1 has x_max as right face
--               so no ghost cells accounted for here
-- Output:
--     x location on a non-uniform mesh
__demand(__inline)
task nonuniform_cell_center(x_min : double,
                            x_max : double,
                            Nx    : uint64,
                            i     : uint64) : double
  var x_non_uniform_neg_face = nonuniform_cell_neg_face(x_min, x_max, Nx, i)
  var x_non_uniform_pos_face = nonuniform_cell_pos_face(x_min, x_max, Nx, i)

  var x_non_uniform_center = 0.5*(x_non_uniform_neg_face + x_non_uniform_pos_face)

  return x_non_uniform_center
end

-- Description:
--     Generate the cell width of a nonuniform mesh
-- Input:
--     x_min = domain minimum
--     x_max = domain maximum
--     Nx = number of cells between x_min and x_max
--     i  = cell index between x_min and x_max
--         Note: i = 0 has x_min as left face and
--               i = Nx-1 has x_max as right face
--               so no ghost cells accounted for here
-- Output:
--     width of the non-uniform mesh cell
__demand(__inline)
task nonuniform_cell_width(x_min : double,
                           x_max : double,
                           Nx    : uint64,
                           i     : uint64) : double
  var x_non_uniform_neg_face = nonuniform_cell_neg_face(x_min, x_max, Nx, i)
  var x_non_uniform_pos_face = nonuniform_cell_pos_face(x_min, x_max, Nx, i)

  var x_non_uniform_dx = x_non_uniform_pos_face - x_non_uniform_neg_face

  return x_non_uniform_dx
end

-------------------------------------------------------------------------------
-- OTHER ROUTINES
-------------------------------------------------------------------------------

__demand(__inline)
task locate(pos : double[3],
            Grid_xType : SCHEMA.GridType, Grid_yType : SCHEMA.GridType, Grid_zType : SCHEMA.GridType,
            Grid_xBnum : int32, Grid_xNum : int32, Grid_xOrigin : double, Grid_xWidth : double,
            Grid_yBnum : int32, Grid_yNum : int32, Grid_yOrigin : double, Grid_yWidth : double,
            Grid_zBnum : int32, Grid_zNum : int32, Grid_zOrigin : double, Grid_zWidth : double)

  var xidx = 0
  if (Grid_xType == SCHEMA.GridType_Uniform) then
    var xcw = Grid_xWidth/Grid_xNum
    var xro = Grid_xOrigin-Grid_xBnum*xcw
    var xpos = int(floor((pos[0]-xro)/xcw))
    var xrnum = Grid_xNum+2*Grid_xBnum
    xidx = max(0, min(xrnum-1, xpos))
  elseif (Grid_xType == SCHEMA.GridType_Stretched) then
    var x_neg_boundary = Grid_xOrigin
    var x_pos_boundary = Grid_xOrigin + Grid_xWidth
    -- assumes only one ghost cell in each direction
    if pos[0] < x_neg_boundary then
      xidx = 0
    elseif pos[0] > x_pos_boundary then
      xidx = (Grid_xNum + 2*Grid_xBnum) - 1
    else
      -- brute force search in interior cells
      for xidx_interior=0,(Grid_xNum+Grid_xBnum) do
        if (pos[0]>nonuniform_cell_neg_face(x_neg_boundary, x_pos_boundary, Grid_xNum, xidx_interior)) and
           (pos[0]<nonuniform_cell_pos_face(x_neg_boundary, x_pos_boundary, Grid_xNum, xidx_interior)) then
          -- Convert from interior index back to global index (accounts for ghost cells)
          xidx = (xidx_interior + Grid_xBnum)
          break
        end
      end
    end
  end

  var yidx = 0
  if (Grid_yType == SCHEMA.GridType_Uniform) then
    var ycw = Grid_yWidth/Grid_yNum
    var yro = Grid_yOrigin-Grid_yBnum*ycw
    var ypos = int(floor((pos[1]-yro)/ycw))
    var yrnum = Grid_yNum+2*Grid_yBnum
    yidx = max(0, min(yrnum-1, ypos))
  elseif (Grid_yType == SCHEMA.GridType_Stretched) then
    var y_neg_boundary = Grid_yOrigin
    var y_pos_boundary = Grid_yOrigin + Grid_yWidth
    -- assumes only one ghost cell in each direction
    if pos[1] < y_neg_boundary then
      yidx = 0
    elseif pos[1] > y_pos_boundary then
      yidx = (Grid_yNum + 2*Grid_yBnum) - 1
    else
      -- brute force search in interior cells
      for yidx_interior=0,(Grid_yNum+Grid_yBnum) do
        if (pos[1]>nonuniform_cell_neg_face(y_neg_boundary, y_pos_boundary, Grid_yNum, yidx_interior)) and
           (pos[1]<nonuniform_cell_pos_face(y_neg_boundary, y_pos_boundary, Grid_yNum, yidx_interior)) then
          -- Convert from interior index back to global index (accounts for ghost cells)
          yidx = (yidx_interior + Grid_yBnum)
          break
        end
      end
    end
  end

  var zidx = 0
  if (Grid_zType == SCHEMA.GridType_Uniform) then
    var zcw = Grid_zWidth/Grid_zNum
    var zro = Grid_zOrigin-Grid_zBnum*zcw
    var zpos = int(floor((pos[2]-zro)/zcw))
    var zrnum = Grid_zNum+2*Grid_zBnum
    zidx = max(0, min(zrnum-1, zpos))
  elseif (Grid_zType == SCHEMA.GridType_Stretched) then
    var z_neg_boundary = Grid_zOrigin
    var z_pos_boundary = Grid_zOrigin + Grid_zWidth
    -- assumes only one ghost cell in each direction
    if pos[2] < z_neg_boundary then
      zidx = 0
    elseif pos[2] > z_pos_boundary then
      zidx = (Grid_zNum + 2*Grid_zBnum) - 1
    else
      -- brute force search in interior cells
      for zidx_interior=0,(Grid_zNum+Grid_zBnum) do
        if (pos[2]>nonuniform_cell_neg_face(z_neg_boundary, z_pos_boundary, Grid_zNum, zidx_interior)) and
           (pos[2]<nonuniform_cell_pos_face(z_neg_boundary, z_pos_boundary, Grid_zNum, zidx_interior)) then
          -- Convert from interior index back to global index (accounts for ghost cells)
          zidx = (zidx_interior + Grid_zBnum)
          break
        end
      end
    end
  end

  return int3d{xidx, yidx, zidx}
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
                                xyz000 : double[3],
                                xyz100 : double[3],
                                xyz010 : double[3],
                                xyz110 : double[3],
                                xyz001 : double[3],
                                xyz101 : double[3],
                                xyz011 : double[3],
                                xyz111 : double[3])

  var dX = (xyz[0] - xyz000[0])/(xyz100[0] - xyz000[0])
  var dY = (xyz[1] - xyz000[1])/(xyz010[1] - xyz000[1])
  var dZ = (xyz[2] - xyz000[2])/(xyz001[2] - xyz000[2])

  var oneMinusdX = (1.0-dX)
  var oneMinusdY = (1.0-dY)
  var oneMinusdZ = (1.0-dZ)

  var weight00 = vv_add(vs_mul(c000, oneMinusdX), vs_mul(c100, dX))
  var weight10 = vv_add(vs_mul(c010, oneMinusdX), vs_mul(c110, dX))
  var weight01 = vv_add(vs_mul(c001, oneMinusdX), vs_mul(c101, dX))
  var weight11 = vv_add(vs_mul(c011, oneMinusdX), vs_mul(c111, dX))

  var weight0 =  vv_add(vs_mul(weight00, oneMinusdY), vs_mul(weight10, dY))
  var weight1 =  vv_add(vs_mul(weight01, oneMinusdY), vs_mul(weight11, dY))

  return vv_add(vs_mul(weight0, oneMinusdZ), vs_mul(weight1, dZ))
end

__demand(__inline)
task InterpolateTriVelocity(c : int3d,
                            xyz : double[3],
                            Fluid : region(ispace(int3d), Fluid_columns))
where
  reads(Fluid.{centerCoordinates, velocity})
do

  var velocity000 = Fluid[c].velocity
  var velocity001 = Fluid[(c+{ 0, 0, 1})%Fluid.bounds].velocity
  var velocity00_ = Fluid[(c+{ 0, 0,-1})%Fluid.bounds].velocity
  var velocity010 = Fluid[(c+{ 0, 1, 0})%Fluid.bounds].velocity
  var velocity011 = Fluid[(c+{ 0, 1, 1})%Fluid.bounds].velocity
  var velocity01_ = Fluid[(c+{ 0, 1,-1})%Fluid.bounds].velocity
  var velocity0_0 = Fluid[(c+{ 0,-1, 0})%Fluid.bounds].velocity
  var velocity0_1 = Fluid[(c+{ 0,-1, 1})%Fluid.bounds].velocity
  var velocity0__ = Fluid[(c+{ 0,-1,-1})%Fluid.bounds].velocity

  var velocity100 = Fluid[(c+{ 1, 0, 0})%Fluid.bounds].velocity
  var velocity101 = Fluid[(c+{ 1, 0, 1})%Fluid.bounds].velocity
  var velocity10_ = Fluid[(c+{ 1, 0,-1})%Fluid.bounds].velocity
  var velocity110 = Fluid[(c+{ 1, 1, 0})%Fluid.bounds].velocity
  var velocity111 = Fluid[(c+{ 1, 1, 1})%Fluid.bounds].velocity
  var velocity11_ = Fluid[(c+{ 1, 1,-1})%Fluid.bounds].velocity
  var velocity1_0 = Fluid[(c+{ 1,-1, 0})%Fluid.bounds].velocity
  var velocity1_1 = Fluid[(c+{ 1,-1, 1})%Fluid.bounds].velocity
  var velocity1__ = Fluid[(c+{ 1,-1,-1})%Fluid.bounds].velocity

  var velocity_00 = Fluid[(c+{-1, 0, 0})%Fluid.bounds].velocity
  var velocity_01 = Fluid[(c+{-1, 0, 1})%Fluid.bounds].velocity
  var velocity_0_ = Fluid[(c+{-1, 0,-1})%Fluid.bounds].velocity
  var velocity_10 = Fluid[(c+{-1, 1, 0})%Fluid.bounds].velocity
  var velocity_11 = Fluid[(c+{-1, 1, 1})%Fluid.bounds].velocity
  var velocity_1_ = Fluid[(c+{-1, 1,-1})%Fluid.bounds].velocity
  var velocity__0 = Fluid[(c+{-1,-1, 0})%Fluid.bounds].velocity
  var velocity__1 = Fluid[(c+{-1,-1, 1})%Fluid.bounds].velocity
  var velocity___ = Fluid[(c+{-1,-1,-1})%Fluid.bounds].velocity


  var centerCoordinates000 = Fluid[c].centerCoordinates
  var centerCoordinates001 = Fluid[(c+{ 0, 0, 1})%Fluid.bounds].centerCoordinates
  var centerCoordinates00_ = Fluid[(c+{ 0, 0,-1})%Fluid.bounds].centerCoordinates
  var centerCoordinates010 = Fluid[(c+{ 0, 1, 0})%Fluid.bounds].centerCoordinates
  var centerCoordinates011 = Fluid[(c+{ 0, 1, 1})%Fluid.bounds].centerCoordinates
  var centerCoordinates01_ = Fluid[(c+{ 0, 1,-1})%Fluid.bounds].centerCoordinates
  var centerCoordinates0_0 = Fluid[(c+{ 0,-1, 0})%Fluid.bounds].centerCoordinates
  var centerCoordinates0_1 = Fluid[(c+{ 0,-1, 1})%Fluid.bounds].centerCoordinates
  var centerCoordinates0__ = Fluid[(c+{ 0,-1,-1})%Fluid.bounds].centerCoordinates

  var centerCoordinates100 = Fluid[(c+{ 1, 0, 0})%Fluid.bounds].centerCoordinates
  var centerCoordinates101 = Fluid[(c+{ 1, 0, 1})%Fluid.bounds].centerCoordinates
  var centerCoordinates10_ = Fluid[(c+{ 1, 0,-1})%Fluid.bounds].centerCoordinates
  var centerCoordinates110 = Fluid[(c+{ 1, 1, 0})%Fluid.bounds].centerCoordinates
  var centerCoordinates111 = Fluid[(c+{ 1, 1, 1})%Fluid.bounds].centerCoordinates
  var centerCoordinates11_ = Fluid[(c+{ 1, 1,-1})%Fluid.bounds].centerCoordinates
  var centerCoordinates1_0 = Fluid[(c+{ 1,-1, 0})%Fluid.bounds].centerCoordinates
  var centerCoordinates1_1 = Fluid[(c+{ 1,-1, 1})%Fluid.bounds].centerCoordinates
  var centerCoordinates1__ = Fluid[(c+{ 1,-1,-1})%Fluid.bounds].centerCoordinates

  var centerCoordinates_00 = Fluid[(c+{-1, 0, 0})%Fluid.bounds].centerCoordinates
  var centerCoordinates_01 = Fluid[(c+{-1, 0, 1})%Fluid.bounds].centerCoordinates
  var centerCoordinates_0_ = Fluid[(c+{-1, 0,-1})%Fluid.bounds].centerCoordinates
  var centerCoordinates_10 = Fluid[(c+{-1, 1, 0})%Fluid.bounds].centerCoordinates
  var centerCoordinates_11 = Fluid[(c+{-1, 1, 1})%Fluid.bounds].centerCoordinates
  var centerCoordinates_1_ = Fluid[(c+{-1, 1,-1})%Fluid.bounds].centerCoordinates
  var centerCoordinates__0 = Fluid[(c+{-1,-1, 0})%Fluid.bounds].centerCoordinates
  var centerCoordinates__1 = Fluid[(c+{-1,-1, 1})%Fluid.bounds].centerCoordinates
  var centerCoordinates___ = Fluid[(c+{-1,-1,-1})%Fluid.bounds].centerCoordinates

  var v000 = array(0.0, 0.0, 0.0)
  var v001 = array(0.0, 0.0, 0.0)
  var v010 = array(0.0, 0.0, 0.0)
  var v011 = array(0.0, 0.0, 0.0)
  var v100 = array(0.0, 0.0, 0.0)
  var v101 = array(0.0, 0.0, 0.0)
  var v110 = array(0.0, 0.0, 0.0)
  var v111 = array(0.0, 0.0, 0.0)

  var xyz000 = array(0.0, 0.0, 0.0)
  var xyz001 = array(0.0, 0.0, 0.0)
  var xyz010 = array(0.0, 0.0, 0.0)
  var xyz011 = array(0.0, 0.0, 0.0)
  var xyz100 = array(0.0, 0.0, 0.0)
  var xyz101 = array(0.0, 0.0, 0.0)
  var xyz110 = array(0.0, 0.0, 0.0)
  var xyz111 = array(0.0, 0.0, 0.0)

  if (xyz[0]>Fluid[c].centerCoordinates[0]) then
    if (xyz[1]>Fluid[c].centerCoordinates[1]) then
      if (xyz[2]>Fluid[c].centerCoordinates[2]) then
        v000 = velocity000
        v001 = velocity001
        v010 = velocity010
        v011 = velocity011
        v100 = velocity100
        v101 = velocity101
        v110 = velocity110
        v111 = velocity111

        xyz000 = centerCoordinates000
        xyz001 = centerCoordinates001
        xyz010 = centerCoordinates010
        xyz011 = centerCoordinates011
        xyz100 = centerCoordinates100
        xyz101 = centerCoordinates101
        xyz110 = centerCoordinates110
        xyz111 = centerCoordinates111
      else
        v000 = velocity00_
        v001 = velocity000
        v010 = velocity01_
        v011 = velocity010
        v100 = velocity10_
        v101 = velocity100
        v110 = velocity11_
        v111 = velocity110

        xyz000 = centerCoordinates00_
        xyz001 = centerCoordinates000
        xyz010 = centerCoordinates01_
        xyz011 = centerCoordinates010
        xyz100 = centerCoordinates10_
        xyz101 = centerCoordinates100
        xyz110 = centerCoordinates11_
        xyz111 = centerCoordinates110
      end
    else
      if (xyz[2]>Fluid[c].centerCoordinates[2]) then
        v000 = velocity0_0
        v001 = velocity0_1
        v010 = velocity000
        v011 = velocity001
        v100 = velocity1_0
        v101 = velocity1_1
        v110 = velocity100
        v111 = velocity101

        xyz000 = centerCoordinates0_0
        xyz001 = centerCoordinates0_1
        xyz010 = centerCoordinates000
        xyz011 = centerCoordinates001
        xyz100 = centerCoordinates1_0
        xyz101 = centerCoordinates1_1
        xyz110 = centerCoordinates100
        xyz111 = centerCoordinates101
      else
        v000 = velocity0__
        v001 = velocity0_0
        v010 = velocity00_
        v011 = velocity000
        v100 = velocity1__
        v101 = velocity1_0
        v110 = velocity10_
        v111 = velocity100

        xyz000 = centerCoordinates0__
        xyz001 = centerCoordinates0_0
        xyz010 = centerCoordinates00_
        xyz011 = centerCoordinates000
        xyz100 = centerCoordinates1__
        xyz101 = centerCoordinates1_0
        xyz110 = centerCoordinates10_
        xyz111 = centerCoordinates100
      end
    end
  else
    if (xyz[1]>Fluid[c].centerCoordinates[1]) then
      if (xyz[2]>Fluid[c].centerCoordinates[2]) then
        v000 = velocity_00
        v001 = velocity_01
        v010 = velocity_10
        v011 = velocity_11
        v100 = velocity000
        v101 = velocity001
        v110 = velocity010
        v111 = velocity011

        xyz000 = centerCoordinates_00
        xyz001 = centerCoordinates_01
        xyz010 = centerCoordinates_10
        xyz011 = centerCoordinates_11
        xyz100 = centerCoordinates000
        xyz101 = centerCoordinates001
        xyz110 = centerCoordinates010
        xyz111 = centerCoordinates011
      else
        v000 = velocity_0_
        v001 = velocity_00
        v010 = velocity_1_
        v011 = velocity_10
        v100 = velocity00_
        v101 = velocity000
        v110 = velocity01_
        v111 = velocity010

        xyz000 = centerCoordinates_0_
        xyz001 = centerCoordinates_00
        xyz010 = centerCoordinates_1_
        xyz011 = centerCoordinates_10
        xyz100 = centerCoordinates00_
        xyz101 = centerCoordinates000
        xyz110 = centerCoordinates01_
        xyz111 = centerCoordinates010
      end
    else
      if (xyz[2]>Fluid[c].centerCoordinates[2]) then
        v000 = velocity__0
        v001 = velocity__1
        v010 = velocity_00
        v011 = velocity_01
        v100 = velocity0_0
        v101 = velocity0_1
        v110 = velocity000
        v111 = velocity001

        xyz000 = centerCoordinates__0
        xyz001 = centerCoordinates__1
        xyz010 = centerCoordinates_00
        xyz011 = centerCoordinates_01
        xyz100 = centerCoordinates0_0
        xyz101 = centerCoordinates0_1
        xyz110 = centerCoordinates000
        xyz111 = centerCoordinates001
      else
        v000 = velocity___
        v001 = velocity__0
        v010 = velocity_0_
        v011 = velocity_00
        v100 = velocity0__
        v101 = velocity0_0
        v110 = velocity00_
        v111 = velocity000

        xyz000 = centerCoordinates___
        xyz001 = centerCoordinates__0
        xyz010 = centerCoordinates_0_
        xyz011 = centerCoordinates_00
        xyz100 = centerCoordinates0__
        xyz101 = centerCoordinates0_0
        xyz110 = centerCoordinates00_
        xyz111 = centerCoordinates000
      end
    end
  end

  return TrilinearInterpolateVelocity(xyz, v000, v100, v010, v110, v001, v101, v011, v111, xyz000, xyz100, xyz010, xyz110, xyz001, xyz101, xyz011, xyz111)
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
                              xyz000 : double[3],
                              xyz100 : double[3],
                              xyz010 : double[3],
                              xyz110 : double[3],
                              xyz001 : double[3],
                              xyz101 : double[3],
                              xyz011 : double[3],
                              xyz111 : double[3])

  var dX = (xyz[0] - xyz000[0])/(xyz100[0] - xyz000[0])
  var dY = (xyz[1] - xyz000[1])/(xyz010[1] - xyz000[1])
  var dZ = (xyz[2] - xyz000[2])/(xyz001[2] - xyz000[2])

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
                        Fluid : region(ispace(int3d), Fluid_columns))
where
  reads(Fluid.{centerCoordinates, temperature})
do
  var temperature000 = Fluid[c].temperature
  var temperature001 = Fluid[((c+{ 0, 0, 1})%Fluid.bounds)].temperature
  var temperature00_ = Fluid[((c+{ 0, 0,-1})%Fluid.bounds)].temperature
  var temperature010 = Fluid[((c+{ 0, 1, 0})%Fluid.bounds)].temperature
  var temperature011 = Fluid[((c+{ 0, 1, 1})%Fluid.bounds)].temperature
  var temperature01_ = Fluid[((c+{ 0, 1,-1})%Fluid.bounds)].temperature
  var temperature0_0 = Fluid[((c+{ 0,-1, 0})%Fluid.bounds)].temperature
  var temperature0_1 = Fluid[((c+{ 0,-1, 1})%Fluid.bounds)].temperature
  var temperature0__ = Fluid[((c+{ 0,-1,-1})%Fluid.bounds)].temperature

  var temperature100 = Fluid[((c+{ 1, 0, 0})%Fluid.bounds)].temperature
  var temperature101 = Fluid[((c+{ 1, 0, 1})%Fluid.bounds)].temperature
  var temperature10_ = Fluid[((c+{ 1, 0,-1})%Fluid.bounds)].temperature
  var temperature110 = Fluid[((c+{ 1, 1, 0})%Fluid.bounds)].temperature
  var temperature111 = Fluid[((c+{ 1, 1, 1})%Fluid.bounds)].temperature
  var temperature11_ = Fluid[((c+{ 1, 1,-1})%Fluid.bounds)].temperature
  var temperature1_0 = Fluid[((c+{ 1,-1, 0})%Fluid.bounds)].temperature
  var temperature1_1 = Fluid[((c+{ 1,-1, 1})%Fluid.bounds)].temperature
  var temperature1__ = Fluid[((c+{ 1,-1,-1})%Fluid.bounds)].temperature

  var temperature_00 = Fluid[((c+{-1, 0, 0})%Fluid.bounds)].temperature
  var temperature_01 = Fluid[((c+{-1, 0, 1})%Fluid.bounds)].temperature
  var temperature_0_ = Fluid[((c+{-1, 0,-1})%Fluid.bounds)].temperature
  var temperature_10 = Fluid[((c+{-1, 1, 0})%Fluid.bounds)].temperature
  var temperature_11 = Fluid[((c+{-1, 1, 1})%Fluid.bounds)].temperature
  var temperature_1_ = Fluid[((c+{-1, 1,-1})%Fluid.bounds)].temperature
  var temperature__0 = Fluid[((c+{-1,-1, 0})%Fluid.bounds)].temperature
  var temperature__1 = Fluid[((c+{-1,-1, 1})%Fluid.bounds)].temperature
  var temperature___ = Fluid[((c+{-1,-1,-1})%Fluid.bounds)].temperature

  var centerCoordinates000 = Fluid[c].centerCoordinates
  var centerCoordinates001 = Fluid[((c+{ 0, 0, 1})%Fluid.bounds)].centerCoordinates
  var centerCoordinates00_ = Fluid[((c+{ 0, 0,-1})%Fluid.bounds)].centerCoordinates
  var centerCoordinates010 = Fluid[((c+{ 0, 1, 0})%Fluid.bounds)].centerCoordinates
  var centerCoordinates011 = Fluid[((c+{ 0, 1, 1})%Fluid.bounds)].centerCoordinates
  var centerCoordinates01_ = Fluid[((c+{ 0, 1,-1})%Fluid.bounds)].centerCoordinates
  var centerCoordinates0_0 = Fluid[((c+{ 0,-1, 0})%Fluid.bounds)].centerCoordinates
  var centerCoordinates0_1 = Fluid[((c+{ 0,-1, 1})%Fluid.bounds)].centerCoordinates
  var centerCoordinates0__ = Fluid[((c+{ 0,-1,-1})%Fluid.bounds)].centerCoordinates

  var centerCoordinates100 = Fluid[((c+{ 1, 0, 0})%Fluid.bounds)].centerCoordinates
  var centerCoordinates101 = Fluid[((c+{ 1, 0, 1})%Fluid.bounds)].centerCoordinates
  var centerCoordinates10_ = Fluid[((c+{ 1, 0,-1})%Fluid.bounds)].centerCoordinates
  var centerCoordinates110 = Fluid[((c+{ 1, 1, 0})%Fluid.bounds)].centerCoordinates
  var centerCoordinates111 = Fluid[((c+{ 1, 1, 1})%Fluid.bounds)].centerCoordinates
  var centerCoordinates11_ = Fluid[((c+{ 1, 1,-1})%Fluid.bounds)].centerCoordinates
  var centerCoordinates1_0 = Fluid[((c+{ 1,-1, 0})%Fluid.bounds)].centerCoordinates
  var centerCoordinates1_1 = Fluid[((c+{ 1,-1, 1})%Fluid.bounds)].centerCoordinates
  var centerCoordinates1__ = Fluid[((c+{ 1,-1,-1})%Fluid.bounds)].centerCoordinates

  var centerCoordinates_00 = Fluid[((c+{-1, 0, 0})%Fluid.bounds)].centerCoordinates
  var centerCoordinates_01 = Fluid[((c+{-1, 0, 1})%Fluid.bounds)].centerCoordinates
  var centerCoordinates_0_ = Fluid[((c+{-1, 0,-1})%Fluid.bounds)].centerCoordinates
  var centerCoordinates_10 = Fluid[((c+{-1, 1, 0})%Fluid.bounds)].centerCoordinates
  var centerCoordinates_11 = Fluid[((c+{-1, 1, 1})%Fluid.bounds)].centerCoordinates
  var centerCoordinates_1_ = Fluid[((c+{-1, 1,-1})%Fluid.bounds)].centerCoordinates
  var centerCoordinates__0 = Fluid[((c+{-1,-1, 0})%Fluid.bounds)].centerCoordinates
  var centerCoordinates__1 = Fluid[((c+{-1,-1, 1})%Fluid.bounds)].centerCoordinates
  var centerCoordinates___ = Fluid[((c+{-1,-1,-1})%Fluid.bounds)].centerCoordinates

  var v000 = 0.0
  var v001 = 0.0
  var v010 = 0.0
  var v011 = 0.0
  var v100 = 0.0
  var v101 = 0.0
  var v110 = 0.0
  var v111 = 0.0

  var xyz000 = array(0.0, 0.0, 0.0)
  var xyz001 = array(0.0, 0.0, 0.0)
  var xyz010 = array(0.0, 0.0, 0.0)
  var xyz011 = array(0.0, 0.0, 0.0)
  var xyz100 = array(0.0, 0.0, 0.0)
  var xyz101 = array(0.0, 0.0, 0.0)
  var xyz110 = array(0.0, 0.0, 0.0)
  var xyz111 = array(0.0, 0.0, 0.0)

  if (xyz[0]>Fluid[c].centerCoordinates[0]) then
    if (xyz[1]>Fluid[c].centerCoordinates[1]) then
      if (xyz[2]>Fluid[c].centerCoordinates[2]) then
        v000 = temperature000
        v001 = temperature001
        v010 = temperature010
        v011 = temperature011
        v100 = temperature100
        v101 = temperature101
        v110 = temperature110
        v111 = temperature111

        xyz000 = centerCoordinates000
        xyz001 = centerCoordinates001
        xyz010 = centerCoordinates010
        xyz011 = centerCoordinates011
        xyz100 = centerCoordinates100
        xyz101 = centerCoordinates101
        xyz110 = centerCoordinates110
        xyz111 = centerCoordinates111
      else
        v000 = temperature00_
        v001 = temperature000
        v010 = temperature01_
        v011 = temperature010
        v100 = temperature10_
        v101 = temperature100
        v110 = temperature11_
        v111 = temperature110

        xyz000 = centerCoordinates00_
        xyz001 = centerCoordinates000
        xyz010 = centerCoordinates01_
        xyz011 = centerCoordinates010
        xyz100 = centerCoordinates10_
        xyz101 = centerCoordinates100
        xyz110 = centerCoordinates11_
        xyz111 = centerCoordinates110
      end
    else
      if (xyz[2]>Fluid[c].centerCoordinates[2]) then
        v000 = temperature0_0
        v001 = temperature0_1
        v010 = temperature000
        v011 = temperature001
        v100 = temperature1_0
        v101 = temperature1_1
        v110 = temperature100
        v111 = temperature101

        xyz000 = centerCoordinates0_0
        xyz001 = centerCoordinates0_1
        xyz010 = centerCoordinates000
        xyz011 = centerCoordinates001
        xyz100 = centerCoordinates1_0
        xyz101 = centerCoordinates1_1
        xyz110 = centerCoordinates100
        xyz111 = centerCoordinates101
      else
        v000 = temperature0__
        v001 = temperature0_0
        v010 = temperature00_
        v011 = temperature000
        v100 = temperature1__
        v101 = temperature1_0
        v110 = temperature10_
        v111 = temperature100

        xyz000 = centerCoordinates0__
        xyz001 = centerCoordinates0_0
        xyz010 = centerCoordinates00_
        xyz011 = centerCoordinates000
        xyz100 = centerCoordinates1__
        xyz101 = centerCoordinates1_0
        xyz110 = centerCoordinates10_
        xyz111 = centerCoordinates100
      end
    end
  else
    if (xyz[1]>Fluid[c].centerCoordinates[1]) then
      if (xyz[2]>Fluid[c].centerCoordinates[2]) then
        v000 = temperature_00
        v001 = temperature_01
        v010 = temperature_10
        v011 = temperature_11
        v100 = temperature000
        v101 = temperature001
        v110 = temperature010
        v111 = temperature011

        xyz000 = centerCoordinates_00
        xyz001 = centerCoordinates_01
        xyz010 = centerCoordinates_10
        xyz011 = centerCoordinates_11
        xyz100 = centerCoordinates000
        xyz101 = centerCoordinates001
        xyz110 = centerCoordinates010
        xyz111 = centerCoordinates011
      else
        v000 = temperature_0_
        v001 = temperature_00
        v010 = temperature_1_
        v011 = temperature_10
        v100 = temperature00_
        v101 = temperature000
        v110 = temperature01_
        v111 = temperature010

        xyz000 = centerCoordinates_0_
        xyz001 = centerCoordinates_00
        xyz010 = centerCoordinates_1_
        xyz011 = centerCoordinates_10
        xyz100 = centerCoordinates00_
        xyz101 = centerCoordinates000
        xyz110 = centerCoordinates01_
        xyz111 = centerCoordinates010
      end
    else
      if (xyz[2]>Fluid[c].centerCoordinates[2]) then
        v000 = temperature__0
        v001 = temperature__1
        v010 = temperature_00
        v011 = temperature_01
        v100 = temperature0_0
        v101 = temperature0_1
        v110 = temperature000
        v111 = temperature001

        xyz000 = centerCoordinates__0
        xyz001 = centerCoordinates__1
        xyz010 = centerCoordinates_00
        xyz011 = centerCoordinates_01
        xyz100 = centerCoordinates0_0
        xyz101 = centerCoordinates0_1
        xyz110 = centerCoordinates000
        xyz111 = centerCoordinates001
      else
        v000 = temperature___
        v001 = temperature__0
        v010 = temperature_0_
        v011 = temperature_00
        v100 = temperature0__
        v101 = temperature0_0
        v110 = temperature00_
        v111 = temperature000

        xyz000 = centerCoordinates___
        xyz001 = centerCoordinates__0
        xyz010 = centerCoordinates_0_
        xyz011 = centerCoordinates_00
        xyz100 = centerCoordinates0__
        xyz101 = centerCoordinates0_0
        xyz110 = centerCoordinates00_
        xyz111 = centerCoordinates000
      end
    end
  end

  return TrilinearInterpolateTemp(xyz, v000, v100, v010, v110, v001, v101, v011, v111, xyz000, xyz100, xyz010, xyz110, xyz001, xyz101, xyz011, xyz111)
end

__demand(__inline)
task GetDynamicViscosity(temperature : double,
                         Flow_constantVisc : double,
                         Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                         Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                         Flow_viscosityModel : SCHEMA.ViscosityModel)
  var viscosity = 0.0
  if Flow_viscosityModel == SCHEMA.ViscosityModel_Constant then
    viscosity = Flow_constantVisc
  elseif Flow_viscosityModel == SCHEMA.ViscosityModel_PowerLaw then
    viscosity = Flow_powerlawViscRef*pow(temperature/Flow_powerlawTempRef, 0.75)
  else -- Flow_viscosityModel == SCHEMA.ViscosityModel_Sutherland
    viscosity =
      Flow_sutherlandViscRef
      * pow(temperature/Flow_sutherlandTempRef, 1.5)
      * (Flow_sutherlandTempRef + Flow_sutherlandSRef)
      / (temperature + Flow_sutherlandSRef)
  end
  return viscosity
end

-- XXX: This task needs the parallelizer ghost regions to do the interpolation,
-- but the parallelizer can't handle this task currently (even if we move the
-- RNG to a separate task). Therefore, at this point this task will only work
-- on a single tile.
__demand(__leaf) -- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task Particles_InitializeRandom(color : int3d,
                                Particles : region(ispace(int1d), Particles_columns),
                                Fluid : region(ispace(int3d), Fluid_columns),
                                config : Config,
                                Grid_volume : double,
                                Grid_xBnum : int, Grid_yBnum : int, Grid_zBnum : int)
where
  reads(Fluid.{centerCoordinates, cellWidth, velocity}),
  writes(Particles.{__valid, cell, position, velocity, density, temperature, diameter})
do
  -- Grid geometry
  var Grid_xNum = config.Grid.xNum
  var Grid_yNum = config.Grid.yNum
  var Grid_zNum = config.Grid.zNum
  var Grid_xWidth = config.Grid.xWidth
  var Grid_yWidth = config.Grid.yWidth
  var Grid_zWidth = config.Grid.zWidth
  var Grid_xOrigin = config.Grid.origin[0]
  var Grid_yOrigin = config.Grid.origin[1]
  var Grid_zOrigin = config.Grid.origin[2]

  -- Tile geometry
  var Tile_xOrigin = config.Grid.origin[0] + config.Grid.xWidth
  for c in Fluid do
    Tile_xOrigin min= Fluid[c].centerCoordinates[0] - 0.5*Fluid[c].cellWidth[0]
  end
  var Tile_yOrigin = config.Grid.origin[1] + config.Grid.yWidth
  for c in Fluid do
    Tile_yOrigin min= Fluid[c].centerCoordinates[1] - 0.5*Fluid[c].cellWidth[1]
  end
  var Tile_zOrigin = config.Grid.origin[3] + config.Grid.zWidth
  for c in Fluid do
    Tile_zOrigin min= Fluid[c].centerCoordinates[2] - 0.5*Fluid[c].cellWidth[2]
  end
  -- Do not initialize particles in ghost cells
  if Tile_xOrigin < config.Grid.origin[0] then Tile_xOrigin = config.Grid.origin[0] end
  if Tile_yOrigin < config.Grid.origin[1] then Tile_yOrigin = config.Grid.origin[1] end
  if Tile_zOrigin < config.Grid.origin[2] then Tile_zOrigin = config.Grid.origin[2] end

  var Tile_xMax = config.Grid.origin[0]
  for c in Fluid do
    Tile_xMax max= Fluid[c].centerCoordinates[0] + 0.5*Fluid[c].cellWidth[0]
  end
  var Tile_yMax = config.Grid.origin[1]
  for c in Fluid do
    Tile_yMax max= Fluid[c].centerCoordinates[1] + 0.5*Fluid[c].cellWidth[1]
  end
  var Tile_zMax = config.Grid.origin[2]
  for c in Fluid do
    Tile_zMax max= Fluid[c].centerCoordinates[2] + 0.5*Fluid[c].cellWidth[2]
  end
  -- Do not initialize particles in ghost cells
  if Tile_xMax < (config.Grid.origin[0] + config.Grid.xWidth) then Tile_xMax = config.Grid.origin[0] + config.Grid.xWidth end
  if Tile_yMax < (config.Grid.origin[1] + config.Grid.yWidth) then Tile_yMax = config.Grid.origin[1] + config.Grid.yWidth end
  if Tile_zMax < (config.Grid.origin[2] + config.Grid.zWidth) then Tile_zMax = config.Grid.origin[2] + config.Grid.zWidth end

  var Tile_xWidth = Tile_xMax - Tile_xOrigin
  var Tile_yWidth = Tile_yMax - Tile_yOrigin
  var Tile_zWidth = Tile_zMax - Tile_zOrigin

  var Tile_volume = Tile_xWidth*Tile_yWidth*Tile_zWidth 

  -- Particle values
  var pBase = Particles.bounds.lo
  var numTiles = config.Mapping.tiles[0]*config.Mapping.tiles[1]*config.Mapping.tiles[2]
  var particlesPerTile = config.Particles.initNum / config.Particles.parcelSize * (Tile_volume/Grid_volume)
  var Particles_density = config.Particles.density
  var Particles_initTemperature = config.Particles.initTemperature
  var Particles_diameterMean = config.Particles.diameterMean
  -- RNG state
  var rngState : C.drand48_data
  C.srand48_r(C.legion_get_current_time_in_nanos(), &rngState)
  -- Fill loop
  for p in Particles do
    var relIdx = int64(p - pBase)
    if relIdx < particlesPerTile then
      -- Pick a random position within the current tile, ignoring boundary cells
      var rx = 0.0; repeat rx = drand48_r(&rngState) until rx ~= 0.0
      var ry = 0.0; repeat ry = drand48_r(&rngState) until ry ~= 0.0
      var rz = 0.0; repeat rz = drand48_r(&rngState) until rz ~= 0.0
      var pos = array( Tile_xOrigin + Tile_xWidth * rx,
                       Tile_yOrigin + Tile_yWidth * ry,
                       Tile_zOrigin + Tile_zWidth * rz )
      var c = locate(pos,
                     config.Grid.xType, config.Grid.yType, config.Grid.zType,
                     Grid_xBnum, Grid_xNum, Grid_xOrigin, Grid_xWidth,
                     Grid_yBnum, Grid_yNum, Grid_yOrigin, Grid_yWidth,
                     Grid_zBnum, Grid_zNum, Grid_zOrigin, Grid_zWidth)
      var flowVelocity = InterpolateTriVelocity(c,
                                                pos,
                                                Fluid)
      Particles[p].__valid = true
      Particles[p].cell = c
      Particles[p].position = pos
      Particles[p].velocity = flowVelocity
      Particles[p].density = Particles_density
      Particles[p].temperature = Particles_initTemperature
      Particles[p].diameter = Particles_diameterMean
    end
  end
end

__demand(__leaf, __cuda) -- MANUALLY PARALLELIZED
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
  var numTiles = config.Mapping.tiles[0]*config.Mapping.tiles[1]*config.Mapping.tiles[2]
  var particlesPerTile = config.Particles.initNum / config.Particles.parcelSize / numTiles
  var Particles_density = config.Particles.density
  var Particles_initTemperature = config.Particles.initTemperature
  var Particles_diameterMean = config.Particles.diameterMean
  __demand(__openmp)
  for p in Particles do
    var relIdx = int64(p - pBase)
    if relIdx < particlesPerTile then
      Particles[p].__valid = true
      var c = lo + int3d{relIdx%xSize, relIdx/xSize%ySize, relIdx/xSize/ySize%zSize}
      Particles[p].cell = c
      Particles[p].position = Fluid[c].centerCoordinates
      Particles[p].velocity = Fluid[c].velocity
      Particles[p].density = Particles_density
      Particles[p].temperature = Particles_initTemperature
      Particles[p].diameter = Particles_diameterMean
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Particles_AbsorbRadiationAlgebraic(Particles : region(ispace(int1d), Particles_columns),
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
    if Particles[p].__valid then
      var crossSectionArea = PI*pow(Particles[p].diameter,2.0)/4.0
      var mass = PI*pow(Particles[p].diameter,3.0)/6.0*Particles[p].density
      var absorbedRadiationIntensity = absorptivity*intensity*crossSectionArea
      Particles[p].temperature_t += absorbedRadiationIntensity/(mass*heatCapacity)
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Particles_initValidField(Particles : region(ispace(int1d), Particles_columns))
where
  writes(Particles.__valid)
do
  __demand(__openmp)
  for p in Particles do
    Particles[p].__valid = false
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_SetCoarseningField(Fluid : region(ispace(int3d), Fluid_columns),
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
  for c in Fluid do
    if in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum) then
      Fluid[c].to_Radiation = int3d{(int3d(c).x-Grid_xBnum)/xFactor,
                                    (int3d(c).y-Grid_yBnum)/yFactor,
                                    (int3d(c).z-Grid_zBnum)/zFactor}
    else
      Fluid[c].to_Radiation = int3d({uint64(0ULL), uint64(0ULL), uint64(0ULL)})
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_InitializeCell(Fluid : region(ispace(int3d), Fluid_columns))
where
  writes(Fluid.centerCoordinates),
  writes(Fluid.cellWidth),
  writes(Fluid.dissipation),
  writes(Fluid.dissipationFlux),
  writes(Fluid.pressure),
  writes(Fluid.rho),
  writes(Fluid.rhoEnergy),
  writes(Fluid.{rhoEnergyFluxX, rhoEnergyFluxY, rhoEnergyFluxZ}),
  writes(Fluid.rhoEnergy_new),
  writes(Fluid.rhoEnergy_old),
  writes(Fluid.rhoEnergy_t),
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
  writes(Fluid.temperatureGradient),
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
    Fluid[c].cellWidth = array(0.0, 0.0, 0.0)
    Fluid[c].velocityGradientX = array(0.0, 0.0, 0.0)
    Fluid[c].velocityGradientY = array(0.0, 0.0, 0.0)
    Fluid[c].velocityGradientZ = array(0.0, 0.0, 0.0)
    Fluid[c].temperature = 0.0
    Fluid[c].temperatureGradient = array(0.0, 0.0, 0.0)
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

-- NOTE: It is safe to not pass the ghost regions to this task, because we
-- always group ghost cells with their neighboring interior cells.
__demand(__leaf, __cuda) -- MANUALLY PARALLELIZED
task Flow_InitializeGeometry(Fluid : region(ispace(int3d), Fluid_columns),
                             Grid_xType : SCHEMA.GridType, Grid_yType : SCHEMA.GridType, Grid_zType : SCHEMA.GridType,
                             Grid_xBnum : int32, Grid_xNum : int32, Grid_xOrigin : double, Grid_xWidth : double,
                             Grid_yBnum : int32, Grid_yNum : int32, Grid_yOrigin : double, Grid_yWidth : double,
                             Grid_zBnum : int32, Grid_zNum : int32, Grid_zOrigin : double, Grid_zWidth : double)
where
  reads writes(Fluid.{centerCoordinates, cellWidth})
do
  -- Find cell center coordinates and cell width of interior cells
  __demand(__openmp)
  for cell in Fluid do

    var xNegGhost = is_xNegGhost(cell, Grid_xBnum)
    var xPosGhost = is_xPosGhost(cell, Grid_xBnum, Grid_xNum)
    var yNegGhost = is_yNegGhost(cell, Grid_yBnum)
    var yPosGhost = is_yPosGhost(cell, Grid_yBnum, Grid_yNum)
    var zNegGhost = is_zNegGhost(cell, Grid_zBnum)
    var zPosGhost = is_zPosGhost(cell, Grid_zBnum, Grid_zNum)

    if not (xNegGhost or xPosGhost) then
      var x_neg_boundary = Grid_xOrigin
      var x_pos_boundary = Grid_xOrigin + Grid_xWidth
      var x_idx_interior = cell.x - Grid_xBnum
      if (Grid_xType == SCHEMA.GridType_Stretched) then
        Fluid[cell].centerCoordinates[0] = nonuniform_cell_center(x_neg_boundary, x_pos_boundary, Grid_xNum, x_idx_interior)
        Fluid[cell].cellWidth[0]         = nonuniform_cell_width( x_neg_boundary, x_pos_boundary, Grid_xNum, x_idx_interior)
      elseif (Grid_xType == SCHEMA.GridType_Uniform) then
        Fluid[cell].centerCoordinates[0] = uniform_cell_center(x_neg_boundary, x_pos_boundary, Grid_xNum, x_idx_interior)
        Fluid[cell].cellWidth[0]         = uniform_cell_width( x_neg_boundary, x_pos_boundary, Grid_xNum)
      end
    end

    if not (yNegGhost or yPosGhost) then
      var y_neg_boundary = Grid_yOrigin
      var y_pos_boundary = Grid_yOrigin + Grid_yWidth
      var y_idx_interior = cell.y - Grid_yBnum
      if (Grid_yType == SCHEMA.GridType_Stretched) then
        Fluid[cell].centerCoordinates[1] = nonuniform_cell_center(y_neg_boundary, y_pos_boundary, Grid_yNum, y_idx_interior)
        Fluid[cell].cellWidth[1]         = nonuniform_cell_width( y_neg_boundary, y_pos_boundary, Grid_yNum, y_idx_interior)
      elseif (Grid_yType == SCHEMA.GridType_Uniform) then
        Fluid[cell].centerCoordinates[1] = uniform_cell_center(y_neg_boundary, y_pos_boundary, Grid_yNum, y_idx_interior)
        Fluid[cell].cellWidth[1]         = uniform_cell_width( y_neg_boundary, y_pos_boundary, Grid_yNum)
      end
    end

    if not (zNegGhost or zPosGhost) then
      var z_neg_boundary = Grid_zOrigin
      var z_pos_boundary = Grid_zOrigin + Grid_zWidth
      var z_idx_interior = cell.z - Grid_zBnum
      if (Grid_zType == SCHEMA.GridType_Stretched) then
        Fluid[cell].centerCoordinates[2] = nonuniform_cell_center(z_neg_boundary, z_pos_boundary, Grid_zNum, z_idx_interior)
        Fluid[cell].cellWidth[2]         = nonuniform_cell_width( z_neg_boundary, z_pos_boundary, Grid_zNum, z_idx_interior)
      elseif (Grid_zType == SCHEMA.GridType_Uniform) then
        Fluid[cell].centerCoordinates[2] = uniform_cell_center(z_neg_boundary, z_pos_boundary, Grid_zNum, z_idx_interior)
        Fluid[cell].cellWidth[2]         = uniform_cell_width( z_neg_boundary, z_pos_boundary, Grid_zNum)
      end
    end

  end

  -- Find cell center coordinates and cell width in ghost cells using the interior cells geometry
  __demand(__openmp)
  for cell in Fluid do

    var xNegGhost = is_xNegGhost(cell, Grid_xBnum)
    var xPosGhost = is_xPosGhost(cell, Grid_xBnum, Grid_xNum)
    var yNegGhost = is_yNegGhost(cell, Grid_yBnum)
    var yPosGhost = is_yPosGhost(cell, Grid_yBnum, Grid_yNum)
    var zNegGhost = is_zNegGhost(cell, Grid_zBnum)
    var zPosGhost = is_zPosGhost(cell, Grid_zBnum, Grid_zNum)

    if xNegGhost then
      Fluid[cell].centerCoordinates[0] = Fluid[cell+{1,0,0}].centerCoordinates[0] - Fluid[cell+{1,0,0}].cellWidth[0]
      Fluid[cell].cellWidth[0] = Fluid[cell+{1,0,0}].cellWidth[0]
    elseif xPosGhost then
      Fluid[cell].centerCoordinates[0] = Fluid[cell-{1,0,0}].centerCoordinates[0] + Fluid[cell-{1,0,0}].cellWidth[0]
      Fluid[cell].cellWidth[0] = Fluid[cell-{1,0,0}].cellWidth[0]
    end

    if yNegGhost then
      Fluid[cell].centerCoordinates[1] = Fluid[cell+{0,1,0}].centerCoordinates[1] - Fluid[cell+{0,1,0}].cellWidth[1]
      Fluid[cell].cellWidth[1] = Fluid[cell+{0,1,0}].cellWidth[1]
    elseif yPosGhost then
      Fluid[cell].centerCoordinates[1] = Fluid[cell-{0,1,0}].centerCoordinates[1] + Fluid[cell-{0,1,0}].cellWidth[1]
      Fluid[cell].cellWidth[1] = Fluid[cell-{0,1,0}].cellWidth[1]
    end

    if zNegGhost then
      Fluid[cell].centerCoordinates[2] = Fluid[cell+{0,0,1}].centerCoordinates[2] - Fluid[cell+{0,0,1}].cellWidth[2]
      Fluid[cell].cellWidth[2] = Fluid[cell+{0,0,1}].cellWidth[2]
    elseif zPosGhost then
      Fluid[cell].centerCoordinates[2] = Fluid[cell-{0,0,1}].centerCoordinates[2] + Fluid[cell-{0,0,1}].cellWidth[2]
      Fluid[cell].cellWidth[2] = Fluid[cell-{0,0,1}].cellWidth[2]
    end

  end

end

__demand(__leaf, __cuda) -- MANUALLY PARALLELIZED
task Radiation_InitializeGeometry(Radiation : region(ispace(int3d), Radiation_columns),
                                  Fluid     : region(ispace(int3d), Fluid_columns),
                                  Grid_xNum : int32,
                                  Grid_yNum : int32,
                                  Grid_zNum : int32,
                                  Radiation_xNum : int32,
                                  Radiation_yNum : int32,
                                  Radiation_zNum : int32)
where
  reads writes(Radiation.{centerCoordinates, cellWidth}),
  reads (Fluid.{centerCoordinates, cellWidth})
do
  -- How many radiation cells per fluid cells
  var xFactor = (Grid_xNum/Radiation_xNum)
  var yFactor = (Grid_yNum/Radiation_yNum)
  var zFactor = (Grid_zNum/Radiation_zNum)

  __demand(__openmp)
  for rad_cell in Radiation do
    var rad_index = int3d(rad_cell)

    var fluid_index_x_neg = (rad_index.x  )*xFactor
    var fluid_index_x_pos = (rad_index.x+1)*xFactor - 1
    var fluid_index_y_neg = (rad_index.y  )*yFactor
    var fluid_index_y_pos = (rad_index.y+1)*yFactor - 1
    var fluid_index_z_neg = (rad_index.z  )*zFactor
    var fluid_index_z_pos = (rad_index.z+1)*zFactor - 1
    var fluid_index_neg = int3d({fluid_index_x_neg, fluid_index_y_neg, fluid_index_z_neg})
    var fluid_index_pos = int3d({fluid_index_x_pos, fluid_index_y_pos, fluid_index_z_pos})

    var x_neg_face = Fluid[fluid_index_neg].centerCoordinates[0]-Fluid[fluid_index_neg].cellWidth[0]
    var x_pos_face = Fluid[fluid_index_pos].centerCoordinates[0]+Fluid[fluid_index_pos].cellWidth[0]
    var y_neg_face = Fluid[fluid_index_neg].centerCoordinates[1]-Fluid[fluid_index_neg].cellWidth[1]
    var y_pos_face = Fluid[fluid_index_pos].centerCoordinates[1]+Fluid[fluid_index_pos].cellWidth[1]
    var z_neg_face = Fluid[fluid_index_neg].centerCoordinates[2]-Fluid[fluid_index_neg].cellWidth[2]
    var z_pos_face = Fluid[fluid_index_pos].centerCoordinates[2]+Fluid[fluid_index_pos].cellWidth[2]

    Radiation[rad_cell].cellWidth = array(x_pos_face - x_neg_face,
                                          y_pos_face - y_neg_face,
                                          z_pos_face - z_neg_face)

    Radiation[rad_cell].centerCoordinates =  array((x_neg_face + x_pos_face)/2.0,
                                                   (y_neg_face + y_pos_face)/2.0,
                                                   (z_neg_face + z_pos_face)/2.0)

  end

end

__demand(__leaf, __parallel, __cuda)
task Flow_FindRealOriginX(Fluid : region(ispace(int3d), Fluid_columns),
                          Grid_xOrigin : double)
where
  reads(Fluid.{centerCoordinates, cellWidth})
do
  var x_min = Grid_xOrigin
  __demand(__openmp)
  for c in Fluid do
    x_min min= Fluid[c].centerCoordinates[0] - 0.5*Fluid[c].cellWidth[0]
  end
  return x_min
end

__demand(__leaf, __parallel, __cuda)
task Flow_FindRealOriginY(Fluid : region(ispace(int3d), Fluid_columns),
                          Grid_yOrigin : double)
where
  reads(Fluid.{centerCoordinates, cellWidth})
do
  var y_min = Grid_yOrigin
  __demand(__openmp)
  for c in Fluid do
    y_min min= Fluid[c].centerCoordinates[1] - 0.5*Fluid[c].cellWidth[1]
  end
  return y_min
end

__demand(__leaf, __parallel, __cuda)
task Flow_FindRealOriginZ(Fluid : region(ispace(int3d), Fluid_columns),
                          Grid_zOrigin : double)
where
  reads(Fluid.{centerCoordinates, cellWidth})
do
  var z_min = Grid_zOrigin
  __demand(__openmp)
  for c in Fluid do
    z_min min= Fluid[c].centerCoordinates[2] - 0.5*Fluid[c].cellWidth[2]
  end
  return z_min
end

__demand(__leaf, __parallel, __cuda)
task Flow_FindRealMaxX(Fluid : region(ispace(int3d), Fluid_columns),
                       Grid_xMax : double)
where
  reads(Fluid.{centerCoordinates, cellWidth})
do
  var x_max = Grid_xMax
  __demand(__openmp)
  for c in Fluid do
    x_max max= Fluid[c].centerCoordinates[0] + 0.5*Fluid[c].cellWidth[0]
  end
  return x_max
end

__demand(__leaf, __parallel, __cuda)
task Flow_FindRealMaxY(Fluid : region(ispace(int3d), Fluid_columns),
                       Grid_yMax : double)
where
  reads(Fluid.{centerCoordinates, cellWidth})
do
  var y_max = Grid_yMax
  __demand(__openmp)
  for c in Fluid do
    y_max max= Fluid[c].centerCoordinates[1] + 0.5*Fluid[c].cellWidth[1]
  end
  return y_max
end

__demand(__leaf, __parallel, __cuda)
task Flow_FindRealMaxZ(Fluid : region(ispace(int3d), Fluid_columns),
                       Grid_zMax : double)
where
  reads(Fluid.{centerCoordinates, cellWidth})
do
  var z_max = Grid_zMax
  __demand(__openmp)
  for c in Fluid do
    z_max max= Fluid[c].centerCoordinates[2] + 0.5*Fluid[c].cellWidth[2]
  end
  return z_max
end

__demand(__leaf, __parallel, __cuda)
task Flow_FindCellOriginX(Fluid : region(ispace(int3d), Fluid_columns),
                          cell : int[3])
where
  reads(Fluid.{centerCoordinates, cellWidth})
do
  var x_origin = 0.0
  __demand(__openmp)
  for c in Fluid do
    if c.x == cell[0] and c.y == cell[1] and c.z == cell[2] then
      x_origin += Fluid[c].centerCoordinates[0] - 0.5*Fluid[c].cellWidth[0]
    end
  end
  return x_origin
end

__demand(__leaf, __parallel, __cuda)
task Flow_FindCellOriginY(Fluid : region(ispace(int3d), Fluid_columns),
                          cell : int[3])
where
  reads(Fluid.{centerCoordinates, cellWidth})
do
  var y_origin = 0.0
  __demand(__openmp)
  for c in Fluid do
    if c.x == cell[0] and c.y == cell[1] and c.z == cell[2] then
      y_origin += Fluid[c].centerCoordinates[1] - 0.5*Fluid[c].cellWidth[1]
    end
  end
  return y_origin
end

__demand(__leaf, __parallel, __cuda)
task Flow_FindCellOriginZ(Fluid : region(ispace(int3d), Fluid_columns),
                          cell : int[3])
where
  reads(Fluid.{centerCoordinates, cellWidth})
do
  var z_origin = 0.0
  __demand(__openmp)
  for c in Fluid do
    if c.x == cell[0] and c.y == cell[1] and c.z == cell[2] then
      z_origin += Fluid[c].centerCoordinates[2] - 0.5*Fluid[c].cellWidth[2]
    end
  end
  return z_origin
end

__demand(__leaf, __parallel, __cuda)
task Flow_InitializeUniform(Fluid : region(ispace(int3d), Fluid_columns), Flow_initParams : double[6])
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

__demand(__leaf) -- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task Flow_InitializeRandom(Fluid : region(ispace(int3d), Fluid_columns),
                           Flow_initParams : double[6])
where
  writes(Fluid.{rho, pressure, velocity})
do
  var magnitude = Flow_initParams[2]
  var rngState : C.drand48_data
  C.srand48_r(C.legion_get_current_time_in_nanos(), &rngState)
  for c in Fluid do
    Fluid[c].rho = Flow_initParams[0]
    Fluid[c].pressure = Flow_initParams[1]
    Fluid[c].velocity = array(2 * (drand48_r(&rngState) - 0.5) * magnitude,
                              2 * (drand48_r(&rngState) - 0.5) * magnitude,
                              2 * (drand48_r(&rngState) - 0.5) * magnitude)
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_InitializeTaylorGreen2D(Fluid : region(ispace(int3d), Fluid_columns),
                                  Flow_initParams : double[6],
                                  Grid_xBnum : int32, Grid_xNum : int32, Grid_xOrigin : double, Grid_xWidth : double,
                                  Grid_yBnum : int32, Grid_yNum : int32, Grid_yOrigin : double, Grid_yWidth : double,
                                  Grid_zBnum : int32, Grid_zNum : int32, Grid_zOrigin : double, Grid_zWidth : double)
where
  reads(Fluid.centerCoordinates),
  writes(Fluid.{rho, velocity, pressure})
do
  __demand(__openmp)
  for c in Fluid do
    var taylorGreenDensity = Flow_initParams[0]
    var taylorGreenPressure = Flow_initParams[1]
    var taylorGreenVelocity = Flow_initParams[2]
    var xy = Fluid[c].centerCoordinates
    var coorZ = 0
    Fluid[c].rho = taylorGreenDensity
    Fluid[c].velocity = vs_mul([double[3]](array(((sin(xy[0])*cos(xy[1]))*cos(coorZ)), (((-cos(xy[0]))*sin(xy[1]))*cos(coorZ)), 0.0)), taylorGreenVelocity)
    var factorA = (cos((2.0*double(coorZ)))+2.0)
    var factorB = (cos((2.0*xy[0]))+cos((2.0*xy[1])))
    Fluid[c].pressure = (taylorGreenPressure+((((taylorGreenDensity*pow(taylorGreenVelocity, 2.0))/16.0)*factorA)*factorB))
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_InitializeTaylorGreen3D(Fluid : region(ispace(int3d), Fluid_columns),
                                  Flow_initParams : double[6],
                                  Grid_xBnum : int32, Grid_xNum : int32, Grid_xOrigin : double, Grid_xWidth : double,
                                  Grid_yBnum : int32, Grid_yNum : int32, Grid_yOrigin : double, Grid_yWidth : double,
                                  Grid_zBnum : int32, Grid_zNum : int32, Grid_zOrigin : double, Grid_zWidth : double)
where
  reads(Fluid.centerCoordinates),
  writes(Fluid.{rho, velocity, pressure})
do
  __demand(__openmp)
  for c in Fluid do
    var taylorGreenDensity = Flow_initParams[0]
    var taylorGreenPressure = Flow_initParams[1]
    var taylorGreenVelocity = Flow_initParams[2]
    var xy = Fluid[c].centerCoordinates
    Fluid[c].rho = taylorGreenDensity
    Fluid[c].velocity = vs_mul([double[3]](array(((sin(xy[0])*cos(xy[1]))*cos(xy[2])), (((-cos(xy[0]))*sin(xy[1]))*cos(xy[2])), 0.0)), taylorGreenVelocity)
    var factorA = (cos((2.0*xy[2]))+2.0)
    var factorB = (cos((2.0*xy[0]))+cos((2.0*xy[1])))
    Fluid[c].pressure = (taylorGreenPressure+((((taylorGreenDensity*pow(taylorGreenVelocity, 2.0))/16.0)*factorA)*factorB))
  end
end

__demand(__leaf) -- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task Flow_InitializePerturbed(Fluid : region(ispace(int3d), Fluid_columns),
                              Flow_initParams : double[6])
where
  writes(Fluid.{rho, pressure, velocity})
do
  var magnitude = Flow_initParams[5]
  var rngState : C.drand48_data
  C.srand48_r(C.legion_get_current_time_in_nanos(), &rngState)
  for c in Fluid do
    Fluid[c].rho = Flow_initParams[0]
    Fluid[c].pressure = Flow_initParams[1]
    Fluid[c].velocity = array(Flow_initParams[2] + 2 * (drand48_r(&rngState) - 0.5) * magnitude,
                              Flow_initParams[3] + 2 * (drand48_r(&rngState) - 0.5) * magnitude,
                              Flow_initParams[4] + 2 * (drand48_r(&rngState) - 0.5) * magnitude)
  end
end

-- Modify initial rho, pressure, velocity to be consistent with NSCBC inflow outflow condition, also set up stuff for RHS of inflow
-- NOTE: It is safe to not pass the ghost regions to this task, because we
-- always group ghost cells with their neighboring interior cells.
__demand(__leaf, __cuda) -- MANUALLY PARALLELIZED
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
  -- Inflow values
  var BC_xBCLeftHeat_type = config.BC.xBCLeftHeat.type
  var BC_xBCLeftHeat_Constant_temperature = config.BC.xBCLeftHeat.u.Constant.temperature
  var BC_xBCLeftInflowProfile_type = config.BC.xBCLeftInflowProfile.type
  var BC_xBCLeftInflowProfile_Constant_velocity = config.BC.xBCLeftInflowProfile.u.Constant.velocity
  var BC_xBCLeftInflowProfile_Duct_meanVelocity = config.BC.xBCLeftInflowProfile.u.Duct.meanVelocity
  var BC_xBCLeftInflowProfile_Incoming_addedVelocity = config.BC.xBCLeftInflowProfile.u.Incoming.addedVelocity
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = is_xNegGhost(c, Grid_xBnum)
    var xPosGhost = is_xPosGhost(c, Grid_xBnum, Grid_xNum)
    var yNegGhost = is_yNegGhost(c, Grid_yBnum)
    var yPosGhost = is_yPosGhost(c, Grid_yBnum, Grid_yNum)
    var zNegGhost = is_zNegGhost(c, Grid_zBnum)
    var zPosGhost = is_zPosGhost(c, Grid_zBnum, Grid_zNum)

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
        velocity_inc[0] += -BC_xBCLeftInflowProfile_Incoming_addedVelocity
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

__demand(__leaf, __parallel, __cuda)
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
    if in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum) then
      var tmpTemperature = (Fluid[c].pressure/(Flow_gasConstant*Fluid[c].rho))
      var velocity = Fluid[c].velocity
      Fluid[c].rhoVelocity = vs_mul(Fluid[c].velocity, Fluid[c].rho)
      var cv = (Flow_gasConstant/(Flow_gamma-1.0))
      Fluid[c].rhoEnergy = Fluid[c].rho*((cv*tmpTemperature)+(0.5*dot(velocity, velocity)))
    end
  end
end

__demand(__leaf, __parallel, __cuda)
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
    var xNegGhost = is_xNegGhost(c, Grid_xBnum)
    var xPosGhost = is_xPosGhost(c, Grid_xBnum, Grid_xNum)
    var yNegGhost = is_yNegGhost(c, Grid_yBnum)
    var yPosGhost = is_yPosGhost(c, Grid_yBnum, Grid_yNum)
    var zNegGhost = is_zNegGhost(c, Grid_zBnum)
    var zPosGhost = is_zPosGhost(c, Grid_zBnum, Grid_zNum)

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

local function mkFlow_CalculateAverageVelocity(dim)
  local I = dim == 'X' and 0 or
            dim == 'Y' and 1 or
            dim == 'Z' and 2 or
            assert(false)
  local __demand(__leaf, __parallel, __cuda)
  task Flow_CalculateAverageVelocity(Fluid : region(ispace(int3d), Fluid_columns),
                                     Grid_xBnum : int32, Grid_xNum : int32,
                                     Grid_yBnum : int32, Grid_yNum : int32,
                                     Grid_zBnum : int32, Grid_zNum : int32)
  where
    reads(Fluid.{rho, rhoVelocity, cellWidth})
  do
    var acc = 0.0
    __demand(__openmp)
    for c in Fluid do
      if in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum) then
        var cellVolume = Fluid[c].cellWidth[0]*Fluid[c].cellWidth[1]*Fluid[c].cellWidth[2]
        acc += (Fluid[c].rhoVelocity[I]/Fluid[c].rho)*cellVolume
      end
    end
    return acc
  end
  local name = 'Flow_CalculateAverageVelocity'..dim
  Flow_CalculateAverageVelocity:set_name(name)
  Flow_CalculateAverageVelocity:get_primary_variant():get_ast().name[1] = name
  return Flow_CalculateAverageVelocity
end
local Flow_CalculateAverageVelocityX = mkFlow_CalculateAverageVelocity('X')
local Flow_CalculateAverageVelocityY = mkFlow_CalculateAverageVelocity('Y')
local Flow_CalculateAverageVelocityZ = mkFlow_CalculateAverageVelocity('Z')

__demand(__leaf, __parallel, __cuda)
task Flow_AdjustAverageVelocity(Fluid : region(ispace(int3d), Fluid_columns),
                                config : Config,
                                Flow_averageVelocityX : double,
                                Flow_averageVelocityY : double,
                                Flow_averageVelocityZ : double,
                                Grid_xBnum : int32, Grid_xNum : int32,
                                Grid_yBnum : int32, Grid_yNum : int32,
                                Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.rho),
  reads writes(Fluid.rhoVelocity)
do
  var adjustmentX = config.Flow.turbForcing.u.HIT.meanVelocity[0] - Flow_averageVelocityX
  var adjustmentY = config.Flow.turbForcing.u.HIT.meanVelocity[1] - Flow_averageVelocityY
  var adjustmentZ = config.Flow.turbForcing.u.HIT.meanVelocity[2] - Flow_averageVelocityZ
  __demand(__openmp)
  for c in Fluid do
    if in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum) then
      Fluid[c].rhoVelocity[0] += adjustmentX*Fluid[c].rho
      Fluid[c].rhoVelocity[1] += adjustmentY*Fluid[c].rho
      Fluid[c].rhoVelocity[2] += adjustmentZ*Fluid[c].rho
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_UpdateAuxiliaryVelocity(Fluid : region(ispace(int3d), Fluid_columns),
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
  -- Domain Size
  var Grid_xWidth = config.Grid.xWidth
  var Grid_yWidth = config.Grid.yWidth
  var Grid_zWidth = config.Grid.zWidth
  -- Inflow values
  var BC_xBCLeftInflowProfile_type = config.BC.xBCLeftInflowProfile.type
  var BC_xBCLeftInflowProfile_Constant_velocity = config.BC.xBCLeftInflowProfile.u.Constant.velocity
  var BC_xBCLeftInflowProfile_Duct_meanVelocity = config.BC.xBCLeftInflowProfile.u.Duct.meanVelocity
  var BC_xBCLeftInflowProfile_Incoming_addedVelocity = config.BC.xBCLeftInflowProfile.u.Incoming.addedVelocity
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = is_xNegGhost(c, Grid_xBnum)
    var xPosGhost = is_xPosGhost(c, Grid_xBnum, Grid_xNum)
    var yNegGhost = is_yNegGhost(c, Grid_yBnum)
    var yPosGhost = is_yPosGhost(c, Grid_yBnum, Grid_yNum)
    var zNegGhost = is_zNegGhost(c, Grid_zBnum)
    var zPosGhost = is_zPosGhost(c, Grid_zBnum, Grid_zNum)
    var interior = in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if NSCBC_inflow_cell then
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
        var mu = GetDynamicViscosity(Fluid[c].temperature,
                                     Flow_constantVisc,
                                     Flow_powerlawTempRef, Flow_powerlawViscRef,
                                     Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef,
                                     Flow_viscosityModel)
        var Re = Fluid[c].rho*meanVelocity*Grid_yWidth / mu
        var n = -1.7 + 1.8*log(Re)
        velocity[0] = meanVelocity*pow((d/d_max), (1.0/n))
      else -- BC_xBCLeftInflowProfile_type == SCHEMA.InflowProfile_Incoming
        velocity = Fluid[c].velocity_inc
        velocity[0] += BC_xBCLeftInflowProfile_Incoming_addedVelocity
      end
      Fluid[c].velocity = velocity
    end
    if interior or NSCBC_outflow_cell then
      Fluid[c].velocity = vs_div(Fluid[c].rhoVelocity, Fluid[c].rho)
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_UpdateGhostConserved(Fluid : region(ispace(int3d), Fluid_columns),
                               config : Config,
                               Grid_xBnum : int32, Grid_xNum : int32,
                               Grid_yBnum : int32, Grid_yNum : int32,
                               Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{rho, velocity, pressure, temperature}),
  writes(Fluid.{rho, rhoEnergy, rhoVelocity})
do
  var BC_xBCLeft  = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  var Flow_gasConstant = config.Flow.gasConstant
  var Flow_gamma = config.Flow.gamma
  var cv = (Flow_gasConstant/(Flow_gamma-1.0))

  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = is_xNegGhost(c, Grid_xBnum)
    var xPosGhost = is_xPosGhost(c, Grid_xBnum, Grid_xNum)
    var yNegGhost = is_yNegGhost(c, Grid_yBnum)
    var yPosGhost = is_yPosGhost(c, Grid_yBnum, Grid_yNum)
    var zNegGhost = is_zNegGhost(c, Grid_zBnum)
    var zPosGhost = is_zPosGhost(c, Grid_zBnum, Grid_zNum)
    var ghost_cell = (xNegGhost or xPosGhost or
                      yNegGhost or yPosGhost or
                      zNegGhost or zPosGhost )
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if ghost_cell then
      var c_bnd = int3d(c)
      var velocity = Fluid[c_bnd].velocity
      var temperature = Fluid[c_bnd].temperature
      var rho : double
      if NSCBC_inflow_cell or NSCBC_outflow_cell then
        rho = Fluid[c_bnd].rho
      else
        rho = Fluid[c_bnd].pressure/(Flow_gasConstant*temperature)
      end
      Fluid[c_bnd].rho = rho
      Fluid[c_bnd].rhoVelocity = vs_mul(velocity, rho)
      Fluid[c_bnd].rhoEnergy = rho*((cv*temperature)+(0.5*dot(velocity, velocity)))
    end
  end
end

-- NOTE: It is safe to not pass the ghost regions to this task, because we
-- always group ghost cells with their neighboring interior cells.
__demand(__leaf, __cuda) -- MANUALLY PARALLELIZED
task Flow_UpdateGhostVelocity(Fluid : region(ispace(int3d), Fluid_columns),
                              config : Config,
                              BC_xNegVelocity : double[3], BC_xPosVelocity : double[3], BC_xNegSign : double[3], BC_xPosSign : double[3],
                              BC_yNegVelocity : double[3], BC_yPosVelocity : double[3], BC_yNegSign : double[3], BC_yPosSign : double[3],
                              BC_zNegVelocity : double[3], BC_zPosVelocity : double[3], BC_zNegSign : double[3], BC_zPosSign : double[3],
                              Grid_xBnum : int32, Grid_xNum : int32,
                              Grid_yBnum : int32, Grid_yNum : int32,
                              Grid_zBnum : int32, Grid_zNum : int32)
where
  reads writes(Fluid.velocity)
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = is_xNegGhost(c, Grid_xBnum)
    var xPosGhost = is_xPosGhost(c, Grid_xBnum, Grid_xNum)
    var yNegGhost = is_yNegGhost(c, Grid_yBnum)
    var yPosGhost = is_yPosGhost(c, Grid_yBnum, Grid_yNum)
    var zNegGhost = is_zNegGhost(c, Grid_zBnum)
    var zPosGhost = is_zPosGhost(c, Grid_zBnum, Grid_zNum)
    if xNegGhost and not BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow then
      var c_int = ((c+{1, 0, 0})%Fluid.bounds)
      var sign = BC_xNegSign
      var bnd_velocity = BC_xNegVelocity
      Fluid[c].velocity = vv_add(vv_mul(Fluid[c_int].velocity, sign), bnd_velocity)
    end
    if xPosGhost and not BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow then
      var c_int = ((c+{-1, 0, 0})%Fluid.bounds)
      var sign = BC_xPosSign
      var bnd_velocity = BC_xPosVelocity
      Fluid[c].velocity = vv_add(vv_mul(Fluid[c_int].velocity, sign), bnd_velocity)
    end
    if yNegGhost then
      var c_int = ((c+{0, 1, 0})%Fluid.bounds)
      var sign = BC_yNegSign
      var bnd_velocity = BC_yNegVelocity
      Fluid[c].velocity = vv_add(vv_mul(Fluid[c_int].velocity, sign), bnd_velocity)
    end
    if yPosGhost then
      var c_int = ((c+{0, -1, 0})%Fluid.bounds)
      var sign = BC_yPosSign
      var bnd_velocity = BC_yPosVelocity
      Fluid[c].velocity = vv_add(vv_mul(Fluid[c_int].velocity, sign), bnd_velocity)
    end
    if zNegGhost then
      var c_int = ((c+{0, 0, 1})%Fluid.bounds)
      var sign = BC_zNegSign
      var bnd_velocity = BC_zNegVelocity
      Fluid[c].velocity = vv_add(vv_mul(Fluid[c_int].velocity, sign), bnd_velocity)
    end
    if zPosGhost then
      var c_int = ((c+{0, 0, -1})%Fluid.bounds)
      var sign = BC_zPosSign
      var bnd_velocity = BC_zPosVelocity
      Fluid[c].velocity = vv_add(vv_mul(Fluid[c_int].velocity, sign), bnd_velocity)
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_ComputeGradients(Fluid : region(ispace(int3d), Fluid_columns),
                           config : Config,
                           Grid_xBnum : int32, Grid_xNum : int32,
                           Grid_yBnum : int32, Grid_yNum : int32,
                           Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{cellWidth, velocity, temperature, temperatureGradient}),
  writes(Fluid.{velocityGradientX, velocityGradientY, velocityGradientZ, temperatureGradient})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = is_xNegGhost(c, Grid_xBnum)
    var xPosGhost = is_xPosGhost(c, Grid_xBnum, Grid_xNum)
    var yNegGhost = is_yNegGhost(c, Grid_yBnum)
    var yPosGhost = is_yPosGhost(c, Grid_yBnum, Grid_yNum)
    var zNegGhost = is_zNegGhost(c, Grid_zBnum)
    var zPosGhost = is_zPosGhost(c, Grid_zBnum, Grid_zNum)
    var interior = in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    var v100 = Fluid[(c+{ 1,  0,  0}) % Fluid.bounds].velocity
    var v010 = Fluid[(c+{ 0,  1,  0}) % Fluid.bounds].velocity
    var v001 = Fluid[(c+{ 0,  0,  1}) % Fluid.bounds].velocity
    var v_00 = Fluid[(c+{-1,  0,  0}) % Fluid.bounds].velocity
    var v0_0 = Fluid[(c+{ 0, -1,  0}) % Fluid.bounds].velocity
    var v00_ = Fluid[(c+{ 0,  0, -1}) % Fluid.bounds].velocity

    var T100 = Fluid[(c+{ 1,  0,  0}) % Fluid.bounds].temperature
    var T010 = Fluid[(c+{ 0,  1,  0}) % Fluid.bounds].temperature
    var T001 = Fluid[(c+{ 0,  0,  1}) % Fluid.bounds].temperature
    var T_00 = Fluid[(c+{-1,  0,  0}) % Fluid.bounds].temperature
    var T0_0 = Fluid[(c+{ 0, -1,  0}) % Fluid.bounds].temperature
    var T00_ = Fluid[(c+{ 0,  0, -1}) % Fluid.bounds].temperature

    var dx100 = Fluid[(c+{ 1,  0,  0}) % Fluid.bounds].cellWidth[0]
    var dx_00 = Fluid[(c+{-1,  0,  0}) % Fluid.bounds].cellWidth[0]
    var dy010 = Fluid[(c+{ 0,  1,  0}) % Fluid.bounds].cellWidth[1]
    var dy0_0 = Fluid[(c+{ 0, -1,  0}) % Fluid.bounds].cellWidth[1]
    var dz001 = Fluid[(c+{ 0,  0,  1}) % Fluid.bounds].cellWidth[2]
    var dz00_ = Fluid[(c+{ 0,  0, -1}) % Fluid.bounds].cellWidth[2]

    if interior then
      Fluid[c].velocityGradientX = vs_div(vv_sub(v100, v_00), 0.5*dx100 + Fluid[c].cellWidth[0] + 0.5*dx_00)
      Fluid[c].velocityGradientY = vs_div(vv_sub(v010, v0_0), 0.5*dy010 + Fluid[c].cellWidth[1] + 0.5*dy0_0)
      Fluid[c].velocityGradientZ = vs_div(vv_sub(v001, v00_), 0.5*dz001 + Fluid[c].cellWidth[2] + 0.5*dz00_)
      Fluid[c].temperatureGradient[0] = (T100-T_00)/(0.5*dx100 + Fluid[c].cellWidth[0] + 0.5*dx_00)
      Fluid[c].temperatureGradient[1] = (T010-T0_0)/(0.5*dy010 + Fluid[c].cellWidth[1] + 0.5*dy0_0)
      Fluid[c].temperatureGradient[2] = (T001-T00_)/(0.5*dz001 + Fluid[c].cellWidth[2] + 0.5*dz00_)
    end
    if NSCBC_inflow_cell  then
      -- forward one sided difference
      Fluid[c].velocityGradientX = vs_div(vv_sub(v100, Fluid[c].velocity), 0.5*(dx100 + Fluid[c].cellWidth[0]))
      Fluid[c].temperatureGradient[0] = (T100-Fluid[c].temperature)/(0.5*(dx100 + Fluid[c].cellWidth[0]))
      -- central difference
      Fluid[c].velocityGradientY = vs_div(vv_sub(v010, v0_0), 0.5*dy010 + Fluid[c].cellWidth[1] + 0.5*dy0_0)
      Fluid[c].velocityGradientZ = vs_div(vv_sub(v001, v00_), 0.5*dz001 + Fluid[c].cellWidth[2] + 0.5*dz00_)
      Fluid[c].temperatureGradient[1] = (T010-T0_0)/(0.5*dy010 + Fluid[c].cellWidth[1] + 0.5*dy0_0)
      Fluid[c].temperatureGradient[2] = (T001-T00_)/(0.5*dz001 + Fluid[c].cellWidth[2] + 0.5*dz00_)
    end
    if NSCBC_outflow_cell  then
      -- backward one sided difference
      Fluid[c].velocityGradientX = vs_div(vv_sub(Fluid[c].velocity, v_00), 0.5*(Fluid[c].cellWidth[0] + dx_00))
      Fluid[c].temperatureGradient[0] = (Fluid[c].temperature-T_00)/(0.5*(dx_00 + Fluid[c].cellWidth[0]))
      -- central difference
      Fluid[c].velocityGradientY = vs_div(vv_sub(v010, v0_0), 0.5*dy010 + Fluid[c].cellWidth[1] + 0.5*dy0_0)
      Fluid[c].velocityGradientZ = vs_div(vv_sub(v001, v00_), 0.5*dz001 + Fluid[c].cellWidth[2] + 0.5*dz00_)
      Fluid[c].temperatureGradient[1] = (T010-T0_0)/(0.5*dy010 + Fluid[c].cellWidth[1] + 0.5*dy0_0)
      Fluid[c].temperatureGradient[2] = (T001-T00_)/(0.5*dz001 + Fluid[c].cellWidth[2] + 0.5*dz00_)
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_UpdateAuxiliaryThermodynamics(Fluid : region(ispace(int3d), Fluid_columns),
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
    var xNegGhost = is_xNegGhost(c, Grid_xBnum)
    var xPosGhost = is_xPosGhost(c, Grid_xBnum, Grid_xNum)
    var yNegGhost = is_yNegGhost(c, Grid_yBnum)
    var yPosGhost = is_yPosGhost(c, Grid_yBnum, Grid_yNum)
    var zNegGhost = is_zNegGhost(c, Grid_zBnum)
    var zPosGhost = is_zPosGhost(c, Grid_zBnum, Grid_zNum)
    var interior = in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if NSCBC_inflow_cell  then
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
    if interior or NSCBC_outflow_cell then
      var kineticEnergy = ((0.5*Fluid[c].rho)*dot(Fluid[c].velocity, Fluid[c].velocity))
      var pressure = ((Flow_gamma-1.0)*(Fluid[c].rhoEnergy-kineticEnergy))
      Fluid[c].pressure = pressure
      Fluid[c].temperature = (pressure/(Flow_gasConstant*Fluid[c].rho))
    end
  end
end

-- NOTE: It is safe to not pass the ghost regions to this task, because we
-- always group ghost cells with their neighboring interior cells.
__demand(__leaf, __cuda) -- MANUALLY PARALLELIZED
task Flow_UpdateGhostThermodynamics(Fluid : region(ispace(int3d), Fluid_columns),
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
  writes(Fluid.{pressure, temperature})
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

    var xNegGhost = is_xNegGhost(c, Grid_xBnum)
    var xPosGhost = is_xPosGhost(c, Grid_xBnum, Grid_xNum)
    var yNegGhost = is_yNegGhost(c, Grid_yBnum)
    var yPosGhost = is_yPosGhost(c, Grid_yBnum, Grid_yNum)
    var zNegGhost = is_zNegGhost(c, Grid_zBnum)
    var zPosGhost = is_zPosGhost(c, Grid_zBnum, Grid_zNum)

    if xNegGhost and BC_xBCLeft ~= SCHEMA.FlowBC_NSCBC_SubsonicInflow then
      var c_bnd = int3d(c)
      var c_int = ((c+{1, 0, 0})%Fluid.bounds)
      var bnd_temperature  = BC_xNegTemperature
      var temp_wall = 0.0
      var temperature = 0.0
      temp_wall = Fluid[c_int].temperature
      if (bnd_temperature>0.0) then
        temp_wall = bnd_temperature
      end
      temperature = (2.0*temp_wall)-Fluid[c_int].temperature
      Fluid[c_bnd].pressure = Fluid[c_int].pressure
      Fluid[c_bnd].temperature = temperature
    end

    if xPosGhost and BC_xBCRight ~= SCHEMA.FlowBC_NSCBC_SubsonicOutflow then
      var c_bnd = int3d(c)
      var c_int = ((c+{-1, 0, 0})%Fluid.bounds)
      var bnd_temperature = BC_xPosTemperature
      var temp_wall = 0.0
      var temperature = 0.0
      temp_wall = Fluid[c_int].temperature
      if (bnd_temperature>0.0) then
        temp_wall = bnd_temperature
      end
      temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
      Fluid[c_bnd].pressure = Fluid[c_int].pressure
      Fluid[c_bnd].temperature = temperature
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
        Fluid[c_bnd].pressure = Fluid[c_int].pressure
        Fluid[c_bnd].temperature = temperature
      else
        var bnd_temperature = BC_yNegTemperature
        var temp_wall = 0.0
        var temperature = 0.0
        temp_wall = Fluid[c_int].temperature
        if (bnd_temperature>0.0) then
          temp_wall = bnd_temperature
        end
        temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
        Fluid[c_bnd].pressure = Fluid[c_int].pressure
        Fluid[c_bnd].temperature = temperature
      end
    end

    if yPosGhost then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, -1, 0})%Fluid.bounds)
      if (BC_yBCRight == SCHEMA.FlowBC_NonUniformTemperatureWall) then
        var c_1 = 2.0/(Grid_xWidth*Grid_xWidth)*( (BC_yBCRightHeat_T_right - BC_yBCRightHeat_T_left) - 2.0*(BC_yBCRightHeat_T_mid - BC_yBCRightHeat_T_left))
        var c_2 = 4.0/(Grid_xWidth)*((BC_yBCRightHeat_T_mid - BC_yBCRightHeat_T_left) - 1.0/4.0*(BC_yBCRightHeat_T_right - BC_yBCRightHeat_T_left))
        var c_3 = BC_yBCRightHeat_T_left
        var wall_temperature = c_1*Fluid[c_bnd].centerCoordinates[0]*Fluid[c_bnd].centerCoordinates[0] + c_2*Fluid[c_bnd].centerCoordinates[0] + c_3
        if wall_temperature < 0.0 then --unphysical.... set wall themperature to zero
          wall_temperature = 0.0
        end
        var temperature = ((2.0*wall_temperature)-Fluid[c_int].temperature)
        Fluid[c_bnd].pressure = Fluid[c_int].pressure
        Fluid[c_bnd].temperature = temperature
      else
        var bnd_temperature = BC_yPosTemperature
        var temp_wall = 0.0
        var temperature = 0.0
        temp_wall = Fluid[c_int].temperature
        if (bnd_temperature>0.0) then
          temp_wall = bnd_temperature
        end
        temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
        Fluid[c_bnd].pressure = Fluid[c_int].pressure
        Fluid[c_bnd].temperature = temperature
      end
    end

    if zNegGhost then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, 0, 1})%Fluid.bounds)
      if (BC_zBCLeft == SCHEMA.FlowBC_NonUniformTemperatureWall) then
        var c_1 = 2.0/(Grid_xWidth*Grid_xWidth)*( (BC_zBCLeftHeat_T_right - BC_zBCLeftHeat_T_left) - 2.0*(BC_zBCLeftHeat_T_mid - BC_zBCLeftHeat_T_left))
        var c_2 = 4.0/(Grid_xWidth)*((BC_zBCLeftHeat_T_mid - BC_zBCLeftHeat_T_left) - 1.0/4.0*(BC_zBCLeftHeat_T_right - BC_zBCLeftHeat_T_left))
        var c_3 = BC_zBCLeftHeat_T_left
        var wall_temperature = c_1*Fluid[c_bnd].centerCoordinates[0]*Fluid[c_bnd].centerCoordinates[0] + c_2*Fluid[c_bnd].centerCoordinates[0] + c_3
        if wall_temperature < 0.0 then --unphysical.... set wall themperature to zero
          wall_temperature = 0.0
        end
        var temperature = ((2.0*wall_temperature)-Fluid[c_int].temperature)
        Fluid[c_bnd].pressure = Fluid[c_int].pressure
        Fluid[c_bnd].temperature = temperature
      else
        var bnd_temperature = BC_zNegTemperature
        var temp_wall = 0.0
        var temperature = 0.0
        temp_wall = Fluid[c_int].temperature
        if (bnd_temperature>0.0) then
          temp_wall = bnd_temperature
        end
        temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
        Fluid[c_bnd].pressure = Fluid[c_int].pressure
        Fluid[c_bnd].temperature = temperature
      end
    end

    if zPosGhost then
      var c_bnd = int3d(c)
      var c_int = ((c+{0, 0, -1})%Fluid.bounds)
      if (BC_zBCRight == SCHEMA.FlowBC_NonUniformTemperatureWall) then
        var c_1 = 2.0/(Grid_xWidth*Grid_xWidth)*( (BC_zBCRightHeat_T_right - BC_zBCRightHeat_T_left) - 2.0*(BC_zBCRightHeat_T_mid - BC_zBCRightHeat_T_left))
        var c_2 = 4.0/(Grid_xWidth)*((BC_zBCRightHeat_T_mid - BC_zBCRightHeat_T_left) - 1.0/4.0*(BC_zBCRightHeat_T_right - BC_zBCRightHeat_T_left))
        var c_3 = BC_zBCRightHeat_T_left
        var wall_temperature = c_1*Fluid[c_bnd].centerCoordinates[0]*Fluid[c_bnd].centerCoordinates[0] + c_2*Fluid[c_bnd].centerCoordinates[0] + c_3
        if wall_temperature < 0.0 then --unphysical.... set wall themperature to zero
          wall_temperature = 0.0
        end
        var temperature = ((2.0*wall_temperature)-Fluid[c_int].temperature)
        Fluid[c_bnd].pressure = Fluid[c_int].pressure
        Fluid[c_bnd].temperature = temperature
      else
        var bnd_temperature = BC_zPosTemperature
        var temp_wall = 0.0
        var temperature = 0.0
        temp_wall = Fluid[c_int].temperature
        if (bnd_temperature>0.0) then
          temp_wall = bnd_temperature
        end
        temperature = ((2.0*temp_wall)-Fluid[c_int].temperature)
        Fluid[c_bnd].pressure = Fluid[c_int].pressure
        Fluid[c_bnd].temperature = temperature
      end
    end

  end
end

__demand(__leaf, __parallel, __cuda)
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

__demand(__leaf, __parallel, __cuda)
task Flow_CalculateInteriorVolume(Fluid : region(ispace(int3d), Fluid_columns),
                                  Grid_xBnum : int32, Grid_xNum : int32,
                                  Grid_yBnum : int32, Grid_yNum : int32,
                                  Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.cellWidth)
do
  var acc = 0.0
  __demand(__openmp)
  for c in Fluid do
    if in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum) then
      acc += Fluid[c].cellWidth[0]*Fluid[c].cellWidth[1]*Fluid[c].cellWidth[2]
    end
  end
  return acc
end

__demand(__leaf, __parallel, __cuda)
task Flow_CalculateAveragePressure(Fluid : region(ispace(int3d), Fluid_columns),
                                   Grid_xBnum : int32, Grid_xNum : int32,
                                   Grid_yBnum : int32, Grid_yNum : int32,
                                   Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{cellWidth,pressure})
do
  var acc = 0.0
  __demand(__openmp)
  for c in Fluid do
    if in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum) then
      var cellVolume = Fluid[c].cellWidth[0]*Fluid[c].cellWidth[1]*Fluid[c].cellWidth[2]
      acc += Fluid[c].pressure*cellVolume
    end
  end
  return acc
end

__demand(__leaf, __parallel, __cuda)
task Flow_CalculateAverageTemperature(Fluid : region(ispace(int3d), Fluid_columns),
                                      Grid_xBnum : int32, Grid_xNum : int32,
                                      Grid_yBnum : int32, Grid_yNum : int32,
                                      Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{cellWidth,temperature})
do
  var acc = 0.0
  __demand(__openmp)
  for c in Fluid do
    if in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum) then
      var cellVolume = Fluid[c].cellWidth[0]*Fluid[c].cellWidth[1]*Fluid[c].cellWidth[2]
      acc += Fluid[c].temperature*cellVolume
    end
  end
  return acc
end

__demand(__leaf, __parallel, __cuda)
task Flow_CalculateAverageKineticEnergy(Fluid : region(ispace(int3d), Fluid_columns),
                                        Grid_xBnum : int32, Grid_xNum : int32,
                                        Grid_yBnum : int32, Grid_yNum : int32,
                                        Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{cellWidth, rho, velocity})
do
  var acc = 0.0
  __demand(__openmp)
  for c in Fluid do
    if in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum) then
      var cellVolume = Fluid[c].cellWidth[0]*Fluid[c].cellWidth[1]*Fluid[c].cellWidth[2]
      var kineticEnergy = 0.5*Fluid[c].rho*dot(Fluid[c].velocity, Fluid[c].velocity)
      acc += kineticEnergy*cellVolume
    end
  end
  return acc
end

__demand(__leaf, __parallel, __cuda)
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

__demand(__leaf, __parallel, __cuda)
task Flow_CalculateMaxMachNumber(Fluid : region(ispace(int3d), Fluid_columns),
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
    var xNegGhost = is_xNegGhost(c, Grid_xBnum)
    var xPosGhost = is_xPosGhost(c, Grid_xBnum, Grid_xNum)
    var yNegGhost = is_yNegGhost(c, Grid_yBnum)
    var yPosGhost = is_yPosGhost(c, Grid_yBnum, Grid_yNum)
    var zNegGhost = is_zNegGhost(c, Grid_zBnum)
    var zPosGhost = is_zPosGhost(c, Grid_zBnum, Grid_zNum)
    var ghost_cell = (xNegGhost or xPosGhost or
                      yNegGhost or yPosGhost or
                      zNegGhost or zPosGhost )
    var interior_cell = not (ghost_cell)
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    if interior_cell or NSCBC_inflow_cell or NSCBC_outflow_cell then
      var c_sound = GetSoundSpeed(Fluid[c].temperature, Flow_gamma, Flow_gasConstant)
      acc max= sqrt(dot(Fluid[c].velocity,Fluid[c].velocity))/c_sound
    end
  end
  return acc
end

__demand(__leaf, __parallel, __cuda)
task Flow_CalculateConvectiveSpectralRadius(Fluid : region(ispace(int3d), Fluid_columns),
                                            Flow_gamma : double,
                                            Flow_gasConstant : double)
where
  reads(Fluid.{cellWidth, velocity, temperature})
do
  var acc = -math.huge
  __demand(__openmp)
  for c in Fluid do
    var Grid_dXYZInverseSquare =
      1.0/Fluid[c].cellWidth[0]/Fluid[c].cellWidth[0] +
      1.0/Fluid[c].cellWidth[1]/Fluid[c].cellWidth[1] +
      1.0/Fluid[c].cellWidth[2]/Fluid[c].cellWidth[2]
    acc max= ((((fabs(Fluid[c].velocity[0])/Fluid[c].cellWidth[0])+(fabs(Fluid[c].velocity[1])/Fluid[c].cellWidth[1]))+(fabs(Fluid[c].velocity[2])/Fluid[c].cellWidth[2]))+(GetSoundSpeed(Fluid[c].temperature, Flow_gamma, Flow_gasConstant)*sqrt(Grid_dXYZInverseSquare)))
  end
  return acc
end

__demand(__leaf, __parallel, __cuda)
task Flow_CalculateViscousSpectralRadius(Fluid : region(ispace(int3d), Fluid_columns),
                                         Flow_constantVisc : double,
                                         Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                                         Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                                         Flow_viscosityModel : SCHEMA.ViscosityModel)
where
  reads(Fluid.{cellWidth, rho, temperature})
do
  var acc = -math.huge
  __demand(__openmp)
  for c in Fluid do
    var Grid_dXYZInverseSquare =
      1.0/Fluid[c].cellWidth[0]/Fluid[c].cellWidth[0] +
      1.0/Fluid[c].cellWidth[1]/Fluid[c].cellWidth[1] +
      1.0/Fluid[c].cellWidth[2]/Fluid[c].cellWidth[2]
    var dynamicViscosity = GetDynamicViscosity(Fluid[c].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)
    acc max= ((((2.0*dynamicViscosity)/Fluid[c].rho)*Grid_dXYZInverseSquare)*4.0)
  end
  return acc
end

__demand(__leaf, __parallel, __cuda)
task Flow_CalculateHeatConductionSpectralRadius(Fluid : region(ispace(int3d), Fluid_columns),
                                                Flow_constantVisc : double,
                                                Flow_gamma : double,
                                                Flow_gasConstant : double,
                                                Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                                                Flow_prandtl : double,
                                                Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                                                Flow_viscosityModel : SCHEMA.ViscosityModel)
where
  reads(Fluid.{cellWidth, rho, temperature})
do
  var acc = -math.huge
  __demand(__openmp)
  for c in Fluid do
    var Grid_dXYZInverseSquare =
      1.0/Fluid[c].cellWidth[0]/Fluid[c].cellWidth[0] +
      1.0/Fluid[c].cellWidth[1]/Fluid[c].cellWidth[1] +
      1.0/Fluid[c].cellWidth[2]/Fluid[c].cellWidth[2]
    var dynamicViscosity = GetDynamicViscosity(Fluid[c].temperature, Flow_constantVisc, Flow_powerlawTempRef, Flow_powerlawViscRef, Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef, Flow_viscosityModel)
    var cv = (Flow_gasConstant/(Flow_gamma-1.0))
    var cp = (Flow_gamma*cv)
    var kappa = ((cp/Flow_prandtl)*dynamicViscosity)
    acc max= (((kappa/(cv*Fluid[c].rho))*Grid_dXYZInverseSquare)*4.0)
  end
  return acc
end

__demand(__leaf, __parallel, __cuda)
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

__demand(__leaf, __parallel, __cuda)
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

__demand(__leaf, __parallel, __cuda)
task Flow_InitializeTimeDerivatives(Fluid : region(ispace(int3d), Fluid_columns))
where
  writes(Fluid.{rho_t, rhoVelocity_t, rhoEnergy_t})
do
  __demand(__openmp)
  for c in Fluid do
    Fluid[c].rho_t = 0.0
    Fluid[c].rhoVelocity_t = array(0.0, 0.0, 0.0)
    Fluid[c].rhoEnergy_t = 0.0
  end
end

-- NOTE: It is safe to not pass the ghost regions to this task, because we
-- always group ghost cells with their neighboring interior cells.
__demand(__leaf, __cuda) -- MANUALLY PARALLELIZED
task Flow_UpdateGhostGradients(Fluid : region(ispace(int3d), Fluid_columns),
                               config : Config,
                               BC_xNegSign : double[3], BC_yNegSign : double[3], BC_zNegSign : double[3],
                               BC_xPosSign : double[3], BC_yPosSign : double[3], BC_zPosSign : double[3],
                               Grid_xBnum : int32, Grid_xNum : int32,
                               Grid_yBnum : int32, Grid_yNum : int32,
                               Grid_zBnum : int32, Grid_zNum : int32)
where
  reads writes(Fluid.{velocityGradientX, velocityGradientY, velocityGradientZ, temperatureGradient})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = is_xNegGhost(c, Grid_xBnum)
    var xPosGhost = is_xPosGhost(c, Grid_xBnum, Grid_xNum)
    var yNegGhost = is_yNegGhost(c, Grid_yBnum)
    var yPosGhost = is_yPosGhost(c, Grid_yBnum, Grid_yNum)
    var zNegGhost = is_zNegGhost(c, Grid_zBnum)
    var zPosGhost = is_zPosGhost(c, Grid_zBnum, Grid_zNum)
    if xNegGhost and BC_xBCLeft ~= SCHEMA.FlowBC_NSCBC_SubsonicInflow then
      Fluid[c].velocityGradientX   = vv_mul(BC_xNegSign, Fluid[(c+{1, 0, 0})%Fluid.bounds].velocityGradientX)
      Fluid[c].velocityGradientY   = vv_mul(BC_xNegSign, Fluid[(c+{1, 0, 0})%Fluid.bounds].velocityGradientY)
      Fluid[c].velocityGradientZ   = vv_mul(BC_xNegSign, Fluid[(c+{1, 0, 0})%Fluid.bounds].velocityGradientZ)
      Fluid[c].temperatureGradient = vv_mul(BC_xNegSign, Fluid[(c+{1, 0, 0})%Fluid.bounds].temperatureGradient)
    end
    if xPosGhost and BC_xBCRight ~= SCHEMA.FlowBC_NSCBC_SubsonicOutflow then
      Fluid[c].velocityGradientX   = vv_mul(BC_xPosSign, Fluid[(c+{-1, 0, 0})%Fluid.bounds].velocityGradientX)
      Fluid[c].velocityGradientY   = vv_mul(BC_xPosSign, Fluid[(c+{-1, 0, 0})%Fluid.bounds].velocityGradientY)
      Fluid[c].velocityGradientZ   = vv_mul(BC_xPosSign, Fluid[(c+{-1, 0, 0})%Fluid.bounds].velocityGradientZ)
      Fluid[c].temperatureGradient = vv_mul(BC_xPosSign, Fluid[(c+{-1, 0, 0})%Fluid.bounds].temperatureGradient)
    end
    if yNegGhost then
      Fluid[c].velocityGradientX   = vv_mul(BC_yNegSign, Fluid[(c+{0, 1, 0})%Fluid.bounds].velocityGradientX)
      Fluid[c].velocityGradientY   = vv_mul(BC_yNegSign, Fluid[(c+{0, 1, 0})%Fluid.bounds].velocityGradientY)
      Fluid[c].velocityGradientZ   = vv_mul(BC_yNegSign, Fluid[(c+{0, 1, 0})%Fluid.bounds].velocityGradientZ)
      Fluid[c].temperatureGradient = vv_mul(BC_yNegSign, Fluid[(c+{0, 1, 0})%Fluid.bounds].temperatureGradient)
    end
    if yPosGhost then
      Fluid[c].velocityGradientX   = vv_mul(BC_yPosSign, Fluid[(c+{0, -1, 0})%Fluid.bounds].velocityGradientX)
      Fluid[c].velocityGradientY   = vv_mul(BC_yPosSign, Fluid[(c+{0, -1, 0})%Fluid.bounds].velocityGradientY)
      Fluid[c].velocityGradientZ   = vv_mul(BC_yPosSign, Fluid[(c+{0, -1, 0})%Fluid.bounds].velocityGradientZ)
      Fluid[c].temperatureGradient = vv_mul(BC_yPosSign, Fluid[(c+{0, -1, 0})%Fluid.bounds].temperatureGradient)
    end
    if zNegGhost then
      Fluid[c].velocityGradientX   = vv_mul(BC_zNegSign, Fluid[(c+{0, 0, 1})%Fluid.bounds].velocityGradientX)
      Fluid[c].velocityGradientY   = vv_mul(BC_zNegSign, Fluid[(c+{0, 0, 1})%Fluid.bounds].velocityGradientY)
      Fluid[c].velocityGradientZ   = vv_mul(BC_zNegSign, Fluid[(c+{0, 0, 1})%Fluid.bounds].velocityGradientZ)
      Fluid[c].temperatureGradient = vv_mul(BC_zNegSign, Fluid[(c+{0, 0, 1})%Fluid.bounds].temperatureGradient)
    end
    if zPosGhost then
      Fluid[c].velocityGradientX   = vv_mul(BC_zPosSign, Fluid[(c+{0, 0, -1})%Fluid.bounds].velocityGradientX)
      Fluid[c].velocityGradientY   = vv_mul(BC_zPosSign, Fluid[(c+{0, 0, -1})%Fluid.bounds].velocityGradientY)
      Fluid[c].velocityGradientZ   = vv_mul(BC_zPosSign, Fluid[(c+{0, 0, -1})%Fluid.bounds].velocityGradientZ)
      Fluid[c].temperatureGradient = vv_mul(BC_zPosSign, Fluid[(c+{0, 0, -1})%Fluid.bounds].temperatureGradient)
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_GetFluxX(Fluid : region(ispace(int3d), Fluid_columns),
                   config : Config,
                   Flow_constantVisc : double,
                   Flow_gamma : double,
                   Flow_gasConstant : double,
                   Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                   Flow_prandtl : double,
                   Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                   Flow_viscosityModel : SCHEMA.ViscosityModel,
                   Grid_xBnum : int32, Grid_xNum : int32,
                   Grid_yBnum : int32, Grid_yNum : int32,
                   Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{centerCoordinates, cellWidth}),
  reads(Fluid.{rho, pressure, velocity, rhoVelocity, rhoEnergy, temperature}),
  reads(Fluid.{velocityGradientX, velocityGradientY, velocityGradientZ, temperatureGradient}),
  writes(Fluid.{rhoEnergyFluxX, rhoFluxX, rhoVelocityFluxX})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight

  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = is_xNegGhost(c, Grid_xBnum)
    var xPosGhost = is_xPosGhost(c, Grid_xBnum, Grid_xNum)
    var yNegGhost = is_yNegGhost(c, Grid_yBnum)
    var yPosGhost = is_yPosGhost(c, Grid_yBnum, Grid_yNum)
    var zNegGhost = is_zNegGhost(c, Grid_zBnum)
    var zPosGhost = is_zPosGhost(c, Grid_zBnum, Grid_zNum)
    var interior = in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if interior or xNegGhost  then
      var rho = Fluid[c].rho
      var pressure = Fluid[c].pressure
      var rhoVelocity = Fluid[c].rhoVelocity
      var rhoEnergy = Fluid[c].rhoEnergy
      var temperature = Fluid[c].temperature
      var temperatureGradient = Fluid[c].temperatureGradient
      var velocity = Fluid[c].velocity
      var velocityGradientX = Fluid[c].velocityGradientX
      var velocityGradientY = Fluid[c].velocityGradientY
      var velocityGradientZ = Fluid[c].velocityGradientZ
      var mu = GetDynamicViscosity(temperature,
                                   Flow_constantVisc,
                                   Flow_powerlawTempRef, Flow_powerlawViscRef,
                                   Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef,
                                   Flow_viscosityModel)
      var xCellWidth = Fluid[c].cellWidth[0]

      var stencil = (c+{1, 0, 0}) % Fluid.bounds
      var rho_stencil = Fluid[stencil].rho
      var pressure_stencil = Fluid[stencil].pressure
      var rhoVelocity_stencil = Fluid[stencil].rhoVelocity
      var rhoEnergy_stencil = Fluid[stencil].rhoEnergy
      var temperature_stencil = Fluid[stencil].temperature
      var temperatureGradient_stencil = Fluid[stencil].temperatureGradient
      var velocity_stencil = Fluid[stencil].velocity
      var velocityGradientX_stencil = Fluid[stencil].velocityGradientX
      var velocityGradientY_stencil = Fluid[stencil].velocityGradientY
      var velocityGradientZ_stencil = Fluid[stencil].velocityGradientZ
      var mu_stencil = GetDynamicViscosity(temperature_stencil,
                                           Flow_constantVisc,
                                           Flow_powerlawTempRef, Flow_powerlawViscRef,
                                           Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef,
                                           Flow_viscosityModel)
      var xCellWidth_stencil = Fluid[stencil].cellWidth[0]


      var muFace = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                        Fluid[c].centerCoordinates[0],
                                        Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                        mu,
                                        mu_stencil)

      var rhoFace = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                         Fluid[c].centerCoordinates[0],
                                         Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                         rho,
                                         rho_stencil)

      var pressureFace = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                              Fluid[c].centerCoordinates[0],
                                              Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                              pressure,
                                              pressure_stencil)

      var rhoEnergyFace = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                               Fluid[c].centerCoordinates[0],
                                               Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                               rhoEnergy,
                                               rhoEnergy_stencil)

      var velocityFace = array(0.0, 0.0, 0.0)
      velocityFace[0] = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                             Fluid[c].centerCoordinates[0],
                                             Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                             velocity[0],
                                             velocity_stencil[0])
      velocityFace[1] = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                             Fluid[c].centerCoordinates[0],
                                             Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                             velocity[1],
                                             velocity_stencil[1])
      velocityFace[2] = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                             Fluid[c].centerCoordinates[0],
                                             Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                             velocity[2],
                                             velocity_stencil[2])


      var velocityX_YFace = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[0],
                                                 Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                                 velocityGradientY[0],
                                                 velocityGradientY_stencil[0])
      var velocityX_ZFace = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[0],
                                                 Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                                 velocityGradientZ[0],
                                                 velocityGradientZ_stencil[0])
      var velocityY_YFace = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[0],
                                                 Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                                 velocityGradientY[1],
                                                 velocityGradientY_stencil[1])
      var velocityZ_ZFace = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[0],
                                                 Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                                 velocityGradientZ[2],
                                                 velocityGradientZ_stencil[2])

      var velocityX_XFace   = 0.0
      var velocityY_XFace   = 0.0
      var velocityZ_XFace   = 0.0
      var temperature_XFace = 0.0
      if xNegGhost or is_xPosGhost(stencil, Grid_xBnum, Grid_xNum) then
        -- do not use velocity gradients in ghost cells for computing the fluxes
        velocityX_XFace   = (velocity_stencil[0] - velocity[0]) / xCellWidth
        velocityY_XFace   = (velocity_stencil[1] - velocity[1]) / xCellWidth
        velocityZ_XFace   = (velocity_stencil[2] - velocity[2]) / xCellWidth
        temperature_XFace = (temperature_stencil - temperature) / xCellWidth
      else
        velocityX_XFace = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                               Fluid[c].centerCoordinates[0],
                                               Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                               velocityGradientX[0],
                                               velocityGradientX_stencil[0])
        velocityY_XFace = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                               Fluid[c].centerCoordinates[0],
                                               Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                               velocityGradientX[1],
                                               velocityGradientX_stencil[1])
        velocityZ_XFace = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                               Fluid[c].centerCoordinates[0],
                                               Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                               velocityGradientX[2],
                                               velocityGradientX_stencil[2])
        temperature_XFace = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[0],
                                                 Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                                 temperatureGradient[0],
                                                 temperatureGradient_stencil[0])
      end

      var sigmaXX = muFace*(4.0*velocityX_XFace-2.0*velocityY_YFace-2.0*velocityZ_ZFace)/3.0
      var sigmaYX = muFace*(velocityY_XFace+velocityX_YFace)
      var sigmaZX = muFace*(velocityZ_XFace+velocityX_ZFace)

      var usigma = velocityFace[0]*sigmaXX + velocityFace[1]*sigmaYX + velocityFace[2]*sigmaZX

      var cp = Flow_gamma * Flow_gasConstant / (Flow_gamma-1.0)
      var heatFlux = (-(cp*muFace/Flow_prandtl))*temperature_XFace

      -- Density Flux
      Fluid[c].rhoFluxX = rhoFace*velocityFace[0]

      -- Momentum Flux
      var rhoVelocityFluxX = vs_mul(vs_mul(velocityFace,rhoFace), velocityFace[0])
      rhoVelocityFluxX[0] += pressureFace
      Fluid[c].rhoVelocityFluxX = vv_sub(rhoVelocityFluxX, array(sigmaXX,sigmaYX,sigmaZX))

      -- Energy Flux Flux
      Fluid[c].rhoEnergyFluxX = (rhoEnergyFace + pressureFace)*velocityFace[0] - (usigma-heatFlux)
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_GetFluxY(Fluid : region(ispace(int3d), Fluid_columns),
                   config : Config,
                   Flow_constantVisc : double,
                   Flow_gamma : double,
                   Flow_gasConstant : double,
                   Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                   Flow_prandtl : double,
                   Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                   Flow_viscosityModel : SCHEMA.ViscosityModel,
                   Grid_xBnum : int32, Grid_xNum : int32,
                   Grid_yBnum : int32, Grid_yNum : int32,
                   Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{centerCoordinates, cellWidth}),
  reads(Fluid.{rho, pressure, velocity, rhoVelocity, rhoEnergy, temperature}),
  reads(Fluid.{velocityGradientX, velocityGradientY, velocityGradientZ, temperatureGradient}),
  writes(Fluid.{rhoEnergyFluxY, rhoFluxY, rhoVelocityFluxY})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight

  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = is_xNegGhost(c, Grid_xBnum)
    var xPosGhost = is_xPosGhost(c, Grid_xBnum, Grid_xNum)
    var yNegGhost = is_yNegGhost(c, Grid_yBnum)
    var yPosGhost = is_yPosGhost(c, Grid_yBnum, Grid_yNum)
    var zNegGhost = is_zNegGhost(c, Grid_zBnum)
    var zPosGhost = is_zPosGhost(c, Grid_zBnum, Grid_zNum)
    var interior = in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if interior or yNegGhost or NSCBC_inflow_cell or NSCBC_outflow_cell then
      var rho = Fluid[c].rho
      var pressure = Fluid[c].pressure
      var rhoVelocity = Fluid[c].rhoVelocity
      var rhoEnergy = Fluid[c].rhoEnergy
      var temperature = Fluid[c].temperature
      var temperatureGradient = Fluid[c].temperatureGradient
      var velocity = Fluid[c].velocity
      var velocityGradientX = Fluid[c].velocityGradientX
      var velocityGradientY = Fluid[c].velocityGradientY
      var velocityGradientZ = Fluid[c].velocityGradientZ
      var mu = GetDynamicViscosity(temperature,
                                   Flow_constantVisc,
                                   Flow_powerlawTempRef, Flow_powerlawViscRef,
                                   Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef,
                                   Flow_viscosityModel)
      var yCellWidth = Fluid[c].cellWidth[1]

      var stencil = (c+{0, 1, 0}) % Fluid.bounds
      var rho_stencil = Fluid[stencil].rho
      var pressure_stencil = Fluid[stencil].pressure
      var rhoVelocity_stencil = Fluid[stencil].rhoVelocity
      var rhoEnergy_stencil = Fluid[stencil].rhoEnergy
      var temperature_stencil = Fluid[stencil].temperature
      var temperatureGradient_stencil = Fluid[stencil].temperatureGradient
      var velocity_stencil = Fluid[stencil].velocity
      var velocityGradientX_stencil = Fluid[stencil].velocityGradientX
      var velocityGradientY_stencil = Fluid[stencil].velocityGradientY
      var velocityGradientZ_stencil = Fluid[stencil].velocityGradientZ
      var mu_stencil = GetDynamicViscosity(temperature_stencil,
                                           Flow_constantVisc,
                                           Flow_powerlawTempRef, Flow_powerlawViscRef,
                                           Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef,
                                           Flow_viscosityModel)
      var yCellWidth_stencil = Fluid[stencil].cellWidth[1]

      var muFace = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                        Fluid[c].centerCoordinates[1],
                                        Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                        mu,
                                        mu_stencil)

      var rhoFace = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                         Fluid[c].centerCoordinates[1],
                                         Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                         rho,
                                         rho_stencil)

      var pressureFace = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                              Fluid[c].centerCoordinates[1],
                                              Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                              pressure,
                                              pressure_stencil)

     var rhoEnergyFace = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                               Fluid[c].centerCoordinates[1],
                                               Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                               rhoEnergy,
                                               rhoEnergy_stencil)

      var velocityFace = array(0.0, 0.0, 0.0)
      velocityFace[0] = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                             Fluid[c].centerCoordinates[1],
                                             Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                             velocity[0],
                                             velocity_stencil[0])
      velocityFace[1] = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                             Fluid[c].centerCoordinates[1],
                                             Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                             velocity[1],
                                             velocity_stencil[1])
      velocityFace[2] = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                             Fluid[c].centerCoordinates[1],
                                             Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                             velocity[2],
                                             velocity_stencil[2])

      var velocityX_XFace = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[1],
                                                 Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                                 velocityGradientX[0],
                                                 velocityGradientX_stencil[0])
      var velocityY_XFace = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[1],
                                                 Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                                 velocityGradientX[1],
                                                 velocityGradientX_stencil[1])
      var velocityY_ZFace = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[1],
                                                 Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                                 velocityGradientZ[1],
                                                 velocityGradientZ_stencil[1])
      var velocityZ_ZFace = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[1],
                                                 Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                                 velocityGradientZ[2],
                                                 velocityGradientZ_stencil[2])

      var velocityX_YFace   = 0.0
      var velocityY_YFace   = 0.0
      var velocityZ_YFace   = 0.0
      var temperature_YFace = 0.0
      if yNegGhost or is_yPosGhost(stencil, Grid_yBnum, Grid_yNum) then
        -- do not use velocity gradients in ghost cells for computing the fluxes
        velocityX_YFace   = (velocity_stencil[0] - velocity[0]) / yCellWidth
        velocityY_YFace   = (velocity_stencil[1] - velocity[1]) / yCellWidth
        velocityZ_YFace   = (velocity_stencil[2] - velocity[2]) / yCellWidth
        temperature_YFace = (temperature_stencil - temperature) / yCellWidth
      else
        velocityX_YFace = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                               Fluid[c].centerCoordinates[1],
                                               Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                               velocityGradientY[0],
                                               velocityGradientY_stencil[0])
        velocityY_YFace = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                               Fluid[c].centerCoordinates[1],
                                               Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                               velocityGradientY[1],
                                               velocityGradientY_stencil[1])
        velocityZ_YFace = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                               Fluid[c].centerCoordinates[1],
                                               Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                               velocityGradientY[2],
                                               velocityGradientY_stencil[2])
        temperature_YFace = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[1],
                                                 Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                                 temperatureGradient[1],
                                                 temperatureGradient_stencil[1])
      end

      var sigmaXY = muFace*(velocityX_YFace+velocityY_XFace)
      var sigmaYY = muFace*(4.0*velocityY_YFace-2.0*velocityX_XFace-2.0*velocityZ_ZFace)/3.0
      var sigmaZY = muFace*(velocityZ_YFace+velocityY_ZFace)

      var usigma = velocityFace[0]*sigmaXY + velocityFace[1]*sigmaYY + velocityFace[2]*sigmaZY

      var cp = Flow_gamma * Flow_gasConstant / (Flow_gamma-1.0)
      var heatFlux = (-(cp*muFace/Flow_prandtl))*temperature_YFace

      -- Density Flux
      Fluid[c].rhoFluxY = rhoFace*velocityFace[1]

      -- Momentum Flux
      var rhoVelocityFluxY = vs_mul(vs_mul(velocityFace, rhoFace), velocityFace[1])
      rhoVelocityFluxY[1] += pressureFace
      Fluid[c].rhoVelocityFluxY = vv_sub(rhoVelocityFluxY, array(sigmaXY,sigmaYY,sigmaZY))

      -- Energy Flux Flux
      Fluid[c].rhoEnergyFluxY = (rhoEnergyFace + pressureFace)*velocityFace[1] - (usigma-heatFlux)
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_GetFluxZ(Fluid : region(ispace(int3d), Fluid_columns),
                   config : Config,
                   Flow_constantVisc : double,
                   Flow_gamma : double,
                   Flow_gasConstant : double,
                   Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                   Flow_prandtl : double,
                   Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                   Flow_viscosityModel : SCHEMA.ViscosityModel,
                   Grid_xBnum : int32, Grid_xNum : int32,
                   Grid_yBnum : int32, Grid_yNum : int32,
                   Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{centerCoordinates, cellWidth}),
  reads(Fluid.{rho, pressure, velocity, rhoVelocity, rhoEnergy, temperature}),
  reads(Fluid.{velocityGradientX, velocityGradientY, velocityGradientZ, temperatureGradient}),
  writes(Fluid.{rhoEnergyFluxZ, rhoFluxZ, rhoVelocityFluxZ})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight

  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = is_xNegGhost(c, Grid_xBnum)
    var xPosGhost = is_xPosGhost(c, Grid_xBnum, Grid_xNum)
    var yNegGhost = is_yNegGhost(c, Grid_yBnum)
    var yPosGhost = is_yPosGhost(c, Grid_yBnum, Grid_yNum)
    var zNegGhost = is_zNegGhost(c, Grid_zBnum)
    var zPosGhost = is_zPosGhost(c, Grid_zBnum, Grid_zNum)
    var interior = in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if interior or zNegGhost or NSCBC_inflow_cell or NSCBC_outflow_cell then
      var rho = Fluid[c].rho
      var pressure = Fluid[c].pressure
      var rhoVelocity = Fluid[c].rhoVelocity
      var rhoEnergy = Fluid[c].rhoEnergy
      var temperature = Fluid[c].temperature
      var temperatureGradient = Fluid[c].temperatureGradient
      var velocity = Fluid[c].velocity
      var velocityGradientX = Fluid[c].velocityGradientX
      var velocityGradientY = Fluid[c].velocityGradientY
      var velocityGradientZ = Fluid[c].velocityGradientZ
      var mu = GetDynamicViscosity(temperature,
                                   Flow_constantVisc,
                                   Flow_powerlawTempRef, Flow_powerlawViscRef,
                                   Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef,
                                   Flow_viscosityModel)
      var zCellWidth = Fluid[c].cellWidth[2]



      var stencil = (c+{0, 0, 1}) % Fluid.bounds
      var rho_stencil = Fluid[stencil].rho
      var pressure_stencil = Fluid[stencil].pressure
      var rhoVelocity_stencil = Fluid[stencil].rhoVelocity
      var rhoEnergy_stencil = Fluid[stencil].rhoEnergy
      var temperature_stencil = Fluid[stencil].temperature
      var temperatureGradient_stencil = Fluid[stencil].temperatureGradient
      var velocity_stencil = Fluid[stencil].velocity
      var velocityGradientX_stencil = Fluid[stencil].velocityGradientX
      var velocityGradientY_stencil = Fluid[stencil].velocityGradientY
      var velocityGradientZ_stencil = Fluid[stencil].velocityGradientZ
      var mu_stencil = GetDynamicViscosity(temperature_stencil,
                                           Flow_constantVisc,
                                           Flow_powerlawTempRef, Flow_powerlawViscRef,
                                           Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef,
                                           Flow_viscosityModel)
      var zCellWidth_stencil = Fluid[stencil].cellWidth[2]

      var muFace = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                        Fluid[c].centerCoordinates[2],
                                        Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                        mu,
                                        mu_stencil)

      var rhoFace = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                         Fluid[c].centerCoordinates[2],
                                         Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                         rho,
                                         rho_stencil)

      var pressureFace = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                              Fluid[c].centerCoordinates[2],
                                              Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                              pressure,
                                              pressure_stencil)

      var rhoEnergyFace = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                               Fluid[c].centerCoordinates[2],
                                               Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                               rhoEnergy,
                                               rhoEnergy_stencil)

      var velocityFace = array(0.0, 0.0, 0.0)
      velocityFace[0] = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                             Fluid[c].centerCoordinates[2],
                                             Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                             velocity[0],
                                             velocity_stencil[0])
      velocityFace[1] = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                             Fluid[c].centerCoordinates[2],
                                             Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                             velocity[1],
                                             velocity_stencil[1])
      velocityFace[2] = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                             Fluid[c].centerCoordinates[2],
                                             Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                             velocity[2],
                                             velocity_stencil[2])

      var velocityX_XFace = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[2],
                                                 Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                                 velocityGradientX[0],
                                                 velocityGradientX_stencil[0])
      var velocityY_YFace = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[2],
                                                 Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                                 velocityGradientY[1],
                                                 velocityGradientY_stencil[1])
      var velocityZ_XFace = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[2],
                                                 Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                                 velocityGradientX[2],
                                                 velocityGradientX_stencil[2])
      var velocityZ_YFace = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[2],
                                                 Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                                 velocityGradientY[2],
                                                 velocityGradientY_stencil[2])

      var velocityX_ZFace   = 0.0
      var velocityY_ZFace   = 0.0
      var velocityZ_ZFace   = 0.0
      var temperature_ZFace = 0.0
      if zNegGhost or is_zPosGhost(stencil, Grid_zBnum, Grid_zNum) then
        -- do not use velocity gradients in ghost cells for computing the fluxes
        velocityX_ZFace   = (velocity_stencil[0] - velocity[0]) / zCellWidth
        velocityY_ZFace   = (velocity_stencil[1] - velocity[1]) / zCellWidth
        velocityZ_ZFace   = (velocity_stencil[2] - velocity[2]) / zCellWidth
        temperature_ZFace = (temperature_stencil - temperature) / zCellWidth
      else
        velocityX_ZFace = linear_interpolation(Fluid[c].centerCoordinates[1] + zCellWidth/2.0,
                                               Fluid[c].centerCoordinates[1],
                                               Fluid[c].centerCoordinates[1] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                               velocityGradientZ[0],
                                               velocityGradientZ_stencil[0])
        velocityY_ZFace = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                               Fluid[c].centerCoordinates[2],
                                               Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                               velocityGradientZ[1],
                                               velocityGradientZ_stencil[1])
        velocityZ_ZFace = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                               Fluid[c].centerCoordinates[2],
                                               Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                               velocityGradientZ[2],
                                               velocityGradientZ_stencil[2])
        temperature_ZFace = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[2],
                                                 Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                                 temperatureGradient[2],
                                                 temperatureGradient_stencil[2])
      end

      var sigmaXZ = muFace*(velocityX_ZFace+velocityZ_XFace)
      var sigmaYZ = muFace*(velocityY_ZFace+velocityZ_YFace)
      var sigmaZZ = muFace*(4.0*velocityZ_ZFace-2.0*velocityX_XFace-2.0*velocityY_YFace)/3.0

      var usigma = velocityFace[0]*sigmaXZ + velocityFace[1]*sigmaYZ + velocityFace[2]*sigmaZZ

      var cp = Flow_gamma * Flow_gasConstant / (Flow_gamma-1.0)
      var heatFlux = (-(cp*muFace/Flow_prandtl))*temperature_ZFace

      -- Density Flux
      Fluid[c].rhoFluxZ = rhoFace*velocityFace[2]

      -- Momentum Flux
      var rhoVelocityFluxZ = vs_mul(vs_mul(velocityFace, rhoFace), velocityFace[2])
      rhoVelocityFluxZ[2] += pressureFace
      Fluid[c].rhoVelocityFluxZ = vv_sub(rhoVelocityFluxZ, array(sigmaXZ,sigmaYZ,sigmaZZ))

      -- Energy Flux Flux
      Fluid[c].rhoEnergyFluxZ = (rhoEnergyFace + pressureFace)*velocityFace[2] - (usigma-heatFlux)
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_UpdateUsingFluxX(Fluid : region(ispace(int3d), Fluid_columns),
                           config : Config,
                           Grid_xBnum : int32, Grid_xNum : int32,
                           Grid_yBnum : int32, Grid_yNum : int32,
                           Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{cellWidth, rhoFluxX, rhoVelocityFluxX, rhoEnergyFluxX}),
  reads writes(Fluid.{rho_t, rhoVelocity_t, rhoEnergy_t})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = is_xNegGhost(c, Grid_xBnum)
    var xPosGhost = is_xPosGhost(c, Grid_xBnum, Grid_xNum)
    var yNegGhost = is_yNegGhost(c, Grid_yBnum)
    var yPosGhost = is_yPosGhost(c, Grid_yBnum, Grid_yNum)
    var zNegGhost = is_zNegGhost(c, Grid_zBnum)
    var zPosGhost = is_zPosGhost(c, Grid_zBnum, Grid_zNum)
    var interior = in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if interior then
      var xCellWidth = Fluid[c].cellWidth[0]
      Fluid[c].rho_t += ((-(Fluid[c].rhoFluxX-Fluid[(c+{-1, 0, 0})%Fluid.bounds].rhoFluxX))/xCellWidth);
      [UTIL.emitArrayReduce(3, '+',
         rexpr Fluid[c].rhoVelocity_t end,
         rexpr vs_div(vs_mul(vv_sub(Fluid[c].rhoVelocityFluxX, Fluid[(c+{-1, 0, 0})%Fluid.bounds].rhoVelocityFluxX), double((-1))), xCellWidth) end)];
      Fluid[c].rhoEnergy_t += ((-(Fluid[c].rhoEnergyFluxX-Fluid[(c+{-1, 0, 0})%Fluid.bounds].rhoEnergyFluxX))/xCellWidth)
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_UpdateUsingFluxY(Fluid : region(ispace(int3d), Fluid_columns),
                           config : Config,
                           Grid_xBnum : int32, Grid_xNum : int32,
                           Grid_yBnum : int32, Grid_yNum : int32,
                           Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{cellWidth, rhoFluxY, rhoVelocityFluxY, rhoEnergyFluxY}),
  reads writes(Fluid.{rho_t, rhoVelocity_t, rhoEnergy_t})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = is_xNegGhost(c, Grid_xBnum)
    var xPosGhost = is_xPosGhost(c, Grid_xBnum, Grid_xNum)
    var yNegGhost = is_yNegGhost(c, Grid_yBnum)
    var yPosGhost = is_yPosGhost(c, Grid_yBnum, Grid_yNum)
    var zNegGhost = is_zNegGhost(c, Grid_zBnum)
    var zPosGhost = is_zPosGhost(c, Grid_zBnum, Grid_zNum)
    var interior = in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if interior or NSCBC_inflow_cell or NSCBC_outflow_cell then
      var yCellWidth = Fluid[c].cellWidth[1]
      Fluid[c].rho_t += (-(Fluid[c].rhoFluxY-Fluid[(c+{0, -1, 0})%Fluid.bounds].rhoFluxY))/yCellWidth;
      [UTIL.emitArrayReduce(3, '+',
         rexpr Fluid[c].rhoVelocity_t end,
         rexpr vs_div(vs_mul(vv_sub(Fluid[c].rhoVelocityFluxY, Fluid[(c+{0, -1, 0})%Fluid.bounds].rhoVelocityFluxY), double((-1))), yCellWidth) end)];
      Fluid[c].rhoEnergy_t += (-(Fluid[c].rhoEnergyFluxY-Fluid[(c+{0, -1, 0})%Fluid.bounds].rhoEnergyFluxY))/yCellWidth
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_UpdateUsingFluxZ(Fluid : region(ispace(int3d), Fluid_columns),
                           config : Config,
                           Grid_xBnum : int32, Grid_xNum : int32,
                           Grid_yBnum : int32, Grid_yNum : int32,
                           Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{cellWidth, rhoFluxZ, rhoVelocityFluxZ, rhoEnergyFluxZ}),
  reads writes(Fluid.{rho_t, rhoVelocity_t, rhoEnergy_t})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = is_xNegGhost(c, Grid_xBnum)
    var xPosGhost = is_xPosGhost(c, Grid_xBnum, Grid_xNum)
    var yNegGhost = is_yNegGhost(c, Grid_yBnum)
    var yPosGhost = is_yPosGhost(c, Grid_yBnum, Grid_yNum)
    var zNegGhost = is_zNegGhost(c, Grid_zBnum)
    var zPosGhost = is_zPosGhost(c, Grid_zBnum, Grid_zNum)
    var interior = in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if interior or NSCBC_inflow_cell or NSCBC_outflow_cell then
      var zCellWidth = Fluid[c].cellWidth[2]
      Fluid[c].rho_t += (-(Fluid[c].rhoFluxZ-Fluid[(c+{0, 0, -1})%Fluid.bounds].rhoFluxZ))/zCellWidth;
      [UTIL.emitArrayReduce(3, '+',
         rexpr Fluid[c].rhoVelocity_t end,
         rexpr vs_div(vs_mul(vv_sub(Fluid[c].rhoVelocityFluxZ, Fluid[(c+{0, 0, -1})%Fluid.bounds].rhoVelocityFluxZ), double((-1))), zCellWidth) end)];
      Fluid[c].rhoEnergy_t += (-(Fluid[c].rhoEnergyFluxZ-Fluid[(c+{0, 0, -1})%Fluid.bounds].rhoEnergyFluxZ))/zCellWidth
    end
  end
end

-- NOTE: It is safe to not pass the ghost regions to this task, because we
-- always group ghost cells with their neighboring interior cells.
__demand(__leaf, __cuda) -- MANUALLY PARALLELIZED
task Flow_UpdateUsingFluxGhostNSCBC(Fluid : region(ispace(int3d), Fluid_columns),
                                    config : Config,
                                    Flow_gamma : double, Flow_gasConstant : double,
                                    Flow_prandtl : double,
                                    Flow_maxMach : double,
                                    Flow_lengthScale : double,
                                    Flow_constantVisc : double,
                                    Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                                    Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                                    Flow_viscosityModel : SCHEMA.ViscosityModel,
                                    BC_xPosP_inf : double,
                                    Grid_xBnum : int32, Grid_xNum : int32,
                                    Grid_yBnum : int32, Grid_yNum : int32,
                                    Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{cellWidth, rho, velocity, pressure, temperature, rhoVelocity, dudtBoundary, dTdtBoundary}),
  reads(Fluid.{velocityGradientX, velocityGradientY, velocityGradientZ}),
  reads writes(Fluid.{rho_t, rhoVelocity_t, rhoEnergy_t})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = is_xNegGhost(c, Grid_xBnum)
    var xPosGhost = is_xPosGhost(c, Grid_xBnum, Grid_xNum)
    var yNegGhost = is_yNegGhost(c, Grid_yBnum)
    var yPosGhost = is_yPosGhost(c, Grid_yBnum, Grid_yNum)
    var zNegGhost = is_zNegGhost(c, Grid_zBnum)
    var zPosGhost = is_zPosGhost(c, Grid_zBnum, Grid_zNum)
    var ghost_cell = (xNegGhost or xPosGhost or
                      yNegGhost or yPosGhost or
                      zNegGhost or zPosGhost )
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if ghost_cell then
      if NSCBC_inflow_cell then
        -- Add in the x flux using NSCBC
        var c_bnd = int3d(c)
        var c_int = ((c+{1, 0, 0})%Fluid.bounds)

        -- compute amplitudes of waves leaving the domain
        var c_sound = GetSoundSpeed(Fluid[c_bnd].temperature, Flow_gamma, Flow_gasConstant)
        var lambda_1 = Fluid[c_bnd].velocity[0] - c_sound

        var dP_dx = (Fluid[c_int].pressure    - Fluid[c_bnd].pressure)    /  (0.5*Fluid[c_int].cellWidth[0] + 0.5*Fluid[c_bnd].cellWidth[0])
        var du_dx = (Fluid[c_int].velocity[0] - Fluid[c_bnd].velocity[0]) /  (0.5*Fluid[c_int].cellWidth[0] + 0.5*Fluid[c_bnd].cellWidth[0])
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
        -- Add in the x fluxes using NSCBC for outflow
        var c_bnd = int3d(c)
        var c_int = ((c+{-1, 0, 0})%Fluid.bounds)

        var sigma = 0.25 -- Specified constant
        var c_sound = GetSoundSpeed(Fluid[c_bnd].temperature, Flow_gamma, Flow_gasConstant)
        var K = sigma*(1.0-Flow_maxMach*Flow_maxMach)*c_sound/Flow_lengthScale

        var L1 = K*(Fluid[c_bnd].pressure - BC_xPosP_inf)

        var lambda_2 = Fluid[c_bnd].velocity[0]
        var lambda_3 = Fluid[c_bnd].velocity[0]
        var lambda_4 = Fluid[c_bnd].velocity[0]
        var lambda_5 = Fluid[c_bnd].velocity[0] + c_sound

        var drho_dx = (Fluid[c_bnd].rho         - Fluid[c_int].rho        ) / (0.5*Fluid[c_int].cellWidth[0] + 0.5*Fluid[c_bnd].cellWidth[0])
        var dp_dx   = (Fluid[c_bnd].pressure    - Fluid[c_int].pressure   ) / (0.5*Fluid[c_int].cellWidth[0] + 0.5*Fluid[c_bnd].cellWidth[0])
        var du_dx   = (Fluid[c_bnd].velocity[0] - Fluid[c_int].velocity[0]) / (0.5*Fluid[c_int].cellWidth[0] + 0.5*Fluid[c_bnd].cellWidth[0])
        var dv_dx   = (Fluid[c_bnd].velocity[1] - Fluid[c_int].velocity[1]) / (0.5*Fluid[c_int].cellWidth[0] + 0.5*Fluid[c_bnd].cellWidth[0])
        var dw_dx   = (Fluid[c_bnd].velocity[2] - Fluid[c_int].velocity[2]) / (0.5*Fluid[c_int].cellWidth[0] + 0.5*Fluid[c_bnd].cellWidth[0])

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
        var dtau11_dx = (tau11_pos - tau11_neg) / (0.5*Fluid[c_int].cellWidth[0] + 0.5*Fluid[c_bnd].cellWidth[0])
        var dtau21_dx = (tau21_pos - tau21_neg) / (0.5*Fluid[c_int].cellWidth[0] + 0.5*Fluid[c_bnd].cellWidth[0])
        var dtau31_dx = (tau31_pos - tau31_neg) / (0.5*Fluid[c_int].cellWidth[0] + 0.5*Fluid[c_bnd].cellWidth[0])

        -- Stuff for energy equation
        var mu = GetDynamicViscosity(Fluid[c_bnd].temperature,
                                     Flow_constantVisc,
                                     Flow_powerlawTempRef, Flow_powerlawViscRef,
                                     Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef,
                                     Flow_viscosityModel)
        var tau_12 =  mu*( Fluid[c_bnd].velocityGradientY[0] + Fluid[c_bnd].velocityGradientX[1] )
        var tau_13 =  mu*( Fluid[c_bnd].velocityGradientZ[0] + Fluid[c_bnd].velocityGradientX[2] )
        var energy_term_x = (Fluid[c_bnd].velocity[0]*tau11_pos - Fluid[c_int].velocity[0]*tau11_neg) / (0.5*Fluid[c_int].cellWidth[0] + 0.5*Fluid[c_bnd].cellWidth[0]) + Fluid[c_bnd].velocityGradientX[1]*tau_12 + Fluid[c_bnd].velocityGradientX[2]*tau_13

        -- Update the RHS of conservation equations with x fluxes
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
__demand(__leaf, __parallel, __cuda)
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
    var xNegGhost = is_xNegGhost(c, Grid_xBnum)
    var xPosGhost = is_xPosGhost(c, Grid_xBnum, Grid_xNum)
    var yNegGhost = is_yNegGhost(c, Grid_yBnum)
    var yPosGhost = is_yPosGhost(c, Grid_yBnum, Grid_yNum)
    var zNegGhost = is_zNegGhost(c, Grid_zBnum)
    var zPosGhost = is_zPosGhost(c, Grid_zBnum, Grid_zNum)
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

__demand(__leaf, __parallel, __cuda)
task Flow_AddBodyForces(Fluid : region(ispace(int3d), Fluid_columns),
                        config : Config,
                        Grid_xBnum : int32, Grid_xNum : int32,
                        Grid_yBnum : int32, Grid_yNum : int32,
                        Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{rho, velocity}),
  reads writes(Fluid.{rhoEnergy_t, rhoVelocity_t})
do
  var BC_xBCLeft = config.BC.xBCLeft
  var BC_xBCRight = config.BC.xBCRight
  var Flow_bodyForce = config.Flow.bodyForce
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = is_xNegGhost(c, Grid_xBnum)
    var xPosGhost = is_xPosGhost(c, Grid_xBnum, Grid_xNum)
    var yNegGhost = is_yNegGhost(c, Grid_yBnum)
    var yPosGhost = is_yPosGhost(c, Grid_yBnum, Grid_yNum)
    var zNegGhost = is_zNegGhost(c, Grid_zBnum)
    var zPosGhost = is_zPosGhost(c, Grid_zBnum, Grid_zNum)
    var interior = in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)
    var NSCBC_inflow_cell  = ((BC_xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow)   and xNegGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))
    var NSCBC_outflow_cell = ((BC_xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow) and xPosGhost and not (yNegGhost or yPosGhost or zNegGhost or zPosGhost))

    if interior or NSCBC_inflow_cell or NSCBC_outflow_cell then
      [UTIL.emitArrayReduce(3, '+',
         rexpr Fluid[c].rhoVelocity_t end,
         rexpr vs_mul(Flow_bodyForce, Fluid[c].rho) end)];
      Fluid[c].rhoEnergy_t += (Fluid[c].rho*dot(Flow_bodyForce, Fluid[c].velocity))
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_AddVelocity(Fluid : region(ispace(int3d), Fluid_columns),
                      velocity : double[3],
                      Grid_xBnum : int32, Grid_xNum : int32,
                      Grid_yBnum : int32, Grid_yNum : int32,
                      Grid_zBnum : int32, Grid_zNum : int32)
where
  reads writes(Fluid.velocity)
do
  __demand(__openmp)
  for c in Fluid do
    if in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum) then
      [UTIL.emitArrayReduce(3, '+',
         rexpr Fluid[c].velocity end,
         rexpr velocity end)];
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_CalculateAveragePD(Fluid : region(ispace(int3d), Fluid_columns),
                             Grid_xBnum : int32, Grid_xNum : int32,
                             Grid_yBnum : int32, Grid_yNum : int32,
                             Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{cellWidth, pressure, velocityGradientX, velocityGradientY, velocityGradientZ})
do
  var acc = 0.0
  __demand(__openmp)
  for c in Fluid do
    if in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum) then
      var cellVolume = Fluid[c].cellWidth[0]*Fluid[c].cellWidth[1]*Fluid[c].cellWidth[2]
      var divU = Fluid[c].velocityGradientX[0] + Fluid[c].velocityGradientY[1] + Fluid[c].velocityGradientZ[2]
      acc += (divU * Fluid[c].pressure)*cellVolume
    end
  end
  return acc
end

__demand(__leaf, __parallel, __cuda)
task Flow_ResetDissipation(Fluid : region(ispace(int3d), Fluid_columns))
where
  writes(Fluid.dissipation)
do
  __demand(__openmp)
  for c in Fluid do
    Fluid[c].dissipation = 0.0
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_ComputeDissipationX(Fluid : region(ispace(int3d), Fluid_columns),
                              Flow_constantVisc : double,
                              Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                              Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                              Flow_viscosityModel : SCHEMA.ViscosityModel,
                              Grid_xBnum : int32, Grid_xNum : int32,
                              Grid_yBnum : int32, Grid_yNum : int32,
                              Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{centerCoordinates, cellWidth}),
  reads(Fluid.{velocity, temperature, velocityGradientX, velocityGradientY, velocityGradientZ}),
  writes(Fluid.dissipationFlux)
do
  __demand(__openmp)
  for c in Fluid do
    var xNegGhost = is_xNegGhost(c, Grid_xBnum)
    var interior = in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)

    if (interior or xNegGhost) then

      var temperature = Fluid[c].temperature
      var velocity = Fluid[c].velocity
      var velocityGradientX = Fluid[c].velocityGradientX
      var velocityGradientY = Fluid[c].velocityGradientY
      var velocityGradientZ = Fluid[c].velocityGradientZ
      var mu = GetDynamicViscosity(temperature,
                                   Flow_constantVisc,
                                   Flow_powerlawTempRef, Flow_powerlawViscRef,
                                   Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef,
                                   Flow_viscosityModel)
      var xCellWidth = Fluid[c].cellWidth[0]

      var stencil = (c+{1, 0, 0}) % Fluid.bounds
      var temperature_stencil = Fluid[stencil].temperature
      var velocity_stencil = Fluid[stencil].velocity
      var velocityGradientX_stencil = Fluid[stencil].velocityGradientX
      var velocityGradientY_stencil = Fluid[stencil].velocityGradientY
      var velocityGradientZ_stencil = Fluid[stencil].velocityGradientZ
      var mu_stencil = GetDynamicViscosity(temperature_stencil,
                                           Flow_constantVisc,
                                           Flow_powerlawTempRef, Flow_powerlawViscRef,
                                           Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef,
                                           Flow_viscosityModel)
      var xCellWidth_stencil = Fluid[stencil].cellWidth[0]


      var muFace = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                        Fluid[c].centerCoordinates[0],
                                        Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                        mu,
                                        mu_stencil)

      var velocityFace = array(0.0, 0.0, 0.0)
      velocityFace[0] = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                             Fluid[c].centerCoordinates[0],
                                             Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                             velocity[0],
                                             velocity_stencil[0])
      velocityFace[1] = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                             Fluid[c].centerCoordinates[0],
                                             Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                             velocity[1],
                                             velocity_stencil[1])
      velocityFace[2] = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                             Fluid[c].centerCoordinates[0],
                                             Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                             velocity[2],
                                             velocity_stencil[2])


      var velocityX_YFace = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[0],
                                                 Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                                 velocityGradientY[0],
                                                 velocityGradientY_stencil[0])
      var velocityX_ZFace = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[0],
                                                 Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                                 velocityGradientZ[0],
                                                 velocityGradientZ_stencil[0])
      var velocityY_YFace = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[0],
                                                 Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                                 velocityGradientY[1],
                                                 velocityGradientY_stencil[1])
      var velocityZ_ZFace = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[0],
                                                 Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                                 velocityGradientZ[2],
                                                 velocityGradientZ_stencil[2])

      var velocityX_XFace   = 0.0
      var velocityY_XFace   = 0.0
      var velocityZ_XFace   = 0.0
      if xNegGhost or is_xPosGhost(stencil, Grid_xBnum, Grid_xNum)then
        -- do not use velocity gradients in ghost cells for computing the fluxes
        velocityX_XFace   = (velocity_stencil[0] - velocity[0]) / xCellWidth
        velocityY_XFace   = (velocity_stencil[1] - velocity[1]) / xCellWidth
        velocityZ_XFace   = (velocity_stencil[2] - velocity[2]) / xCellWidth
      else
        velocityX_XFace = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                               Fluid[c].centerCoordinates[0],
                                               Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                               velocityGradientX[0],
                                               velocityGradientX_stencil[0])
        velocityY_XFace = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                               Fluid[c].centerCoordinates[0],
                                               Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                               velocityGradientX[1],
                                               velocityGradientX_stencil[1])
        velocityZ_XFace = linear_interpolation(Fluid[c].centerCoordinates[0] + xCellWidth/2.0,
                                               Fluid[c].centerCoordinates[0],
                                               Fluid[c].centerCoordinates[0] + xCellWidth/2.0 + xCellWidth_stencil/2.0,
                                               velocityGradientX[2],
                                               velocityGradientX_stencil[2])
      end

      var sigmaXX = muFace*(4.0*velocityX_XFace-2.0*velocityY_YFace-2.0*velocityZ_ZFace)/3.0
      var sigmaYX = muFace*(velocityY_XFace+velocityX_YFace)
      var sigmaZX = muFace*(velocityZ_XFace+velocityX_ZFace)

      var usigma = velocityFace[0]*sigmaXX + velocityFace[1]*sigmaYX + velocityFace[2]*sigmaZX
      Fluid[c].dissipationFlux = usigma
    end
  end
end

-- CHANGE to reduces
__demand(__leaf, __parallel, __cuda)
task Flow_UpdateDissipationX(Fluid : region(ispace(int3d), Fluid_columns),
                             Grid_xBnum : int32, Grid_xNum : int32,
                             Grid_yBnum : int32, Grid_yNum : int32,
                             Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{cellWidth, dissipationFlux}),
  reads writes(Fluid.dissipation)
do
  __demand(__openmp)
  for c in Fluid do
    if in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum) then
      Fluid[c].dissipation += (Fluid[c].dissipationFlux-Fluid[((c+{-1, 0, 0})%Fluid.bounds)].dissipationFlux)/Fluid[c].cellWidth[0]
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_ComputeDissipationY(Fluid : region(ispace(int3d), Fluid_columns),
                              Flow_constantVisc : double,
                              Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                              Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                              Flow_viscosityModel : SCHEMA.ViscosityModel,
                              Grid_xBnum : int32, Grid_xNum : int32,
                              Grid_yBnum : int32, Grid_yNum : int32,
                              Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{centerCoordinates, cellWidth}),
  reads(Fluid.{velocity, temperature, velocityGradientX, velocityGradientY, velocityGradientZ}),
  writes(Fluid.dissipationFlux)
do
  __demand(__openmp)
  for c in Fluid do
    var yNegGhost = is_yNegGhost(c, Grid_yBnum)
    var interior = in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)

    if (interior or yNegGhost) then

      var temperature = Fluid[c].temperature
      var velocity = Fluid[c].velocity
      var velocityGradientX = Fluid[c].velocityGradientX
      var velocityGradientY = Fluid[c].velocityGradientY
      var velocityGradientZ = Fluid[c].velocityGradientZ
      var mu = GetDynamicViscosity(temperature,
                                   Flow_constantVisc,
                                   Flow_powerlawTempRef, Flow_powerlawViscRef,
                                   Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef,
                                   Flow_viscosityModel)
      var yCellWidth = Fluid[c].cellWidth[1]

      var stencil = (c+{0, 1, 0}) % Fluid.bounds
      var temperature_stencil = Fluid[stencil].temperature
      var velocity_stencil = Fluid[stencil].velocity
      var velocityGradientX_stencil = Fluid[stencil].velocityGradientX
      var velocityGradientY_stencil = Fluid[stencil].velocityGradientY
      var velocityGradientZ_stencil = Fluid[stencil].velocityGradientZ
      var mu_stencil = GetDynamicViscosity(temperature_stencil,
                                           Flow_constantVisc,
                                           Flow_powerlawTempRef, Flow_powerlawViscRef,
                                           Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef,
                                           Flow_viscosityModel)
      var yCellWidth_stencil = Fluid[stencil].cellWidth[1]

      var muFace = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                        Fluid[c].centerCoordinates[1],
                                        Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                        mu,
                                        mu_stencil)

      var velocityFace = array(0.0, 0.0, 0.0)
      velocityFace[0] = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                             Fluid[c].centerCoordinates[1],
                                             Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                             velocity[0],
                                             velocity_stencil[0])
      velocityFace[1] = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                             Fluid[c].centerCoordinates[1],
                                             Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                             velocity[1],
                                             velocity_stencil[1])
      velocityFace[2] = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                             Fluid[c].centerCoordinates[1],
                                             Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                             velocity[2],
                                             velocity_stencil[2])

      var velocityX_XFace = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[1],
                                                 Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                                 velocityGradientX[0],
                                                 velocityGradientX_stencil[0])
      var velocityY_XFace = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[1],
                                                 Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                                 velocityGradientX[1],
                                                 velocityGradientX_stencil[1])
      var velocityY_ZFace = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[1],
                                                 Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                                 velocityGradientZ[1],
                                                 velocityGradientZ_stencil[1])
      var velocityZ_ZFace = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[1],
                                                 Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                                 velocityGradientZ[2],
                                                 velocityGradientZ_stencil[2])

      var velocityX_YFace   = 0.0
      var velocityY_YFace   = 0.0
      var velocityZ_YFace   = 0.0
      if yNegGhost or is_yPosGhost(stencil, Grid_yBnum, Grid_yNum) then
        velocityX_YFace   = (velocity_stencil[0] - velocity[0]) / yCellWidth
        velocityY_YFace   = (velocity_stencil[1] - velocity[1]) / yCellWidth
        velocityZ_YFace   = (velocity_stencil[2] - velocity[2]) / yCellWidth
      else
        velocityX_YFace = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                               Fluid[c].centerCoordinates[1],
                                               Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                               velocityGradientY[0],
                                               velocityGradientY_stencil[0])
        velocityY_YFace = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                               Fluid[c].centerCoordinates[1],
                                               Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                               velocityGradientY[1],
                                               velocityGradientY_stencil[1])
        velocityZ_YFace = linear_interpolation(Fluid[c].centerCoordinates[1] + yCellWidth/2.0,
                                               Fluid[c].centerCoordinates[1],
                                               Fluid[c].centerCoordinates[1] + yCellWidth/2.0 + yCellWidth_stencil/2.0,
                                               velocityGradientY[2],
                                               velocityGradientY_stencil[2])
      end

      var sigmaXY = muFace*(velocityX_YFace+velocityY_XFace)
      var sigmaYY = muFace*(4.0*velocityY_YFace-2.0*velocityX_XFace-2.0*velocityZ_ZFace)/3.0
      var sigmaZY = muFace*(velocityZ_YFace+velocityY_ZFace)

      var usigma = velocityFace[0]*sigmaXY + velocityFace[1]*sigmaYY + velocityFace[2]*sigmaZY

      Fluid[c].dissipationFlux = usigma
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_UpdateDissipationY(Fluid : region(ispace(int3d), Fluid_columns),
                             Grid_xBnum : int32, Grid_xNum : int32,
                             Grid_yBnum : int32, Grid_yNum : int32,
                             Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{cellWidth, dissipationFlux}),
  reads writes(Fluid.dissipation)
do
  __demand(__openmp)
  for c in Fluid do
    if in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum) then
      Fluid[c].dissipation += (Fluid[c].dissipationFlux-Fluid[((c+{0, -1, 0})%Fluid.bounds)].dissipationFlux)/Fluid[c].cellWidth[1]
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_ComputeDissipationZ(Fluid : region(ispace(int3d), Fluid_columns),
                              Flow_constantVisc : double,
                              Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                              Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                              Flow_viscosityModel : SCHEMA.ViscosityModel,
                              Grid_xBnum : int32, Grid_xNum : int32,
                              Grid_yBnum : int32, Grid_yNum : int32,
                              Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{centerCoordinates, cellWidth}),
  reads(Fluid.{velocity, temperature, velocityGradientX, velocityGradientY, velocityGradientZ}),
  writes(Fluid.dissipationFlux)
do
  __demand(__openmp)
  for c in Fluid do
    var zNegGhost = is_zNegGhost(c, Grid_zBnum)
    var interior = in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum)

    if (interior or zNegGhost) then

      var temperature = Fluid[c].temperature
      var velocity = Fluid[c].velocity
      var velocityGradientX = Fluid[c].velocityGradientX
      var velocityGradientY = Fluid[c].velocityGradientY
      var velocityGradientZ = Fluid[c].velocityGradientZ
      var mu = GetDynamicViscosity(temperature,
                                   Flow_constantVisc,
                                   Flow_powerlawTempRef, Flow_powerlawViscRef,
                                   Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef,
                                   Flow_viscosityModel)
      var zCellWidth = Fluid[c].cellWidth[2]



      var stencil = (c+{0, 0, 1}) % Fluid.bounds
      var temperature_stencil = Fluid[stencil].temperature
      var velocity_stencil = Fluid[stencil].velocity
      var velocityGradientX_stencil = Fluid[stencil].velocityGradientX
      var velocityGradientY_stencil = Fluid[stencil].velocityGradientY
      var velocityGradientZ_stencil = Fluid[stencil].velocityGradientZ
      var mu_stencil = GetDynamicViscosity(temperature_stencil,
                                           Flow_constantVisc,
                                           Flow_powerlawTempRef, Flow_powerlawViscRef,
                                           Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef,
                                           Flow_viscosityModel)
      var zCellWidth_stencil = Fluid[stencil].cellWidth[2]

      var muFace = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                        Fluid[c].centerCoordinates[2],
                                        Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                        mu,
                                        mu_stencil)

      var velocityFace = array(0.0, 0.0, 0.0)
      velocityFace[0] = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                             Fluid[c].centerCoordinates[2],
                                             Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                             velocity[0],
                                             velocity_stencil[0])
      velocityFace[1] = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                             Fluid[c].centerCoordinates[2],
                                             Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                             velocity[1],
                                             velocity_stencil[1])
      velocityFace[2] = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                             Fluid[c].centerCoordinates[2],
                                             Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                             velocity[2],
                                             velocity_stencil[2])

      var velocityX_XFace = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[2],
                                                 Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                                 velocityGradientX[0],
                                                 velocityGradientX_stencil[0])
      var velocityY_YFace = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[2],
                                                 Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                                 velocityGradientY[1],
                                                 velocityGradientY_stencil[1])
      var velocityZ_XFace = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[2],
                                                 Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                                 velocityGradientX[2],
                                                 velocityGradientX_stencil[2])
      var velocityZ_YFace = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                                 Fluid[c].centerCoordinates[2],
                                                 Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                                 velocityGradientY[2],
                                                 velocityGradientY_stencil[2])

      var velocityX_ZFace   = 0.0
      var velocityY_ZFace   = 0.0
      var velocityZ_ZFace   = 0.0
      if zNegGhost or is_zPosGhost(stencil, Grid_zBnum, Grid_zNum) then
        velocityX_ZFace   = (velocity_stencil[0] - velocity[0]) / zCellWidth
        velocityY_ZFace   = (velocity_stencil[1] - velocity[1]) / zCellWidth
        velocityZ_ZFace   = (velocity_stencil[2] - velocity[2]) / zCellWidth
      else
        velocityX_ZFace = linear_interpolation(Fluid[c].centerCoordinates[1] + zCellWidth/2.0,
                                               Fluid[c].centerCoordinates[1],
                                               Fluid[c].centerCoordinates[1] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                               velocityGradientZ[0],
                                               velocityGradientZ_stencil[0])
        velocityY_ZFace = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                               Fluid[c].centerCoordinates[2],
                                               Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                               velocityGradientZ[1],
                                               velocityGradientZ_stencil[1])
        velocityZ_ZFace = linear_interpolation(Fluid[c].centerCoordinates[2] + zCellWidth/2.0,
                                               Fluid[c].centerCoordinates[2],
                                               Fluid[c].centerCoordinates[2] + zCellWidth/2.0 + zCellWidth_stencil/2.0,
                                               velocityGradientZ[2],
                                               velocityGradientZ_stencil[2])
      end

      var sigmaXZ = muFace*(velocityX_ZFace+velocityZ_XFace)
      var sigmaYZ = muFace*(velocityY_ZFace+velocityZ_YFace)
      var sigmaZZ = muFace*(4.0*velocityZ_ZFace-2.0*velocityX_XFace-2.0*velocityY_YFace)/3.0

      var usigma = velocityFace[0]*sigmaXZ + velocityFace[1]*sigmaYZ + velocityFace[2]*sigmaZZ

      Fluid[c].dissipationFlux = usigma
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_UpdateDissipationZ(Fluid : region(ispace(int3d), Fluid_columns),
                             Grid_xBnum : int32, Grid_xNum : int32,
                             Grid_yBnum : int32, Grid_yNum : int32,
                             Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{cellWidth, dissipationFlux}),
  reads writes(Fluid.dissipation)
do
  __demand(__openmp)
  for c in Fluid do
    if in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum) then
      Fluid[c].dissipation += (Fluid[c].dissipationFlux-Fluid[((c+{0, 0, -1})%Fluid.bounds)].dissipationFlux)/Fluid[c].cellWidth[2]
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_CalculateAverageDissipation(Fluid : region(ispace(int3d), Fluid_columns),
                                      Grid_xBnum : int32, Grid_xNum : int32,
                                      Grid_yBnum : int32, Grid_yNum : int32,
                                      Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{cellWidth, dissipation})
do
  var acc = 0.0
  __demand(__openmp)
  for c in Fluid do
    if in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum) then
      var cellVolume = Fluid[c].cellWidth[0]*Fluid[c].cellWidth[1]*Fluid[c].cellWidth[2]
      acc += Fluid[c].dissipation*cellVolume
    end
  end
  return acc
end

__demand(__leaf, __parallel, __cuda)
task Flow_CalculateAverageK(Fluid : region(ispace(int3d), Fluid_columns),
                            Grid_xBnum : int32, Grid_xNum : int32,
                            Grid_yBnum : int32, Grid_yNum : int32,
                            Grid_zBnum : int32, Grid_zNum : int32)
where
  reads(Fluid.{cellWidth, rho, velocity})
do
  var acc = 0.0
  __demand(__openmp)
  for c in Fluid do
    if in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum) then
      var cellVolume = Fluid[c].cellWidth[0]*Fluid[c].cellWidth[1]*Fluid[c].cellWidth[2]
      acc += 0.5*Fluid[c].rho*dot(Fluid[c].velocity, Fluid[c].velocity)*cellVolume
    end
  end
  return acc
end

__demand(__leaf, __parallel, __cuda)
task Flow_AddTurbulentSource(Fluid : region(ispace(int3d), Fluid_columns),
                             Flow_averageDissipation : double,
                             Flow_averageK : double,
                             Flow_averagePD : double,
                             Grid_xBnum : int32, Grid_xNum : int32,
                             Grid_yBnum : int32, Grid_yNum : int32,
                             Grid_zBnum : int32, Grid_zNum : int32,
                             config : Config)
where
  reads(Fluid.{cellWidth, rho, velocity}),
  reads writes(Fluid.{rhoVelocity_t, rhoEnergy_t})
do
  var W = Flow_averagePD + Flow_averageDissipation
  var G   = config.Flow.turbForcing.u.HIT.G
  var t_o = config.Flow.turbForcing.u.HIT.t_o
  var K_o = config.Flow.turbForcing.u.HIT.K_o
  var A = (-W-G*(Flow_averageK-K_o)/t_o) / (2.0*Flow_averageK)
  var acc = 0.0
  __demand(__openmp)
  for c in Fluid do
    if in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum) then
      var force = vs_mul(Fluid[c].velocity, Fluid[c].rho*A);
      [UTIL.emitArrayReduce(3, '+',
         rexpr Fluid[c].rhoVelocity_t end,
         rexpr force end)];
      Fluid[c].rhoEnergy_t += dot(force, Fluid[c].velocity)

      var cellVolume = Fluid[c].cellWidth[0]*Fluid[c].cellWidth[1]*Fluid[c].cellWidth[2]
      acc += dot(force, Fluid[c].velocity) * cellVolume
    end
  end
  return acc
end

-- CHANGE to reduces+?
__demand(__leaf, __parallel, __cuda)
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
    if in_interior(c, Grid_xBnum, Grid_xNum, Grid_yBnum, Grid_yNum, Grid_zBnum, Grid_zNum) then
      Fluid[c].rhoEnergy_t += (-Flow_averageFe)
    end
  end
end

-------------------------------------------------------------------------------
-- PARTICLE MOVEMENT
-------------------------------------------------------------------------------

__demand(__leaf, __cuda) -- MANUALLY PARALLELIZED
task Particles_LocateInCells(Particles : region(ispace(int1d), Particles_columns),
                             Grid_xType : SCHEMA.GridType, Grid_yType : SCHEMA.GridType, Grid_zType : SCHEMA.GridType,
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
                                 Grid_xType, Grid_yType, Grid_zType,
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

__demand(__leaf, __cuda) -- MANUALLY PARALLELIZED
task Particles_CheckPartitioning(color : int3d,
                                 Particles : region(ispace(int1d), Particles_columns),
                                 Grid_xBnum : int32, Grid_xNum : int32, NX : int32,
                                 Grid_yBnum : int32, Grid_yNum : int32, NY : int32,
                                 Grid_zBnum : int32, Grid_zNum : int32, NZ : int32)
where
  reads(Particles.{cell, __valid})
do
  var num_invalid = int64(0)
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid and
       color ~= Fluid_elemColor(Particles[p].cell,
                                Grid_xBnum, Grid_xNum, NX,
                                Grid_yBnum, Grid_yNum, NY,
                                Grid_zBnum, Grid_zNum, NZ) then
      num_invalid += 1
    end
  end
  regentlib.assert(num_invalid == 0, 'Invalid particle partitioning')
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

__demand(__leaf, __cuda) -- MANUALLY PARALLELIZED
task TradeQueue_push(partColor : int3d,
                     Particles : region(ispace(int1d), Particles_columns),
                     [tradeQueues],
                     config : Config,
                     Grid_xBnum : int32, Grid_xNum : int32, NX : int32,
                     Grid_yBnum : int32, Grid_yNum : int32, NY : int32,
                     Grid_zBnum : int32, Grid_zNum : int32, NZ : int32)
where
  reads(Particles.[Particles_subStepConserved]),
  reads writes(Particles.{__valid, __xfer_dir, __xfer_slot}),
  [tradeQueues:map(function(queue)
     return Particles_subStepConserved:map(function(fld)
       return regentlib.privilege(regentlib.writes, queue, fld)
     end)
   end):flatten()]
do
  -- Fill in movement direction
  var toTransfer = int64(0)
  __demand(__openmp)
  for i in Particles do
    Particles[i].__xfer_dir = 0
    if Particles[i].__valid then
      var elemColor = Fluid_elemColor(Particles[i].cell,
                                      Grid_xBnum, Grid_xNum, NX,
                                      Grid_yBnum, Grid_yNum, NY,
                                      Grid_zBnum, Grid_zNum, NZ)
      if elemColor ~= partColor then
        toTransfer += 1;
        @ESCAPE for k = 1,26 do @EMIT
          if Particles[i].__xfer_dir == 0 and
             elemColor == (partColor + [colorOffsets[k]] + {NX,NY,NZ}) % {NX,NY,NZ} then
            Particles[i].__xfer_dir = k
            toTransfer += -1
          end
        @TIME end @EPACSE
      end
    end
  end
  [UTIL.emitAssert(
     rexpr toTransfer == 0 end,
     'Sample %d: %ld particle(s) moved past expected stencil',
     rexpr config.Mapping.sampleId end,
     rexpr toTransfer end)];
  var total_xfers = int64(0);
  -- For each movement direction...
  @ESCAPE for k = 1,26 do local queue = tradeQueues[k] @EMIT
    -- Clear the transfer queue
    __demand(__openmp)
    for j in queue do
      queue[j].__valid = false
    end
    -- Assign slots on the transfer queue for moving particles
    var transferred = int64(0)
    __demand(__openmp)
    for i in Particles do
      if Particles[i].__xfer_dir == k then
        Particles[i].__xfer_slot = 1
        transferred += 1
      else
        Particles[i].__xfer_slot = 0
      end
    end
    total_xfers += transferred
    __parallel_prefix(Particles.__xfer_slot, Particles.__xfer_slot, +, 1);
    -- Check that there's enough space in the transfer queue
    [UTIL.emitAssert(
       rexpr transferred <= int64(queue.bounds.hi - queue.bounds.lo + 1) end,
       'Sample %d: Ran out of space in transfer queue',
       rexpr config.Mapping.sampleId end)];
    -- Copy moving particles to the transfer queue
    __demand(__openmp)
    for i in Particles do
      if Particles[i].__xfer_dir == k then
        var j = Particles[i].__xfer_slot - 1 + queue.bounds.lo;
        @ESCAPE for _,fld in ipairs(Particles_subStepConserved) do @EMIT
          queue[j].[fld] = Particles[i].[fld]
        @TIME end @EPACSE
        Particles[i].__valid = false
      end
    end
  @TIME end @EPACSE
  return total_xfers
end

__demand(__leaf, __cuda) -- MANUALLY PARALLELIZED
task TradeQueue_pull(Particles : region(ispace(int1d), Particles_columns),
                     [tradeQueues],
                     config : Config)
where
  reads(Particles.__valid),
  writes(Particles.[Particles_subStepConserved]),
  reads writes(Particles.__xfer_slot),
  [tradeQueues:map(function(queue)
     return Particles_subStepConserved:map(function(fld)
       return regentlib.privilege(regentlib.reads, queue, fld)
     end)
   end):flatten()]
do
  -- Count number of particles coming in from each transfer queue
  var xfer_bounds : int64[27]
  xfer_bounds[0] = 0
  var total_xfers = int64(0);
  @ESCAPE for k = 1,26 do local queue = tradeQueues[k] @EMIT
    __demand(__openmp)
    for j in queue do
      if queue[j].__valid then
        total_xfers += 1
      end
    end
    xfer_bounds[k] = total_xfers
  @TIME end @EPACSE
  -- Number all empty slots on particles sub-region
  var avail_slots = int64(0)
  __demand(__openmp)
  for i in Particles do
    if Particles[i].__valid then
      Particles[i].__xfer_slot = 0
    else
      Particles[i].__xfer_slot = 1
      avail_slots += 1
    end
  end
  __parallel_prefix(Particles.__xfer_slot, Particles.__xfer_slot, +, 1);
  -- Check that there's enough space in the particles sub-region
  [UTIL.emitAssert(
     rexpr total_xfers <= avail_slots end,
     'Sample %d: Not enough space in sub-region for incoming particles',
     rexpr config.Mapping.sampleId end)];
  -- Copy moving particles from the transfer queues
  -- NOTE: This part assumes that transfer queues are filled contiguously.
  @ESCAPE for k = 1,26 do local queue = tradeQueues[k] @EMIT
    var lo = xfer_bounds[k-1]
    var hi = xfer_bounds[k]
    __demand(__openmp)
    for i in Particles do
      if not Particles[i].__valid then
        var j_off = Particles[i].__xfer_slot - 1
        if j_off >= lo and j_off < hi then
          var j = j_off - lo + queue.bounds.lo;
          @ESCAPE for _,fld in ipairs(Particles_subStepConserved) do @EMIT
            Particles[i].[fld] = queue[j].[fld]
          @TIME end @EPACSE
        end
      end
    end
  @TIME end @EPACSE
  return total_xfers
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
  return int64(ceil(
    copiedCells
    / [double](totalCells)
    * (config.Particles.maxNum / config.Particles.parcelSize)
    * config.Particles.maxSkew
  ))
end

__demand(__leaf) -- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task CopyQueue_push(Particles : region(ispace(int1d), Particles_columns),
                    CopyQueue : region(ispace(int1d), CopyQueue_columns),
                    config : Config,
                    copySrc : SCHEMA.Volume,
                    copySrcOrigin : double[3], copyTgtOrigin : double[3])
where
  reads(Particles.[Particles_primitives], Particles.cell),
  writes(CopyQueue.[Particles_primitives])
do
  var p2 = CopyQueue.bounds.lo
  for p1 in Particles do
    if Particles[p1].__valid then
      var cell = Particles[p1].cell
      if  copySrc.fromCell[0] <= cell.x and cell.x <= copySrc.uptoCell[0]
      and copySrc.fromCell[1] <= cell.y and cell.y <= copySrc.uptoCell[1]
      and copySrc.fromCell[2] <= cell.z and cell.z <= copySrc.uptoCell[2] then
        [UTIL.emitAssert(
           rexpr p2 <= CopyQueue.bounds.hi end,
           'Sample %d: Ran out of space in cross-section particles copy queue',
           rexpr config.Mapping.sampleId end)];
        CopyQueue[p2].position = vv_add(copyTgtOrigin, vv_sub(Particles[p1].position, copySrcOrigin))
        CopyQueue[p2].velocity = Particles[p1].velocity
        CopyQueue[p2].temperature = Particles[p1].temperature
        CopyQueue[p2].diameter = Particles[p1].diameter
        CopyQueue[p2].density = Particles[p1].density
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
__demand(__leaf) -- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task CopyQueue_pull(partColor : int3d,
                    Particles : region(ispace(int1d), Particles_columns),
                    CopyQueue : region(ispace(int1d), CopyQueue_columns),
                    config : Config,
                    Grid_xBnum : int32, Grid_yBnum : int32, Grid_zBnum : int32)
where
  reads(CopyQueue.[Particles_primitives], Particles.__valid),
  writes(Particles.[Particles_primitives], Particles.cell)
do
  var acc = int64(0)
  var addedVelocity = config.Particles.feeding.u.Incoming.addedVelocity
  var p1 = Particles.bounds.lo
  for p2 in CopyQueue do
    if CopyQueue[p2].__valid then
      var cell = locate(CopyQueue[p2].position,
                        config.Grid.xType, config.Grid.yType, config.Grid.zType,
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
        [UTIL.emitAssert(
           rexpr p1 <= Particles.bounds.hi end,
           'Sample %d: Ran out of space while copying particles from other section',
           rexpr config.Mapping.sampleId end)];
        Particles[p1].cell = cell
        Particles[p1].position = CopyQueue[p2].position
        Particles[p1].velocity = vv_add(CopyQueue[p2].velocity, addedVelocity)
        Particles[p1].temperature = CopyQueue[p2].temperature
        Particles[p1].diameter = CopyQueue[p2].diameter
        Particles[p1].density = CopyQueue[p2].density
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

__demand(__leaf, __parallel, __cuda)
task Particles_CalcDeltaTerms(Particles : region(ispace(int1d), Particles_columns),
                              Fluid : region(ispace(int3d), Fluid_columns),
                              Flow_constantVisc : double,
                              Flow_powerlawTempRef : double, Flow_powerlawViscRef : double,
                              Flow_sutherlandSRef : double, Flow_sutherlandTempRef : double, Flow_sutherlandViscRef : double,
                              Flow_viscosityModel : SCHEMA.ViscosityModel,
                              Particles_convectiveCoeff : double)
where
  reads(Fluid.{centerCoordinates, velocity, temperature}),
  reads(Particles.{cell, position, velocity, diameter, density, temperature, __valid}),
  writes(Particles.{deltaTemperatureTerm, deltaVelocityOverRelaxationTime})
do
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      var flowVelocity = InterpolateTriVelocity(Particles[p].cell,
                                                Particles[p].position,
                                                Fluid)
      var flowTemperature = InterpolateTriTemp(Particles[p].cell,
                                               Particles[p].position,
                                               Fluid)
      var flowDynamicViscosity = GetDynamicViscosity(flowTemperature,
                                                     Flow_constantVisc,
                                                     Flow_powerlawTempRef, Flow_powerlawViscRef,
                                                     Flow_sutherlandSRef, Flow_sutherlandTempRef, Flow_sutherlandViscRef,
                                                     Flow_viscosityModel)
      var relaxationTime = Particles[p].density * pow(Particles[p].diameter,2.0) / (18.0 * flowDynamicViscosity)
      Particles[p].deltaVelocityOverRelaxationTime = vs_div(vv_sub(flowVelocity, Particles[p].velocity), relaxationTime)
      Particles[p].deltaTemperatureTerm = PI * pow(Particles[p].diameter,2.0) * Particles_convectiveCoeff * (flowTemperature-Particles[p].temperature)
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Particles_AddFlowCoupling(Particles : region(ispace(int1d), Particles_columns),
                               Particles_heatCapacity : double)
where
  reads(Particles.{diameter, density, deltaTemperatureTerm, deltaVelocityOverRelaxationTime, __valid}),
  writes(Particles.{velocity_t, temperature_t})
do
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      Particles[p].velocity_t = Particles[p].deltaVelocityOverRelaxationTime
      Particles[p].temperature_t = Particles[p].deltaTemperatureTerm/(PI*pow(Particles[p].diameter,3.0)/6.0*Particles[p].density*Particles_heatCapacity)
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Particles_AddBodyForces(Particles : region(ispace(int1d), Particles_columns),
                             Particles_bodyForce : double[3])
where
  reads(Particles.__valid),
  reads writes(Particles.velocity_t)
do
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      [UTIL.emitArrayReduce(3, '+',
         rexpr Particles[p].velocity_t end,
         rexpr Particles_bodyForce end)];
    end
  end
end

__demand(__leaf, __cuda) -- MANUALLY PARALLELIZED
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

__demand(__leaf, __parallel, __cuda)
task Radiation_UpdateFieldValues(Radiation : region(ispace(int3d), Radiation_columns),
                                 config : Config,
                                 Radiation_qa : double,
                                 Radiation_qs : double)
where
  reads(Radiation.{cellWidth, acc_d2, acc_d2t4}),
  writes(Radiation.{Ib, sigma})
do
  var Particles_parcelSize = config.Particles.parcelSize
  __demand(__openmp)
  for c in Radiation do
    Radiation[c].sigma = Radiation[c].acc_d2*PI*Particles_parcelSize*(Radiation_qa+Radiation_qs)/(4.0*Radiation[c].cellWidth[0]*Radiation[c].cellWidth[1]*Radiation[c].cellWidth[2])
    if Radiation[c].acc_d2 == 0.0 then
      Radiation[c].Ib = 0.0
    else
      Radiation[c].Ib = (SB*Radiation[c].acc_d2t4)/(PI*Radiation[c].acc_d2)
    end
  end
end

__demand(__leaf, __cuda) -- MANUALLY PARALLELIZED
task Particles_AbsorbRadiationDOM(Particles : region(ispace(int1d), Particles_columns),
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
      var mass = PI*pow(Particles[p].diameter,3.0)/6.0*Particles[p].density
      var t4 = pow(Particles[p].temperature, 4.0)
      var alpha = PI*Radiation_qa*pow(Particles[p].diameter, 2.0)*(Radiation[Fluid[Particles[p].cell].to_Radiation].G-4.0*SB*t4)/4.0
      Particles[p].temperature_t += alpha/(mass*Particles_heatCapacity)
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_AddParticlesCoupling(Particles : region(ispace(int1d), Particles_columns),
                               Fluid : region(ispace(int3d), Fluid_columns),
                               config : Config)
where
  reads(Fluid.{cellWidth}),
  reads(Particles.{cell, diameter, density, deltaTemperatureTerm, deltaVelocityOverRelaxationTime, __valid}),
  reads writes(Fluid.{rhoVelocity_t, rhoEnergy_t})
do
  var Particles_parcelSize = config.Particles.parcelSize
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      var cellVolume = Fluid[Particles[p].cell].cellWidth[0]*Fluid[Particles[p].cell].cellWidth[1]*Fluid[Particles[p].cell].cellWidth[2]
      var mass = PI*pow(Particles[p].diameter,3.0)/6.0*Particles[p].density;
      [UTIL.emitArrayReduce(3, '+',
         rexpr Fluid[Particles[p].cell].rhoVelocity_t end,
         rexpr vs_mul(Particles[p].deltaVelocityOverRelaxationTime, -mass*Particles_parcelSize/cellVolume) end)];
      Fluid[Particles[p].cell].rhoEnergy_t += -Particles_parcelSize*Particles[p].deltaTemperatureTerm/cellVolume
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Flow_UpdateVars(Fluid : region(ispace(int3d), Fluid_columns),
                     Integrator_deltaTime : double,
                     Integrator_stage : int32,
                     config : Config)
where
  reads(Fluid.{rho_old, rhoEnergy_old, rhoVelocity_old}),
  reads(Fluid.{rho_t, rhoEnergy_t, rhoVelocity_t}),
  writes(Fluid.{rho, rhoEnergy, rhoVelocity}),
  reads writes(Fluid.{rho_new, rhoEnergy_new, rhoVelocity_new})
do
  var dt = Integrator_deltaTime;
  @ESCAPE for ORDER = RK_MIN_ORDER,RK_MAX_ORDER do @EMIT
    if config.Integrator.rkOrder == ORDER then
      @ESCAPE for STAGE = 1,ORDER do @EMIT
        if Integrator_stage == STAGE then
          __demand(__openmp)
          for c in Fluid do
            -- Accumulate intermediate values into final values
            Fluid[c].rho_new +=
              Fluid[c].rho_t * [RK_B[ORDER][STAGE]] * dt;
            [UTIL.emitArrayReduce(3, '+',
               rexpr Fluid[c].rhoVelocity_new end,
               rexpr vs_mul(Fluid[c].rhoVelocity_t, [RK_B[ORDER][STAGE]] * dt) end)];
            Fluid[c].rhoEnergy_new +=
              Fluid[c].rhoEnergy_t * [RK_B[ORDER][STAGE]] * dt;
            @ESCAPE if STAGE == ORDER then @EMIT
              -- Set final values
              Fluid[c].rho = Fluid[c].rho_new
              Fluid[c].rhoVelocity = Fluid[c].rhoVelocity_new
              Fluid[c].rhoEnergy = Fluid[c].rhoEnergy_new
            @TIME else @EMIT
              -- Set values for next substep
              Fluid[c].rho = Fluid[c].rho_old +
                Fluid[c].rho_t * [RK_C[ORDER][STAGE]] * dt
              Fluid[c].rhoVelocity = vv_add(Fluid[c].rhoVelocity_old,
                vs_mul(Fluid[c].rhoVelocity_t, [RK_C[ORDER][STAGE]] * dt))
              Fluid[c].rhoEnergy = Fluid[c].rhoEnergy_old +
                Fluid[c].rhoEnergy_t * [RK_C[ORDER][STAGE]] * dt
            @TIME end @EPACSE
          end
        end
      @TIME end @EPACSE
    end
  @TIME end @EPACSE
end

__demand(__leaf, __parallel, __cuda)
task Particles_UpdateVars(Particles : region(ispace(int1d), Particles_columns),
                          Particles_deltaTime : double,
                          Integrator_stage : int32,
                          config : Config)
where
  reads(Particles.{position_old, velocity_old, temperature_old}),
  reads(Particles.{velocity, velocity_t, temperature_t}),
  reads(Particles.__valid),
  writes(Particles.{position, temperature, velocity}),
  reads writes(Particles.{position_new, temperature_new, velocity_new})
do
  var dt = Particles_deltaTime;
  @ESCAPE for ORDER = RK_MIN_ORDER,RK_MAX_ORDER do @EMIT
    if config.Integrator.rkOrder == ORDER then
      @ESCAPE for STAGE = 1,ORDER do @EMIT
        if Integrator_stage == STAGE then
          __demand(__openmp)
          for p in Particles do
            if Particles[p].__valid then
              -- Accumulate intermediate values into final values
              [UTIL.emitArrayReduce(3, '+',
                 rexpr Particles[p].position_new end,
                 rexpr vs_mul(Particles[p].velocity, [RK_B[ORDER][STAGE]] * dt) end)];
              [UTIL.emitArrayReduce(3, '+',
                 rexpr Particles[p].velocity_new end,
                 rexpr vs_mul(Particles[p].velocity_t, [RK_B[ORDER][STAGE]] * dt) end)];
              Particles[p].temperature_new +=
                Particles[p].temperature_t * [RK_B[ORDER][STAGE]] * dt;
              @ESCAPE if STAGE == ORDER then @EMIT
                -- Set final values
                Particles[p].position = Particles[p].position_new
                Particles[p].velocity = Particles[p].velocity_new
                Particles[p].temperature = Particles[p].temperature_new
              @TIME else @EMIT
                -- Set values for next substep
                Particles[p].position = vv_add(Particles[p].position_old,
                  vs_mul(Particles[p].velocity, [RK_C[ORDER][STAGE]] * dt))
                Particles[p].velocity = vv_add(Particles[p].velocity_old,
                  vs_mul(Particles[p].velocity_t, [RK_C[ORDER][STAGE]] * dt))
                Particles[p].temperature = Particles[p].temperature_old +
                  Particles[p].temperature_t * [RK_C[ORDER][STAGE]] * dt
              @TIME end @EPACSE
            end
          end
        end
      @TIME end @EPACSE
    end
  @TIME end @EPACSE
end

__demand(__leaf) -- MANUALLY PARALLELIZED, NO CUDA, NO OPENMP
task Particles_HandleCollisions(Particles : region(ispace(int1d), Particles_columns),
                                config : Config,
                                Particles_deltaTime : double,
                                Particles_restitutionCoeff : double )
-- This is an adaption of collisionPrt routine of the Soleil-MPI version
-- TODO: search box implementation
where
  reads(Particles.{position_old, diameter, density, __valid}),
  reads writes(Particles.{position, velocity})
do
  var Particles_parcelSize = config.Particles.parcelSize
  for p1 in Particles do
    if Particles[p1].__valid then
      for p2 in Particles do
        if Particles[p2].__valid and p1 < p2 then

          -- Relative position of particles
          var x = Particles[p2].position[0] - Particles[p1].position[0]
          var y = Particles[p2].position[1] - Particles[p1].position[1]
          var z = Particles[p2].position[2] - Particles[p1].position[2]

          -- Old relative position of particles
          var xold = Particles[p2].position_old[0] - Particles[p1].position_old[0]
          var yold = Particles[p2].position_old[1] - Particles[p1].position_old[1]
          var zold = Particles[p2].position_old[2] - Particles[p1].position_old[2]


          -- Relative velocity
          var ux = (x-xold)/Particles_deltaTime
          var uy = (y-yold)/Particles_deltaTime
          var uz = (z-zold)/Particles_deltaTime

          -- Relevant scalar products
          var x_scal_u = xold*ux + yold*uy + zold*uz
          var x_scal_x = xold*xold + yold*yold + zold*zold
          var u_scal_u = ux*ux + uy*uy + uz*uz

          -- Critical distance
          var dcrit = 0.5 * sqrt(Particles_parcelSize) * ( Particles[p1].diameter + Particles[p2].diameter )

          -- Checking if particles are getting away from each other
          if x_scal_u<0.0 then

            -- Checking if particles are in a collision path
            var det = x_scal_u*x_scal_u - u_scal_u*(x_scal_x - dcrit*dcrit)
            if det>0.0 then

              -- Checking if collision occurs in this time step
              var timecol = ( -x_scal_u - sqrt(det) ) / u_scal_u
              if (timecol>0.0 and timecol<Particles_deltaTime) then

                -- We do have a collision


                -- Mass ratio of particles
                var mr = (Particles[p2].density * Particles[p2].diameter * Particles[p2].diameter * Particles[p2].diameter)
                mr = mr/ (Particles[p1].density * Particles[p1].diameter * Particles[p1].diameter * Particles[p1].diameter)

                -- Change of velocity and particle location after impact
                -- Note: for now particle restitution coeff is the same for all particles ?
                var du = ( 1.0 + min( Particles_restitutionCoeff, Particles_restitutionCoeff ) ) / (1.0 + mr)*x_scal_u/x_scal_x
                var dx = du * ( Particles_deltaTime - timecol )

                -- Update velocities

                Particles[p1].velocity[0] = Particles[p1].velocity[0] + du*xold*mr
                Particles[p1].velocity[1] = Particles[p1].velocity[1] + du*yold*mr
                Particles[p1].velocity[2] = Particles[p1].velocity[2] + du*zold*mr

                Particles[p2].velocity[0] = Particles[p2].velocity[0] - du*xold
                Particles[p2].velocity[1] = Particles[p2].velocity[1] - du*yold
                Particles[p2].velocity[2] = Particles[p2].velocity[2] - du*zold

                -- Update positions
                Particles[p1].position[0] = Particles[p1].position[0] + dx*xold*mr
                Particles[p1].position[1] = Particles[p1].position[1] + dx*yold*mr
                Particles[p1].position[2] = Particles[p1].position[2] + dx*zold*mr

                Particles[p2].position[0] = Particles[p2].position[0] - dx*xold
                Particles[p2].position[1] = Particles[p2].position[1] - dx*yold
                Particles[p2].position[2] = Particles[p2].position[2] - dx*zold

              end

            end
          end

        end
      end
    end
  end
end

__demand(__leaf, __parallel, __cuda)
task Particles_UpdateAuxiliary(Particles : region(ispace(int1d), Particles_columns),
                               BC_xBCParticles : SCHEMA.ParticlesBC,
                               BC_yBCParticles : SCHEMA.ParticlesBC,
                               BC_zBCParticles : SCHEMA.ParticlesBC,
                               Grid_xOrigin : double, Grid_xWidth : double,
                               Grid_yOrigin : double, Grid_yWidth : double,
                               Grid_zOrigin : double, Grid_zWidth : double,
                               Particles_restitutionCoeff : double)
where
  reads(Particles.__valid),
  reads writes(Particles.{position, velocity, velocity_t})
do
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      if (Particles[p].position[0]<Grid_xOrigin) then
        if BC_xBCParticles == SCHEMA.ParticlesBC_Periodic then
          Particles[p].position[0] += Grid_xWidth
        elseif BC_xBCParticles == SCHEMA.ParticlesBC_Bounce then
          Particles[p].position[0] = Grid_xOrigin
          var impulse = ((-(1.0+Particles_restitutionCoeff))*Particles[p].velocity[0])
          if (impulse<=0.0) then
            Particles[p].velocity[0] += impulse
          end
          Particles[p].velocity_t[0] max= 0.0
        else -- BC_xBCParticles == SCHEMA.ParticlesBC_Disappear
          -- Do nothing, let out-of-bounds particles get deleted
        end
      end
      if (Particles[p].position[0]>(Grid_xOrigin+Grid_xWidth)) then
        if BC_xBCParticles == SCHEMA.ParticlesBC_Periodic then
          Particles[p].position[0] += -Grid_xWidth
        elseif BC_xBCParticles == SCHEMA.ParticlesBC_Bounce then
          Particles[p].position[0] = (Grid_xOrigin+Grid_xWidth)
          var impulse = ((-(1.0+Particles_restitutionCoeff))*Particles[p].velocity[0])
          if (impulse>=0.0) then
            Particles[p].velocity[0] += impulse
          end
          Particles[p].velocity_t[0] min= 0.0
        else -- BC_xBCParticles == SCHEMA.ParticlesBC_Disappear
          -- Do nothing, let out-of-bounds particles get deleted
        end
      end
      if (Particles[p].position[1]<Grid_yOrigin) then
        if BC_yBCParticles == SCHEMA.ParticlesBC_Periodic then
          Particles[p].position[1] += Grid_yWidth
        elseif BC_yBCParticles == SCHEMA.ParticlesBC_Bounce then
          Particles[p].position[1] = Grid_yOrigin
          var impulse = ((-(1.0+Particles_restitutionCoeff))*Particles[p].velocity[1])
          if (impulse<=0.0) then
            Particles[p].velocity[1] += impulse
          end
          Particles[p].velocity_t[1] max= 0.0
        else -- BC_yBCParticles == SCHEMA.ParticlesBC_Disappear
          -- Do nothing, let out-of-bounds particles get deleted
        end
      end
      if (Particles[p].position[1]>(Grid_yOrigin+Grid_yWidth)) then
        if BC_yBCParticles == SCHEMA.ParticlesBC_Periodic then
          Particles[p].position[1] += -Grid_yWidth
        elseif BC_yBCParticles == SCHEMA.ParticlesBC_Bounce then
          Particles[p].position[1] = (Grid_yOrigin+Grid_yWidth)
          var impulse = ((-(1.0+Particles_restitutionCoeff))*Particles[p].velocity[1])
          if (impulse>=0.0) then
            Particles[p].velocity[1] += impulse
          end
          Particles[p].velocity_t[1] min= 0.0
        else -- BC_yBCParticles == SCHEMA.ParticlesBC_Disappear
          -- Do nothing, let out-of-bounds particles get deleted
        end
      end
      if (Particles[p].position[2]<Grid_zOrigin) then
        if BC_zBCParticles == SCHEMA.ParticlesBC_Periodic then
          Particles[p].position[2] += Grid_zWidth
        elseif BC_zBCParticles == SCHEMA.ParticlesBC_Bounce then
          Particles[p].position[2] = Grid_zOrigin
          var impulse = ((-(1.0+Particles_restitutionCoeff))*Particles[p].velocity[2])
          if (impulse<=0.0) then
            Particles[p].velocity[2] += impulse
          end
          Particles[p].velocity_t[2] max= 0.0
        else -- BC_zBCParticles == SCHEMA.ParticlesBC_Disappear
          -- Do nothing, let out-of-bounds particles get deleted
        end
      end
      if (Particles[p].position[2]>(Grid_zOrigin+Grid_zWidth)) then
        if BC_zBCParticles == SCHEMA.ParticlesBC_Periodic then
          Particles[p].position[2] += -Grid_zWidth
        elseif BC_zBCParticles == SCHEMA.ParticlesBC_Bounce then
          Particles[p].position[2] = (Grid_zOrigin+Grid_zWidth)
          var impulse = ((-(1.0+Particles_restitutionCoeff))*Particles[p].velocity[2])
          if (impulse>=0.0) then
            Particles[p].velocity[2] += impulse
          end
          Particles[p].velocity_t[2] min= 0.0
        else -- BC_zBCParticles == SCHEMA.ParticlesBC_Disappear
          -- Do nothing, let out-of-bounds particles get deleted
        end
      end
    end
  end
end

__demand(__leaf, __cuda) -- MANUALLY PARALLELIZED
task Particles_DeleteEscapingParticles(Particles : region(ispace(int1d), Particles_columns),
                                       x_min : double, x_max : double,
                                       y_min : double, y_max : double,
                                       z_min : double, z_max : double)
where
  reads(Particles.position),
  reads writes(Particles.__valid)
do
  var acc = int64(0)
  __demand(__openmp)
  for p in Particles do
    if Particles[p].__valid then
      var pos = Particles[p].position
      if pos[0] < x_min
      or pos[1] < y_min
      or pos[2] < z_min
      or pos[0] > x_max
      or pos[1] > y_max
      or pos[2] > z_max then
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

local function mkInstance(config) local INSTANCE = {}

  -----------------------------------------------------------------------------
  -- Symbols shared between quotes
  -----------------------------------------------------------------------------

  local SYMBOLS = UTIL.SymbolFactory()

  local DEBUG_COPYING = SYMBOLS:scalarVar(bool)
  local startTime = SYMBOLS:scalarVar(int)
  local Grid = {
    volume = SYMBOLS:scalarVar(double),
    xBnum = SYMBOLS:scalarVar(int),
    yBnum = SYMBOLS:scalarVar(int),
    zBnum = SYMBOLS:scalarVar(int),
    xRealOrigin = SYMBOLS:scalarVar(double),
    yRealOrigin = SYMBOLS:scalarVar(double),
    zRealOrigin = SYMBOLS:scalarVar(double),
    xRealMax = SYMBOLS:scalarVar(double),
    yRealMax = SYMBOLS:scalarVar(double),
    zRealMax = SYMBOLS:scalarVar(double),
    particleCopyOrigin = SYMBOLS:scalarVar(double[3]),
  }
  local BC = {
    xPosSign = SYMBOLS:scalarVar(double[3]),
    xNegSign = SYMBOLS:scalarVar(double[3]),
    xPosVelocity = SYMBOLS:scalarVar(double[3]),
    xNegVelocity = SYMBOLS:scalarVar(double[3]),
    xPosTemperature = SYMBOLS:scalarVar(double),
    xNegTemperature = SYMBOLS:scalarVar(double),
    yPosSign = SYMBOLS:scalarVar(double[3]),
    yNegSign = SYMBOLS:scalarVar(double[3]),
    yPosVelocity = SYMBOLS:scalarVar(double[3]),
    yNegVelocity = SYMBOLS:scalarVar(double[3]),
    yPosTemperature = SYMBOLS:scalarVar(double),
    yNegTemperature = SYMBOLS:scalarVar(double),
    zPosSign = SYMBOLS:scalarVar(double[3]),
    zNegSign = SYMBOLS:scalarVar(double[3]),
    zPosVelocity = SYMBOLS:scalarVar(double[3]),
    zNegVelocity = SYMBOLS:scalarVar(double[3]),
    zPosTemperature = SYMBOLS:scalarVar(double),
    zNegTemperature = SYMBOLS:scalarVar(double),
    xBCParticles = SYMBOLS:scalarVar(SCHEMA.ParticlesBC),
    yBCParticles = SYMBOLS:scalarVar(SCHEMA.ParticlesBC),
    zBCParticles = SYMBOLS:scalarVar(SCHEMA.ParticlesBC),
  }
  local numTiles = SYMBOLS:scalarVar(int)
  local Integrator = {
    deltaTime = SYMBOLS:scalarVar(double),
    simTime = SYMBOLS:scalarVar(double),
    timeStep = SYMBOLS:scalarVar(int),
    exitCond = SYMBOLS:scalarVar(bool),
  }
  local Particles_number = SYMBOLS:scalarVar(int64)
  local Flow_averagePressure = SYMBOLS:scalarVar(double)
  local Flow_averageTemperature = SYMBOLS:scalarVar(double)
  local Flow_averageKineticEnergy = SYMBOLS:scalarVar(double)
  local Particles_averageTemperature = SYMBOLS:scalarVar(double)

  local is_Fluid = SYMBOLS:ispaceVar(int3d)
  local Fluid = SYMBOLS:regionVar(is_Fluid, Fluid_columns)
  local Fluid_copy = SYMBOLS:regionVar(is_Fluid, Fluid_columns)
  local is_Particles = SYMBOLS:ispaceVar(int1d)
  local Particles = SYMBOLS:regionVar(is_Particles, Particles_columns)
  local Particles_copy = SYMBOLS:regionVar(is_Particles, Particles_columns)
  local is_TradeQueue = UTIL.generate(26, function()
    return SYMBOLS:ispaceVar(int1d)
  end)
  local TradeQueue = is_TradeQueue:map(function(is)
    return SYMBOLS:regionVar(is, TradeQueue_columns)
  end)
  local is_Radiation = SYMBOLS:ispaceVar(int3d)
  local Radiation = SYMBOLS:regionVar(is_Radiation, Radiation_columns)

  local tiles = SYMBOLS:ispaceVar(int3d)
  local p_Fluid = SYMBOLS:partitionVar(Fluid, tiles)
  local p_Fluid_copy = SYMBOLS:partitionVar(Fluid_copy, tiles)
  local p_Particles = SYMBOLS:partitionVar(Particles, tiles)
  local p_Particles_copy = SYMBOLS:partitionVar(Particles_copy, tiles)
  local p_TradeQueue_bySrc = TradeQueue:map(function(r)
    return SYMBOLS:partitionVar(r, tiles)
  end)
  local p_TradeQueue_byDst = TradeQueue:map(function(r)
    return SYMBOLS:partitionVar(r, tiles)
  end)
  local p_Radiation = SYMBOLS:partitionVar(Radiation, tiles)

  -----------------------------------------------------------------------------
  -- Exported symbols
  -----------------------------------------------------------------------------

  INSTANCE.DEBUG_COPYING = DEBUG_COPYING
  INSTANCE.Grid = Grid
  INSTANCE.Integrator = Integrator
  INSTANCE.Flow_averagePressure = Flow_averagePressure
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
  -- Imported instances
  -----------------------------------------------------------------------------

  local DOM_INST = DOM.mkInstance(SYMBOLS, config, tiles, p_Radiation)

  -----------------------------------------------------------------------------
  -- Symbol declaration & initialization
  -----------------------------------------------------------------------------

  function INSTANCE.DeclSymbols() return rquote

    ---------------------------------------------------------------------------
    -- Environment options
    ---------------------------------------------------------------------------

    var [DEBUG_COPYING] = false
    if C.getenv('DEBUG_COPYING') ~= [&int8](0) and
       C.strcmp(C.getenv('DEBUG_COPYING'), '1') == 0 then
      DEBUG_COPYING = true
    end

    ---------------------------------------------------------------------------
    -- Preparation
    ---------------------------------------------------------------------------

    -- Start timer
    var [startTime] = C.legion_get_current_time_in_micros() / 1000;

    -- Write console header
    Console_WriteHeader(0, config)

    -- Write probe file headers
    var probeId = 0
    while probeId < config.IO.probes.length do
      Probe_WriteHeader(0, config, probeId)
      probeId += 1
    end

    ---------------------------------------------------------------------------
    -- Declare & initialize state variables
    ---------------------------------------------------------------------------

    var [Grid.volume] = 0.0

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

    var [Grid.xRealOrigin] = config.Grid.origin[0]
    var [Grid.yRealOrigin] = config.Grid.origin[1]
    var [Grid.zRealOrigin] = config.Grid.origin[2]

    var [Grid.xRealMax] = config.Grid.origin[0] + config.Grid.xWidth
    var [Grid.yRealMax] = config.Grid.origin[1] + config.Grid.yWidth
    var [Grid.zRealMax] = config.Grid.origin[2] + config.Grid.zWidth

    var [Grid.particleCopyOrigin] = array(0.0, 0.0, 0.0)

    var [numTiles] = config.Mapping.tiles[0] * config.Mapping.tiles[1] * config.Mapping.tiles[2]

    var [Integrator.exitCond] = true
    var [Integrator.simTime] = config.Integrator.startTime
    var [Integrator.timeStep] = config.Integrator.startIter
    var [Integrator.deltaTime] = config.Integrator.fixedDeltaTime
    regentlib.assert(
      config.Integrator.rkOrder >= RK_MIN_ORDER and
      config.Integrator.rkOrder <= RK_MAX_ORDER,
      'Unsupported RK integration scheme')

    var [Particles_number] = int64(0)

    var [Flow_averagePressure] = 0.0
    var [Flow_averageTemperature] = 0.0
    var [Flow_averageKineticEnergy] = 0.0
    var [Particles_averageTemperature] = 0.0

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
    var [is_Fluid] = ispace(int3d, {x = config.Grid.xNum + 2*Grid.xBnum,
                                    y = config.Grid.yNum + 2*Grid.yBnum,
                                    z = config.Grid.zNum + 2*Grid.zBnum})
    var [Fluid] = region(is_Fluid, Fluid_columns);
    [UTIL.emitRegionTagAttach(Fluid, MAPPER.SAMPLE_ID_TAG, sampleId, int)];
    var [Fluid_copy] = region(is_Fluid, Fluid_columns);
    [UTIL.emitRegionTagAttach(Fluid_copy, MAPPER.SAMPLE_ID_TAG, sampleId, int)];

    -- Create Particles Regions
    regentlib.assert((config.Particles.maxNum / config.Particles.parcelSize) % numTiles == 0,
                     'Uneven partitioning of particles')
    var maxParticlesPerTile = config.Particles.maxNum / config.Particles.parcelSize / numTiles
    if numTiles > 1 then
      maxParticlesPerTile =
        int64(ceil(maxParticlesPerTile * config.Particles.maxSkew))
    end
    var [is_Particles] = ispace(int1d, maxParticlesPerTile * numTiles)
    var [Particles] = region(is_Particles, Particles_columns);
    [UTIL.emitRegionTagAttach(Particles, MAPPER.SAMPLE_ID_TAG, sampleId, int)];
    var [Particles_copy] = region(is_Particles, Particles_columns);
    [UTIL.emitRegionTagAttach(Particles_copy, MAPPER.SAMPLE_ID_TAG, sampleId, int)];
    @ESCAPE for k = 1,26 do @EMIT
      -- Make tradequeues smaller for diagonal movement
      var off = [colorOffsets[k]]
      var num_dirs = off.x*off.x + off.y*off.y + off.z*off.z
      var escapeRatio = 1.0
      for i = 0, num_dirs do
        escapeRatio *= config.Particles.escapeRatioPerDir
      end
      var [is_TradeQueue[k]] = ispace(int1d, int64(ceil(escapeRatio * maxParticlesPerTile) * numTiles))
      var [TradeQueue[k]] = region([is_TradeQueue[k]], TradeQueue_columns);
      [UTIL.emitRegionTagAttach(TradeQueue[k], MAPPER.SAMPLE_ID_TAG, sampleId, int)];
    @TIME end @EPACSE

    -- Create Radiation Regions
    var rad_x = config.Mapping.tiles[0]
    var rad_y = config.Mapping.tiles[1]
    var rad_z = config.Mapping.tiles[2]
    if config.Radiation.type == SCHEMA.RadiationModel_DOM then
      rad_x = config.Radiation.u.DOM.xNum
      rad_y = config.Radiation.u.DOM.yNum
      rad_z = config.Radiation.u.DOM.zNum
    end
    var [is_Radiation] = ispace(int3d, {x = rad_x, y = rad_y, z = rad_z})
    var [Radiation] = region(is_Radiation, Radiation_columns);
    [UTIL.emitRegionTagAttach(Radiation, MAPPER.SAMPLE_ID_TAG, sampleId, int)];

    -- Partitioning domain
    var [tiles] = ispace(int3d, {config.Mapping.tiles[0],config.Mapping.tiles[1],config.Mapping.tiles[2]})

    -- Fluid Partitioning
    var [p_Fluid] =
      [UTIL.mkPartitionByTile(int3d, int3d, Fluid_columns)]
      (Fluid, tiles, int3d{Grid.xBnum,Grid.yBnum,Grid.zBnum}, int3d{0,0,0})
    var [p_Fluid_copy] =
      [UTIL.mkPartitionByTile(int3d, int3d, Fluid_columns)]
      (Fluid_copy, tiles, int3d{Grid.xBnum,Grid.yBnum,Grid.zBnum}, int3d{0,0,0})

    -- Particles Partitioning
    var [p_Particles] =
      [UTIL.mkPartitionByTile(int1d, int3d, Particles_columns)]
      (Particles, tiles, 0, int3d{0,0,0})
    var [p_Particles_copy] =
      [UTIL.mkPartitionByTile(int1d, int3d, Particles_columns)]
      (Particles_copy, tiles, 0, int3d{0,0,0});
    @ESCAPE for k = 1,26 do @EMIT
      var [p_TradeQueue_bySrc[k]] =
        [UTIL.mkPartitionByTile(int1d, int3d, TradeQueue_columns)]
        ([TradeQueue[k]], tiles, 0, int3d{0,0,0});
      var [p_TradeQueue_byDst[k]] =
        [UTIL.mkPartitionByTile(int1d, int3d, TradeQueue_columns)]
        ([TradeQueue[k]], tiles, 0, [colorOffsets[k]]);
    @TIME end @EPACSE

    -- Radiation Partitioning
    var [p_Radiation] =
      [UTIL.mkPartitionByTile(int3d, int3d, Radiation_columns)]
      (Radiation, tiles, int3d{0,0,0}, int3d{0,0,0});

    ---------------------------------------------------------------------------
    -- DOM code declarations
    ---------------------------------------------------------------------------

    [DOM_INST.DeclSymbols()];

  end end -- DeclSymbols

  -----------------------------------------------------------------------------
  -- Region initialization
  -----------------------------------------------------------------------------

  local function SyncConservedPrimitive() return rquote

    -- Use the interior conserved values (and BC settings) to update primitives everywhere
    Flow_UpdateAuxiliaryVelocity(Fluid,
                                 config,
                                 config.Flow.constantVisc,
                                 config.Flow.powerlawTempRef, config.Flow.powerlawViscRef,
                                 config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef,
                                 config.Flow.viscosityModel,
                                 Grid.xBnum, config.Grid.xNum,
                                 Grid.yBnum, config.Grid.yNum,
                                 Grid.zBnum, config.Grid.zNum)
    for c in tiles do
      Flow_UpdateGhostVelocity(p_Fluid[c],
                               config,
                               BC.xNegVelocity, BC.xPosVelocity, BC.xNegSign, BC.xPosSign,
                               BC.yNegVelocity, BC.yPosVelocity, BC.yNegSign, BC.yPosSign,
                               BC.zNegVelocity, BC.zPosVelocity, BC.zNegSign, BC.zPosSign,
                               Grid.xBnum, config.Grid.xNum,
                               Grid.yBnum, config.Grid.yNum,
                               Grid.zBnum, config.Grid.zNum)
    end
    Flow_UpdateAuxiliaryThermodynamics(Fluid,
                                       config,
                                       config.Flow.gamma,
                                       config.Flow.gasConstant,
                                       Grid.xBnum, config.Grid.xNum,
                                       Grid.yBnum, config.Grid.yNum,
                                       Grid.zBnum, config.Grid.zNum)
    for c in tiles do
      Flow_UpdateGhostThermodynamics(p_Fluid[c],
                                     config,
                                     config.Flow.gamma,
                                     config.Flow.gasConstant,
                                     BC.xNegTemperature, BC.xPosTemperature,
                                     BC.yNegTemperature, BC.yPosTemperature,
                                     BC.zNegTemperature, BC.zPosTemperature,
                                     Grid.xBnum, config.Grid.xNum,
                                     Grid.yBnum, config.Grid.yNum,
                                     Grid.zBnum, config.Grid.zNum)
    end

    -- Compute the conserved values in the ghost cells
    Flow_UpdateGhostConserved(Fluid,
                              config,
                              Grid.xBnum, config.Grid.xNum,
                              Grid.yBnum, config.Grid.yNum,
                              Grid.zBnum, config.Grid.zNum)

  end end -- SyncConservedPrimitive

  function INSTANCE.InitRegions() return rquote

    -- initialize base region contents (dummy values & connectivity info)
    if config.Particles.maxNum > 0 then
      Particles_initValidField(Particles)
    end
    if config.Radiation.type == SCHEMA.RadiationModel_DOM then
      Flow_SetCoarseningField(Fluid,
                              Grid.xBnum, config.Grid.xNum,
                              Grid.yBnum, config.Grid.yNum,
                              Grid.zBnum, config.Grid.zNum,
                              config.Radiation.u.DOM.xNum,
                              config.Radiation.u.DOM.yNum,
                              config.Radiation.u.DOM.zNum)
    end
    Flow_InitializeCell(Fluid)
    for c in tiles do
      Flow_InitializeGeometry(p_Fluid[c],
                              config.Grid.xType, config.Grid.yType, config.Grid.zType,
                              Grid.xBnum, config.Grid.xNum, config.Grid.origin[0], config.Grid.xWidth,
                              Grid.yBnum, config.Grid.yNum, config.Grid.origin[1], config.Grid.yWidth,
                              Grid.zBnum, config.Grid.zNum, config.Grid.origin[2], config.Grid.zWidth)
    end
    Grid.volume += Flow_CalculateInteriorVolume(Fluid,
                                                Grid.xBnum, config.Grid.xNum,
                                                Grid.yBnum, config.Grid.yNum,
                                                Grid.zBnum, config.Grid.zNum)
    Grid.xRealOrigin min= Flow_FindRealOriginX(Fluid, config.Grid.origin[0])
    Grid.yRealOrigin min= Flow_FindRealOriginY(Fluid, config.Grid.origin[1])
    Grid.zRealOrigin min= Flow_FindRealOriginZ(Fluid, config.Grid.origin[2])
    Grid.xRealMax max= Flow_FindRealMaxX(Fluid, config.Grid.origin[0] + config.Grid.xWidth)
    Grid.yRealMax max= Flow_FindRealMaxY(Fluid, config.Grid.origin[1] + config.Grid.yWidth)
    Grid.zRealMax max= Flow_FindRealMaxZ(Fluid, config.Grid.origin[2] + config.Grid.zWidth)

    -- initialize fluid proper
    if config.Flow.initCase == SCHEMA.FlowInitCase_Uniform then
      Flow_InitializeUniform(Fluid, config.Flow.initParams)
    elseif config.Flow.initCase == SCHEMA.FlowInitCase_Random then
      for c in tiles do
        Flow_InitializeRandom(p_Fluid[c], config.Flow.initParams)
      end
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
      for c in tiles do
        Flow_InitializePerturbed(p_Fluid[c], config.Flow.initParams)
      end
    elseif config.Flow.initCase == SCHEMA.FlowInitCase_Restart then
      HDF_FLUID.load(0, tiles, config.Flow.restartDir, Fluid, Fluid_copy, p_Fluid, p_Fluid_copy)
    else regentlib.assert(false, 'Unhandled case in switch') end

    -- initialize ghost cells to their specified values in NSCBC case
    if ((config.BC.xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow) and (config.BC.xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow)) then
      for c in tiles do
        Flow_InitializeGhostNSCBC(p_Fluid[c],
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
    end

    -- update all cells based on initialized primitive values
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
    [SyncConservedPrimitive(config)];

    -- Initialize particles
    if config.Particles.maxNum > 0 then
      if config.Particles.initCase == SCHEMA.ParticlesInitCase_Random then
        regentlib.assert(numTiles == 1, 'Random particle initialization will only work on 1 tile')
        regentlib.assert((config.Particles.initNum / config.Particles.parcelSize) % numTiles == 0,
                         'Uneven partitioning of particles')
        regentlib.assert(config.Particles.initNum <= config.Particles.maxNum,
                         'Not enough space for initial number of particles')
        for c in tiles do
          Particles_InitializeRandom(c,
                                     p_Particles[c],
                                     p_Fluid[c],
                                     config,
                                     Grid.volume,
                                     Grid.xBnum, Grid.yBnum, Grid.zBnum)
        end
      elseif config.Particles.initCase == SCHEMA.ParticlesInitCase_Restart then
        HDF_PARTICLES.load(0, tiles, config.Particles.restartDir, Particles, Particles_copy, p_Particles, p_Particles_copy)
        for c in tiles do
          Particles_LocateInCells(p_Particles[c],
                                  config.Grid.xType, config.Grid.yType, config.Grid.zType,
                                  Grid.xBnum, config.Grid.xNum, config.Grid.origin[0], config.Grid.xWidth,
                                  Grid.yBnum, config.Grid.yNum, config.Grid.origin[1], config.Grid.yWidth,
                                  Grid.zBnum, config.Grid.zNum, config.Grid.origin[2], config.Grid.zWidth)
        end
      elseif config.Particles.initCase == SCHEMA.ParticlesInitCase_Uniform then
        regentlib.assert((config.Particles.initNum / config.Particles.parcelSize) % numTiles == 0,
                         'Uneven partitioning of particles')
        regentlib.assert(config.Particles.initNum <= config.Particles.maxNum,
                         'Not enough space for initial number of particles')
        for c in tiles do
          Particles_InitializeUniform(p_Particles[c],
                                      p_Fluid[c],
                                      config,
                                      Grid.xBnum, Grid.yBnum, Grid.zBnum)
        end
      else regentlib.assert(false, 'Unhandled case in switch') end
      for c in tiles do
        Particles_CheckPartitioning(c,
                                    p_Particles[c],
                                    Grid.xBnum, config.Grid.xNum, config.Mapping.tiles[0],
                                    Grid.yBnum, config.Grid.yNum, config.Mapping.tiles[1],
                                    Grid.zBnum, config.Grid.zNum, config.Mapping.tiles[2])
      end
      Particles_number += Particles_CalculateNumber(Particles)
    end

    -- Initialize radiation
    if config.Radiation.type == SCHEMA.RadiationModel_OFF then
      -- Do nothing
    elseif config.Radiation.type == SCHEMA.RadiationModel_Algebraic then
      -- Do nothing
    elseif config.Radiation.type == SCHEMA.RadiationModel_DOM then
      [DOM_INST.InitRegions()];
      for c in tiles do
        Radiation_InitializeGeometry(p_Radiation[c],
                                     p_Fluid[c],
                                     config.Grid.xNum,
                                     config.Grid.yNum,
                                     config.Grid.zNum,
                                     config.Radiation.u.DOM.xNum,
                                     config.Radiation.u.DOM.yNum,
                                     config.Radiation.u.DOM.zNum)
      end
    else regentlib.assert(false, 'Unhandled case in switch') end

  end end -- InitRegions

  -----------------------------------------------------------------------------
  -- Init for copy particles in 2 section simulations
  -----------------------------------------------------------------------------

  function INSTANCE.InitParticleCopyOrigin(fromCell) return rquote

    var copyOriginX = 0.0
    copyOriginX += Flow_FindCellOriginX(Fluid, fromCell)
    var copyOriginY = 0.0
    copyOriginY += Flow_FindCellOriginY(Fluid, fromCell)
    var copyOriginZ = 0.0
    copyOriginZ += Flow_FindCellOriginZ(Fluid, fromCell)

    Grid.particleCopyOrigin = array(copyOriginX, copyOriginY, copyOriginZ)

  end end -- InitParticleCopyOrigin

  -----------------------------------------------------------------------------
  -- Main time-step loop header
  -----------------------------------------------------------------------------

  function INSTANCE.MainLoopHeader() return rquote

    -- Calculate exit condition
    Integrator.exitCond =
      Integrator.timeStep >= config.Integrator.maxIter

    -- Determine time step size
    if config.Integrator.cfl > 0.0 then
      var Integrator_maxConvectiveSpectralRadius = 0.0
      var Integrator_maxViscousSpectralRadius = 0.0
      var Integrator_maxHeatConductionSpectralRadius = 0.0
      Integrator_maxConvectiveSpectralRadius max=
        Flow_CalculateConvectiveSpectralRadius(Fluid,
                                               config.Flow.gamma, config.Flow.gasConstant)
      Integrator_maxViscousSpectralRadius max=
        Flow_CalculateViscousSpectralRadius(Fluid,
                                            config.Flow.constantVisc,
                                            config.Flow.powerlawTempRef, config.Flow.powerlawViscRef,
                                            config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef,
                                            config.Flow.viscosityModel)
      Integrator_maxHeatConductionSpectralRadius max=
        Flow_CalculateHeatConductionSpectralRadius(Fluid,
                                                   config.Flow.constantVisc,
                                                   config.Flow.gamma, config.Flow.gasConstant,
                                                   config.Flow.powerlawTempRef, config.Flow.powerlawViscRef,
                                                   config.Flow.prandtl,
                                                   config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef,
                                                   config.Flow.viscosityModel)
      Integrator.deltaTime = (config.Integrator.cfl/max(Integrator_maxConvectiveSpectralRadius, max(Integrator_maxViscousSpectralRadius, Integrator_maxHeatConductionSpectralRadius)))
    end

  end end -- MainLoopHeader

  -----------------------------------------------------------------------------
  -- Per-time-step I/O
  -----------------------------------------------------------------------------

  function INSTANCE.DumpHDF(nameFmt, ...) local args = terralib.newlist{...} return rquote

    var dirname = [&int8](C.malloc(256))
    C.snprintf(dirname, 256, ['%s/fluid_'..nameFmt], config.Mapping.outDir, [args])
    var _1 = IO_CreateDir(0, dirname)
    _1 = HDF_FLUID.dump(_1, tiles, dirname, Fluid, Fluid_copy, p_Fluid, p_Fluid_copy)
    _1 = HDF_FLUID.write.timeStep(_1, tiles, dirname, Fluid, p_Fluid, Integrator.timeStep)
    _1 = HDF_FLUID.write.simTime(_1, tiles, dirname, Fluid, p_Fluid, Integrator.simTime)
    C.snprintf(dirname, 256, ['%s/particles_'..nameFmt], config.Mapping.outDir, [args])
    var _2 = IO_CreateDir(0, dirname)
    _2 = HDF_PARTICLES.dump(_2, tiles, dirname, Particles, Particles_copy, p_Particles, p_Particles_copy)
    _2 = HDF_PARTICLES.write.timeStep(_2, tiles, dirname, Particles, p_Particles, Integrator.timeStep)
    _2 = HDF_PARTICLES.write.simTime(_2, tiles, dirname, Particles, p_Particles, Integrator.simTime)
    C.free(dirname)

  end end -- DumpHDF

  function INSTANCE.PerformIO() return rquote

    -- Write to console
    Flow_averagePressure = 0.0
    Flow_averageTemperature = 0.0
    Flow_averageKineticEnergy = 0.0
    Particles_averageTemperature = 0.0
    Flow_averagePressure += Flow_CalculateAveragePressure(Fluid,
                                                          Grid.xBnum, config.Grid.xNum,
                                                          Grid.yBnum, config.Grid.yNum,
                                                          Grid.zBnum, config.Grid.zNum)
    Flow_averageTemperature += Flow_CalculateAverageTemperature(Fluid,
                                                                Grid.xBnum, config.Grid.xNum,
                                                                Grid.yBnum, config.Grid.yNum,
                                                                Grid.zBnum, config.Grid.zNum)
    Flow_averageKineticEnergy += Flow_CalculateAverageKineticEnergy(Fluid,
                                                                    Grid.xBnum, config.Grid.xNum,
                                                                    Grid.yBnum, config.Grid.yNum,
                                                                    Grid.zBnum, config.Grid.zNum)
    if config.Particles.maxNum > 0 then
      Particles_averageTemperature += Particles_IntegrateQuantities(Particles)
    end
    Flow_averagePressure /= Grid.volume
    Flow_averageTemperature /= Grid.volume
    Flow_averageKineticEnergy /= Grid.volume
    Particles_averageTemperature /= Particles_number
    Console_Write(config,
                  Integrator.timeStep,
                  Integrator.simTime,
                  startTime,
                  Integrator.deltaTime,
                  Flow_averagePressure,
                  Flow_averageTemperature,
                  Flow_averageKineticEnergy,
                  Particles_number,
                  Particles_averageTemperature)

    -- Write probe files
    for i = 0,config.IO.probes.length do
      var probe = config.IO.probes.values[i]
      var totalCells =
        (probe.uptoCell[0] - probe.fromCell[0] + 1) *
        (probe.uptoCell[1] - probe.fromCell[1] + 1) *
        (probe.uptoCell[2] - probe.fromCell[2] + 1)
      var avgFluidT = 0.0
      avgFluidT += Probe_AvgFluidT(Fluid, probe, totalCells)
      var totalParticles = int64(0)
      var avgParticleT = 0.0
      var avgCellOfParticleT = 0.0
      if config.Particles.maxNum > 0 then
        totalParticles += Probe_CountParticles(Particles, probe)
        avgParticleT += Probe_AvgParticleT(Particles, probe, totalParticles)
        avgCellOfParticleT += Probe_AvgCellOfParticleT(Fluid, Particles, probe, totalParticles)
      end
      Probe_Write(0, config, i, Integrator.timeStep, avgFluidT, avgParticleT, avgCellOfParticleT)
    end

    -- Dump restart files
    if config.IO.wrtRestart then
      if Integrator.exitCond or Integrator.timeStep % config.IO.restartEveryTimeSteps == 0 then
        [INSTANCE.DumpHDF('iter%010d', Integrator.timeStep)];
      end
    end

  end end -- PerformIO

  -----------------------------------------------------------------------------
  -- Main time-step loop body
  -----------------------------------------------------------------------------

  function INSTANCE.MainLoopBody(incoming, CopyQueue) return rquote

    -- Process incoming values from other section
    if incoming then
      if DEBUG_COPYING then
        [INSTANCE.DumpHDF('precopy%010d', Integrator.timeStep)];
      end
      -- Feed fluid
      [SyncConservedPrimitive(config)];
      -- Feed particles
      if config.Particles.maxNum > 0 then
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
      end
      if DEBUG_COPYING then
        [INSTANCE.DumpHDF('postcopy%010d', Integrator.timeStep)];
      end
    end

    -- Set iteration-specific fields that persist across RK sub-steps
    Flow_InitializeTemporaries(Fluid)
    if config.Particles.maxNum > 0 and Integrator.timeStep % config.Particles.staggerFactor == 0 then
      Particles_InitializeTemporaries(Particles)
    end

    -- RK sub-time-stepping loop
    var Integrator_time_old = Integrator.simTime
    for Integrator_stage = 1,config.Integrator.rkOrder+1 do

      -- Compute velocity gradients
      Flow_ComputeGradients(Fluid,
                            config,
                            Grid.xBnum, config.Grid.xNum,
                            Grid.yBnum, config.Grid.yNum,
                            Grid.zBnum, config.Grid.zNum)

      for c in tiles do
        Flow_UpdateGhostGradients(p_Fluid[c],
                                  config,
                                  BC.xNegSign, BC.yNegSign, BC.zNegSign,
                                  BC.xPosSign, BC.yPosSign, BC.zPosSign,
                                  Grid.xBnum, config.Grid.xNum,
                                  Grid.yBnum, config.Grid.yNum,
                                  Grid.zBnum, config.Grid.zNum)
      end

      -- Compute fluxes
      Flow_GetFluxX(Fluid,
                    config,
                    config.Flow.constantVisc,
                    config.Flow.gamma, config.Flow.gasConstant,
                    config.Flow.powerlawTempRef, config.Flow.powerlawViscRef,
                    config.Flow.prandtl,
                    config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef,
                    config.Flow.viscosityModel,
                    Grid.xBnum, config.Grid.xNum,
                    Grid.yBnum, config.Grid.yNum,
                    Grid.zBnum, config.Grid.zNum)
      Flow_GetFluxY(Fluid,
                    config,
                    config.Flow.constantVisc,
                    config.Flow.gamma, config.Flow.gasConstant,
                    config.Flow.powerlawTempRef, config.Flow.powerlawViscRef,
                    config.Flow.prandtl,
                    config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef,
                    config.Flow.viscosityModel,
                    Grid.xBnum, config.Grid.xNum,
                    Grid.yBnum, config.Grid.yNum,
                    Grid.zBnum, config.Grid.zNum)
      Flow_GetFluxZ(Fluid,
                    config,
                    config.Flow.constantVisc,
                    config.Flow.gamma, config.Flow.gasConstant,
                    config.Flow.powerlawTempRef, config.Flow.powerlawViscRef,
                    config.Flow.prandtl,
                    config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef,
                    config.Flow.viscosityModel,
                    Grid.xBnum, config.Grid.xNum,
                    Grid.yBnum, config.Grid.yNum,
                    Grid.zBnum, config.Grid.zNum)

      -- Initialize conserved derivatives to 0
      Flow_InitializeTimeDerivatives(Fluid)

      -- Add body forces
      Flow_AddBodyForces(Fluid,
                         config,
                         Grid.xBnum, config.Grid.xNum,
                         Grid.yBnum, config.Grid.yNum,
                         Grid.zBnum, config.Grid.zNum)

      -- Add turbulent forcing
      if config.Flow.turbForcing.type == SCHEMA.TurbForcingModel_HIT then
        Flow_AddVelocity(Fluid,
                         vs_mul(config.Flow.turbForcing.u.HIT.meanVelocity, -1.0),
                         Grid.xBnum, config.Grid.xNum,
                         Grid.yBnum, config.Grid.yNum,
                         Grid.zBnum, config.Grid.zNum)
        var Flow_averageDissipation = 0.0
        var Flow_averageFe = 0.0
        var Flow_averageK = 0.0
        var Flow_averagePD = 0.0

        Flow_averagePD += Flow_CalculateAveragePD(Fluid,
                                                  Grid.xBnum, config.Grid.xNum,
                                                  Grid.yBnum, config.Grid.yNum,
                                                  Grid.zBnum, config.Grid.zNum)
        Flow_averagePD /= Grid.volume

        Flow_ResetDissipation(Fluid)

        Flow_ComputeDissipationX(Fluid,
                                 config.Flow.constantVisc,
                                 config.Flow.powerlawTempRef, config.Flow.powerlawViscRef,
                                 config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef,
                                 config.Flow.viscosityModel,
                                 Grid.xBnum, config.Grid.xNum,
                                 Grid.yBnum, config.Grid.yNum,
                                 Grid.zBnum, config.Grid.zNum)
        Flow_UpdateDissipationX(Fluid,
                                Grid.xBnum, config.Grid.xNum,
                                Grid.yBnum, config.Grid.yNum,
                                Grid.zBnum, config.Grid.zNum)

        Flow_ComputeDissipationY(Fluid,
                                 config.Flow.constantVisc,
                                 config.Flow.powerlawTempRef, config.Flow.powerlawViscRef,
                                 config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef,
                                 config.Flow.viscosityModel,
                                 Grid.xBnum, config.Grid.xNum,
                                 Grid.yBnum, config.Grid.yNum,
                                 Grid.zBnum, config.Grid.zNum)
        Flow_UpdateDissipationY(Fluid,
                                Grid.xBnum, config.Grid.xNum,
                                Grid.yBnum, config.Grid.yNum,
                                Grid.zBnum, config.Grid.zNum)

        Flow_ComputeDissipationZ(Fluid,
                                 config.Flow.constantVisc,
                                 config.Flow.powerlawTempRef, config.Flow.powerlawViscRef,
                                 config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef,
                                 config.Flow.viscosityModel,
                                 Grid.xBnum, config.Grid.xNum,
                                 Grid.yBnum, config.Grid.yNum,
                                 Grid.zBnum, config.Grid.zNum)
        Flow_UpdateDissipationZ(Fluid,
                                Grid.xBnum, config.Grid.xNum,
                                Grid.yBnum, config.Grid.yNum,
                                Grid.zBnum, config.Grid.zNum)

        Flow_averageDissipation += Flow_CalculateAverageDissipation(Fluid,
                                                                    Grid.xBnum, config.Grid.xNum,
                                                                    Grid.yBnum, config.Grid.yNum,
                                                                    Grid.zBnum, config.Grid.zNum)
        Flow_averageDissipation /= Grid.volume

        Flow_averageK += Flow_CalculateAverageK(Fluid,
                                                Grid.xBnum, config.Grid.xNum,
                                                Grid.yBnum, config.Grid.yNum,
                                                Grid.zBnum, config.Grid.zNum)
        Flow_averageK /= Grid.volume

        Flow_averageFe += Flow_AddTurbulentSource(Fluid,
                                                  Flow_averageDissipation,
                                                  Flow_averageK,
                                                  Flow_averagePD,
                                                  Grid.xBnum, config.Grid.xNum,
                                                  Grid.yBnum, config.Grid.yNum,
                                                  Grid.zBnum, config.Grid.zNum,
                                                  config)
        Flow_averageFe /= Grid.volume

        Flow_AdjustTurbulentSource(Fluid,
                                   Flow_averageFe,
                                   Grid.xBnum, config.Grid.xNum,
                                   Grid.yBnum, config.Grid.yNum,
                                   Grid.zBnum, config.Grid.zNum)
        Flow_AddVelocity(Fluid,
                         config.Flow.turbForcing.u.HIT.meanVelocity,
                         Grid.xBnum, config.Grid.xNum,
                         Grid.yBnum, config.Grid.yNum,
                         Grid.zBnum, config.Grid.zNum)
      end

      -- Particles & radiation solve
      if config.Particles.maxNum > 0 and (Integrator.timeStep % config.Particles.staggerFactor == 0 or Integrator.timeStep == config.Integrator.startIter) then
        Particles_CalcDeltaTerms(Particles,
                                 Fluid,
                                 config.Flow.constantVisc,
                                 config.Flow.powerlawTempRef, config.Flow.powerlawViscRef,
                                 config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef,
                                 config.Flow.viscosityModel,
                                 config.Particles.convectiveCoeff)
      end
      if config.Particles.maxNum > 0 and Integrator.timeStep % config.Particles.staggerFactor == 0 then
        -- Add fluid forces to particles
        Particles_AddFlowCoupling(Particles, config.Particles.heatCapacity)
        Particles_AddBodyForces(Particles, config.Particles.bodyForce)
        -- Add radiation
        if config.Radiation.type == SCHEMA.RadiationModel_OFF then
          -- Do nothing
        elseif config.Radiation.type == SCHEMA.RadiationModel_Algebraic then
          Particles_AbsorbRadiationAlgebraic(Particles, config)
        elseif config.Radiation.type == SCHEMA.RadiationModel_DOM then
          fill(Radiation.acc_d2, 0.0)
          fill(Radiation.acc_d2t4, 0.0)
          for c in tiles do
            Radiation_AccumulateParticleValues(p_Particles[c], p_Fluid[c], p_Radiation[c])
          end
          Radiation_UpdateFieldValues(Radiation,
                                      config,
                                      config.Radiation.u.DOM.qa,
                                      config.Radiation.u.DOM.qs);
          [DOM_INST.ComputeRadiationField()];
          for c in tiles do
            Particles_AbsorbRadiationDOM(p_Particles[c],
                                         p_Fluid[c],
                                         p_Radiation[c],
                                         config.Particles.heatCapacity,
                                         config.Radiation.u.DOM.qa)
          end
        else regentlib.assert(false, 'Unhandled case in switch') end
      end

      -- Add particle forces to fluid
      if config.Particles.maxNum > 0 and config.Particles.twoWayCoupled then
        Flow_AddParticlesCoupling(Particles, Fluid, config)
      end

      -- Use fluxes to update conserved value derivatives
      Flow_UpdateUsingFluxZ(Fluid,
                            config,
                            Grid.xBnum, config.Grid.xNum,
                            Grid.yBnum, config.Grid.yNum,
                            Grid.zBnum, config.Grid.zNum)
      Flow_UpdateUsingFluxY(Fluid,
                            config,
                            Grid.xBnum, config.Grid.xNum,
                            Grid.yBnum, config.Grid.yNum,
                            Grid.zBnum, config.Grid.zNum)
      Flow_UpdateUsingFluxX(Fluid,
                            config,
                            Grid.xBnum, config.Grid.xNum,
                            Grid.yBnum, config.Grid.yNum,
                            Grid.zBnum, config.Grid.zNum)
      if ((config.BC.xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow) and (config.BC.xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow)) then
        var Flow_maxMach = -math.huge
        Flow_maxMach max= Flow_CalculateMaxMachNumber(Fluid,
                                                      config,
                                                      config.Flow.gamma, config.Flow.gasConstant,
                                                      Grid.xBnum, config.Grid.xNum,
                                                      Grid.yBnum, config.Grid.yNum,
                                                      Grid.zBnum, config.Grid.zNum)
        -- XXX: Assumes that inflow out/flow only in the x direction
        var Flow_lengthScale = config.Grid.xWidth
        for c in tiles do
          Flow_UpdateUsingFluxGhostNSCBC(p_Fluid[c],
                                         config,
                                         config.Flow.gamma, config.Flow.gasConstant,
                                         config.Flow.prandtl,
                                         Flow_maxMach,
                                         Flow_lengthScale,
                                         config.Flow.constantVisc,
                                         config.Flow.powerlawTempRef, config.Flow.powerlawViscRef,
                                         config.Flow.sutherlandSRef, config.Flow.sutherlandTempRef, config.Flow.sutherlandViscRef,
                                         config.Flow.viscosityModel,
                                         config.BC.xBCRightP_inf,
                                         Grid.xBnum, config.Grid.xNum,
                                         Grid.yBnum, config.Grid.yNum,
                                         Grid.zBnum, config.Grid.zNum)
        end
      end

      -- Time step
      Flow_UpdateVars(Fluid, Integrator.deltaTime, Integrator_stage, config)
      if config.Particles.maxNum > 0 and Integrator.timeStep % config.Particles.staggerFactor == 0 then
        Particles_UpdateVars(Particles,
                             Integrator.deltaTime * config.Particles.staggerFactor,
                             Integrator_stage,
                             config)
      end

      -- Impose desired mean velocity
      if config.Flow.turbForcing.type == SCHEMA.TurbForcingModel_HIT then
        var Flow_averageVelocityX = 0.0
        var Flow_averageVelocityY = 0.0
        var Flow_averageVelocityZ = 0.0
        Flow_averageVelocityX += Flow_CalculateAverageVelocityX(Fluid,
                                                                Grid.xBnum, config.Grid.xNum,
                                                                Grid.yBnum, config.Grid.yNum,
                                                                Grid.zBnum, config.Grid.zNum)
        Flow_averageVelocityY += Flow_CalculateAverageVelocityY(Fluid,
                                                                Grid.xBnum, config.Grid.xNum,
                                                                Grid.yBnum, config.Grid.yNum,
                                                                Grid.zBnum, config.Grid.zNum)
        Flow_averageVelocityZ += Flow_CalculateAverageVelocityZ(Fluid,
                                                                Grid.xBnum, config.Grid.xNum,
                                                                Grid.yBnum, config.Grid.yNum,
                                                                Grid.zBnum, config.Grid.zNum)
        Flow_averageVelocityX /= Grid.volume
        Flow_averageVelocityY /= Grid.volume
        Flow_averageVelocityZ /= Grid.volume
        Flow_AdjustAverageVelocity(Fluid,
                                   config,
                                   Flow_averageVelocityX,
                                   Flow_averageVelocityY,
                                   Flow_averageVelocityZ,
                                   Grid.xBnum, config.Grid.xNum,
                                   Grid.yBnum, config.Grid.yNum,
                                   Grid.zBnum, config.Grid.zNum)
      end

      -- Update all cell values (conserved & primitive) based on updated interior conserved
      [SyncConservedPrimitive(config)];

      -- Particle movement post-processing
      if config.Particles.maxNum > 0 and Integrator.timeStep % config.Particles.staggerFactor == 0 then
        -- Handle particle collisions
        -- TODO: Collisions across tiles are not handled.
        if config.Particles.collisions and Integrator_stage == config.Integrator.rkOrder then
          for c in tiles do
            Particles_HandleCollisions(p_Particles[c],
                                       config,
                                       Integrator.deltaTime * config.Particles.staggerFactor,
                                       config.Particles.restitutionCoeff)
          end
        end
        -- Handle particle boundary conditions
        Particles_UpdateAuxiliary(Particles,
                                  BC.xBCParticles,
                                  BC.yBCParticles,
                                  BC.zBCParticles,
                                  config.Grid.origin[0], config.Grid.xWidth,
                                  config.Grid.origin[1], config.Grid.yWidth,
                                  config.Grid.origin[2], config.Grid.zWidth,
                                  config.Particles.restitutionCoeff)
        for c in tiles do
          Particles_number +=
            Particles_DeleteEscapingParticles(p_Particles[c],
                                              Grid.xRealOrigin, Grid.xRealMax,
                                              Grid.yRealOrigin, Grid.yRealMax,
                                              Grid.zRealOrigin, Grid.zRealMax)
        end

        -- Move particles to new partitions
        for c in tiles do
          Particles_LocateInCells(p_Particles[c],
                                  config.Grid.xType, config.Grid.yType, config.Grid.zType,
                                  Grid.xBnum, config.Grid.xNum, config.Grid.origin[0], config.Grid.xWidth,
                                  Grid.yBnum, config.Grid.yNum, config.Grid.origin[1], config.Grid.yWidth,
                                  Grid.zBnum, config.Grid.zNum, config.Grid.origin[2], config.Grid.zWidth)
        end
        if numTiles > 1 then
          var totalPushed = int64(0)
          for c in tiles do
            totalPushed +=
              TradeQueue_push(c,
                              p_Particles[c],
                              [UTIL.range(1,26):map(function(k) return rexpr
                                 [p_TradeQueue_bySrc[k]][c]
                               end end)],
                              config,
                              Grid.xBnum, config.Grid.xNum, config.Mapping.tiles[0],
                              Grid.yBnum, config.Grid.yNum, config.Mapping.tiles[1],
                              Grid.zBnum, config.Grid.zNum, config.Mapping.tiles[2])
          end
          var totalPulled = int64(0)
          for c in tiles do
            totalPulled +=
              TradeQueue_pull(p_Particles[c],
                              [UTIL.range(1,26):map(function(k) return rexpr
                                 [p_TradeQueue_byDst[k]][c]
                               end end)],
                              config)
          end
          regentlib.assert(totalPushed == totalPulled, 'Internal error in particle trading')
        end
      end

      -- Advance the time for the next sub-step
      @ESCAPE for ORDER = RK_MIN_ORDER,RK_MAX_ORDER do @EMIT
        if config.Integrator.rkOrder == ORDER then
          @ESCAPE for STAGE = 1,ORDER do @EMIT
            if Integrator_stage == STAGE then
              @ESCAPE if STAGE == ORDER then @EMIT
                Integrator.simTime = Integrator_time_old + Integrator.deltaTime
              @TIME else @EMIT
                Integrator.simTime = Integrator_time_old + [RK_B[ORDER][STAGE]] * Integrator.deltaTime
              @TIME end @EPACSE
            end
          @TIME end @EPACSE
        end
      @TIME end @EPACSE

    end -- RK sub-time-stepping

    -- Update time derivatives at boundary for NSCBC
    if config.BC.xBCLeft == SCHEMA.FlowBC_NSCBC_SubsonicInflow and config.BC.xBCRight == SCHEMA.FlowBC_NSCBC_SubsonicOutflow then
      Flow_UpdateNSCBCGhostCellTimeDerivatives(Fluid,
                                               config,
                                               Grid.xBnum, config.Grid.xNum,
                                               Grid.yBnum, config.Grid.yNum,
                                               Grid.zBnum, config.Grid.zNum,
                                               Integrator.deltaTime)
    end

    if incoming then
      if DEBUG_COPYING then
        [INSTANCE.DumpHDF('postiter%010d', Integrator.timeStep)];
      end
    end

    Integrator.timeStep += 1

  end end -- MainLoopBody

  -----------------------------------------------------------------------------
  -- Cleanup code
  -----------------------------------------------------------------------------

  function INSTANCE.Cleanup() return rquote

    -- Wait for everything above to finish
    __fence(__execution, __block)

  end end -- Cleanup

return INSTANCE end -- mkInstance

-------------------------------------------------------------------------------
-- TOP-LEVEL INTERFACE
-------------------------------------------------------------------------------

local function parallelizeFor(sim, stmts)
  return rquote
    __parallelize_with
      sim.p_Fluid, sim.p_Particles, sim.p_Radiation, sim.tiles,
      image(sim.Fluid, sim.p_Particles, [sim.Particles].cell) <= sim.p_Fluid
    do [stmts] end
  end
end

__forbid(__optimize) __demand(__inner, __replicable)
task workSingle(config : Config)
  @ESCAPE local SIM = mkInstance(config) @EMIT
    [SIM.DeclSymbols()];
    var is_FakeCopyQueue = ispace(int1d, 0)
    var FakeCopyQueue = region(is_FakeCopyQueue, CopyQueue_columns);
    [UTIL.emitRegionTagAttach(FakeCopyQueue, MAPPER.SAMPLE_ID_TAG, -1, int)];
    [parallelizeFor(SIM, rquote
      [SIM.InitRegions()];
      while true do
        [SIM.MainLoopHeader()];
        [SIM.PerformIO()];
        if SIM.Integrator.exitCond then
          break
        end
        [SIM.MainLoopBody(rexpr false end, FakeCopyQueue)];
      end
    end)];
    [SIM.Cleanup()];
  @TIME @EPACSE
end

__forbid(__optimize) __demand(__inner, __replicable)
task workDual(mc : MultiConfig)
  -- Create simulation instances
  @ESCAPE
    local SIM0 = mkInstance(rexpr mc.configs[0] end)
    local SIM1 = mkInstance(rexpr mc.configs[1] end)
    @EMIT
      -- Declare symbols
      [SIM0.DeclSymbols()];
      [SIM1.DeclSymbols()];
      var is_FakeCopyQueue = ispace(int1d, 0)
      var FakeCopyQueue = region(is_FakeCopyQueue, CopyQueue_columns);
      [UTIL.emitRegionTagAttach(FakeCopyQueue, MAPPER.SAMPLE_ID_TAG, -1, int)];
      var CopyQueue_size = int64(0)
      var coloring_CopyQueue = C.legion_domain_point_coloring_create()
      for c in SIM0.tiles do
        var partSize = CopyQueue_partSize(SIM0.p_Fluid[c].bounds,
                                          mc.configs[0],
                                          mc.copySrc)
        C.legion_domain_point_coloring_color_domain(
          coloring_CopyQueue, c, rect1d{CopyQueue_size,CopyQueue_size+partSize-1})
        CopyQueue_size += partSize
      end
      var is_CopyQueue = ispace(int1d, CopyQueue_size)
      var CopyQueue = region(is_CopyQueue, CopyQueue_columns);
      [UTIL.emitRegionTagAttach(CopyQueue, MAPPER.SAMPLE_ID_TAG, rexpr mc.configs[0].Mapping.sampleId end, int)];
      var p_CopyQueue = partition(disjoint, CopyQueue, coloring_CopyQueue, SIM0.tiles)
      C.legion_domain_point_coloring_destroy(coloring_CopyQueue)
      -- Check 2-section configuration
      regentlib.assert(
        SIM0.Integrator.simTime == SIM1.Integrator.simTime and
        SIM0.Integrator.timeStep == SIM1.Integrator.timeStep,
        'Coupled sections disagree on starting time')
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
      [parallelizeFor(SIM0, SIM0.InitRegions())];
      [parallelizeFor(SIM1, SIM1.InitRegions())];
      [parallelizeFor(SIM0, SIM0.InitParticleCopyOrigin(rexpr mc.copySrc.fromCell end))];
      [parallelizeFor(SIM1, SIM1.InitParticleCopyOrigin(rexpr mc.copyTgt.fromCell end))];
      var srcOrigin = int3d{mc.copySrc.fromCell[0], mc.copySrc.fromCell[1], mc.copySrc.fromCell[2]}
      var tgtOrigin = int3d{mc.copyTgt.fromCell[0], mc.copyTgt.fromCell[1], mc.copyTgt.fromCell[2]}
      var srcColoring = C.legion_domain_point_coloring_create()
      for c in SIM1.tiles do
        var tgtRect = intersection(SIM1.p_Fluid[c].bounds, mc.copyTgt)
        if rectSize(tgtRect) > 0 then
          var srcRect = rect3d{lo = tgtRect.lo - tgtOrigin + srcOrigin,
                               hi = tgtRect.hi - tgtOrigin + srcOrigin}
          C.legion_domain_point_coloring_color_domain(srcColoring, c, srcRect)
        end
      end
      var p_Fluid0_src = partition(disjoint, SIM0.Fluid, srcColoring, SIM1.tiles)
      C.legion_domain_point_coloring_destroy(srcColoring)
      var tgtColoring = C.legion_domain_point_coloring_create()
      C.legion_domain_point_coloring_color_domain(tgtColoring, int1d(0),
        rect3d{lo = int3d{mc.copyTgt.fromCell[0], mc.copyTgt.fromCell[1], mc.copyTgt.fromCell[2]},
               hi = int3d{mc.copyTgt.uptoCell[0], mc.copyTgt.uptoCell[1], mc.copyTgt.uptoCell[2]}})
      var p_Fluid1_isCopied = partition(disjoint, SIM1.Fluid, tgtColoring, ispace(int1d,1))
      C.legion_domain_point_coloring_destroy(tgtColoring)
      var p_Fluid1_tgt = cross_product(SIM1.p_Fluid, p_Fluid1_isCopied)
      -- Main simulation loop
      while true do
        var Integrator_timeStep = SIM0.Integrator.timeStep;
        -- Perform preliminary actions before each timestep
        [parallelizeFor(SIM0, SIM0.MainLoopHeader())];
        [parallelizeFor(SIM1, SIM1.MainLoopHeader())];
        -- Make sure both simulations are using the same timestep
        SIM0.Integrator.deltaTime = min(SIM0.Integrator.deltaTime, SIM1.Integrator.deltaTime)
        SIM1.Integrator.deltaTime = min(SIM0.Integrator.deltaTime, SIM1.Integrator.deltaTime);
        [parallelizeFor(SIM0, SIM0.PerformIO())];
        [parallelizeFor(SIM1, SIM1.PerformIO())];
        if SIM0.Integrator.exitCond or SIM1.Integrator.exitCond then
          break
        end
        -- Run one iteration of first section
        [parallelizeFor(SIM0, SIM0.MainLoopBody(rexpr false end, FakeCopyQueue))];
        -- Copy fluid & particles to second section
        var incoming = Integrator_timeStep % mc.copyEveryTimeSteps == 0
        if incoming then
          if SIM0.DEBUG_COPYING then
            [SIM0.DumpHDF('copysrc%010d', Integrator_timeStep)];
          end
          for c in SIM1.tiles do
            var src = p_Fluid0_src[c]
            var tgt = p_Fluid1_tgt[c][0]
            copy(src.temperature, tgt.temperature_inc)
            copy(src.velocity, tgt.velocity_inc)
          end
          fill(CopyQueue.__valid, false)
          fill(CopyQueue.position, array(-1.0, -1.0, -1.0))
          fill(CopyQueue.velocity, array(-1.0, -1.0, -1.0))
          fill(CopyQueue.temperature, -1.0)
          fill(CopyQueue.diameter, -1.0)
          fill(CopyQueue.density, -1.0)
          if CopyQueue_size > 0 and C.finite(SIM1.Flow_averagePressure) == 1 then
            for c in SIM0.tiles do
              CopyQueue_push(SIM0.p_Particles[c],
                             p_CopyQueue[c],
                             mc.configs[0],
                             mc.copySrc,
                             SIM0.Grid.particleCopyOrigin, SIM1.Grid.particleCopyOrigin)
            end
          end
        end
        -- Run one iteration of second section
        [parallelizeFor(SIM1, SIM1.MainLoopBody(incoming, CopyQueue))];
      end
      -- Cleanups
      [SIM0.Cleanup()];
      [SIM1.Cleanup()];
    @TIME
  @EPACSE
end

__demand(__inline)
task initSingle(config : &Config, launched : int, outDirBase : &int8)
  config.Mapping.sampleId = launched
  C.snprintf([&int8](config.Mapping.outDir), 256, "%s/sample%d", outDirBase, launched)
  UTIL.createDir(config.Mapping.outDir)
end

__demand(__inline)
task initDual(mc : &MultiConfig, launched : int, outDirBase : &int8)
  initSingle([&Config](mc.configs), launched, outDirBase)
  initSingle([&Config](mc.configs) + 1, launched + 1, outDirBase)
  -- Check 2-section configuration
  regentlib.assert(
    mc.configs[0].Particles.staggerFactor <= mc.copyEveryTimeSteps and
    mc.configs[1].Particles.staggerFactor <= mc.copyEveryTimeSteps and
    mc.copyEveryTimeSteps % mc.configs[0].Particles.staggerFactor == 0 and
    mc.copyEveryTimeSteps % mc.configs[1].Particles.staggerFactor == 0,
    'Invalid stagger factor configuration')
end

__forbid(__optimize) __demand(__inner)
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
      var config : Config
      SCHEMA.parse_Config(&config, args.argv[i+1])
      initSingle(&config, launched, outDirBase)
      launched += 1
      workSingle(config)
    elseif C.strcmp(args.argv[i], '-m') == 0 and i < args.argc-1 then
      var mc : MultiConfig
      SCHEMA.parse_MultiConfig(&mc, args.argv[i+1])
      initDual(&mc, launched, outDirBase)
      launched += 2
      workDual(mc)
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
