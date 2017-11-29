import "regent"

local NUM_ANGLES=14

local JSON = terralib.includec("json.h")
local C = terralib.includecstring[[
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
]]

local c = regentlib.c

local json_dir = os.getenv("JSON_DIR")
assert(json_dir)

local root_dir = arg[0]:match(".*/") or "./"
local cmapper
do
  assert(os.getenv('LG_RT_DIR'))
  local runtime_dir = os.getenv('LG_RT_DIR') .. "/"
  local legion_dir = runtime_dir .. "legion/"
  local mapper_dir = runtime_dir .. "mappers/"
  local realm_dir = runtime_dir .. "realm/"
  local mapper_cc = root_dir .. "ensemble_mapper.cc"
  local mapper_so = root_dir .. "libensemble_mapper.so"
  local cxx = os.getenv('CXX') or 'c++'

  local cxx_flags = os.getenv('CC_FLAGS') or ''
  cxx_flags = cxx_flags .. " -O2 -Wall -Werror"
  if os.execute('test "$(uname)" = Darwin') == 0 then
    cxx_flags =
      (cxx_flags ..
         " -dynamiclib -single_module -undefined dynamic_lookup -fPIC")
  else
    cxx_flags = cxx_flags .. " -shared -fPIC"
  end

  local cmd = (cxx .. " " .. cxx_flags .. " -I " .. runtime_dir .. " " ..
                 " -I " .. mapper_dir .. " " .. " -I " .. legion_dir .. " " ..
                 " -I " .. realm_dir .. " " ..
                 " -I " .. json_dir .. "/include/json-parser " ..
                 mapper_cc .. " -o " .. mapper_so)
  if os.execute(cmd) ~= 0 then
    print("Error: failed to compile " .. mapper_cc)
    assert(false)
  end
  if os.getenv("SAVE_MAPPER_ONLY") == "1" then os.exit() end
  terralib.linklibrary(mapper_so)
  cmapper = terralib.includec("ensemble_mapper.h", {"-I", root_dir, "-I", runtime_dir,
                                                   "-I", mapper_dir, "-I", legion_dir,
                                                   "-I", realm_dir})
end
local SOLEIL_TYPES = terralib.includec("soleil_types.h", { "-I", root_dir })

struct particles_columns {
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

struct Fluid_columns {
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

local IO = SOLEIL_TYPES.IO
local Integrator = SOLEIL_TYPES.Integrator
local Flow = SOLEIL_TYPES.Flow
local Particles = SOLEIL_TYPES.Particles
local BC = SOLEIL_TYPES.BC
local Grid = SOLEIL_TYPES.Grid
local Radiation = SOLEIL_TYPES.Radiation
local Config = SOLEIL_TYPES.Config

terra Grid:grid_volume()
  return self.xNum * self.yNum * self.zNum
end

__demand(__parallel)
task Particles_InitializeDensity(v5574 : region(ispace(int1d), particles_columns), v5576 : double)
-- leaf (false), inner (false), idempotent (false)
where
  reads(v5574.density), writes(v5574.density), reads(v5574.__valid)
do
  for v5578 in v5574 do
    if v5574[v5578].__valid then
      v5574[v5578].density = v5576
    else
    end
  end
end

__demand(__parallel)
task Particles_CalculateNumber(v5584 : region(ispace(int1d), particles_columns)) : int64
-- leaf (false), inner (false), idempotent (false)
where
  reads(v5584.__valid)
do
  var v5589 = int64(0)
  for v5590 in v5584 do
    if v5584[v5590].__valid then
      v5589 += int64(int32(1))
    else
    end
  end
  return v5589
end

__demand(__parallel)
task InitParticlesUniform(v5 : region(ispace(int1d), particles_columns), v7 : region(ispace(int3d), Fluid_columns), v9 : Config, v10 : int32, v11 : int32, v12 : int32)
where
  reads(v5), writes(v5), reads(v7.velocity), reads(v7.centerCoordinates)
do
  var v13 = 0
  for v14 in v5 do
    v13 = int32(v14)
    break
  end
  var v15 = v7.bounds.lo
  v15.x = max(v15.x, v10)
  v15.y = max(v15.y, v11)
  v15.z = max(v15.z, v12)
  var v16 = v7.bounds.hi
  v16.x = min(v16.x, ((v9.grid.xNum+v10)-1))
  v16.y = min(v16.y, ((v9.grid.yNum+v11)-1))
  v16.z = min(v16.z, ((v9.grid.zNum+v12)-1))
  var v17 = ((v16.x-v15.x)+1)
  var v18 = ((v16.y-v15.y)+1)
  var v19 = (v9.particles.initNum/((v9.grid.xTiles*v9.grid.yTiles)*v9.grid.zTiles))
  for v20 in v5 do
    if ((int32(v20)-v13)<v19) then
      v20.__valid = true
      var v21 = (int32(v20)-v13)
      var v22 = int3d({(v15.x+(v21%v17)), (v15.y+((v21/v17)%v18)), (v15.z+((v21/v17)/v18))})
      v20.cell = v22
      v20.position = v7[v20.cell].centerCoordinates
      v20.velocity = v7[v20.cell].velocity
      v20.density = v9.particles.density
      v20.temperature = v9.particles.initTemperature
      v20.diameter = v9.particles.diameterMean
    else
    end
  end
end

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
         terra read_val(f : &C._IO_FILE,val : &double) : int32
             return C.fscanf(f, "%lf\n", &val[0])
         end
struct angle {
  xi : double;
  eta : double;
  mu : double;
  w : double;
}
struct face {
  I_1 : double[NUM_ANGLES];
  I_2 : double[NUM_ANGLES];
  I_3 : double[NUM_ANGLES];
  I_4 : double[NUM_ANGLES];
  I_5 : double[NUM_ANGLES];
  I_6 : double[NUM_ANGLES];
  I_7 : double[NUM_ANGLES];
  I_8 : double[NUM_ANGLES];
}
task initialize_angles(vangles : region(ispace(int1d), angle))
-- leaf (false), inner (false), idempotent (false)
where
  reads(vangles.xi), reads(vangles.eta), reads(vangles.mu), reads(vangles.w), writes(vangles.xi), writes(vangles.eta), writes(vangles.mu), writes(vangles.w)
do
  var vval : double[1]
  var vf = C.fopen(["LMquads/"..tostring(NUM_ANGLES)..".txt"], "rb")
  read_val(vf, [&double](vval))
  for va in vangles do
    read_val(vf, [&double](vval))
    va.xi = vval[0]
  end
  for va in vangles do
    read_val(vf, [&double](vval))
    va.eta = vval[0]
  end
  for va in vangles do
    read_val(vf, [&double](vval))
    va.mu = vval[0]
  end
  for va in vangles do
    read_val(vf, [&double](vval))
    va.w = vval[0]
  end
  C.fclose(vf)
end

task make_interior_partition_x_hi(vfaces : region(ispace(int3d), face), vx_tiles : ispace(int3d), vNx : int32, vNy : int32, vNz : int32, vntx : int32, vnty : int32, vntz : int32)
-- leaf (false), inner (false), idempotent (false)
  var vcoloring = regentlib.c.legion_domain_point_coloring_create()
  for vtile in vx_tiles do
    var vval = (-1)
    if (vtile.x==(vx_tiles.bounds.hi.x-1)) then
      vval = 0
    else
    end
    var vlo = int3d({x = ((vtile.x*vNx)/vntx), y = ((vtile.y*vNy)/vnty), z = ((vtile.z*vNz)/vntz)})
    var vhi = int3d({x = ((((vtile.x+1)*vNx)/vntx)+vval), y = ((((vtile.y+1)*vNy)/vnty)-1), z = ((((vtile.z+1)*vNz)/vntz)-1)})
    if (vhi.x>=(vNx+1)) then
      vlo.x = int64(1)
      vhi.x = int64(0)
    else
    end
    var vrect = rect3d({lo = vlo, hi = vhi})
    regentlib.c.legion_domain_point_coloring_color_domain(vcoloring, regentlib.c.legion_domain_point_t(vtile), regentlib.c.legion_domain_t(vrect))
  end
  var vp = partition(disjoint, vfaces, vcoloring, vx_tiles)
  regentlib.c.legion_domain_point_coloring_destroy(vcoloring)
  return vp
end

task make_interior_partition_x_lo(vfaces : region(ispace(int3d), face), vx_tiles : ispace(int3d), vNx : int32, vNy : int32, vNz : int32, vntx : int32, vnty : int32, vntz : int32)
-- leaf (false), inner (false), idempotent (false)
  var vcoloring = regentlib.c.legion_domain_point_coloring_create()
  for vtile in vx_tiles do
    var vval = 1
    if (vtile.x==1) then
      vval = 0
    else
    end
    var vlo = int3d({x = ((((vtile.x-1)*vNx)/vntx)+vval), y = ((vtile.y*vNy)/vnty), z = ((vtile.z*vNz)/vntz)})
    var vhi = int3d({x = ((vtile.x*vNx)/vntx), y = ((((vtile.y+1)*vNy)/vnty)-1), z = ((((vtile.z+1)*vNz)/vntz)-1)})
    if (vlo.x<0) then
      vlo.x = int64(1)
      vhi.x = int64(0)
    else
    end
    var vrect = rect3d({lo = vlo, hi = vhi})
    regentlib.c.legion_domain_point_coloring_color_domain(vcoloring, regentlib.c.legion_domain_point_t(vtile), regentlib.c.legion_domain_t(vrect))
  end
  var vp = partition(disjoint, vfaces, vcoloring, vx_tiles)
  regentlib.c.legion_domain_point_coloring_destroy(vcoloring)
  return vp
end

task make_interior_partition_y_hi(vfaces : region(ispace(int3d), face), vy_tiles : ispace(int3d), vNx : int32, vNy : int32, vNz : int32, vntx : int32, vnty : int32, vntz : int32)
-- leaf (false), inner (false), idempotent (false)
  var vcoloring = regentlib.c.legion_domain_point_coloring_create()
  for vtile in vy_tiles do
    var vval = (-1)
    if (vtile.y==(vy_tiles.bounds.hi.y-1)) then
      vval = 0
    else
    end
    var vlo = int3d({x = ((vtile.x*vNx)/vntx), y = ((vtile.y*vNy)/vnty), z = ((vtile.z*vNz)/vntz)})
    var vhi = int3d({x = ((((vtile.x+1)*vNx)/vntx)-1), y = ((((vtile.y+1)*vNy)/vnty)+vval), z = ((((vtile.z+1)*vNz)/vntz)-1)})
    if (vhi.y>=(vNy+1)) then
      vlo.y = int64(1)
      vhi.y = int64(0)
    else
    end
    var vrect = rect3d({lo = vlo, hi = vhi})
    regentlib.c.legion_domain_point_coloring_color_domain(vcoloring, regentlib.c.legion_domain_point_t(vtile), regentlib.c.legion_domain_t(vrect))
  end
  var vp = partition(disjoint, vfaces, vcoloring, vy_tiles)
  regentlib.c.legion_domain_point_coloring_destroy(vcoloring)
  return vp
end

task make_interior_partition_y_lo(vfaces : region(ispace(int3d), face), vy_tiles : ispace(int3d), vNx : int32, vNy : int32, vNz : int32, vntx : int32, vnty : int32, vntz : int32)
-- leaf (false), inner (false), idempotent (false)
  var vcoloring = regentlib.c.legion_domain_point_coloring_create()
  for vtile in vy_tiles do
    var vval = 1
    if (vtile.y==1) then
      vval = 0
    else
    end
    var vlo = int3d({x = ((vtile.x*vNx)/vntx), y = ((((vtile.y-1)*vNy)/vnty)+vval), z = ((vtile.z*vNz)/vntz)})
    var vhi = int3d({x = ((((vtile.x+1)*vNx)/vntx)-1), y = ((vtile.y*vNy)/vnty), z = ((((vtile.z+1)*vNz)/vntz)-1)})
    if (vlo.y<0) then
      vlo.y = int64(1)
      vhi.y = int64(0)
    else
    end
    var vrect = rect3d({lo = vlo, hi = vhi})
    regentlib.c.legion_domain_point_coloring_color_domain(vcoloring, regentlib.c.legion_domain_point_t(vtile), regentlib.c.legion_domain_t(vrect))
  end
  var vp = partition(disjoint, vfaces, vcoloring, vy_tiles)
  regentlib.c.legion_domain_point_coloring_destroy(vcoloring)
  return vp
end

task make_interior_partition_z_hi(vfaces : region(ispace(int3d), face), vz_tiles : ispace(int3d), vNx : int32, vNy : int32, vNz : int32, vntx : int32, vnty : int32, vntz : int32)
-- leaf (false), inner (false), idempotent (false)
  var vcoloring = regentlib.c.legion_domain_point_coloring_create()
  for vtile in vz_tiles do
    var vval = (-1)
    if (vtile.z==(vz_tiles.bounds.hi.z-1)) then
      vval = 0
    else
    end
    var vlo = int3d({x = ((vtile.x*vNx)/vntx), y = ((vtile.y*vNy)/vnty), z = ((vtile.z*vNz)/vntz)})
    var vhi = int3d({x = ((((vtile.x+1)*vNx)/vntx)-1), y = ((((vtile.y+1)*vNy)/vnty)-1), z = ((((vtile.z+1)*vNz)/vntz)+vval)})
    if (vhi.z>=(vNz+1)) then
      vlo.z = int64(1)
      vhi.z = int64(0)
    else
    end
    var vrect = rect3d({lo = vlo, hi = vhi})
    regentlib.c.legion_domain_point_coloring_color_domain(vcoloring, regentlib.c.legion_domain_point_t(vtile), regentlib.c.legion_domain_t(vrect))
  end
  var vp = partition(disjoint, vfaces, vcoloring, vz_tiles)
  regentlib.c.legion_domain_point_coloring_destroy(vcoloring)
  return vp
end

task make_interior_partition_z_lo(vfaces : region(ispace(int3d), face), vz_tiles : ispace(int3d), vNx : int32, vNy : int32, vNz : int32, vntx : int32, vnty : int32, vntz : int32)
-- leaf (false), inner (false), idempotent (false)
  var vcoloring = regentlib.c.legion_domain_point_coloring_create()
  for vtile in vz_tiles do
    var vval = 1
    if (vtile.z==1) then
      vval = 0
    else
    end
    var vlo = int3d({x = ((vtile.x*vNx)/vntx), y = ((vtile.y*vNy)/vnty), z = ((((vtile.z-1)*vNz)/vntz)+vval)})
    var vhi = int3d({x = ((((vtile.x+1)*vNx)/vntx)-1), y = ((((vtile.y+1)*vNy)/vnty)-1), z = ((vtile.z*vNz)/vntz)})
    if (vlo.z<0) then
      vlo.z = int64(1)
      vhi.z = int64(0)
    else
    end
    var vrect = rect3d({lo = vlo, hi = vhi})
    regentlib.c.legion_domain_point_coloring_color_domain(vcoloring, regentlib.c.legion_domain_point_t(vtile), regentlib.c.legion_domain_t(vrect))
  end
  var vp = partition(disjoint, vfaces, vcoloring, vz_tiles)
  regentlib.c.legion_domain_point_coloring_destroy(vcoloring)
  return vp
end

task source_term(vpoints : region(ispace(int3d), Radiation_columns), vangles : region(ispace(int1d), angle), vomega : double)
-- leaf (false), inner (false), idempotent (false)
where
  reads(vpoints.Iiter_1), reads(vpoints.Iiter_2), reads(vpoints.Iiter_3), reads(vpoints.Iiter_4), reads(vpoints.Iiter_5), reads(vpoints.Iiter_6), reads(vpoints.Iiter_7), reads(vpoints.Iiter_8), reads(vpoints.Ib), reads(vpoints.sigma), reads(vangles.w), reads(vpoints.S), writes(vpoints.S)
do
  for vp in vpoints do
    vp.S = (((1-vomega)*vp.sigma)*vp.Ib)
    for vm : int32 = 0, NUM_ANGLES do
      vp.S += (((((((((((vomega*vp.sigma)/(4*3.1415926535898))*vangles[vm].w)*vp.Iiter_1[vm])+((((vomega*vp.sigma)/(4*3.1415926535898))*vangles[vm].w)*vp.Iiter_2[vm]))+((((vomega*vp.sigma)/(4*3.1415926535898))*vangles[vm].w)*vp.Iiter_3[vm]))+((((vomega*vp.sigma)/(4*3.1415926535898))*vangles[vm].w)*vp.Iiter_4[vm]))+((((vomega*vp.sigma)/(4*3.1415926535898))*vangles[vm].w)*vp.Iiter_5[vm]))+((((vomega*vp.sigma)/(4*3.1415926535898))*vangles[vm].w)*vp.Iiter_6[vm]))+((((vomega*vp.sigma)/(4*3.1415926535898))*vangles[vm].w)*vp.Iiter_7[vm]))+((((vomega*vp.sigma)/(4*3.1415926535898))*vangles[vm].w)*vp.Iiter_8[vm]))
    end
  end
end

task west_bound(vfaces_1 : region(ispace(int3d), face), vfaces_2 : region(ispace(int3d), face), vfaces_3 : region(ispace(int3d), face), vfaces_4 : region(ispace(int3d), face), vfaces_5 : region(ispace(int3d), face), vfaces_6 : region(ispace(int3d), face), vfaces_7 : region(ispace(int3d), face), vfaces_8 : region(ispace(int3d), face), vangles : region(ispace(int1d), angle), vemissWest : double, vtempWest : double)
-- leaf (false), inner (false), idempotent (false)
where
  reads(vangles.w), reads(vangles.xi), reads(vangles.eta), reads(vangles.mu), reads(vfaces_1.I_1), reads(vfaces_2.I_2), reads(vfaces_3.I_3), reads(vfaces_4.I_4), reads(vfaces_5.I_5), reads(vfaces_6.I_6), reads(vfaces_7.I_7), reads(vfaces_8.I_8), writes(vfaces_1.I_1), writes(vfaces_2.I_2), writes(vfaces_3.I_3), writes(vfaces_4.I_4), writes(vfaces_5.I_5), writes(vfaces_6.I_6), writes(vfaces_7.I_7), writes(vfaces_8.I_8)
do
  var vlimits = vfaces_1.bounds
  var vreflect = 0
  var vepsw = vemissWest
  var vTw = vtempWest
  for vj : int64 = vlimits.lo.y, (vlimits.hi.y+1) do
    for vk : int64 = vlimits.lo.z, (vlimits.hi.z+1) do
      vreflect = double(0)
      for vm : int32 = 0, NUM_ANGLES do
        if (vangles[vm].xi<0) then
          var vface_value = 0
          if ((vangles[vm].eta>0) and (vangles[vm].mu>0)) then
            vface_value = vfaces_5[{vlimits.lo.x, vj, vk}].I_5[vm]
          elseif ((vangles[vm].eta>0) and (vangles[vm].mu<0)) then
            vface_value = vfaces_6[{vlimits.lo.x, vj, vk}].I_6[vm]
          elseif ((vangles[vm].eta<0) and (vangles[vm].mu>0)) then
            vface_value = vfaces_7[{vlimits.lo.x, vj, vk}].I_7[vm]
          else
            vface_value = vfaces_8[{vlimits.lo.x, vj, vk}].I_8[vm]
          end
          vreflect += (((((1-vepsw)/3.1415926535898)*vangles[vm].w)*C.fabs(vangles[vm].xi))*vface_value)
        else
        end
      end
      for vm : int32 = 0, NUM_ANGLES do
        if (vangles[vm].xi>0) then
          var vvalue = ((((vepsw*5.67e-08)*C.pow(vTw, 4))/3.1415926535898)+vreflect)
          if ((vangles[vm].eta>0) and (vangles[vm].mu>0)) then
            vfaces_1[{vlimits.lo.x, vj, vk}].I_1[vm] = vvalue
          elseif ((vangles[vm].eta>0) and (vangles[vm].mu<0)) then
            vfaces_2[{vlimits.lo.x, vj, vk}].I_2[vm] = vvalue
          elseif ((vangles[vm].eta<0) and (vangles[vm].mu>0)) then
            vfaces_3[{vlimits.lo.x, vj, vk}].I_3[vm] = vvalue
          else
            vfaces_4[{vlimits.lo.x, vj, vk}].I_4[vm] = vvalue
          end
        else
        end
      end
    end
  end
end

task east_bound(vfaces_1 : region(ispace(int3d), face), vfaces_2 : region(ispace(int3d), face), vfaces_3 : region(ispace(int3d), face), vfaces_4 : region(ispace(int3d), face), vfaces_5 : region(ispace(int3d), face), vfaces_6 : region(ispace(int3d), face), vfaces_7 : region(ispace(int3d), face), vfaces_8 : region(ispace(int3d), face), vangles : region(ispace(int1d), angle), vemissEast : double, vtempEast : double)
-- leaf (false), inner (false), idempotent (false)
where
  reads(vangles.w), reads(vangles.xi), reads(vangles.eta), reads(vangles.mu), reads(vfaces_1.I_1), reads(vfaces_2.I_2), reads(vfaces_3.I_3), reads(vfaces_4.I_4), reads(vfaces_5.I_5), reads(vfaces_6.I_6), reads(vfaces_7.I_7), reads(vfaces_8.I_8), writes(vfaces_1.I_1), writes(vfaces_2.I_2), writes(vfaces_3.I_3), writes(vfaces_4.I_4), writes(vfaces_5.I_5), writes(vfaces_6.I_6), writes(vfaces_7.I_7), writes(vfaces_8.I_8)
do
  var vlimits = vfaces_1.bounds
  var vreflect = 0
  var vepsw = vemissEast
  var vTw = vtempEast
  for vj : int64 = vlimits.lo.y, (vlimits.hi.y+1) do
    for vk : int64 = vlimits.lo.z, (vlimits.hi.z+1) do
      vreflect = double(0)
      for vm : int32 = 0, NUM_ANGLES do
        if (vangles[vm].xi>0) then
          var vface_value = 0
          if ((vangles[vm].eta>0) and (vangles[vm].mu>0)) then
            vface_value = vfaces_1[{vlimits.hi.x, vj, vk}].I_1[vm]
          elseif ((vangles[vm].eta>0) and (vangles[vm].mu<0)) then
            vface_value = vfaces_2[{vlimits.hi.x, vj, vk}].I_2[vm]
          elseif ((vangles[vm].eta<0) and (vangles[vm].mu>0)) then
            vface_value = vfaces_3[{vlimits.hi.x, vj, vk}].I_3[vm]
          else
            vface_value = vfaces_4[{vlimits.hi.x, vj, vk}].I_4[vm]
          end
          vreflect += (((((1-vepsw)/3.1415926535898)*vangles[vm].w)*vangles[vm].xi)*vface_value)
        else
        end
      end
      for vm : int32 = 0, NUM_ANGLES do
        if (vangles[vm].xi<0) then
          var vvalue = ((((vepsw*5.67e-08)*C.pow(vTw, 4))/3.1415926535898)+vreflect)
          if ((vangles[vm].eta>0) and (vangles[vm].mu>0)) then
            vfaces_5[{vlimits.hi.x, vj, vk}].I_5[vm] = vvalue
          elseif ((vangles[vm].eta>0) and (vangles[vm].mu<0)) then
            vfaces_6[{vlimits.hi.x, vj, vk}].I_6[vm] = vvalue
          elseif ((vangles[vm].eta<0) and (vangles[vm].mu>0)) then
            vfaces_7[{vlimits.hi.x, vj, vk}].I_7[vm] = vvalue
          else
            vfaces_8[{vlimits.hi.x, vj, vk}].I_8[vm] = vvalue
          end
        else
        end
      end
    end
  end
end

task north_bound(vfaces_1 : region(ispace(int3d), face), vfaces_2 : region(ispace(int3d), face), vfaces_3 : region(ispace(int3d), face), vfaces_4 : region(ispace(int3d), face), vfaces_5 : region(ispace(int3d), face), vfaces_6 : region(ispace(int3d), face), vfaces_7 : region(ispace(int3d), face), vfaces_8 : region(ispace(int3d), face), vangles : region(ispace(int1d), angle), vemissNorth : double, vtempNorth : double)
-- leaf (false), inner (false), idempotent (false)
where
  reads(vangles.w), reads(vangles.xi), reads(vangles.eta), reads(vangles.mu), reads(vfaces_1.I_1), reads(vfaces_2.I_2), reads(vfaces_3.I_3), reads(vfaces_4.I_4), reads(vfaces_5.I_5), reads(vfaces_6.I_6), reads(vfaces_7.I_7), reads(vfaces_8.I_8), writes(vfaces_1.I_1), writes(vfaces_2.I_2), writes(vfaces_3.I_3), writes(vfaces_4.I_4), writes(vfaces_5.I_5), writes(vfaces_6.I_6), writes(vfaces_7.I_7), writes(vfaces_8.I_8)
do
  var vlimits = vfaces_1.bounds
  var vreflect = 0
  var vepsw = vemissNorth
  var vTw = vtempNorth
  for vi : int64 = vlimits.lo.x, (vlimits.hi.x+1) do
    for vk : int64 = vlimits.lo.z, (vlimits.hi.z+1) do
      vreflect = double(0)
      for vm : int32 = 0, NUM_ANGLES do
        if (vangles[vm].eta>0) then
          var vface_value = 0
          if ((vangles[vm].xi>0) and (vangles[vm].mu>0)) then
            vface_value = vfaces_1[{vi, vlimits.hi.y, vk}].I_1[vm]
          elseif ((vangles[vm].xi>0) and (vangles[vm].mu<0)) then
            vface_value = vfaces_2[{vi, vlimits.hi.y, vk}].I_2[vm]
          elseif ((vangles[vm].xi<0) and (vangles[vm].mu>0)) then
            vface_value = vfaces_5[{vi, vlimits.hi.y, vk}].I_5[vm]
          else
            vface_value = vfaces_6[{vi, vlimits.hi.y, vk}].I_6[vm]
          end
          vreflect += (((((1-vepsw)/3.1415926535898)*vangles[vm].w)*vangles[vm].eta)*vface_value)
        else
        end
      end
      for vm : int32 = 0, NUM_ANGLES do
        if (vangles[vm].eta<0) then
          var vvalue = ((((vepsw*5.67e-08)*C.pow(vTw, 4))/3.1415926535898)+vreflect)
          if ((vangles[vm].xi>0) and (vangles[vm].mu>0)) then
            vfaces_3[{vi, vlimits.hi.y, vk}].I_3[vm] = vvalue
          elseif ((vangles[vm].xi>0) and (vangles[vm].mu<0)) then
            vfaces_4[{vi, vlimits.hi.y, vk}].I_4[vm] = vvalue
          elseif ((vangles[vm].xi<0) and (vangles[vm].mu>0)) then
            vfaces_7[{vi, vlimits.hi.y, vk}].I_7[vm] = vvalue
          else
            vfaces_8[{vi, vlimits.hi.y, vk}].I_8[vm] = vvalue
          end
        else
        end
      end
    end
  end
end

task south_bound(vfaces_1 : region(ispace(int3d), face), vfaces_2 : region(ispace(int3d), face), vfaces_3 : region(ispace(int3d), face), vfaces_4 : region(ispace(int3d), face), vfaces_5 : region(ispace(int3d), face), vfaces_6 : region(ispace(int3d), face), vfaces_7 : region(ispace(int3d), face), vfaces_8 : region(ispace(int3d), face), vangles : region(ispace(int1d), angle), vemissSouth : double, vtempSouth : double)
-- leaf (false), inner (false), idempotent (false)
where
  reads(vangles.w), reads(vangles.xi), reads(vangles.eta), reads(vangles.mu), reads(vfaces_1.I_1), reads(vfaces_2.I_2), reads(vfaces_3.I_3), reads(vfaces_4.I_4), reads(vfaces_5.I_5), reads(vfaces_6.I_6), reads(vfaces_7.I_7), reads(vfaces_8.I_8), writes(vfaces_1.I_1), writes(vfaces_2.I_2), writes(vfaces_3.I_3), writes(vfaces_4.I_4), writes(vfaces_5.I_5), writes(vfaces_6.I_6), writes(vfaces_7.I_7), writes(vfaces_8.I_8)
do
  var vlimits = vfaces_1.bounds
  var vreflect = 0
  var vepsw = vemissSouth
  var vTw = vtempSouth
  for vi : int64 = vlimits.lo.x, (vlimits.hi.x+1) do
    for vk : int64 = vlimits.lo.z, (vlimits.hi.z+1) do
      vreflect = double(0)
      for vm : int32 = 0, NUM_ANGLES do
        if (vangles[vm].eta<0) then
          var vface_value = 0
          if ((vangles[vm].xi>0) and (vangles[vm].mu>0)) then
            vface_value = vfaces_3[{vi, vlimits.lo.y, vk}].I_3[vm]
          elseif ((vangles[vm].xi>0) and (vangles[vm].mu<0)) then
            vface_value = vfaces_4[{vi, vlimits.lo.y, vk}].I_4[vm]
          elseif ((vangles[vm].xi<0) and (vangles[vm].mu>0)) then
            vface_value = vfaces_7[{vi, vlimits.lo.y, vk}].I_7[vm]
          else
            vface_value = vfaces_8[{vi, vlimits.lo.y, vk}].I_8[vm]
          end
          vreflect += (((((1-vepsw)/3.1415926535898)*vangles[vm].w)*C.fabs(vangles[vm].eta))*vface_value)
        else
        end
      end
      for vm : int32 = 0, NUM_ANGLES do
        if (vangles[vm].eta>0) then
          var vvalue = ((((vepsw*5.67e-08)*C.pow(vTw, 4))/3.1415926535898)+vreflect)
          if ((vangles[vm].xi>0) and (vangles[vm].mu>0)) then
            vfaces_1[{vi, vlimits.lo.y, vk}].I_1[vm] = vvalue
          elseif ((vangles[vm].xi>0) and (vangles[vm].mu<0)) then
            vfaces_2[{vi, vlimits.lo.y, vk}].I_2[vm] = vvalue
          elseif ((vangles[vm].xi<0) and (vangles[vm].mu>0)) then
            vfaces_5[{vi, vlimits.lo.y, vk}].I_5[vm] = vvalue
          else
            vfaces_6[{vi, vlimits.lo.y, vk}].I_6[vm] = vvalue
          end
        else
        end
      end
    end
  end
end

task up_bound(vfaces_1 : region(ispace(int3d), face), vfaces_2 : region(ispace(int3d), face), vfaces_3 : region(ispace(int3d), face), vfaces_4 : region(ispace(int3d), face), vfaces_5 : region(ispace(int3d), face), vfaces_6 : region(ispace(int3d), face), vfaces_7 : region(ispace(int3d), face), vfaces_8 : region(ispace(int3d), face), vangles : region(ispace(int1d), angle), vemissUp : double, vtempUp : double)
-- leaf (false), inner (false), idempotent (false)
where
  reads(vangles.w), reads(vangles.xi), reads(vangles.eta), reads(vangles.mu), reads(vfaces_1.I_1), reads(vfaces_2.I_2), reads(vfaces_3.I_3), reads(vfaces_4.I_4), reads(vfaces_5.I_5), reads(vfaces_6.I_6), reads(vfaces_7.I_7), reads(vfaces_8.I_8), writes(vfaces_1.I_1), writes(vfaces_2.I_2), writes(vfaces_3.I_3), writes(vfaces_4.I_4), writes(vfaces_5.I_5), writes(vfaces_6.I_6), writes(vfaces_7.I_7), writes(vfaces_8.I_8)
do
  var vlimits = vfaces_1.bounds
  var vreflect = 0
  var vepsw = vemissUp
  var vTw = vtempUp
  for vi : int64 = vlimits.lo.x, (vlimits.hi.x+1) do
    for vj : int64 = vlimits.lo.y, (vlimits.hi.y+1) do
      vreflect = double(0)
      for vm : int32 = 0, NUM_ANGLES do
        if (vangles[vm].mu<0) then
          var vface_value = 0
          if ((vangles[vm].xi>0) and (vangles[vm].eta>0)) then
            vface_value = vfaces_2[{vi, vj, vlimits.lo.z}].I_2[vm]
          elseif ((vangles[vm].xi>0) and (vangles[vm].eta<0)) then
            vface_value = vfaces_4[{vi, vj, vlimits.lo.z}].I_4[vm]
          elseif ((vangles[vm].xi<0) and (vangles[vm].eta>0)) then
            vface_value = vfaces_6[{vi, vj, vlimits.lo.z}].I_6[vm]
          else
            vface_value = vfaces_8[{vi, vj, vlimits.lo.z}].I_8[vm]
          end
          vreflect += (((((1-vepsw)/3.1415926535898)*vangles[vm].w)*C.fabs(vangles[vm].mu))*vface_value)
        else
        end
      end
      for vm : int32 = 0, NUM_ANGLES do
        if (vangles[vm].mu>0) then
          var vvalue = ((((vepsw*5.67e-08)*C.pow(vTw, 4))/3.1415926535898)+vreflect)
          if ((vangles[vm].xi>0) and (vangles[vm].eta>0)) then
            vfaces_1[{vi, vj, vlimits.lo.z}].I_1[vm] = vvalue
          elseif ((vangles[vm].xi>0) and (vangles[vm].eta<0)) then
            vfaces_3[{vi, vj, vlimits.lo.z}].I_3[vm] = vvalue
          elseif ((vangles[vm].xi<0) and (vangles[vm].eta>0)) then
            vfaces_5[{vi, vj, vlimits.lo.z}].I_5[vm] = vvalue
          else
            vfaces_7[{vi, vj, vlimits.lo.z}].I_7[vm] = vvalue
          end
        else
        end
      end
    end
  end
end

task down_bound(vfaces_1 : region(ispace(int3d), face), vfaces_2 : region(ispace(int3d), face), vfaces_3 : region(ispace(int3d), face), vfaces_4 : region(ispace(int3d), face), vfaces_5 : region(ispace(int3d), face), vfaces_6 : region(ispace(int3d), face), vfaces_7 : region(ispace(int3d), face), vfaces_8 : region(ispace(int3d), face), vangles : region(ispace(int1d), angle), vemissDown : double, vtempDown : double)
-- leaf (false), inner (false), idempotent (false)
where
  reads(vangles.w), reads(vangles.xi), reads(vangles.eta), reads(vangles.mu), reads(vfaces_1.I_1), reads(vfaces_2.I_2), reads(vfaces_3.I_3), reads(vfaces_4.I_4), reads(vfaces_5.I_5), reads(vfaces_6.I_6), reads(vfaces_7.I_7), reads(vfaces_8.I_8), writes(vfaces_1.I_1), writes(vfaces_2.I_2), writes(vfaces_3.I_3), writes(vfaces_4.I_4), writes(vfaces_5.I_5), writes(vfaces_6.I_6), writes(vfaces_7.I_7), writes(vfaces_8.I_8)
do
  var vlimits = vfaces_1.bounds
  var vreflect = 0
  var vepsw = vemissDown
  var vTw = vtempDown
  for vi : int64 = vlimits.lo.x, (vlimits.hi.x+1) do
    for vj : int64 = vlimits.lo.y, (vlimits.hi.y+1) do
      vreflect = double(0)
      for vm : int32 = 0, NUM_ANGLES do
        if (vangles[vm].mu>0) then
          var vface_value = 0
          if ((vangles[vm].xi>0) and (vangles[vm].eta>0)) then
            vface_value = vfaces_1[{vi, vj, vlimits.hi.z}].I_1[vm]
          elseif ((vangles[vm].xi>0) and (vangles[vm].eta<0)) then
            vface_value = vfaces_3[{vi, vj, vlimits.hi.z}].I_3[vm]
          elseif ((vangles[vm].xi<0) and (vangles[vm].eta>0)) then
            vface_value = vfaces_5[{vi, vj, vlimits.hi.z}].I_5[vm]
          else
            vface_value = vfaces_7[{vi, vj, vlimits.hi.z}].I_7[vm]
          end
          vreflect += (((((1-vepsw)/3.1415926535898)*vangles[vm].w)*vangles[vm].mu)*vface_value)
        else
        end
      end
      for vm : int32 = 0, NUM_ANGLES do
        if (vangles[vm].mu<0) then
          var vvalue = ((((vepsw*5.67e-08)*C.pow(vTw, 4))/3.1415926535898)+vreflect)
          if ((vangles[vm].xi>0) and (vangles[vm].eta>0)) then
            vfaces_2[{vi, vj, vlimits.hi.z}].I_2[vm] = vvalue
          elseif ((vangles[vm].xi>0) and (vangles[vm].eta<0)) then
            vfaces_4[{vi, vj, vlimits.hi.z}].I_4[vm] = vvalue
          elseif ((vangles[vm].xi<0) and (vangles[vm].eta>0)) then
            vfaces_6[{vi, vj, vlimits.hi.z}].I_6[vm] = vvalue
          else
            vfaces_8[{vi, vj, vlimits.hi.z}].I_8[vm] = vvalue
          end
        else
        end
      end
    end
  end
end

task sweep_1(vpoints : region(ispace(int3d), Radiation_columns), vx_faces : region(ispace(int3d), face), vy_faces : region(ispace(int3d), face), vz_faces : region(ispace(int3d), face), vghost_x_faces : region(ispace(int3d), face), vghost_y_faces : region(ispace(int3d), face), vghost_z_faces : region(ispace(int3d), face), vangles : region(ispace(int1d), angle), vxi : int64, veta : int64, vmu : int64, vdx : double, vdy : double, vdz : double)
-- leaf (false), inner (false), idempotent (false)
where
  reads(vangles.xi), reads(vangles.eta), reads(vangles.mu), reads(vpoints.S), reads(vpoints.sigma), reads(vghost_x_faces.I_1), reads(vghost_y_faces.I_1), reads(vghost_z_faces.I_1), reads(vpoints.I_1), reads(vx_faces.I_1), reads(vy_faces.I_1), reads(vz_faces.I_1), writes(vpoints.I_1), writes(vx_faces.I_1), writes(vy_faces.I_1), writes(vz_faces.I_1)
do
  var vdAx = (vdy*vdz)
  var vdAy = (vdx*vdz)
  var vdAz = (vdx*vdy)
  var vdV = ((vdx*vdy)*vdz)
  var vlimits = vpoints.bounds
  var vdindx = int64(1)
  var vstartx = vlimits.lo.x
  var vendx = (vlimits.hi.x+1)
  var vdindy = int64(1)
  var vstarty = vlimits.lo.y
  var vendy = (vlimits.hi.y+1)
  var vdindz = int64(1)
  var vstartz = vlimits.lo.z
  var vendz = (vlimits.hi.z+1)
  if (vxi<0) then
    vdindx = int64((-1))
    vstartx = vlimits.hi.x
    vendx = (vlimits.lo.x-1)
  else
  end
  if (veta<0) then
    vdindy = int64((-1))
    vstarty = vlimits.hi.y
    vendy = (vlimits.lo.y-1)
  else
  end
  if (vmu<0) then
    vdindz = int64((-1))
    vstartz = vlimits.hi.z
    vendz = (vlimits.lo.z-1)
  else
  end
  for vm : int32 = 0, NUM_ANGLES do
    if (((((vangles[vm].xi*vxi)>0) or ((vangles[vm].xi==0) and (vxi<0))) and (((vangles[vm].eta*veta)>0) or ((vangles[vm].eta==0) and (veta<0)))) and (((vangles[vm].mu*vmu)>0) or ((vangles[vm].mu==0) and (vmu<0)))) then
      for vk : int64 = vstartz, vendz, vdindz do
        for vj : int64 = vstarty, vendy, vdindy do
          for vi : int64 = vstartx, vendx, vdindx do
            var vindx = (vi-min(vdindx, 0))
            var vindy = (vj-min(vdindy, 0))
            var vindz = (vk-min(vdindz, 0))
            var vghost_x_limits = vghost_x_faces.bounds
            var vupwind_x_value = 0
            if (vindx<vx_faces.bounds.lo.x) then
              vupwind_x_value = vghost_x_faces[{vghost_x_limits.hi.x, vj, vk}].I_1[vm]
            elseif (vindx>vx_faces.bounds.hi.x) then
              vupwind_x_value = vghost_x_faces[{vghost_x_limits.lo.x, vj, vk}].I_1[vm]
            else
              vupwind_x_value = vx_faces[{vindx, vj, vk}].I_1[vm]
            end
            var vghost_y_limits = vghost_y_faces.bounds
            var vupwind_y_value = 0
            if (vindy<vy_faces.bounds.lo.y) then
              vupwind_y_value = vghost_y_faces[{vi, vghost_y_limits.hi.y, vk}].I_1[vm]
            elseif (vindy>vy_faces.bounds.hi.y) then
              vupwind_y_value = vghost_y_faces[{vi, vghost_y_limits.lo.y, vk}].I_1[vm]
            else
              vupwind_y_value = vy_faces[{vi, vindy, vk}].I_1[vm]
            end
            var vghost_z_limits = vghost_z_faces.bounds
            var vupwind_z_value = 0
            if (vindz<vz_faces.bounds.lo.z) then
              vupwind_z_value = vghost_z_faces[{vi, vj, vghost_z_limits.hi.z}].I_1[vm]
            elseif (vindz>vz_faces.bounds.hi.z) then
              vupwind_z_value = vghost_z_faces[{vi, vj, vghost_z_limits.lo.z}].I_1[vm]
            else
              vupwind_z_value = vz_faces[{vi, vj, vindz}].I_1[vm]
            end
            vpoints[{vi, vj, vk}].I_1[vm] = (((((vpoints[{vi, vj, vk}].S*vdV)+(((C.fabs(vangles[vm].xi)*vdAx)*vupwind_x_value)/0.5))+(((C.fabs(vangles[vm].eta)*vdAy)*vupwind_y_value)/0.5))+(((C.fabs(vangles[vm].mu)*vdAz)*vupwind_z_value)/0.5))/((((vpoints[{vi, vj, vk}].sigma*vdV)+((C.fabs(vangles[vm].xi)*vdAx)/0.5))+((C.fabs(vangles[vm].eta)*vdAy)/0.5))+((C.fabs(vangles[vm].mu)*vdAz)/0.5)))
            vx_faces[{(vindx+vdindx), vj, vk}].I_1[vm] = ((vpoints[{vi, vj, vk}].I_1[vm]-((1-0.5)*vupwind_x_value))/0.5)
            vy_faces[{vi, (vindy+vdindy), vk}].I_1[vm] = ((vpoints[{vi, vj, vk}].I_1[vm]-((1-0.5)*vupwind_y_value))/0.5)
            vz_faces[{vi, vj, (vindz+vdindz)}].I_1[vm] = ((vpoints[{vi, vj, vk}].I_1[vm]-((1-0.5)*vupwind_z_value))/0.5)
          end
        end
      end
    else
    end
  end
end

task sweep_2(vpoints : region(ispace(int3d), Radiation_columns), vx_faces : region(ispace(int3d), face), vy_faces : region(ispace(int3d), face), vz_faces : region(ispace(int3d), face), vghost_x_faces : region(ispace(int3d), face), vghost_y_faces : region(ispace(int3d), face), vghost_z_faces : region(ispace(int3d), face), vangles : region(ispace(int1d), angle), vxi : int64, veta : int64, vmu : int64, vdx : double, vdy : double, vdz : double)
-- leaf (false), inner (false), idempotent (false)
where
  reads(vangles.xi), reads(vangles.eta), reads(vangles.mu), reads(vpoints.S), reads(vpoints.sigma), reads(vghost_x_faces.I_2), reads(vghost_y_faces.I_2), reads(vghost_z_faces.I_2), reads(vpoints.I_2), reads(vx_faces.I_2), reads(vy_faces.I_2), reads(vz_faces.I_2), writes(vpoints.I_2), writes(vx_faces.I_2), writes(vy_faces.I_2), writes(vz_faces.I_2)
do
  var vdAx = (vdy*vdz)
  var vdAy = (vdx*vdz)
  var vdAz = (vdx*vdy)
  var vdV = ((vdx*vdy)*vdz)
  var vlimits = vpoints.bounds
  var vdindx = int64(1)
  var vstartx = vlimits.lo.x
  var vendx = (vlimits.hi.x+1)
  var vdindy = int64(1)
  var vstarty = vlimits.lo.y
  var vendy = (vlimits.hi.y+1)
  var vdindz = int64(1)
  var vstartz = vlimits.lo.z
  var vendz = (vlimits.hi.z+1)
  if (vxi<0) then
    vdindx = int64((-1))
    vstartx = vlimits.hi.x
    vendx = (vlimits.lo.x-1)
  else
  end
  if (veta<0) then
    vdindy = int64((-1))
    vstarty = vlimits.hi.y
    vendy = (vlimits.lo.y-1)
  else
  end
  if (vmu<0) then
    vdindz = int64((-1))
    vstartz = vlimits.hi.z
    vendz = (vlimits.lo.z-1)
  else
  end
  for vm : int32 = 0, NUM_ANGLES do
    if (((((vangles[vm].xi*vxi)>0) or ((vangles[vm].xi==0) and (vxi<0))) and (((vangles[vm].eta*veta)>0) or ((vangles[vm].eta==0) and (veta<0)))) and (((vangles[vm].mu*vmu)>0) or ((vangles[vm].mu==0) and (vmu<0)))) then
      for vk : int64 = vstartz, vendz, vdindz do
        for vj : int64 = vstarty, vendy, vdindy do
          for vi : int64 = vstartx, vendx, vdindx do
            var vindx = (vi-min(vdindx, 0))
            var vindy = (vj-min(vdindy, 0))
            var vindz = (vk-min(vdindz, 0))
            var vghost_x_limits = vghost_x_faces.bounds
            var vupwind_x_value = 0
            if (vindx<vx_faces.bounds.lo.x) then
              vupwind_x_value = vghost_x_faces[{vghost_x_limits.hi.x, vj, vk}].I_2[vm]
            elseif (vindx>vx_faces.bounds.hi.x) then
              vupwind_x_value = vghost_x_faces[{vghost_x_limits.lo.x, vj, vk}].I_2[vm]
            else
              vupwind_x_value = vx_faces[{vindx, vj, vk}].I_2[vm]
            end
            var vghost_y_limits = vghost_y_faces.bounds
            var vupwind_y_value = 0
            if (vindy<vy_faces.bounds.lo.y) then
              vupwind_y_value = vghost_y_faces[{vi, vghost_y_limits.hi.y, vk}].I_2[vm]
            elseif (vindy>vy_faces.bounds.hi.y) then
              vupwind_y_value = vghost_y_faces[{vi, vghost_y_limits.lo.y, vk}].I_2[vm]
            else
              vupwind_y_value = vy_faces[{vi, vindy, vk}].I_2[vm]
            end
            var vghost_z_limits = vghost_z_faces.bounds
            var vupwind_z_value = 0
            if (vindz<vz_faces.bounds.lo.z) then
              vupwind_z_value = vghost_z_faces[{vi, vj, vghost_z_limits.hi.z}].I_2[vm]
            elseif (vindz>vz_faces.bounds.hi.z) then
              vupwind_z_value = vghost_z_faces[{vi, vj, vghost_z_limits.lo.z}].I_2[vm]
            else
              vupwind_z_value = vz_faces[{vi, vj, vindz}].I_2[vm]
            end
            vpoints[{vi, vj, vk}].I_2[vm] = (((((vpoints[{vi, vj, vk}].S*vdV)+(((C.fabs(vangles[vm].xi)*vdAx)*vupwind_x_value)/0.5))+(((C.fabs(vangles[vm].eta)*vdAy)*vupwind_y_value)/0.5))+(((C.fabs(vangles[vm].mu)*vdAz)*vupwind_z_value)/0.5))/((((vpoints[{vi, vj, vk}].sigma*vdV)+((C.fabs(vangles[vm].xi)*vdAx)/0.5))+((C.fabs(vangles[vm].eta)*vdAy)/0.5))+((C.fabs(vangles[vm].mu)*vdAz)/0.5)))
            vx_faces[{(vindx+vdindx), vj, vk}].I_2[vm] = ((vpoints[{vi, vj, vk}].I_2[vm]-((1-0.5)*vupwind_x_value))/0.5)
            vy_faces[{vi, (vindy+vdindy), vk}].I_2[vm] = ((vpoints[{vi, vj, vk}].I_2[vm]-((1-0.5)*vupwind_y_value))/0.5)
            vz_faces[{vi, vj, (vindz+vdindz)}].I_2[vm] = ((vpoints[{vi, vj, vk}].I_2[vm]-((1-0.5)*vupwind_z_value))/0.5)
          end
        end
      end
    else
    end
  end
end

task sweep_3(vpoints : region(ispace(int3d), Radiation_columns), vx_faces : region(ispace(int3d), face), vy_faces : region(ispace(int3d), face), vz_faces : region(ispace(int3d), face), vghost_x_faces : region(ispace(int3d), face), vghost_y_faces : region(ispace(int3d), face), vghost_z_faces : region(ispace(int3d), face), vangles : region(ispace(int1d), angle), vxi : int64, veta : int64, vmu : int64, vdx : double, vdy : double, vdz : double)
-- leaf (false), inner (false), idempotent (false)
where
  reads(vangles.xi), reads(vangles.eta), reads(vangles.mu), reads(vpoints.S), reads(vpoints.sigma), reads(vghost_x_faces.I_3), reads(vghost_y_faces.I_3), reads(vghost_z_faces.I_3), reads(vpoints.I_3), reads(vx_faces.I_3), reads(vy_faces.I_3), reads(vz_faces.I_3), writes(vpoints.I_3), writes(vx_faces.I_3), writes(vy_faces.I_3), writes(vz_faces.I_3)
do
  var vdAx = (vdy*vdz)
  var vdAy = (vdx*vdz)
  var vdAz = (vdx*vdy)
  var vdV = ((vdx*vdy)*vdz)
  var vlimits = vpoints.bounds
  var vdindx = int64(1)
  var vstartx = vlimits.lo.x
  var vendx = (vlimits.hi.x+1)
  var vdindy = int64(1)
  var vstarty = vlimits.lo.y
  var vendy = (vlimits.hi.y+1)
  var vdindz = int64(1)
  var vstartz = vlimits.lo.z
  var vendz = (vlimits.hi.z+1)
  if (vxi<0) then
    vdindx = int64((-1))
    vstartx = vlimits.hi.x
    vendx = (vlimits.lo.x-1)
  else
  end
  if (veta<0) then
    vdindy = int64((-1))
    vstarty = vlimits.hi.y
    vendy = (vlimits.lo.y-1)
  else
  end
  if (vmu<0) then
    vdindz = int64((-1))
    vstartz = vlimits.hi.z
    vendz = (vlimits.lo.z-1)
  else
  end
  for vm : int32 = 0, NUM_ANGLES do
    if (((((vangles[vm].xi*vxi)>0) or ((vangles[vm].xi==0) and (vxi<0))) and (((vangles[vm].eta*veta)>0) or ((vangles[vm].eta==0) and (veta<0)))) and (((vangles[vm].mu*vmu)>0) or ((vangles[vm].mu==0) and (vmu<0)))) then
      for vk : int64 = vstartz, vendz, vdindz do
        for vj : int64 = vstarty, vendy, vdindy do
          for vi : int64 = vstartx, vendx, vdindx do
            var vindx = (vi-min(vdindx, 0))
            var vindy = (vj-min(vdindy, 0))
            var vindz = (vk-min(vdindz, 0))
            var vghost_x_limits = vghost_x_faces.bounds
            var vupwind_x_value = 0
            if (vindx<vx_faces.bounds.lo.x) then
              vupwind_x_value = vghost_x_faces[{vghost_x_limits.hi.x, vj, vk}].I_3[vm]
            elseif (vindx>vx_faces.bounds.hi.x) then
              vupwind_x_value = vghost_x_faces[{vghost_x_limits.lo.x, vj, vk}].I_3[vm]
            else
              vupwind_x_value = vx_faces[{vindx, vj, vk}].I_3[vm]
            end
            var vghost_y_limits = vghost_y_faces.bounds
            var vupwind_y_value = 0
            if (vindy<vy_faces.bounds.lo.y) then
              vupwind_y_value = vghost_y_faces[{vi, vghost_y_limits.hi.y, vk}].I_3[vm]
            elseif (vindy>vy_faces.bounds.hi.y) then
              vupwind_y_value = vghost_y_faces[{vi, vghost_y_limits.lo.y, vk}].I_3[vm]
            else
              vupwind_y_value = vy_faces[{vi, vindy, vk}].I_3[vm]
            end
            var vghost_z_limits = vghost_z_faces.bounds
            var vupwind_z_value = 0
            if (vindz<vz_faces.bounds.lo.z) then
              vupwind_z_value = vghost_z_faces[{vi, vj, vghost_z_limits.hi.z}].I_3[vm]
            elseif (vindz>vz_faces.bounds.hi.z) then
              vupwind_z_value = vghost_z_faces[{vi, vj, vghost_z_limits.lo.z}].I_3[vm]
            else
              vupwind_z_value = vz_faces[{vi, vj, vindz}].I_3[vm]
            end
            vpoints[{vi, vj, vk}].I_3[vm] = (((((vpoints[{vi, vj, vk}].S*vdV)+(((C.fabs(vangles[vm].xi)*vdAx)*vupwind_x_value)/0.5))+(((C.fabs(vangles[vm].eta)*vdAy)*vupwind_y_value)/0.5))+(((C.fabs(vangles[vm].mu)*vdAz)*vupwind_z_value)/0.5))/((((vpoints[{vi, vj, vk}].sigma*vdV)+((C.fabs(vangles[vm].xi)*vdAx)/0.5))+((C.fabs(vangles[vm].eta)*vdAy)/0.5))+((C.fabs(vangles[vm].mu)*vdAz)/0.5)))
            vx_faces[{(vindx+vdindx), vj, vk}].I_3[vm] = ((vpoints[{vi, vj, vk}].I_3[vm]-((1-0.5)*vupwind_x_value))/0.5)
            vy_faces[{vi, (vindy+vdindy), vk}].I_3[vm] = ((vpoints[{vi, vj, vk}].I_3[vm]-((1-0.5)*vupwind_y_value))/0.5)
            vz_faces[{vi, vj, (vindz+vdindz)}].I_3[vm] = ((vpoints[{vi, vj, vk}].I_3[vm]-((1-0.5)*vupwind_z_value))/0.5)
          end
        end
      end
    else
    end
  end
end

task sweep_4(vpoints : region(ispace(int3d), Radiation_columns), vx_faces : region(ispace(int3d), face), vy_faces : region(ispace(int3d), face), vz_faces : region(ispace(int3d), face), vghost_x_faces : region(ispace(int3d), face), vghost_y_faces : region(ispace(int3d), face), vghost_z_faces : region(ispace(int3d), face), vangles : region(ispace(int1d), angle), vxi : int64, veta : int64, vmu : int64, vdx : double, vdy : double, vdz : double)
-- leaf (false), inner (false), idempotent (false)
where
  reads(vangles.xi), reads(vangles.eta), reads(vangles.mu), reads(vpoints.S), reads(vpoints.sigma), reads(vghost_x_faces.I_4), reads(vghost_y_faces.I_4), reads(vghost_z_faces.I_4), reads(vpoints.I_4), reads(vx_faces.I_4), reads(vy_faces.I_4), reads(vz_faces.I_4), writes(vpoints.I_4), writes(vx_faces.I_4), writes(vy_faces.I_4), writes(vz_faces.I_4)
do
  var vdAx = (vdy*vdz)
  var vdAy = (vdx*vdz)
  var vdAz = (vdx*vdy)
  var vdV = ((vdx*vdy)*vdz)
  var vlimits = vpoints.bounds
  var vdindx = int64(1)
  var vstartx = vlimits.lo.x
  var vendx = (vlimits.hi.x+1)
  var vdindy = int64(1)
  var vstarty = vlimits.lo.y
  var vendy = (vlimits.hi.y+1)
  var vdindz = int64(1)
  var vstartz = vlimits.lo.z
  var vendz = (vlimits.hi.z+1)
  if (vxi<0) then
    vdindx = int64((-1))
    vstartx = vlimits.hi.x
    vendx = (vlimits.lo.x-1)
  else
  end
  if (veta<0) then
    vdindy = int64((-1))
    vstarty = vlimits.hi.y
    vendy = (vlimits.lo.y-1)
  else
  end
  if (vmu<0) then
    vdindz = int64((-1))
    vstartz = vlimits.hi.z
    vendz = (vlimits.lo.z-1)
  else
  end
  for vm : int32 = 0, NUM_ANGLES do
    if (((((vangles[vm].xi*vxi)>0) or ((vangles[vm].xi==0) and (vxi<0))) and (((vangles[vm].eta*veta)>0) or ((vangles[vm].eta==0) and (veta<0)))) and (((vangles[vm].mu*vmu)>0) or ((vangles[vm].mu==0) and (vmu<0)))) then
      for vk : int64 = vstartz, vendz, vdindz do
        for vj : int64 = vstarty, vendy, vdindy do
          for vi : int64 = vstartx, vendx, vdindx do
            var vindx = (vi-min(vdindx, 0))
            var vindy = (vj-min(vdindy, 0))
            var vindz = (vk-min(vdindz, 0))
            var vghost_x_limits = vghost_x_faces.bounds
            var vupwind_x_value = 0
            if (vindx<vx_faces.bounds.lo.x) then
              vupwind_x_value = vghost_x_faces[{vghost_x_limits.hi.x, vj, vk}].I_4[vm]
            elseif (vindx>vx_faces.bounds.hi.x) then
              vupwind_x_value = vghost_x_faces[{vghost_x_limits.lo.x, vj, vk}].I_4[vm]
            else
              vupwind_x_value = vx_faces[{vindx, vj, vk}].I_4[vm]
            end
            var vghost_y_limits = vghost_y_faces.bounds
            var vupwind_y_value = 0
            if (vindy<vy_faces.bounds.lo.y) then
              vupwind_y_value = vghost_y_faces[{vi, vghost_y_limits.hi.y, vk}].I_4[vm]
            elseif (vindy>vy_faces.bounds.hi.y) then
              vupwind_y_value = vghost_y_faces[{vi, vghost_y_limits.lo.y, vk}].I_4[vm]
            else
              vupwind_y_value = vy_faces[{vi, vindy, vk}].I_4[vm]
            end
            var vghost_z_limits = vghost_z_faces.bounds
            var vupwind_z_value = 0
            if (vindz<vz_faces.bounds.lo.z) then
              vupwind_z_value = vghost_z_faces[{vi, vj, vghost_z_limits.hi.z}].I_4[vm]
            elseif (vindz>vz_faces.bounds.hi.z) then
              vupwind_z_value = vghost_z_faces[{vi, vj, vghost_z_limits.lo.z}].I_4[vm]
            else
              vupwind_z_value = vz_faces[{vi, vj, vindz}].I_4[vm]
            end
            vpoints[{vi, vj, vk}].I_4[vm] = (((((vpoints[{vi, vj, vk}].S*vdV)+(((C.fabs(vangles[vm].xi)*vdAx)*vupwind_x_value)/0.5))+(((C.fabs(vangles[vm].eta)*vdAy)*vupwind_y_value)/0.5))+(((C.fabs(vangles[vm].mu)*vdAz)*vupwind_z_value)/0.5))/((((vpoints[{vi, vj, vk}].sigma*vdV)+((C.fabs(vangles[vm].xi)*vdAx)/0.5))+((C.fabs(vangles[vm].eta)*vdAy)/0.5))+((C.fabs(vangles[vm].mu)*vdAz)/0.5)))
            vx_faces[{(vindx+vdindx), vj, vk}].I_4[vm] = ((vpoints[{vi, vj, vk}].I_4[vm]-((1-0.5)*vupwind_x_value))/0.5)
            vy_faces[{vi, (vindy+vdindy), vk}].I_4[vm] = ((vpoints[{vi, vj, vk}].I_4[vm]-((1-0.5)*vupwind_y_value))/0.5)
            vz_faces[{vi, vj, (vindz+vdindz)}].I_4[vm] = ((vpoints[{vi, vj, vk}].I_4[vm]-((1-0.5)*vupwind_z_value))/0.5)
          end
        end
      end
    else
    end
  end
end

task sweep_5(vpoints : region(ispace(int3d), Radiation_columns), vx_faces : region(ispace(int3d), face), vy_faces : region(ispace(int3d), face), vz_faces : region(ispace(int3d), face), vghost_x_faces : region(ispace(int3d), face), vghost_y_faces : region(ispace(int3d), face), vghost_z_faces : region(ispace(int3d), face), vangles : region(ispace(int1d), angle), vxi : int64, veta : int64, vmu : int64, vdx : double, vdy : double, vdz : double)
-- leaf (false), inner (false), idempotent (false)
where
  reads(vangles.xi), reads(vangles.eta), reads(vangles.mu), reads(vpoints.S), reads(vpoints.sigma), reads(vghost_x_faces.I_5), reads(vghost_y_faces.I_5), reads(vghost_z_faces.I_5), reads(vpoints.I_5), reads(vx_faces.I_5), reads(vy_faces.I_5), reads(vz_faces.I_5), writes(vpoints.I_5), writes(vx_faces.I_5), writes(vy_faces.I_5), writes(vz_faces.I_5)
do
  var vdAx = (vdy*vdz)
  var vdAy = (vdx*vdz)
  var vdAz = (vdx*vdy)
  var vdV = ((vdx*vdy)*vdz)
  var vlimits = vpoints.bounds
  var vdindx = int64(1)
  var vstartx = vlimits.lo.x
  var vendx = (vlimits.hi.x+1)
  var vdindy = int64(1)
  var vstarty = vlimits.lo.y
  var vendy = (vlimits.hi.y+1)
  var vdindz = int64(1)
  var vstartz = vlimits.lo.z
  var vendz = (vlimits.hi.z+1)
  if (vxi<0) then
    vdindx = int64((-1))
    vstartx = vlimits.hi.x
    vendx = (vlimits.lo.x-1)
  else
  end
  if (veta<0) then
    vdindy = int64((-1))
    vstarty = vlimits.hi.y
    vendy = (vlimits.lo.y-1)
  else
  end
  if (vmu<0) then
    vdindz = int64((-1))
    vstartz = vlimits.hi.z
    vendz = (vlimits.lo.z-1)
  else
  end
  for vm : int32 = 0, NUM_ANGLES do
    if (((((vangles[vm].xi*vxi)>0) or ((vangles[vm].xi==0) and (vxi<0))) and (((vangles[vm].eta*veta)>0) or ((vangles[vm].eta==0) and (veta<0)))) and (((vangles[vm].mu*vmu)>0) or ((vangles[vm].mu==0) and (vmu<0)))) then
      for vk : int64 = vstartz, vendz, vdindz do
        for vj : int64 = vstarty, vendy, vdindy do
          for vi : int64 = vstartx, vendx, vdindx do
            var vindx = (vi-min(vdindx, 0))
            var vindy = (vj-min(vdindy, 0))
            var vindz = (vk-min(vdindz, 0))
            var vghost_x_limits = vghost_x_faces.bounds
            var vupwind_x_value = 0
            if (vindx<vx_faces.bounds.lo.x) then
              vupwind_x_value = vghost_x_faces[{vghost_x_limits.hi.x, vj, vk}].I_5[vm]
            elseif (vindx>vx_faces.bounds.hi.x) then
              vupwind_x_value = vghost_x_faces[{vghost_x_limits.lo.x, vj, vk}].I_5[vm]
            else
              vupwind_x_value = vx_faces[{vindx, vj, vk}].I_5[vm]
            end
            var vghost_y_limits = vghost_y_faces.bounds
            var vupwind_y_value = 0
            if (vindy<vy_faces.bounds.lo.y) then
              vupwind_y_value = vghost_y_faces[{vi, vghost_y_limits.hi.y, vk}].I_5[vm]
            elseif (vindy>vy_faces.bounds.hi.y) then
              vupwind_y_value = vghost_y_faces[{vi, vghost_y_limits.lo.y, vk}].I_5[vm]
            else
              vupwind_y_value = vy_faces[{vi, vindy, vk}].I_5[vm]
            end
            var vghost_z_limits = vghost_z_faces.bounds
            var vupwind_z_value = 0
            if (vindz<vz_faces.bounds.lo.z) then
              vupwind_z_value = vghost_z_faces[{vi, vj, vghost_z_limits.hi.z}].I_5[vm]
            elseif (vindz>vz_faces.bounds.hi.z) then
              vupwind_z_value = vghost_z_faces[{vi, vj, vghost_z_limits.lo.z}].I_5[vm]
            else
              vupwind_z_value = vz_faces[{vi, vj, vindz}].I_5[vm]
            end
            vpoints[{vi, vj, vk}].I_5[vm] = (((((vpoints[{vi, vj, vk}].S*vdV)+(((C.fabs(vangles[vm].xi)*vdAx)*vupwind_x_value)/0.5))+(((C.fabs(vangles[vm].eta)*vdAy)*vupwind_y_value)/0.5))+(((C.fabs(vangles[vm].mu)*vdAz)*vupwind_z_value)/0.5))/((((vpoints[{vi, vj, vk}].sigma*vdV)+((C.fabs(vangles[vm].xi)*vdAx)/0.5))+((C.fabs(vangles[vm].eta)*vdAy)/0.5))+((C.fabs(vangles[vm].mu)*vdAz)/0.5)))
            vx_faces[{(vindx+vdindx), vj, vk}].I_5[vm] = ((vpoints[{vi, vj, vk}].I_5[vm]-((1-0.5)*vupwind_x_value))/0.5)
            vy_faces[{vi, (vindy+vdindy), vk}].I_5[vm] = ((vpoints[{vi, vj, vk}].I_5[vm]-((1-0.5)*vupwind_y_value))/0.5)
            vz_faces[{vi, vj, (vindz+vdindz)}].I_5[vm] = ((vpoints[{vi, vj, vk}].I_5[vm]-((1-0.5)*vupwind_z_value))/0.5)
          end
        end
      end
    else
    end
  end
end

task sweep_6(vpoints : region(ispace(int3d), Radiation_columns), vx_faces : region(ispace(int3d), face), vy_faces : region(ispace(int3d), face), vz_faces : region(ispace(int3d), face), vghost_x_faces : region(ispace(int3d), face), vghost_y_faces : region(ispace(int3d), face), vghost_z_faces : region(ispace(int3d), face), vangles : region(ispace(int1d), angle), vxi : int64, veta : int64, vmu : int64, vdx : double, vdy : double, vdz : double)
-- leaf (false), inner (false), idempotent (false)
where
  reads(vangles.xi), reads(vangles.eta), reads(vangles.mu), reads(vpoints.S), reads(vpoints.sigma), reads(vghost_x_faces.I_6), reads(vghost_y_faces.I_6), reads(vghost_z_faces.I_6), reads(vpoints.I_6), reads(vx_faces.I_6), reads(vy_faces.I_6), reads(vz_faces.I_6), writes(vpoints.I_6), writes(vx_faces.I_6), writes(vy_faces.I_6), writes(vz_faces.I_6)
do
  var vdAx = (vdy*vdz)
  var vdAy = (vdx*vdz)
  var vdAz = (vdx*vdy)
  var vdV = ((vdx*vdy)*vdz)
  var vlimits = vpoints.bounds
  var vdindx = int64(1)
  var vstartx = vlimits.lo.x
  var vendx = (vlimits.hi.x+1)
  var vdindy = int64(1)
  var vstarty = vlimits.lo.y
  var vendy = (vlimits.hi.y+1)
  var vdindz = int64(1)
  var vstartz = vlimits.lo.z
  var vendz = (vlimits.hi.z+1)
  if (vxi<0) then
    vdindx = int64((-1))
    vstartx = vlimits.hi.x
    vendx = (vlimits.lo.x-1)
  else
  end
  if (veta<0) then
    vdindy = int64((-1))
    vstarty = vlimits.hi.y
    vendy = (vlimits.lo.y-1)
  else
  end
  if (vmu<0) then
    vdindz = int64((-1))
    vstartz = vlimits.hi.z
    vendz = (vlimits.lo.z-1)
  else
  end
  for vm : int32 = 0, NUM_ANGLES do
    if (((((vangles[vm].xi*vxi)>0) or ((vangles[vm].xi==0) and (vxi<0))) and (((vangles[vm].eta*veta)>0) or ((vangles[vm].eta==0) and (veta<0)))) and (((vangles[vm].mu*vmu)>0) or ((vangles[vm].mu==0) and (vmu<0)))) then
      for vk : int64 = vstartz, vendz, vdindz do
        for vj : int64 = vstarty, vendy, vdindy do
          for vi : int64 = vstartx, vendx, vdindx do
            var vindx = (vi-min(vdindx, 0))
            var vindy = (vj-min(vdindy, 0))
            var vindz = (vk-min(vdindz, 0))
            var vghost_x_limits = vghost_x_faces.bounds
            var vupwind_x_value = 0
            if (vindx<vx_faces.bounds.lo.x) then
              vupwind_x_value = vghost_x_faces[{vghost_x_limits.hi.x, vj, vk}].I_6[vm]
            elseif (vindx>vx_faces.bounds.hi.x) then
              vupwind_x_value = vghost_x_faces[{vghost_x_limits.lo.x, vj, vk}].I_6[vm]
            else
              vupwind_x_value = vx_faces[{vindx, vj, vk}].I_6[vm]
            end
            var vghost_y_limits = vghost_y_faces.bounds
            var vupwind_y_value = 0
            if (vindy<vy_faces.bounds.lo.y) then
              vupwind_y_value = vghost_y_faces[{vi, vghost_y_limits.hi.y, vk}].I_6[vm]
            elseif (vindy>vy_faces.bounds.hi.y) then
              vupwind_y_value = vghost_y_faces[{vi, vghost_y_limits.lo.y, vk}].I_6[vm]
            else
              vupwind_y_value = vy_faces[{vi, vindy, vk}].I_6[vm]
            end
            var vghost_z_limits = vghost_z_faces.bounds
            var vupwind_z_value = 0
            if (vindz<vz_faces.bounds.lo.z) then
              vupwind_z_value = vghost_z_faces[{vi, vj, vghost_z_limits.hi.z}].I_6[vm]
            elseif (vindz>vz_faces.bounds.hi.z) then
              vupwind_z_value = vghost_z_faces[{vi, vj, vghost_z_limits.lo.z}].I_6[vm]
            else
              vupwind_z_value = vz_faces[{vi, vj, vindz}].I_6[vm]
            end
            vpoints[{vi, vj, vk}].I_6[vm] = (((((vpoints[{vi, vj, vk}].S*vdV)+(((C.fabs(vangles[vm].xi)*vdAx)*vupwind_x_value)/0.5))+(((C.fabs(vangles[vm].eta)*vdAy)*vupwind_y_value)/0.5))+(((C.fabs(vangles[vm].mu)*vdAz)*vupwind_z_value)/0.5))/((((vpoints[{vi, vj, vk}].sigma*vdV)+((C.fabs(vangles[vm].xi)*vdAx)/0.5))+((C.fabs(vangles[vm].eta)*vdAy)/0.5))+((C.fabs(vangles[vm].mu)*vdAz)/0.5)))
            vx_faces[{(vindx+vdindx), vj, vk}].I_6[vm] = ((vpoints[{vi, vj, vk}].I_6[vm]-((1-0.5)*vupwind_x_value))/0.5)
            vy_faces[{vi, (vindy+vdindy), vk}].I_6[vm] = ((vpoints[{vi, vj, vk}].I_6[vm]-((1-0.5)*vupwind_y_value))/0.5)
            vz_faces[{vi, vj, (vindz+vdindz)}].I_6[vm] = ((vpoints[{vi, vj, vk}].I_6[vm]-((1-0.5)*vupwind_z_value))/0.5)
          end
        end
      end
    else
    end
  end
end

task sweep_7(vpoints : region(ispace(int3d), Radiation_columns), vx_faces : region(ispace(int3d), face), vy_faces : region(ispace(int3d), face), vz_faces : region(ispace(int3d), face), vghost_x_faces : region(ispace(int3d), face), vghost_y_faces : region(ispace(int3d), face), vghost_z_faces : region(ispace(int3d), face), vangles : region(ispace(int1d), angle), vxi : int64, veta : int64, vmu : int64, vdx : double, vdy : double, vdz : double)
-- leaf (false), inner (false), idempotent (false)
where
  reads(vangles.xi), reads(vangles.eta), reads(vangles.mu), reads(vpoints.S), reads(vpoints.sigma), reads(vghost_x_faces.I_7), reads(vghost_y_faces.I_7), reads(vghost_z_faces.I_7), reads(vpoints.I_7), reads(vx_faces.I_7), reads(vy_faces.I_7), reads(vz_faces.I_7), writes(vpoints.I_7), writes(vx_faces.I_7), writes(vy_faces.I_7), writes(vz_faces.I_7)
do
  var vdAx = (vdy*vdz)
  var vdAy = (vdx*vdz)
  var vdAz = (vdx*vdy)
  var vdV = ((vdx*vdy)*vdz)
  var vlimits = vpoints.bounds
  var vdindx = int64(1)
  var vstartx = vlimits.lo.x
  var vendx = (vlimits.hi.x+1)
  var vdindy = int64(1)
  var vstarty = vlimits.lo.y
  var vendy = (vlimits.hi.y+1)
  var vdindz = int64(1)
  var vstartz = vlimits.lo.z
  var vendz = (vlimits.hi.z+1)
  if (vxi<0) then
    vdindx = int64((-1))
    vstartx = vlimits.hi.x
    vendx = (vlimits.lo.x-1)
  else
  end
  if (veta<0) then
    vdindy = int64((-1))
    vstarty = vlimits.hi.y
    vendy = (vlimits.lo.y-1)
  else
  end
  if (vmu<0) then
    vdindz = int64((-1))
    vstartz = vlimits.hi.z
    vendz = (vlimits.lo.z-1)
  else
  end
  for vm : int32 = 0, NUM_ANGLES do
    if (((((vangles[vm].xi*vxi)>0) or ((vangles[vm].xi==0) and (vxi<0))) and (((vangles[vm].eta*veta)>0) or ((vangles[vm].eta==0) and (veta<0)))) and (((vangles[vm].mu*vmu)>0) or ((vangles[vm].mu==0) and (vmu<0)))) then
      for vk : int64 = vstartz, vendz, vdindz do
        for vj : int64 = vstarty, vendy, vdindy do
          for vi : int64 = vstartx, vendx, vdindx do
            var vindx = (vi-min(vdindx, 0))
            var vindy = (vj-min(vdindy, 0))
            var vindz = (vk-min(vdindz, 0))
            var vghost_x_limits = vghost_x_faces.bounds
            var vupwind_x_value = 0
            if (vindx<vx_faces.bounds.lo.x) then
              vupwind_x_value = vghost_x_faces[{vghost_x_limits.hi.x, vj, vk}].I_7[vm]
            elseif (vindx>vx_faces.bounds.hi.x) then
              vupwind_x_value = vghost_x_faces[{vghost_x_limits.lo.x, vj, vk}].I_7[vm]
            else
              vupwind_x_value = vx_faces[{vindx, vj, vk}].I_7[vm]
            end
            var vghost_y_limits = vghost_y_faces.bounds
            var vupwind_y_value = 0
            if (vindy<vy_faces.bounds.lo.y) then
              vupwind_y_value = vghost_y_faces[{vi, vghost_y_limits.hi.y, vk}].I_7[vm]
            elseif (vindy>vy_faces.bounds.hi.y) then
              vupwind_y_value = vghost_y_faces[{vi, vghost_y_limits.lo.y, vk}].I_7[vm]
            else
              vupwind_y_value = vy_faces[{vi, vindy, vk}].I_7[vm]
            end
            var vghost_z_limits = vghost_z_faces.bounds
            var vupwind_z_value = 0
            if (vindz<vz_faces.bounds.lo.z) then
              vupwind_z_value = vghost_z_faces[{vi, vj, vghost_z_limits.hi.z}].I_7[vm]
            elseif (vindz>vz_faces.bounds.hi.z) then
              vupwind_z_value = vghost_z_faces[{vi, vj, vghost_z_limits.lo.z}].I_7[vm]
            else
              vupwind_z_value = vz_faces[{vi, vj, vindz}].I_7[vm]
            end
            vpoints[{vi, vj, vk}].I_7[vm] = (((((vpoints[{vi, vj, vk}].S*vdV)+(((C.fabs(vangles[vm].xi)*vdAx)*vupwind_x_value)/0.5))+(((C.fabs(vangles[vm].eta)*vdAy)*vupwind_y_value)/0.5))+(((C.fabs(vangles[vm].mu)*vdAz)*vupwind_z_value)/0.5))/((((vpoints[{vi, vj, vk}].sigma*vdV)+((C.fabs(vangles[vm].xi)*vdAx)/0.5))+((C.fabs(vangles[vm].eta)*vdAy)/0.5))+((C.fabs(vangles[vm].mu)*vdAz)/0.5)))
            vx_faces[{(vindx+vdindx), vj, vk}].I_7[vm] = ((vpoints[{vi, vj, vk}].I_7[vm]-((1-0.5)*vupwind_x_value))/0.5)
            vy_faces[{vi, (vindy+vdindy), vk}].I_7[vm] = ((vpoints[{vi, vj, vk}].I_7[vm]-((1-0.5)*vupwind_y_value))/0.5)
            vz_faces[{vi, vj, (vindz+vdindz)}].I_7[vm] = ((vpoints[{vi, vj, vk}].I_7[vm]-((1-0.5)*vupwind_z_value))/0.5)
          end
        end
      end
    else
    end
  end
end

task sweep_8(vpoints : region(ispace(int3d), Radiation_columns), vx_faces : region(ispace(int3d), face), vy_faces : region(ispace(int3d), face), vz_faces : region(ispace(int3d), face), vghost_x_faces : region(ispace(int3d), face), vghost_y_faces : region(ispace(int3d), face), vghost_z_faces : region(ispace(int3d), face), vangles : region(ispace(int1d), angle), vxi : int64, veta : int64, vmu : int64, vdx : double, vdy : double, vdz : double)
-- leaf (false), inner (false), idempotent (false)
where
  reads(vangles.xi), reads(vangles.eta), reads(vangles.mu), reads(vpoints.S), reads(vpoints.sigma), reads(vghost_x_faces.I_8), reads(vghost_y_faces.I_8), reads(vghost_z_faces.I_8), reads(vpoints.I_8), reads(vx_faces.I_8), reads(vy_faces.I_8), reads(vz_faces.I_8), writes(vpoints.I_8), writes(vx_faces.I_8), writes(vy_faces.I_8), writes(vz_faces.I_8)
do
  var vdAx = (vdy*vdz)
  var vdAy = (vdx*vdz)
  var vdAz = (vdx*vdy)
  var vdV = ((vdx*vdy)*vdz)
  var vlimits = vpoints.bounds
  var vdindx = int64(1)
  var vstartx = vlimits.lo.x
  var vendx = (vlimits.hi.x+1)
  var vdindy = int64(1)
  var vstarty = vlimits.lo.y
  var vendy = (vlimits.hi.y+1)
  var vdindz = int64(1)
  var vstartz = vlimits.lo.z
  var vendz = (vlimits.hi.z+1)
  if (vxi<0) then
    vdindx = int64((-1))
    vstartx = vlimits.hi.x
    vendx = (vlimits.lo.x-1)
  else
  end
  if (veta<0) then
    vdindy = int64((-1))
    vstarty = vlimits.hi.y
    vendy = (vlimits.lo.y-1)
  else
  end
  if (vmu<0) then
    vdindz = int64((-1))
    vstartz = vlimits.hi.z
    vendz = (vlimits.lo.z-1)
  else
  end
  for vm : int32 = 0, NUM_ANGLES do
    if (((((vangles[vm].xi*vxi)>0) or ((vangles[vm].xi==0) and (vxi<0))) and (((vangles[vm].eta*veta)>0) or ((vangles[vm].eta==0) and (veta<0)))) and (((vangles[vm].mu*vmu)>0) or ((vangles[vm].mu==0) and (vmu<0)))) then
      for vk : int64 = vstartz, vendz, vdindz do
        for vj : int64 = vstarty, vendy, vdindy do
          for vi : int64 = vstartx, vendx, vdindx do
            var vindx = (vi-min(vdindx, 0))
            var vindy = (vj-min(vdindy, 0))
            var vindz = (vk-min(vdindz, 0))
            var vghost_x_limits = vghost_x_faces.bounds
            var vupwind_x_value = 0
            if (vindx<vx_faces.bounds.lo.x) then
              vupwind_x_value = vghost_x_faces[{vghost_x_limits.hi.x, vj, vk}].I_8[vm]
            elseif (vindx>vx_faces.bounds.hi.x) then
              vupwind_x_value = vghost_x_faces[{vghost_x_limits.lo.x, vj, vk}].I_8[vm]
            else
              vupwind_x_value = vx_faces[{vindx, vj, vk}].I_8[vm]
            end
            var vghost_y_limits = vghost_y_faces.bounds
            var vupwind_y_value = 0
            if (vindy<vy_faces.bounds.lo.y) then
              vupwind_y_value = vghost_y_faces[{vi, vghost_y_limits.hi.y, vk}].I_8[vm]
            elseif (vindy>vy_faces.bounds.hi.y) then
              vupwind_y_value = vghost_y_faces[{vi, vghost_y_limits.lo.y, vk}].I_8[vm]
            else
              vupwind_y_value = vy_faces[{vi, vindy, vk}].I_8[vm]
            end
            var vghost_z_limits = vghost_z_faces.bounds
            var vupwind_z_value = 0
            if (vindz<vz_faces.bounds.lo.z) then
              vupwind_z_value = vghost_z_faces[{vi, vj, vghost_z_limits.hi.z}].I_8[vm]
            elseif (vindz>vz_faces.bounds.hi.z) then
              vupwind_z_value = vghost_z_faces[{vi, vj, vghost_z_limits.lo.z}].I_8[vm]
            else
              vupwind_z_value = vz_faces[{vi, vj, vindz}].I_8[vm]
            end
            vpoints[{vi, vj, vk}].I_8[vm] = (((((vpoints[{vi, vj, vk}].S*vdV)+(((C.fabs(vangles[vm].xi)*vdAx)*vupwind_x_value)/0.5))+(((C.fabs(vangles[vm].eta)*vdAy)*vupwind_y_value)/0.5))+(((C.fabs(vangles[vm].mu)*vdAz)*vupwind_z_value)/0.5))/((((vpoints[{vi, vj, vk}].sigma*vdV)+((C.fabs(vangles[vm].xi)*vdAx)/0.5))+((C.fabs(vangles[vm].eta)*vdAy)/0.5))+((C.fabs(vangles[vm].mu)*vdAz)/0.5)))
            vx_faces[{(vindx+vdindx), vj, vk}].I_8[vm] = ((vpoints[{vi, vj, vk}].I_8[vm]-((1-0.5)*vupwind_x_value))/0.5)
            vy_faces[{vi, (vindy+vdindy), vk}].I_8[vm] = ((vpoints[{vi, vj, vk}].I_8[vm]-((1-0.5)*vupwind_y_value))/0.5)
            vz_faces[{vi, vj, (vindz+vdindz)}].I_8[vm] = ((vpoints[{vi, vj, vk}].I_8[vm]-((1-0.5)*vupwind_z_value))/0.5)
          end
        end
      end
    else
    end
  end
end

task residual(vpoints : region(ispace(int3d), Radiation_columns), vNx : int32, vNy : int32, vNz : int32) : double
-- leaf (false), inner (false), idempotent (false)
where
  reads(vpoints.I_1), reads(vpoints.I_2), reads(vpoints.I_3), reads(vpoints.I_4), reads(vpoints.I_5), reads(vpoints.I_6), reads(vpoints.I_7), reads(vpoints.I_8), reads(vpoints.Iiter_1), reads(vpoints.Iiter_2), reads(vpoints.Iiter_3), reads(vpoints.Iiter_4), reads(vpoints.Iiter_5), reads(vpoints.Iiter_6), reads(vpoints.Iiter_7), reads(vpoints.Iiter_8)
do
  var vres = 0
  for vp in vpoints do
    for vm : int32 = 0, NUM_ANGLES do
      if (vp.I_1[vm]>0) then
        vres += (((1/(((vNx*vNy)*vNz)*NUM_ANGLES))*C.pow((vp.I_1[vm]-vp.Iiter_1[vm]), 2))/C.pow(vp.I_1[vm], 2))
      else
      end
      if (vp.I_2[vm]>0) then
        vres += (((1/(((vNx*vNy)*vNz)*NUM_ANGLES))*C.pow((vp.I_2[vm]-vp.Iiter_2[vm]), 2))/C.pow(vp.I_2[vm], 2))
      else
      end
      if (vp.I_3[vm]>0) then
        vres += (((1/(((vNx*vNy)*vNz)*NUM_ANGLES))*C.pow((vp.I_3[vm]-vp.Iiter_3[vm]), 2))/C.pow(vp.I_3[vm], 2))
      else
      end
      if (vp.I_4[vm]>0) then
        vres += (((1/(((vNx*vNy)*vNz)*NUM_ANGLES))*C.pow((vp.I_4[vm]-vp.Iiter_4[vm]), 2))/C.pow(vp.I_4[vm], 2))
      else
      end
      if (vp.I_5[vm]>0) then
        vres += (((1/(((vNx*vNy)*vNz)*NUM_ANGLES))*C.pow((vp.I_5[vm]-vp.Iiter_5[vm]), 2))/C.pow(vp.I_5[vm], 2))
      else
      end
      if (vp.I_6[vm]>0) then
        vres += (((1/(((vNx*vNy)*vNz)*NUM_ANGLES))*C.pow((vp.I_6[vm]-vp.Iiter_6[vm]), 2))/C.pow(vp.I_6[vm], 2))
      else
      end
      if (vp.I_7[vm]>0) then
        vres += (((1/(((vNx*vNy)*vNz)*NUM_ANGLES))*C.pow((vp.I_7[vm]-vp.Iiter_7[vm]), 2))/C.pow(vp.I_7[vm], 2))
      else
      end
      if (vp.I_8[vm]>0) then
        vres += (((1/(((vNx*vNy)*vNz)*NUM_ANGLES))*C.pow((vp.I_8[vm]-vp.Iiter_8[vm]), 2))/C.pow(vp.I_8[vm], 2))
      else
      end
    end
  end
  return vres
end

task update(vpoints : region(ispace(int3d), Radiation_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads(vpoints.I_1), reads(vpoints.I_2), reads(vpoints.I_3), reads(vpoints.I_4), reads(vpoints.I_5), reads(vpoints.I_6), reads(vpoints.I_7), reads(vpoints.I_8), reads(vpoints.Iiter_1), reads(vpoints.Iiter_2), reads(vpoints.Iiter_3), reads(vpoints.Iiter_4), reads(vpoints.Iiter_5), reads(vpoints.Iiter_6), reads(vpoints.Iiter_7), reads(vpoints.Iiter_8), writes(vpoints.Iiter_1), writes(vpoints.Iiter_2), writes(vpoints.Iiter_3), writes(vpoints.Iiter_4), writes(vpoints.Iiter_5), writes(vpoints.Iiter_6), writes(vpoints.Iiter_7), writes(vpoints.Iiter_8)
do
  for vp in vpoints do
    for vm : int32 = 0, NUM_ANGLES do
      vp.Iiter_1[vm] = vp.I_1[vm]
      vp.Iiter_2[vm] = vp.I_2[vm]
      vp.Iiter_3[vm] = vp.I_3[vm]
      vp.Iiter_4[vm] = vp.I_4[vm]
      vp.Iiter_5[vm] = vp.I_5[vm]
      vp.Iiter_6[vm] = vp.I_6[vm]
      vp.Iiter_7[vm] = vp.I_7[vm]
      vp.Iiter_8[vm] = vp.I_8[vm]
    end
  end
end

task reduce_intensity(vpoints : region(ispace(int3d), Radiation_columns), vangles : region(ispace(int1d), angle))
-- leaf (false), inner (false), idempotent (false)
where
  reads(vpoints.I_1), reads(vpoints.I_2), reads(vpoints.I_3), reads(vpoints.I_4), reads(vpoints.I_5), reads(vpoints.I_6), reads(vpoints.I_7), reads(vpoints.I_8), reads(vangles.w), reads(vpoints.G), writes(vpoints.G)
do
  for vp in vpoints do
    for vm : int32 = 0, NUM_ANGLES do
      vp.G += ((((((((vangles[vm].w*vp.I_1[vm])+(vangles[vm].w*vp.I_2[vm]))+(vangles[vm].w*vp.I_3[vm]))+(vangles[vm].w*vp.I_4[vm]))+(vangles[vm].w*vp.I_5[vm]))+(vangles[vm].w*vp.I_6[vm]))+(vangles[vm].w*vp.I_7[vm]))+(vangles[vm].w*vp.I_8[vm]))
    end
  end
end

task particles_initValidField(v1019 : region(ispace(int1d), particles_columns))
where
  writes(v1019.__valid)
do
  for v1021 in v1019 do
    v1021.__valid = false
  end
end

__demand(__parallel)
task SetCoarseningField(vFluid : region(ispace(int3d), Fluid_columns), vGrid_xBnum : int32, vGrid_xNum : int32, vGrid_yBnum : int32, vGrid_yNum : int32, vGrid_zBnum : int32, vGrid_zNum : int32, vRadiation_xBnum : int32, vRadiation_xNum : int32, vRadiation_yBnum : int32, vRadiation_yNum : int32, vRadiation_zBnum : int32, vRadiation_zNum : int32)
-- leaf (false), inner (false), idempotent (false)
where
  reads(vFluid.to_Radiation), writes(vFluid.to_Radiation)
do
  for vf in vFluid do
    var vxFactor = (vGrid_xNum/vRadiation_xNum)
    var vyFactor = (vGrid_yNum/vRadiation_yNum)
    var vzFactor = (vRadiation_zNum/vRadiation_zNum)
    if (not ((((((max(int32((uint64(vGrid_xBnum)-int3d(vf).x)), int32(0))>int32(0)) or (max(int32((int3d(vf).x-uint64(((vGrid_xNum+vGrid_xBnum)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(vGrid_yBnum)-int3d(vf).y)), int32(0))>int32(0))) or (max(int32((int3d(vf).y-uint64(((vGrid_yNum+vGrid_yBnum)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(vGrid_zBnum)-int3d(vf).z)), int32(0))>int32(0))) or (max(int32((int3d(vf).z-uint64(((vGrid_zNum+vGrid_zBnum)-int32(1))))), int32(0))>int32(0)))) then
      vFluid[vf].to_Radiation = int3d({(((int3d(vf).x-uint64(vGrid_xBnum))/uint64(vxFactor))+uint64(vRadiation_xBnum)), (((int3d(vf).y-uint64(vGrid_yBnum))/uint64(vyFactor))+uint64(vRadiation_yBnum)), (((int3d(vf).z-uint64(vGrid_zBnum))/uint64(vzFactor))+uint64(vRadiation_zBnum))})
    else
      vFluid[vf].to_Radiation = int3d({uint64(0ULL), uint64(0ULL), uint64(0ULL)})
    end
  end
end

__demand(__parallel)
task Flow_InitializeCell(v1154 : region(ispace(int3d), Fluid_columns))
where
  reads(v1154.PD), writes(v1154.PD), reads(v1154.centerCoordinates), writes(v1154.centerCoordinates), reads(v1154.convectiveSpectralRadius), writes(v1154.convectiveSpectralRadius), reads(v1154.dissipation), writes(v1154.dissipation), reads(v1154.dissipationFlux), writes(v1154.dissipationFlux), reads(v1154.heatConductionSpectralRadius), writes(v1154.heatConductionSpectralRadius), reads(v1154.kineticEnergy), writes(v1154.kineticEnergy), reads(v1154.pressure), writes(v1154.pressure), reads(v1154.pressureBoundary), writes(v1154.pressureBoundary), reads(v1154.rho), writes(v1154.rho), reads(v1154.rhoBoundary), writes(v1154.rhoBoundary), reads(v1154.rhoEnergy), writes(v1154.rhoEnergy), reads(v1154.rhoEnergyBoundary), writes(v1154.rhoEnergyBoundary), reads(v1154.rhoEnergyFluxX), writes(v1154.rhoEnergyFluxX), reads(v1154.rhoEnergyFluxY), writes(v1154.rhoEnergyFluxY), reads(v1154.rhoEnergyFluxZ), writes(v1154.rhoEnergyFluxZ), reads(v1154.rhoEnergy_new), writes(v1154.rhoEnergy_new), reads(v1154.rhoEnergy_old), writes(v1154.rhoEnergy_old), reads(v1154.rhoEnergy_t), writes(v1154.rhoEnergy_t), reads(v1154.rhoEnthalpy), writes(v1154.rhoEnthalpy), reads(v1154.rhoFluxX), writes(v1154.rhoFluxX), reads(v1154.rhoFluxY), writes(v1154.rhoFluxY), reads(v1154.rhoFluxZ), writes(v1154.rhoFluxZ), reads(v1154.rhoVelocity), writes(v1154.rhoVelocity), reads(v1154.rhoVelocityBoundary), writes(v1154.rhoVelocityBoundary), reads(v1154.rhoVelocityFluxX), writes(v1154.rhoVelocityFluxX), reads(v1154.rhoVelocityFluxY), writes(v1154.rhoVelocityFluxY), reads(v1154.rhoVelocityFluxZ), writes(v1154.rhoVelocityFluxZ), reads(v1154.rhoVelocity_new), writes(v1154.rhoVelocity_new), reads(v1154.rhoVelocity_old), writes(v1154.rhoVelocity_old), reads(v1154.rhoVelocity_t), writes(v1154.rhoVelocity_t), reads(v1154.rho_new), writes(v1154.rho_new), reads(v1154.rho_old), writes(v1154.rho_old), reads(v1154.rho_t), writes(v1154.rho_t), reads(v1154.sgsEddyKappa), writes(v1154.sgsEddyKappa), reads(v1154.sgsEddyViscosity), writes(v1154.sgsEddyViscosity), reads(v1154.sgsEnergy), writes(v1154.sgsEnergy), reads(v1154.temperature), writes(v1154.temperature), reads(v1154.temperatureBoundary), writes(v1154.temperatureBoundary), reads(v1154.velocity), writes(v1154.velocity), reads(v1154.velocityBoundary), writes(v1154.velocityBoundary), reads(v1154.velocityGradientX), writes(v1154.velocityGradientX), reads(v1154.velocityGradientXBoundary), writes(v1154.velocityGradientXBoundary), reads(v1154.velocityGradientY), writes(v1154.velocityGradientY), reads(v1154.velocityGradientYBoundary), writes(v1154.velocityGradientYBoundary), reads(v1154.velocityGradientZ), writes(v1154.velocityGradientZ), reads(v1154.velocityGradientZBoundary), writes(v1154.velocityGradientZBoundary), reads(v1154.viscousSpectralRadius), writes(v1154.viscousSpectralRadius)
do
  for v1157 in v1154 do
    v1154[v1157].rho = double(0)
    v1154[v1157].pressure = double(0)
    v1154[v1157].velocity = [double[3]](array(double(0), double(0), double(0)))
    v1154[v1157].centerCoordinates = [double[3]](array(double(0), double(0), double(0)))
    v1154[v1157].velocityGradientX = [double[3]](array(double(0), double(0), double(0)))
    v1154[v1157].velocityGradientY = [double[3]](array(double(0), double(0), double(0)))
    v1154[v1157].velocityGradientZ = [double[3]](array(double(0), double(0), double(0)))
    v1154[v1157].temperature = double(0)
    v1154[v1157].rhoEnthalpy = double(0)
    v1154[v1157].kineticEnergy = double(0)
    v1154[v1157].sgsEnergy = double(0)
    v1154[v1157].sgsEddyViscosity = double(0)
    v1154[v1157].sgsEddyKappa = double(0)
    v1154[v1157].convectiveSpectralRadius = double(0)
    v1154[v1157].viscousSpectralRadius = double(0)
    v1154[v1157].heatConductionSpectralRadius = double(0)
    v1154[v1157].rhoVelocity = [double[3]](array(double(0), double(0), double(0)))
    v1154[v1157].rhoEnergy = double(0)
    v1154[v1157].rhoBoundary = double(0)
    v1154[v1157].rhoVelocityBoundary = [double[3]](array(double(0), double(0), double(0)))
    v1154[v1157].rhoEnergyBoundary = double(0)
    v1154[v1157].velocityBoundary = [double[3]](array(double(0), double(0), double(0)))
    v1154[v1157].pressureBoundary = double(0)
    v1154[v1157].temperatureBoundary = double(0)
    v1154[v1157].velocityGradientXBoundary = [double[3]](array(double(0), double(0), double(0)))
    v1154[v1157].velocityGradientYBoundary = [double[3]](array(double(0), double(0), double(0)))
    v1154[v1157].velocityGradientZBoundary = [double[3]](array(double(0), double(0), double(0)))
    v1154[v1157].rho_old = double(0)
    v1154[v1157].rhoVelocity_old = [double[3]](array(double(0), double(0), double(0)))
    v1154[v1157].rhoEnergy_old = double(0)
    v1154[v1157].rho_new = double(0)
    v1154[v1157].rhoVelocity_new = [double[3]](array(double(0), double(0), double(0)))
    v1154[v1157].rhoEnergy_new = double(0)
    v1154[v1157].rho_t = double(0)
    v1154[v1157].rhoVelocity_t = [double[3]](array(double(0), double(0), double(0)))
    v1154[v1157].rhoEnergy_t = double(0)
    v1154[v1157].rhoFluxX = double(0)
    v1154[v1157].rhoVelocityFluxX = [double[3]](array(double(0), double(0), double(0)))
    v1154[v1157].rhoEnergyFluxX = double(0)
    v1154[v1157].rhoFluxY = double(0)
    v1154[v1157].rhoVelocityFluxY = [double[3]](array(double(0), double(0), double(0)))
    v1154[v1157].rhoEnergyFluxY = double(0)
    v1154[v1157].rhoFluxZ = double(0)
    v1154[v1157].rhoVelocityFluxZ = [double[3]](array(double(0), double(0), double(0)))
    v1154[v1157].rhoEnergyFluxZ = double(0)
    v1154[v1157].PD = double(0)
    v1154[v1157].dissipation = double(0)
    v1154[v1157].dissipationFlux = double(0)
  end
end

__demand(__parallel)
task Flow_InitializeCenterCoordinates(v1265 : region(ispace(int3d), Fluid_columns), v1267 : int32, v1268 : int32, v1269 : double, v1270 : double, v1271 : int32, v1272 : int32, v1273 : double, v1274 : double, v1275 : int32, v1276 : int32, v1277 : double, v1278 : double)
where
  reads(v1265.centerCoordinates), writes(v1265.centerCoordinates)
do
  for v1284 in v1265 do
    var v1285 = [double[3]](array((v1269+((v1270/double(v1268))*(double((int3d(v1284).x-uint64(v1267)))+double(0.5)))), (v1273+((v1274/double(v1272))*(double((int3d(v1284).y-uint64(v1271)))+double(0.5)))), (v1277+((v1278/double(v1276))*(double((int3d(v1284).z-uint64(v1275)))+double(0.5))))))
    v1265[v1284].centerCoordinates = [double[3]](array(double(v1285[int32(0)]), double(v1285[int32(1)]), double(v1285[int32(2)])))
  end
end

__demand(__parallel)
task Flow_InitializeUniform(v1305 : region(ispace(int3d), Fluid_columns), v1307 : double[5])
where
  reads(v1305.pressure), writes(v1305.pressure), reads(v1305.rho), writes(v1305.rho), reads(v1305.velocity), writes(v1305.velocity)
do
  for v1309 in v1305 do
    v1305[v1309].rho = v1307[int32(0)]
    v1305[v1309].pressure = v1307[int32(1)]
    v1305[v1309].velocity[int32(0)] = v1307[int32(2)]
    v1305[v1309].velocity[int32(1)] = v1307[int32(3)]
    v1305[v1309].velocity[int32(2)] = v1307[int32(4)]
  end
end
terra vs_mul_double_3(va : double[3],vb : double) : double[3]
    return array([&double](va)[0] * vb, [&double](va)[1] * vb, [&double](va)[2] * vb)
end

task Flow_InitializeTaylorGreen2D(v1317 : region(ispace(int3d), Fluid_columns), v1319 : double[5], v1320 : int32, v1321 : int32, v1322 : double, v1323 : double, v1324 : int32, v1325 : int32, v1326 : double, v1327 : double, v1328 : int32, v1329 : int32, v1330 : double, v1331 : double)
where
  reads(v1317.pressure), writes(v1317.pressure), reads(v1317.rho), writes(v1317.rho), reads(v1317.velocity), writes(v1317.velocity)
do
  for v1361 in v1317 do
    var v1362 = v1319[int32(0)]
    var v1363 = v1319[int32(1)]
    var v1364 = v1319[int32(2)]
    var v1365 = [double[3]](array((v1322+((v1323/double(v1321))*(double((int3d(v1361).x-uint64(v1320)))+double(0.5)))), (v1326+((v1327/double(v1325))*(double((int3d(v1361).y-uint64(v1324)))+double(0.5)))), (v1330+((v1331/double(v1329))*(double((int3d(v1361).z-uint64(v1328)))+double(0.5))))))
    var v1366 = int32(0)
    v1317[v1361].rho = v1362
    v1317[v1361].velocity = vs_mul_double_3([double[3]](array((([regentlib.sin(double)](v1365[int32(0)])*[regentlib.cos(double)](v1365[int32(1)]))*[regentlib.cos(double)](v1366)), (((-[regentlib.cos(double)](v1365[int32(0)]))*[regentlib.sin(double)](v1365[int32(1)]))*[regentlib.cos(double)](v1366)), double(int32(0)))), v1364)
    var v1367 = ([regentlib.cos(double)]((double(2)*double(v1366)))+double(2))
    var v1368 = ([regentlib.cos(double)]((double(2)*v1365[int32(0)]))+[regentlib.cos(double)]((double(2)*v1365[int32(1)])))
    v1317[v1361].pressure = (v1363+((((v1362*C.pow(v1364, double(int32(2))))/double(int32(16)))*v1367)*v1368))
  end
end

__demand(__parallel)
task Flow_InitializeTaylorGreen3D(v1417 : region(ispace(int3d), Fluid_columns), v1419 : double[5], v1420 : int32, v1421 : int32, v1422 : double, v1423 : double, v1424 : int32, v1425 : int32, v1426 : double, v1427 : double, v1428 : int32, v1429 : int32, v1430 : double, v1431 : double)
where
  reads(v1417.pressure), writes(v1417.pressure), reads(v1417.rho), writes(v1417.rho), reads(v1417.velocity), writes(v1417.velocity)
do
  for v1457 in v1417 do
    var v1458 = v1419[int32(0)]
    var v1459 = v1419[int32(1)]
    var v1460 = v1419[int32(2)]
    var v1461 = [double[3]](array((v1422+((v1423/double(v1421))*(double((int3d(v1457).x-uint64(v1420)))+double(0.5)))), (v1426+((v1427/double(v1425))*(double((int3d(v1457).y-uint64(v1424)))+double(0.5)))), (v1430+((v1431/double(v1429))*(double((int3d(v1457).z-uint64(v1428)))+double(0.5))))))
    v1417[v1457].rho = v1458
    v1417[v1457].velocity = vs_mul_double_3([double[3]](array((([regentlib.sin(double)](v1461[int32(0)])*[regentlib.cos(double)](v1461[int32(1)]))*[regentlib.cos(double)](v1461[int32(2)])), (((-[regentlib.cos(double)](v1461[int32(0)]))*[regentlib.sin(double)](v1461[int32(1)]))*[regentlib.cos(double)](v1461[int32(2)])), double(int32(0)))), v1460)
    var v1462 = ([regentlib.cos(double)]((double(2)*v1461[int32(2)]))+double(2))
    var v1463 = ([regentlib.cos(double)]((double(2)*v1461[int32(0)]))+[regentlib.cos(double)]((double(2)*v1461[int32(1)])))
    v1417[v1457].pressure = (v1459+((((v1458*C.pow(v1460, double(int32(2))))/double(int32(16)))*v1462)*v1463))
  end
end

task Flow_InitializePerturbed(v1513 : region(ispace(int3d), Fluid_columns), v1515 : double[5])
where
  reads(v1513.pressure), writes(v1513.pressure), reads(v1513.rho), writes(v1513.rho), reads(v1513.velocity), writes(v1513.velocity)
do
  for v1517 in v1513 do
    v1513[v1517].rho = v1515[int32(0)]
    v1513[v1517].pressure = v1515[int32(1)]
    v1513[v1517].velocity[int32(0)] = (v1515[int32(2)]+(((double(C.rand())/2147483647)-double(0.5))*double(10)))
    v1513[v1517].velocity[int32(1)] = (v1515[int32(3)]+(((double(C.rand())/2147483647)-double(0.5))*double(10)))
    v1513[v1517].velocity[int32(2)] = (v1515[int32(4)]+(((double(C.rand())/2147483647)-double(0.5))*double(10)))
  end
end
terra dot_double_3(va : double[3],vb : double[3]) : double
    return [&double](va)[0] * [&double](vb)[0] + [&double](va)[1] * [&double](vb)[1] + [&double](va)[2] * [&double](vb)[2]
end

__demand(__parallel)
task Flow_UpdateConservedFromPrimitive(v1524 : region(ispace(int3d), Fluid_columns), v1526 : double, v1527 : double, v1528 : int32, v1529 : int32, v1530 : int32, v1531 : int32, v1532 : int32, v1533 : int32)
where
  reads(v1524.pressure), reads(v1524.rho), reads(v1524.rhoEnergy), writes(v1524.rhoEnergy), reads(v1524.rhoVelocity), writes(v1524.rhoVelocity), reads(v1524.sgsEnergy), reads(v1524.velocity)
do
  for v1553 in v1524 do
    if (not ((((((max(int32((uint64(v1528)-int3d(v1553).x)), int32(0))>int32(0)) or (max(int32((int3d(v1553).x-uint64(((v1529+v1528)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v1530)-int3d(v1553).y)), int32(0))>int32(0))) or (max(int32((int3d(v1553).y-uint64(((v1531+v1530)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v1532)-int3d(v1553).z)), int32(0))>int32(0))) or (max(int32((int3d(v1553).z-uint64(((v1533+v1532)-int32(1))))), int32(0))>int32(0)))) then
      var v1554 = (v1524[v1553].pressure/(v1527*v1524[v1553].rho))
      var v1555 = v1524[v1553].velocity
      v1524[v1553].rhoVelocity = vs_mul_double_3(v1524[v1553].velocity, v1524[v1553].rho)
      var v1556 = (v1527/(v1526-double(1)))
      v1524[v1553].rhoEnergy = ((v1524[v1553].rho*((v1556*v1554)+(double(0.5)*dot_double_3(v1555, v1555))))+v1524[v1553].sgsEnergy)
    else
    end
  end
end
terra vs_div_double_3(va : double[3],vb : double) : double[3]
    return array([&double](va)[0] / vb, [&double](va)[1] / vb, [&double](va)[2] / vb)
end

__demand(__parallel)
task Flow_UpdateAuxiliaryVelocity(v1571 : region(ispace(int3d), Fluid_columns), v1573 : int32, v1574 : int32, v1575 : int32, v1576 : int32, v1577 : int32, v1578 : int32)
where
  reads(v1571.kineticEnergy), writes(v1571.kineticEnergy), reads(v1571.rho), reads(v1571.rhoVelocity), reads(v1571.velocity), writes(v1571.velocity)
do
  for v1586 in v1571 do
    if (not ((((((max(int32((uint64(v1573)-int3d(v1586).x)), int32(0))>int32(0)) or (max(int32((int3d(v1586).x-uint64(((v1574+v1573)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v1575)-int3d(v1586).y)), int32(0))>int32(0))) or (max(int32((int3d(v1586).y-uint64(((v1576+v1575)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v1577)-int3d(v1586).z)), int32(0))>int32(0))) or (max(int32((int3d(v1586).z-uint64(((v1578+v1577)-int32(1))))), int32(0))>int32(0)))) then
      var v1587 = vs_div_double_3(v1571[v1586].rhoVelocity, v1571[v1586].rho)
      v1571[v1586].velocity = v1587
      v1571[v1586].kineticEnergy = ((double(0.5)*v1571[v1586].rho)*dot_double_3(v1587, v1587))
    else
    end
  end
end
terra vv_mul_double_3(va : double[3],vb : double[3]) : double[3]
    return array([&double](va)[0] * [&double](vb)[0], [&double](va)[1] * [&double](vb)[1], [&double](va)[2] * [&double](vb)[2])
end
terra vv_add_double_3(va : double[3],vb : double[3]) : double[3]
    return array([&double](va)[0] + [&double](vb)[0], [&double](va)[1] + [&double](vb)[1], [&double](va)[2] + [&double](vb)[2])
end

__demand(__parallel)
task Flow_UpdateGhostConservedStep1(v1600 : region(ispace(int3d), Fluid_columns), v1602 : double, v1603 : double[3], v1604 : double, v1605 : double[3], v1606 : double[3], v1607 : double, v1608 : double[3], v1609 : double, v1610 : double[3], v1611 : double[3], v1612 : double, v1613 : double[3], v1614 : double, v1615 : double[3], v1616 : double[3], v1617 : double, v1618 : double, v1619 : int32, v1620 : int32, v1621 : int32, v1622 : int32, v1623 : int32, v1624 : int32)
where
  reads(v1600.pressure), reads(v1600.rho), reads(v1600.rhoBoundary), writes(v1600.rhoBoundary), reads(v1600.rhoEnergyBoundary), writes(v1600.rhoEnergyBoundary), reads(v1600.rhoVelocity), reads(v1600.rhoVelocityBoundary), writes(v1600.rhoVelocityBoundary), reads(v1600.temperature)
do
  for v2088 in v1600 do
    if (max(int32((uint64(v1619)-int3d(v2088).x)), int32(0))>int32(0)) then
      var v2089 = int3d(v2088)
      var v2090 = ((v2088+{1, 0, 0})%v1600.bounds)
      var v2091 = v1606
      var v2092 = v1603
      var v2093 = v1602
      var v2094 = double(double(0))
      var v2095 = double(double(0))
      var v2096 = double(double(0))
      var v2097 = [double[3]](array(double(0), double(0), double(0)))
      var v2098 = (v1618/(v1617-double(1)))
      var v2099 = [double[3]](array(double(0), double(0), double(0)))
      v2099 = vv_add_double_3(vv_mul_double_3(vs_div_double_3(v1600[v2090].rhoVelocity, v1600[v2090].rho), v2091), v2092)
      v2095 = v1600[v2090].temperature
      if (v2093>double(0)) then
        v2095 = v2093
      else
      end
      v2096 = ((double(2)*v2095)-v1600[v2090].temperature)
      v2094 = (v1600[v2090].pressure/(v1618*v2096))
      v1600[v2089].rhoBoundary = v2094
      v1600[v2089].rhoVelocityBoundary = vs_mul_double_3(v2099, v2094)
      v1600[v2089].rhoEnergyBoundary = (v2094*((v2098*v2096)+(double(0.5)*dot_double_3(v2099, v2099))))
    else
    end
    if (max(int32((int3d(v2088).x-uint64(((v1620+v1619)-int32(1))))), int32(0))>int32(0)) then
      var v2100 = int3d(v2088)
      var v2101 = ((v2088+{-1, 0, 0})%v1600.bounds)
      var v2102 = v1606
      var v2103 = v1605
      var v2104 = v1604
      var v2105 = double(double(0))
      var v2106 = double(double(0))
      var v2107 = double(double(0))
      var v2108 = [double[3]](array(double(0), double(0), double(0)))
      var v2109 = (v1618/(v1617-double(1)))
      var v2110 = [double[3]](array(double(0), double(0), double(0)))
      v2110 = vv_add_double_3(vv_mul_double_3(vs_div_double_3(v1600[v2101].rhoVelocity, v1600[v2101].rho), v2102), v2103)
      v2106 = v1600[v2101].temperature
      if (v2104>double(0)) then
        v2106 = v2104
      else
      end
      v2107 = ((double(2)*v2106)-v1600[v2101].temperature)
      v2105 = (v1600[v2101].pressure/(v1618*v2107))
      v1600[v2100].rhoBoundary = v2105
      v1600[v2100].rhoVelocityBoundary = vs_mul_double_3(v2110, v2105)
      v1600[v2100].rhoEnergyBoundary = (v2105*((v2109*v2107)+(double(0.5)*dot_double_3(v2110, v2110))))
    else
    end
    if (max(int32((uint64(v1621)-int3d(v2088).y)), int32(0))>int32(0)) then
      var v2111 = int3d(v2088)
      var v2112 = ((v2088+{0, 1, 0})%v1600.bounds)
      var v2113 = v1611
      var v2114 = v1608
      var v2115 = v1607
      var v2116 = double(double(0))
      var v2117 = double(double(0))
      var v2118 = double(double(0))
      var v2119 = [double[3]](array(double(0), double(0), double(0)))
      var v2120 = (v1618/(v1617-double(1)))
      var v2121 = [double[3]](array(double(0), double(0), double(0)))
      v2121 = vv_add_double_3(vv_mul_double_3(vs_div_double_3(v1600[v2112].rhoVelocity, v1600[v2112].rho), v2113), v2114)
      v2117 = v1600[v2112].temperature
      if (v2115>double(0)) then
        v2117 = v2115
      else
      end
      v2118 = ((double(2)*v2117)-v1600[v2112].temperature)
      v2116 = (v1600[v2112].pressure/(v1618*v2118))
      v1600[v2111].rhoBoundary = v2116
      v1600[v2111].rhoVelocityBoundary = vs_mul_double_3(v2121, v2116)
      v1600[v2111].rhoEnergyBoundary = (v2116*((v2120*v2118)+(double(0.5)*dot_double_3(v2121, v2121))))
    else
    end
    if (max(int32((int3d(v2088).y-uint64(((v1622+v1621)-int32(1))))), int32(0))>int32(0)) then
      var v2122 = int3d(v2088)
      var v2123 = ((v2088+{0, -1, 0})%v1600.bounds)
      var v2124 = v1611
      var v2125 = v1610
      var v2126 = v1609
      var v2127 = double(double(0))
      var v2128 = double(double(0))
      var v2129 = double(double(0))
      var v2130 = [double[3]](array(double(0), double(0), double(0)))
      var v2131 = (v1618/(v1617-double(1)))
      var v2132 = [double[3]](array(double(0), double(0), double(0)))
      v2132 = vv_add_double_3(vv_mul_double_3(vs_div_double_3(v1600[v2123].rhoVelocity, v1600[v2123].rho), v2124), v2125)
      v2128 = v1600[v2123].temperature
      if (v2126>double(0)) then
        v2128 = v2126
      else
      end
      v2129 = ((double(2)*v2128)-v1600[v2123].temperature)
      v2127 = (v1600[v2123].pressure/(v1618*v2129))
      v1600[v2122].rhoBoundary = v2127
      v1600[v2122].rhoVelocityBoundary = vs_mul_double_3(v2132, v2127)
      v1600[v2122].rhoEnergyBoundary = (v2127*((v2131*v2129)+(double(0.5)*dot_double_3(v2132, v2132))))
    else
    end
    if (max(int32((uint64(v1623)-int3d(v2088).z)), int32(0))>int32(0)) then
      var v2133 = int3d(v2088)
      var v2134 = ((v2088+{0, 0, 1})%v1600.bounds)
      var v2135 = v1616
      var v2136 = v1613
      var v2137 = v1612
      var v2138 = double(double(0))
      var v2139 = double(double(0))
      var v2140 = double(double(0))
      var v2141 = [double[3]](array(double(0), double(0), double(0)))
      var v2142 = (v1618/(v1617-double(1)))
      var v2143 = [double[3]](array(double(0), double(0), double(0)))
      v2143 = vv_add_double_3(vv_mul_double_3(vs_div_double_3(v1600[v2134].rhoVelocity, v1600[v2134].rho), v2135), v2136)
      v2139 = v1600[v2134].temperature
      if (v2137>double(0)) then
        v2139 = v2137
      else
      end
      v2140 = ((double(2)*v2139)-v1600[v2134].temperature)
      v2138 = (v1600[v2134].pressure/(v1618*v2140))
      v1600[v2133].rhoBoundary = v2138
      v1600[v2133].rhoVelocityBoundary = vs_mul_double_3(v2143, v2138)
      v1600[v2133].rhoEnergyBoundary = (v2138*((v2142*v2140)+(double(0.5)*dot_double_3(v2143, v2143))))
    else
    end
    if (max(int32((int3d(v2088).z-uint64(((v1624+v1623)-int32(1))))), int32(0))>int32(0)) then
      var v2144 = int3d(v2088)
      var v2145 = ((v2088+{0, 0, -1})%v1600.bounds)
      var v2146 = v1616
      var v2147 = v1615
      var v2148 = v1614
      var v2149 = double(double(0))
      var v2150 = double(double(0))
      var v2151 = double(double(0))
      var v2152 = [double[3]](array(double(0), double(0), double(0)))
      var v2153 = (v1618/(v1617-double(1)))
      var v2154 = [double[3]](array(double(0), double(0), double(0)))
      v2154 = vv_add_double_3(vv_mul_double_3(vs_div_double_3(v1600[v2145].rhoVelocity, v1600[v2145].rho), v2146), v2147)
      v2150 = v1600[v2145].temperature
      if (v2148>double(0)) then
        v2150 = v2148
      else
      end
      v2151 = ((double(2)*v2150)-v1600[v2145].temperature)
      v2149 = (v1600[v2145].pressure/(v1618*v2151))
      v1600[v2144].rhoBoundary = v2149
      v1600[v2144].rhoVelocityBoundary = vs_mul_double_3(v2154, v2149)
      v1600[v2144].rhoEnergyBoundary = (v2149*((v2153*v2151)+(double(0.5)*dot_double_3(v2154, v2154))))
    else
    end
  end
end

__demand(__parallel)
task Flow_UpdateGhostConservedStep2(v2454 : region(ispace(int3d), Fluid_columns), v2456 : int32, v2457 : int32, v2458 : int32, v2459 : int32, v2460 : int32, v2461 : int32)
where
  reads(v2454.rho), writes(v2454.rho), reads(v2454.rhoBoundary), reads(v2454.rhoEnergy), writes(v2454.rhoEnergy), reads(v2454.rhoEnergyBoundary), reads(v2454.rhoVelocity), writes(v2454.rhoVelocity), reads(v2454.rhoVelocityBoundary)
do
  for v2463 in v2454 do
    if ((((((max(int32((uint64(v2456)-int3d(v2463).x)), int32(0))>int32(0)) or (max(int32((int3d(v2463).x-uint64(((v2457+v2456)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v2458)-int3d(v2463).y)), int32(0))>int32(0))) or (max(int32((int3d(v2463).y-uint64(((v2459+v2458)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v2460)-int3d(v2463).z)), int32(0))>int32(0))) or (max(int32((int3d(v2463).z-uint64(((v2461+v2460)-int32(1))))), int32(0))>int32(0))) then
      v2454[v2463].rho = v2454[v2463].rhoBoundary
      v2454[v2463].rhoVelocity = v2454[v2463].rhoVelocityBoundary
      v2454[v2463].rhoEnergy = v2454[v2463].rhoEnergyBoundary
    else
    end
  end
end

__demand(__parallel)
task Flow_UpdateGhostVelocityStep1(v2469 : region(ispace(int3d), Fluid_columns), v2471 : double[3], v2472 : double[3], v2473 : double[3], v2474 : double[3], v2475 : double[3], v2476 : double[3], v2477 : double[3], v2478 : double[3], v2479 : double[3], v2480 : int32, v2481 : int32, v2482 : int32, v2483 : int32, v2484 : int32, v2485 : int32)
where
  reads(v2469.velocity), reads(v2469.velocityBoundary), writes(v2469.velocityBoundary)
do
  for v2655 in v2469 do
    if (max(int32((uint64(v2480)-int3d(v2655).x)), int32(0))>int32(0)) then
      var v2656 = int3d(v2655)
      var v2657 = ((v2655+{1, 0, 0})%v2469.bounds)
      var v2658 = v2473
      var v2659 = v2471
      v2469[v2656].velocityBoundary = vv_add_double_3(vv_mul_double_3(v2469[v2657].velocity, v2658), v2659)
    else
    end
    if (max(int32((int3d(v2655).x-uint64(((v2481+v2480)-int32(1))))), int32(0))>int32(0)) then
      var v2660 = int3d(v2655)
      var v2661 = ((v2655+{-1, 0, 0})%v2469.bounds)
      var v2662 = v2473
      var v2663 = v2472
      v2469[v2660].velocityBoundary = vv_add_double_3(vv_mul_double_3(v2469[v2661].velocity, v2662), v2663)
    else
    end
    if (max(int32((uint64(v2482)-int3d(v2655).y)), int32(0))>int32(0)) then
      var v2664 = int3d(v2655)
      var v2665 = ((v2655+{0, 1, 0})%v2469.bounds)
      var v2666 = v2476
      var v2667 = v2474
      v2469[v2664].velocityBoundary = vv_add_double_3(vv_mul_double_3(v2469[v2665].velocity, v2666), v2667)
    else
    end
    if (max(int32((int3d(v2655).y-uint64(((v2483+v2482)-int32(1))))), int32(0))>int32(0)) then
      var v2668 = int3d(v2655)
      var v2669 = ((v2655+{0, -1, 0})%v2469.bounds)
      var v2670 = v2476
      var v2671 = v2475
      v2469[v2668].velocityBoundary = vv_add_double_3(vv_mul_double_3(v2469[v2669].velocity, v2670), v2671)
    else
    end
    if (max(int32((uint64(v2484)-int3d(v2655).z)), int32(0))>int32(0)) then
      var v2672 = int3d(v2655)
      var v2673 = ((v2655+{0, 0, 1})%v2469.bounds)
      var v2674 = v2479
      var v2675 = v2477
      v2469[v2672].velocityBoundary = vv_add_double_3(vv_mul_double_3(v2469[v2673].velocity, v2674), v2675)
    else
    end
    if (max(int32((int3d(v2655).z-uint64(((v2485+v2484)-int32(1))))), int32(0))>int32(0)) then
      var v2676 = int3d(v2655)
      var v2677 = ((v2655+{0, 0, -1})%v2469.bounds)
      var v2678 = v2479
      var v2679 = v2478
      v2469[v2676].velocityBoundary = vv_add_double_3(vv_mul_double_3(v2469[v2677].velocity, v2678), v2679)
    else
    end
  end
end

__demand(__parallel)
task Flow_UpdateGhostVelocityStep2(v2769 : region(ispace(int3d), Fluid_columns), v2771 : int32, v2772 : int32, v2773 : int32, v2774 : int32, v2775 : int32, v2776 : int32)
where
  reads(v2769.velocity), writes(v2769.velocity), reads(v2769.velocityBoundary)
do
  for v2778 in v2769 do
    if ((((((max(int32((uint64(v2771)-int3d(v2778).x)), int32(0))>int32(0)) or (max(int32((int3d(v2778).x-uint64(((v2772+v2771)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v2773)-int3d(v2778).y)), int32(0))>int32(0))) or (max(int32((int3d(v2778).y-uint64(((v2774+v2773)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v2775)-int3d(v2778).z)), int32(0))>int32(0))) or (max(int32((int3d(v2778).z-uint64(((v2776+v2775)-int32(1))))), int32(0))>int32(0))) then
      v2769[v2778].velocity = v2769[v2778].velocityBoundary
    else
    end
  end
end
terra vv_sub_double_3(va : double[3],vb : double[3]) : double[3]
    return array([&double](va)[0] - [&double](vb)[0], [&double](va)[1] - [&double](vb)[1], [&double](va)[2] - [&double](vb)[2])
end

__demand(__parallel)
task Flow_ComputeVelocityGradientAll(v2784 : region(ispace(int3d), Fluid_columns), v2786 : int32, v2787 : double, v2788 : int32, v2789 : int32, v2790 : double, v2791 : int32, v2792 : int32, v2793 : double, v2794 : int32)
where
  reads(v2784.velocity), reads(v2784.velocityGradientX), writes(v2784.velocityGradientX), reads(v2784.velocityGradientY), writes(v2784.velocityGradientY), reads(v2784.velocityGradientZ), writes(v2784.velocityGradientZ)
do
  for v2796 in v2784 do
    if (not ((((((max(int32((uint64(v2786)-int3d(v2796).x)), int32(0))>int32(0)) or (max(int32((int3d(v2796).x-uint64(((v2788+v2786)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v2789)-int3d(v2796).y)), int32(0))>int32(0))) or (max(int32((int3d(v2796).y-uint64(((v2791+v2789)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v2792)-int3d(v2796).z)), int32(0))>int32(0))) or (max(int32((int3d(v2796).z-uint64(((v2794+v2792)-int32(1))))), int32(0))>int32(0)))) then
      v2784[v2796].velocityGradientX = vs_div_double_3(vs_mul_double_3(vv_sub_double_3(v2784[((v2796+{1, 0, 0})%v2784.bounds)].velocity, v2784[((v2796+{-1, 0, 0})%v2784.bounds)].velocity), double(0.5)), v2787)
      v2784[v2796].velocityGradientY = vs_div_double_3(vs_mul_double_3(vv_sub_double_3(v2784[((v2796+{0, 1, 0})%v2784.bounds)].velocity, v2784[((v2796+{0, -1, 0})%v2784.bounds)].velocity), double(0.5)), v2790)
      v2784[v2796].velocityGradientZ = vs_div_double_3(vs_mul_double_3(vv_sub_double_3(v2784[((v2796+{0, 0, 1})%v2784.bounds)].velocity, v2784[((v2796+{0, 0, -1})%v2784.bounds)].velocity), double(0.5)), v2793)
    else
    end
  end
end

__demand(__parallel)
task Flow_UpdateAuxiliaryThermodynamics(v2865 : region(ispace(int3d), Fluid_columns), v2867 : double, v2868 : double, v2869 : int32, v2870 : int32, v2871 : int32, v2872 : int32, v2873 : int32, v2874 : int32)
where
  reads(v2865.pressure), writes(v2865.pressure), reads(v2865.rho), reads(v2865.rhoEnergy), reads(v2865.temperature), writes(v2865.temperature), reads(v2865.velocity)
do
  for v2888 in v2865 do
    if (not ((((((max(int32((uint64(v2869)-int3d(v2888).x)), int32(0))>int32(0)) or (max(int32((int3d(v2888).x-uint64(((v2870+v2869)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v2871)-int3d(v2888).y)), int32(0))>int32(0))) or (max(int32((int3d(v2888).y-uint64(((v2872+v2871)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v2873)-int3d(v2888).z)), int32(0))>int32(0))) or (max(int32((int3d(v2888).z-uint64(((v2874+v2873)-int32(1))))), int32(0))>int32(0)))) then
      var v2889 = ((double(0.5)*v2865[v2888].rho)*dot_double_3(v2865[v2888].velocity, v2865[v2888].velocity))
      var v2890 = ((v2867-double(1))*(v2865[v2888].rhoEnergy-v2889))
      v2865[v2888].pressure = v2890
      v2865[v2888].temperature = (v2890/(v2868*v2865[v2888].rho))
    else
    end
  end
end

__demand(__parallel)
task Flow_UpdateGhostThermodynamicsStep1(v2902 : region(ispace(int3d), Fluid_columns), v2904 : double, v2905 : double, v2906 : double, v2907 : double, v2908 : double, v2909 : double, v2910 : int32, v2911 : int32, v2912 : int32, v2913 : int32, v2914 : int32, v2915 : int32)
where
  reads(v2902.pressure), reads(v2902.pressureBoundary), writes(v2902.pressureBoundary), reads(v2902.temperature), reads(v2902.temperatureBoundary), writes(v2902.temperatureBoundary)
do
  for v3127 in v2902 do
    if (max(int32((uint64(v2910)-int3d(v3127).x)), int32(0))>int32(0)) then
      var v3128 = int3d(v3127)
      var v3129 = ((v3127+{1, 0, 0})%v2902.bounds)
      var v3130 = v2904
      var v3131 = double(double(0))
      var v3132 = double(double(0))
      v3131 = v2902[v3129].temperature
      if (v3130>double(0)) then
        v3131 = v3130
      else
      end
      v3132 = ((double(2)*v3131)-v2902[v3129].temperature)
      v2902[v3128].pressureBoundary = v2902[v3129].pressure
      v2902[v3128].temperatureBoundary = v3132
    else
    end
    if (max(int32((int3d(v3127).x-uint64(((v2911+v2910)-int32(1))))), int32(0))>int32(0)) then
      var v3133 = int3d(v3127)
      var v3134 = ((v3127+{-1, 0, 0})%v2902.bounds)
      var v3135 = v2905
      var v3136 = double(double(0))
      var v3137 = double(double(0))
      v3136 = v2902[v3134].temperature
      if (v3135>double(0)) then
        v3136 = v3135
      else
      end
      v3137 = ((double(2)*v3136)-v2902[v3134].temperature)
      v2902[v3133].pressureBoundary = v2902[v3134].pressure
      v2902[v3133].temperatureBoundary = v3137
    else
    end
    if (max(int32((uint64(v2912)-int3d(v3127).y)), int32(0))>int32(0)) then
      var v3138 = int3d(v3127)
      var v3139 = ((v3127+{0, 1, 0})%v2902.bounds)
      var v3140 = v2906
      var v3141 = double(double(0))
      var v3142 = double(double(0))
      v3141 = v2902[v3139].temperature
      if (v3140>double(0)) then
        v3141 = v3140
      else
      end
      v3142 = ((double(2)*v3141)-v2902[v3139].temperature)
      v2902[v3138].pressureBoundary = v2902[v3139].pressure
      v2902[v3138].temperatureBoundary = v3142
    else
    end
    if (max(int32((int3d(v3127).y-uint64(((v2913+v2912)-int32(1))))), int32(0))>int32(0)) then
      var v3143 = int3d(v3127)
      var v3144 = ((v3127+{0, -1, 0})%v2902.bounds)
      var v3145 = v2907
      var v3146 = double(double(0))
      var v3147 = double(double(0))
      v3146 = v2902[v3144].temperature
      if (v3145>double(0)) then
        v3146 = v3145
      else
      end
      v3147 = ((double(2)*v3146)-v2902[v3144].temperature)
      v2902[v3143].pressureBoundary = v2902[v3144].pressure
      v2902[v3143].temperatureBoundary = v3147
    else
    end
    if (max(int32((uint64(v2914)-int3d(v3127).z)), int32(0))>int32(0)) then
      var v3148 = int3d(v3127)
      var v3149 = ((v3127+{0, 0, 1})%v2902.bounds)
      var v3150 = v2908
      var v3151 = double(double(0))
      var v3152 = double(double(0))
      v3151 = v2902[v3149].temperature
      if (v3150>double(0)) then
        v3151 = v3150
      else
      end
      v3152 = ((double(2)*v3151)-v2902[v3149].temperature)
      v2902[v3148].pressureBoundary = v2902[v3149].pressure
      v2902[v3148].temperatureBoundary = v3152
    else
    end
    if (max(int32((int3d(v3127).z-uint64(((v2915+v2914)-int32(1))))), int32(0))>int32(0)) then
      var v3153 = int3d(v3127)
      var v3154 = ((v3127+{0, 0, -1})%v2902.bounds)
      var v3155 = v2909
      var v3156 = double(double(0))
      var v3157 = double(double(0))
      v3156 = v2902[v3154].temperature
      if (v3155>double(0)) then
        v3156 = v3155
      else
      end
      v3157 = ((double(2)*v3156)-v2902[v3154].temperature)
      v2902[v3153].pressureBoundary = v2902[v3154].pressure
      v2902[v3153].temperatureBoundary = v3157
    else
    end
  end
end

__demand(__parallel)
task Flow_UpdateGhostThermodynamicsStep2(v3241 : region(ispace(int3d), Fluid_columns), v3243 : int32, v3244 : int32, v3245 : int32, v3246 : int32, v3247 : int32, v3248 : int32)
where
  reads(v3241.pressure), writes(v3241.pressure), reads(v3241.pressureBoundary), reads(v3241.temperature), writes(v3241.temperature), reads(v3241.temperatureBoundary)
do
  for v3250 in v3241 do
    if ((((((max(int32((uint64(v3243)-int3d(v3250).x)), int32(0))>int32(0)) or (max(int32((int3d(v3250).x-uint64(((v3244+v3243)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v3245)-int3d(v3250).y)), int32(0))>int32(0))) or (max(int32((int3d(v3250).y-uint64(((v3246+v3245)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v3247)-int3d(v3250).z)), int32(0))>int32(0))) or (max(int32((int3d(v3250).z-uint64(((v3248+v3247)-int32(1))))), int32(0))>int32(0))) then
      v3241[v3250].pressure = v3241[v3250].pressureBoundary
      v3241[v3250].temperature = v3241[v3250].temperatureBoundary
    else
    end
  end
end

__demand(__parallel)
task Flow_UpdateGhostFieldsStep1(v3256 : region(ispace(int3d), Fluid_columns), v3258 : double, v3259 : double[3], v3260 : double, v3261 : double[3], v3262 : double[3], v3263 : double, v3264 : double[3], v3265 : double, v3266 : double[3], v3267 : double[3], v3268 : double, v3269 : double[3], v3270 : double, v3271 : double[3], v3272 : double[3], v3273 : double, v3274 : double, v3275 : int32, v3276 : int32, v3277 : int32, v3278 : int32, v3279 : int32, v3280 : int32)
where
  reads(v3256.pressure), reads(v3256.pressureBoundary), writes(v3256.pressureBoundary), reads(v3256.rho), reads(v3256.rhoBoundary), writes(v3256.rhoBoundary), reads(v3256.rhoEnergyBoundary), writes(v3256.rhoEnergyBoundary), reads(v3256.rhoVelocity), reads(v3256.rhoVelocityBoundary), writes(v3256.rhoVelocityBoundary), reads(v3256.temperature), reads(v3256.temperatureBoundary), writes(v3256.temperatureBoundary), reads(v3256.velocityBoundary), writes(v3256.velocityBoundary)
do
  for v3702 in v3256 do
    if (max(int32((uint64(v3275)-int3d(v3702).x)), int32(0))>int32(0)) then
      var v3703 = int3d(v3702)
      var v3704 = ((v3702+{1, 0, 0})%v3256.bounds)
      var v3705 = v3262
      var v3706 = v3259
      var v3707 = v3258
      var v3708 = double(double(0))
      var v3709 = double(double(0))
      var v3710 = double(double(0))
      var v3711 = [double[3]](array(double(0), double(0), double(0)))
      var v3712 = (v3274/(v3273-double(1)))
      v3711 = vv_add_double_3(vv_mul_double_3(vs_div_double_3(v3256[v3704].rhoVelocity, v3256[v3704].rho), v3705), v3706)
      v3709 = v3256[v3704].temperature
      if (v3707>double(0)) then
        v3709 = v3707
      else
      end
      v3710 = ((double(2)*v3709)-v3256[v3704].temperature)
      v3708 = (v3256[v3704].pressure/(v3274*v3710))
      v3256[v3703].rhoBoundary = v3708
      v3256[v3703].rhoVelocityBoundary = vs_mul_double_3(v3711, v3708)
      v3256[v3703].rhoEnergyBoundary = (v3708*((v3712*v3710)+(double(0.5)*dot_double_3(v3711, v3711))))
      v3256[v3703].velocityBoundary = v3711
      v3256[v3703].pressureBoundary = v3256[v3704].pressure
      v3256[v3703].temperatureBoundary = v3710
    else
    end
    if (max(int32((int3d(v3702).x-uint64(((v3276+v3275)-int32(1))))), int32(0))>int32(0)) then
      var v3713 = int3d(v3702)
      var v3714 = ((v3702+{-1, 0, 0})%v3256.bounds)
      var v3715 = v3262
      var v3716 = v3261
      var v3717 = v3260
      var v3718 = double(double(0))
      var v3719 = double(double(0))
      var v3720 = double(double(0))
      var v3721 = [double[3]](array(double(0), double(0), double(0)))
      var v3722 = (v3274/(v3273-double(1)))
      v3721 = vv_add_double_3(vv_mul_double_3(vs_div_double_3(v3256[v3714].rhoVelocity, v3256[v3714].rho), v3715), v3716)
      v3719 = v3256[v3714].temperature
      if (v3717>double(0)) then
        v3719 = v3717
      else
      end
      v3720 = ((double(2)*v3719)-v3256[v3714].temperature)
      v3718 = (v3256[v3714].pressure/(v3274*v3720))
      v3256[v3713].rhoBoundary = v3718
      v3256[v3713].rhoVelocityBoundary = vs_mul_double_3(v3721, v3718)
      v3256[v3713].rhoEnergyBoundary = (v3718*((v3722*v3720)+(double(0.5)*dot_double_3(v3721, v3721))))
      v3256[v3713].velocityBoundary = v3721
      v3256[v3713].pressureBoundary = v3256[v3714].pressure
      v3256[v3713].temperatureBoundary = v3720
    else
    end
    if (max(int32((uint64(v3277)-int3d(v3702).y)), int32(0))>int32(0)) then
      var v3723 = int3d(v3702)
      var v3724 = ((v3702+{0, 1, 0})%v3256.bounds)
      var v3725 = v3267
      var v3726 = v3264
      var v3727 = v3263
      var v3728 = double(double(0))
      var v3729 = double(double(0))
      var v3730 = double(double(0))
      var v3731 = [double[3]](array(double(0), double(0), double(0)))
      var v3732 = (v3274/(v3273-double(1)))
      v3731 = vv_add_double_3(vv_mul_double_3(vs_div_double_3(v3256[v3724].rhoVelocity, v3256[v3724].rho), v3725), v3726)
      v3729 = v3256[v3724].temperature
      if (v3727>double(0)) then
        v3729 = v3727
      else
      end
      v3730 = ((double(2)*v3729)-v3256[v3724].temperature)
      v3728 = (v3256[v3724].pressure/(v3274*v3730))
      v3256[v3723].rhoBoundary = v3728
      v3256[v3723].rhoVelocityBoundary = vs_mul_double_3(v3731, v3728)
      v3256[v3723].rhoEnergyBoundary = (v3728*((v3732*v3730)+(double(0.5)*dot_double_3(v3731, v3731))))
      v3256[v3723].velocityBoundary = v3731
      v3256[v3723].pressureBoundary = v3256[v3724].pressure
      v3256[v3723].temperatureBoundary = v3730
    else
    end
    if (max(int32((int3d(v3702).y-uint64(((v3278+v3277)-int32(1))))), int32(0))>int32(0)) then
      var v3733 = int3d(v3702)
      var v3734 = ((v3702+{0, -1, 0})%v3256.bounds)
      var v3735 = v3267
      var v3736 = v3266
      var v3737 = v3265
      var v3738 = double(double(0))
      var v3739 = double(double(0))
      var v3740 = double(double(0))
      var v3741 = [double[3]](array(double(0), double(0), double(0)))
      var v3742 = (v3274/(v3273-double(1)))
      v3741 = vv_add_double_3(vv_mul_double_3(vs_div_double_3(v3256[v3734].rhoVelocity, v3256[v3734].rho), v3735), v3736)
      v3739 = v3256[v3734].temperature
      if (v3737>double(0)) then
        v3739 = v3737
      else
      end
      v3740 = ((double(2)*v3739)-v3256[v3734].temperature)
      v3738 = (v3256[v3734].pressure/(v3274*v3740))
      v3256[v3733].rhoBoundary = v3738
      v3256[v3733].rhoVelocityBoundary = vs_mul_double_3(v3741, v3738)
      v3256[v3733].rhoEnergyBoundary = (v3738*((v3742*v3740)+(double(0.5)*dot_double_3(v3741, v3741))))
      v3256[v3733].velocityBoundary = v3741
      v3256[v3733].pressureBoundary = v3256[v3734].pressure
      v3256[v3733].temperatureBoundary = v3740
    else
    end
    if (max(int32((uint64(v3279)-int3d(v3702).z)), int32(0))>int32(0)) then
      var v3743 = int3d(v3702)
      var v3744 = ((v3702+{0, 0, 1})%v3256.bounds)
      var v3745 = v3272
      var v3746 = v3269
      var v3747 = v3268
      var v3748 = double(double(0))
      var v3749 = double(double(0))
      var v3750 = double(double(0))
      var v3751 = [double[3]](array(double(0), double(0), double(0)))
      var v3752 = (v3274/(v3273-double(1)))
      v3751 = vv_add_double_3(vv_mul_double_3(vs_div_double_3(v3256[v3744].rhoVelocity, v3256[v3744].rho), v3745), v3746)
      v3749 = v3256[v3744].temperature
      if (v3747>double(0)) then
        v3749 = v3747
      else
      end
      v3750 = ((double(2)*v3749)-v3256[v3744].temperature)
      v3748 = (v3256[v3744].pressure/(v3274*v3750))
      v3256[v3743].rhoBoundary = v3748
      v3256[v3743].rhoVelocityBoundary = vs_mul_double_3(v3751, v3748)
      v3256[v3743].rhoEnergyBoundary = (v3748*((v3752*v3750)+(double(0.5)*dot_double_3(v3751, v3751))))
      v3256[v3743].velocityBoundary = v3751
      v3256[v3743].pressureBoundary = v3256[v3744].pressure
      v3256[v3743].temperatureBoundary = v3750
    else
    end
    if (max(int32((int3d(v3702).z-uint64(((v3280+v3279)-int32(1))))), int32(0))>int32(0)) then
      var v3753 = int3d(v3702)
      var v3754 = ((v3702+{0, 0, -1})%v3256.bounds)
      var v3755 = v3272
      var v3756 = v3271
      var v3757 = v3270
      var v3758 = double(double(0))
      var v3759 = double(double(0))
      var v3760 = double(double(0))
      var v3761 = [double[3]](array(double(0), double(0), double(0)))
      var v3762 = (v3274/(v3273-double(1)))
      v3761 = vv_add_double_3(vv_mul_double_3(vs_div_double_3(v3256[v3754].rhoVelocity, v3256[v3754].rho), v3755), v3756)
      v3759 = v3256[v3754].temperature
      if (v3757>double(0)) then
        v3759 = v3757
      else
      end
      v3760 = ((double(2)*v3759)-v3256[v3754].temperature)
      v3758 = (v3256[v3754].pressure/(v3274*v3760))
      v3256[v3753].rhoBoundary = v3758
      v3256[v3753].rhoVelocityBoundary = vs_mul_double_3(v3761, v3758)
      v3256[v3753].rhoEnergyBoundary = (v3758*((v3762*v3760)+(double(0.5)*dot_double_3(v3761, v3761))))
      v3256[v3753].velocityBoundary = v3761
      v3256[v3753].pressureBoundary = v3256[v3754].pressure
      v3256[v3753].temperatureBoundary = v3760
    else
    end
  end
end

__demand(__parallel)
task Flow_UpdateGhostFieldsStep2(v4032 : region(ispace(int3d), Fluid_columns), v4034 : int32, v4035 : int32, v4036 : int32, v4037 : int32, v4038 : int32, v4039 : int32)
where
  reads(v4032.pressure), writes(v4032.pressure), reads(v4032.pressureBoundary), reads(v4032.rho), writes(v4032.rho), reads(v4032.rhoBoundary), reads(v4032.rhoEnergy), writes(v4032.rhoEnergy), reads(v4032.rhoEnergyBoundary), reads(v4032.rhoVelocity), writes(v4032.rhoVelocity), reads(v4032.rhoVelocityBoundary), reads(v4032.temperature), writes(v4032.temperature), reads(v4032.temperatureBoundary)
do
  for v4041 in v4032 do
    if ((((((max(int32((uint64(v4034)-int3d(v4041).x)), int32(0))>int32(0)) or (max(int32((int3d(v4041).x-uint64(((v4035+v4034)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v4036)-int3d(v4041).y)), int32(0))>int32(0))) or (max(int32((int3d(v4041).y-uint64(((v4037+v4036)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v4038)-int3d(v4041).z)), int32(0))>int32(0))) or (max(int32((int3d(v4041).z-uint64(((v4039+v4038)-int32(1))))), int32(0))>int32(0))) then
      v4032[v4041].rho = v4032[v4041].rhoBoundary
      v4032[v4041].rhoVelocity = v4032[v4041].rhoVelocityBoundary
      v4032[v4041].rhoEnergy = v4032[v4041].rhoEnergyBoundary
      v4032[v4041].pressure = v4032[v4041].pressureBoundary
      v4032[v4041].temperature = v4032[v4041].temperatureBoundary
    else
    end
  end
end

__demand(__parallel)
task CalculateAveragePressure(v4051 : region(ispace(int3d), Fluid_columns), v4054 : double, v4055 : int32, v4056 : int32, v4057 : int32, v4058 : int32, v4059 : int32, v4060 : int32) : double
where
  reads(v4051.pressure)
do
  var v4063 = double(0)
  for v4064 in v4051 do
    if (not ((((((max(int32((uint64(v4055)-int3d(v4064).x)), int32(0))>int32(0)) or (max(int32((int3d(v4064).x-uint64(((v4056+v4055)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v4057)-int3d(v4064).y)), int32(0))>int32(0))) or (max(int32((int3d(v4064).y-uint64(((v4058+v4057)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v4059)-int3d(v4064).z)), int32(0))>int32(0))) or (max(int32((int3d(v4064).z-uint64(((v4060+v4059)-int32(1))))), int32(0))>int32(0)))) then
      v4063 += (v4051[v4064].pressure*v4054)
    else
    end
  end
  return v4063
end

__demand(__parallel)
task CalculateAverageTemperature(v4071 : region(ispace(int3d), Fluid_columns), v4074 : double, v4075 : int32, v4076 : int32, v4077 : int32, v4078 : int32, v4079 : int32, v4080 : int32) : double
where
  reads(v4071.temperature)
do
  var v4083 = double(0)
  for v4084 in v4071 do
    if (not ((((((max(int32((uint64(v4075)-int3d(v4084).x)), int32(0))>int32(0)) or (max(int32((int3d(v4084).x-uint64(((v4076+v4075)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v4077)-int3d(v4084).y)), int32(0))>int32(0))) or (max(int32((int3d(v4084).y-uint64(((v4078+v4077)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v4079)-int3d(v4084).z)), int32(0))>int32(0))) or (max(int32((int3d(v4084).z-uint64(((v4080+v4079)-int32(1))))), int32(0))>int32(0)))) then
      v4083 += (v4071[v4084].temperature*v4074)
    else
    end
  end
  return v4083
end

__demand(__parallel)
task CalculateAverageKineticEnergy(v4091 : region(ispace(int3d), Fluid_columns), v4094 : double, v4095 : int32, v4096 : int32, v4097 : int32, v4098 : int32, v4099 : int32, v4100 : int32) : double
where
  reads(v4091.kineticEnergy)
do
  var v4103 = double(0)
  for v4104 in v4091 do
    if (not ((((((max(int32((uint64(v4095)-int3d(v4104).x)), int32(0))>int32(0)) or (max(int32((int3d(v4104).x-uint64(((v4096+v4095)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v4097)-int3d(v4104).y)), int32(0))>int32(0))) or (max(int32((int3d(v4104).y-uint64(((v4098+v4097)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v4099)-int3d(v4104).z)), int32(0))>int32(0))) or (max(int32((int3d(v4104).z-uint64(((v4100+v4099)-int32(1))))), int32(0))>int32(0)))) then
      v4103 += (v4091[v4104].kineticEnergy*v4094)
    else
    end
  end
  return v4103
end

__demand(__parallel)
task CalculateMinTemperature(v4111 : region(ispace(int3d), Fluid_columns), v4114 : int32, v4115 : int32, v4116 : int32, v4117 : int32, v4118 : int32, v4119 : int32) : double
where
  reads(v4111.temperature)
do
  var v4122 = math.huge
  for v4123 in v4111 do
    if (not ((((((max(int32((uint64(v4114)-int3d(v4123).x)), int32(0))>int32(0)) or (max(int32((int3d(v4123).x-uint64(((v4115+v4114)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v4116)-int3d(v4123).y)), int32(0))>int32(0))) or (max(int32((int3d(v4123).y-uint64(((v4117+v4116)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v4118)-int3d(v4123).z)), int32(0))>int32(0))) or (max(int32((int3d(v4123).z-uint64(((v4119+v4118)-int32(1))))), int32(0))>int32(0)))) then
      v4122 min= v4111[v4123].temperature
    else
    end
  end
  return v4122
end

__demand(__parallel)
task CalculateMaxTemperature(v4130 : region(ispace(int3d), Fluid_columns), v4133 : int32, v4134 : int32, v4135 : int32, v4136 : int32, v4137 : int32, v4138 : int32) : double
where
  reads(v4130.temperature)
do
  var v4141 = -math.huge
  for v4142 in v4130 do
    if (not ((((((max(int32((uint64(v4133)-int3d(v4142).x)), int32(0))>int32(0)) or (max(int32((int3d(v4142).x-uint64(((v4134+v4133)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v4135)-int3d(v4142).y)), int32(0))>int32(0))) or (max(int32((int3d(v4142).y-uint64(((v4136+v4135)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v4137)-int3d(v4142).z)), int32(0))>int32(0))) or (max(int32((int3d(v4142).z-uint64(((v4138+v4137)-int32(1))))), int32(0))>int32(0)))) then
      v4141 max= v4130[v4142].temperature
    else
    end
  end
  return v4141
end

__demand(__parallel)
task Particles_IntegrateQuantities(v4149 : region(ispace(int1d), particles_columns)) : double
where
  reads(v4149.temperature), reads(v4149.__valid)
do
  var v4154 = double(0)
  for v4155 in v4149 do
    if v4149[v4155].__valid then
      v4154 += v4149[v4155].temperature
    else
    end
  end
  return v4154
end


__demand(__parallel)
task Radiation_InitializeCell(vRadiation : region(ispace(int3d), Radiation_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads(vRadiation.G), writes(vRadiation.G), reads(vRadiation.I_1), writes(vRadiation.I_1), reads(vRadiation.I_2), writes(vRadiation.I_2), reads(vRadiation.I_3), writes(vRadiation.I_3), reads(vRadiation.I_4), writes(vRadiation.I_4), reads(vRadiation.I_5), writes(vRadiation.I_5), reads(vRadiation.I_6), writes(vRadiation.I_6), reads(vRadiation.I_7), writes(vRadiation.I_7), reads(vRadiation.I_8), writes(vRadiation.I_8), reads(vRadiation.Iiter_1), writes(vRadiation.Iiter_1), reads(vRadiation.Iiter_2), writes(vRadiation.Iiter_2), reads(vRadiation.Iiter_3), writes(vRadiation.Iiter_3), reads(vRadiation.Iiter_4), writes(vRadiation.Iiter_4), reads(vRadiation.Iiter_5), writes(vRadiation.Iiter_5), reads(vRadiation.Iiter_6), writes(vRadiation.Iiter_6), reads(vRadiation.Iiter_7), writes(vRadiation.Iiter_7), reads(vRadiation.Iiter_8), writes(vRadiation.Iiter_8), reads(vRadiation.S), writes(vRadiation.S)
do
  for vc in vRadiation do
    for vm : int32 = int32(0), int32(NUM_ANGLES) do
      vRadiation[vc].I_1[vm] = double(0)
      vRadiation[vc].I_2[vm] = double(0)
      vRadiation[vc].I_3[vm] = double(0)
      vRadiation[vc].I_4[vm] = double(0)
      vRadiation[vc].I_5[vm] = double(0)
      vRadiation[vc].I_6[vm] = double(0)
      vRadiation[vc].I_7[vm] = double(0)
      vRadiation[vc].I_8[vm] = double(0)
      vRadiation[vc].Iiter_1[vm] = double(0)
      vRadiation[vc].Iiter_2[vm] = double(0)
      vRadiation[vc].Iiter_3[vm] = double(0)
      vRadiation[vc].Iiter_4[vm] = double(0)
      vRadiation[vc].Iiter_5[vm] = double(0)
      vRadiation[vc].Iiter_6[vm] = double(0)
      vRadiation[vc].Iiter_7[vm] = double(0)
      vRadiation[vc].Iiter_8[vm] = double(0)
    end
    vRadiation[vc].G = double(0)
    vRadiation[vc].S = double(0)
  end
end

task print__________________(v4163 : double)
  C.printf("\n Current time step: %2.6e s.\n", v4163)
end

task print_(v4166 : double, v4167 : double)
  C.printf(" Min Flow Temp: %11.6f K. Max Flow Temp: %11.6f K.\n", v4166, v4167)
end

task print__(v4170 : int64)
  C.printf(" Current number of particles: %d.\n", v4170)
end

task print___()
  C.printf("\n")
end

task print____()
  C.printf("    Iter     Time(s)   Avg Press    Avg Temp      Avg KE  Particle T\n")
end

task print_____(v4178 : int32, v4179 : double, v4180 : double, v4181 : double, v4182 : double, v4183 : double)
  C.printf("%8d %11.6f %11.6f %11.6f %11.6f %11.6f\n", v4178, v4179, v4180, v4181, v4182, v4183)
end

task GetSoundSpeed(v4202 : double, v4203 : double, v4204 : double) : double
  return [regentlib.sqrt(double)](((v4203*v4204)*v4202))
end

__demand(__parallel)
task CalculateConvectiveSpectralRadius(v4192 : region(ispace(int3d), Fluid_columns), v4194 : double, v4195 : double, v4196 : double, v4197 : double, v4198 : double, v4199 : double) : double
where
  reads(v4192.convectiveSpectralRadius), writes(v4192.convectiveSpectralRadius), reads(v4192.temperature), reads(v4192.velocity)
do
  var v4208 = -math.huge
  for v4209 in v4192 do
    var v4223 : double
    do
      var v4220 = v4192[v4209].temperature
      var v4221 = v4194
      var v4222 = v4195
      v4223 = [regentlib.sqrt(double)](((v4221*v4222)*v4220))
    end
    var v4224 = 0
    v4192[v4209].convectiveSpectralRadius = (((([regentlib.fabs(double)](v4192[v4209].velocity[int32(0)])/v4197)+([regentlib.fabs(double)](v4192[v4209].velocity[int32(1)])/v4198))+([regentlib.fabs(double)](v4192[v4209].velocity[int32(2)])/v4199))+(v4223*[regentlib.sqrt(double)](v4196)))
    v4208 max= v4192[v4209].convectiveSpectralRadius
  end
  return v4208
end

__demand(__parallel)
task CalculateViscousSpectralRadius(v4229 : region(ispace(int3d), Fluid_columns), v4231 : double, v4232 : double, v4233 : double, v4234 : double, v4235 : double, v4236 : double, v4237 : int32, v4238 : double) : double
where
  reads(v4229.rho), reads(v4229.sgsEddyViscosity), reads(v4229.temperature), reads(v4229.viscousSpectralRadius), writes(v4229.viscousSpectralRadius)
do
  var v4275 = -math.huge
  for v4276 in v4229 do
    var v4292 : double
    do
      var v4284 = v4229[v4276].temperature
      var v4285 = v4231
      var v4286 = v4232
      var v4287 = v4233
      var v4288 = v4234
      var v4289 = v4235
      var v4290 = v4236
      var v4291 = v4237
      var v4253 = double(double(0))
      if (v4291==int32(0)) then
        v4253 = v4285
      else
        if (v4291==int32(1)) then
          v4253 = (v4287*C.pow((v4284/v4286), double(0.75)))
        else
          if (v4291==int32(2)) then
            v4253 = ((v4290*C.pow((v4284/v4289), (double(3)/double(2))))*((v4289+v4288)/(v4284+v4288)))
          else
            regentlib.assert(false, "(Liszt assertion)")
          end
        end
      end
      v4292 = v4253
    end
    var v4293 = 0
    var v4277 = v4292
    var v4278 = v4229[v4276].sgsEddyViscosity
    v4229[v4276].viscousSpectralRadius = ((((double(2)*(v4277+v4278))/v4229[v4276].rho)*v4238)*double(4))
    v4275 max= v4229[v4276].viscousSpectralRadius
  end
  return v4275
end

__demand(__parallel)
task CalculateHeatConductionSpectralRadius(v4298 : region(ispace(int3d), Fluid_columns), v4300 : double, v4301 : double, v4302 : double, v4303 : double, v4304 : double, v4305 : double, v4306 : double, v4307 : double, v4308 : double, v4309 : int32, v4310 : double) : double
where
  reads(v4298.heatConductionSpectralRadius), writes(v4298.heatConductionSpectralRadius), reads(v4298.rho), reads(v4298.sgsEddyKappa), reads(v4298.temperature)
do
  var v4330 = -math.huge
  for v4331 in v4298 do
    var v4351 : double
    do
      var v4343 = v4298[v4331].temperature
      var v4344 = v4300
      var v4345 = v4303
      var v4346 = v4304
      var v4347 = v4306
      var v4348 = v4307
      var v4349 = v4308
      var v4350 = v4309
      var v4253 = double(double(0))
      if (v4350==int32(0)) then
        v4253 = v4344
      else
        if (v4350==int32(1)) then
          v4253 = (v4346*C.pow((v4343/v4345), double(0.75)))
        else
          if (v4350==int32(2)) then
            v4253 = ((v4349*C.pow((v4343/v4348), (double(3)/double(2))))*((v4348+v4347)/(v4343+v4347)))
          else
            regentlib.assert(false, "(Liszt assertion)")
          end
        end
      end
      v4351 = v4253
    end
    var v4352 = 0
    var v4332 = v4351
    var v4333 = (v4302/(v4301-double(1)))
    var v4334 = (v4301*v4333)
    var v4335 = ((v4334/v4305)*v4332)
    v4298[v4331].heatConductionSpectralRadius = ((((v4335+v4298[v4331].sgsEddyKappa)/(v4333*v4298[v4331].rho))*v4310)*double(4))
    v4330 max= v4298[v4331].heatConductionSpectralRadius
  end
  return v4330
end

__demand(__parallel)
task Flow_InitializeTemporaries(v4359 : region(ispace(int3d), Fluid_columns))
where
  reads(v4359.rho), reads(v4359.rhoEnergy), reads(v4359.rhoEnergy_new), writes(v4359.rhoEnergy_new), reads(v4359.rhoEnergy_old), writes(v4359.rhoEnergy_old), reads(v4359.rhoVelocity), reads(v4359.rhoVelocity_new), writes(v4359.rhoVelocity_new), reads(v4359.rhoVelocity_old), writes(v4359.rhoVelocity_old), reads(v4359.rho_new), writes(v4359.rho_new), reads(v4359.rho_old), writes(v4359.rho_old)
do
  for v4362 in v4359 do
    v4359[v4362].rho_old = v4359[v4362].rho
    v4359[v4362].rhoVelocity_old = v4359[v4362].rhoVelocity
    v4359[v4362].rhoEnergy_old = v4359[v4362].rhoEnergy
    v4359[v4362].rho_new = v4359[v4362].rho
    v4359[v4362].rhoVelocity_new = v4359[v4362].rhoVelocity
    v4359[v4362].rhoEnergy_new = v4359[v4362].rhoEnergy
  end
end

__demand(__parallel)
task Particles_InitializeTemporaries(v4368 : region(ispace(int1d), particles_columns))
where
  reads(v4368.position), reads(v4368.position_new), writes(v4368.position_new), reads(v4368.position_old), writes(v4368.position_old), reads(v4368.temperature), reads(v4368.temperature_new), writes(v4368.temperature_new), reads(v4368.temperature_old), writes(v4368.temperature_old), reads(v4368.velocity), reads(v4368.velocity_new), writes(v4368.velocity_new), reads(v4368.velocity_old), writes(v4368.velocity_old), reads(v4368.__valid)
do
  for v4371 in v4368 do
    if v4368[v4371].__valid then
      v4368[v4371].position_old = v4368[v4371].position
      v4368[v4371].velocity_old = v4368[v4371].velocity
      v4368[v4371].temperature_old = v4368[v4371].temperature
      v4368[v4371].position_new = v4368[v4371].position
      v4368[v4371].velocity_new = v4368[v4371].velocity
      v4368[v4371].temperature_new = v4368[v4371].temperature
    else
    end
  end
end

__demand(__parallel)
task Flow_InitializeTimeDerivatives(v4377 : region(ispace(int3d), Fluid_columns))
where
  reads(v4377.pressure), reads(v4377.rhoEnergy), reads(v4377.rhoEnergy_t), writes(v4377.rhoEnergy_t), reads(v4377.rhoEnthalpy), writes(v4377.rhoEnthalpy), reads(v4377.rhoVelocity_t), writes(v4377.rhoVelocity_t), reads(v4377.rho_t), writes(v4377.rho_t)
do
  for v4380 in v4377 do
    v4377[v4380].rho_t = double(double(0))
    v4377[v4380].rhoVelocity_t = [double[3]](array(double(0), double(0), double(0)))
    v4377[v4380].rhoEnergy_t = double(double(0))
    v4377[v4380].rhoEnthalpy = (v4377[v4380].rhoEnergy+v4377[v4380].pressure)
  end
end

__demand(__parallel)
task Particles_InitializeTimeDerivatives(v4392 : region(ispace(int1d), particles_columns))
where
  reads(v4392.position_t), writes(v4392.position_t), reads(v4392.temperature_t), writes(v4392.temperature_t), reads(v4392.velocity_t), writes(v4392.velocity_t), reads(v4392.__valid)
do
  for v4395 in v4392 do
    if v4392[v4395].__valid then
      v4392[v4395].position_t = [double[3]](array(double(0), double(0), double(0)))
      v4392[v4395].velocity_t = [double[3]](array(double(0), double(0), double(0)))
      v4392[v4395].temperature_t = double(int32(0))
    else
    end
  end
end

__demand(__parallel)
task Flow_UpdateGhostVelocityGradientStep1(v4413 : region(ispace(int3d), Fluid_columns), v4415 : double[3], v4416 : double[3], v4417 : double[3], v4418 : int32, v4419 : int32, v4420 : int32, v4421 : int32, v4422 : int32, v4423 : int32)
where
  reads(v4413.velocityGradientX), reads(v4413.velocityGradientXBoundary), writes(v4413.velocityGradientXBoundary), reads(v4413.velocityGradientY), reads(v4413.velocityGradientYBoundary), writes(v4413.velocityGradientYBoundary), reads(v4413.velocityGradientZ), reads(v4413.velocityGradientZBoundary), writes(v4413.velocityGradientZBoundary)
do
  for v4551 in v4413 do
    if (max(int32((uint64(v4418)-int3d(v4551).x)), int32(0))>int32(0)) then
      var v4552 = int3d(v4551)
      var v4553 = ((v4551+{1, 0, 0})%v4413.bounds)
      var v4554 = v4415
      v4413[v4552].velocityGradientXBoundary = vv_mul_double_3(v4554, v4413[v4553].velocityGradientX)
      v4413[v4552].velocityGradientYBoundary = vv_mul_double_3(v4554, v4413[v4553].velocityGradientY)
      v4413[v4552].velocityGradientZBoundary = vv_mul_double_3(v4554, v4413[v4553].velocityGradientZ)
    else
    end
    if (max(int32((int3d(v4551).x-uint64(((v4419+v4418)-int32(1))))), int32(0))>int32(0)) then
      var v4555 = int3d(v4551)
      var v4556 = ((v4551+{-1, 0, 0})%v4413.bounds)
      var v4557 = v4415
      v4413[v4555].velocityGradientXBoundary = vv_mul_double_3(v4557, v4413[v4556].velocityGradientX)
      v4413[v4555].velocityGradientYBoundary = vv_mul_double_3(v4557, v4413[v4556].velocityGradientY)
      v4413[v4555].velocityGradientZBoundary = vv_mul_double_3(v4557, v4413[v4556].velocityGradientZ)
    else
    end
    if (max(int32((uint64(v4420)-int3d(v4551).y)), int32(0))>int32(0)) then
      var v4558 = int3d(v4551)
      var v4559 = ((v4551+{0, 1, 0})%v4413.bounds)
      var v4560 = v4416
      v4413[v4558].velocityGradientXBoundary = vv_mul_double_3(v4560, v4413[v4559].velocityGradientX)
      v4413[v4558].velocityGradientYBoundary = vv_mul_double_3(v4560, v4413[v4559].velocityGradientY)
      v4413[v4558].velocityGradientZBoundary = vv_mul_double_3(v4560, v4413[v4559].velocityGradientZ)
    else
    end
    if (max(int32((int3d(v4551).y-uint64(((v4421+v4420)-int32(1))))), int32(0))>int32(0)) then
      var v4561 = int3d(v4551)
      var v4562 = ((v4551+{0, -1, 0})%v4413.bounds)
      var v4563 = v4416
      v4413[v4561].velocityGradientXBoundary = vv_mul_double_3(v4563, v4413[v4562].velocityGradientX)
      v4413[v4561].velocityGradientYBoundary = vv_mul_double_3(v4563, v4413[v4562].velocityGradientY)
      v4413[v4561].velocityGradientZBoundary = vv_mul_double_3(v4563, v4413[v4562].velocityGradientZ)
    else
    end
    if (max(int32((uint64(v4422)-int3d(v4551).z)), int32(0))>int32(0)) then
      var v4564 = int3d(v4551)
      var v4565 = ((v4551+{0, 0, 1})%v4413.bounds)
      var v4566 = v4417
      v4413[v4564].velocityGradientXBoundary = vv_mul_double_3(v4566, v4413[v4565].velocityGradientX)
      v4413[v4564].velocityGradientYBoundary = vv_mul_double_3(v4566, v4413[v4565].velocityGradientY)
      v4413[v4564].velocityGradientZBoundary = vv_mul_double_3(v4566, v4413[v4565].velocityGradientZ)
    else
    end
    if (max(int32((int3d(v4551).z-uint64(((v4423+v4422)-int32(1))))), int32(0))>int32(0)) then
      var v4567 = int3d(v4551)
      var v4568 = ((v4551+{0, 0, -1})%v4413.bounds)
      var v4569 = v4417
      v4413[v4567].velocityGradientXBoundary = vv_mul_double_3(v4569, v4413[v4568].velocityGradientX)
      v4413[v4567].velocityGradientYBoundary = vv_mul_double_3(v4569, v4413[v4568].velocityGradientY)
      v4413[v4567].velocityGradientZBoundary = vv_mul_double_3(v4569, v4413[v4568].velocityGradientZ)
    else
    end
  end
end

__demand(__parallel)
task Flow_UpdateGhostVelocityGradientStep2(v4695 : region(ispace(int3d), Fluid_columns), v4697 : int32, v4698 : int32, v4699 : int32, v4700 : int32, v4701 : int32, v4702 : int32)
where
  reads(v4695.velocityGradientX), writes(v4695.velocityGradientX), reads(v4695.velocityGradientXBoundary), reads(v4695.velocityGradientY), writes(v4695.velocityGradientY), reads(v4695.velocityGradientYBoundary), reads(v4695.velocityGradientZ), writes(v4695.velocityGradientZ), reads(v4695.velocityGradientZBoundary)
do
  for v4704 in v4695 do
    if ((((((max(int32((uint64(v4697)-int3d(v4704).x)), int32(0))>int32(0)) or (max(int32((int3d(v4704).x-uint64(((v4698+v4697)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v4699)-int3d(v4704).y)), int32(0))>int32(0))) or (max(int32((int3d(v4704).y-uint64(((v4700+v4699)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v4701)-int3d(v4704).z)), int32(0))>int32(0))) or (max(int32((int3d(v4704).z-uint64(((v4702+v4701)-int32(1))))), int32(0))>int32(0))) then
      v4695[v4704].velocityGradientX = v4695[v4704].velocityGradientXBoundary
      v4695[v4704].velocityGradientY = v4695[v4704].velocityGradientYBoundary
      v4695[v4704].velocityGradientZ = v4695[v4704].velocityGradientZBoundary
    else
    end
  end
end

__demand(__parallel)
task Flow_AddGetFlux(v4710 : region(ispace(int3d), Fluid_columns), v4712 : double, v4713 : double, v4714 : double, v4715 : double, v4716 : double, v4717 : double, v4718 : double, v4719 : double, v4720 : double, v4721 : int32, v4722 : int32, v4723 : double, v4724 : int32, v4725 : int32, v4726 : double, v4727 : int32, v4728 : int32, v4729 : double, v4730 : int32)
where
  reads(v4710.pressure), reads(v4710.rho), reads(v4710.rhoEnergyFluxX), writes(v4710.rhoEnergyFluxX), reads(v4710.rhoEnergyFluxX), writes(v4710.rhoEnergyFluxX), reads(v4710.rhoEnergyFluxY), writes(v4710.rhoEnergyFluxY), reads(v4710.rhoEnergyFluxY), writes(v4710.rhoEnergyFluxY), reads(v4710.rhoEnergyFluxZ), writes(v4710.rhoEnergyFluxZ), reads(v4710.rhoEnergyFluxZ), writes(v4710.rhoEnergyFluxZ), reads(v4710.rhoEnthalpy), reads(v4710.rhoFluxX), writes(v4710.rhoFluxX), reads(v4710.rhoFluxY), writes(v4710.rhoFluxY), reads(v4710.rhoFluxZ), writes(v4710.rhoFluxZ), reads(v4710.rhoVelocity), reads(v4710.rhoVelocityFluxX), writes(v4710.rhoVelocityFluxX), reads(v4710.rhoVelocityFluxX), writes(v4710.rhoVelocityFluxX), reads(v4710.rhoVelocityFluxY), writes(v4710.rhoVelocityFluxY), reads(v4710.rhoVelocityFluxY), writes(v4710.rhoVelocityFluxY), reads(v4710.rhoVelocityFluxZ), writes(v4710.rhoVelocityFluxZ), reads(v4710.rhoVelocityFluxZ), writes(v4710.rhoVelocityFluxZ), reads(v4710.temperature), reads(v4710.velocity), reads(v4710.velocityGradientX), reads(v4710.velocityGradientY), reads(v4710.velocityGradientZ)
do
  for v5477 in v4710 do
    if ((not ((((((max(int32((uint64(v4722)-int3d(v5477).x)), int32(0))>int32(0)) or (max(int32((int3d(v5477).x-uint64(((v4724+v4722)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v4725)-int3d(v5477).y)), int32(0))>int32(0))) or (max(int32((int3d(v5477).y-uint64(((v4727+v4725)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v4728)-int3d(v5477).z)), int32(0))>int32(0))) or (max(int32((int3d(v5477).z-uint64(((v4730+v4728)-int32(1))))), int32(0))>int32(0)))) or (max(int32((uint64(v4722)-int3d(v5477).x)), int32(0))==int32(1))) then
      var v5660 : double[5]
      do
        var v5658 = int3d(v5477)
        var v5659 = ((v5477+{1, 0, 0})%v4710.bounds)
        var v4788 = double(double(0))
        var v4789 = [double[3]](array(double(0), double(0), double(0)))
        var v4790 = double(double(0))
        var v4791 = double(double(0))
        v4788 = (double(0.5)*((v4710[v5658].rho*v4710[v5658].velocity[int32(0)])+(v4710[v5659].rho*v4710[v5659].velocity[int32(0)])))
        v4789 = vs_mul_double_3(vv_add_double_3(vs_mul_double_3(v4710[v5658].rhoVelocity, v4710[v5658].velocity[int32(0)]), vs_mul_double_3(v4710[v5659].rhoVelocity, v4710[v5659].velocity[int32(0)])), double(0.5))
        v4790 = (double(0.5)*((v4710[v5658].rhoEnthalpy*v4710[v5658].velocity[int32(0)])+(v4710[v5659].rhoEnthalpy*v4710[v5659].velocity[int32(0)])))
        v4791 += (double(0.5)*(v4710[v5658].pressure+v4710[v5659].pressure))
        var v4792 = double(double(0))
        var v4793 = [double[3]](array(double(0), double(0), double(0)))
        var v4794 = double(double(0))
        var v4795 = double(double(0))
        v4795 = (double(0.5)*v4710[v5659].velocity[int32(0)])
        v4792 += (v4710[v5658].rho*v4795)
        var v4796 = vs_mul_double_3(v4710[v5658].rhoVelocity, v4795)
        var v4797 = v4793
        v4797[0] += v4796[0]
        v4797[1] += v4796[1]
        v4797[2] += v4796[2]
        v4793 = v4797
        v4794 += (v4710[v5658].rhoEnthalpy*v4795)
        v4795 = (double(0.5)*v4710[v5658].velocity[int32(0)])
        v4792 += (v4710[v5659].rho*v4795)
        var v4798 = vs_mul_double_3(v4710[v5659].rhoVelocity, v4795)
        var v4799 = v4793
        v4799[0] += v4798[0]
        v4799[1] += v4798[1]
        v4799[2] += v4798[2]
        v4793 = v4799
        v4794 += (v4710[v5659].rhoEnthalpy*v4795)
        var v4800 = double(0.5)
        var v4801 = ((v4800*v4788)+((double(int32(1))-v4800)*v4792))
        var v4802 = vv_add_double_3(vs_mul_double_3(v4789, v4800), vs_mul_double_3(v4793, (double(int32(1))-v4800)))
        var v4803 = ((v4800*v4790)+((double(int32(1))-v4800)*v4794))
        v4802[int32(0)] += v4791
        v5660 = array(v4801, v4802[int32(0)], v4802[int32(1)], v4802[int32(2)], v4803)
      end
      var v5661 = 0
      var v5478 = v5660
      v4710[v5477].rhoFluxX = v5478[int32(0)]
      v4710[v5477].rhoVelocityFluxX = array(v5478[int32(1)], v5478[int32(2)], v5478[int32(3)])
      v4710[v5477].rhoEnergyFluxX = v5478[int32(4)]
      var v5670 : double
      do
        var v5662 = v4710[((v5477+{1, 0, 0})%v4710.bounds)].temperature
        var v5663 = v4712
        var v5664 = v4715
        var v5665 = v4716
        var v5666 = v4718
        var v5667 = v4719
        var v5668 = v4720
        var v5669 = v4721
        var v4253 = double(double(0))
        if (v5669==int32(0)) then
          v4253 = v5663
        else
          if (v5669==int32(1)) then
            v4253 = (v5665*C.pow((v5662/v5664), double(0.75)))
          else
            if (v5669==int32(2)) then
              v4253 = ((v5668*C.pow((v5662/v5667), (double(3)/double(2))))*((v5667+v5666)/(v5662+v5666)))
            else
              regentlib.assert(false, "(Liszt assertion)")
            end
          end
        end
        v5670 = v4253
      end
      var v5671 = 0
      var v5680 : double
      do
        var v5672 = v4710[v5477].temperature
        var v5673 = v4712
        var v5674 = v4715
        var v5675 = v4716
        var v5676 = v4718
        var v5677 = v4719
        var v5678 = v4720
        var v5679 = v4721
        var v4253 = double(double(0))
        if (v5679==int32(0)) then
          v4253 = v5673
        else
          if (v5679==int32(1)) then
            v4253 = (v5675*C.pow((v5672/v5674), double(0.75)))
          else
            if (v5679==int32(2)) then
              v4253 = ((v5678*C.pow((v5672/v5677), (double(3)/double(2))))*((v5677+v5676)/(v5672+v5676)))
            else
              regentlib.assert(false, "(Liszt assertion)")
            end
          end
        end
        v5680 = v4253
      end
      var v5681 = 0
      var v5479 = (double(0.5)*(v5680+v5670))
      var v5480 = [double[3]](array(double(0), double(0), double(0)))
      var v5481 = double(double(0))
      var v5482 = double(double(0))
      var v5483 = double(double(0))
      var v5484 = double(double(0))
      v5480 = vs_mul_double_3(vv_add_double_3(v4710[v5477].velocity, v4710[((v5477+{1, 0, 0})%v4710.bounds)].velocity), double(0.5))
      v5481 = (double(0.5)*(v4710[v5477].velocityGradientY[int32(0)]+v4710[((v5477+{1, 0, 0})%v4710.bounds)].velocityGradientY[int32(0)]))
      v5482 = (double(0.5)*(v4710[v5477].velocityGradientZ[int32(0)]+v4710[((v5477+{1, 0, 0})%v4710.bounds)].velocityGradientZ[int32(0)]))
      v5483 = (double(0.5)*(v4710[v5477].velocityGradientY[int32(1)]+v4710[((v5477+{1, 0, 0})%v4710.bounds)].velocityGradientY[int32(1)]))
      v5484 = (double(0.5)*(v4710[v5477].velocityGradientZ[int32(2)]+v4710[((v5477+{1, 0, 0})%v4710.bounds)].velocityGradientZ[int32(2)]))
      var v5485 = double(double(0))
      var v5486 = double(double(0))
      var v5487 = double(double(0))
      var v5488 = double(double(0))
      v5485 = (double(0.5)*(v4710[((v5477+{1, 0, 0})%v4710.bounds)].velocity[int32(0)]-v4710[v5477].velocity[int32(0)]))
      v5486 = (double(0.5)*(v4710[((v5477+{1, 0, 0})%v4710.bounds)].velocity[int32(1)]-v4710[v5477].velocity[int32(1)]))
      v5487 = (double(0.5)*(v4710[((v5477+{1, 0, 0})%v4710.bounds)].velocity[int32(2)]-v4710[v5477].velocity[int32(2)]))
      v5488 = (double(0.5)*(v4710[((v5477+{1, 0, 0})%v4710.bounds)].temperature-v4710[v5477].temperature))
      v5485 *= (1/(v4723*double(0.5)))
      v5486 *= (1/(v4723*double(0.5)))
      v5487 *= (1/(v4723*double(0.5)))
      v5488 *= (1/(v4723*double(0.5)))
      var v5489 = ((v5479*(((double(4)*v5485)-(double(2)*v5483))-(double(2)*v5484)))/double(3))
      var v5490 = (v5479*(v5486+v5481))
      var v5491 = (v5479*(v5487+v5482))
      var v5492 = (((v5480[int32(0)]*v5489)+(v5480[int32(1)]*v5490))+(v5480[int32(2)]*v5491))
      var v5493 = ((v4713*v4714)/(v4713-double(1)))
      var v5494 = ((-((v5493*v5479)/v4717))*v5488)
      v4710[v5477].rhoVelocityFluxX[int32(0)] += (-v5489)
      v4710[v5477].rhoVelocityFluxX[int32(1)] += (-v5490)
      v4710[v5477].rhoVelocityFluxX[int32(2)] += (-v5491)
      v4710[v5477].rhoEnergyFluxX += (-(v5492-v5494))
    else
    end
    if ((not ((((((max(int32((uint64(v4722)-int3d(v5477).x)), int32(0))>int32(0)) or (max(int32((int3d(v5477).x-uint64(((v4724+v4722)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v4725)-int3d(v5477).y)), int32(0))>int32(0))) or (max(int32((int3d(v5477).y-uint64(((v4727+v4725)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v4728)-int3d(v5477).z)), int32(0))>int32(0))) or (max(int32((int3d(v5477).z-uint64(((v4730+v4728)-int32(1))))), int32(0))>int32(0)))) or (max(int32((uint64(v4725)-int3d(v5477).y)), int32(0))==int32(1))) then
      var v5684 : double[5]
      do
        var v5682 = int3d(v5477)
        var v5683 = ((v5477+{0, 1, 0})%v4710.bounds)
        var v4995 = double(double(0))
        var v4996 = [double[3]](array(double(0), double(0), double(0)))
        var v4997 = double(double(0))
        var v4998 = double(double(0))
        v4995 = (double(0.5)*((v4710[v5682].rho*v4710[v5682].velocity[int32(1)])+(v4710[v5683].rho*v4710[v5683].velocity[int32(1)])))
        v4996 = vs_mul_double_3(vv_add_double_3(vs_mul_double_3(v4710[v5682].rhoVelocity, v4710[v5682].velocity[int32(1)]), vs_mul_double_3(v4710[v5683].rhoVelocity, v4710[v5683].velocity[int32(1)])), double(0.5))
        v4997 = (double(0.5)*((v4710[v5682].rhoEnthalpy*v4710[v5682].velocity[int32(1)])+(v4710[v5683].rhoEnthalpy*v4710[v5683].velocity[int32(1)])))
        v4998 += (double(0.5)*(v4710[v5682].pressure+v4710[v5683].pressure))
        var v4999 = double(double(0))
        var v5000 = [double[3]](array(double(0), double(0), double(0)))
        var v5001 = double(double(0))
        var v5002 = double(double(0))
        v5002 = (double(0.5)*v4710[v5683].velocity[int32(1)])
        v4999 += (v4710[v5682].rho*v5002)
        var v5003 = vs_mul_double_3(v4710[v5682].rhoVelocity, v5002)
        var v5004 = v5000
        v5004[0] += v5003[0]
        v5004[1] += v5003[1]
        v5004[2] += v5003[2]
        v5000 = v5004
        v5001 += (v4710[v5682].rhoEnthalpy*v5002)
        v5002 = (double(0.5)*v4710[v5682].velocity[int32(1)])
        v4999 += (v4710[v5683].rho*v5002)
        var v5005 = vs_mul_double_3(v4710[v5683].rhoVelocity, v5002)
        var v5006 = v5000
        v5006[0] += v5005[0]
        v5006[1] += v5005[1]
        v5006[2] += v5005[2]
        v5000 = v5006
        v5001 += (v4710[v5683].rhoEnthalpy*v5002)
        var v5007 = double(0.5)
        var v5008 = ((v5007*v4995)+((double(int32(1))-v5007)*v4999))
        var v5009 = vv_add_double_3(vs_mul_double_3(v4996, v5007), vs_mul_double_3(v5000, (double(int32(1))-v5007)))
        var v5010 = ((v5007*v4997)+((double(int32(1))-v5007)*v5001))
        v5009[int32(1)] += v4998
        v5684 = array(v5008, v5009[int32(0)], v5009[int32(1)], v5009[int32(2)], v5010)
      end
      var v5685 = 0
      var v5495 = v5684
      v4710[v5477].rhoFluxY = v5495[int32(0)]
      v4710[v5477].rhoVelocityFluxY = array(v5495[int32(1)], v5495[int32(2)], v5495[int32(3)])
      v4710[v5477].rhoEnergyFluxY = v5495[int32(4)]
      var v5694 : double
      do
        var v5686 = v4710[((v5477+{0, 1, 0})%v4710.bounds)].temperature
        var v5687 = v4712
        var v5688 = v4715
        var v5689 = v4716
        var v5690 = v4718
        var v5691 = v4719
        var v5692 = v4720
        var v5693 = v4721
        var v4253 = double(double(0))
        if (v5693==int32(0)) then
          v4253 = v5687
        else
          if (v5693==int32(1)) then
            v4253 = (v5689*C.pow((v5686/v5688), double(0.75)))
          else
            if (v5693==int32(2)) then
              v4253 = ((v5692*C.pow((v5686/v5691), (double(3)/double(2))))*((v5691+v5690)/(v5686+v5690)))
            else
              regentlib.assert(false, "(Liszt assertion)")
            end
          end
        end
        v5694 = v4253
      end
      var v5695 = 0
      var v5704 : double
      do
        var v5696 = v4710[v5477].temperature
        var v5697 = v4712
        var v5698 = v4715
        var v5699 = v4716
        var v5700 = v4718
        var v5701 = v4719
        var v5702 = v4720
        var v5703 = v4721
        var v4253 = double(double(0))
        if (v5703==int32(0)) then
          v4253 = v5697
        else
          if (v5703==int32(1)) then
            v4253 = (v5699*C.pow((v5696/v5698), double(0.75)))
          else
            if (v5703==int32(2)) then
              v4253 = ((v5702*C.pow((v5696/v5701), (double(3)/double(2))))*((v5701+v5700)/(v5696+v5700)))
            else
              regentlib.assert(false, "(Liszt assertion)")
            end
          end
        end
        v5704 = v4253
      end
      var v5705 = 0
      var v5496 = (double(0.5)*(v5704+v5694))
      var v5497 = [double[3]](array(double(0), double(0), double(0)))
      var v5498 = double(double(0))
      var v5499 = double(double(0))
      var v5500 = double(double(0))
      var v5501 = double(double(0))
      v5497 = vs_mul_double_3(vv_add_double_3(v4710[v5477].velocity, v4710[((v5477+{0, 1, 0})%v4710.bounds)].velocity), double(0.5))
      v5498 = (double(0.5)*(v4710[v5477].velocityGradientX[int32(1)]+v4710[((v5477+{0, 1, 0})%v4710.bounds)].velocityGradientX[int32(1)]))
      v5499 = (double(0.5)*(v4710[v5477].velocityGradientZ[int32(1)]+v4710[((v5477+{0, 1, 0})%v4710.bounds)].velocityGradientZ[int32(1)]))
      v5500 = (double(0.5)*(v4710[v5477].velocityGradientX[int32(0)]+v4710[((v5477+{0, 1, 0})%v4710.bounds)].velocityGradientX[int32(0)]))
      v5501 = (double(0.5)*(v4710[v5477].velocityGradientZ[int32(2)]+v4710[((v5477+{0, 1, 0})%v4710.bounds)].velocityGradientZ[int32(2)]))
      var v5502 = double(double(0))
      var v5503 = double(double(0))
      var v5504 = double(double(0))
      var v5505 = double(double(0))
      v5502 = (double(0.5)*(v4710[((v5477+{0, 1, 0})%v4710.bounds)].velocity[int32(0)]-v4710[v5477].velocity[int32(0)]))
      v5503 = (double(0.5)*(v4710[((v5477+{0, 1, 0})%v4710.bounds)].velocity[int32(1)]-v4710[v5477].velocity[int32(1)]))
      v5504 = (double(0.5)*(v4710[((v5477+{0, 1, 0})%v4710.bounds)].velocity[int32(2)]-v4710[v5477].velocity[int32(2)]))
      v5505 = (double(0.5)*(v4710[((v5477+{0, 1, 0})%v4710.bounds)].temperature-v4710[v5477].temperature))
      v5502 *= (1/(v4726*double(0.5)))
      v5503 *= (1/(v4726*double(0.5)))
      v5504 *= (1/(v4726*double(0.5)))
      v5505 *= (1/(v4726*double(0.5)))
      var v5506 = (v5496*(v5502+v5498))
      var v5507 = ((v5496*(((double(4)*v5503)-(double(2)*v5500))-(double(2)*v5501)))/double(3))
      var v5508 = (v5496*(v5504+v5499))
      var v5509 = (((v5497[int32(0)]*v5506)+(v5497[int32(1)]*v5507))+(v5497[int32(2)]*v5508))
      var v5510 = ((v4713*v4714)/(v4713-double(1)))
      var v5511 = ((-((v5510*v5496)/v4717))*v5505)
      v4710[v5477].rhoVelocityFluxY[int32(0)] += (-v5506)
      v4710[v5477].rhoVelocityFluxY[int32(1)] += (-v5507)
      v4710[v5477].rhoVelocityFluxY[int32(2)] += (-v5508)
      v4710[v5477].rhoEnergyFluxY += (-(v5509-v5511))
    else
    end
    if ((not ((((((max(int32((uint64(v4722)-int3d(v5477).x)), int32(0))>int32(0)) or (max(int32((int3d(v5477).x-uint64(((v4724+v4722)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v4725)-int3d(v5477).y)), int32(0))>int32(0))) or (max(int32((int3d(v5477).y-uint64(((v4727+v4725)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v4728)-int3d(v5477).z)), int32(0))>int32(0))) or (max(int32((int3d(v5477).z-uint64(((v4730+v4728)-int32(1))))), int32(0))>int32(0)))) or (max(int32((uint64(v4728)-int3d(v5477).z)), int32(0))==int32(1))) then
      var v5708 : double[5]
      do
        var v5706 = int3d(v5477)
        var v5707 = ((v5477+{0, 0, 1})%v4710.bounds)
        var v5202 = double(double(0))
        var v5203 = [double[3]](array(double(0), double(0), double(0)))
        var v5204 = double(double(0))
        var v5205 = double(double(0))
        v5202 = (double(0.5)*((v4710[v5706].rho*v4710[v5706].velocity[int32(2)])+(v4710[v5707].rho*v4710[v5707].velocity[int32(2)])))
        v5203 = vs_mul_double_3(vv_add_double_3(vs_mul_double_3(v4710[v5706].rhoVelocity, v4710[v5706].velocity[int32(2)]), vs_mul_double_3(v4710[v5707].rhoVelocity, v4710[v5707].velocity[int32(2)])), double(0.5))
        v5204 = (double(0.5)*((v4710[v5706].rhoEnthalpy*v4710[v5706].velocity[int32(2)])+(v4710[v5707].rhoEnthalpy*v4710[v5707].velocity[int32(2)])))
        v5205 += (double(0.5)*(v4710[v5706].pressure+v4710[v5707].pressure))
        var v5206 = double(double(0))
        var v5207 = [double[3]](array(double(0), double(0), double(0)))
        var v5208 = double(double(0))
        var v5209 = double(double(0))
        v5209 = (double(0.5)*v4710[v5707].velocity[int32(2)])
        v5206 += (v4710[v5706].rho*v5209)
        var v5210 = vs_mul_double_3(v4710[v5706].rhoVelocity, v5209)
        var v5211 = v5207
        v5211[0] += v5210[0]
        v5211[1] += v5210[1]
        v5211[2] += v5210[2]
        v5207 = v5211
        v5208 += (v4710[v5706].rhoEnthalpy*v5209)
        v5209 = (double(0.5)*v4710[v5706].velocity[int32(2)])
        v5206 += (v4710[v5707].rho*v5209)
        var v5212 = vs_mul_double_3(v4710[v5707].rhoVelocity, v5209)
        var v5213 = v5207
        v5213[0] += v5212[0]
        v5213[1] += v5212[1]
        v5213[2] += v5212[2]
        v5207 = v5213
        v5208 += (v4710[v5707].rhoEnthalpy*v5209)
        var v5214 = double(0.5)
        var v5215 = ((v5214*v5202)+((double(int32(1))-v5214)*v5206))
        var v5216 = vv_add_double_3(vs_mul_double_3(v5203, v5214), vs_mul_double_3(v5207, (double(int32(1))-v5214)))
        var v5217 = ((v5214*v5204)+((double(int32(1))-v5214)*v5208))
        v5216[int32(2)] += v5205
        v5708 = array(v5215, v5216[int32(0)], v5216[int32(1)], v5216[int32(2)], v5217)
      end
      var v5709 = 0
      var v5512 = v5708
      v4710[v5477].rhoFluxZ = v5512[int32(0)]
      v4710[v5477].rhoVelocityFluxZ = array(v5512[int32(1)], v5512[int32(2)], v5512[int32(3)])
      v4710[v5477].rhoEnergyFluxZ = v5512[int32(4)]
      var v5718 : double
      do
        var v5710 = v4710[((v5477+{0, 0, 1})%v4710.bounds)].temperature
        var v5711 = v4712
        var v5712 = v4715
        var v5713 = v4716
        var v5714 = v4718
        var v5715 = v4719
        var v5716 = v4720
        var v5717 = v4721
        var v4253 = double(double(0))
        if (v5717==int32(0)) then
          v4253 = v5711
        else
          if (v5717==int32(1)) then
            v4253 = (v5713*C.pow((v5710/v5712), double(0.75)))
          else
            if (v5717==int32(2)) then
              v4253 = ((v5716*C.pow((v5710/v5715), (double(3)/double(2))))*((v5715+v5714)/(v5710+v5714)))
            else
              regentlib.assert(false, "(Liszt assertion)")
            end
          end
        end
        v5718 = v4253
      end
      var v5719 = 0
      var v5728 : double
      do
        var v5720 = v4710[v5477].temperature
        var v5721 = v4712
        var v5722 = v4715
        var v5723 = v4716
        var v5724 = v4718
        var v5725 = v4719
        var v5726 = v4720
        var v5727 = v4721
        var v4253 = double(double(0))
        if (v5727==int32(0)) then
          v4253 = v5721
        else
          if (v5727==int32(1)) then
            v4253 = (v5723*C.pow((v5720/v5722), double(0.75)))
          else
            if (v5727==int32(2)) then
              v4253 = ((v5726*C.pow((v5720/v5725), (double(3)/double(2))))*((v5725+v5724)/(v5720+v5724)))
            else
              regentlib.assert(false, "(Liszt assertion)")
            end
          end
        end
        v5728 = v4253
      end
      var v5729 = 0
      var v5513 = (double(0.5)*(v5728+v5718))
      var v5514 = [double[3]](array(double(0), double(0), double(0)))
      var v5515 = double(double(0))
      var v5516 = double(double(0))
      var v5517 = double(double(0))
      var v5518 = double(double(0))
      v5514 = vs_mul_double_3(vv_add_double_3(v4710[v5477].velocity, v4710[((v5477+{0, 0, 1})%v4710.bounds)].velocity), double(0.5))
      v5515 = (double(0.5)*(v4710[v5477].velocityGradientX[int32(2)]+v4710[((v5477+{0, 0, 1})%v4710.bounds)].velocityGradientX[int32(2)]))
      v5516 = (double(0.5)*(v4710[v5477].velocityGradientY[int32(2)]+v4710[((v5477+{0, 0, 1})%v4710.bounds)].velocityGradientY[int32(2)]))
      v5517 = (double(0.5)*(v4710[v5477].velocityGradientX[int32(0)]+v4710[((v5477+{0, 0, 1})%v4710.bounds)].velocityGradientX[int32(0)]))
      v5518 = (double(0.5)*(v4710[v5477].velocityGradientY[int32(1)]+v4710[((v5477+{0, 0, 1})%v4710.bounds)].velocityGradientY[int32(1)]))
      var v5519 = double(double(0))
      var v5520 = double(double(0))
      var v5521 = double(double(0))
      var v5522 = double(double(0))
      v5519 = (double(0.5)*(v4710[((v5477+{0, 0, 1})%v4710.bounds)].velocity[int32(0)]-v4710[v5477].velocity[int32(0)]))
      v5520 = (double(0.5)*(v4710[((v5477+{0, 0, 1})%v4710.bounds)].velocity[int32(1)]-v4710[v5477].velocity[int32(1)]))
      v5521 = (double(0.5)*(v4710[((v5477+{0, 0, 1})%v4710.bounds)].velocity[int32(2)]-v4710[v5477].velocity[int32(2)]))
      v5522 = (double(0.5)*(v4710[((v5477+{0, 0, 1})%v4710.bounds)].temperature-v4710[v5477].temperature))
      v5519 *= (1/(v4729*double(0.5)))
      v5520 *= (1/(v4729*double(0.5)))
      v5521 *= (1/(v4729*double(0.5)))
      v5522 *= (1/(v4729*double(0.5)))
      var v5523 = (v5513*(v5519+v5515))
      var v5524 = (v5513*(v5520+v5516))
      var v5525 = ((v5513*(((double(4)*v5521)-(double(2)*v5517))-(double(2)*v5518)))/double(3))
      var v5526 = (((v5514[int32(0)]*v5523)+(v5514[int32(1)]*v5524))+(v5514[int32(2)]*v5525))
      var v5527 = ((v4713*v4714)/(v4713-double(1)))
      var v5528 = ((-((v5527*v5513)/v4717))*v5522)
      v4710[v5477].rhoVelocityFluxZ[int32(0)] += (-v5523)
      v4710[v5477].rhoVelocityFluxZ[int32(1)] += (-v5524)
      v4710[v5477].rhoVelocityFluxZ[int32(2)] += (-v5525)
      v4710[v5477].rhoEnergyFluxZ += (-(v5526-v5528))
    else
    end
    var v5529 = int32(0)
    if (v5529==int32(1)) then
      var v5530 = v4710[((v5477+{-1, 0, 0})%v4710.bounds)].velocity[int32(0)]
      var v5531 = v4710[((v5477+{0, -1, 0})%v4710.bounds)].velocity[int32(1)]
      var v5532 = v4710[((v5477+{0, 0, -1})%v4710.bounds)].velocity[int32(2)]
    else
    end
  end
end

__demand(__parallel)
task Flow_AddUpdateUsingFlux(v5911 : region(ispace(int3d), Fluid_columns), v5913 : int32, v5914 : double, v5915 : int32, v5916 : int32, v5917 : double, v5918 : int32, v5919 : int32, v5920 : double, v5921 : int32)
where
  reads(v5911.rhoEnergyFluxX), reads(v5911.rhoEnergyFluxY), reads(v5911.rhoEnergyFluxZ), reads(v5911.rhoEnergy_t), writes(v5911.rhoEnergy_t), reads(v5911.rhoFluxX), reads(v5911.rhoFluxY), reads(v5911.rhoFluxZ), reads(v5911.rhoVelocityFluxX), reads(v5911.rhoVelocityFluxY), reads(v5911.rhoVelocityFluxZ), reads(v5911.rhoVelocity_t), writes(v5911.rhoVelocity_t), reads(v5911.rho_t), writes(v5911.rho_t)
do
  for v5965 in v5911 do
    if (not ((((((max(int32((uint64(v5913)-int3d(v5965).x)), int32(0))>int32(0)) or (max(int32((int3d(v5965).x-uint64(((v5915+v5913)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v5916)-int3d(v5965).y)), int32(0))>int32(0))) or (max(int32((int3d(v5965).y-uint64(((v5918+v5916)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v5919)-int3d(v5965).z)), int32(0))>int32(0))) or (max(int32((int3d(v5965).z-uint64(((v5921+v5919)-int32(1))))), int32(0))>int32(0)))) then
      v5911[v5965].rho_t += ((-(v5911[v5965].rhoFluxX-v5911[((v5965+{-1, 0, 0})%v5911.bounds)].rhoFluxX))/v5914)
      var v5966 = vs_div_double_3(vs_mul_double_3(vv_sub_double_3(v5911[v5965].rhoVelocityFluxX, v5911[((v5965+{-1, 0, 0})%v5911.bounds)].rhoVelocityFluxX), double((-1))), v5914)
      var v5967 = v5911[v5965].rhoVelocity_t
      v5967[0] += v5966[0]
      v5967[1] += v5966[1]
      v5967[2] += v5966[2]
      v5911[v5965].rhoVelocity_t = v5967
      v5911[v5965].rhoEnergy_t += ((-(v5911[v5965].rhoEnergyFluxX-v5911[((v5965+{-1, 0, 0})%v5911.bounds)].rhoEnergyFluxX))/v5914)
      v5911[v5965].rho_t += ((-(v5911[v5965].rhoFluxY-v5911[((v5965+{0, -1, 0})%v5911.bounds)].rhoFluxY))/v5917)
      var v5968 = vs_div_double_3(vs_mul_double_3(vv_sub_double_3(v5911[v5965].rhoVelocityFluxY, v5911[((v5965+{0, -1, 0})%v5911.bounds)].rhoVelocityFluxY), double((-1))), v5917)
      var v5969 = v5911[v5965].rhoVelocity_t
      v5969[0] += v5968[0]
      v5969[1] += v5968[1]
      v5969[2] += v5968[2]
      v5911[v5965].rhoVelocity_t = v5969
      v5911[v5965].rhoEnergy_t += ((-(v5911[v5965].rhoEnergyFluxY-v5911[((v5965+{0, -1, 0})%v5911.bounds)].rhoEnergyFluxY))/v5917)
      v5911[v5965].rho_t += ((-(v5911[v5965].rhoFluxZ-v5911[((v5965+{0, 0, -1})%v5911.bounds)].rhoFluxZ))/v5920)
      var v5970 = vs_div_double_3(vs_mul_double_3(vv_sub_double_3(v5911[v5965].rhoVelocityFluxZ, v5911[((v5965+{0, 0, -1})%v5911.bounds)].rhoVelocityFluxZ), double((-1))), v5920)
      var v5971 = v5911[v5965].rhoVelocity_t
      v5971[0] += v5970[0]
      v5971[1] += v5970[1]
      v5971[2] += v5970[2]
      v5911[v5965].rhoVelocity_t = v5971
      v5911[v5965].rhoEnergy_t += ((-(v5911[v5965].rhoEnergyFluxZ-v5911[((v5965+{0, 0, -1})%v5911.bounds)].rhoEnergyFluxZ))/v5920)
    else
    end
  end
end

__demand(__parallel)
task Flow_AddBodyForces(v6049 : region(ispace(int3d), Fluid_columns), v6051 : double[3], v6052 : int32, v6053 : int32, v6054 : int32, v6055 : int32, v6056 : int32, v6057 : int32)
where
  reads(v6049.rho), reads(v6049.rhoEnergy_t), writes(v6049.rhoEnergy_t), reads(v6049.rhoVelocity_t), writes(v6049.rhoVelocity_t), reads(v6049.velocity)
do
  for v6073 in v6049 do
    if (not ((((((max(int32((uint64(v6052)-int3d(v6073).x)), int32(0))>int32(0)) or (max(int32((int3d(v6073).x-uint64(((v6053+v6052)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6054)-int3d(v6073).y)), int32(0))>int32(0))) or (max(int32((int3d(v6073).y-uint64(((v6055+v6054)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6056)-int3d(v6073).z)), int32(0))>int32(0))) or (max(int32((int3d(v6073).z-uint64(((v6057+v6056)-int32(1))))), int32(0))>int32(0)))) then
      var v6074 = vs_mul_double_3(v6051, v6049[v6073].rho)
      var v6075 = v6049[v6073].rhoVelocity_t
      v6075[0] += v6074[0]
      v6075[1] += v6074[1]
      v6075[2] += v6074[2]
      v6049[v6073].rhoVelocity_t = v6075
      v6049[v6073].rhoEnergy_t += (v6049[v6073].rho*dot_double_3(v6051, v6049[v6073].velocity))
    else
    end
  end
end

-- XXX: Turn off turbulent forcing
-- task Flow_UpdatePD(v6090 : region(ispace(int3d), Fluid_columns), v6092 : int32, v6093 : int32, v6094 : int32, v6095 : int32, v6096 : int32, v6097 : int32)
-- where
--   reads(v6090.PD), writes(v6090.PD), reads(v6090.pressure), reads(v6090.velocityGradientX), reads(v6090.velocityGradientY), reads(v6090.velocityGradientZ)
-- do
--   for v6105 in v6090 do
--     if (not ((((((max(int32((uint64(v6092)-int3d(v6105).x)), int32(0))>int32(0)) or (max(int32((int3d(v6105).x-uint64(((v6093+v6092)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6094)-int3d(v6105).y)), int32(0))>int32(0))) or (max(int32((int3d(v6105).y-uint64(((v6095+v6094)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6096)-int3d(v6105).z)), int32(0))>int32(0))) or (max(int32((int3d(v6105).z-uint64(((v6097+v6096)-int32(1))))), int32(0))>int32(0)))) then
--       var v6106 = double(double(0))
--       v6106 = ((v6090[v6105].velocityGradientX[int32(0)]+v6090[v6105].velocityGradientY[int32(1)])+v6090[v6105].velocityGradientZ[int32(2)])
--       v6090[v6105].PD = (v6106*v6090[v6105].pressure)
--     else
--     end
--   end
-- end

-- task Flow_ResetDissipation(v6113 : region(ispace(int3d), Fluid_columns))
-- where
--   reads(v6113.dissipation), writes(v6113.dissipation)
-- do
--   for v6116 in v6113 do
--     v6113[v6116].dissipation = double(0)
--   end
-- end

-- task Flow_ComputeDissipationX(v6122 : region(ispace(int3d), Fluid_columns), v6124 : double, v6125 : double, v6126 : double, v6127 : double, v6128 : double, v6129 : double, v6130 : int32, v6131 : int32, v6132 : double, v6133 : int32, v6134 : int32, v6135 : int32, v6136 : int32, v6137 : int32)
-- where
--   reads(v6122.dissipationFlux), writes(v6122.dissipationFlux), reads(v6122.temperature), reads(v6122.velocity), reads(v6122.velocityGradientY), reads(v6122.velocityGradientZ)
-- do
--   for v6223 in v6122 do
--     if ((not ((((((max(int32((uint64(v6131)-int3d(v6223).x)), int32(0))>int32(0)) or (max(int32((int3d(v6223).x-uint64(((v6133+v6131)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6134)-int3d(v6223).y)), int32(0))>int32(0))) or (max(int32((int3d(v6223).y-uint64(((v6135+v6134)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6136)-int3d(v6223).z)), int32(0))>int32(0))) or (max(int32((int3d(v6223).z-uint64(((v6137+v6136)-int32(1))))), int32(0))>int32(0)))) or (max(int32((uint64(v6131)-int3d(v6223).x)), int32(0))==int32(1))) then
--       var v6277 : double
--       do
--         var v6269 = v6122[((v6223+{1, 0, 0})%v6122.bounds)].temperature
--         var v6270 = v6124
--         var v6271 = v6125
--         var v6272 = v6126
--         var v6273 = v6127
--         var v6274 = v6128
--         var v6275 = v6129
--         var v6276 = v6130
--         var v4253 = double(double(0))
--         if (v6276==int32(0)) then
--           v4253 = v6270
--         else
--           if (v6276==int32(1)) then
--             v4253 = (v6272*C.pow((v6269/v6271), double(0.75)))
--           else
--             if (v6276==int32(2)) then
--               v4253 = ((v6275*C.pow((v6269/v6274), (double(3)/double(2))))*((v6274+v6273)/(v6269+v6273)))
--             else
--               regentlib.assert(false, "(Liszt assertion)")
--             end
--           end
--         end
--         v6277 = v4253
--       end
--       var v6278 = 0
--       var v6287 : double
--       do
--         var v6279 = v6122[v6223].temperature
--         var v6280 = v6124
--         var v6281 = v6125
--         var v6282 = v6126
--         var v6283 = v6127
--         var v6284 = v6128
--         var v6285 = v6129
--         var v6286 = v6130
--         var v4253 = double(double(0))
--         if (v6286==int32(0)) then
--           v4253 = v6280
--         else
--           if (v6286==int32(1)) then
--             v4253 = (v6282*C.pow((v6279/v6281), double(0.75)))
--           else
--             if (v6286==int32(2)) then
--               v4253 = ((v6285*C.pow((v6279/v6284), (double(3)/double(2))))*((v6284+v6283)/(v6279+v6283)))
--             else
--               regentlib.assert(false, "(Liszt assertion)")
--             end
--           end
--         end
--         v6287 = v4253
--       end
--       var v6288 = 0
--       var v6224 = (double(0.5)*(v6287+v6277))
--       var v6225 = [double[3]](array(double(0), double(0), double(0)))
--       var v6226 = double(double(0))
--       var v6227 = double(double(0))
--       var v6228 = double(double(0))
--       var v6229 = double(double(0))
--       v6225 = vs_mul_double_3(vv_add_double_3(v6122[v6223].velocity, v6122[((v6223+{1, 0, 0})%v6122.bounds)].velocity), double(0.5))
--       v6226 = (double(0.5)*(v6122[v6223].velocityGradientY[int32(0)]+v6122[((v6223+{1, 0, 0})%v6122.bounds)].velocityGradientY[int32(0)]))
--       v6227 = (double(0.5)*(v6122[v6223].velocityGradientZ[int32(0)]+v6122[((v6223+{1, 0, 0})%v6122.bounds)].velocityGradientZ[int32(0)]))
--       v6228 = (double(0.5)*(v6122[v6223].velocityGradientY[int32(1)]+v6122[((v6223+{1, 0, 0})%v6122.bounds)].velocityGradientY[int32(1)]))
--       v6229 = (double(0.5)*(v6122[v6223].velocityGradientZ[int32(2)]+v6122[((v6223+{1, 0, 0})%v6122.bounds)].velocityGradientZ[int32(2)]))
--       var v6230 = double(double(0))
--       var v6231 = double(double(0))
--       var v6232 = double(double(0))
--       var v6233 = double(double(0))
--       v6230 = (double(0.5)*(v6122[((v6223+{1, 0, 0})%v6122.bounds)].velocity[int32(0)]-v6122[v6223].velocity[int32(0)]))
--       v6231 = (double(0.5)*(v6122[((v6223+{1, 0, 0})%v6122.bounds)].velocity[int32(1)]-v6122[v6223].velocity[int32(1)]))
--       v6232 = (double(0.5)*(v6122[((v6223+{1, 0, 0})%v6122.bounds)].velocity[int32(2)]-v6122[v6223].velocity[int32(2)]))
--       v6233 = (double(0.5)*(v6122[((v6223+{1, 0, 0})%v6122.bounds)].temperature-v6122[v6223].temperature))
--       v6230 *= (1/(v6132*double(0.5)))
--       v6231 *= (1/(v6132*double(0.5)))
--       v6232 *= (1/(v6132*double(0.5)))
--       v6233 *= (1/(v6132*double(0.5)))
--       var v6234 = ((v6224*(((double(4)*v6230)-(double(2)*v6228))-(double(2)*v6229)))/double(3))
--       var v6235 = (v6224*(v6231+v6226))
--       var v6236 = (v6224*(v6232+v6227))
--       var v6237 = (((v6225[int32(0)]*v6234)+(v6225[int32(1)]*v6235))+(v6225[int32(2)]*v6236))
--       v6122[v6223].dissipationFlux = v6237
--     else
--     end
--   end
-- end
--
-- task Flow_UpdateDissipationX(v6325 : region(ispace(int3d), Fluid_columns), v6327 : int32, v6328 : double, v6329 : int32, v6330 : int32, v6331 : int32, v6332 : int32, v6333 : int32)
-- where
--   reads(v6325.dissipation), writes(v6325.dissipation), reads(v6325.dissipationFlux)
-- do
--   for v6335 in v6325 do
--     if (not ((((((max(int32((uint64(v6327)-int3d(v6335).x)), int32(0))>int32(0)) or (max(int32((int3d(v6335).x-uint64(((v6329+v6327)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6330)-int3d(v6335).y)), int32(0))>int32(0))) or (max(int32((int3d(v6335).y-uint64(((v6331+v6330)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6332)-int3d(v6335).z)), int32(0))>int32(0))) or (max(int32((int3d(v6335).z-uint64(((v6333+v6332)-int32(1))))), int32(0))>int32(0)))) then
--       v6325[v6335].dissipation += ((v6325[v6335].dissipationFlux-v6325[((v6335+{-1, 0, 0})%v6325.bounds)].dissipationFlux)/v6328)
--     else
--     end
--   end
-- end
--
-- task Flow_ComputeDissipationY(v6346 : region(ispace(int3d), Fluid_columns), v6348 : double, v6349 : double, v6350 : double, v6351 : double, v6352 : double, v6353 : double, v6354 : int32, v6355 : int32, v6356 : int32, v6357 : int32, v6358 : double, v6359 : int32, v6360 : int32, v6361 : int32)
-- where
--   reads(v6346.dissipationFlux), writes(v6346.dissipationFlux), reads(v6346.temperature), reads(v6346.velocity), reads(v6346.velocityGradientX), reads(v6346.velocityGradientZ)
-- do
--   for v6447 in v6346 do
--     if ((not ((((((max(int32((uint64(v6355)-int3d(v6447).x)), int32(0))>int32(0)) or (max(int32((int3d(v6447).x-uint64(((v6356+v6355)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6357)-int3d(v6447).y)), int32(0))>int32(0))) or (max(int32((int3d(v6447).y-uint64(((v6359+v6357)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6360)-int3d(v6447).z)), int32(0))>int32(0))) or (max(int32((int3d(v6447).z-uint64(((v6361+v6360)-int32(1))))), int32(0))>int32(0)))) or (max(int32((uint64(v6357)-int3d(v6447).y)), int32(0))==int32(1))) then
--       var v6501 : double
--       do
--         var v6493 = v6346[((v6447+{0, 1, 0})%v6346.bounds)].temperature
--         var v6494 = v6348
--         var v6495 = v6349
--         var v6496 = v6350
--         var v6497 = v6351
--         var v6498 = v6352
--         var v6499 = v6353
--         var v6500 = v6354
--         var v4253 = double(double(0))
--         if (v6500==int32(0)) then
--           v4253 = v6494
--         else
--           if (v6500==int32(1)) then
--             v4253 = (v6496*C.pow((v6493/v6495), double(0.75)))
--           else
--             if (v6500==int32(2)) then
--               v4253 = ((v6499*C.pow((v6493/v6498), (double(3)/double(2))))*((v6498+v6497)/(v6493+v6497)))
--             else
--               regentlib.assert(false, "(Liszt assertion)")
--             end
--           end
--         end
--         v6501 = v4253
--       end
--       var v6502 = 0
--       var v6511 : double
--       do
--         var v6503 = v6346[v6447].temperature
--         var v6504 = v6348
--         var v6505 = v6349
--         var v6506 = v6350
--         var v6507 = v6351
--         var v6508 = v6352
--         var v6509 = v6353
--         var v6510 = v6354
--         var v4253 = double(double(0))
--         if (v6510==int32(0)) then
--           v4253 = v6504
--         else
--           if (v6510==int32(1)) then
--             v4253 = (v6506*C.pow((v6503/v6505), double(0.75)))
--           else
--             if (v6510==int32(2)) then
--               v4253 = ((v6509*C.pow((v6503/v6508), (double(3)/double(2))))*((v6508+v6507)/(v6503+v6507)))
--             else
--               regentlib.assert(false, "(Liszt assertion)")
--             end
--           end
--         end
--         v6511 = v4253
--       end
--       var v6512 = 0
--       var v6448 = (double(0.5)*(v6511+v6501))
--       var v6449 = [double[3]](array(double(0), double(0), double(0)))
--       var v6450 = double(double(0))
--       var v6451 = double(double(0))
--       var v6452 = double(double(0))
--       var v6453 = double(double(0))
--       v6449 = vs_mul_double_3(vv_add_double_3(v6346[v6447].velocity, v6346[((v6447+{0, 1, 0})%v6346.bounds)].velocity), double(0.5))
--       v6450 = (double(0.5)*(v6346[v6447].velocityGradientX[int32(1)]+v6346[((v6447+{0, 1, 0})%v6346.bounds)].velocityGradientX[int32(1)]))
--       v6451 = (double(0.5)*(v6346[v6447].velocityGradientZ[int32(1)]+v6346[((v6447+{0, 1, 0})%v6346.bounds)].velocityGradientZ[int32(1)]))
--       v6452 = (double(0.5)*(v6346[v6447].velocityGradientX[int32(0)]+v6346[((v6447+{0, 1, 0})%v6346.bounds)].velocityGradientX[int32(0)]))
--       v6453 = (double(0.5)*(v6346[v6447].velocityGradientZ[int32(2)]+v6346[((v6447+{0, 1, 0})%v6346.bounds)].velocityGradientZ[int32(2)]))
--       var v6454 = double(double(0))
--       var v6455 = double(double(0))
--       var v6456 = double(double(0))
--       var v6457 = double(double(0))
--       v6454 = (double(0.5)*(v6346[((v6447+{0, 1, 0})%v6346.bounds)].velocity[int32(0)]-v6346[v6447].velocity[int32(0)]))
--       v6455 = (double(0.5)*(v6346[((v6447+{0, 1, 0})%v6346.bounds)].velocity[int32(1)]-v6346[v6447].velocity[int32(1)]))
--       v6456 = (double(0.5)*(v6346[((v6447+{0, 1, 0})%v6346.bounds)].velocity[int32(2)]-v6346[v6447].velocity[int32(2)]))
--       v6457 = (double(0.5)*(v6346[((v6447+{0, 1, 0})%v6346.bounds)].temperature-v6346[v6447].temperature))
--       v6454 *= (1/(v6358*double(0.5)))
--       v6455 *= (1/(v6358*double(0.5)))
--       v6456 *= (1/(v6358*double(0.5)))
--       v6457 *= (1/(v6358*double(0.5)))
--       var v6458 = (v6448*(v6454+v6450))
--       var v6459 = ((v6448*(((double(4)*v6455)-(double(2)*v6452))-(double(2)*v6453)))/double(3))
--       var v6460 = (v6448*(v6456+v6451))
--       var v6461 = (((v6449[int32(0)]*v6458)+(v6449[int32(1)]*v6459))+(v6449[int32(2)]*v6460))
--       v6346[v6447].dissipationFlux = v6461
--     else
--     end
--   end
-- end
--
-- task Flow_UpdateDissipationY(v6549 : region(ispace(int3d), Fluid_columns), v6551 : int32, v6552 : int32, v6553 : int32, v6554 : double, v6555 : int32, v6556 : int32, v6557 : int32)
-- where
--   reads(v6549.dissipation), writes(v6549.dissipation), reads(v6549.dissipationFlux)
-- do
--   for v6559 in v6549 do
--     if (not ((((((max(int32((uint64(v6551)-int3d(v6559).x)), int32(0))>int32(0)) or (max(int32((int3d(v6559).x-uint64(((v6552+v6551)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6553)-int3d(v6559).y)), int32(0))>int32(0))) or (max(int32((int3d(v6559).y-uint64(((v6555+v6553)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6556)-int3d(v6559).z)), int32(0))>int32(0))) or (max(int32((int3d(v6559).z-uint64(((v6557+v6556)-int32(1))))), int32(0))>int32(0)))) then
--       v6549[v6559].dissipation += ((v6549[v6559].dissipationFlux-v6549[((v6559+{0, -1, 0})%v6549.bounds)].dissipationFlux)/v6554)
--     else
--     end
--   end
-- end
--
-- task Flow_ComputeDissipationZ(v6570 : region(ispace(int3d), Fluid_columns), v6572 : double, v6573 : double, v6574 : double, v6575 : double, v6576 : double, v6577 : double, v6578 : int32, v6579 : int32, v6580 : int32, v6581 : int32, v6582 : int32, v6583 : int32, v6584 : double, v6585 : int32)
-- where
--   reads(v6570.dissipationFlux), writes(v6570.dissipationFlux), reads(v6570.temperature), reads(v6570.velocity), reads(v6570.velocityGradientX), reads(v6570.velocityGradientY)
-- do
--   for v6671 in v6570 do
--     if ((not ((((((max(int32((uint64(v6579)-int3d(v6671).x)), int32(0))>int32(0)) or (max(int32((int3d(v6671).x-uint64(((v6580+v6579)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6581)-int3d(v6671).y)), int32(0))>int32(0))) or (max(int32((int3d(v6671).y-uint64(((v6582+v6581)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6583)-int3d(v6671).z)), int32(0))>int32(0))) or (max(int32((int3d(v6671).z-uint64(((v6585+v6583)-int32(1))))), int32(0))>int32(0)))) or (max(int32((uint64(v6583)-int3d(v6671).z)), int32(0))==int32(1))) then
--       var v6725 : double
--       do
--         var v6717 = v6570[((v6671+{0, 0, 1})%v6570.bounds)].temperature
--         var v6718 = v6572
--         var v6719 = v6573
--         var v6720 = v6574
--         var v6721 = v6575
--         var v6722 = v6576
--         var v6723 = v6577
--         var v6724 = v6578
--         var v4253 = double(double(0))
--         if (v6724==int32(0)) then
--           v4253 = v6718
--         else
--           if (v6724==int32(1)) then
--             v4253 = (v6720*C.pow((v6717/v6719), double(0.75)))
--           else
--             if (v6724==int32(2)) then
--               v4253 = ((v6723*C.pow((v6717/v6722), (double(3)/double(2))))*((v6722+v6721)/(v6717+v6721)))
--             else
--               regentlib.assert(false, "(Liszt assertion)")
--             end
--           end
--         end
--         v6725 = v4253
--       end
--       var v6726 = 0
--       var v6735 : double
--       do
--         var v6727 = v6570[v6671].temperature
--         var v6728 = v6572
--         var v6729 = v6573
--         var v6730 = v6574
--         var v6731 = v6575
--         var v6732 = v6576
--         var v6733 = v6577
--         var v6734 = v6578
--         var v4253 = double(double(0))
--         if (v6734==int32(0)) then
--           v4253 = v6728
--         else
--           if (v6734==int32(1)) then
--             v4253 = (v6730*C.pow((v6727/v6729), double(0.75)))
--           else
--             if (v6734==int32(2)) then
--               v4253 = ((v6733*C.pow((v6727/v6732), (double(3)/double(2))))*((v6732+v6731)/(v6727+v6731)))
--             else
--               regentlib.assert(false, "(Liszt assertion)")
--             end
--           end
--         end
--         v6735 = v4253
--       end
--       var v6736 = 0
--       var v6672 = (double(0.5)*(v6735+v6725))
--       var v6673 = [double[3]](array(double(0), double(0), double(0)))
--       var v6674 = double(double(0))
--       var v6675 = double(double(0))
--       var v6676 = double(double(0))
--       var v6677 = double(double(0))
--       v6673 = vs_mul_double_3(vv_add_double_3(v6570[v6671].velocity, v6570[((v6671+{0, 0, 1})%v6570.bounds)].velocity), double(0.5))
--       v6674 = (double(0.5)*(v6570[v6671].velocityGradientX[int32(2)]+v6570[((v6671+{0, 0, 1})%v6570.bounds)].velocityGradientX[int32(2)]))
--       v6675 = (double(0.5)*(v6570[v6671].velocityGradientY[int32(2)]+v6570[((v6671+{0, 0, 1})%v6570.bounds)].velocityGradientY[int32(2)]))
--       v6676 = (double(0.5)*(v6570[v6671].velocityGradientX[int32(0)]+v6570[((v6671+{0, 0, 1})%v6570.bounds)].velocityGradientX[int32(0)]))
--       v6677 = (double(0.5)*(v6570[v6671].velocityGradientY[int32(1)]+v6570[((v6671+{0, 0, 1})%v6570.bounds)].velocityGradientY[int32(1)]))
--       var v6678 = double(double(0))
--       var v6679 = double(double(0))
--       var v6680 = double(double(0))
--       var v6681 = double(double(0))
--       v6678 = (double(0.5)*(v6570[((v6671+{0, 0, 1})%v6570.bounds)].velocity[int32(0)]-v6570[v6671].velocity[int32(0)]))
--       v6679 = (double(0.5)*(v6570[((v6671+{0, 0, 1})%v6570.bounds)].velocity[int32(1)]-v6570[v6671].velocity[int32(1)]))
--       v6680 = (double(0.5)*(v6570[((v6671+{0, 0, 1})%v6570.bounds)].velocity[int32(2)]-v6570[v6671].velocity[int32(2)]))
--       v6681 = (double(0.5)*(v6570[((v6671+{0, 0, 1})%v6570.bounds)].temperature-v6570[v6671].temperature))
--       v6678 *= (1/(v6584*double(0.5)))
--       v6679 *= (1/(v6584*double(0.5)))
--       v6680 *= (1/(v6584*double(0.5)))
--       v6681 *= (1/(v6584*double(0.5)))
--       var v6682 = (v6672*(v6678+v6674))
--       var v6683 = (v6672*(v6679+v6675))
--       var v6684 = ((v6672*(((double(4)*v6680)-(double(2)*v6676))-(double(2)*v6677)))/double(3))
--       var v6685 = (((v6673[int32(0)]*v6682)+(v6673[int32(1)]*v6683))+(v6673[int32(2)]*v6684))
--       v6570[v6671].dissipationFlux = v6685
--     else
--     end
--   end
-- end
--
-- task Flow_UpdateDissipationZ(v6773 : region(ispace(int3d), Fluid_columns), v6775 : int32, v6776 : int32, v6777 : int32, v6778 : int32, v6779 : int32, v6780 : double, v6781 : int32)
-- where
--   reads(v6773.dissipation), writes(v6773.dissipation), reads(v6773.dissipationFlux)
-- do
--   for v6783 in v6773 do
--     if (not ((((((max(int32((uint64(v6775)-int3d(v6783).x)), int32(0))>int32(0)) or (max(int32((int3d(v6783).x-uint64(((v6776+v6775)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6777)-int3d(v6783).y)), int32(0))>int32(0))) or (max(int32((int3d(v6783).y-uint64(((v6778+v6777)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6779)-int3d(v6783).z)), int32(0))>int32(0))) or (max(int32((int3d(v6783).z-uint64(((v6781+v6779)-int32(1))))), int32(0))>int32(0)))) then
--       v6773[v6783].dissipation += ((v6773[v6783].dissipationFlux-v6773[((v6783+{0, 0, -1})%v6773.bounds)].dissipationFlux)/v6780)
--     else
--     end
--   end
-- end
--
-- task CalculateAveragePD(v6794 : region(ispace(int3d), Fluid_columns), v6797 : double, v6798 : int32, v6799 : int32, v6800 : int32, v6801 : int32, v6802 : int32, v6803 : int32) : double
-- where
--   reads(v6794.PD)
-- do
--   var v6806 = double(0)
--   for v6807 in v6794 do
--     if (not ((((((max(int32((uint64(v6798)-int3d(v6807).x)), int32(0))>int32(0)) or (max(int32((int3d(v6807).x-uint64(((v6799+v6798)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6800)-int3d(v6807).y)), int32(0))>int32(0))) or (max(int32((int3d(v6807).y-uint64(((v6801+v6800)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6802)-int3d(v6807).z)), int32(0))>int32(0))) or (max(int32((int3d(v6807).z-uint64(((v6803+v6802)-int32(1))))), int32(0))>int32(0)))) then
--       v6806 += (v6794[v6807].PD*v6797)
--     else
--     end
--   end
--   return v6806
-- end
--
-- task CalculateAverageDissipation(v6814 : region(ispace(int3d), Fluid_columns), v6817 : double, v6818 : int32, v6819 : int32, v6820 : int32, v6821 : int32, v6822 : int32, v6823 : int32) : double
-- where
--   reads(v6814.dissipation)
-- do
--   var v6826 = double(0)
--   for v6827 in v6814 do
--     if (not ((((((max(int32((uint64(v6818)-int3d(v6827).x)), int32(0))>int32(0)) or (max(int32((int3d(v6827).x-uint64(((v6819+v6818)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6820)-int3d(v6827).y)), int32(0))>int32(0))) or (max(int32((int3d(v6827).y-uint64(((v6821+v6820)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6822)-int3d(v6827).z)), int32(0))>int32(0))) or (max(int32((int3d(v6827).z-uint64(((v6823+v6822)-int32(1))))), int32(0))>int32(0)))) then
--       v6826 += (v6814[v6827].dissipation*v6817)
--     else
--     end
--   end
--   return v6826
-- end
--
-- task CalculateAverageK(v6834 : region(ispace(int3d), Fluid_columns), v6837 : double, v6838 : int32, v6839 : int32, v6840 : int32, v6841 : int32, v6842 : int32, v6843 : int32) : double
-- where
--   reads(v6834.rho), reads(v6834.velocity)
-- do
--   var v6846 = double(0)
--   for v6847 in v6834 do
--     if (not ((((((max(int32((uint64(v6838)-int3d(v6847).x)), int32(0))>int32(0)) or (max(int32((int3d(v6847).x-uint64(((v6839+v6838)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6840)-int3d(v6847).y)), int32(0))>int32(0))) or (max(int32((int3d(v6847).y-uint64(((v6841+v6840)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6842)-int3d(v6847).z)), int32(0))>int32(0))) or (max(int32((int3d(v6847).z-uint64(((v6843+v6842)-int32(1))))), int32(0))>int32(0)))) then
--       v6846 += (((double(0.5)*v6834[v6847].rho)*dot_double_3(v6834[v6847].velocity, v6834[v6847].velocity))*v6837)
--     else
--     end
--   end
--   return v6846
-- end
--
-- task Flow_AddTurbulentSource(v6858 : region(ispace(int3d), Fluid_columns), v6860 : double, v6862 : double, v6863 : double, v6864 : double, v6865 : int32, v6866 : int32, v6867 : int32, v6868 : int32, v6869 : int32, v6870 : int32) : double
-- where
--   reads(v6858.rho), reads(v6858.rhoEnergy_t), writes(v6858.rhoEnergy_t), reads(v6858.rhoVelocity_t), writes(v6858.rhoVelocity_t), reads(v6858.velocity)
-- do
--   var v6923 = double(0)
--   for v6924 in v6858 do
--     if (not ((((((max(int32((uint64(v6865)-int3d(v6924).x)), int32(0))>int32(0)) or (max(int32((int3d(v6924).x-uint64(((v6866+v6865)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6867)-int3d(v6924).y)), int32(0))>int32(0))) or (max(int32((int3d(v6924).y-uint64(((v6868+v6867)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6869)-int3d(v6924).z)), int32(0))>int32(0))) or (max(int32((int3d(v6924).z-uint64(((v6870+v6869)-int32(1))))), int32(0))>int32(0)))) then
--       var v6925 = double(double(0))
--       var v6926 = double(double(0))
--       var v6927 = double(double(0))
--       var v6928 = double(double(0))
--       var v6929 = double(double(0))
--       var v6930 = [double[3]](array(double(0), double(0), double(0)))
--       v6925 = (v6863+v6860)
--       v6927 = double(300)
--       v6928 = double(3.00889e-06)
--       v6929 = double(66.27348)
--       v6926 = (((-v6925)-((v6927*(v6862-v6929))/v6928))/(double(2)*v6862))
--       v6930 = vs_mul_double_3(v6858[v6924].velocity, (v6858[v6924].rho*v6926))
--       var v6931 = v6930
--       var v6932 = v6858[v6924].rhoVelocity_t
--       v6932[0] += v6931[0]
--       v6932[1] += v6931[1]
--       v6932[2] += v6931[2]
--       v6858[v6924].rhoVelocity_t = v6932
--       v6858[v6924].rhoEnergy_t += dot_double_3(v6930, v6858[v6924].velocity)
--       v6923 += (dot_double_3(v6930, v6858[v6924].velocity)*v6864)
--     else
--     end
--   end
--   return v6923
-- end
--
-- task Flow_AdjustTurbulentSource(v6963 : region(ispace(int3d), Fluid_columns), v6965 : double, v6966 : int32, v6967 : int32, v6968 : int32, v6969 : int32, v6970 : int32, v6971 : int32)
-- where
--   reads(v6963.rhoEnergy_t), writes(v6963.rhoEnergy_t)
-- do
--   for v6973 in v6963 do
--     if (not ((((((max(int32((uint64(v6966)-int3d(v6973).x)), int32(0))>int32(0)) or (max(int32((int3d(v6973).x-uint64(((v6967+v6966)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6968)-int3d(v6973).y)), int32(0))>int32(0))) or (max(int32((int3d(v6973).y-uint64(((v6969+v6968)-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(v6970)-int3d(v6973).z)), int32(0))>int32(0))) or (max(int32((int3d(v6973).z-uint64(((v6971+v6970)-int32(1))))), int32(0))>int32(0)))) then
--       v6963[v6973].rhoEnergy_t += (-v6965)
--     else
--     end
--   end
-- end
--
task Particles_LocateInCells(v6980 : region(ispace(int1d), particles_columns), v6982 : bool, v6983 : bool, v6984 : bool, v6985 : int32, v6986 : int32, v6987 : double, v6988 : double, v6989 : int32, v6990 : int32, v6991 : double, v6992 : double, v6993 : int32, v6994 : int32, v6995 : double, v6996 : double)
where
  reads(v6980.cell), writes(v6980.cell), reads(v6980.position), reads(v6980.__valid)
do
  for v7095 in v6980 do
    if v6980[v7095].__valid then
      var v7114 : int3d
      do
        var v7098 = v6980[v7095].position
        var v7099 = v6982
        var v7100 = v6983
        var v7101 = v6984
        var v7102 = v6985
        var v7103 = v6986
        var v7104 = v6987
        var v7105 = v6988
        var v7106 = v6989
        var v7107 = v6990
        var v7108 = v6991
        var v7109 = v6992
        var v7110 = v6993
        var v7111 = v6994
        var v7112 = v6995
        var v7113 = v6996
        var v7058 = (v7105/double(v7103))
        var v7059 = (v7104-(double(v7102)*v7058))
        var v7060 = ((v7098[int32(0)]-v7059)/v7058)
        var v7061 = (v7103+(int32(2)*v7102))
        var v7062 : uint64
        if v7099 then
          v7062 = (uint64((C.fmod(v7060, double(v7061))+double(v7061)))%uint64(v7061))
        else
          v7062 = uint64(max(double(0), min(double((v7061-int32(1))), v7060)))
        end
        var v7063 = (v7109/double(v7107))
        var v7064 = (v7108-(double(v7106)*v7063))
        var v7065 = ((v7098[int32(1)]-v7064)/v7063)
        var v7066 = (v7107+(int32(2)*v7106))
        var v7067 : uint64
        if v7100 then
          v7067 = (uint64((C.fmod(v7065, double(v7066))+double(v7066)))%uint64(v7066))
        else
          v7067 = uint64(max(double(0), min(double((v7066-int32(1))), v7065)))
        end
        var v7068 = (v7113/double(v7111))
        var v7069 = (v7112-(double(v7110)*v7068))
        var v7070 = ((v7098[int32(2)]-v7069)/v7068)
        var v7071 = (v7111+(int32(2)*v7110))
        var v7072 : uint64
        if v7101 then
          v7072 = (uint64((C.fmod(v7070, double(v7071))+double(v7071)))%uint64(v7071))
        else
          v7072 = uint64(max(double(0), min(double((v7071-int32(1))), v7070)))
        end
        v7114 = int3d({v7062, v7067, v7072})
      end
      var v7115 = 0
      v6980[v7095].cell = v7114
    else
    end
  end
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

task particles_pushAll(v7134 : int3d, v7132 : region(ispace(int1d), particles_columns), v7138 : region(ispace(int1d), int8[376]), v7150 : region(ispace(int1d), int8[376]), v7162 : region(ispace(int1d), int8[376]), v7174 : region(ispace(int1d), int8[376]), v7186 : region(ispace(int1d), int8[376]), v7198 : region(ispace(int1d), int8[376]), v7210 : region(ispace(int1d), int8[376]), v7222 : region(ispace(int1d), int8[376]), v7234 : region(ispace(int1d), int8[376]), v7246 : region(ispace(int1d), int8[376]), v7258 : region(ispace(int1d), int8[376]), v7270 : region(ispace(int1d), int8[376]), v7282 : region(ispace(int1d), int8[376]), v7294 : region(ispace(int1d), int8[376]), v7306 : region(ispace(int1d), int8[376]), v7318 : region(ispace(int1d), int8[376]), v7330 : region(ispace(int1d), int8[376]), v7342 : region(ispace(int1d), int8[376]), v7354 : region(ispace(int1d), int8[376]), v7366 : region(ispace(int1d), int8[376]), v7378 : region(ispace(int1d), int8[376]), v7390 : region(ispace(int1d), int8[376]), v7402 : region(ispace(int1d), int8[376]), v7414 : region(ispace(int1d), int8[376]), v7426 : region(ispace(int1d), int8[376]), v7438 : region(ispace(int1d), int8[376]), v7449 : int32, v7450 : int32, v7451 : int32, v7452 : int32, v7453 : int32, v7454 : int32, v7118 : int32, v7119 : int32, v7120 : int32)
where
  reads(v7132), writes(v7132.__valid), reads(v7138), writes(v7138), reads(v7150), writes(v7150), reads(v7162), writes(v7162), reads(v7174), writes(v7174), reads(v7186), writes(v7186), reads(v7198), writes(v7198), reads(v7210), writes(v7210), reads(v7222), writes(v7222), reads(v7234), writes(v7234), reads(v7246), writes(v7246), reads(v7258), writes(v7258), reads(v7270), writes(v7270), reads(v7282), writes(v7282), reads(v7294), writes(v7294), reads(v7306), writes(v7306), reads(v7318), writes(v7318), reads(v7330), writes(v7330), reads(v7342), writes(v7342), reads(v7354), writes(v7354), reads(v7366), writes(v7366), reads(v7378), writes(v7378), reads(v7390), writes(v7390), reads(v7402), writes(v7402), reads(v7414), writes(v7414), reads(v7426), writes(v7426), reads(v7438), writes(v7438)
do
  for v7455 in v7138 do
    v7138[v7455][368LL] = int8(false)
  end
  var v7456 = particles_getBasePointer(__physical(v7138)[0], __fields(v7138)[0], __runtime())
  for v7457 in v7150 do
    v7150[v7457][368LL] = int8(false)
  end
  var v7458 = particles_getBasePointer(__physical(v7150)[0], __fields(v7150)[0], __runtime())
  for v7459 in v7162 do
    v7162[v7459][368LL] = int8(false)
  end
  var v7460 = particles_getBasePointer(__physical(v7162)[0], __fields(v7162)[0], __runtime())
  for v7461 in v7174 do
    v7174[v7461][368LL] = int8(false)
  end
  var v7462 = particles_getBasePointer(__physical(v7174)[0], __fields(v7174)[0], __runtime())
  for v7463 in v7186 do
    v7186[v7463][368LL] = int8(false)
  end
  var v7464 = particles_getBasePointer(__physical(v7186)[0], __fields(v7186)[0], __runtime())
  for v7465 in v7198 do
    v7198[v7465][368LL] = int8(false)
  end
  var v7466 = particles_getBasePointer(__physical(v7198)[0], __fields(v7198)[0], __runtime())
  for v7467 in v7210 do
    v7210[v7467][368LL] = int8(false)
  end
  var v7468 = particles_getBasePointer(__physical(v7210)[0], __fields(v7210)[0], __runtime())
  for v7469 in v7222 do
    v7222[v7469][368LL] = int8(false)
  end
  var v7470 = particles_getBasePointer(__physical(v7222)[0], __fields(v7222)[0], __runtime())
  for v7471 in v7234 do
    v7234[v7471][368LL] = int8(false)
  end
  var v7472 = particles_getBasePointer(__physical(v7234)[0], __fields(v7234)[0], __runtime())
  for v7473 in v7246 do
    v7246[v7473][368LL] = int8(false)
  end
  var v7474 = particles_getBasePointer(__physical(v7246)[0], __fields(v7246)[0], __runtime())
  for v7475 in v7258 do
    v7258[v7475][368LL] = int8(false)
  end
  var v7476 = particles_getBasePointer(__physical(v7258)[0], __fields(v7258)[0], __runtime())
  for v7477 in v7270 do
    v7270[v7477][368LL] = int8(false)
  end
  var v7478 = particles_getBasePointer(__physical(v7270)[0], __fields(v7270)[0], __runtime())
  for v7479 in v7282 do
    v7282[v7479][368LL] = int8(false)
  end
  var v7480 = particles_getBasePointer(__physical(v7282)[0], __fields(v7282)[0], __runtime())
  for v7481 in v7294 do
    v7294[v7481][368LL] = int8(false)
  end
  var v7482 = particles_getBasePointer(__physical(v7294)[0], __fields(v7294)[0], __runtime())
  for v7483 in v7306 do
    v7306[v7483][368LL] = int8(false)
  end
  var v7484 = particles_getBasePointer(__physical(v7306)[0], __fields(v7306)[0], __runtime())
  for v7485 in v7318 do
    v7318[v7485][368LL] = int8(false)
  end
  var v7486 = particles_getBasePointer(__physical(v7318)[0], __fields(v7318)[0], __runtime())
  for v7487 in v7330 do
    v7330[v7487][368LL] = int8(false)
  end
  var v7488 = particles_getBasePointer(__physical(v7330)[0], __fields(v7330)[0], __runtime())
  for v7489 in v7342 do
    v7342[v7489][368LL] = int8(false)
  end
  var v7490 = particles_getBasePointer(__physical(v7342)[0], __fields(v7342)[0], __runtime())
  for v7491 in v7354 do
    v7354[v7491][368LL] = int8(false)
  end
  var v7492 = particles_getBasePointer(__physical(v7354)[0], __fields(v7354)[0], __runtime())
  for v7493 in v7366 do
    v7366[v7493][368LL] = int8(false)
  end
  var v7494 = particles_getBasePointer(__physical(v7366)[0], __fields(v7366)[0], __runtime())
  for v7495 in v7378 do
    v7378[v7495][368LL] = int8(false)
  end
  var v7496 = particles_getBasePointer(__physical(v7378)[0], __fields(v7378)[0], __runtime())
  for v7497 in v7390 do
    v7390[v7497][368LL] = int8(false)
  end
  var v7498 = particles_getBasePointer(__physical(v7390)[0], __fields(v7390)[0], __runtime())
  for v7499 in v7402 do
    v7402[v7499][368LL] = int8(false)
  end
  var v7500 = particles_getBasePointer(__physical(v7402)[0], __fields(v7402)[0], __runtime())
  for v7501 in v7414 do
    v7414[v7501][368LL] = int8(false)
  end
  var v7502 = particles_getBasePointer(__physical(v7414)[0], __fields(v7414)[0], __runtime())
  for v7503 in v7426 do
    v7426[v7503][368LL] = int8(false)
  end
  var v7504 = particles_getBasePointer(__physical(v7426)[0], __fields(v7426)[0], __runtime())
  for v7505 in v7438 do
    v7438[v7505][368LL] = int8(false)
  end
  var v7506 = particles_getBasePointer(__physical(v7438)[0], __fields(v7438)[0], __runtime())
  for v7507 in v7132 do
    if v7507.__valid then
      var v8202 : int3d
      do
        var v8192 = v7507.cell
        var v8193 = v7449
        var v8194 = v7450
        var v8195 = v7451
        var v8196 = v7452
        var v8197 = v7453
        var v8198 = v7454
        var v8199 = v7118
        var v8200 = v7119
        var v8201 = v7120
        v8192.x = min(max(v8192.x, v8196), ((v8193+v8196)-1))
        v8192.y = min(max(v8192.y, v8197), ((v8194+v8197)-1))
        v8192.z = min(max(v8192.z, v8198), ((v8195+v8198)-1))
        v8202 = int3d({((v8192.x-v8196)/(v8193/v8199)), ((v8192.y-v8197)/(v8194/v8200)), ((v8192.z-v8198)/(v8195/v8201))})
      end
      var v8203 = 0
      var v7508 = v8202
      if (v7508~=v7134) then
        do
          var v7509 = int3d({0, 0, 1})
          if (v7507.__valid and (v7508==(((v7134+v7509)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7510 = 0
            for v7511 in v7138 do
              if (not bool(v7138[v7511][368LL])) then
                particles_pushElement(v7456, v7510, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7138[v7511][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7510 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7512 = int3d({0, 0, -1})
          if (v7507.__valid and (v7508==(((v7134+v7512)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7513 = 0
            for v7514 in v7150 do
              if (not bool(v7150[v7514][368LL])) then
                particles_pushElement(v7458, v7513, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7150[v7514][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7513 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7515 = int3d({0, 1, 0})
          if (v7507.__valid and (v7508==(((v7134+v7515)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7516 = 0
            for v7517 in v7162 do
              if (not bool(v7162[v7517][368LL])) then
                particles_pushElement(v7460, v7516, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7162[v7517][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7516 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7518 = int3d({0, 1, 1})
          if (v7507.__valid and (v7508==(((v7134+v7518)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7519 = 0
            for v7520 in v7174 do
              if (not bool(v7174[v7520][368LL])) then
                particles_pushElement(v7462, v7519, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7174[v7520][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7519 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7521 = int3d({0, 1, -1})
          if (v7507.__valid and (v7508==(((v7134+v7521)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7522 = 0
            for v7523 in v7186 do
              if (not bool(v7186[v7523][368LL])) then
                particles_pushElement(v7464, v7522, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7186[v7523][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7522 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7524 = int3d({0, -1, 0})
          if (v7507.__valid and (v7508==(((v7134+v7524)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7525 = 0
            for v7526 in v7198 do
              if (not bool(v7198[v7526][368LL])) then
                particles_pushElement(v7466, v7525, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7198[v7526][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7525 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7527 = int3d({0, -1, 1})
          if (v7507.__valid and (v7508==(((v7134+v7527)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7528 = 0
            for v7529 in v7210 do
              if (not bool(v7210[v7529][368LL])) then
                particles_pushElement(v7468, v7528, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7210[v7529][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7528 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7530 = int3d({0, -1, -1})
          if (v7507.__valid and (v7508==(((v7134+v7530)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7531 = 0
            for v7532 in v7222 do
              if (not bool(v7222[v7532][368LL])) then
                particles_pushElement(v7470, v7531, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7222[v7532][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7531 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7533 = int3d({1, 0, 0})
          if (v7507.__valid and (v7508==(((v7134+v7533)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7534 = 0
            for v7535 in v7234 do
              if (not bool(v7234[v7535][368LL])) then
                particles_pushElement(v7472, v7534, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7234[v7535][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7534 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7536 = int3d({1, 0, 1})
          if (v7507.__valid and (v7508==(((v7134+v7536)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7537 = 0
            for v7538 in v7246 do
              if (not bool(v7246[v7538][368LL])) then
                particles_pushElement(v7474, v7537, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7246[v7538][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7537 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7539 = int3d({1, 0, -1})
          if (v7507.__valid and (v7508==(((v7134+v7539)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7540 = 0
            for v7541 in v7258 do
              if (not bool(v7258[v7541][368LL])) then
                particles_pushElement(v7476, v7540, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7258[v7541][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7540 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7542 = int3d({1, 1, 0})
          if (v7507.__valid and (v7508==(((v7134+v7542)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7543 = 0
            for v7544 in v7270 do
              if (not bool(v7270[v7544][368LL])) then
                particles_pushElement(v7478, v7543, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7270[v7544][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7543 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7545 = int3d({1, 1, 1})
          if (v7507.__valid and (v7508==(((v7134+v7545)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7546 = 0
            for v7547 in v7282 do
              if (not bool(v7282[v7547][368LL])) then
                particles_pushElement(v7480, v7546, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7282[v7547][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7546 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7548 = int3d({1, 1, -1})
          if (v7507.__valid and (v7508==(((v7134+v7548)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7549 = 0
            for v7550 in v7294 do
              if (not bool(v7294[v7550][368LL])) then
                particles_pushElement(v7482, v7549, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7294[v7550][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7549 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7551 = int3d({1, -1, 0})
          if (v7507.__valid and (v7508==(((v7134+v7551)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7552 = 0
            for v7553 in v7306 do
              if (not bool(v7306[v7553][368LL])) then
                particles_pushElement(v7484, v7552, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7306[v7553][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7552 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7554 = int3d({1, -1, 1})
          if (v7507.__valid and (v7508==(((v7134+v7554)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7555 = 0
            for v7556 in v7318 do
              if (not bool(v7318[v7556][368LL])) then
                particles_pushElement(v7486, v7555, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7318[v7556][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7555 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7557 = int3d({1, -1, -1})
          if (v7507.__valid and (v7508==(((v7134+v7557)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7558 = 0
            for v7559 in v7330 do
              if (not bool(v7330[v7559][368LL])) then
                particles_pushElement(v7488, v7558, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7330[v7559][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7558 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7560 = int3d({-1, 0, 0})
          if (v7507.__valid and (v7508==(((v7134+v7560)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7561 = 0
            for v7562 in v7342 do
              if (not bool(v7342[v7562][368LL])) then
                particles_pushElement(v7490, v7561, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7342[v7562][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7561 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7563 = int3d({-1, 0, 1})
          if (v7507.__valid and (v7508==(((v7134+v7563)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7564 = 0
            for v7565 in v7354 do
              if (not bool(v7354[v7565][368LL])) then
                particles_pushElement(v7492, v7564, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7354[v7565][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7564 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7566 = int3d({-1, 0, -1})
          if (v7507.__valid and (v7508==(((v7134+v7566)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7567 = 0
            for v7568 in v7366 do
              if (not bool(v7366[v7568][368LL])) then
                particles_pushElement(v7494, v7567, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7366[v7568][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7567 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7569 = int3d({-1, 1, 0})
          if (v7507.__valid and (v7508==(((v7134+v7569)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7570 = 0
            for v7571 in v7378 do
              if (not bool(v7378[v7571][368LL])) then
                particles_pushElement(v7496, v7570, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7378[v7571][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7570 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7572 = int3d({-1, 1, 1})
          if (v7507.__valid and (v7508==(((v7134+v7572)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7573 = 0
            for v7574 in v7390 do
              if (not bool(v7390[v7574][368LL])) then
                particles_pushElement(v7498, v7573, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7390[v7574][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7573 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7575 = int3d({-1, 1, -1})
          if (v7507.__valid and (v7508==(((v7134+v7575)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7576 = 0
            for v7577 in v7402 do
              if (not bool(v7402[v7577][368LL])) then
                particles_pushElement(v7500, v7576, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7402[v7577][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7576 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7578 = int3d({-1, -1, 0})
          if (v7507.__valid and (v7508==(((v7134+v7578)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7579 = 0
            for v7580 in v7414 do
              if (not bool(v7414[v7580][368LL])) then
                particles_pushElement(v7502, v7579, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7414[v7580][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7579 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7581 = int3d({-1, -1, 1})
          if (v7507.__valid and (v7508==(((v7134+v7581)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7582 = 0
            for v7583 in v7426 do
              if (not bool(v7426[v7583][368LL])) then
                particles_pushElement(v7504, v7582, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7426[v7583][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7582 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var v7584 = int3d({-1, -1, -1})
          if (v7507.__valid and (v7508==(((v7134+v7584)+{v7118, v7119, v7120})%{v7118, v7119, v7120}))) then
            var v7585 = 0
            for v7586 in v7438 do
              if (not bool(v7438[v7586][368LL])) then
                particles_pushElement(v7506, v7585, v7132[v7507])
                v7507.__valid = false
                regentlib.assert(bool(v7438[v7586][368LL]), "Element did not get copied properly")
                break
              else
              end
              v7585 += 1
            end
            regentlib.assert((not v7507.__valid), "Transfer queue ran out of space")
          else
          end
        end
        regentlib.assert((not v7507.__valid), "Element moved past predicted stencil")
      else
      end
    else
    end
  end
end

terra particles_pullElement(src : &int8) : particles_columns
    var dst : particles_columns
    C.memcpy([&opaque](&dst), [&opaque](src), [uint64](376))
    return dst
end

task particles_pullAll(v8259 : int3d, v8260 : region(ispace(int1d), particles_columns), v8205 : region(ispace(int1d), int8[376]), v8207 : region(ispace(int1d), int8[376]), v8209 : region(ispace(int1d), int8[376]), v8211 : region(ispace(int1d), int8[376]), v8213 : region(ispace(int1d), int8[376]), v8215 : region(ispace(int1d), int8[376]), v8217 : region(ispace(int1d), int8[376]), v8219 : region(ispace(int1d), int8[376]), v8221 : region(ispace(int1d), int8[376]), v8223 : region(ispace(int1d), int8[376]), v8225 : region(ispace(int1d), int8[376]), v8227 : region(ispace(int1d), int8[376]), v8229 : region(ispace(int1d), int8[376]), v8231 : region(ispace(int1d), int8[376]), v8233 : region(ispace(int1d), int8[376]), v8235 : region(ispace(int1d), int8[376]), v8237 : region(ispace(int1d), int8[376]), v8239 : region(ispace(int1d), int8[376]), v8241 : region(ispace(int1d), int8[376]), v8243 : region(ispace(int1d), int8[376]), v8245 : region(ispace(int1d), int8[376]), v8247 : region(ispace(int1d), int8[376]), v8249 : region(ispace(int1d), int8[376]), v8251 : region(ispace(int1d), int8[376]), v8253 : region(ispace(int1d), int8[376]), v8255 : region(ispace(int1d), int8[376]))
where
  reads(v8260), writes(v8260), reads(v8205), reads(v8207), reads(v8209), reads(v8211), reads(v8213), reads(v8215), reads(v8217), reads(v8219), reads(v8221), reads(v8223), reads(v8225), reads(v8227), reads(v8229), reads(v8231), reads(v8233), reads(v8235), reads(v8237), reads(v8239), reads(v8241), reads(v8243), reads(v8245), reads(v8247), reads(v8249), reads(v8251), reads(v8253), reads(v8255)
do
  for v8418 in v8205 do
    if bool(v8205[v8418][368LL]) then
      var v8419 = false
      for v8420 in v8260 do
        if (not v8420.__valid) then
          v8260[v8420] = particles_pullElement([&int8](v8205[v8418]))
          v8419 = true
          regentlib.assert(v8260[v8420].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8419, "Ran out of space on sub-partition")
    else
    end
  end
  for v8421 in v8207 do
    if bool(v8207[v8421][368LL]) then
      var v8422 = false
      for v8423 in v8260 do
        if (not v8423.__valid) then
          v8260[v8423] = particles_pullElement([&int8](v8207[v8421]))
          v8422 = true
          regentlib.assert(v8260[v8423].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8422, "Ran out of space on sub-partition")
    else
    end
  end
  for v8424 in v8209 do
    if bool(v8209[v8424][368LL]) then
      var v8425 = false
      for v8426 in v8260 do
        if (not v8426.__valid) then
          v8260[v8426] = particles_pullElement([&int8](v8209[v8424]))
          v8425 = true
          regentlib.assert(v8260[v8426].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8425, "Ran out of space on sub-partition")
    else
    end
  end
  for v8427 in v8211 do
    if bool(v8211[v8427][368LL]) then
      var v8428 = false
      for v8429 in v8260 do
        if (not v8429.__valid) then
          v8260[v8429] = particles_pullElement([&int8](v8211[v8427]))
          v8428 = true
          regentlib.assert(v8260[v8429].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8428, "Ran out of space on sub-partition")
    else
    end
  end
  for v8430 in v8213 do
    if bool(v8213[v8430][368LL]) then
      var v8431 = false
      for v8432 in v8260 do
        if (not v8432.__valid) then
          v8260[v8432] = particles_pullElement([&int8](v8213[v8430]))
          v8431 = true
          regentlib.assert(v8260[v8432].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8431, "Ran out of space on sub-partition")
    else
    end
  end
  for v8433 in v8215 do
    if bool(v8215[v8433][368LL]) then
      var v8434 = false
      for v8435 in v8260 do
        if (not v8435.__valid) then
          v8260[v8435] = particles_pullElement([&int8](v8215[v8433]))
          v8434 = true
          regentlib.assert(v8260[v8435].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8434, "Ran out of space on sub-partition")
    else
    end
  end
  for v8436 in v8217 do
    if bool(v8217[v8436][368LL]) then
      var v8437 = false
      for v8438 in v8260 do
        if (not v8438.__valid) then
          v8260[v8438] = particles_pullElement([&int8](v8217[v8436]))
          v8437 = true
          regentlib.assert(v8260[v8438].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8437, "Ran out of space on sub-partition")
    else
    end
  end
  for v8439 in v8219 do
    if bool(v8219[v8439][368LL]) then
      var v8440 = false
      for v8441 in v8260 do
        if (not v8441.__valid) then
          v8260[v8441] = particles_pullElement([&int8](v8219[v8439]))
          v8440 = true
          regentlib.assert(v8260[v8441].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8440, "Ran out of space on sub-partition")
    else
    end
  end
  for v8442 in v8221 do
    if bool(v8221[v8442][368LL]) then
      var v8443 = false
      for v8444 in v8260 do
        if (not v8444.__valid) then
          v8260[v8444] = particles_pullElement([&int8](v8221[v8442]))
          v8443 = true
          regentlib.assert(v8260[v8444].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8443, "Ran out of space on sub-partition")
    else
    end
  end
  for v8445 in v8223 do
    if bool(v8223[v8445][368LL]) then
      var v8446 = false
      for v8447 in v8260 do
        if (not v8447.__valid) then
          v8260[v8447] = particles_pullElement([&int8](v8223[v8445]))
          v8446 = true
          regentlib.assert(v8260[v8447].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8446, "Ran out of space on sub-partition")
    else
    end
  end
  for v8448 in v8225 do
    if bool(v8225[v8448][368LL]) then
      var v8449 = false
      for v8450 in v8260 do
        if (not v8450.__valid) then
          v8260[v8450] = particles_pullElement([&int8](v8225[v8448]))
          v8449 = true
          regentlib.assert(v8260[v8450].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8449, "Ran out of space on sub-partition")
    else
    end
  end
  for v8451 in v8227 do
    if bool(v8227[v8451][368LL]) then
      var v8452 = false
      for v8453 in v8260 do
        if (not v8453.__valid) then
          v8260[v8453] = particles_pullElement([&int8](v8227[v8451]))
          v8452 = true
          regentlib.assert(v8260[v8453].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8452, "Ran out of space on sub-partition")
    else
    end
  end
  for v8454 in v8229 do
    if bool(v8229[v8454][368LL]) then
      var v8455 = false
      for v8456 in v8260 do
        if (not v8456.__valid) then
          v8260[v8456] = particles_pullElement([&int8](v8229[v8454]))
          v8455 = true
          regentlib.assert(v8260[v8456].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8455, "Ran out of space on sub-partition")
    else
    end
  end
  for v8457 in v8231 do
    if bool(v8231[v8457][368LL]) then
      var v8458 = false
      for v8459 in v8260 do
        if (not v8459.__valid) then
          v8260[v8459] = particles_pullElement([&int8](v8231[v8457]))
          v8458 = true
          regentlib.assert(v8260[v8459].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8458, "Ran out of space on sub-partition")
    else
    end
  end
  for v8460 in v8233 do
    if bool(v8233[v8460][368LL]) then
      var v8461 = false
      for v8462 in v8260 do
        if (not v8462.__valid) then
          v8260[v8462] = particles_pullElement([&int8](v8233[v8460]))
          v8461 = true
          regentlib.assert(v8260[v8462].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8461, "Ran out of space on sub-partition")
    else
    end
  end
  for v8463 in v8235 do
    if bool(v8235[v8463][368LL]) then
      var v8464 = false
      for v8465 in v8260 do
        if (not v8465.__valid) then
          v8260[v8465] = particles_pullElement([&int8](v8235[v8463]))
          v8464 = true
          regentlib.assert(v8260[v8465].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8464, "Ran out of space on sub-partition")
    else
    end
  end
  for v8466 in v8237 do
    if bool(v8237[v8466][368LL]) then
      var v8467 = false
      for v8468 in v8260 do
        if (not v8468.__valid) then
          v8260[v8468] = particles_pullElement([&int8](v8237[v8466]))
          v8467 = true
          regentlib.assert(v8260[v8468].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8467, "Ran out of space on sub-partition")
    else
    end
  end
  for v8469 in v8239 do
    if bool(v8239[v8469][368LL]) then
      var v8470 = false
      for v8471 in v8260 do
        if (not v8471.__valid) then
          v8260[v8471] = particles_pullElement([&int8](v8239[v8469]))
          v8470 = true
          regentlib.assert(v8260[v8471].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8470, "Ran out of space on sub-partition")
    else
    end
  end
  for v8472 in v8241 do
    if bool(v8241[v8472][368LL]) then
      var v8473 = false
      for v8474 in v8260 do
        if (not v8474.__valid) then
          v8260[v8474] = particles_pullElement([&int8](v8241[v8472]))
          v8473 = true
          regentlib.assert(v8260[v8474].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8473, "Ran out of space on sub-partition")
    else
    end
  end
  for v8475 in v8243 do
    if bool(v8243[v8475][368LL]) then
      var v8476 = false
      for v8477 in v8260 do
        if (not v8477.__valid) then
          v8260[v8477] = particles_pullElement([&int8](v8243[v8475]))
          v8476 = true
          regentlib.assert(v8260[v8477].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8476, "Ran out of space on sub-partition")
    else
    end
  end
  for v8478 in v8245 do
    if bool(v8245[v8478][368LL]) then
      var v8479 = false
      for v8480 in v8260 do
        if (not v8480.__valid) then
          v8260[v8480] = particles_pullElement([&int8](v8245[v8478]))
          v8479 = true
          regentlib.assert(v8260[v8480].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8479, "Ran out of space on sub-partition")
    else
    end
  end
  for v8481 in v8247 do
    if bool(v8247[v8481][368LL]) then
      var v8482 = false
      for v8483 in v8260 do
        if (not v8483.__valid) then
          v8260[v8483] = particles_pullElement([&int8](v8247[v8481]))
          v8482 = true
          regentlib.assert(v8260[v8483].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8482, "Ran out of space on sub-partition")
    else
    end
  end
  for v8484 in v8249 do
    if bool(v8249[v8484][368LL]) then
      var v8485 = false
      for v8486 in v8260 do
        if (not v8486.__valid) then
          v8260[v8486] = particles_pullElement([&int8](v8249[v8484]))
          v8485 = true
          regentlib.assert(v8260[v8486].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8485, "Ran out of space on sub-partition")
    else
    end
  end
  for v8487 in v8251 do
    if bool(v8251[v8487][368LL]) then
      var v8488 = false
      for v8489 in v8260 do
        if (not v8489.__valid) then
          v8260[v8489] = particles_pullElement([&int8](v8251[v8487]))
          v8488 = true
          regentlib.assert(v8260[v8489].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8488, "Ran out of space on sub-partition")
    else
    end
  end
  for v8490 in v8253 do
    if bool(v8253[v8490][368LL]) then
      var v8491 = false
      for v8492 in v8260 do
        if (not v8492.__valid) then
          v8260[v8492] = particles_pullElement([&int8](v8253[v8490]))
          v8491 = true
          regentlib.assert(v8260[v8492].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8491, "Ran out of space on sub-partition")
    else
    end
  end
  for v8493 in v8255 do
    if bool(v8255[v8493][368LL]) then
      var v8494 = false
      for v8495 in v8260 do
        if (not v8495.__valid) then
          v8260[v8495] = particles_pullElement([&int8](v8255[v8493]))
          v8494 = true
          regentlib.assert(v8260[v8495].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      regentlib.assert(v8494, "Ran out of space on sub-partition")
    else
    end
  end
end

__demand(__parallel)
task Particles_AddFlowCoupling(v8816 : region(ispace(int1d), particles_columns), v8819 : region(ispace(int3d), Fluid_columns), v8820 : double, v8821 : double, v8822 : double, v8823 : double, v8824 : double, v8825 : double, v8826 : int32, v8827 : double, v8828 : double, v8829 : double, v8830 : double, v8831 : double, v8832 : double, v8833 : double, v8834 : double)
where
  reads(v8819.centerCoordinates), reads(v8819.temperature), reads(v8819.velocity), reads(v8816.cell), reads(v8816.deltaTemperatureTerm), writes(v8816.deltaTemperatureTerm), reads(v8816.deltaVelocityOverRelaxationTime), writes(v8816.deltaVelocityOverRelaxationTime), reads(v8816.density), reads(v8816.diameter), reads(v8816.position), reads(v8816.position_t), writes(v8816.position_t), reads(v8816.temperature), reads(v8816.temperature_t), writes(v8816.temperature_t), reads(v8816.velocity), reads(v8816.velocity_t), writes(v8816.velocity_t), reads(v8816.__valid)
do
  for v9479 in v8816 do
    if v8816[v9479].__valid then
      var v9532 : double[3]
      do
        var v9524 = v8816[v9479].cell
        var v9525 = v8816[v9479].position
        var v9526 = v8827
        var v9527 = v8828
        var v9528 = v8829
        var v9529 = v8830
        var v9530 = v8831
        var v9531 = v8832
        var v9082 = [double[3]](array(double(0), double(0), double(0)))
        var v9083 = [double[3]](array(double(0), double(0), double(0)))
        var v9084 = [double[3]](array(double(0), double(0), double(0)))
        var v9085 = [double[3]](array(double(0), double(0), double(0)))
        var v9086 = [double[3]](array(double(0), double(0), double(0)))
        var v9087 = [double[3]](array(double(0), double(0), double(0)))
        var v9088 = [double[3]](array(double(0), double(0), double(0)))
        var v9089 = [double[3]](array(double(0), double(0), double(0)))
        var v9090 = v8819[v9524].velocity
        if (v9525[int32(0)]>v8819[v9524].centerCoordinates[int32(0)]) then
          var v9091 = v8819[((v9524+{1, 0, 0})%v8819.bounds)].velocity
          if (v9525[int32(1)]>v8819[v9524].centerCoordinates[int32(1)]) then
            var v9092 = v8819[((v9524+{0, 1, 0})%v8819.bounds)].velocity
            var v9093 = v8819[((v9524+{1, 1, 0})%v8819.bounds)].velocity
            if (v9525[int32(2)]>v8819[v9524].centerCoordinates[int32(2)]) then
              v9082 = v9090
              v9083 = v9091
              v9084 = v9092
              v9085 = v9093
              v9086 = v8819[((v9524+{0, 0, 1})%v8819.bounds)].velocity
              v9087 = v8819[((v9524+{1, 0, 1})%v8819.bounds)].velocity
              v9088 = v8819[((v9524+{0, 1, 1})%v8819.bounds)].velocity
              v9089 = v8819[((v9524+{1, 1, 1})%v8819.bounds)].velocity
            else
              v9082 = v8819[((v9524+{0, 0, -1})%v8819.bounds)].velocity
              v9083 = v8819[((v9524+{1, 0, -1})%v8819.bounds)].velocity
              v9084 = v8819[((v9524+{0, 1, -1})%v8819.bounds)].velocity
              v9085 = v8819[((v9524+{1, 1, -1})%v8819.bounds)].velocity
              v9086 = v9090
              v9087 = v9091
              v9088 = v9092
              v9089 = v9093
            end
          else
            var v9094 = v8819[((v9524+{0, -1, 0})%v8819.bounds)].velocity
            var v9095 = v8819[((v9524+{1, -1, 0})%v8819.bounds)].velocity
            if (v9525[int32(2)]>v8819[v9524].centerCoordinates[int32(2)]) then
              v9082 = v9094
              v9083 = v9095
              v9084 = v9090
              v9085 = v9091
              v9086 = v8819[((v9524+{0, -1, 1})%v8819.bounds)].velocity
              v9087 = v8819[((v9524+{1, -1, 1})%v8819.bounds)].velocity
              v9088 = v8819[((v9524+{0, 0, 1})%v8819.bounds)].velocity
              v9089 = v8819[((v9524+{1, 0, 1})%v8819.bounds)].velocity
            else
              v9082 = v8819[((v9524+{0, -1, -1})%v8819.bounds)].velocity
              v9083 = v8819[((v9524+{1, -1, -1})%v8819.bounds)].velocity
              v9084 = v8819[((v9524+{0, 0, -1})%v8819.bounds)].velocity
              v9085 = v8819[((v9524+{1, 0, -1})%v8819.bounds)].velocity
              v9086 = v9094
              v9087 = v9095
              v9088 = v9090
              v9089 = v9091
            end
          end
        else
          var v9096 = v8819[((v9524+{-1, 0, 0})%v8819.bounds)].velocity
          if (v9525[int32(1)]>v8819[v9524].centerCoordinates[int32(1)]) then
            var v9097 = v8819[((v9524+{-1, 1, 0})%v8819.bounds)].velocity
            var v9098 = v8819[((v9524+{0, 1, 0})%v8819.bounds)].velocity
            if (v9525[int32(2)]>v8819[v9524].centerCoordinates[int32(2)]) then
              v9082 = v9096
              v9083 = v9090
              v9084 = v9097
              v9085 = v9098
              v9086 = v8819[((v9524+{-1, 0, 1})%v8819.bounds)].velocity
              v9087 = v8819[((v9524+{0, 0, 1})%v8819.bounds)].velocity
              v9088 = v8819[((v9524+{-1, 1, 1})%v8819.bounds)].velocity
              v9089 = v8819[((v9524+{0, 1, 1})%v8819.bounds)].velocity
            else
              v9082 = v8819[((v9524+{-1, 0, -1})%v8819.bounds)].velocity
              v9083 = v8819[((v9524+{0, 0, -1})%v8819.bounds)].velocity
              v9084 = v8819[((v9524+{-1, 1, -1})%v8819.bounds)].velocity
              v9085 = v8819[((v9524+{0, 1, -1})%v8819.bounds)].velocity
              v9086 = v9096
              v9087 = v9090
              v9088 = v9097
              v9089 = v9098
            end
          else
            var v9099 = v8819[((v9524+{-1, -1, 0})%v8819.bounds)].velocity
            var v9100 = v8819[((v9524+{0, -1, 0})%v8819.bounds)].velocity
            if (v9525[int32(2)]>v8819[v9524].centerCoordinates[int32(2)]) then
              v9082 = v9099
              v9083 = v9100
              v9084 = v9096
              v9085 = v9090
              v9086 = v8819[((v9524+{-1, -1, 1})%v8819.bounds)].velocity
              v9087 = v8819[((v9524+{0, -1, 1})%v8819.bounds)].velocity
              v9088 = v8819[((v9524+{-1, 0, 1})%v8819.bounds)].velocity
              v9089 = v8819[((v9524+{0, 0, 1})%v8819.bounds)].velocity
            else
              v9082 = v8819[((v9524+{-1, -1, -1})%v8819.bounds)].velocity
              v9083 = v8819[((v9524+{0, -1, -1})%v8819.bounds)].velocity
              v9084 = v8819[((v9524+{-1, 0, -1})%v8819.bounds)].velocity
              v9085 = v8819[((v9524+{0, 0, -1})%v8819.bounds)].velocity
              v9086 = v9099
              v9087 = v9100
              v9088 = v9096
              v9089 = v9090
            end
          end
        end
        var v9183 : double[3]
        do
          var v9168 = v9525
          var v9169 = v9082
          var v9170 = v9083
          var v9171 = v9084
          var v9172 = v9085
          var v9173 = v9086
          var v9174 = v9087
          var v9175 = v9088
          var v9176 = v9089
          var v9177 = v9526
          var v9178 = v9527
          var v9179 = v9528
          var v9180 = v9529
          var v9181 = v9530
          var v9182 = v9531
          var v8990 = C.fmod((((v9168[int32(0)]-v9178)/v9177)+double(0.5)), double(1))
          var v8991 = C.fmod((((v9168[int32(1)]-v9180)/v9179)+double(0.5)), double(1))
          var v8992 = C.fmod((((v9168[int32(2)]-v9182)/v9181)+double(0.5)), double(1))
          var v8993 = (double(1)-v8990)
          var v8994 = (double(1)-v8991)
          var v8995 = (double(1)-v8992)
          var v8996 = vv_add_double_3(vs_mul_double_3(v9169, v8993), vs_mul_double_3(v9170, v8990))
          var v8997 = vv_add_double_3(vs_mul_double_3(v9171, v8993), vs_mul_double_3(v9172, v8990))
          var v8998 = vv_add_double_3(vs_mul_double_3(v9173, v8993), vs_mul_double_3(v9174, v8990))
          var v8999 = vv_add_double_3(vs_mul_double_3(v9175, v8993), vs_mul_double_3(v9176, v8990))
          var v9000 = vv_add_double_3(vs_mul_double_3(v8996, v8994), vs_mul_double_3(v8997, v8991))
          var v9001 = vv_add_double_3(vs_mul_double_3(v8998, v8994), vs_mul_double_3(v8999, v8991))
          v9183 = vv_add_double_3(vs_mul_double_3(v9000, v8995), vs_mul_double_3(v9001, v8992))
        end
        var v9184 = 0
        v9532 = v9183
      end
      var v9533 = 0
      var v9480 = v9532
      var v9542 : double
      do
        var v9534 = v8816[v9479].cell
        var v9535 = v8816[v9479].position
        var v9536 = v8827
        var v9537 = v8828
        var v9538 = v8829
        var v9539 = v8830
        var v9540 = v8831
        var v9541 = v8832
        var v9377 = double(double(0))
        var v9378 = double(double(0))
        var v9379 = double(double(0))
        var v9380 = double(double(0))
        var v9381 = double(double(0))
        var v9382 = double(double(0))
        var v9383 = double(double(0))
        var v9384 = double(double(0))
        var v9385 = v8819[v9534].temperature
        if (v9535[int32(0)]>v8819[v9534].centerCoordinates[int32(0)]) then
          var v9386 = v8819[((v9534+{1, 0, 0})%v8819.bounds)].temperature
          if (v9535[int32(1)]>v8819[v9534].centerCoordinates[int32(1)]) then
            var v9387 = v8819[((v9534+{0, 1, 0})%v8819.bounds)].temperature
            var v9388 = v8819[((v9534+{1, 1, 0})%v8819.bounds)].temperature
            if (v9535[int32(2)]>v8819[v9534].centerCoordinates[int32(2)]) then
              v9377 = v9385
              v9378 = v9386
              v9379 = v9387
              v9380 = v9388
              v9381 = v8819[((v9534+{0, 0, 1})%v8819.bounds)].temperature
              v9382 = v8819[((v9534+{1, 0, 1})%v8819.bounds)].temperature
              v9383 = v8819[((v9534+{0, 1, 1})%v8819.bounds)].temperature
              v9384 = v8819[((v9534+{1, 1, 1})%v8819.bounds)].temperature
            else
              v9377 = v8819[((v9534+{0, 0, -1})%v8819.bounds)].temperature
              v9378 = v8819[((v9534+{1, 0, -1})%v8819.bounds)].temperature
              v9379 = v8819[((v9534+{0, 1, -1})%v8819.bounds)].temperature
              v9380 = v8819[((v9534+{1, 1, -1})%v8819.bounds)].temperature
              v9381 = v9385
              v9382 = v9386
              v9383 = v9387
              v9384 = v9388
            end
          else
            var v9389 = v8819[((v9534+{0, -1, 0})%v8819.bounds)].temperature
            var v9390 = v8819[((v9534+{1, -1, 0})%v8819.bounds)].temperature
            if (v9535[int32(2)]>v8819[v9534].centerCoordinates[int32(2)]) then
              v9377 = v9389
              v9378 = v9390
              v9379 = v9385
              v9380 = v9386
              v9381 = v8819[((v9534+{0, -1, 1})%v8819.bounds)].temperature
              v9382 = v8819[((v9534+{1, -1, 1})%v8819.bounds)].temperature
              v9383 = v8819[((v9534+{0, 0, 1})%v8819.bounds)].temperature
              v9384 = v8819[((v9534+{1, 0, 1})%v8819.bounds)].temperature
            else
              v9377 = v8819[((v9534+{0, -1, -1})%v8819.bounds)].temperature
              v9378 = v8819[((v9534+{1, -1, -1})%v8819.bounds)].temperature
              v9379 = v8819[((v9534+{0, 0, -1})%v8819.bounds)].temperature
              v9380 = v8819[((v9534+{1, 0, -1})%v8819.bounds)].temperature
              v9381 = v9389
              v9382 = v9390
              v9383 = v9385
              v9384 = v9386
            end
          end
        else
          var v9391 = v8819[((v9534+{-1, 0, 0})%v8819.bounds)].temperature
          if (v9535[int32(1)]>v8819[v9534].centerCoordinates[int32(1)]) then
            var v9392 = v8819[((v9534+{-1, 1, 0})%v8819.bounds)].temperature
            var v9393 = v8819[((v9534+{0, 1, 0})%v8819.bounds)].temperature
            if (v9535[int32(2)]>v8819[v9534].centerCoordinates[int32(2)]) then
              v9377 = v9391
              v9378 = v9385
              v9379 = v9392
              v9380 = v9393
              v9381 = v8819[((v9534+{-1, 0, 1})%v8819.bounds)].temperature
              v9382 = v8819[((v9534+{0, 0, 1})%v8819.bounds)].temperature
              v9383 = v8819[((v9534+{-1, 1, 1})%v8819.bounds)].temperature
              v9384 = v8819[((v9534+{0, 1, 1})%v8819.bounds)].temperature
            else
              v9377 = v8819[((v9534+{-1, 0, -1})%v8819.bounds)].temperature
              v9378 = v8819[((v9534+{0, 0, -1})%v8819.bounds)].temperature
              v9379 = v8819[((v9534+{-1, 1, -1})%v8819.bounds)].temperature
              v9380 = v8819[((v9534+{0, 1, -1})%v8819.bounds)].temperature
              v9381 = v9391
              v9382 = v9385
              v9383 = v9392
              v9384 = v9393
            end
          else
            var v9394 = v8819[((v9534+{-1, -1, 0})%v8819.bounds)].temperature
            var v9395 = v8819[((v9534+{0, -1, 0})%v8819.bounds)].temperature
            if (v9535[int32(2)]>v8819[v9534].centerCoordinates[int32(2)]) then
              v9377 = v9394
              v9378 = v9395
              v9379 = v9391
              v9380 = v9385
              v9381 = v8819[((v9534+{-1, -1, 1})%v8819.bounds)].temperature
              v9382 = v8819[((v9534+{0, -1, 1})%v8819.bounds)].temperature
              v9383 = v8819[((v9534+{-1, 0, 1})%v8819.bounds)].temperature
              v9384 = v8819[((v9534+{0, 0, 1})%v8819.bounds)].temperature
            else
              v9377 = v8819[((v9534+{-1, -1, -1})%v8819.bounds)].temperature
              v9378 = v8819[((v9534+{0, -1, -1})%v8819.bounds)].temperature
              v9379 = v8819[((v9534+{-1, 0, -1})%v8819.bounds)].temperature
              v9380 = v8819[((v9534+{0, 0, -1})%v8819.bounds)].temperature
              v9381 = v9394
              v9382 = v9395
              v9383 = v9391
              v9384 = v9385
            end
          end
        end
        var v9430 : double
        do
          var v9415 = v9535
          var v9416 = v9377
          var v9417 = v9378
          var v9418 = v9379
          var v9419 = v9380
          var v9420 = v9381
          var v9421 = v9382
          var v9422 = v9383
          var v9423 = v9384
          var v9424 = v9536
          var v9425 = v9537
          var v9426 = v9538
          var v9427 = v9539
          var v9428 = v9540
          var v9429 = v9541
          var v9341 = C.fmod((((v9415[int32(0)]-v9425)/v9424)+double(0.5)), double(1))
          var v9342 = C.fmod((((v9415[int32(1)]-v9427)/v9426)+double(0.5)), double(1))
          var v9343 = C.fmod((((v9415[int32(2)]-v9429)/v9428)+double(0.5)), double(1))
          var v9344 = (double(1)-v9341)
          var v9345 = (double(1)-v9342)
          var v9346 = (double(1)-v9343)
          var v9347 = ((v9416*v9344)+(v9417*v9341))
          var v9348 = ((v9418*v9344)+(v9419*v9341))
          var v9349 = ((v9420*v9344)+(v9421*v9341))
          var v9350 = ((v9422*v9344)+(v9423*v9341))
          var v9351 = ((v9347*v9345)+(v9348*v9342))
          var v9352 = ((v9349*v9345)+(v9350*v9342))
          v9430 = ((v9351*v9346)+(v9352*v9343))
        end
        var v9431 = 0
        v9542 = v9430
      end
      var v9543 = 0
      var v9481 = v9542
      var v9552 : double
      do
        var v9544 = v9481
        var v9545 = v8820
        var v9546 = v8821
        var v9547 = v8822
        var v9548 = v8823
        var v9549 = v8824
        var v9550 = v8825
        var v9551 = v8826
        var v4253 = double(double(0))
        if (v9551==int32(0)) then
          v4253 = v9545
        else
          if (v9551==int32(1)) then
            v4253 = (v9547*C.pow((v9544/v9546), double(0.75)))
          else
            if (v9551==int32(2)) then
              v4253 = ((v9550*C.pow((v9544/v9549), (double(3)/double(2))))*((v9549+v9548)/(v9544+v9548)))
            else
              regentlib.assert(false, "(Liszt assertion)")
            end
          end
        end
        v9552 = v4253
      end
      var v9553 = 0
      var v9482 = v9552
      var v9483 = v8816[v9479].velocity
      var v9484 = v8816[v9479].position_t
      v9484[0] += v9483[0]
      v9484[1] += v9483[1]
      v9484[2] += v9483[2]
      v8816[v9479].position_t = v9484
      var v9485 = double(0)
      var v9486 = (((v8816[v9479].density*C.pow(v8816[v9479].diameter, double(int32(2))))/(double(18)*v9482))/(double(1)+(double(0.15)*C.pow(v9485, double(0.687)))))
      v8816[v9479].deltaVelocityOverRelaxationTime = vs_div_double_3(vv_sub_double_3(v9480, v8816[v9479].velocity), v9486)
      v8816[v9479].deltaTemperatureTerm = (((double(3.1415926535898)*C.pow(v8816[v9479].diameter, double(int32(2))))*v8833)*(v9481-v8816[v9479].temperature))
      var v9487 = v8816[v9479].deltaVelocityOverRelaxationTime
      var v9488 = v8816[v9479].velocity_t
      v9488[0] += v9487[0]
      v9488[1] += v9487[1]
      v9488[2] += v9487[2]
      v8816[v9479].velocity_t = v9488
      v8816[v9479].temperature_t += (v8816[v9479].deltaTemperatureTerm/((((double(3.1415926535898)*C.pow(v8816[v9479].diameter, double(int32(3))))/double(6))*v8816[v9479].density)*v8834))
    else
    end
  end
end

__demand(__parallel)
task Particles_AddBodyForces(v9953 : region(ispace(int1d), particles_columns), v9955 : double[3])
where
  reads(v9953.velocity_t), writes(v9953.velocity_t), reads(v9953.__valid)
do
  for v9969 in v9953 do
    if v9953[v9969].__valid then
      var v9970 = v9955
      var v9971 = v9953[v9969].velocity_t
      v9971[0] += v9970[0]
      v9971[1] += v9970[1]
      v9971[2] += v9970[2]
      v9953[v9969].velocity_t = v9971
    else
    end
  end
end


__demand(__parallel)
task Radiation_ClearAccumulators(vRadiation : region(ispace(int3d), Radiation_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads(vRadiation.acc_d2), writes(vRadiation.acc_d2), reads(vRadiation.acc_d2t4), writes(vRadiation.acc_d2t4)
do
  for vc in vRadiation do
    vRadiation[vc].acc_d2 = double(0)
    vRadiation[vc].acc_d2t4 = double(0)
  end
end

task Radiation_AccumulateParticleValues(vparticles : region(ispace(int1d), particles_columns), vFluid : region(ispace(int3d), Fluid_columns), vRadiation : region(ispace(int3d), Radiation_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads(vFluid.to_Radiation), reads(vRadiation.acc_d2), writes(vRadiation.acc_d2), reads(vRadiation.acc_d2t4), writes(vRadiation.acc_d2t4), reads(vparticles.cell), reads(vparticles.diameter), reads(vparticles.temperature), reads(vparticles.__valid)
do
  for vp in vparticles do
    if vparticles[vp].__valid then
      vRadiation[vFluid[vparticles[vp].cell].to_Radiation].acc_d2 += C.pow(vparticles[vp].diameter, double(2))
      vRadiation[vFluid[vparticles[vp].cell].to_Radiation].acc_d2t4 += (C.pow(vparticles[vp].diameter, double(2))*C.pow(vparticles[vp].temperature, double(4)))
    else
    end
  end
end

__demand(__parallel)
task Radiation_UpdateFieldValues(vRadiation : region(ispace(int3d), Radiation_columns), vRadiation_cellVolume : double, vRadiation_qa : double, vRadiation_qs : double)
-- leaf (false), inner (false), idempotent (false)
where
  reads(vRadiation.Ib), writes(vRadiation.Ib), reads(vRadiation.acc_d2), reads(vRadiation.acc_d2t4), reads(vRadiation.sigma), writes(vRadiation.sigma)
do
  for vc in vRadiation do
    vRadiation[vc].sigma = (((vRadiation[vc].acc_d2*double(3.1415926535898))*(vRadiation_qa+vRadiation_qs))/(double(4)*vRadiation_cellVolume))
    if (vRadiation[vc].acc_d2==double(0)) then
      vRadiation[vc].Ib = double(0)
    else
      vRadiation[vc].Ib = ((double(5.67e-08)*vRadiation[vc].acc_d2t4)/(double(3.1415926535898)*vRadiation[vc].acc_d2))
    end
  end
end

task Particles_AbsorbRadiation(vparticles : region(ispace(int1d), particles_columns), vFluid : region(ispace(int3d), Fluid_columns), vRadiation : region(ispace(int3d), Radiation_columns), vParticles_heatCapacity : double, vRadiation_qa : double)
-- leaf (false), inner (false), idempotent (false)
where
  reads(vFluid.to_Radiation), reads(vRadiation.G), reads(vparticles.cell), reads(vparticles.density), reads(vparticles.diameter), reads(vparticles.temperature), reads(vparticles.temperature_t), writes(vparticles.temperature_t), reads(vparticles.__valid)
do
  for vp in vparticles do
    if vparticles[vp].__valid then
      var vt4 = C.pow(vparticles[vp].temperature, double(4))
      var valpha = ((((double(3.1415926535898)*vRadiation_qa)*C.pow(vparticles[vp].diameter, double(2)))*(vRadiation[vFluid[vparticles[vp].cell].to_Radiation].G-((double(4)*double(5.67e-08))*vt4)))/double(4))
      vparticles[vp].temperature_t += (valpha/((((double(3.1415926535898)*C.pow(vparticles[vp].diameter, double(int32(3))))/double(6))*vparticles[vp].density)*vParticles_heatCapacity))
    else
    end
  end
end

-- XXX: Turn off two-way coupling
-- task Flow_AddParticlesCoupling(v9979 : region(ispace(int1d), particles_columns), v9982 : region(ispace(int3d), Fluid_columns), v9983 : double)
-- where
--   reads(v9982.rhoEnergy_t), writes(v9982.rhoEnergy_t), reads(v9982.rhoVelocity_t), writes(v9982.rhoVelocity_t), reads(v9979.cell), reads(v9979.deltaTemperatureTerm), reads(v9979.deltaVelocityOverRelaxationTime), reads(v9979.density), reads(v9979.diameter), reads(v9979.__valid)
-- do
--   for v9997 in v9979 do
--     if v9979[v9997].__valid then
--       var v9998 = vs_div_double_3(vs_mul_double_3(v9979[v9997].deltaVelocityOverRelaxationTime, (-(((double(3.1415926535898)*C.pow(v9979[v9997].diameter, double(int32(3))))/double(6))*v9979[v9997].density))), v9983)
--       var v9999 = v9982[v9979[v9997].cell].rhoVelocity_t
--       v9999[0] += v9998[0]
--       v9999[1] += v9998[1]
--       v9999[2] += v9998[2]
--       v9982[v9979[v9997].cell].rhoVelocity_t = v9999
--       v9982[v9979[v9997].cell].rhoEnergy_t += ((-v9979[v9997].deltaTemperatureTerm)/v9983)
--     else
--     end
--   end
-- end

__demand(__parallel)
task Flow_UpdateVars(v10022 : region(ispace(int3d), Fluid_columns), v10024 : double, v10025 : int32)
where
  reads(v10022.rho), writes(v10022.rho), reads(v10022.rhoEnergy), writes(v10022.rhoEnergy), reads(v10022.rhoEnergy_new), reads(v10022.rhoEnergy_new), writes(v10022.rhoEnergy_new), reads(v10022.rhoEnergy_old), reads(v10022.rhoEnergy_t), reads(v10022.rhoVelocity), writes(v10022.rhoVelocity), reads(v10022.rhoVelocity_new), reads(v10022.rhoVelocity_new), writes(v10022.rhoVelocity_new), reads(v10022.rhoVelocity_old), reads(v10022.rhoVelocity_t), reads(v10022.rho_new), reads(v10022.rho_new), writes(v10022.rho_new), reads(v10022.rho_old), reads(v10022.rho_t)
do
  for v10079 in v10022 do
    var v10080 = v10024
    if (v10025==int32(1)) then
      v10022[v10079].rho_new += (((double(1)/double(6))*v10080)*v10022[v10079].rho_t)
      v10022[v10079].rho = (v10022[v10079].rho_old+((double(0.5)*v10080)*v10022[v10079].rho_t))
      var v10081 = vs_mul_double_3(v10022[v10079].rhoVelocity_t, ((double(1)/double(6))*v10080))
      var v10082 = v10022[v10079].rhoVelocity_new
      v10082[0] += v10081[0]
      v10082[1] += v10081[1]
      v10082[2] += v10081[2]
      v10022[v10079].rhoVelocity_new = v10082
      v10022[v10079].rhoVelocity = vv_add_double_3(v10022[v10079].rhoVelocity_old, vs_mul_double_3(v10022[v10079].rhoVelocity_t, (double(0.5)*v10080)))
      v10022[v10079].rhoEnergy_new += (((double(1)/double(6))*v10080)*v10022[v10079].rhoEnergy_t)
      v10022[v10079].rhoEnergy = (v10022[v10079].rhoEnergy_old+((double(0.5)*v10080)*v10022[v10079].rhoEnergy_t))
    else
      if (v10025==int32(2)) then
        v10022[v10079].rho_new += (((double(1)/double(3))*v10080)*v10022[v10079].rho_t)
        v10022[v10079].rho = (v10022[v10079].rho_old+((double(0.5)*v10080)*v10022[v10079].rho_t))
        var v10083 = vs_mul_double_3(v10022[v10079].rhoVelocity_t, ((double(1)/double(3))*v10080))
        var v10084 = v10022[v10079].rhoVelocity_new
        v10084[0] += v10083[0]
        v10084[1] += v10083[1]
        v10084[2] += v10083[2]
        v10022[v10079].rhoVelocity_new = v10084
        v10022[v10079].rhoVelocity = vv_add_double_3(v10022[v10079].rhoVelocity_old, vs_mul_double_3(v10022[v10079].rhoVelocity_t, (double(0.5)*v10080)))
        v10022[v10079].rhoEnergy_new += (((double(1)/double(3))*v10080)*v10022[v10079].rhoEnergy_t)
        v10022[v10079].rhoEnergy = (v10022[v10079].rhoEnergy_old+((double(0.5)*v10080)*v10022[v10079].rhoEnergy_t))
      else
        if (v10025==int32(3)) then
          v10022[v10079].rho_new += (((double(1)/double(3))*v10080)*v10022[v10079].rho_t)
          v10022[v10079].rho = (v10022[v10079].rho_old+((double(1)*v10080)*v10022[v10079].rho_t))
          var v10085 = vs_mul_double_3(v10022[v10079].rhoVelocity_t, ((double(1)/double(3))*v10080))
          var v10086 = v10022[v10079].rhoVelocity_new
          v10086[0] += v10085[0]
          v10086[1] += v10085[1]
          v10086[2] += v10085[2]
          v10022[v10079].rhoVelocity_new = v10086
          v10022[v10079].rhoVelocity = vv_add_double_3(v10022[v10079].rhoVelocity_old, vs_mul_double_3(v10022[v10079].rhoVelocity_t, (double(1)*v10080)))
          v10022[v10079].rhoEnergy_new += (((double(1)/double(3))*v10080)*v10022[v10079].rhoEnergy_t)
          v10022[v10079].rhoEnergy = (v10022[v10079].rhoEnergy_old+((double(1)*v10080)*v10022[v10079].rhoEnergy_t))
        else
          v10022[v10079].rho = (v10022[v10079].rho_new+(((double(1)/double(6))*v10080)*v10022[v10079].rho_t))
          v10022[v10079].rhoVelocity = vv_add_double_3(v10022[v10079].rhoVelocity_new, vs_mul_double_3(v10022[v10079].rhoVelocity_t, ((double(1)/double(6))*v10080)))
          v10022[v10079].rhoEnergy = (v10022[v10079].rhoEnergy_new+(((double(1)/double(6))*v10080)*v10022[v10079].rhoEnergy_t))
        end
      end
    end
  end
end

__demand(__parallel)
task Particles_UpdateVars(v10143 : region(ispace(int1d), particles_columns), v10145 : double, v10146 : int32)
where
  reads(v10143.position), writes(v10143.position), reads(v10143.position_new), reads(v10143.position_new), writes(v10143.position_new), reads(v10143.position_old), reads(v10143.position_t), reads(v10143.temperature), writes(v10143.temperature), reads(v10143.temperature_new), reads(v10143.temperature_new), writes(v10143.temperature_new), reads(v10143.temperature_old), reads(v10143.temperature_t), reads(v10143.velocity), writes(v10143.velocity), reads(v10143.velocity_new), reads(v10143.velocity_new), writes(v10143.velocity_new), reads(v10143.velocity_old), reads(v10143.velocity_t), reads(v10143.__valid)
do
  for v10261 in v10143 do
    if v10143[v10261].__valid then
      var v10262 = v10145
      if (v10146==int32(1)) then
        var v10263 = vs_mul_double_3(v10143[v10261].position_t, ((double(1)/double(6))*v10262))
        var v10264 = v10143[v10261].position_new
        v10264[0] += v10263[0]
        v10264[1] += v10263[1]
        v10264[2] += v10263[2]
        v10143[v10261].position_new = v10264
        v10143[v10261].position = vv_add_double_3(v10143[v10261].position_old, vs_mul_double_3(v10143[v10261].position_t, (double(0.5)*v10262)))
        var v10265 = vs_mul_double_3(v10143[v10261].velocity_t, ((double(1)/double(6))*v10262))
        var v10266 = v10143[v10261].velocity_new
        v10266[0] += v10265[0]
        v10266[1] += v10265[1]
        v10266[2] += v10265[2]
        v10143[v10261].velocity_new = v10266
        v10143[v10261].velocity = vv_add_double_3(v10143[v10261].velocity_old, vs_mul_double_3(v10143[v10261].velocity_t, (double(0.5)*v10262)))
        v10143[v10261].temperature_new += (((double(1)/double(6))*v10262)*v10143[v10261].temperature_t)
        v10143[v10261].temperature = (v10143[v10261].temperature_old+((double(0.5)*v10262)*v10143[v10261].temperature_t))
      else
        if (v10146==int32(2)) then
          var v10267 = vs_mul_double_3(v10143[v10261].position_t, ((double(1)/double(3))*v10262))
          var v10268 = v10143[v10261].position_new
          v10268[0] += v10267[0]
          v10268[1] += v10267[1]
          v10268[2] += v10267[2]
          v10143[v10261].position_new = v10268
          v10143[v10261].position = vv_add_double_3(v10143[v10261].position_old, vs_mul_double_3(v10143[v10261].position_t, (double(0.5)*v10262)))
          var v10269 = vs_mul_double_3(v10143[v10261].velocity_t, ((double(1)/double(3))*v10262))
          var v10270 = v10143[v10261].velocity_new
          v10270[0] += v10269[0]
          v10270[1] += v10269[1]
          v10270[2] += v10269[2]
          v10143[v10261].velocity_new = v10270
          v10143[v10261].velocity = vv_add_double_3(v10143[v10261].velocity_old, vs_mul_double_3(v10143[v10261].velocity_t, (double(0.5)*v10262)))
          v10143[v10261].temperature_new += (((double(1)/double(3))*v10262)*v10143[v10261].temperature_t)
          v10143[v10261].temperature = (v10143[v10261].temperature_old+((double(0.5)*v10262)*v10143[v10261].temperature_t))
        else
          if (v10146==int32(3)) then
            var v10271 = vs_mul_double_3(v10143[v10261].position_t, ((double(1)/double(3))*v10262))
            var v10272 = v10143[v10261].position_new
            v10272[0] += v10271[0]
            v10272[1] += v10271[1]
            v10272[2] += v10271[2]
            v10143[v10261].position_new = v10272
            v10143[v10261].position = vv_add_double_3(v10143[v10261].position_old, vs_mul_double_3(v10143[v10261].position_t, (double(1)*v10262)))
            var v10273 = vs_mul_double_3(v10143[v10261].velocity_t, ((double(1)/double(3))*v10262))
            var v10274 = v10143[v10261].velocity_new
            v10274[0] += v10273[0]
            v10274[1] += v10273[1]
            v10274[2] += v10273[2]
            v10143[v10261].velocity_new = v10274
            v10143[v10261].velocity = vv_add_double_3(v10143[v10261].velocity_old, vs_mul_double_3(v10143[v10261].velocity_t, (double(1)*v10262)))
            v10143[v10261].temperature_new += (((double(1)/double(3))*v10262)*v10143[v10261].temperature_t)
            v10143[v10261].temperature = (v10143[v10261].temperature_old+((double(1)*v10262)*v10143[v10261].temperature_t))
          else
            v10143[v10261].position = vv_add_double_3(v10143[v10261].position_new, vs_mul_double_3(v10143[v10261].position_t, ((double(1)/double(6))*v10262)))
            v10143[v10261].velocity = vv_add_double_3(v10143[v10261].velocity_new, vs_mul_double_3(v10143[v10261].velocity_t, ((double(1)/double(6))*v10262)))
            v10143[v10261].temperature = (v10143[v10261].temperature_new+(((double(1)/double(6))*v10262)*v10143[v10261].temperature_t))
          end
        end
      end
    else
    end
  end
end

__demand(__parallel)
task Particles_UpdateAuxiliaryStep1(v10390 : region(ispace(int1d), particles_columns), v10392 : int32, v10393 : int32, v10394 : int32, v10395 : int32, v10396 : int32, v10397 : int32, v10398 : double, v10399 : double, v10400 : double, v10401 : double, v10402 : double, v10403 : double, v10404 : double)
where
  reads(v10390.position), reads(v10390.position_ghost), writes(v10390.position_ghost), reads(v10390.velocity), reads(v10390.velocity_ghost), writes(v10390.velocity_ghost), reads(v10390.velocity_ghost), writes(v10390.velocity_ghost), reads(v10390.velocity_t), reads(v10390.velocity_t_ghost), writes(v10390.velocity_t_ghost), reads(v10390.velocity_t_ghost), writes(v10390.velocity_t_ghost), reads(v10390.__valid)
do
  for v10526 in v10390 do
    if v10390[v10526].__valid then
      v10390[v10526].position_ghost[int32(0)] = v10390[v10526].position[int32(0)]
      v10390[v10526].position_ghost[int32(1)] = v10390[v10526].position[int32(1)]
      v10390[v10526].position_ghost[int32(2)] = v10390[v10526].position[int32(2)]
      v10390[v10526].velocity_ghost[int32(0)] = v10390[v10526].velocity[int32(0)]
      v10390[v10526].velocity_ghost[int32(1)] = v10390[v10526].velocity[int32(1)]
      v10390[v10526].velocity_ghost[int32(2)] = v10390[v10526].velocity[int32(2)]
      v10390[v10526].velocity_t_ghost[int32(0)] = v10390[v10526].velocity_t[int32(0)]
      v10390[v10526].velocity_t_ghost[int32(1)] = v10390[v10526].velocity_t[int32(1)]
      v10390[v10526].velocity_t_ghost[int32(2)] = v10390[v10526].velocity_t[int32(2)]
      if (v10390[v10526].position[int32(0)]<v10398) then
        if (v10392==int32(0)) then
          v10390[v10526].position_ghost[int32(0)] = (v10390[v10526].position[int32(0)]+v10399)
        else
          if (v10392==int32(1)) then
            v10390[v10526].position_ghost[int32(0)] = v10398
            var v10527 = ((-(double(1)+v10404))*v10390[v10526].velocity[int32(0)])
            if (v10527<=double(int32(0))) then
              v10390[v10526].velocity_ghost[int32(0)] += v10527
            else
            end
            var v10528 = (double(-1)*v10390[v10526].velocity_t[int32(0)])
            if (v10528>double(int32(0))) then
              v10390[v10526].velocity_t_ghost[int32(0)] += v10528
            else
            end
          else
            regentlib.assert(false, "(Liszt assertion)")
          end
        end
      else
      end
      if (v10390[v10526].position[int32(0)]>(v10398+v10399)) then
        if (v10393==int32(0)) then
          v10390[v10526].position_ghost[int32(0)] = (v10390[v10526].position[int32(0)]-v10399)
        else
          if (v10393==int32(1)) then
            v10390[v10526].position_ghost[int32(0)] = (v10398+v10399)
            var v10529 = ((-(double(1)+v10404))*v10390[v10526].velocity[int32(0)])
            if (v10529>=double(int32(0))) then
              v10390[v10526].velocity_ghost[int32(0)] += v10529
            else
            end
            var v10530 = (double(-1)*v10390[v10526].velocity_t[int32(0)])
            if (v10530<double(int32(0))) then
              v10390[v10526].velocity_t_ghost[int32(0)] += v10530
            else
            end
          else
            regentlib.assert(false, "(Liszt assertion)")
          end
        end
      else
      end
      if (v10390[v10526].position[int32(1)]<v10400) then
        if (v10394==int32(0)) then
          v10390[v10526].position_ghost[int32(1)] = (v10390[v10526].position[int32(1)]+v10401)
        else
          if (v10394==int32(1)) then
            v10390[v10526].position_ghost[int32(1)] = v10400
            var v10531 = ((-(double(1)+v10404))*v10390[v10526].velocity[int32(1)])
            if (v10531<=double(int32(0))) then
              v10390[v10526].velocity_ghost[int32(1)] += v10531
            else
            end
            var v10532 = (double(-1)*v10390[v10526].velocity_t[int32(1)])
            if (v10532>double(int32(0))) then
              v10390[v10526].velocity_t_ghost[int32(1)] += v10532
            else
            end
          else
            regentlib.assert(false, "(Liszt assertion)")
          end
        end
      else
      end
      if (v10390[v10526].position[int32(1)]>(v10400+v10401)) then
        if (v10395==int32(0)) then
          v10390[v10526].position_ghost[int32(1)] = (v10390[v10526].position[int32(1)]-v10401)
        else
          if (v10395==int32(1)) then
            v10390[v10526].position_ghost[int32(1)] = (v10400+v10401)
            var v10533 = ((-(double(1)+v10404))*v10390[v10526].velocity[int32(1)])
            if (v10533>=double(int32(0))) then
              v10390[v10526].velocity_ghost[int32(1)] += v10533
            else
            end
            var v10534 = (double(-1)*v10390[v10526].velocity_t[int32(1)])
            if (v10534<double(int32(0))) then
              v10390[v10526].velocity_t_ghost[int32(1)] += v10534
            else
            end
          else
            regentlib.assert(false, "(Liszt assertion)")
          end
        end
      else
      end
      if (v10390[v10526].position[int32(2)]<v10402) then
        if (v10396==int32(0)) then
          v10390[v10526].position_ghost[int32(2)] = (v10390[v10526].position[int32(2)]+v10403)
        else
          if (v10396==int32(1)) then
            v10390[v10526].position_ghost[int32(2)] = v10402
            var v10535 = ((-(double(1)+v10404))*v10390[v10526].velocity[int32(2)])
            if (v10535<=double(int32(0))) then
              v10390[v10526].velocity_ghost[int32(2)] += v10535
            else
            end
            var v10536 = (double(-1)*v10390[v10526].velocity_t[int32(2)])
            if (v10536>double(int32(0))) then
              v10390[v10526].velocity_t_ghost[int32(2)] += v10536
            else
            end
          else
            regentlib.assert(false, "(Liszt assertion)")
          end
        end
      else
      end
      if (v10390[v10526].position[int32(2)]>(v10402+v10403)) then
        if (v10397==int32(0)) then
          v10390[v10526].position_ghost[int32(2)] = (v10390[v10526].position[int32(2)]-v10403)
        else
          if (v10397==int32(1)) then
            v10390[v10526].position_ghost[int32(2)] = (v10402+v10403)
            var v10537 = ((-(double(1)+v10404))*v10390[v10526].velocity[int32(2)])
            if (v10537>=double(int32(0))) then
              v10390[v10526].velocity_ghost[int32(2)] += v10537
            else
            end
            var v10538 = (double(-1)*v10390[v10526].velocity_t[int32(2)])
            if (v10538<double(int32(0))) then
              v10390[v10526].velocity_t_ghost[int32(2)] += v10538
            else
            end
          else
            regentlib.assert(false, "(Liszt assertion)")
          end
        end
      else
      end
    else
    end
  end
end

__demand(__parallel)
task Particles_UpdateAuxiliaryStep2(v10580 : region(ispace(int1d), particles_columns))
where
  reads(v10580.position), writes(v10580.position), reads(v10580.position_ghost), reads(v10580.velocity), writes(v10580.velocity), reads(v10580.velocity_ghost), reads(v10580.velocity_t), writes(v10580.velocity_t), reads(v10580.velocity_t_ghost), reads(v10580.__valid)
do
  for v10583 in v10580 do
    if v10580[v10583].__valid then
      v10580[v10583].position = v10580[v10583].position_ghost
      v10580[v10583].velocity = v10580[v10583].velocity_ghost
      v10580[v10583].velocity_t = v10580[v10583].velocity_t_ghost
    else
    end
  end
end

task Particles_DeleteEscapingParticles(v10589 : region(ispace(int1d), particles_columns), v10591 : double, v10592 : double, v10593 : double, v10594 : double, v10595 : double, v10596 : double) : int64
where
  reads(v10589.position), reads(v10589.__valid), writes(v10589.__valid), reads(v10589.__valid)
do
  var v10635 = int64(0)
  for v10636 in v10589 do
    if v10589[v10636].__valid then
      var v10637 = v10591
      var v10638 = (v10591+v10592)
      var v10639 = v10593
      var v10640 = (v10593+v10594)
      var v10641 = v10595
      var v10642 = (v10595+v10596)
      var v10643 = v10589[v10636].position
      if ((((((v10643[int32(0)]>v10638) or (v10643[int32(0)]<v10637)) or (v10643[int32(1)]>v10640)) or (v10643[int32(1)]<v10639)) or (v10643[int32(2)]>v10642)) or (v10643[int32(2)]<v10641)) then
        v10589[v10636].__valid = false
        v10635 += (-int64(int32(1)))
      else
      end
    else
    end
  end
  return v10635
end

-- task print______(v10675 : double)
--   C.printf("\n Current time step: %2.6e s.\n", v10675)
-- end
--
-- task print_______(v10678 : double, v10679 : double)
--   C.printf(" Min Flow Temp: %11.6f K. Max Flow Temp: %11.6f K.\n", v10678, v10679)
-- end
--
-- task print________(v10682 : int64)
--   C.printf(" Current number of particles: %d.\n", v10682)
-- end
--
-- task print_________()
--   C.printf("\n")
-- end
--
-- task print__________()
--   C.printf("    Iter     Time(s)   Avg Press    Avg Temp      Avg KE  Particle T\n")
-- end
--
-- task print___________(v10690 : int32, v10691 : double, v10692 : double, v10693 : double, v10694 : double, v10695 : double)
--   C.printf("%8d %11.6f %11.6f %11.6f %11.6f %11.6f\n", v10690, v10691, v10692, v10693, v10694, v10695)
-- end
--
-- task print____________(v10728 : double)
--   C.printf("\n Current time step: %2.6e s.\n", v10728)
-- end
--
-- task print_____________(v10731 : double, v10732 : double)
--   C.printf(" Min Flow Temp: %11.6f K. Max Flow Temp: %11.6f K.\n", v10731, v10732)
-- end
--
-- task print______________(v10735 : int64)
--   C.printf(" Current number of particles: %d.\n", v10735)
-- end
--
-- task print_______________()
--   C.printf("\n")
-- end
--
-- task print________________()
--   C.printf("    Iter     Time(s)   Avg Press    Avg Temp      Avg KE  Particle T\n")
-- end
--
-- task print_________________(v10743 : int32, v10744 : double, v10745 : double, v10746 : double, v10747 : double, v10748 : double)
--   C.printf("%8d %11.6f %11.6f %11.6f %11.6f %11.6f\n", v10743, v10744, v10745, v10746, v10747, v10748)
-- end

task work(v46 : Config)
-- leaf (false), inner (false), idempotent (false)
  var v12725 = v46.grid.xTiles
  var v12726 = v46.grid.yTiles
  var v12727 = v46.grid.zTiles
  var v12728 = double(int32(0))
  var v12729 = double(int32(1))
  var v12730 = v46.grid.xNum
  var v12731 = v46.grid.yNum
  var v12732 = v46.grid.zNum
  var v12733 = v46.grid.origin[0]
  var v12734 = v46.grid.origin[1]
  var v12735 = v46.grid.origin[2]
  var v12736 = v46.grid.xWidth
  var v12737 = v46.grid.yWidth
  var v12738 = v46.grid.zWidth
  var v12739 = (v12736/v12730)
  var v12740 = (v12737/v12731)
  var v12741 = (v12738/v12732)
  var v12742 = ((v12739*v12740)*v12741)
  var v12743 = (((((int32(1)/v12739)*int32(1))/v12739)+(((int32(1)/v12740)*int32(1))/v12740))+(((int32(1)/v12741)*int32(1))/v12741))
  var v12744 = (v46.bc.xBCLeft==int32(0))
  var v12745 = array(double(0.1), double(0.1), double(0.1))
  var v12746 = array(double(0.1), double(0.1), double(0.1))
  var v12747 = array(double(0.1), double(0.1), double(0.1))
  var v12748 = double(int32(0))
  var v12749 = double(int32(0))
  var v12750 = (v46.bc.yBCLeft==int32(0))
  var v12751 = array(double(0.1), double(0.1), double(0.1))
  var v12752 = array(double(0.1), double(0.1), double(0.1))
  var v12753 = array(double(0.1), double(0.1), double(0.1))
  var v12754 = double(int32(0))
  var v12755 = double(int32(0))
  var v12756 = (v46.bc.zBCLeft==int32(0))
  var v12757 = array(double(0.1), double(0.1), double(0.1))
  var v12758 = array(double(0.1), double(0.1), double(0.1))
  var v12759 = array(double(0.1), double(0.1), double(0.1))
  var v12760 = double(int32(0))
  var v12761 = double(int32(0))
  var v12762 = int32(-1)
  var v12763 = int32(-1)
  var v12764 = int32(-1)
  var v12765 = int32(-1)
  var v12766 = int32(-1)
  var v12767 = int32(-1)
  var v12768 = min(int32(1), ((v46.bc.xBCLeft-int32(0))*(v46.bc.xBCLeft-int32(0))))
  var v12769 = min(int32(1), ((v46.bc.yBCLeft-int32(0))*(v46.bc.yBCLeft-int32(0))))
  var v12770 = min(int32(1), ((v46.bc.zBCLeft-int32(0))*(v46.bc.zBCLeft-int32(0))))
  var v12771 = (v12733-(v12739*v12768))
  var v12772 = (v12734-(v12740*v12769))
  var v12773 = (v12735-(v12741*v12770))
  var v12774 = (v12736+(int32(2)*(v12739*v12768)))
  var v12775 = (v12737+(int32(2)*(v12740*v12769)))
  var v12776 = (v12738+(int32(2)*(v12741*v12770)))
  var v12777 = double(int32(0))
  var v12778 = double(int32(0))
  var v12779 = int32(0)
  var v12780 = double(0.0001)
  var v12781 = int32(0)
  var v12782 = double(int32(0))
  var v12783 = double(int32(0))
  var v12784 = double(int32(0))
  var v12785 = v46.flow.gasConstant
  var v12786 = v46.flow.gamma
  var v12787 = v46.flow.prandtl
  var v12788 = v46.flow.viscosityModel
  var v12789 = v46.flow.constantVisc
  var v12790 = v46.flow.powerlawViscRef
  var v12791 = v46.flow.powerlawTempRef
  var v12792 = v46.flow.sutherlandViscRef
  var v12793 = v46.flow.sutherlandTempRef
  var v12794 = v46.flow.sutherlandSRef
  var v12795 = v46.flow.initParams
  var v12796 = v46.flow.bodyForce
  var v12797 = double(int32(0))
  var v12798 = double(int32(0))
  var v12799 = double(int32(0))
  var v12800 = double(int32(0))
  var v12801 = double(int32(0))
  var v12802 = double(int32(0))
  var v12803 = double(int32(0))
  var v12804 = double(int32(0))
  var v12805 = double(int32(0))
  var v12806 = v46.particles.maxNum
  var v12807 = v46.particles.restitutionCoeff
  var v12808 = v46.particles.convectiveCoeff
  var v12809 = v46.particles.heatCapacity
  var v12810 = v46.particles.density
  var v12811 = v46.particles.bodyForce
  var v12812 = v46.particles.maxSkew
  var v12813 = v46.particles.maxXferNum
  var v12814 = double(int32(0))
  var v12815 = int64(int32(0))
  var v12816 = v46.radiation.qa
  var v12817 = v46.radiation.qs
  var v12818 = v46.radiation.xNum
  var v12819 = v46.radiation.yNum
  var v12820 = v46.radiation.zNum
  var v12821 = int32(0)
  var v12822 = int32(0)
  var v12823 = int32(0)
  var v12824 = false
  var v12825 = false
  var v12826 = false
  var v12827 = (v12736/v12818)
  var v12828 = (v12737/v12819)
  var v12829 = (v12738/v12820)
  var v12830 = ((v12827*v12828)*v12829)
  var v12831 = ispace(int3d, int3d({x = (v12730+(2*v12768)), y = (v12731+(2*v12769)), z = (v12732+(2*v12770))}))
  var v12832 = region(v12831, Fluid_columns)
  var v12833 = region(v12831, Fluid_columns)
  var v12834 = ispace(int1d, int1d((C.ceil(((v12806/((v12725*v12726)*v12727))*v12812))*((v12725*v12726)*v12727))))
  var v12835 = region(v12834, particles_columns)
  var v12836 = region(v12834, particles_columns)
  var v12837 = ispace(int3d, int3d({x = (v12818+(2*v12821)), y = (v12819+(2*v12822)), z = (v12820+(2*v12823))}))
  var v12838 = region(v12837, Radiation_columns)
  var v12839 = region(v12837, Radiation_columns)
  var v12840 = ispace(int3d, int3d({v12725, v12726, v12727}))
  regentlib.assert(((v12730%v12725)==0), "Uneven partitioning")
  regentlib.assert(((v12731%v12726)==0), "Uneven partitioning")
  regentlib.assert(((v12732%v12727)==0), "Uneven partitioning")
  var v12841 = regentlib.c.legion_domain_point_coloring_create()
  for v12842 in v12840 do
    var v12843 = rect3d({lo = int3d({x = (v12768+((v12730/v12725)*v12842.x)), y = (v12769+((v12731/v12726)*v12842.y)), z = (v12770+((v12732/v12727)*v12842.z))}), hi = int3d({x = ((v12768+((v12730/v12725)*(v12842.x+1)))-1), y = ((v12769+((v12731/v12726)*(v12842.y+1)))-1), z = ((v12770+((v12732/v12727)*(v12842.z+1)))-1)})})
    if (v12842.x==0) then
      v12843.lo.x -= v12768
    else
    end
    if (v12842.x==(v12725-1)) then
      v12843.hi.x += v12768
    else
    end
    if (v12842.y==0) then
      v12843.lo.y -= v12769
    else
    end
    if (v12842.y==(v12726-1)) then
      v12843.hi.y += v12769
    else
    end
    if (v12842.z==0) then
      v12843.lo.z -= v12770
    else
    end
    if (v12842.z==(v12727-1)) then
      v12843.hi.z += v12770
    else
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v12841, regentlib.c.legion_domain_point_t(v12842), regentlib.c.legion_domain_t(v12843))
  end
  var v12844 = partition(disjoint, v12832, v12841, v12840)
  var v12845 = partition(disjoint, v12833, v12841, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v12841)
  regentlib.assert(((v12806%((v12725*v12726)*v12727))==0), "Uneven partitioning")
  var v12846 = regentlib.c.legion_domain_point_coloring_create()
  for v12847 : int32 = 0, v12727 do
    for v12848 : int32 = 0, v12726 do
      for v12849 : int32 = 0, v12725 do
        var v12850 : int64
        for v12851 in v12835 do
          v12850 = int64((v12851+(((((v12847*v12725)*v12726)+(v12848*v12725))+v12849)*C.ceil(((v12806/((v12725*v12726)*v12727))*v12812)))))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v12846, regentlib.c.legion_domain_point_t(int3d({v12849, v12848, v12847})), regentlib.c.legion_domain_t(rect1d({v12850, ((v12850+C.ceil(((v12806/((v12725*v12726)*v12727))*v12812)))-1)})))
      end
    end
  end
  var v12852 = partition(disjoint, v12835, v12846, v12840)
  var v12853 = partition(disjoint, v12836, v12846, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v12846)
  var v12854 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v12855 = regentlib.c.legion_domain_point_coloring_create()
  for v12856 : int32 = 0, v12727 do
    for v12857 : int32 = 0, v12726 do
      for v12858 : int32 = 0, v12725 do
        var v12859 : int64
        for v12860 in v12854 do
          v12859 = int64((v12860+(((((v12856*v12725)*v12726)+(v12857*v12725))+v12858)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v12855, regentlib.c.legion_domain_point_t(int3d({v12858, v12857, v12856})), regentlib.c.legion_domain_t(rect1d({v12859, ((v12859+v12813)-1)})))
      end
    end
  end
  var v12861 = partition(disjoint, v12854, v12855, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v12855)
  var v12862 = regentlib.c.legion_domain_point_coloring_create()
  var v12863 = int3d({0, 0, 1})
  for v12864 in v12840 do
    var v12865 : int64
    for v12866 in v12861[(((v12864-v12863)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v12865 = int64(int1d(v12866))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v12862, regentlib.c.legion_domain_point_t(v12864), regentlib.c.legion_domain_t(rect1d({v12865, ((v12865+v12813)-1)})))
  end
  var v12867 = partition(aliased, v12854, v12862, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v12862)
  var v12868 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v12869 = regentlib.c.legion_domain_point_coloring_create()
  for v12870 : int32 = 0, v12727 do
    for v12871 : int32 = 0, v12726 do
      for v12872 : int32 = 0, v12725 do
        var v12873 : int64
        for v12874 in v12868 do
          v12873 = int64((v12874+(((((v12870*v12725)*v12726)+(v12871*v12725))+v12872)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v12869, regentlib.c.legion_domain_point_t(int3d({v12872, v12871, v12870})), regentlib.c.legion_domain_t(rect1d({v12873, ((v12873+v12813)-1)})))
      end
    end
  end
  var v12875 = partition(disjoint, v12868, v12869, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v12869)
  var v12876 = regentlib.c.legion_domain_point_coloring_create()
  var v12877 = int3d({0, 0, -1})
  for v12878 in v12840 do
    var v12879 : int64
    for v12880 in v12875[(((v12878-v12877)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v12879 = int64(int1d(v12880))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v12876, regentlib.c.legion_domain_point_t(v12878), regentlib.c.legion_domain_t(rect1d({v12879, ((v12879+v12813)-1)})))
  end
  var v12881 = partition(aliased, v12868, v12876, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v12876)
  var v12882 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v12883 = regentlib.c.legion_domain_point_coloring_create()
  for v12884 : int32 = 0, v12727 do
    for v12885 : int32 = 0, v12726 do
      for v12886 : int32 = 0, v12725 do
        var v12887 : int64
        for v12888 in v12882 do
          v12887 = int64((v12888+(((((v12884*v12725)*v12726)+(v12885*v12725))+v12886)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v12883, regentlib.c.legion_domain_point_t(int3d({v12886, v12885, v12884})), regentlib.c.legion_domain_t(rect1d({v12887, ((v12887+v12813)-1)})))
      end
    end
  end
  var v12889 = partition(disjoint, v12882, v12883, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v12883)
  var v12890 = regentlib.c.legion_domain_point_coloring_create()
  var v12891 = int3d({0, 1, 0})
  for v12892 in v12840 do
    var v12893 : int64
    for v12894 in v12889[(((v12892-v12891)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v12893 = int64(int1d(v12894))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v12890, regentlib.c.legion_domain_point_t(v12892), regentlib.c.legion_domain_t(rect1d({v12893, ((v12893+v12813)-1)})))
  end
  var v12895 = partition(aliased, v12882, v12890, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v12890)
  var v12896 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v12897 = regentlib.c.legion_domain_point_coloring_create()
  for v12898 : int32 = 0, v12727 do
    for v12899 : int32 = 0, v12726 do
      for v12900 : int32 = 0, v12725 do
        var v12901 : int64
        for v12902 in v12896 do
          v12901 = int64((v12902+(((((v12898*v12725)*v12726)+(v12899*v12725))+v12900)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v12897, regentlib.c.legion_domain_point_t(int3d({v12900, v12899, v12898})), regentlib.c.legion_domain_t(rect1d({v12901, ((v12901+v12813)-1)})))
      end
    end
  end
  var v12903 = partition(disjoint, v12896, v12897, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v12897)
  var v12904 = regentlib.c.legion_domain_point_coloring_create()
  var v12905 = int3d({0, 1, 1})
  for v12906 in v12840 do
    var v12907 : int64
    for v12908 in v12903[(((v12906-v12905)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v12907 = int64(int1d(v12908))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v12904, regentlib.c.legion_domain_point_t(v12906), regentlib.c.legion_domain_t(rect1d({v12907, ((v12907+v12813)-1)})))
  end
  var v12909 = partition(aliased, v12896, v12904, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v12904)
  var v12910 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v12911 = regentlib.c.legion_domain_point_coloring_create()
  for v12912 : int32 = 0, v12727 do
    for v12913 : int32 = 0, v12726 do
      for v12914 : int32 = 0, v12725 do
        var v12915 : int64
        for v12916 in v12910 do
          v12915 = int64((v12916+(((((v12912*v12725)*v12726)+(v12913*v12725))+v12914)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v12911, regentlib.c.legion_domain_point_t(int3d({v12914, v12913, v12912})), regentlib.c.legion_domain_t(rect1d({v12915, ((v12915+v12813)-1)})))
      end
    end
  end
  var v12917 = partition(disjoint, v12910, v12911, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v12911)
  var v12918 = regentlib.c.legion_domain_point_coloring_create()
  var v12919 = int3d({0, 1, -1})
  for v12920 in v12840 do
    var v12921 : int64
    for v12922 in v12917[(((v12920-v12919)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v12921 = int64(int1d(v12922))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v12918, regentlib.c.legion_domain_point_t(v12920), regentlib.c.legion_domain_t(rect1d({v12921, ((v12921+v12813)-1)})))
  end
  var v12923 = partition(aliased, v12910, v12918, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v12918)
  var v12924 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v12925 = regentlib.c.legion_domain_point_coloring_create()
  for v12926 : int32 = 0, v12727 do
    for v12927 : int32 = 0, v12726 do
      for v12928 : int32 = 0, v12725 do
        var v12929 : int64
        for v12930 in v12924 do
          v12929 = int64((v12930+(((((v12926*v12725)*v12726)+(v12927*v12725))+v12928)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v12925, regentlib.c.legion_domain_point_t(int3d({v12928, v12927, v12926})), regentlib.c.legion_domain_t(rect1d({v12929, ((v12929+v12813)-1)})))
      end
    end
  end
  var v12931 = partition(disjoint, v12924, v12925, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v12925)
  var v12932 = regentlib.c.legion_domain_point_coloring_create()
  var v12933 = int3d({0, -1, 0})
  for v12934 in v12840 do
    var v12935 : int64
    for v12936 in v12931[(((v12934-v12933)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v12935 = int64(int1d(v12936))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v12932, regentlib.c.legion_domain_point_t(v12934), regentlib.c.legion_domain_t(rect1d({v12935, ((v12935+v12813)-1)})))
  end
  var v12937 = partition(aliased, v12924, v12932, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v12932)
  var v12938 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v12939 = regentlib.c.legion_domain_point_coloring_create()
  for v12940 : int32 = 0, v12727 do
    for v12941 : int32 = 0, v12726 do
      for v12942 : int32 = 0, v12725 do
        var v12943 : int64
        for v12944 in v12938 do
          v12943 = int64((v12944+(((((v12940*v12725)*v12726)+(v12941*v12725))+v12942)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v12939, regentlib.c.legion_domain_point_t(int3d({v12942, v12941, v12940})), regentlib.c.legion_domain_t(rect1d({v12943, ((v12943+v12813)-1)})))
      end
    end
  end
  var v12945 = partition(disjoint, v12938, v12939, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v12939)
  var v12946 = regentlib.c.legion_domain_point_coloring_create()
  var v12947 = int3d({0, -1, 1})
  for v12948 in v12840 do
    var v12949 : int64
    for v12950 in v12945[(((v12948-v12947)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v12949 = int64(int1d(v12950))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v12946, regentlib.c.legion_domain_point_t(v12948), regentlib.c.legion_domain_t(rect1d({v12949, ((v12949+v12813)-1)})))
  end
  var v12951 = partition(aliased, v12938, v12946, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v12946)
  var v12952 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v12953 = regentlib.c.legion_domain_point_coloring_create()
  for v12954 : int32 = 0, v12727 do
    for v12955 : int32 = 0, v12726 do
      for v12956 : int32 = 0, v12725 do
        var v12957 : int64
        for v12958 in v12952 do
          v12957 = int64((v12958+(((((v12954*v12725)*v12726)+(v12955*v12725))+v12956)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v12953, regentlib.c.legion_domain_point_t(int3d({v12956, v12955, v12954})), regentlib.c.legion_domain_t(rect1d({v12957, ((v12957+v12813)-1)})))
      end
    end
  end
  var v12959 = partition(disjoint, v12952, v12953, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v12953)
  var v12960 = regentlib.c.legion_domain_point_coloring_create()
  var v12961 = int3d({0, -1, -1})
  for v12962 in v12840 do
    var v12963 : int64
    for v12964 in v12959[(((v12962-v12961)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v12963 = int64(int1d(v12964))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v12960, regentlib.c.legion_domain_point_t(v12962), regentlib.c.legion_domain_t(rect1d({v12963, ((v12963+v12813)-1)})))
  end
  var v12965 = partition(aliased, v12952, v12960, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v12960)
  var v12966 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v12967 = regentlib.c.legion_domain_point_coloring_create()
  for v12968 : int32 = 0, v12727 do
    for v12969 : int32 = 0, v12726 do
      for v12970 : int32 = 0, v12725 do
        var v12971 : int64
        for v12972 in v12966 do
          v12971 = int64((v12972+(((((v12968*v12725)*v12726)+(v12969*v12725))+v12970)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v12967, regentlib.c.legion_domain_point_t(int3d({v12970, v12969, v12968})), regentlib.c.legion_domain_t(rect1d({v12971, ((v12971+v12813)-1)})))
      end
    end
  end
  var v12973 = partition(disjoint, v12966, v12967, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v12967)
  var v12974 = regentlib.c.legion_domain_point_coloring_create()
  var v12975 = int3d({1, 0, 0})
  for v12976 in v12840 do
    var v12977 : int64
    for v12978 in v12973[(((v12976-v12975)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v12977 = int64(int1d(v12978))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v12974, regentlib.c.legion_domain_point_t(v12976), regentlib.c.legion_domain_t(rect1d({v12977, ((v12977+v12813)-1)})))
  end
  var v12979 = partition(aliased, v12966, v12974, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v12974)
  var v12980 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v12981 = regentlib.c.legion_domain_point_coloring_create()
  for v12982 : int32 = 0, v12727 do
    for v12983 : int32 = 0, v12726 do
      for v12984 : int32 = 0, v12725 do
        var v12985 : int64
        for v12986 in v12980 do
          v12985 = int64((v12986+(((((v12982*v12725)*v12726)+(v12983*v12725))+v12984)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v12981, regentlib.c.legion_domain_point_t(int3d({v12984, v12983, v12982})), regentlib.c.legion_domain_t(rect1d({v12985, ((v12985+v12813)-1)})))
      end
    end
  end
  var v12987 = partition(disjoint, v12980, v12981, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v12981)
  var v12988 = regentlib.c.legion_domain_point_coloring_create()
  var v12989 = int3d({1, 0, 1})
  for v12990 in v12840 do
    var v12991 : int64
    for v12992 in v12987[(((v12990-v12989)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v12991 = int64(int1d(v12992))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v12988, regentlib.c.legion_domain_point_t(v12990), regentlib.c.legion_domain_t(rect1d({v12991, ((v12991+v12813)-1)})))
  end
  var v12993 = partition(aliased, v12980, v12988, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v12988)
  var v12994 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v12995 = regentlib.c.legion_domain_point_coloring_create()
  for v12996 : int32 = 0, v12727 do
    for v12997 : int32 = 0, v12726 do
      for v12998 : int32 = 0, v12725 do
        var v12999 : int64
        for v13000 in v12994 do
          v12999 = int64((v13000+(((((v12996*v12725)*v12726)+(v12997*v12725))+v12998)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v12995, regentlib.c.legion_domain_point_t(int3d({v12998, v12997, v12996})), regentlib.c.legion_domain_t(rect1d({v12999, ((v12999+v12813)-1)})))
      end
    end
  end
  var v13001 = partition(disjoint, v12994, v12995, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v12995)
  var v13002 = regentlib.c.legion_domain_point_coloring_create()
  var v13003 = int3d({1, 0, -1})
  for v13004 in v12840 do
    var v13005 : int64
    for v13006 in v13001[(((v13004-v13003)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v13005 = int64(int1d(v13006))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v13002, regentlib.c.legion_domain_point_t(v13004), regentlib.c.legion_domain_t(rect1d({v13005, ((v13005+v12813)-1)})))
  end
  var v13007 = partition(aliased, v12994, v13002, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13002)
  var v13008 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v13009 = regentlib.c.legion_domain_point_coloring_create()
  for v13010 : int32 = 0, v12727 do
    for v13011 : int32 = 0, v12726 do
      for v13012 : int32 = 0, v12725 do
        var v13013 : int64
        for v13014 in v13008 do
          v13013 = int64((v13014+(((((v13010*v12725)*v12726)+(v13011*v12725))+v13012)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v13009, regentlib.c.legion_domain_point_t(int3d({v13012, v13011, v13010})), regentlib.c.legion_domain_t(rect1d({v13013, ((v13013+v12813)-1)})))
      end
    end
  end
  var v13015 = partition(disjoint, v13008, v13009, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13009)
  var v13016 = regentlib.c.legion_domain_point_coloring_create()
  var v13017 = int3d({1, 1, 0})
  for v13018 in v12840 do
    var v13019 : int64
    for v13020 in v13015[(((v13018-v13017)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v13019 = int64(int1d(v13020))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v13016, regentlib.c.legion_domain_point_t(v13018), regentlib.c.legion_domain_t(rect1d({v13019, ((v13019+v12813)-1)})))
  end
  var v13021 = partition(aliased, v13008, v13016, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13016)
  var v13022 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v13023 = regentlib.c.legion_domain_point_coloring_create()
  for v13024 : int32 = 0, v12727 do
    for v13025 : int32 = 0, v12726 do
      for v13026 : int32 = 0, v12725 do
        var v13027 : int64
        for v13028 in v13022 do
          v13027 = int64((v13028+(((((v13024*v12725)*v12726)+(v13025*v12725))+v13026)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v13023, regentlib.c.legion_domain_point_t(int3d({v13026, v13025, v13024})), regentlib.c.legion_domain_t(rect1d({v13027, ((v13027+v12813)-1)})))
      end
    end
  end
  var v13029 = partition(disjoint, v13022, v13023, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13023)
  var v13030 = regentlib.c.legion_domain_point_coloring_create()
  var v13031 = int3d({1, 1, 1})
  for v13032 in v12840 do
    var v13033 : int64
    for v13034 in v13029[(((v13032-v13031)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v13033 = int64(int1d(v13034))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v13030, regentlib.c.legion_domain_point_t(v13032), regentlib.c.legion_domain_t(rect1d({v13033, ((v13033+v12813)-1)})))
  end
  var v13035 = partition(aliased, v13022, v13030, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13030)
  var v13036 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v13037 = regentlib.c.legion_domain_point_coloring_create()
  for v13038 : int32 = 0, v12727 do
    for v13039 : int32 = 0, v12726 do
      for v13040 : int32 = 0, v12725 do
        var v13041 : int64
        for v13042 in v13036 do
          v13041 = int64((v13042+(((((v13038*v12725)*v12726)+(v13039*v12725))+v13040)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v13037, regentlib.c.legion_domain_point_t(int3d({v13040, v13039, v13038})), regentlib.c.legion_domain_t(rect1d({v13041, ((v13041+v12813)-1)})))
      end
    end
  end
  var v13043 = partition(disjoint, v13036, v13037, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13037)
  var v13044 = regentlib.c.legion_domain_point_coloring_create()
  var v13045 = int3d({1, 1, -1})
  for v13046 in v12840 do
    var v13047 : int64
    for v13048 in v13043[(((v13046-v13045)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v13047 = int64(int1d(v13048))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v13044, regentlib.c.legion_domain_point_t(v13046), regentlib.c.legion_domain_t(rect1d({v13047, ((v13047+v12813)-1)})))
  end
  var v13049 = partition(aliased, v13036, v13044, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13044)
  var v13050 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v13051 = regentlib.c.legion_domain_point_coloring_create()
  for v13052 : int32 = 0, v12727 do
    for v13053 : int32 = 0, v12726 do
      for v13054 : int32 = 0, v12725 do
        var v13055 : int64
        for v13056 in v13050 do
          v13055 = int64((v13056+(((((v13052*v12725)*v12726)+(v13053*v12725))+v13054)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v13051, regentlib.c.legion_domain_point_t(int3d({v13054, v13053, v13052})), regentlib.c.legion_domain_t(rect1d({v13055, ((v13055+v12813)-1)})))
      end
    end
  end
  var v13057 = partition(disjoint, v13050, v13051, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13051)
  var v13058 = regentlib.c.legion_domain_point_coloring_create()
  var v13059 = int3d({1, -1, 0})
  for v13060 in v12840 do
    var v13061 : int64
    for v13062 in v13057[(((v13060-v13059)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v13061 = int64(int1d(v13062))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v13058, regentlib.c.legion_domain_point_t(v13060), regentlib.c.legion_domain_t(rect1d({v13061, ((v13061+v12813)-1)})))
  end
  var v13063 = partition(aliased, v13050, v13058, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13058)
  var v13064 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v13065 = regentlib.c.legion_domain_point_coloring_create()
  for v13066 : int32 = 0, v12727 do
    for v13067 : int32 = 0, v12726 do
      for v13068 : int32 = 0, v12725 do
        var v13069 : int64
        for v13070 in v13064 do
          v13069 = int64((v13070+(((((v13066*v12725)*v12726)+(v13067*v12725))+v13068)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v13065, regentlib.c.legion_domain_point_t(int3d({v13068, v13067, v13066})), regentlib.c.legion_domain_t(rect1d({v13069, ((v13069+v12813)-1)})))
      end
    end
  end
  var v13071 = partition(disjoint, v13064, v13065, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13065)
  var v13072 = regentlib.c.legion_domain_point_coloring_create()
  var v13073 = int3d({1, -1, 1})
  for v13074 in v12840 do
    var v13075 : int64
    for v13076 in v13071[(((v13074-v13073)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v13075 = int64(int1d(v13076))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v13072, regentlib.c.legion_domain_point_t(v13074), regentlib.c.legion_domain_t(rect1d({v13075, ((v13075+v12813)-1)})))
  end
  var v13077 = partition(aliased, v13064, v13072, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13072)
  var v13078 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v13079 = regentlib.c.legion_domain_point_coloring_create()
  for v13080 : int32 = 0, v12727 do
    for v13081 : int32 = 0, v12726 do
      for v13082 : int32 = 0, v12725 do
        var v13083 : int64
        for v13084 in v13078 do
          v13083 = int64((v13084+(((((v13080*v12725)*v12726)+(v13081*v12725))+v13082)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v13079, regentlib.c.legion_domain_point_t(int3d({v13082, v13081, v13080})), regentlib.c.legion_domain_t(rect1d({v13083, ((v13083+v12813)-1)})))
      end
    end
  end
  var v13085 = partition(disjoint, v13078, v13079, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13079)
  var v13086 = regentlib.c.legion_domain_point_coloring_create()
  var v13087 = int3d({1, -1, -1})
  for v13088 in v12840 do
    var v13089 : int64
    for v13090 in v13085[(((v13088-v13087)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v13089 = int64(int1d(v13090))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v13086, regentlib.c.legion_domain_point_t(v13088), regentlib.c.legion_domain_t(rect1d({v13089, ((v13089+v12813)-1)})))
  end
  var v13091 = partition(aliased, v13078, v13086, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13086)
  var v13092 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v13093 = regentlib.c.legion_domain_point_coloring_create()
  for v13094 : int32 = 0, v12727 do
    for v13095 : int32 = 0, v12726 do
      for v13096 : int32 = 0, v12725 do
        var v13097 : int64
        for v13098 in v13092 do
          v13097 = int64((v13098+(((((v13094*v12725)*v12726)+(v13095*v12725))+v13096)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v13093, regentlib.c.legion_domain_point_t(int3d({v13096, v13095, v13094})), regentlib.c.legion_domain_t(rect1d({v13097, ((v13097+v12813)-1)})))
      end
    end
  end
  var v13099 = partition(disjoint, v13092, v13093, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13093)
  var v13100 = regentlib.c.legion_domain_point_coloring_create()
  var v13101 = int3d({-1, 0, 0})
  for v13102 in v12840 do
    var v13103 : int64
    for v13104 in v13099[(((v13102-v13101)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v13103 = int64(int1d(v13104))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v13100, regentlib.c.legion_domain_point_t(v13102), regentlib.c.legion_domain_t(rect1d({v13103, ((v13103+v12813)-1)})))
  end
  var v13105 = partition(aliased, v13092, v13100, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13100)
  var v13106 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v13107 = regentlib.c.legion_domain_point_coloring_create()
  for v13108 : int32 = 0, v12727 do
    for v13109 : int32 = 0, v12726 do
      for v13110 : int32 = 0, v12725 do
        var v13111 : int64
        for v13112 in v13106 do
          v13111 = int64((v13112+(((((v13108*v12725)*v12726)+(v13109*v12725))+v13110)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v13107, regentlib.c.legion_domain_point_t(int3d({v13110, v13109, v13108})), regentlib.c.legion_domain_t(rect1d({v13111, ((v13111+v12813)-1)})))
      end
    end
  end
  var v13113 = partition(disjoint, v13106, v13107, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13107)
  var v13114 = regentlib.c.legion_domain_point_coloring_create()
  var v13115 = int3d({-1, 0, 1})
  for v13116 in v12840 do
    var v13117 : int64
    for v13118 in v13113[(((v13116-v13115)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v13117 = int64(int1d(v13118))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v13114, regentlib.c.legion_domain_point_t(v13116), regentlib.c.legion_domain_t(rect1d({v13117, ((v13117+v12813)-1)})))
  end
  var v13119 = partition(aliased, v13106, v13114, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13114)
  var v13120 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v13121 = regentlib.c.legion_domain_point_coloring_create()
  for v13122 : int32 = 0, v12727 do
    for v13123 : int32 = 0, v12726 do
      for v13124 : int32 = 0, v12725 do
        var v13125 : int64
        for v13126 in v13120 do
          v13125 = int64((v13126+(((((v13122*v12725)*v12726)+(v13123*v12725))+v13124)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v13121, regentlib.c.legion_domain_point_t(int3d({v13124, v13123, v13122})), regentlib.c.legion_domain_t(rect1d({v13125, ((v13125+v12813)-1)})))
      end
    end
  end
  var v13127 = partition(disjoint, v13120, v13121, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13121)
  var v13128 = regentlib.c.legion_domain_point_coloring_create()
  var v13129 = int3d({-1, 0, -1})
  for v13130 in v12840 do
    var v13131 : int64
    for v13132 in v13127[(((v13130-v13129)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v13131 = int64(int1d(v13132))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v13128, regentlib.c.legion_domain_point_t(v13130), regentlib.c.legion_domain_t(rect1d({v13131, ((v13131+v12813)-1)})))
  end
  var v13133 = partition(aliased, v13120, v13128, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13128)
  var v13134 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v13135 = regentlib.c.legion_domain_point_coloring_create()
  for v13136 : int32 = 0, v12727 do
    for v13137 : int32 = 0, v12726 do
      for v13138 : int32 = 0, v12725 do
        var v13139 : int64
        for v13140 in v13134 do
          v13139 = int64((v13140+(((((v13136*v12725)*v12726)+(v13137*v12725))+v13138)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v13135, regentlib.c.legion_domain_point_t(int3d({v13138, v13137, v13136})), regentlib.c.legion_domain_t(rect1d({v13139, ((v13139+v12813)-1)})))
      end
    end
  end
  var v13141 = partition(disjoint, v13134, v13135, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13135)
  var v13142 = regentlib.c.legion_domain_point_coloring_create()
  var v13143 = int3d({-1, 1, 0})
  for v13144 in v12840 do
    var v13145 : int64
    for v13146 in v13141[(((v13144-v13143)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v13145 = int64(int1d(v13146))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v13142, regentlib.c.legion_domain_point_t(v13144), regentlib.c.legion_domain_t(rect1d({v13145, ((v13145+v12813)-1)})))
  end
  var v13147 = partition(aliased, v13134, v13142, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13142)
  var v13148 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v13149 = regentlib.c.legion_domain_point_coloring_create()
  for v13150 : int32 = 0, v12727 do
    for v13151 : int32 = 0, v12726 do
      for v13152 : int32 = 0, v12725 do
        var v13153 : int64
        for v13154 in v13148 do
          v13153 = int64((v13154+(((((v13150*v12725)*v12726)+(v13151*v12725))+v13152)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v13149, regentlib.c.legion_domain_point_t(int3d({v13152, v13151, v13150})), regentlib.c.legion_domain_t(rect1d({v13153, ((v13153+v12813)-1)})))
      end
    end
  end
  var v13155 = partition(disjoint, v13148, v13149, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13149)
  var v13156 = regentlib.c.legion_domain_point_coloring_create()
  var v13157 = int3d({-1, 1, 1})
  for v13158 in v12840 do
    var v13159 : int64
    for v13160 in v13155[(((v13158-v13157)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v13159 = int64(int1d(v13160))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v13156, regentlib.c.legion_domain_point_t(v13158), regentlib.c.legion_domain_t(rect1d({v13159, ((v13159+v12813)-1)})))
  end
  var v13161 = partition(aliased, v13148, v13156, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13156)
  var v13162 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v13163 = regentlib.c.legion_domain_point_coloring_create()
  for v13164 : int32 = 0, v12727 do
    for v13165 : int32 = 0, v12726 do
      for v13166 : int32 = 0, v12725 do
        var v13167 : int64
        for v13168 in v13162 do
          v13167 = int64((v13168+(((((v13164*v12725)*v12726)+(v13165*v12725))+v13166)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v13163, regentlib.c.legion_domain_point_t(int3d({v13166, v13165, v13164})), regentlib.c.legion_domain_t(rect1d({v13167, ((v13167+v12813)-1)})))
      end
    end
  end
  var v13169 = partition(disjoint, v13162, v13163, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13163)
  var v13170 = regentlib.c.legion_domain_point_coloring_create()
  var v13171 = int3d({-1, 1, -1})
  for v13172 in v12840 do
    var v13173 : int64
    for v13174 in v13169[(((v13172-v13171)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v13173 = int64(int1d(v13174))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v13170, regentlib.c.legion_domain_point_t(v13172), regentlib.c.legion_domain_t(rect1d({v13173, ((v13173+v12813)-1)})))
  end
  var v13175 = partition(aliased, v13162, v13170, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13170)
  var v13176 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v13177 = regentlib.c.legion_domain_point_coloring_create()
  for v13178 : int32 = 0, v12727 do
    for v13179 : int32 = 0, v12726 do
      for v13180 : int32 = 0, v12725 do
        var v13181 : int64
        for v13182 in v13176 do
          v13181 = int64((v13182+(((((v13178*v12725)*v12726)+(v13179*v12725))+v13180)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v13177, regentlib.c.legion_domain_point_t(int3d({v13180, v13179, v13178})), regentlib.c.legion_domain_t(rect1d({v13181, ((v13181+v12813)-1)})))
      end
    end
  end
  var v13183 = partition(disjoint, v13176, v13177, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13177)
  var v13184 = regentlib.c.legion_domain_point_coloring_create()
  var v13185 = int3d({-1, -1, 0})
  for v13186 in v12840 do
    var v13187 : int64
    for v13188 in v13183[(((v13186-v13185)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v13187 = int64(int1d(v13188))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v13184, regentlib.c.legion_domain_point_t(v13186), regentlib.c.legion_domain_t(rect1d({v13187, ((v13187+v12813)-1)})))
  end
  var v13189 = partition(aliased, v13176, v13184, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13184)
  var v13190 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v13191 = regentlib.c.legion_domain_point_coloring_create()
  for v13192 : int32 = 0, v12727 do
    for v13193 : int32 = 0, v12726 do
      for v13194 : int32 = 0, v12725 do
        var v13195 : int64
        for v13196 in v13190 do
          v13195 = int64((v13196+(((((v13192*v12725)*v12726)+(v13193*v12725))+v13194)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v13191, regentlib.c.legion_domain_point_t(int3d({v13194, v13193, v13192})), regentlib.c.legion_domain_t(rect1d({v13195, ((v13195+v12813)-1)})))
      end
    end
  end
  var v13197 = partition(disjoint, v13190, v13191, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13191)
  var v13198 = regentlib.c.legion_domain_point_coloring_create()
  var v13199 = int3d({-1, -1, 1})
  for v13200 in v12840 do
    var v13201 : int64
    for v13202 in v13197[(((v13200-v13199)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v13201 = int64(int1d(v13202))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v13198, regentlib.c.legion_domain_point_t(v13200), regentlib.c.legion_domain_t(rect1d({v13201, ((v13201+v12813)-1)})))
  end
  var v13203 = partition(aliased, v13190, v13198, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13198)
  var v13204 = region(ispace(int1d, int1d((v12813*((v12725*v12726)*v12727)))), int8[376])
  var v13205 = regentlib.c.legion_domain_point_coloring_create()
  for v13206 : int32 = 0, v12727 do
    for v13207 : int32 = 0, v12726 do
      for v13208 : int32 = 0, v12725 do
        var v13209 : int64
        for v13210 in v13204 do
          v13209 = int64((v13210+(((((v13206*v12725)*v12726)+(v13207*v12725))+v13208)*v12813)))
          break
        end
        regentlib.c.legion_domain_point_coloring_color_domain(v13205, regentlib.c.legion_domain_point_t(int3d({v13208, v13207, v13206})), regentlib.c.legion_domain_t(rect1d({v13209, ((v13209+v12813)-1)})))
      end
    end
  end
  var v13211 = partition(disjoint, v13204, v13205, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13205)
  var v13212 = regentlib.c.legion_domain_point_coloring_create()
  var v13213 = int3d({-1, -1, -1})
  for v13214 in v12840 do
    var v13215 : int64
    for v13216 in v13211[(((v13214-v13213)+{v12725, v12726, v12727})%{v12725, v12726, v12727})] do
      v13215 = int64(int1d(v13216))
      break
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v13212, regentlib.c.legion_domain_point_t(v13214), regentlib.c.legion_domain_t(rect1d({v13215, ((v13215+v12813)-1)})))
  end
  var v13217 = partition(aliased, v13204, v13212, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13212)
  regentlib.assert(((v12818%v12725)==0), "Uneven partitioning")
  regentlib.assert(((v12819%v12726)==0), "Uneven partitioning")
  regentlib.assert(((v12820%v12727)==0), "Uneven partitioning")
  var v13218 = regentlib.c.legion_domain_point_coloring_create()
  for v13219 in v12840 do
    var v13220 = rect3d({lo = int3d({x = (v12821+((v12818/v12725)*v13219.x)), y = (v12822+((v12819/v12726)*v13219.y)), z = (v12823+((v12820/v12727)*v13219.z))}), hi = int3d({x = ((v12821+((v12818/v12725)*(v13219.x+1)))-1), y = ((v12822+((v12819/v12726)*(v13219.y+1)))-1), z = ((v12823+((v12820/v12727)*(v13219.z+1)))-1)})})
    if (v13219.x==0) then
      v13220.lo.x -= v12821
    else
    end
    if (v13219.x==(v12725-1)) then
      v13220.hi.x += v12821
    else
    end
    if (v13219.y==0) then
      v13220.lo.y -= v12822
    else
    end
    if (v13219.y==(v12726-1)) then
      v13220.hi.y += v12822
    else
    end
    if (v13219.z==0) then
      v13220.lo.z -= v12823
    else
    end
    if (v13219.z==(v12727-1)) then
      v13220.hi.z += v12823
    else
    end
    regentlib.c.legion_domain_point_coloring_color_domain(v13218, regentlib.c.legion_domain_point_t(v13219), regentlib.c.legion_domain_t(v13220))
  end
  var v13221 = partition(disjoint, v12838, v13218, v12840)
  var v13222 = partition(disjoint, v12839, v13218, v12840)
  regentlib.c.legion_domain_point_coloring_destroy(v13218)
  __parallelize_with v12844, v12852, v13221, v12840, (image(v12832, v12852, v12835.cell)<=v12844) do
    particles_initValidField(v12835)
    var v13223 = int32(((not ((v12730%v12818)==int32(0))) or ((not ((v12731%v12819)==int32(0))) or (not ((v12732%v12820)==int32(0))))))
    while (v13223>0) do
      regentlib.assert(false, "Inexact coarsening factor")
      v13223 -= 1
    end
    SetCoarseningField(v12832, v12768, v12730, v12769, v12731, v12770, v12732, v12821, v12818, v12822, v12819, v12823, v12820)
    var v13224 = int32(((v46.bc.xBCLeft==int32(0)) and (v46.bc.xBCRight==int32(0))))
    var v13225 = (1-v13224)
    while (v13224>0) do
      v12745 = array(v12729, v12729, v12729)
      v12746 = array(v12728, v12728, v12728)
      v12747 = array(v12728, v12728, v12728)
      v12748 = (-v12729)
      v12749 = (-v12729)
      v12762 = int32(0)
      v12763 = int32(0)
      v13224 -= 1
    end
    while (v13225>0) do
      var v13226 = int32(((v46.bc.xBCLeft==int32(1)) and (v46.bc.xBCRight==int32(1))))
      var v13227 = (1-v13226)
      while (v13226>0) do
        v12745 = array((-v12729), v12729, v12729)
        v12746 = array(v12728, v12728, v12728)
        v12747 = array(v12728, v12728, v12728)
        v12748 = (-v12729)
        v12749 = (-v12729)
        v12762 = int32(1)
        v12763 = int32(1)
        v13226 -= 1
      end
      while (v13227>0) do
        var v13228 = int32(((v46.bc.xBCLeft==int32(2)) and (v46.bc.xBCRight==int32(2))))
        var v13229 = (1-v13228)
        while (v13228>0) do
          v12745 = array((-v12729), (-v12729), (-v12729))
          v12746 = array((int32(2)*v46.bc.xBCRightVel[0]), (int32(2)*v46.bc.xBCRightVel[1]), (int32(2)*v46.bc.xBCRightVel[2]))
          v12747 = array((int32(2)*v46.bc.xBCLeftVel[0]), (int32(2)*v46.bc.xBCLeftVel[1]), (int32(2)*v46.bc.xBCLeftVel[2]))
          v12748 = (-v12729)
          v12749 = (-v12729)
          v12762 = int32(1)
          v12763 = int32(1)
          v13228 -= 1
        end
        while (v13229>0) do
          var v13230 = int32(((v46.bc.xBCLeft==int32(3)) and (v46.bc.xBCRight==int32(3))))
          var v13231 = (1-v13230)
          while (v13230>0) do
            v12745 = array((-v12729), (-v12729), (-v12729))
            v12746 = array((int32(2)*v46.bc.xBCRightVel[0]), (int32(2)*v46.bc.xBCRightVel[1]), (int32(2)*v46.bc.xBCRightVel[2]))
            v12747 = array((int32(2)*v46.bc.xBCLeftVel[0]), (int32(2)*v46.bc.xBCLeftVel[1]), (int32(2)*v46.bc.xBCLeftVel[2]))
            v12748 = v46.bc.xBCRightTemp
            v12749 = v46.bc.xBCLeftTemp
            v12762 = int32(1)
            v12763 = int32(1)
            v13230 -= 1
          end
          while (v13231>0) do
            regentlib.assert(false, "Boundary conditions in x not implemented")
            v13231 -= 1
          end
          v13229 -= 1
        end
        v13227 -= 1
      end
      v13225 -= 1
    end
    var v13232 = int32(((v46.bc.yBCLeft==int32(0)) and (v46.bc.yBCRight==int32(0))))
    var v13233 = (1-v13232)
    while (v13232>0) do
      v12751 = array(v12729, v12729, v12729)
      v12752 = array(v12728, v12728, v12728)
      v12753 = array(v12728, v12728, v12728)
      v12754 = (-v12729)
      v12755 = (-v12729)
      v12764 = int32(0)
      v12765 = int32(0)
      v13232 -= 1
    end
    while (v13233>0) do
      var v13234 = int32(((v46.bc.yBCLeft==int32(1)) and (v46.bc.yBCRight==int32(1))))
      var v13235 = (1-v13234)
      while (v13234>0) do
        v12751 = array(v12729, (-v12729), v12729)
        v12752 = array(v12728, v12728, v12728)
        v12753 = array(v12728, v12728, v12728)
        v12754 = (-v12729)
        v12755 = (-v12729)
        v12764 = int32(1)
        v12765 = int32(1)
        v13234 -= 1
      end
      while (v13235>0) do
        var v13236 = int32(((v46.bc.yBCLeft==int32(2)) and (v46.bc.yBCRight==int32(2))))
        var v13237 = (1-v13236)
        while (v13236>0) do
          v12751 = array((-v12729), (-v12729), (-v12729))
          v12752 = array((int32(2)*v46.bc.yBCRightVel[0]), (int32(2)*v46.bc.yBCRightVel[1]), (int32(2)*v46.bc.yBCRightVel[2]))
          v12753 = array((int32(2)*v46.bc.yBCLeftVel[0]), (int32(2)*v46.bc.yBCLeftVel[1]), (int32(2)*v46.bc.yBCLeftVel[2]))
          v12754 = (-v12729)
          v12755 = (-v12729)
          v12764 = int32(1)
          v12765 = int32(1)
          v13236 -= 1
        end
        while (v13237>0) do
          var v13238 = int32(((v46.bc.yBCLeft==int32(3)) and (v46.bc.yBCRight==int32(3))))
          var v13239 = (1-v13238)
          while (v13238>0) do
            v12751 = array((-v12729), (-v12729), (-v12729))
            v12752 = array((int32(2)*v46.bc.yBCRightVel[0]), (int32(2)*v46.bc.yBCRightVel[1]), (int32(2)*v46.bc.yBCRightVel[2]))
            v12753 = array((int32(2)*v46.bc.yBCLeftVel[0]), (int32(2)*v46.bc.yBCLeftVel[1]), (int32(2)*v46.bc.yBCLeftVel[2]))
            v12754 = v46.bc.yBCRightTemp
            v12755 = v46.bc.yBCLeftTemp
            v12764 = int32(1)
            v12765 = int32(1)
            v13238 -= 1
          end
          while (v13239>0) do
            regentlib.assert(false, "Boundary conditions in y not implemented")
            v13239 -= 1
          end
          v13237 -= 1
        end
        v13235 -= 1
      end
      v13233 -= 1
    end
    var v13240 = int32(((v46.bc.zBCLeft==int32(0)) and (v46.bc.zBCRight==int32(0))))
    var v13241 = (1-v13240)
    while (v13240>0) do
      v12757 = array(v12729, v12729, v12729)
      v12758 = array(v12728, v12728, v12728)
      v12759 = array(v12728, v12728, v12728)
      v12760 = (-v12729)
      v12761 = (-v12729)
      v12766 = int32(0)
      v12767 = int32(0)
      v13240 -= 1
    end
    while (v13241>0) do
      var v13242 = int32(((v46.bc.zBCLeft==int32(1)) and (v46.bc.zBCRight==int32(1))))
      var v13243 = (1-v13242)
      while (v13242>0) do
        v12757 = array(v12729, v12729, (-v12729))
        v12758 = array(v12728, v12728, v12728)
        v12759 = array(v12728, v12728, v12728)
        v12760 = (-v12729)
        v12761 = (-v12729)
        v12766 = int32(1)
        v12767 = int32(1)
        v13242 -= 1
      end
      while (v13243>0) do
        var v13244 = int32(((v46.bc.zBCLeft==int32(2)) and (v46.bc.zBCRight==int32(2))))
        var v13245 = (1-v13244)
        while (v13244>0) do
          v12757 = array((-v12729), (-v12729), (-v12729))
          v12758 = array((int32(2)*v46.bc.zBCRightVel[0]), (int32(2)*v46.bc.zBCRightVel[1]), (int32(2)*v46.bc.zBCRightVel[2]))
          v12759 = array((int32(2)*v46.bc.zBCLeftVel[0]), (int32(2)*v46.bc.zBCLeftVel[1]), (int32(2)*v46.bc.zBCLeftVel[2]))
          v12760 = (-v12729)
          v12761 = (-v12729)
          v12766 = int32(1)
          v12767 = int32(1)
          v13244 -= 1
        end
        while (v13245>0) do
          var v13246 = int32(((v46.bc.zBCLeft==int32(3)) and (v46.bc.zBCRight==int32(3))))
          var v13247 = (1-v13246)
          while (v13246>0) do
            v12757 = array((-v12729), (-v12729), (-v12729))
            v12758 = array((int32(2)*v46.bc.zBCRightVel[0]), (int32(2)*v46.bc.zBCRightVel[1]), (int32(2)*v46.bc.zBCRightVel[2]))
            v12759 = array((int32(2)*v46.bc.zBCLeftVel[0]), (int32(2)*v46.bc.zBCLeftVel[1]), (int32(2)*v46.bc.zBCLeftVel[2]))
            v12760 = v46.bc.zBCRightTemp
            v12761 = v46.bc.zBCLeftTemp
            v12766 = int32(1)
            v12767 = int32(1)
            v13246 -= 1
          end
          while (v13247>0) do
            regentlib.assert(false, "Boundary conditions in z not implemented")
            v13247 -= 1
          end
          v13245 -= 1
        end
        v13243 -= 1
      end
      v13241 -= 1
    end
    var v13248 = int32((not (((v46.bc.xBCLeft==int32(0)) and (v46.bc.xBCRight==int32(0))) or ((not (v46.bc.xBCLeft==int32(0))) and (not (v46.bc.xBCRight==int32(0)))))))
    while (v13248>0) do
      regentlib.assert(false, "Boundary conditions in x should match for periodicity")
      v13248 -= 1
    end
    var v13249 = int32((not (((v46.bc.yBCLeft==int32(0)) and (v46.bc.yBCRight==int32(0))) or ((not (v46.bc.yBCLeft==int32(0))) and (not (v46.bc.yBCRight==int32(0)))))))
    while (v13249>0) do
      regentlib.assert(false, "Boundary conditions in y should match for periodicity")
      v13249 -= 1
    end
    var v13250 = int32((not (((v46.bc.zBCLeft==int32(0)) and (v46.bc.zBCRight==int32(0))) or ((not (v46.bc.zBCLeft==int32(0))) and (not (v46.bc.zBCRight==int32(0)))))))
    while (v13250>0) do
      regentlib.assert(false, "Boundary conditions in z should match for periodicity")
      v13250 -= 1
    end
    var v13251 = int32((v46.flow.initCase==int32(1)))
    while (v13251>0) do
      v12779 = v46.integrator.restartIter
      v13251 -= 1
    end
    Flow_InitializeCell(v12832)
    Flow_InitializeCenterCoordinates(v12832, v12768, v12730, v12733, v12736, v12769, v12731, v12734, v12737, v12770, v12732, v12735, v12738)
    var v13252 = int32((v46.flow.initCase==int32(0)))
    while (v13252>0) do
      Flow_InitializeUniform(v12832, v12795)
      v13252 -= 1
    end
    var v13253 = int32((v46.flow.initCase==int32(3)))
    while (v13253>0) do
      Flow_InitializeTaylorGreen2D(v12832, v12795, v12768, v12730, v12733, v12736, v12769, v12731, v12734, v12737, v12770, v12732, v12735, v12738)
      v13253 -= 1
    end
    var v13254 = int32((v46.flow.initCase==int32(4)))
    while (v13254>0) do
      Flow_InitializeTaylorGreen3D(v12832, v12795, v12768, v12730, v12733, v12736, v12769, v12731, v12734, v12737, v12770, v12732, v12735, v12738)
      v13254 -= 1
    end
    var v13255 = int32((v46.flow.initCase==int32(2)))
    while (v13255>0) do
      Flow_InitializePerturbed(v12832, v12795)
      v13255 -= 1
    end
    var v13256 = int32((v46.flow.initCase==int32(1)))
    while (v13256>0) do
      var v13257 = [&int8](C.malloc(uint64(256)))
      C.snprintf(v13257, uint64(256), "restart_fluid_%d.hdf", v46.integrator.restartIter)
      attach(hdf5, v12833.{rho, pressure, velocity}, v13257, uint32(0))
      for v13258 in v12840 do
        var v13259 = v12844[int3d(v13258)]
        var v13260 = v12845[int3d(v13258)]
        acquire(v13260.{rho, pressure, velocity})
        copy(v13260.{rho, pressure, velocity}, v13259.{rho, pressure, velocity})
        release(v13260.{rho, pressure, velocity})
      end
      detach(hdf5, v12833.{rho, pressure, velocity})
      C.free([&opaque](v13257))
      v13256 -= 1
    end
    Flow_UpdateConservedFromPrimitive(v12832, v12786, v12785, v12768, v12730, v12769, v12731, v12770, v12732)
    Flow_UpdateAuxiliaryVelocity(v12832, v12768, v12730, v12769, v12731, v12770, v12732)
    Flow_UpdateGhostConservedStep1(v12832, v12749, v12747, v12748, v12746, v12745, v12755, v12753, v12754, v12752, v12751, v12761, v12759, v12760, v12758, v12757, v12786, v12785, v12768, v12730, v12769, v12731, v12770, v12732)
    Flow_UpdateGhostConservedStep2(v12832, v12768, v12730, v12769, v12731, v12770, v12732)
    Flow_UpdateGhostVelocityStep1(v12832, v12747, v12746, v12745, v12753, v12752, v12751, v12759, v12758, v12757, v12768, v12730, v12769, v12731, v12770, v12732)
    Flow_UpdateGhostVelocityStep2(v12832, v12768, v12730, v12769, v12731, v12770, v12732)
    Flow_ComputeVelocityGradientAll(v12832, v12768, v12739, v12730, v12769, v12740, v12731, v12770, v12741, v12732)
    Flow_UpdateAuxiliaryThermodynamics(v12832, v12786, v12785, v12768, v12730, v12769, v12731, v12770, v12732)
    Flow_UpdateGhostThermodynamicsStep1(v12832, v12749, v12748, v12755, v12754, v12761, v12760, v12768, v12730, v12769, v12731, v12770, v12732)
    Flow_UpdateGhostThermodynamicsStep2(v12832, v12768, v12730, v12769, v12731, v12770, v12732)
    Flow_UpdateGhostFieldsStep1(v12832, v12749, v12747, v12748, v12746, v12745, v12755, v12753, v12754, v12752, v12751, v12761, v12759, v12760, v12758, v12757, v12786, v12785, v12768, v12730, v12769, v12731, v12770, v12732)
    Flow_UpdateGhostFieldsStep2(v12832, v12768, v12730, v12769, v12731, v12770, v12732)
    var v13261 = int32((v46.particles.initCase==int32(0)))
    while (v13261>0) do
      regentlib.assert(false, "Random particle initialization is disabled")
      v13261 -= 1
    end
    var v13267 = int32((v46.particles.initCase==int32(2)))
    while (v13267>0) do
      InitParticlesUniform(v12835, v12832, v46, v12768, v12769, v12770)
      v12815 = int64(((v46.particles.initNum/((v46.grid.xTiles*v46.grid.yTiles)*v46.grid.zTiles))*((v46.grid.xTiles*v46.grid.yTiles)*v46.grid.zTiles)))
      v13267 -= 1
    end
    v12797 = double(int32(0))
    v12798 = double(int32(0))
    v12799 = double(int32(0))
    v12800 = double(int32(math.huge))
    v12801 = double(int32(-math.huge))
    v12802 = double(int32(0))
    v12803 = double(int32(0))
    v12814 = double(int32(0))
    v12797 += CalculateAveragePressure(v12832, v12742, v12768, v12730, v12769, v12731, v12770, v12732)
    v12798 += CalculateAverageTemperature(v12832, v12742, v12768, v12730, v12769, v12731, v12770, v12732)
    v12799 += CalculateAverageKineticEnergy(v12832, v12742, v12768, v12730, v12769, v12731, v12770, v12732)
    v12800 min= CalculateMinTemperature(v12832, v12768, v12730, v12769, v12731, v12770, v12732)
    v12801 max= CalculateMaxTemperature(v12832, v12768, v12730, v12769, v12731, v12770, v12732)
    v12814 += Particles_IntegrateQuantities(v12835)
    v12797 = (v12797/(((v12730*v12731)*v12732)*v12742))
    v12798 = (v12798/(((v12730*v12731)*v12732)*v12742))
    v12799 = (v12799/(((v12730*v12731)*v12732)*v12742))
    v12814 = (v12814/v12815)
    Radiation_InitializeCell(v12838)
    var v13268 = v46.grid.xTiles
    var v13269 = v46.grid.yTiles
    var v13270 = v46.grid.zTiles
    var v13271 = ispace(int3d, int3d({x = (v12818+1), y = v12819, z = v12820}))
    var v13272 = ispace(int3d, int3d({x = v12818, y = (v12819+1), z = v12820}))
    var v13273 = ispace(int3d, int3d({x = v12818, y = v12819, z = (v12820+1)}))
    var v13274 = region(v13271, face)
    var v13275 = region(v13271, face)
    var v13276 = region(v13271, face)
    var v13277 = region(v13271, face)
    var v13278 = region(v13271, face)
    var v13279 = region(v13271, face)
    var v13280 = region(v13271, face)
    var v13281 = region(v13271, face)
    var v13282 = region(v13272, face)
    var v13283 = region(v13272, face)
    var v13284 = region(v13272, face)
    var v13285 = region(v13272, face)
    var v13286 = region(v13272, face)
    var v13287 = region(v13272, face)
    var v13288 = region(v13272, face)
    var v13289 = region(v13272, face)
    var v13290 = region(v13273, face)
    var v13291 = region(v13273, face)
    var v13292 = region(v13273, face)
    var v13293 = region(v13273, face)
    var v13294 = region(v13273, face)
    var v13295 = region(v13273, face)
    var v13296 = region(v13273, face)
    var v13297 = region(v13273, face)
    var v13298 = ispace(int1d, int1d(14))
    var v13299 = region(v13298, angle)
    var v13300 = ispace(int3d, int3d({x = (v13268+1), y = v13269, z = v13270}))
    var v13301 = ispace(int3d, int3d({x = v13268, y = (v13269+1), z = v13270}))
    var v13302 = ispace(int3d, int3d({x = v13268, y = v13269, z = (v13270+1)}))
    var v13303 = make_interior_partition_x_lo(v13274, v13300, v12818, v12819, v12820, v13268, v13269, v13270)
    var v13304 = make_interior_partition_y_lo(v13282, v13301, v12818, v12819, v12820, v13268, v13269, v13270)
    var v13305 = make_interior_partition_z_lo(v13290, v13302, v12818, v12819, v12820, v13268, v13269, v13270)
    var v13306 = make_interior_partition_x_lo(v13275, v13300, v12818, v12819, v12820, v13268, v13269, v13270)
    var v13307 = make_interior_partition_y_lo(v13283, v13301, v12818, v12819, v12820, v13268, v13269, v13270)
    var v13308 = make_interior_partition_z_hi(v13291, v13302, v12818, v12819, v12820, v13268, v13269, v13270)
    var v13309 = make_interior_partition_x_lo(v13276, v13300, v12818, v12819, v12820, v13268, v13269, v13270)
    var v13310 = make_interior_partition_y_hi(v13284, v13301, v12818, v12819, v12820, v13268, v13269, v13270)
    var v13311 = make_interior_partition_z_lo(v13292, v13302, v12818, v12819, v12820, v13268, v13269, v13270)
    var v13312 = make_interior_partition_x_lo(v13277, v13300, v12818, v12819, v12820, v13268, v13269, v13270)
    var v13313 = make_interior_partition_y_hi(v13285, v13301, v12818, v12819, v12820, v13268, v13269, v13270)
    var v13314 = make_interior_partition_z_hi(v13293, v13302, v12818, v12819, v12820, v13268, v13269, v13270)
    var v13315 = make_interior_partition_x_hi(v13278, v13300, v12818, v12819, v12820, v13268, v13269, v13270)
    var v13316 = make_interior_partition_y_lo(v13286, v13301, v12818, v12819, v12820, v13268, v13269, v13270)
    var v13317 = make_interior_partition_z_lo(v13294, v13302, v12818, v12819, v12820, v13268, v13269, v13270)
    var v13318 = make_interior_partition_x_hi(v13279, v13300, v12818, v12819, v12820, v13268, v13269, v13270)
    var v13319 = make_interior_partition_y_lo(v13287, v13301, v12818, v12819, v12820, v13268, v13269, v13270)
    var v13320 = make_interior_partition_z_hi(v13295, v13302, v12818, v12819, v12820, v13268, v13269, v13270)
    var v13321 = make_interior_partition_x_hi(v13280, v13300, v12818, v12819, v12820, v13268, v13269, v13270)
    var v13322 = make_interior_partition_y_hi(v13288, v13301, v12818, v12819, v12820, v13268, v13269, v13270)
    var v13323 = make_interior_partition_z_lo(v13296, v13302, v12818, v12819, v12820, v13268, v13269, v13270)
    var v13324 = make_interior_partition_x_hi(v13281, v13300, v12818, v12819, v12820, v13268, v13269, v13270)
    var v13325 = make_interior_partition_y_hi(v13289, v13301, v12818, v12819, v12820, v13268, v13269, v13270)
    var v13326 = make_interior_partition_z_hi(v13297, v13302, v12818, v12819, v12820, v13268, v13269, v13270)
    initialize_angles(v13299)
    var v13327 = 1
    __demand(__spmd)
    while ((v12777<v46.integrator.finalTime) and (v12779<v46.integrator.maxIter)) do
      var v13341 = int32((v46.integrator.cfl<int32(0)))
      var v13342 = (1-v13341)
      while (v13341>0) do
        v12780 = v46.integrator.fixedDeltaTime
        v13341 -= 1
      end
      while (v13342>0) do
        v12782 max= CalculateConvectiveSpectralRadius(v12832, v12786, v12785, v12743, v12739, v12740, v12741)
        v12783 max= CalculateViscousSpectralRadius(v12832, v12789, v12791, v12790, v12794, v12793, v12792, v12788, v12743)
        v12784 max= CalculateHeatConductionSpectralRadius(v12832, v12789, v12786, v12785, v12791, v12790, v12787, v12794, v12793, v12792, v12788, v12743)
        v12780 = (v46.integrator.cfl/max(v12782, max(v12783, v12784)))
        v13342 -= 1
      end
      Flow_InitializeTemporaries(v12832)
      Particles_InitializeTemporaries(v12835)
      v12778 = v12777
      v12781 = int32(1)
      while (v12781<int32(5)) do
        Flow_InitializeTimeDerivatives(v12832)
        Particles_InitializeTimeDerivatives(v12835)
        Flow_UpdateGhostVelocityGradientStep1(v12832, v12745, v12751, v12757, v12768, v12730, v12769, v12731, v12770, v12732)
        Flow_UpdateGhostVelocityGradientStep2(v12832, v12768, v12730, v12769, v12731, v12770, v12732)
        Flow_AddGetFlux(v12832, v12789, v12786, v12785, v12791, v12790, v12787, v12794, v12793, v12792, v12788, v12768, v12739, v12730, v12769, v12740, v12731, v12770, v12741, v12732)
        Flow_AddUpdateUsingFlux(v12832, v12768, v12739, v12730, v12769, v12740, v12731, v12770, v12741, v12732)
        Flow_AddBodyForces(v12832, v12796, v12768, v12730, v12769, v12731, v12770, v12732)
        for v13344 in v12840 do
          Particles_LocateInCells(v12852[int3d(v13344)], v12744, v12750, v12756, v12768, v12730, v12733, v12736, v12769, v12731, v12734, v12737, v12770, v12732, v12735, v12738)
        end
        for v13345 in v12840 do
          particles_pushAll(int3d(v13345), v12852[int3d(v13345)], v12861[int3d(v13345)], v12875[int3d(v13345)], v12889[int3d(v13345)], v12903[int3d(v13345)], v12917[int3d(v13345)], v12931[int3d(v13345)], v12945[int3d(v13345)], v12959[int3d(v13345)], v12973[int3d(v13345)], v12987[int3d(v13345)], v13001[int3d(v13345)], v13015[int3d(v13345)], v13029[int3d(v13345)], v13043[int3d(v13345)], v13057[int3d(v13345)], v13071[int3d(v13345)], v13085[int3d(v13345)], v13099[int3d(v13345)], v13113[int3d(v13345)], v13127[int3d(v13345)], v13141[int3d(v13345)], v13155[int3d(v13345)], v13169[int3d(v13345)], v13183[int3d(v13345)], v13197[int3d(v13345)], v13211[int3d(v13345)], v12730, v12731, v12732, v12768, v12769, v12770, v12725, v12726, v12727)
        end
        for v13346 in v12840 do
          particles_pullAll(int3d(v13346), v12852[int3d(v13346)], v12867[int3d(v13346)], v12881[int3d(v13346)], v12895[int3d(v13346)], v12909[int3d(v13346)], v12923[int3d(v13346)], v12937[int3d(v13346)], v12951[int3d(v13346)], v12965[int3d(v13346)], v12979[int3d(v13346)], v12993[int3d(v13346)], v13007[int3d(v13346)], v13021[int3d(v13346)], v13035[int3d(v13346)], v13049[int3d(v13346)], v13063[int3d(v13346)], v13077[int3d(v13346)], v13091[int3d(v13346)], v13105[int3d(v13346)], v13119[int3d(v13346)], v13133[int3d(v13346)], v13147[int3d(v13346)], v13161[int3d(v13346)], v13175[int3d(v13346)], v13189[int3d(v13346)], v13203[int3d(v13346)], v13217[int3d(v13346)])
        end
        Particles_AddFlowCoupling(v12835, v12832, v12789, v12791, v12790, v12794, v12793, v12792, v12788, v12739, v12771, v12740, v12772, v12741, v12773, v12808, v12809)
        Particles_AddBodyForces(v12835, v12811)
        Radiation_ClearAccumulators(v12838)
        for v13347 in v12840 do
          Radiation_AccumulateParticleValues(v12852[int3d(v13347)], v12844[int3d(v13347)], v13221[int3d(v13347)])
        end
        Radiation_UpdateFieldValues(v12838, v12830, v12816, v12817)
        --var v13348 = int64(1)
        var v13349 = (v46.radiation.qs/(v46.radiation.qa+v46.radiation.qs))
        var v13327 = 1
        while (v13327>1e-06) do
          for v13350 in v12840 do
            source_term(v13221[int3d(v13350)], v13299, v13349)
          end
          for v13351 : int32 = 0, v13269 do
            for v13352 : int32 = 0, v13270 do
              west_bound(v13303[int3d({1, v13351, v13352})], v13306[int3d({1, v13351, v13352})], v13309[int3d({1, v13351, v13352})], v13312[int3d({1, v13351, v13352})], v13315[int3d({0, v13351, v13352})], v13318[int3d({0, v13351, v13352})], v13321[int3d({0, v13351, v13352})], v13324[int3d({0, v13351, v13352})], v13299, v46.radiation.emissWest, v46.radiation.tempWest)
              east_bound(v13303[int3d({v13268, v13351, v13352})], v13306[int3d({v13268, v13351, v13352})], v13309[int3d({v13268, v13351, v13352})], v13312[int3d({v13268, v13351, v13352})], v13315[int3d({(v13268-1), v13351, v13352})], v13318[int3d({(v13268-1), v13351, v13352})], v13321[int3d({(v13268-1), v13351, v13352})], v13324[int3d({(v13268-1), v13351, v13352})], v13299, v46.radiation.emissEast, v46.radiation.tempEast)
            end
          end
          for v13353 : int32 = 0, v13268 do
            for v13354 : int32 = 0, v13270 do
              south_bound(v13304[int3d({v13353, 1, v13354})], v13307[int3d({v13353, 1, v13354})], v13310[int3d({v13353, 0, v13354})], v13313[int3d({v13353, 0, v13354})], v13316[int3d({v13353, 1, v13354})], v13319[int3d({v13353, 1, v13354})], v13322[int3d({v13353, 0, v13354})], v13325[int3d({v13353, 0, v13354})], v13299, v46.radiation.emissSouth, v46.radiation.tempSouth)
              north_bound(v13304[int3d({v13353, v13269, v13354})], v13307[int3d({v13353, v13269, v13354})], v13310[int3d({v13353, (v13269-1), v13354})], v13313[int3d({v13353, (v13269-1), v13354})], v13316[int3d({v13353, v13269, v13354})], v13319[int3d({v13353, v13269, v13354})], v13322[int3d({v13353, (v13269-1), v13354})], v13325[int3d({v13353, (v13269-1), v13354})], v13299, v46.radiation.emissNorth, v46.radiation.tempNorth)
            end
          end
          for v13355 : int32 = 0, v13268 do
            for v13356 : int32 = 0, v13269 do
              up_bound(v13305[int3d({v13355, v13356, 1})], v13308[int3d({v13355, v13356, 0})], v13311[int3d({v13355, v13356, 1})], v13314[int3d({v13355, v13356, 0})], v13317[int3d({v13355, v13356, 1})], v13320[int3d({v13355, v13356, 0})], v13323[int3d({v13355, v13356, 1})], v13326[int3d({v13355, v13356, 0})], v13299, v46.radiation.emissUp, v46.radiation.tempUp)
              down_bound(v13305[int3d({v13355, v13356, v13270})], v13308[int3d({v13355, v13356, (v13270-1)})], v13311[int3d({v13355, v13356, v13270})], v13314[int3d({v13355, v13356, (v13270-1)})], v13317[int3d({v13355, v13356, v13270})], v13320[int3d({v13355, v13356, (v13270-1)})], v13323[int3d({v13355, v13356, v13270})], v13326[int3d({v13355, v13356, (v13270-1)})], v13299, v46.radiation.emissDown, v46.radiation.tempDown)
            end
          end
          for v13357 : int32 = 0, v13268 do
            for v13358 : int32 = 0, v13269 do
              for v13359 : int32 = 0, v13270 do
                sweep_1(v13221[int3d({v13357, v13358, v13359})], v13303[int3d({(v13357+1), v13358, v13359})], v13304[int3d({v13357, (v13358+1), v13359})], v13305[int3d({v13357, v13358, (v13359+1)})], v13303[int3d({v13357, v13358, v13359})], v13304[int3d({v13357, v13358, v13359})], v13305[int3d({v13357, v13358, v13359})], v13299, int64(1), int64(1), int64(1), v12827, v12828, v12829)
              end
            end
          end
          for v13360 : int32 = 0, v13268 do
            for v13361 : int32 = 0, v13269 do
              for v13362 : int32 = (v13270-1), (-1), (-1) do
                sweep_2(v13221[int3d({v13360, v13361, v13362})], v13306[int3d({(v13360+1), v13361, v13362})], v13307[int3d({v13360, (v13361+1), v13362})], v13308[int3d({v13360, v13361, v13362})], v13306[int3d({v13360, v13361, v13362})], v13307[int3d({v13360, v13361, v13362})], v13308[int3d({v13360, v13361, (v13362+1)})], v13299, int64(1), int64(1), int64((-1)), v12827, v12828, v12829)
              end
            end
          end
          for v13363 : int32 = 0, v13268 do
            for v13364 : int32 = (v13269-1), (-1), (-1) do
              for v13365 : int32 = 0, v13270 do
                sweep_3(v13221[int3d({v13363, v13364, v13365})], v13309[int3d({(v13363+1), v13364, v13365})], v13310[int3d({v13363, v13364, v13365})], v13311[int3d({v13363, v13364, (v13365+1)})], v13309[int3d({v13363, v13364, v13365})], v13310[int3d({v13363, (v13364+1), v13365})], v13311[int3d({v13363, v13364, v13365})], v13299, int64(1), int64((-1)), int64(1), v12827, v12828, v12829)
              end
            end
          end
          for v13366 : int32 = 0, v13268 do
            for v13367 : int32 = (v13269-1), (-1), (-1) do
              for v13368 : int32 = (v13270-1), (-1), (-1) do
                sweep_4(v13221[int3d({v13366, v13367, v13368})], v13312[int3d({(v13366+1), v13367, v13368})], v13313[int3d({v13366, v13367, v13368})], v13314[int3d({v13366, v13367, v13368})], v13312[int3d({v13366, v13367, v13368})], v13313[int3d({v13366, (v13367+1), v13368})], v13314[int3d({v13366, v13367, (v13368+1)})], v13299, int64(1), int64((-1)), int64((-1)), v12827, v12828, v12829)
              end
            end
          end
          for v13369 : int32 = (v13268-1), (-1), (-1) do
            for v13370 : int32 = 0, v13269 do
              for v13371 : int32 = 0, v13270 do
                sweep_5(v13221[int3d({v13369, v13370, v13371})], v13315[int3d({v13369, v13370, v13371})], v13316[int3d({v13369, (v13370+1), v13371})], v13317[int3d({v13369, v13370, (v13371+1)})], v13315[int3d({(v13369+1), v13370, v13371})], v13316[int3d({v13369, v13370, v13371})], v13317[int3d({v13369, v13370, v13371})], v13299, int64((-1)), int64(1), int64(1), v12827, v12828, v12829)
              end
            end
          end
          for v13372 : int32 = (v13268-1), (-1), (-1) do
            for v13373 : int32 = 0, v13269 do
              for v13374 : int32 = (v13270-1), (-1), (-1) do
                sweep_6(v13221[int3d({v13372, v13373, v13374})], v13318[int3d({v13372, v13373, v13374})], v13319[int3d({v13372, (v13373+1), v13374})], v13320[int3d({v13372, v13373, v13374})], v13318[int3d({(v13372+1), v13373, v13374})], v13319[int3d({v13372, v13373, v13374})], v13320[int3d({v13372, v13373, (v13374+1)})], v13299, int64((-1)), int64(1), int64((-1)), v12827, v12828, v12829)
              end
            end
          end
          for v13375 : int32 = (v13268-1), (-1), (-1) do
            for v13376 : int32 = (v13269-1), (-1), (-1) do
              for v13377 : int32 = 0, v13270 do
                sweep_7(v13221[int3d({v13375, v13376, v13377})], v13321[int3d({v13375, v13376, v13377})], v13322[int3d({v13375, v13376, v13377})], v13323[int3d({v13375, v13376, (v13377+1)})], v13321[int3d({(v13375+1), v13376, v13377})], v13322[int3d({v13375, (v13376+1), v13377})], v13323[int3d({v13375, v13376, v13377})], v13299, int64((-1)), int64((-1)), int64(1), v12827, v12828, v12829)
              end
            end
          end
          for v13378 : int32 = (v13268-1), (-1), (-1) do
            for v13379 : int32 = (v13269-1), (-1), (-1) do
              for v13380 : int32 = (v13270-1), (-1), (-1) do
                sweep_8(v13221[int3d({v13378, v13379, v13380})], v13324[int3d({v13378, v13379, v13380})], v13325[int3d({v13378, v13379, v13380})], v13326[int3d({v13378, v13379, v13380})], v13324[int3d({(v13378+1), v13379, v13380})], v13325[int3d({v13378, (v13379+1), v13380})], v13326[int3d({v13378, v13379, (v13380+1)})], v13299, int64((-1)), int64((-1)), int64((-1)), v12827, v12828, v12829)
              end
            end
          end
          v13327 = 0
          for v13381 in v12840 do
            v13327 += residual(v13221[int3d(v13381)], v12818, v12819, v12820)
          end
          v13327 = C.sqrt(v13327)
          for v13382 in v12840 do
            update(v13221[int3d(v13382)])
          end
          --if (v13348==1) then
          --  C.printf("\n")
          --  C.printf(" Iteration     Residual         \n")
          --  C.printf(" ------------------------------ \n")
          --else
          --end
          --C.printf("   %3d    %.15e \n", v13348, v13327)
          --v13348 = (v13348+1)
        end
        for v13383 in v12840 do
          reduce_intensity(v13221[int3d(v13383)], v13299)
        end
        for v13384 in v12840 do
          Particles_AbsorbRadiation(v12852[int3d(v13384)], v12844[int3d(v13384)], v13221[int3d(v13384)], v12809, v12816)
        end
        -- Flow_AddParticlesCoupling(v12835, v12832, v12742)
        Flow_UpdateVars(v12832, v12780, v12781)
        Particles_UpdateVars(v12835, v12780, v12781)
        Flow_UpdateAuxiliaryVelocity(v12832, v12768, v12730, v12769, v12731, v12770, v12732)
        Flow_UpdateGhostConservedStep1(v12832, v12749, v12747, v12748, v12746, v12745, v12755, v12753, v12754, v12752, v12751, v12761, v12759, v12760, v12758, v12757, v12786, v12785, v12768, v12730, v12769, v12731, v12770, v12732)
        Flow_UpdateGhostConservedStep2(v12832, v12768, v12730, v12769, v12731, v12770, v12732)
        Flow_UpdateGhostVelocityStep1(v12832, v12747, v12746, v12745, v12753, v12752, v12751, v12759, v12758, v12757, v12768, v12730, v12769, v12731, v12770, v12732)
        Flow_UpdateGhostVelocityStep2(v12832, v12768, v12730, v12769, v12731, v12770, v12732)
        Flow_ComputeVelocityGradientAll(v12832, v12768, v12739, v12730, v12769, v12740, v12731, v12770, v12741, v12732)
        Flow_UpdateAuxiliaryThermodynamics(v12832, v12786, v12785, v12768, v12730, v12769, v12731, v12770, v12732)
        Flow_UpdateGhostThermodynamicsStep1(v12832, v12749, v12748, v12755, v12754, v12761, v12760, v12768, v12730, v12769, v12731, v12770, v12732)
        Flow_UpdateGhostThermodynamicsStep2(v12832, v12768, v12730, v12769, v12731, v12770, v12732)
        Particles_UpdateAuxiliaryStep1(v12835, v12762, v12763, v12764, v12765, v12766, v12767, v12733, v12736, v12734, v12737, v12735, v12738, v12807)
        Particles_UpdateAuxiliaryStep2(v12835)
        v12777 = (v12778+((double(0.5)*(int32(1)+(v12781/int32(3))))*v12780))
        v12781 = (v12781+int32(1))
        for v13385 in v12840 do
          v12815 += Particles_DeleteEscapingParticles(v12852[int3d(v13385)], v12771, v12774, v12772, v12775, v12773, v12776)
        end
      end
      v12779 = (v12779+int32(1))
    end
    v12797 = double(int32(0))
    v12798 = double(int32(0))
    v12799 = double(int32(0))
    v12800 = double(int32(math.huge))
    v12801 = double(int32(-math.huge))
    v12802 = double(int32(0))
    v12803 = double(int32(0))
    v12814 = double(int32(0))
    v12797 += CalculateAveragePressure(v12832, v12742, v12768, v12730, v12769, v12731, v12770, v12732)
    v12798 += CalculateAverageTemperature(v12832, v12742, v12768, v12730, v12769, v12731, v12770, v12732)
    v12799 += CalculateAverageKineticEnergy(v12832, v12742, v12768, v12730, v12769, v12731, v12770, v12732)
    v12800 min= CalculateMinTemperature(v12832, v12768, v12730, v12769, v12731, v12770, v12732)
    v12801 max= CalculateMaxTemperature(v12832, v12768, v12730, v12769, v12731, v12770, v12732)
    v12814 += Particles_IntegrateQuantities(v12835)
    v12797 = (v12797/(((v12730*v12731)*v12732)*v12742))
    v12798 = (v12798/(((v12730*v12731)*v12732)*v12742))
    v12799 = (v12799/(((v12730*v12731)*v12732)*v12742))
    v12814 = (v12814/v12815)
    C.printf(["Test Case: %s\nCurrent time step: %2.6e s.\n" ..
              " Min Flow Temp: %11.6f K. Max Flow Temp: %11.6f K.\n" ..
              " Current number of particles: %d.\n" ..
              "\n" ..
              "    Iter     Time(s)   Avg Press    Avg Temp      Avg KE  Particle T\n" ..
              "%8d %11.6f %11.6f %11.6f %11.6f %11.6f\n"],
             v46.filename, v12780, v12800, v12801, v12815,
             v12779, v12777, v12797, v12798, v12799, v12814)
    __delete(v12832)
    __delete(v12833)
    __delete(v12835)
    __delete(v12836)
    __delete(v12838)
    __delete(v12839)
    __delete(v12854)
    __delete(v12868)
    __delete(v12882)
    __delete(v12896)
    __delete(v12910)
    __delete(v12924)
    __delete(v12938)
    __delete(v12952)
    __delete(v12966)
    __delete(v12980)
    __delete(v12994)
    __delete(v13008)
    __delete(v13022)
    __delete(v13036)
    __delete(v13050)
    __delete(v13064)
    __delete(v13078)
    __delete(v13092)
    __delete(v13106)
    __delete(v13120)
    __delete(v13134)
    __delete(v13148)
    __delete(v13162)
    __delete(v13176)
    __delete(v13190)
    __delete(v13204)
    __delete(v13274)
    __delete(v13275)
    __delete(v13276)
    __delete(v13277)
    __delete(v13278)
    __delete(v13279)
    __delete(v13280)
    __delete(v13281)
    __delete(v13282)
    __delete(v13283)
    __delete(v13284)
    __delete(v13285)
    __delete(v13286)
    __delete(v13287)
    __delete(v13288)
    __delete(v13289)
    __delete(v13290)
    __delete(v13291)
    __delete(v13292)
    __delete(v13293)
    __delete(v13294)
    __delete(v13295)
    __delete(v13296)
    __delete(v13297)
    __delete(v13299)
  end
end

terra get_configs_from_future(configs_from_mapper : regentlib.c.legion_future_t)
  var p : &int8 = [&int8](c.legion_future_get_untyped_pointer(configs_from_mapper))
  var sz : int64 = c.legion_future_get_untyped_size(configs_from_mapper)
  var num_configs = @[&uint32](p)
  p = p + sizeof(uint32)
  regentlib.assert(sz == num_configs * [sizeof(Config)] + [sizeof(uint32)],
                   "incomplete data from the mapper")
  var configs = [&Config](C.malloc([sizeof(Config)] * num_configs))
  for i = 0, num_configs do
    C.memcpy(&configs[i], [&Config](p), [sizeof(Config)])
    p = p + [sizeof(Config)]
  end
  return { configs = configs, num_configs = num_configs }
end

task main()
  var configs_from_mapper =
    c.legion_runtime_select_tunable_value(__runtime(), __context(),
                                          SOLEIL_TYPES.TUNABLE_CONFIG, 0, 0)
  var result = get_configs_from_future(configs_from_mapper)
  var configs : &Config, num_configs : uint32 = result.configs, result.num_configs
  c.legion_future_destroy(configs_from_mapper)

  for i = 0, num_configs do
    C.printf("conf %3d: %4d   x %4d   x %4d / %3d   x %3d   x %3d (# iter: %d)\n", i,
        configs[i].grid.xNum, configs[i].grid.yNum,
        configs[i].grid.zNum, configs[i].grid.xTiles,
        configs[i].grid.yTiles, configs[i].grid.zTiles,
        configs[i].integrator.maxIter)
    work(configs[i])
  end
end

regentlib.saveobj(main, "soleil_tgv_dom", "executable",
                  cmapper.register_mappers,
                  { "-lm", "-lensemble_mapper", "-L" .. root_dir,
                    "-ljsonparser", "-L" .. json_dir .. "/lib" })
