import "regent"

-------------------------------------------------------------------------------
-- IMPORTS
-------------------------------------------------------------------------------

local C = regentlib.c
local MAPPER = terralib.includec("soleil_mapper.h")

local struct Particles_columns {
  id : int64,
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

-------------------------------------------------------------------------------
-- Visualization
-------------------------------------------------------------------------------

local root_dir = arg[0]:match(".*/") or "./"
assert(os.getenv('LG_RT_DIR'), "LG_RT_DIR should be set!")
local runtime_dir = os.getenv('LG_RT_DIR') .. "/"
local legion_dir = runtime_dir .. "legion/"
local mapper_dir = runtime_dir .. "mappers/"
local realm_dir = runtime_dir .. "realm/"

render = terralib.includec("render.h",
{"-I", root_dir,
"-I", runtime_dir,
"-I", mapper_dir,
"-I", legion_dir,
"-I", realm_dir,
})



task partitionByTile(r : region(ispace(int1d), Particles_columns),
                     cs : ispace(int3d),
                     halo : int64,
                     offset : int3d)
  var N = int64(r.bounds.hi - 2*halo + 1)
  var ntx = cs.bounds.hi.x + 1
  var nty = cs.bounds.hi.y + 1
  var ntz = cs.bounds.hi.z + 1
  regentlib.assert(int64(r.bounds.lo) == 0, "Can only partition root region")
  regentlib.assert(N % (ntx*nty*ntz) == 0, "Uneven partitioning")
  regentlib.assert(-ntx <= offset.x and offset.x <= ntx, "offset.x too large")
  regentlib.assert(-nty <= offset.y and offset.y <= nty, "offset.y too large")
  regentlib.assert(-ntz <= offset.z and offset.z <= ntz, "offset.z too large")
  var coloring = regentlib.c.legion_domain_point_coloring_create()
  for c_real in cs do
    var c = (c_real - offset + {ntx,nty,ntz}) % {ntx,nty,ntz}
    var rect = rect1d{
      lo = halo + (N/ntx/nty/ntz)*(c.x*nty*ntz+c.y*ntz+c.z),
      hi = halo + (N/ntx/nty/ntz)*(c.x*nty*ntz+c.y*ntz+c.z+1) - 1}
    if c.x == 0 and
       c.y == 0 and
       c.z == 0 then rect.lo -= halo end
    if c.x == ntx-1 and
       c.y == nty-1 and
       c.z == ntz-1 then rect.hi += halo end
    regentlib.c.legion_domain_point_coloring_color_domain(coloring, c_real, rect)
  end
  var p = partition(disjoint, r, coloring, cs)
  regentlib.c.legion_domain_point_coloring_destroy(coloring)
  return p
end

task Particles_InitializeUniform(Particles : region(ispace(int1d), Particles_columns))
where
writes(Particles.{__valid, id, cell, position, velocity, density, temperature, diameter})
do
  for p in Particles do
    Particles[p].__valid = true
    Particles[p].id = 0
    Particles[p].cell = {0,0,0}
    var zero3 : double[3]
    Particles[p].position = zero3
    Particles[p].velocity = zero3
    Particles[p].density = 0.0
    Particles[p].temperature = 0.0
    Particles[p].diameter = 0.0
  end
end
-------------------------------------------------------------------------------
-- MAIN SIMULATION
-------------------------------------------------------------------------------

local function mkInstance() local INSTANCE = {}

  -----------------------------------------------------------------------------
  -- Symbols shared between quotes
  -----------------------------------------------------------------------------

  local Particles = regentlib.newsymbol()
  local p_Particles = regentlib.newsymbol()

  -----------------------------------------------------------------------------
  -- Exported symbols
  -----------------------------------------------------------------------------

  INSTANCE.Particles = Particles
  INSTANCE.p_Particles = p_Particles

  -----------------------------------------------------------------------------
  -- Symbol declaration & initialization
  -----------------------------------------------------------------------------

  function INSTANCE.DeclSymbols(config) return rquote

    var is_Particles = ispace(int1d, 100000)
    var [Particles] = region(is_Particles, Particles_columns);

    -- Partitioning domain
    var [tiles] = ispace(int3d, {1,1,5})

    -- Particles Partitioning
    var [p_Particles] = partitionByTile(Particles, tiles, 0, int3d{0,0,0})

  end end -- DeclSymbols

  -----------------------------------------------------------------------------
  -- Region initialization
  -----------------------------------------------------------------------------


  function INSTANCE.InitRegions(config) return rquote

  for c in tiles do
    Particles_InitializeUniform(p_Particles[c])

  end end -- InitRegions
end end


__forbid(__inner)
task initializeVisualization(
        Particles : region(ispace(int1d), Particles_columns),
        p_Particles : partition(disjoint, Particles, ispace(int3d))
)
where
  reads(Particles.{id, position, temperature, density})
do
  render.cxx_initialize(__runtime(), __context(),
    __raw(Particles),
    __raw(p_Particles),
    __fields([Particles].{id, position, temperature, density}),
    4,
    1000)
end

local SIM = mkInstance()

--__forbid(__optimize) __demand(__inner, __replicable)
__forbid(__optimize) __demand(__inner)
task workSingle()

  [SIM.DeclSymbols()];
  [SIM.InitRegions()];
  __fence(__execution, __block)
  var cameraFromAtUp : double[9]
  var colorScale  : double[2]
  initializeVisualization(SIM.Particles, SIM.p_Particles);
  __fence(__execution, __block)
      render.cxx_render(__runtime(), __context(),
        cameraFromAtUp,
        colorScale)
      __fence(__execution, __block)
      render.cxx_reduce(__context(), cameraFromAtUp)
      __fence(__execution, __block)
      render.cxx_saveImage(__runtime(), __context(), ".")
end


__forbid(__optimize) __demand(__inner)
task main()
      workSingle()
end

-------------------------------------------------------------------------------
-- COMPILATION CALL
-------------------------------------------------------------------------------

regentlib.saveobj(main, "soleil.o", "object", MAPPER.register_mappers)
