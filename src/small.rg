import "regent"

-------------------------------------------------------------------------------
-- IMPORTS
-------------------------------------------------------------------------------

local C = regentlib.c
local MAPPER = terralib.includec("soleil_mapper.h")
local SCHEMA = terralib.includec("config_schema.h")
local Config = SCHEMA.Config
local MultiConfig = SCHEMA.MultiConfig


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

task Particles_InitializeUniform(
  Particles : region(ispace(int1d), Particles_columns),
  tile_id : int64)
where
reads(Particles.{__valid, id, cell, position, velocity, density, temperature, diameter}),
writes(Particles.{__valid, id, cell, position, velocity, density, temperature, diameter})
do
  var id : int64 = 0
  var shifted_tile_id : int64 = tile_id * C.pow(2.0, 32)
  for p in Particles do
    Particles[p].__valid = true
    Particles[p].id = id + shifted_tile_id
    Particles[p].cell = {0,0,0}
    var zero3 : double[3]
    Particles[p].position[0] = 0.1
    Particles[p].position[1] = 0.2
    Particles[p].position[2] = 0.3
    Particles[p].velocity = zero3
    Particles[p].density = 0.0
    Particles[p].temperature = 0.0
    Particles[p].diameter = 0.0
    id = id + 1
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
  local tiles = regentlib.newsymbol()
  local NX = regentlib.newsymbol()
  local NY = regentlib.newsymbol()
  local NZ = regentlib.newsymbol()

  -----------------------------------------------------------------------------
  -- Exported symbols
  -----------------------------------------------------------------------------

  INSTANCE.Particles = Particles
  INSTANCE.p_Particles = p_Particles

  -----------------------------------------------------------------------------
  -- Symbol declaration & initialization
  -----------------------------------------------------------------------------

  function INSTANCE.DeclSymbols(config) return rquote

    -- Partitioning domain
    var [NX] = config.Mapping.tiles[0]
    var [NY] = config.Mapping.tiles[1]
    var [NZ] = config.Mapping.tiles[2]
    var numTiles = NX * NY * NZ

    var maxParticlesPerTile = config.Particles.maxNum / config.Particles.parcelSize / numTiles
    if numTiles > 1 then
      maxParticlesPerTile =
        int64(C.ceil(maxParticlesPerTile * config.Particles.maxSkew))
    end
    var is_Particles = ispace(int1d, maxParticlesPerTile * numTiles)
    var [Particles] = region(is_Particles, Particles_columns);

    var sampleId = config.Mapping.sampleId
    var info : int = sampleId
    regentlib.c.legion_logical_region_attach_semantic_information(
      __runtime(), __raw(Particles), MAPPER.SAMPLE_ID_TAG, &info, [sizeof(int)], false)

    var [tiles] = ispace(int3d, {NX,NY,NZ})

    -- Particles Partitioning
    var [p_Particles] = partitionByTile(Particles, tiles, 0, int3d{0,0,0})

  end end -- DeclSymbols

  -----------------------------------------------------------------------------
  -- Region initialization
  -----------------------------------------------------------------------------


  function INSTANCE.InitRegions(config) return rquote

  var tile_id : int64 = 0
  for c in [tiles] do
    C.printf("initialie particles for one tile\n");C.fflush(C.stdout);
    Particles_InitializeUniform(p_Particles[c], tile_id)
    tile_id = tile_id + 1
  end

  end end -- InitRegions

return INSTANCE end -- mkInstance


__forbid(__inner)
task initializeVisualization(
        config : Config,
        Particles : region(ispace(int1d), Particles_columns),
        p_Particles : partition(disjoint, Particles, ispace(int3d))
)
where
  reads(Particles.{id, position, temperature, density, __valid})
do

  render.cxx_initialize(__runtime(), __context(),
    __raw(Particles),
    __raw(p_Particles),
    __fields([Particles].{id, position, temperature, density, __valid}),
    5,
    1000,
    config.Mapping.sampleId,
    MAPPER.SAMPLE_ID_TAG,
    config.Mapping.tiles)

  C.printf("initializeVisualization\n");C.fflush(C.stdout);
  for p in Particles do
    if p.__valid then
      C.printf("particle id %ld position %g %g %g temperature %g density %g valid %d\n", 
        p.id, p.position[0], p.position[1], p.position[2], p.temperature, p.density, p.__valid);C.fflush(C.stdout);
    end
  end
end

local SIM = mkInstance()

--__forbid(__optimize) __demand(__inner, __replicable)
__forbid(__optimize) __demand(__inner)
task workSingle(config : Config)
  [SIM.DeclSymbols(config)];
  [SIM.InitRegions(rexpr config end)];
  initializeVisualization(config, SIM.Particles, SIM.p_Particles);
  __fence(__execution, __block)
    -- Visualization
      render.cxx_render(__runtime(), __context(),
        config.Visualization.cameraFromAtUp,
        config.Visualization.colorScale)
      __fence(__execution, __block)
      render.cxx_reduce(__context(), config.Visualization.cameraFromAtUp)
      __fence(__execution, __block)
      render.cxx_saveImage(__runtime(), __context(), ".")
end


local SIM0 = mkInstance()
local SIM1 = mkInstance()

--__forbid(__optimize) __demand(__inner, __replicable)
__forbid(__optimize) __demand(__inner)
task workDual(mc : MultiConfig)

  C.printf("workDual calls DeclSymbols\n");C.fflush(C.stdout);
  [SIM0.DeclSymbols(rexpr mc.configs[0] end)];
  [SIM1.DeclSymbols(rexpr mc.configs[1] end)];

  C.printf("workDual calls InitRegions\n");C.fflush(C.stdout);
  [SIM0.InitRegions(rexpr mc.configs[0] end)];
  [SIM1.InitRegions(rexpr mc.configs[1] end)];
  C.printf("workDual back from InitRegions\n");C.fflush(C.stdout);
  __fence(__execution, __block)
  C.printf("workDual calls initializeVisualization\n");C.fflush(C.stdout);
  var cameraFromAtUp : double[9]
  var colorScale  : double[2]
  initializeVisualization(mc.configs[1], SIM1.Particles, SIM1.p_Particles);
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
  var args = regentlib.c.legion_runtime_get_input_args()
  for i = 1, args.argc do
    if C.strcmp(args.argv[i], '-m') == 0 and i < args.argc-1 then
      var mc : MultiConfig
      SCHEMA.parse_MultiConfig(&mc, args.argv[i+1])
      mc.configs[0].Mapping.sampleId = 0
      mc.configs[1].Mapping.sampleId = 1
      workDual(mc)
    end
    if C.strcmp(args.argv[i], '-i') == 0 and i < args.argc-1 then
      var config : Config
      SCHEMA.parse_Config(&config, args.argv[i+1])
      config.Mapping.sampleId = 0
      workSingle(config)
    end
  end
end

-------------------------------------------------------------------------------
-- COMPILATION CALL
-------------------------------------------------------------------------------

regentlib.saveobj(main, "small.o", "object", MAPPER.register_mappers)
