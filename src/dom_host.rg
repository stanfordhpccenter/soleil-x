-- Runs dom.rg standalone.
-- Reads configuration options in the same format as main simulation.
-- Uses default values for Ib and sigma.

-------------------------------------------------------------------------------
-- Imports
-------------------------------------------------------------------------------

import 'regent'

local C = regentlib.c
local SCHEMA = terralib.includec("config_schema.h")
local UTIL = require 'util'

-------------------------------------------------------------------------------
-- Compile-time configuration options
-------------------------------------------------------------------------------

local MAX_ANGLES_PER_QUAD = 44

-------------------------------------------------------------------------------
-- Proxy radiation grid
-------------------------------------------------------------------------------

struct Point {
    I_1 : regentlib.array(double, MAX_ANGLES_PER_QUAD);
    I_2 : regentlib.array(double, MAX_ANGLES_PER_QUAD);
    I_3 : regentlib.array(double, MAX_ANGLES_PER_QUAD);
    I_4 : regentlib.array(double, MAX_ANGLES_PER_QUAD);
    I_5 : regentlib.array(double, MAX_ANGLES_PER_QUAD);
    I_6 : regentlib.array(double, MAX_ANGLES_PER_QUAD);
    I_7 : regentlib.array(double, MAX_ANGLES_PER_QUAD);
    I_8 : regentlib.array(double, MAX_ANGLES_PER_QUAD);
    G : double;
    S : double;
    Ib : double;
    sigma : double;
}

-------------------------------------------------------------------------------
-- Import DOM module
-------------------------------------------------------------------------------

local DOM = (require 'dom-desugared')(MAX_ANGLES_PER_QUAD, Point, SCHEMA.Config)
local DOM_INST = DOM.mkInstance()

-------------------------------------------------------------------------------
-- Proxy tasks
-------------------------------------------------------------------------------

local SB = 5.67e-8
local PI = 3.1415926535898
local pow = regentlib.pow(double)

local task InitPoints(points : region(ispace(int3d),Point))
where
  writes(points.{G, S, Ib, sigma}),
  reads writes(points.{I_1, I_2, I_3, I_4, I_5, I_6, I_7, I_8})
do
  for p in points do
    for m = 0, MAX_ANGLES_PER_QUAD do
      p.I_1[m] = 0.0
      p.I_2[m] = 0.0
      p.I_3[m] = 0.0
      p.I_4[m] = 0.0
      p.I_5[m] = 0.0
      p.I_6[m] = 0.0
      p.I_7[m] = 0.0
      p.I_8[m] = 0.0
    end
    p.G = 0.0
    p.S = 0.0
    p.Ib = (SB/PI) * pow(1000.0, 4.0)
    p.sigma = 5.0
  end
end

local task writeIntensity(points : region(ispace(int3d), Point))
where
  reads(points.G)
do
  var limits = points.bounds
  var f = UTIL.openFile("intensity.dat", "w")
  for i = limits.lo.x, limits.hi.x+1 do
    for j = limits.lo.y, limits.hi.y+1 do
      for k = limits.lo.z, limits.hi.z+1 do
        C.fprintf(f,' %.15e \n', points[{i,j,k}].G)
      end
      C.fprintf(f,'\n')
    end
    C.fprintf(f,'\n')
  end
  C.fclose(f)
end

-------------------------------------------------------------------------------
-- Proxy main
-------------------------------------------------------------------------------

local __forbid(__optimize)
task main()
  -- Read configuration
  var args = C.legion_runtime_get_input_args()
  var stderr = C.fdopen(2, 'w')
  if args.argc < 2 then
    C.fprintf(stderr, "Usage: %s config.json\n", args.argv[0])
    C.fflush(stderr)
    C.exit(1)
  end
  var config : SCHEMA.Config[1]
  SCHEMA.parse_Config([&SCHEMA.Config](config), args.argv[1])
  var Nx = config[0].Radiation.xNum
  var Ny = config[0].Radiation.yNum
  var Nz = config[0].Radiation.zNum
  var ntx = config[0].Mapping.tiles[0]
  var nty = config[0].Mapping.tiles[1]
  var ntz = config[0].Mapping.tiles[2]
  regentlib.assert(Nx % ntx == 0, "Uneven partitioning of radiation grid on x")
  regentlib.assert(Ny % nty == 0, "Uneven partitioning of radiation grid on y")
  regentlib.assert(Nz % ntz == 0, "Uneven partitioning of radiation grid on z")
  var is = ispace(int3d, {Nx,Ny,Nz})
  var points = region(is, Point)
  var tiles = ispace(int3d, {ntx,nty,ntz})
  var coloring = regentlib.c.legion_domain_point_coloring_create()
  for c in tiles do
    var rect = rect3d{lo = int3d{(Nx/ntx)*c.x,       (Ny/nty)*c.y,       (Nz/ntz)*c.z      },
                      hi = int3d{(Nx/ntx)*(c.x+1)-1, (Ny/nty)*(c.y+1)-1, (Nz/ntz)*(c.z+1)-1}}
    regentlib.c.legion_domain_point_coloring_color_domain(coloring, c, rect)
  end
  var p_points = partition(disjoint, points, coloring, tiles)
  regentlib.c.legion_domain_point_coloring_destroy(coloring);
  -- Inline quotes from external module
  [DOM_INST.DeclSymbols(rexpr config[0] end)];
  [DOM_INST.InitRegions(rexpr config[0] end)];
  for c in tiles do
    InitPoints(p_points[c])
  end
  [DOM_INST.ComputeRadiationField(rexpr config[0] end, tiles, p_points)];
  writeIntensity(points)
end

regentlib.saveobj(main, 'dom_host.o', 'object')
