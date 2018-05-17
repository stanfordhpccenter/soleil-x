-- Runs dom.rg standalone.
-- Reads configuration options in the same format as main simulation.
-- Uses default values for Ib and sigma.

-------------------------------------------------------------------------------
-- Imports
-------------------------------------------------------------------------------

import 'regent'

local C = regentlib.c
local SCHEMA = terralib.includec("config_schema.h")

-------------------------------------------------------------------------------
-- Compile-time configuration options
-------------------------------------------------------------------------------

local NUM_ANGLES = 14

-------------------------------------------------------------------------------
-- Proxy radiation grid
-------------------------------------------------------------------------------

struct Point {
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
}

-------------------------------------------------------------------------------
-- Import DOM module
-------------------------------------------------------------------------------

local DOM = (require 'dom')(NUM_ANGLES, Point, SCHEMA.Config)

-------------------------------------------------------------------------------
-- Proxy tasks
-------------------------------------------------------------------------------

local SB = 5.67e-8
local pi = 3.1415926535898
local pow = regentlib.pow(double)

local task InitPoints(points : region(ispace(int3d),Point))
where reads writes(points.{I_1, I_2, I_3, I_4, I_5, I_6, I_7, I_8,
                           Iiter_1, Iiter_2, Iiter_3, Iiter_4,
                           Iiter_5, Iiter_6, Iiter_7, Iiter_8,
                           G, S, Ib, sigma}) do
  for p in points do
    for m = 0,NUM_ANGLES do
      p.I_1[m]     = 0.0
      p.I_2[m]     = 0.0
      p.I_3[m]     = 0.0
      p.I_4[m]     = 0.0
      p.I_5[m]     = 0.0
      p.I_6[m]     = 0.0
      p.I_7[m]     = 0.0
      p.I_8[m]     = 0.0
      p.Iiter_1[m] = 0.0
      p.Iiter_2[m] = 0.0
      p.Iiter_3[m] = 0.0
      p.Iiter_4[m] = 0.0
      p.Iiter_5[m] = 0.0
      p.Iiter_6[m] = 0.0
      p.Iiter_7[m] = 0.0
      p.Iiter_8[m] = 0.0
    end
    p.G = 0.0
    p.S = 0.0
    p.Ib = (SB/pi) * pow(1000.0, 4.0)
    p.sigma = 5.0
  end
end

local task writeIntensity(points : region(ispace(int3d), Point))
where
  reads (points.G)
do
  var limits = points.bounds
  var f = C.fopen("dom_intensities/intensity.dat", "w")
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

local task main()
  -- Read configuration
  var args = C.legion_runtime_get_input_args()
  var stderr = C.fdopen(2, 'w')
  if args.argc < 2 then
    C.fprintf(stderr, "Usage: %s config.json\n", args.argv[0])
    C.fflush(stderr)
    C.exit(1)
  end
  var config = SCHEMA.parse_config(args.argv[1])

  if config.Radiation.angles ~= NUM_ANGLES then
    C.fprintf(stderr, "angles in config file (%d) must match NUM_ANGLES in dom_host.rg (%d)\n",
      config.Radiation.angles, NUM_ANGLES)
    C.fflush(stderr)
    C.exit(1)
  end
  -- Initialize symbols
  var is = ispace(int3d, {config.Radiation.xNum, config.Radiation.yNum, config.Radiation.zNum})
  var points = region(is, Point)
  var colors = ispace(int3d, {config.Mapping.xTiles, config.Mapping.yTiles, config.Mapping.zTiles})
  var p_points = partition(equal, points, colors);
  -- Inline quotes from external module
  [DOM.DeclSymbols(config)];
  [DOM.InitRegions()];
  for color in colors do
    InitPoints(p_points[color])
  end
  [DOM.ComputeRadiationField(config, colors, p_points)];
  writeIntensity(points)
end

regentlib.saveobj(main, 'dom_host.o', 'object')
