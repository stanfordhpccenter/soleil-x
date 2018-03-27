-- Runs dom.rg standalone.
-- Uses static values for DOM configuration options.
-- Uses default values for Ib and sigma.

-------------------------------------------------------------------------------
-- Imports
-------------------------------------------------------------------------------

import 'regent'

-------------------------------------------------------------------------------
-- Proxy runtime configuration object
-------------------------------------------------------------------------------

struct Config {
  Mapping : struct {
    xTiles : int,
    yTiles : int,
    zTiles : int,
  },
  Grid : struct {
    xWidth : double,
    yWidth : double,
    zWidth : double,
  },
  Radiation : struct {
    xNum : int,
    yNum : int,
    zNum : int,
    qa : double,
    qs : double,
    emissWest  : double,
    emissEast  : double,
    emissSouth : double,
    emissNorth : double,
    emissUp    : double,
    emissDown  : double,
    tempWest  : double,
    tempEast  : double,
    tempSouth : double,
    tempNorth : double,
    tempUp    : double,
    tempDown  : double,
  },
}

-------------------------------------------------------------------------------
-- Set configuration options statically
-------------------------------------------------------------------------------

local terra readConfig() : Config
  var config : Config
  config.Mapping.xTiles = 2
  config.Mapping.yTiles = 2
  config.Mapping.zTiles = 2
  config.Grid.xWidth = 1.0
  config.Grid.yWidth = 1.0
  config.Grid.zWidth = 1.0
  config.Radiation.xNum = 6
  config.Radiation.yNum = 6
  config.Radiation.zNum = 6
  config.Radiation.qa = 0.5
  config.Radiation.qs = 0.5
  config.Radiation.emissWest  = 1.0
  config.Radiation.emissEast  = 1.0
  config.Radiation.emissSouth = 1.0
  config.Radiation.emissNorth = 1.0
  config.Radiation.emissUp    = 1.0
  config.Radiation.emissDown  = 1.0
  config.Radiation.tempWest  = 2000.0
  config.Radiation.tempEast  = 300.0
  config.Radiation.tempSouth = 300.0
  config.Radiation.tempNorth = 300.0
  config.Radiation.tempUp    = 300.0
  config.Radiation.tempDown  = 300.0
  return config
end

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

local domMod = (require 'dom')(NUM_ANGLES, Point)

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
  var f = regentlib.c.fopen("intensity.dat", "w")
  for i = limits.lo.x, limits.hi.x+1 do
    for j = limits.lo.y, limits.hi.y+1 do
      for k = limits.lo.z, limits.hi.z+1 do
        regentlib.c.fprintf(f,' %.6e ', points[{i,j,k}].G)
      end
      regentlib.c.fprintf(f,'\n')
    end
    regentlib.c.fprintf(f,'\n')
  end
  regentlib.c.fclose(f)
end

-------------------------------------------------------------------------------
-- Proxy main
-------------------------------------------------------------------------------

local task main()
  -- "Read" configuration
  var config = readConfig()
  -- Initialize symbols
  var is = ispace(int3d, {config.Radiation.xNum, config.Radiation.yNum, config.Radiation.zNum})
  var points = region(is, Point)
  var colors = ispace(int3d, {config.Mapping.xTiles, config.Mapping.yTiles, config.Mapping.zTiles})
  var p_points = partition(equal, points, colors);
  -- Inline quotes from external module
  [domMod.DeclSymbols(config)];
  [domMod.InitRegions()];
  InitPoints(points);
  [domMod.ComputeRadiationField(config, colors, p_points)];
  writeIntensity(points)
end

print('Saving standalone DOM executable to a.out')
regentlib.saveobj(main, 'a.out', 'executable', nil, {'-lm'})
