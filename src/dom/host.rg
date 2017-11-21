-- Runs dom.rg standalone.
-- Compile using: regent.py host.rg -i <config>.lua
-- Run binary using: LD_LIBRARY_PATH=$LEGION_PATH/bindings/terra/ ./a.out

-- Reads the following configuration options from a Soleil-X config file
-- (supplied using the '-i' option):
-- * xnum,ynum,znum : number of cells in fluid grid
-- * coarsenFactor : how coarser the rad. grid is vs the fluid, e.g. {2,2,2}
-- * xWidth,yWidth,zWidth : physical size of the fluid grid
-- * [emiss_west/south/north/up/down : wall emissivity]
-- * [T_west/south/north/up/down : wall temperature]
-- * qa,qs : albedo
-- * numAngles : number of DOM angles

-- Uses default values for:
-- * Ib
-- * sigma

-------------------------------------------------------------------------------
-- Imports
-------------------------------------------------------------------------------

import 'regent'

local A = require 'admiral'

local cmath = terralib.includec("math.h")

-------------------------------------------------------------------------------
-- Read configuration
-------------------------------------------------------------------------------

local config
local i = 1
while i <= #arg do
  if arg[i] == '-i' and i < #arg then
    config = loadfile(arg[i+1])()
    break
  end
  i = i + 1
end
if not config then
  print('config file required (-i <file> option)')
  os.exit(1)
end

-------------------------------------------------------------------------------
-- Proxy radiation grid
-------------------------------------------------------------------------------

struct Point {
    I_1 : double[config.numAngles];
    I_2 : double[config.numAngles];
    I_3 : double[config.numAngles];
    I_4 : double[config.numAngles];
    I_5 : double[config.numAngles];
    I_6 : double[config.numAngles];
    I_7 : double[config.numAngles];
    I_8 : double[config.numAngles];
    Iiter_1 : double[config.numAngles];
    Iiter_2 : double[config.numAngles];
    Iiter_3 : double[config.numAngles];
    Iiter_4 : double[config.numAngles];
    Iiter_5 : double[config.numAngles];
    Iiter_6 : double[config.numAngles];
    Iiter_7 : double[config.numAngles];
    Iiter_8 : double[config.numAngles];
    G : double;
    S : double;
    Ib : double;
    sigma : double;
}

local pointsRel = {}

pointsRel.regionSymbol = terralib.memoize(function(self)
  return regentlib.newsymbol()
end)

pointsRel.primPartSymbol = terralib.memoize(function(self)
  return regentlib.newsymbol()
end)

function pointsRel:regionType()
  return region(ispace(int3d),Point)
end

function pointsRel:Dims()
  return {config.xnum / config.coarsenFactor[1],
          config.ynum / config.coarsenFactor[2],
          config.znum / config.coarsenFactor[3]}
end

-------------------------------------------------------------------------------
-- Import DOM module
-------------------------------------------------------------------------------

local domMod = (require 'dom')(pointsRel)

-------------------------------------------------------------------------------
-- Proxy tasks
-------------------------------------------------------------------------------

local SB = 5.67e-8
local pi = 2.0*cmath.acos(0.0)

local task InitPoints(points : region(ispace(int3d),Point))
where reads writes(points.{I_1, I_2, I_3, I_4, I_5, I_6, I_7, I_8,
                           Iiter_1, Iiter_2, Iiter_3, Iiter_4,
                           Iiter_5, Iiter_6, Iiter_7, Iiter_8,
                           G, S, Ib, sigma}) do
  for p in points do
    for m = 0,config.numAngles do
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
    p.Ib = (SB/pi) * cmath.pow(1000.0, 4.0)
    p.sigma = 5.0
  end
end

local dims = pointsRel:Dims()
local nx,ny,nz = dims[1],dims[2],dims[3]
local ntx,nty,ntz = A.primPartDims()

local points = pointsRel:regionSymbol()
local colors = A.primColors()
local p_points = pointsRel:primPartSymbol()

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

local task main()
  var is = ispace(int3d, {nx,ny,nz})
  var [points] = region(is, Point)
  var [colors] = ispace(int3d, {ntx,nty,ntz})
  var [p_points] = partition(equal, points, colors);
  [domMod.InitModule]
  InitPoints(points);
  [domMod.ComputeRadiationField]
  writeIntensity(points)
end

print('Saving standalone DOM executable to a.out')
regentlib.saveobj(main, 'a.out', 'executable')
