/* Copyright 2017 Stanford University
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef __ENSEMBLE_TYPES_H__
#define __ENSEMBLE_TYPES_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef struct IO
{
  int     headerFrequency;
  int     consoleFrequency;
  int     restartEveryTimeSteps;
  int     wrtRestart;
} IO;

typedef struct Integrator
{
  double  finalTime;
  int     restartIter;
  double  fixedDeltaTime;
  double  cfl;
  int     maxIter;
} Integrator;

typedef struct Flow
{
  int     initCase;
  double  powerlawViscRef;
  double  sutherlandTempRef;
  int     turbForcing;
  double  prandtl;
  double  sutherlandSRef;
  double  initParams[5];
  double  bodyForce[3];
  double  sutherlandViscRef;
  int     viscosityModel;
  double  powerlawTempRef;
  double  constantVisc;
  double  gasConstant;
  double  gamma;
} Flow;

typedef struct Particles
{
  double  convectiveCoeff;
  double  diameterMean;
  double  absorptivity;
  double  restitutionCoeff;
  double  density;
  int     maxNum;
  double  bodyForce[3];
  int     maxXferNum;
  double  maxSkew;
  int     initNum;
  int     initCase;
  double  initTemperature;
  double  heatCapacity;
} Particles;

typedef struct BC
{
  double  xBCLeftVel[3];
  double  yBCLeftVel[3];
  double  xBCRightVel[3];
  int     xBCRight;
  double  xBCRightTemp;
  int     zBCRight;
  int     yBCRight;
  double  zBCLeftTemp;
  int     xBCLeft;
  double  xBCLeftTemp;
  double  zBCRightTemp;
  int     yBCLeft;
  double  yBCRightVel[3];
  double  zBCRightVel[3];
  double  yBCLeftTemp;
  double  yBCRightTemp;
  double  zBCLeftVel[3];
  int     zBCLeft;
} BC;

typedef struct Grid
{
  int     xNum;
  int     yNum;
  int     zNum;
  int     yTiles;
  int     zTiles;
  int     xTiles;
  double  origin[3];
  double  xWidth;
  double  yWidth;
  double  zWidth;
} Grid;

typedef struct Radiation
{
  double intensity;
  double emissSouth;
  double tempSouth;
  double qa;
  int yNum;
  double tempNorth;
  double emissWest;
  double emissEast;
  double tempWest;
  int xNum;
  double emissDown;
  double emissUp;
  double tempDown;
  int zNum;
  double emissNorth;
  double tempEast;
  double qs;
  double tempUp;
} Radiation;

typedef struct Config
{
  int        unique_id;
  Grid       grid;
  BC         bc;
  Integrator integrator;
  Flow       flow;
  Particles  particles;
  Radiation  radiation;
  IO         io;
  char       filename[512];
} Config;

enum Tunables
{
  TUNABLE_CONFIG = 0xCAFE,
};

#ifdef __cplusplus
}
#endif

#endif // __ENSEMBLE_TYPES_H__
