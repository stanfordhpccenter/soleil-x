#include <stdbool.h>
#include <stdint.h>

typedef int ParticlesInitCase;
#define ParticlesInitCase_Random 0
#define ParticlesInitCase_Restart 1
#define ParticlesInitCase_Uniform 2

typedef int TurbForcingModel;
#define TurbForcingModel_OFF 0
#define TurbForcingModel_HIT 1

typedef int FlowBC;
#define FlowBC_IsothermalWall 3
#define FlowBC_NonUniformTemperatureWall 6
#define FlowBC_Periodic 0
#define FlowBC_NSCBC_SubsonicOutflow 5
#define FlowBC_AdiabaticWall 2
#define FlowBC_Symmetry 1
#define FlowBC_NSCBC_SubsonicInflow 4

typedef int ParticlesBC;
#define ParticlesBC_Bounce 1
#define ParticlesBC_Disappear 2
#define ParticlesBC_Periodic 0

typedef int FeedModel;
#define FeedModel_OFF 0
#define FeedModel_Incoming 1

typedef int ViscosityModel;
#define ViscosityModel_PowerLaw 1
#define ViscosityModel_Constant 0
#define ViscosityModel_Sutherland 2

typedef int FlowInitCase;
#define FlowInitCase_Perturbed 3
#define FlowInitCase_Uniform 0
#define FlowInitCase_Random 1
#define FlowInitCase_TaylorGreen3DVortex 5
#define FlowInitCase_TaylorGreen2DVortex 4
#define FlowInitCase_Restart 2

typedef int TempProfile;
#define TempProfile_Incoming 0
#define TempProfile_Constant 1
#define TempProfile_Parabola 2

typedef int InflowProfile;
#define InflowProfile_Incoming 0
#define InflowProfile_Duct 1
#define InflowProfile_Constant 2

typedef int RadiationModel;
#define RadiationModel_DOM 0
#define RadiationModel_Algebraic 1
#define RadiationModel_OFF 2

struct Volume {
  int32_t uptoCell[3];
  int32_t fromCell[3];
};

struct Window {
  int32_t uptoCell[2];
  int32_t fromCell[2];
};

struct Config {
  struct  {
    bool collisions;
    double convectiveCoeff;
    double diameterMean;
    double escapeRatioPerDir;
    double restitutionCoeff;
    double density;
    struct  {
      int32_t type;
      union  {
        struct  {
        } OFF;
        struct  {
          double addedVelocity[3];
        } Incoming;
      } u;
    } feeding;
    int32_t maxNum;
    double bodyForce[3];
    double maxSkew;
    int8_t restartDir[256];
    int32_t initNum;
    int32_t initCase;
    double initTemperature;
    double heatCapacity;
  } Particles;
  struct  {
    struct  {
      uint32_t length;
      struct Volume values[5];
    } probes;
    int32_t restartEveryTimeSteps;
    bool wrtRestart;
  } IO;
  struct  {
    int32_t restartIter;
    double restartTime;
    double fixedDeltaTime;
    double cfl;
    int32_t maxIter;
  } Integrator;
  struct  {
    int32_t initCase;
    double powerlawViscRef;
    double sutherlandTempRef;
    int8_t restartDir[256];
    struct  {
      int32_t type;
      union  {
        struct  {
        } OFF;
        struct  {
          double K_o;
          double t_o;
          double G;
        } HIT;
      } u;
    } turbForcing;
    double prandtl;
    double sutherlandSRef;
    double initParams[5];
    double bodyForce[3];
    double sutherlandViscRef;
    int32_t viscosityModel;
    double powerlawTempRef;
    double constantVisc;
    double gasConstant;
    double gamma;
  } Flow;
  struct  {
    int32_t tiles[3];
    int32_t tilesPerRank[3];
    int8_t outDir[256];
    int32_t sampleId;
    int32_t wallTime;
  } Mapping;
  struct  {
    double xBCLeftVel[3];
    double yBCLeftVel[3];
    struct  {
      int32_t type;
      union  {
        struct  {
          double addedVelocity;
        } Incoming;
        struct  {
          double meanVelocity;
        } Duct;
        struct  {
          double velocity;
        } Constant;
      } u;
    } xBCLeftInflowProfile;
    double xBCRightVel[3];
    int32_t xBCRight;
    int32_t zBCRight;
    int32_t yBCRight;
    struct  {
      int32_t type;
      union  {
        struct  {
        } Incoming;
        struct  {
          double temperature;
        } Constant;
        struct  {
          double T_right;
          double T_left;
          double T_mid;
        } Parabola;
      } u;
    } zBCRightHeat;
    int32_t zBCLeft;
    struct  {
      int32_t type;
      union  {
        struct  {
        } Incoming;
        struct  {
          double temperature;
        } Constant;
        struct  {
          double T_right;
          double T_left;
          double T_mid;
        } Parabola;
      } u;
    } yBCRightHeat;
    int32_t xBCLeft;
    struct  {
      int32_t type;
      union  {
        struct  {
        } Incoming;
        struct  {
          double temperature;
        } Constant;
        struct  {
          double T_right;
          double T_left;
          double T_mid;
        } Parabola;
      } u;
    } zBCLeftHeat;
    double xBCRightP_inf;
    struct  {
      int32_t type;
      union  {
        struct  {
        } Incoming;
        struct  {
          double temperature;
        } Constant;
        struct  {
          double T_right;
          double T_left;
          double T_mid;
        } Parabola;
      } u;
    } xBCLeftHeat;
    double yBCRightVel[3];
    double zBCRightVel[3];
    struct  {
      int32_t type;
      union  {
        struct  {
        } Incoming;
        struct  {
          double temperature;
        } Constant;
        struct  {
          double T_right;
          double T_left;
          double T_mid;
        } Parabola;
      } u;
    } yBCLeftHeat;
    int32_t yBCLeft;
    double zBCLeftVel[3];
    struct  {
      int32_t type;
      union  {
        struct  {
        } Incoming;
        struct  {
          double temperature;
        } Constant;
        struct  {
          double T_right;
          double T_left;
          double T_mid;
        } Parabola;
      } u;
    } xBCRightHeat;
  } BC;
  struct  {
    double xWidth;
    double zWidth;
    double origin[3];
    int32_t yNum;
    double yWidth;
    int32_t xNum;
    int32_t zNum;
  } Grid;
  struct  {
    int32_t type;
    union  {
      struct  {
        double yLoTemp;
        struct Window zLoWindow;
        double qa;
        double zLoIntensity;
        double yLoIntensity;
        struct Window zHiWindow;
        struct Window yLoWindow;
        struct Window yHiWindow;
        struct Window xHiWindow;
        int32_t yNum;
        double yHiTemp;
        double yLoEmiss;
        double xHiTemp;
        struct Window xLoWindow;
        double zHiIntensity;
        double yHiIntensity;
        int32_t angles;
        int32_t xNum;
        double xLoIntensity;
        double xLoTemp;
        double zLoTemp;
        double zHiTemp;
        int32_t zNum;
        double xHiEmiss;
        double xLoEmiss;
        double zHiEmiss;
        double yHiEmiss;
        double zLoEmiss;
        double qs;
        double xHiIntensity;
      } DOM;
      struct  {
        double absorptivity;
        double intensity;
      } Algebraic;
      struct  {
      } OFF;
    } u;
  } Radiation;
};
