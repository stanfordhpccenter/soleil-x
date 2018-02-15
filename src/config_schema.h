#include "enum.h"

MK_ENUM(FlowBC, Periodic, Symmetry, AdiabaticWall, IsothermalWall);
MK_ENUM(ParticleBC, Permeable, Solid);
MK_ENUM(ViscosityModel, Constant, PowerLaw, Sutherland);
MK_ENUM(FlowInitCase, Uniform, Restart, Perturbed, TaylorGreen2DVortex, TaylorGreen3DVortex);
MK_ENUM(OnOrOff, OFF, ON);
MK_ENUM(ParticlesInitCase, Random, Restart, Uniform);
MK_ENUM(RadiationType, OFF, Algebraic, DOM);

struct Config {
    struct {
        int xNum;
        int yNum;
        int zNum;
        int xTiles;
        int yTiles;
        int zTiles;
        double origin[3];
        double xWidth;
        double yWidth;
        double zWidth;
    } Grid;
    struct {
        FlowBC xBCLeft;
        double xBCLeftVel[3];
        double xBCLeftTemp;
        FlowBC xBCRight;
        double xBCRightVel[3];
        double xBCRightTemp;
        FlowBC yBCLeft;
        double yBCLeftVel[3];
        double yBCLeftTemp;
        FlowBC yBCRight;
        double yBCRightVel[3];
        double yBCRightTemp;
        FlowBC zBCLeft;
        double zBCLeftVel[3];
        double zBCLeftTemp;
        FlowBC zBCRight;
        double zBCRightVel[3];
        double zBCRightTemp;
    } BC;
    struct {
        double finalTime;
        int restartIter;
        int maxIter;
        double cfl;
        double fixedDeltaTime;
    } Integrator;
    struct {
        double gasConstant;
        double gamma;
        double prandtl;
        ViscosityModel viscosityModel;
        double constantVisc;
        double powerlawViscRef;
        double powerlawTempRef;
        double sutherlandViscRef;
        double sutherlandTempRef;
        double sutherlandSRef;
        FlowInitCase initCase;
        double initParams[5];
        double bodyForce[3];
        OnOrOff turbForcing;
    } Flow;
    struct {
        ParticlesInitCase initCase;
        int initNum;
        int maxNum;
        double restitutionCoeff;
        double convectiveCoeff;
        double absorptivity;
        double heatCapacity;
        double initTemperature;
        double density;
        double diameterMean;
        double bodyForce[3];
        double maxSkew;
        int maxXferNum;
    } Particles;
    struct {
        RadiationType type;
        double intensity;
        double qa;
        double qs;
        int xNum;
        int yNum;
        int zNum;
        double emissEast;
        double emissWest;
        double emissSouth;
        double emissNorth;
        double emissUp;
        double emissDown;
        double tempEast;
        double tempWest;
        double tempSouth;
        double tempNorth;
        double tempUp;
        double tempDown;
    } Radiation;
    struct {
        OnOrOff wrtRestart;
        int restartEveryTimeSteps;
        int consoleFrequency;
        int headerFrequency;
    } IO;
};
