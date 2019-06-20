
#ifndef __renderimage_h__
#define __renderimage_h__

#include "render_standalone.h"

typedef double FieldData;
typedef struct {
  FieldData x[3];
} FieldData3;

typedef enum {
  rhoField = 0,
  pressureField,
  velocityField,
  temperatureField
} VisualizationField;

void setCameraPosition(FieldData domainMin[3], FieldData domainMax[3]);

void renderInitialize(FieldData domainMin[3], FieldData domainMax[3],
                      OSMesaContext& mesaCtx, GLubyte*& rgbaBuffer, GLfloat*& depthBuffer);
void renderImage(int numFluidX,
                 int numFluidY,
                 int numFluidZ,
                 const FieldData* rho,
                 const FieldData* pressure,
                 const FieldData3* velocity,
                 const FieldData3* centerCoordinates,
                 const FieldData* temperature,
                 FieldData domainMin[3],
                 FieldData domainMax[3],
                 VisualizationField visualizationField,
                 FieldData targetValue,
                 FieldData isosurfaceScale[2],
                 int numParticles,
                 const long int* particlesID,
                 const FieldData3* particlesPosition,
                 const FieldData* particlesTemperature,
                 const FieldData* particlesDensity,
                 long int* particlesToDraw,
                 int numParticlesToDraw);
void renderTerminate(OSMesaContext mesaCtx, GLubyte*& rgbaBuffer, GLfloat*& depthBuffer);


void write_targa(const char *filename, const GLubyte *rgba, int width, int height);
void write_ppm(const char *filename, const GLubyte *rgba, int width, int height);

#endif // __renderimage_h__

