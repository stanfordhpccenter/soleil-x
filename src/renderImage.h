
#ifndef __renderimage_h__
#define __renderimage_h__

#include "render_standalone.h"


typedef enum {
  rhoField = 1,
  pressureField,
  velocityField,
  temperatureField
} VisualizationField;

void setCameraPosition(FieldData domainMin[3], FieldData domainMax[3]);

void renderInitialize(FieldData domainMin[3], FieldData domainMax[3]);
void renderImage(int numFluidX,
                 int numFluidY,
                 int numFluidZ,
                 FieldData* rho,
                 FieldData* pressure,
                 FieldData* velocity,
                 FieldData* centerCoordinates,
                 FieldData* temperature,
                 FieldData domainMin[3],
                 FieldData domainMax[3],
                 VisualizationField visualizationField);


void write_targa(const char *filename, const GLubyte *rgba, int width, int height);
void write_ppm(const char *filename, const GLubyte *rgba, int width, int height);

#endif // __renderimage_h__

