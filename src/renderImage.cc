
#include <assert.h>
#include <stdio.h>
#include "renderImage.h"

void vDrawScene(int numFluidX,
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


void setCameraPosition(FieldData domainMin[3], FieldData domainMax[3]) {
  // TODO add camera motion
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  FieldData scale[3];
  for(unsigned i = 0; i < 3; ++i) scale[i] = domainMax[i] - domainMin[i];
  FieldData scaleOffset = 0.1;
  glOrtho(domainMin[0], domainMax[0], domainMin[1], domainMax[1], domainMin[2], domainMax[2]);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  gluLookAt(domainMin[0] - scale[0] * scaleOffset,
            domainMax[1] + scale[1] * scaleOffset,
            domainMin[2] - scale[2] * scaleOffset,
            domainMin[0] + domainMax[0] * 0.5,
            domainMax[1],
            domainMin[2] + domainMax[2] * 0.5,
            0.0, 0.0, 1.0);
}


void getLightPosition(FieldData domainMin[3], FieldData domainMax[3], GLfloat lightPosition[3]) {
  lightPosition[0] = 0.5 * (domainMax[0] - domainMin[0]);
  lightPosition[1] = domainMax[1] * 2;
  lightPosition[3] = 0.5 * (domainMax[2] - domainMin[2]);
}


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
                 VisualizationField visualizationField) {
  
  GLfloat lightPosition[3] = { 0 };
  getLightPosition(domainMin, domainMax, lightPosition);
  setCameraPosition(domainMin, domainMax);
  vDrawScene(numFluidX, numFluidY, numFluidZ, rho, pressure, velocity, centerCoordinates, temperature, domainMin, domainMax, visualizationField);
}


void
write_ppm(const char *filename, const GLfloat *rgba, int width, int height)
{
  FILE *f = fopen( filename, "w" );
  if (f) {
    printf ("writing ppm file %s\n", filename);
    int x, y;
    fprintf(f,"P6\n");
    fprintf(f,"# binary ppm file %s\n", filename);
    fprintf(f,"%i %i\n", width,height);
    fprintf(f,"255\n");
    fclose(f);
    
    f = fopen( filename, "ab" );  /* reopen in binary append mode */
    assert(f != NULL);
    for (y = height - 1; y >= 0; y--) {
      unsigned char outputBuffer[width * 3];
      unsigned char *outputPtr = outputBuffer;
      GLfloat* rgbaPtr = (GLfloat*)rgba + (y * width * 4);
      for (x = 0; x < width; x++) {
        int r, g, b;
        r = (int) (rgbaPtr[0] * 255.0);
        g = (int) (rgbaPtr[1] * 255.0);
        b = (int) (rgbaPtr[2] * 255.0);
        if (r > 255) r = 255;
        if (g > 255) g = 255;
        if (b > 255) b = 255;
        outputPtr[0] = r;
        outputPtr[1] = g;
        outputPtr[2] = b;
        outputPtr += 3;
        rgbaPtr += 4;
      }
      fwrite(outputBuffer, 3 * sizeof(unsigned char), width, f);
    }
    fclose(f);
    printf("successfully wrote %s\n", filename);
  } else {
    printf("could not write %s\n", filename);
  }
}
