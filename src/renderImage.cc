

#include "renderImage.h"

void vDrawScene();


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


void renderImage(int numLines,
                 FieldData* rho,
                 FieldData* pressure,
                 FieldData* velocity,
                 FieldData* centerCoordinates,
                 FieldData* temperature,
                 FieldData domainMin[3],
                 FieldData domainMax[3]) {
  
  GLfloat lightPosition[3] = { 0 };
  getLightPosition(domainMin, domainMax, lightPosition);
  setCameraPosition(domainMin, domainMax);
  vDrawScene();
}
