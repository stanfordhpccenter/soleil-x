
#define USE_SOFTWARE_OPENGL

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include "OpenGL/glu.h"
#else
#ifdef USE_SOFTWARE_OPENGL
#include "GL/osmesa.h"
#else
#include "vtk_glew.h"

#endif
#include "GL/glu.h"
#endif

#include "renderImage.h"
#include <stdio.h>
#include <iostream>

typedef double FieldData;


static void temperatureToColor(float temperature,
                               float color[4]) {
  
  // https://stackoverflow.com/questions/7229895/display-temperature-as-a-color-with-c
  
  float x = (float)(temperature / 1000.0);
  float x2 = x * x;
  float x3 = x2 * x;
  float x4 = x3 * x;
  float x5 = x4 * x;
  
  float R, G, B = 0.0f;
  
  // red
  if (temperature <= 6600)
    R = 1.0f;
  else
    R = 0.0002889f * x5 - 0.01258f * x4 + 0.2148f * x3 - 1.776f * x2 + 6.907f * x - 8.723f;
  
  // green
  if (temperature <= 6600)
    G = -4.593e-05f * x5 + 0.001424f * x4 - 0.01489f * x3 + 0.0498f * x2 + 0.1669f * x - 0.1653f;
  else
    G = -1.308e-07f * x5 + 1.745e-05f * x4 - 0.0009116f * x3 + 0.02348f * x2 - 0.3048f * x + 2.159f;
  
  // blue
  if (temperature <= 2000)
    B = 0.0f;
  else if (temperature < 6600)
    B = 1.764e-05f * x5 + 0.0003575f * x4 - 0.01554f * x3 + 0.1549f * x2 - 0.3682f * x + 0.2386f;
  else
    B = 1.0f;
  
  color[0] = R;
  color[1] = G;
  color[2] = B;
  color[3] = 1.0f;
}


static void scaledTemperatureToColor(float temperature,
                                     float color[4],
                                     FieldData isosurfaceScale[2]) {
  const float min = isosurfaceScale[0];
  const float max = isosurfaceScale[1];
  const float Kmin = 0.0f;
  const float Kmax = 10000.0f;
  // stretch it on a scale of Kmin...Kmax
  float scaledTemperature = (temperature - min) * ((Kmax - Kmin) / (max - min));
  return temperatureToColor(scaledTemperature, color);
}


void drawParticle(GLUquadricObj* qobj, const FieldData3* position, FieldData density, FieldData particleTemperature, float particleSize, FieldData isosurfaceScale[2]) {
  
  GLfloat t = particleTemperature;
  GLfloat color[4];
  scaledTemperatureToColor(t, color, isosurfaceScale);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, color);
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color);
  
#if 0
  std::cout << "particle at " << position->x[0] << " " << position->x[1] << " " << position->x[2] << std::endl;
#endif
  glPushMatrix();
  glTranslatef(position->x[0], position->x[1], position->x[2]);
  gluSphere(qobj, particleSize, 7, 7);
  glPopMatrix();
}


static bool drawThis(long int id, int numParticlesToDraw, long int* particlesToDraw) {
  int first = 0;
  int last = numParticlesToDraw - 1;
  int middle = (first + last) / 2;
  
  while (first <= last) {
    if (particlesToDraw[middle] < id) first = middle + 1;
    else if (particlesToDraw[middle] == id) return true;
    else last = middle - 1;
    middle = (first + last) / 2;
  }
  if (first > last) return false;
  return true;
}


static float particleSize(float systemScale, FieldData density) {
  // 0.0002 at density 8900 systemscale .023 numParticles 50 should be 500
  // systemScale*100/density
  return systemScale * 100.0 / density;
}


void renderParticles(int numParticles,
                     const long int* particlesID,
                     const FieldData3* particlesPosition,
                     const FieldData* particlesTemperature,
                     const FieldData* particlesDensity,
                     int numParticlesToDraw,
                     long int* particlesToDraw,
                     float systemScale,
                     FieldData isosurfaceScale[2]) {
  GLUquadricObj *qobj = gluNewQuadric();
  
  unsigned drawnCount = 0;
  for(int i = 0; i < numParticles; ++i) {
    if(drawThis(particlesID[i], numParticlesToDraw, particlesToDraw)) {
      if(particlesDensity[i] > 0) {
        drawParticle(qobj, particlesPosition + i, particlesDensity[i], particlesTemperature[i], particleSize(systemScale, particlesDensity[i]), isosurfaceScale);
        drawnCount++;
      }
    }
  }
  printf("drew %d particles\n", drawnCount);
  gluDeleteQuadric(qobj);
  
}
