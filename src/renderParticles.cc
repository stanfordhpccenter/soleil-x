
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


void drawParticle(GLUquadricObj* qobj, const FieldData3* position, FieldData density, FieldData particleTemperature, float particleSize, FieldData colorScale[2]) {
  
  GLfloat t = particleTemperature;
  GLfloat color[4];
  scaledTemperatureToColor(t, color, colorScale);
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
  return systemScale * 25.0 / density;
}


void renderParticles(int numParticles,
                     const long int* particlesID,
                     const FieldData3* particlesPosition,
                     const FieldData* particlesTemperature,
                     const FieldData* particlesDensity,
                     int numParticlesToDraw,
                     long int* particlesToDraw,
                     float systemScale,
                     FieldData colorScale[2]) {
  GLUquadricObj *qobj = gluNewQuadric();
  
  unsigned drawnCount = 0;
  for(int i = 0; i < numParticles; ++i) {
    if(drawThis(particlesID[i], numParticlesToDraw, particlesToDraw)) {
      if(particlesDensity[i] > 0) {
        drawParticle(qobj, particlesPosition + i, particlesDensity[i], particlesTemperature[i], particleSize(systemScale, particlesDensity[i]), colorScale);
        drawnCount++;
      }
    }
  }
  printf("drew %d particles\n", drawnCount);
  gluDeleteQuadric(qobj);
  
}
