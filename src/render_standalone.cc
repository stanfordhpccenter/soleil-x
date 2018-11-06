
//
// This is a standalone renderer, used to develop the visualization algorithm offline.
// Keep this file for future development.

#include <iostream>
#include <stdio.h>

#include "render_standalone.h"
#include "renderImage.h"


void initializeMarchingCubes();


void createGraphicsContext(OSMesaContext &mesaCtx,
                           GLfloat* &rgbaBuffer,
                           GLfloat* &depthBuffer) {
#if OSMESA_MAJOR_VERSION * 100 + OSMESA_MINOR_VERSION >= 305
  /* specify Z, stencil, accum sizes */
  mesaCtx = OSMesaCreateContextExt(GL_RGBA, 32, 0, 0, NULL);
#else
  mesaCtx = OSMesaCreateContext(GL_RGBA, NULL);
#endif
  if (!mesaCtx) {
    printf("OSMesaCreateContext failed!\n");
    return;
  }
  /* Allocate the image buffer */
  const int fieldsPerPixel = 4;
  rgbaBuffer = (GLfloat *) malloc(WIDTH * HEIGHT * fieldsPerPixel * sizeof(GLfloat));
  if (!rgbaBuffer) {
    printf("Alloc image buffer failed!\n");
    return;
  }
  
  /* Bind the buffer to the context and make it current */
  if (!OSMesaMakeCurrent(mesaCtx, rgbaBuffer, GL_FLOAT, WIDTH, HEIGHT)) {
    printf("OSMesaMakeCurrent failed!\n");
    return;
  }
  
  /* Allocate the depth buffer. */
  depthBuffer = (GLfloat*)malloc(WIDTH * HEIGHT * sizeof(GLfloat));
  if (!depthBuffer) {
    printf("Alloc depth buffer failed!\n");
    return;
  }

}


static void destroyGraphicsContext(OSMesaContext mesaCtx) {
  /* destroy the context */
  OSMesaDestroyContext(mesaCtx);
}


void loadFluidData(char* fileName,
                   int numLines,
                   FieldData* rho,
                   FieldData* pressure,
                   FieldData* velocity,
                   FieldData* centerCoordinates,
                   FieldData* temperature,
                   FieldData domainMin[3],
                   FieldData domainMax[3]) {
  FILE* fluidIn = fopen(fileName, "r");
  for(int i = 0; i < numLines; ++i) {
    int ret = fscanf(fluidIn, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
                     rho + i, pressure + i,
                     velocity + 3 * i, velocity + 3 * i + 1, velocity + 3 * i + 2,
                     centerCoordinates + 3 * i, centerCoordinates + 3 * i + 1, centerCoordinates + 3 * i + 2,
                     temperature + i);
    if(ret == EOF) {
      std::cerr << "error reading fluid data file" << std::endl;
    }
    for(int j = 0; j < 3; ++j) {
      domainMin[j] = std::min(domainMin[j], centerCoordinates[3 * i + j]);
      domainMax[j] = std::max(domainMax[j], centerCoordinates[3 * i + j]);
    }
  }
}


void retrieveRenderedImage(GLfloat* rgbaBuffer,
                           GLfloat* depthBuffer) {
  glReadPixels(0, 0, WIDTH, HEIGHT, GL_RGBA, GL_FLOAT, rgbaBuffer);
  glReadPixels(0, 0, WIDTH, HEIGHT, GL_DEPTH_COMPONENT, GL_FLOAT, depthBuffer);
}


void saveImageToFile(char* fluidFileName,
                     GLfloat* rgbaBuffer) {
  // base output file name on fluidFileName
  char filename[256] = { 0 };
  sprintf(filename, "%s.ppm", fluidFileName);
  write_ppm(filename, rgbaBuffer, WIDTH, HEIGHT);
}


int main(int argc, char **argv) {
  
  if(argc < 2) {
    std::cerr << "missing number of points in fluid file string in the form XxYxZ where X,Y,Z are integers" << std::endl;
    return -1;
  }
  int numFluidLines = 0;
  int numFluidX = 0;
  int numFluidY = 0;
  int numFluidZ = 0;
  sscanf(argv[1], "%dx%dx%d", &numFluidX, &numFluidY, &numFluidZ);
  numFluidLines = numFluidX * numFluidY * numFluidZ;
  
  if(argc < 3) {
    std::cerr << "missing name of fluid file" << std::endl;
  }
  char* fluidFileName = argv[2];
  
  FieldData* rho = new FieldData[numFluidLines];
  FieldData* pressure = new FieldData[numFluidLines];
  FieldData* velocity = new FieldData[numFluidLines * 3];
  FieldData* centerCoordinates = new FieldData[numFluidLines * 3];
  FieldData* temperature = new FieldData[numFluidLines];
  FieldData domainMin[3] = { 99999, 99999, 99999 };
  FieldData domainMax[3] = { 0 };
  
  loadFluidData(fluidFileName, numFluidLines, rho, pressure, velocity, centerCoordinates,  temperature, domainMin, domainMax);
  
  OSMesaContext mesaCtx;
  GLfloat* rgbaBuffer;
  GLfloat* depthBuffer;
  createGraphicsContext(mesaCtx, rgbaBuffer, depthBuffer);
  
  initializeMarchingCubes();
  
  renderImage(numFluidX, numFluidY, numFluidZ, rho, pressure, velocity, centerCoordinates, temperature, domainMin, domainMax, temperatureField);
  
  retrieveRenderedImage(rgbaBuffer, depthBuffer);
  
  saveImageToFile(fluidFileName, rgbaBuffer);
  
  destroyGraphicsContext(mesaCtx);
}


