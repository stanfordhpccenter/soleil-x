
#include <iostream>
#include <stdio.h>

#include "render_standalone.h"
#include "renderImage.h"





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
  fclose(fluidIn);
}


void loadParticlesData(char* filename,
                       int numParticles,
                       long int* particlesID,
                       FieldData* particlesPosition,
                       FieldData* particlesTemperature,
                       FieldData* particlesDensity) {
  FieldData particlesMin[3] = { 9999, 9999, 9999 };
  FieldData particlesMax[3] = { -9999, -9999, -9999 };
  FILE* particlesIn = fopen(filename, "r");
  for(int i = 0; i < numParticles; ++i) {
    int ret = fscanf(particlesIn, "%ld %lf %lf %lf %lf %lf",
                     particlesID + i,
                     particlesPosition + 3 * i,
                     particlesPosition + 3 * i + 1,
                     particlesPosition + 3 * i + 2,
                     particlesTemperature + i,
                     particlesDensity + i);
    for(int j = 0; j < 3; ++j) {
      particlesMin[j] = std::min(particlesMin[j], particlesPosition[3 * i + j]);
      particlesMax[j] = std::max(particlesMax[j], particlesPosition[3 * i + j]);
    }
    if(ret == EOF) {
      std::cerr << "error reading particles file" << std::endl;
    }
  }
  fclose(particlesIn);
  std::cout << "particles min " << particlesMin[0] << " " << particlesMin[1] << " " << particlesMin[2] << std::endl;
  std::cout << "particles max " << particlesMax[0] << " " << particlesMax[1] << " " << particlesMax[2] << std::endl;
}



void saveImageToFile(char* fluidFileName,
                     GLubyte* rgbaBuffer) {
  // base output file name on fluidFileName
  char filename[256] = { 0 };
  sprintf(filename, "%s.tga", fluidFileName);
  write_targa(filename, rgbaBuffer, WIDTH, HEIGHT);
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
    return -1;
  }
  char* fluidFileName = argv[2];
  
  if(argc < 4) {
    std::cerr << "missing number of particles" << std::endl;
    return -1;
  }
  int numParticles = 0;
  sscanf(argv[3], "%d", &numParticles);
  
  if(argc < 5) {
    std::cerr << "missing name of particles file" << std::endl;
    return -1;
  }
  char* particlesFileName = argv[4];
  
  if(argc < 6) {
    std::cerr << "missing lower bounds eg XxYxZ" << std::endl;
    return -1;
  }
  FieldData lowerBound[3];
  sscanf(argv[5], "%lfx%lfx%lf", lowerBound, lowerBound + 1, lowerBound + 2);
  
  if(argc < 7) {
    std::cerr << "missing upper bounds eg XxYxZ" << std::endl;
    return -1;
  }
  FieldData upperBound[3];
  sscanf(argv[6], "%lfx%lfx%lf", upperBound, upperBound + 1, upperBound + 2);

  FieldData* rho = new FieldData[numFluidLines];
  FieldData* pressure = new FieldData[numFluidLines];
  FieldData* velocity = new FieldData[numFluidLines * 3];
  FieldData* centerCoordinates = new FieldData[numFluidLines * 3];
  FieldData* temperature = new FieldData[numFluidLines];
  FieldData domainMin[3] = { 99999, 99999, 99999 };
  FieldData domainMax[3] = { 0 };
  
  loadFluidData(fluidFileName, numFluidLines, rho, pressure, velocity, centerCoordinates,  temperature, domainMin, domainMax);
  
  std::cout << "domainMin " << domainMin[0] << " " << domainMin[1] << " " << domainMin[2] << std::endl;
  std::cout << "domainMax " << domainMax[0] << " " << domainMax[1] << " " << domainMax[2] << std::endl;

  long int* particlesID = new long int[numParticles];
  FieldData* particlesPosition = new FieldData[numParticles * 3];
  FieldData* particlesTemperature = new FieldData[numParticles];
  FieldData* particlesDensity = new FieldData[numParticles];

  loadParticlesData(particlesFileName, numParticles, particlesID, particlesPosition, particlesTemperature, particlesDensity);
  const int numParticlesToDraw = 50;
  long int particlesToDraw[numParticlesToDraw] = { 0 };
  for(int i = 0; i < numParticlesToDraw; ++i) particlesToDraw[i] = i;
  
  OSMesaContext mesaCtx;
  GLubyte* rgbaBuffer;
  GLfloat* depthBuffer;

  renderInitialize(lowerBound, upperBound, mesaCtx, rgbaBuffer, depthBuffer);
  renderImage(numFluidX, numFluidY, numFluidZ, rho, pressure, velocity, centerCoordinates, temperature, lowerBound, upperBound, temperatureField, 4.88675,
              numParticles, particlesID, particlesPosition, particlesTemperature, particlesDensity,
              particlesToDraw, numParticlesToDraw);

  saveImageToFile(fluidFileName, rgbaBuffer);
  
  renderTerminate(mesaCtx, rgbaBuffer, depthBuffer);
  
  delete [] rho;
  delete [] pressure;
  delete [] velocity;
  delete [] centerCoordinates;
  delete [] temperature;

}

