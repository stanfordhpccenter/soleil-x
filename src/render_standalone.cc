
//
// This is a standalone renderer, used to develop the visualization algorithm offline.
// Keep this file for future development.

#include <iostream>
#include <stdio.h>

#include "render_standalone.h"
#include "renderImage.h"


void initializeMarchingCubes();



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



void saveImageToFile(char* fluidFileName,
                     int numLines,
                     FieldData* rho,
                     FieldData* pressure,
                     FieldData* velocity,
                     FieldData* centerCoordinates,
                     FieldData* temperature) {
  // base output file name on fluidFileName
}


int main(int argc, char **argv) {
  
  if(argc < 2) {
    std::cerr << "missing number of lines in fluid file" << std::endl;
    return -1;
  }
  int numFluidLines = 0;
  sscanf(argv[1], "%d", &numFluidLines);
  
  if(argc < 3) {
    std::cerr << "missing number of lines in particles file" << std::endl;
    return -1;
  }
  int numParticlesLines = 0;
  sscanf(argv[2], "%d", &numParticlesLines);
  
  if(argc < 4) {
    std::cerr << "missing name of fluid file" << std::endl;
  }
  char* fluidFileName = argv[3];
  
  FieldData* rho = new FieldData[numFluidLines];
  FieldData* pressure = new FieldData[numFluidLines];
  FieldData* velocity = new FieldData[numFluidLines * 3];
  FieldData* centerCoordinates = new FieldData[numFluidLines * 3];
  FieldData* temperature = new FieldData[numFluidLines];
  FieldData domainMin[3] = { 99999, 99999, 99999 };
  FieldData domainMax[3] = { 0 };
  
  loadFluidData(fluidFileName, numFluidLines, rho, pressure, velocity, centerCoordinates,  temperature, domainMin, domainMax);
  
  initializeMarchingCubes();
  
  renderImage(numFluidLines, rho, pressure, velocity, centerCoordinates, temperature, domainMin, domainMax);
  
  saveImageToFile(fluidFileName, numFluidLines, rho, pressure, velocity, centerCoordinates, temperature);
}


