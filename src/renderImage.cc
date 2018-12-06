
#include <assert.h>
#include <stdio.h>
#include <iostream>
#include <math.h>

#include "renderImage.h"

void vDrawScene(int numFluidX,
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
                FieldData targetValue);

void renderParticles(int numParticles, const long int* particlesID, const FieldData3* particlesPosition, const FieldData* particlesTemperature, const FieldData* particlesDensity, int numParticlesToDraw, long int* particlesToDraw);

void initializeMarchingCubes(GLfloat lightPosition[4]);


void setCameraPosition(FieldData domainMin[3], FieldData domainMax[3]) {
  // TODO add camera motion
  glViewport( 0, 0, WIDTH, HEIGHT );
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  float depthMax = sqrt((domainMax[2] - domainMin[2]) * (domainMax[2] - domainMin[2])
                        + (domainMax[1] - domainMin[1]) * (domainMax[1] - domainMin[1])
                        + (domainMax[0] - domainMin[0]) * (domainMax[0] - domainMin[0]));
  glOrtho(-depthMax * 0.67, depthMax * 0.67, -depthMax * 0.3, depthMax * 1.3, 0, depthMax);
  
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  
#define FIX_CAMERA 0

#if !FIX_CAMERA
  FieldData scale[3];
  for(unsigned i = 0; i < 3; ++i) scale[i] = domainMax[i] - domainMin[i];
  FieldData scaleOffset = 0.25;
#endif
  
  GLfloat from[] =
#if FIX_CAMERA
  { 3, .1, -.5 };
#else
  { (GLfloat)(domainMin[0] - scale[0] * scaleOffset),
    (GLfloat)(domainMin[1] - (scale[1] * scaleOffset)),
    (GLfloat)(domainMin[2] - scale[2] * scaleOffset) };
#endif
  GLfloat at[] =
#if FIX_CAMERA
  { 4, .1, .1 };
#else
  { (GLfloat)((domainMin[0] + domainMax[0]) * 0.5),
    (GLfloat)((domainMin[1] + domainMax[1]) * 0.1),
    (GLfloat)((domainMin[2] + domainMax[2]) * 0.5) };
#endif
  GLfloat up[] = { 0, 1, 0 };
  std::cout << "camera from " << from[0] << "," << from[1] << "," << from[2] << std::endl;
  std::cout << "camera at   " << at[0] << "," << at[1] << "," << at[2] << std::endl;
  std::cout << "camera up   " << up[0] << "," << up[1] << "," << up[2] << std::endl;
  gluLookAt(from[0], from[1], from[2], at[0], at[1], at[2], up[0], up[1], up[2]);
}

void createGraphicsContext(OSMesaContext &mesaCtx,
                           GLubyte* &rgbaBuffer,
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
  rgbaBuffer = new GLubyte[WIDTH * HEIGHT * fieldsPerPixel];
  if (!rgbaBuffer) {
    printf("Alloc image buffer failed!\n");
    return;
  }

  /* Bind the buffer to the context and make it current */
  if (!OSMesaMakeCurrent(mesaCtx, rgbaBuffer, GL_UNSIGNED_BYTE, WIDTH, HEIGHT)) {
    printf("OSMesaMakeCurrent failed!\n");
    return;
  }

  {
    int z, s, a;
    glGetIntegerv(GL_DEPTH_BITS, &z);
    glGetIntegerv(GL_STENCIL_BITS, &s);
    glGetIntegerv(GL_ACCUM_RED_BITS, &a);
    printf("Depth=%d Stencil=%d Accum=%d\n", z, s, a);
  }

  /* Allocate the depth buffer. */
  depthBuffer = new GLfloat[WIDTH * HEIGHT];
  if (!depthBuffer) {
    printf("Alloc depth buffer failed!\n");
    return;
  }
}

static void destroyGraphicsContext(OSMesaContext mesaCtx) {
  /* destroy the context */
  OSMesaDestroyContext(mesaCtx);
}


void renderInitialize(FieldData domainMin[3], FieldData domainMax[3],
                      OSMesaContext& mesaCtx, GLubyte*& rgbaBuffer, GLfloat*& depthBuffer) {
  createGraphicsContext(mesaCtx, rgbaBuffer, depthBuffer);
  GLfloat lightPosition[4];
  lightPosition[0] = 0.5 * (domainMax[0] - domainMin[0]);
  lightPosition[1] = domainMax[1] * 1.5;
  lightPosition[2] = 0.5 * (domainMax[2] - domainMin[2]);
  lightPosition[3] = 1.0;
  initializeMarchingCubes(lightPosition);
  std::cout << "light position " << lightPosition[0] << "," << lightPosition[1] << "," << lightPosition[2] << "," << lightPosition[3] << std::endl;
}

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
                 int numParticles,
                 const long int* particlesID,
                 const FieldData3* particlesPosition,
                 const FieldData* particlesTemperature,
                 const FieldData* particlesDensity,
                 long int* particlesToDraw,
                 int numParticlesToDraw) {
  
  std::cout << "domain min " << domainMin[0] << "," << domainMin[1] << "," << domainMin[2] << std::endl;
  std::cout << "domain max " << domainMax[0] << "," << domainMax[1] << "," << domainMax[2] << std::endl;
  
  setCameraPosition(domainMin, domainMax);
  vDrawScene(numFluidX, numFluidY, numFluidZ, rho, pressure, velocity, centerCoordinates, temperature, domainMin, domainMax, visualizationField, targetValue);
  renderParticles(numParticles, particlesID, particlesPosition, particlesTemperature, particlesDensity, numParticlesToDraw, particlesToDraw);
}


void renderTerminate(OSMesaContext mesaCtx, GLubyte*& rgbaBuffer, GLfloat*& depthBuffer) {
  destroyGraphicsContext(mesaCtx);
  delete [] rgbaBuffer;
  delete [] depthBuffer;
}


void
write_targa(const char *filename, const GLubyte *buffer, int width, int height)
{
  FILE *f = fopen( filename, "w" );
  if (f) {
    int i, x, y;
    const GLubyte *ptr = buffer;
    fputc (0x00, f);  /* ID Length, 0 => No ID   */
    fputc (0x00, f);  /* Color Map Type, 0 => No color map included   */
    fputc (0x02, f);  /* Image Type, 2 => Uncompressed, True-color Image */
    fputc (0x00, f);  /* Next five bytes are about the color map entries */
    fputc (0x00, f);  /* 2 bytes Index, 2 bytes length, 1 byte size */
    fputc (0x00, f);
    fputc (0x00, f);
    fputc (0x00, f);
    fputc (0x00, f);  /* X-origin of Image */
    fputc (0x00, f);
    fputc (0x00, f);  /* Y-origin of Image */
    fputc (0x00, f);
    fputc (WIDTH & 0xff, f);      /* Image Width */
    fputc ((WIDTH>>8) & 0xff, f);
    fputc (HEIGHT & 0xff, f);     /* Image Height   */
    fputc ((HEIGHT>>8) & 0xff, f);
    fputc (0x18, f);     /* Pixel Depth, 0x18 => 24 Bits  */
    fputc (0x20, f);     /* Image Descriptor  */
    fclose(f);
    f = fopen( filename, "ab" );  /* reopen in binary append mode */
    for (y=HEIGHT-1; y>=0; y--) {
      for (x=0; x<WIDTH; x++) {
        i = (y*WIDTH + x) * 4;
        fputc(ptr[i+2], f); /* write blue */
        fputc(ptr[i+1], f); /* write green */
        fputc(ptr[i], f);   /* write red */
      }
    }
    fclose(f);
    std::cout << "wrote targa file " << filename << std::endl;
  }
}


void
write_ppm(const char *filename, const GLubyte *rgba, int width, int height)
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
      GLubyte* rgbaPtr = (GLubyte*)rgba + (y * width * 4);
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
