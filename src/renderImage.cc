
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
                VisualizationField isosurfaceField,
                FieldData targetValue,
                FieldData colorScale[2],
                GLfloat cameraLookAt[3]);

void renderParticles(int numParticles, const long int* particlesID, const FieldData3* particlesPosition, const FieldData* particlesTemperature, const FieldData* particlesDensity, int numParticlesToDraw, long int* particlesToDraw, float systemScale,
                     FieldData colorScale[2]);

void initializeMarchingCubes(GLfloat lightPosition[4]);


void setupRender(FieldData domainMin[3], FieldData domainMax[3], float* depthMax, GLfloat cameraLookAt[3]) {
  // TODO add camera motion
  glViewport( 0, 0, WIDTH, HEIGHT );
  std::cout << "viewport width, height " << WIDTH << "," << HEIGHT << std::endl;

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  *depthMax = 1.0 / 3.0 * sqrt((domainMax[2] - domainMin[2]) * (domainMax[2] - domainMin[2])
                        + (domainMax[1] - domainMin[1]) * (domainMax[1] - domainMin[1])
                        + (domainMax[0] - domainMin[0]) * (domainMax[0] - domainMin[0]));
  
#define ORTHO 0
#if ORTHO
  GLfloat left = - *depthMax;
  GLfloat right = *depthMax;
  GLfloat bottom = -*depthMax * 0.5;
  GLfloat top = *depthMax * 2.0;
  GLfloat near = 0.0;
  GLfloat far = *depthMax * 3.0;
  glOrtho(left, right, bottom, top, near, far);
#else
  GLfloat fovy = 20;
  GLfloat aspect = (GLfloat)WIDTH / (GLfloat)HEIGHT;
  GLfloat near = 0.0;
  GLfloat far = *depthMax * 2.0;
  gluPerspective(fovy, aspect, near, far);
#endif
  
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  
  // view the particles coming toward the camera
  GLfloat from[] =
  { (GLfloat)(domainMax[0] * 2.0),
    (GLfloat)(domainMax[1] * 4.0),
    (GLfloat)(domainMin[2] + domainMax[2] * 0.5) };

  GLfloat at[] =
  { (GLfloat)((domainMin[0] + domainMax[0]) * 0.5),
    (GLfloat)((domainMin[1] + domainMax[1]) * 0.5),
    (GLfloat)((domainMin[2] + domainMax[2]) * 0.5) };

  GLfloat up[] = { 0, 1, 0 };

#if 1
#if ORTHO
  std::cout << "glOrtho " << left << "," << right << "  " << bottom << "," << top << "  " << near << "," << far << std::endl;
#else
  std::cout << "gluPersepctive " << fovy << " " << aspect << " " << near << " " << far << std::endl;
#endif
  std::cout << "camera from " << from[0] << "," << from[1] << "," << from[2] << std::endl;
  std::cout << "camera at   " << at[0] << "," << at[1] << "," << at[2] << std::endl;
  std::cout << "camera up   " << up[0] << "," << up[1] << "," << up[2] << std::endl;
  std::cout << "domainMin   " << domainMin[0] << "," << domainMin[1] << "," << domainMin[2] << std::endl;
  std::cout << "domainMax   " << domainMax[0] << "," << domainMax[1] << "," << domainMax[2] << std::endl;
  std::cout << "systemScale " << *depthMax << std::endl;
#endif
  gluLookAt(from[0], from[1], from[2], at[0], at[1], at[2], up[0], up[1], up[2]);
  for(unsigned i = 0; i < 3; ++i) cameraLookAt[i] = at[i] - from[i];
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
}

#define DRAW_COORDINATES 0

#if DRAW_COORDINATES

static void drawParticle(GLUquadricObj* qobj, const FieldData3* position, FieldData density, FieldData particleTemperature) {

  GLfloat color[4];
  color[0] = color[1] = color[2] = color[3] = 1;//debug
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, color);
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color);

  glPushMatrix();
  GLfloat x = position->x[0];
  GLfloat y = position->x[1];
  GLfloat z = position->x[2];
  glTranslatef(x, y, z);
  const GLfloat densityScale = 1.0e-5;
  gluSphere(qobj, density * densityScale, 7, 7);
  glPopMatrix();
}


static void drawCoordinates(int nx, int ny, int nz, const FieldData3* coordinates) {
  GLUquadricObj *qobj = gluNewQuadric();
   for(int i = 0; i < nx * ny * nz; ++i) {
     drawParticle(qobj, coordinates + i, 2000, 0);
   }
  gluDeleteQuadric(qobj);
}

#endif

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
                 VisualizationField isosurfaceField,
                 FieldData targetValue,
                 FieldData colorScale[2],
                 int numParticles,
                 const long int* particlesID,
                 const FieldData3* particlesPosition,
                 const FieldData* particlesTemperature,
                 const FieldData* particlesDensity,
                 long int* particlesToDraw,
                 int numParticlesToDraw) {
  
//  std::cout << "domain min " << domainMin[0] << "," << domainMin[1] << "," << domainMin[2] << std::endl;
//  std::cout << "domain max " << domainMax[0] << "," << domainMax[1] << "," << domainMax[2] << std::endl;
#if 1
  switch(isosurfaceField) {
    case noneField:
      std::cout << "no isosurface" << std::endl;
      break;
    case temperatureField:
      std::cout << "isosurface of temperature = " << targetValue << std::endl;
      break;
    case rhoField:
      std::cout << "isosurface of rho = " << targetValue << std::endl;
      break;
    case pressureField:
      std::cout << "isosurface of pressure = " << targetValue << std::endl;
      break;
    case velocityField:
      std::cerr << "error, velocity field not selectable for isosurface" << std::endl;
      break;
  }
#endif
  
  float systemScale;
  GLfloat cameraLookAt[3];
  setupRender(domainMin, domainMax, &systemScale, cameraLookAt);
std::cout << "isosurfaceField " << isosurfaceField << std::endl;
  vDrawScene(numFluidX, numFluidY, numFluidZ, rho, pressure, velocity, centerCoordinates, temperature, domainMin, domainMax, isosurfaceField, targetValue, colorScale, cameraLookAt);
  renderParticles(numParticles, particlesID, particlesPosition, particlesTemperature, particlesDensity, numParticlesToDraw, particlesToDraw, systemScale, colorScale);

#if 1
  double rhoSum = 0;
  double pressureSum = 0;
  double temperatureSum = 0;
  unsigned numCells = numFluidX * numFluidY * numFluidZ;
  for(unsigned i = 0; i < numCells; ++i) {
    rhoSum += rho[i];
    pressureSum += pressure[i];
    temperatureSum += temperature[i];
  }
  std::cout << "mean rho " << (rhoSum / numCells) << std::endl;
  std::cout << "mean pressure " << (pressureSum / numCells) << std::endl;
  std::cout << "mean temperature " << (temperatureSum / numCells) << std::endl;
#endif

#if DRAW_COORDINATES
  drawCoordinates(numFluidX, numFluidY, numFluidZ, centerCoordinates);
#endif
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

#define DEBUG_COLOR 0

void temperatureToColorOLD(float temperature, float color[4])
{
  
  // https://stackoverflow.com/questions/7229895/display-temperature-as-a-color-with-c
  
  float x = (float)(temperature / 1000.0);
  float x2 = x * x;
  float x3 = x2 * x;
  float x4 = x3 * x;
  float x5 = x4 * x;
  
  float R, G, B = 0.0f;

#if DEBUG_COLOR
  std::cout << "scaled temperature (0..10000) " << temperature << std::endl;
#endif
  
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
  
#if 1
  color[0] = R;
  color[1] = G;
  color[2] = B;
  color[3] = 1.0f;
#else
  color[0] = B;
  color[1] = G;
  color[2] = R;
  color[3] = 1.0f;
#endif
}


// reference http://www.tannerhelland.com/4435/convert-temperature-rgb-algorithm-code/

void temperatureToColor(float temperature, float color[4])
{
  float R, G, B;
  float t = temperature / 100.0;
  
#if DEBUG_COLOR
  std::cout << "scaled temperature " << temperature << std::endl;
#endif

  if(t <= 66) {
    B = 255.0;
#if DEBUG_COLOR
    std::cout << "t<=66 " << t << " blue full on" << std::endl;
#endif
  } else {
    if(t > 113) {
      B = 0.0;
#if DEBUG_COLOR
      std::cout << "t>113 " << t << " blue off" << std::endl;
#endif
    } else {
      B = t - 60.0;
      B = 329.698727446 * pow(B, -0.1332047592);
      if(B < 0) B = 0.0;
      if(B > 255.0) B = 255.0;
#if DEBUG_COLOR
      std::cout << "t>66 " << t << " blue " << B << std::endl;
#endif
    }
  }

  if(t <= 66) {
    G = t;
    G = 99.4708025861 * log(G) - 161.1195681661;
    if(G < 0) G = 0.0;
    if(G > 255.0) G = 255.0;
#if DEBUG_COLOR
    std::cout << "t<=66 " << t << " green " << G << std::endl;
#endif
  } else {
    if(t > 113) {
      G = 0.0;
#if DEBUG_COLOR
      std::cout << "t>113 " << t << " green off" << std::endl;
#endif
    } else {
      G = t - 60.0;
      G = 288.1221695283 * pow(G, -0.0755148492);
      if(G < 0) G = 0.0;
      if(G > 255.0) G = 255.0;
#if DEBUG_COLOR
      std::cout << "t>66 " << t << " green " << G << std::endl;
#endif
    }
  }

  if(t >= 66) {
    R = 255.0;
#if DEBUG_COLOR
    std::cout << "t>=66 " << t << " red full on" << std::endl;
#endif
  } else {
    if(t <= 19) {
      R = 0.0;
#if DEBUG_COLOR
      std::cout << "t<=19 " << t << " red zero" << std::endl;
#endif
    } else {
      R = t - 10.0;
      R = 138.5177312231 * log(R) - 305.0447927307;
      if(R < 0) R = 0.0;
      if(R > 255.0) R = 255.0;
#if DEBUG_COLOR
      std::cout << "t>19 " << t << " red " << R << std::endl;
#endif
    }
  }

  color[0] = R / 255.0;
  color[1] = G / 255.0;
  color[2] = B / 255.0;
  color[3] = 1.0;
}


void scaledTemperatureToColor(float temperature, float color[4], FieldData colorScale[2])
{
  const float min = colorScale[0];
  const float max = colorScale[1];
  const float Kmin = 1500.0f;
  const float Kmax = 15000.0f;
  // stretch it on a scale of Kmin...Kmax
  float scaledTemperature = Kmin + (temperature - min) * ((Kmax - Kmin) / (max - min));
#if DEBUG_COLOR
  std::cout << "raw temperature " << temperature << " scaled " << scaledTemperature << " colorScale " << colorScale[0] << " " << colorScale[1] << std::endl;
#endif
  return temperatureToColor(scaledTemperature, color);
}



