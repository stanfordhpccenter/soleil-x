//
// render.cc
//

#include "render.h"
#include "legion_visualization.h"
#include "GL/osmesa.h"
#include "GL/glu.h"
#include <math.h>

  using namespace Legion;

  // global data

  static Visualization::ImageReduction* gImageCompositor = nullptr;
  static int gRenderTaskID = 0;
  static int gSaveImageTaskID = 0;
  static int gFrameNumber = 0;
  static int gImageWidth = 1280;
  static int gImageHeight = 720;
  static int gNumParticlesToDraw = 0;
  static long int* gParticlesToDraw = nullptr;
  static legion_field_id_t* gParticlesFields;
  static int gNumParticlesFields;
  static MapperID gImageReductionMapperID = 0;
  static int gInitialRenderTaskID = 0;
  static int gRenderProjectionFunctorID = 99;//TODO system assigns this

  typedef double FieldData;

      class RenderProjectionFunctor : public ProjectionFunctor {
      public:
        RenderProjectionFunctor(int tiles[3]) {
          memcpy(mTiles, tiles, sizeof(mTiles));
          mNx = tiles[0];
          mNy = tiles[1];
          mNz = tiles[2];
        }

        virtual LogicalRegion project(const Mappable *mappable, unsigned index,
                                      LogicalPartition upperBound,
                                      const DomainPoint &point) {
          long int z = point[0] * mNy * mNz + point[1] * mNz + point[2];
          Legion::Point<3> remappedPoint;
          remappedPoint[0] = 0;
          remappedPoint[1] = 0;
          remappedPoint[2] = z;
          DomainPoint remappedDomainPoint(remappedPoint);

          LogicalRegion result = Legion::Runtime::get_runtime()->get_logical_subregion_by_color(upperBound, remappedDomainPoint);
          return result;
        }

        virtual LogicalRegion project(LogicalPartition upper_bound, const DomainPoint &point, const Domain &launch_domain) {
          assert(false);
          LogicalRegion result;
          return result;
        }

        virtual LogicalRegion project(const Mappable *mappable, unsigned index, LogicalRegion upper_bound, const DomainPoint &point) {
          assert(false);
          LogicalRegion result;
          return result;
        }

        virtual LogicalRegion project(LogicalRegion upper_bound, const DomainPoint &point, const Domain &launch_domain) {
          assert(false);
          LogicalRegion result;
          return result;
        }


        virtual bool is_exclusive(void) const{ return true; }
        virtual unsigned get_depth(void) const{ return 0; }
        virtual bool is_functional(void) const { return true; }

      private:
        int mTiles[3];
        long int mNx;
        long int mNy;
        long int mNz;
	};

  static RenderProjectionFunctor* mRenderProjectionFunctor = nullptr;


  static void createGraphicsContext(OSMesaContext &mesaCtx,
                             GLubyte* &rgbaBuffer,
                             GLfloat* &depthBuffer,
                             int width,
                             int height)
  {
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
        rgbaBuffer = new GLubyte[width * height * fieldsPerPixel];
        if (!rgbaBuffer) {
          printf("Alloc image buffer failed!\n");
          return;
        }

        /* Bind the buffer to the context and make it current */
        if (!OSMesaMakeCurrent(mesaCtx, rgbaBuffer, GL_UNSIGNED_BYTE, width, height)) {
          printf("OSMesaMakeCurrent failed!\n");
          return;
        }

        /* Allocate the depth buffer. */
        depthBuffer = new GLfloat[width * height];
        if (!depthBuffer) {
          printf("Alloc depth buffer failed!\n");
          return;
        }
    }

    static void initializeRender(Camera* camera, int width, int height) {
      GLfloat afPropertiesAmbient [] = {1.00, 1.00, 1.00, 0.5};
      GLfloat afPropertiesDiffuse [] = {1.00, 1.00, 1.00, 0.5};
      GLfloat afPropertiesSpecular[] = {1.00, 1.00, 1.00, 0.5};
      glClearColor( 0.3, 0.3, 0.3, 0 );
      glEnable(GL_DEPTH_TEST);
      glDepthFunc(GL_LESS);
      //glEnable(GL_BLEND);
      //glBlendFunc(GL_SRC_ALPHA, GL_DST_ALPHA);
      glEnable(GL_LIGHTING);
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      glLightfv(GL_LIGHT0, GL_AMBIENT,  afPropertiesAmbient);
      glLightfv(GL_LIGHT0, GL_DIFFUSE,  afPropertiesDiffuse);
      glLightfv(GL_LIGHT0, GL_SPECULAR, afPropertiesSpecular);
      GLfloat lightPosition[] = { 1, 4, 1, 1 };
      glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
      glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0);
      glEnable(GL_LIGHT0);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      glViewport(0, 0, width, height);
      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      glLoadIdentity();
      GLfloat fovy = 20;
      GLfloat aspect = (GLfloat)width / (GLfloat)height;
      GLfloat near = 0.0;
      GLfloat far = 50.0;
      gluPerspective(fovy, aspect, near, far);

      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();
      glLoadIdentity();
      gluLookAt(camera->from[0], camera->from[1], camera->from[2],
                camera->at[0], camera->at[1], camera->at[2],
                camera->up[0], camera->up[1], camera->up[2]);
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

    static void temperatureToColor(float temperature, float color[4])
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


    static void scaledTemperatureToColor(float temperature, float color[4], double colorScale[2])
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


  static void initial_render_task(const Task *task,
                          const std::vector<PhysicalRegion> &regions,
                          Context ctx, HighLevelRuntime *runtime) {
    long int* args = (long int*)task->args;
    int tiles[3];
    tiles[0] = args[0];
    tiles[1] = args[1];
    tiles[2] = args[2];
    mRenderProjectionFunctor = new RenderProjectionFunctor(tiles);
    runtime->register_projection_functor(gRenderProjectionFunctorID, mRenderProjectionFunctor);
    if(gParticlesToDraw == nullptr) {
      gNumParticlesToDraw = args[3];
#if 1
{
char b[256];
gethostname(b, sizeof(b));
sprintf(b + strlen(b), " gNumParticlesToDraw = %d\n", gNumParticlesToDraw);
}
#endif
      gParticlesToDraw = new long int[gNumParticlesToDraw];
      memcpy(gParticlesToDraw, args + 4, gNumParticlesToDraw * sizeof(long int));
    }
  }

  static int comparPoints(const void* p0, const void* p1) {
    GLfloat* f0 = (GLfloat*)p0;
    GLfloat* f1 = (GLfloat*)p1;
    if(f0[3] < f1[3]) return -1;
    if(f0[3] > f1[3]) return +1;
    return 0;
  }



  static void transformPoints(AccessorRO<FieldData3, 1> particlesPosition,
                              AccessorRO<bool, 1> particlesValid,
                              AccessorRO<long int, 1> particlesID,
                              long int numParticlesToDraw,
                              long int *particlesToDraw,
                              int lo,
                              int hi,
                              GLfloat points[],
                              long int& numValid)
{
  numValid = 0;
  int drawn = 0;
  GLfloat modelView[16];
  glGetFloatv(GL_MODELVIEW_MATRIX, modelView);

  for(long long int i = lo; i < hi; ++i) {
    FieldData3 position = particlesPosition[i];
    bool valid = particlesValid[i];
    if(valid) numValid++;
    if(valid && drawThis(particlesID[i], numParticlesToDraw, particlesToDraw)) {
      drawn++;
      GLfloat x = position.x[0];
      GLfloat y = position.x[1];
      GLfloat z = position.x[2];
      GLfloat zTransformed = x * modelView[2] + y * modelView[6] + z * modelView[10] + modelView[14];
      points[numValid * 5] = x;
      points[numValid * 5 + 1] = y;
      points[numValid * 5 + 2] = z;
      points[numValid * 5 + 3] = zTransformed;
      points[numValid * 5 + 4] = i;
    }
  }
  qsort(points, numValid, 5 * sizeof(GLfloat), comparPoints);
#if 1
{
char b[128];
gethostname(b, sizeof(b));
sprintf(b + strlen(b), " numValid %ld drawn %d\n", numValid, drawn);
std::cout<<b;
}
#endif
}


  static void render_task(const Task *task,
                          const std::vector<PhysicalRegion> &regions,
                          Context ctx, HighLevelRuntime *runtime) {
__TRACE
    PhysicalRegion particles = regions[0];
    PhysicalRegion image = regions[1];
    char* argsPtr = (char*)task->args;

    ImageDescriptor* imageDescriptor = (ImageDescriptor*)(argsPtr);
    argsPtr += sizeof(ImageDescriptor);
    Camera* camera = (Camera*)argsPtr;
    argsPtr += sizeof(Camera);
    double* colorScale = (double*)argsPtr;
    argsPtr += 2 * sizeof(double);
    legion_field_id_t* particlesFields = (legion_field_id_t*)argsPtr;

    OSMesaContext mesaCtx;
    unsigned char* rgbaBuffer = nullptr;
    float* depthBuffer = nullptr;
    createGraphicsContext(mesaCtx, rgbaBuffer, depthBuffer,
      imageDescriptor->width, imageDescriptor->height);
    initializeRender(camera, imageDescriptor->width, imageDescriptor->height);

    AccessorRO<long int, 1> particlesID(particles, particlesFields[0]);
    AccessorRO<FieldData3, 1> particlesPosition(particles, particlesFields[1]);
    AccessorRO<FieldData, 1> particlesTemperature(particles, particlesFields[2]);
    AccessorRO<FieldData, 1> particlesDensity(particles, particlesFields[3]);
    AccessorRO<bool, 1> particlesValid(particles, particlesFields[4]);

    Domain particlesDomain = runtime->get_index_space_domain(
      particles.get_logical_region().get_index_space());

    long long int lo = particlesDomain.bounds<1, long long int>().lo[0];
    long long int hi = particlesDomain.bounds<1, long long int>().hi[0];
    GLfloat* points = new GLfloat[5 * (hi - lo)];
    long int numValid;
    transformPoints(particlesPosition, particlesValid, particlesID,
      gNumParticlesToDraw, gParticlesToDraw, lo, hi, points, numValid);

    GLfloat min[3] = { 999, 999, 999 };
    GLfloat max[3] = { -999, -999, -999 };
    glBegin(GL_POINTS);
    for(long int i = 0; i < numValid; ++i) {
        long long int index = (long long int)points[5 * i + 4];
        if(particlesDensity[index] > 0) {
          GLfloat t = particlesTemperature[index];
          GLfloat color[4];
          scaledTemperatureToColor(t, color, colorScale);

          glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, color);
          glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color);
          glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color);

          glVertex3f(points[5 * i],  points[5 * i + 1], points[5 * i + 2]);

          for(unsigned j = 0; j < 3; j++) {
            if(min[j] > points[5 * i + j]) min[j] = points[5 * i + j];
            if(max[j] < points[5 * i + j]) max[j] = points[5 * i + j];
          }
        }
    }
    glEnd();

#if 1
{
__TRACE
char b[256];
gethostname(b, sizeof(b));
sprintf(b + strlen(b), " drawn min %g %g %g max %g %g %g\n",
min[0], min[1], min[2], max[0], max[1], max[2]);
std::cout<<b;
}
#endif

    delete [] points;

    // now copy the image data into the image logical region

    glReadPixels(0, 0, imageDescriptor->width, imageDescriptor->height,
      GL_DEPTH_COMPONENT, GL_FLOAT, depthBuffer);

    std::vector<legion_field_id_t> imageFields;
    image.get_fields(imageFields);
    AccessorWO<ImageReduction::PixelField, 3> r(image, imageFields[0]);
    AccessorWO<ImageReduction::PixelField, 3> g(image, imageFields[1]);
    AccessorWO<ImageReduction::PixelField, 3> b(image, imageFields[2]);
    AccessorWO<ImageReduction::PixelField, 3> a(image, imageFields[3]);
    AccessorWO<ImageReduction::PixelField, 3> z(image, imageFields[4]);
    AccessorWO<ImageReduction::PixelField, 3> u(image, imageFields[5]);

    IndexSpace saveIndexSpace = image.get_logical_region().get_index_space();
    Legion::Rect<3> saveRect = runtime->get_index_space_domain(ctx, saveIndexSpace);

    int index = 0;
    for(PointInRectIterator<3> pir(saveRect); pir(); pir++) {
      r[*pir] = rgbaBuffer[index * 4] / 255.0;
      g[*pir] = rgbaBuffer[index * 4 + 1] / 255.0;
      b[*pir] = rgbaBuffer[index * 4 + 2] / 255.0;
      a[*pir] = rgbaBuffer[index * 4 + 3] / 255.0;
      z[*pir] = depthBuffer[index];
      u[*pir] = 0; // user defined channel, unused
      index++;
    }

    OSMesaDestroyContext(mesaCtx);
    delete [] rgbaBuffer;
    delete [] depthBuffer;
__TRACE
  }



  static int save_image_task(const Task *task,
                             const std::vector<PhysicalRegion> &regions,
                             Context ctx, HighLevelRuntime *runtime) {
    ImageDescriptor* imageDescriptor = (ImageDescriptor*)(task->args);
    unsigned char* outDir = (unsigned char*)(task->args) + sizeof(ImageDescriptor);
    PhysicalRegion image = regions[0];
    IndexSpace indexSpace = image.get_logical_region().get_index_space();
    Legion::Rect<3> bounds = runtime->get_index_space_domain(ctx, indexSpace);
    std::vector<legion_field_id_t> imageFields;
    image.get_fields(imageFields);

    AccessorRO<ImageReduction::PixelField, 3> r(image, imageFields[0]);
    AccessorRO<ImageReduction::PixelField, 3> g(image, imageFields[1]);
    AccessorRO<ImageReduction::PixelField, 3> b(image, imageFields[2]);
    AccessorRO<ImageReduction::PixelField, 3> a(image, imageFields[3]);
    AccessorRO<ImageReduction::PixelField, 3> z(image, imageFields[4]);

    char filename[1024];
    sprintf(filename, "%s/image_%d_%d_%d.%05d.tga", outDir, (
      int)bounds.lo.x, (int)bounds.lo.y, (int)bounds.lo.z, gFrameNumber++);
    FILE* f = fopen(filename, "w");
    if(f == nullptr) {
      std::cerr << "could not create file " << filename << std::endl;
      return -1;
    }
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
    fputc (imageDescriptor->width & 0xff, f);      /* Image Width */
    fputc ((imageDescriptor->width>>8) & 0xff, f);
    fputc (imageDescriptor->height & 0xff, f);     /* Image Height   */
    fputc ((imageDescriptor->height>>8) & 0xff, f);
    fputc (0x18, f);     /* Pixel Depth, 0x18 => 24 Bits  */
    fputc (0x20, f);     /* Image Descriptor  */
    fclose(f);

    f = fopen(filename, "ab");  /* reopen in binary append mode */
    IndexSpace saveIndexSpace = image.get_logical_region().get_index_space();
    Legion::Rect<3> saveRect = runtime->get_index_space_domain(ctx, saveIndexSpace);
    PointInRectIterator<3> pir(saveRect);
    ImageReduction::PixelField* BB = (ImageReduction::PixelField*)b.ptr(*pir);
    ImageReduction::PixelField* GG = (ImageReduction::PixelField*)g.ptr(*pir);
    ImageReduction::PixelField* RR = (ImageReduction::PixelField*)r.ptr(*pir);

    for(int y = imageDescriptor->height - 1; y >= 0; y--) {
      for(int x = 0; x < imageDescriptor->width; ++x) {
        int index = x + y * imageDescriptor->width;
        GLubyte b_ = BB[index] * 255;
        fputc(b_, f); /* write blue */
        GLubyte g_ = GG[index] * 255;
        fputc(g_, f); /* write green */
        GLubyte r_ = RR[index] * 255;
        fputc(r_, f);   /* write red */
      }
    }
    fclose(f);
    std::cout << "wrote image " << filename << std::endl;
    return 0;
  }








  // Called from mapper before runtime has started
  void cxx_preinitialize(MapperID mapperID) {
    gImageReductionMapperID = mapperID;

    Visualization::ImageReduction::preinitializeBeforeRuntimeStarts();
    // allocate physical regions contiguously in memory
    LayoutConstraintRegistrar layout_registrar(FieldSpace::NO_SPACE, "SOA layout");
    std::vector<DimensionKind> dim_order(4);
    dim_order[0] = DIM_X;
    dim_order[1] = DIM_Y;
    dim_order[2] = DIM_Z;
    dim_order[3] = DIM_F; // fields go last for SOA
    layout_registrar.add_constraint(OrderingConstraint(dim_order, true/*contig*/));
    LayoutConstraintID soa_layout_id = Runtime::preregister_layout(layout_registrar);

    // preregister initial_render_task
    gInitialRenderTaskID = Legion::HighLevelRuntime::generate_static_task_id();
    TaskVariantRegistrar registrarCPFT(gInitialRenderTaskID, "initial_render_task");
    registrarCPFT.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
    Runtime::preregister_task_variant<initial_render_task>(registrarCPFT, "initial_render_task");

    // preregister render task
    gRenderTaskID = Legion::HighLevelRuntime::generate_static_task_id();
    TaskVariantRegistrar registrarRender(gRenderTaskID, "render_task");
    registrarRender.add_constraint(ProcessorConstraint(Processor::LOC_PROC))
    .add_layout_constraint_set(0/*index*/, soa_layout_id)
    .add_layout_constraint_set(1/*index*/, soa_layout_id)
    .add_layout_constraint_set(2/*index*/, soa_layout_id);
    Runtime::preregister_task_variant<render_task>(registrarRender, "render_task");

    // preregister save image task
    gSaveImageTaskID = Legion::HighLevelRuntime::generate_static_task_id();
    TaskVariantRegistrar registrarSaveImage(gSaveImageTaskID, "save_image_task");
    registrarSaveImage.add_constraint(ProcessorConstraint(Processor::LOC_PROC))
    .add_layout_constraint_set(0/*index*/, soa_layout_id);
    Runtime::preregister_task_variant<int, save_image_task>(registrarSaveImage,
      "save_image_task");
  }





static int compar(const void* p1, const void* p2) {
  int i1 = *(int*)p1;
  int i2 = *(int*)p2;
  if(i1 < i2) return -1;
  if(i1 > i2) return +1;
  return 0;
}




  // this entry point is called once from the main task
  void cxx_initialize(
                     legion_runtime_t runtime_,
                     legion_context_t ctx_,
                     legion_logical_region_t region_,
                     legion_logical_partition_t partition_,
                     legion_field_id_t pFields[],
                     int numPFields,
                     int numParticlesToDraw_,
                     int sampleId,
                     int tag,
                     int tiles[3],
                     long int numParticles
                     )
  {
    Runtime *runtime = CObjectWrapper::unwrap(runtime_);
    Context ctx = CObjectWrapper::unwrap(ctx_)->context();
    LogicalRegion region = CObjectWrapper::unwrap(region_);
    LogicalPartition partition = CObjectWrapper::unwrap(partition_);

    Visualization::ImageDescriptor imageDescriptor = { gImageWidth, gImageHeight, 1 };

    gImageCompositor = new Visualization::ImageReduction(region,
      partition, pFields, numPFields, imageDescriptor, ctx, runtime, gImageReductionMapperID);

    // pack the tiles into the first three positions of the arguments
    gNumParticlesToDraw = numParticlesToDraw_;
    const int numTiles = 3;
    const int numNum = 1;
    long int* particleBuffer = new long int[gNumParticlesToDraw + numTiles + numNum];
    particleBuffer[0] = tiles[0];
    particleBuffer[1] = tiles[1];
    particleBuffer[2] = tiles[2];
    particleBuffer[3] = gNumParticlesToDraw;
    srandom(0);
    for(int i = 0; i < gNumParticlesToDraw; ++i) {
      long double r = (long double)random();
      long int id = numParticles * r / RAND_MAX;
      particleBuffer[i + numTiles + numNum] = id;
    }
    qsort(particleBuffer + numTiles + numNum, gNumParticlesToDraw, sizeof(particleBuffer[0]), compar);

    gNumParticlesFields = numPFields;
    gParticlesFields = new legion_field_id_t[numPFields];
    for(int i = 0; i < numPFields; ++i) {
      gParticlesFields[i] = pFields[i];
    }
    gParticlesToDraw = particleBuffer + numTiles + numNum;

    size_t argLen = (gNumParticlesToDraw + numTiles + numNum) * sizeof(long int);
    gImageCompositor->initializeRenderNodes(runtime, ctx, gInitialRenderTaskID, (char*)particleBuffer, argLen);
  }







  // this entry point is called once from the main task
  void cxx_render(legion_runtime_t runtime_,
                  legion_context_t ctx_,
                  double cameraFromAtUp[9],
                  double colorScale[2]
                  ) {
__TRACE
    // Unwrap objects

    Runtime *runtime = CObjectWrapper::unwrap(runtime_);
    Context ctx = CObjectWrapper::unwrap(ctx_)->context();
    Camera camera;
    camera.from[0] = cameraFromAtUp[0];
    camera.from[1] = cameraFromAtUp[1];
    camera.from[2] = cameraFromAtUp[2];
    camera.at[0] = cameraFromAtUp[3];
    camera.at[1] = cameraFromAtUp[4];
    camera.at[2] = cameraFromAtUp[5];
    camera.up[0] = cameraFromAtUp[6];
    camera.up[1] = cameraFromAtUp[7];
    camera.up[2] = cameraFromAtUp[8];
    Visualization::ImageReduction* compositor = gImageCompositor;

    // Setup the render task launch with region requirements
    ArgumentMap argMap;
    ImageDescriptor imageDescriptor = compositor->imageDescriptor();
    size_t argSize = sizeof(ImageDescriptor) + sizeof(camera) + 2 * sizeof(double) + gNumParticlesFields * sizeof(int);
    char args[argSize];
    char *argsPtr = args;
    memcpy(argsPtr, (char*)&imageDescriptor, sizeof(ImageDescriptor));
    argsPtr += sizeof(ImageDescriptor);
    memcpy(argsPtr, (char*)&camera, sizeof(camera));
    argsPtr += sizeof(camera);
    memcpy(argsPtr, (char*)colorScale, 2 * sizeof(double));
    argsPtr += 2 * sizeof(double);
    memcpy(argsPtr, gParticlesFields, gNumParticlesFields * sizeof(legion_field_id_t));
    argsPtr += gNumParticlesFields * sizeof(legion_field_id_t);

    IndexTaskLauncher renderLauncher(gRenderTaskID, compositor->renderImageDomain(), TaskArgument(args, argSize),
                                     argMap, Predicate::TRUE_PRED, false, gImageReductionMapperID);

    RegionRequirement req0(imageDescriptor.simulationLogicalPartition, 0, READ_ONLY, EXCLUSIVE,
      imageDescriptor.simulationLogicalRegion, gImageReductionMapperID);
    for(int i = 0; i < imageDescriptor.numPFields; ++i) {
      req0.add_field(imageDescriptor.pFields[i]);
    }
    renderLauncher.add_region_requirement(req0);

    RegionRequirement req1(compositor->compositeImagePartition(), gRenderProjectionFunctorID, WRITE_DISCARD, EXCLUSIVE,
      compositor->sourceImage(), gImageReductionMapperID);
    req1.add_field(Visualization::ImageReduction::FID_FIELD_R);
    req1.add_field(Visualization::ImageReduction::FID_FIELD_G);
    req1.add_field(Visualization::ImageReduction::FID_FIELD_B);
    req1.add_field(Visualization::ImageReduction::FID_FIELD_A);
    req1.add_field(Visualization::ImageReduction::FID_FIELD_Z);
    req1.add_field(Visualization::ImageReduction::FID_FIELD_USERDATA);
    renderLauncher.add_region_requirement(req1);

    FutureMap futures = runtime->execute_index_space(ctx, renderLauncher);
    futures.wait_all_results();
__TRACE
  }






  void cxx_reduce(legion_context_t ctx_, double cameraFromAtUp[9]
) {
    Context ctx = CObjectWrapper::unwrap(ctx_)->context();
    Camera camera;
    camera.from[0] = cameraFromAtUp[0];
    camera.from[1] = cameraFromAtUp[1];
    camera.from[2] = cameraFromAtUp[2];
    camera.at[0] = cameraFromAtUp[3];
    camera.at[1] = cameraFromAtUp[4];
    camera.at[2] = cameraFromAtUp[5];
    camera.up[0] = cameraFromAtUp[6];
    camera.up[1] = cameraFromAtUp[7];
    camera.up[2] = cameraFromAtUp[8];
    Visualization::ImageReduction* compositor = gImageCompositor;
    compositor->set_blend_func(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    compositor->set_blend_equation(GL_FUNC_ADD);
    float cameraDirection[] = {
      (float)(camera.at[0] - camera.from[0]),
      (float)(camera.at[1] - camera.from[1]),
      (float)(camera.at[2] - camera.from[2])
    };
    FutureMap futures = compositor->reduceImages(ctx, cameraDirection);
    futures.wait_all_results();
  }




  void cxx_saveImage(legion_runtime_t runtime_,
                     legion_context_t ctx_,
                     const char* outDir
                     ) {
    Runtime *runtime = CObjectWrapper::unwrap(runtime_);
    Context ctx = CObjectWrapper::unwrap(ctx_)->context();

    // save the image
    Visualization::ImageReduction* compositor = gImageCompositor;
    ImageDescriptor imageDescriptor = compositor->imageDescriptor();
    size_t argLen = sizeof(ImageDescriptor) + strlen(outDir) + 1;
    char args[argLen];
    memcpy(args, &imageDescriptor, sizeof(imageDescriptor));
    strcpy(args + sizeof(imageDescriptor), outDir);
    TaskLauncher saveImageLauncher(gSaveImageTaskID, TaskArgument(args, argLen),
      Predicate::TRUE_PRED, gImageReductionMapperID);
    DomainPoint slice0 = Legion::Point<3>::ZEROES();
    LogicalRegion imageSlice0 = runtime->get_logical_subregion_by_color(compositor->compositeImagePartition(), slice0);
    RegionRequirement req(imageSlice0, READ_ONLY, EXCLUSIVE, compositor->sourceImage(),
      gImageReductionMapperID);
    saveImageLauncher.add_region_requirement(req);
    saveImageLauncher.add_field(0/*idx*/, Visualization::ImageReduction::FID_FIELD_R);
    saveImageLauncher.add_field(0/*idx*/, Visualization::ImageReduction::FID_FIELD_G);
    saveImageLauncher.add_field(0/*idx*/, Visualization::ImageReduction::FID_FIELD_B);
    saveImageLauncher.add_field(0/*idx*/, Visualization::ImageReduction::FID_FIELD_A);
    saveImageLauncher.add_field(0/*idx*/, Visualization::ImageReduction::FID_FIELD_Z);
    Future result = runtime->execute_task(ctx, saveImageLauncher);
    result.get_result<int>();
  }


  // this is provided for debugging the renderer

  void cxx_saveIndividualImages(legion_runtime_t runtime_,
                                legion_context_t ctx_,
                                const char* outDir
                                ) {
    Runtime *runtime = CObjectWrapper::unwrap(runtime_);
    Context ctx = CObjectWrapper::unwrap(ctx_)->context();

    // save the image
    Visualization::ImageReduction* compositor = gImageCompositor;
    ImageDescriptor imageDescriptor = compositor->imageDescriptor();
    ArgumentMap argMap;
    size_t argLen = sizeof(ImageDescriptor) + strlen(outDir) + 1;
    char args[argLen];
    memcpy(args, &imageDescriptor, sizeof(imageDescriptor));
    strcpy(args + sizeof(imageDescriptor), outDir);
    IndexTaskLauncher saveImageLauncher(gSaveImageTaskID, compositor->compositeImageDomain(), TaskArgument(args, argLen),
                                     argMap, Predicate::TRUE_PRED, false);
    RegionRequirement req(compositor->compositeImagePartition(), 0, READ_ONLY, EXCLUSIVE, compositor->sourceImage());
    saveImageLauncher.add_region_requirement(req);
    saveImageLauncher.add_field(0/*idx*/, Visualization::ImageReduction::FID_FIELD_R);
    saveImageLauncher.add_field(0/*idx*/, Visualization::ImageReduction::FID_FIELD_G);
    saveImageLauncher.add_field(0/*idx*/, Visualization::ImageReduction::FID_FIELD_B);
    saveImageLauncher.add_field(0/*idx*/, Visualization::ImageReduction::FID_FIELD_A);
    saveImageLauncher.add_field(0/*idx*/, Visualization::ImageReduction::FID_FIELD_Z);
    FutureMap futures = runtime->execute_index_space(ctx, saveImageLauncher);
    futures.wait_all_results();
  }


  void cxx_terminate() {
    delete gImageCompositor;
  }
