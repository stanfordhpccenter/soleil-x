/* Copyright 2016 Stanford University
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "render.h"
#include "renderImage.h"
#include "legion_visualization.h"
#include "legion_c_util.h"
#include "image_reduction_mapper.h"
#include <unistd.h>

using namespace Legion;
using namespace LegionRuntime::Accessor;


#define _T {std::cout<<__FUNCTION__<<":"<<__LINE__<<std::endl;}


#ifdef __cplusplus
extern "C" {
#endif
  
  
  
  static map<int, Visualization::ImageReduction*> gImageCompositors;
  static MapperID gImageReductionMapperID = 0;
  static int gRenderTaskID = 0;
  static int gSaveImageTaskID = 0;

  typedef double FieldData;
  
#define SAVE_RENDER_DATA 0
  
  
  static
  void create_field_pointer_3D(PhysicalRegion region,
                               FieldData* &field,
                               int fieldID,
                               ByteOffset stride[],
                               Runtime* runtime) {
    
    const int dim = 3;
    Domain indexSpaceDomain = runtime->get_index_space_domain(region.get_logical_region().get_index_space());
    LegionRuntime::Arrays::Rect<dim> bounds = indexSpaceDomain.get_rect<dim>();
    RegionAccessor<AccessorType::Generic, FieldData> acc = region.get_field_accessor(fieldID).typeify<FieldData>();
    LegionRuntime::Arrays::Rect<dim> tempBounds;
    field = acc.raw_rect_ptr<dim>(bounds, tempBounds, stride);
    assert(bounds == tempBounds);
  }
  
  static
  void create_field_pointer_1D(PhysicalRegion region,
                               FieldData* &field,
                               int fieldID,
                               ByteOffset stride[],
                               Runtime* runtime) {
    
    const int dim = 1;
    Domain indexSpaceDomain = runtime->get_index_space_domain(region.get_logical_region().get_index_space());
    LegionRuntime::Arrays::Rect<dim> bounds = indexSpaceDomain.get_rect<dim>();
    RegionAccessor<AccessorType::Generic, FieldData> acc = region.get_field_accessor(fieldID).typeify<FieldData>();
    LegionRuntime::Arrays::Rect<dim> tempBounds;
    field = acc.raw_rect_ptr<dim>(bounds, tempBounds, stride);
    assert(bounds == tempBounds);
  }
  
  static
  void create_int_pointer_1D(PhysicalRegion region,
                             long int* &field,
                             int fieldID,
                             ByteOffset stride[],
                             Runtime* runtime) {
    
    const int dim = 1;
    Domain indexSpaceDomain = runtime->get_index_space_domain(region.get_logical_region().get_index_space());
    LegionRuntime::Arrays::Rect<dim> bounds = indexSpaceDomain.get_rect<dim>();
    RegionAccessor<AccessorType::Generic, long int> acc = region.get_field_accessor(fieldID).typeify<long int>();
    LegionRuntime::Arrays::Rect<dim> tempBounds;
    field = acc.raw_rect_ptr<dim>(bounds, tempBounds, stride);
    assert(bounds == tempBounds);
  }
  
#if SAVE_RENDER_DATA

  static void saveFluidRenderData(Context ctx,
                                  HighLevelRuntime *runtime,
                                  const Task* task,
                                  PhysicalRegion& fluid,
                                  std::vector<legion_field_id_t> fluidFields) {
    
    char filename[256] = "render.fluid.";
    sprintf(filename + strlen(filename), "%lld", task->get_unique_id());
    FILE *fluidOut = fopen(filename, "w");
    
    FieldData* rho;
    FieldData* pressure;
    FieldData* velocity;
    FieldData* centerCoordinates;
    FieldData* temperature;
    
    ByteOffset rhoStride[3];
    ByteOffset pressureStride[3];
    ByteOffset velocityStride[3];
    ByteOffset centerCoordinatesStride[3];
    ByteOffset temperatureStride[3];
    
    create_field_pointer_3D(fluid, rho, fluidFields[0], rhoStride, runtime);
    create_field_pointer_3D(fluid, pressure, fluidFields[1], pressureStride, runtime);
    create_field_pointer_3D(fluid, velocity, fluidFields[2], velocityStride, runtime);
    create_field_pointer_3D(fluid, centerCoordinates, fluidFields[3], centerCoordinatesStride, runtime);
    create_field_pointer_3D(fluid, temperature, fluidFields[7], temperatureStride, runtime);
    
    IndexSpace indexSpace = fluid.get_logical_region().get_index_space();
    Domain domain = runtime->get_index_space_domain(ctx, indexSpace);
    Rect<3> rect = domain;
    
    for (PointInRectIterator<3> pir(rect); pir(); pir++) {
      fprintf(fluidOut, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
              *rho, *pressure,
              velocity[0], velocity[1], velocity[2],
              centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
              *temperature);
      rho += rhoStride[0].offset / sizeof(*rho);
      pressure += pressureStride[0].offset / sizeof(*pressure);
      velocity += velocityStride[0].offset / sizeof(*velocity);
      centerCoordinates += centerCoordinatesStride[0].offset / sizeof(*centerCoordinates);
      temperature += temperatureStride[0].offset / sizeof(*temperature);
    }
    
    fclose(fluidOut);
  }
  
  static void saveParticlesRenderData(Context ctx,
                                      HighLevelRuntime *runtime,
                                      const Task* task,
                                      PhysicalRegion& particles,
                                      std::vector<legion_field_id_t> particlesFields) {
    
    char filename[256] = "render.particles.";
    sprintf(filename + strlen(filename), "%lld", task->get_unique_id());
    FILE *particlesOut = fopen(filename, "w");
    
    long int* id;
    FieldData* position;
    FieldData* temperature;
    FieldData* density;
    
    ByteOffset idStride[1];
    ByteOffset positionStride[1];
    ByteOffset temperatureStride[1];
    ByteOffset densityStride[1];
    
    create_int_pointer_1D(particles, id, particlesFields[0], idStride, runtime);
    create_field_pointer_1D(particles, position, particlesFields[2], positionStride, runtime);
    create_field_pointer_1D(particles, temperature, particlesFields[4], temperatureStride, runtime);
    create_field_pointer_1D(particles, density, particlesFields[6], densityStride, runtime);
    
    IndexSpace indexSpace = particles.get_logical_region().get_index_space();
    Domain domain = runtime->get_index_space_domain(ctx, indexSpace);
    Rect<1> rect = domain;
    
    for (PointInRectIterator<1> pir(rect); pir(); pir++) {
      fprintf(particlesOut, "%ld\t%g\t%g\t%g\t%g\t%g\n",
              *id, position[0], position[1], position[2],
              *temperature, *density);
      id += idStride[0].offset / sizeof(*id);
      position += positionStride[0].offset / sizeof(*position);
      temperature += temperatureStride[0].offset / sizeof(*temperature);
      density += densityStride[0].offset / sizeof(*density);
    }
    
    fclose(particlesOut);
  }
#endif
  
  
  static void render_task(const Task *task,
                          const std::vector<PhysicalRegion> &regions,
                          Context ctx, HighLevelRuntime *runtime) {
    char hostname[128];
    gethostname(hostname, sizeof(hostname));
    std::cout << "in render_task " << task->task_id << " " << task->get_unique_id() << " pid " << getpid() << " " << hostname << std::endl;
    
    PhysicalRegion fluid = regions[0];
    PhysicalRegion particles = regions[1];
    PhysicalRegion image = regions[2];
    
    std::vector<legion_field_id_t> fluidFields;
    fluid.get_fields(fluidFields);
    
    std::vector<legion_field_id_t> particlesFields;
    particles.get_fields(particlesFields);
    
    std::vector<legion_field_id_t> imageFields;
    image.get_fields(imageFields);
    
#if SAVE_RENDER_DATA
    saveFluidRenderData(ctx, runtime, task, fluid, fluidFields);
    saveParticlesRenderData(ctx, runtime, task, particles, particlesFields);
#else
    ImageDescriptor* imageDescriptor = (ImageDescriptor*)(task->args);
    FieldData* lowerBound = (FieldData*)((char*)task->args + sizeof(ImageDescriptor));
    FieldData* upperBound = lowerBound + 3;
    int numParticlesToDraw = *((int*)(upperBound + 3));
    long int* particlesToDraw = (long int*)((char*)task->args + sizeof(ImageDescriptor) + 6 * sizeof(FieldData) + sizeof(int));
    
    OSMesaContext mesaCtx;
    GLubyte* rgbaBuffer;
    GLfloat* depthBuffer;

    renderInitialize(lowerBound, upperBound, mesaCtx, rgbaBuffer, depthBuffer);
    
    FieldData* rho;
    FieldData* pressure;
    FieldData* velocity;
    FieldData* centerCoordinates;
    FieldData* temperature;
    
    ByteOffset rhoStride[3];
    ByteOffset pressureStride[3];
    ByteOffset velocityStride[3];
    ByteOffset centerCoordinatesStride[3];
    ByteOffset temperatureStride[3];
    
    create_field_pointer_3D(fluid, rho, fluidFields[0], rhoStride, runtime);
    create_field_pointer_3D(fluid, pressure, fluidFields[1], pressureStride, runtime);
    create_field_pointer_3D(fluid, velocity, fluidFields[2], velocityStride, runtime);
    create_field_pointer_3D(fluid, centerCoordinates, fluidFields[3], centerCoordinatesStride, runtime);
    create_field_pointer_3D(fluid, temperature, fluidFields[7], temperatureStride, runtime);
    
    long int* id;
    FieldData* particlesPosition;
    FieldData* particlesTemperature;
    FieldData* particlesDensity;
    
    ByteOffset idStride[1];
    ByteOffset positionStride[1];
    ByteOffset particlesTemperatureStride[1];
    ByteOffset densityStride[1];
    
    create_int_pointer_1D(particles, id, particlesFields[0], idStride, runtime);
    create_field_pointer_1D(particles, particlesPosition, particlesFields[2], positionStride, runtime);
    create_field_pointer_1D(particles, particlesTemperature, particlesFields[4], particlesTemperatureStride, runtime);
    create_field_pointer_1D(particles, particlesDensity, particlesFields[6], densityStride, runtime);
    
    IndexSpace indexSpace = fluid.get_logical_region().get_index_space();
    Rect<3> bounds = runtime->get_index_space_domain(ctx, indexSpace);
    long int num[3] = {
      bounds.hi[0] - bounds.lo[0], bounds.hi[1] - bounds.lo[1], bounds.hi[2] - bounds.lo[2]
    };
    IndexSpace particlesIndexSpace = particles.get_logical_region().get_index_space();
    Rect<1> particlesBounds = runtime->get_index_space_domain(ctx, particlesIndexSpace);
    long int numParticles = particlesBounds.hi[0] - particlesBounds.lo[0];
    
    renderImage(num[0], num[1], num[2], rho, pressure, velocity, centerCoordinates, temperature, lowerBound, upperBound, temperatureField, 4.88675,
                numParticles, id, particlesPosition, particlesTemperature, particlesDensity,
                particlesToDraw, numParticlesToDraw);
    
    // copy rendered pixels into source image region
    glReadPixels(0, 0, imageDescriptor->width, imageDescriptor->height, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, depthBuffer);
    FieldData* r;
    FieldData* g;
    FieldData* b;
    FieldData* a;
    FieldData* z;
    
    ByteOffset rStride[3];
    ByteOffset gStride[3];
    ByteOffset bStride[3];
    ByteOffset aStride[3];
    ByteOffset zStride[3];
    
    create_field_pointer_3D(image, r, imageFields[0], rStride, runtime);
    create_field_pointer_3D(image, g, imageFields[1], gStride, runtime);
    create_field_pointer_3D(image, b, imageFields[2], bStride, runtime);
    create_field_pointer_3D(image, a, imageFields[3], aStride, runtime);
    create_field_pointer_3D(image, z, imageFields[4], zStride, runtime);

    for(int i = 0; i < imageDescriptor->width * imageDescriptor->height; ++i) {
      *r = rgbaBuffer[i * 4];
      *g = rgbaBuffer[i * 4 + 1];
      *b = rgbaBuffer[i * 4 + 2];
      *a = rgbaBuffer[i * 4 + 3];
      *z = depthBuffer[i];
      r += rStride[0].offset / sizeof(*r);
      g += gStride[0].offset / sizeof(*g);
      b += bStride[0].offset / sizeof(*b);
      a += aStride[0].offset / sizeof(*a);
      z += zStride[0].offset / sizeof(*z);
    }
    
    renderTerminate(mesaCtx, rgbaBuffer, depthBuffer);
#endif
    
  }
  
  static int gFrameNumber = 0;
  
  static void save_image_task(const Task *task,
                          const std::vector<PhysicalRegion> &regions,
                          Context ctx, HighLevelRuntime *runtime) {
    ImageDescriptor* imageDescriptor = (ImageDescriptor*)(task->args);
    PhysicalRegion image = regions[0];
    std::vector<legion_field_id_t> imageFields;
    image.get_fields(imageFields);

    FieldData* r;
    FieldData* g;
    FieldData* b;
    FieldData* a;
    FieldData* z;
    
    ByteOffset rStride[3];
    ByteOffset gStride[3];
    ByteOffset bStride[3];
    ByteOffset aStride[3];
    ByteOffset zStride[3];
    
    create_field_pointer_3D(image, r, imageFields[0], rStride, runtime);
    create_field_pointer_3D(image, g, imageFields[1], gStride, runtime);
    create_field_pointer_3D(image, b, imageFields[2], bStride, runtime);
    create_field_pointer_3D(image, a, imageFields[3], aStride, runtime);
    create_field_pointer_3D(image, z, imageFields[4], zStride, runtime);

    char filename[128];
    sprintf(filename, "image.%05d.tga", gFrameNumber++);
    FILE* f = fopen(filename, "w");
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
    int x, y;
    for (y=imageDescriptor->height-1; y>=0; y--) {
      for (x=0; x<imageDescriptor->width; x++) {
        GLubyte b_ = *b;
        fputc(b_, f); /* write blue */
        GLubyte g_ = *g;
        fputc(g_, f); /* write green */
        GLubyte r_ = *r;
        fputc(r_, f);   /* write red */
        r += rStride[0].offset / sizeof(*r);
        g += gStride[0].offset / sizeof(*g);
        b += bStride[0].offset / sizeof(*b);
      }
    }
    fclose(f);
    std::cout << "wrote image " << filename << std::endl;
  }
  
  
  // Called from mapper before runtime has started
  void cxx_preinitialize(MapperID mapperID)
  {
    Visualization::ImageReduction::preinitializeBeforeRuntimeStarts();
    
    // allocate physical regions contiguously in memory
    LayoutConstraintRegistrar layout_registrar(FieldSpace::NO_SPACE, "SOA layout");
    std::vector<DimensionKind> dim_order(2);
    dim_order[0] = DIM_X;
    dim_order[1] = DIM_F; // fields go last for SOA
    layout_registrar.add_constraint(OrderingConstraint(dim_order, true/*contig*/));
    LayoutConstraintID soa_layout_id = Runtime::preregister_layout(layout_registrar);
    
    // preregister render task
    gImageReductionMapperID = mapperID;
    gRenderTaskID = Legion::HighLevelRuntime::generate_static_task_id();
    TaskVariantRegistrar registrar(gRenderTaskID, "render_task");
    registrar.add_constraint(ProcessorConstraint(Processor::LOC_PROC))
    .add_layout_constraint_set(0/*index*/, soa_layout_id)
    .add_layout_constraint_set(1/*index*/, soa_layout_id)
    .add_layout_constraint_set(2/*index*/, soa_layout_id);
    Runtime::preregister_task_variant<render_task>(registrar, "render_task");
    
    // preregister save image task
    gSaveImageTaskID = Legion::HighLevelRuntime::generate_static_task_id();
    TaskVariantRegistrar registrarSaveImage(gSaveImageTaskID, "save_image_task");
    registrarSaveImage.add_constraint(ProcessorConstraint(Processor::LOC_PROC))
    .add_layout_constraint_set(0/*index*/, soa_layout_id);
    Runtime::preregister_task_variant<save_image_task>(registrarSaveImage, "save_image_task");

  }
  
  
  void cxx_render(legion_runtime_t runtime_,
                  legion_context_t ctx_,
                  legion_mapper_id_t sampleId,
                  legion_physical_region_t* fluid_,
                  FieldID fluidFields[],
                  int numFluidFields,
                  legion_physical_region_t* particles_,
                  FieldID particlesFields[],
                  int numParticlesFields,
                  legion_index_space_t tiles_,
                  legion_logical_partition_t fluidPartition_,
                  legion_logical_partition_t particlesPartition_,
                  int numParticlesToDraw,
                  legion_physical_region_t* particlesToDraw_,
                  FieldData lowerBound[3],
                  FieldData upperBound[3]
                  )
  {
    std::cout << __FUNCTION__ << " pid " << getpid() << std::endl;
    
    PhysicalRegion* fluid = CObjectWrapper::unwrap(fluid_[0]);
    std::vector<legion_field_id_t> fluid_fields;
    fluid->get_fields(fluid_fields);
    const FieldAccessor<READ_ONLY, float, 3> fluid_temperature_acc(*fluid, fluid_fields[0]);
    const FieldAccessor<READ_ONLY, float, 3> fluid_pressure_acc(*fluid, fluid_fields[1]);
    const FieldAccessor<READ_ONLY, float, 3> fluid_velocity_acc(*fluid, fluid_fields[2]);
    
    Runtime *runtime = CObjectWrapper::unwrap(runtime_);
    Context ctx = CObjectWrapper::unwrap(ctx_)->context();
    LogicalPartition logicalPartition = CObjectWrapper::unwrap(fluidPartition_);
    Visualization::ImageDescriptor imageDescriptor = { imageWidth, imageHeight, 1, 1 };
    
    if(gImageCompositors.find(sampleId) == gImageCompositors.end()) {
      gImageCompositors[sampleId] = new Visualization::ImageReduction(logicalPartition, imageDescriptor, ctx, runtime, gImageReductionMapperID);
      ImageReductionMapper::registerRenderTaskName("render_task");
    }
    
    Visualization::ImageReduction* compositor = gImageCompositors[sampleId];
    imageDescriptor = compositor->imageDescriptor();
    
    PhysicalRegion* particles = CObjectWrapper::unwrap(particles_[0]);
    std::vector<legion_field_id_t> particles_fields;
    particles->get_fields(particles_fields);
    const FieldAccessor<READ_ONLY, float, 1> particles_position_acc(*particles, particles_fields[0]);
    const FieldAccessor<READ_ONLY, float, 1> particles_position_history_acc(*particles, particles_fields[1]);
    
    PhysicalRegion* particlesToDraw = CObjectWrapper::unwrap(particlesToDraw_[0]);
    std::vector<legion_field_id_t> particlesToDrawFields;
    particlesToDraw->get_fields(particlesToDrawFields);
    long int* particlesToDrawInt;
    ByteOffset particlesToDrawStride[1];
    create_int_pointer_1D(*particlesToDraw, particlesToDrawInt, particlesToDrawFields[0], particlesToDrawStride, runtime);
    
    ArgumentMap argMap;
    size_t argSize = sizeof(imageDescriptor) + 6 * sizeof(FieldData) + sizeof(int) + numParticlesToDraw * sizeof(long int);
    char args[argSize] = { 0 };
    memcpy(args, &imageDescriptor, sizeof(imageDescriptor));
    memcpy(args + sizeof(imageDescriptor), lowerBound, 3 * sizeof(FieldData));
    memcpy(args + sizeof(imageDescriptor) + 3 * sizeof(FieldData), upperBound, 3 * sizeof(FieldData));
    memcpy(args + sizeof(imageDescriptor) + 6 * sizeof(FieldData), &numParticlesToDraw, sizeof(int));
    memcpy(args + sizeof(imageDescriptor) + 6 * sizeof(FieldData) + sizeof(int), particlesToDrawInt, numParticlesToDraw * sizeof(long int));

    IndexTaskLauncher renderLauncher(gRenderTaskID, compositor->everywhereDomain(), TaskArgument(args, sizeof(args)), argMap, Predicate::TRUE_PRED, false, gImageReductionMapperID);
    
    LogicalPartition fluidPartition = CObjectWrapper::unwrap(fluidPartition_);
    ImageReductionProjectionFunctor functor0(compositor->everywhereDomain(), fluidPartition);
    runtime->register_projection_functor(1, &functor0);
    RegionRequirement req0(fluidPartition, 1, READ_ONLY, SIMULTANEOUS, fluid->get_logical_region(), gImageReductionMapperID);
    for(int i = 0; i < numFluidFields; ++i) req0.add_field(fluidFields[i]);
    renderLauncher.add_region_requirement(req0);
    
    LogicalPartition particlesPartition = CObjectWrapper::unwrap(particlesPartition_);
    ImageReductionProjectionFunctor functor1(compositor->everywhereDomain(), particlesPartition);
    runtime->register_projection_functor(2, &functor1);
    RegionRequirement req1(particlesPartition, 2, READ_ONLY, SIMULTANEOUS, particles->get_logical_region(), gImageReductionMapperID);
    for(int i = 0; i < numParticlesFields; ++i) req1.add_field(particlesFields[i]);
    renderLauncher.add_region_requirement(req1);
    
    RegionRequirement req2(compositor->depthPartition(), 0, READ_WRITE, EXCLUSIVE, compositor->sourceImage(), gImageReductionMapperID);
    ImageReductionProjectionFunctor functor2(compositor->everywhereDomain(), compositor->depthPartition());
    runtime->register_projection_functor(3, &functor2);
    req2.add_field(Visualization::ImageReduction::FID_FIELD_R);
    req2.add_field(Visualization::ImageReduction::FID_FIELD_G);
    req2.add_field(Visualization::ImageReduction::FID_FIELD_B);
    req2.add_field(Visualization::ImageReduction::FID_FIELD_A);
    req2.add_field(Visualization::ImageReduction::FID_FIELD_Z);
    renderLauncher.add_region_requirement(req2);
    
    ImageReductionMapper::clearPlacement(logicalPartition);
    
    FutureMap futures = runtime->execute_index_space(ctx, renderLauncher);
    futures.wait_all_results();
  }
  
  
  
  
  void cxx_reduce(legion_runtime_t runtime_,
                  legion_context_t ctx_,
                  legion_mapper_id_t sampleId
                  )
  {
    Runtime *runtime = CObjectWrapper::unwrap(runtime_);
    Context ctx = CObjectWrapper::unwrap(ctx_)->context();

    Visualization::ImageReduction* compositor = gImageCompositors[sampleId];
    compositor->set_depth_func(GL_LESS);
    FutureMap futures = compositor->reduce_associative_commutative();
    futures.wait_all_results();
    
    // save the image
    ImageDescriptor imageDescriptor = compositor->imageDescriptor();
    TaskLauncher saveImageLauncher(gSaveImageTaskID, TaskArgument(&imageDescriptor, sizeof(ImageDescriptor)), Predicate::TRUE_PRED, gImageReductionMapperID);
    DomainPoint slice0 = Point<3>::ZEROES();
    LogicalRegion imageSlice0 = runtime->get_logical_subregion_by_color(compositor->depthPartition(), slice0);
    RegionRequirement req(imageSlice0, READ_ONLY, EXCLUSIVE, compositor->sourceImage());
    saveImageLauncher.add_region_requirement(req);
    saveImageLauncher.add_field(0/*idx*/, Visualization::ImageReduction::FID_FIELD_R);
    saveImageLauncher.add_field(0/*idx*/, Visualization::ImageReduction::FID_FIELD_G);
    saveImageLauncher.add_field(0/*idx*/, Visualization::ImageReduction::FID_FIELD_B);
    saveImageLauncher.add_field(0/*idx*/, Visualization::ImageReduction::FID_FIELD_A);
    saveImageLauncher.add_field(0/*idx*/, Visualization::ImageReduction::FID_FIELD_Z);
    Future result = runtime->execute_task(ctx, saveImageLauncher);
    result.get_result<int>();
  }
  
#ifdef __cplusplus
}
#endif


