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
#include "legion_visualization.h"
#include "legion_c_util.h"
#include "image_reduction_mapper.h"
#include <unistd.h>

using namespace Legion;
using namespace LegionRuntime::Accessor;


#ifdef __cplusplus
extern "C" {
#endif
  
  
  
  static map<int, Visualization::ImageReduction*> gImageCompositors;
  static MapperID gImageReductionMapperID = 0;
  static int gRenderTaskID = 0;
  
  typedef double FieldData;
  
#define SAVE_RENDER_DATA 0
  
#if SAVE_RENDER_DATA
  
  static
  void create_field_pointer(PhysicalRegion region,
                            FieldData* &field,
                            int fieldID,
                            ByteOffset stride[3],
                            Runtime* runtime) {
    
    Domain indexSpaceDomain = runtime->get_index_space_domain(region.get_logical_region().get_index_space());
    LegionRuntime::Arrays::Rect<3> bounds = indexSpaceDomain.get_rect<3>();
    RegionAccessor<AccessorType::Generic, FieldData> acc = region.get_field_accessor(fieldID).typeify<FieldData>();
    LegionRuntime::Arrays::Rect<3> tempBounds;
    field = acc.raw_rect_ptr<3>(bounds, tempBounds, stride);
    assert(bounds == tempBounds);
  }
  
  static
  void create_int_pointer(PhysicalRegion region,
                          long int* &field,
                          int fieldID,
                          ByteOffset stride[3],
                          Runtime* runtime) {
    
    Domain indexSpaceDomain = runtime->get_index_space_domain(region.get_logical_region().get_index_space());
    LegionRuntime::Arrays::Rect<3> bounds = indexSpaceDomain.get_rect<3>();
    RegionAccessor<AccessorType::Generic, long int> acc = region.get_field_accessor(fieldID).typeify<long int>();
    LegionRuntime::Arrays::Rect<3> tempBounds;
    field = acc.raw_rect_ptr<3>(bounds, tempBounds, stride);
    assert(bounds == tempBounds);
  }
  
  static void saveFluidRenderData(Context ctx,
                                  HighLevelRuntime *runtime,
                                  PhysicalRegion& fluid,
                                  std::vector<legion_field_id_t> fluidFields) {
    
    char filename[256] = "render.fluid.";
    gethostname(filename + strlen(filename), sizeof(filename) - strlen(filename));
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
    
    create_field_pointer(fluid, rho, fluidFields[0], rhoStride, runtime);
    create_field_pointer(fluid, pressure, fluidFields[1], pressureStride, runtime);
    create_field_pointer(fluid, velocity, fluidFields[2], velocityStride, runtime);
    create_field_pointer(fluid, centerCoordinates, fluidFields[3], centerCoordinatesStride, runtime);
    create_field_pointer(fluid, temperature, fluidFields[7], temperatureStride, runtime);
    
    IndexSpace indexSpace = fluid.get_logical_region().get_index_space();
    Domain domain = runtime->get_index_space_domain(ctx, indexSpace);
    Rect<3> rect = domain;
    std::cout << "renderTask subdomain " << rect << std::endl;
    std::cout << "renderTask indexSpace " << indexSpace << std::endl;
    
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
                                      PhysicalRegion& particles,
                                      std::vector<legion_field_id_t> particlesFields) {
    
    char filename[256] = "render.particles.";
    gethostname(filename + strlen(filename), sizeof(filename) - strlen(filename));
    FILE *particlesOut = fopen(filename, "w");
    
    long int* id;
    FieldData* position;
    FieldData* temperature;
    FieldData* density;
    
    ByteOffset idStride[1];
    ByteOffset positionStride[1];
    ByteOffset temperatureStride[1];
    ByteOffset densityStride[1];
    
    create_int_pointer(particles, id, particlesFields[0], idStride, runtime);
    create_field_pointer(particles, position, particlesFields[2], positionStride, runtime);
    create_field_pointer(particles, temperature, particlesFields[4], temperatureStride, runtime);
    create_field_pointer(particles, density, particlesFields[6], densityStride, runtime);
    
    IndexSpace indexSpace = particles.get_logical_region().get_index_space();
    Domain domain = runtime->get_index_space_domain(ctx, indexSpace);
    Rect<1> rect = domain;
    std::cout << "renderTask subdomain " << rect << std::endl;
    std::cout << "renderTask indexSpace " << indexSpace << std::endl;
    
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
    
#if 0
    PhysicalRegion fluid = regions[0];
    PhysicalRegion particles = regions[1];
    PhysicalRegion image = regions[2];
    
    std::vector<legion_field_id_t> fluidFields;
    fluid.get_fields(fluidFields);
    
    std::vector<legion_field_id_t> particlesFields;
    particles.get_fields(particlesFields);
    
    std::vector<legion_field_id_t> imageFields;
    image.get_fields(imageFields);
#endif
    
#if SAVE_RENDER_DATA
    saveFluidRenderData(ctx, runtime, fluid, fluidFields);
    saveParticlesRenderData(ctx, runtime, particles, particlesFields);
#else
    
#endif
    
  }
  
  
  // Called from mapper before runtime has started
  void cxx_preinitialize(MapperID mapperID)
  {
    Visualization::ImageReduction::preinitializeBeforeRuntimeStarts();
    gImageReductionMapperID = mapperID;
    gRenderTaskID = Legion::HighLevelRuntime::generate_static_task_id();
    TaskVariantRegistrar registrar(gRenderTaskID, "render_task");
    registrar.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
    Runtime::preregister_task_variant<render_task>(registrar, "render_task");
    
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
                  legion_logical_partition_t particlesPartition_
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
    
    ArgumentMap argMap;
std::cout << __FUNCTION__ << " everywhereDomain " << compositor->everywhereDomain() << std::endl;
    IndexTaskLauncher renderLauncher(gRenderTaskID, compositor->everywhereDomain(), TaskArgument(&imageDescriptor, sizeof(imageDescriptor)), argMap, Predicate::TRUE_PRED, false, gImageReductionMapperID);
    
    LogicalPartition fluidPartition = CObjectWrapper::unwrap(fluidPartition_);
    std::cout << "fluid partition " << fluidPartition << std::endl;
    ImageReductionProjectionFunctor functor0(compositor->everywhereDomain(), fluidPartition);
    runtime->register_projection_functor(1, &functor0);
    RegionRequirement req0(fluidPartition, 1, READ_ONLY, SIMULTANEOUS, fluid->get_logical_region(), gImageReductionMapperID);
    for(int i = 0; i < numFluidFields; ++i) req0.add_field(fluidFields[i]);
    renderLauncher.add_region_requirement(req0);
    
#if 0
    LogicalPartition particlesPartition = CObjectWrapper::unwrap(particlesPartition_);
    RegionRequirement req1(particlesPartition, 0, READ_ONLY, SIMULTANEOUS, particles->get_logical_region(), gImageReductionMapperID);
    for(int i = 0; i < numParticlesFields; ++i) req1.add_field(particlesFields[i]);
    renderLauncher.add_region_requirement(req1);
    
    RegionRequirement req2(fluidPartition, 0, READ_WRITE, EXCLUSIVE, compositor->sourceImage(), gImageReductionMapperID);
    req2.add_field(Visualization::ImageReduction::FID_FIELD_R);
    req2.add_field(Visualization::ImageReduction::FID_FIELD_G);
    req2.add_field(Visualization::ImageReduction::FID_FIELD_B);
    req2.add_field(Visualization::ImageReduction::FID_FIELD_A);
    req2.add_field(Visualization::ImageReduction::FID_FIELD_Z);
    renderLauncher.add_region_requirement(req2);
#endif
    
    
    ImageReductionMapper::clearPlacement(logicalPartition);
    
    FutureMap futures = runtime->execute_index_space(ctx, renderLauncher);
    futures.wait_all_results();
  }
  
  
  
  
  void cxx_reduce(legion_runtime_t runtime_,
                  legion_context_t ctx_,
                  legion_mapper_id_t sampleId
                  )
  {
    Visualization::ImageReduction* compositor = gImageCompositors[sampleId];
    compositor->set_depth_func(GL_LESS);
    FutureMap futures = compositor->reduce_associative_commutative();
    futures.wait_all_results();
  }
  
#ifdef __cplusplus
}
#endif


