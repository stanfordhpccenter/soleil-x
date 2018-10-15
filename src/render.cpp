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

using namespace Legion;

#ifdef __cplusplus
extern "C" {
#endif
  
  
  
  static map<int, Visualization::ImageReduction*> gImageCompositors;
  static MapperID gImageReductionMapperID = 0;
  static int gRenderTaskID = 0;
  
  static void render_task(const Task *task,
                          const std::vector<PhysicalRegion> &regions,
                          Context ctx, HighLevelRuntime *runtime) {
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
                  FieldID* fluidFields,
                  legion_physical_region_t* particles_,
                  FieldID* particlesFields,
                  legion_index_space_t tiles_,
                  legion_logical_partition_t fluidPartition_,
                  legion_logical_partition_t particlesPartition_
                  )
  {
    
    PhysicalRegion* fluid = CObjectWrapper::unwrap(fluid_[0]);
    std::vector<legion_field_id_t> fluid_fields;
    fluid->get_fields(fluid_fields);
    const FieldAccessor<READ_ONLY, float, 3> fluid_temperature_acc(*fluid, fluid_fields[0]);
    const FieldAccessor<READ_ONLY, float, 3> fluid_pressure_acc(*fluid, fluid_fields[1]);
    const FieldAccessor<READ_ONLY, float, 3> fluid_velocity_acc(*fluid, fluid_fields[2]);
    
    Runtime *runtime = CObjectWrapper::unwrap(runtime_);
    Context ctx = CObjectWrapper::unwrap(ctx_)->context();
    LogicalPartition logicalPartition = CObjectWrapper::unwrap(fluidPartition_);
    Visualization::ImageDescriptor imageDescriptor = { imageWidth, imageHeight, 1, 1, logicalPartition };
    
    if(gImageCompositors.find(sampleId) == gImageCompositors.end()) {
      gImageCompositors[sampleId] = new Visualization::ImageReduction(logicalPartition, imageDescriptor, ctx, runtime, gImageReductionMapperID);
    }
    
    Visualization::ImageReduction* compositor = gImageCompositors[sampleId];
    
    PhysicalRegion* particles = CObjectWrapper::unwrap(particles_[0]);
    std::vector<legion_field_id_t> particles_fields;
    particles->get_fields(particles_fields);
    const FieldAccessor<READ_ONLY, float, 1> particles_position_acc(*particles, particles_fields[0]);
    const FieldAccessor<READ_ONLY, float, 1> particles_position_history_acc(*particles, particles_fields[1]);
    
//    IndexSpace tiles = CObjectWrapper::unwrap(tiles_);
//    LogicalPartition fluidPartition = CObjectWrapper::unwrap(fluidPartition_);
//    LogicalPartition particlesPartition = CObjectWrapper::unwrap(particlesPartition_);

    ArgumentMap argMap;
    IndexTaskLauncher renderLauncher(gRenderTaskID, compositor->everywhereDomain(), TaskArgument(&imageDescriptor, sizeof(imageDescriptor)), argMap, Predicate::TRUE_PRED, false, gImageReductionMapperID);
    
    RegionRequirement req0(fluid->get_logical_region(), 0, READ_ONLY, SIMULTANEOUS, fluid->get_logical_region(), gImageReductionMapperID);
    unsigned numFluidFields = sizeof(fluidFields) / sizeof(fluidFields[0]);
    for(unsigned i = 0; i < numFluidFields; ++i) req0.add_field(fluidFields[i]);
    renderLauncher.add_region_requirement(req0);
    
    RegionRequirement req1(particles->get_logical_region(), 0, READ_ONLY, SIMULTANEOUS, particles->get_logical_region(), gImageReductionMapperID);
    unsigned numParticlesFields = sizeof(particlesFields) / sizeof(particlesFields[0]);
    for(unsigned i = 0; i < numParticlesFields; ++i) req1.add_field(particlesFields[i]);
    renderLauncher.add_region_requirement(req1);
    
    RegionRequirement req2(compositor->sourceImage(), 0, READ_ONLY, EXCLUSIVE, compositor->sourceImage(), gImageReductionMapperID);
    req2.add_field(Visualization::ImageReduction::FID_FIELD_R);
    req2.add_field(Visualization::ImageReduction::FID_FIELD_G);
    req2.add_field(Visualization::ImageReduction::FID_FIELD_B);
    req2.add_field(Visualization::ImageReduction::FID_FIELD_A);
    req2.add_field(Visualization::ImageReduction::FID_FIELD_Z);
    renderLauncher.add_region_requirement(req2);

  }
  
  
  
  
  void cxx_reduce(legion_runtime_t runtime_,
                  legion_context_t ctx_,
                  legion_mapper_id_t sampleId
                  )
  {
    std::cout << "in cxx_reduce" << std::endl;
        
#if 0
    Visualization::ImageReduction* compositor = gImageCompositors[sampleId];
    Runtime *runtime = CObjectWrapper::unwrap(runtime_);
    Context ctx = CObjectWrapper::unwrap(ctx_)->context();
#endif
  }
  
#ifdef __cplusplus
}
#endif

