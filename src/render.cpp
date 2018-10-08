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



static bool imageCompositorInitialized = false;
static int numDomainNodes = -1;

    // Called from mapper to set number of nodes
void cxx_preinitialize1(int numDomainNodes_)
{
std::cout << __FUNCTION__ << " " << numDomainNodes_ << std::endl;
  numDomainNodes = numDomainNodes_;
std::cout << "return from " << __FUNCTION__ << std::endl;
}

    // Called from mapper before runtime has started
void cxx_preinitialize2()
{
  Visualization::ImageReduction::initialize();
}

    // Called from application after runtime has started
static void initialize(legion_runtime_t runtime_,
                       legion_context_t ctx_) 
{
std::cout << __FUNCTION__ << " " << numDomainNodes << std::endl;
  if(!imageCompositorInitialized) {
    Runtime *runtime = CObjectWrapper::unwrap(runtime_);
    Context ctx = CObjectWrapper::unwrap(ctx_)->context();
    Visualization::ImageSize imageSize = { imageWidth, imageHeight, numDomainNodes, 1 };
    Visualization::ImageReduction imageReduction(imageSize, ctx, runtime);
    imageCompositorInitialized = true;
  } 
}


void cxx_render(legion_runtime_t runtime_,
                legion_context_t ctx_,
                legion_physical_region_t* fluid_,
                legion_physical_region_t* particles_,
                legion_index_space_t tiles_,
                legion_logical_partition_t fluidPartition_,
                legion_logical_partition_t particlesPartition_
                )
{
  initialize(runtime_, ctx_);

  PhysicalRegion* fluid = CObjectWrapper::unwrap(fluid_[0]);
  std::vector<legion_field_id_t> fluid_fields;
  fluid->get_fields(fluid_fields);
  const FieldAccessor<READ_ONLY, float, 3> fluid_temperature_acc(*fluid, fluid_fields[0]);
  const FieldAccessor<READ_ONLY, float, 3> fluid_pressure_acc(*fluid, fluid_fields[1]);
  const FieldAccessor<READ_ONLY, float, 3> fluid_velocity_acc(*fluid, fluid_fields[2]);
  
  PhysicalRegion* particles = CObjectWrapper::unwrap(particles_[0]);
  std::vector<legion_field_id_t> particles_fields;
  particles->get_fields(particles_fields);
  const FieldAccessor<READ_ONLY, float, 1> particles_position_acc(*particles, particles_fields[0]);
  const FieldAccessor<READ_ONLY, float, 1> particles_position_history_acc(*particles, particles_fields[1]);
  
#if 0 // notused yet
  IndexSpace tiles = CObjectWrapper::unwrap(tiles_);
  LogicalPartition fluidPartition = CObjectWrapper::unwrap(fluidPartition_);
  LogicalPartition particlesPartition = CObjectWrapper::unwrap(particlesPartition_);
#endif
  
#ifdef RENDER_IMAGE
  Runtime *runtime = CObjectWrapper::unwrap(runtime_);
  Context ctx = CObjectWrapper::unwrap(ctx_)->context();
  // create OpenGL context on first occurrence
  // glClear
  // execute marching cubes
  // visualize particles with tails
  // write result to image pr
#endif
  

}




void cxx_reduce(legion_runtime_t runtime_,
                legion_context_t ctx_
                )
{
  std::cout << "in cxx_reduce" << std::endl;
  initialize(runtime_, ctx_);
#if 0
  Runtime *runtime = CObjectWrapper::unwrap(runtime_);
  Context ctx = CObjectWrapper::unwrap(ctx_)->context();
#endif
}

#ifdef __cplusplus
}
#endif


