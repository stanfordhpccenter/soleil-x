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

#ifndef render_hpp
#define render_hpp

#include <stdio.h>
#include "legion.h"
#include "legion_c.h"

#ifdef __cplusplus
extern "C" {
#endif
  
  void cxx_preinitialize(legion_mapper_id_t mapperID);
  
  typedef struct {
    legion_index_space_t indexSpace;
    legion_logical_region_t imageX;
    legion_index_space_t colorSpace;
    legion_logical_partition_t p_Image;
    legion_field_id_t imageFields[6];
  } RegionPartition;
  
  RegionPartition cxx_initialize(
                  legion_runtime_t runtime_,
                  legion_context_t ctx_,
                  legion_mapper_id_t sampleId,
                  legion_logical_partition_t fluidPartition_);
  
  
  void cxx_render(legion_runtime_t runtime_,
                  legion_context_t ctx_,
                  legion_mapper_id_t sampleId,
                  legion_physical_region_t *fluid_,
                  legion_field_id_t* fluidFields_,
                  int numFluidFields,
                  legion_physical_region_t *particles_,
                  legion_field_id_t* particlesFields_,
                  int numParticlesFields,
                  legion_physical_region_t* image_,
                  legion_field_id_t* imageFields_,
                  int numImageFields,
                  legion_logical_partition_t fluidPartition_,
                  legion_logical_partition_t particlesPartition_,
                  int numParticlesToDraw,
                  int isosurfaceField,
                  double isosurfaceValue,
                  legion_physical_region_t *particlesToDraw_,
                  double lowerBound[3],
                  double upperBound[3]
                  );
  
  void cxx_reduce(legion_runtime_t runtime_,
                  legion_context_t ctx_,
                  legion_mapper_id_t sampleId
                  );
  
  
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
template<typename FT, int N, typename T = long long>
using AccessorRO = Legion::FieldAccessor<READ_ONLY,FT,N,T,Realm::AffineAccessor<FT,N,T> >;
template<typename FT, int N, typename T = long long>
using AccessorWO = Legion::FieldAccessor<WRITE_DISCARD,FT,N,T,Realm::AffineAccessor<FT,N,T> >;
template<typename FT, int N, typename T = long long>
using AccessorRW = Legion::FieldAccessor<READ_WRITE,FT,N,T,Realm::AffineAccessor<FT,N,T> >;
#endif


#endif /* render_hpp */
