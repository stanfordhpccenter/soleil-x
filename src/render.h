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
  
  const unsigned imageWidth = 1280;
  const unsigned imageHeight = 720;

  void cxx_preinitialize1(int numDomainNodes);
  void cxx_preinitialize2();
  
  void cxx_render(legion_runtime_t runtime_,
                  legion_context_t ctx_,
                  legion_physical_region_t *fluid_,
                  legion_physical_region_t *particles_,
                  legion_index_space_t tiles,
                  legion_logical_partition_t p_fluid,
                  legion_logical_partition_t p_particles
                  );
  
  void cxx_reduce(legion_runtime_t runtime_,
                  legion_context_t ctx_
                  );
  
  
#ifdef __cplusplus
}
#endif


#endif /* render_hpp */
