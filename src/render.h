//
// render.h
//

#include "legion.h"
#include "legion/legion_c.h"


#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
  double from[3];
  double at[3];
  double up[3];
} Camera;

typedef double FieldData;
typedef struct {
  FieldData x[3];
} FieldData3;

  void cxx_preinitialize(legion_mapper_id_t mapperID);

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
                     );

  void cxx_render(legion_runtime_t runtime_, legion_context_t ctx_,
   double camera[9], double colorScale[2]);

  void cxx_reduce(legion_context_t ctx_, double camera[9]);

  void cxx_saveImage(legion_runtime_t runtime_, legion_context_t ctx_, const char* outDir);

  void cxx_saveIndividualImages(legion_runtime_t runtime_, legion_context_t ctx_, const char* outDir);

  void cxx_terminate();

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
