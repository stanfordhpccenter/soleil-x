#ifndef __DOM_MAPPER_H__
#define __DOM_MAPPER_H__

#ifdef __cplusplus
extern "C" {
#endif

void register_mappers();

void register_sharding_functor(void *runtime, void *config);

#ifdef __cplusplus
}
#endif

#endif // __DOM_MAPPER_H__
