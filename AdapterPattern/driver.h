/* driver.h : a C-friendly interface                                       */
#pragma once
#ifdef __cplusplus // if compiler is C++, switch to C linkage
extern "C" {
#endif

typedef int (*compute_fn)(void* ctx, int); // function pointer type
compute_fn make_driver(int parameter, void** out_ctx); //equiv to model_driver_create

#ifdef __cplusplus
}  /* extern "C" */
#endif
