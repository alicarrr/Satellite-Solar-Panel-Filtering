#ifndef ALLOCND_H
#define ALLOCND_H

#if defined(__cplusplus)
extern "C" {
#endif

#include <stdlib.h>

#define UVSPEC_FREE(p)                                                                                                             \
  free (p);                                                                                                                        \
  p = NULL

#define ALLOC_ND(type, name)                                                                                                       \
  type*      malloc_##name##_1D (size_t sx, const char* varname);                                                                  \
  type**     malloc_##name##_2D (size_t sx, size_t sy, const char* varname);                                                       \
  type***    malloc_##name##_3D (size_t sx, size_t sy, size_t sz, const char* varname);                                            \
  type****   malloc_##name##_4D (size_t sx, size_t sy, size_t sz, size_t dim4, const char* varname);                               \
  type*****  malloc_##name##_5D (size_t sx, size_t sy, size_t sz, size_t dim4, size_t dim5, const char* varname);                  \
  type****** malloc_##name##_6D (size_t sx, size_t sy, size_t sz, size_t dim4, size_t dim5, size_t dim6, const char* varname);     \
  type*******                                                                                                                      \
    malloc_##name##_7D (size_t sx, size_t sy, size_t sz, size_t dim4, size_t dim5, size_t dim6, size_t dim7, const char* varname); \
                                                                                                                                   \
  type*      calloc_##name##_1D (size_t sx, const char* varname);                                                                  \
  type**     calloc_##name##_2D (size_t sx, size_t sy, const char* varname);                                                       \
  type***    calloc_##name##_3D (size_t sx, size_t sy, size_t sz, const char* varname);                                            \
  type****   calloc_##name##_4D (size_t sx, size_t sy, size_t sz, size_t dim4, const char* varname);                               \
  type*****  calloc_##name##_5D (size_t sx, size_t sy, size_t sz, size_t dim4, size_t dim5, const char* varname);                  \
  type****** calloc_##name##_6D (size_t sx, size_t sy, size_t sz, size_t dim4, size_t dim5, size_t dim6, const char* varname);     \
  type*******                                                                                                                      \
    calloc_##name##_7D (size_t sx, size_t sy, size_t sz, size_t dim4, size_t dim5, size_t dim6, size_t dim7, const char* varname); \
                                                                                                                                   \
  void free_##name##_1D (type* data);                                                                                              \
  void free_##name##_2D (type** data);                                                                                             \
  void free_##name##_3D (type*** data);                                                                                            \
  void free_##name##_4D (type**** data);                                                                                           \
  void free_##name##_5D (type***** data);                                                                                          \
  void free_##name##_6D (type****** data);                                                                                         \
  void free_##name##_7D (type******* data);                                                                                        \
                                                                                                                                   \
  type***** calloc_##name##_5D_restricted (const size_t s1,                                                                        \
                                           const size_t s2,                                                                        \
                                           const size_t s3,                                                                        \
                                           const size_t s4_lower,                                                                  \
                                           const size_t s4_upper,                                                                  \
                                           const int*   s5,                                                                        \
                                           const char*  varname);                                                                   \
  void      free_##name##_5D_restricted (type***** data);

ALLOC_ND (char, char)
ALLOC_ND (unsigned char, uchar)
ALLOC_ND (int, int)
ALLOC_ND (size_t, size_t)
ALLOC_ND (float, float)
ALLOC_ND (double, double)
#undef ALLOC_ND

float*** calloc_float_3D_sparse (size_t nz, size_t nx, size_t ny, int* allocz, const char* varname);

#define APPLY_ARRAY_OFFSETS_ND(type, name)                                                                                         \
  int apply_array_offsets_##name##_1D (const size_t* dimlen, const int* offsets, type** arr);                                      \
  int apply_array_offsets_##name##_2D (const size_t* dimlen, const int* offsets, type*** arr);                                     \
  int apply_array_offsets_##name##_3D (const size_t* dimlen, const int* offsets, type**** arr);                                    \
  int apply_array_offsets_##name##_4D (const size_t* dimlen, const int* offsets, type***** arr);                                   \
  int apply_array_offsets_##name##_5D (const size_t* dimlen, const int* offsets, type****** arr);                                  \
  int apply_array_offsets_##name##_6D (const size_t* dimlen, const int* offsets, type******* arr);                                 \
  int apply_array_offsets_##name##_7D (const size_t* dimlen, const int* offsets, type******** arr);

APPLY_ARRAY_OFFSETS_ND (int, int)
APPLY_ARRAY_OFFSETS_ND (size_t, size_t)
APPLY_ARRAY_OFFSETS_ND (char, char)
APPLY_ARRAY_OFFSETS_ND (unsigned char, uchar)
APPLY_ARRAY_OFFSETS_ND (float, float)
APPLY_ARRAY_OFFSETS_ND (double, double)
#undef APPLY_ARRAY_OFFSETS_ND

int alloc_array_with_offsets_4d (const size_t* dims, const int* offsets, const char* varname, float***** var);
int alloc_array_with_offsets_5d (const size_t* dims, const int* offsets, const char* varname, float****** var);
int alloc_array_with_offsets_6d (const size_t* dims, const int* offsets, const char* varname, float******* var);
int alloc_array_with_offsets_7d (const size_t* dims, const int* offsets, const char* varname, float******** var);

#if defined(__cplusplus)
}
#endif

#endif
