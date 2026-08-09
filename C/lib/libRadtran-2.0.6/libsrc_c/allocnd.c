#include <math.h>
#include <stdio.h>

#include "allocnd.h"
#include "common_math.h"
#include "errors.h"

/** ALLOC_ND routines to allocate ndimensional arrays
 *
 * The data block is allocated contiguously
 * and we reserve a block in front of it for pointers into the data block
 */
static void* alloc_nd_internal (const size_t Ndim, const size_t* dims, const size_t dtype_size, int zero_entries) {
  size_t dims_r[Ndim];
  size_t offsets_data[Ndim + 1];
  size_t offsets_ptr[Ndim];
  int    ierr;
  size_t Nptr; /* number of entries in prefixed pointer block */

  for (size_t i = 0; i < Ndim; ++i) {
    dims_r[i] = dims[Ndim - 1 - i];
  }
  offsets_data[0] = 1;
  ierr            = vec_cum_prod_size_t (Ndim, dims_r, offsets_data + 1);
  if (ierr)
    VALIDPTROUT (NULL, "failed setting up offsets");

  ierr = vec_cum_prod_size_t (Ndim, dims, offsets_ptr);
  if (ierr)
    VALIDPTROUT (NULL, "failed setting up offsets");

  ierr = vec_sum_size_t (Ndim - 1, offsets_ptr, &Nptr);
  if (ierr)
    VALIDPTROUT (NULL, "failed setting up Nr of pointers");

  const size_t Ndata = offsets_data[Ndim];                                // number of data values
  const size_t Npp = ceil ((float) sizeof (void*) / (float) dtype_size);  // number of data vals to have the same space as a pointer
  const size_t Nentries = Nptr * Npp + Ndata;

  // put enough asterisks there to accomodate pointer arithmetic down below, otherwise the notation is just really ugly
  void******* arr;
  if (zero_entries) {
    arr = calloc (Nentries, dtype_size);
  } else {
    arr = malloc (Nentries * dtype_size);
  }
  VALIDPTROUT (arr, "failed allocating mem");

  if (Ndim == 1) {
    return (void*) arr;
  }

  void* data = arr + Nptr;
  {
    if (Ndim == 2) {
      for (size_t i = 0; i < dims[0]; ++i) {
        const size_t idx[2] = {0, i};
        arr[i]              = data + idx1d (Ndim, idx, offsets_data) * dtype_size;
      }
      return (void*) arr;
    }
  }
  {
    void* ptr1 = arr + dims[0];
    for (size_t i = 0; i < dims[0]; ++i) {
      arr[i] = ptr1 + i * dims[1] * sizeof (void*);
    }

    if (Ndim == 3) {
      for (size_t i = 0; i < dims[0]; ++i) {
        for (size_t j = 0; j < dims[1]; ++j) {
          const size_t idx[3] = {0, j, i};
          arr[i][j]           = data + idx1d (Ndim, idx, offsets_data) * dtype_size;
        }
      }
      return (void*) arr;
    }
  }
  {
    void* ptr2 = arr + dims[0] + dims[0] * dims[1];
    for (size_t i = 0; i < dims[0]; ++i) {
      for (size_t j = 0; j < dims[1]; ++j) {
        const size_t ptr_off = i * dims[1] * dims[2] + j * dims[2];
        arr[i][j]            = ptr2 + ptr_off * sizeof (void*);
      }
    }

    if (Ndim == 4) {
      for (size_t i = 0; i < dims[0]; ++i) {
        for (size_t j = 0; j < dims[1]; ++j) {
          for (size_t k = 0; k < dims[2]; ++k) {
            const size_t idx[4] = {0, k, j, i};
            arr[i][j][k]        = data + idx1d (Ndim, idx, offsets_data) * dtype_size;
          }
        }
      }
      return (void*) arr;
    }
  }
  {
    void* ptr3 = arr + dims[0] + dims[0] * dims[1] + dims[0] * dims[1] * dims[2];
    for (size_t i = 0; i < dims[0]; ++i) {
      for (size_t j = 0; j < dims[1]; ++j) {
        for (size_t k = 0; k < dims[2]; ++k) {
          const size_t ptr_off = i * dims[1] * dims[2] * dims[3] + j * dims[2] * dims[3] + k * dims[3];
          arr[i][j][k]         = ptr3 + ptr_off * sizeof (void*);
        }
      }
    }

    if (Ndim == 5) {
      for (size_t i = 0; i < dims[0]; ++i) {
        for (size_t j = 0; j < dims[1]; ++j) {
          for (size_t k = 0; k < dims[2]; ++k) {
            for (size_t l = 0; l < dims[3]; ++l) {
              const size_t idx[5] = {0, l, k, j, i};
              arr[i][j][k][l]     = data + idx1d (Ndim, idx, offsets_data) * dtype_size;
            }
          }
        }
      }
      return (void*) arr;
    }
  }

  {
    void* ptr4 = arr + dims[0] + dims[0] * dims[1] + dims[0] * dims[1] * dims[2] + dims[0] * dims[1] * dims[2] * dims[3];
    for (size_t i = 0; i < dims[0]; ++i) {
      for (size_t j = 0; j < dims[1]; ++j) {
        for (size_t k = 0; k < dims[2]; ++k) {
          for (size_t l = 0; l < dims[3]; ++l) {
            const size_t ptr_off =
              i * dims[1] * dims[2] * dims[3] * dims[4] + j * dims[2] * dims[3] * dims[4] + k * dims[3] * dims[4] + l * dims[4];
            arr[i][j][k][l] = ptr4 + ptr_off * sizeof (void*);
          }
        }
      }
    }

    if (Ndim == 6) {
      for (size_t i = 0; i < dims[0]; ++i) {
        for (size_t j = 0; j < dims[1]; ++j) {
          for (size_t k = 0; k < dims[2]; ++k) {
            for (size_t l = 0; l < dims[3]; ++l) {
              for (size_t m = 0; m < dims[4]; ++m) {
                const size_t idx[6] = {0, m, l, k, j, i};
                arr[i][j][k][l][m]  = data + idx1d (Ndim, idx, offsets_data) * dtype_size;
              }
            }
          }
        }
      }
      return (void*) arr;
    }
  }

  {
    void* ptr5 = arr + dims[0] + dims[0] * dims[1] + dims[0] * dims[1] * dims[2] + dims[0] * dims[1] * dims[2] * dims[3] +
                 dims[0] * dims[1] * dims[2] * dims[3] * dims[4];
    for (size_t i = 0; i < dims[0]; ++i) {
      for (size_t j = 0; j < dims[1]; ++j) {
        for (size_t k = 0; k < dims[2]; ++k) {
          for (size_t l = 0; l < dims[3]; ++l) {
            for (size_t m = 0; m < dims[4]; ++m) {
              const size_t ptr_off = i * dims[1] * dims[2] * dims[3] * dims[4] * dims[5] +
                                     j * dims[2] * dims[3] * dims[4] * dims[5] + k * dims[3] * dims[4] * dims[5] +
                                     l * dims[4] * dims[5] + m * dims[5];
              arr[i][j][k][l][m] = ptr5 + ptr_off * sizeof (void*);
            }
          }
        }
      }
    }

    if (Ndim == 7) {
      for (size_t i = 0; i < dims[0]; ++i) {
        for (size_t j = 0; j < dims[1]; ++j) {
          for (size_t k = 0; k < dims[2]; ++k) {
            for (size_t l = 0; l < dims[3]; ++l) {
              for (size_t m = 0; m < dims[4]; ++m) {
                for (size_t n = 0; n < dims[5]; ++n) {
                  const size_t idx[7]   = {0, n, m, l, k, j, i};
                  arr[i][j][k][l][m][n] = data + idx1d (Ndim, idx, offsets_data) * dtype_size;
                }
              }
            }
          }
        }
      }
      return (void*) arr;
    }
  }

  VALIDPTROUT (NULL, "not implemented for higher dimensions");
  return NULL;
}

#define ALLOC_1D(type, name)                                                                                                       \
  type* malloc_##name##_1D (size_t s1, const char* varname) {                                                                      \
    const size_t dims[1] = {s1};                                                                                                   \
    type*        arr     = (type*) alloc_nd_internal (1, dims, sizeof (type), 0);                                                  \
    VALIDPTROUT (arr, varname);                                                                                                    \
    return arr;                                                                                                                    \
  }                                                                                                                                \
  type* calloc_##name##_1D (size_t s1, const char* varname) {                                                                      \
    const size_t dims[1] = {s1};                                                                                                   \
    type*        arr     = (type*) alloc_nd_internal (1, dims, sizeof (type), 1);                                                  \
    VALIDPTROUT (arr, varname);                                                                                                    \
    return arr;                                                                                                                    \
  }                                                                                                                                \
  void free_##name##_1D (type* data) { free (data); }

#define ALLOC_2D(type, name)                                                                                                       \
  type** malloc_##name##_2D (size_t s1, size_t s2, const char* varname) {                                                          \
    const size_t dims[2] = {s1, s2};                                                                                               \
    type**       arr     = (type**) alloc_nd_internal (2, dims, sizeof (type), 0);                                                 \
    VALIDPTROUT (arr, varname);                                                                                                    \
    return arr;                                                                                                                    \
  }                                                                                                                                \
  type** calloc_##name##_2D (size_t s1, size_t s2, const char* varname) {                                                          \
    const size_t dims[2] = {s1, s2};                                                                                               \
    type**       arr     = (type**) alloc_nd_internal (2, dims, sizeof (type), 1);                                                 \
    VALIDPTROUT (arr, varname);                                                                                                    \
    return arr;                                                                                                                    \
  }                                                                                                                                \
  void free_##name##_2D (type** data) { free (data); }

#define ALLOC_3D(type, name)                                                                                                       \
  type*** malloc_##name##_3D (size_t s1, size_t s2, const size_t s3, const char* varname) {                                        \
    const size_t dims[3] = {s1, s2, s3};                                                                                           \
    type***      arr     = (type***) alloc_nd_internal (3, dims, sizeof (type), 0);                                                \
    VALIDPTROUT (arr, varname);                                                                                                    \
    return arr;                                                                                                                    \
  }                                                                                                                                \
  type*** calloc_##name##_3D (size_t s1, size_t s2, const size_t s3, const char* varname) {                                        \
    const size_t dims[3] = {s1, s2, s3};                                                                                           \
    type***      arr     = (type***) alloc_nd_internal (3, dims, sizeof (type), 1);                                                \
    VALIDPTROUT (arr, varname);                                                                                                    \
    return arr;                                                                                                                    \
  }                                                                                                                                \
  void free_##name##_3D (type*** data) { free (data); }

#define ALLOC_4D(type, name)                                                                                                       \
  type**** malloc_##name##_4D (size_t s1, size_t s2, size_t s3, size_t s4, const char* varname) {                                  \
    const size_t dims[4] = {s1, s2, s3, s4};                                                                                       \
    type****     arr     = (type****) alloc_nd_internal (4, dims, sizeof (type), 0);                                               \
    VALIDPTROUT (arr, varname);                                                                                                    \
    return arr;                                                                                                                    \
  }                                                                                                                                \
  type**** calloc_##name##_4D (size_t s1, size_t s2, size_t s3, size_t s4, const char* varname) {                                  \
    const size_t dims[4] = {s1, s2, s3, s4};                                                                                       \
    type****     arr     = (type****) alloc_nd_internal (4, dims, sizeof (type), 1);                                               \
    VALIDPTROUT (arr, varname);                                                                                                    \
    return arr;                                                                                                                    \
  }                                                                                                                                \
  void free_##name##_4D (type**** data) { free (data); }

#define ALLOC_5D(type, name)                                                                                                       \
  type***** malloc_##name##_5D (size_t s1, size_t s2, size_t s3, size_t s4, size_t s5, const char* varname) {                      \
    const size_t dims[5] = {s1, s2, s3, s4, s5};                                                                                   \
    type*****    arr     = (type*****) alloc_nd_internal (5, dims, sizeof (type), 0);                                              \
    VALIDPTROUT (arr, varname);                                                                                                    \
    return arr;                                                                                                                    \
  }                                                                                                                                \
  type***** calloc_##name##_5D (size_t s1, size_t s2, size_t s3, size_t s4, size_t s5, const char* varname) {                      \
    const size_t dims[5] = {s1, s2, s3, s4, s5};                                                                                   \
    type*****    arr     = (type*****) alloc_nd_internal (5, dims, sizeof (type), 1);                                              \
    VALIDPTROUT (arr, varname);                                                                                                    \
    return arr;                                                                                                                    \
  }                                                                                                                                \
  void free_##name##_5D (type***** data) { free (data); }

#define ALLOC_6D(type, name)                                                                                                       \
  type****** malloc_##name##_6D (size_t s1, size_t s2, size_t s3, size_t s4, size_t s5, size_t s6, const char* varname) {          \
    const size_t dims[6] = {s1, s2, s3, s4, s5, s6};                                                                               \
    type******   arr     = (type******) alloc_nd_internal (6, dims, sizeof (type), 0);                                             \
    VALIDPTROUT (arr, varname);                                                                                                    \
    return arr;                                                                                                                    \
  }                                                                                                                                \
  type****** calloc_##name##_6D (size_t s1, size_t s2, size_t s3, size_t s4, size_t s5, size_t s6, const char* varname) {          \
    const size_t dims[6] = {s1, s2, s3, s4, s5, s6};                                                                               \
    type******   arr     = (type******) alloc_nd_internal (6, dims, sizeof (type), 1);                                             \
    VALIDPTROUT (arr, varname);                                                                                                    \
    return arr;                                                                                                                    \
  }                                                                                                                                \
  void free_##name##_6D (type****** data) { free (data); }

#define ALLOC_7D(type, name)                                                                                                       \
  type*******                                                                                                                      \
    malloc_##name##_7D (size_t s1, size_t s2, size_t s3, size_t s4, size_t s5, size_t s6, size_t s7, const char* varname) {        \
    const size_t dims[7] = {s1, s2, s3, s4, s5, s6, s7};                                                                           \
    type*******  arr     = (type*******) alloc_nd_internal (7, dims, sizeof (type), 0);                                            \
    VALIDPTROUT (arr, varname);                                                                                                    \
    return arr;                                                                                                                    \
  }                                                                                                                                \
  type*******                                                                                                                      \
    calloc_##name##_7D (size_t s1, size_t s2, size_t s3, size_t s4, size_t s5, size_t s6, size_t s7, const char* varname) {        \
    const size_t dims[7] = {s1, s2, s3, s4, s5, s6, s7};                                                                           \
    type*******  arr     = (type*******) alloc_nd_internal (7, dims, sizeof (type), 1);                                            \
    VALIDPTROUT (arr, varname);                                                                                                    \
    return arr;                                                                                                                    \
  }                                                                                                                                \
  void free_##name##_7D (type******* data) { free (data); }

/* special case for a ragged array in the last dimension with offsets in 4th */
#define CALLOC_5D_RESTRICTED(type, name)                                                                                           \
  type***** calloc_##name##_5D_restricted (const size_t s1,                                                                        \
                                           const size_t s2,                                                                        \
                                           const size_t s3,                                                                        \
                                           const size_t s4_lower,                                                                  \
                                           const size_t s4_upper,                                                                  \
                                           const int*   s5,                                                                        \
                                           const char*  varname) {                                                                  \
    const size_t s4     = s4_upper - s4_lower + 1;                                                                                 \
    size_t       s5_max = 0;                                                                                                       \
    for (size_t k = s4_lower; k <= s4_upper; ++k) {                                                                                \
      s5_max = MAX (s5[k], s5_max);                                                                                                \
    }                                                                                                                              \
    type***** p = malloc_##name##_5D (s1, s2, s3, s4, s5_max, varname);                                                            \
    VALIDPTR (p);                                                                                                                  \
    const size_t dimlen[5]  = {s1, s2, s3, s4, s5_max};                                                                            \
    const int    offsets[5] = {0, 0, 0, -s4_lower, 0};                                                                             \
    const int    ierr       = apply_array_offsets_##name##_5D (dimlen, offsets, &p);                                               \
    if (ierr) {                                                                                                                    \
      p = NULL;                                                                                                                    \
    }                                                                                                                              \
    VALIDPTR (p);                                                                                                                  \
    /* zeroing valid entries */                                                                                                    \
    for (size_t i1 = 0; i1 < s1; ++i1) {                                                                                           \
      for (size_t i2 = 0; i2 < s2; ++i2) {                                                                                         \
        for (size_t i3 = 0; i3 < s3; ++i3) {                                                                                       \
          for (size_t i4 = s4_lower; i4 <= s4_upper; ++i4) {                                                                       \
            for (size_t i5 = 0; i5 < s5[i4]; ++i5) {                                                                               \
              p[i1][i2][i3][i4][i5] = 0;                                                                                           \
            }                                                                                                                      \
          }                                                                                                                        \
        }                                                                                                                          \
      }                                                                                                                            \
    }                                                                                                                              \
    return p;                                                                                                                      \
  }                                                                                                                                \
  void free_##name##_5D_restricted (type***** data) { free_##name##_5D (data); }

#define ALLOC_ND(type, name)                                                                                                       \
  ALLOC_1D (type, name)                                                                                                            \
  ALLOC_2D (type, name)                                                                                                            \
  ALLOC_3D (type, name)                                                                                                            \
  ALLOC_4D (type, name)                                                                                                            \
  ALLOC_5D (type, name)                                                                                                            \
  ALLOC_6D (type, name)                                                                                                            \
  ALLOC_7D (type, name)                                                                                                            \
  CALLOC_5D_RESTRICTED (type, name)

ALLOC_ND (char, char);
ALLOC_ND (unsigned char, uchar);
ALLOC_ND (int, int);
ALLOC_ND (size_t, size_t);
ALLOC_ND (float, float);
ALLOC_ND (double, double);

float*** calloc_float_3D_sparse (size_t nz, size_t nx, size_t ny, int* allocz, const char* varname) {
  size_t sz = 0;
  for (size_t i = 0; i < nz; ++i) {
    if (allocz[i])
      ++sz;
  }
  float** p1 = calloc_float_2D (sz * nx, ny, varname);
  VALIDPTR (p1);
  float*** p2 = malloc (nz * sizeof (float**));
  if (!p2) {
    free_float_2D (p1);
    fprintf (stderr, "Error allocating memory for 2st pointer level of %s\n", varname);
    VALIDPTR (p2);
  }
  size_t i_sparse = 0;
  for (size_t i = 0; i < nz; ++i) {
    if (allocz[i]) {
      p2[i] = p1 + i_sparse * nx;
      ++i_sparse;
    } else {
      p2[i] = NULL;
    }
  }
  return p2;
}

/** Routines to add offsets to N-dimensional arrays
 * Funtion: apply_array_offsets_*_*D
 *
 * Description:
 *   Offset the adresses by a given amount in each dimension
 *   to allow for negative indices
 *   Note that the offsets have to be reverted before freeing the memory!
 *
 * Parameters:
 *   * dimlen  - array of each dimension size, itself has to be length (Ndim-1)
 *   * offsets - array with offset length in each dimension size
 *   * arr     - pointer to an address construct
 *
 * Example:
 *   given a 1D array with size N
 *
 *   int * arr = calloc(N, sizeof(int);
 *
 *   we can access the pointer beginning at [0] upt to [N-1]
 *
 *   The routines apply_array_offsets_* allow to shift the adresses
 *   so that the valid range is a different one.
 *   e.g. after:
 *
 *   apply_array_offsets_int_1D({N}, {-2}, &arr);
 *   the valid range is move to [2] up to [N+1]
 *   i.e. we shifted the entries by two
 *
 */
#define APPLY_ARRAY_OFFSETS_1D(type, name)                                                                                         \
  int apply_array_offsets_##name##_1D (const size_t* dimlen, const int* offsets, type** arr) {                                     \
    CHKPOINTER (*arr);                                                                                                             \
    (*arr) += offsets[0];                                                                                                          \
    return 0;                                                                                                                      \
  }

#define APPLY_ARRAY_OFFSETS_2D(type, name)                                                                                         \
  int apply_array_offsets_##name##_2D (const size_t* dimlen, const int* offsets, type*** arr) {                                    \
    CHKPOINTER (*arr);                                                                                                             \
    for (size_t i = 0; i < dimlen[0]; ++i) {                                                                                       \
      CHKPOINTER ((*arr)[i]);                                                                                                      \
      (*arr)[i] += offsets[1];                                                                                                     \
    }                                                                                                                              \
    (*arr) += offsets[0];                                                                                                          \
    return 0;                                                                                                                      \
  }

#define APPLY_ARRAY_OFFSETS_3D(type, name)                                                                                         \
  int apply_array_offsets_##name##_3D (const size_t* dimlen, const int* offsets, type**** arr) {                                   \
    CHKPOINTER (*arr);                                                                                                             \
    for (size_t i = 0; i < dimlen[0]; ++i) {                                                                                       \
      CHKPOINTER ((*arr)[i]);                                                                                                      \
      for (size_t j = 0; j < dimlen[1]; ++j) {                                                                                     \
        CHKPOINTER ((*arr)[i][j]);                                                                                                 \
        (*arr)[i][j] += offsets[2];                                                                                                \
      }                                                                                                                            \
      (*arr)[i] += offsets[1];                                                                                                     \
    }                                                                                                                              \
    (*arr) += offsets[0];                                                                                                          \
    return 0;                                                                                                                      \
  }

#define APPLY_ARRAY_OFFSETS_4D(type, name)                                                                                         \
  int apply_array_offsets_##name##_4D (const size_t* dimlen, const int* offsets, type***** arr) {                                  \
    CHKPOINTER (*arr);                                                                                                             \
    for (size_t i = 0; i < dimlen[0]; ++i) {                                                                                       \
      CHKPOINTER ((*arr)[i]);                                                                                                      \
      for (size_t j = 0; j < dimlen[1]; ++j) {                                                                                     \
        CHKPOINTER ((*arr)[i][j]);                                                                                                 \
        for (size_t k = 0; k < dimlen[2]; ++k) {                                                                                   \
          CHKPOINTER ((*arr)[i][j][k]);                                                                                            \
          (*arr)[i][j][k] += offsets[3];                                                                                           \
        }                                                                                                                          \
        (*arr)[i][j] += offsets[2];                                                                                                \
      }                                                                                                                            \
      (*arr)[i] += offsets[1];                                                                                                     \
    }                                                                                                                              \
    (*arr) += offsets[0];                                                                                                          \
    return 0;                                                                                                                      \
  }

#define APPLY_ARRAY_OFFSETS_5D(type, name)                                                                                         \
  int apply_array_offsets_##name##_5D (const size_t* dimlen, const int* offsets, type****** arr) {                                 \
    CHKPOINTER (*arr);                                                                                                             \
    for (size_t i = 0; i < dimlen[0]; ++i) {                                                                                       \
      CHKPOINTER ((*arr)[i]);                                                                                                      \
      for (size_t j = 0; j < dimlen[1]; ++j) {                                                                                     \
        CHKPOINTER ((*arr)[i][j]);                                                                                                 \
        for (size_t k = 0; k < dimlen[2]; ++k) {                                                                                   \
          CHKPOINTER ((*arr)[i][j][k]);                                                                                            \
          for (size_t l = 0; l < dimlen[3]; ++l) {                                                                                 \
            CHKPOINTER ((*arr)[i][j][k][l]);                                                                                       \
            (*arr)[i][j][k][l] += offsets[4];                                                                                      \
          }                                                                                                                        \
          (*arr)[i][j][k] += offsets[3];                                                                                           \
        }                                                                                                                          \
        (*arr)[i][j] += offsets[2];                                                                                                \
      }                                                                                                                            \
      (*arr)[i] += offsets[1];                                                                                                     \
    }                                                                                                                              \
    (*arr) += offsets[0];                                                                                                          \
    return 0;                                                                                                                      \
  }

#define APPLY_ARRAY_OFFSETS_6D(type, name)                                                                                         \
  int apply_array_offsets_##name##_6D (const size_t* dimlen, const int* offsets, type******* arr) {                                \
    CHKPOINTER (*arr);                                                                                                             \
    for (size_t i = 0; i < dimlen[0]; ++i) {                                                                                       \
      CHKPOINTER ((*arr)[i]);                                                                                                      \
      for (size_t j = 0; j < dimlen[1]; ++j) {                                                                                     \
        CHKPOINTER ((*arr)[i][j]);                                                                                                 \
        for (size_t k = 0; k < dimlen[2]; ++k) {                                                                                   \
          CHKPOINTER ((*arr)[i][j][k]);                                                                                            \
          for (size_t l = 0; l < dimlen[3]; ++l) {                                                                                 \
            CHKPOINTER ((*arr)[i][j][k][l]);                                                                                       \
            for (size_t m = 0; m < dimlen[4]; ++m) {                                                                               \
              CHKPOINTER ((*arr)[i][j][k][l][m]);                                                                                  \
              (*arr)[i][j][k][l][m] += offsets[5];                                                                                 \
            }                                                                                                                      \
            (*arr)[i][j][k][l] += offsets[4];                                                                                      \
          }                                                                                                                        \
          (*arr)[i][j][k] += offsets[3];                                                                                           \
        }                                                                                                                          \
        (*arr)[i][j] += offsets[2];                                                                                                \
      }                                                                                                                            \
      (*arr)[i] += offsets[1];                                                                                                     \
    }                                                                                                                              \
    (*arr) += offsets[0];                                                                                                          \
    return 0;                                                                                                                      \
  }

#define APPLY_ARRAY_OFFSETS_7D(type, name)                                                                                         \
  int apply_array_offsets_##name##_7D (const size_t* dimlen, const int* offsets, type******** arr) {                               \
    CHKPOINTER (*arr);                                                                                                             \
    for (size_t i = 0; i < dimlen[0]; ++i) {                                                                                       \
      CHKPOINTER ((*arr)[i]);                                                                                                      \
      for (size_t j = 0; j < dimlen[1]; ++j) {                                                                                     \
        CHKPOINTER ((*arr)[i][j]);                                                                                                 \
        for (size_t k = 0; k < dimlen[2]; ++k) {                                                                                   \
          CHKPOINTER ((*arr)[i][j][k]);                                                                                            \
          for (size_t l = 0; l < dimlen[3]; ++l) {                                                                                 \
            CHKPOINTER ((*arr)[i][j][k][l]);                                                                                       \
            for (size_t m = 0; m < dimlen[4]; ++m) {                                                                               \
              CHKPOINTER ((*arr)[i][j][k][l][m]);                                                                                  \
              for (size_t n = 0; n < dimlen[5]; ++n) {                                                                             \
                CHKPOINTER ((*arr)[i][j][k][l][m][n]);                                                                             \
                (*arr)[i][j][k][l][m][n] += offsets[6];                                                                            \
              }                                                                                                                    \
              (*arr)[i][j][k][l][m] += offsets[5];                                                                                 \
            }                                                                                                                      \
            (*arr)[i][j][k][l] += offsets[4];                                                                                      \
          }                                                                                                                        \
          (*arr)[i][j][k] += offsets[3];                                                                                           \
        }                                                                                                                          \
        (*arr)[i][j] += offsets[2];                                                                                                \
      }                                                                                                                            \
      (*arr)[i] += offsets[1];                                                                                                     \
    }                                                                                                                              \
    (*arr) += offsets[0];                                                                                                          \
    return 0;                                                                                                                      \
  }

#define APPLY_ARRAY_OFFSETS_ND(type, name)                                                                                         \
  APPLY_ARRAY_OFFSETS_1D (type, name)                                                                                              \
  APPLY_ARRAY_OFFSETS_2D (type, name)                                                                                              \
  APPLY_ARRAY_OFFSETS_3D (type, name)                                                                                              \
  APPLY_ARRAY_OFFSETS_4D (type, name)                                                                                              \
  APPLY_ARRAY_OFFSETS_5D (type, name)                                                                                              \
  APPLY_ARRAY_OFFSETS_6D (type, name)                                                                                              \
  APPLY_ARRAY_OFFSETS_7D (type, name)

APPLY_ARRAY_OFFSETS_ND (char, char);
APPLY_ARRAY_OFFSETS_ND (unsigned char, uchar);
APPLY_ARRAY_OFFSETS_ND (int, int);
APPLY_ARRAY_OFFSETS_ND (size_t, size_t);
APPLY_ARRAY_OFFSETS_ND (float, float);
APPLY_ARRAY_OFFSETS_ND (double, double);

/**
 * Description:
 *   Local convenience functions to allocate memory for 3D results like output->rfldir3d_r
 *
 *   One key aspect here is that we shift the indices to the left so that arrays which have dimension:
 *
 *     [wl.nlambda_r, output->mc.sample.Nx, output->mc.sample.Ny, output->atm.nzout]
 *
 *   get offsetted by:
 *     [0, out->islower, out->jslower, wl.iv_lower]
 *
 *   and then have valid index ranges:
 *    [0               , out->atm.nzout]
 *    [out->islower    , out->isupper]
 *    [output->jslower , output->jsupper]
 *    [out->wl.iv_lower, out->wl.iv_upper]
 */
int alloc_array_with_offsets_4d (const size_t* dims, const int* offsets, const char* varname, float***** var) {
  *var = calloc_float_4D (dims[0], dims[1], dims[2], dims[3], varname);
  CHKPOINTER (*var);
  int status = apply_array_offsets_float_4D (dims, offsets, var);
  CHKERR (status);
  return 0;
}
int alloc_array_with_offsets_5d (const size_t* dims, const int* offsets, const char* varname, float****** var) {
  *var = calloc_float_5D (dims[0], dims[1], dims[2], dims[3], dims[4], varname);
  CHKPOINTER (*var);
  int status = apply_array_offsets_float_5D (dims, offsets, var);
  CHKERR (status);
  return 0;
}
int alloc_array_with_offsets_6d (const size_t* dims, const int* offsets, const char* varname, float******* var) {
  *var = calloc_float_6D (dims[0], dims[1], dims[2], dims[3], dims[4], dims[5], varname);
  CHKPOINTER (*var);
  int status = apply_array_offsets_float_6D (dims, offsets, var);
  CHKERR (status);
  return 0;
}
int alloc_array_with_offsets_7d (const size_t* dims, const int* offsets, const char* varname, float******** var) {
  *var = calloc_float_7D (dims[0], dims[1], dims[2], dims[3], dims[4], dims[5], dims[6], varname);
  CHKPOINTER (*var);
  int status = apply_array_offsets_float_7D (dims, offsets, var);
  CHKERR (status);
  return 0;
}
