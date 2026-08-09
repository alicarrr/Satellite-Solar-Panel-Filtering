#include <stdio.h>
#include <float.h>
#include <math.h>

#include "common_math.h"
#include "errors.h"

/***********************************************************************************/
/* Functions: v_mult, v_sub, vn_mult                                      @62_30i@ */
/* Description:                                                                    */
/*   short cuts for vector*vector -> number, vector - vector -> vector,            */
/*                  vector*number -> vector                                        */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/
int v_cross_product_norm (const double a[3], const double b[3], double c[3]) {
  double n = 0.0;
  int    i = 0;

  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];

  /* Normalization*/
  for (i = 0; i < 3; i++)
    n += c[i] * c[i];

  /* a and b are parallel */
  if (n == 0.0)
    return 1;

  n = sqrt (n);

  for (i = 0; i < 3; i++)
    c[i] /= n;

  return 0;
}

/* const problem in C for a, see header file for details! */
void mat_v_mult (size_t N, double** a, const double* b, double* c) {
  int i = 0, j = 0;

  for (i = 0; i < N; i++) {
    c[i] = 0.0;
    for (j = 0; j < N; j++)
      c[i] += a[i][j] * b[j];
  }
}

void v_mult (const double a[3], const double b[3], double* c) {
  *c = vec3_dot_product (a, b);
}

void v_sub (const double a[3], const double b[3], double c[3]) {
  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
  c[2] = a[2] - b[2];
}

void v_add (const double a[3], const double b[3], double c[3]) {
  c[0] = a[0] + b[0];
  c[1] = a[1] + b[1];
  c[2] = a[2] + b[2];
}

void vn_mult (const double a[3], const double* b, double c[3]) {
  c[0] = a[0] * *b;
  c[1] = a[1] * *b;
  c[2] = a[2] * *b;
}

// Cross product, right hand rule, a(thumb), b(pointing finger)
void v_cross_product (const double a[3], const double b[3], double c[3]) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

int v_mult_mu (const double a[3], const double b[3], double* c) {
  *c = vec3_dot_product (a, b);

  if (*c > 1.0) {
    if (*c > 1.0 + DBL_EPSILON)
      *c = NOT_A_NUMBER;
    *c = 1.0;
  }

  if (*c < -1.0) {
    if (*c < -1.0 - DBL_EPSILON)
      *c = NOT_A_NUMBER;
    *c = -1.0;
  }

  if (isnan (*c)) {
    fprintf (stderr, "Error! v_mult_mu([%f %f %f] [%f %f %f]) has unphysical result %e \n", a[0], a[1], a[2], b[0], b[1], b[2], *c);
    CHKERR (-1);
  }
  return 0;
}

void v_neg (double* a) {
  a[0] = -a[0];
  a[1] = -a[1];
  a[2] = -a[2];
}

// constrains angle in degress to range [0,360]
inline double constrainAngle (double x) {
  x = fmod (x, 360);
  if (x < 0)
    x += 360;
  return x;
}

inline double cosd (double angle) {
  const double deg2rad = M_PI / 180.0;
  return cos (angle * deg2rad);
}

inline double sind (double angle) {
  const double deg2rad = M_PI / 180.0;
  return sin (angle * deg2rad);
}

// perform a optimized computation of cosd and sind at the same time (this is about 30% percent faster than calling them separately)
inline int cosd_sind (double angle, double* cos_val, double* sin_val) {
  *sin_val             = sind (angle);
  const double c_angle = constrainAngle (angle);
  if (c_angle >= 90 && c_angle <= 270) {
    *cos_val = -sqrt (1.0 - (*sin_val) * (*sin_val));
  } else {
    *cos_val = sqrt (1.0 - (*sin_val) * (*sin_val));
  }
  return 0;
}

inline double tand (double angle) {
  // tan is extremely slow in current gcc version, CE 180824
  // return tan ( angle * PI / 180.0 );
  return sind (angle) / cosd (angle);
}

double acosd (double rad) {
  return acos (rad) * 180.0 / M_PI;
}

double asind (double rad) {
  return asin (rad) * 180.0 / M_PI;
}

double atand (double rad) {
  return atan (rad) * 180.0 / M_PI;
}

/* return -1 or +1 according to the sign of floating point variable x */
int sign (double x) {
  if (x < 0)
    return -1;
  return 1;
}

/* applies out[i] = fabs(inp[i]) on all elements of inp[N] */
void vec_fabs (const size_t N, const double* inp, double* out) {
  for (size_t i = 0; i < N; ++i)
    out[i] = fabs (inp[i]);
}

/* Dot product of two vecs a, b with size N */
int vec_dot_product (const size_t N, const double* a, const double* b, double* c) {
  *c = 0;
  for (size_t i = 0; i < N; ++i)
    *c += a[i] * b[i];
  return 0;
}

double vec3_dot_product (const double a[3], const double b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/* Compute 2-Norm of vector */
int vec_norm2 (const size_t N, const double* v, double* norm) {
  int ierr = vec_dot_product (N, v, v, norm);
  CHKERR (ierr);
  CHKERR (*norm < 0);
  *norm = sqrt (*norm);
  return 0;
}
double vec3_norm2 (const double v[3]) {
  return sqrt (vec3_dot_product (v, v));
}

/* Computes y = alpha y */
int vec3_scale (const double alpha, double y[3]) {
  if (alpha == 1)
    return 0;
  y[0] = alpha * y[0];
  y[1] = alpha * y[1];
  y[2] = alpha * y[2];
  return 0;
}

int vec3_normalize (double x[3]) {
  double norm = vec3_norm2 (x);
  if (norm <= DBL_MIN)
    return -1;
  return vec3_scale (1. / norm, x);
}

int vec_scale (size_t N, const double alpha, double* y) {
  if (alpha == 1)
    return 0;
  for (size_t i = 0; i < N; ++i)
    y[i] = alpha * y[i];
  return 0;
}

#define TYPED_VEC_FUNCS(type)                                                                                                      \
  int vec_sum_##type (const size_t N, const type* y, type* sum) {                                                                  \
    *sum = 0;                                                                                                                      \
    for (size_t i = 0; i < N; ++i)                                                                                                 \
      *sum += y[i];                                                                                                                \
    return 0;                                                                                                                      \
  }                                                                                                                                \
  int vec_prod_sum_##type (const size_t N, const type* y, type* sum) {                                                             \
    if (N < 1) {                                                                                                                   \
      *sum = 0;                                                                                                                    \
      return 0;                                                                                                                    \
    }                                                                                                                              \
    *sum = y[0];                                                                                                                   \
    for (size_t i = 1; i < N; ++i)                                                                                                 \
      *sum *= y[i];                                                                                                                \
    return 0;                                                                                                                      \
  }                                                                                                                                \
  int vec_cum_sum_##type (const size_t N, const type* x, type* y) {                                                                \
    if (N < 1)                                                                                                                     \
      return 0;                                                                                                                    \
    y[0] = x[0];                                                                                                                   \
    for (size_t i = 1; i < N; ++i)                                                                                                 \
      y[i] = y[i - 1] + x[i];                                                                                                      \
    return 0;                                                                                                                      \
  }                                                                                                                                \
  int vec_cum_prod_##type (const size_t N, const type* x, type* y) {                                                               \
    if (N < 1)                                                                                                                     \
      return 0;                                                                                                                    \
    y[0] = x[0];                                                                                                                   \
    for (size_t i = 1; i < N; ++i)                                                                                                 \
      y[i] = y[i - 1] * x[i];                                                                                                      \
    return 0;                                                                                                                      \
  }

TYPED_VEC_FUNCS (int)
TYPED_VEC_FUNCS (size_t)
TYPED_VEC_FUNCS (float)
TYPED_VEC_FUNCS (double)
#undef TYPED_VEC_FUNCS

/* Computes y = alpha x + y */
int vec_AXPY (size_t N, const double alpha, const double* x, double* y) {
  if (alpha == 0)
    return 0;
  for (size_t i = 0; i < N; ++i)
    y[i] = alpha * x[i] + y[i];
  return 0;
}

/* Computes y = alpha x + beta y */
int vec_AXPBY (size_t N, const double alpha, const double* x, const double beta, double* y) {
  if (beta == 1)
    return vec_AXPY (N, alpha, x, y);
  if (alpha == 0)
    return vec_scale (N, beta, y);
  for (size_t i = 0; i < N; ++i)
    y[i] = alpha * x[i] + beta * y[i];
  return 0;
}

/* Computes z = alpha x + beta y + gamma z
     The implementation is optimized for alpha of 1 and gamma of 1 or 0
*/
int vec_AXPBYPCZ (size_t        N,
                  const double  alpha,
                  const double* x,
                  const double  beta,
                  const double* y,
                  const double  gamma,
                  double*       z) {
  if (alpha == 1) {
    if (gamma == 0) {
      for (size_t i = 0; i < N; ++i)
        z[i] = x[i] + beta * y[i];
    } else if (gamma == 1) {
      for (size_t i = 0; i < N; ++i)
        z[i] = x[i] + beta * y[i] + z[i];
    } else {
      for (size_t i = 0; i < N; ++i)
        z[i] = x[i] + beta * y[i] + gamma * z[i];
    }
  } else {
    if (gamma == 0) {
      for (size_t i = 0; i < N; ++i)
        z[i] = alpha * x[i] + beta * y[i];
    } else if (gamma == 1) {
      for (size_t i = 0; i < N; ++i)
        z[i] = alpha * x[i] + beta * y[i] + z[i];
    } else {
      for (size_t i = 0; i < N; ++i)
        z[i] = alpha * x[i] + beta * y[i] + gamma * z[i];
    }
  }
  return 0;
}

/* Compares 2 floats if they are equal within a specified range, if precision==0, use the machine precision */
int approx_rel (const double a, const double b, const double precision) {
  double epsilon = precision == 0 ? DBL_EPSILON : precision;
  return fabs (a - b) <= MAX (a, b) * epsilon;
}
int approx_abs (const double a, const double b, const double precision) {
  double epsilon = precision == 0 ? DBL_EPSILON : precision;
  return fabs (a - b) <= epsilon;
}

/* Searches for index of max value in array
    - use ints to allow for negative search indices
    - returns error 1 if no maximum is found or -1 if ambiguous results
*/
int max_loc (const double* arr, int start, int end, int* ierr) {
  *ierr      = 1;
  int maxloc = start;
  if (arr[start] >= arr[end])
    *ierr = 0;
  for (int i = start + 1; i <= end; ++i) {
    if (arr[i] >= arr[maxloc]) {
      *ierr = 0;
      if (arr[i] == arr[maxloc])
        *ierr = -1;
      maxloc = i;
    }
  }
  return maxloc;
}

/* Inplace reverse array entries */
int reverse_1d_size_t (size_t N, size_t* arr) {
  // dont do stack allocations if temp array is bigger than approx 1MB
  if ((N / 2) * sizeof (size_t) > 1e6)
    CHKERROUT (1, "array too large for this implementation");

  size_t tmp[N / 2];
  size_t i;
  for (i = 0; i < N / 2; i++)
    tmp[i] = arr[i];

  for (i = 0; i < N / 2; i++)
    arr[i] = arr[N - i - 1];

  for (i = 0; i < N / 2; i++)
    arr[N - i - 1] = tmp[i];

  return 0;
}
int reverse_1d_f4 (size_t N, float* arr) {
  // dont do stack allocations if temp array is bigger than approx 1MB
  if ((N / 2) * sizeof (float) > 1e6)
    CHKERROUT (1, "array too large for this implementation");

  float  tmp[N / 2];
  size_t i;
  for (i = 0; i < N / 2; i++)
    tmp[i] = arr[i];

  for (i = 0; i < N / 2; i++)
    arr[i] = arr[N - i - 1];

  for (i = 0; i < N / 2; i++)
    arr[N - i - 1] = tmp[i];

  return 0;
}
int reverse_1d_f8 (size_t N, double* arr) {
  // dont do stack allocations if temp array is bigger than approx 1MB
  if ((N / 2) * sizeof (double) > 1e6)
    CHKERROUT (1, "array too large for this implementation");

  double tmp[N / 2];
  size_t i;
  for (i = 0; i < N / 2; i++)
    tmp[i] = arr[i];

  for (i = 0; i < N / 2; i++)
    arr[i] = arr[N - i - 1];

  for (i = 0; i < N / 2; i++)
    arr[N - i - 1] = tmp[i];

  return 0;
}

// count all true/false values in arr[N]
size_t vec_count_true (const size_t N, const int* arr) {
  size_t cnt = 0;
  for (size_t i = 0; i < N; ++i)
    if (arr[i])
      ++cnt;
  return cnt;
}
size_t vec_count_false (const size_t N, const int* arr) {
  size_t cnt = 0;
  for (size_t i = 0; i < N; ++i)
    if (!arr[i])
      ++cnt;
  return cnt;
}

// Determine distance/edge length between two points with N dimensions
double vec_distance_between_two_points (const size_t N, const double* p1, const double* p2) {
  double distance = 0;
  for (size_t i = 0; i < N; ++i)
    distance += (p2[i] - p1[i]) * (p2[i] - p1[i]);
  return sqrt (distance);
}

// Use Herons Formula to determine the area of a triangle given the 3 edge lengths
double triangle_area_by_edgelengths (const double e1, const double e2, const double e3) {
  const double p = (e1 + e2 + e3) / 2;
  return sqrt (p * (p - e1) * (p - e2) * (p - e3));
}

double triangle_area_by_vertices (const double p1[3], const double p2[3], const double p3[3]) {
  const double e1 = vec_distance_between_two_points (3, p2, p3);
  const double e2 = vec_distance_between_two_points (3, p3, p1);
  const double e3 = vec_distance_between_two_points (3, p1, p2);
  return triangle_area_by_edgelengths (e1, e2, e3);
}

// compute offsets to compute indices into a flattend multidimensional array with a dot product, i.e. the offsets that are used in conjunction with the idx1d routine
// Example Input:
//   Ndim = 3
//   dimsizes = [4, 2, 3]
//
// Output:
//   offsets -> [2*3, 3, 1]
//
// So that idx1d(Ndim, [2, 1, 0], offsets) yields 15
int nd_offsets (const size_t Ndim, const size_t* dimsizes, size_t* offsets) {
  offsets[Ndim - 1] = 1;
  for (size_t idim = Ndim - 2; idim < Ndim; --idim) {
    offsets[idim] = offsets[idim + 1] * dimsizes[idim + 1];
  }
  return 0;
}

size_t idx1d (const size_t Ndim, const size_t* idx_nd, const size_t* offsets) {
  size_t idx = 0;
  for (size_t i = 0; i < Ndim; ++i)
    idx += idx_nd[i] * offsets[i];
  return idx;
}

inline int is_between (const double x, const double a, const double b) {
  return ((x - a) * (x - b)) < 0;
}
inline int is_inrange (const double x, const double a, const double b) {
  return ((x - a) * (x - b)) <= 0;
}

long factorial (const int n) {
  long fact = 1;
  for (int i = 1; i <= n; i++)
    fact *= i;
  return fact;
}

inline double rad2deg (double rad) {
  const double f = 180 / M_PI;
  return rad * f;
}

inline double deg2rad (double deg) {
  const double f = M_PI / 180;
  return deg * f;
}

//https://www.maplesoft.com/support/help/maple/view.aspx?path=MathApps%2FProjectionOfVectorOntoPlane
int vec3_proj_on_plane (const double v[3], const double plane_normal[3], double proj[3]) {
  double dp      = vec3_dot_product (v, plane_normal);
  double pn_norm = vec3_norm2 (plane_normal);
  pn_norm        = pn_norm * pn_norm;
  proj[0]        = v[0] - dp * plane_normal[0] / pn_norm;
  proj[1]        = v[1] - dp * plane_normal[1] / pn_norm;
  proj[2]        = v[2] - dp * plane_normal[2] / pn_norm;
  return 0;
}

double angle_between_two_normed_vec (const double v1[3], const double v2[3]) {
  const double dp = vec3_dot_product (v1, v2);
  return acos (dp);
}

// given cartesian direction [u,v,w] return azimuth phi and zenith theta
// w == 1, i.e. upward -> theta == 180
// vector going toward minus y, i.e. from north to south -> azimuth 180
// vector going toward plus x, i.e. from west to east -> azimuth 90
int cartesian_2_spherical (const double u, const double v, const double w, double* phi, double* theta) {
  const double e_x[3] = {1, 0, 0};
  const double e_y[3] = {0, 1, 0};
  const double e_z[3] = {0, 0, -1};

  double dir[3] = {u, v, w};
  int    status = vec3_normalize (dir);
  CHKERR (status);

  *theta = rad2deg (angle_between_two_normed_vec (dir, e_z));

  double proj[3];
  status = vec3_proj_on_plane (dir, e_z, proj);
  CHKERR (status);

  vec3_normalize (proj);  // delibaretly ignore error in case of sza == 0

  double azi = angle_between_two_normed_vec (proj, e_y);
  azi *= sign (vec3_dot_product (proj, e_x));

  *phi = rad2deg (azi);

  return 0;
}
