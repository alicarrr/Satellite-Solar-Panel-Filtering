/*--------------------------------------------------------------------
 * $Id$
 *
 * This file is part of libRadtran.
 * Copyright (c) 1997-2012 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras
 *
 * ######### Contact info: http://www.libradtran.org #########
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 *--------------------------------------------------------------------*/

#ifndef __common_math_h
#define __common_math_h

#if defined(__cplusplus)
extern "C" {
#endif

#ifndef NOT_A_NUMBER
#ifdef NAN
#define NOT_A_NUMBER NAN
#else
#define NOT_A_NUMBER 0.0 / 0.0
#endif
#endif

#ifndef POS_INFINITY
#ifdef INFINITY
#define POS_INFINITY INFINITY
#else
#define POS_INFINITY 1.0 / 0.0
#endif
#endif

#define MIN(a, b)                                                                                                                  \
  ({                                                                                                                               \
    __typeof__ (a) _a = (a);                                                                                                       \
    __typeof__ (b) _b = (b);                                                                                                       \
    (_a) < (_b) ? (_a) : (_b);                                                                                                     \
  })

#define MAX(a, b)                                                                                                                  \
  ({                                                                                                                               \
    __typeof__ (a) _a = (a);                                                                                                       \
    __typeof__ (b) _b = (b);                                                                                                       \
    (_a) > (_b) ? (_a) : (_b);                                                                                                     \
  })

/* the first argument of mat_mult should really be:
 * const double * const * a
 * but C specification does not allow it :-(
 * see: https://stackoverflow.com/a/29553017
 * if converted to C++, this would be valid.
 */
void mat_v_mult (size_t N, double** a, const double* b, double* c);
void v_mult (const double* a, const double* b, double* c);
void v_sub (const double* a, const double* b, double* c);
void v_add (const double* a, const double* b, double* c);
void vn_mult (const double* a, const double* b, double* c);
void v_cross_product (const double* a, const double* b, double* c);
int  v_cross_product_norm (const double* a, const double* b, double* c);
int  v_mult_mu (const double* a, const double* b, double* c);
void v_neg (double* a);

double constrainAngle (double x);

double cosd (double angle);
double sind (double angle);
double tand (double angle);
double acosd (double rad);
double asind (double rad);
double atand (double rad);

int cosd_sind (double angle, double* sin_val, double* cos_val);

int sign (double x);

void vec_fabs (const size_t N, const double* a, double* c);

int    vec_dot_product (const size_t N, const double* a, const double* b, double* c);
double vec3_dot_product (const double a[3], const double b[3]);

int    vec_norm2 (const size_t N, const double* v, double* norm);
double vec3_norm2 (const double v[3]);
int    vec3_normalize (double x[3]);

int vec3_scale (const double alpha, double y[3]);
int vec_scale (size_t N, const double alpha, double* y);

#define TYPED_VEC_FUNCS(type)                                                                                                      \
  int vec_sum_##type (const size_t N, const type* y, type* sum);                                                                   \
  int vec_prod_sum_##type (const size_t N, const type* y, type* sum);                                                              \
  int vec_cum_sum_##type (const size_t N, const type* x, type* y);                                                                 \
  int vec_cum_prod_##type (const size_t N, const type* x, type* y);

TYPED_VEC_FUNCS (int)
TYPED_VEC_FUNCS (size_t)
TYPED_VEC_FUNCS (float)
TYPED_VEC_FUNCS (double)
#undef TYPED_VEC_FUNCS

int vec_AXPY (size_t N, const double alpha, const double* x, double* y);
int vec_AXPBY (size_t N, const double alpha, const double* x, const double beta, double* y);
int vec_AXPBYPCZ (size_t N, const double alpha, const double* x, const double beta, const double* y, const double gamma, double* z);

int approx_rel (const double a, const double b, const double precision);
int approx_abs (const double a, const double b, const double precision);

int max_loc (const double* arr, int start, int end, int* ierr);

int reverse_1d_size_t (size_t N, size_t* arr);
int reverse_1d_f4 (size_t N, float* arr);
int reverse_1d_f8 (size_t N, double* arr);

size_t vec_count_true (const size_t N, const int* arr);
size_t vec_count_false (const size_t N, const int* arr);

double vec_distance_between_two_points (const size_t N, const double* p1, const double* p2);
double triangle_area_by_edgelengths (const double e1, const double e2, const double e3);
double triangle_area_by_vertices (const double p1[3], const double p2[3], const double p3[3]);

int    nd_offsets (const size_t Ndim, const size_t* dimsizes, size_t* offsets);
size_t idx1d (const size_t Ndim, const size_t* idx_nd, const size_t* offsets);

int is_between (const double x, const double a, const double b);
int is_inrange (const double x, const double a, const double b);

long factorial (const int n);

double rad2deg (double rad);
double deg2rad (double deg);

int vec3_proj_on_plane (const double v[3], const double plane_normal[3], double proj[3]);

double angle_between_two_normed_vec (const double v1[3], const double v2[3]);

int cartesian_2_spherical (const double u, const double v, const double w, double* phi, double* theta);

#if defined(__cplusplus)
}
#endif

#endif
