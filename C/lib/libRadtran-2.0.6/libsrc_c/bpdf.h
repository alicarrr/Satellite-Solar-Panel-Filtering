/************************************************************************
 * $Id$
 *
 * MYSTIC - Monte Carlo code for the physically correct tracing of
 *          photons in cloudy atmospheres.
 *
 * Copyright (c) 2000-2012 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras
 *
 * Correspondence: bernhard.mayer@lmu.de
 *
 ************************************************************************/

/** @file 
 * The files `bpdf.h`, `bpdf.c` and parts of `mystic.c` contain code related to
 * land surface bpdfs.
 * @author Cedric Wehrum
 * @date April 1 2023
 */

#ifndef __bpdf_h
#define __bpdf_h

#if defined(__cplusplus)
extern "C" {
#endif

#include <complex.h>
#include <stdlib.h>

#define BPDF_NONE 0
#define BPDF_TSANG 1
#define BPDF_LITVINOV 2
#define BPDF_MAIGN 3

#define LITVINOV_CONSTANT 0
#define LITVINOV_FROM_FILE 1

#define MAIGN_CONSTANT 2
#define MAIGN_FROM_FILE 3

typedef struct {
  float u10; /* wind speed for bpdf_tsang */
} bpdf_ocean_struct;

typedef struct {
  /* parameters for Maignan Breon BPDF*/
  /* single empirical parameter */
  float c_maign;
  /* ARVI or NDVI */
  /* defaults to 0.0 */
  float arvi;
  /* real and imaginary part of the refractive index */
  /* defaults to 1.5 + 0i */
  float refr_index;
  float refr_index_i;
} bpdf_maign_struct;

typedef struct {
  /* parameters for Litvinov BPDF */
  /* black sky albedo */
  float albedo;
  /* root mean square of slopes in first scale */
  float rms_slope;
  /* portion of surface elements that can provide Fresnel scattering */
  float fresnel;
  /* shadow parameter */
  float shadowing;
  /* real and imaginary part of the refractive index */
  /* defaults to 1.5 + 0i */
  float refr_index;
  float refr_index_i;
  /* root mean square of slopes in second scale */
  float rms_slope2;
} bpdf_land_struct;

/* Input structure for polarized reflectance function*/
typedef struct {
  int               type;
  bpdf_ocean_struct bpdf_ocean;
  bpdf_land_struct  bpdf_land;
  bpdf_maign_struct bpdf_maign;
} bpdf_struct;

//TODO remove debugging operations
// debug Litvinonv and Maignan Breon matrices
// can be removed when Litvinov and Maignan are stable
#define DEBUG_REFL_PROB_MATRIX

#ifndef M_PI
// taken from GCC's math.h
#define M_E 2.7182818284590452354         /* e */
#define M_LOG2E 1.4426950408889634074     /* log_2 e */
#define M_LOG10E 0.43429448190325182765   /* log_10 e */
#define M_LN2 0.69314718055994530942      /* log_e 2 */
#define M_LN10 2.30258509299404568402     /* log_e 10 */
#define M_PI 3.14159265358979323846       /* pi */
#define M_PI_2 1.57079632679489661923     /* pi/2 */
#define M_PI_4 0.78539816339744830962     /* pi/4 */
#define M_1_PI 0.31830988618379067154     /* 1/pi */
#define M_2_PI 0.63661977236758134308     /* 2/pi */
#define M_2_SQRTPI 1.12837916709551257390 /* 2/sqrt(pi) */
#define M_SQRT2 1.41421356237309504880    /* sqrt(2) */
#define M_SQRT1_2 0.70710678118654752440  /* 1/sqrt(2) */
#endif

// 'matIndexByType' provides access to a 2D-array that is stored inside a 1D-array
// by changing 'matIndex' one can switch from row-major to column-major order and vice versa
// this was used to adapt some code that was originally implemented in FORTRAN
//TODO Use 2D C arrays everywhere and assume row-major order wherever relevant
//TODO only MATRIX4 is used in MYSTIC, all others can be removed
typedef unsigned (*MatIndexFunc) (unsigned, unsigned);

enum MatrixType { SCALAR = 1, MATRIX3 = 3, MATRIX4 = 4 };

static inline unsigned matIndex (enum MatrixType matrixType, unsigned i, unsigned j) {
  return i + matrixType * j;
}
static unsigned matIndex1 (unsigned i, unsigned j) {
  return matIndex (SCALAR, i, j);
}
static unsigned matIndex3 (unsigned i, unsigned j) {
  return matIndex (MATRIX3, i, j);
}
static unsigned matIndex4 (unsigned i, unsigned j) {
  return matIndex (MATRIX4, i, j);
}

static const MatIndexFunc matIndexByType[5] = {
  NULL,
  &matIndex1,
  NULL,
  &matIndex3,
  &matIndex4,
};

//TODO get rid of 'matrixType' and only support 4x4 matrices
/*     The implementation of this model was based Litvinov's FORTRAN that he kindly shared with us.
 *
 *     Physical model for land bpdf based on: (Litvinov et all. JQSRT (2012), 113, 2023-2039. doi:10.1016/j.jqsrt.2012.06.027.)
 *     matrixType: calculate SCALAR, MATRIX3 (3x3) or MATRIX4 (4x4)
 *     solar_za - solar zenith angle
 *     view_za - viewing zenith angle
 *     solar_aa  - solar azimuth angle
 *     view_aa - viewing azimuth angle
 *     numParams - number of parameters given for params
 *     params - input parameters
 *     In this version the meaning of parameters is follows:
 *        params[0] is the albedo of a surface element
 *        params[1] is the rms of surface slope (in scale 1)
 *        params[2] is the fraction of surface which can provide Fresnel reflection
 *        params[3] is a parameter for the shadowing function
 *        params[4] is the rms of surface slope (in scale 2)
 *     An example of retrieved values of the parameters for soil (Mongu):
 *        params[0]=0.042 (440 nm)
 *                  =0.072 (490 nm)
 *                  =0.115 (550 nm)
 *                  =0.155 (670 nm)
 *                  =0.351 (865 nm)
 *                  =0.434 (1020 nm)
 *        params[1]=0.203
 *        params[2]=0.635
 *        params[3]=0.17
 *     An example of retrieved values of the parameters for Banizoumbou:
 *        params[0]=0.072 (440 nm)
 *                  =0.126 (490 nm)
 *                  =0.229 (550 nm)
 *                  =0.371(670 nm)
 *                  =0.510 (865 nm)
 *                  =0.590 (1020 nm)
 *        params[1]=0.115
 *        params[2]=0.725
 *        params[3]=0.1
 *     Examples for sigma_2:
 *        params[4]=0.005    for ocean
 *        params[4]=0.9      for land surfaces
 *     rMat: calculated reflection matrix
 */
int reflection_matrix_litvinov (enum MatrixType      matrixType,
                                double               solar_za,
                                double               view_za,
                                double               solar_aa,
                                double               view_aa,
                                const unsigned short numParams,
                                const double         params[5],
                                double _Complex m_ind,
                                double* rMat);

/* An abstraction function that uses MYSTIC's matrix format and the same notation for the angles
 * TODO can be removed and merged with the function above
 */
int reflection_matrix_litvinov_mystic (double mu_esc,
                                       double mu_inc,
                                       double phi_inc,
                                       double phi_esc,
                                       double _Complex refr_index,
                                       double    albedo,
                                       double    rms_slope,
                                       double    fresnel,
                                       double    shadowing,
                                       double    rms_slope2,
                                       int       backward,
                                       double*** refl_mat);

/* An abstraction function that uses MYSTIC's matrix format and the same notation for the angles
 * TODO can be removed and merged with the function below
 */
int reflection_matrix_maignan_breon_mystic (double mu_esc,
                                            double mu_inc,
                                            double phi_inc,
                                            double phi_esc,
                                            double _Complex refr_index,
                                            double    c_maign,
                                            double    arvi,
                                            double    rpv_rho0,
                                            double    rpv_k,
                                            double    rpv_theta,
                                            double    rpv_scale,
                                            double    rpv_sigma,
                                            double    rpv_t1,
                                            double    rpv_t2,
                                            int       backward,
                                            double*** refl_mat);

/*     The implementation of this model was based Litvinov's FORTRAN that he kindly shared with us.
 *
 *     Physical model for land bpdf based on:
 *       Maignan et al.: Polarized reflectances of natural surfaces:
 *                       Spaceborne measurements and analytical modeling,
 *                       Remote Sensing of Environment, 113, 2642–2650, doi:10.1016/j.rse.2009.07.022
 *     matrixType: calculate SCALAR, MATRIX3 (3x3) or MATRIX4 (4x4)
 *     solar_za - solar zenith angle
 *     view_za - viewing zenith angle
 *     phi - relative azimuth angle
 *     m_ind - complex refractive index
 *     c_mgn - single empirical parameter
 *     arvi - ARVI or NDVI
 *       examples:
 *       surface (IGBP class)     c_mgn   ndvi
 *       evergreen needleleaf (1) 4.98    0.8
 *       evergreen broadleaf (2)  7.36    0.8
 *       grassland (10)           6.11    0.7
 *       barren (16)              7.29    0.0
 *     rMat: calculated reflection matrix
 *
 *     maign_breon_xh is based on a slightly different code that is more clear to understand than maign_breon
 *     Both functions yield identical results
 *     maign_breon may be faster though?
 *     TODO remove one of maign_breon or maign_breon_xh
 */

int maign_breon (double solar_za, double view_za, double phi, double _Complex m_ind, double c_mgn, double arvi, double* rMat);

int maign_breon_xh (enum MatrixType matrixType,
                    double          solar_za,
                    double          view_za,
                    double          phi,
                    double _Complex m_ind,
                    const double c_mgn,
                    const double arvi,
                    double*      rMat);

/* The following are helper function that are also based on Litvinov's FORTRAN code
 * TODO these can likely be removed from the header file
 */

// CB for infinite medium
double CB_peak_new (double zenithAngle);

// Shadowing function
double shadowFunc (double mu0, double muv, double sigma2);

// calculates a Fresnel reflection matrix
double* reflFresnelMy (enum MatrixType matrixType,
                       double          mu01,
                       double          muv1,
                       double          phi,
                       double _Complex m,
                       double* rMat_fr,
                       double* alpha1);

/***********************************************************************
* Given the lower and upper limits of integration a and b,
* and given the number of Gauss-Legendre points numGaussLegPoints,
* this routine returns through array abscissas the abscissas and through
* array weights the weights of the Gauss-Legendre quadrature formula.
* eps is the desired accuracy of the abscissas.
* This routine is documented further in:
*   W.H. Press et al. 'Numerical Recipes' Cambridge Univ. Pr. (1987)
*   page 125 ISBN 0-521-30811-9
**********************************************************************/
double*
gaussLegendre (unsigned numDimensions, unsigned numGaussLegPoints, double a, double b, double* abscissas, double* weights  //,
               //struct DebugValue* debugValues,
               //unsigned* debugCount
);

// PDF used for Litvinov's model
double pdfFresnel (double  alpha,
                   double  beta,
                   double  sigma2,
                   double* dmui1,
                   double  phii,
                   double* dmur1,
                   double  phir,
                   double  mu1r,
                   double  mu2i);

#if defined(__cplusplus)
}
#endif

#endif  //__bpdf_h
