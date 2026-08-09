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

/* @44c@ */
/* @c44@ */

/********************************************************************************/
/* In order to use the functions provided by the mie library,          @44_10c@ */
/* #include <miecalc.h> in your source code and link with libRadtran_c.a.       */
/*                                                                              */
/* @strong{Example:}                                                            */
/* Example for a source file:                                                   */
/* @example                                                                     */
/*                                                                              */
/*   ...                                                                        */
/*   #include "../src_c/miecalc.h"                                              */
/*   ...                                                                        */
/*                                                                              */
/* @end example                                                                 */
/*                                                                              */
/* Linking of the executable, using the GNU compiler gcc:                       */
/* @example                                                                     */
/*                                                                              */
/*   gcc -o test test.c -lRadtran_c                                             */
/*                                                                              */
/* @end example                                                        @c44_10@ */
/********************************************************************************/

/****************************************************************************/
/* The Mie library provides functions for Mie calculations,       @44_20c@  */
/* interfacing the MIEV0 and BHMIE codes by Warren Wiscombe                 */
/* (ftp://climate.gsfc.nasa.gov/pub/wiscombe) and Bohren and Huffman        */
/* (ftp://astro.princeton.edu/draine/scat/bhmie/). Functions for            */
/* evaluating the phase function and the integrated phase function are      */
/* also provided.                                                 @c44_20@  */
/****************************************************************************/

#include <math.h>
#include <float.h>
#include <string.h>

#include "common_math.h"
#include "errors.h"
#include "f77-uscore.h"
#include "fortran_and_c.h"
#include "miecalc.h"
#include "numeric.h"

#if HAVE_NETCDF4
#include <netcdf.h>
#endif

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif

/* This is a fix needed at MIM, need to test whether this also works on DLR */
#if COMPILING_CONDOR_AT_MIM
#define fscanf __fscanf
#endif

int refwat_iapws (float wavelength, float temperature, float* reff_re, float* reff_im);

/***********************************************************************************/
/* Function: mie_calc                                                     @44_30i@ */
/* Description:                                                                    */
/*  Mie calculations, using the MIEV0 and BHMIE codes by Warren Wiscombe           */
/*  (ftp://climate.gsfc.nasa.gov/pub/wiscombe) and Bohren and Huffman              */
/*  (ftp://astro.princeton.edu/draine/scat/bhmie/).                                */
/*                                                                                 */
/* Parameters:                                                                     */
/*  mie_inp_struct input:    mie input structure (see src_c/miecalc.h)             */
/*  mie_out_struct *output:  mie output structure (see src_c/miecalc.h)            */
/*  int program:             MIEV0 or BHMIE                                        */
/*  int medium:              WATER, ICE, or USER; if USER or AEROSOL,              */
/*                           the refractive index is read from crefin              */
/*  mie_complex crefin:      Complex refractive index (both numbers positive)      */
/*  float wavelength:        Wavelength [micron]                                   */
/*  float radius:            Droplet radius [micron]                               */
/*  float temperature:       Temperature [K]                                       */
/*  mie_complex *ref:        Complex index of refraction                           */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/*  Syntax and parameters of this function are subject to change.                  */
/* Author:                                                                         */
/*                                                                        @i44_30@ */
/***********************************************************************************/
int mie_calc (mie_inp_struct  input,
              mie_out_struct* output,
              int             program,
              int             medium,
              int             water_properties,
              mie_complex     crefin,
              double          wavelength,
              double          radius,
              float           temperature,
              int             nstokes,
              int             numang,
              mie_complex*    ref) {
  void F77_FUNC (wrefice, WREFICE) (float* wl, float* temp, float* re, float* im);

  void F77_FUNC (wrefwat, WREFWAT) (float* wl, float* temp, float* re, float* im);

  void F77_FUNC (bhmie, BHMIE) (float*       xxf,
                                mie_complex* crefin,
                                int*         numang,
                                mie_complex* s1,
                                mie_complex* s2,
                                float*       qext,
                                float*       qsca,
                                float*       qback,
                                float*       gsca);

  void F77_FUNC (miev0, MIEV0) (double*      xxd,
                                mie_complex* crefin,
                                int*         perfct,
                                float*       mimcut,
                                int*         anyang,
                                int*         numang,
                                double*      xmu,
                                int*         nmom,
                                int*         ipolzn,
                                int*         momdim,
                                int*         prnt,
                                float*       qext,
                                float*       qsca,
                                float*       gqsc,
                                float*       pmom,
                                mie_complex* sforw,
                                mie_complex* sback,
                                mie_complex* s1,
                                mie_complex* s2,
                                mie_complex* tforw,
                                mie_complex* tback,
                                float*       spike);

  /* input.nmom+2*4, because miev0 starts counting Legendre polynomials
     from 0 and fortran adds one element at end of string*/
  float*  tmp_pmom = calloc ((input.nmom + 2) * nstokes, sizeof (float));
  float** pmom_array;

  float wavelength_f = (float) wavelength;
  float xxf;

  int status;

  static int first = 1;

  pmom_array = calloc (nstokes, sizeof (float*));
  for (int ip = 0; ip < nstokes; ip++)
    pmom_array[ip] = calloc (input.nmom + 2, sizeof (float));

  if (nstokes > 1 && program != MIEV0) {
    CHKERROUT (-1, "Error: For nstokes > 1 only the program MIEV0 can be used.");
  }

  switch (medium) {
  case ICE:
    switch (water_properties) {
    case WATPROP_IAPWS:
      fprintf (stderr,
               "Warning: water_properties iapws together with refrac ice uses the default values for the refractive index and not "
               "the IAPWS equation which is only valid for liquid water.\n");
    }
    F77_FUNC (wrefice, WREFICE) (&wavelength_f, &temperature, &(ref->re), &(ref->im));
    if (ref->re == 0.0 && ref->im == 0.0) {
      CHKERROUT (-1, "encountered invalid re/im values from wrefice");
    }
    break;

  case WATER:
    switch (water_properties) {
    case WATPROP_SEGELSTEIN:
      F77_FUNC (wrefwat, WREFWAT) (&wavelength_f, &temperature, &(ref->re), &(ref->im));
      if (ref->re == 0.0 && ref->im == 0.0) {
        CHKERR (-1);
      }
      break;
    case WATPROP_IAPWS:
      status = refwat_iapws (wavelength_f, temperature, &(ref->re), &(ref->im));
      CHKERR (status);
      break;
    default:
      fprintf (stderr, "Error, unknown water_properties %d\n", water_properties);
      CHKERR (-1);
    }
    break;

  case USER:
  case AEROSOL:
    ref->re = crefin.re;
    ref->im = crefin.im;
    break;

  default:
    CHKERROUT (-1, "unknown medium");
  }
  output->crefin.re = ref->re;
  output->crefin.im = ref->im;
  output->xxd       = 2.0 * PI * radius / wavelength;
  xxf               = (float) output->xxd;

  switch (program) {
  case BHMIE:;

    int numang_bhmie = (numang + 1) / 2;
    F77_FUNC (bhmie, BHMIE)
    (&xxf, &output->crefin, &numang_bhmie, output->s1, output->s2, &output->qext, &output->qsca, &output->qback, &output->gsca);

    if (numang > 1 && first) {
      first       = 0;
      FILE*  f    = fopen ("bhmie.amplitudes", "w");
      double dang = PI / (double) (numang - 1);
      for (int iang = 0; iang < numang; iang++) {
        double theta      = iang * dang;
        output->xmu[iang] = cos (theta);
        fprintf (f,
                 "%11.6f %14.6e %14.6e  %14.6e %14.6e\n",
                 theta * 180.0 / PI,
                 output->s1[iang].re,
                 output->s1[iang].im,
                 output->s2[iang].re,
                 output->s2[iang].im);
      }
      fclose (f);
    }
    break;

  case MIEV0:

    /* Loop over number of Stokes components
       description of *ipolzn* in Mie documentation. */
    if (nstokes == 1)
      input.ipolzn = 0;
    else if (nstokes == 2)
      input.ipolzn = 12;
    else if (nstokes == 3)
      input.ipolzn = 123;
    else if (nstokes == 4)
      input.ipolzn = 1234;
    else {
      CHKERROUT (-1, "Number of Stokes components (nstokes) not between 1 and 4");
    }

    /* This is needed to avoid segmentation fault in case miev0 is used
       without calculating Legendre polynomials (miev0 gives segmentation
       fault if input.momdim=0.)*/
    if (input.nmom == 0)
      input.momdim = 1;

    F77_FUNC (miev0, MIEV0)
    (&output->xxd,
     &output->crefin,
     &input.perfct,
     &input.mimcut,
     &input.anyang,
     &numang,
     output->xmu,
     &input.nmom,
     &input.ipolzn,
     &input.momdim,
     input.prnt,
     &output->qext,
     &output->qsca,
     &output->gqsc,
     tmp_pmom,
     &output->sforw,
     &output->sback,
     output->s1,
     output->s2,
     output->tforw,
     output->tback,
     &output->spike);

    output->gsca = output->gqsc / output->qsca;

    if (input.nmom > 0) {
      fortran2c_2D_float_ary_noalloc (nstokes, input.nmom + 2, tmp_pmom, pmom_array);

      for (int im = 0; im <= input.nmom; im++)
        for (int ip = 0; ip < nstokes; ip++)
          output->pmom[ip][im] = pmom_array[ip][im];
    }

    break;
  default:
    CHKERROUT (-1, "unknown Mie program");
  }

  free (tmp_pmom);

  for (int ip = 0; ip < nstokes; ip++)
    free (pmom_array[ip]);
  free (pmom_array);

  return 0;
}

/***********************************************************************************/
/* Function: mie_calc_sizedist                                            @44_30i@ */
/* Description:                                                                    */
/*  Mie calculations, using the MIEV0 and BHMIE codes by Warren Wiscombe           */
/*  (ftp://climate.gsfc.nasa.gov/pub/wiscombe) and Bohren and Huffman              */
/*  (ftp://astro.princeton.edu/draine/scat/bhmie/).                                */
/*                                                                                 */
/* Parameters:                                                                     */
/*  mie_inp_struct input:    mie input structure (see src_c/miecalc.h)             */
/*  mie_out_struct *output:  mie output structure (see src_c/miecalc.h)            */
/*  int program:             MIEV0 or BHMIE                                        */
/*  int medium:              WATER, ICE, or USER; if USER, the refractive index    */
/*                           is read from crefin                                   */
/*  mie_complex crefin:      Complex refractive index (both numbers positive)      */
/*  float wavelength:        wavelength [micron]                                   */
/*  float temperature:       temperature                                           */
/*  double *x_size:          size distribution, radius [um]                        */
/*  double *y_size:          size distribution, n(r)                               */
/*  int n_size:              size distribution, number of radii                    */
/*  double *beta:            extinction coefficient [km-1] per unit                */
/*                           liquid water content (returned)                       */
/*  double *omega:           Single scattering albedo (returned)                   */
/*  double *g:               Asymmetry parameter (returned)                        */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/*  Syntax and parameters of this function are subject to change.                  */
/* Author:                                                                         */
/*                                                                        @i44_30@ */
/***********************************************************************************/

int mie_calc_sizedist (mie_inp_struct  input,
                       mie_out_struct* output,
                       int             program,
                       int             medium,
                       int             water_properties,
                       mie_complex     crefin,
                       double          wavelength,
                       float           temperature,
                       int             nstokes,
                       double*         x_size,
                       double*         y_size,
                       int             n_size,
                       double*         beta,
                       double*         omega,
                       double*         g,
                       mie_complex*    ref) {
  int status = 0, is = 0, ip = 0, im = 0;

  double xsquared = 0, xcubed = 0, norm = 0;

  double ext = 0, sca = 0, asy = 0, vol = 0;

  /* allocate memory for fields */
  double* extinction = calloc (n_size, sizeof (double));
  double* scattering = calloc (n_size, sizeof (double));
  double* asymmetry  = calloc (n_size, sizeof (double));
  double* volume     = calloc (n_size, sizeof (double));

  double*** pmom = NULL;

  if (input.nmom > 0) {
    pmom = calloc (4, sizeof (double**));

    for (ip = 0; ip < 4; ip++) {
      pmom[ip] = calloc (input.nmom + 1, sizeof (double*));

      for (im = 0; im <= input.nmom; im++)
        pmom[ip][im] = calloc (n_size, sizeof (double));
    }
  }

  // allocating one entry for xmu because TESTMI in miev0.f expects that
  if (input.numang < 1)
    output->xmu = calloc (1, sizeof (double));

  for (is = 0; is < n_size; is++) {

    status = mie_calc (input,
                       output,
                       program,
                       medium,
                       water_properties,
                       crefin,
                       wavelength,
                       x_size[is],
                       temperature,
                       nstokes,
                       input.numang,
                       ref);
    CHKERR (status);

    xsquared = x_size[is] * x_size[is] * 1e-12;
    xcubed   = xsquared * x_size[is] * 1e-6;

    extinction[is] = output->qext * xsquared * y_size[is];
    scattering[is] = output->qsca * xsquared * y_size[is];
    asymmetry[is]  = output->gsca * scattering[is];
    volume[is]     = xcubed * y_size[is];

    /* quote from MIEV.doc, concerning the normalization of pmom[]:   */
    /*	                                                              */
    /*	The normalized moments are                                    */
    /*      4 / ( XX**2 * QSCA ) * PMOM,  but it is PMOM itself, not  */
    /*      these normalized moments, which should be integrated over */
    /*      a size distribution.                                      */

    if (input.nmom > 0)
      for (ip = 0; ip < 4; ip++)
        for (im = 0; im <= input.nmom; im++)
          pmom[ip][im][is] = output->pmom[ip][im] * y_size[is];
  }

  norm = integrate (x_size, y_size, n_size);

  ext = integrate (x_size, extinction, n_size) / norm;
  sca = integrate (x_size, scattering, n_size) / norm;
  asy = integrate (x_size, asymmetry, n_size) / norm;
  vol = integrate (x_size, volume, n_size) / norm;

  if (input.nmom > 0)
    for (ip = 0; ip < 4; ip++)
      for (im = 0; im <= input.nmom; im++)
        output->pmom[ip][im] = integrate (x_size, pmom[ip][im], n_size) / norm;

  *beta  = 3.0 * ext / (4.0 * vol) / 1000.0;
  *omega = sca / ext;
  *g     = asy / sca;

  /* free memory */
  free (extinction);
  free (scattering);
  free (asymmetry);
  free (volume);

  if (input.nmom > 0) {
    for (ip = 0; ip < 4; ip++) {
      for (im = 0; im <= input.nmom; im++)
        free (pmom[ip][im]);

      free (pmom[ip]);
    }

    free (pmom);
  }

  // freeing xmu again
  if (input.numang < 1)
    free (output->xmu);

  return 0;
}

/************************************************************************************/
/* Function: phase_function                                               @44_30i@  */
/* Description:                                                                     */
/*  Calculate the phase function from its moments and output to a file "phase.dat". */
/*  The phase function                                                              */
/*    @iftex                                                                        */
/*    @tex                                                                          */
/*     $p(\mu)$ is @*                                                               */
/*      $$ p (\mu) = \sum_{m=0}^{\infty} (2m+1) \cdot k_m \cdot P_m (\mu) $$        */
/*      where $k_m$ is the m'th moment and $P_m (\mu)$ is the m'th                  */
/*      Legendre polynomial.                                                        */
/*     @end tex                                                                     */
/*    @end iftex                                                                    */
/*    @ifinfo                                                                       */
/*     p(mu) is @*                                                                  */
/*     p (mu) = sum (m=0 to infinity) (2m+1) * k(m) * Pm (mu) @*                    */
/*     where k(m) is the m'th moment and Pm (mu) is the m'th Legendre polynomial.   */
/*    @end ifinfo                                                                   */
/*                                                                                  */
/* Parameters:                                                                      */
/*  float *moment: moments of the phase function, f[0, ..., L-1]                    */
/*  int L:         number of moments                                                */
/*                                                                                  */
/* Return value:                                                                    */
/*  0  if o.k., <0 if error                                                         */
/*                                                                                  */
/* Example:                                                                         */
/* Files:                                                                           */
/* Known bugs:                                                                      */
/*  Syntax and parameters of this function are subject to change.                   */
/* Author:                                                                          */
/*                                                                        @i44_30@  */
/************************************************************************************/

int phase_function (float* moment, int L) {
  FILE* outfile = NULL;

  outfile = fopen ("phase.dat", "w");
  CHKPOINTEROUT (outfile, "failed opening phase.dat for writing");

  const int    Nt  = 100001;
  const double dmu = 2.0 / (double) (Nt - 1);

  const double factor = 1.0 / 4.0 / 3.141692653590 / moment[0];

  for (int it = 0; it < Nt; it++) {

    const double mu = -1.0 + (double) it * dmu;

    double p0 = 1;
    double p1 = mu;

    /* CE: The following can be merged with mom2phase */

    double       sum  = moment[0] * p0 + 3.0 * moment[1] * p1;
    const double sum2 = 0;

    for (int l = 2; l <= L; l++) {
      const double p2 = (double) (2 * l - 1) / (double) l * mu * p1 - (double) (l - 1) / (double) l * p0;

      sum += (double) (2 * l + 1) * moment[l] * p2;
      p0 = p1;
      p1 = p2;
    }

    fprintf (outfile, "%g  %g  %g\n", mu, sum * factor, sum2 * factor);
  }

  fclose (outfile);
  return 0;
}

/**
 * mom2phase
 *
 * calculate phase function from phase function moments
 *
 * @param x cos(theta)
 * @param f Legendre moment vector
 * @param L length of f
 *
 * @return phase function at x
 *
 * @author Bernhard Mayer
 * @date   2009-06-29 Moved to miecalc.c from phase.c by Claudia Emde
 */
double mom2phase (double x, double* f, int L) {

  double plm2 = 1.0;
  double plm1 = x;

  double sum = plm2 * f[0] + plm1 * f[1];

  for (int l = 2; l < L; l++) {
    const double pl = ((double) (2 * l - 1) * x * plm1 - (double) (l - 1) * plm2) / (double) l;

    sum += f[l] * pl;

    plm2 = plm1;
    plm1 = pl;
  }

  return sum;
}

/***********************************************************************************/
/* Function: read_mie_table                                               @44_30i@ */
/* Description:                                                                    */
/*  Read Frank Evans' Mie table, similar to his subroutine READ_MIE_TABLE          */
/*  provided in plotmietab.f; the Mie table can be either a single-column file     */
/*  with the moments of the phase function, or a file complying with Frank Evans'  */
/*  'plotmietab' routine, containing moments for different droplet sizes.          */
/*  IMPORTANT: First, read_mie_table() looks for a file 'filename'.cdf which,      */
/*  when available, is interpreted as the netCDF version of the file.              */
/*  For the netCDF format of the Mie table file, have a look at                    */
/*  src/cloudprp2cdf.sh which converts the ASCII to the netCDF version.            */
/*  ATTENTION: Frank Evans' stores Legendre coefficients f(l), not moments p(l):   */
/*             f = p * (2*l + 1)                                                   */
/*                                                                                 */
/* Parameters:                                                                     */
/*  char *filename:   filename where data is stored                                */
/*  double *r0:       smallest effective radius found in filename                  */
/*  double *dr:       effective radius step                                        */
/*  int *nreff:       number of effective radii                                    */
/*  double *reff:     effective radius grid                                        */
/*  double *wavelen:  wavelength [nm]                                              */
/*  double *nre:      real part of the refractive index                            */
/*  double *nim:      imaginary part of the refractive index                       */
/*  double **extinc:  array[0 ... nreff-1] of extinction coefficients              */
/*  double **albedo:  array[0 ... nreff-1] of single scattering albedos            */
/*  double **f:       array[0 ... nreff-1] of delta scaling factors                */
/*  int **nleg:       array[0 ... nreff-1] storing the number of Legendre coeffc.  */
/*  float ***legen:   array[0 ... nreff-1][0...ip][0 ... nleg]                     */
/*                                                     of Legendre coefficients    */
/*  size_t nphamat:   number of phase matrix elements                              */
/*  int ***ntheta:    array[0 ... nreff-1][0...ip] storing the number of           */
/*                                                                    theta values */
/*  float ****theta:  array[0 ... nreff-1][0...ip][0 ... ntheta] of theta values   */
/*  double ****mu:    array[0 ... nreff-1][0...ip][0 ... ntheta] of mu values      */
/*  float ****phase:  array[0 ... nreff-1][0...ip][0 ... ntheta] of phase function */
/*                                                                          values */
/*  int *alloc_moments: return logical whether moments have been read/allocated    */
/*  int *alloc_explicit: return logical whether phase function have been read/all. */
/*  int aerosol:      aerosol flag                                                 */
/*  int nstokes:      number of Stokes components in RTE calculation               */
/*  int nstrmax:                                                                   */
/*  int read_scattering_phase_function:                                            */
/*  int iws:                                                                       */
/*  int niw:                                                                       */
/*  int quiet:        'Shut up!' flag                                              */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i44_30@ */
/***********************************************************************************/

int read_mie_table (char*      filename,
                    double*    r0,
                    double*    dr,
                    size_t*    nreff,
                    double**   reff,
                    double*    wavelen,
                    double*    nre,
                    double*    nim,
                    double**   extinc,
                    double**   gg,
                    double**   albedo,
                    double**   f,
                    int**      nleg,
                    float****  legen,
                    size_t*    nphamat,
                    int***     ntheta,
                    float****  theta,
                    double**** mu,
                    float****  phase,
                    int*       alloc_moments,
                    int*       alloc_explicit,
                    int        aerosol,
                    int        nstokes,
                    int        nstrmax,
                    int        read_scattering_phase_function,
                    int        iws,
                    int        niw,
                    int        quiet) {
  int i = 0, ip = 0, im = 0, status = 0;
  int rows = 0, min_columns = 0, max_columns = 0, max_length = 0, temp = 0;

  double reffmin = 0, reffmax = 0;

  char dummy[1024];

  FILE* file = NULL;

#if HAVE_NETCDF4
  /* try to open the CDF file */

  char cdf_filename[FILENAME_MAX] = "";
  int  ncid = 0, iw = 0;

  int idd_nreff = 0, idd_nmommax = 0, idd_nphamat = 0, idd_nthetamax = 0;
  int id_reff = 0, id_ext = 0, id_ssa = 0, id_f = 0, id_nmom = 0, id_pmom = 0;
  int id_ntheta = 0, id_theta = 0, id_phase = 0;
  int id_wavelen    = 0;
  int new_format    = 2;
  int delta_scaling = 1;
  int found_cdf     = 0;

  size_t nolam = 0;

  size_t nmommax = 0, nthetamax = 0;

  size_t start[5] = {0, 0, 0, 0, 0};
  size_t count[5] = {0, 0, 0, 0, 0};

  long   version;
  int    complete_moments = 0;
  double Ftot             = 0.0;

  strcpy (cdf_filename, filename);
  strcat (cdf_filename, ".cdf");

  /* open netcdf file */
  status = nc_open (cdf_filename, NC_NOWRITE, &ncid);
  if (status == NC_NOERR)
    found_cdf = 1;
  else {
    strcpy (cdf_filename, filename);
    status = nc_open (cdf_filename, NC_NOWRITE, &ncid);
    if (status == NC_NOERR)
      found_cdf = 1;
  }

  if (found_cdf) {
    /* read file */

    /* read version number. If not defined, old version with only one wavelength per file */
    status = nc_get_att_long (ncid, NC_GLOBAL, "version", &version);
    if (status != NC_NOERR) {
      nolam = 1;
      if (!quiet)
        fprintf (stderr, "No version number existent. Assuming old file version\n");
      version = 0;
    } else
      nolam = 0;

    /* CE new version with only one theta grid for all phamat elements, not yet finished */
    /*  if(version==20110916) */
    /*        new_format=3; */

    /*read complete_moments, which is "1" if all moments to generate */
    /*the phase matrix are stored and "0" if not"  */
    /*   if(new_format==3 ) */
    /*       status = nc_get_att_int (ncid, NC_GLOBAL, "complete_moments", &complete_moments); */

    /**********************/
    /* 1. read dimensions */
    /**********************/

    /***********/
    /* nmommax */
    /***********/

    /* get dimension id for "nmommax" */
    status = nc_inq_dimid_err (ncid, "nmommax", &idd_nmommax, cdf_filename, 1);
    if (status) /* old name was "maxmom" */
      status = nc_inq_dimid_err (ncid, "maxmom", &idd_nmommax, cdf_filename, 0);
    CHKERR (status);

    /* get dimension length for "nmommax" */
    status = nc_inq_dimlen_err (ncid, "nmommax", idd_nmommax, cdf_filename, &nmommax);
    CHKERR (status);

    /*********/
    /* nreff */
    /*********/

    /* get dimension id for "nreff" */
    status = nc_inq_dimid_err (ncid, "nreff", &idd_nreff, cdf_filename, 1);
    if (status) /* old name was "reff" */
      status = nc_inq_dimid_err (ncid, "reff", &idd_nreff, cdf_filename, 1);
    if (status) /* In case of aerosol, the effective radius dimension is replaced by relative humidity */
      status = nc_inq_dimid_err (ncid, "nhum", &idd_nreff, cdf_filename, 0);
    if (status) /* old name was hum */
      status = nc_inq_dimid_err (ncid, "hum", &idd_nreff, cdf_filename, 0);

    CHKERR (status);

    /* get dimension length for "nreff" */
    status = nc_inq_dimlen_err (ncid, "nreff/hum", idd_nreff, cdf_filename, nreff);
    CHKERR (status);

    for (iw = 1; iw < niw; iw++)
      *(nreff + iw) = *nreff;

    /***********/
    /* nphamat */
    /***********/

    /* get dimension id for "nphamat" == "nphase" */
    status = nc_inq_dimid_err (ncid, "nphamat", &idd_nphamat, cdf_filename, 1);
    if (status) /* old name was "nphase" */
      status = nc_inq_dimid_err (ncid, "nphase", &idd_nphamat, cdf_filename, 1);
    if (status) {
      if (!quiet) { /* XXX && nstokes > 1){ */
        fprintf (stderr, "Warning: nphamat not found in netcdf file %s. Assuming old \n", cdf_filename);
        fprintf (stderr, "         dataformat, where only phase function moments are included. \n");
      }
      *nphamat   = 1;
      new_format = 0;
    } else {
      /* get dimension length for "nphamat" */
      status = nc_inq_dimlen_err (ncid, "nphamat", idd_nphamat, cdf_filename, nphamat);
      CHKERR (status);
    }

    if (nstokes > 1 && *nphamat == 1) {
      CHKERROUT (-1,
                 "Error in function read_mie_table: For polarized radiative transfer calculations the mie table has to include the "
                 "phase matrix, your file includes only the phase function!");
    }

    /*************/
    /* nthetamax */
    /*************/

    /* get dimension id for "nthetamax" */
    status = nc_inq_dimid_err (ncid, "nthetamax", &idd_nthetamax, cdf_filename, 1);
    if (status) {
      if (new_format) {
        if (!quiet) {
          fprintf (stderr, "Warning: nthetamax not found in netcdf file %s. Assuming old \n", cdf_filename);
          fprintf (stderr, "         dataformat, where only phase function moments are included. \n");
        }
        nthetamax  = 0;
        new_format = 1;
      }
    } else {
      if (!new_format) {
        fprintf (stderr, "Error %d reading nthetamax from %s\n", status, cdf_filename);
        fprintf (stderr, "Error in function read_mie_table: nthetamax is defined, but nphamat not!\n");
        CHKERR (status);
      }
    }

    if (new_format == 2) {
      /* get dimension length for "nthetamax" */
      status = nc_inq_dimlen_err (ncid, "nthetamax", idd_nthetamax, cdf_filename, &nthetamax);
      CHKERR (status);
    }

    /******************************************************/
    /* 2. check what needs to be read depending on solver */
    /******************************************************/

    /* set whether phases should be read and how many streams should be read */
    switch (read_scattering_phase_function) {
    case READ_SPF_MOMENTS:
      /* this case makes sense for simple nstream-methods */
      /* do not read phases, read nstrmax moments */
      break;
    case READ_SPF_BOTH:
      /* this case makes sense for DISORT2 */
      if (new_format != 2) {
        /* no phases defined, error message */
        CHKERROUT (-1, "no phase function defined, you need to use the option 'disort_icm moments'");
      }
      /* else read phases, read nstrmax moments */
      read_scattering_phase_function = READ_SPF_PHASES;
      break;
    case READ_SPF_PHASES:
      /* this case makes sense for MYSTIC */
      if (new_format != 2) {
        /* no phases defined, no not read them, read all moments */
        read_scattering_phase_function = READ_SPF_MOMENTS;
        nstrmax                        = -1;
      } else {
        /* else read phases only */
        read_scattering_phase_function = READ_SPF_PHASES;
        nstrmax                        = 0;
      }
      break;
    case READ_SPF_ALL_MOMENTS:
      /* this case makes sense for DISORT_ICM MOMENTS */
      if (new_format == 2) {
        if (complete_moments == 0) {
          CHKERROUT (-1,
                     "you need all moments of phase functions if you want to use original intensity correction in fdisort/cdisort! "
                     "... else use the option 'disort_icm phase'!");
        }
      }
      /* do not read phases, read all moments */
      read_scattering_phase_function = READ_SPF_MOMENTS;
      nstrmax                        = -1;
      break;
    default:
      fprintf (stderr, "Error, read_scattering_phase_function = %d no allowed!\n", read_scattering_phase_function);
      CHKERR (-1);
    }

    /***************************/
    /* 3. find id's for fields */
    /***************************/

    /* get variable id for "wavelen" */
    status = nc_inq_varid_err (ncid, "wavelen", &id_wavelen, cdf_filename, 0);
    CHKERR (status);

    /* get variable id for "reff" */
    status = nc_inq_varid_err (ncid, "reff", &id_reff, cdf_filename, 1);
    if (status) /* In case of aerosol, the effective radius dimension is replaced by relative humidity */
      status = nc_inq_varid_err (ncid, "hum", &id_reff, cdf_filename, 0);
    CHKERR (status);

    /* get variable id for "ext" */
    status = nc_inq_varid_err (ncid, "ext", &id_ext, cdf_filename, 0);
    CHKERR (status);

    /* get variable id for "ssa" */
    status = nc_inq_varid_err (ncid, "ssa", &id_ssa, cdf_filename, 0);
    CHKERR (status);

    /* get variable id for "f" (delta scaling factor) */
    status = nc_inq_varid (ncid, "sf", &id_f);
    if (status != NC_NOERR)
      /* Do nothing, because delta scaling factor is only included
         in "Hu-fitted". Set delta_scaling flag to 0. */
      delta_scaling = 0;

    /* get variable id for "nmom" */
    status = nc_inq_varid_err (ncid, "nmom", &id_nmom, cdf_filename, 0);
    CHKERR (status);

    /* get variable id for "pmom" */
    if (nstrmax) {
      status = nc_inq_varid_err (ncid, "pmom", &id_pmom, cdf_filename, 0);
      CHKERR (status);
    }

    if (read_scattering_phase_function == READ_SPF_PHASES) {
      /* get variable id for "ntheta" */
      status = nc_inq_varid_err (ncid, "ntheta", &id_ntheta, cdf_filename, 0);
      CHKERR (status);

      /* get variable id for "theta" */
      status = nc_inq_varid_err (ncid, "theta", &id_theta, cdf_filename, 0);
      CHKERR (status);

      /* get variable id for "phase" */
      status = nc_inq_varid_err (ncid, "phase", &id_phase, cdf_filename, 0);
      CHKERR (status);
    }

    /* (CE) The refractive index is not used in uvspec. Commented this
       because for aerosols the refractive index varies with
       relative humidity. */
    /* get variable id for "nre" */
    /*
      status = nc_inq_varid_err (ncid, "nre", &id_nre, cdf_filename, 0);
      if (status) return status;
    */

    /* get variable id for "nim" */
    /*
      status = nc_inq_varid_err (ncid, "nim", &id_nim, cdf_filename, 0);
      if (status) return status;
    */

    /******************/
    /* 4. read fields */
    /******************/

    for (iw = 0; iw < niw; iw++) {

      /* allocate memory for fields */
      *(reff + iw)   = calloc (*nreff, sizeof (double));
      *(extinc + iw) = calloc (*nreff, sizeof (double));
      *(albedo + iw) = calloc (*nreff, sizeof (double));
      *(gg + iw)     = calloc (*nreff, sizeof (double));
      *(f + iw)      = calloc (*nreff, sizeof (double));
      *(nleg + iw)   = calloc (*nreff, sizeof (int));

      /* set array index for wavelength. Is ignored if nolam=1. */
      start[0] = iw + iws;
      count[0] = 1;

      /* read "wavelen" */
      status = nc_get_vara_double_err (ncid, "wavelen", id_wavelen, start, count, cdf_filename, wavelen + iw);
      CHKERR (status);

      /* convert wavelength from micron to nm */
      *(wavelen + iw) *= 1000.0;

      /* read "nre" */
      /*
	status = nc_get_vara_double_err (ncid, "nre", id_nre, start+nolam, count+nolam, cdf_filename, nre);
	if (status) return status;
      */

      /* read "nim" */
      /*
	status = nc_get_vara_double_err (ncid, "nim", id_nim, start+nolam, count+nolam, cdf_filename, nim);
	if (status) return status;
      */

      /* set array index for reff */
      start[1] = 0;
      count[1] = *nreff;

      /* read "reff"; never depends on wavelen, hence "start+1" */
      status = nc_get_vara_double_err (ncid, "reff", id_reff, start + 1, count + 1, cdf_filename, *(reff + iw));
      CHKERR (status);

      /* read "ext" */
      status = nc_get_vara_double_err (ncid, "ext", id_ext, start + nolam, count + nolam, cdf_filename, *(extinc + iw));
      CHKERR (status);

      /* read "ssa" */
      status = nc_get_vara_double_err (ncid, "ssa", id_ssa, start + nolam, count + nolam, cdf_filename, *(albedo + iw));
      CHKERR (status);

      if (delta_scaling == 1) {
        /* read "sf" (delta_scaling factor)*/
        status = nc_get_vara_double_err (ncid, "sf", id_f, start + nolam, count + nolam, cdf_filename, *(f + iw));
        CHKERR (status);
      }

      /* BEGIN: Paul Ockenfuss: read asymmetry parameter g from mie file for use in veiled core */
      /* read "nmom" */
      status = nc_get_vara_int_err (ncid, "nmom", id_nmom, start + nolam, count + nolam, cdf_filename, *(nleg + iw));
      CHKERR (status);

      if (*(nleg[iw]) > 1) {
        start[1] = 0;
        start[2] = 0;
        start[3] = 0;
        count[1] = *nreff;
        count[2] = 1;
        count[3] = 1;

        if (new_format)
          start[3] = 1;
        else
          start[2] = 1;

        /* get variable id for "pmom" */
        status = nc_inq_varid_err (ncid, "pmom", &id_pmom, cdf_filename, 0);
        CHKERR (status);

        status = nc_get_vara_double_err (ncid, "pmom", id_pmom, start + nolam, count + nolam, cdf_filename, *(gg + iw));
        CHKERR (status);

        /* convert 1st moment to asymmetry parameter */
        for (i = 0; i < *nreff; i++) {
          (*(gg + iw))[i] /= 3.0;
        }
      }
      /* END Paul Ockenfuss */

      if (nstrmax != 0) { /* only read if needed (not needed for MC) */

        /* CE, the following two lines can be deleted... */
        start[2] = 0;
        count[2] = 1;

        /* read "nmom" */
        status = nc_get_vara_int_err (ncid, "nmom", id_nmom, start + nolam, count + nolam, cdf_filename, *(nleg + iw));
        CHKERR (status);

        /* test whether moments are sufficient */
        if (nmommax <= nstrmax || nstrmax == -1)
          for (i = 0; i < *nreff; i++)
            if ((*(nleg + iw))[i] >= nmommax) {
              fprintf (
                stderr,
                "Error! You have specified an nstr which is larger than the number of moments (%d) available in the file %s!\n",
                (int) nmommax,
                cdf_filename);
              fprintf (stderr,
                       "  ... Note that the phase functions described by the file would need more legendre moments to be fully "
                       "described.\n");
              fprintf (stderr, "  ...  please use a complete set of legendre moments, or use less streams.\n");
              CHKERR (-1);
            }

        if (nstrmax != -1) {
          /* reduce number of read moments if limited by nstrmax */
          if (nmommax > nstrmax + 1)
            nmommax = nstrmax + 1;
          for (i = 0; i < *nreff; i++)
            if ((*(nleg + iw))[i] > nmommax)
              (*(nleg + iw))[i] = nmommax;
        }

        /* allocate Legendre polynomials */
        status = ASCII_calloc_float_3D (legen + iw, *nreff, *nphamat, nmommax);
        CHKERR (status);

        *alloc_moments = 1;

        start[1] = 0;
        start[2] = 0;
        start[3] = 0;
        count[1] = 1;
        count[2] = 1;
        count[3] = nmommax;

        for (i = 0; i < *nreff; i++) {
          start[1] = i;
          if (new_format)
            count[3] = (*(nleg + iw))[i];
          else
            count[2] = (*(nleg + iw))[i];

          for (ip = 0; ip < *nphamat; ip++) {
            /* read Legendre polynomials */
            start[2] = ip;

            status =
              nc_get_vara_float_err (ncid, "pmom", id_pmom, start + nolam, count + nolam, cdf_filename, (*(legen + iw))[i][ip]);
            CHKERR (status);

            for (im = 0; im < (*(nleg + iw))[i]; im++)
              (*(legen + iw))[i][ip][im] /= (double) (2 * im + 1);

            /*	    fprintf(stderr, "(*legen)[i][ip][im] %d %d %d  %g \n",i, ip, im, (*legen)[i][ip][im]);  */
          }
        }
      }

      if (read_scattering_phase_function == READ_SPF_PHASES) {
        /* allocate ntheta, theta, mu, phase */
        status = ASCII_calloc_int (ntheta + iw, *nreff, *nphamat);
        CHKERR (status);
        status = ASCII_calloc_float_3D (theta + iw, *nreff, *nphamat, nthetamax);
        CHKERR (status);
        status = ASCII_calloc_double_3D (mu + iw, *nreff, *nphamat, nthetamax);
        CHKERR (status);
        status = ASCII_calloc_float_3D (phase + iw, *nreff, *nphamat, nthetamax);
        CHKERR (status);

        *alloc_explicit = 1;

        start[1] = 0;
        start[2] = 0;
        count[1] = 1;
        count[2] = *nphamat;

        for (i = 0; i < *nreff; i++) {
          start[1] = i;
          /* read "ntheta" */
          status = nc_get_vara_int_err (ncid, "ntheta", id_ntheta, start + nolam, count + nolam, cdf_filename, (*(ntheta + iw))[i]);
          CHKERR (status);
        }

        start[1] = 0;
        start[2] = 0;
        start[3] = 0;
        count[1] = 1;
        count[2] = 1;
        count[3] = nthetamax;

        for (i = 0; i < *nreff; i++) {
          start[1] = i;
          for (ip = 0; ip < *nphamat; ip++) {
            start[2] = ip;
            count[3] = (*(ntheta + iw))[i][ip];

            /* read "theta" */
            status =
              nc_get_vara_float_err (ncid, "theta", id_theta, start + nolam, count + nolam, cdf_filename, (*(theta + iw))[i][ip]);
            CHKERR (status);

            if ((*(theta + iw))[i][ip][0] != 180.0) {
              fprintf (stderr,
                       "Error, theta must begin with 180 degrees in file %s. The value is %g.!\n",
                       cdf_filename,
                       (*(theta + iw))[i][ip][0]);
              CHKERR (-1);
            }
            (*(theta + iw))[i][ip][0] = PI;
            (*(mu + iw))[i][ip][0]    = -1.0;

            for (im = 1; im < (*(ntheta + iw))[i][ip]; im++) {
              (*(theta + iw))[i][ip][im] *= PI / 180.0;
              (*(mu + iw))[i][ip][im] = cos ((double) (*(theta + iw))[i][ip][im]);
            }

            /* read "phase" */
            status =
              nc_get_vara_float_err (ncid, "phase", id_phase, start + nolam, count + nolam, cdf_filename, (*(phase + iw))[i][ip]);
            CHKERR (status);
          }
        }
      } /* endif read phases */

    } /* iw loop */

    nc_close (ncid);

    /*********************/
    /* 5. do some checks */
    /*********************/

    /* Check number of phase matrix elements. */
    if (*nphamat != 1 && *nphamat != 6 && *nphamat != 4) {
      fprintf (stderr, "Error, number of phase matrix elements in %s must be either 1 or 4 or 6\n", filename);
      fprintf (stderr, "but the file includes the polynomials for %zu phase matrix elements. \n", *nphamat);
      CHKERR (-1);
    }

    if (!aerosol) {
      /* check if reff[i] is an equidistant grid */

      if (*nreff > 1) {
        *r0 = (*reff)[0];
        *dr = (*reff)[1] - (*reff)[0];
        for (i = 1; i < *nreff; i++)
          if (!approx_rel ((*reff)[i - 1] + *dr, (*reff)[i], FLT_EPSILON)) {
            fprintf (stderr, "Error, effective radii in %s are not equidistant\n", filename);
            fprintf (stderr, "reff[1] - reff[0] = %g, reff[i]-reff[i-1] = %g\n", *dr, (*reff)[i] - (*reff)[i - 1]);
            CHKERR (-1);
          }
      } else {
        *r0 = (*reff)[0];
        *dr = 0.0;
      }
    }

    /* check if phase function is normalized; otherwise break */
    if (nstrmax != 0) { /* only read if needed (not needed for MC) */
      for (iw = 0; iw < niw; iw++) {
        for (i = 0; i < *nreff; i++)
          if (!approx_rel ((*(legen + iw))[i][0][0], 1, FLT_EPSILON)) {
            fprintf (stderr, "Phase function %s not normalized.\n", cdf_filename);
            fprintf (stderr, "%d %d %d %d %d %e \n", nstrmax, iw, niw, i, (int) (*nreff), (*(legen + iw))[i][0][0]);
            CHKERR (-1);
          }
      }
    }

    if (read_scattering_phase_function != READ_SPF_MOMENTS) {
      for (iw = 0; iw < niw; iw++) {
        for (i = 0; i < *nreff; i++) {
          Ftot = 0.0;
          for (im = 1; im < (*(ntheta + iw))[i][0]; im++) {
            Ftot += ((*(mu + iw))[i][0][im] - (*(mu + iw))[i][0][im - 1]) *
                    ((*(phase + iw))[i][0][im] + (*(phase + iw))[i][0][im - 1]) / 2.0;
          }
          if (Ftot != 2.0) {
            /*  fprintf(stderr, "Normalize phase function with factor %g\n", 2.0/Ftot); */
            for (ip = 0; ip < *nphamat; ip++)
              for (im = 1; im < (*(ntheta + iw))[i][ip]; im++) {
                (*(phase + iw))[i][ip][im] *= 2.0 / Ftot;
              }
          }
        }
      }
    }

    return 0;
  } /* endif found_cdf */

  /*******************************************************************/
  /* The following has not been changed during the major revision of */
  /* meteoric properties done by RPB. There is no guarantee that it  */
  /* still works...                                                  */
  /*******************************************************************/

  if (!quiet) {
    fprintf (stderr, "*** Cannot open netCDF file %s, using\n", cdf_filename);
    fprintf (stderr, "*** ASCII file %s instead.\n", filename);
    fprintf (stderr, "*** It is recommended to use the netCDF format because access is much faster.\n");
  }

#endif

  /* check if the moments file is a 1-column file */
  status = ASCII_checkfile (filename, &rows, &min_columns, &max_columns, &max_length);
  CHKERR (status);

  if (min_columns == 1 && max_columns == 1) {

    *wavelen = -999;

    *nre   = -999;
    *nim   = -999;
    *nreff = 1;

    /* allocate memory for fields */
    *extinc = calloc ((size_t) *nreff, sizeof (double));
    *gg     = calloc ((size_t) *nreff, sizeof (double));
    *albedo = calloc ((size_t) *nreff, sizeof (double));
    *f      = calloc ((size_t) *nreff, sizeof (double));
    *nleg   = calloc ((size_t) *nreff, sizeof (int));

    (*nleg)[0] = rows;

    *r0 = -999;
    *dr = 0;

    *legen      = calloc ((size_t) *nreff, sizeof (float**));
    (*legen)[0] = calloc (1, sizeof (float*));
    /* ??? do we really need the following ??? */
    (*legen)[0][0] = calloc ((size_t) (*nleg)[0], sizeof (float));
    *alloc_moments = 1;

    fprintf (stderr, " ... reading one column of data from %s\n", filename);

    status = read_1c_file_float (filename, &((*legen)[0][0]), &temp);
    CHKERR (status);

    /* final check, just to make sure ... */
    if (temp != rows) {
      fprintf (stderr, "FATAL error reading %s:\n", filename);
      fprintf (stderr, "  found %d rows earlier, %d rows now!\n", rows, temp);
      CHKERR (temp - rows);
    }
    if (rows >= 2)
      (*gg)[0] = (*legen)[0][0][1];
    else {
      (*gg)[0] = 0.0;
      fprintf (stderr, "Warning: found less than two moments in %s!\n", filename);
    }

    fprintf (stderr, " ... read %d coefficients from %s\n", rows, filename);
  } else {

    if ((file = fopen (filename, "r")) == NULL) {
      fprintf (stderr, "Error opening %s for reading\n", filename);
      CHKERR (-1);
    }

    /* first line ignored          */
    if (fgets (dummy, 1024, file) == NULL) {
      fprintf (stderr, "Error ignoring first line of %s\n", filename);
      CHKERR (-1);
    };

    if (!fscanf (file, "%lf", wavelen)) { /* 2nd line, wavelength        */
      fprintf (stderr, "Error reading file %s\n", filename);
      CHKERR (-1);
    }
    /* ignore rest of line         */
    if (fgets (dummy, 1024, file) == NULL) {
      fprintf (stderr, "Error ignoring rest of line of %s\n", filename);
      CHKERR (-1);
    };

    if (!fscanf (file, "%lf %lf", nre, nim)) { /* 3rd line, refractive index  */
      fprintf (stderr, "Error reading file %s\n", filename);
      CHKERR (-1);
    }
    /* ignore rest of line         */
    if (fgets (dummy, 1024, file) == NULL) {
      fprintf (stderr, "Error ignoring rest of line of %s\n", filename);
      CHKERR (-1);
    };

    /* 4th line, distflag, ignored */
    if (fgets (dummy, 1024, file) == NULL) {
      fprintf (stderr, "Error ignoring 4th line, distflag, of %s\n", filename);
      CHKERR (-1);
    };

    /* 5th line, alpha, ignored    */
    if (fgets (dummy, 1024, file) == NULL) {
      fprintf (stderr, "Error ignoring 5th line, alpha, of %s\n", filename);
      CHKERR (-1);
    };

    /* 6th line, number of radii, max and min radius */
    if (!fscanf (file, "%zu %lf %lf", nreff, &reffmin, &reffmax)) {
      fprintf (stderr, "Error reading file %s\n", filename);
      CHKERR (-1);
    }

    /* ignore rest of line         */
    if (fgets (dummy, 1024, file) == NULL) {
      fprintf (stderr, "Error ignoring rest of line of %s\n", filename);
      CHKERR (-1);
    };

    /* convert wavelength from micron to nm */
    *wavelen *= 1000.0;

    /* allocate memory for fields */
    *reff   = calloc (*nreff, sizeof (double));
    *extinc = calloc (*nreff, sizeof (double));
    *gg     = calloc (*nreff, sizeof (double));
    *albedo = calloc (*nreff, sizeof (double));
    *f      = calloc (*nreff, sizeof (double));

    *nleg = calloc (*nreff, sizeof (int));

    *legen = calloc ((size_t) *nreff, sizeof (float**));

    for (i = 0; i < *nreff; i++) {
      if (!fscanf (file, "%lf %lf %lf %d %zu", &((*reff)[i]), &((*extinc)[i]), &((*albedo)[i]), &((*nleg)[i]), nphamat)) {
        fprintf (stderr, "Error reading file %s\n", filename);
        CHKERR (-1);
      }

      /* ignore rest of line */
      if (fgets (dummy, 1024, file) == NULL) {
        fprintf (stderr, "Error ignoring rest of line of %s\n", filename);
        CHKERR (-1);
      };

      /* If nphamat is neither 1 of 6, it is probably the old dataformat and only
         the phase function polynomials are read.*/
      if (*nphamat != 1 && *nphamat != 6)
        *nphamat = 1;

      (*nleg)[i]++;

      /* allocate array for legendre polynomials*/
      (*legen)[i] = calloc (*nphamat, sizeof (float*));
      for (ip = 0; ip < *nphamat; ip++) {
        (*legen)[i][ip] = calloc ((*nleg)[i], sizeof (float));
      }
      *alloc_moments = 1;

      for (ip = 0; ip < *nphamat; ip++) {
        for (im = 0; im < (*nleg)[i]; im++) {
          if (!fscanf (file, "%f", &((*legen)[i][ip][im]))) {
            fprintf (stderr, "Error reading file %s\n", filename);
            CHKERR (-1);
          }
          (*legen)[i][ip][im] /= (double) (2 * im + 1);
        }
      }
      (*gg)[i] = (*legen)[i][0][1];
    }

    /* check if reff[i] is an equidistant grid */
    if (*nreff > 1) {
      *r0 = (*reff)[0];
      *dr = (*reff)[1] - (*reff)[0];
      for (i = 1; i < *nreff; i++)
        if (!approx_rel ((*reff)[i - 1] + *dr, (*reff)[i], FLT_EPSILON)) {
          fprintf (stderr, "Error, effective radii in %s are not equidistant\n", filename);
          fprintf (stderr, "reff[1] - reff[0] = %g, reff[i]-reff[i-1] = %g\n", *dr, (*reff)[i] - (*reff)[i - 1]);
          CHKERR (-1);
        }
    } else {
      *r0 = (*reff)[0];
      *dr = 0.0;
    }
  }

  /* check if phase function is normalized; otherwise break */
  for (i = 0; i < *nreff; i++)
    if (!approx_rel ((*legen)[i][0][0], 1, FLT_EPSILON)) {
      fprintf (stderr, "Phase function %s not normalized :: %g", filename, (*legen)[i][0][0]);
      CHKERR (-1);
    }

  return 0;
}

/***********************************************************************************/
/* Function: read_mie_table_lambda                                        @44_30i@ */
/* Description:                                                                    */
/*  Similar to read_mie_table, but reads only wavelength                           */
/*                                                                                 */
/* Parameters:                                                                     */
/*  char *filename:   filename where data is stored                                */
/*  double *wavelen:  wavelength [nm]                                              */
/*  int quiet:        'Shut up!' flag                                              */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i44_30@ */
/***********************************************************************************/

int read_mie_table_lambda (char* filename, size_t* nlam, double** wavelen, int quiet) {
  int status = 0;
  int rows = 0, min_columns = 0, max_columns = 0, max_length = 0;

  char dummy[1024];

  FILE* f = NULL;

#if HAVE_NETCDF4
  int found_cdf = 0;
  /* try to open the CDF file */

  char cdf_filename[FILENAME_MAX] = "";
  int  ncid = 0, iw = 0;

  int idd_nlam = 0, id_wavelen = 0;

  size_t start[1] = {0};
  size_t count[1] = {0};

  strcpy (cdf_filename, filename);
  strcat (cdf_filename, ".cdf");

  /* open netcdf file */
  status = nc_open (cdf_filename, NC_NOWRITE, &ncid);
  if (status == NC_NOERR)
    found_cdf = 1;
  else {
    strcpy (cdf_filename, filename);
    status = nc_open (cdf_filename, NC_NOWRITE, &ncid);
    if (status == NC_NOERR)
      found_cdf = 1;
  }

  if (found_cdf) {

    /**********************/
    /* 1. read dimensions */
    /**********************/

    /********/
    /* nlam */
    /********/

    /* get dimension id for "nlam" */
    status = nc_inq_dimid_err (ncid, "nlam", &idd_nlam, cdf_filename, 1);
    if (status) { /* old version, only one wavelength per file */
      *nlam    = 1;
      count[0] = 1;
    } else {
      /* get dimension length for "nlam" */
      status = nc_inq_dimlen_err (ncid, "nlam", idd_nlam, cdf_filename, nlam);
      if (status)
        return status;
      count[0] = *nlam;
    }

    /***************************/
    /* 3. find id's for fields */
    /***************************/

    /* get variable id for "wavelen" */
    status = nc_inq_varid_err (ncid, "wavelen", &id_wavelen, cdf_filename, 0);
    if (status)
      return status;

    /******************/
    /* 4. read fields */
    /******************/

    *wavelen = calloc (*nlam, sizeof (double));

    status = nc_get_vara_double_err (ncid, "wavelen", id_wavelen, start, count, cdf_filename, *wavelen);
    if (status)
      return status;

    for (iw = 0; iw < *nlam; iw++) /* convert wavelength from micron to nm */
      (*wavelen)[iw] *= 1000.0;

    nc_close (ncid);

    return 0;
  } /* endif found_cdf */

  if (!quiet) {
    fprintf (stderr, "*** Cannot open netCDF file %s, using\n", cdf_filename);
    fprintf (stderr, "*** ASCII file %s instead.\n", filename);
    fprintf (stderr, "*** It is recommended to use the netCDF format because access is much faster.\n");
  }

#endif

  /* check if the moments file is a 1-column file */
  status = ASCII_checkfile (filename, &rows, &min_columns, &max_columns, &max_length);

  if (status != 0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }

  if (min_columns == 1 && max_columns == 1) {
    **wavelen = -999;
    fprintf (stderr, " ... one-column Mie file, cannot determine wavelength from %s\n", filename);
    return -1;
  } else {

    if ((f = fopen (filename, "r")) == NULL) {
      fprintf (stderr, "Error opening %s for reading\n", filename);
      return -1;
    }

    /* first line ignored          */
    if (fgets (dummy, 1024, f) == NULL) {
      fprintf (stderr, "Error ignoring first line of %s\n", filename);
      return -1;
    };

    if (!fscanf (f, "%lf", *wavelen)) { /* 2nd line, wavelength        */
      fprintf (stderr, "Error reading file %s\n", filename);
      return -1;
    }
    /* ignore rest of line         */
    if (fgets (dummy, 1024, f) == NULL) {
      fprintf (stderr, "Error ignoring rest of line of %s\n", filename);
      return -1;
    };

    /* convert wavelength from micron to nm */
    **wavelen *= 1000.0;
  }

  return 0;
}

/***********************************************************************************/
/* Function: nc_inq_varid_err                                             @44_30i@ */
/* Description:                                                                    */
/*  Calls nc_inq_varid_er and outputs error message if wanted                      */
/*                                                                                 */
/* Parameters:                                                                     */
/*                                                                                 */
/* Return value:                                                                   */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i44_30@ */
/***********************************************************************************/

#if HAVE_NETCDF4
int nc_inq_dimid_err (int ncid, char* varstr, int* id_var, char* cdf_filename, int quiet) {
  int status = 0;

  status = nc_inq_dimid (ncid, varstr, id_var);
  if (status != NC_NOERR) {
    if (!(quiet && status == NC_EBADDIM))
      fprintf (stderr, "Error %d reading %s id in file %s\n", status, varstr, cdf_filename);
    return 1;
  }

  return 0;
}

int nc_inq_varid_err (int ncid, char* varstr, int* id_var, char* cdf_filename, int quiet) {
  int status = 0;

  status = nc_inq_varid (ncid, varstr, id_var);
  if (status != NC_NOERR) {
    if (!(quiet && status == NC_ENOTVAR))
      fprintf (stderr, "Error %d reading %s id in file %s\n", status, varstr, cdf_filename);
    return 1;
  }

  return 0;
}

int nc_inq_dimlen_err (int ncid, char* varstr, int id_var, char* cdf_filename, size_t* var) {
  int status = 0;

  status = nc_inq_dimlen (ncid, id_var, var);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d reading %s dimension in file %s\n", status, varstr, cdf_filename);
    return 1;
  }

  return 0;
}

int nc_get_vara_int_err (int ncid, char* varstr, int id_var, size_t* start, size_t* count, char* cdf_filename, int* var) {
  int status = 0;

  status = nc_get_vara_int (ncid, id_var, start, count, var);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d reading %s variable in file %s\n", status, varstr, cdf_filename);
    return 1;
  }

  return 0;
}

int nc_get_vara_float_err (int ncid, char* varstr, int id_var, size_t* start, size_t* count, char* cdf_filename, float* var) {
  int status = 0;

  status = nc_get_vara_float (ncid, id_var, start, count, var);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d reading %s variable in file %s\n", status, varstr, cdf_filename);
    return 1;
  }

  return 0;
}

int nc_get_vara_double_err (int ncid, char* varstr, int id_var, size_t* start, size_t* count, char* cdf_filename, double* var) {
  int status = 0;

  status = nc_get_vara_double (ncid, id_var, start, count, var);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d reading %s variable in file %s\n", status, varstr, cdf_filename);
    return 1;
  }

  return 0;
}
#endif

/**
 * optimize theta_grid
 *
 * Function to optimize the angular grid of the phase matrix.
 * As input the phase matrix is required on a fine grid with exactly
 * two or three decimals (18001/180001 points). The grid is optimized in such a way,
 * that first, grid points with zero decimals are included; second, if
 * the required accuracy is not reached, points with one decimal are
 * included; and if the accuracy is still not reached, points with two
 * decimals are included.
 *
 * @author Claudia Emde
 * @date   2009-06-29
 */
int optimize_theta_grid (double** grid_opt,
                         double** data_opt,
                         int*     Nopt,
                         double*  grid_fine,
                         double*  data_fine,
                         int      Npoints,
                         double   acc,
                         int      nthetamax,
                         int*     warnings_ntheta)

{
  int    i = 0, j = 0;
  int    indices[Npoints];
  float  data_approx[Npoints];
  int*   ind_tmp = NULL;
  double diff = 0, max_diff = 100;
  int    Nind = 2, index_acc = 0;

  /* The function can be made more general if higher accuracy is*/
  /* required but currently it is hardcoded, because of the*/
  /* requirement that the optimized grids should have as many common*/
  /* points as possible.*/
  if (Npoints != 180001 && Npoints != 18001) {
    fprintf (stderr, "Function is implemented only for a theta grid with 2 or 3 decimals accuracy.\n");
    fprintf (stderr, "Npoints must be equal to 18001 or 180001.\n");
    return -1;
  }

  indices[0] = 0;
  indices[1] = Npoints - 1;

  while (max_diff > acc && Nind < nthetamax) {
    max_diff = 0.;

    /* piecewise linear interpolation */
    for (i = 0; i < Nind - 1; i++) {
      for (j = indices[i]; j < indices[i + 1]; j++) {

        data_approx[j] = data_fine[indices[i]] + (data_fine[indices[i + 1]] - data_fine[indices[i]]) /
                                                   (grid_fine[indices[i + 1]] - grid_fine[indices[i]]) *
                                                   (grid_fine[j] - grid_fine[indices[i]]);

        diff = fabs ((data_fine[j] - data_approx[j]) / data_fine[j]);

        if (diff > max_diff) {
          max_diff  = diff;
          index_acc = j;
          /* If this routine would be used as a general grid */
          /* optimization routine without selecting prefered grid */
          /* points, the following line should be used. */
          /* indices[Nind]=j; */
        }
      }
    }

    /* Select prefered grid points: */
    /* find index for theta with 0 decimals */
    if (!find_opt_index (index_acc, indices, Nind, 0, Npoints)) /* if not found, find index for theta with 1 decimal */
      if (!find_opt_index (index_acc, indices, Nind, 1, Npoints))
        /* if not found and optimization for 3 digits,
           find index for theta with 2 decimal */
        if (!find_opt_index (index_acc, indices, Nind, 2, Npoints) && Npoints == 180001)
          /* if again not found, use the accurate index for a theta with 2 */
          /* decimals (corresponding to fine grid)*/
          indices[Nind] = index_acc;

    if (Nind > 2)
      free (ind_tmp);
    Nind += 1;

    ind_tmp = calloc (Nind, sizeof (int));

    /* sort indices */
    for (i = 0; i < Nind; i++)
      ind_tmp[i] = indices[i];

    qsort (ind_tmp, Nind, sizeof (int), int_cmp);
    for (i = 0; i < Nind; i++)
      indices[i] = ind_tmp[i];
  }

  for (i = 0; i < Nind; i++) {
    (*data_opt)[i] = data_fine[indices[i]];
    (*grid_opt)[i] = grid_fine[indices[i]];
  }

  if (Nind == nthetamax) {
    fprintf (stderr, "Warning: accuracy of phase function only %.2f percent.\n", max_diff * 100);
    *warnings_ntheta = 1;
  }

  *Nopt = Nind;

  free (ind_tmp);

  return 0;
}

/**
 * optimize_theta_grid_phamat
 *
 * Function to optimize the angular grid of the phase matrix.
 * As input the phase matrix is required on a fine grid with exactly
 * two or three decimals (18001/180001 points). The grid is optimized in such a way,
 * that first, grid points with zero decimals are included; second, if
 * the required accuracy is not reached, points with one decimal are
 * included; and if the accuracy is still not reached, points with two
 * decimals are included.
 *
 * This function is very similar to *optimize_theta_grid*, the only
 * difference is that it takes the full phase matrix as input and
 * optimizes all elements at once resulting in a common theta grid for
 * the whole phase matrix. *optimize_theta_grid* will be removed
 * later when it is not needed anymore.
 *
 * @author Claudia Emde
 * @date   2011-09-21
 */
int optimize_theta_grid_phamat (double**  grid_opt,
                                double*** data_opt,
                                int*      Nopt,
                                double*   grid_fine,
                                double**  data_fine,
                                int       Npoints,
                                double    acc,
                                int       nthetamax,
                                int       nstokes,
                                int*      warnings_ntheta)

{
  int    i = 0, j = 0, ip = 0;
  int    indices[Npoints];
  float  data_approx[Npoints];
  int*   ind_tmp = NULL;
  double diff = 0, max_diff = 100;
  int    Nind = 2, index_acc = 0;

  /* The function can be made more general if higher accuracy is*/
  /* required but currently it is hardcoded, because of the*/
  /* requirement that the optimized grids should have as many common*/
  /* points as possible.*/
  if (Npoints != 180001 && Npoints != 18001) {
    fprintf (stderr, "Function is implemented only for a theta grid with 2 or 3 decimals accuracy.\n");
    fprintf (stderr, "Npoints must be equal to 18001 or 180001.\n");
    return -1;
  }

  indices[0] = 0;
  indices[1] = Npoints - 1;

  for (ip = 0; ip < nstokes; ip++) {
    max_diff = 100.;
    //fprintf(stderr, "here !!\n");
    //fprintf(stderr, "max_diff %g acc %g Nind %d nthetamax %d\n", max_diff, acc, Nind, nthetamax);
    while (max_diff > acc && Nind < nthetamax) {
      max_diff = 0;

      /* piecewise linear interpolation */
      for (i = 0; i < Nind - 1; i++) {
        for (j = indices[i]; j < indices[i + 1]; j++) {

          data_approx[j] = data_fine[ip][indices[i]] + (data_fine[ip][indices[i + 1]] - data_fine[ip][indices[i]]) /
                                                         (grid_fine[indices[i + 1]] - grid_fine[indices[i]]) *
                                                         (grid_fine[j] - grid_fine[indices[i]]);

          diff = fabs ((data_fine[ip][j] - data_approx[j]) / data_fine[ip][j]);

          if (diff > max_diff) {
            max_diff  = diff;
            index_acc = j;
            //fprintf (stderr, "ip %d Nind %d maxdiff %g j %d\n", ip, Nind, max_diff, j);
            /* If this routine would be used as a general grid */
            /* optimization routine without selecting prefered grid */
            /* points, the following line should be used. */
            /* indices[Nind]=j; */
          }
        }
      }

      /* Select prefered grid points: */
      /* find index for theta with 0 decimals */
      if (!find_opt_index (index_acc, indices, Nind, 0, Npoints)) /* if not found, find index for theta with 1 decimal */
        if (!find_opt_index (index_acc, indices, Nind, 1, Npoints))
          /* if not found and optimization for 3 digits,
             find index for theta with 2 decimal */
          if (!find_opt_index (index_acc, indices, Nind, 2, Npoints) && Npoints == 180001)
            /* if again not found, use the accurate index for a theta with 2 */
            /* decimals (corresponding to fine grid)*/
            indices[Nind] = index_acc;

      if (Nind > 2)
        free (ind_tmp);
      Nind += 1;

      ind_tmp = calloc (Nind, sizeof (int));

      /* sort indices */
      for (i = 0; i < Nind; i++)
        ind_tmp[i] = indices[i];

      qsort (ind_tmp, Nind, sizeof (int), int_cmp);
      for (i = 0; i < Nind; i++) {
        indices[i] = ind_tmp[i];
      }
    }
  }

  for (i = 0; i < Nind; i++) {
    for (ip = 0; ip < nstokes; ip++)
      (*data_opt)[ip][i] = data_fine[ip][indices[i]];
    (*grid_opt)[i] = grid_fine[indices[i]];
  }

  if (Nind == nthetamax) {
    fprintf (stderr, "Warning: accuracy of phase function only %.2f percent.\n", max_diff * 100);
    *warnings_ntheta = 1;
  }

  *Nopt = Nind;

  free (ind_tmp);

  return 0;
}

/**
 * find_opt_index
 *
 * Find index for a grid point with given decimals (used in
 * optimize_theta_grid).
 *
 * @author Claudia Emde
 * @date   2009-06-29
 */
int find_opt_index (int index_acc, int* indices, int Nind, int decimal, int Npoints) {
  int index_low = 0, index_up = 0, lower = 0, upper = 0;
  int factor = 0;

  /* Settings for 180001 points*/
  if (decimal == 0)
    factor = 1000;
  else if (decimal == 1)
    factor = 100;
  else if (decimal == 2)
    factor = 10;
  else
    fprintf (stderr, "Error: find_opt_index() is only implemented for decimal 0, 1 or 2\n");

  if (Npoints == 18001)
    factor /= 10;

  /* Indices of round angles below and above j*/
  index_low = floor (index_acc / factor) * factor;
  index_up  = index_low + factor;

  /* Check whether indices are already in optimized grid */
  lower = element (index_low, indices, Nind);
  upper = element (index_up, indices, Nind);

  if (!lower || !upper) {
    if (index_low - index_acc == 0)
      indices[Nind] = index_low;
    else if (index_up - index_acc == 0)
      indices[Nind] = index_up;
    else if ((index_up - index_acc) >= 0.5 * factor && !(lower))
      indices[Nind] = index_low;
    else if (!upper)
      indices[Nind] = index_up;
    else if (!upper)
      indices[Nind] = index_up;
    else if (!lower)
      indices[Nind] = index_low;
  }
  /* No index found */
  else
    return 0;

  /* Index has been found and added to indices vector */
  return 1;
}

/**
 * element
 *
 * This function checks whether a is an element of the vector b.
 *
 * @return 1 if a in b, else 0
 *
 * @author Claudia Emde
 * @date   2009-06-29
*/
int element (int a, int* b, int L) {
  int i      = 0;
  int status = 0;

  while (i < L) {
    if (a == b[i]) {
      status = 1;
      break;
    }
    i++;
  }
  return status;
}

/**
 * int_cmp
 *
 * qsort int comparison function
 *
 * Taken from http://www.anyexample.com/programming/c/.
 *
 * @author included by Claudia Emde
 * @date   2009-06-29
*/
int int_cmp (const void* a, const void* b) {
  const int* ia = (const int*) a; /* casting pointer types*/
  const int* ib = (const int*) b;
  return *ia - *ib;
  /* integer comparison: returns negative if b > a
     and positive if a > b */
}

/* refractive index according to IAPWS */
/* Veronika Poertge, 14.4.2021         */
int refwat_iapws (float wavelength, float temperature, float* ref_re, float* ref_im) {
  /**
  *
  * @author Veronika Poertge
  * @date   2020-10-14
  * Taken from http://www.iapws.org/relguide/rindex.pdf
  */

  char function_name[] = "refwat_iapws";
  char file_name[]     = "src_c/miecalc.c";

  void F77_FUNC (wrefwat, WREFWAT) (float* wl, float* temp, float* re, float* im);

  /* first calculcate standard Warren et al refractive index because IAPWS does not provide imaginary part */
  F77_FUNC (wrefwat, WREFWAT) (&wavelength, &temperature, ref_re, ref_im);
  if (*ref_re == 0.0 && *ref_im == 0.0) {
    fprintf (stderr, "Error, during execution of 'wrefwat' in %s (%s)\n", function_name, file_name);
    CHKERR (-1);
  }

  /* now over-write real part of the refractive index with IAPWS data */
  if (temperature <= 233.15 || temperature > 360) {
    /* valid range of temperature: -12 deg C < t ≤ 86.85 deg C*/
    fprintf (stderr,
             "Error: Temperature outside range of tabulated and extrapolated IAPWS95 "
             "relationship between temperature, pressure and density. Valid range of temperature: "
             "-40 deg C ≤ t ≤ 86.85 deg C!\n");
    CHKERR (-1);
  }

  if (temperature >= 233.15 && temperature <= 261.15) {
    /* range of extrapolated temperature: -40 deg C < t ≤ -12 deg C*/
    fprintf (stderr,
             "Temperature inside range of extrapolated IAPWS95 relationship between "
             "temperature, pressure and density. (Extrapolated) range of temperature: -40 deg C ≤ t ≤ -12 deg C!\n"
             "Result of the refractive index could be slightly incorrect.\n");
  }

  if (wavelength < 0.2 || wavelength > 1.1) {
    /* valid range of wavelength: 200 nm ≤ wavelength ≤ 1100 nm */
    fprintf (stderr,
             "Error: Wavelength outside range of IAPWS formula:"
             "valid range of wavelength: 200 nm ≤ wavelength ≤ 1100 nm!\n");
    CHKERR (-1);
  }

  const double a_0 = 0.244257733, a_1 = 9.74634476E-3, a_2 = -3.73234996E-3;
  const double a_3 = 2.68678472E-4, a_4 = 1.58920570E-3, a_5 = 2.45934259E-3;
  const double a_6 = 0.900704920, a_7 = -1.66626219E-2;
  const double lambda_UV = 0.2292020, lambda_IR = 5.432937;
  const double rho_prime            = 1000 /* kg/m**3 */;
  const double T_prime              = 273.15 /* K */;
  const double wavelength_prime     = 0.589 /* micrometers  */;
  const double temp_ratio           = temperature / T_prime;
  const double wavelength_ratio     = wavelength / wavelength_prime;
  double       density_interpolated = 0;
  int          status;
  const int    number_points = 37;

  /*
  * calculate density from temperature (IAPWS95 formulation)
  * density_array and temperature_array come from
  * https://aip.scitation.org/doi/10.1063/1.1461829 Page 497 (table for 0.101 325 MPa)
  * and extrapolated to -12 deg C, compared with Kell 1975, see Figure 7.49 in IAPWS.
  * d0 array and d1_array were defined via the linear_coeffc_fast function
  */
  double temperature_array[] = {233.15, 234.15, 235.15, 237.15, 239.15, 241.15, 243.15, 245.15, 247.15,  249.15,
                                251.15, 253.15, 255.15, 257.15, 259.15, 261.15, 265,    270,    273.153, 275,
                                280,    285,    290,    295,    300,    305,    310,    315,    320,     325,
                                330,    335,    340,    345,    350,    355,    360};

  double d0_array[] = {959.241, 964.256, 967.987, 973.204, 977.341, 980.839, 983.832, 986.402, 988.612, 990.515, 992.156, 993.57,
                       994.789, 995.836, 996.731, 997.490, 998.620, 999.536, 999.843, 999.938, 999.911, 999.517, 998.804, 997.807,
                       996.557, 995.076, 993.384, 991.496, 989.427, 987.187, 984.787, 982.234, 979.536, 976.699, 973.728, 970.629};

  double d1_array[] = {5.015,    3.731,   2.6085,  2.0685,  1.749,   1.4965,  1.285,    1.105,    0.9515,
                       0.8205,   0.707,   0.6095,  0.5235,  0.4475,  0.3795,  0.293506, 0.183200, 0.097368,
                       0.051435, -0.0054, -0.0788, -0.1426, -0.1994, -0.25,   -0.2962,  -0.3384,  -0.3776,
                       -0.4138,  -0.448,  -0.48,   -0.5106, -0.5396, -0.5674, -0.5942,  -0.6198,  -0.645};

  //uncomment
  //double *d0 = NULL, *d1 = NULL;

  //recomment this:
  double* d0 = d0_array;
  double* d1 = d1_array;

  /*
  double density_array[] = {959.241, 964.256, 967.987, 973.204, 977.341, 980.839, 983.832, 986.402, 988.612, 990.515,
                            992.156, 993.57,  994.789, 995.836, 996.731, 997.490, 998.620, 999.536, 999.843, 999.938,
                            999.911, 999.517, 998.804, 997.807, 996.557, 995.076, 993.384, 991.496, 989.427, 987.187,
                            984.787, 982.234, 979.536, 976.699, 973.728, 970.629, 967.404};

  status = linear_coeffc_fast (temperature_array, density_array, number_points, &d0, &d1);
  for (int ip = 0; ip < number_points-1; ip++)
      printf("d0 %f d1 %f \n", d0[ip], d1[ip]);
  */

  status = calc_linear_value (temperature, &density_interpolated, temperature_array, number_points, d0, d1);
  CHKERR (status);

  // Calculate refractive index from IAPWS formula
  const double rho_ratio  = density_interpolated / rho_prime;
  const double pow_lambda = pow (wavelength_ratio, 2);
  const double part0 = a_0, part1 = a_1 * rho_ratio, part2 = a_2 * temp_ratio;
  const double part3 = a_3 * pow_lambda * temp_ratio, part4 = a_4 / pow_lambda;
  const double part5 = a_5 / (pow_lambda - pow (lambda_UV, 2));
  const double part6 = a_6 / (pow_lambda - pow (lambda_IR, 2)), part7 = a_7 * pow (rho_ratio, 2);
  const double y = (part0 + part1 + part2 + part3 + part4 + part5 + part6 + part7) / (1 / rho_ratio);
  *ref_re        = sqrt (2 * y + 1) / sqrt (1 - y);

  return 0;
}
