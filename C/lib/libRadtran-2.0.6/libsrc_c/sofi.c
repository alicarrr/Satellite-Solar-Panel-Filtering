/************************************************************************/
/* $Id$                 */
/*                                                                      */
/* Functions to calculate radiation in the umbral shadow during a       */
/* total eclipse using MYSTIC.                                          */
/*                                                                      */
/* Authors: Claudia Emde, Paul Ockenfuss                                */
/*                                                                      */
/* Correspondence: claudia.emde@lmu.de                                  */
/*                                                                      */
/************************************************************************/
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ascii.h"
#include "common_math.h"
#include "sofi.h"
#include "uvspecrandom.h"

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif

int    limbdarkening_Koepke (double lambda, double Nr, double* Sl);
int    limbdarkening_Pierce (double lambda, double sdist, int Nr, double* Sl);
int    limbdarkening_Neckel (double lambda, double sdist, int Nr, double* Sl);
double calcPsi (double X, double sdist);

/***********************************************************************************/
/* Function: sample_photons_sofi                                                   */
/* Description:                                                                    */
/*  Generate a probability distribution functions which corresponds to the         */
/*  irradiance distribution at the top pf the atmosphere during a the total        */
/*  eclipse of the sun.                                                            */
/*  The integrated probability is returned.                                        */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                                 */
/***********************************************************************************/
int sample_photons_sofi (float lambda, double dx, double ratio, double sdist, int Nx, int limb_darkening, double* pd) {

  int     status = 0;
  double *x, *r;
  double *Ic, *Sl, *alpha;
  int     i = 0, j = 0;
  double  I_int                 = 0.;
  FILE*   file                  = NULL;
  int     write_solar_intensity = 0;

  /* Number of sun radius grid points */
  int Nr = 5001;

  /*extraterrestrial solar irradiance */
  double E_0 = 1368.0;

  /* Intensity in centre of solar disk */
  double I_0;

  /* This is not needed anymore, it is calculated in sofi main
     function based on sun and moon positions */
  /* FIX??  ratio between moon diameter and sun diameter */
  /*double ratio=1.0494;  29-03-2006 in Turkey.*/
  /*double ratio=0.965;  Ratio for solar eclipse over south pole */
  /* calculated from data from Nasa report, corresponds to approx 1km */
  // CE FIX double dx=0.0006;
  /*double dx = 5.56e-4; */

  /* Sun disk radius normalized to x */
  double r_M = ratio; /* PO: Alle Größen in Einheiten von r_S*/

  double epsilon = 1.e-8;

  /* Distance between centre of moon and centre of sun*/
  x = calloc (Nx + 1, sizeof (double));

  /* Irradiance as a function of x*/
  Ic = calloc (Nx + 1, sizeof (double));

  /* Initialize x. This is *not* X in Koepke! Normalisation not included. 
     X=x/(r_M+r_S) */
  x[0] = 0.;
  for (j = 1; j < Nx + 1; j++)
    x[j] = x[j - 1] + dx;

  /* Grid for radius on solar disk, needed for integration of apparant irradiance. 
     Also in units of apparent moon radius.*/
  /* FIXCE PO: nicht eig. sun radius? */
  r = calloc (Nr + 1, sizeof (double));
  for (i = 0; i < Nr; i++) /*FIXCE?? <Nr+1 in old version */
  {
    r[i] = ((double) i) / (Nr - 1);
  }

  /* Angle for integration of solar disk.*/
  alpha = calloc (Nr + 1, sizeof (double));

  /* Solar limb darkening */
  /* Gamma in Koepke Paper: relative Helligkeit zum Zentrum in
     Abhängigkeit von der Position auf der Sonne */
  Sl = calloc (Nr + 1, sizeof (double));

  switch (limb_darkening) {
  case 1:
    status = limbdarkening_Pierce (lambda, sdist, Nr, Sl);
    break;
  case 2:
    status = limbdarkening_Neckel (lambda, sdist, Nr, Sl);
    break;
  case 3:
    status = limbdarkening_Koepke (lambda, Nr, Sl);
    break;
  default:
    fprintf (stderr, "Error, unknown limb darkening function. \n");
    return -1;
  }

  /* Print solar irradiance of uncovered solar disk */
  /*   for (i=0; i<Nr; i++)
       fprintf(stderr, " %g %g \n",  r[i], Sl[i]);*/

  /* Integrate over solar disk to get I_0, take numerical integration
     for consistency */
  for (i = 1; i < Nr; i++) /* I_int+=2*PI*(Sl[i]*r[i]+Sl[i+1]*r[i+1])*0.5*(r[i+1]-r[i]); FIXCE: Why changed? */
    I_int += 2 * PI * (Sl[i] * r[i]) * (r[1] - r[0]); /* Integral 2PI*Sl*r*dr */

  /* Intensity in centre of solar disk */
  I_0 = E_0 / I_int;

  /* Calculate irradiance as a function of X*/
  for (j = 0; j < Nx; j++) {

    /* Moon has passed the sun, so we get the total intensity */
    if (x[j] >= r_M + 1.0)
      Ic[j] = E_0;

    /* Moon in front of sun, integration over partly covered disk*/
    else {

      for (i = 0; i < Nr; i++) {

        if (x[j] < r_M + 1.0) {

          if ((r[i] < (x[j] - r_M) && x[j] > r_M) || (r[i] > (r_M + x[j]) && x[j] < 1.0 - r_M))
            /* second part applies for ring eclipse */

            alpha[i] = PI;

          else if (r[i] < (r_M - x[j]) + epsilon && x[j] < r_M + epsilon)
            alpha[i] = 0.0;

          else
            alpha[i] = acos ((r_M * r_M - r[i] * r[i] - x[j] * x[j]) / (2 * r[i] * x[j]));

        } else
          fprintf (stderr, "Forgotten case x %g r %g \n. This should not happen!!!", x[j], r[i]);
      }

      /* Integrate uncovered part */
      for (i = 0; i < Nr; i++) {

        /* Print alpha for x[j] */
        /* if( j==1 ) */
        /*   fprintf(stderr, " %g %g %g \n", x[j], r[i], alpha[i]); */

        /* Ic[j] += 2.0 * I_0 * ( ( alpha[i] * Sl[i] * r[i] ) +  */
        /*                        ( alpha[i+1] * Sl[i+1] * r[i+1] ) ) * */
        /*   0.5 * ( r[i+1] - r[i] );  FIXCE: Why changed?? */

        Ic[j] += 2.0 * I_0 * (alpha[i] * Sl[i] * r[i]) * (r[1] - r[0]);
      }
    }
  }

  /* Photon density function */
  for (j = 0; j < Nx; j++)
    pd[j] = Ic[j] / E_0;

  /*  pd[i] = (Ic[i]+1.1e-7*I_0)/E_0;
      1e-7: Upper limit of corona intensity.*/

  /* Print solar intensity (weighting function) as function of distance to file */
  if (write_solar_intensity) {
    if ((file = fopen ("solar_intensity_pd.dat", "w")) == NULL)
      return ASCIIFILE_NOT_FOUND;

    for (j = 0; j < Nx; j++)
      fprintf (file, "%g %g %g \n", x[j], Ic[j], pd[j]);
    fclose (file);
  }

  /* Free variables.*/
  free (r);
  free (alpha);
  free (Sl);

  return status;
}

/***********************************************************************************/
/* Function: generate_photon_sofi                                         @62_30i@ */
/* Description: Generate photons to simulate solar eclipse (forward
   Monte Carlo). This function has not been used in solar eclipse
   simulations because it is too slow.  */
/* */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void generate_photon_sofi (atmosphere_struct* atmos, photon_struct* p, double sza, double phi, double* pd) {
  double sinphi = 0, cosphi = 0, cossza = 0, sinsza = 0;
  double alpha = 0, x = 0, y = 0, r = 0, dx = 0, dy = 0;
  int    i_re = 0;

  /* sine and cosine of solar zenith */
  cossza = cosd (sza);
  sinsza = sind (sza);

  /* Different definition of azimuth angle */
  phi = phi - 180.;

  /* sine and cosine of solar azimuth */
  cosphi = cosd (phi);
  sinphi = sind (phi);

  /* Center of umbral shadow at TOA */
  r  = tand (sza) * atmos->Z[atmos->Nz];
  dx = r * sinphi;
  dy = r * cosphi;

  /* x and y coordinates, (x,y)=0 corresponds to center of ellipse */
  x = p->x[0] - 200e3 + dx;
  y = p->x[1] - 200e3 + dy;

  /* Starting point in polar coordinates */
  r     = sqrt (x * x + y * y);
  alpha = atand (x / y);

  /* Equivalent radius, if the ellipsoidal shadow is transformed back to a circle, 
     for which the probability density function is given.*/
  i_re = (int) fabs (sqrt (pow (r, 2) / (pow (cosd (alpha - phi), 2) / pow (cossza, 2) + pow (sind (alpha - phi), 2)))) / 1000.;

  /* photon weight */
  p->weight *= pd[i_re];

  /* initialize direction using the specified solar zenith and azimuth angles */
  /* and copy direction to "initial direction" struct                         */
  init_direction (sinsza, cossza, sinphi, cosphi, &(p->dir));
  cp_direction (&(p->dir0), &(p->dir));

  /* vertical start position: TOA */
  /* Make always sure to call set_photon_z() AFTER the new photon     */
  /* direction has been assigned because the start index and position */
  /* might depend on direction                                        */
  set_photon_z (atmos->Z[atmos->Nz], atmos, p);
}

/***********************************************************************************/
/* Function: limbdarkening_Koepke                                                  */
/*    Description:                                                                 */
/*       Limb darkening of sun calculated following                                */
/*       Koepke, P., Reuder, J., and Schween, J.: Spectral variation of the so-    */
/*            lar radiation during an eclipse, Meteorologische Zeitschrift, 10,    */
/*            179–186, 2001.                                                       */
/*                                                                                 */
/* Input:   lambda - wavelength in nm                                              */
/*          Nr     - Number of points where Sl will be calculated                  */
/*                    (equidistant X[0,1] from center of solar disk to limb)       */
/* Output:  Sl     - Solar limb darkening                                          */
/***********************************************************************************/
int limbdarkening_Koepke (double lambda, double Nr, double* Sl) {
  /* Definition of constants:*/
  double planck         = 6.63e-34;
  double speed_of_light = 2.998e8;
  double boltzmann      = 1.38e-23;
  double T_S            = 5740.0;
  double beta           = 3.0 * planck * speed_of_light * sqrt (sqrt (2.0)) / (8.0 * boltzmann * lambda * 1e-9 * T_S);
  double X              = 0;

  for (int i = 0; i < Nr; i++) {
    X     = ((double) i) / (Nr - 1);
    Sl[i] = (1 + beta * sqrt (fabs (1. - X * X))) / (1 + beta);
  }

  return 0;
}

/***********************************************************************************/
/* Function: limbdarkening_Neckel                                                  */
/* Description:                                                                    */
/*       Limb darkening of sun calculated following                                */
/*            Neckel, H. (2005). “Analytical reference functions F(λ) for the      */
/*	      sun’s limb darkening and its absolute continuum intensities”.        */
/*	      In: Solar Physics 229, pp. 13–33.                                    */
/*                                                                                 */
/*                                                                                 */
/* Input:   lambda - wavelength in nm                                              */
/*          Nr     - Number of points where Sl will be calculated                  */
/*                    (equidistant X[0,1] from center of solar disk to limb)       */
/* Output:  Sl     - Solar limb darkening                                          */
/***********************************************************************************/
int limbdarkening_Neckel (double lambda, double sdist, int Nr, double* Sl) {
  double a00, a01, a10, a11, a15, a20, a25, a30, a35, a40, a45, a50, a55;
  double lambda_1 = 1. / lambda * 1000.;
  double lambda_5 = 1. / pow ((lambda * 1e-3), 5);
  if (lambda >= 300 && lambda < 372.98) {
    a00 = 0.35601;
    a01 = -0.085217;
    a10 = +1.11529;
    a11 = +0.085217;
    a15 = -0.001871;
    a20 = -0.67237;
    a25 = +0.003589;
    a30 = +0.18696;
    a35 = -0.002415;
    a40 = +0.00038;
    a45 = +0.000897;
    a50 = +0.01373;
    a55 = -0.000200;
  } else if (lambda >= 385 && lambda < 422.57) {
    a00 = +0.09900;
    a01 = +0.010643;
    a10 = +1.96884;
    a11 = -0.010643;
    a15 = -0.009166;
    a20 = -2.80511;
    a25 = +0.024873;
    a30 = +3.32780;
    a35 = -0.029317;
    a40 = -2.17878;
    a45 = +0.018273;
    a50 = +0.58825;
    a55 = -0.004663;
  } else if (lambda >= 422.57 && lambda < 1100) {
    a00 = +0.75267;
    a01 = -0.265577;
    a10 = +0.93874;
    a11 = +0.265577;
    a15 = -0.004095;
    a20 = -1.89287;
    a25 = +0.012582;
    a30 = +2.42234;
    a35 = -0.017117;
    a40 = -1.71150;
    a45 = +0.011977;
    a50 = +0.49062;
    a55 = -0.003347;
  } else {
    fprintf (stdout, "Error in Neckel: no data available for this wavelength.");
    return -1;
  }
  double a0 = a00 + a01 * lambda_1;
  double a1 = a10 + a11 * lambda_1 + a15 * lambda_5;
  double a2 = a20 + a25 * lambda_5;
  double a3 = a30 + a35 * lambda_5;
  double a4 = a40 + a45 * lambda_5;
  double a5 = a50 + a55 * lambda_5;

  double X  = 0;
  double mu = 0;
  for (int i = 0; i < Nr; i++) {
    X     = ((double) i) / (Nr - 1);
    mu    = cos (calcPsi (X, sdist));
    Sl[i] = a0 + a1 * mu + a2 * mu * mu + a3 * pow (mu, 3) + a4 * pow (mu, 4) + a5 * pow (mu, 5);
  }

  return 0;
}

/***********************************************************************************/
/* Function: limbdarkening_Pierce                                                  */
/* Description:                                                                    */
/*       Limb darkening of sun calculated following                                */
/*         Pierce, A. K. and C. D. Slaughter (1977). "Solar limb                   */
/*         darkening". In: Solar Physics 51, pp. 25–41.                            */
/*                                                                                 */
/*                                                                                 */
/* Input:   lambda - wavelength in nm                                              */
/*          sdist  - Distance Sun-Earth in km                                      */
/*          Nr     - Number of points where Sl will be calculated                  */
/*                    (equidistant X[0,1] from center of solar disk to limb)       */
/* Output:  Sl     - Solar limb darkening                                          */
/***********************************************************************************/
int limbdarkening_Pierce (double lambda, double sdist, int Nr, double* Sl) {
  double* wvl    = NULL;
  double* number = NULL;
  double* a      = NULL;
  double* b      = NULL;
  double* c      = NULL;
  double* d      = NULL;
  double* e      = NULL;
  double* f      = NULL;
  double* pe_x   = NULL;
  int     L      = 0;
  int     status = 0;

  status = read_9c_file ("PierceCoefficients.csv", &wvl, &number, &a, &b, &c, &d, &e, &f, &pe_x, &L);
  if (status != 0) {
    printf ("Error: Table with coefficients for Pierce parametrization not found.");
    return status;
  }
  if (lambda < 304 || lambda >= 1046) {
    printf ("Error: Pierce parametrization not implemented for %.3fnm", lambda);
    return -1;
  }
  int ind = 0;
  for (int i = 1; i < L; i++) {
    if (lambda < wvl[i]) {
      ind = i - 1;
      break;
    }
  }

  double mu = 0;
  double X  = 0;
  double pstart, pend;
  for (int i = 0; i < Nr; i++) {
    X      = ((double) i) / (Nr - 1);
    mu     = cos (calcPsi (X, sdist));
    pstart = a[ind] + b[ind] * mu + c[ind] * mu * mu + d[ind] * pow (mu, 3) + e[ind] * pow (mu, 4) + f[ind] * pow (mu, 5);
    pend   = a[ind + 1] + b[ind + 1] * mu + c[ind + 1] * mu * mu + d[ind + 1] * pow (mu, 3) + e[ind + 1] * pow (mu, 4) +
           f[ind + 1] * pow (mu, 5);
    Sl[i] = pstart + (pend - pstart) / (wvl[ind + 1] - wvl[ind]) * (lambda - wvl[ind]);
  }
  free (wvl);
  free (number);
  free (a);
  free (b);
  free (c);
  free (d);
  free (e);
  free (f);
  free (pe_x);
  return 0;
}

/***********************************************************************************/
/* Function: calcPsi                                                               */
/* Description: Calculate Psi from Neckel Limb Darkening                           */
/*                                                                                 */
/* Return: Psi                                                                     */
/* Input:  X - Position on sun disk in units of apparent Sun Radius                */
/*         sdist - sun Distance im km                                              */
/*                                                                                 */
/***********************************************************************************/
double calcPsi (double X, double sdist) {
  /* Weg Nr 1: Exakt */
  double Rs   = 695700.;                      /* Sunradius in km */
  double RadS = Rs * sin (acos (Rs / sdist)); /* Apparent Sun Radius in km */
  double r    = X * RadS;                     /* Position on Sun disk in km */
  double a    = sdist - sqrt (Rs * Rs - r * r);
  return PI - atan (a / r) - acos (r / Rs);
  /* Weg Nr 2: Approximation für sdist=>infinity (nicht vollständig der Fall) */
  /* return PI/2.-acos(X); */
}
