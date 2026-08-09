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

/* full solution of the multi-layer twostream equation solar + thermal */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "equation.h"
#include "twostrebe.h"
#include "cdisort.h"

static int adding (double*  dtau_org,
                   double*  omega0_org,
                   double*  g_org,
                   int      nlev,
                   double   S0,
                   double   mu0,
                   double   Ag,
                   int      delta,
                   double** Edir,
                   double** Edn,
                   double** Eup,
                   double** Lavg);

static int twostream (double*  dtau_org,
                      double*  omega0_org,
                      double*  g_org,
                      int      nlev,
                      double   S0,
                      double   mu0,
                      double   Ag,
                      double   Bg,
                      double*  B,
                      int      delta,
                      double** Edir,
                      double** Edn,
                      double** Eup,
                      double** Lavg);

int twostrebe (float* dtau_org,
               float* omega0_org,
               float* g_org,
               int    nlev,
               double S0,
               double mu0,
               double Ag,
               int    planck,
               int    delta,
               int    nzout,
               float* zd,
               float* temper,
               float  btemp,
               float  wvnmlo,
               float  wvnmhi,
               float* zout,
               float* fldn,
               float* flup,
               float* fldir,
               float* uavg) {
  int status = 0;

  double* dtau_org_d   = calloc (nlev - 1, sizeof (double));
  double* omega0_org_d = calloc (nlev - 1, sizeof (double));
  double* g_org_d      = calloc (nlev - 1, sizeof (double));
  double* B            = calloc (nlev, sizeof (double));
  double *Edir = NULL, *Edn = NULL, *Eup = NULL, *Lavg = NULL;

  double taumax = 100.0;
  double Bg     = 0;

  if (planck) {
    for (int ilev = 0; ilev < nlev; ilev++) /* level temperatures and Planck functions */
      B[ilev] = c_planck_func1 (wvnmlo, wvnmhi, temper[ilev]);

    /* surface temperature and Planck function */
    Bg = c_planck_func1 (wvnmlo, wvnmhi, btemp);
  }

  /* copy float to double arrays */
  for (int ilyr = 0; ilyr < nlev - 1; ilyr++) {
    dtau_org_d[ilyr] = dtau_org[ilyr];

    /* restrict layer optical thickness to 100 */
    if (dtau_org_d[ilyr] > taumax)
      dtau_org_d[ilyr] = taumax;

    omega0_org_d[ilyr] = omega0_org[ilyr];
    g_org_d[ilyr]      = g_org[ilyr];
  }

  /* call twostream code: original solution with 5-diagonal matrix for thermal and adding() for solar */
  if (planck)
    status = twostream (dtau_org_d, omega0_org_d, g_org_d, nlev, S0, mu0, Ag, Bg, B, delta, &Edir, &Edn, &Eup, &Lavg);
  else
    status = adding (dtau_org_d, omega0_org_d, g_org_d, nlev, S0, mu0, Ag, delta, &Edir, &Edn, &Eup, &Lavg);

  if (status != 0) {
    fprintf (stderr, "Error %d returned by twostream()\n", status);
    return status;
  }

  /* copy results to final fields */
  for (int ilev = 0; ilev < nlev; ilev++)
    for (int lu = 0; lu < nzout; lu++)
      // avoid roundoff errors
      if (fabs (zout[lu] - zd[ilev]) <= 1E-6 * zout[lu]) {
        fldn[lu]  = Edn[ilev];
        flup[lu]  = Eup[ilev];
        fldir[lu] = Edir[ilev];
        uavg[lu]  = Lavg[ilev];
      }

  /* free memory */
  free (Edir);
  free (Edn);
  free (Eup);
  free (Lavg);
  free (dtau_org_d);
  free (omega0_org_d);
  free (g_org_d);
  free (B);

  return 0;
}

/* Bernhard's twostream routine - copied from FPdA */
static int twostream (double*  dtau_org,
                      double*  omega0_org,
                      double*  g_org,
                      int      nlev,
                      double   S0,
                      double   mu0,
                      double   Ag,
                      double   Bg,
                      double*  B,
                      int      delta,
                      double** Edir,
                      double** Edn,
                      double** Eup,
                      double** Lavg) {
  int status = 0;

  double* res = NULL;

  double alpha1 = 0, alpha2 = 0, alpha3 = 0, alpha4 = 0, alpha5 = 0, alpha6 = 0, alpha7 = 0, alpha8 = 0, alpha9 = 0, alpha10 = 0;
  double a11 = 0, a12 = 0, a13 = 0, a23 = 0, a33 = 0, a33_unscaled = 0;
  double lambda = 0, kappa = 0, b = 0, A = 0;
  double denom = 0;

  double dtau = 0, omega0 = 0, g = 0;

  double b1 = 0, b2 = 0;

  double B0 = 0, B1 = 0;

  double** AA = NULL;
  double*  bb = NULL;
  double * S = NULL, *S_unscaled = NULL;

  /* allocate memory for equation system */
  AA = calloc (2 * nlev, sizeof (double*));
  for (int ilev = 0; ilev < 2 * nlev; ilev++)
    if ((AA[ilev] = calloc (5, sizeof (double))) == NULL) {
      fprintf (stderr, "Error allocating memory for AA[%d]\n", ilev);
      return -1;
    }

  bb         = calloc (2 * nlev, sizeof (double));
  S          = calloc (2 * nlev, sizeof (double));
  S_unscaled = calloc (2 * nlev, sizeof (double));

  S[0]          = S0;
  S_unscaled[0] = S0;

  /* delta scaling */

  for (int ilyr = 0; ilyr < nlev - 1; ilyr++) {

    if (delta) {
      delta_scale_hg (dtau_org[ilyr], omega0_org[ilyr], g_org[ilyr], &dtau, &omega0, &g);
    } else {
      dtau   = dtau_org[ilyr];
      omega0 = omega0_org[ilyr];
      g      = g_org[ilyr];
    }

    /* eddington */

    /* restrict omega0 to 0.999999 */
    omega0 = (omega0 < 0.999999 ? omega0 : 0.999999);

    if (omega0 == 0) {
      alpha1 = 1.75;
      alpha2 = -0.25;
    } else {
      alpha1 = (1.0 - omega0) + 0.75 * (1.0 - omega0 * g);
      alpha2 = -(1.0 - omega0) + 0.75 * (1.0 - omega0 * g);
    }

    lambda = sqrt (alpha1 * alpha1 - alpha2 * alpha2);

    if (alpha2 != 0) {
      A = 1.0 / (alpha2 / (alpha1 - lambda) * exp (lambda * dtau) - alpha2 / (alpha1 + lambda) * exp (-lambda * dtau));

      a11 = A * 2.0 * lambda / alpha2;
      a12 = A * (exp (lambda * dtau) - exp (-lambda * dtau));
    } else {  // special solution for alpha2 == 0 which causes A = 0/0 == nan
      a11 = exp (-lambda * dtau);
      a12 = 0;
    }

    /* calculate B0, B1 */
    if (dtau > 0.01) {
      B1 = (B[ilyr + 1] - B[ilyr]) / dtau;
      B0 = B[ilyr];
    } else {
      B1 = 0;
      B0 = (B[ilyr + 1] + B[ilyr]) * 0.5;
    }

    /* set thermal source */
    kappa  = 2.0 * M_PI * (1.0 - omega0);
    alpha7 = kappa * (B0 / (alpha1 - alpha2) + B1 / lambda / lambda);
    alpha9 = kappa * (B0 / (alpha1 - alpha2) - B1 / lambda / lambda);

    alpha8  = kappa * B1 / (alpha1 - alpha2);
    alpha10 = alpha8;

    b1 = -a11 * (alpha7 + alpha8 * dtau) - a12 * (alpha9) + alpha7;
    b2 = -a12 * (alpha7 + alpha8 * dtau) - a11 * (alpha9) + alpha9 + alpha10 * dtau;

    AA[2 * ilyr][3]     = -a12;
    AA[2 * ilyr][4]     = -a11;
    AA[2 * ilyr + 3][0] = -a11;
    AA[2 * ilyr + 3][1] = -a12;

    bb[2 * ilyr]     = b1;
    bb[2 * ilyr + 3] = b2;

    /* add solar source */
    if (S0 > 0) { /* solar component */
      b      = 0.5 - 0.75 * g * mu0;
      alpha3 = -omega0 * b;
      alpha4 = omega0 * (1 - b);

      denom  = (1.0 / mu0 / mu0 - lambda * lambda);
      alpha5 = ((alpha1 - 1.0 / mu0) * alpha3 - alpha2 * alpha4) / denom;
      alpha6 = (alpha2 * alpha3 - (alpha1 + 1.0 / mu0) * alpha4) / denom;

      a33          = exp (-dtau / mu0);
      a33_unscaled = exp (-dtau_org[ilyr] / mu0);

      a13 = alpha5 * (1.0 - a11 * a33) - alpha6 * a12;
      a23 = -a12 * alpha5 * a33 + alpha6 * (a33 - a11);

      bb[2 * ilyr] += a13 * S[ilyr];
      bb[2 * ilyr + 3] += a23 * S[ilyr];

      S[ilyr + 1]          = S[ilyr] * a33;
      S_unscaled[ilyr + 1] = S_unscaled[ilyr] * a33_unscaled;
    }
  }

  /* diagonal elements */
  for (int ilev = 0; ilev < 2 * nlev; ilev++)
    AA[ilev][2] = 1.0;

  /* boundary condition: surface albedo */
  AA[2 * nlev - 2][3] = -Ag;
  bb[2 * nlev - 2]    = (1.0 - Ag) * M_PI * Bg;
  if (S0 > 0)
    bb[2 * nlev - 2] += S[nlev - 1] * mu0 * Ag;

  status = solve_five_ms (AA, bb, 2 * nlev, &res);
  if (status != 0) {
    fprintf (stderr, "Error %d solving equation system\n", status);
    return status;
  }

  /* allocate memory for result */
  *Edir = calloc (nlev, sizeof (double));
  *Edn  = calloc (nlev, sizeof (double));
  *Eup  = calloc (nlev, sizeof (double));
  *Lavg = calloc (nlev, sizeof (double));

  for (int ilev = 0; ilev < nlev; ilev++) {
    /* direct irradiance with scaled optical thickness */
    /* (*Edir)[ilev]          = S[ilev]*mu0;           */
    /* (*Edn)[ilev]           = res[2*ilev+1];         */

    /* direct irradiance with unscaled optical thickness */
    (*Edir)[ilev] = S_unscaled[ilev] * mu0;
    (*Edn)[ilev]  = res[2 * ilev + 1] + (S[ilev] - S_unscaled[ilev]) * mu0;

    (*Eup)[ilev] = res[2 * ilev];

    (*Lavg)[ilev] = S[ilev] / 4.0 / M_PI + (res[2 * ilev] + res[2 * ilev + 1]) / 2.0 / M_PI;
    /* we see large differences compared to disort and twostr for the averaged radiance; */
    /* need to check if that is inherent to the delta-Eddington method  (probably yes)   */
    /* or due to an error in the implementation; so long, we provide NAN instead!        */
    (*Lavg)[ilev] = NAN;
  }

  free (res);
  for (int ilev = 0; ilev < 2 * nlev; ilev++)
    free (AA[ilev]);

  free (AA);
  free (bb);
  free (S);
  free (S_unscaled);

  return 0;
}

void delta_scale_hg (double tau, double ssa, double g, double* tauscale, double* ssascale, double* gscale) {
  double f = g * g;

  *tauscale = (1.0 - ssa * f) * tau;
  *ssascale = (1.0 - f) * ssa / (1.0 - ssa * f);
  *gscale   = (g - f) / (1.0 - f);
}

/* calculate layer properties t, r, rdir, sdir, and tdir from            */
/* layer optical thickness dtau, asymmetry parameter g,                  */
/* single scattering albedo omega0, and cosine of solar zenith angle mu0 */

void eddington_v2 (double  dtau,
                   double  g,
                   double  omega0,
                   double  mu0,
                   double* t,
                   double* r,
                   double* rdir,
                   double* sdir,
                   double* tdir) {
  double alpha1 = 0, alpha2 = 0, alpha3 = 0, alpha4 = 0, alpha5 = 0, alpha6 = 0;
  double a11 = 0, a12 = 0, a13 = 0, a23 = 0, a33 = 0;
  double lambda = 0, b = 0, A = 0;
  double denom = 0;

  /* first, avoid omega0=1 because of instability */
  if (omega0 > 0.999999)
    omega0 = 0.999999;

  alpha1 = (1.0 - omega0) + 0.75 * (1.0 - omega0 * g);
  alpha2 = -(1.0 - omega0) + 0.75 * (1.0 - omega0 * g);

  lambda = sqrt (alpha1 * alpha1 - alpha2 * alpha2);

  A = 1.0 / (alpha2 / (alpha1 - lambda) * exp (lambda * dtau) - alpha2 / (alpha1 + lambda) * exp (-lambda * dtau));

  a11 = A * 2.0 * lambda / alpha2;
  a12 = A * (exp (lambda * dtau) - exp (-lambda * dtau));

  b      = 0.5 - 0.75 * g * mu0;
  alpha3 = -omega0 * b;
  alpha4 = omega0 * (1 - b);

  denom  = (1.0 / mu0 / mu0 - lambda * lambda);
  alpha5 = ((alpha1 - 1.0 / mu0) * alpha3 - alpha2 * alpha4) / denom;
  alpha6 = (alpha2 * alpha3 - (alpha1 + 1.0 / mu0) * alpha4) / denom;

  a33 = exp (-dtau / mu0);

  a13 = alpha5 * (1.0 - a11 * a33) - alpha6 * a12;
  a23 = -(a12) *alpha5 * a33 + alpha6 * (a33 - a11);

  *t    = a11;
  *r    = a12;
  *tdir = a33;
  *rdir = a13 / mu0;
  *sdir = a23 / mu0;
}

static int adding (double*  dtau_org,
                   double*  omega0_org,
                   double*  g_org,
                   int      nlev,
                   double   S0,
                   double   mu0,
                   double   Ag,
                   int      delta,
                   double** Edir,
                   double** Edn,
                   double** Eup,
                   double** Lavg) {
  double dtau = 0, omega0 = 0, g = 0;

  double *t = NULL, *r = NULL, *tdir = NULL, *rdir = NULL, *sdir = NULL, *tdir_unscaled = NULL;
  double *T = NULL, *R = NULL, *Tdir = NULL, *Sdir = NULL, *Tdir_unscaled = NULL;
  double* Edir_unscaled = NULL;

  /* allocate memory for reflectivity / transmission */
  t    = calloc (nlev - 1, sizeof (double));
  r    = calloc (nlev - 1, sizeof (double));
  tdir = calloc (nlev - 1, sizeof (double));
  rdir = calloc (nlev - 1, sizeof (double));
  sdir = calloc (nlev - 1, sizeof (double));

  T    = calloc (nlev - 1, sizeof (double));
  R    = calloc (nlev - 1, sizeof (double));
  Tdir = calloc (nlev - 1, sizeof (double));
  Sdir = calloc (nlev - 1, sizeof (double));

  tdir_unscaled = calloc (nlev - 1, sizeof (double));
  Tdir_unscaled = calloc (nlev - 1, sizeof (double));

  /* allocate memory for result arrays */
  *Edir = calloc (nlev, sizeof (double));
  *Edn  = calloc (nlev, sizeof (double));
  *Eup  = calloc (nlev, sizeof (double));
  *Lavg = calloc (nlev, sizeof (double));

  Edir_unscaled = calloc (nlev, sizeof (double));

  if (mu0 <= 0) { /* no need to do calculation if neither solar nor thermal source */
    for (int ilev = 0; ilev < nlev; ilev++) {
      (*Edir)[ilev] = 0;
      (*Edn)[ilev]  = 0;
      (*Eup)[ilev]  = 0;
      (*Lavg)[ilev] = 0;
    }
    return 0;
  }

  if (mu0 > 0)
    (*Edir)[0] = S0 * mu0;
  else
    (*Edir)[0] = 0;

  Edir_unscaled[0] = (*Edir)[0];

  for (int ilyr = 0; ilyr < nlev - 1; ilyr++) {

    if (delta)
      delta_scale_hg (dtau_org[ilyr], omega0_org[ilyr], g_org[ilyr], &dtau, &omega0, &g);
    else {
      dtau   = dtau_org[ilyr];
      omega0 = omega0_org[ilyr];
      g      = g_org[ilyr];
    }

    /* limit dtau to 100 to avoid numerical overflow */
    if (dtau > 100)
      dtau = 100;

    eddington_v2 (dtau, g, omega0, mu0, &(t[ilyr]), &(r[ilyr]), &(rdir[ilyr]), &(sdir[ilyr]), &(tdir[ilyr]));

    tdir_unscaled[ilyr] = exp (-dtau_org[ilyr] / mu0);
  }

  R[0]    = r[0];
  T[0]    = t[0];
  Tdir[0] = tdir[0];
  Sdir[0] = sdir[0];

  Tdir_unscaled[0] = tdir[0];

  for (int ilev = 0; ilev < nlev - 2; ilev++) {
    R[ilev + 1] = r[ilev + 1] + (R[ilev] * t[ilev + 1] * t[ilev + 1]) / (1.0 - R[ilev] * r[ilev + 1]);
    T[ilev + 1] = T[ilev] * t[ilev + 1] / (1.0 - R[ilev] * r[ilev + 1]);

    Tdir[ilev + 1] = Tdir[ilev] * tdir[ilev + 1];
    Sdir[ilev + 1] =
      (t[ilev + 1] * Sdir[ilev] + Tdir[ilev] * rdir[ilev + 1] * R[ilev] * t[ilev + 1]) / (1.0 - R[ilev] * r[ilev + 1]) +
      Tdir[ilev] * sdir[ilev + 1];

    Tdir_unscaled[ilev + 1] = Tdir_unscaled[ilev] * tdir_unscaled[ilev + 1];
  }

  for (int ilev = nlev - 1; ilev >= 1; ilev--) {
    (*Edir)[ilev]       = Tdir[ilev - 1] * (*Edir)[0];
    Edir_unscaled[ilev] = Tdir_unscaled[ilev - 1] * Edir_unscaled[0];
  }

  (*Edn)[nlev - 1] = (Sdir[nlev - 2] + Tdir[nlev - 2] * R[nlev - 2] * Ag) / (1.0 - R[nlev - 2] * Ag) * (*Edir)[0];
  (*Eup)[nlev - 1] = Ag * ((*Edn)[nlev - 1] + (*Edir)[nlev - 1]);

  for (int ilev = nlev - 1; ilev >= 2; ilev--) {

    (*Edn)[ilev - 1] =
      (R[ilev - 2] * t[ilev - 1] * (*Eup)[ilev] + (*Edir)[0] * Sdir[ilev - 2] + (*Edir)[ilev - 1] * rdir[ilev - 1] * R[ilev - 2]) /
      (1 - R[ilev - 2] * r[ilev - 1]);

    (*Eup)[ilev - 1] =
      (t[ilev - 1] * (*Eup)[ilev] + (*Edir)[0] * Sdir[ilev - 2] * r[ilev - 1] + (*Edir)[ilev - 1] * rdir[ilev - 1]) /
      (1 - R[ilev - 2] * r[ilev - 1]);
  }

  (*Eup)[0] = t[0] * (*Eup)[1] + rdir[0] * (*Edir)[0];

  /* replace direct irradiance with non-delta-scaled value and correct diffuse irradiance accordingly */
  for (int ilev = 0; ilev < nlev; ilev++) {
    (*Edn)[ilev] += (*Edir)[ilev] - Edir_unscaled[ilev];
    (*Edir)[ilev] = Edir_unscaled[ilev];
  }

  free (t);
  free (r);
  free (tdir);
  free (rdir);
  free (sdir);

  free (T);
  free (R);
  free (Tdir);
  free (Sdir);

  free (tdir_unscaled);
  free (Tdir_unscaled);

  return 0;
}
