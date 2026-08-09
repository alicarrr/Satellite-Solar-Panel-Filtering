
/*----------------------------------------------------------------------
Full solution of the multi-layer two-stream equation (solar, thermal)
with maximum-random overlap assumption for partial cloudiness.
Author: Nina Črnivec, nina.crnivec@physik.uni-muenchen.de
----------------------------------------------------------------------*/

/*-------------------------------------------------------------------
 * This file is part of libRadtran.
 * Copyright (c) 1997-2018 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras,
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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "allocnd.h"
#include "bandec.h"
#include "common_math.h"
#include "equation.h"
#include "errors.h"
#include "solver.h"
#include "twomaxrnd.h"
#include "twostrebe.h"

#define NFLUX 4  // number of diffuse fluxes defined at each level (Edn_c, Edn_f, Eup_c, Eup_f)

//======================
// FUNCTION calcp1p2p3p4
//======================
// INPUT:  nlev, cf;
// OUTPUT: ar_p1, ar_p2, ar_p3, ar_p4;

/*
* From Zdunkowski (pages 180-183):
* p1 = (1 - max(cf[i], cf[i-1])) / (1 - cf[i-1])
* p2 = (1 - max(cf[i], cf[i+1]) / (1 - cf[i+1])
* p3 = min(cf[i], cf[i-1]) / cf[i-1]
* p4 = min(cf[i], cf[i+1]) / cf[i+1]

* Remark from Zdunkowski:
* It is noteworthy that a particular coefficient pj (j=1,2,3,4) is set equal to 1
* if an undetermined expression 0/0 occurs.
* This follows from physical reasoning or from applying l'Hopital's rule.
*/

static int calcp1p2p3p4 (int nlev, double* cf, double* ar_p1, double* ar_p2, double* ar_p3, double* ar_p4) {

  const int nlyr = nlev - 1;

  /*
  * FORMULAS:
  * ar_p1[ilyr] = (1.0 - max(cf[ilyr], cf[ilyr-1])) / (1.0 - cf[ilyr-1]);
  * ar_p2[ilyr] = (1.0 - max(cf[ilyr], cf[ilyr+1]) / (1.0 - cf[ilyr+1]);
  * ar_p3[ilyr] = min(cf[ilyr], cf[ilyr-1]) / cf[ilyr-1];
  * ar_p4[ilyr] = min(cf[ilyr], cf[ilyr+1]) / cf[ilyr+1];
  */

  // Calculate vertical profiles of p1, p3:
  // Special case: ilyr=0;
  for (int ilyr = 0; ilyr < nlyr; ilyr++) {
    if (ilyr == 0) {
      ar_p1[ilyr] = 1.0;
      ar_p3[ilyr] = 1.0;
    } else {

      double max_p1;  // = max(cf[ilyr],cf[ilyr-1]); appears in the expression for p1;
      // Find max_p1:
      if (cf[ilyr] > cf[ilyr - 1]) {
        max_p1 = cf[ilyr];
      } else {
        max_p1 = cf[ilyr - 1];
      }  //e-if

      // Find_min_p3:
      double min_p3;  // = min(cf[ilyr],cf[ilyr-1]); appears in the expression for p3;
      if (cf[ilyr] < cf[ilyr - 1]) {
        min_p3 = cf[ilyr];
      } else {
        min_p3 = cf[ilyr - 1];
      }  //e-if

      if (cf[ilyr - 1] == 1.0)
        ar_p1[ilyr] = 1.0;
      else
        ar_p1[ilyr] = (1.0 - max_p1) / (1.0 - cf[ilyr - 1]);

      if (cf[ilyr - 1] == 0.0)
        ar_p3[ilyr] = 1.0;
      else
        ar_p3[ilyr] = min_p3 / cf[ilyr - 1];

    }  //e-if
  }    //e-for over ilyr;

  // Calculate vertical profiles of p2 and p4:
  // Special case: ilyr = nlyr-1 = nlev-2;
  for (int ilyr = 0; ilyr < nlyr; ilyr++) {

    if (ilyr == (nlyr - 1)) {
      ar_p2[ilyr] = 1.0;
      ar_p4[ilyr] = 1.0;
    } else {

      // Find max_p2:
      double max_p2;  // = max(cf[ilyr],cf[ilyr+1]); appears in the expression for p2;
      if (cf[ilyr] > cf[ilyr + 1]) {
        max_p2 = cf[ilyr];
      } else {
        max_p2 = cf[ilyr + 1];
      }  //e-if

      // Find_min_p4:
      double min_p4;  // = min(cf[ilyr],cf[ilyr+1]); appears in the expression for p4;
      if (cf[ilyr] < cf[ilyr + 1]) {
        min_p4 = cf[ilyr];
      } else {
        min_p4 = cf[ilyr + 1];
      }  //e-if

      if (cf[ilyr + 1] == 1.0)
        ar_p2[ilyr] = 1.0;
      else
        ar_p2[ilyr] = (1.0 - max_p2) / (1.0 - cf[ilyr + 1]);

      if (cf[ilyr + 1] == 0.0)
        ar_p4[ilyr] = 1.0;
      else
        ar_p4[ilyr] = min_p4 / cf[ilyr + 1];
    }  //e-if
  }    //e-for

  return 0;
}  //e-calcp1p2p3p4

//===========================
// FUNCTION: eddington_coeffc
//===========================

/*
Calculate Eddington coefficients a11, a12, a13, a23, and a33
from layer optical thickness dtau, asymmetry parameter g,
single scattering albedo omega0, and cosine of solar zenith angle mu0.
*/

void eddington_coeffc (double  dtau,
                       double  g,
                       double  omega0,
                       double  mu0,
                       double* a11,
                       double* a12,
                       double* a13,
                       double* a23,
                       double* a33) {
  double alpha1 = 1.75;
  double alpha2 = -0.25;

  if (omega0 != 0) {
    alpha1 = (1.0 - omega0) + 0.75 * (1.0 - omega0 * g);
    alpha2 = -(1.0 - omega0) + 0.75 * (1.0 - omega0 * g);
  }  //e-if

  const double lambda = sqrt (alpha1 * alpha1 - alpha2 * alpha2);

  const double A = 1.0 / (alpha2 / (alpha1 - lambda) * exp (lambda * dtau) - alpha2 / (alpha1 + lambda) * exp (-lambda * dtau));

  *a11 = A * 2.0 * lambda / alpha2;
  *a12 = A * (exp (lambda * dtau) - exp (-lambda * dtau));

  const double b      = 0.5 - 0.75 * g * mu0;
  double       alpha3 = -omega0 * b;
  double       alpha4 = omega0 * (1 - b);

  const double denom  = (1.0 / mu0 / mu0 - lambda * lambda);
  const double alpha5 = ((alpha1 - 1.0 / mu0) * alpha3 - alpha2 * alpha4) / denom;
  const double alpha6 = (alpha2 * alpha3 - (alpha1 + 1.0 / mu0) * alpha4) / denom;

  *a33 = exp (-dtau / mu0);

  *a13 = alpha5 * (1.0 - (*a11) * (*a33)) - alpha6 * (*a12);
  *a23 = -(*a12) * alpha5 * (*a33) + alpha6 * ((*a33) - (*a11));
}  //e-eddington_coeffc

//================================
// FUNCTION: calcThermalComponents
//================================
// INPUT:  ilyr, B, dtau, omega0, g;
// OUTPUT: theComp1, theComp2;

// theComp1 = Planckian emission term of a given layer (upward direction);
// theComp2 = Planckian emission term of a given layer (downward direction);

void calcThermalComponents (int ilyr, double* B, double dtau, double omega0, double g, double* theComp1, double* theComp2) {
  double B0;
  double B1;
  double alpha1;
  double alpha2;
  double lambda;
  double A;
  double a11;
  double a12;
  double kappa;
  double alpha7;
  double alpha8;
  double alpha9;
  double alpha10;

  // Calculate B0, B1:
  if (dtau > 0.01) {
    B1 = (B[ilyr + 1] - B[ilyr]) / dtau;
    B0 = B[ilyr];
  } else {
    B1 = 0;
    B0 = (B[ilyr + 1] + B[ilyr]) * 0.5;
  }  //e-if

  if (omega0 == 0) {
    alpha1 = 1.75;
    alpha2 = -0.25;
  } else {
    alpha1 = (1.0 - omega0) + 0.75 * (1.0 - omega0 * g);
    alpha2 = -(1.0 - omega0) + 0.75 * (1.0 - omega0 * g);
  }  //e-if

  lambda = sqrt (alpha1 * alpha1 - alpha2 * alpha2);

  A = 1.0 / (alpha2 / (alpha1 - lambda) * exp (lambda * dtau) - alpha2 / (alpha1 + lambda) * exp (-lambda * dtau));

  a11 = A * 2.0 * lambda / alpha2;
  a12 = A * (exp (lambda * dtau) - exp (-lambda * dtau));

  kappa  = 2.0 * M_PI * (1.0 - omega0);
  alpha7 = kappa * (B0 / (alpha1 - alpha2) + B1 / lambda / lambda);
  alpha9 = kappa * (B0 / (alpha1 - alpha2) - B1 / lambda / lambda);

  alpha8  = kappa * B1 / (alpha1 - alpha2);
  alpha10 = alpha8;

  *theComp1 = -a11 * (alpha7 + alpha8 * dtau) - a12 * (alpha9) + alpha7;
  *theComp2 = -a12 * (alpha7 + alpha8 * dtau) - a11 * (alpha9) + alpha9 + alpha10 * dtau;

}  //e-calcThermalComponents

//=======================
// FUNCTION: buildMatrixA
//=======================
// INPUT:  nlev, Ag, ar_a11_c, ar_a11_f, ar_a12_c, ar_a12_f, ar_p1, ar_p2, ar_p3, ar_p4;
// OUTPUT: matrixA;

static int buildMatrixA (int      nlev,
                         double   Ag,
                         double*  ar_a11_c,
                         double*  ar_a11_f,
                         double*  ar_a12_c,
                         double*  ar_a12_f,
                         double*  ar_p1,
                         double*  ar_p2,
                         double*  ar_p3,
                         double*  ar_p4,
                         double** matrixA) {
  // Set values for initial four rows of matrix A:
  // First row:
  matrixA[0][2] = ar_a12_f[0] * ar_p1[0];
  matrixA[0][3] = ar_a12_f[0] * (1.0 - ar_p3[0]);
  matrixA[0][4] = ar_a11_f[0] * ar_p2[0];
  matrixA[0][5] = ar_a11_f[0] * (1.0 - ar_p4[0]);

  // Second row:
  matrixA[1][2] = ar_a12_c[0] * (1.0 - ar_p1[0]);
  matrixA[1][3] = ar_a12_c[0] * ar_p3[0];
  matrixA[1][4] = ar_a11_c[0] * (1.0 - ar_p2[0]);
  matrixA[1][5] = ar_a11_c[0] * ar_p4[0];
  // Third row is already zero; (needs to be zero due to upper boundary condition);
  // Forth row is already zero; (needs to be zero due to upper boundary condition);

  int nextColEup = 6;  // index of column for Eup(level1)
  int nextColEdw = 2;  // index of column for Edn(level1)

  for (int i = 1; i < nlev - 1; i++) {
    const int iRow = NFLUX * i;

    matrixA[iRow][nextColEup]     = ar_a12_f[i] * ar_p1[i];
    matrixA[iRow][nextColEup + 1] = ar_a12_f[i] * (1.0 - ar_p3[i]);
    matrixA[iRow][nextColEup + 2] = ar_a11_f[i] * ar_p2[i];
    matrixA[iRow][nextColEup + 3] = ar_a11_f[i] * (1.0 - ar_p4[i]);

    matrixA[iRow + 1][nextColEup]     = ar_a12_c[i] * (1.0 - ar_p1[i]);
    matrixA[iRow + 1][nextColEup + 1] = ar_a12_c[i] * ar_p3[i];
    matrixA[iRow + 1][nextColEup + 2] = ar_a11_c[i] * (1.0 - ar_p2[i]);
    matrixA[iRow + 1][nextColEup + 3] = ar_a11_c[i] * ar_p4[i];
    nextColEup += NFLUX;

    matrixA[iRow + 2][nextColEdw]     = ar_a11_f[i - 1] * ar_p1[i - 1];
    matrixA[iRow + 2][nextColEdw + 1] = ar_a11_f[i - 1] * (1.0 - ar_p3[i - 1]);
    matrixA[iRow + 2][nextColEdw + 2] = ar_a12_f[i - 1] * ar_p2[i - 1];
    matrixA[iRow + 2][nextColEdw + 3] = ar_a12_f[i - 1] * (1.0 - ar_p4[i - 1]);

    matrixA[iRow + 3][nextColEdw]     = ar_a11_c[i - 1] * (1.0 - ar_p1[i - 1]);
    matrixA[iRow + 3][nextColEdw + 1] = ar_a11_c[i - 1] * ar_p3[i - 1];
    matrixA[iRow + 3][nextColEdw + 2] = ar_a12_c[i - 1] * (1.0 - ar_p2[i - 1]);
    matrixA[iRow + 3][nextColEdw + 3] = ar_a12_c[i - 1] * ar_p4[i - 1];

    nextColEdw += NFLUX;
  }  //e-for

  {  // lower boundary condition for the forth and third row from bottom up of matrix A;
    const int i                       = (nlev - 1);
    const int iRow                    = NFLUX * i;
    matrixA[iRow][nextColEup]         = Ag;
    matrixA[iRow + 1][nextColEup + 1] = Ag;

    matrixA[iRow + 2][nextColEdw]     = ar_a11_f[i - 1] * ar_p1[i - 1];
    matrixA[iRow + 2][nextColEdw + 1] = ar_a11_f[i - 1] * (1.0 - ar_p3[i - 1]);
    matrixA[iRow + 2][nextColEdw + 2] = ar_a12_f[i - 1] * ar_p2[i - 1];
    matrixA[iRow + 2][nextColEdw + 3] = ar_a12_f[i - 1] * (1.0 - ar_p4[i - 1]);

    matrixA[iRow + 3][nextColEdw]     = ar_a11_c[i - 1] * (1.0 - ar_p1[i - 1]);
    matrixA[iRow + 3][nextColEdw + 1] = ar_a11_c[i - 1] * ar_p3[i - 1];
    matrixA[iRow + 3][nextColEdw + 2] = ar_a12_c[i - 1] * (1.0 - ar_p2[i - 1]);
    matrixA[iRow + 3][nextColEdw + 3] = ar_a12_c[i - 1] * ar_p4[i - 1];
  }

  return 0;
}  //e-buildMatrixA

//=========================
// FUNCTION buildVectorBsol
//=========================
// INPUT:  nlev, Ag, mu0, ar_a13_c, ar_a13_f, ar_a23_c, ar_a23_f, ar_S_c, ar_S_f, ar_p1, ar_p3;
// OUTPUT: vectB;

static int buildVectorBsol (int     nlev,
                            double  Ag,
                            double  mu0,
                            double* ar_a13_c,
                            double* ar_a13_f,
                            double* ar_a23_c,
                            double* ar_a23_f,
                            double* ar_S_c,
                            double* ar_S_f,
                            double* ar_p1,
                            double* ar_p3,
                            double* vectB) {
  // Set initial four values:
  vectB[0] = ar_a13_f[0] * ar_S_f[0];
  vectB[1] = ar_a13_c[0] * ar_S_c[0];
  vectB[2] = 0.0;  // upper boundary condition
  vectB[3] = 0.0;  // upper boundary condition

  for (int i = 1; i < (nlev - 1); i++) {
    const int j  = NFLUX * i;
    vectB[j]     = ar_a13_f[i] * (ar_p1[i] * ar_S_f[i] + (1.0 - ar_p3[i]) * ar_S_c[i]);
    vectB[j + 1] = ar_a13_c[i] * ((1.0 - ar_p1[i]) * ar_S_f[i] + ar_p3[i] * ar_S_c[i]);
    vectB[j + 2] = ar_a23_f[i - 1] * (ar_p1[i - 1] * ar_S_f[i - 1] + (1.0 - ar_p3[i - 1]) * ar_S_c[i - 1]);
    vectB[j + 3] = ar_a23_c[i - 1] * ((1.0 - ar_p1[i - 1]) * ar_S_f[i - 1] + ar_p3[i - 1] * ar_S_c[i - 1]);
  }  //e-for

  // Treat last four values seperately:
  const int j  = NFLUX * (nlev - 1);
  vectB[j]     = Ag * mu0 * ar_S_f[nlev - 1];  // lower boundary condition
  vectB[j + 1] = Ag * mu0 * ar_S_c[nlev - 1];  // lower boundary condition
  vectB[j + 2] = ar_a23_f[nlev - 2] * (ar_p1[nlev - 2] * ar_S_f[nlev - 2] + (1.0 - ar_p3[nlev - 2]) * ar_S_c[nlev - 2]);
  vectB[j + 3] = ar_a23_c[nlev - 2] * ((1.0 - ar_p1[nlev - 2]) * ar_S_f[nlev - 2] + ar_p3[nlev - 2] * ar_S_c[nlev - 2]);

  return 0;
}  //e-buildVectorBsol

//=========================
// FUNCTION buildVectorBthe
//=========================
// INPUT:  nlev, Ag, Bg, ar_theComp1_c, ar_theComp1_f, ar_theComp2_c, ar_theComp2_f, cf;
// OUTPUT: bb_the;

// thermal component 1 = upward emission
// thermal component 2 = downward emission

static int buildVectorBthe (int     nlev,
                            double  Ag,
                            double  Bg,
                            double* ar_theComp1_c,
                            double* ar_theComp1_f,
                            double* ar_theComp2_c,
                            double* ar_theComp2_f,
                            double* cf,
                            double* vectB) {
  // Set initial four values:
  vectB[0] = (1.0 - cf[0]) * ar_theComp1_f[0];
  vectB[1] = cf[0] * ar_theComp1_c[0];
  vectB[2] = 0.0;  // upper boundary condition
  vectB[3] = 0.0;  // upper boundary condition

  for (int i = 1; i < (nlev - 1); i++) {
    const int j  = NFLUX * i;
    vectB[j]     = (1.0 - cf[i]) * ar_theComp1_f[i];
    vectB[j + 1] = cf[i] * ar_theComp1_c[i];
    vectB[j + 2] = (1.0 - cf[i - 1]) * ar_theComp2_f[i - 1];
    vectB[j + 3] = cf[i - 1] * ar_theComp2_c[i - 1];
  }  //e-for

  // Treat last four values separately:
  // i=nlev-1; // bottom (ground) level;
  const int j  = NFLUX * (nlev - 1);
  vectB[j]     = (1.0 - cf[nlev - 2]) * (1.0 - Ag) * M_PI * Bg;  // lower boundary condition
  vectB[j + 1] = cf[nlev - 2] * (1.0 - Ag) * M_PI * Bg;          // lower boundary condition
  vectB[j + 2] = (1.0 - cf[nlev - 2]) * ar_theComp2_f[nlev - 2];
  vectB[j + 3] = cf[nlev - 2] * ar_theComp2_c[nlev - 2];

  return 0;
}  //e-buildVectorBthe

//===========================
// FUNCTION twostream_maxrand
//===========================

/*
OUTPUT:
Edir = direct irradiance;
Edn  = diffuse downward irradiance;
Eup  = diffuse upward irradiance;
Lavg = radiance (set to NaN);
*/

static int twostream_maxrand (double* dtau_org_c,
                              double* omega0_org_c,
                              double* g_org_c,  //"_c" for cloudy region
                              double* dtau_org_f,
                              double* omega0_org_f,
                              double* g_org_f,  //"_f" for cloud-free region
                              double* cf,
                              int     nlev,
                              double  S0,
                              double  mu0,
                              double  Ag,
                              double  Bg,
                              double* B,
                              int     delta,
                              int     flagSolar,
                              int     flagThermal,
                              double* Edir,
                              double* Edn,
                              double* Eup,
                              double* Lavg) {
  const int nlyr = nlev - 1;
  int       status;

  /*
  Eddington coefficients:
  a11 = transmission coefficient for diffuse radiation;
  a12 = reflection coefficient for diffuse radiation;
  a13 = reflection coefficient for the primary scattered parallel solar radiation;
  a23 = transmission coefficient for the primary scattered parallel solar radiation;
  a33 = transmission coefficient for the direct parallel solar radiation;
  */

  double a11_c[nlyr];
  double a12_c[nlyr];
  double a13_c[nlyr];
  double a23_c[nlyr];
  double a33_c[nlyr];

  double a11_f[nlyr];
  double a12_f[nlyr];
  double a13_f[nlyr];
  double a23_f[nlyr];
  double a33_f[nlyr];

  double theComp1_c[nlev];
  double theComp1_f[nlev];
  double theComp2_c[nlev];
  double theComp2_f[nlev];

  // Parameters related to cloud cover of two adjacent layers:
  // in Zdunkowski (pages: 180-183) denoted as: b1,b2,b3,b4;
  // here: p1, p2, p3, p4;
  double p1[nlyr];
  double p2[nlyr];
  double p3[nlyr];
  double p4[nlyr];

  double bb_sol[NFLUX * nlev];
  double bb_the[NFLUX * nlev];
  double bb[NFLUX * nlev];
  double xx[NFLUX * nlev];

  double S_c[nlev];
  double S_f[nlev];

  double Edir_c[nlev];
  double Edir_f[nlev];

  double Eup_c[nlev];
  double Eup_f[nlev];

  double Edn_c[nlev];
  double Edn_f[nlev];

  double** AA = calloc_double_2D (NFLUX * nlev, NFLUX * nlev, "AA");
  CHKPOINTER (AA);

  // At the moment it is only possible to calculate solar OR thermal RT, but not both simultaneously;
  if ((flagSolar && flagThermal) || (!flagSolar && !flagThermal)) {
    CHKERROUT (-2, "Error - invalid input parameters - use flagSolar and flagThermal alternatingly");
  }  //e-if

  // Calculate vertical profiles of p1, p2, p3 and p4 from vertical profile of cloud fraction:
  status = calcp1p2p3p4 (nlev, cf, p1, p2, p3, p4);
  CHKERROUT (status, "Error calculating vertical profiles of p1, p2, p3 and p4 from cloud cover");

  // Calculate vertical profiles of Eddington coefficients
  // and vertical profiles of thermal components:
  for (int ilyr = 0; ilyr < nlyr; ilyr++) {

    double dtau_c   = dtau_org_c[ilyr];
    double omega0_c = omega0_org_c[ilyr];
    double g_c      = g_org_c[ilyr];

    // Delta scaling of optical properties of cloudy regions:
    if (delta) {
      delta_scale_hg (dtau_org_c[ilyr], omega0_org_c[ilyr], g_org_c[ilyr], &dtau_c, &omega0_c, &g_c);
    }  //e-if
    // omega0 should not be 1 (avoiding singularity problem), restrict omega0 to 0.999999:
    omega0_c = MIN (omega0_c, 0.999999);

    // Optical properties of cloud-free regions:
    const double dtau_f   = dtau_org_f[ilyr];
    const double omega0_f = MIN (omega0_org_f[ilyr], 0.999999);
    const double g_f      = g_org_f[ilyr];

    eddington_coeffc (dtau_c, g_c, omega0_c, mu0, &(a11_c[ilyr]), &(a12_c[ilyr]), &(a13_c[ilyr]), &(a23_c[ilyr]), &(a33_c[ilyr]));
    eddington_coeffc (dtau_f, g_f, omega0_f, mu0, &(a11_f[ilyr]), &(a12_f[ilyr]), &(a13_f[ilyr]), &(a23_f[ilyr]), &(a33_f[ilyr]));

    if (flagThermal) {  // Components of vector B in the thermal spectral range:
      calcThermalComponents (ilyr, B, dtau_c, omega0_c, g_c, &(theComp1_c[ilyr]), &(theComp2_c[ilyr]));
      calcThermalComponents (ilyr, B, dtau_f, omega0_f, g_f, &(theComp1_f[ilyr]), &(theComp2_f[ilyr]));
    }  //e-if
  }    //e-for

  // Initialize vectors S_c and S_f:
  // S[0]= S0 = S_c[0] + S_f[0] = cf[0]*S0 + (1.0-cf[0])*S0;
  if (flagSolar) {
    S_c[0] = cf[0] * S0;
    S_f[0] = (1.0 - cf[0]) * S0;

    for (int ilev = 1; ilev < nlev; ilev++) {
      S_c[ilev] = a33_c[ilev - 1] * ((1.0 - p1[ilev - 1]) * S_f[ilev - 1] + p3[ilev - 1] * S_c[ilev - 1]);
      S_f[ilev] = a33_f[ilev - 1] * (p1[ilev - 1] * S_f[ilev - 1] + (1.0 - p3[ilev - 1]) * S_c[ilev - 1]);
    }  //e-for
  }    //e-if

  // Equation system has the following form: xx = AA*xx + bb;

  // Build matrix AA:
  status = buildMatrixA (nlev, Ag, a11_c, a11_f, a12_c, a12_f, p1, p2, p3, p4, AA);
  CHKERR (status);

  // Make matrix A2 = AA - IdentityMatrix and save it to the same matrix AA,
  // since the equation system in the form Ax=b must be passed to solve_gauss;
  // xx = AA*xx + bb;
  // (AA-II)*xx = -bb;
  for (int iRow = 0; iRow < nlev * NFLUX; iRow++)
    AA[iRow][iRow] = -1.0;

  // Build vector bb_sol:
  if (flagSolar) {
    status = buildVectorBsol (nlev, Ag, mu0, a13_c, a13_f, a23_c, a23_f, S_c, S_f, p1, p3, bb_sol);
    CHKERR (status);
  }  //e-if

  // Build vector bb_the:
  if (flagThermal) {
    status = buildVectorBthe (nlev, Ag, Bg, theComp1_c, theComp1_f, theComp2_c, theComp2_f, cf, bb_the);
    CHKERR (status);
  }  //e-if

  // Assign bb_sol or bb_the to the final vector bb:
  for (int ilev = 0; ilev < NFLUX * nlev; ilev++) {
    bb[ilev] = 0;

    if (flagSolar) {
      bb[ilev] += bb_sol[ilev];
    }  //e-if

    if (flagThermal) {
      bb[ilev] += bb_the[ilev];
    }  //e-if

  }  //e-for

  // Make vector -bb and save it to the same vector bb:
  status = vec_scale (nlev * NFLUX, -1, bb);
  CHKERR (status);

  // Solve AA*xx=bb to obtain xx:

  const long nn = NFLUX * nlev;
  const long m1 = 5;
  const long m2 = 5;  // n=2*input.nlyr+2;

  double** matr = malloc_double_2D (nn, m1 + m2 + 1, "matr");
  double** mlu  = malloc_double_2D (nn, m2, "mlu");
  long     indx[nn];

  matrix_to_band_form (AA, nn, m1, m2, matr);  // Transform matrix AA to band form matr

  free_double_2D (AA);

  bandec_M (matr, nn, m1, m2, mlu, indx);
  banbks_M (matr, nn, m1, m2, mlu, indx, bb);

  for (int i = 0; i < nn; i++) {
    xx[i] = bb[i];
  }

  free_double_2D (matr);
  free_double_2D (mlu);

  for (int ilev = 0; ilev < nlev; ilev++) {
    Eup_f[ilev] = xx[NFLUX * ilev];
    Eup_c[ilev] = xx[NFLUX * ilev + 1];
    Edn_f[ilev] = xx[NFLUX * ilev + 2];
    Edn_c[ilev] = xx[NFLUX * ilev + 3];
  }  //e-for

  for (int ilev = 0; ilev < nlev; ilev++) {
    Edir_c[ilev] = S_c[ilev] * mu0;
    Edir_f[ilev] = S_f[ilev] * mu0;
  }  //e-for

  // Sum up the irradiances for cloudy and cloud-free regions to obtain the final result:
  for (int ilev = 0; ilev < nlev; ilev++) {
    Edir[ilev] = Edir_c[ilev] + Edir_f[ilev];
    Eup[ilev]  = Eup_c[ilev] + Eup_f[ilev];
    Edn[ilev]  = Edn_c[ilev] + Edn_f[ilev];
    Lavg[ilev] = NAN;
  }  //e-for

  return 0;
}  //e-twostream_maxrand

int twomaxrnd (float* dtau_org,
               float* omega0_org,
               float* g_org,
               float* dtau_clr_org,
               float* omega0_clr_org,
               float* g_clr_org,
               float* cf,
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

  double dtau_clr_org_d[nlev - 1];
  double omega0_clr_org_d[nlev - 1];
  double g_clr_org_d[nlev - 1];

  double dtau_org_d[nlev - 1];
  double omega0_org_d[nlev - 1];
  double g_org_d[nlev - 1];

  double cf_d[nlev - 1];
  double B[nlev];
  double Edir[nlev];
  double Edn[nlev];
  double Eup[nlev];
  double Lavg[nlev];

  const double taumax = 100.0;
  double       Bg     = 0.0;  // Used only for thermal RT, unused (0.0) for solar RT;

  const int flagSolar = S0 > 0.0;

  if (planck) {
    float plkavg;
    for (int ilev = 0; ilev < nlev; ilev++) {  // level temperatures and Planck functions
      F77_FUNC (cplkavg, CPLKAVG) (&wvnmlo, &wvnmhi, &temper[ilev], &plkavg);
      B[ilev] = plkavg;
    }  //e-for

    // surface temperature and Planck function
    F77_FUNC (cplkavg, CPLKAVG) (&wvnmlo, &wvnmhi, &btemp, &plkavg);
    Bg = plkavg;
  }  //e-if

  // copy float arrays to double arrays
  for (int ilyr = 0; ilyr < nlev - 1; ilyr++) {
    dtau_org_d[ilyr]     = dtau_org[ilyr];
    dtau_clr_org_d[ilyr] = dtau_clr_org[ilyr];

    // restrict layer optical thickness to 100
    if (dtau_org_d[ilyr] > taumax)
      dtau_org_d[ilyr] = taumax;

    if (dtau_clr_org_d[ilyr] > taumax)
      dtau_clr_org_d[ilyr] = taumax;

    omega0_org_d[ilyr]     = omega0_org[ilyr];
    omega0_clr_org_d[ilyr] = omega0_clr_org[ilyr];

    g_org_d[ilyr]     = g_org[ilyr];
    g_clr_org_d[ilyr] = g_clr_org[ilyr];

    cf_d[ilyr] = cf[ilyr];
  }  //e-for

  // Call twostream_maxrand code:
  int status = twostream_maxrand (dtau_org_d,
                                  omega0_org_d,
                                  g_org_d,
                                  dtau_clr_org_d,
                                  omega0_clr_org_d,
                                  g_clr_org_d,
                                  cf_d,
                                  nlev,
                                  S0,
                                  mu0,
                                  Ag,
                                  Bg,
                                  B,
                                  delta,
                                  flagSolar,
                                  planck,  // flagThermal = planck
                                  Edir,
                                  Edn,
                                  Eup,
                                  Lavg);
  CHKERR (status);

  // Copy results to final fields:
  for (int ilev = 0; ilev < nlev; ilev++) {
    for (int lu = 0; lu < nzout; lu++) {
      if (zout[lu] == zd[ilev]) {
        fldn[lu]  = Edn[ilev];
        flup[lu]  = Eup[ilev];
        fldir[lu] = Edir[ilev];
        uavg[lu]  = Lavg[ilev];
      }  //e-if
    }    //e-for
  }      //e-for

  return 0;
}
