
/*----------------------------------------------------------------------
Full solution of the multi-layer twostream equation (solar, thermal),
where each vertical layer is divided into three regions ("Tripleclouds"):
i)   one region of optically thick cloud ("ck"),
ii)  one region of optically thin cloud  ("cn"),
iii) one cloud-free region ("f");
The maximum-random overlap is assumed for regions of optically thick cloud
as well as for cloudy regions as a whole (i.e., ck + cn regions).
Author: Nina Črnivec, nina.crnivec@physik.uni-muenchen.de
----------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 * $Id: twomaxrnd3C.c 2623 2011-12-23 10:52:38Z bernhard.mayer $
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
#include "twomaxrnd3C.h"
#include "twomaxrnd.h"
#include "twostrebe.h"

#define NFLUX 6  // number of diffuse fluxes defined at each level (Edn_ck, Edn_cn, Edn_f, Eup_ck, Eup_cn, Eup_f)

//=========================================
//FUNCTION: calc_pdn_transport_coefficients
//=========================================
//INPUT:  nlev, cf, cf_ck, cf_cn;
//OUTPUT: pdn_CKtCK, pdn_CKtCN, pdn_CKtF, pdn_CNtCK, pdn_CNtCN, pdn_CNtF, pdn_FtCK, pdn_FtCN, pdn_FtF;

static int calc_pdn_transport_coefficients_twomaxrnd3C (int     nlev,
                                                        double* cf,
                                                        double* cf_ck,
                                                        double* cf_cn,
                                                        double* pdn_CKtCK,
                                                        double* pdn_CKtCN,
                                                        double* pdn_CKtF,
                                                        double* pdn_CNtCK,
                                                        double* pdn_CNtCN,
                                                        double* pdn_CNtF,
                                                        double* pdn_FtCK,
                                                        double* pdn_FtCN,
                                                        double* pdn_FtF) {
  const int nlyr = nlev - 1;
  {
    int ilyr = 0;  //SPECIAL CASE

    pdn_CKtCK[ilyr] = 1.0;
    pdn_CKtCN[ilyr] = 0.0;
    pdn_CKtF[ilyr]  = 0.0;

    pdn_CNtCK[ilyr] = 0.0;
    pdn_CNtCN[ilyr] = 1.0;
    pdn_CNtF[ilyr]  = 0.0;

    pdn_FtCK[ilyr] = 0.0;
    pdn_FtCN[ilyr] = 0.0;
    pdn_FtF[ilyr]  = 1.0;
  }

  for (int ilyr = 1; ilyr < nlyr; ilyr++) {
    const double a = cf_ck[ilyr - 1];     //upper-layer cloud thick;
    const double b = cf_cn[ilyr - 1];     //upper-layer cloud thin;
    const double c = 1.0 - cf[ilyr - 1];  //upper-layer clear-sky;
    const double x = cf_ck[ilyr];         //this-layer cloud thick;
    const double y = cf_cn[ilyr];         //this-layer cloud thin;
    const double z = 1.0 - cf[ilyr];      //this-layer clear-sky;

    if (cf[ilyr - 1] < 0.0000000001) {  //SPECIAL CASE: In the layer above entirely clear-sky (i.e., a=0, b=0, c=1)
      pdn_CKtCK[ilyr] = 1.0;
      pdn_CKtCN[ilyr] = 0.0;
      pdn_CKtF[ilyr]  = 0.0;

      pdn_CNtCK[ilyr] = 0.0;
      pdn_CNtCN[ilyr] = 1.0;
      pdn_CNtF[ilyr]  = 0.0;

      pdn_FtCK[ilyr] = x;
      pdn_FtCN[ilyr] = y;
      pdn_FtF[ilyr]  = z;

      continue;
    }  //-eif

    if (cf[ilyr] < 0.0000000001) {  // SPECIAL CASE: In this layer entirely clear-sky (i.e., x=0, y=0, z=1);
      pdn_CKtCK[ilyr] = 0.0;
      pdn_CKtCN[ilyr] = 0.0;
      pdn_CKtF[ilyr]  = 1.0;

      pdn_CNtCK[ilyr] = 0.0;
      pdn_CNtCN[ilyr] = 0.0;
      pdn_CNtF[ilyr]  = 1.0;

      pdn_FtCK[ilyr] = 0.0;  //x
      pdn_FtCN[ilyr] = 0.0;  //y
      pdn_FtF[ilyr]  = 1.0;  //z

      continue;
    }  //e-if

    if (cf[ilyr] > cf[ilyr - 1]) {  // Case 1: Cloud in this layer larger than cloud above;
                                    //       = Thick cloud in this layer larger than thick cloud above;
                                    // (implies also: c not equal to 0.0, but: c > 0.0);

      if (x < (a + b)) {  // Case 1a: Thick cloud in this layer smaller than the entire cloud above;
        pdn_CKtCK[ilyr] = 1.0;
        pdn_CKtCN[ilyr] = 0.0;
        pdn_CKtF[ilyr]  = 0.0;

        if (b > 0.0) {
          pdn_CNtCK[ilyr] = (x - a) / b;
          pdn_CNtCN[ilyr] = (a + b - x) / b;
        } else {
          fprintf (stderr, "pdn, Case 1a: This should NOT happen: b = 0.0 !!!, ilyr = %d \n", ilyr);
          return -101;
        }  //e-if
        pdn_CNtF[ilyr] = 0.0;

        pdn_FtCK[ilyr] = 0.0;
        if (c > 0.0) {
          pdn_FtCN[ilyr] = (x + y - a - b) / c;
          pdn_FtF[ilyr]  = z / c;
        } else {
          fprintf (stderr, "pdn, Case 1a: This should NOT happen: c = 0.0 !!!, ilyr = %d \n", ilyr);
          return -102;
        }  //e-if

      } else {  // Case 1b: Thick cloud in this layer larger than or equal to the entire cloud above;
        pdn_CKtCK[ilyr] = 1.0;
        pdn_CKtCN[ilyr] = 0.0;
        pdn_CKtF[ilyr]  = 0.0;

        pdn_CNtCK[ilyr] = 1.0;
        pdn_CNtCN[ilyr] = 0.0;
        pdn_CNtF[ilyr]  = 0.0;

        if (c > 0.0) {
          pdn_FtCK[ilyr] = (x - a - b) / c;
          pdn_FtCN[ilyr] = y / c;
          pdn_FtF[ilyr]  = z / c;
        } else {
          fprintf (stderr, "pdn, Case 1b: This should NOT happen: c = 0.0 !!!, ilyr = %d \n", ilyr);
          return -103;
        }  //e-if

      }  //e-if-Case1a1b

    } else {  // Case 2: Cloud in this layer smaller than or equal to the cloud above;
              //       = Thick cloud in this layer smaller than or equal to the thick cloud above;

      if (a < (x + y)) {  // Case 2a: Thick cloud above smaller than the entire cloud in this layer;
        if (a > 0.0) {
          pdn_CKtCK[ilyr] = x / a;
          pdn_CKtCN[ilyr] = (a - x) / a;
        } else {
          fprintf (stderr, "pdn, Case 2a: This should NOT happen: a = 0.0 !!!, ilyr = %d \n", ilyr);
          return -104;
        }  //e-if
        pdn_CKtF[ilyr] = 0.0;

        pdn_CNtCK[ilyr] = 0.0;
        if (b > 0.0) {
          pdn_CNtCN[ilyr] = (x + y - a) / b;
          pdn_CNtF[ilyr]  = (a + b - x - y) / b;
        } else {
          fprintf (stderr, "pdn, Case 2a: This should NOT happen: b = 0.0 !!!, ilyr = %d \n", ilyr);
          return -105;
        }  //e-if

        pdn_FtCK[ilyr] = 0.0;
        pdn_FtCN[ilyr] = 0.0;
        pdn_FtF[ilyr]  = 1.0;

      } else {  // Case 2b: Thick cloud above larger than or equal to the entire cloud in this layer;
        if (a > 0.0) {
          pdn_CKtCK[ilyr] = x / a;
          pdn_CKtCN[ilyr] = y / a;
          pdn_CKtF[ilyr]  = (a - x - y) / a;
        } else {
          fprintf (stderr, "pdn, Case 2b: This should NOT happen: a = 0.0 !!!, ilyr = %d \n", ilyr);
          return -106;
        }  //e-if

        pdn_CNtCK[ilyr] = 0.0;
        pdn_CNtCN[ilyr] = 0.0;
        pdn_CNtF[ilyr]  = 1.0;

        pdn_FtCK[ilyr] = 0.0;
        pdn_FtCN[ilyr] = 0.0;
        pdn_FtF[ilyr]  = 1.0;

      }  //e-if-Case2a2b
    }    //e-if-Case12

  }  //e-for over ilyr;

  return 0;
}  //e-calc_pdn_transport_coefficients

//=========================================
//FUNCTION: calc_pup_transport_coefficients
//=========================================
//INPUT:  nlev, cf, cf_ck, cf_cn;
//OUTPUT: pup_CKtCK, pup_CKtCN, pup_CKtF, pup_CNtCK, pup_CNtCN, pup_CNtF, pup_FtCK, pup_FtCN, pup_FtF;

static int calc_pup_transport_coefficients_twomaxrnd3C (int     nlev,
                                                        double* cf,
                                                        double* cf_ck,
                                                        double* cf_cn,
                                                        double* pup_CKtCK,
                                                        double* pup_CKtCN,
                                                        double* pup_CKtF,
                                                        double* pup_CNtCK,
                                                        double* pup_CNtCN,
                                                        double* pup_CNtF,
                                                        double* pup_FtCK,
                                                        double* pup_FtCN,
                                                        double* pup_FtF) {
  const int nlyr = nlev - 1;

  {
    int ilyr = (nlyr - 1);  // SPECIAL CASE

    pup_CKtCK[ilyr] = 1.0;
    pup_CKtCN[ilyr] = 0.0;
    pup_CKtF[ilyr]  = 0.0;

    pup_CNtCK[ilyr] = 0.0;
    pup_CNtCN[ilyr] = 1.0;
    pup_CNtF[ilyr]  = 0.0;

    pup_FtCK[ilyr] = 0.0;
    pup_FtCN[ilyr] = 0.0;
    pup_FtF[ilyr]  = 1.0;
  }

  for (int ilyr = 0; ilyr < (nlyr - 1); ilyr++) {

    const double x = cf_ck[ilyr];         //this-layer cloud thick;
    const double y = cf_cn[ilyr];         //this-layer cloud thin;
    const double z = 1.0 - cf[ilyr];      //this-layer clear-sky;
    const double a = cf_ck[ilyr + 1];     //bottom-layer cloud thick;
    const double b = cf_cn[ilyr + 1];     //bottom-layer cloud thin;
    const double c = 1.0 - cf[ilyr + 1];  //bottom-layer clear-sky;

    if (cf[ilyr + 1] < 0.0000000001) {  // SPECIAL CASE: In the layer below entirely clear-sky (i.e., a=0, b=0, c=1);
      pup_CKtCK[ilyr] = 1.0;
      pup_CKtCN[ilyr] = 0.0;
      pup_CKtF[ilyr]  = 0.0;

      pup_CNtCK[ilyr] = 0.0;
      pup_CNtCN[ilyr] = 1.0;
      pup_CNtF[ilyr]  = 0.0;

      pup_FtCK[ilyr] = x;
      pup_FtCN[ilyr] = y;
      pup_FtF[ilyr]  = z;

      continue;
    }  //e-if

    if (cf[ilyr] < 0.0000000001) {  // SPECIAL CASE: In this layer entirely clear-sky (i.e., x=0, y=0, z=1);
      pup_CKtCK[ilyr] = 0.0;
      pup_CKtCN[ilyr] = 0.0;
      pup_CKtF[ilyr]  = 1.0;

      pup_CNtCK[ilyr] = 0.0;
      pup_CNtCN[ilyr] = 0.0;
      pup_CNtF[ilyr]  = 1.0;

      pup_FtCK[ilyr] = 0.0;  //x
      pup_FtCN[ilyr] = 0.0;  //y
      pup_FtF[ilyr]  = 1.0;  //z

      continue;
    }  //e-if

    if (cf[ilyr] > cf[ilyr + 1]) {  // Case 1: Cloud in this layer larger than cloud below;
                                    //       = thick cloud in this layer larger than thick cloud below;
                                    // (implies also: c not equal to 0.0, but: c > 0.0);

      if (x < (a + b)) {  // Case 1a: Thick cloud in this layer smaller than the entire cloud below;
        pup_CKtCK[ilyr] = 1.0;
        pup_CKtCN[ilyr] = 0.0;
        pup_CKtF[ilyr]  = 0.0;

        if (b > 0.0) {
          pup_CNtCK[ilyr] = (x - a) / b;
          pup_CNtCN[ilyr] = (a + b - x) / b;
        } else {
          fprintf (stderr, "pup, Case 1a: This should NOT happen: b = 0.0 !!!, ilyr = %d \n", ilyr);
          return -201;
        }  //e-if
        pup_CNtF[ilyr] = 0.0;

        pup_FtCK[ilyr] = 0.0;
        if (c > 0.0) {
          pup_FtCN[ilyr] = (x + y - a - b) / c;
          pup_FtF[ilyr]  = z / c;
        } else {
          fprintf (stderr, "pup, Case 1a: This should NOT happen: c = 0.0 !!!, ilyr = %d \n", ilyr);
          return -202;
        }  //e-if

      } else {  // Case 1b: Thick cloud in this layer larger than or equal to the entire cloud below;
        pup_CKtCK[ilyr] = 1.0;
        pup_CKtCN[ilyr] = 0.0;
        pup_CKtF[ilyr]  = 0.0;

        pup_CNtCK[ilyr] = 1.0;
        pup_CNtCN[ilyr] = 0.0;
        pup_CNtF[ilyr]  = 0.0;

        if (c > 0.0) {
          pup_FtCK[ilyr] = (x - a - b) / c;
          pup_FtCN[ilyr] = y / c;
          pup_FtF[ilyr]  = z / c;
        } else {
          fprintf (stderr, "pup, Case 1b: This should NOT happen: c = 0.0 !!!, ilyr = %d \n", ilyr);
          return -203;
        }  //e-if

      }  //e-if-Case1a1b

    } else {  // Case 2: Cloud in this layer smaller than or equal to the cloud below;
              //       = thick cloud in this layer smaller than or equal to the thick cloud below;

      if (a < (x + y)) {  // Case 2a: Thick cloud below smaller than the entire cloud in this layer;
        if (a > 0.0) {
          pup_CKtCK[ilyr] = x / a;
          pup_CKtCN[ilyr] = (a - x) / a;
        } else {
          fprintf (stderr, "pup, Case 2a: This should NOT happen: a = 0.0 !!!, ilyr = %d \n", ilyr);
          return -204;
        }  //e-if
        pup_CKtF[ilyr] = 0.0;

        pup_CNtCK[ilyr] = 0.0;
        if (b > 0.0) {
          pup_CNtCN[ilyr] = (x + y - a) / b;
          pup_CNtF[ilyr]  = (a + b - x - y) / b;
        } else {
          fprintf (stderr, "pup, Case 2a: This should NOT happen: b = 0.0 !!!, ilyr = %d \n", ilyr);
          return -205;
        }  //e-if

        pup_FtCK[ilyr] = 0.0;
        pup_FtCN[ilyr] = 0.0;
        pup_FtF[ilyr]  = 1.0;

      } else {  // Case 2b: Thick cloud below larger than or equal to the entire cloud in this layer;
        if (a > 0.0) {
          pup_CKtCK[ilyr] = x / a;
          pup_CKtCN[ilyr] = y / a;
          pup_CKtF[ilyr]  = (a - x - y) / a;
        } else {
          fprintf (stderr, "pup, Case 2b: This should NOT happen: a = 0.0 !!!, ilyr = %d \n", ilyr);
          return -206;
        }  //e-if

        pup_CNtCK[ilyr] = 0.0;
        pup_CNtCN[ilyr] = 0.0;
        pup_CNtF[ilyr]  = 1.0;

        pup_FtCK[ilyr] = 0.0;
        pup_FtCN[ilyr] = 0.0;
        pup_FtF[ilyr]  = 1.0;

      }  //e-if-Case2a2b
    }    //e-if-Case12

  }  //e-for over ilyr;

  return 0;
}  //e-calc_pup_transport_coefficient

//=======================
// FUNCTION: buildMatrixA
//=======================
// INPUT:  nlev, Ag, ar_a11_ck, ar_a11_cn, ar_a11_f, ar_a12_ck, ar_a12_cn, ar_a12_f,
//         pdn_CKtCK, pdn_CKtCN, pdn_CKtF, pdn_CNtCK, pdn_CNtCN, pdn_CNtF, pdn_FtCK, pdn_FtCN, pdn_FtF,
//         pup_CKtCK, pup_CKtCN, pup_CKtF, pup_CNtCK, pup_CNtCN, pup_CNtF, pup_FtCK, pup_FtCN, pup_FtF;
// OUTPUT: matrixA;

static int buildMatrixA_twomaxrnd3C (int      nlev,
                                     double   Ag,
                                     double*  ar_a11_ck,
                                     double*  ar_a11_cn,
                                     double*  ar_a11_f,
                                     double*  ar_a12_ck,
                                     double*  ar_a12_cn,
                                     double*  ar_a12_f,
                                     double*  pdn_CKtCK,
                                     double*  pdn_CKtCN,
                                     double*  pdn_CKtF,
                                     double*  pdn_CNtCK,
                                     double*  pdn_CNtCN,
                                     double*  pdn_CNtF,
                                     double*  pdn_FtCK,
                                     double*  pdn_FtCN,
                                     double*  pdn_FtF,
                                     double*  pup_CKtCK,
                                     double*  pup_CKtCN,
                                     double*  pup_CKtF,
                                     double*  pup_CNtCK,
                                     double*  pup_CNtCN,
                                     double*  pup_CNtF,
                                     double*  pup_FtCK,
                                     double*  pup_FtCN,
                                     double*  pup_FtF,
                                     double** matrixA) {
  int i;     // position in levels (e.g. for nlev=21, i in range: 0 - 20)
  int iRow;  // row position in matrix A (e.g. for nlev=21, iRow in range: 0 - 125)

  int nextColEup;  // index of column for next Eup data
  int nextColEdw;  // index of column for next Edn data; nextColEdn = nextColEup - 6

  // Set values for initial six rows of matrix A:
  // First row:
  matrixA[0][3] = ar_a12_f[0] * pdn_FtF[0];
  matrixA[0][4] = ar_a12_f[0] * pdn_CKtF[0];
  matrixA[0][5] = ar_a12_f[0] * pdn_CNtF[0];

  matrixA[0][6] = ar_a11_f[0] * pup_FtF[0];
  matrixA[0][7] = ar_a11_f[0] * pup_CKtF[0];
  matrixA[0][8] = ar_a11_f[0] * pup_CNtF[0];

  // Second row:
  matrixA[1][3] = ar_a12_ck[0] * pdn_FtCK[0];
  matrixA[1][4] = ar_a12_ck[0] * pdn_CKtCK[0];
  matrixA[1][5] = ar_a12_ck[0] * pdn_CNtCK[0];

  matrixA[1][6] = ar_a11_ck[0] * pup_FtCK[0];
  matrixA[1][7] = ar_a11_ck[0] * pup_CKtCK[0];
  matrixA[1][8] = ar_a11_ck[0] * pup_CNtCK[0];

  // Third row:
  matrixA[2][3] = ar_a12_cn[0] * pdn_FtCN[0];
  matrixA[2][4] = ar_a12_cn[0] * pdn_CKtCN[0];
  matrixA[2][5] = ar_a12_cn[0] * pdn_CNtCN[0];

  matrixA[2][6] = ar_a11_cn[0] * pup_FtCN[0];
  matrixA[2][7] = ar_a11_cn[0] * pup_CKtCN[0];
  matrixA[2][8] = ar_a11_cn[0] * pup_CNtCN[0];

  // 4th row is already zero; (needs to be zero due to upper boundary condition);
  // 5th row is already zero; (needs to be zero due to upper boundary condition);
  // 6th row is already zero; (needs to be zero due to upper boundary condition);

  nextColEup = 9;  // index of column for Eup(level1)
  nextColEdw = 3;  // index of column for Edn(level1)

  for (i = 1; i < nlev; i++) {
    iRow = NFLUX * i;

    if (i == (nlev - 1)) {  // lower boundary condition for the forth and third row from bottom up of matrix A;
      matrixA[iRow][nextColEup]         = Ag;
      matrixA[iRow + 1][nextColEup + 1] = Ag;
      matrixA[iRow + 2][nextColEup + 2] = Ag;
    } else {
      matrixA[iRow][nextColEup]     = ar_a12_f[i] * pdn_FtF[i];
      matrixA[iRow][nextColEup + 1] = ar_a12_f[i] * pdn_CKtF[i];
      matrixA[iRow][nextColEup + 2] = ar_a12_f[i] * pdn_CNtF[i];

      matrixA[iRow][nextColEup + 3] = ar_a11_f[i] * pup_FtF[i];
      matrixA[iRow][nextColEup + 4] = ar_a11_f[i] * pup_CKtF[i];
      matrixA[iRow][nextColEup + 5] = ar_a11_f[i] * pup_CNtF[i];

      matrixA[iRow + 1][nextColEup]     = ar_a12_ck[i] * pdn_FtCK[i];
      matrixA[iRow + 1][nextColEup + 1] = ar_a12_ck[i] * pdn_CKtCK[i];
      matrixA[iRow + 1][nextColEup + 2] = ar_a12_ck[i] * pdn_CNtCK[i];

      matrixA[iRow + 1][nextColEup + 3] = ar_a11_ck[i] * pup_FtCK[i];
      matrixA[iRow + 1][nextColEup + 4] = ar_a11_ck[i] * pup_CKtCK[i];
      matrixA[iRow + 1][nextColEup + 5] = ar_a11_ck[i] * pup_CNtCK[i];

      matrixA[iRow + 2][nextColEup]     = ar_a12_cn[i] * pdn_FtCN[i];
      matrixA[iRow + 2][nextColEup + 1] = ar_a12_cn[i] * pdn_CKtCN[i];
      matrixA[iRow + 2][nextColEup + 2] = ar_a12_cn[i] * pdn_CNtCN[i];

      matrixA[iRow + 2][nextColEup + 3] = ar_a11_cn[i] * pup_FtCN[i];
      matrixA[iRow + 2][nextColEup + 4] = ar_a11_cn[i] * pup_CKtCN[i];
      matrixA[iRow + 2][nextColEup + 5] = ar_a11_cn[i] * pup_CNtCN[i];
    }  //e-if

    nextColEup += NFLUX;

    matrixA[iRow + 3][nextColEdw]     = ar_a11_f[i - 1] * pdn_FtF[i - 1];
    matrixA[iRow + 3][nextColEdw + 1] = ar_a11_f[i - 1] * pdn_CKtF[i - 1];
    matrixA[iRow + 3][nextColEdw + 2] = ar_a11_f[i - 1] * pdn_CNtF[i - 1];

    matrixA[iRow + 3][nextColEdw + 3] = ar_a12_f[i - 1] * pup_FtF[i - 1];
    matrixA[iRow + 3][nextColEdw + 4] = ar_a12_f[i - 1] * pup_CKtF[i - 1];
    matrixA[iRow + 3][nextColEdw + 5] = ar_a12_f[i - 1] * pup_CNtF[i - 1];

    matrixA[iRow + 4][nextColEdw]     = ar_a11_ck[i - 1] * pdn_FtCK[i - 1];
    matrixA[iRow + 4][nextColEdw + 1] = ar_a11_ck[i - 1] * pdn_CKtCK[i - 1];
    matrixA[iRow + 4][nextColEdw + 2] = ar_a11_ck[i - 1] * pdn_CNtCK[i - 1];

    matrixA[iRow + 4][nextColEdw + 3] = ar_a12_ck[i - 1] * pup_FtCK[i - 1];
    matrixA[iRow + 4][nextColEdw + 4] = ar_a12_ck[i - 1] * pup_CKtCK[i - 1];
    matrixA[iRow + 4][nextColEdw + 5] = ar_a12_ck[i - 1] * pup_CNtCK[i - 1];

    matrixA[iRow + 5][nextColEdw]     = ar_a11_cn[i - 1] * pdn_FtCN[i - 1];
    matrixA[iRow + 5][nextColEdw + 1] = ar_a11_cn[i - 1] * pdn_CKtCN[i - 1];
    matrixA[iRow + 5][nextColEdw + 2] = ar_a11_cn[i - 1] * pdn_CNtCN[i - 1];

    matrixA[iRow + 5][nextColEdw + 3] = ar_a12_cn[i - 1] * pup_FtCN[i - 1];
    matrixA[iRow + 5][nextColEdw + 4] = ar_a12_cn[i - 1] * pup_CKtCN[i - 1];
    matrixA[iRow + 5][nextColEdw + 5] = ar_a12_cn[i - 1] * pup_CNtCN[i - 1];

    nextColEdw += NFLUX;
  }  //e-for

  return 0;
}  //e-buildMatrixA

//=========================
// FUNCTION buildVectorBsol
//=========================
// INPUT:  nlev, Ag, mu0, ar_a13_ck, ar_a13_cn, ar_a13_f, ar_a23_ck, ar_a23_cn, ar_a23_f, ar_S_ck, ar_S_cn, ar_S_f,
//         pdn_CKtCK, pdn_CKtCN, pdn_CKtF, pdn_CNtCK, pdn_CNtCN, pdn_CNtF, pdn_FtCK,  pdn_FtCN,  pdn_FtF;
// OUTPUT: vectB;

static int buildVectorBsol_twomaxrnd3C (int     nlev,
                                        double  Ag,
                                        double  mu0,
                                        double* ar_a13_ck,
                                        double* ar_a13_cn,
                                        double* ar_a13_f,
                                        double* ar_a23_ck,
                                        double* ar_a23_cn,
                                        double* ar_a23_f,
                                        double* ar_S_ck,
                                        double* ar_S_cn,
                                        double* ar_S_f,
                                        double* pdn_CKtCK,
                                        double* pdn_CKtCN,
                                        double* pdn_CKtF,
                                        double* pdn_CNtCK,
                                        double* pdn_CNtCN,
                                        double* pdn_CNtF,
                                        double* pdn_FtCK,
                                        double* pdn_FtCN,
                                        double* pdn_FtF,
                                        double* vectB) {
  // Set initial six values:
  vectB[0] = ar_a13_f[0] * ar_S_f[0];
  vectB[1] = ar_a13_ck[0] * ar_S_ck[0];
  vectB[2] = ar_a13_cn[0] * ar_S_cn[0];
  vectB[3] = 0.0;  // upper boundary condition
  vectB[4] = 0.0;  // upper boundary condition
  vectB[5] = 0.0;  // upper boundary condition

  for (int i = 1; i < (nlev - 1); i++) {
    const int j  = NFLUX * i;
    vectB[j]     = ar_a13_f[i] * (pdn_FtF[i] * ar_S_f[i] + pdn_CKtF[i] * ar_S_ck[i] + pdn_CNtF[i] * ar_S_cn[i]);
    vectB[j + 1] = ar_a13_ck[i] * (pdn_FtCK[i] * ar_S_f[i] + pdn_CKtCK[i] * ar_S_ck[i] + pdn_CNtCK[i] * ar_S_cn[i]);
    vectB[j + 2] = ar_a13_cn[i] * (pdn_FtCN[i] * ar_S_f[i] + pdn_CKtCN[i] * ar_S_ck[i] + pdn_CNtCN[i] * ar_S_cn[i]);
    vectB[j + 3] =
      ar_a23_f[i - 1] * (pdn_FtF[i - 1] * ar_S_f[i - 1] + pdn_CKtF[i - 1] * ar_S_ck[i - 1] + pdn_CNtF[i - 1] * ar_S_cn[i - 1]);
    vectB[j + 4] =
      ar_a23_ck[i - 1] * (pdn_FtCK[i - 1] * ar_S_f[i - 1] + pdn_CKtCK[i - 1] * ar_S_ck[i - 1] + pdn_CNtCK[i - 1] * ar_S_cn[i - 1]);
    vectB[j + 5] =
      ar_a23_cn[i - 1] * (pdn_FtCN[i - 1] * ar_S_f[i - 1] + pdn_CKtCN[i - 1] * ar_S_ck[i - 1] + pdn_CNtCN[i - 1] * ar_S_cn[i - 1]);
  }  //e-for

  // Treat last six values separately:
  const int j  = NFLUX * (nlev - 1);
  vectB[j]     = Ag * mu0 * ar_S_f[nlev - 1];   // lower boundary condition
  vectB[j + 1] = Ag * mu0 * ar_S_ck[nlev - 1];  // lower boundary condition
  vectB[j + 2] = Ag * mu0 * ar_S_cn[nlev - 1];  // lower boundary condition
  vectB[j + 3] = ar_a23_f[nlev - 2] * (pdn_FtF[nlev - 2] * ar_S_f[nlev - 2] + pdn_CKtF[nlev - 2] * ar_S_ck[nlev - 2] +
                                       pdn_CNtF[nlev - 2] * ar_S_cn[nlev - 2]);
  vectB[j + 4] = ar_a23_ck[nlev - 2] * (pdn_FtCK[nlev - 2] * ar_S_f[nlev - 2] + pdn_CKtCK[nlev - 2] * ar_S_ck[nlev - 2] +
                                        pdn_CNtCK[nlev - 2] * ar_S_cn[nlev - 2]);
  vectB[j + 5] = ar_a23_cn[nlev - 2] * (pdn_FtCN[nlev - 2] * ar_S_f[nlev - 2] + pdn_CKtCN[nlev - 2] * ar_S_ck[nlev - 2] +
                                        pdn_CNtCN[nlev - 2] * ar_S_cn[nlev - 2]);

  return 0;
}  //e-buildVectorBsol

//=========================
// FUNCTION buildVectorBthe
//=========================
// INPUT:  nlev, Ag, Bg, ar_theComp1_ck, ar_theComp1_cn, ar_theComp1_f,
//                       ar_theComp2_ck, ar_theComp2_cn, ar_theComp2_f, cf, cf_ck, cf_cn;
// OUTPUT: bb_the;

// thermal component 1 = upward emission
// thermal component 2 = downward emission

static int buildVectorBthe_twomaxrnd3C (int     nlev,
                                        double  Ag,
                                        double  Bg,
                                        double* ar_theComp1_ck,
                                        double* ar_theComp1_cn,
                                        double* ar_theComp1_f,
                                        double* ar_theComp2_ck,
                                        double* ar_theComp2_cn,
                                        double* ar_theComp2_f,
                                        double* cf,
                                        double* cf_ck,
                                        double* cf_cn,
                                        double* vectB) {
  // Set initial six values:
  vectB[0] = (1.0 - cf[0]) * ar_theComp1_f[0];
  vectB[1] = cf_ck[0] * ar_theComp1_ck[0];
  vectB[2] = cf_cn[0] * ar_theComp1_cn[0];
  vectB[3] = 0.0;  // upper boundary condition
  vectB[4] = 0.0;  // upper boundary condition
  vectB[5] = 0.0;  // upper boundary condition

  for (int i = 1; i < (nlev - 1); i++) {
    const int j  = NFLUX * i;
    vectB[j]     = (1.0 - cf[i]) * ar_theComp1_f[i];
    vectB[j + 1] = cf_ck[i] * ar_theComp1_ck[i];
    vectB[j + 2] = cf_cn[i] * ar_theComp1_cn[i];
    vectB[j + 3] = (1.0 - cf[i - 1]) * ar_theComp2_f[i - 1];
    vectB[j + 4] = cf_ck[i - 1] * ar_theComp2_ck[i - 1];
    vectB[j + 5] = cf_cn[i - 1] * ar_theComp2_cn[i - 1];
  }  //e-for

  // Treat last six values seperately:
  // i=nlev-1; // bottom (ground) level;
  const int j  = NFLUX * (nlev - 1);
  vectB[j]     = (1.0 - cf[nlev - 2]) * (1.0 - Ag) * M_PI * Bg;  // lower boundary condition
  vectB[j + 1] = cf_ck[nlev - 2] * (1.0 - Ag) * M_PI * Bg;       // lower boundary condition
  vectB[j + 2] = cf_cn[nlev - 2] * (1.0 - Ag) * M_PI * Bg;       // lower boundary condition
  vectB[j + 3] = (1.0 - cf[nlev - 2]) * ar_theComp2_f[nlev - 2];
  vectB[j + 4] = cf_ck[nlev - 2] * ar_theComp2_ck[nlev - 2];
  vectB[j + 5] = cf_cn[nlev - 2] * ar_theComp2_cn[nlev - 2];

  return 0;
}  //e-buildVectorBthe

//=============================
// FUNCTION twostream_maxrand3C
//=============================

/*
OUTPUT:
Edir = direct irradiance;
Edn  = diffuse downward irradiance;
Eup  = diffuse upward irradiance;
Lavg = radiance (set to NaN);
*/

static int twostream_maxrand3C (double* dtau_org_ck,
                                double* omega0_org_ck,
                                double* g_org_ck,  // "_ck" for optically thick cloudy region
                                double* dtau_org_cn,
                                double* omega0_org_cn,
                                double* g_org_cn,  // "_cn" for optically thin cloudy region
                                double* dtau_org_f,
                                double* omega0_org_f,
                                double* g_org_f,  // "_f"  for cloud-free region
                                double* cf,
                                double  scale_cf,
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

  double a11_ck[nlyr];
  double a12_ck[nlyr];
  double a13_ck[nlyr];
  double a23_ck[nlyr];
  double a33_ck[nlyr];

  double a11_cn[nlyr];
  double a12_cn[nlyr];
  double a13_cn[nlyr];
  double a23_cn[nlyr];
  double a33_cn[nlyr];

  double a11_f[nlyr];
  double a12_f[nlyr];
  double a13_f[nlyr];
  double a23_f[nlyr];
  double a33_f[nlyr];

  // Components of vector B in the thermal spectral range:
  double theComp1_ck[nlev];
  double theComp1_cn[nlev];
  double theComp1_f[nlev];

  double theComp2_ck[nlev];
  double theComp2_cn[nlev];
  double theComp2_f[nlev];

  // Transport coefficients:
  double pdn_CKtCK[nlyr];
  double pdn_CKtCN[nlyr];
  double pdn_CKtF[nlyr];

  double pdn_CNtCK[nlyr];
  double pdn_CNtCN[nlyr];
  double pdn_CNtF[nlyr];

  double pdn_FtCK[nlyr];
  double pdn_FtCN[nlyr];
  double pdn_FtF[nlyr];

  double pup_CKtCK[nlyr];
  double pup_CKtCN[nlyr];
  double pup_CKtF[nlyr];

  double pup_CNtCK[nlyr];
  double pup_CNtCN[nlyr];
  double pup_CNtF[nlyr];

  double pup_FtCK[nlyr];
  double pup_FtCN[nlyr];
  double pup_FtF[nlyr];

  double S_ck[nlev];
  double S_cn[nlev];
  double S_f[nlev];

  double Edir_ck[nlev];
  double Edir_cn[nlev];
  double Edir_f[nlev];

  double Eup_ck[nlev];
  double Eup_cn[nlev];
  double Eup_f[nlev];

  double Edn_ck[nlev];
  double Edn_cn[nlev];
  double Edn_f[nlev];

  double cf_ck[nlyr];
  double cf_cn[nlyr];

  const long Nrow = NFLUX * nlev;

  double** AA = calloc_double_2D (Nrow, Nrow, "AA");
  CHKPOINTER (AA);
  double bb_sol[Nrow];
  double bb_the[Nrow];
  double bb[Nrow];
  double xx[Nrow];  // result vector

  // Check if scale_cf is set properly (0.0 <= scale_cf >= 1.0):
  if (scale_cf < 0.0 || scale_cf > 1.0) {
    CHKERROUT (-2, "Error - invalid value of input parameter: scale_cf should be between 0.0 and 1.0 !!!");
  }  //e-if

  // Calculate cf_ck and cf_cn:
  for (int ilyr = 0; ilyr < nlyr; ilyr++) {
    cf_ck[ilyr] = scale_cf * cf[ilyr];
    cf_cn[ilyr] = (1.0 - scale_cf) * cf[ilyr];
  }  //e-for

  // At the moment it is only possible to calculate solar OR thermal RT, but not both simultaneously;
  if ((flagSolar && flagThermal) || (!flagSolar && !flagThermal)) {
    CHKERROUT (-3, "Error - invalid input parameters - use flagSolar and flagThermal alternatingly");
  }

  // Calculate vertical profiles of transport coefficients pdn_... from vertical profiles of partial cloudiness:
  status = calc_pdn_transport_coefficients_twomaxrnd3C (nlev,
                                                        cf,
                                                        cf_ck,
                                                        cf_cn,
                                                        pdn_CKtCK,
                                                        pdn_CKtCN,
                                                        pdn_CKtF,
                                                        pdn_CNtCK,
                                                        pdn_CNtCN,
                                                        pdn_CNtF,
                                                        pdn_FtCK,
                                                        pdn_FtCN,
                                                        pdn_FtF);
  CHKERR (status);

  // Calculate vertical profiles of transport coefficients pup_... from vertical profiles of partial cloudiness:
  status = calc_pup_transport_coefficients_twomaxrnd3C (nlev,
                                                        cf,
                                                        cf_ck,
                                                        cf_cn,
                                                        pup_CKtCK,
                                                        pup_CKtCN,
                                                        pup_CKtF,
                                                        pup_CNtCK,
                                                        pup_CNtCN,
                                                        pup_CNtF,
                                                        pup_FtCK,
                                                        pup_FtCN,
                                                        pup_FtF);
  CHKERR (status);

  // Calculate vertical profiles of Eddington coefficients
  // and vertical profiles of thermal components:
  for (int ilyr = 0; ilyr < nlyr; ilyr++) {

    double dtau_ck;
    double omega0_ck;
    double g_ck;
    double dtau_cn;
    double omega0_cn;
    double g_cn;

    // Delta scaling of optical properties of cloudy regions:
    if (delta) {
      delta_scale_hg (dtau_org_ck[ilyr], omega0_org_ck[ilyr], g_org_ck[ilyr], &dtau_ck, &omega0_ck, &g_ck);
      delta_scale_hg (dtau_org_cn[ilyr], omega0_org_cn[ilyr], g_org_cn[ilyr], &dtau_cn, &omega0_cn, &g_cn);
    } else {
      dtau_ck   = dtau_org_ck[ilyr];
      omega0_ck = omega0_org_ck[ilyr];
      g_ck      = g_org_ck[ilyr];

      dtau_cn   = dtau_org_cn[ilyr];
      omega0_cn = omega0_org_cn[ilyr];
      g_cn      = g_org_cn[ilyr];
    }  //e-if
    // omega0 should not be 1 (avoiding singularity problem), restrict omega0 to 0.999999:
    omega0_ck = MIN (omega0_ck, 0.999999);
    omega0_cn = MIN (omega0_cn, 0.999999);

    // Optical properties of cloud-free regions:
    const double dtau_f   = dtau_org_f[ilyr];
    const double omega0_f = MIN (omega0_org_f[ilyr], 0.999999);
    const double g_f      = g_org_f[ilyr];

    eddington_coeffc (dtau_ck,
                      g_ck,
                      omega0_ck,
                      mu0,
                      &(a11_ck[ilyr]),
                      &(a12_ck[ilyr]),
                      &(a13_ck[ilyr]),
                      &(a23_ck[ilyr]),
                      &(a33_ck[ilyr]));
    eddington_coeffc (dtau_cn,
                      g_cn,
                      omega0_cn,
                      mu0,
                      &(a11_cn[ilyr]),
                      &(a12_cn[ilyr]),
                      &(a13_cn[ilyr]),
                      &(a23_cn[ilyr]),
                      &(a33_cn[ilyr]));
    eddington_coeffc (dtau_f, g_f, omega0_f, mu0, &(a11_f[ilyr]), &(a12_f[ilyr]), &(a13_f[ilyr]), &(a23_f[ilyr]), &(a33_f[ilyr]));

    if (flagThermal) {
      calcThermalComponents (ilyr, B, dtau_ck, omega0_ck, g_ck, &(theComp1_ck[ilyr]), &(theComp2_ck[ilyr]));
      calcThermalComponents (ilyr, B, dtau_cn, omega0_cn, g_cn, &(theComp1_cn[ilyr]), &(theComp2_cn[ilyr]));
      calcThermalComponents (ilyr, B, dtau_f, omega0_f, g_f, &(theComp1_f[ilyr]), &(theComp2_f[ilyr]));
    }  //e-if
  }    //e-for

  // Initialize vectors S_ck, S_cn and S_f:
  // S[0]= S0 = S_ck[0] + S_cn[0] + S_f[0] = cf[0]*S0 + (1.0-cf[0])*S0 = cf_ck[0]*S0 + cf_cn[0]*S0 + (1.0-cf[0])*S0
  if (flagSolar) {
    S_ck[0] = cf_ck[0] * S0;
    S_cn[0] = cf_cn[0] * S0;
    S_f[0]  = (1.0 - cf[0]) * S0;

    for (int ilev = 1; ilev < nlev; ilev++) {
      S_ck[ilev] = a33_ck[ilev - 1] * (pdn_FtCK[ilev - 1] * S_f[ilev - 1] + pdn_CKtCK[ilev - 1] * S_ck[ilev - 1] +
                                       pdn_CNtCK[ilev - 1] * S_cn[ilev - 1]);
      S_cn[ilev] = a33_cn[ilev - 1] * (pdn_FtCN[ilev - 1] * S_f[ilev - 1] + pdn_CKtCN[ilev - 1] * S_ck[ilev - 1] +
                                       pdn_CNtCN[ilev - 1] * S_cn[ilev - 1]);
      S_f[ilev]  = a33_f[ilev - 1] *
                  (pdn_FtF[ilev - 1] * S_f[ilev - 1] + pdn_CKtF[ilev - 1] * S_ck[ilev - 1] + pdn_CNtF[ilev - 1] * S_cn[ilev - 1]);
    }  //e-for
  }    //e-if

  // Equation system has the following form: xx = AA*xx + bb;

  // Build matrix AA:
  status = buildMatrixA_twomaxrnd3C (nlev,
                                     Ag,
                                     a11_ck,
                                     a11_cn,
                                     a11_f,
                                     a12_ck,
                                     a12_cn,
                                     a12_f,
                                     pdn_CKtCK,
                                     pdn_CKtCN,
                                     pdn_CKtF,
                                     pdn_CNtCK,
                                     pdn_CNtCN,
                                     pdn_CNtF,
                                     pdn_FtCK,
                                     pdn_FtCN,
                                     pdn_FtF,
                                     pup_CKtCK,
                                     pup_CKtCN,
                                     pup_CKtF,
                                     pup_CNtCK,
                                     pup_CNtCN,
                                     pup_CNtF,
                                     pup_FtCK,
                                     pup_FtCN,
                                     pup_FtF,
                                     AA);
  CHKERR (status);

  // Make matrix A2 = AA - IdentityMatrix and save it to the same matrix AA,
  // since the equation system in the form Ax=b must be passed to solve_gauss;
  // xx = AA*xx + bb;
  // (AA-II)*xx = -bb;
  for (int iRow = 0; iRow < Nrow; iRow++)
    AA[iRow][iRow] = -1.0;

  // Build vector bb_sol:
  if (flagSolar) {
    status = buildVectorBsol_twomaxrnd3C (nlev,
                                          Ag,
                                          mu0,
                                          a13_ck,
                                          a13_cn,
                                          a13_f,
                                          a23_ck,
                                          a23_cn,
                                          a23_f,
                                          S_ck,
                                          S_cn,
                                          S_f,
                                          pdn_CKtCK,
                                          pdn_CKtCN,
                                          pdn_CKtF,
                                          pdn_CNtCK,
                                          pdn_CNtCN,
                                          pdn_CNtF,
                                          pdn_FtCK,
                                          pdn_FtCN,
                                          pdn_FtF,
                                          bb_sol);
    CHKERR (status);
  }  //e-if

  // Build vector bb_the:
  if (flagThermal) {
    status = buildVectorBthe_twomaxrnd3C (nlev,
                                          Ag,
                                          Bg,
                                          theComp1_ck,
                                          theComp1_cn,
                                          theComp1_f,
                                          theComp2_ck,
                                          theComp2_cn,
                                          theComp2_f,
                                          cf,
                                          cf_ck,
                                          cf_cn,
                                          bb_the);
    CHKERR (status);
  }  //e-if

  // Assign bb_sol or bb_the to the final vector bb:
  for (int ilev = 0; ilev < Nrow; ilev++) {
    bb[ilev] = 0;

    if (flagSolar) {
      bb[ilev] += bb_sol[ilev];
    }  //e-if

    if (flagThermal) {
      bb[ilev] += bb_the[ilev];
    }  //e-if

  }  //e-for

  // Make vector -bb and save it to the same vector bb:
  status = vec_scale (Nrow, -1, bb);
  CHKERR (status);

  const long m1   = 8;
  const long m2   = 8;  // n=2*input.nlyr+2;
  double**   matr = malloc_double_2D (Nrow, m1 + m2 + 1, "matr");

  matrix_to_band_form (AA, Nrow, m1, m2, matr);  // Transform matrix AA to band form matr
  free_double_2D (AA);

  {
    long     indx[Nrow];
    double** mlu = malloc_double_2D (Nrow, m2, "mlu");
    bandec_M (matr, Nrow, m1, m2, mlu, indx);
    banbks_M (matr, Nrow, m1, m2, mlu, indx, bb);
    free_double_2D (mlu);
  }
  free_double_2D (matr);

  for (int i = 0; i < Nrow; i++) {
    xx[i] = bb[i];
  }

  for (int ilev = 0; ilev < nlev; ilev++) {
    Eup_f[ilev]  = xx[NFLUX * ilev];
    Eup_ck[ilev] = xx[NFLUX * ilev + 1];
    Eup_cn[ilev] = xx[NFLUX * ilev + 2];

    Edn_f[ilev]  = xx[NFLUX * ilev + 3];
    Edn_ck[ilev] = xx[NFLUX * ilev + 4];
    Edn_cn[ilev] = xx[NFLUX * ilev + 5];
  }  //e-for

  for (int ilev = 0; ilev < nlev; ilev++) {
    Edir_ck[ilev] = S_ck[ilev] * mu0;
    Edir_cn[ilev] = S_cn[ilev] * mu0;
    Edir_f[ilev]  = S_f[ilev] * mu0;
  }  //e-for

  // Sum up the irradiances of the three regions to obtain the final result:
  for (int ilev = 0; ilev < nlev; ilev++) {
    Edir[ilev] = Edir_ck[ilev] + Edir_cn[ilev] + Edir_f[ilev];
    Eup[ilev]  = Eup_ck[ilev] + Eup_cn[ilev] + Eup_f[ilev];
    Edn[ilev]  = Edn_ck[ilev] + Edn_cn[ilev] + Edn_f[ilev];
    Lavg[ilev] = NAN;
  }  //e-for

  return 0;
}  //e-twostream_maxrand3C

int twomaxrnd3C (float* dtau_cldk,
                 float* omega0_cldk,
                 float* g_cldk,
                 float* dtau_cldn,
                 float* omega0_cldn,
                 float* g_cldn,
                 float* dtau_clr,
                 float* omega0_clr,
                 float* g_clr,
                 float* cf,
                 float  scale_cf,
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

  double dtau_cldk_org_d[nlev - 1];
  double omega0_cldk_org_d[nlev - 1];
  double g_cldk_org_d[nlev - 1];

  double dtau_cldn_org_d[nlev - 1];
  double omega0_cldn_org_d[nlev - 1];
  double g_cldn_org_d[nlev - 1];

  double cf_d[nlev - 1];
  double B[nlev];

  double Edir[nlev];
  double Edn[nlev];
  double Eup[nlev];
  double Lavg[nlev];

  const double taumax = 100.0;
  double       Bg     = 0;  // Used only for thermal RT; unused (0.0) for solar RT;

  const int flagSolar = S0 > 0.0;

  if (planck) {
    float plkavg;

    for (int ilev = 0; ilev < nlev; ilev++) {  // level temperatures and Planck functions
      F77_FUNC (cplkavg, CPLKAVG) (&wvnmlo, &wvnmhi, &temper[ilev], &plkavg);
      B[ilev] = plkavg;
    }

    // surface temperature and Planck function
    F77_FUNC (cplkavg, CPLKAVG) (&wvnmlo, &wvnmhi, &btemp, &plkavg);
    Bg = plkavg;
  }

  // copy float arrays to double arrays
  for (int ilyr = 0; ilyr < nlev - 1; ilyr++) {
    // restrict layer optical thickness to 100
    dtau_cldk_org_d[ilyr] = MIN (dtau_cldk[ilyr], taumax);
    dtau_cldn_org_d[ilyr] = MIN (dtau_cldn[ilyr], taumax);
    dtau_clr_org_d[ilyr]  = MIN (dtau_clr[ilyr], taumax);

    omega0_cldk_org_d[ilyr] = omega0_cldk[ilyr];
    g_cldk_org_d[ilyr]      = g_cldk[ilyr];

    omega0_cldn_org_d[ilyr] = omega0_cldn[ilyr];
    g_cldn_org_d[ilyr]      = g_cldn[ilyr];

    omega0_clr_org_d[ilyr] = omega0_clr[ilyr];
    g_clr_org_d[ilyr]      = g_clr[ilyr];

    cf_d[ilyr] = cf[ilyr];
  }  //e-for

  // Call twostream_maxrand3C code:
  int status = twostream_maxrand3C (dtau_cldk_org_d,
                                    omega0_cldk_org_d,
                                    g_cldk_org_d,
                                    dtau_cldn_org_d,
                                    omega0_cldn_org_d,
                                    g_cldn_org_d,
                                    dtau_clr_org_d,
                                    omega0_clr_org_d,
                                    g_clr_org_d,
                                    cf_d,
                                    scale_cf,
                                    nlev,
                                    S0,
                                    mu0,
                                    Ag,
                                    Bg,
                                    B,
                                    delta,
                                    flagSolar,
                                    planck,
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
}  //e-twomaxrnd3C
