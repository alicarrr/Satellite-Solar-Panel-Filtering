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

/****************************************************************************/
/* This file is an interface for the Fortran ocean BRDF routine oceabrdf.   */
/* It turns out that DISORT2 requires many calls to this routine, something */
/* like 140,000 calls for an ordinary 16 streams disort run. As the         */
/* oceabrdf routine is quite complicated, the calculation of the BRDF       */
/* causes the major part of the computational cost of the disort call.      */
/*                                                                          */
/* In many applications, disort is called several times for the same        */
/* wavelength during one uvspec run (all correlated_k applications).        */
/* In this case, the same 140,000 oceabrdf calls are repeated for each      */
/* individual wavelength.                                                   */
/*                                                                          */
/* ocean stores all quadruplets (mu, mup, dphi, bdref) that are calculated  */
/* in a search tree (code courtesy of Kearnighan and Ritchie).              */
/* A search tree is an efficient way of storing and extracting large        */
/* amounts of data the amount of which is not known a priori. Data are      */
/* added with addtree() and may be extracted again with searchtree().       */
/* Hence, for each now input triple (mu, mup, dphi) it is checked if the    */
/* value is already available in the tree. If yes, it is got from there;    */
/* if not, it is calculated and stored in the tree.                         */
/*                                                                          */
/* The bad news is that such trees are most efficient if the data to be     */
/* stored in the tree come in completely random (unsorted) order. This is   */
/* unfortunately not the case here because DISORT2 calls oceabrdf for       */
/* well sorted mu, mup, and dphi values. Hence, in a first step, the data   */
/* are gathered in a different tree (rtree) during the first DISORT2 call   */
/* for each wavelength (DISORT2 has been modified to pass the number of     */
/* the DISORT2 run since since the start of uvspec to oceanc). A function   */
/* rtree2tree then copies these data to their final destination in fully    */
/* randomized order.                                                        */
/*                                                                          */
/* A first test showed a speed improvement of a factor 3-4 for a            */
/* 16 streams DISORT2 run for the 32 KATO wavelength bands.                 */
/*                                                                          */
/* Removed this after nearly 15 years. It doesn't help at all with REPTRAN  */
/* where each wavelength is calculated only once and BRDFs may not be       */
/* recycled.                                                                */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "ocean.h"
#include "uvspecrandom.h"
#include "f77-uscore.h"
#include "fortran_and_c.h"
#include <math.h>
#include <complex.h>
#include "common_math.h"
#include "mystic.h"

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif

/* search tree node; this will become a 3D tree: "value" stores the  */
/* polar angle mu while "next" points to another tree where "value"  */
/* stores mup and "next" points to a third tree where "value" stores */
/* dphi and bdref finally holds the BRDF.                            */

struct node {
  float value;
  float bdref;

  struct node* left;
  struct node* right;
  struct node* next;
};

/* randomizing tree node, similar to the node above but one-dimensional */
struct rnode {
  float data[4];

  struct rnode* left;
  struct rnode* right;
};

/* prototypes */

static struct node* addtree (struct node* p, float x, float y, float z, float bdref);
static struct node* addtree2 (struct node* p, float y, float z, float bdref);
static struct node* addtree3 (struct node* p, float z, float bdref);

static float searchtree (struct node* p, float x, float y, float z);
static float searchtree2 (struct node* p, float y, float z);
static float searchtree3 (struct node* p, float z);

static void treefree (struct node* p);
static void treefree2 (struct node* p);
static void treefree3 (struct node* p);

struct rnode* addrtree (struct rnode* p, float mu, float mup, float dphi, float bdref);
static void   rtree2tree (struct rnode* source, struct node** dest);
static void   rtreefree (struct rnode* p);

/* prototype of the Fortran oceabrdf function */
void F77_FUNC (oceabrdf, OCEABRDF) (float* wvnmlo,
                                    float* wvnmhi,
                                    float* mu,
                                    float* mup,
                                    float* dphi,
                                    float* u10,
                                    float* pcl,
                                    float* xsal,
                                    float* paw,
                                    float* bdref);

void F77_FUNC (indwat, INDWAT) (float* lambda, float* dummy_sal, float* refre, float* refim);
void F77_FUNC (wrefice, WREFICE) (float* wl, float* temp, float* re, float* im);
void F77_FUNC (morcasiwat, MORCASIWAT) (float* wl, float* pcl, float* R);
void F77_FUNC (tdstdv, TDSTDV) (float* tetas, float* tetav, float* u10, float* tds, float* tdv);
void F77_FUNC (whitecaps, WHITECAPS) (float* wl, float* Ref);

/* simple C wrapper function for the standard oceabrdf call */
double oceabrdfc (float wvnmlo, float wvnmhi, float mu, float mup, float dphi, float u10, float pcl, float xsal, float paw) {
  float bdref = 0;

  F77_FUNC (oceabrdf, OCEABRDF) (&wvnmlo, &wvnmhi, &mu, &mup, &dphi, &u10, &pcl, &xsal, &paw, &bdref);
  return (double) bdref;
}

/* add a node to the search tree, see Kernighan and Ritchie */
static struct node* addtree (struct node* p, float x, float y, float z, float bdref) {
  if (p == NULL) {
    p        = calloc (1, sizeof (struct node));
    p->value = x;
    p->next  = addtree2 (p->next, y, z, bdref);
  } else {
    if (x == p->value) {
      p->next = addtree2 (p->next, y, z, bdref);
    } else {
      if (x < p->value)
        p->left = addtree (p->left, x, y, z, bdref);
      else
        p->right = addtree (p->right, x, y, z, bdref);
    }
  }
  return p;
}

/* add a node to the first subtree; required by addtree */
static struct node* addtree2 (struct node* p, float y, float z, float bdref) {
  if (p == NULL) {
    p        = calloc (1, sizeof (struct node));
    p->value = y;
    p->bdref = bdref;
    p->next  = addtree3 (p->next, z, bdref);
  } else {
    if (y == p->value) {
      p->next = addtree3 (p->next, z, bdref);
    } else {
      if (y < p->value)
        p->left = addtree2 (p->left, y, z, bdref);
      else
        p->right = addtree2 (p->right, y, z, bdref);
    }
  }
  return p;
}

/* add a node to the second subtree; required by addtree */
static struct node* addtree3 (struct node* p, float z, float bdref) {
  if (p == NULL) {
    p        = calloc (1, sizeof (struct node));
    p->value = z;
    p->bdref = bdref;
  } else {
    if (z < p->value)
      p->left = addtree3 (p->left, z, bdref);

    if (z > p->value)
      p->right = addtree3 (p->right, z, bdref);
  }
  return p;
}

/* search the search tree, see Kernighan and Ritchie */
static float searchtree (struct node* p, float x, float y, float z) {
  if (p == NULL) {
    return -999;
  } else {
    if (x == p->value) {
      return searchtree2 (p->next, y, z);
    } else {
      if (x < p->value)
        return searchtree (p->left, x, y, z);
      else
        return searchtree (p->right, x, y, z);
    }
  }
}

/* search the first subtree, required by searchtree */
static float searchtree2 (struct node* p, float y, float z) {
  if (p == NULL) {
    return -999;
  } else {
    if (y == p->value) {
      return searchtree3 (p->next, z);
    } else {
      if (y < p->value)
        return searchtree2 (p->left, y, z);
      else
        return searchtree2 (p->right, y, z);
    }
  }
}

/* search the second subtree, required by searchtree */
static float searchtree3 (struct node* p, float z) {
  if (p == NULL) {
    return -999;
  } else {
    if (z == p->value) {
      return p->bdref;
    } else {
      if (z < p->value)
        return searchtree3 (p->left, z);
      else
        return searchtree3 (p->right, z);
    }
  }
}

/* free randomization tree */
static void rtreefree (struct rnode* p) {
  if (p != NULL) {
    rtreefree (p->left);
    rtreefree (p->right);
    free (p);
  }
}

/* free search tree */
static void treefree (struct node* p) {
  if (p != NULL) {
    treefree (p->left);
    treefree (p->right);

    treefree2 (p->next);
    free (p);
  }
}

/* free first subtree, required by treefree */
static void treefree2 (struct node* p) {
  if (p != NULL) {
    treefree2 (p->left);
    treefree2 (p->right);

    treefree3 (p->next);
    free (p);
  }
}

/* free second subtree, required by treefree */
static void treefree3 (struct node* p) {
  if (p != NULL) {
    treefree3 (p->left);
    treefree3 (p->right);
    free (p);
  }
}

/* copy randomization tree to search tree */
static void rtree2tree (struct rnode* source, struct node** dest) {
  if (source != NULL) {
    rtree2tree (source->left, dest);
    *dest = addtree (*dest, source->data[0], source->data[1], source->data[2], source->data[3]);
    rtree2tree (source->right, dest);
  }
}

/* add node to the radomization tree; new nodes are added at random branches; */
/* that is, the tree is completely unsorted which is what we want             */
struct rnode* addrtree (struct rnode* p, float mu, float mup, float dphi, float bdref) {
  if (p == NULL) {
    p          = calloc (1, sizeof (struct rnode));
    p->data[0] = mu;
    p->data[1] = mup;
    p->data[2] = dphi;
    p->data[3] = bdref;
  } else {
    if (uvspec_random() < 0.5)
      p->left = addrtree (p->left, mu, mup, dphi, bdref);
    else
      p->right = addrtree (p->right, mu, mup, dphi, bdref);
  }
  return p;
}

/* prototype of the Fortran oceabrdf function */
void F77_FUNC (oceabrdf, OCEABRDF) (float* wvnmlo,
                                    float* wvnmhi,
                                    float* mu,
                                    float* mup,
                                    float* dphi,
                                    float* u10,
                                    float* pcl,
                                    float* xsal,
                                    float* paw,
                                    float* bdref);

/* new version without search tree, after November 30, 2016 */
float ocean_brdf (float wvnmlo, float wvnmhi, float mu, float mup, float dphi, float u10, float pcl, float xsal, int firstcall) {
  float paw;
  float bdref = 0.0;

  /* NOTE! paw = 0.0 means that wind comes from sun! This cannot be fixed! BCA */
  paw = 0.0;

  /* if anything (wavelength, salinity etc has changed, */
  /* we need to create a completely new search tree     */

  F77_FUNC (oceabrdf, OCEABRDF) (&wvnmlo, &wvnmhi, &mu, &mup, &dphi, &u10, &pcl, &xsal, &paw, &bdref);
  return bdref;
}

/* old version with search tree, before November 30, 2016 */
float ocean_brdf_old (float wvnmlo,
                      float wvnmhi,
                      float mu,
                      float mup,
                      float dphi,
                      float u10,
                      float pcl,
                      float xsal,
                      int   firstcall) {
  static struct node*  root    = NULL;
  static struct rnode* rtree   = NULL;
  static float         swvnmlo = 0, swvnmhi = 0, su10 = 0, spcl = 0, sxsal = 0;
  static int           randomize = 1, randomized = 0, randomrun = 0;
  static float         paw;
  float                bdref = 0.0;

  /* NOTE! paw = 0.0 means that wind comes from sun! This cannot be fixed! BCA */
  paw = 0.0;

  /* if anything (wavelength, salinity etc has changed, */
  /* we need to create a completely new search tree     */

  if (swvnmlo != wvnmlo || swvnmhi != wvnmhi || su10 != u10 || spcl != pcl || sxsal != xsal) {

    swvnmlo = wvnmlo;
    swvnmhi = wvnmhi;
    su10    = u10;
    spcl    = pcl;
    sxsal   = xsal;

    randomize  = 1;
    randomized = 0;
    randomrun  = firstcall;

    /* free memory of the search and randomization trees */
    rtreefree (rtree);
    rtree = NULL;

    treefree (root);
    root = NULL;
  }

  if (firstcall != randomrun)
    randomize = 0;

  /* add data to randomization tree during the first DISORT2 call */
  if (randomize) {
    F77_FUNC (oceabrdf, OCEABRDF) (&wvnmlo, &wvnmhi, &mu, &mup, &dphi, &u10, &pcl, &xsal, &paw, &bdref);

    rtree = addrtree (rtree, mu, mup, dphi, bdref);
    return bdref;
  }

  /* after creating the randomization tree, copy the data to the search tree */
  if (!randomized) {
    rtree2tree (rtree, &root);

    randomized = 1;
  }

  /* this is the main function of the program: */
  /* look if the data exist in the tree ...    */
  bdref = searchtree (root, mu, mup, dphi);

  /* ... if not, calculate them and store them in the search tree */
  if (bdref == -999) {

    F77_FUNC (oceabrdf, OCEABRDF) (&wvnmlo, &wvnmhi, &mu, &mup, &dphi, &u10, &pcl, &xsal, &paw, &bdref);

    root = addtree (root, mu, mup, dphi, bdref);
  }

  return bdref;
}

/***********************************************************************************/
/* Function: bpdf_tsang                                                            */
/* Description:                                                                    */
/*      Calculate reflection matrix (Tsang, 1985) using Mishchenko's code.         */
/*      The reflection matrix is for an ocean surface and depends on the wind      */
/*      speed. Shadowing of the waves is taken into account.                       */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Claudia Emde                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/
int bpdf_tsang (double    wind_speed,
                float     wavelength,
                double    mu_inc,
                double    phi_inc,
                double    mu_esc,
                double    phi_esc,
                double*** refl_mat) {

  int status = 0;

  ref_complex cn1;
  ref_complex cn2;

  double sigma2 = 0.0;

  float dummy_sal = -1;
  float refre     = 0.0;
  float refim     = 0.0;
  /*float T=288.200; surface temperature, only needed in Warren's */
  /* refractive index is used. */

  /* mean surface slope, should be calculated from wind speed 
     according to Eq. (18) Mishchenko and Travis (1997) */
  /* Numbers correspond to original Cox and Munk (1954) definition */
  sigma2 = 0.5 * (0.003 + 0.00512 * wind_speed);

  /* To be very accurate, cn1 should be set according to refractive*/
  /* index profile */
  cn1.re = 1.0;
  cn1.im = 0.0;

  /* wavelength dependent refractive index of water */

  /* this is used in Cox and Munk BRDF  */
  F77_FUNC (indwat, INDWAT) (&wavelength, &dummy_sal, &refre, &refim);

  /* This function uses the Warren database */
  /* F77_FUNC (indwat, WREFWAT) (&wavelength, &T, &refre, &refim); */
  cn2.re = refre;
  cn2.im = refim;

  /* Set values for refractive index, uncomment these lines if you want to do this */
  /* cn2.re=1.34; */
  /* cn2.im=0.0; */

  /*convention is clockwise in MYSTIC, anticlockwise in Mishchenkos program!*/
  /* phi must be replaced by -phi */
  phi_inc *= -1.;
  phi_esc *= -1.;

  /* to use the c translation */
  status = bpdf_tsang_rmatr (cn1, cn2, sigma2, mu_inc, phi_inc, mu_esc, phi_esc, 1, refl_mat);

  return status;
}

void F77_FUNC (oceanfort, OCEANFORT) (float* wvnmlo,
                                      float* wvnmhi,
                                      float* mu,
                                      float* mup,
                                      float* dphi,
                                      float* u10,
                                      float* pcl,
                                      float* xsal,
                                      int*   firstcall,
                                      float* bdref) {
  *bdref = (float) ocean_brdf (*wvnmlo, *wvnmhi, *mu, *mup, *dphi, *u10, *pcl, *xsal, *firstcall);
}

/***********************************************************************************/
/* Function: bpdf_tsang_rmatr                                                      */
/* Description:                                                                    */
/*      Calculate reflection matrix (Tsang, 1985). Translation of Mishchenko's     */
/*      fortran code. The reflection matrix is for an ocean surface and depends    */
/*      on the wind speed. Shadowing of the waves is taken into account.           */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Anna Weber                                                              */
/*                                                                        @i62_30@ */
/***********************************************************************************/
int bpdf_tsang_rmatr (ref_complex CN1,
                      ref_complex CN2,
                      double      sigma2,
                      double      mu_inc,
                      double      phi_inc,
                      double      mu_esc,
                      double      phi_esc,
                      int         shadowing,
                      double***   refl_mat) {

  int status = 0;
  int i = 0, j = 0;

  double dir_inc[3], dir_esc[3];
  double k_d[3];
  double h_inc[3], v_inc[3], h_esc[3], v_esc[3];

  double complex R_p, R_s;
  double complex f_hh, f_hv, f_vh, f_vv;

  double coeff, k_d_norm2, shadow;

  float complex cn1;
  float complex cn2;
  cn1 = CN1.re + CN1.im * I;
  cn2 = CN2.re + CN2.im * I;

  /* Calculate Stokes reflection matrix, translated from Mishchenko's fortran code in ocean.shadow.f */
  if (mu_inc > 0.999999999999)
    mu_inc = 0.999999999999;
  if (mu_esc > 0.999999999999)
    mu_esc = 0.999999999999;

  /* unit vectors of incident and reflected rays */
  double sitheta_inc = sqrt (1.0 - mu_inc * mu_inc);
  double sitheta_esc = sqrt (1.0 - mu_esc * mu_esc);

  dir_inc[0] = sitheta_inc * cos (phi_inc);
  dir_inc[1] = sitheta_inc * sin (phi_inc);
  dir_inc[2] = -mu_inc;

  dir_esc[0] = sitheta_esc * cos (phi_esc);
  dir_esc[1] = sitheta_esc * sin (phi_esc);
  dir_esc[2] = mu_esc;

  /* local surface normal k_d*/
  for (i = 0; i < 3; i++)
    k_d[i] = dir_inc[i] - dir_esc[i];
  k_d_norm2 = (k_d[0] * k_d[0] + k_d[1] * k_d[1] + k_d[2] * k_d[2]);

  /* fresnel reflection coefficients */
  double mu_inc_l = 1.0 / sqrt (k_d_norm2) *
                    (k_d[0] * dir_inc[0] + k_d[1] * dir_inc[1] + k_d[2] * dir_inc[2]);  //incident angle wrt local surface normal
  double complex mu_refr_l = sqrt (1.0 - (1.0 - mu_inc_l * mu_inc_l) * cn1 * cn1 / (cn2 * cn2));
  R_s                      = (cn1 * mu_inc_l - cn2 * mu_refr_l) / (cn1 * mu_inc_l + cn2 * mu_refr_l);
  R_p                      = (cn2 * mu_inc_l - cn1 * mu_refr_l) / (cn2 * mu_inc_l + cn1 * mu_refr_l);

  /* polarization vectors (h = phi, v = theta in Mishchenko's paper)*/
  h_inc[0] = -sin (phi_inc);
  h_inc[1] = cos (phi_inc);
  h_inc[2] = 0.0;

  v_inc[0] = -mu_inc * cos (phi_inc);
  v_inc[1] = -mu_inc * sin (phi_inc);
  v_inc[2] = -sitheta_inc;

  h_esc[0] = -sin (phi_esc);
  h_esc[1] = cos (phi_esc);
  h_esc[2] = 0.0;

  v_esc[0] = mu_esc * cos (phi_esc);
  v_esc[1] = mu_esc * sin (phi_esc);
  v_esc[2] = -sitheta_esc;

  /* amplitude transmission matrix */
  double h_inc_dir_esc = 0, h_esc_dir_inc = 0, v_inc_dir_esc = 0, v_esc_dir_inc = 0;
  for (i = 0; i < 3; i++) {
    h_inc_dir_esc += h_inc[i] * dir_esc[i];
    h_esc_dir_inc += h_esc[i] * dir_inc[i];
    v_inc_dir_esc += v_inc[i] * dir_esc[i];
    v_esc_dir_inc += v_esc[i] * dir_inc[i];
  }

  f_vv = h_inc_dir_esc * h_esc_dir_inc * R_s + v_inc_dir_esc * v_esc_dir_inc * R_p;
  f_vh = -v_inc_dir_esc * h_esc_dir_inc * R_s + h_inc_dir_esc * v_esc_dir_inc * R_p;
  f_hv = -h_inc_dir_esc * v_esc_dir_inc * R_s + v_inc_dir_esc * h_esc_dir_inc * R_p;
  f_hh = v_inc_dir_esc * v_esc_dir_inc * R_s + h_inc_dir_esc * h_esc_dir_inc * R_p;

  /* stokes transmission matrix */
  double dir_inc_cross_dir_esc[3];
  dir_inc_cross_dir_esc[0] = dir_inc[1] * dir_esc[2] - dir_inc[2] * dir_esc[1];
  dir_inc_cross_dir_esc[1] = dir_inc[2] * dir_esc[0] - dir_inc[0] * dir_esc[2];
  dir_inc_cross_dir_esc[2] = dir_inc[0] * dir_esc[1] - dir_inc[1] * dir_esc[0];
  double norm2             = 0;
  for (i = 0; i < 3; i++)
    norm2 += dir_inc_cross_dir_esc[i] * dir_inc_cross_dir_esc[i];

  coeff = (k_d_norm2 * k_d_norm2) / (norm2 * norm2 * k_d[2] * k_d[2] * k_d[2] * k_d[2]);
  coeff *= 1.0 / (4.0 * mu_inc * mu_esc);
  coeff *= 1.0 / (2.0 * sigma2) *
           exp (-(k_d[0] * k_d[0] + k_d[1] * k_d[1]) / (k_d[2] * k_d[2] * 2.0 * sigma2));  // slope probability distribution

  (*refl_mat)[0][0] =
    0.5 * coeff * (cabs (f_vv) * cabs (f_vv) + cabs (f_vh) * cabs (f_vh) + cabs (f_hv) * cabs (f_hv) + cabs (f_hh) * cabs (f_hh));
  (*refl_mat)[0][1] =
    0.5 * coeff * (cabs (f_vv) * cabs (f_vv) - cabs (f_vh) * cabs (f_vh) + cabs (f_hv) * cabs (f_hv) - cabs (f_hh) * cabs (f_hh));
  (*refl_mat)[1][0] =
    0.5 * coeff * (cabs (f_vv) * cabs (f_vv) + cabs (f_vh) * cabs (f_vh) - cabs (f_hv) * cabs (f_hv) - cabs (f_hh) * cabs (f_hh));
  (*refl_mat)[1][1] =
    0.5 * coeff * (cabs (f_vv) * cabs (f_vv) - cabs (f_vh) * cabs (f_vh) - cabs (f_hv) * cabs (f_hv) + cabs (f_hh) * cabs (f_hh));
  (*refl_mat)[0][2] = -coeff * creal (f_vv * conj (f_vh) + f_hv * conj (f_hh));
  (*refl_mat)[0][3] = -coeff * cimag (f_vv * conj (f_vh) + f_hv * conj (f_hh));
  (*refl_mat)[1][2] = -coeff * creal (f_vv * conj (f_vh) - f_hv * conj (f_hh));
  (*refl_mat)[1][3] = -coeff * cimag (f_vv * conj (f_vh) - f_hv * conj (f_hh));
  (*refl_mat)[2][0] = -coeff * creal (f_vv * conj (f_hv) + f_vh * conj (f_hh));
  (*refl_mat)[2][1] = -coeff * creal (f_vv * conj (f_hv) - f_vh * conj (f_hh));
  (*refl_mat)[2][2] = coeff * creal (f_vv * conj (f_hh) + f_vh * conj (f_hv));
  (*refl_mat)[2][3] = coeff * cimag (f_vv * conj (f_hh) - f_vh * conj (f_hv));
  (*refl_mat)[3][0] = coeff * cimag (f_vv * conj (f_hv) + f_vh * conj (f_hh));
  (*refl_mat)[3][1] = coeff * cimag (f_vv * conj (f_hv) - f_vh * conj (f_hh));
  (*refl_mat)[3][2] = -coeff * cimag (f_vv * conj (f_hh) + f_vh * conj (f_hv));
  (*refl_mat)[3][3] = coeff * creal (f_vv * conj (f_hh) - f_vh * conj (f_hv));

  /* shadowing */
  if (shadowing) {
    shadow = 1.0 / (1.0 + gamma_shadowing (mu_inc, sigma2) + gamma_shadowing (mu_esc, sigma2));
    for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++)
        (*refl_mat)[i][j] *= shadow;
    }
  }

  return status;
}

double gamma_shadowing (double mu, double sigma2) {

  double cot = mu / sqrt (1 - mu * mu);
  return 0.5 * (sqrt (2.0 * (1.0 - mu * mu) / PI) * sqrt (sigma2) / mu * exp (-cot * cot / (2.0 * sigma2)) -
                erfc (cot / sqrt (2 * sigma2)));
}

/***********************************************************************************/
/* Function: fresnel_rmatr                                                         */
/* Description:                                                                    */
/*      Calculate Fresnel reflection matrix (Tsang, 1985). Translation of          */
/*      Mishchenko's fortran code.                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Anna Weber                                                              */
/*                                                                        @i62_30@ */
/***********************************************************************************/
int fresnel_rmatr (ref_complex CN1,
                   ref_complex CN2,
                   double      mu_inc,
                   double      phi_inc,
                   double      mu_esc,
                   double      phi_esc,
                   double***   refl_mat) {

  int status = 0;
  int i      = 0;

  double dir_inc[3], dir_esc[3];
  double k_d[3];
  double h_inc[3], v_inc[3], h_esc[3], v_esc[3];

  double complex R_p, R_s;
  double complex f_hh, f_hv, f_vh, f_vv;

  double coeff, k_d_norm2;

  float complex cn1;
  float complex cn2;
  cn1 = CN1.re + CN1.im * I;
  cn2 = CN2.re + CN2.im * I;

  /* Calculate Stokes reflection matrix, translated from Mishchenko's fortran code in ocean.shadow.f */
  if (mu_inc > 0.999999999999)
    mu_inc = 0.999999999999;
  if (mu_esc > 0.999999999999)
    mu_esc = 0.999999999999;

  /* unit vectors of incident and reflected rays */
  double sitheta_inc = sqrt (1.0 - mu_inc * mu_inc);
  double sitheta_esc = sqrt (1.0 - mu_esc * mu_esc);

  dir_inc[0] = sitheta_inc * cos (phi_inc);
  dir_inc[1] = sitheta_inc * sin (phi_inc);
  dir_inc[2] = -mu_inc;

  dir_esc[0] = sitheta_esc * cos (phi_esc);
  dir_esc[1] = sitheta_esc * sin (phi_esc);
  dir_esc[2] = mu_esc;

  /* local surface normal k_d*/
  for (i = 0; i < 3; i++)
    k_d[i] = dir_inc[i] - dir_esc[i];
  k_d_norm2 = (k_d[0] * k_d[0] + k_d[1] * k_d[1] + k_d[2] * k_d[2]);

  /* fresnel reflection coefficients */
  double mu_inc_l = 1.0 / sqrt (k_d_norm2) *
                    (k_d[0] * dir_inc[0] + k_d[1] * dir_inc[1] + k_d[2] * dir_inc[2]);  //incident angle wrt local surface normal
  double complex mu_refr_l = sqrt (1.0 - (1.0 - mu_inc_l * mu_inc_l) * cn1 * cn1 / (cn2 * cn2));
  R_s                      = (cn1 * mu_inc_l - cn2 * mu_refr_l) / (cn1 * mu_inc_l + cn2 * mu_refr_l);
  R_p                      = (cn2 * mu_inc_l - cn1 * mu_refr_l) / (cn2 * mu_inc_l + cn1 * mu_refr_l);

  /* polarization vectors (h = phi, v = theta in Mishchenko's paper)*/
  h_inc[0] = -sin (phi_inc);
  h_inc[1] = cos (phi_inc);
  h_inc[2] = 0.0;

  v_inc[0] = -mu_inc * cos (phi_inc);
  v_inc[1] = -mu_inc * sin (phi_inc);
  v_inc[2] = -sitheta_inc;

  h_esc[0] = -sin (phi_esc);
  h_esc[1] = cos (phi_esc);
  h_esc[2] = 0.0;

  v_esc[0] = mu_esc * cos (phi_esc);
  v_esc[1] = mu_esc * sin (phi_esc);
  v_esc[2] = -sitheta_esc;

  /* amplitude reflection matrix */
  double h_inc_dir_esc = 0, h_esc_dir_inc = 0, v_inc_dir_esc = 0, v_esc_dir_inc = 0;
  for (i = 0; i < 3; i++) {
    h_inc_dir_esc += h_inc[i] * dir_esc[i];
    h_esc_dir_inc += h_esc[i] * dir_inc[i];
    v_inc_dir_esc += v_inc[i] * dir_esc[i];
    v_esc_dir_inc += v_esc[i] * dir_inc[i];
  }

  f_vv = h_inc_dir_esc * h_esc_dir_inc * R_s + v_inc_dir_esc * v_esc_dir_inc * R_p;
  f_vh = -v_inc_dir_esc * h_esc_dir_inc * R_s + h_inc_dir_esc * v_esc_dir_inc * R_p;
  f_hv = -h_inc_dir_esc * v_esc_dir_inc * R_s + v_inc_dir_esc * h_esc_dir_inc * R_p;
  f_hh = v_inc_dir_esc * v_esc_dir_inc * R_s + h_inc_dir_esc * h_esc_dir_inc * R_p;

  /* stokes reflection matrix */
  double dir_inc_cross_dir_esc[3];
  dir_inc_cross_dir_esc[0] = dir_inc[1] * dir_esc[2] - dir_inc[2] * dir_esc[1];
  dir_inc_cross_dir_esc[1] = dir_inc[2] * dir_esc[0] - dir_inc[0] * dir_esc[2];
  dir_inc_cross_dir_esc[2] = dir_inc[0] * dir_esc[1] - dir_inc[1] * dir_esc[0];
  double norm2             = 0;
  for (i = 0; i < 3; i++)
    norm2 += dir_inc_cross_dir_esc[i] * dir_inc_cross_dir_esc[i];

  coeff = 1 / (norm2 * norm2);  //fresnel

  (*refl_mat)[0][0] =
    0.5 * coeff * (cabs (f_vv) * cabs (f_vv) + cabs (f_vh) * cabs (f_vh) + cabs (f_hv) * cabs (f_hv) + cabs (f_hh) * cabs (f_hh));
  (*refl_mat)[0][1] =
    0.5 * coeff * (cabs (f_vv) * cabs (f_vv) - cabs (f_vh) * cabs (f_vh) + cabs (f_hv) * cabs (f_hv) - cabs (f_hh) * cabs (f_hh));
  (*refl_mat)[1][0] =
    0.5 * coeff * (cabs (f_vv) * cabs (f_vv) + cabs (f_vh) * cabs (f_vh) - cabs (f_hv) * cabs (f_hv) - cabs (f_hh) * cabs (f_hh));
  (*refl_mat)[1][1] =
    0.5 * coeff * (cabs (f_vv) * cabs (f_vv) - cabs (f_vh) * cabs (f_vh) - cabs (f_hv) * cabs (f_hv) + cabs (f_hh) * cabs (f_hh));
  (*refl_mat)[0][2] = -coeff * creal (f_vv * conj (f_vh) + f_hv * conj (f_hh));
  (*refl_mat)[0][3] = -coeff * cimag (f_vv * conj (f_vh) + f_hv * conj (f_hh));
  (*refl_mat)[1][2] = -coeff * creal (f_vv * conj (f_vh) - f_hv * conj (f_hh));
  (*refl_mat)[1][3] = -coeff * cimag (f_vv * conj (f_vh) - f_hv * conj (f_hh));
  (*refl_mat)[2][0] = -coeff * creal (f_vv * conj (f_hv) + f_vh * conj (f_hh));
  (*refl_mat)[2][1] = -coeff * creal (f_vv * conj (f_hv) - f_vh * conj (f_hh));
  (*refl_mat)[2][2] = coeff * creal (f_vv * conj (f_hh) + f_vh * conj (f_hv));
  (*refl_mat)[2][3] = coeff * cimag (f_vv * conj (f_hh) - f_vh * conj (f_hv));
  (*refl_mat)[3][0] = coeff * cimag (f_vv * conj (f_hv) + f_vh * conj (f_hh));
  (*refl_mat)[3][1] = coeff * cimag (f_vv * conj (f_hv) - f_vh * conj (f_hh));
  (*refl_mat)[3][2] = -coeff * cimag (f_vv * conj (f_hh) + f_vh * conj (f_hv));
  (*refl_mat)[3][3] = coeff * creal (f_vv * conj (f_hh) - f_vh * conj (f_hv));

  return status;
}

/***********************************************************************************/
/* Function: cox_and_munk_polarized                                                */
/* Description:                                                                    */
/*      C translation of the Cox and Munk BRDF (oceabrdf.f),                       */
/*      extension to polarization (Fresnel coefficient reflaced by Fresnel matrix) */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Anna Weber                                                              */
/*                                                                        @i62_30@ */
/***********************************************************************************/
int cox_and_munk_polarized (float     wavelength,  // wavelength in micrometer
                            float     mu_i,        // cos(incident angle)
                            float     mu_r,        // cos(reflected angle)
                            float     phi_i,       // incident azimuth in radians
                            float     phi_r,       // reflected azimuth in radians
                            float     u10,         // wind speed
                            float     pcl,         // pigment concentration
                            float     xsal,        // salinity
                            float     azw,         // wind direction in degree
                            int       backward,
                            double*** refl_mat) {
  int status = 0;
  int i, j;

  ref_complex cn1;
  ref_complex cn2;

  float refre = 0.0;
  float refim = 0.0;

  float sitheta_i, sitheta_r;
  float sigmaC, sigmaU, C21, C03, C40, C22, C04;
  float phi, zx, zy, tilt, xe, xn, coef, proba;
  float W, R_wc, a, R_w, tds, tdv;

  if (u10 < 0.1)
    u10 = 0.1;

  azw *= PI / 180;

  if (wavelength < 0.4)
    wavelength = 0.4;
  else if (wavelength > 4.)
    wavelength = 4.0;

  /* To be very accurate, cn1 should be set according to refractive*/
  /* index profile */
  cn1.re = 1.0;
  cn1.im = 0.0;

  /* wavelength dependent refractive index of water */
  /* this is used in Cox and Munk BRDF  */
  F77_FUNC (indwat, INDWAT) (&wavelength, &xsal, &refre, &refim);

  cn2.re = refre;
  cn2.im = refim;

  /* sunglint contribution */
  /* anisotropic Gaussian probability distribution according to Cox and Munk 1954*/
  sigmaC = 0.003 + 0.00192 * u10;
  sigmaU = 0.00316 * u10;
  C21    = 0.01 - 0.0086 * u10;
  C03    = 0.04 - 0.033 * u10;
  C40    = 0.40;
  C22    = 0.12;
  C04    = 0.23;

  phi       = PI - (phi_r - phi_i);
  sitheta_i = sqrt (1.0 - mu_i * mu_i);
  sitheta_r = sqrt (1.0 - mu_r * mu_r);

  zx   = -sitheta_r * sin (phi) / (mu_i + mu_r);
  zy   = (sitheta_i + sitheta_r * cos (phi)) / (mu_i + mu_r);
  tilt = atan (sqrt (zx * zx + zy * zy));

  xe = (cos (azw) * zx + sin (azw) * zy) / sqrt (sigmaC);
  xn = (-sin (azw) * zx + cos (azw) * zy) / sqrt (sigmaU);

  coef = 1 - C21 / 2. * (xe * xe - 1) * xn - C03 / 6. * (xn * xn - 3) * xn;
  coef += C40 / 24. * (xe * xe * xe * xe - 6 * xe * xe + 3);
  coef += C04 / 24. * (xn * xn * xn * xn - 6 * xn * xn + 3);
  coef += C22 / 4. * (xe * xe - 1) * (xn * xn - 1);
  proba = coef / (2. * PI * sqrt (sigmaU) * sqrt (sigmaC)) * exp (-(xe * xe + xn * xn) / 2.);

  /* fresnel reflection matrix*/
  if (backward)
    status = fresnel_rmatr (cn1,
                            cn2,
                            mu_r,
                            phi_r * (-1),
                            mu_i,
                            phi_i * (-1),
                            refl_mat);  // convention for phi in mystic different from convention for Fresnel matrix
  else
    status = fresnel_rmatr (cn1,
                            cn2,
                            mu_i,
                            phi_i * (-1),
                            mu_r,
                            phi_r * (-1),
                            refl_mat);  // convention for phi in mystic different from convention for Fresnel matrix

  /* reflectance of the sunglint */
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++)
      (*refl_mat)[i][j] *= proba / (cos (tilt) * cos (tilt) * cos (tilt) * cos (tilt)) * 1 / (4 * mu_i * mu_r) * PI;
  }

  /* whitecaps reflectance (Koepke 1984)*/
  W = 2.95e-06 * pow (u10, 3.52);
  F77_FUNC (whitecaps, WHITECAPS) (&wavelength, &R_wc);
  R_wc *= W;

  /* backscattered sea water reflectance from below the sea surface (Morel 1988) */
  F77_FUNC (morcasiwat, MORCASIWAT) (&wavelength, &pcl, &R_w);

  float theta_i = acos (mu_i) * 180 / PI;
  float theta_r = acos (mu_r) * 180 / PI;
  F77_FUNC (tdstdv, TDSTDV) (&theta_i, &theta_r, &u10, &tds, &tdv);

  /* water reflectance above the sea surface, see oceabrdf.f */
  a   = 0.485;
  R_w = (1 / (cn2.re * cn2.re + cn2.im * cn2.im)) * tds * tdv * R_w / (1 - a * R_w);

  /* total reflectance (assume that whitecaps and from below sea surface contributions are unpolarized) */
  (*refl_mat)[0][0] = R_wc + (1 - W) * (*refl_mat)[0][0] + (1 - R_wc) * R_w;

  return status;
}

/***********************************************************************************/
/* Function: bpdf_tsang_ice                                                        */
/* Description:                                                                    */
/*      Calculate reflection matrix (Tsang, 1985) using Mishchenko's code.         */
/*      The reflection matrix depends on the width of the angular distribution.    */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/
int bpdf_tsang_ice (double      angdist_width,
                    ref_complex CN1,
                    ref_complex CN2,
                    double      mu_inc,
                    double      phi_inc,
                    double      mu_esc,
                    double      phi_esc,
                    double***   refl_mat) {

  int    status = 0;
  double sigma2 = 0.0;

  /* Calculate mean squared surface slope from angular distribution width*/
  sigma2 = tan (angdist_width / 180 * PI) * tan (angdist_width / 180 * PI);

  /*convention is clockwise in MYSTIC, anticlockwise in Mishchenkos program!*/
  /* phi must be replaced by -phi */
  phi_inc *= -1.;
  phi_esc *= -1.;

  status = bpdf_tsang_rmatr (CN1, CN2, sigma2, mu_inc, phi_inc, mu_esc, phi_esc, 0, refl_mat);

  return status;
}

/***********************************************************************************/
/* Function: bpdf_tsang_tmatr                                                      */
/* Description:                                                                    */
/*      Calculate stokes transmission matrix (Tsang, 1985) for a horizontal        */
/*      surface in the Kirchhoff approximation (without wave contribution)         */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Anna Weber                                                              */
/*                                                                        @i62_30@ */
/***********************************************************************************/
int bpdf_tsang_tmatr (ref_complex CN1,
                      ref_complex CN2,
                      double      mu_inc,
                      double      phi_inc,
                      double      mu_esc,
                      double      phi_esc,
                      double***   t_mat) {
  int status = 0;

  int i = 0;

  double dir_inc[3], dir_esc[3];
  double k_d[3];
  double h_inc[3], v_inc[3], h_esc[3], v_esc[3];

  double complex T_p, T_s;
  double complex f_hh, f_hv, f_vh, f_vv;

  double coeff, k_d_norm2;

  float complex cn1;
  float complex cn2;
  cn1 = CN1.re + CN1.im * I;
  cn2 = CN2.re + CN2.im * I;

  /*convention is clockwise in MYSTIC, anticlockwise in Mishchenkos program!*/
  /* phi must be replaced by -phi */
  phi_inc *= -1.;
  phi_esc *= -1.;

  /* Calculate Stokes transmission matrix analogous to reflection matrix in Mishchenkos program but without wave contribution */
  if (mu_inc > 0.9999999)
    mu_inc = 0.9999999;
  else if (mu_esc > 0.9999999)
    mu_esc = 0.9999999;

  /* unit vectors of incident and transmitted rays */
  double sitheta_inc = sqrt (1.0 - mu_inc * mu_inc);
  double sitheta_esc = sqrt (1.0 - mu_esc * mu_esc);

  dir_inc[0] = sitheta_inc * cos (phi_inc);
  dir_inc[1] = sitheta_inc * sin (phi_inc);
  dir_inc[2] = -mu_inc;

  dir_esc[0] = sitheta_esc * cos (phi_esc);
  dir_esc[1] = sitheta_esc * sin (phi_esc);
  dir_esc[2] = -mu_esc;

  /* local surface normal */
  for (i = 0; i < 3; i++) {
    if (creal (cn1) < creal (cn2))
      k_d[i] = (creal (cn1) * dir_inc[i] - creal (cn2) * dir_esc[i]);
    else
      k_d[i] = -(creal (cn1) * dir_inc[i] - creal (cn2) * dir_esc[i]);
  }
  k_d_norm2 = (k_d[0] * k_d[0] + k_d[1] * k_d[1] + k_d[2] * k_d[2]);

  /* fresnel transmission coefficients */
  double mu_inc_l =
    -1.0 / sqrt (k_d_norm2) *
    (k_d[0] * dir_inc[0] + k_d[1] * dir_inc[1] +
     k_d[2] * dir_inc[2]);  //incident angle wrt local surface normal, for horizontal surface/specular transmission = mu_inc
  double mu_esc_l = -1.0 / sqrt (k_d_norm2) *
                    (k_d[0] * dir_esc[0] + k_d[1] * dir_esc[1] +
                     k_d[2] * dir_esc[2]);  //escape angle wrt local surface normal, for horizontal surface = mu_esc

  T_s = 2.0 * cn1 * mu_inc_l / (cn1 * mu_inc_l + cn2 * mu_esc_l);
  T_p = 2.0 * cn1 * mu_inc_l / (cn2 * mu_inc_l + cn1 * mu_esc_l);

  /* total reflexion */
  if ((creal (cn1) / creal (cn2) * sqrt (1.0 - mu_inc_l * mu_inc_l)) > 1.0) {
    T_s = 0.0;
    T_p = 0.0;
  }

  /* polarization vectors */
  h_inc[0] = -sin (phi_inc);
  h_inc[1] = cos (phi_inc);
  h_inc[2] = 0.0;

  v_inc[0] = -mu_inc * cos (phi_inc);
  v_inc[1] = -mu_inc * sin (phi_inc);
  v_inc[2] = -sitheta_inc;

  h_esc[0] = -sin (phi_esc);
  h_esc[1] = cos (phi_esc);
  h_esc[2] = 0.0;

  v_esc[0] = -mu_esc * cos (phi_esc);
  v_esc[1] = -mu_esc * sin (phi_esc);
  v_esc[2] = -sitheta_esc;

  /* amplitude transmission matrix */
  double h_inc_dir_esc = 0, h_esc_dir_inc = 0, v_inc_dir_esc = 0, v_esc_dir_inc = 0;
  for (i = 0; i < 3; i++) {
    h_inc_dir_esc += h_inc[i] * dir_esc[i];
    h_esc_dir_inc += h_esc[i] * dir_inc[i];
    v_inc_dir_esc += v_inc[i] * dir_esc[i];
    v_esc_dir_inc += v_esc[i] * dir_inc[i];
  }

  f_vv = h_inc_dir_esc * h_esc_dir_inc * T_s + v_inc_dir_esc * v_esc_dir_inc * T_p;
  f_vh = -v_inc_dir_esc * h_esc_dir_inc * T_s + h_inc_dir_esc * v_esc_dir_inc * T_p;
  f_hv = -h_inc_dir_esc * v_esc_dir_inc * T_s + v_inc_dir_esc * h_esc_dir_inc * T_p;
  f_hh = v_inc_dir_esc * v_esc_dir_inc * T_s + h_inc_dir_esc * h_esc_dir_inc * T_p;

  /* stokes transmission matrix */
  double dir_inc_cross_dir_esc[3];
  dir_inc_cross_dir_esc[0] = dir_inc[1] * dir_esc[2] - dir_inc[2] * dir_esc[1];
  dir_inc_cross_dir_esc[1] = dir_inc[2] * dir_esc[0] - dir_inc[0] * dir_esc[2];
  dir_inc_cross_dir_esc[2] = dir_inc[0] * dir_esc[1] - dir_inc[1] * dir_esc[0];
  double norm2             = 0;
  for (i = 0; i < 3; i++)
    norm2 += dir_inc_cross_dir_esc[i] * dir_inc_cross_dir_esc[i];

  coeff = 1 / (norm2 * norm2) * creal (cn2) / creal (cn1) * mu_esc_l / mu_inc_l;  // fresnel
  //coeff *= 1.0/(4.0*mu_inc*mu_esc); // for radiance/BTDF
  /* to include slope probability distribution: */
  //coeff *= 1.0 / (2.0 * sigma2) * exp(-(k_d[0] * k_d[0] + k_d[1] * k_d[1]) / (k_d[2] * k_d[2] * 2.0 * sigma2) * (k_d_norm2 * k_d_norm2) / ( k_d[2] * k_d[2] * k_d[2] * k_d[2]);

  (*t_mat)[0][0] =
    0.5 * coeff * (cabs (f_vv) * cabs (f_vv) + cabs (f_vh) * cabs (f_vh) + cabs (f_hv) * cabs (f_hv) + cabs (f_hh) * cabs (f_hh));
  (*t_mat)[0][1] =
    0.5 * coeff * (cabs (f_vv) * cabs (f_vv) - cabs (f_vh) * cabs (f_vh) + cabs (f_hv) * cabs (f_hv) - cabs (f_hh) * cabs (f_hh));
  (*t_mat)[1][0] =
    0.5 * coeff * (cabs (f_vv) * cabs (f_vv) + cabs (f_vh) * cabs (f_vh) - cabs (f_hv) * cabs (f_hv) - cabs (f_hh) * cabs (f_hh));
  (*t_mat)[1][1] =
    0.5 * coeff * (cabs (f_vv) * cabs (f_vv) - cabs (f_vh) * cabs (f_vh) - cabs (f_hv) * cabs (f_hv) + cabs (f_hh) * cabs (f_hh));
  (*t_mat)[0][2] = -coeff * creal (f_vv * conj (f_vh) + f_hv * conj (f_hh));
  (*t_mat)[0][3] = -coeff * cimag (f_vv * conj (f_vh) + f_hv * conj (f_hh));
  (*t_mat)[1][2] = -coeff * creal (f_vv * conj (f_vh) - f_hv * conj (f_hh));
  (*t_mat)[1][3] = -coeff * cimag (f_vv * conj (f_vh) - f_hv * conj (f_hh));
  (*t_mat)[2][0] = -coeff * creal (f_vv * conj (f_hv) + f_vh * conj (f_hh));
  (*t_mat)[2][1] = -coeff * creal (f_vv * conj (f_hv) - f_vh * conj (f_hh));
  (*t_mat)[2][2] = coeff * creal (f_vv * conj (f_hh) + f_vh * conj (f_hv));
  (*t_mat)[2][3] = coeff * cimag (f_vv * conj (f_hh) - f_vh * conj (f_hv));
  (*t_mat)[3][0] = coeff * cimag (f_vv * conj (f_hv) + f_vh * conj (f_hh));
  (*t_mat)[3][1] = coeff * cimag (f_vv * conj (f_hv) - f_vh * conj (f_hh));
  (*t_mat)[3][2] = -coeff * cimag (f_vv * conj (f_hh) + f_vh * conj (f_hv));
  (*t_mat)[3][3] = coeff * creal (f_vv * conj (f_hh) - f_vh * conj (f_hv));

  return status;
}

/***********************************************************************************/
/* Function: bpdf_tsang_ice_internal_reflection                                    */
/* Description:                                                                    */
/*      Calculate a combined reflection matrix for transmission, internal          */
/*      reflection and transmission at two parallel interfaces of e.g. an ice      */
/*      crystal plate in the Kirchhoff approximation similar to the bpdf function  */
/*      by Tsang.                                                                  */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Anna Weber                                                              */
/*                                                                        @i62_30@ */
/***********************************************************************************/
int bpdf_tsang_ice_internal_reflection (double      angdist_width,
                                        ref_complex CN1,
                                        ref_complex CN2,
                                        double      mu_inc,
                                        double      phi_inc,
                                        double      mu_esc,
                                        double      phi_esc,
                                        int         backward,
                                        double***   refl_mat) {

  int status = 0;
  int i = 0, j = 0, k = 0;

  double dir_inc[3], dir_esc[3], dir_refr_1[3], dir_refr_2[3];
  double k_d[3];
  double k_d_norm2;
  double n_z[3] = {0, 0, 1};

  double mu_inc_l, mu_refr_l, mu_refr_1, phi_refr_1, mu_refr_2, phi_refr_2;
  double rotation_angle;
  double rotation_axis[3];

  double** t_mat             = NULL;
  double** internal_refl_mat = NULL;
  double** t_r_mat           = NULL;

  /* Calculate a combined Stokes reflection matrix for transmission, internal reflection and transmission at two parallel interfaces of e.g. an ice crystal plate */
  if (mu_inc > 0.999999999999)
    mu_inc = 0.999999999999;
  if (mu_esc > 0.999999999999)
    mu_esc = 0.999999999999;

  /* unit vectors of incident and reflected rays (MYSTIC convention) */
  double sitheta_inc = sqrt (1.0 - mu_inc * mu_inc);
  double sitheta_esc = sqrt (1.0 - mu_esc * mu_esc);

  dir_inc[0] = sitheta_inc * sin (phi_inc);
  dir_inc[1] = sitheta_inc * cos (phi_inc);
  dir_inc[2] = -mu_inc;

  dir_esc[0] = sitheta_esc * sin (phi_esc);
  dir_esc[1] = sitheta_esc * cos (phi_esc);
  dir_esc[2] = mu_esc;

  /* local surface normal k_d*/
  for (i = 0; i < 3; i++)
    k_d[i] = dir_inc[i] - dir_esc[i];
  k_d_norm2 = (k_d[0] * k_d[0] + k_d[1] * k_d[1] + k_d[2] * k_d[2]);

  /************* 1. transmission *****************/
  /* local incident angle and local angle of refraction */
  mu_inc_l  = 1.0 / sqrt (k_d_norm2) * (k_d[0] * dir_inc[0] + k_d[1] * dir_inc[1] + k_d[2] * dir_inc[2]);
  mu_refr_l = sqrt (1.0 - (1.0 - mu_inc_l * mu_inc_l) * CN1.re * CN1.re / (CN2.re * CN2.re));

  /* refracted photon direction */
  rotation_angle = acos (mu_inc_l) - acos (mu_refr_l);
  status         = v_cross_product_norm (dir_inc, dir_esc, rotation_axis);
  status         = rodrigues (dir_inc, rotation_axis, -rotation_angle, dir_refr_1);

  status     = v_mult_mu (n_z, dir_refr_1, &mu_refr_1);
  mu_refr_1  = -mu_refr_1;
  phi_refr_1 = calc_phi_horz (dir_refr_1, n_z) / 180 * PI;

  /* transmission matrix */
  ASCII_calloc_double (&t_mat, 4, 4);
  status = bpdf_tsang_tmatr (CN1, CN2, mu_inc, phi_inc, mu_refr_1, phi_refr_1, &t_mat);

  /*********** internal reflection ****************/
  /* reflected photon direction */
  status = rodrigues (dir_esc, rotation_axis, rotation_angle, dir_refr_2);

  status     = v_mult_mu (n_z, dir_refr_2, &mu_refr_2);
  phi_refr_2 = calc_phi_horz (dir_refr_2, n_z) / 180 * PI;

  /* reflection matrix (includes surface slope probability) */
  ASCII_calloc_double (&internal_refl_mat, 4, 4);
  status = bpdf_tsang_ice (angdist_width, CN2, CN1, mu_refr_1, phi_refr_1, mu_refr_2, phi_refr_2, &internal_refl_mat);

  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++)
      internal_refl_mat[i][j] *= (mu_refr_1 * mu_refr_2) / (mu_inc * mu_esc);  // correct the prefactor for combined matrix
  }

  /* matrix multiplication */
  ASCII_calloc_double (&t_r_mat, 4, 4);
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      t_r_mat[i][j] = 0.0;
      for (k = 0; k < 4; k++) {
        if (backward)
          t_r_mat[i][j] += t_mat[i][k] * internal_refl_mat[k][j];
        else
          t_r_mat[i][j] += internal_refl_mat[i][k] * t_mat[k][j];
      }
    }
  }

  /************* 2. transmission ******************/
  /* transmission matrix */
  status = bpdf_tsang_tmatr (CN2, CN1, mu_refr_2, phi_refr_2, mu_esc, phi_esc, &t_mat);

  /* matrix multiplication */
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      (*refl_mat)[i][j] = 0.0;
      for (k = 0; k < 4; k++) {
        if (backward)
          (*refl_mat)[i][j] += t_r_mat[i][k] * t_mat[k][j];
        else
          (*refl_mat)[i][j] += t_mat[i][k] * t_r_mat[k][j];
      }
    }
  }

  ASCII_free_double (t_mat, 4);
  ASCII_free_double (internal_refl_mat, 4);
  ASCII_free_double (t_r_mat, 4);

  return status;
}

int rodrigues (const double vector[3], const double axis[3], double angle, double rotated[3]) {

  /* rotate vector around axis by angle */
  double axis_cross_vector[3];
  double axis_dot_vector;
  int    i = 0;

  v_cross_product (axis, vector, axis_cross_vector);
  axis_dot_vector = vec3_dot_product (axis, vector);

  for (i = 0; i < 3; i++)
    rotated[i] = vector[i] * cos (angle) + axis_cross_vector[i] * sin (angle) + axis[i] * axis_dot_vector * (1 - cos (angle));

  return 0;
}

/***********************************************************************************/
/* Function: bpdf_tsang_ice_transmission                                           */
/* Description:                                                                    */
/*      Calculate a transmission matrix for transmission through two parallel and  */
/*      horizontally oriented interfaces, e.g. a horizontally oriented ice crystal */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Anna Weber                                                              */
/*                                                                        @i62_30@ */
/***********************************************************************************/
int bpdf_tsang_ice_transmission (ref_complex cn1, ref_complex cn2, double mu_inc, double phi_inc, int backward, double*** t_mat) {

  double   mu_refr;
  double** t_mat1 = NULL;
  double** t_mat2 = NULL;
  int      i = 0, j = 0, k = 0;
  int      status = 0;

  /* angle of refraction */
  mu_refr = sqrt (1 - ((cn1.re * cn1.re + cn1.im * cn1.im) / (cn2.re * cn2.re + cn2.im * cn2.im) * (1 - mu_inc * mu_inc)));
  if (mu_inc < 0)
    mu_refr *= -1.0;

  /* 1. transmission */
  ASCII_calloc_double (&t_mat1, 4, 4);
  status = bpdf_tsang_tmatr (cn1, cn2, mu_inc, phi_inc, mu_refr, phi_inc, &t_mat1);

  /* 2. transmission */
  ASCII_calloc_double (&t_mat2, 4, 4);
  status = bpdf_tsang_tmatr (cn2, cn1, mu_refr, phi_inc, mu_inc, phi_inc, &t_mat2);

  /* matrix multiplication */
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      (*t_mat)[i][j] = 0.0;
      for (k = 0; k < 4; k++) {
        if (backward)
          (*t_mat)[i][j] += t_mat1[i][k] * t_mat2[k][j];
        else
          (*t_mat)[i][j] += t_mat2[i][k] * t_mat1[k][j];
      }
    }
  }

  ASCII_free_double (t_mat1, 4);
  ASCII_free_double (t_mat2, 4);

  return status;
}
