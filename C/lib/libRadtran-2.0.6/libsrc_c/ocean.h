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

#ifndef __ocean_h
#define __ocean_h

#if defined(__cplusplus)
extern "C" {
#endif

#include "f77-uscore.h"

/* prototypes */

float ocean_brdf (float wvnmlo, float wvnmhi, float mu, float mup, float dphi, float u10, float pcl, float xsal, int firstcall);

/* Fortran-type oceanc function which builds a search tree etc. */
void F77_FCN (oceanc) (float* wvnmlo,
                       float* wvnmhi,
                       float* mu,
                       float* mup,
                       float* dphi,
                       float* u10,
                       float* pcl,
                       float* xsal,
                       int*   firstcall,
                       float* bdref);

/* simple C wrapper for the standard oceabrdf() function, without search tree etc */
double oceabrdfc (float wvnmlo, float wvnmhi, float mu, float mup, float dphi, float u10, float pcl, float xsal, float paw);

/* simple C wrapper for the bpdf function by Mishchenko */
int bpdf_tsang (double    wind_speed,
                float     wavelength,
                double    mu_inc,
                double    phi_inc,
                double    mu_esc,
                double    phi_esc,
                double*** refl_mat);

typedef struct {
  double re;
  double im;
} ref_complex;

/* C translation of bpdf function by Mishchenko */
int bpdf_tsang_rmatr (ref_complex CN1,
                      ref_complex CN2,
                      double      sigma2,
                      double      mu_inc,
                      double      phi_inc,
                      double      mu_esc,
                      double      phi_esc,
                      int         shadowing,
                      double***   refl_mat);

double gamma_shadowing (double mu, double sigma2);

/* polarized version of the BDRF by Cox and Munk */
int cox_and_munk_polarized (float     wavelength,
                            float     mu_i,
                            float     mu_r,
                            float     phi_i,
                            float     phi_r,
                            float     u10,
                            float     pcl,
                            float     xsal,
                            float     azw,
                            int       backward,
                            double*** refl_mat);

/* bpdf function by Mishchenko for oriented ice crystals*/
int bpdf_tsang_ice (double      angdist_width,
                    ref_complex CN1,
                    ref_complex CN2,
                    double      mu_inc,
                    double      phi_inc,
                    double      mu_esc,
                    double      phi_esc,
                    double***   refl_mat);

/* stokes transmission matrix for reflection on a horizontal surface */
/*int bpdf_tsang_tmatr (ref_complex CN1,
                      ref_complex CN2,
                      double   mu_inc,
                      double   phi_inc,
                      double   mu_esc,
                      double   phi_esc,
                      double*** t_mat);*/

/* combined stokes reflection matrix for transmission and internal reflection in Kirchhoff approximation */
int bpdf_tsang_ice_internal_reflection (double      angdist_width,
                                        ref_complex CN1,
                                        ref_complex CN2,
                                        double      mu_inc,
                                        double      phi_inc,
                                        double      mu_esc,
                                        double      phi_esc,
                                        int         backward,
                                        double***   refl_mat);

int rodrigues (const double vector[3], const double axis[3], double angle, double rotated[3]);

/* combined stokes transmission matrix for transmission through a horizontally oriented crystal*/
int bpdf_tsang_ice_transmission (ref_complex cn1, ref_complex cn2, double mu_inc, double phi_inc, int backward, double*** t_mat);

#if defined(__cplusplus)
}
#endif

#endif
