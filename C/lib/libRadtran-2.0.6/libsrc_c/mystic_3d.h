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

#ifndef _MYSTIC_3D_H
#define _MYSTIC_3D_H 1

#include "mystic.h"

#if defined(__cplusplus)
extern "C" {
#endif

int step3D (photon_struct* p, atmosphere_struct* atmos, double step, int bcond, int photonpath, int spherical3D, int visualize);

int intersection3D (photon_struct* p, atmosphere_struct* atmos, double tau, double tausca, double* step);

int intersection3D_spherical (photon_struct* p, atmosphere_struct* atmos, int spherical3D_scene, int* first, double* step);

int coord_spherical3D (photon_struct* p,
                       int            Nx,
                       double*        X,
                       int            Ny,
                       double*        Y,
                       int            Nz,
                       float*         Z,
                       int            Nyg,
                       double*        Yg,
                       double         r_earth,
                       int*           ic,
                       int*           jc,
                       int*           kc,
                       int*           jcg,
                       int            warn);

int setup_spherical3D (const char* type,
                       int         Nx,
                       int         Ny,
                       int         spherical3D_scene,
                       double      spherical3D_scene_lon_min,
                       double      spherical3D_scene_lon_max,
                       double      spherical3D_scene_lat_min,
                       double      spherical3D_scene_lat_max,
                       double*     delX,
                       double*     delY,
                       double*     xmax,
                       double*     ymax,
                       double**    X,
                       double**    Y,
                       int*        Nyg,
                       double**    Yg,
                       int**       realY,
                       double**    sin2_lat,
                       double**    tan2_lat,
                       int*        jc_eq,
                       int         quiet);

void cart_to_spher (const double* x, const double* dx, double* r, double* costheta, double* sintheta, double* e_r);

int move_satellite_photon_to_toa (photon_struct* p, atmosphere_struct* atmos);

int derive_cphi (double* dx1, double* dx2, double stheta, double* cphi);

#if defined(__cplusplus)
}
#endif

#endif /* _MYSTIC_3D_H */
