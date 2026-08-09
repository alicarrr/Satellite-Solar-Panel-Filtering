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
 * Authors: Fabian Jakub
 ************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

int triBoxOverlap (const double boxcenter[3], const double boxhalfsize[3], const double A[3], const double B[3], const double C[3]);

int check_aabb_intersection (const double a_min[3], const double a_max[3], const double b_min[3], const double b_max[3]);
#ifdef __cplusplus
}
#endif
