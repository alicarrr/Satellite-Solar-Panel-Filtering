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

#ifndef _ALIS_H
#define _ALIS_H 1

#include "mystic.h"

#if defined(__cplusplus)
extern "C" {
#endif

int spectral_is_weight (double** totweight_spectral, photon_struct* p, photon_struct* photon, atmosphere_struct* atmos);

int concentration_is_weight (double** totweight_concentration, photon_struct* p, photon_struct* photon, atmosphere_struct* atmos);

jacobian_result_field* calloc_jacobian_result_field (int std, int Nx, int Ny, int Nz, int n_caoth, int nstokes);

void write_jacobian_result (sample_struct* sample, atmosphere_struct* atmos, jacobian_result_field* jacobian, FILE* backjacfile);

void free_jacobian_result_field (jacobian_result_field* res);

int calloc_photon_jacobian_part (photon_struct* p, sample_struct* sample, int Nx, int Ny, int Nz, int source, int n_caoth);

#if defined(__cplusplus)
}
#endif

#endif /* _ALIS_H */
