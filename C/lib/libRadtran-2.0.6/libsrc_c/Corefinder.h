/************************************************************************
 * $Id$
 *
 * Corefinder - Identify veiled core; authored by Paul Ockenfuss
 *
 * Copyright (c) 2000-2019 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Paul Ockenfuss
 *
 * Correspondence: bernhard.mayer@lmu.de
 *
 ************************************************************************/

#ifndef COREFINDER_H
#define COREFINDER_H

#include "LinArray.h"
#include "Stack.h"

int  find_core (arr3d_f*      k_ext,
                float         dx,
                float         dy,
                float*        dz,
                const size_t* xstart,
                const size_t* ystart,
                const size_t* zstart,
                const size_t  nstart,
                float         threshold,
                arr3d_i*      core);
int  black_white_filter (arr3d_f* image, float threshold, arr3d_i* binary);
int  get_distance (arr3d_f* tau, float dx, float dy, float* dz, int xstart, int ystart, int zstart, int nrep, arr3d_f* distances);
int  compare_images (arr3d_f* image1, arr3d_f* image2, int* equal);
int  fill_strel4d (arr4d_f* strel, const float dx, const float dy, const float* dz);
int  fill_strel3d (arr4d_f* strel, const float dx, const float dy, const float dz, const int index_4);
int  minimum_filter_enhanced (arr3d_f* distances, arr3d_f* k_ext, arr3d_f* out_distances, arr4d_f* weight_mask);
int  minimum_filter_enhanced2 (arr3d_f* distances,
                               arr3d_f* k_ext,
                               arr3d_f* out_distances,
                               arr3d_i* originx,
                               arr3d_i* originy,
                               arr3d_i* originz,
                               arr4d_f* weight_mask);
int  add_images (arr3d_f* image1, arr3d_f* image2);
void handle_point (my_stack_t*   stack,
                   arr3d_f*      best_guess,
                   arr3d_f*      best_guess_surface,
                   const arr3d_f k_ext,
                   arr3d_f*      distances,
                   const arr4d_f strel);
arr3d_f* get_distances2 (const arr3d_f k_ext,
                         const float   dx,
                         const float   dy,
                         const float*  dz,
                         const size_t* xstart,
                         const size_t* ystart,
                         const size_t* zstart,
                         const size_t  nstart,
                         const float   threshold);

#endif /*COREFINDER_H*/
