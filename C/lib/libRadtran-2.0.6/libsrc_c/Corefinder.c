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

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include "Corefinder.h"
#include "Stack.h"
#include "LinArray.h"

/**
 * @brief Locate veiled core
 * 
 * @param tau 3D array with optical thickness values
 * @param nx number or grid boxes in x direction
 * @param ny number or grid boxes in y direction
 * @param nz number or grid boxes in z direction
 * @param xstart x index of starting grid box
 * @param ystart y index of starting grid box
 * @param zstart y index of starting grid box
 * @param nrep number of repetitions of the algorithm. If <1, algorithm runs until total convergence.
 * @param threshold Threshold value of distance to startpoint, where the veiled core begins
 * @param core Binary return array with 1 for every pixel belonging to the veiled core
 * @return int 
 */
int find_core (arr3d_f*      k_ext,
               float         dx,
               float         dy,
               float*        dz,
               const size_t* xstart,
               const size_t* ystart,
               const size_t* zstart,
               const size_t  nstart,
               float         threshold,
               arr3d_i*      core) {
  //Verbose
  // for (size_t i = 0; i < nstart; i++)
  // {
  //  fprintf(stdout, "start: (%zu %zu %zu)\n", xstart[i], ystart[i], zstart[i]);
  // }
  // fprintf(stdout, "threshold: %.4f\n", threshold);
  // fprintf(stdout, "dx: %.4f dy: %.4f\n", dx, dy);
  // for (size_t i = 0; i < k_ext->nz; i++)
  // {
  //     fprintf(stdout, "k_ext[%zu]=%.4f dz[%zu]=%.4f\n", i,get3d_f(k_ext,0,0,i), i,dz[i]);
  // }
  // fprintf(stdout, "\n");
  arr3d_f* distances = get_distances2 (*k_ext, dx, dy, dz, xstart, ystart, zstart, nstart, threshold);
  black_white_filter (distances, threshold, core);
  free3d_float (distances);
  return 0;
}

int black_white_filter (arr3d_f* image, float threshold, arr3d_i* binary) {
  int size  = (image->nx) * (image->ny) * (image->nz);
  int count = 0;
  while (count < size) {
    if (image->data[count] > threshold)
      binary->data[count] = 1;
    else
      binary->data[count] = 0;
    count++;
  }
  return 0;
}

int get_distance (arr3d_f* k_ext, float dx, float dy, float* dz, int xstart, int ystart, int zstart, int nrep, arr3d_f* distances) {
  arr3d_f* out_image = distances;
  // copy3d_f_metadata(distances, out_image);
  // out_image->data=distances->data;//out_image points to same data as distances
  int size = (distances->nx) * (distances->ny) * (distances->nz);
  for (int count = 0; count < size; count++)
    distances->data[count] = FLT_MAX;

  arr3d_f* out_image2 = calloc3d_float (k_ext->nx, k_ext->ny, k_ext->nz);
  arr4d_f* strel      = calloc4d_float (k_ext->nz, 3, 3, 3);
  fill_strel4d (strel, dx, dy, dz);
  arr3d_f* dummy;
  int      equal = 0;
  int      count = 0;
  while (1) {
    count++;
    set3d_f (out_image, xstart, ystart, zstart, 0.0);
    minimum_filter_enhanced (out_image, k_ext, out_image2, strel);
    dummy      = out_image;  //swap arrays (notebook entry 16)
    out_image  = out_image2;
    out_image2 = dummy;
    set3d_f (out_image, xstart, ystart, zstart, 0.0);  //start has zero distance to itself
    if (nrep < 1) {
      compare_images (out_image, out_image2, &equal);
      if (equal == 1)
        break;
    } else if (count == nrep) {
      break;
    }
    if (count % 10 == 0)
      printf ("%d\n", count);
  }
  printf ("Iterations: %d\n", count);
  set3d_f (out_image, xstart, ystart, zstart, 0.0);
  if (
    distances !=
    out_image)  //out_image always points to the final array. If this is not the original one, but the one we defined internally, we have to copy the data. (notebook entry 16)
  {
    if (distances == out_image2)  //Another confirmation
    {
      copy3d_f2 (out_image, distances);
      dummy      = out_image;  //do another swap
      out_image  = out_image2;
      out_image2 = dummy;
    } else
      return -1;
  }
  free3d_float (out_image2);
  free4d_float (strel);
  return 0;
}

/**
 * @brief Compare two images and return equal=1 if they are the same
 * 
 * @param image1 
 * @param image2 
 * @param nx 
 * @param ny 
 * @param nz 
 * @param equal 
 * @return int 
 */
int compare_images (arr3d_f* image1, arr3d_f* image2, int* equal) {
  *equal    = 1;
  int count = 0;
  int size  = (image1->nx) * (image1->ny) * (image1->nz);
  while (count < size) {
    if (image1->data[count] != image2->data[count])  //Traverse Memory linearly!
    {
      *equal = 0;
      break;
    }
    count++;
  }
  return 0;
}

int fill_strel4d (arr4d_f* strel, const float dx, const float dy, const float* dz) {
  for (int i = 0; i < strel->na; i++) {
    fill_strel3d (strel, dx, dy, dz[i], i);
  }
  return 0;
}

/**
 * @brief Create a cubic structuring element that takes into account a non-cubic atmosphere grid (notebook entry 15)
 * The values in the "strel" contain all the different diagonals in a cube (divided by 2)
 * @param strel nzx3x3x3 array
 * @param dx grid spacing in x
 * @param dy grid spacing in y
 * @param dz grid spacing in z
 * @param index_4 The index of the first dimension to be filled
 * @return int 
 */
int fill_strel3d (arr4d_f* strel, const float dx, const float dy, const float dz, const int index_4) {
  //Conversion from old array style: strel\[(\d)\]\[(\d)\]\[(\d)\]=([^;]+);
  //set4d_f(strel, index_4, $1,$2,$3,$4);
  //cube center
  set4d_f (strel, index_4, 1, 1, 1, FLT_MAX);
  //Faces centers
  set4d_f (strel, index_4, 1, 1, 0, dz / 2.0);
  set4d_f (strel, index_4, 1, 1, 2, dz / 2.0);

  set4d_f (strel, index_4, 1, 0, 1, dy / 2.0);
  set4d_f (strel, index_4, 1, 2, 1, dy / 2.0);

  set4d_f (strel, index_4, 0, 1, 1, dx / 2.0);
  set4d_f (strel, index_4, 2, 1, 1, dx / 2.0);
  //Edges
  //x=1
  set4d_f (strel, index_4, 1, 0, 2, sqrt (dy * dy + dz * dz) / 2.0);
  set4d_f (strel, index_4, 1, 2, 0, sqrt (dy * dy + dz * dz) / 2.0);
  set4d_f (strel, index_4, 1, 2, 2, sqrt (dy * dy + dz * dz) / 2.0);
  set4d_f (strel, index_4, 1, 0, 0, sqrt (dy * dy + dz * dz) / 2.0);
  //y=1
  set4d_f (strel, index_4, 0, 1, 0, sqrt (dx * dx + dz * dz) / 2.0);
  set4d_f (strel, index_4, 0, 1, 2, sqrt (dx * dx + dz * dz) / 2.0);
  set4d_f (strel, index_4, 2, 1, 0, sqrt (dx * dx + dz * dz) / 2.0);
  set4d_f (strel, index_4, 2, 1, 2, sqrt (dx * dx + dz * dz) / 2.0);
  //z=1
  set4d_f (strel, index_4, 0, 0, 1, sqrt (dx * dx + dy * dy) / 2.0);
  set4d_f (strel, index_4, 0, 2, 1, sqrt (dx * dx + dy * dy) / 2.0);
  set4d_f (strel, index_4, 2, 0, 1, sqrt (dx * dx + dy * dy) / 2.0);
  set4d_f (strel, index_4, 2, 2, 1, sqrt (dx * dx + dy * dy) / 2.0);
  //Corners
  set4d_f (strel, index_4, 0, 0, 0, sqrt (dx * dx + dy * dy + dz * dz) / 2.0);
  set4d_f (strel, index_4, 0, 0, 2, sqrt (dx * dx + dy * dy + dz * dz) / 2.0);
  set4d_f (strel, index_4, 0, 2, 0, sqrt (dx * dx + dy * dy + dz * dz) / 2.0);
  set4d_f (strel, index_4, 0, 2, 2, sqrt (dx * dx + dy * dy + dz * dz) / 2.0);
  set4d_f (strel, index_4, 2, 0, 0, sqrt (dx * dx + dy * dy + dz * dz) / 2.0);
  set4d_f (strel, index_4, 2, 0, 2, sqrt (dx * dx + dy * dy + dz * dz) / 2.0);
  set4d_f (strel, index_4, 2, 2, 0, sqrt (dx * dx + dy * dy + dz * dz) / 2.0);
  set4d_f (strel, index_4, 2, 2, 2, sqrt (dx * dx + dy * dy + dz * dz) / 2.0);
  return 0;
}

/**
 * @brief Given a guess of minimal distances from the start and the underlying extinction field, try to optimize the distances.
 * 
 * @param distances Distance field at the recent step. Set to infinity for starting.
 * @param k_ext Extinction coefficient
 * @param out_distances Resulting, optimized field
 * @param nx Number of grid boxes in x
 * @param ny Number of grid boxes in y
 * @param nz Number of grid boxes in z
 * @param weight_mask 
 * @return int 
 */
int minimum_filter_enhanced (arr3d_f* distances, arr3d_f* k_ext, arr3d_f* out_distances, arr4d_f* weight_mask) {
  float path;
  int   i, j, k, a, b, c;
  for (i = 0; i < k_ext->nx; i++) {
    for (j = 0; j < k_ext->ny; j++) {
      for (k = 0; k < k_ext->nz; k++) {
        set3d_f (out_distances, i, j, k, FLT_MAX);
        for (a = -1; a < 2; a++) {
          for (b = -1; b < 2; b++) {
            for (c = -1; c < 2; c++) {
              if (i + a > -1 && i + a < k_ext->nx && j + b > -1 && j + b < k_ext->ny && k + c > -1 && k + c < k_ext->nz) {
                path =
                  get3d_f (distances, i + a, j + b, k + c) +
                  get4d_f (weight_mask, k, a + 1, b + 1, c + 1) * (get3d_f (k_ext, i, j, k) + get3d_f (k_ext, i + a, j + b, k + c));
                if (path < get3d_f (out_distances, i, j, k))
                  set3d_f (out_distances, i, j, k, path);
              }
            } /* ends loop over c*/
          }   /* ends loop over b*/
        }     /* ends loop over a*/
      }       /* ends loop over k*/
    }         /* ends loop over j */
  }           /* ends loop over i */
  return 0;
}
/**
 * @brief like minimum_filter_enhanced(), but prints out origin of photons
 * 
 * @param distances 
 * @param k_ext 
 * @param out_distances 
 * @param nx 
 * @param ny 
 * @param nz 
 * @param weight_mask 
 * @return int 
 */
int minimum_filter_enhanced2 (arr3d_f* distances,
                              arr3d_f* k_ext,
                              arr3d_f* out_distances,
                              arr3d_i* originx,
                              arr3d_i* originy,
                              arr3d_i* originz,
                              arr4d_f* weight_mask) {
  float path;
  int   i, j, k, a, b, c;
  for (i = 0; i < k_ext->nx; i++) {
    for (j = 0; j < k_ext->ny; j++) {
      for (k = 0; k < k_ext->nz; k++) {
        set3d_f (out_distances, i, j, k, FLT_MAX);
        for (a = -1; a < 2; a++) {
          for (b = -1; b < 2; b++) {
            for (c = -1; c < 2; c++) {
              if (i + a > -1 && i + a < k_ext->nx && j + b > -1 && j + b < k_ext->ny && k + c > -1 && k + c < k_ext->nz) {
                path =
                  get3d_f (distances, i + a, j + b, k + c) +
                  get4d_f (weight_mask, k, a + 1, b + 1, c + 1) * (get3d_f (k_ext, i, j, k) + get3d_f (k_ext, i + a, j + b, k + c));
                if (path < get3d_f (out_distances, i, j, k)) {
                  set3d_f (out_distances, i, j, k, path);
                  set3d_i (originx, i, j, k, a);
                  set3d_i (originy, i, j, k, b);
                  set3d_i (originz, i, j, k, c);
                }
              }
            } /* ends loop over c*/
          }   /* ends loop over b*/
        }     /* ends loop over a*/
      }       /* ends loop over k*/
    }         /* ends loop over j */
  }           /* ends loop over i */
  return 0;
}

/**
 * @brief Add image2 to image1
 * 
 * @param image1 
 * @param image2 
 * @return int 
 */
int add_images (arr3d_f* image1, arr3d_f* image2) {
  int size  = (image1->nx) * (image1->ny) * (image1->nz);
  int count = 0;
  while (count < size) {
    image1->data[count] += image2->data[count];
    count++;
  }
  return 0;
}

//=========================================================================================================================================
void handle_point (my_stack_t*   stack,
                   arr3d_f*      best_guess,
                   arr3d_f*      best_guess_surface,
                   const arr3d_f k_ext,
                   arr3d_f*      distances,
                   const arr4d_f strel) {
  point_t p = stack_pop (stack);
  if (get3d_f (distances, p.x, p.y, p.z) != FLT_MAX)  //already reached this point
  {
    return;
  }
  set3d_f (distances, p.x, p.y, p.z, p.cost);  //Set this point and add all unset neighbours to stack
  assert (get3d_f (distances, p.x, p.y, p.z) == get3d_f (best_guess, p.x, p.y, p.z));
  // fprintf(stderr, "%d %d %d %.2f\n", p.x,p.y,p.z,p.cost);
  int   a, b, c;
  float path;
  float path_surface;
  for (a = -1; a < 2; a++) {
    for (b = -1; b < 2; b++) {
      for (c = -1; c < 2; c++) {
        if (p.x + a > -1 && p.x + a < k_ext.nx && p.y + b > -1 && p.y + b < k_ext.ny && p.z + c > -1 &&
            p.z + c < k_ext.nz)  //Check boundaries
        {
          if (a != 0 || b != 0 || c != 0)  //Avoid center
          {
            path_surface = get3d_f (distances, p.x, p.y, p.z) +
                           get4d_f (&strel, p.z, a + 1, b + 1, c + 1) *
                             get3d_f (&k_ext, p.x, p.y, p.z);  //Path to the surface of the neighbouring cube
            path = path_surface + get4d_f (&strel, p.z + c, a + 1, b + 1, c + 1) *
                                    get3d_f (&k_ext, p.x + a, p.y + b, p.z + c);  //Path to the center of the neighbouring cube
            assert (get3d_f (best_guess_surface, p.x, p.y, p.z) <= get3d_f (best_guess, p.x, p.y, p.z));
            assert (path >= p.cost);
            if (get3d_f (distances, p.x + a, p.y + b, p.z + c) == FLT_MAX)  //if neighbour is not yet set finally
            {
              if (path < get3d_f (best_guess, p.x + a, p.y + b, p.z + c))  //and if the actual guess is better
              {
                point_t p_changed = {p.x + a, p.y + b, p.z + c, path};
                set3d_f (best_guess, p_changed.x, p_changed.y, p_changed.z, p_changed.cost);  //update best guess
                stack_push (stack, p_changed);                                                //and add to stack
                // fprintf(stderr, "Stack push: %d %d %d %.4f",p_changed.x, p_changed.y, p_changed.z, p_changed.cost);
              }
              if (path_surface < get3d_f (best_guess_surface, p.x + a, p.y + b, p.z + c)) {
                set3d_f (best_guess_surface, p.x + a, p.y + b, p.z + c, path_surface);
              }
            }
          }
        }
      } /* ends loop over c*/
    }   /* ends loop over b*/
  }     /* ends loop over a*/
}

arr3d_f* get_distances2 (const arr3d_f k_ext,
                         const float   dx,
                         const float   dy,
                         const float*  dz,
                         const size_t* xstart,
                         const size_t* ystart,
                         const size_t* zstart,
                         const size_t  nstart,
                         const float   threshold) {
  arr3d_f* distances          = calloc3d_float (k_ext.nx, k_ext.ny, k_ext.nz);
  arr3d_f* best_guess         = calloc3d_float (k_ext.nx, k_ext.ny, k_ext.nz);
  arr3d_f* best_guess_surface = calloc3d_float (k_ext.nx, k_ext.ny, k_ext.nz);
  size_t   size               = (k_ext.nx) * (k_ext.ny) * (k_ext.nz);
  for (size_t i = 0; i < size; i++) {
    distances->data[i]          = FLT_MAX;
    best_guess->data[i]         = FLT_MAX;
    best_guess_surface->data[i] = FLT_MAX;
  }
  // set3d_f(distances, xstart, ystart, zstart,0.0);
  arr4d_f* strel = calloc4d_float (k_ext.nz, 3, 3, 3);
  fill_strel4d (strel, dx, dy, dz);
  my_stack_t stack = create_stack (nstart);
  for (size_t i = 0; i < nstart; i++) {
    point_t start_p = {xstart[i], ystart[i], zstart[i], 0.0};
    stack_push (&stack, start_p);
    set3d_f (best_guess, start_p.x, start_p.y, start_p.z, 0.0);
    set3d_f (best_guess_surface, start_p.x, start_p.y, start_p.z, 0.0);
  }

  // int count=0;
  while (!stack_empty (&stack)) {
    if (stack.points[0].cost > threshold) {
      break;
    }
    handle_point (&stack, best_guess, best_guess_surface, k_ext, distances, *strel);
    // if(count%1==0)
    // {
    //     fprintf(stdout, "count %d, size=%d\n", count, stack.size);
    // }
    // count++;
  }
  free3d_float (best_guess);
  free3d_float (distances);
  destroy_stack (&stack);
  free4d_float (strel);
  return best_guess_surface;
}
