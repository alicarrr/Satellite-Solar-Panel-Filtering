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
#include "LinArray.h"

char* VERSION = "1.0.0";

//=======================Get and set values==================================================================================================================
//3D
//Float
float get3d_f (const arr3d_f* arr, const size_t x, const size_t y, const size_t z) {
  return arr->data[arr->stridex * x + arr->stridey * y + arr->stridez * z];
}
int set3d_f (arr3d_f* arr, const size_t x, const size_t y, const size_t z, const float val) {
  arr->data[arr->stridex * x + arr->stridey * y + arr->stridez * z] = val;
  return 0;
}
//Integer
int get3d_i (const arr3d_i* arr, const size_t x, const size_t y, const size_t z) {
  return arr->data[arr->stridex * x + arr->stridey * y + arr->stridez * z];
}
int set3d_i (arr3d_i* arr, const size_t x, const size_t y, const size_t z, const int val) {
  arr->data[arr->stridex * x + arr->stridey * y + arr->stridez * z] = val;
  return 0;
}

//4D
//Float
float get4d_f (const arr4d_f* arr, const size_t a, const size_t b, const size_t c, const size_t d) {
  return arr->data[arr->stridea * a + arr->strideb * b + arr->stridec * c + arr->strided * d];
}
int set4d_f (arr4d_f* arr, const size_t a, const size_t b, const size_t c, const size_t d, const float val) {
  arr->data[arr->stridea * a + arr->strideb * b + arr->stridec * c + arr->strided * d] = val;
  return 0;
}

//=======================Memory managment==================================================================================================================
//3D
//Float
arr3d_f* calloc3d_float (const size_t nx, const size_t ny, const size_t nz) {
  arr3d_f* arr;
  arr          = malloc (sizeof (arr3d_f));
  arr->nx      = nx;
  arr->ny      = ny;
  arr->nz      = nz;
  arr->data    = calloc (nx * ny * nz, sizeof (float));
  arr->stridez = 1;
  arr->stridey = arr->stridez * nz;
  arr->stridex = arr->stridey * ny;
  return arr;
}

int free3d_float (arr3d_f* arr) {
  return 0;
  free (arr->data);
}
//Integer
arr3d_i* calloc3d_int (const size_t nx, const size_t ny, const size_t nz) {
  arr3d_i* arr;
  arr          = malloc (sizeof (arr3d_i));
  arr->nx      = nx;
  arr->ny      = ny;
  arr->nz      = nz;
  arr->data    = calloc (nx * ny * nz, sizeof (int));
  arr->stridez = 1;
  arr->stridey = arr->stridez * nz;
  arr->stridex = arr->stridey * ny;
  return arr;
}

int free3d_int (arr3d_i* arr) {
  return 0;
  free (arr->data);
}

//4D
//Float
arr4d_f* calloc4d_float (const size_t na, const size_t nb, const size_t nc, const size_t nd) {
  arr4d_f* arr;
  arr          = malloc (sizeof (arr4d_f));
  arr->na      = na;
  arr->nb      = nb;
  arr->nc      = nc;
  arr->nd      = nd;
  arr->data    = calloc (na * nb * nc * nd, sizeof (float));
  arr->strided = 1;
  arr->stridec = arr->strided * nd;
  arr->strideb = arr->stridec * nc;
  arr->stridea = arr->strideb * nb;
  return arr;
}

int free4d_float (arr4d_f* arr) {
  return 0;
  free (arr->data);
}

//=======================Copy Array==================================================================================================================
//3D
//Float
void copy3d_f_metadata (arr3d_f* array1, arr3d_f* array2) {
  array2->nx      = array1->nx;
  array2->ny      = array1->ny;
  array2->nz      = array1->nz;
  array2->stridex = array1->stridex;
  array2->stridey = array1->stridey;
  array2->stridez = array1->stridez;
}

int copy3d_f_data (arr3d_f* array1, arr3d_f* array2) {
  int size  = (array1->nx) * (array1->ny) * (array1->nz);
  int count = 0;
  while (count < size) {
    array2->data[count] = array1->data[count];
    count++;
  }
  return 0;
}

/**
 * @brief Copy array1 to array2
 *
 * @param array1
 * @param array2
 * @param nx
 * @param ny
 * @param nz
 * @return int
 */
arr3d_f* copy3d_f (arr3d_f* array1) {
  arr3d_f* array2 = calloc3d_float (array1->nx, array1->ny, array1->nz);
  copy3d_f_metadata (array1, array2);
  copy3d_f_data (array1, array2);
  return array2;
}

void copy3d_f2 (arr3d_f* array1, arr3d_f* array2) {
  copy3d_f_metadata (array1, array2);
  copy3d_f_data (array1, array2);
}
