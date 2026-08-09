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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#if HAVE_NETCDF4
#include <netcdf.h>
#define MAX_NAME_LEN NC_MAX_NAME
#else
#define MAX_NAME_LEN 128
#endif

#include "allocnd.h"
#include "common_math.h"
#include "errors.h"
#include "netCDF_functions.h"

#define NCCHKERR(e)                                                                                                                \
  {                                                                                                                                \
    if (e) {                                                                                                                       \
      printf ("Netcdf Error: %s\n", nc_strerror (e));                                                                              \
    }                                                                                                                              \
    CHKERR (e);                                                                                                                    \
  }

int print_have_netcdf() {
  fprintf (stderr, " ***********************************************************************\n");
  fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
  fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
  fprintf (stderr, " ***********************************************************************\n");
  CHKERR (-1);
}

/***********************************************************************************/
/* Function: write_netCDF_1D_float                                                 */
/* Description:                                                                    */
/*  Write a 1D float (like latitude array), according attributes and dimension     */
/*  (number of array entries) to an opened netCDF result file.                     */
/*  The function checks, if this variables already                                 */
/*  exists in the netCDF header. If not, it is created.                            */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the netCDF input file                       */
/*     char *filename     filename of the result file                              */
/*     char *data_name    short name of the data set, e.g. "lat"                   */
/*     int *id_data       id number of the data (is created, of data is new)       */
/*     char *long_name    long name of the data set. e.g."latitude"                */
/*     char *units        units of the data set, e.g. "degrees east"               */
/*     float *data        the data set itself                                      */
/*     char *dim_name     name of the according dimension                          */
/*                        (in general the same as data_name), e.g. "lat"           */
/*     int id_dim         id number of the dimension                               */
/*     int n              number of array entries (dimension size)                 */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int deprecate_write_netCDF_1D_float (int    ncid,
                                     char*  filename,
                                     char*  data_name,
                                     int*   id_data,
                                     char*  long_name,
                                     char*  units,
                                     float* data,
                                     char*  dim_name,
                                     int*   id_dim,
                                     int    n) {
  DEPRECATE (2023);
#if HAVE_NETCDF4

  int    status = NC_NOERR;
  size_t iv     = 0;

  /* get variable id for data_name */

  status = nc_inq_varid (ncid, data_name, id_data);
  if (status != NC_NOERR) {
    NCCHKERR (status);

    /* variable does not exist, create it! */
    status = nc_redef (ncid); /* put in define mode */
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d putting netCDF file %s into define mode\n", status, filename);
      NCCHKERR (status);
    }

    /* create dimension */
    status = nc_def_dim (ncid, dim_name, n, id_dim);
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d defining data '%s' in netCDF file %s\n", status, data_name, filename);
      NCCHKERR (status);
    }

    status = nc_def_var (ncid, data_name, NC_FLOAT, 1, id_dim, id_data);
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d defining data %s in netCDF file %s\n", status, data_name, filename);
      NCCHKERR (status);
    }

    status = nc_enddef (ncid); /* end define mode */
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d leaving the define mode for %s \n", status, filename);
      NCCHKERR (status);
      ;
    }

    /* write data into the variable */
    for (iv = 0; iv < n; iv++) {
      /* WRITE ONE VARIABLE */
      status = nc_put_var1_float (ncid, (*id_data), &iv, &(data[iv]));
      if (status != NC_NOERR) {
        fprintf (stderr, "Error %d writing %s in netCDF file %s\n", status, data_name, filename);
        NCCHKERR (status);
      }
    }
  } else {
    /* yes, check if they are identical */
  }

#else
  CHKERR (print_have_netcdf());
#endif
  return 0;
}

/***********************************************************************************/
/* Function: write_netCDF_3D_float                                                 */
/* Description:                                                                    */
/*  Write one result (float) (time, lat, lon) into a 3D array                      */
/*  The function checks, if this data set already exits  in  the netCDF header.    */
/*  If not it will be created with according attributes, units ... and             */
/*  initialiesed to 0.                                                             */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the netCDF input file                       */
/*     char *filename     filename of the result file                              */
/*     char *data_name    short name of the data set, e.g. "lat"                   */
/*     char *long_name    long name of the data set. e.g."latitude"                */
/*     char *units        units of the data set, e.g. "degrees east"               */
/*     float data         data to write                                            */
/*     int *id_dim_array  array of dimension id numbers                            */
/*     size_t *index_3D   where to write the data                                  */
/*     int ndim           number of dimenision (here always 3)                     */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Jan 2008   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int deprecate_write_netCDF_3D_float (int     ncid,
                                     char*   filename,
                                     char*   data_name,
                                     char*   long_name,
                                     char*   units,
                                     float   data,
                                     int*    id_dim_array,
                                     size_t* index,
                                     int     ndim) {
  DEPRECATE (2023);
#if HAVE_NETCDF4

  int status  = NC_NOERR;
  int id_data = NOT_DEFINED_INTEGER;

  status = deprecate_create_float_variable_or_get_id (ncid, filename, data_name, &id_data, ndim, id_dim_array, long_name, units);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d creating float variable %s in netCDF file %s\n", status, long_name, filename);
    NCCHKERR (status);
  }

  /* WRITE ONE VARIABLE */
  status = nc_put_var1_float (ncid, id_data, index, &(data));
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d writing %s in netCDF file %s\n", status, data_name, filename);
    NCCHKERR (status);
    ;
  }

#else
  CHKERR (print_have_netcdf());
#endif
  return 0;
}

/***********************************************************************************/
/* Function: write_netCDF_4D_float                                                 */
/* Description:                                                                    */
/*  Write one column of results (mostly z-column) into a 4D array                  */
/*  (time, z, lat, lon). The function checks, if this data set already exits in    */
/*  the netCDF header. If not it will be created with according attributes,        */
/*  units ...                                                                      */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the netCDF input file                       */
/*     char *filename     filename of the result file                              */
/*     char *data_name    short name of the data set, e.g. "lat"                   */
/*     int *id_data       id number of the data (is created, of data is new)       */
/*     char *long_name    long name of the data set. e.g."latitude"                */
/*     char *units        units of the data set, e.g. "degrees east"               */
/*     float *data        the data set itself                                      */
/*     int *id_dim_array  array of dimension id numbers                            */
/*     size_t *index      specifies the position, where to write the entry         */
/*                        into the 4D file                                         */
/*     int ndim           number of dimenision (here always 4)                     */
/*     size_t *n_end      size of the dimensions                                   */
/*     int *zout_index    index of the zout levels                                 */
/*                        (nessesary, as some zout-levels, which are below the     */
/*                        ground altitude are obmitted)                            */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int deprecate_write_netCDF_4D_float (int     ncid,
                                     char*   filename,
                                     char*   data_name,
                                     char*   long_name,
                                     char*   units,
                                     float*  data,
                                     int*    id_dim_array,
                                     size_t* index,
                                     int     ndim,
                                     size_t* n_end,
                                     int*    zout_index,
                                     int*    simulation_done) {
  DEPRECATE (2023);
#if HAVE_NETCDF4

  int    status  = NC_NOERR;
  int    id_data = NOT_DEFINED_INTEGER;
  size_t lev     = 0;

  /* check if this result is valid */
  for (lev = 0; lev < n_end[POSITION_ZOUT_4D]; lev++) {
    /* replace infinite numbers with NaN (Not a number) */
    if (isinf (data[lev]))
      data[lev] = 0.0 / 0.0;
    /* check for NaN (Not a Number) */
    if (isnan (data[lev])) {
      *simulation_done = FALSE;
      fprintf (stderr, " *** Warning, unresonable result for '%s' (lev=%d)\n", long_name, (int) lev);
    }
  }

  status = deprecate_create_float_variable_or_get_id (ncid, filename, data_name, &id_data, ndim, id_dim_array, long_name, units);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d creating float variable %s in netCDF file %s\n", status, long_name, filename);
    NCCHKERR (status);
  }

  /* write data into the variable */
  for (lev = 0; lev < n_end[POSITION_ZOUT_4D]; lev++) {
    index[POSITION_ZOUT_4D] = zout_index[lev];

    /* WRITE ONE VARIABLE */
    status = nc_put_var1_float (ncid, id_data, index, &(data[lev]));
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d writing %s in netCDF file %s\n", status, data_name, filename);
      NCCHKERR (status);
    }
  }

#else
  CHKERR (print_have_netcdf());
#endif
  return 0;
}

/***********************************************************************************/
/* Function: write_netCDF_5D_float                                                 */
/* Description:                                                                    */
/*  Write one column of results (mostly z-column) into a 5D array                  */
/*  (time, z, sza, lat, lon). The function checks, if this data set already exits  */
/*   in  the netCDF header. If not it will be created with according attributes,   */
/*  units ...                                                                      */
/*  If result is not valid (infinity or Not-a-number) simulation done will be      */
/*  set to FALSE == not done. infinity is replaced by Not-a-number.                */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the netCDF input file                       */
/*     char *filename     filename of the result file                              */
/*     char *data_name    short name of the data set, e.g. "lat"                   */
/*     int *id_data       id number of the data (is created, of data is new)       */
/*     char *long_name    long name of the data set. e.g."latitude"                */
/*     char *units        units of the data set, e.g. "degrees east"               */
/*     float *data        the data set itself                                      */
/*     int *id_dim_array  array of dimension id numbers                            */
/*     size_t *index      specifies the position, where to write the entry         */
/*                        into the 4D file                                         */
/*     int ndim           number of dimenision (here always 5)                     */
/*     size_t *n_end      size of the dimensions                                   */
/*     int *zout_index    index of the zout levels                                 */
/*                        (nessesary, as some zout-levels, which are below the     */
/*                        ground altitude are obmitted)                            */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int deprecate_write_netCDF_5D_float (int     ncid,
                                     char*   filename,
                                     char*   data_name,
                                     char*   long_name,
                                     char*   units,
                                     float** data,
                                     int*    id_dim_array,
                                     size_t* index,
                                     int     ndim,
                                     size_t* n_end,
                                     int*    zout_index,
                                     int*    simulation_done) {
  DEPRECATE (2023);
#if HAVE_NETCDF4

  int    status  = NC_NOERR;
  int    id_data = NOT_DEFINED_INTEGER;
  size_t lev     = 0;
  size_t iv      = 0;

  /* check if this result is valid */
  for (lev = 0; lev < n_end[POSITION_ZOUT_5D]; lev++)
    for (iv = 0; iv < n_end[POSITION_LAMBDA_5D]; iv++) {
      /* replace infinite numbers with NaN (Not a number) */
      if (isinf (data[lev][iv]))
        data[lev][iv] = 0.0 / 0.0;
      /* check for NaN (Not a Number) */
      if (isnan (data[lev][iv])) {
        *simulation_done = FALSE;
        fprintf (stderr, " *** Warning, unresonable result for '%s' (iv=%d, lev=%d)\n", long_name, (int) iv, (int) lev);
      }
    }

  status = deprecate_create_float_variable_or_get_id (ncid, filename, data_name, &id_data, ndim, id_dim_array, long_name, units);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d creating float variable %s in netCDF file %s\n", status, long_name, filename);
    NCCHKERR (status);
  }

  for (lev = 0; lev < n_end[POSITION_ZOUT_5D]; lev++) {
    index[POSITION_ZOUT_5D] = zout_index[lev];

    for (iv = 0; iv < n_end[POSITION_LAMBDA_5D]; iv++) {
      index[POSITION_LAMBDA_5D] = iv;

      /* WRITE ONE VARIABLE */
      status = nc_put_var1_float (ncid, id_data, index, &(data[lev][iv]));
      if (status != NC_NOERR) {
        fprintf (stderr, "Error %d writing %s in netCDF file %s\n", status, data_name, filename);
        NCCHKERR (status);
      }
    }
  }

#else
  CHKERR (print_have_netcdf());
#endif
  return 0;
}

/***********************************************************************************/
/* Function: write_netCDF_6D_float                                                 */
/* Description:                                                                    */
/*  Write one column of results (mostly z-column) into a 6D array                  */
/*  (time, z, sza, lat, lon). The function checks, if this data set already exits  */
/*   in  the netCDF header. If not it will be created with according attributes,   */
/*  units ...                                                                      */
/*  If result is not valid (infinity or Not-a-number) simulation done will be      */
/*  set to FALSE == not done. infinity is replaced by Not-a-number.                */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the netCDF input file                       */
/*     char *filename     filename of the result file                              */
/*     char *data_name    short name of the data set, e.g. "lat"                   */
/*     int *id_data       id number of the data (is created, of data is new)       */
/*     char *long_name    long name of the data set. e.g."latitude"                */
/*     char *units        units of the data set, e.g. "degrees east"               */
/*     float *data        the data set itself                                      */
/*     int *id_dim_array  array of dimension id numbers                            */
/*     size_t *index      specifies the position, where to write the entry         */
/*                        into the 4D file                                         */
/*     int ndim           number of dimenision (here always 6)                     */
/*     size_t *n_end      size of the dimensions                                   */
/*     int *zout_index    index of the zout levels                                 */
/*                        (nessesary, as some zout-levels, which are below the     */
/*                        ground altitude are obmitted)                            */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int deprecate_write_netCDF_6D_float (int      ncid,
                                     char*    filename,
                                     char*    data_name,
                                     char*    long_name,
                                     char*    units,
                                     float*** data,
                                     int*     id_dim_array,
                                     size_t*  index,
                                     int      ndim,
                                     size_t*  n_end,
                                     int*     zout_index,
                                     int*     simulation_done) {
  DEPRECATE (2023);
#if HAVE_NETCDF4

  int    status  = NC_NOERR;
  int    id_data = NOT_DEFINED_INTEGER;
  size_t lev     = 0;
  size_t iv      = 0;
  size_t iu      = 0;

  /* check if this result is valid */
  for (lev = 0; lev < n_end[POSITION_ZOUT_6D]; lev++)
    for (iu = 0; iu < n_end[POSITION_UMU_6D]; iu++)
      for (iv = 0; iv < n_end[POSITION_LAMBDA_6D]; iv++) {
        /* replace infinite numbers with NaN (Not a number) */
        if (isinf (data[lev][iu][iv]))
          data[lev][iu][iv] = 0.0 / 0.0;
        /* check for NaN (Not a Number) */
        if (isnan (data[lev][iu][iv])) {
          *simulation_done = FALSE;
          fprintf (stderr,
                   " *** Warning, unresonable result for '%s' (iv=%d, lev=%d, iu=%d)\n",
                   long_name,
                   (int) iv,
                   (int) lev,
                   (int) iu);
        }
      }

  status = deprecate_create_float_variable_or_get_id (ncid, filename, data_name, &id_data, ndim, id_dim_array, long_name, units);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d creating float variable %s in netCDF file %s\n", status, long_name, filename);
    NCCHKERR (status);
  }

  for (lev = 0; lev < n_end[POSITION_ZOUT_6D]; lev++) {
    index[POSITION_ZOUT_6D] = zout_index[lev];

    for (iv = 0; iv < n_end[POSITION_LAMBDA_6D]; iv++) {
      index[POSITION_LAMBDA_6D] = iv;

      for (iu = 0; iu < n_end[POSITION_UMU_6D]; iu++) {
        index[POSITION_UMU_6D] = iu;

        /* WRITE ONE VARIABLE */
        status = nc_put_var1_float (ncid, id_data, index, &(data[lev][iu][iv]));
        if (status != NC_NOERR) {
          fprintf (stderr, "Error %d writing %s in netCDF file %s\n", status, data_name, filename);
          NCCHKERR (status);
        }
      }
    }
  }

#else
  CHKERR (print_have_netcdf());
#endif
  return 0;
}

/***********************************************************************************/
/* Function: write_netCDF_7D_float                                                 */
/* Description:                                                                    */
/*  Write one column of results (mostly z-column) into a 7D array                  */
/*  (time, z, sza, lat, lon). The function checks, if this data set already exits  */
/*   in  the netCDF header. If not it will be created with according attributes,   */
/*  units ...                                                                      */
/*  If result is not valid (infinity or Not-a-number) simulation done will be      */
/*  set to FALSE == not done. infinity is replaced by Not-a-number.                */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the netCDF input file                       */
/*     char *filename     filename of the result file                              */
/*     char *data_name    short name of the data set, e.g. "lat"                   */
/*     int *id_data       id number of the data (is created, of data is new)       */
/*     char *long_name    long name of the data set. e.g."latitude"                */
/*     char *units        units of the data set, e.g. "degrees east"               */
/*     float *data        the data set itself                                      */
/*     int *id_dim_array  array of dimension id numbers                            */
/*     size_t *index      specifies the position, where to write the entry         */
/*                        into the 4D file                                         */
/*     int ndim           number of dimenision (here always 7)                     */
/*     size_t *n_end      size of the dimensions                                   */
/*     int *zout_index    index of the zout levels                                 */
/*                        (nessesary, as some zout-levels, which are below the     */
/*                        ground altitude are obmitted)                            */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int deprecate_write_netCDF_7D_float (int       ncid,
                                     char*     filename,
                                     char*     data_name,
                                     char*     long_name,
                                     char*     units,
                                     float**** data,
                                     int*      id_dim_array,
                                     size_t*   index,
                                     int       ndim,
                                     size_t*   n_end,
                                     int*      zout_index,
                                     int*      simulation_done) {
  DEPRECATE (2023);
#if HAVE_NETCDF4

  int    status  = NC_NOERR;
  int    id_data = NOT_DEFINED_INTEGER;
  size_t lev = 0, iv = 0, iu = 0, ip = 0;

  /* check if this result is valid */
  for (lev = 0; lev < n_end[POSITION_ZOUT_6D]; lev++)
    for (ip = 0; ip < n_end[POSITION_PHI_7D]; ip++)
      for (iu = 0; iu < n_end[POSITION_UMU_6D]; iu++)
        for (iv = 0; iv < n_end[POSITION_LAMBDA_6D]; iv++) {
          /* replace infinite numbers with NaN (Not a number) */
          if (isinf (data[lev][ip][iu][iv]))
            data[lev][ip][iu][iv] = 0.0 / 0.0;
          /* check for NaN (Not a Number) */
          if (isnan (data[lev][ip][iu][iv])) {
            *simulation_done = FALSE;
            fprintf (stderr,
                     " *** Warning, unresonable result for '%s' (iv=%d, lev=%d, iu=%d, ip=%d)\n",
                     long_name,
                     (int) iv,
                     (int) lev,
                     (int) iu,
                     (int) ip);
          }
        }

  status = deprecate_create_float_variable_or_get_id (ncid, filename, data_name, &id_data, ndim, id_dim_array, long_name, units);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d creating float variable %s in netCDF file %s\n", status, long_name, filename);
    NCCHKERR (status);
  }

  for (lev = 0; lev < n_end[POSITION_ZOUT_7D]; lev++) {
    index[POSITION_ZOUT_7D] = zout_index[lev];

    for (iv = 0; iv < n_end[POSITION_LAMBDA_7D]; iv++) {
      index[POSITION_LAMBDA_7D] = iv;

      for (iu = 0; iu < n_end[POSITION_UMU_7D]; iu++) {
        index[POSITION_UMU_7D] = iu;

        for (ip = 0; ip < n_end[POSITION_PHI_7D]; ip++) {
          index[POSITION_PHI_7D] = ip;

          /* WRITE ONE VARIABLE */
          status = nc_put_var1_float (ncid, id_data, index, &(data[lev][ip][iu][iv]));
          if (status != NC_NOERR) {
            fprintf (stderr, "Error %d writing %s in netCDF file %s\n", status, data_name, filename);
            NCCHKERR (status);
          }
        }
      }
    }
  }
#else
  CHKERR (print_have_netcdf());
#endif
  return 0;
}

/***********************************************************************************/
/* Function: create_float_variable_or_get_id                                       */
/* Description:                                                                    */
/*  The function checks, if this data set already exits                            */
/*  in  the netCDF header. If not it will be created with according attributes,    */
/*  units ...                                                                      */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the netCDF input file                       */
/*     char *filename     filename of the result file                              */
/*     char *data_name    short name of the data set, e.g. "lat"                   */
/*     int *id_data       id number of the data (is created, of data is new)       */
/*     char *long_name    long name of the data set. e.g."latitude"                */
/*     char *units        units of the data set, e.g. "degrees east"               */
/*     int *id_dim_array  array of dimension id numbers                            */
/*     int ndim           number of dimenision                                     */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int deprecate_create_float_variable_or_get_id (int   ncid,
                                               char* filename,
                                               char* data_name,
                                               int*  id_data,
                                               int   ndim,
                                               int*  id_dim_array,
                                               char* long_name,
                                               char* units) {
  DEPRECATE (2023);
#if HAVE_NETCDF4

  int   status        = NC_NOERR;
  float missing_value = -1.E10;

  /* get variable id for data_name */
  status = nc_inq_varid (ncid, data_name, id_data);

  if (status != NC_NOERR) {

    /* variable does not exist, create it! */
    status = nc_redef (ncid); /* put in define mode */
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d putting netCDF file %s into define mode\n", status, filename);
      NCCHKERR (status);
    }

    /* create variable */
    status = nc_def_var (ncid, data_name, NC_FLOAT, ndim, id_dim_array, id_data);
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d defining data '%s' in netCDF file %s\n", status, data_name, filename);
      NCCHKERR (status);
    }

    status = nc_put_att_text (ncid, *id_data, "long_name", strlen (long_name), long_name);
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d put attribute long_name to '%s' in netCDF file %s\n", status, data_name, filename);
      NCCHKERR (status);
    }

    status = nc_put_att_text (ncid, *id_data, "units", strlen (units), units);
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d put attribute units to '%s' in netCDF file %s\n", status, data_name, filename);
      NCCHKERR (status);
    }

    status = nc_put_att_float (ncid, *id_data, "missing_value", NC_FLOAT, 1, &(missing_value));
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d defining missing argument in netCDF file %s\n", status, filename);
      NCCHKERR (status);
    }

    status = nc_put_att_float (ncid, *id_data, "_FillValue", NC_FLOAT, 1, &(missing_value));
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d defining missing argument in netCDF file %s\n", status, filename);
      NCCHKERR (status);
    }

    status = nc_enddef (ncid); /* end define mode */
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d leaving the define mode for %s \n", status, filename);
      NCCHKERR (status);
      ;
    }
  }

  return NC_NOERR;

#else
  CHKERR (print_have_netcdf());
#endif
}

/***********************************************************************************/
/* Function: write_netCDF_3D_byte                                                  */
/* Description:                                                                    */
/*  Write one number of results (TRUE or FLASE) into a 3D array (simulation done)  */
/*  The function checks, if this data set already exits  in  the netCDF header.    */
/*  If not it will be created with according attributes, units ... and             */
/*  initialiesed to 0.                                                             */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the netCDF input file                       */
/*     char *filename     filename of the result file                              */
/*     char *data_name    short name of the data set, e.g. "lat"                   */
/*     char *long_name    long name of the data set. e.g."latitude"                */
/*     char *units        units of the data set, e.g. "degrees east"               */
/*     short data         data to write                                            */
/*     int *id_dim_array  array of dimension id numbers                            */
/*     size_t *index_3D   where to write the data                                  */
/*     int ndim           number of dimenision (here always 3)                     */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int deprecate_write_netCDF_3D_byte (int     ncid,
                                    char*   filename,
                                    char*   data_name,
                                    char*   long_name,
                                    char*   units,
                                    short   data,
                                    int*    id_dim_array,
                                    size_t* index_3D,
                                    int     ndim) {
  DEPRECATE (2023);
#if HAVE_NETCDF4

  int status  = NC_NOERR;
  int id_data = NOT_DEFINED_INTEGER;

  short        missing_value = 0;
  size_t       ntime = 0, nlat = 0, nlon = 0;
  int          n_all    = 0;
  int          t        = -999;
  signed char* tmp_char = NULL;

  /* get variable id for data_name */
  status = nc_inq_varid (ncid, data_name, &id_data);
  NCCHKERR (status);

  if (status != NC_NOERR) {

    /* variable does not exist, create it! */
    status = nc_redef (ncid); /* put in define mode */
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d putting netCDF file %s into define mode\n", status, filename);
      NCCHKERR (status);
    }

    /* create variable */
    status = nc_def_var (ncid, data_name, NC_BYTE, ndim, id_dim_array, &id_data);
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d defining data '%s' in netCDF file %s\n", status, data_name, filename);
      NCCHKERR (status);
    }

    status = nc_put_att_text (ncid, id_data, "long_name", strlen (long_name), long_name);
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d put attribute long_name to '%s' in netCDF file %s\n", status, data_name, filename);
      NCCHKERR (status);
    }

    status = nc_put_att_text (ncid, id_data, "units", strlen (units), units);
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d put attribute units to '%s' in netCDF file %s\n", status, data_name, filename);
      NCCHKERR (status);
    }

    status = nc_put_att_short (ncid, id_data, "missing_value", NC_BYTE, 1, &(missing_value));
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d defining missing argument in netCDF file %s\n", status, filename);
      NCCHKERR (status);
    }

    status = nc_enddef (ncid); /* end define mode */
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d leaving the define mode for %s) \n", status, filename);
      NCCHKERR (status);
    }

    /* get dimension length for "ntime" */
    status = nc_inq_dimlen (ncid, id_dim_array[POSITION_TIME_3D], &ntime);
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d reading ntime from %s\n", status, filename);
      NCCHKERR (status);
    }

    /* get dimension length for "nlat" */
    status = nc_inq_dimlen (ncid, id_dim_array[POSITION_LAT_3D], &nlat);
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d reading nlat from %s\n", status, filename);
      NCCHKERR (status);
    }

    /* get dimension length for "nlon" */
    status = nc_inq_dimlen (ncid, id_dim_array[POSITION_LON_3D], &nlon);
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d reading nlon from %s\n", status, filename);
      NCCHKERR (status);
    }

    n_all = ntime * nlat * nlon;

    tmp_char = malloc (n_all * sizeof (signed char));
    CHKPOINTER (tmp_char);

    for (t = 0; t < n_all; t++)
      tmp_char[t] = 0;

    /* WRITE ONE VARIABLE */
    /* does not work on the cluster with unsigned char, that why I use signed chars here !?! */
    status = nc_put_var_schar (ncid, id_data, tmp_char);
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d writing %s in netCDF file %s\n", status, data_name, filename);
      NCCHKERR (status);
    }

    free (tmp_char);
  }

  /* WRITE ONE VARIABLE */
  status = nc_put_var1_short (ncid, id_data, index_3D, &(data));
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d writing %s to netCDF file %s\n", status, data_name, filename);
    NCCHKERR (status);
  }

#else
  CHKERR (print_have_netcdf());
#endif
  return 0;
}

/****************************************************************/
/* alloc_and_read_netCDF_1D_double                              */
/* allocate and                                                 */
/* read an 1D array direct from netCDF file                     */
/* ncid input   input   id_number of the file                   */
/* dimension    input   name of dimension number in netCDF file */
/* dim          output  number of elements                      */
/* variable     input   name of data in netCDF file             */
/* data         output  pointer to data (double)                */
/* January 2007 Ulrich Hamann                                   */
/****************************************************************/

int alloc_and_read_netCDF_1D_double (int ncid, char* dimension, size_t* dim, char* variable, double** data)

{
  int status = 0;

#if HAVE_NETCDF4
  int id_dim = 0;
  int id_var = 0;

  /* get dimension id for "dimension" */
  status = nc_inq_dimid (ncid, dimension, &id_dim);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error '%s' locating '%s' from netCDF file\n", nc_strerror (status), dimension);
    CHKERR (status);
  }

  /* get dimension length for "dimension" */
  status = nc_inq_dimlen (ncid, id_dim, dim);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error '%s' while reading '%s' from netCDF file\n", nc_strerror (status), dimension);
    CHKERR (status);
  }

  /* fprintf (stderr, "reading %s, dim=%u \n", variable, *dim); */

  /* allocate data */
  *data = calloc (*(dim), sizeof (double));
  CHKPOINTER (*data);

  /* get id number for variable */
  status = nc_inq_varid (ncid, variable, &id_var);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error '%s' while getting id for '%s'\n", nc_strerror (status), variable);
    CHKERR (status);
  }

  /* read variable (double) */
  status = nc_get_var_double (ncid, id_var, *(data));
  if (status != NC_NOERR) {
    fprintf (stderr, "Error '%s' while reading '%s' from netCDF file\n", nc_strerror (status), variable);
    CHKERR (status);
  }

#endif

  return status;
}

/***********************************************************************************/
/* Function: read_3d_float                                                         */
/* Description:                                                                    */
/*  Reads a 3D data set from an netCDF file and stores it in an 3D float array,    */
/*  which will be automatically be allocated                                       */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the read netCDF file                        */
/*     char *file         filename of the read netCDF file                         */
/*     char *var_name     short name of the data set, e.g. "lat"                   */
/*     long n1            number of entries in the 1st dimension                   */
/*     long n2            number of entries in the 2nd dimension                   */
/*     long n3            number of entries in the 3rd dimension                   */
/* Return value:                                                                   */
/*     float ****float3d  read data                                                */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Jun 2008   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int deprecate_read_3d_float (int       ncid,
                             char*     file,
                             char*     var_name,
                             float**** float3d,
                             long      n1,
                             long      n2,
                             long      n3,
                             char***   done,
                             int       verbose) {
  DEPRECATE (2023);
#if HAVE_NETCDF4

  int    status      = 0;
  int    id_float3d  = 0;
  float* tmp_float3d = NULL;

  if (verbose)
    fprintf (stderr, " ... read_float_3d: var = %-10s, data_type = %-8s, ndim = %3d \n", var_name, "float", 3);

  /* allocate space for "float3d" */
  *float3d = malloc (n1 * sizeof (float**));
  CHKPOINTER (*float3d);

  for (size_t t = 0; t < n1; t++) {
    (*float3d)[t] = malloc (n2 * sizeof (float*));
    CHKPOINTER ((*float3d)[t]);

    for (size_t j = 0; j < n2; j++) {
      (*float3d)[t][j] = malloc (n3 * sizeof (float));
      CHKPOINTER ((*float3d)[t][j]);
    }
  }

  /* get id for float3d */
  status = nc_inq_varid (ncid, var_name, &id_float3d);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error, get id for '%s' field in the netCDF file %s\n", var_name, file);
    NCCHKERR (status);
  }

  /* allocate space for "tmp_float3d" */
  tmp_float3d = malloc (n1 * n2 * n3 * sizeof (float));
  CHKPOINTER (tmp_float3d);

  /* read values from netCDF variable */
  status = nc_get_var_float (ncid, id_float3d, tmp_float3d);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d reading '%s' from %s\n", status, var_name, file);
    free (tmp_float3d);
    NCCHKERR (status);
  }

  for (size_t t = 0; t < n1; t++)
    for (size_t j = 0; j < n2; j++)
      for (size_t i = 0; i < n3; i++) {
        (*float3d)[t][j][i] = tmp_float3d[t * n2 * n3 + j * n3 + i];
        /* fprintf (stderr, "t=%4d, z=%4d, y=%4d, x=%4d, float3d=%8.2f \n",t,c2,j,i,(*float3d)[t][c2][j][i]); */
        /* fflush(stderr); */
      }

  free (tmp_float3d);

#else
  CHKERR (print_have_netcdf());
#endif
  return 0;
}

/***********************************************************************************/
/* Function: read_4d_float                                                         */
/* Description:                                                                    */
/*  Reads a 4D data set from an netCDF file and stores it in an 4D float array     */
/*  will will be accordingly allocated                                             */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the read netCDF file                        */
/*     char *file         filename of the read netCDF file                         */
/*     char *var_name     short name of the data set, e.g. "lat"                   */
/*     long n1            number of entries in the 1st dimension                   */
/*     long n2            number of entries in the 2nd dimension                   */
/*     long n3            number of entries in the 3rd dimension                   */
/*     long n4            number of entries in the 4th dimension                   */
/* Return value:                                                                   */
/*     float *****float4d read data                                                */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int deprecate_read_4d_float (int        ncid,
                             char*      file,
                             char*      var_name,
                             float***** float4d,
                             long       n1,
                             long       n2,
                             long       n3,
                             long       n4,
                             char***    done,
                             int        verbose) {
  DEPRECATE (2023);
#if HAVE_NETCDF4

  int    status = 0;
  int    t = 0, c2 = 0, j = 0, i = 0;
  int    id_float4d  = 0;
  float* tmp_float4d = NULL;

  if (verbose)
    fprintf (stderr, " ... read_float_4d: var = %-10s, data_type = %-8s, ndim = %3d \n", var_name, "float", 4);

  /* allocate space for "float4d" */
  *float4d = malloc (n1 * sizeof (float***));
  CHKPOINTER (*float4d);

  for (t = 0; t < n1; t++) {
    (*float4d)[t] = malloc (n2 * sizeof (float**));
    CHKPOINTER ((*float4d)[t]);

    for (c2 = 0; c2 < n2; c2++) {
    }

    for (j = 0; j < n3; j++) {
      (*float4d)[t][c2][j] = malloc (n4 * sizeof (float));
      CHKPOINTER ((*float4d)[t][c2][j]);
    }
  }

  /* allocate space for "tmp_float4d" */
  tmp_float4d = malloc (n1 * n2 * n3 * n4 * sizeof (float));
  CHKPOINTER (tmp_float4d);

  /* get id for float4d */
  status = nc_inq_varid (ncid, var_name, &id_float4d);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error, get id for '%s' field in the netCDF result file %s \n", var_name, file);
    NCCHKERR (-1);
  }

  /* read values from netCDF variable */
  status = nc_get_var_float (ncid, id_float4d, tmp_float4d);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d reading '%s' from %s\n", status, var_name, file);
    NCCHKERR (status);
    ;
  }

  for (t = 0; t < n1; t++)
    for (c2 = 0; c2 < n2; c2++)
      for (j = 0; j < n3; j++)
        for (i = 0; i < n4; i++) {
          (*float4d)[t][c2][j][i] = tmp_float4d[t * n2 * n3 * n4 + c2 * n3 * n4 + j * n4 + i];
          /* fprintf (stderr, "t=%4d, z=%4d, y=%4d, x=%4d, float4d=%8.2f \n",t,c2,j,i,(*float4d)[t][c2][j][i]); */
          /* fflush(stderr); */
        }

  free (tmp_float4d);
#else
  CHKERR (print_have_netcdf());
#endif
  return 0;
}

/***********************************************************************************/
/* Function: read_5d_float                                                         */
/* Description:                                                                    */
/*  Reads a 5D data set from an netCDF file and stores it in an 5D float array     */
/*  will will be accordingly allocated                                             */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the read netCDF file                        */
/*     char *file         filename of the read netCDF file                         */
/*     char *var_name     short name of the data set, e.g. "lat"                   */
/*     long n1             number of entries in the 1st dimension                  */
/*     long n2             number of entries in the 2nd dimension                  */
/*     long n3             number of entries in the 3rd dimension                  */
/*     long n4             number of entries in the 4th dimension                  */
/*     long n5             number of entries in the 5th dimension                  */
/* Return value:                                                                   */
/*     float ******float5d read data                                               */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/*                                                                                 */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int deprecate_read_5d_float (int         ncid,
                             char*       file,
                             char*       var_name,
                             float****** float5d,
                             long        n1,
                             long        n2,
                             long        n3,
                             long        n4,
                             long        n5,
                             char***     done,
                             int         verbose) {
  DEPRECATE (2023);
#if HAVE_NETCDF4

  int    status = 0;
  int    t = 0, c2 = 0, c3 = 0, j = 0, i = 0;
  int    id_float5d  = 0;
  float* tmp_float5d = NULL;

  if (verbose)
    fprintf (stderr, " ... read_float_5d: var = %-10s, data_type = %-8s, ndim = %3d \n", var_name, "float", 5);

  /* allocate space for "float5d" */
  *float5d = malloc (n1 * sizeof (float****));
  CHKPOINTER (*float5d);

  for (t = 0; t < n1; t++) {
    (*float5d)[t] = malloc (n2 * sizeof (float***));
    CHKPOINTER ((*float5d)[t]);

    for (c2 = 0; c2 < n2; c2++) {
      (*float5d)[t][c2] = malloc (n3 * sizeof (float**));
      CHKPOINTER ((*float5d)[t][c2]);

      for (c3 = 0; c3 < n3; c3++) {
        (*float5d)[t][c2][c3] = malloc (n4 * sizeof (float*));
        CHKPOINTER ((*float5d)[t][c2][c3]);

        for (j = 0; j < n4; j++) {
          (*float5d)[t][c2][c3][j] = malloc (n5 * sizeof (float));
          CHKPOINTER ((*float5d)[t][c2][c3][j]);
        }
      }
    }
  }

  /* get id for float5d */
  status = nc_inq_varid (ncid, var_name, &id_float5d);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error, get id for '%s' field in the netCDF result file %s \n", var_name, file);
    NCCHKERR (-1);
  }

  /* allocate space for "tmp_float5d" */
  tmp_float5d = malloc (n1 * n2 * n3 * n4 * n5 * sizeof (float));
  CHKPOINTER (tmp_float5d);

  /* read values from netCDF variable */
  status = nc_get_var_float (ncid, id_float5d, tmp_float5d);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d reading '%s' from %s\n", status, var_name, file);
    free (tmp_float5d);
    NCCHKERR (status);
    ;
  }

  for (t = 0; t < n1; t++)
    for (c2 = 0; c2 < n2; c2++)
      for (c3 = 0; c3 < n3; c3++)
        for (j = 0; j < n4; j++)
          for (i = 0; i < n5; i++) {
            (*float5d)[t][c2][c3][j][i] = tmp_float5d[t * n2 * n3 * n4 * n5 + c2 * n3 * n4 * n5 + c3 * n4 * n5 + j * n5 + i];
            /* fprintf (stderr, "t=%4d, lambda=%4d, z=%4d, y=%4d, x=%4d, float5d=%8.2f \n",t,c2,c3,j,i,(*float5d)[t][c2][c3][j][i]); */
            /* fflush(stderr); */
          }

  free (tmp_float5d);

#else
  CHKERR (print_have_netcdf());
#endif
  return 0;
}

/***********************************************************************************/
/* Function: read_6d_float                                                         */
/* Description:                                                                    */
/*  Reads a 6D data set from an netCDF file and stores it in an 6D float array     */
/*  will will be accordingly allocated                                             */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the read netCDF file                        */
/*     char *file         filename of the read netCDF file                         */
/*     char *var_name     short name of the data set, e.g. "lat"                   */
/*     long n1              number of entries in the 1st dimension                 */
/*     long n2              number of entries in the 2nd dimension                 */
/*     long n3              number of entries in the 3rd dimension                 */
/*     long n4              number of entries in the 4th dimension                 */
/*     long n5              number of entries in the 5th dimension                 */
/*     long n6              number of entries in the 6th dimension                 */
/* Return value:                                                                   */
/*     float *******float6d read data                                              */
/*     int status         == 0, if everthing is OK                                 */
/*                                                                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int deprecate_read_6d_float (int          ncid,
                             char*        file,
                             char*        var_name,
                             float******* float6d,
                             long         n1,
                             long         n2,
                             long         n3,
                             long         n4,
                             long         n5,
                             long         n6,
                             char***      done,
                             int          verbose) {
  DEPRECATE (2023);
#if HAVE_NETCDF4

  int    status = 0;
  int    t = 0, c2 = 0, c3 = 0, c4 = 0, j = 0, i = 0;
  int    id_float6d  = 0;
  float* tmp_float6d = NULL;

  if (verbose)
    fprintf (stderr, " ... read_float_6d: var = %-10s, data_type = %-8s, ndim = %3d \n", var_name, "float", 6);

  /* allocate space for "float6d" */
  *float6d = malloc (n1 * sizeof (float*****));
  CHKPOINTER (*float6d);

  for (t = 0; t < n1; t++) {
    (*float6d)[t] = malloc (n2 * sizeof (float****));
    CHKPOINTER ((*float6d)[t]);

    for (c2 = 0; c2 < n2; c2++) {
      (*float6d)[t][c2] = malloc (n3 * sizeof (float***));
      CHKPOINTER ((*float6d)[t][c2]);

      for (c3 = 0; c3 < n3; c3++) {
        (*float6d)[t][c2][c3] = malloc (n4 * sizeof (float**));
        CHKPOINTER ((*float6d)[t][c2][c3]);

        for (c4 = 0; c4 < n4; c4++) {
          (*float6d)[t][c2][c3][c4] = malloc (n5 * sizeof (float*));
          CHKPOINTER ((*float6d)[t][c2][c3][c4]);

          for (j = 0; j < n5; j++) {
            (*float6d)[t][c2][c3][c4][j] = malloc (n6 * sizeof (float));
            CHKPOINTER ((*float6d)[t][c2][c3][c4][j]);
          }
        }
      }
    }
  }

  /* get id for float6d */
  status = nc_inq_varid (ncid, var_name, &id_float6d);
  NCCHKERR (status - NC_NOERR);

  /* allocate space for "tmp_float6d" */
  tmp_float6d = malloc (n1 * n2 * n3 * n4 * n5 * n6 * sizeof (float));
  CHKPOINTER (tmp_float6d);

  /* read values from netCDF variable */
  status = nc_get_var_float (ncid, id_float6d, tmp_float6d);
  if (status - NC_NOERR) {
    free (tmp_float6d);
    NCCHKERR (status - NC_NOERR);
  }

  for (t = 0; t < n1; t++)
    for (c2 = 0; c2 < n2; c2++)
      for (c3 = 0; c3 < n3; c3++)
        for (c4 = 0; c4 < n4; c4++)
          for (j = 0; j < n5; j++)
            for (i = 0; i < n6; i++) {
              /* if (done[t][j][i]) { */
              (*float6d)[t][c2][c3][c4][j][i] =
                tmp_float6d[t * n2 * n3 * n4 * n5 * n6 + c2 * n3 * n4 * n5 * n6 + c3 * n4 * n5 * n6 + c4 * n5 * n6 + j * n6 + i];
              /* fprintf (stderr, "t=%4d, lambda=%4d, z=%4d, umu=%4d, y=%4d, x=%4d, float6d=%8.2f \n",t,c2,c3,c4,j,i,(*float6d)[t][c2][c3][c4][j][i]); */
              /* fflush(stderr); */
              /* } */
            }

  free (tmp_float6d);

#else
  CHKERR (print_have_netcdf());
#endif
  return 0;
}

/***********************************************************************************/
/* Function: read_7d_float                                                         */
/* Description:                                                                    */
/*  Reads a 7D data set from an netCDF file and stores it in an 7D float array     */
/*  will will be accordingly allocated                                             */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the read netCDF file                        */
/*     char *file         filename of the read netCDF file                         */
/*     char *var_name     short name of the data set, e.g. "lat"                   */
/*     long n1              number of entries in the 1st dimension                 */
/*     long n2              number of entries in the 2nd dimension                 */
/*     long n3              number of entries in the 3rd dimension                 */
/*     long n4              number of entries in the 4th dimension                 */
/*     long n5              number of entries in the 5th dimension                 */
/*     long n6              number of entries in the 6th dimension                 */
/*     long n7              number of entries in the 7th dimension                 */
/* Return value:                                                                   */
/*     float *******float6d read data                                              */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/*                                                                                 */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int deprecate_read_7d_float (int           ncid,
                             char*         file,
                             char*         var_name,
                             float******** float7d,
                             long          n1,
                             long          n2,
                             long          n3,
                             long          n4,
                             long          n5,
                             long          n6,
                             long          n7,
                             char***       done,
                             int           verbose) {
  DEPRECATE (2023);
#if HAVE_NETCDF4

  int    status = 0;
  int    t = 0, c2 = 0, c3 = 0, c4 = 0, c5 = 0, j = 0, i = 0;
  int    id_float7d  = 0;
  float* tmp_float7d = NULL;

  if (verbose)
    fprintf (stderr, " ... read_float_7d: var = %-10s, data_type = %-8s, ndim = %3d \n", var_name, "float", 7);

  /* allocate space for "float7d" */
  *float7d = malloc (n1 * sizeof (float******));
  CHKPOINTER (*float7d);

  for (t = 0; t < n1; t++) {
    (*float7d)[t] = malloc (n2 * sizeof (float*****));
    CHKPOINTER ((*float7d)[t]);

    for (c2 = 0; c2 < n2; c2++) {
      (*float7d)[t][c2] = malloc (n3 * sizeof (float****));
      CHKPOINTER ((*float7d)[t][c2]);

      for (c3 = 0; c3 < n3; c3++) {
        (*float7d)[t][c2][c3] = malloc (n4 * sizeof (float***));
        CHKPOINTER ((*float7d)[t][c2][c3]);

        for (c4 = 0; c4 < n4; c4++) {
          (*float7d)[t][c2][c3][c4] = malloc (n5 * sizeof (float**));
          CHKPOINTER ((*float7d)[t][c2][c3][c4]);

          for (c5 = 0; c5 < n5; c5++) {
            (*float7d)[t][c2][c3][c4][c5] = malloc (n6 * sizeof (float*));
            CHKPOINTER ((*float7d)[t][c2][c3][c4][c5]);

            for (j = 0; j < n6; j++) {
              (*float7d)[t][c2][c3][c4][c5][j] = malloc (n7 * sizeof (float));
              CHKPOINTER ((*float7d)[t][c2][c3][c4][c5][j]);
            }
          }
        }
      }
    }
  }

  /* get id for float7d */
  status = nc_inq_varid (ncid, var_name, &id_float7d);
  NCCHKERR (status - NC_NOERR);

  /* allocate space for "tmp_float7d" */
  tmp_float7d = malloc (n1 * n2 * n3 * n4 * n5 * n6 * n7 * sizeof (float));
  CHKPOINTER (tmp_float7d);

  /* read values from netCDF variable */
  status = nc_get_var_float (ncid, id_float7d, tmp_float7d);
  if (status - NC_NOERR) {
    free (tmp_float7d);
    NCCHKERR (status - NC_NOERR);
  }

  for (t = 0; t < n1; t++)
    for (c2 = 0; c2 < n2; c2++)
      for (c3 = 0; c3 < n3; c3++)
        for (c4 = 0; c4 < n4; c4++)
          for (c5 = 0; c5 < n5; c5++)
            for (j = 0; j < n6; j++)
              for (i = 0; i < n7; i++) {
                /* if (done[t][j][i]) { */
                (*float7d)[t][c2][c3][c4][c5][j][i] =
                  tmp_float7d[(((((t * n2 + c2) * n3 + c3) * n4 + c4) * n5 + c5) * n6 + j) * n7 + i];
                /* fprintf (stderr, "t=%4d, lambda=%4d, z=%4d, umu=%4d, y=%4d, x=%4d, float7d=%8.2f \n",t,c2,c3,c4,c5,j,(*float7d)[t][c2][c3][c4][c5][j][i]); */
                /* fflush(stderr); */
                /* } */
              }

  free (tmp_float7d);

#else
  CHKERR (print_have_netcdf());
#endif
  return 0;
}

/***********************************************************************************/
/* Function: read_3d_char                                                          */
/* Description:                                                                    */
/*  Reads a 3D data set from an netCDF file and stores it in an 3D char array      */
/*  will will be accordingly allocated.                                            */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the read netCDF file                        */
/*     char *file         filename of the read netCDF file                         */
/*     char *var_name     short name of the data set, e.g. "lat"                   */
/*     long ntime         number of entries in the 1st dimension                   */
/*     long nlat          number of entries in the 2nd dimension                   */
/*     long nlon          number of entries in the 3rd dimension                   */
/*                                                                                 */
/* Return value:                                                                   */
/*     char ****done      read data                                                */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int deprecate_read_3d_char (int ncid, char* file, char* var_name, char**** done, long ntime, long nlat, long nlon, int verbose) {
  DEPRECATE (2023);
#if HAVE_NETCDF4

  int          status = 0;
  long         t = 0, j = 0, i = 0;
  int          id_done  = 0;
  signed char* tmp_done = NULL;

  if (verbose)
    fprintf (stderr, " ... read_char_3d: var = %-10s, data_type = %-8s, ndim = %3d \n", var_name, "char", 3);

  /* allocate space for "done" */
  *done = malloc ((ntime) * sizeof (char**));
  CHKPOINTER (*done);

  for (t = 0; t < (ntime); t++) {
    (*done)[t] = malloc ((nlat) * sizeof (char*));
    CHKPOINTER ((*done)[t]);

    for (j = 0; j < (nlat); j++) {
      (*done)[t][j] = malloc ((nlon) * sizeof (char));
      CHKPOINTER ((*done)[t][j]);
    }
  }

  /* get id for done */
  status = nc_inq_varid (ncid, "done", &id_done);
  if (status != NC_NOERR) {
    /* this is in most of the cases OK, that's the start of an worldloop job */
    /* fprintf (stderr,"Error, get id for 'done' field in the netCDF result file %s\n", file); */
    return status;
  }

  /* allocate space for "tmp_done" */
  tmp_done = malloc ((ntime) * (nlon) * (nlat) * sizeof (signed char));
  CHKPOINTER (tmp_done);

  /* read values from netCDF variable */
  status = nc_get_var_schar (ncid, id_done, tmp_done);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error '%s' while reading '%s' from %s\n", nc_strerror (status), var_name, file);
    NCCHKERR (status);
  }

  for (t = 0; t < (ntime); t++)
    for (j = 0; j < (nlat); j++)
      for (i = 0; i < (nlon); i++) {
        (*done)[t][j][i] = tmp_done[t * (nlon) * (nlat) + j * (nlon) + i];
        /* if ((*done)[t][j][i]) { */
        /*   fprintf (stderr, "read_3d_char: t=%3ld j=%3ld i=%3ld, done=%d \n",t,j,i,(*done)[t][j][i]); */
        /*   fflush(stderr); } */
      }

  free (tmp_done);

#else
  CHKERR (print_have_netcdf());
#endif
  return 0;
}

// create a new netcdf file, will overwrite an existing file.
// uses netcdf4 format if possible, otherwise generate classic 64bit format
int create_new_netcdf (const char* filename, int* ncid) {
#if HAVE_NETCDF4
  int status = nc_create (filename, NC_CLOBBER | NC_NETCDF4, ncid);
  if (status) {
    status = nc_create (filename, NC_CLOBBER | NC_64BIT_OFFSET, ncid);
  }
  NCCHKERR (status);
  status = nc_enddef (*ncid);
  NCCHKERR (status);
#else
  CHKERR (print_have_netcdf());
#endif
  return 0;
}

// check if file is a valid netcdf file
// returns 0 if valid
int is_valid_netcdf (const char* filename) {
#if HAVE_NETCDF4
  int ncid;
  int status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status)
    return 1;
  status = close_netcdf (ncid);
  if (status)
    return 2;
#else
  CHKERR (print_have_netcdf());
#endif
  return 0;
}

// open netcdf file in read only mode
int open_netcdf (const char* filename, int* ncid) {
#if HAVE_NETCDF4
  int status = nc_open (filename, NC_NOWRITE, ncid);
  NCCHKERR (status);
#else
  CHKERR (print_have_netcdf());
#endif
  return 0;
}

// wrapper around nc_close
int close_netcdf (int ncid) {
#if HAVE_NETCDF4
  int status = nc_close (ncid);
  NCCHKERR (status);
#else
  CHKERR (print_have_netcdf());
#endif
  return 0;
}

/***********************************************************************************/
/* Function: write_netcdf_macros                                                   */
/* Description: Write a variable to netCDF file; variable might be either          */
/*        int, size_t, float or double with specific dimension.                    */
/*        normal write routines expect the data to be in a                         */
/*        pointer of pointer of pointer etc. format.                               */
/*        In these cases the memory is copied in contiguous mem.                   */
/*        It is possible to apply a elemental function on each entry.              */
/*        The `..._sliced` routines allow to define slices of dimensions,          */
/*        i.e. as x-slice = (start, end, step)                                     */
/*        The `..._contiguous` functions are the backbone where the netcdf         */
/*        is actually prepared and written.                                        */
/* Parameters:                                                                     */
/*  int    ncid:          Id of netCDF file give by create_new_netcdf()            */
/*  dtype* data:          Pointer to a variable                                    */
/*  size_t size N:        Size of Nth dimension                                    */
/*  char*  var_name:      Name of variable                                         */
/*  char*  dim N:         Name on Nth dimension                                    */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                         Bettina Richter */
/***********************************************************************************/

#if HAVE_NETCDF4
// write a n-dimensional variable to netcdf. The data needs to be laid out consecutively.
#define WRITE_NETCDF_CONTIGUOUS_ND(name, type, ncfunc, nctype, ncdtypecast)                                                        \
  int write_netcdf_contiguous_ND_##name (int const         ncid,                                                                   \
                                         type* const       data,                                                                   \
                                         const char* const var_name,                                                               \
                                         int const         ndims,                                                                  \
                                         size_t* const     dimsize,                                                                \
                                         char** const      dimname) {                                                                   \
    const int idebug = 0;                                                                                                          \
    CHKERROUT (!ncid, "Invalid netcdf id 0");                                                                                      \
                                                                                                                                   \
    int status = nc_redef (ncid);                                                                                                  \
    if (status != NC_EINDEFINE) {                                                                                                  \
      NCCHKERR (status);                                                                                                           \
    }                                                                                                                              \
                                                                                                                                   \
    int dimids[ndims];                                                                                                             \
    for (size_t idim = 0; idim < ndims; ++idim) {                                                                                  \
      if (idebug) {                                                                                                                \
        fprintf (stderr, "netcdf var <%s> dim %lu <%s> size %lu\n", var_name, idim, dimname[idim], dimsize[idim]);                 \
      }                                                                                                                            \
      if ((nc_inq_dimid (ncid, dimname[idim], &(dimids[idim])) != NC_NOERR)) {                                                     \
        status = nc_def_dim (ncid, dimname[idim], dimsize[idim], &(dimids[idim]));                                                 \
        NCCHKERR (status);                                                                                                         \
      }                                                                                                                            \
    }                                                                                                                              \
                                                                                                                                   \
    int varid;                                                                                                                     \
    status = nc_def_var (ncid, var_name, nctype, ndims, dimids, &varid);                                                           \
    NCCHKERR (status);                                                                                                             \
    status = nc_enddef (ncid);                                                                                                     \
    NCCHKERR (status);                                                                                                             \
    status = nc_put_var_##ncfunc (ncid, varid, (ncdtypecast*) data);                                                               \
    NCCHKERR (status);                                                                                                             \
    return 0;                                                                                                                      \
  }
#else
#define WRITE_NETCDF_CONTIGUOUS_ND(name, type, ncfunc, nctype, ncdtypecast)                                                        \
  int write_netcdf_contiguous_ND_##name (int const         ncid,                                                                   \
                                         type* const       data,                                                                   \
                                         const char* const var_name,                                                               \
                                         int const         ndims,                                                                  \
                                         size_t* const     dimsize,                                                                \
                                         char** const      dimname) {                                                                   \
    CHKERR (print_have_netcdf());                                                                                                  \
    return -1;                                                                                                                     \
  }
#endif

#if HAVE_NETCDF4
#define WRITE_NETCDF_0D(name, type, ncfunc, nctype, ncdtypecast)                                                                   \
  int write_netcdf_0D_##name (const int ncid, const type value, const char* var_name) {                                            \
    CHKERROUT (!ncid, "Invalid netcdf id 0");                                                                                      \
    int varid, status;                                                                                                             \
    status = nc_redef (ncid);                                                                                                      \
    NCCHKERR (status);                                                                                                             \
    status = nc_def_var (ncid, var_name, nctype, 0, 0, &varid);                                                                    \
    NCCHKERR (status);                                                                                                             \
    status = nc_enddef (ncid);                                                                                                     \
    NCCHKERR (status);                                                                                                             \
    status = nc_put_var_##ncfunc (ncid, varid, (ncdtypecast*) &value);                                                             \
    NCCHKERR (status);                                                                                                             \
    return 0;                                                                                                                      \
  }
#else
#define WRITE_NETCDF_0D(name, type, ncfunc, nctype, ncdtypecast)                                                                   \
  int write_netcdf_0D_##name (const int ncid, const type value, const char* var_name) {                                            \
    CHKERR (print_have_netcdf());                                                                                                  \
    return -1;                                                                                                                     \
  }
#endif

#define WRITE_NETCDF_1D(name, type, ncfunc, nctype, ncdtypecast)                                                                   \
  int write_netcdf_1D_##name (int ncid, type* data, size_t N1, char* var_name, char* dim1, type (*elemental_func) (type)) {        \
    size_t slice1[3] = {0, N1, 1};                                                                                                 \
                                                                                                                                   \
    const int status = write_netcdf_1D_##name##_sliced (ncid, data, slice1, var_name, dim1, elemental_func);                       \
    CHKERR (status);                                                                                                               \
                                                                                                                                   \
    return 0;                                                                                                                      \
  }

#define WRITE_NETCDF_1D_SLICED(name, type, ncfunc, nctype, ncdtypecast)                                                            \
  int write_netcdf_1D_##name##_sliced (int    ncid,                                                                                \
                                       type*  value,                                                                               \
                                       size_t N1[3],                                                                               \
                                       char*  var_name,                                                                            \
                                       char*  dim1,                                                                                \
                                       type (*elemental_func) (type)) {                                                            \
                                                                                                                                   \
    const int NDIMS       = 1;                                                                                                     \
    size_t    dimsizes[1] = {0};                                                                                                   \
    for (size_t i = N1[0]; i < N1[1]; i += N1[2], ++dimsizes[0])                                                                   \
      ;                                                                                                                            \
    size_t Nvals;                                                                                                                  \
    int    status = vec_prod_sum_size_t (NDIMS, dimsizes, &Nvals);                                                                 \
    CHKERR (status);                                                                                                               \
                                                                                                                                   \
    type* data = malloc (Nvals * sizeof (type));                                                                                   \
    CHKPOINTER (data);                                                                                                             \
    size_t cnt = 0;                                                                                                                \
    for (size_t i = N1[0]; i < N1[1]; i += N1[2]) {                                                                                \
      data[cnt] = value[i];                                                                                                        \
      ++cnt;                                                                                                                       \
    }                                                                                                                              \
    if (elemental_func) {                                                                                                          \
      for (size_t i = 0; i < cnt; ++i) {                                                                                           \
        data[i] = (*elemental_func) (data[i]);                                                                                     \
      }                                                                                                                            \
    }                                                                                                                              \
                                                                                                                                   \
    char** dimnames = malloc_char_2D (NDIMS, MAX_NAME_LEN, "dim names for netcdf");                                                \
    strcpy (dimnames[0], dim1);                                                                                                    \
    status = write_netcdf_contiguous_ND_##name (ncid, data, var_name, NDIMS, dimsizes, dimnames);                                  \
    CHKERR (status);                                                                                                               \
    free_char_2D (dimnames);                                                                                                       \
    free (data);                                                                                                                   \
                                                                                                                                   \
    return 0;                                                                                                                      \
  }

#define WRITE_NETCDF_2D(name, type, ncfunc, nctype, ncdtypecast)                                                                   \
  int write_netcdf_2D_##name (int    ncid,                                                                                         \
                              type** data,                                                                                         \
                              size_t N1,                                                                                           \
                              size_t N2,                                                                                           \
                              char*  var_name,                                                                                     \
                              char*  dim1,                                                                                         \
                              char*  dim2,                                                                                         \
                              type (*elemental_func) (type)) {                                                                     \
    size_t slice1[3] = {0, N1, 1};                                                                                                 \
    size_t slice2[3] = {0, N2, 1};                                                                                                 \
                                                                                                                                   \
    const int status = write_netcdf_2D_##name##_sliced (ncid, data, slice1, slice2, var_name, dim1, dim2, elemental_func);         \
    CHKERR (status);                                                                                                               \
                                                                                                                                   \
    return 0;                                                                                                                      \
  }

#define WRITE_NETCDF_2D_SLICED(name, type, ncfunc, nctype, ncdtypecast)                                                            \
  int write_netcdf_2D_##name##_sliced (int    ncid,                                                                                \
                                       type** value,                                                                               \
                                       size_t N1[3],                                                                               \
                                       size_t N2[3],                                                                               \
                                       char*  var_name,                                                                            \
                                       char*  dim1,                                                                                \
                                       char*  dim2,                                                                                \
                                       type (*elemental_func) (type)) {                                                            \
                                                                                                                                   \
    const int NDIMS       = 2;                                                                                                     \
    size_t    dimsizes[2] = {0};                                                                                                   \
    for (size_t i = N1[0]; i < N1[1]; i += N1[2], ++dimsizes[0])                                                                   \
      ;                                                                                                                            \
    for (size_t j = N2[0]; j < N2[1]; j += N2[2], ++dimsizes[1])                                                                   \
      ;                                                                                                                            \
    size_t Nvals;                                                                                                                  \
    int    status = vec_prod_sum_size_t (NDIMS, dimsizes, &Nvals);                                                                 \
    CHKERR (status);                                                                                                               \
                                                                                                                                   \
    type* data = malloc (Nvals * sizeof (type));                                                                                   \
    CHKPOINTER (data);                                                                                                             \
    size_t cnt = 0;                                                                                                                \
    for (size_t i = N1[0]; i < N1[1]; i += N1[2]) {                                                                                \
      for (size_t j = N2[0]; j < N2[1]; j += N2[2]) {                                                                              \
        data[cnt] = value[i][j];                                                                                                   \
        ++cnt;                                                                                                                     \
      }                                                                                                                            \
    }                                                                                                                              \
    if (elemental_func) {                                                                                                          \
      for (size_t i = 0; i < cnt; ++i) {                                                                                           \
        data[i] = (*elemental_func) (data[i]);                                                                                     \
      }                                                                                                                            \
    }                                                                                                                              \
                                                                                                                                   \
    char** dimnames = malloc_char_2D (NDIMS, MAX_NAME_LEN, "dim names for netcdf");                                                \
    strcpy (dimnames[0], dim1);                                                                                                    \
    strcpy (dimnames[1], dim2);                                                                                                    \
    status = write_netcdf_contiguous_ND_##name (ncid, data, var_name, NDIMS, dimsizes, dimnames);                                  \
    CHKERR (status);                                                                                                               \
    free_char_2D (dimnames);                                                                                                       \
    free (data);                                                                                                                   \
                                                                                                                                   \
    return 0;                                                                                                                      \
  }

#define WRITE_NETCDF_3D(name, type, ncfunc, nctype, ncdtypecast)                                                                   \
  int write_netcdf_3D_##name (int     ncid,                                                                                        \
                              type*** data,                                                                                        \
                              size_t  N1,                                                                                          \
                              size_t  N2,                                                                                          \
                              size_t  N3,                                                                                          \
                              char*   var_name,                                                                                    \
                              char*   dim1,                                                                                        \
                              char*   dim2,                                                                                        \
                              char*   dim3,                                                                                        \
                              type (*elemental_func) (type)) {                                                                     \
    size_t slice1[3] = {0, N1, 1};                                                                                                 \
    size_t slice2[3] = {0, N2, 1};                                                                                                 \
    size_t slice3[3] = {0, N3, 1};                                                                                                 \
                                                                                                                                   \
    const int status =                                                                                                             \
      write_netcdf_3D_##name##_sliced (ncid, data, slice1, slice2, slice3, var_name, dim1, dim2, dim3, elemental_func);            \
    CHKERR (status);                                                                                                               \
                                                                                                                                   \
    return 0;                                                                                                                      \
  }

#define WRITE_NETCDF_3D_SLICED(name, type, ncfunc, nctype, ncdtypecast)                                                            \
  int write_netcdf_3D_##name##_sliced (int     ncid,                                                                               \
                                       type*** value,                                                                              \
                                       size_t  N1[3],                                                                              \
                                       size_t  N2[3],                                                                              \
                                       size_t  N3[3],                                                                              \
                                       char*   var_name,                                                                           \
                                       char*   dim1,                                                                               \
                                       char*   dim2,                                                                               \
                                       char*   dim3,                                                                               \
                                       type (*elemental_func) (type)) {                                                            \
                                                                                                                                   \
    const int NDIMS       = 3;                                                                                                     \
    size_t    dimsizes[3] = {0};                                                                                                   \
    for (size_t i = N1[0]; i < N1[1]; i += N1[2], ++dimsizes[0])                                                                   \
      ;                                                                                                                            \
    for (size_t j = N2[0]; j < N2[1]; j += N2[2], ++dimsizes[1])                                                                   \
      ;                                                                                                                            \
    for (size_t k = N3[0]; k < N3[1]; k += N3[2], ++dimsizes[2])                                                                   \
      ;                                                                                                                            \
    size_t Nvals;                                                                                                                  \
    int    status = vec_prod_sum_size_t (NDIMS, dimsizes, &Nvals);                                                                 \
    CHKERR (status);                                                                                                               \
                                                                                                                                   \
    type* data = malloc (Nvals * sizeof (type));                                                                                   \
    CHKPOINTER (data);                                                                                                             \
    size_t cnt = 0;                                                                                                                \
    for (size_t i = N1[0]; i < N1[1]; i += N1[2]) {                                                                                \
      for (size_t j = N2[0]; j < N2[1]; j += N2[2]) {                                                                              \
        for (size_t k = N3[0]; k < N3[1]; k += N3[2]) {                                                                            \
          data[cnt] = value[i][j][k];                                                                                              \
          ++cnt;                                                                                                                   \
        }                                                                                                                          \
      }                                                                                                                            \
    }                                                                                                                              \
    if (elemental_func) {                                                                                                          \
      for (size_t i = 0; i < cnt; ++i) {                                                                                           \
        data[i] = (*elemental_func) (data[i]);                                                                                     \
      }                                                                                                                            \
    }                                                                                                                              \
                                                                                                                                   \
    char** dimnames = malloc_char_2D (NDIMS, MAX_NAME_LEN, "dim names for netcdf");                                                \
    strcpy (dimnames[0], dim1);                                                                                                    \
    strcpy (dimnames[1], dim2);                                                                                                    \
    strcpy (dimnames[2], dim3);                                                                                                    \
    status = write_netcdf_contiguous_ND_##name (ncid, data, var_name, NDIMS, dimsizes, dimnames);                                  \
    CHKERR (status);                                                                                                               \
    free_char_2D (dimnames);                                                                                                       \
    free (data);                                                                                                                   \
                                                                                                                                   \
    return 0;                                                                                                                      \
  }

#define WRITE_NETCDF_4D(name, type, ncfunc, nctype, ncdtypecast)                                                                   \
  int write_netcdf_4D_##name (int      ncid,                                                                                       \
                              type**** data,                                                                                       \
                              size_t   N1,                                                                                         \
                              size_t   N2,                                                                                         \
                              size_t   N3,                                                                                         \
                              size_t   N4,                                                                                         \
                              char*    var_name,                                                                                   \
                              char*    dim1,                                                                                       \
                              char*    dim2,                                                                                       \
                              char*    dim3,                                                                                       \
                              char*    dim4,                                                                                       \
                              type (*elemental_func) (type)) {                                                                     \
    size_t slice1[3] = {0, N1, 1};                                                                                                 \
    size_t slice2[3] = {0, N2, 1};                                                                                                 \
    size_t slice3[3] = {0, N3, 1};                                                                                                 \
    size_t slice4[3] = {0, N4, 1};                                                                                                 \
                                                                                                                                   \
    const int status = write_netcdf_4D_##name##_sliced (ncid,                                                                      \
                                                        data,                                                                      \
                                                        slice1,                                                                    \
                                                        slice2,                                                                    \
                                                        slice3,                                                                    \
                                                        slice4,                                                                    \
                                                        var_name,                                                                  \
                                                        dim1,                                                                      \
                                                        dim2,                                                                      \
                                                        dim3,                                                                      \
                                                        dim4,                                                                      \
                                                        elemental_func);                                                           \
    CHKERR (status);                                                                                                               \
                                                                                                                                   \
    return 0;                                                                                                                      \
  }

#define WRITE_NETCDF_4D_SLICED(name, type, ncfunc, nctype, ncdtypecast)                                                            \
  int write_netcdf_4D_##name##_sliced (int      ncid,                                                                              \
                                       type**** value,                                                                             \
                                       size_t   N1[3],                                                                             \
                                       size_t   N2[3],                                                                             \
                                       size_t   N3[3],                                                                             \
                                       size_t   N4[3],                                                                             \
                                       char*    var_name,                                                                          \
                                       char*    dim1,                                                                              \
                                       char*    dim2,                                                                              \
                                       char*    dim3,                                                                              \
                                       char*    dim4,                                                                              \
                                       type (*elemental_func) (type)) {                                                            \
                                                                                                                                   \
    const int NDIMS       = 4;                                                                                                     \
    size_t    dimsizes[4] = {0};                                                                                                   \
    for (size_t i = N1[0]; i < N1[1]; i += N1[2], ++dimsizes[0])                                                                   \
      ;                                                                                                                            \
    for (size_t j = N2[0]; j < N2[1]; j += N2[2], ++dimsizes[1])                                                                   \
      ;                                                                                                                            \
    for (size_t k = N3[0]; k < N3[1]; k += N3[2], ++dimsizes[2])                                                                   \
      ;                                                                                                                            \
    for (size_t l = N4[0]; l < N4[1]; l += N4[2], ++dimsizes[3])                                                                   \
      ;                                                                                                                            \
    size_t Nvals;                                                                                                                  \
    int    status = vec_prod_sum_size_t (NDIMS, dimsizes, &Nvals);                                                                 \
    CHKERR (status);                                                                                                               \
                                                                                                                                   \
    type* data = malloc (Nvals * sizeof (type));                                                                                   \
    CHKPOINTER (data);                                                                                                             \
    size_t cnt = 0;                                                                                                                \
    for (size_t i = N1[0]; i < N1[1]; i += N1[2]) {                                                                                \
      for (size_t j = N2[0]; j < N2[1]; j += N2[2]) {                                                                              \
        for (size_t k = N3[0]; k < N3[1]; k += N3[2]) {                                                                            \
          for (size_t l = N4[0]; l < N4[1]; l += N4[2]) {                                                                          \
            data[cnt] = value[i][j][k][l];                                                                                         \
            ++cnt;                                                                                                                 \
          }                                                                                                                        \
        }                                                                                                                          \
      }                                                                                                                            \
    }                                                                                                                              \
    if (elemental_func) {                                                                                                          \
      for (size_t i = 0; i < cnt; ++i) {                                                                                           \
        data[i] = (*elemental_func) (data[i]);                                                                                     \
      }                                                                                                                            \
    }                                                                                                                              \
                                                                                                                                   \
    char** dimnames = malloc_char_2D (NDIMS, MAX_NAME_LEN, "dim names for netcdf");                                                \
    strcpy (dimnames[0], dim1);                                                                                                    \
    strcpy (dimnames[1], dim2);                                                                                                    \
    strcpy (dimnames[2], dim3);                                                                                                    \
    strcpy (dimnames[3], dim4);                                                                                                    \
    status = write_netcdf_contiguous_ND_##name (ncid, data, var_name, NDIMS, dimsizes, dimnames);                                  \
    CHKERR (status);                                                                                                               \
    free_char_2D (dimnames);                                                                                                       \
    free (data);                                                                                                                   \
                                                                                                                                   \
    return 0;                                                                                                                      \
  }

#define WRITE_NETCDF_5D(name, type, ncfunc, nctype, ncdtypecast)                                                                   \
  int write_netcdf_5D_##name (int       ncid,                                                                                      \
                              type***** data,                                                                                      \
                              size_t    N1,                                                                                        \
                              size_t    N2,                                                                                        \
                              size_t    N3,                                                                                        \
                              size_t    N4,                                                                                        \
                              size_t    N5,                                                                                        \
                              char*     var_name,                                                                                  \
                              char*     dim1,                                                                                      \
                              char*     dim2,                                                                                      \
                              char*     dim3,                                                                                      \
                              char*     dim4,                                                                                      \
                              char*     dim5,                                                                                      \
                              type (*elemental_func) (type)) {                                                                     \
    size_t slice1[3] = {0, N1, 1};                                                                                                 \
    size_t slice2[3] = {0, N2, 1};                                                                                                 \
    size_t slice3[3] = {0, N3, 1};                                                                                                 \
    size_t slice4[3] = {0, N4, 1};                                                                                                 \
    size_t slice5[3] = {0, N5, 1};                                                                                                 \
                                                                                                                                   \
    const int status = write_netcdf_5D_##name##_sliced (ncid,                                                                      \
                                                        data,                                                                      \
                                                        slice1,                                                                    \
                                                        slice2,                                                                    \
                                                        slice3,                                                                    \
                                                        slice4,                                                                    \
                                                        slice5,                                                                    \
                                                        var_name,                                                                  \
                                                        dim1,                                                                      \
                                                        dim2,                                                                      \
                                                        dim3,                                                                      \
                                                        dim4,                                                                      \
                                                        dim5,                                                                      \
                                                        elemental_func);                                                           \
    CHKERR (status);                                                                                                               \
                                                                                                                                   \
    return 0;                                                                                                                      \
  }

#define WRITE_NETCDF_5D_SLICED(name, type, ncfunc, nctype, ncdtypecast)                                                            \
  int write_netcdf_5D_##name##_sliced (int       ncid,                                                                             \
                                       type***** value,                                                                            \
                                       size_t    N1[3],                                                                            \
                                       size_t    N2[3],                                                                            \
                                       size_t    N3[3],                                                                            \
                                       size_t    N4[3],                                                                            \
                                       size_t    N5[3],                                                                            \
                                       char*     var_name,                                                                         \
                                       char*     dim1,                                                                             \
                                       char*     dim2,                                                                             \
                                       char*     dim3,                                                                             \
                                       char*     dim4,                                                                             \
                                       char*     dim5,                                                                             \
                                       type (*elemental_func) (type)) {                                                            \
                                                                                                                                   \
    const int NDIMS       = 5;                                                                                                     \
    size_t    dimsizes[5] = {0};                                                                                                   \
    for (size_t i = N1[0]; i < N1[1]; i += N1[2], ++dimsizes[0])                                                                   \
      ;                                                                                                                            \
    for (size_t j = N2[0]; j < N2[1]; j += N2[2], ++dimsizes[1])                                                                   \
      ;                                                                                                                            \
    for (size_t k = N3[0]; k < N3[1]; k += N3[2], ++dimsizes[2])                                                                   \
      ;                                                                                                                            \
    for (size_t l = N4[0]; l < N4[1]; l += N4[2], ++dimsizes[3])                                                                   \
      ;                                                                                                                            \
    for (size_t m = N5[0]; m < N5[1]; m += N5[2], ++dimsizes[4])                                                                   \
      ;                                                                                                                            \
    ;                                                                                                                              \
    size_t Nvals;                                                                                                                  \
    int    status = vec_prod_sum_size_t (NDIMS, dimsizes, &Nvals);                                                                 \
    CHKERR (status);                                                                                                               \
                                                                                                                                   \
    type* data = malloc (Nvals * sizeof (type));                                                                                   \
    CHKPOINTER (data);                                                                                                             \
    size_t cnt = 0;                                                                                                                \
    for (size_t i = N1[0]; i < N1[1]; i += N1[2]) {                                                                                \
      for (size_t j = N2[0]; j < N2[1]; j += N2[2]) {                                                                              \
        for (size_t k = N3[0]; k < N3[1]; k += N3[2]) {                                                                            \
          for (size_t l = N4[0]; l < N4[1]; l += N4[2]) {                                                                          \
            for (size_t m = N5[0]; m < N5[1]; m += N5[2]) {                                                                        \
              data[cnt] = value[i][j][k][l][m];                                                                                    \
              ++cnt;                                                                                                               \
            }                                                                                                                      \
          }                                                                                                                        \
        }                                                                                                                          \
      }                                                                                                                            \
    }                                                                                                                              \
    if (elemental_func) {                                                                                                          \
      for (size_t i = 0; i < cnt; ++i) {                                                                                           \
        data[i] = (*elemental_func) (data[i]);                                                                                     \
      }                                                                                                                            \
    }                                                                                                                              \
                                                                                                                                   \
    char** dimnames = malloc_char_2D (NDIMS, MAX_NAME_LEN, "dim names for netcdf");                                                \
    strcpy (dimnames[0], dim1);                                                                                                    \
    strcpy (dimnames[1], dim2);                                                                                                    \
    strcpy (dimnames[2], dim3);                                                                                                    \
    strcpy (dimnames[3], dim4);                                                                                                    \
    strcpy (dimnames[4], dim5);                                                                                                    \
    status = write_netcdf_contiguous_ND_##name (ncid, data, var_name, NDIMS, dimsizes, dimnames);                                  \
    CHKERR (status);                                                                                                               \
    free_char_2D (dimnames);                                                                                                       \
    free (data);                                                                                                                   \
                                                                                                                                   \
    return 0;                                                                                                                      \
  }

#define WRITE_NETCDF_6D(name, type, ncfunc, nctype, ncdtypecast)                                                                   \
  int write_netcdf_6D_##name (int        ncid,                                                                                     \
                              type****** data,                                                                                     \
                              size_t     N1,                                                                                       \
                              size_t     N2,                                                                                       \
                              size_t     N3,                                                                                       \
                              size_t     N4,                                                                                       \
                              size_t     N5,                                                                                       \
                              size_t     N6,                                                                                       \
                              char*      var_name,                                                                                 \
                              char*      dim1,                                                                                     \
                              char*      dim2,                                                                                     \
                              char*      dim3,                                                                                     \
                              char*      dim4,                                                                                     \
                              char*      dim5,                                                                                     \
                              char*      dim6,                                                                                     \
                              type (*elemental_func) (type)) {                                                                     \
    size_t slice1[3] = {0, N1, 1};                                                                                                 \
    size_t slice2[3] = {0, N2, 1};                                                                                                 \
    size_t slice3[3] = {0, N3, 1};                                                                                                 \
    size_t slice4[3] = {0, N4, 1};                                                                                                 \
    size_t slice5[3] = {0, N5, 1};                                                                                                 \
    size_t slice6[3] = {0, N6, 1};                                                                                                 \
                                                                                                                                   \
    const int status = write_netcdf_6D_##name##_sliced (ncid,                                                                      \
                                                        data,                                                                      \
                                                        slice1,                                                                    \
                                                        slice2,                                                                    \
                                                        slice3,                                                                    \
                                                        slice4,                                                                    \
                                                        slice5,                                                                    \
                                                        slice6,                                                                    \
                                                        var_name,                                                                  \
                                                        dim1,                                                                      \
                                                        dim2,                                                                      \
                                                        dim3,                                                                      \
                                                        dim4,                                                                      \
                                                        dim5,                                                                      \
                                                        dim6,                                                                      \
                                                        elemental_func);                                                           \
    CHKERR (status);                                                                                                               \
                                                                                                                                   \
    return 0;                                                                                                                      \
  }

#define WRITE_NETCDF_6D_SLICED(name, type, ncfunc, nctype, ncdtypecast)                                                            \
  int write_netcdf_6D_##name##_sliced (int        ncid,                                                                            \
                                       type****** value,                                                                           \
                                       size_t     N1[3],                                                                           \
                                       size_t     N2[3],                                                                           \
                                       size_t     N3[3],                                                                           \
                                       size_t     N4[3],                                                                           \
                                       size_t     N5[3],                                                                           \
                                       size_t     N6[3],                                                                           \
                                       char*      var_name,                                                                        \
                                       char*      dim1,                                                                            \
                                       char*      dim2,                                                                            \
                                       char*      dim3,                                                                            \
                                       char*      dim4,                                                                            \
                                       char*      dim5,                                                                            \
                                       char*      dim6,                                                                            \
                                       type (*elemental_func) (type)) {                                                            \
                                                                                                                                   \
    const int NDIMS       = 6;                                                                                                     \
    size_t    dimsizes[6] = {0};                                                                                                   \
    for (size_t i = N1[0]; i < N1[1]; i += N1[2], ++dimsizes[0])                                                                   \
      ;                                                                                                                            \
    for (size_t j = N2[0]; j < N2[1]; j += N2[2], ++dimsizes[1])                                                                   \
      ;                                                                                                                            \
    for (size_t k = N3[0]; k < N3[1]; k += N3[2], ++dimsizes[2])                                                                   \
      ;                                                                                                                            \
    for (size_t l = N4[0]; l < N4[1]; l += N4[2], ++dimsizes[3])                                                                   \
      ;                                                                                                                            \
    for (size_t m = N5[0]; m < N5[1]; m += N5[2], ++dimsizes[4])                                                                   \
      ;                                                                                                                            \
    for (size_t n = N6[0]; n < N6[1]; n += N6[2], ++dimsizes[5])                                                                   \
      ;                                                                                                                            \
    size_t Nvals;                                                                                                                  \
    int    status = vec_prod_sum_size_t (NDIMS, dimsizes, &Nvals);                                                                 \
    CHKERR (status);                                                                                                               \
                                                                                                                                   \
    type* data = malloc (Nvals * sizeof (type));                                                                                   \
    CHKPOINTER (data);                                                                                                             \
    size_t cnt = 0;                                                                                                                \
    for (size_t i = N1[0]; i < N1[1]; i += N1[2]) {                                                                                \
      for (size_t j = N2[0]; j < N2[1]; j += N2[2]) {                                                                              \
        for (size_t k = N3[0]; k < N3[1]; k += N3[2]) {                                                                            \
          for (size_t l = N4[0]; l < N4[1]; l += N4[2]) {                                                                          \
            for (size_t m = N5[0]; m < N5[1]; m += N5[2]) {                                                                        \
              for (size_t n = N6[0]; n < N6[1]; n += N6[2]) {                                                                      \
                data[cnt] = value[i][j][k][l][m][n];                                                                               \
                ++cnt;                                                                                                             \
              }                                                                                                                    \
            }                                                                                                                      \
          }                                                                                                                        \
        }                                                                                                                          \
      }                                                                                                                            \
    }                                                                                                                              \
    if (elemental_func) {                                                                                                          \
      for (size_t i = 0; i < cnt; ++i) {                                                                                           \
        data[i] = (*elemental_func) (data[i]);                                                                                     \
      }                                                                                                                            \
    }                                                                                                                              \
                                                                                                                                   \
    char** dimnames = malloc_char_2D (NDIMS, MAX_NAME_LEN, "dim names for netcdf");                                                \
    strcpy (dimnames[0], dim1);                                                                                                    \
    strcpy (dimnames[1], dim2);                                                                                                    \
    strcpy (dimnames[2], dim3);                                                                                                    \
    strcpy (dimnames[3], dim4);                                                                                                    \
    strcpy (dimnames[4], dim5);                                                                                                    \
    strcpy (dimnames[5], dim6);                                                                                                    \
    status = write_netcdf_contiguous_ND_##name (ncid, data, var_name, NDIMS, dimsizes, dimnames);                                  \
    CHKERR (status);                                                                                                               \
    free_char_2D (dimnames);                                                                                                       \
    free (data);                                                                                                                   \
                                                                                                                                   \
    return 0;                                                                                                                      \
  }

#if HAVE_NETCDF4
// read a scalar variable
#define READ_NETCDF_SCALAR(name, type, ncfunc, nctype, ncdtypecast)                                                                \
  int read_netcdf_scalar_##name (const int ncid, const char* var_name, type* data) {                                               \
    CHKERROUT (!ncid, "Invalid netcdf id 0");                                                                                      \
    int varid, status;                                                                                                             \
    status = nc_inq_varid (ncid, var_name, &varid);                                                                                \
    NCCHKERR (status);                                                                                                             \
    status = nc_get_var_##ncfunc (ncid, varid, (ncdtypecast*) data);                                                               \
    NCCHKERR (status);                                                                                                             \
    return 0;                                                                                                                      \
  }
#else
#define READ_NETCDF_SCALAR(name, type, ncfunc, nctype, ncdtypecast)                                                                \
  int read_netcdf_scalar_##name (const int ncid, const char* var_name, type* data) {                                               \
    CHKERR (print_have_netcdf());                                                                                                  \
    return -1;                                                                                                                     \
  }
#endif

#if HAVE_NETCDF4
// read a multidimensional variable, if memory is not allocated, we will malloc consecutive mem for the variable.
#define READ_NETCDF(name, type, ncfunc, nctype, ncdtypecast)                                                                       \
  int read_netcdf_##name (const int ncid, const char* var_name, type** data, int* Ndims, size_t* dim_sizes) {                      \
    const int idebug = 0;                                                                                                          \
    CHKERROUT (!ncid, "Invalid netcdf id 0");                                                                                      \
    int varid, status;                                                                                                             \
    status = nc_inq_varid (ncid, var_name, &varid);                                                                                \
    NCCHKERR (status);                                                                                                             \
                                                                                                                                   \
    if (Ndims) {                                                                                                                   \
      status = nc_inq_varndims (ncid, varid, Ndims);                                                                               \
      if (idebug)                                                                                                                  \
        printf ("read_netcdf var <%s> has N=%d dimensions\n", var_name, *Ndims);                                                   \
    }                                                                                                                              \
                                                                                                                                   \
    int  dimids[*Ndims];                                                                                                           \
    char dimname[MAX_NAME_LEN];                                                                                                    \
    status = nc_inq_vardimid (ncid, varid, dimids);                                                                                \
    NCCHKERR (status);                                                                                                             \
    for (size_t idim = 0; idim < *Ndims; ++idim) {                                                                                 \
      status = nc_inq_dim (ncid, dimids[idim], dimname, &(dim_sizes[idim]));                                                       \
      if (idebug)                                                                                                                  \
        printf ("read_netcdf var <%s> investigate dimension %lu <%-20s> size %lu\n", var_name, idim, dimname, dim_sizes[idim]);    \
    }                                                                                                                              \
                                                                                                                                   \
    size_t Nvalues;                                                                                                                \
    status = vec_prod_sum_size_t (*Ndims, dim_sizes, &Nvalues);                                                                    \
    CHKERR (status);                                                                                                               \
    if (!*data) {                                                                                                                  \
      *data = malloc (Nvalues * sizeof (type));                                                                                    \
      if (idebug)                                                                                                                  \
        printf ("read_netcdf var <%s> allocated mem for %lu values\n", var_name, Nvalues);                                         \
    }                                                                                                                              \
                                                                                                                                   \
    status = nc_get_var_##ncfunc (ncid, varid, (ncdtypecast*) (*data));                                                            \
    NCCHKERR (status);                                                                                                             \
    return 0;                                                                                                                      \
  }
#else
#define READ_NETCDF(name, type, ncfunc, nctype, ncdtypecast)                                                                       \
  int read_netcdf_##name (const int ncid, const char* var_name, type** data, int* Ndims, size_t* dim_sizes) {                      \
    CHKERR (print_have_netcdf());                                                                                                  \
    return -1;                                                                                                                     \
  }
#endif

#define NETCDF_ND(name, type, ncfunc, nctype, ncdtypecast)                                                                         \
  WRITE_NETCDF_CONTIGUOUS_ND (name, type, ncfunc, nctype, ncdtypecast)                                                             \
  WRITE_NETCDF_0D (name, type, ncfunc, nctype, ncdtypecast)                                                                        \
  WRITE_NETCDF_1D (name, type, ncfunc, nctype, ncdtypecast)                                                                        \
  WRITE_NETCDF_1D_SLICED (name, type, ncfunc, nctype, ncdtypecast)                                                                 \
  WRITE_NETCDF_2D (name, type, ncfunc, nctype, ncdtypecast)                                                                        \
  WRITE_NETCDF_2D_SLICED (name, type, ncfunc, nctype, ncdtypecast)                                                                 \
  WRITE_NETCDF_3D (name, type, ncfunc, nctype, ncdtypecast)                                                                        \
  WRITE_NETCDF_3D_SLICED (name, type, ncfunc, nctype, ncdtypecast)                                                                 \
  WRITE_NETCDF_4D (name, type, ncfunc, nctype, ncdtypecast)                                                                        \
  WRITE_NETCDF_4D_SLICED (name, type, ncfunc, nctype, ncdtypecast)                                                                 \
  WRITE_NETCDF_5D (name, type, ncfunc, nctype, ncdtypecast)                                                                        \
  WRITE_NETCDF_5D_SLICED (name, type, ncfunc, nctype, ncdtypecast)                                                                 \
  WRITE_NETCDF_6D (name, type, ncfunc, nctype, ncdtypecast)                                                                        \
  WRITE_NETCDF_6D_SLICED (name, type, ncfunc, nctype, ncdtypecast)                                                                 \
  READ_NETCDF_SCALAR (name, type, ncfunc, nctype, ncdtypecast)                                                                     \
  READ_NETCDF (name, type, ncfunc, nctype, ncdtypecast)

NETCDF_ND (int, int, int, NC_INT, int);
NETCDF_ND (size_t, size_t, ulonglong, NC_LONG, unsigned long long);
NETCDF_ND (float, float, float, NC_FLOAT, float);
NETCDF_ND (double, double, double, NC_DOUBLE, double);

#undef WRITE_NETCDF_CONTIGUOUS_ND
#undef WRITE_NETCDF_0D
#undef WRITE_NETCDF_1D
#undef WRITE_NETCDF_1D_SLICED
#undef WRITE_NETCDF_2D
#undef WRITE_NETCDF_2D_SLICED
#undef WRITE_NETCDF_3D
#undef WRITE_NETCDF_3D_SLICED
#undef WRITE_NETCDF_4D
#undef WRITE_NETCDF_4D_SLICED
#undef WRITE_NETCDF_5D
#undef WRITE_NETCDF_5D_SLICED
#undef WRITE_NETCDF_6D
#undef WRITE_NETCDF_6D_SLICED
#undef READ_NETCDF_SCALAR
#undef READ_NETCDF
#undef NETCDF_ND

// special cases where ragged arrays are copied
int write_netcdf_3Dirr_row_float (int      ncid,
                                  float*** value,
                                  int      rows,
                                  int      columns,
                                  int*     length,
                                  char*    var_name,
                                  char*    row_name,
                                  char*    column_name,
                                  char*    length_name) {
#if HAVE_NETCDF4

  CHKERROUT (!ncid, "Invalid netcdf id 0");

  if (rows < 0)
    CHKERROUT (rows, "need to positive number of rows");
  if (columns < 0)
    CHKERROUT (columns, "need to positive number of columns");

  /**** Maximum length of value; in case of regular grid, length_max == length ****/
  int length_max = 0;
  for (size_t row = 0; row < rows; row++) {
    length_max = MAX (length_max, length[row]);
  }

  float* data = malloc (rows * columns * length_max * sizeof (float));
  CHKPOINTER (data);
  size_t k = 0;
  for (size_t row = 0; row < rows; row++) {
    for (size_t column = 0; column < columns; column++) {
      for (size_t il = 0; il < length_max; il++, k++) {
        if (il < length[row])
          data[k] = value[row][column][il];
        else
          data[k] = NOT_DEFINED_FLOAT;
      }
    }
  }

  const int NDIMS       = 3;
  size_t    dimsizes[3] = {rows, columns, length_max};
  char**    dimnames    = malloc_char_2D (NDIMS, MAX_NAME_LEN, "dim names for netcdf");
  strcpy (dimnames[0], row_name);
  strcpy (dimnames[1], column_name);
  strcpy (dimnames[2], length_name);
  int status = write_netcdf_contiguous_ND_float (ncid, data, var_name, NDIMS, dimsizes, dimnames);
  CHKERR (status);
  free_char_2D (dimnames);
  UVSPEC_FREE (data);

#else
  CHKERR (print_have_netcdf());
#endif
  return 0;
}

// special cases where ragged arrays are copied
int write_netcdf_3Dirr_float (int      ncid,
                              float*** value,
                              int      rows,
                              int      columns,
                              int**    length,
                              char*    var_name,
                              char*    row_name,
                              char*    column_name,
                              char*    length_name) {
#if HAVE_NETCDF4

  CHKERROUT (!ncid, "Invalid netcdf id 0");

  /**** Maximum length of value; in case of regular grid, length_max == length ****/
  int length_max = 0;
  for (size_t row = 0; row < rows; row++) {
    for (size_t col = 0; col < columns; col++) {
      length_max = MAX (length_max, length[row][col]);
    }
  }

  if (length_max <= 0) {
    fprintf (stderr,
             "Warning: did not write %s to netcdf_file because max_length of values in last dim was %d\n",
             var_name,
             length_max);
    return 0;
  }

  float* data = malloc (rows * columns * length_max * sizeof (float));
  CHKPOINTER (data);
  size_t k = 0;
  for (size_t row = 0; row < rows; row++) {
    for (size_t col = 0; col < columns; col++) {
      for (size_t il = 0; il < length_max; il++, k++) {
        if (il < length[row][col])
          data[k] = value[row][col][il];
        else
          data[k] = NOT_DEFINED_FLOAT;
      }
    }
  }

  const int NDIMS       = 3;
  size_t    dimsizes[3] = {rows, columns, length_max};
  char**    dimnames    = malloc_char_2D (NDIMS, MAX_NAME_LEN, "dim names for netcdf");
  strcpy (dimnames[0], row_name);
  strcpy (dimnames[1], column_name);
  strcpy (dimnames[2], length_name);
  int status = write_netcdf_contiguous_ND_float (ncid, data, var_name, NDIMS, dimsizes, dimnames);
  CHKERR (status);
  free_char_2D (dimnames);
  UVSPEC_FREE (data);

#else
  CHKERR (print_have_netcdf());
#endif
  return 0;
}

// special cases where ragged arrays are copied
int write_netcdf_3Dirr_double (int       ncid,
                               double*** value,
                               int       rows,
                               int       columns,
                               int**     length,
                               char*     var_name,
                               char*     row_name,
                               char*     column_name,
                               char*     length_name) {
#if HAVE_NETCDF4
  CHKERROUT (!ncid, "Invalid netcdf id 0");

  /**** Maximum length of value; in case of regular grid, length_max == length ****/
  int length_max = 0;
  for (size_t row = 0; row < rows; row++) {
    for (size_t col = 0; col < columns; col++) {
      length_max = MAX (length_max, length[row][col]);
    }
  }
  if (length_max <= 0) {
    fprintf (stderr, "WARNING: Did not write %s to netcdf_file because it is empty\n", var_name);
    return 0;
  }
  double* data = malloc (rows * columns * length_max * sizeof (double));
  CHKPOINTER (data);
  size_t k = 0;
  for (size_t row = 0; row < rows; row++) {
    for (size_t col = 0; col < columns; col++) {
      for (size_t il = 0; il < length_max; il++, k++) {
        if (il < length[row][col])
          data[k] = value[row][col][il];
        else
          data[k] = NOT_DEFINED_FLOAT;
      }
    }
  }

  const int NDIMS       = 3;
  size_t    dimsizes[3] = {rows, columns, length_max};
  char**    dimnames    = malloc_char_2D (NDIMS, MAX_NAME_LEN, "dim names for netcdf");
  strcpy (dimnames[0], row_name);
  strcpy (dimnames[1], column_name);
  strcpy (dimnames[2], length_name);
  int status = write_netcdf_contiguous_ND_double (ncid, data, var_name, NDIMS, dimsizes, dimnames);
  CHKERR (status);
  free_char_2D (dimnames);
  UVSPEC_FREE (data);

#else
  CHKERR (print_have_netcdf());
#endif
  return 0;
}

int write_netcdf_attribute_text (const int ncid, const char* var_name, const char* attr_name, char* attr_value) {
#if HAVE_NETCDF4
  CHKERROUT (!ncid, "Invalid netcdf id 0");
  int status;
  int var_id;
  if (var_name) {
    status = nc_inq_varid (ncid, var_name, &var_id);
    NCCHKERR (status);
  } else {
    var_id = NC_GLOBAL;
  }

  status = nc_redef (ncid);
  NCCHKERR (status);
  status = nc_put_att_text (ncid, var_id, attr_name, strlen (attr_value), attr_value);
  NCCHKERR (status);
  status = nc_enddef (ncid);
  NCCHKERR (status);
#else
  CHKERR (print_have_netcdf());
#endif
  return 0;
}

int read_netcdf_attribute_text (const int ncid, const char* var_name, const char* attr_name, char** attr_value) {
#if HAVE_NETCDF4
  CHKERROUT (!ncid, "Invalid netcdf id 0");
  int status;
  int var_id;
  if (var_name) {
    status = nc_inq_varid (ncid, var_name, &var_id);
    NCCHKERR (status);
  } else {
    var_id = NC_GLOBAL;
  }

  size_t attr_len;
  status = nc_inq_attlen (ncid, var_id, attr_name, &attr_len);
  NCCHKERR (status);

  *attr_value = malloc ((attr_len + 1) * sizeof (char));  // make extra space to append terminator

  status = nc_get_att_text (ncid, var_id, attr_name, *attr_value);
  NCCHKERR (status);

  (*attr_value)[attr_len] = '\0';  // append terminator
#else
  CHKERR (print_have_netcdf());
#endif
  return 0;
}

int read_netcdf_attribute_double (const int ncid, const char* var_name, const char* attr_name, double** attr_value) {
#if HAVE_NETCDF4
  CHKERROUT (!ncid, "Invalid netcdf id 0");
  int status;
  int var_id;
  if (var_name) {
    status = nc_inq_varid (ncid, var_name, &var_id);
    NCCHKERR (status);
  } else {
    var_id = NC_GLOBAL;
  }

  if (!(*attr_value)) {  // allocate space if attr_value is NULL
    size_t attr_len;
    status = nc_inq_attlen (ncid, var_id, attr_name, &attr_len);
    NCCHKERR (status);
    *attr_value = malloc (attr_len * sizeof (double));
  }

  status = nc_get_att_double (ncid, var_id, attr_name, *attr_value);
  NCCHKERR (status);

#else
  CHKERR (print_have_netcdf());
#endif
  return 0;
}
