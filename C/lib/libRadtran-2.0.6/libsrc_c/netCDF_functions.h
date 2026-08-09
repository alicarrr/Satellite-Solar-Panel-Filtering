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

#ifndef __netcdf_functions_h
#define __netcdf_functions_h

#if defined(__cplusplus)
extern "C" {
#endif

#include <stddef.h>

#define FALSE 0
#define TRUE 1

#define POSITION_TIME_3D 0
#define POSITION_LAT_3D 1
#define POSITION_LON_3D 2

#define POSITION_TIME_4D 0
#define POSITION_ZOUT_4D 1
#define POSITION_LAT_4D 2
#define POSITION_LON_4D 3

#define POSITION_TIME_5D 0
#define POSITION_LAMBDA_5D 1
#define POSITION_ZOUT_5D 2
#define POSITION_LAT_5D 3
#define POSITION_LON_5D 4

#define POSITION_TIME_6D 0
#define POSITION_LAMBDA_6D 1
#define POSITION_ZOUT_6D 2
#define POSITION_UMU_6D 3
#define POSITION_LAT_6D 4
#define POSITION_LON_6D 5

#define POSITION_TIME_7D 0
#define POSITION_LAMBDA_7D 1
#define POSITION_ZOUT_7D 2
#define POSITION_UMU_7D 3
#define POSITION_PHI_7D 4
#define POSITION_LAT_7D 5
#define POSITION_LON_7D 6

#define NOT_DEFINED_FLOAT -999.0
#define NOT_DEFINED_INTEGER -999

int deprecate_write_netCDF_1D_float (int    ncid,
                                     char*  filename,
                                     char*  data_name,
                                     int*   id_data,
                                     char*  long_name,
                                     char*  units,
                                     float* data,
                                     char*  dim_name,
                                     int*   id_dim,
                                     int    n);

int deprecate_write_netCDF_3D_float (int     ncid,
                                     char*   filename,
                                     char*   data_name,
                                     char*   long_name,
                                     char*   units,
                                     float   data,
                                     int*    id_dim_array,
                                     size_t* index_3D,
                                     int     ndim);

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
                                     int*    simulation_done);

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
                                     int*    simulation_done);

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
                                     int*     simulation_done);

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
                                     int*      simulation_done);

int deprecate_create_float_variable_or_get_id (int   ncid,
                                               char* filename,
                                               char* data_name,
                                               int*  id_data,
                                               int   ndim,
                                               int*  id_dim_array,
                                               char* long_name,
                                               char* units);

int deprecate_write_netCDF_3D_byte (int     ncid,
                                    char*   filename,
                                    char*   data_name,
                                    char*   long_name,
                                    char*   units,
                                    short   data,
                                    int*    id_dim_array,
                                    size_t* index_3D,
                                    int     ndim);

int deprecate_read_3d_float (int       ncid,
                             char*     file,
                             char*     var_name,
                             float**** float3d,
                             long      n1,
                             long      n2,
                             long      n3,
                             char***   done,
                             int       verbose);

int deprecate_read_4d_float (int        ncid,
                             char*      file,
                             char*      var_name,
                             float***** float4d,
                             long       n1,
                             long       n2,
                             long       n3,
                             long       n4,
                             char***    done,
                             int        verbose);

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
                             int         verbose);

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
                             int          verbose);

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
                             int           verbose);
int deprecate_read_3d_char (int ncid, char* file, char* var_name, char**** done, long ntime, long nlat, long nlon, int verbose);

int alloc_and_read_netCDF_1D_double (int ncid, char* dimension, size_t* dim, char* variable, double** data);

int create_new_netcdf (const char* filename, int* ncid);

int is_valid_netcdf (const char* filename);

int open_netcdf (const char* filename, int* ncid);

int close_netcdf (int ncid);

// macro set for multidimensional and multi datatype netcdf routines
#define NETCDF_ND(name, type, ncfunc, nctype)                                                                                      \
  int write_netcdf_contiguous_ND_##name (int const         ncid,                                                                   \
                                         type* const       data,                                                                   \
                                         const char* const var_name,                                                               \
                                         int const         ndims,                                                                  \
                                         size_t* const     dimsize,                                                                \
                                         char** const      dimname);                                                                    \
  int write_netcdf_0D_##name (const int ncid, const type value, const char* var_name);                                             \
                                                                                                                                   \
  int write_netcdf_1D_##name (int ncid, type* value, size_t N1, char* var_name, char* dim1, type (*elemental_func) (type));        \
                                                                                                                                   \
  int write_netcdf_1D_##name##_sliced (int    ncid,                                                                                \
                                       type*  value,                                                                               \
                                       size_t N1[3],                                                                               \
                                       char*  var_name,                                                                            \
                                       char*  dim1,                                                                                \
                                       type (*elemental_func) (type));                                                             \
                                                                                                                                   \
  int write_netcdf_2D_##name (int    ncid,                                                                                         \
                              type** value,                                                                                        \
                              size_t N1,                                                                                           \
                              size_t N2,                                                                                           \
                              char*  var_name,                                                                                     \
                              char*  dim1,                                                                                         \
                              char*  dim2,                                                                                         \
                              type (*elemental_func) (type));                                                                      \
                                                                                                                                   \
  int write_netcdf_2D_##name##_sliced (int    ncid,                                                                                \
                                       type** value,                                                                               \
                                       size_t N1[3],                                                                               \
                                       size_t N2[3],                                                                               \
                                       char*  var_name,                                                                            \
                                       char*  dim1,                                                                                \
                                       char*  dim2,                                                                                \
                                       type (*elemental_func) (type));                                                             \
                                                                                                                                   \
  int write_netcdf_3D_##name (int     ncid,                                                                                        \
                              type*** value,                                                                                       \
                              size_t  N1,                                                                                          \
                              size_t  N2,                                                                                          \
                              size_t  N3,                                                                                          \
                              char*   var_name,                                                                                    \
                              char*   dim1,                                                                                        \
                              char*   dim2,                                                                                        \
                              char*   dim3,                                                                                        \
                              type (*elemental_func) (type));                                                                      \
                                                                                                                                   \
  int write_netcdf_3D_##name##_sliced (int     ncid,                                                                               \
                                       type*** value,                                                                              \
                                       size_t  N1[3],                                                                              \
                                       size_t  N2[3],                                                                              \
                                       size_t  N3[3],                                                                              \
                                       char*   var_name,                                                                           \
                                       char*   dim1,                                                                               \
                                       char*   dim2,                                                                               \
                                       char*   dim3,                                                                               \
                                       type (*elemental_func) (type));                                                             \
                                                                                                                                   \
  int write_netcdf_4D_##name (int      ncid,                                                                                       \
                              type**** value,                                                                                      \
                              size_t   N1,                                                                                         \
                              size_t   N2,                                                                                         \
                              size_t   N3,                                                                                         \
                              size_t   N4,                                                                                         \
                              char*    var_name,                                                                                   \
                              char*    dim1,                                                                                       \
                              char*    dim2,                                                                                       \
                              char*    dim3,                                                                                       \
                              char*    dim4,                                                                                       \
                              type (*elemental_func) (type));                                                                      \
                                                                                                                                   \
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
                                       type (*elemental_func) (type));                                                             \
                                                                                                                                   \
  int write_netcdf_5D_##name (int       ncid,                                                                                      \
                              type***** value,                                                                                     \
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
                              type (*elemental_func) (type));                                                                      \
                                                                                                                                   \
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
                                       type (*elemental_func) (type));                                                             \
                                                                                                                                   \
  int write_netcdf_6D_##name (int        ncid,                                                                                     \
                              type****** value,                                                                                    \
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
                              type (*elemental_func) (type));                                                                      \
                                                                                                                                   \
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
                                       type (*elemental_func) (type));                                                             \
                                                                                                                                   \
  int read_netcdf_scalar_##name (const int ncid, const char* var_name, type* data);                                                \
  int read_netcdf_##name (const int ncid, const char* var_name, type** data, int* Ndims, size_t* dim_sizes);

// Generate ND routines for various datatypes
NETCDF_ND (int, int, int, NC_INT);
NETCDF_ND (size_t, size_t, ulonglong, NC_LONG);
NETCDF_ND (float, float, float, NC_FLOAT);
NETCDF_ND (double, double, double, NC_DOUBLE);

#undef NETCDF_ND

int write_netcdf_3Dirr_float (int      ncid,
                              float*** value,
                              int      rows,
                              int      columns,
                              int**    length,
                              char*    var_name,
                              char*    row_name,
                              char*    column_name,
                              char*    length_name);

int write_netcdf_3Dirr_row_float (int      ncid,
                                  float*** value,
                                  int      rows,
                                  int      columns,
                                  int*     length,
                                  char*    var_name,
                                  char*    row_name,
                                  char*    column_name,
                                  char*    length_name);

int write_netcdf_3Dirr_double (int       ncid,
                               double*** value,
                               int       rows,
                               int       columns,
                               int**     length,
                               char*     var_name,
                               char*     row_name,
                               char*     column_name,
                               char*     length_name);

int write_netcdf_attribute_text (const int ncid, const char* var_name, const char* attr_name, char* attr_value);
int read_netcdf_attribute_text (const int ncid, const char* var_name, const char* attr_name, char** attr_value);
int read_netcdf_attribute_double (const int ncid, const char* var_name, const char* attr_name, double** attr_value);

#if defined(__cplusplus)
}
#endif

#endif
