#include <netCDF_surface.h>

#if HAVE_NETCDF4
#include <netcdf.h>
#endif

// Needed for colored debugging lines
#if defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#endif

#include <stdbool.h>

// Some error functions similar to those in errors.h
// CHKERROUTV, ASSERTOUTV, INFO and INFOV could be moved to errors.h
// NC_CHKERR and NC_CHKERROUTV are specific to NetCDF errors

// Prefix for debugging or error messages (e.g. a blue "[dbg] ")
// (only on UNIX and only on TTYs)
#define PFX_DEBUG "[dbg] "
#define PFX_ERROR "[err] "
#define PFX_WARN "[warn] "

#define TERM_BOLD(str) "\033[1m" str "\033[0m"
#define TERM_BLUE(str) "\033[34m" str "\033[0m"
#define TERM_GREEN(str) "\033[32m" str "\033[0m"
#define TERM_RED(str) "\033[31m" str "\033[0m"
#define TERM_YELLOW(str) "\033[33m" str "\033[0m"

#define PFX_DEBUG_COLOR TERM_BOLD (TERM_BLUE (PFX_DEBUG))
#define PFX_ERROR_COLOR TERM_BOLD (TERM_RED (PFX_ERROR))
#define PFX_WARN_COLOR TERM_BOLD (TERM_YELLOW (PFX_WARN))

#if defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
#define TERM_PFX(pfx, suffix) (isatty (1) ? (pfx##_COLOR suffix) : (pfx suffix))
#else
#define TERM_PFX(pfx, suffix) (pfx suffix)
#endif

/**
 * @brief Macro to check error status code and print a meaningful error.
 * @author Cedric Wehrum
 * @date Nov 20 2023
 * This is a variadic macro to check error status codes and print a meaningful
 * error with a custom message that follows `printf` formatting rules.
 * 
 * @param[in] int_status  Error status code
 * @param[in] description Format string
 * @param[in] ...         Additional arguments used for `printf`.
 */
#define CHKERROUTV(int_status, description, ...)                                                                                   \
  if (int_status) {                                                                                                                \
    fprintf (stderr,                                                                                                               \
             TERM_PFX (PFX_ERROR, "Code %d in (function '%s', file '%s', line %d) with reason: " description "\n"),                \
             int_status,                                                                                                           \
             __func__,                                                                                                             \
             __FILE__,                                                                                                             \
             __LINE__,                                                                                                             \
             __VA_ARGS__);                                                                                                         \
    return int_status;                                                                                                             \
  }

/**
 * @brief Macro to check assertion and print a meaningful error.
 * @author Cedric Wehrum
 * @date Nov 20 2023
 * This variadic macro checks an assertion and prints a meaningful
 * error with a custom message that follows `printf` formatting rules.
 * 
 * @param[in] assertion   Error status code
 * @param[in] error_code  An error code that describes an assertion failure
 * @param[in] description Format string
 * @param[in] ...         Additional arguments used for `printf`.
 */
#define ASSERTOUTV(assertion, error_code, description, ...)                                                                        \
  if (!(assertion)) {                                                                                                              \
    fprintf (stderr,                                                                                                               \
             TERM_PFX (PFX_ERROR,                                                                                                  \
                       "Assertion failed: " #assertion                                                                             \
                       "\nError %d in (function '%s', file '%s', line %d) with reason: " description "\n"),                        \
             error_code,                                                                                                           \
             __func__,                                                                                                             \
             __FILE__,                                                                                                             \
             __LINE__,                                                                                                             \
             __VA_ARGS__);                                                                                                         \
    return error_code;                                                                                                             \
  }

/**
 * @brief Macro to check assertion and print a meaningful error.
 * @author Cedric Wehrum
 * @date Dez 12 2023
 * This macro checks an assertion and prints a meaningful error with a custom
 * message without any formatting.
 * 
 * @param[in] assertion   Error status code
 * @param[in] error_code  An error code that describes an assertion failure
 * @param[in] description Format string
 * @param[in] ...         Additional arguments used for `printf`.
 */

#define ASSERTOUT(assertion, error_code, description)                                                                              \
  if (!(assertion)) {                                                                                                              \
    fprintf (stderr,                                                                                                               \
             TERM_PFX (PFX_ERROR,                                                                                                  \
                       "Assertion failed: " #assertion                                                                             \
                       "\nError %d in (function '%s', file '%s', line %d) with reason: " description "\n"),                        \
             error_code,                                                                                                           \
             __func__,                                                                                                             \
             __FILE__,                                                                                                             \
             __LINE__);                                                                                                            \
    return error_code;                                                                                                             \
  }

/**
 * @brief Prints a warning.
 * @author Cedric Wehrum
 * @date Nov 20 2023
 * 
 * @param[in] description Warning string
 */
#define WARN(description) fprintf (stderr, TERM_PFX (PFX_WARN, description "\n"));

/**
 * @brief Prints a warning.
 * @author Cedric Wehrum
 * @date Nov 20 2023
 * 
 * @param[in] description Warning string
 * @param[in] ...         Additional arguments used for `printf`.
 */
#define WARNV(description, ...) fprintf (stderr, TERM_PFX (PFX_WARN, description "\n"), __VA_ARGS__);

/**
 * @brief Prints info messages.
 * @author Cedric Wehrum
 * @date Nov 20 2023
 * 
 * @param[in] quiet       Wether logging level allows info messages
 * @param[in] description Info string
 */
#define INFO(quiet, description)                                                                                                   \
  if (!quiet) {                                                                                                                    \
    fprintf (stderr, TERM_PFX (PFX_DEBUG, description "\n"));                                                                      \
  }

/**
 * @brief Prints info messages.
 * @author Cedric Wehrum
 * @date Nov 20 2023
 * 
 * @param[in] quiet       Wether logging level allows info messages
 * @param[in] description Info format string
 * @param[in] ...         Additional arguments used for `printf`.
 */
#define INFOV(quiet, description, ...)                                                                                             \
  if (!quiet) {                                                                                                                    \
    fprintf (stderr, TERM_PFX (PFX_DEBUG, description "\n"), __VA_ARGS__);                                                         \
  }

#if HAVE_NETCDF4

/**
 * @brief Macro to check NetCDF status code and print a meaningful error.
 * @author Cedric Wehrum
 * @date Nov 20 2023
 *
 * @param[in] nc_status NetCDF status code
 */
#define NC_CHKERR(nc_status)                                                                                                       \
  if (nc_status != NC_NOERR) {                                                                                                     \
    fprintf (stderr,                                                                                                               \
             TERM_PFX (PFX_ERROR, "NetCDF Error %d in (function '%s', file '%s', line %d) with reason: %s\n"),                     \
             nc_status,                                                                                                            \
             func,                                                                                                                 \
             filename,                                                                                                             \
             line,                                                                                                                 \
             nc_strerror (nc_status));                                                                                             \
    return NCDF_NETCDF_ERROR;                                                                                                      \
  }

/**
 * @brief Macro to check NetCDF status code and print a meaningful error.
 * @author Cedric Wehrum
 * @date Nov 20 2023
 * This variadic macro checks NetCDF status codes and prints a meaningful
 * error with a custom message that follows `printf` formatting rules.
 * 
 * @param[in] nc_status   NetCDF status code
 * @param[in] description Format string
 * @param[in] ...         Additional arguments used for `printf`.
 */
#define NC_CHKERROUTV(nc_status, description, ...)                                                                                 \
  if (nc_status != NC_NOERR) {                                                                                                     \
    fprintf (stderr,                                                                                                               \
             TERM_PFX (PFX_ERROR, "NetCDF Error %d in (function '%s', file '%s', line %d) with reason: %s\n"),                     \
             nc_status,                                                                                                            \
             __func__,                                                                                                             \
             __FILE__,                                                                                                             \
             __LINE__,                                                                                                             \
             nc_strerror (nc_status));                                                                                             \
    fprintf (stderr, TERM_PFX (PFX_ERROR, "Note: " description "\n"), __VA_ARGS__);                                                \
    return NCDF_NETCDF_ERROR;                                                                                                      \
  }

/**
 * @brief Macro to check NetCDF status code and print a meaningful error.
 * @author Cedric Wehrum
 * @date Nov 20 2023
 * This variadic macro checks NetCDF status codes and prints a meaningful
 * error with a custom message that follows `printf` formatting rules.
 * 
 * @param[in] nc_status   NetCDF status code
 * @param[in] description Format string
 * @param[in] ...         Additional arguments used for `printf`.
 */
#define NC_CHKERROUT(nc_status, description)                                                                                       \
  if (nc_status != NC_NOERR) {                                                                                                     \
    fprintf (stderr,                                                                                                               \
             TERM_PFX (PFX_ERROR, "NetCDF Error %d in (function '%s', file '%s', line %d) with reason: %s\n"),                     \
             nc_status,                                                                                                            \
             __func__,                                                                                                             \
             __FILE__,                                                                                                             \
             __LINE__,                                                                                                             \
             nc_strerror (nc_status));                                                                                             \
    fprintf (stderr, TERM_PFX (PFX_ERROR, "Note: " description "\n"));                                                             \
    return NCDF_NETCDF_ERROR;                                                                                                      \
  }

// defines default units
#define TGT_UNIT_ANY -1
#define TGT_UNIT_DEG_EAST 0
#define TGT_UNIT_DEG_NORTH 1
#define TGT_UNIT_NM 2
// type used for units
typedef int tgt_unit_t;

// defines types (uses) of dimensions
#define DIM_UNKNOWN -1
#define DIM_X 0
#define DIM_Y 1
#define DIM_SPEC 2
#define NUM_DIMENSIONS 3
// type used to define usage of a dimension
typedef int dim_t;

/**
 * @brief Finds the albedo method that corresponds to the given name
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * 
 * @param[in] model_name Name of the albedo method (See `mc_model_str`)
 * @return               Albedo method as is defined in `mystic_definitions.h`
 */
mc_surface_method find_surface_model_id (const char* model_name);

/**
 * @brief Finds the parameter id that corresponds to the given name
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * 
 * @param[in] model_name Name of the albedo method (See `mc_model_str`)
 * @param[in] param_name Name of the parameter (See `mc_model_param_str`)
 * @return               Parameter id as is defined in `mystic_definitions.h`
 */
mc_surface_parameter find_model_param_id (const char* model_name, const char* param_name);

/**
 * @brief Successively calculates the stride parameters for an array
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * Assume a 3d array has the shape a[x=10][y=20][wl=30] and that its values are
 * saved contiguously in a flat array (i.e. `float* data;`) as is the case for
 * NetCDF arrays. Also assume that we have coordinates `i`, `j` and `k` to
 * access this array. Then this array can be accessed via
 * `data[i*y*wl + j*wl + k]`. To simplify the array access we define the stride
 * parameters `stride_x=y*wl`, `stride_y=wl` and `stride_wl=1`. This function
 * provides a convenient way to calculate these strides for an arbitrary order
 * of dimensions. Example:
 * mc_ncdf_array_stride strides;
 * reset_stride (&strides);
 * calc_stride (&strides.x, &strides, 10);
 * calc_stride (&strides.y, &strides, 20);
 * calc_stride (&strides.wl, &strides, 30);
 *
 * @param[in]     stride  Stride parameter to be updated. This *must* point to a
 *                        member of `strides` (The next function argument)
 * @param[in,out] strides The stride struct to be updated
 * @param[in]     length  The length of the dimension associated with `stride`
 */
void calc_stride (size_t* stride, mc_ncdf_array_stride* strides, size_t length);

/**
 * @brief Resets stride parameters
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * Resets a stride struct to be used by `calc_stride`. See `calc_stride` for
 * more information.
 * 
 * @param[out] stride The stride struct to be updated
 */
void reset_stride (mc_ncdf_array_stride* stride);

/**
 * @brief Retrieves parameter from an `mc_ncdf_array`
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * 
 * @param[in] array The array that will be accessed
 * @param[in] i     Longitude or x coordinate
 * @param[in] j     Latitude or y coordinate
 * @param[in] wl    Wavelength coordinate
 * @param[in] param Internal parameter index (see mc_ncdf_array.param_indices)
 * @return          Value at position (i, j, wl, param)
 */
float get_param (mc_ncdf_array* array, size_t i, size_t j, size_t wl, size_t param);

/**
 * @brief Retrieves parameter from a C array using a stride struct
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * 
 * @param[in] stride The stride struct used to access `values`
 * @param[in] values The C array that will be accessed
 * @param[in] i      Longitude or x coordinate
 * @param[in] j      Latitude or y coordinate
 * @param[in] wl     Wavelength coordinate
 * @param[in] param  Internal parameter index (see mc_ncdf_array.param_indices)
 * @return           Value at position (i, j, wl, param)
 */
float get_param_by_stride (mc_ncdf_array_stride* stride, float* values, size_t i, size_t j, size_t wl, size_t param);

/**
 * @brief Sets a parameter from an `mc_ncdf_array`
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * 
 * @param[in,out] array The array that will be accessed
 * @param[in]     i     Longitude or x coordinate
 * @param[in]     j     Latitude or y coordinate
 * @param[in]     wl    Wavelength coordinate
 * @param[in]     param Internal parameter index
 *                      (see mc_ncdf_array.param_indices)
 * @param[in]     value New value in position (i, j, wl, param)
 * @return              NCDF_NOERR (i.e. 0) on success, error code on failure
 */
ncdf_status_t set_param (mc_ncdf_array* array, size_t i, size_t j, size_t wl, size_t param, float value);

/**
 * @brief Checks if a given type belongs to one of the integer types defined by
 *        NetCDF
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * 
 * @param[in] type A NetCDF type
 * @return         true if `type` is an integer
 */
bool is_integer (nc_type type);

/**
 * @brief Checks if a given type belongs to one of the float types defined by
 *        NetCDF
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * 
 * @param[in] type A NetCDF type
 * @return         true if `type` is a floating point type
 */
bool is_float (nc_type type);

/**
 * @brief Checks if a given type belongs to one of the string types defined by
 *        NetCDF
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * 
 * @param[in] type A NetCDF type
 * @return         true if `type` is a string
 */
bool is_string (nc_type type);

/**
 * @brief Frees all resources associated with an `mc_ncdf_array`
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * 
 * @param[in,out] array An `mc_ncdf_array`
 */
void free_array (mc_ncdf_array* array);

/**
 * @brief Frees all raw values in an `mc_ncdf_surface`
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * This frees all resources in a surface that are no longer needed after
 * interpolating to the internal wavelength grid
 * 
 * @param[in,out] surface An `mc_ncdf_surface`
 */
void free_raw_values (mc_ncdf_surface* surface);

/**
 * @brief Adds a parameter to an `mc_ncdf_surface`
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * Given a C array of parameters and a stride struct to define its dimensions
 * this function adds these values to one of the arrays of a surface.
 * If possible it only copies a pointer. In case the values have to be copied
 * because the C array has a different structure than existing arrays `values`
 * will be freed afterwards for consistency. Therefore `values` no longer refers
 * to a valid memory block after calling this function.
 * 
 * @param[in,out] surface   An `mc_ncdf_surface`
 * @param[in]     model     An albedo method as is defined in
 *                          `mystic_definitions.h` TODO redundant?
 * @param[in]     parameter The parameter id; possible values for
 *                          `mc_surface_parameter` are defined in
 *                          `mystic_definitions.h`
 * @param[in,out] values    C Array of values
 * @param[in]     stride    Stride values describing the dimensions in `values`
 * @return                  NCDF_NOERR (i.e. 0) on success, error code on failure
 */
ncdf_status_t add_parameter (mc_ncdf_surface*      surface,
                             mc_surface_method     model,
                             mc_surface_parameter  parameter,
                             float*                values,
                             mc_ncdf_array_stride* stride);

/**
 * @brief Loads a variable from a NetCDF file
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * 
 * @param[in]  ncid           NetCDF file or group id
 * @param[in]  varid          NetCDF variable id
 * @param[out] data           The pointer where all the values of this variable
 *                            should be loaded to
 * @param[in]  target_type    The type of `data` expressed as NetCDF type
 *                            constant. (e.g. NC_INT)
 * @param[in]  accepted_types The types that will be accepted inside the NetCDF
 *                            file. The array has to be terminated by
 *                            `NC_NAT` (not a type). If `NULL` all types will be
 *                            accepted.
 * @param[in]  var_size       Number of values that should be read from the
 *                            variable.
 * @param[in]  target_unit    The unit in which the values in `data` should be
 *                            expressed in. Possible values are defined in
 *                            `ncdf_loader.c` as `TGT_UNIT_*`. Use
 *                            `TGT_UNIT_ANY` for no unit conversion.
 * @return                    NCDF_NOERR (i.e. 0) on success, error code on failure
 */
ncdf_status_t load_variable (int        ncid,
                             int        varid,
                             void**     data,
                             nc_type    target_dtype,
                             nc_type*   accepted_types,
                             size_t     var_size,
                             tgt_unit_t target_unit);

/**
 * @brief Loads packing information for a variable
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * This function loads packing information for a variable that can be used to
 * unpack its values.
 * 
 * @param[in]  ncid           NetCDF file or group id
 * @param[in]  varid          NetCDF variable id
 * @param[in]  name           Name of the variable
 * @param[out] unpacked_dtype The datatype after unpacking
 * @param[out] offset         Offset used for unpacking
 * @param[in]  scale          Scale used for unpacking
 * @return                    NCDF_NOERR (i.e. 0) on success, error code on failure
 */
ncdf_status_t
load_packing_information (int ncid, int varid, const char* name, nc_type* unpacked_dtype, double* offset, double* scale);

/**
 * @brief Indicates the type of dimension
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * Each dimension has to have a corresponding variable that contains the values
 * along its axis and describes the dimension via its attributes. This function
 * tries to find out the type or purpose of a dimension by looking at the
 * corresponding variable.
 * 
 * @param[in]  ncid       NetCDF file or group id
 * @param[in]  varid      NetCDF variable id of dimension values
 * @param[out] dim_use    How this dimension should be used. Possible values are
 *                        defined in `ncdf_loader.c` as `DIM_*`.
 * @param[out] projection Indicated projection. Possible values are
 *                        `MC_PROJ_XY`, `MC_PROJ_LON_LAT` or `MC_PROJ_UNDEFINED`
 * @param[in]  dim_name   Name of the dimension TODO redundant?
 * @return                NCDF_NOERR (i.e. 0) on success, error code on failure
 */
ncdf_status_t get_dimension_purpose (int ncid, int varid, dim_t* dim_use, ncdf_projection_t* projection, const char* dim_name);

/**
 * @brief Calculates the unit scale factor.
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * This function calculates the unit scale factor for converting values of a
 * variable to a certain unit and returns an error if there is no possible
 * conversion. More specifically `unit_scale` describes the factor with which
 * the raw values have to be multiplied with.
 * 
 * @param[in]  ncid        NetCDF file or group id
 * @param[in]  varid       NetCDF variable id
 * @param[in]  target_unit The unit in which the resulting values should be
 *                         expressed in. Possible values are defined in
 *                         `ncdf_loader.c` as `TGT_UNIT_*`. Use `TGT_UNIT_ANY`
 *                         for no unit conversion.
 * @param[out] unit_scale  Scale factor
 * @param[in]  name        Name of the variable TODO redundant?
 * @return                 NCDF_NOERR (i.e. 0) on success, error code on failure
 */
ncdf_status_t get_unit_scale (int ncid, int varid, tgt_unit_t target_unit, double* unit_scale, const char* name);

/**
 * @brief Load dimensions associated with a variable
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * This function retrieves metadata and values from the dimensions that belong
 * to a certain variable. If the dimensions were already previously loaded this
 * function checks if they are identical. It also initializes a stride struct
 * that allows to access values of the variable in a well defined manner and
 * calculates the length of the array. All arrays need to be allocated with
 * `NUM_DIMENSIONS` members and they are accessed using the dimension constants
 * `DIM_*` as they are defined in `ncdf_loader.c`. The arguments
 * `dimension_values` and `dimension_lengths` take an array of pointers that
 * will be initialized inside `load_dimensions`. This allows to directly
 * initialize members of a surface struct. E.g.:
 * double** dimension_values[NUM_DIMENSIONS] = {&s.x, &s.y, &s.spectrum};
 * size_t* dimension_lengths[NUM_DIMENSIONS] = {&s.nx, &s.ny, &s.nspectrum};
 * 
 * @param[in]     global_ncid       NetCDF file id
 * @param[in]     ncid              NetCDF file or group id
 * @param[in]     varid             NetCDF variable id
 * @param[in,out] dimids            Array of NetCDF dimension ids; the array
 *                                  must be allocated and have space for
 *                                  `NUM_DIMENSIONS` ids.
 * @param[in,out] dimension_values  Values along the axis of the dimension; the
 *                                  array must be initialized with
 *                                  `NUM_DIMENSIONS` pointers.
 * @param[in,out] dimension_lengths Lengths of all dimensions; the array must be
 *                                  initialized with `NUM_DIMENSIONS` pointers.
 * @param[in,out] projection        Returns the projection indicated by
 *                                  spatial dimensions.
 * @param[out]    stride            Stride struct that can be used to access an
 *                                  array
 * @return                          NCDF_NOERR (i.e. 0) on success, error code
 *                                  on failure
 */
ncdf_status_t load_dimensions (int                   global_ncid,
                               int                   ncid,
                               int                   varid,
                               int*                  dimids,
                               double***             dimension_values,
                               size_t**              dimension_lengths,
                               ncdf_projection_t*    projection,
                               mc_ncdf_array_stride* stride,
                               const char*           name);

/**
 * @brief Returns the default unit for each dimension
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * 
 * @param[in] dim Dimension constant; possible values are of the form `DIM_*`
 *                and are defined in `ncdf_loader.c`.
 * @return        Unit constant; possible values are of the form `TGT_UNIT_*`
 *                and are defined in `ncdf_loader.c`.
 */
tgt_unit_t tgt_unit_from_dim (int dim);

/**
 * @brief Converts the model array whose values are defined in the NetCDF file
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * This function converts the values of a model array which are defined by flag
 * description inside the NetCDF file to an array whose values are defined by
 * albedo method constants from `mystic_definitions.h` (i.e. `MCALB_*`).
 * 
 * @param[in,out] model_map     Integer array of values that were read from a
 *                              NetCDF file.
 * @param[in]     map_length    Length of `model_map`.
 * @param[in]     flag_meanings Array of strings that describe models.
 *                              (See `mc_model_str`)
 * @param[in]     flag_values   Array of integers that correspond to
 *                              `flag_meanings` and map values in `model_map` to
 *                              albedo methods.
 * @return                      NCDF_NOERR (i.e. 0) on success, error code
 *                              on failure
 */
ncdf_status_t
model_convert (mc_surface_method* model_map, size_t map_length, char** flag_meanings, long* flag_values, size_t num_flags);

/**
 * @brief Loads flag descriptions from a variable with flags
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * This function loads meta data of a variable that contains flags. This
 * metadata can be used in `model_convert` to convert the values in a variable
 * to albedo methods.
 * 
 * @param[in]  ncid          NetCDF file or group id
 * @param[in]  varid         NetCDF variable id
 * @param[in]  name          Name of the variable TODO redundant?
 * @param[out] flag_meanings Uninitialized pointer that will contain an array of
 *                           strings describing the flags (See `mc_model_str`).
 *                           Note: This function is very lazy and only replaces
 *                           spaces with NULL characters. Therefore the string
 *                           is allocated as a whole memory block. So it is very
 *                           important to *only* free the first string in the
 *                           array. Freeing other strings will cause a double
 *                           free.
 * @param[out] flag_values   Uninitialized pointer that will contain an array of
 *                           integers describing the flag ids used in the
 *                           variable. Each entry corresponds to an entry in
 *                           `flag_meanings`.
 * @param[out] num_flags     Number of flags.
 * @return                   NCDF_NOERR (i.e. 0) on success, error code on
 *                           failure.
 */
ncdf_status_t load_flags (int ncid, int varid, const char* name, char*** flag_meanings, long** flag_values, size_t* num_flags);

/**
 * @brief Initializes an interpolated array
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * Interpolates `src` to `interpolated` according to the wavelength grids
 *`spectrum` and `int_spectrum` defined in `mc_ncdf_surface`
 * 
 * @param[in,out] surface      The surface that contains both wavelength grids
 * @param[in]     src          An array that is used as source. It is expected
 *                             to be either `surface->spectral_parameters` or
 *                             `surface->spectral_spatial_parameters`
 * @param[out]    interpolated An array that is used as destination. It is
 *                             expected to be either
 *                             `surface->int_spectral_parameters` or
 *                             `surface->int_spectral_spatial_parameters`
 * @return                     NCDF_NOERR (i.e. 0) on success, error code on
 *                             failure.
 */
ncdf_status_t interpolate_wl (mc_ncdf_surface* surface, mc_ncdf_array* src, mc_ncdf_array* interpolated);

/**
 * @brief Sets the required parameters for each model
 * @author Cedric Wehrum
 * @date Jan 27 2024
 * 
 * @param[in]     method     A surface model. See `mystic_definitions.h` for
 *                           possible values.
 * @param[in,out] parameters An array of booleans that define which parameter
 *                           is required for the given model. The array is
 *                           supposed to be initialized and big enough to
 *                           contain `MCALB_NUM_PARAMS` entries.
 */
void set_requested_parameters_for_model_default (mc_surface_method method, bool* parameters);

/**
 * @brief Sets the required parameters for each model
 * @author Cedric Wehrum
 * @date Jan 27 2024
 * This function allows to define additional parameters that may be necessary
 * for certain models.
 * 
 * @param[in]     method     A surface model. See `mystic_definitions.h` for
 *                           possible values.
 * @param[in,out] parameters An array of booleans that define which parameter
 *                           is required for the given model. The array is
 *                           supposed to be initialized and big enough to
 *                           contain `MCALB_NUM_PARAMS` entries.
 */
void set_requested_parameters_for_model (mc_surface_method method, bool* parameters);

const char* mc_model_str (mc_surface_method model_id);
const char* lrt_nc_type_str (nc_type type);
const char* mc_model_param_str (mc_surface_parameter param_id);

void calc_stride (size_t* stride, mc_ncdf_array_stride* strides, size_t length) {
  strides->param *= length;
  strides->wl *= length;
  strides->x *= length;
  strides->y *= length;
  strides->length *= length;
  *stride = 1;
}

void reset_stride (mc_ncdf_array_stride* stride) {
  stride->x      = 0;
  stride->y      = 0;
  stride->wl     = 0;
  stride->param  = 0;
  stride->length = 1;
}

const char* mc_model_str (mc_surface_method model_id) {
  switch (model_id) {
  case MCALB_LITVINOV:
    return "litvinov";
  case MCALB_TSANG:
    return "tsang";
  case MCALB_MAIGN_BREON:
    return "maignan";
  case MCALB_LAM:
    return "lambertian";
  case MCALB_RPV:
    return "rpv";
  default:
    return "<unknown model>";
  }
}

const char* mc_model_param_str (mc_surface_parameter param_id) {
  switch (param_id) {
  case MCALB_LAM_PARAM_ALBEDO:
    return "lambertian/albedo";
  case MCALB_LITVINOV_PARAM_ALBEDO:
    return "litvinov/albedo";
  case MCALB_LITVINOV_PARAM_RMS_SLOPE:
    return "litvinov/rms_slope";
  case MCALB_LITVINOV_PARAM_FRESNEL:
    return "litvinov/fresnel";
  case MCALB_LITVINOV_PARAM_SHADOWING:
    return "litvinov/shadowing";
  case MCALB_LITVINOV_PARAM_REFR_INDEX:
    return "litvinov/refr_index";
  case MCALB_LITVINOV_PARAM_REFR_INDEX_I:
    return "litvinov/refr_index_i";
  case MCALB_LITVINOV_PARAM_RMS_SLOPE2:
    return "litvinov/rms_slope2";
  case MCALB_MAIGN_BREON_PARAM_C_MAIGN:
    return "maignan/c_maign";
  case MCALB_MAIGN_BREON_PARAM_ARVI:
    return "maignan/arvi";
  case MCALB_MAIGN_BREON_PARAM_REFR_INDEX:
    return "maignan/refr_index";
  case MCALB_MAIGN_BREON_PARAM_REFR_INDEX_I:
    return "maignan/refr_index_i";
  case MCALB_RPV_PARAM_RHO0:
    return "rpv/rho0";
  case MCALB_RPV_PARAM_K:
    return "rpv/k";
  case MCALB_RPV_PARAM_THETA:
    return "rpv/theta";
  case MCALB_RPV_PARAM_SCALE:
    return "rpv/scale";
  case MCALB_RPV_PARAM_SIGMA:
    return "rpv/sigma";
  case MCALB_RPV_PARAM_T1:
    return "rpv/t1";
  case MCALB_RPV_PARAM_T2:
    return "rpv/t2";
  case MCALB_TSANG_PARAM_WIND:
    return "tsang/wind";
  default:
    return "<unknown parameter>";
  }
}

const char* lrt_nc_type_str (nc_type type) {
  switch (type) {
  case NC_UBYTE:
    return "uint8";
  case NC_BYTE:
    return "int8";
  case NC_USHORT:
    return "uint16";
  case NC_SHORT:
    return "int16";
  case NC_UINT:
    return "uint32";
  case NC_INT:
    return "int32";
  case NC_UINT64:
    return "uint64";
  case NC_INT64:
    return "int64";
  case NC_FLOAT:
    return "float";
  case NC_DOUBLE:
    return "double";
  case NC_STRING:
    return "string";
  case NC_CHAR:
    return "char";
  default:
    return "<unknown type>";
  }
}

void set_requested_parameters_for_model_default (mc_surface_method method, bool* parameters) {
  const char* model_name = mc_model_str (method);
  char        model_prefix[strlen (model_name) + 2];
  memset (&model_prefix[0], 0, sizeof (model_prefix));
  strcat (strcat (&model_prefix[0], model_name), "/");
  for (mc_surface_parameter i = 0; i < MCALB_NUM_PARAMS; i++) {
    if (strncmp (mc_model_param_str (i), &model_prefix[0], strlen (&model_prefix[0])) == 0) {
      parameters[i] = true;
    }
  }
}

void set_requested_parameters_for_model (mc_surface_method method, bool* parameters) {
  switch (method) {
  case MCALB_MAIGN_BREON:
    // Maignan needs RPV parameters
    set_requested_parameters_for_model_default (MCALB_MAIGN_BREON, parameters);
    set_requested_parameters_for_model_default (MCALB_RPV, parameters);
    return;
  default:
    set_requested_parameters_for_model_default (method, parameters);
    return;
  }
}

#define NCDF_INTEGERS NC_UBYTE, NC_BYTE, NC_USHORT, NC_SHORT, NC_UINT, NC_INT, NC_UINT64, NC_INT64
#define NCDF_FLOATS NC_FLOAT, NC_DOUBLE

bool is_integer (nc_type type) {
  return type == NC_UBYTE || type == NC_BYTE || type == NC_USHORT || type == NC_SHORT || type == NC_UINT || type == NC_INT ||
         type == NC_UINT64 || type == NC_INT64;
}

bool is_float (nc_type type) {
  return type == NC_FLOAT || type == NC_DOUBLE;
}

bool is_string (nc_type type) {
  return type == NC_CHAR || type == NC_STRING;
}

mc_surface_method find_surface_model_id (const char* model_name) {
  if (strncmp ("model_", model_name, 6) == 0) {
    model_name = &model_name[6];
  }
  for (mc_surface_method i = 0; i < MCALB_NUM_METHODS; i++) {
    if (strcmp (mc_model_str (i), model_name) == 0) {
      return i;
    }
  }
  return MCALB_NONE;
}

mc_surface_parameter find_model_param_id (const char* model_name, const char* param_name) {
  char name[strlen (model_name) + strlen (param_name) + 2];
  sprintf (&name[0], "%s/%s", model_name, param_name);
  for (mc_surface_parameter i = 0; i < MCALB_NUM_PARAMS; i++) {
    if (strcmp (mc_model_param_str (i), name) == 0) {
      return i;
    }
  }
  return MCALB_INVALID_PARAM;
}

void free_array (mc_ncdf_array* array) {
  free (array->data);
  array->data = NULL;
}

void free_raw_values (mc_ncdf_surface* data) {
  free (data->spectrum);
  data->spectrum = NULL;

  free_array (&data->spectral_parameters);
  free_array (&data->spectral_spatial_parameters);
}

void mc_ncdf_free (mc_ncdf_surface* data) {
  free (data->x);
  free (data->y);
  free (data->spectrum);
  free (data->int_spectrum);
  free (data->models);

  data->x            = NULL;
  data->y            = NULL;
  data->spectrum     = NULL;
  data->int_spectrum = NULL;
  data->models       = NULL;

  free_array (&data->parameters);
  free_array (&data->spatial_parameters);
  free_array (&data->spectral_parameters);
  free_array (&data->spectral_spatial_parameters);
  free_array (&data->int_spectral_parameters);
  free_array (&data->int_spectral_spatial_parameters);
}

ncdf_status_t add_parameter (mc_ncdf_surface*      surface,
                             mc_surface_method     model,
                             mc_surface_parameter  parameter,
                             float*                values,
                             mc_ncdf_array_stride* stride) {
  mc_ncdf_array* array;
  mc_ncdf_array* int_array;
  size_t         nx, ny, nspectrum;

  if (stride->wl) {
    nspectrum = surface->nspectrum;
    if (stride->x || stride->y) {
      array     = &surface->spectral_spatial_parameters;
      int_array = &surface->int_spectral_spatial_parameters;
    } else {
      array     = &surface->spectral_parameters;
      int_array = &surface->int_spectral_parameters;
    }
  } else {
    nspectrum = 1;
    if (stride->x || stride->y) {
      array = &surface->spatial_parameters;
    } else {
      array = &surface->parameters;
    }
    int_array = array;
  }
  if (stride->x || stride->y) {
    nx = surface->nx;
    ny = surface->ny;
  } else {
    nx    = 1;
    ny    = 1;
    model = 0;
  }

  ASSERTOUTV (surface->param_link[parameter] == NULL,
              NCDF_PARAMETER_DOUBLE_INIT,
              "Parameter '%s' already initialized!",
              mc_model_param_str (parameter));
  surface->param_link[parameter]     = array;
  surface->int_param_link[parameter] = int_array;
  array->param_indices[parameter]    = array->param_count[model];
  array->param_count[model]++;
  if (array->data == NULL) {
    array->stride          = *stride;
    array->stride.param    = array->stride.length;
    array->data            = values;
    array->max_param_count = array->param_count[model];
  } else {
    if (array->param_count[model] > array->max_param_count) {
      array->stride.length += array->stride.param;
      array->data            = realloc (array->data, array->stride.length * sizeof (array->data[0]));
      array->max_param_count = array->param_count[model];
    }

    for (size_t x = 0; x < nx; x++) {
      for (size_t y = 0; y < ny; y++) {
        if ((nx == 1 && ny == 1) || mc_ncdf_get_model (surface, x, y) == model) {
          for (size_t wl = 0; wl < nspectrum; wl++) {
            set_param (array, x, y, wl, array->param_indices[parameter], get_param_by_stride (stride, values, x, y, wl, 0));
          }
        }
      }
    }
    free (values);
  }
  return NCDF_NOERR;
}

//TODO test without conversion
//TODO round to nearest integer
#define LOAD_VARIABLE_TPL(ncid,                                                                                                    \
                          varid,                                                                                                   \
                          data,                                                                                                    \
                          var_size,                                                                                                \
                          target_dtype,                                                                                            \
                          var_type,                                                                                                \
                          type,                                                                                                    \
                          ncdf_load_function,                                                                                      \
                          name,                                                                                                    \
                          offset,                                                                                                  \
                          scale,                                                                                                   \
                          unit_scale)                                                                                              \
  {                                                                                                                                \
    type*         var;                                                                                                             \
    ncdf_status_t status;                                                                                                          \
    if (target_dtype == NC_INT64) {                                                                                                \
      *data = calloc (var_size, sizeof (int));                                                                                     \
    } else if (target_dtype == NC_FLOAT) {                                                                                         \
      *data = calloc (var_size, sizeof (float));                                                                                   \
    } else if (target_dtype == NC_DOUBLE) {                                                                                        \
      *data = calloc (var_size, sizeof (double));                                                                                  \
    } else if (target_dtype == NC_UBYTE) {                                                                                         \
      *data = calloc (var_size, sizeof (unsigned char));                                                                           \
    } else {                                                                                                                       \
      CHKERROUTV (NCDF_UNSUPPORTED_MAP_TYPE, "Unsupported type (%d) for map %s!", target_dtype, name);                             \
    }                                                                                                                              \
    if (var_type == target_dtype) {                                                                                                \
      var = (type*) *data;                                                                                                         \
    } else {                                                                                                                       \
      var = calloc (var_size, sizeof (type));                                                                                      \
    }                                                                                                                              \
    status = ncdf_load_function (ncid, varid, var);                                                                                \
    NC_CHKERROUTV (status, "Cannot retrieve integer map %s!", name);                                                               \
    if (var_type != target_dtype) {                                                                                                \
      for (size_t i = 0; i < var_size; i++) {                                                                                      \
        if (target_dtype == NC_INT64) {                                                                                            \
          ((long*) *data)[i] = (long) ((var[i] * scale + offset) * unit_scale);                                                    \
        } else if (target_dtype == NC_INT) {                                                                                       \
          ((int*) *data)[i] = (int) ((var[i] * scale + offset) * unit_scale);                                                      \
        } else if (target_dtype == NC_FLOAT) {                                                                                     \
          ((float*) *data)[i] = (float) ((var[i] * scale + offset) * unit_scale);                                                  \
        } else if (target_dtype == NC_DOUBLE) {                                                                                    \
          ((double*) *data)[i] = (double) ((var[i] * scale + offset) * unit_scale);                                                \
        } else if (target_dtype == NC_UBYTE) {                                                                                     \
          ((unsigned char*) *data)[i] = (unsigned char) ((var[i] * scale + offset) * unit_scale);                                  \
        } else {                                                                                                                   \
          CHKERROUTV (NCDF_UNSUPPORTED_MAP_TYPE, "Unsupported type (%d) for map %s!", target_dtype, name);                         \
        }                                                                                                                          \
      }                                                                                                                            \
      free (var);                                                                                                                  \
    }                                                                                                                              \
  }

ncdf_status_t load_variable (int        ncid,
                             int        varid,
                             void**     data,
                             nc_type    target_dtype,
                             nc_type*   accepted_types,
                             size_t     var_size,
                             tgt_unit_t target_unit) {
  nc_type       var_type;
  char          var_name[NC_MAX_NAME + 1];
  nc_type       unpacked_dtype;
  double        offset, scale, unit_scale;
  ncdf_status_t status;

  status = nc_inq_var (ncid, varid, &var_name[0], &var_type, NULL, NULL, NULL);
  NC_CHKERROUTV (status, "Cannot inquiry variable %d!", varid);
  status = load_packing_information (ncid, varid, &var_name[0], &unpacked_dtype, &offset, &scale);
  CHKERROUTV (status, "Cannot load packing information of variable %s!", &var_name[0]);

  if (unpacked_dtype == NC_NAT) {
    unpacked_dtype = var_type;
  }

  if (accepted_types) {
    size_t i;
    for (i = 0; accepted_types[i] != NC_NAT; i++) {
      if (accepted_types[i] == unpacked_dtype) {
        break;
      }
    }
    if (accepted_types[i] == NC_NAT) {
      fprintf (stderr,
               "Variable '%s' has unsupported type %s (%d)!\n",
               &var_name[0],
               lrt_nc_type_str (unpacked_dtype),
               unpacked_dtype);
      fprintf (stderr, "One of the following types were expected: ");
      for (i = 0; accepted_types[i] != NC_NAT; i++) {
        if (i > 0) {
          fprintf (stderr, ", ");
        }
        fprintf (stderr, "%s", lrt_nc_type_str (accepted_types[i]));
      }
      fprintf (stderr, "\n");
      CHKERR (NCDF_TYPE_NOT_ACCEPTABLE);
    }
  }

  status = get_unit_scale (ncid, varid, target_unit, &unit_scale, &var_name[0]);
  CHKERROUTV (status, "Cannot find out units of variable '%s' or units are incompatible", &var_name[0]);

  //TODO add other types
  if (var_type == NC_USHORT) {
    LOAD_VARIABLE_TPL (ncid,
                       varid,
                       data,
                       var_size,
                       target_dtype,
                       var_type,
                       unsigned short,
                       nc_get_var_ushort,
                       &var_name[0],
                       offset,
                       scale,
                       unit_scale);
  } else if (var_type == NC_INT64) {
    LOAD_VARIABLE_TPL (ncid,
                       varid,
                       data,
                       var_size,
                       target_dtype,
                       var_type,
                       int,
                       nc_get_var_int,
                       &var_name[0],
                       offset,
                       scale,
                       unit_scale);
  } else if (var_type == NC_FLOAT) {
    LOAD_VARIABLE_TPL (ncid,
                       varid,
                       data,
                       var_size,
                       target_dtype,
                       var_type,
                       float,
                       nc_get_var_float,
                       &var_name[0],
                       offset,
                       scale,
                       unit_scale);
  } else if (var_type == NC_DOUBLE) {
    LOAD_VARIABLE_TPL (ncid,
                       varid,
                       data,
                       var_size,
                       target_dtype,
                       var_type,
                       double,
                       nc_get_var_double,
                       &var_name[0],
                       offset,
                       scale,
                       unit_scale);
  } else if (var_type == NC_UBYTE) {
    LOAD_VARIABLE_TPL (ncid,
                       varid,
                       data,
                       var_size,
                       target_dtype,
                       var_type,
                       unsigned char,
                       nc_get_var_ubyte,
                       &var_name[0],
                       offset,
                       scale,
                       unit_scale);
  } else {
    CHKERROUTV (NCDF_UNSUPPORTED_VAR_TYPE,
                "Unsupported type '%s' (%d) for map %s!",
                lrt_nc_type_str (var_type),
                var_type,
                &var_name[0]);
  }
  return NCDF_NOERR;
}

ncdf_status_t get_dimension_purpose (int ncid, int varid, dim_t* dim_use, ncdf_projection_t* projection, const char* dim_name) {
  size_t        standard_name_chars, axis_chars;
  char*         standard_name = NULL;
  char*         axis          = NULL;
  nc_type       type          = NC_NAT;
  ncdf_status_t status;

  *dim_use = DIM_UNKNOWN;

  if (nc_inq_atttype (ncid, varid, "standard_name", &type) == NC_NOERR) {
    ASSERTOUTV (is_string (type),
                NCDF_ERROR_STANDARD_NAME,
                "Attribute 'standard_name' associated with the variable '%s' is not a string. Instead got '%s' (%d)",
                dim_name,
                lrt_nc_type_str (type),
                type);
    status = nc_inq_attlen (ncid, varid, "standard_name", &standard_name_chars);
    NC_CHKERROUTV (status,
                   "Cannot inquiry attribute length of the attribute 'standard_name' associated with the variable '%s'!",
                   dim_name);
    standard_name = calloc (standard_name_chars + 1, sizeof (char));
    status        = nc_get_att_text (ncid, varid, "standard_name", standard_name);
    NC_CHKERROUTV (status, "Cannot retrive attribute 'standard_name' for variable %s!", dim_name);
  }

  if (nc_inq_atttype (ncid, varid, "axis", &type) == NC_NOERR) {
    ASSERTOUTV (is_string (type),
                NCDF_ERROR_AXIS,
                "Attribute 'axis' associated with the variable '%s' is not a string. Instead got '%s' (%d)",
                dim_name,
                lrt_nc_type_str (type),
                type);
    status = nc_inq_attlen (ncid, varid, "axis", &axis_chars);
    NC_CHKERROUTV (status, "Cannot inquiry attribute length of the attribute 'axis' associated with the variable '%s'!", dim_name);
    axis   = calloc (axis_chars + 1, sizeof (char));
    status = nc_get_att_text (ncid, varid, "axis", axis);
    NC_CHKERROUTV (status, "Cannot retrive attribute 'axis' for variable %s!", dim_name);
  }

  *projection = MC_PROJ_UNDEFINED;

  if (standard_name && strcasecmp ("grid_longitude", standard_name) == 0) {
    *dim_use    = DIM_X;
    *projection = MC_PROJ_LON_LAT;
  } else if (standard_name && strcasecmp ("grid_latitude", standard_name) == 0) {
    *dim_use    = DIM_Y;
    *projection = MC_PROJ_LON_LAT;
  } else if (standard_name && strcasecmp ("radiation_wavelength", standard_name) == 0) {
    *dim_use = DIM_SPEC;
  } else if (strcasecmp ("lon", dim_name) == 0) {
    *dim_use    = DIM_X;
    *projection = MC_PROJ_LON_LAT;
  } else if (strcasecmp ("lat", dim_name) == 0) {
    *dim_use    = DIM_Y;
    *projection = MC_PROJ_LON_LAT;
  } else if (strcasecmp ("spectral", dim_name) == 0) {
    *dim_use = DIM_SPEC;
  } else if (axis && strcasecmp ("X", axis) == 0) {
    *dim_use    = DIM_X;
    *projection = MC_PROJ_XY;
  } else if (axis && strcasecmp ("Y", axis) == 0) {
    *dim_use    = DIM_Y;
    *projection = MC_PROJ_XY;
  } else {
    CHKERROUTV (NCDF_INVALID_DIMENSION,
                "Dimension %s is of unknown type. Make sure to use a standardized name or the correct attributes.",
                dim_name);
  }
  return NCDF_NOERR;
}

tgt_unit_t tgt_unit_from_dim (dim_t dim) {
  switch (dim) {
  case DIM_X:
    return TGT_UNIT_DEG_EAST;
  case DIM_Y:
    return TGT_UNIT_DEG_NORTH;
  case DIM_SPEC:
    return TGT_UNIT_NM;
  default:
    WARNV ("Unknown dimension %d in 'tgt_unit_from_dim'", dim);
    return TGT_UNIT_ANY;
  }
}

ncdf_status_t get_unit_scale (int ncid, int varid, tgt_unit_t target_unit, double* unit_scale, const char* name) {
  size_t  units_chars;
  char*   units = NULL;
  nc_type type;
  *unit_scale = 1.;
  ncdf_status_t status;

  if (nc_inq_atttype (ncid, varid, "units", &type) == NC_NOERR) {
    ASSERTOUTV (is_string (type),
                NCDF_UNSUPPORTED_DIM_TYPE,
                "Attribute 'units' associated with the variable '%s' is not a string. Instead got '%s' (%d)",
                name,
                lrt_nc_type_str (type),
                type);
    status = nc_inq_attlen (ncid, varid, "units", &units_chars);
    NC_CHKERROUTV (status, "Cannot inquiry attribute length of the attribute 'units' associated with the variable '%s'!", name);
    units  = calloc (units_chars + 1, sizeof (char));
    status = nc_get_att_text (ncid, varid, "units", units);
    NC_CHKERROUTV (status, "Cannot retrive attribute 'units' for variable %s!", name);
  }

  //TODO quick & dirty, use UDUNITS library instead
  switch (target_unit) {
  case TGT_UNIT_NM:
    if (units && strcasecmp ("nm", units) == 0) {
      *unit_scale = 1.;
    } else if (units && strcasecmp ("µm", units) == 0) {
      *unit_scale = 1e3;
    } else if (units && strcasecmp ("mm", units) == 0) {
      *unit_scale = 1e6;
    } else if (units && strcasecmp ("cm", units) == 0) {
      *unit_scale = 1e7;
    } else if (units && strcasecmp ("m", units) == 0) {
      *unit_scale = 1e9;
    } else if (units) {
      CHKERROUTV (NCDF_UNKNOWN_UNIT, "Unknown unit '%s' for variable '%s'", units, name);
    } else {
      *unit_scale = 1.;  // Assuming 'nm'
    }
    break;
  case TGT_UNIT_DEG_EAST:
    if (units && strcasecmp ("degrees_east", units) != 0 && strcasecmp ("degree_east", units) != 0 &&
        strcasecmp ("degree_E", units) != 0 && strcasecmp ("degrees_E", units) != 0 && strcasecmp ("degreeE", units) != 0 &&
        strcasecmp ("degreesE", units) != 0 && strcasecmp ("m", units) != 0) {
      CHKERROUTV (NCDF_UNKNOWN_UNIT, "Unknown unit '%s' for variable '%s'", units, name);
    }
    break;
  case TGT_UNIT_DEG_NORTH:
    if (units && strcasecmp ("degrees_north", units) != 0 && strcasecmp ("degree_north", units) != 0 &&
        strcasecmp ("degree_N", units) != 0 && strcasecmp ("degrees_N", units) != 0 && strcasecmp ("degreeN", units) != 0 &&
        strcasecmp ("degreesN", units) != 0 && strcasecmp ("m", units) != 0) {
      CHKERROUTV (NCDF_UNKNOWN_UNIT, "Unknown unit '%s' for variable '%s'", units, name);
    }
    break;
  case TGT_UNIT_ANY:
    *unit_scale = 1.;
    break;
  default:
    CHKERROUTV (NCDF_UNKNOWN_UNIT, "Unknown target unit %d", target_unit);
  }

  free (units);
  return NCDF_NOERR;
}

ncdf_status_t load_dimensions (int                   global_ncid,
                               int                   ncid,
                               int                   varid,
                               int*                  dimids,
                               double***             dimension_values,
                               size_t**              dimension_lengths,
                               ncdf_projection_t*    spatial_projection,
                               mc_ncdf_array_stride* stride,
                               const char*           name) {
  int           ndims;
  int*          dim_ids;
  int           accepted_types[] = {NC_DOUBLE, NC_FLOAT, NC_NAT};
  ncdf_status_t status;

  status = nc_inq_var (ncid, varid, NULL, NULL, &ndims, NULL, NULL);
  NC_CHKERROUTV (status, "Cannot inquiry number of dimensions of variable %d!", varid);
  dim_ids = calloc (ndims, sizeof (int));
  status  = nc_inq_var (ncid, varid, NULL, NULL, NULL, dim_ids, NULL);
  NC_CHKERROUTV (status, "Cannot inquiry dimension ids of variable '%s' (%d)!", name, varid);
  reset_stride (stride);

  for (size_t i = 0; i < ndims; i++) {
    char              dim_name[NC_MAX_NAME + 1];
    size_t            dim_len;
    int               dim_varid;
    dim_t             dim_use;
    ncdf_projection_t projection;

    status = nc_inq_dim (ncid, dim_ids[i], &dim_name[0], &dim_len);
    NC_CHKERROUTV (status, "Cannot inquiry dimension id %d (%lu)!", dim_ids[i], i);
    status = nc_inq_varid (global_ncid, &dim_name[0], &dim_varid);
    NC_CHKERROUTV (status, "Cannot find values for dimension %s (%d)!", &dim_name[0], dim_ids[i]);
    status = get_dimension_purpose (global_ncid, dim_varid, &dim_use, &projection, &dim_name[0]);
    CHKERROUTV (status, "Cannot find out intended purpose of dimension '%s'", &dim_name[0]);

    if (projection != MC_PROJ_UNDEFINED) {
      if (*spatial_projection == MC_PROJ_UNDEFINED) {
        *spatial_projection = projection;
      } else if (*spatial_projection != projection) {
        CHKERROUTV (NCDF_INVALID_PROJECTION,
                    "Indicated projection of dimension '%s' (%d) conflicts earlier projection (%d)",
                    &dim_name[0],
                    projection,
                    *spatial_projection);
      }
    }

    switch (dim_use) {
    case DIM_X:
      calc_stride (&stride->x, stride, dim_len);
      break;
    case DIM_Y:
      calc_stride (&stride->y, stride, dim_len);
      break;
    case DIM_SPEC:
      calc_stride (&stride->wl, stride, dim_len);
      break;
    default:
      CHKERROUTV (NCDF_INVALID_DIMENSION, "Unknown dimension purpose %d!", dim_use);
    }

    if (*dimension_values[dim_use] == NULL) {
      dimids[dim_use]               = dim_ids[i];
      *(dimension_lengths[dim_use]) = dim_len;

      status = load_variable (global_ncid,
                              dim_varid,
                              (void**) (dimension_values[dim_use]),
                              NC_DOUBLE,
                              &accepted_types[0],
                              dim_len,
                              tgt_unit_from_dim (dim_use));
      CHKERROUTV (status, "Cannot load dimension values of dimension %s", &dim_name[0]);
    } else if (dimids[dim_use] == -1) {
      double* alt_dimension_values;

      ASSERTOUTV (*(dimension_lengths[dim_use]) == dim_len,
                  NCDF_INVALID_DIMENSION,
                  "Dimension '%s' was previously defined in a different file, but has a different length! (%lu != %lu)\n"
                  "Note that the longitudinal, latitudinal and spectral dimension must be identical between all variables",
                  &dim_name[0],
                  *(dimension_lengths[dim_use]),
                  dim_len);

      status = load_variable (global_ncid,
                              dim_varid,
                              (void**) (&alt_dimension_values),
                              NC_DOUBLE,
                              &accepted_types[0],
                              dim_len,
                              tgt_unit_from_dim (dim_use));
      CHKERROUTV (status, "Cannot load dimension values of dimension %s", &dim_name[0]);
      for (size_t j = 0; j < dim_len; j++) {
        ASSERTOUTV (
          (*dimension_values[dim_use])[j] == alt_dimension_values[j],
          NCDF_INVALID_DIMENSION,
          "Dimension '%s' was previously defined in a different file, but its values differ! ((previous[%lu]=%f) != (%s[%lu]=%f))\n"
          "Note that the longitudinal, latitudinal and spectral dimension must be identical between all variables",
          &dim_name[0],
          j,
          (*dimension_values[dim_use])[j],
          &dim_name[0],
          j,
          alt_dimension_values[j]);
      }
      dimids[dim_use] = dim_ids[i];
      free (alt_dimension_values);
    } else if (dimids[dim_use] != dim_ids[i]) {
      CHKERROUTV (NCDF_INVALID_DIMENSION,
                  "Dimension '%s' is already defined in the same file! (id %d != id %d)\n"
                  "Note that the longitudinal, latitudinal and spectral dimension must be the same between all variables",
                  &dim_name[0],
                  dimids[dim_use],
                  dim_ids[i]);
    }
  }
  free (dim_ids);
  return NCDF_NOERR;
}

ncdf_status_t
model_convert (mc_surface_method* model_map, size_t map_length, char** flag_meanings, long* flag_values, size_t num_flags) {
  int dictionary[256] = {-1};
  for (size_t i = 0; i < num_flags; i++) {
    ASSERTOUTV (
      flag_values[i] >= 0 && flag_values[i] < 256,
      NCDF_FLAG_ERROR,
      "The flag values for the model map must be in the range of [0-255]! However the flag '%s' was given with the value '%ld'",
      flag_meanings[i],
      flag_values[i]);
    dictionary[flag_values[i]] = find_surface_model_id (flag_meanings[i]);
  }
  for (size_t i = 0; i < map_length; i++) {
    ASSERTOUTV (dictionary[model_map[i]] != -1,
                NCDF_FLAG_ERROR,
                "The flag with the value %d (position: %lu) in the model map was not defined in the flag_values attribute!",
                model_map[i],
                i);
    model_map[i] = dictionary[model_map[i]];
  }
  return NCDF_NOERR;
}

ncdf_status_t load_flags (int ncid, int varid, const char* name, char*** flag_meanings, long** flag_values, size_t* num_flags) {
  int           flag_values_id, flag_meanings_id;
  char*         raw_flag_meanings;
  nc_type       values_type, meanings_type;
  size_t        num_chars;
  ncdf_status_t status;

  status = nc_inq_attid (ncid, varid, "flag_values", &flag_values_id);
  NC_CHKERROUTV (status, "Cannot find the 'flag_values' attribute associated with the variable '%s'", name);
  status = nc_inq_attid (ncid, varid, "flag_meanings", &flag_meanings_id);
  NC_CHKERROUTV (status, "Cannot find the 'flag_meanings' attribute associated with the variable '%s'", name);
  status = nc_inq_atttype (ncid, varid, "flag_values", &values_type);
  NC_CHKERROUTV (status, "Cannot inquiry attribute type of the attribute 'flag_values' associated with the variable '%s'!", name);
  ASSERTOUTV (is_integer (values_type),
              NCDF_FLAG_ERROR,
              "Attribute 'flag_values' associated with the variable '%s' is not an integer. Instead got '%s' (%d)",
              name,
              lrt_nc_type_str (values_type),
              values_type);
  status = nc_inq_atttype (ncid, varid, "flag_meanings", &meanings_type);
  NC_CHKERROUTV (status, "Cannot inquiry attribute type of the attribute 'flag_meanings' associated with the variable '%s'!", name);
  ASSERTOUTV (is_string (meanings_type),
              NCDF_FLAG_ERROR,
              "Attribute 'flag_meanings' associated with the variable '%s' is not a string. Instead got '%s' (%d)",
              name,
              lrt_nc_type_str (meanings_type),
              meanings_type);
  status = nc_inq_attlen (ncid, varid, "flag_values", num_flags);
  NC_CHKERROUTV (status, "Cannot inquiry attribute length of the attribute 'flag_values' associated with the variable '%s'!", name);
  status = nc_inq_attlen (ncid, varid, "flag_meanings", &num_chars);
  NC_CHKERROUTV (status,
                 "Cannot inquiry attribute length of the attribute 'flag_meanings' associated with the variable '%s'!",
                 name);
  *flag_meanings    = calloc (*num_flags, sizeof (char*));
  raw_flag_meanings = calloc (num_chars + 1, sizeof (char));
  *flag_values      = calloc (*num_flags, sizeof (long));
  status            = nc_get_att_long (ncid, varid, "flag_values", *flag_values);
  NC_CHKERROUTV (status, "Cannot get flag_values for variable %s!", name);
  status = nc_get_att_text (ncid, varid, "flag_meanings", raw_flag_meanings);
  NC_CHKERROUTV (status, "Cannot get flag_meanings for variable %s!", name);
  (*flag_meanings)[0] = &raw_flag_meanings[0];
  for (size_t i = 0, j = 1; i < num_chars; i++) {
    if (raw_flag_meanings[i] == ' ') {
      raw_flag_meanings[i]  = 0;
      (*flag_meanings)[j++] = &raw_flag_meanings[i + 1];
    }
  }
  return NCDF_NOERR;
}

ncdf_status_t
load_packing_information (int ncid, int varid, const char* name, nc_type* target_dtype, double* offset, double* scale) {
  int           n_missing_atts = 0;
  nc_type       scale_type, offset_type;
  ncdf_status_t status;

  if ((status = nc_inq_atttype (ncid, varid, "actual_range", target_dtype))) {
    n_missing_atts++;
  }
  if ((status = nc_inq_atttype (ncid, varid, "scale_factor", &scale_type))) {
    n_missing_atts++;
  }
  if ((status = nc_inq_atttype (ncid, varid, "add_offset", &offset_type))) {
    n_missing_atts++;
  }
  if (n_missing_atts == 3) {
    *target_dtype = NC_NAT;
    *offset       = 0.;
    *scale        = 1.;
    return NCDF_NOERR;
  } else if (n_missing_atts != 0) {
    CHKERROUTV (NCDF_PACK_ERROR,
                "Some packing attributes for the variable '%s' are missing or cannot be loaded! Make sure to include "
                "'actual_range', 'scale_factor' and 'add_offset'\n",
                name);
  }
  ASSERTOUTV (is_float (scale_type),
              NCDF_PACK_ERROR,
              "Attribute 'scale_factor' associated with the variable '%s' is not a floating point number. Instead got '%s' (%d)\n",
              name,
              lrt_nc_type_str (scale_type),
              scale_type);
  ASSERTOUTV (is_float (offset_type),
              NCDF_PACK_ERROR,
              "Attribute 'add_offset' associated with the variable '%s' is not a floating point number. Instead got '%s' (%d)\n",
              name,
              lrt_nc_type_str (offset_type),
              offset_type);
  status = nc_get_att_double (ncid, varid, "scale_factor", scale);
  NC_CHKERROUTV (status, "Cannot get scale_factor for variable %s!\n", name);
  status = nc_get_att_double (ncid, varid, "add_offset", offset);
  NC_CHKERROUTV (status, "Cannot get scale_factor for variable %s!\n", name);
  return NCDF_NOERR;
}

int has_trailing (char* str, char trail) {
  return str[strlen (str) - 1] == trail;
}

ncdf_status_t mc_ncdf_load (char** filenames, mc_ncdf_surface* data, int quiet) {
  size_t n, i;
  char** files = NULL;
  for (n = 0, i = 0; filenames[n]; n++) {
    if (!has_trailing (filenames[n], '/')) {
      files      = realloc (files, (i + 1) * sizeof (char*));
      files[i++] = filenames[n];
    }
  }
  return mc_ncdf_load_n (files, i, data, quiet);
}

bool stride_is_spatial (mc_ncdf_array_stride* stride) {
  return stride->x && stride->y && !stride->wl;
}

bool stride_is_spectral (mc_ncdf_array_stride* stride) {
  return !stride->x && !stride->y && stride->wl;
}

bool stride_is_spatial_spectral (mc_ncdf_array_stride* stride) {
  return stride->x && stride->y && stride->wl;
}

bool stride_is_const (mc_ncdf_array_stride* stride) {
  return !stride->x && !stride->y && !stride->wl;
}

ncdf_status_t mc_ncdf_load_n (char** filenames, size_t n_filenames, mc_ncdf_surface* data, int quiet) {
  int           ncid[n_filenames];
  int           nvars;
  int           n_groups;
  int*          groups;
  char*         model_file = NULL;
  int           n_model_groups;
  int           n_parameters;
  int*          model_groups;
  int           dimension_ids[n_filenames][NUM_DIMENSIONS];
  double**      dimension_values[NUM_DIMENSIONS]  = {&data->x, &data->y, &data->spectrum};
  size_t*       dimension_lengths[NUM_DIMENSIONS] = {&data->nx, &data->ny, &data->nspectrum};
  nc_type       accepted_model_types[]            = {NCDF_INTEGERS, NC_NAT};
  nc_type       accepted_parameter_types[]        = {NCDF_FLOATS, NC_NAT};
  ncdf_status_t status;
  bool          requested_parameters[MCALB_NUM_PARAMS] = {false};

  *data = (const mc_ncdf_surface){0};
  ASSERTOUT (n_filenames > 0, NCDF_INTERNAL_ERROR, "No files provided\n");
  for (size_t i = 0; i < n_filenames; i++) {
    for (size_t j = 0; j < NUM_DIMENSIONS; j++) {
      dimension_ids[i][j] = -1;
    }
  }
  for (size_t i = 0; i < n_filenames; i++) {
    INFOV (quiet, "Loading variables from NetCDF file '%s'", filenames[i]);

    status = nc_open (filenames[i], NC_NOWRITE, &ncid[i]);
    NC_CHKERROUT (status, "Cannot open netcdf file!");
    status = nc_inq (ncid[i], NULL, &nvars, NULL, NULL);
    NC_CHKERROUT (status, "Cannot inquiry global variables!");

    for (int varid = 0; varid < nvars; varid++) {
      char var_name[NC_MAX_NAME + 1];

      status = nc_inq_var (ncid[i], varid, &var_name[0], NULL, NULL, NULL, NULL);
      NC_CHKERROUTV (status, "Cannot inquiry variable %d!", varid);

      if (strcmp ("models", &var_name[0]) == 0) {
        char** model_meanings;
        long*  model_values;
        size_t num_models;

        ASSERTOUTV (model_file == NULL,
                    NCDF_MODELS_DOUBLE_DEFINITION,
                    "'models' map is defined multiple times! First defined in '%s', last defined in '%s'.",
                    model_file,
                    filenames[i]);
        model_file = filenames[i];

        INFO (quiet, "Loading model map");

        status = load_dimensions (ncid[i],
                                  ncid[i],
                                  varid,
                                  &dimension_ids[i][0],
                                  &dimension_values[0],
                                  &dimension_lengths[0],
                                  &data->projection,
                                  &data->models_stride,
                                  &var_name[0]);
        CHKERROUTV (status, "Cannot load dimensions of variable '%s'!", &var_name[0]);
        ASSERTOUTV (stride_is_spatial (&data->models_stride) || stride_is_const (&data->models_stride),
                    NCDF_UNSUPPORTED_DIM_TYPE,
                    "Unsupported map dimensions for '%s': Exactly two spatial dimensions or a constant needs to be provided!",
                    &var_name[0]);
        status = load_variable (ncid[i],
                                varid,
                                (void**) &data->models,
                                NC_UBYTE,
                                &accepted_model_types[0],
                                data->models_stride.length,
                                TGT_UNIT_ANY);
        CHKERROUTV (status, "Cannot load contents of the variable '%s'", &var_name[0]);
        status = load_flags (ncid[i], varid, &var_name[0], &model_meanings, &model_values, &num_models);
        CHKERROUT (status, "Cannot load flag meanings for the models variable!");
        status = model_convert (data->models, data->models_stride.length, model_meanings, model_values, num_models);
        CHKERROUT (status, "Cannot convert models map!");
        for (size_t j = 0; j < num_models; j++) {
          mc_surface_method method = find_surface_model_id (model_meanings[j]);
          set_requested_parameters_for_model (method, &requested_parameters[0]);
        }

        free (model_values);
        free (model_meanings[0]);
        free (model_meanings);
      }
    }
  }

  for (size_t i = 0; i < n_filenames; i++) {
    INFOV (quiet, "Loading groups from NetCDF file '%s'", filenames[i]);

    status = nc_inq_grps (ncid[i], &n_groups, NULL);
    NC_CHKERROUT (status, "Cannot inquiry groups!");
    groups = calloc (n_groups, sizeof (int));
    status = nc_inq_grps (ncid[i], NULL, groups);
    NC_CHKERROUT (status, "Cannot inquiry group ids!");

    for (int j = 0; j < n_groups; j++) {
      char grp_name[NC_MAX_NAME + 1];
      status = nc_inq_grpname (groups[j], &grp_name[0]);
      NC_CHKERROUTV (status, "Cannot inquiry group %d!", groups[j]);
      if (strcmp ("maps", &grp_name[0]) == 0) {
        status = nc_inq_grps (groups[j], &n_model_groups, NULL);
        NC_CHKERROUT (status, "Cannot inquiry subgroups!");
        model_groups = calloc (n_model_groups, sizeof (int));
        status       = nc_inq_grps (groups[j], NULL, model_groups);
        NC_CHKERROUT (status, "Cannot inquiry subgroup ids!");

        for (int k = 0; k < n_model_groups; k++) {
          char model_name[NC_MAX_NAME + 1];
          status = nc_inq_grpname (model_groups[k], &model_name[0]);
          NC_CHKERROUTV (status, "Cannot inquiry subgroup %s/%d!", &grp_name[0], model_groups[k]);
          status = nc_inq (model_groups[k], NULL, &n_parameters, NULL, NULL);
          NC_CHKERROUTV (status, "Cannot inquiry variables of subgroup %s/%s!", &grp_name[0], &model_name[0]);

          for (int varid = 0; varid < n_parameters; varid++) {
            char                 parameter_name[NC_MAX_NAME + 1];
            mc_surface_parameter param_id;
            mc_surface_method    model_id;
            mc_ncdf_array_stride stride;
            float*               values;
            status = nc_inq_varname (model_groups[k], varid, &parameter_name[0]);
            NC_CHKERROUTV (status, "Cannot inquiry variable %d of group %s/%s!", varid, &grp_name[0], &model_name[0]);

            INFOV (quiet, "Loading parameter '%s/%s'", &model_name[0], &parameter_name[0]);

            param_id = find_model_param_id (&model_name[0], &parameter_name[0]);
            ASSERTOUTV (param_id != MCALB_INVALID_PARAM,
                        NCDF_INVALID_MODEL_PARAM,
                        "Invalid parameter %s/%s",
                        &model_name[0],
                        &parameter_name[0]);
            model_id = find_surface_model_id (&model_name[0]);
            ASSERTOUTV (model_id != MCALB_NONE, NCDF_INVALID_MODEL_PARAM, "Invalid model %s", &model_name[0]);
            status = load_dimensions (ncid[i],
                                      model_groups[k],
                                      varid,
                                      &dimension_ids[i][0],
                                      &dimension_values[0],
                                      &dimension_lengths[0],
                                      &data->projection,
                                      &stride,
                                      &parameter_name[0]);
            CHKERROUTV (status, "Cannot load dimensions from parameter '%s/%s'!", &model_name[0], &parameter_name[0]);
            status = load_variable (model_groups[k],
                                    varid,
                                    (void**) &values,
                                    NC_FLOAT,
                                    &accepted_parameter_types[0],
                                    stride.length,
                                    TGT_UNIT_ANY);
            CHKERROUTV (status, "Cannot load values from parameter '%s/%s'!", &model_name[0], &parameter_name[0]);
            status = add_parameter (data, model_id, param_id, values, &stride);
            CHKERROUTV (status, "Cannot add parameter '%s/%s' to surface data", &model_name[0], &parameter_name[0]);

            requested_parameters[param_id] = false;
          }
        }
        free (model_groups);
      }
    }

    free (groups);
    status = nc_close (ncid[i]);
    CHKERROUTV (status, "Cannot close netcdf file %s!", filenames[i]);
  }

  data->int_spectrum                         = NULL;
  data->int_nspectrum                        = 0;
  data->int_spectral_parameters.data         = NULL;
  data->int_spectral_spatial_parameters.data = NULL;
  reset_stride (&data->int_spectral_parameters.stride);
  reset_stride (&data->int_spectral_spatial_parameters.stride);

  for (mc_surface_parameter param = 0; param < MCALB_NUM_PARAMS; param++) {
    ASSERTOUTV (requested_parameters[param] == false,
                NCDF_MISSING_PARAMETER,
                "The parameter %s (%d) was requested by the 'models' field, but never provided in 'maps'\n",
                mc_model_param_str (param),
                param);
  }

  INFO (quiet, "NetCDF files sucessfully loaded");

  return NCDF_NOERR;
}

ncdf_status_t mc_ncdf_int_wl (mc_ncdf_surface* data, double wl_min, double wl_max, size_t count, int quiet) {
  ncdf_status_t status;
  double        delta;
  float*        spectrum = calloc (count, sizeof (float));

  delta = (wl_max - wl_min) / (count - 1);
  for (size_t i = 0; i < count; i++) {
    spectrum[i] = wl_min + i * delta;
  }

  status = mc_ncdf_int_wl_fromspec (data, spectrum, count, quiet);
  free (spectrum);
  return status;
}

ncdf_status_t mc_ncdf_int_wl_fromspec (mc_ncdf_surface* surface, float* spectrum, size_t count, int quiet) {
  ncdf_status_t status;

  INFOV (quiet, "Interpolating to [%f..%f] x %lu", spectrum[0], spectrum[count - 1], count);

  if (!surface->spectral_parameters.data && !surface->spectral_spatial_parameters.data) {
    // Nothing to convert
    return NCDF_NOERR;
  }

  ASSERTOUT (surface->spectrum, NCDF_INTERPOLATION_ERROR, "Cannot find spectrum. This is likely a bug!");
  ASSERTOUTV (surface->spectrum[0] <= spectrum[0],
              NCDF_INTERPOLATION_ERROR,
              "Error, minimum requested wavelength below data range (%f > %f)",
              surface->spectrum[0],
              spectrum[0]);
  ASSERTOUTV (surface->spectrum[surface->nspectrum - 1] >= spectrum[count - 1],
              NCDF_INTERPOLATION_ERROR,
              "Error, maximum requested wavelength above data range (%f < %f)\n",
              surface->spectrum[surface->nspectrum - 1],
              spectrum[count - 1]);

  surface->int_nspectrum = count;
  free (surface->int_spectrum);
  surface->int_spectrum = calloc (count, sizeof (double));
  for (size_t i = 0; i < count; i++) {
    surface->int_spectrum[i] = spectrum[i];
  }

  status = interpolate_wl (surface, &surface->spectral_parameters, &surface->int_spectral_parameters);
  CHKERROUT (status, "Error, cannot interpolate wavelengths of spectral parameters!");
  status = interpolate_wl (surface, &surface->spectral_spatial_parameters, &surface->int_spectral_spatial_parameters);
  CHKERROUT (status, "Error, cannot interpolate wavelengths of spatial parameters!");

  //free_raw_values (surface);

  INFOV (quiet, "Interpolation done [%f..%f] x %lu\n", spectrum[0], spectrum[count - 1], count);

  return NCDF_NOERR;
}

ncdf_status_t interpolate_wl (mc_ncdf_surface* surface, mc_ncdf_array* src, mc_ncdf_array* interpolated) {
  reset_stride (&interpolated->stride);
  calc_stride (&interpolated->stride.x, &interpolated->stride, surface->nx);
  calc_stride (&interpolated->stride.y, &interpolated->stride, surface->ny);
  calc_stride (&interpolated->stride.wl, &interpolated->stride, surface->int_nspectrum);
  calc_stride (&interpolated->stride.param, &interpolated->stride, src->max_param_count);

  memcpy (&interpolated->param_indices, &src->param_indices, sizeof (src->param_indices));
  memcpy (&interpolated->param_count, &src->param_count, sizeof (src->param_count));
  interpolated->max_param_count = src->max_param_count;

  free (interpolated->data);
  interpolated->data = calloc (src->max_param_count * surface->nx * surface->ny * surface->int_nspectrum, sizeof (float));
  for (size_t wl = 0, wlr_old = 1; wl < surface->int_nspectrum; wl++) {
    size_t wlr;
    double target_wl, wl_lower, wl_upper, weight;
    target_wl = surface->int_spectrum[wl];
    for (wlr = wlr_old; wlr < surface->nspectrum; wlr++) {
      if (surface->spectrum[wlr] > target_wl) {
        wlr_old = wlr;
        break;
      }
    }
    wl_lower = surface->spectrum[wlr - 1];
    wl_upper = surface->spectrum[wlr];
    weight   = (target_wl - wl_lower) / (wl_upper - wl_lower);

    for (size_t i = 0; i < surface->nx; i++) {
      for (size_t j = 0; j < surface->ny; j++) {
        for (size_t p = 0; p < src->max_param_count; p++) {
          double param_lower, param_upper, new_param;

          param_lower = get_param (src, i, j, wlr - 1, p);
          param_upper = get_param (src, i, j, wlr, p);
          new_param   = param_upper * weight + param_lower * (1. - weight);
          set_param (interpolated, i, j, wl, p, new_param);
        }
      }
    }
  }
  return NCDF_NOERR;
}

mc_surface_method mc_ncdf_get_model (mc_ncdf_surface* data, size_t i, size_t j) {
  return data->models[i * data->models_stride.x + j * data->models_stride.y];
}

float get_param (mc_ncdf_array* array, size_t i, size_t j, size_t wl, size_t param_index) {
  return array->data[i * array->stride.x + j * array->stride.y + wl * array->stride.wl + param_index * array->stride.param];
}

float get_param_by_stride (mc_ncdf_array_stride* stride, float* values, size_t i, size_t j, size_t wl, size_t param_index) {
  return values[i * stride->x + j * stride->y + wl * stride->wl + param_index * stride->param];
}

ncdf_status_t set_param (mc_ncdf_array* array, size_t i, size_t j, size_t wl, size_t param_index, float value) {
  array->data[i * array->stride.x + j * array->stride.y + wl * array->stride.wl + param_index * array->stride.param] = value;
  return NCDF_NOERR;
}

float mc_ncdf_get_param (mc_ncdf_surface* data, size_t i, size_t j, size_t wl, mc_surface_parameter param) {
  return get_param (data->param_link[param], i, j, wl, data->param_link[param]->param_indices[param]);
}

float mc_ncdf_get_param_int (mc_ncdf_surface* data, size_t i, size_t j, long wl, mc_surface_parameter param) {
  if (wl == -1) {
    return get_param (data->int_param_link[param], i, j, 0, data->int_param_link[param]->param_indices[param]);
  } else {
    return get_param (data->int_param_link[param], i, j, wl, data->int_param_link[param]->param_indices[param]);
  }
}

typedef float (*mc_ncdf_get_param_ptr_t) (mc_ncdf_surface*, size_t, size_t, size_t, mc_surface_parameter);

ncdf_status_t mc_ncdf_print_params_generic (mc_ncdf_surface*        data,
                                            size_t                  i,
                                            size_t                  j,
                                            long                    wl,
                                            mc_ncdf_get_param_ptr_t get_param_generic,
                                            double*                 spectrum) {
  mc_surface_method model = mc_ncdf_get_model (data, i, j);
  double            wavelength;
  if (wl > 0) {
    wavelength = spectrum[wl];
  } else {
    wavelength = spectrum[0];
  }
  switch (model) {
  case MCALB_LITVINOV:
    fprintf (stderr,
             "(%.4g,%.4g,%.4g): %s -> (%f, %f, %f, %f, %f, %f, %f)\n",
             data->x[i],
             data->y[j],
             wavelength,
             mc_model_str (model),
             (*get_param_generic) (data, i, j, wl, MCALB_LITVINOV_PARAM_ALBEDO),
             (*get_param_generic) (data, i, j, wl, MCALB_LITVINOV_PARAM_RMS_SLOPE),
             (*get_param_generic) (data, i, j, wl, MCALB_LITVINOV_PARAM_FRESNEL),
             (*get_param_generic) (data, i, j, wl, MCALB_LITVINOV_PARAM_SHADOWING),
             (*get_param_generic) (data, i, j, wl, MCALB_LITVINOV_PARAM_REFR_INDEX),
             (*get_param_generic) (data, i, j, wl, MCALB_LITVINOV_PARAM_REFR_INDEX_I),
             (*get_param_generic) (data, i, j, wl, MCALB_LITVINOV_PARAM_RMS_SLOPE2));
    break;
  case MCALB_TSANG:
    fprintf (stderr,
             "(%.4g,%.4g,%.4g): %s -> (%f)\n",
             data->x[i],
             data->y[j],
             wavelength,
             mc_model_str (model),
             (*get_param_generic) (data, i, j, wl, MCALB_TSANG_PARAM_WIND));
    break;
  case MCALB_MAIGN_BREON:
    fprintf (stderr,
             "(%.4g,%.4g,%.4g): %s -> (%f, %f, %f, %f), rpv -> (%f, %f, %f, %f, %f, %f, %f)\n",
             data->x[i],
             data->y[j],
             wavelength,
             mc_model_str (model),
             (*get_param_generic) (data, i, j, wl, MCALB_MAIGN_BREON_PARAM_C_MAIGN),
             (*get_param_generic) (data, i, j, wl, MCALB_MAIGN_BREON_PARAM_ARVI),
             (*get_param_generic) (data, i, j, wl, MCALB_MAIGN_BREON_PARAM_REFR_INDEX),
             (*get_param_generic) (data, i, j, wl, MCALB_MAIGN_BREON_PARAM_REFR_INDEX_I),
             (*get_param_generic) (data, i, j, wl, MCALB_RPV_PARAM_RHO0),
             (*get_param_generic) (data, i, j, wl, MCALB_RPV_PARAM_K),
             (*get_param_generic) (data, i, j, wl, MCALB_RPV_PARAM_THETA),
             (*get_param_generic) (data, i, j, wl, MCALB_RPV_PARAM_SCALE),
             (*get_param_generic) (data, i, j, wl, MCALB_RPV_PARAM_SIGMA),
             (*get_param_generic) (data, i, j, wl, MCALB_RPV_PARAM_T1),
             (*get_param_generic) (data, i, j, wl, MCALB_RPV_PARAM_T2));
    break;
  case MCALB_LAM:
    fprintf (stderr,
             "(%.4g,%.4g,%.4g): %s -> (%f)\n",
             data->x[i],
             data->y[j],
             wavelength,
             mc_model_str (model),
             (*get_param_generic) (data, i, j, wl, MCALB_LAM_PARAM_ALBEDO));
    break;
  default:
    fprintf (stderr,
             "(%.4g,%.4g,%.4g): %s[%d] -> no parameters available\n",
             data->x[i],
             data->y[j],
             wavelength,
             mc_model_str (model),
             model);
    return NCDF_PRINT_ERROR;
  }
  return NCDF_NOERR;
}

ncdf_status_t mc_ncdf_print_params (mc_ncdf_surface* data, size_t i, size_t j, size_t wl) {
  return mc_ncdf_print_params_generic (data, i, j, (long) wl, &mc_ncdf_get_param, data->spectrum);
}

float mc_ncdf_get_param_conv_int (mc_ncdf_surface* data, size_t i, size_t j, size_t wl, mc_surface_parameter param) {
  return mc_ncdf_get_param_int (data, i, j, wl, param);
}

ncdf_status_t mc_ncdf_print_params_int (mc_ncdf_surface* data, size_t i, size_t j, long wl) {
  return mc_ncdf_print_params_generic (data, i, j, wl > 0 ? wl : 0, &mc_ncdf_get_param_conv_int, data->int_spectrum);
}

#else  // HAVE_NETCDF4

// Fail all external functions if there was no netCDF 4 support at compile time
#define FAIL_NCDF4_MSG                                                                                                             \
  "Error, mystic has been compiled without netCDF which is\n"                                                                      \
  "required in order to open surface netCDF files.\n"
#define FAIL_NCDF CHKERROUT (NCDF_NO_NETCDF4, FAIL_NCDF4_MSG)
#define FAIL_NCDF_NO_RETURN chkerr_out (NCDF_NO_NETCDF4, FAIL_NCDF4_MSG, ERROR_POSITION)
#define FAIL_NCDF_DEF_RETURN(value) (chkerr_out (NCDF_NO_NETCDF4, FAIL_NCDF4_MSG, ERROR_POSITION), value)

ncdf_status_t mc_ncdf_load (char** filenames, mc_ncdf_surface* surface, int quiet) {
  FAIL_NCDF;
}
ncdf_status_t mc_ncdf_load_n (char** filenames, size_t n_filenames, mc_ncdf_surface* surface, int quiet) {
  FAIL_NCDF;
}
ncdf_status_t mc_ncdf_int_wl (mc_ncdf_surface* surface, double wl_min, double wl_max, size_t count, int quiet) {
  FAIL_NCDF;
}
ncdf_status_t mc_ncdf_int_wl_fromspec (mc_ncdf_surface* surface, float* spectrum, size_t count, int quiet) {
  FAIL_NCDF;
}
void mc_ncdf_free (mc_ncdf_surface* surface) {
  FAIL_NCDF_NO_RETURN;
}
mc_surface_method mc_ncdf_get_model (mc_ncdf_surface* surface, size_t i, size_t j) {
  FAIL_NCDF_DEF_RETURN (0);
}
float mc_ncdf_get_param (mc_ncdf_surface* surface, size_t i, size_t j, size_t wl, mc_surface_parameter param) {
  FAIL_NCDF_DEF_RETURN (0.);
}
float mc_ncdf_get_param_int (mc_ncdf_surface* surface, size_t i, size_t j, long wl, mc_surface_parameter param) {
  FAIL_NCDF_DEF_RETURN (0.);
}
ncdf_status_t mc_ncdf_print_params (mc_ncdf_surface* surface, size_t i, size_t j, size_t wl) {
  FAIL_NCDF;
}
ncdf_status_t mc_ncdf_print_params_int (mc_ncdf_surface* surface, size_t i, size_t j, long wl) {
  FAIL_NCDF;
}

#endif  // HAVE_NETCDF4

int reflection_matrix_mixed_surface_properties (albedo_struct* albedo,
                                                photon_struct* p,
                                                int            nstokes,
                                                int            backward,
                                                double         mu_inc,
                                                double         phi_inc,
                                                double         mu_esc,
                                                double         phi_esc,
                                                int            ia,
                                                int            ja,
                                                int            il,
                                                int            iv,
                                                /* Output */
                                                double*** refl_mat) {

  int    status     = 0;
  double wavelength = 0;

  switch (mc_ncdf_get_model (albedo->mixed_surface, ia, ja)) {
  case MCALB_LITVINOV:
    status = reflection_matrix_litvinov_mystic (
      mu_esc,
      mu_inc,
      phi_inc,
      phi_esc,
      mc_ncdf_get_param_int (albedo->mixed_surface, ia, ja, p->iv, MCALB_LITVINOV_PARAM_REFR_INDEX) +
        I * mc_ncdf_get_param_int (albedo->mixed_surface, ia, ja, p->iv, MCALB_LITVINOV_PARAM_REFR_INDEX_I),
      mc_ncdf_get_param_int (albedo->mixed_surface, ia, ja, p->iv, MCALB_LITVINOV_PARAM_ALBEDO),
      mc_ncdf_get_param_int (albedo->mixed_surface, ia, ja, p->iv, MCALB_LITVINOV_PARAM_RMS_SLOPE),
      mc_ncdf_get_param_int (albedo->mixed_surface, ia, ja, p->iv, MCALB_LITVINOV_PARAM_FRESNEL),
      mc_ncdf_get_param_int (albedo->mixed_surface, ia, ja, p->iv, MCALB_LITVINOV_PARAM_SHADOWING),
      mc_ncdf_get_param_int (albedo->mixed_surface, ia, ja, p->iv, MCALB_LITVINOV_PARAM_RMS_SLOPE2),
      backward,
      refl_mat);
    break;
  case MCALB_MAIGN_BREON:
    status = reflection_matrix_maignan_breon_mystic (
      mu_esc,
      mu_inc,
      phi_inc,
      phi_esc,
      mc_ncdf_get_param_int (albedo->mixed_surface, ia, ja, p->iv, MCALB_MAIGN_BREON_PARAM_REFR_INDEX) +
        I * mc_ncdf_get_param_int (albedo->mixed_surface, ia, ja, p->iv, MCALB_MAIGN_BREON_PARAM_REFR_INDEX_I),
      mc_ncdf_get_param_int (albedo->mixed_surface, ia, ja, p->iv, MCALB_MAIGN_BREON_PARAM_C_MAIGN),
      mc_ncdf_get_param_int (albedo->mixed_surface, ia, ja, p->iv, MCALB_MAIGN_BREON_PARAM_ARVI),
      mc_ncdf_get_param_int (albedo->mixed_surface, ia, ja, p->iv, MCALB_RPV_PARAM_RHO0),
      mc_ncdf_get_param_int (albedo->mixed_surface, ia, ja, p->iv, MCALB_RPV_PARAM_K),
      mc_ncdf_get_param_int (albedo->mixed_surface, ia, ja, p->iv, MCALB_RPV_PARAM_THETA),
      mc_ncdf_get_param_int (albedo->mixed_surface, ia, ja, p->iv, MCALB_RPV_PARAM_SCALE),
      mc_ncdf_get_param_int (albedo->mixed_surface, ia, ja, p->iv, MCALB_RPV_PARAM_SIGMA),
      mc_ncdf_get_param_int (albedo->mixed_surface, ia, ja, p->iv, MCALB_RPV_PARAM_T1),
      mc_ncdf_get_param_int (albedo->mixed_surface, ia, ja, p->iv, MCALB_RPV_PARAM_T2),
      backward,
      refl_mat);
    break;
  case MCALB_LAM:
    (*refl_mat)[0][0] = mc_ncdf_get_param_int (albedo->mixed_surface, ia, ja, p->iv, MCALB_LAM_PARAM_ALBEDO);
    break;
  case MCALB_TSANG:
    /* Wavelength in micron */
    wavelength = 0.5 * 1e4 * (1.0 / p->wvnmlo + 1.0 / p->wvnmhi);
    if (backward)
      status = bpdf_tsang (mc_ncdf_get_param_int (albedo->mixed_surface, ia, ja, p->iv, MCALB_TSANG_PARAM_WIND),
                           wavelength,
                           mu_esc,
                           phi_esc,
                           mu_inc,
                           phi_inc,
                           refl_mat);
    else
      status = bpdf_tsang (mc_ncdf_get_param_int (albedo->mixed_surface, ia, ja, p->iv, MCALB_TSANG_PARAM_WIND),
                           wavelength,
                           mu_inc,
                           phi_inc,
                           mu_esc,
                           phi_esc,
                           refl_mat);
    break;
  }
  return status;
}
