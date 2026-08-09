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
 * Authors of complex srfc: Marc Schwaerzel, Claudia Emde and Fabian Jakub
 ************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

typedef struct t_star_engine {
#ifdef HAVE_STAR_ENGINE
  struct s3d_device*     dev;
  struct s3d_scene_view* view;
  struct s3d_scene*      scene;
#else
  char dummy;  // to prevent empty struct
#endif
} t_star_engine;

typedef struct t_native_engine {
  size_t**  cell2map;       // foreach cell adresses into triangle map
  size_t*   triangle_map;   // mapping of triangles that cross a bounding box
  size_t*** kij2cellid;     // cellid for ragged kij array structure
  size_t*   k2cellid;       // cellid for k-only array, i.e. only use it for 1D layers
  size_t (*cellid2kij)[3];  // indices kij for a particular cellid, size(Ncells,3)
  size_t Ncells;            // Global number of voxels
  int*   layer_is_3D;       // flag if layer is 1D/3D, size(Nz), starting at TOA
} t_native_engine;

typedef struct t_triangle_radiation_field {
  size_t  N_triangles;
  size_t* ndir; /* number of photons per triangle */
  size_t* ndn;
  size_t* nup;

  double* edir; /* irradiance */
  double* edn;
  double* eup;
} t_triangle_radiation_field;

typedef struct {
  size_t N_vertices;
  size_t N_triangles;
  size_t N_materials;
  double (*vertices)[3];
  unsigned int (*triangles)[3];
  unsigned int*           material_of_triangle;
  double*                 material_albedo;
  double*                 temperature_of_triangle;
  double*                 Bplanck_of_triangle;
  struct t_star_engine*   star_engine_tracer;
  struct t_native_engine* native_tracer;
} t_triangular_surface;

int setup_triangular_surface (char* filename, t_triangular_surface* triangular_surface, int quiet);

int init_triangular_surface_result_struct (const size_t N_triangles, t_triangle_radiation_field** result);
int reset_triangular_surface_result_struct (t_triangle_radiation_field* result);

int init_spectral_triangular_surface_result_struct (const size_t                  N_wavelength,
                                                    const size_t                  N_triangles,
                                                    t_triangle_radiation_field*** result);

int free_triangular_surface_result_struct (t_triangle_radiation_field* result);

int add_triangular_surface_result (const double factor, t_triangle_radiation_field* add, t_triangle_radiation_field* result);

int init_native_tracer (int                   Nx,          // number of voxels in x dimension
                        int                   Ny,          // number of voxels in y dimension
                        int                   Nz,          // number of voxels in z dimension, i.e. layers
                        double                delX,        // voxel size in x dimension [m]
                        double                delY,        // voxel size in y dimension [m]
                        float*                heightlvls,  // level heights, starting at TOA [km]
                        int*                  threed,      // flag if a layer is 3D (size Nz), note that this one starts at the srfc
                        t_triangular_surface* srfc,        // will initialize the native engine struct in srfc
                        int                   quiet);

int cellid2boundingbox (int    k,              // voxel index in z dimension, i.e. layers
                        int    i,              // voxel index in x dimension, set to -Nx if k-layer is 1D
                        int    j,              // voxel index in y dimension, set to -Ny if k-layer is 1D
                        int    Nx,             // grid size in x-dimension
                        int    Ny,             // grid size in y-dimension
                        int    threedlayer,    // flag if k-layer is 1D or 3D
                        double delX,           // voxel size in x dimension [m]
                        double delY,           // voxel size in y dimension [m]
                        float* heightlvls,     // level heights, starting at the surface [km]
                        double boxcenter[3],   // returns cartesian center of bounding box
                        double boxhalfsize[3]  // and half distance to the edge of the box
);

int check_triangle_boundingbox_overlap (const size_t          iface,          // triangle index
                                        t_triangular_surface* srfc,           // triangle data
                                        const double          boxcenter[3],   // bounding box center
                                        const double          boxhalfsize[3]  // bounding box extent
);

int get_global_vertices_bounding_box (const t_triangular_surface* srfc, double min[3], double max[3]);

int get_triangle_bounding_box (const double A[3], const double B[3], const double C[3], double min[3], double max[3]);

int wiggle_vertices (const t_triangular_surface* srfc, const double wiggle);

int cross_triangular_surface (const t_triangular_surface* srfc,
                              const double                origin[3],
                              const double                dir[3],
                              const double                range[2],
                              const size_t                kij[3],
                              double*                     distance,
                              size_t*                     primitive_id);

int get_triangle_albedo (const t_triangular_surface* srfc, const size_t primitive_id, double* triangle_albedo);

int get_triangle_normal (const t_triangular_surface* srfc, const size_t primitive_id, double* triangle_normal);

int triangle_intersection_woop (const double origin[3],
                                const double direction[3],
                                const double A[3],
                                const double B[3],
                                const double C[3],
                                const double tnear,
                                const double tfar,
                                double*      distance,
                                double       hit[3]);

#ifdef __cplusplus
}
#endif
