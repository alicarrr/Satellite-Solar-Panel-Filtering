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
 * Authors of triangle srfc: Marc Schwaerzel, Claudia Emde and Fabian Jakub
 ************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if HAVE_NETCDF4
#include <netcdf.h>

#define RUN_NC(command, msg)                                                                                                       \
  do {                                                                                                                             \
    int status = (command);                                                                                                        \
    if (status != NC_NOERR) {                                                                                                      \
      fprintf (stderr, "netCDF error %d: \"%s\" %s from %s\n", status, nc_strerror (status), (msg), filename);                     \
      return status;                                                                                                               \
    }                                                                                                                              \
  } while (0)

#endif

#include "allocnd.h"
#include "common_math.h"
#include "triangle_surface.h"
#include "errors.h"
#include "tribox3.h"

#if HAVE_STAR_ENGINE
#include "rsys/mem_allocator.h"
#include "star/s3daw.h"
#include "star/s3d.h"
#endif

#if HAVE_STAR_ENGINE
/* Initialize Star-Engine objects used for raytracing */
void surf_get_indices (const unsigned primitive_id, unsigned indices[3], void* data) {
  /* The following array lists the indices toward the 3D vertices of each
   * triangle.

   *       ,2----,3
   *     ,'0 \1,'
   *    0----1'
   * e.g. id=1 -> indices={1,3,2}
   */
  const t_triangular_surface* srfc = (t_triangular_surface*) data;
  if (primitive_id >= srfc->N_triangles)
    fprintf (stderr, "surf_get_indices %ud out of bounds %lu\n", primitive_id, srfc->N_triangles);
  indices[0] = srfc->triangles[primitive_id][0];
  indices[1] = srfc->triangles[primitive_id][1];
  indices[2] = srfc->triangles[primitive_id][2];
}

void surf_get_position (const unsigned vertex_id, float position[3], void* data) {
  /* Memory layout of the vertex positions
   *   2.....3
   *  /     /
   * 0-----1
  */
  const t_triangular_surface* srfc = (t_triangular_surface*) data;
  if (vertex_id >= srfc->N_vertices)
    fprintf (stderr, "surf_get_position %ud out of bounds %lu\n", vertex_id, srfc->N_vertices);
  position[0] = srfc->vertices[vertex_id][0];
  position[1] = srfc->vertices[vertex_id][1];
  position[2] = srfc->vertices[vertex_id][2];
}

int s4vs_discard_self_hit (const struct s3d_hit* hit,
                           const float           ray_org[3],
                           const float           ray_dir[3],
                           void*                 ray_data,
                           void*                 filter_data) {
  const struct s3d_primitive* prim_from = ray_data;

  /* Avoid unused variable warn */
  (void) ray_org, (void) ray_dir, (void) filter_data;
  return prim_from ? S3D_PRIMITIVE_EQ (prim_from, &hit->prim) : 0;
}

int set_attach_surf (struct s3d_device* s3d, struct s3d_scene* scene, t_triangular_surface* srfc) {

  int                    ierr;
  struct s3d_shape*      surf;
  struct s3d_vertex_data vertex_attribs[2];

  vertex_attribs[0].usage = S3D_POSITION;
  vertex_attribs[0].type  = S3D_FLOAT3;
  vertex_attribs[0].get   = surf_get_position;
  vertex_attribs[1]       = S3D_VERTEX_DATA_NULL;

  ierr = s3d_shape_create_mesh (s3d, &surf);
  CHKERR (ierr);

  const unsigned ntris  = srfc->N_triangles;
  const unsigned nverts = srfc->N_vertices;

  ierr = s3d_mesh_setup_indexed_vertices (surf, ntris, surf_get_indices, nverts, vertex_attribs, 1, (void*) srfc);
  CHKERR (ierr);

  ierr = s3d_mesh_set_hit_filter_function (surf, s4vs_discard_self_hit, (void*) srfc);
  CHKERR (ierr);

  ierr = s3d_scene_attach_shape (scene, surf);
  CHKERR (ierr);

  return 0;
  ierr = s3d_shape_ref_put (surf);
  CHKERR (ierr);

  return 0;
}

int init_star_engine (t_triangular_surface* triangular_surface, struct t_star_engine** star_engine, int quiet) {
  int                    ierr;
  struct s3d_device*     dev;
  struct s3d_scene_view* view;
  struct s3d_scene*      scene;
  if (!quiet)
    fprintf (stderr, "Initializing Star-Engine-Tracer!\n");

  /*initialize device*/
  ierr = s3d_device_create (NULL, NULL, 0, &dev);
  CHKERR (ierr);

  /*initialize scene*/
  ierr = s3d_scene_create (dev, &scene);
  CHKERR (ierr);

  ierr = set_attach_surf (dev, scene, triangular_surface);
  CHKERR (ierr);

  /*initialize geometry*/
  ierr = s3d_scene_view_create (scene, S3D_TRACE, &view);
  CHKERR (ierr);

  *star_engine          = malloc (sizeof (struct t_star_engine));
  (*star_engine)->dev   = dev;
  (*star_engine)->view  = view;
  (*star_engine)->scene = scene;
  S3D (device_ref_put (dev));
  S3D (scene_ref_put (scene));
  if (!quiet)
    fprintf (stderr, "Initializing Star-Engine-Tracer ... done!\n");
  return 0;
}
#endif

// transforms a cellid to a bounding box specification
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
) {
  if (threedlayer) {
    if (i < 0)
      return err_out ("cannot have negative `i` index in threed layers", i);
    if (j < 0)
      return err_out ("cannot have negative `j` index in threed layers", j);
    boxcenter[0]   = (i + 0.5) * delX;
    boxhalfsize[0] = delX * 0.5 + FLT_EPSILON;
    boxcenter[1]   = (j + 0.5) * delY;
    boxhalfsize[1] = delY * 0.5 + FLT_EPSILON;
  } else {
    boxcenter[0]   = delX * Nx * 0.5;
    boxhalfsize[0] = delX * Nx * 0.5 + FLT_EPSILON;
    boxcenter[1]   = delY * Ny * 0.5;
    boxhalfsize[1] = delY * Ny * 0.5 + FLT_EPSILON;
  }

  boxcenter[2]   = (heightlvls[k + 1] + heightlvls[k]) * .5 * 1e3;
  boxhalfsize[2] = (heightlvls[k + 1] - heightlvls[k]) * .5 * 1e3;

  // increase bounding box size by a small increment to make sure that edge cases get caught
  const double EPS = sqrt (FLT_EPSILON);
  for (size_t k = 0; k < 3; ++k)
    boxhalfsize[k] += EPS;

  return 0;
}

int check_triangle_boundingbox_overlap (const size_t          iface,          // triangle index
                                        t_triangular_surface* srfc,           // triangle data
                                        const double          boxcenter[3],   // bounding box center
                                        const double          boxhalfsize[3]  // bounding box extent
) {
  const int ldebug = 0;  // enable to see a lot of debug output

  const int use_axis_aligned_check = 0;

  // indices of vertices of triangle face
  const size_t vA = srfc->triangles[iface][0];
  const size_t vB = srfc->triangles[iface][1];
  const size_t vC = srfc->triangles[iface][2];
  // coordinates of vertices
  const double A[3] = {srfc->vertices[vA][0], srfc->vertices[vA][1], srfc->vertices[vA][2]};
  const double B[3] = {srfc->vertices[vB][0], srfc->vertices[vB][1], srfc->vertices[vB][2]};
  const double C[3] = {srfc->vertices[vC][0], srfc->vertices[vC][1], srfc->vertices[vC][2]};

  int intersec;
  if (use_axis_aligned_check) {
    const double box_min[3] = {boxcenter[0] - boxhalfsize[0], boxcenter[1] - boxhalfsize[1], boxcenter[2] - boxhalfsize[2]};
    const double box_max[3] = {boxcenter[0] + boxhalfsize[0], boxcenter[1] + boxhalfsize[1], boxcenter[2] + boxhalfsize[2]};
    double       vbox_min[3], vbox_max[3];
    const int    ierr = get_triangle_bounding_box (A, B, C, vbox_min, vbox_max);
    CHKERR (ierr);
    intersec = check_aabb_intersection (vbox_min, vbox_max, box_min, box_max);
  } else {
    intersec = triBoxOverlap (boxcenter, boxhalfsize, A, B, C);
  }

  if (intersec && ldebug) {
    fprintf (stderr, "A %8.2f %8.2f %8.2f\n", A[0], A[1], A[2]);
    fprintf (stderr, "B %8.2f %8.2f %8.2f\n", B[0], B[1], B[2]);
    fprintf (stderr, "C %8.2f %8.2f %8.2f\n", C[0], C[1], C[2]);
    fprintf (stderr, "center %f %f %f\n", boxcenter[0], boxcenter[1], boxcenter[2]);
    fprintf (stderr, "boxhalfsize %f %f %f\n", boxhalfsize[0], boxhalfsize[1], boxhalfsize[2]);
    fprintf (stderr, "face: %lu intersects? %i \n", iface, intersec);
  }

  return intersec;
}

// For testing purposes: if the dimension of all vertices of a triangle are the same number,
// we can wiggle the values a wee bit because we have numerical problems if the triangles coincide
// with mystic layer borders
int wiggle_vertices (const t_triangular_surface* srfc, const double wiggle) {
  for (size_t i = 0; i < srfc->N_triangles; ++i) {
    double* A = srfc->vertices[srfc->triangles[i][0]];
    double* B = srfc->vertices[srfc->triangles[i][1]];
    double* C = srfc->vertices[srfc->triangles[i][2]];
    for (size_t k = 0; k < 3; ++k) {
      if (A[k] == B[k] && B[k] == C[k]) {
        A[k] -= wiggle;
        C[k] += wiggle;
      }
    }
  }
  return 0;
}

/* Building accelerator structure for native tracer engine triangle lookups.
   The idea is to have a mapping for which triangles intersect or lie inside the bounding box of a cell
   For that, we iterate over each cell, determine the intersecting triangles and create a mapping data structure
   which can be used for an efficient lookup
*/
int init_native_tracer (int                   Nx,          // number of voxels in x dimension
                        int                   Ny,          // number of voxels in y dimension
                        int                   Nz,          // number of voxels in z dimension, i.e. layers
                        double                delX,        // voxel size in x dimension [m]
                        double                delY,        // voxel size in y dimension [m]
                        float*                heightlvls,  // level heights, starting at TOA [km]
                        int*                  threed,      // flag if a layer is 3D (size Nz), note that this one starts at the srfc
                        t_triangular_surface* srfc,        // will initialize the native engine struct in srfc
                        int                   quiet) {

  const int ldebug = 0;  // enable to see a lot of debug output

  if (!srfc->N_triangles)
    return 0;
  if (!quiet)
    fprintf (stderr,
             " ... Setup of Triangular surface accelerators -\n     Nx %i Ny %i Nz %i dx,dy (%6.2f %6.2f [m]) domain size (%8.3f "
             "%8.3f [km])\n",
             Nx,
             Ny,
             Nz,
             delX,
             delY,
             Nx * delX * 1e-3,
             Ny * delY * 1e-3);

  if (Nx < 0)
    return err_out ("Nx cannot be negative!", Nx);
  if (Ny < 0)
    return err_out ("Ny cannot be negative!", Ny);
  if (Nz < 0)
    return err_out ("Nz cannot be negative!", Nz);

#if HAVE_STAR_ENGINE
  return 0;
#endif

  // reverse the heightlvls... we stick to mystic conventions here where vectors start at the surface
  float heightlvls_rev[Nz + 1];
  for (size_t k = 0; k < Nz + 1; ++k) {
    size_t rev_k      = Nz - k;
    heightlvls_rev[k] = heightlvls[rev_k];
  }

  if (srfc->native_tracer) {
    fprintf (stderr, "native_engine already allocated?\n");
    CHKERR (1);
  }
  srfc->native_tracer = malloc (sizeof (struct t_native_engine));

  // save threed on our own because we need it later on
  srfc->native_tracer->layer_is_3D = malloc (Nz * sizeof (int));
  for (size_t k = 0; k < Nz; ++k) {
    srfc->native_tracer->layer_is_3D[k] = threed[k];
  }

  {  // First we initialize a mapping from the ragged datastructure of libRadtran to allow for a linearized indexing
    // count the global number of voxels
    srfc->native_tracer->Ncells = 0;
    for (size_t k = 0; k < Nz; ++k) {
      if (srfc->native_tracer->layer_is_3D[k]) {
        for (size_t i = 0; i < Nx; ++i)
          for (size_t j = 0; j < Ny; ++j)
            srfc->native_tracer->Ncells++;
      } else {
        srfc->native_tracer->Ncells++;
      }
    }
    if (ldebug)
      fprintf (stderr, "Global scene consists of %lu voxels\n", srfc->native_tracer->Ncells);
    srfc->native_tracer->cellid2kij = malloc (srfc->native_tracer->Ncells * sizeof (size_t[3]));

    srfc->native_tracer->kij2cellid = malloc (Nz * sizeof (size_t**));
    srfc->native_tracer->k2cellid   = malloc (Nz * sizeof (size_t));
    size_t cellid                   = 0;
    for (size_t k = 0; k < Nz; ++k) {
      size_t NX, NY;
      if (srfc->native_tracer->layer_is_3D[k]) {
        NX = Nx;
        NY = Ny;
      } else {
        NX = 1;
        NY = 1;
      }
      srfc->native_tracer->kij2cellid[k] = malloc (NX * sizeof (size_t*));
      for (size_t i = 0; i < NX; ++i) {
        srfc->native_tracer->kij2cellid[k][i] = malloc (NY * sizeof (size_t));
        for (size_t j = 0; j < NY; ++j) {
          srfc->native_tracer->kij2cellid[k][i][j]   = cellid;
          srfc->native_tracer->k2cellid[k]           = cellid;
          srfc->native_tracer->cellid2kij[cellid][0] = k;
          srfc->native_tracer->cellid2kij[cellid][1] = i;
          srfc->native_tracer->cellid2kij[cellid][2] = j;
          cellid++;
        }
      }
    }

    if (ldebug) {
      for (size_t k = 0; k < Nz; ++k) {
        if (srfc->native_tracer->layer_is_3D[k]) {
          for (size_t i = 0; i < Nx; ++i) {
            for (size_t j = 0; j < Ny; ++j) {
              fprintf (stderr,
                       "3D Layer @ k,i,j (%lu %lu %lu) has cellid %lu\n",
                       k,
                       i,
                       j,
                       srfc->native_tracer->kij2cellid[k][i][j]);
            }
          }
        } else {
          fprintf (stderr, "1D Layer @ k = %lu has cellid %lu\n", k, srfc->native_tracer->kij2cellid[k][0][0]);
        }
      }

      for (size_t cellid = 0; cellid < srfc->native_tracer->Ncells; ++cellid) {
        fprintf (stderr,
                 "CellID %lu has indices: %lu %lu %lu\n",
                 cellid,
                 srfc->native_tracer->cellid2kij[cellid][0],
                 srfc->native_tracer->cellid2kij[cellid][1],
                 srfc->native_tracer->cellid2kij[cellid][2]);
      }
    }
  }

  // Now build accelerator structure
  if (srfc->N_triangles > 1e4) {
    fprintf (stderr,
             " ... ATTENTION: Generating the accelerator structures for the native tracer for %lu triangles will take a long time. "
             "You better use the Star-Engine tracer!\n",
             srfc->N_triangles);
  }
  size_t count                   = 0;
  size_t max_intersects_per_cell = 0;
  {
    if (ldebug) {
      fprintf (stderr, "First count how big the triangle mapping has to be, i.e., do a dry run\n");
    }
    double vbox_min[3], vbox_max[3];
    {
      const int ierr = get_global_vertices_bounding_box (srfc, vbox_min, vbox_max);
      CHKERR (ierr);
      if (ldebug) {
        fprintf (stderr, "vertices bounding box %f %f %f\n", vbox_min[0], vbox_min[1], vbox_min[2]);
        fprintf (stderr, "vertices bounding box %f %f %f\n", vbox_max[0], vbox_max[1], vbox_max[2]);
      }
    }
    for (size_t cellid = 0; cellid < srfc->native_tracer->Ncells; ++cellid) {
      size_t k = srfc->native_tracer->cellid2kij[cellid][0];
      size_t i = srfc->native_tracer->cellid2kij[cellid][1];
      size_t j = srfc->native_tracer->cellid2kij[cellid][2];

      double boxcenter[3];
      double boxhalfsize[3];
      // first check if 1D cell box intersects global vertices bounding box:
      {
        int ierr = cellid2boundingbox (k, -1, -1, Nx, Ny, 0, delX, delY, heightlvls_rev, boxcenter, boxhalfsize);
        CHKERR (ierr);
        const double box_min[3] = {boxcenter[0] - boxhalfsize[0], boxcenter[1] - boxhalfsize[1], boxcenter[2] - boxhalfsize[2]};
        const double box_max[3] = {boxcenter[0] + boxhalfsize[0], boxcenter[1] + boxhalfsize[1], boxcenter[2] + boxhalfsize[2]};
        if (!check_aabb_intersection (vbox_min, vbox_max, box_min, box_max)) {

          if (ldebug) {
            fprintf (stderr, "layer %lu cannot have any intersections because bounding boxes do not overlap\n", k);
          }
          continue;
        }
      }

      int is3D = srfc->native_tracer->layer_is_3D[k];
      int ierr = cellid2boundingbox (k, i, j, Nx, Ny, is3D, delX, delY, heightlvls_rev, boxcenter, boxhalfsize);
      CHKERR (ierr);
      size_t intersects_per_cell = 0;
      for (size_t iface = 0; iface < srfc->N_triangles; ++iface) {
        int intersec = check_triangle_boundingbox_overlap (iface, srfc, boxcenter, boxhalfsize);
        count += intersec;
        intersects_per_cell += intersec;
      }
      max_intersects_per_cell = intersects_per_cell > max_intersects_per_cell ? intersects_per_cell : max_intersects_per_cell;
    }

    if (max_intersects_per_cell > 100) {
      fprintf (stderr, "Using the native triangle intersection code will take forever... \n");
      fprintf (stderr, "I urge you to use a more optimized version!\n");
      fprintf (stderr, "E.g. you can use the StarEngine raytracing code in libRadtran.\n");
      fprintf (stderr, "Consult the documentation or ask Fabian Jakub on how to install and use it.\n");
      CHKERR (max_intersects_per_cell);
    }

    if (!quiet)
      fprintf (stderr,
               " ... in total, the %lu triangles intersect %lu times with cells. (max %lu intersections per cell)\n",
               srfc->N_triangles,
               count,
               max_intersects_per_cell);
  }

  // Then allocate memory for the mapping and start building it
  srfc->native_tracer->triangle_map = malloc ((count) * sizeof (size_t));
  srfc->native_tracer->cell2map     = malloc ((srfc->native_tracer->Ncells + 1) * sizeof (size_t*));

  if (ldebug) {
    fprintf (stderr,
             "Start/End of triangle_map, %p - %p\n",
             srfc->native_tracer->triangle_map,
             srfc->native_tracer->triangle_map + count);
  }

  srfc->native_tracer->cell2map[0] = srfc->native_tracer->triangle_map;
  size_t* current_map_loc          = srfc->native_tracer->triangle_map;
  {
    for (size_t cellid = 0; cellid < srfc->native_tracer->Ncells; ++cellid) {
      size_t k    = srfc->native_tracer->cellid2kij[cellid][0];
      size_t i    = srfc->native_tracer->cellid2kij[cellid][1];
      size_t j    = srfc->native_tracer->cellid2kij[cellid][2];
      int    is3D = srfc->native_tracer->layer_is_3D[k];

      double boxcenter[3];
      double boxhalfsize[3];
      int    ierr = cellid2boundingbox (k, i, j, Nx, Ny, is3D, delX, delY, heightlvls_rev, boxcenter, boxhalfsize);
      CHKERR (ierr);

      for (size_t iface = 0; iface < srfc->N_triangles; ++iface) {
        int intersec = check_triangle_boundingbox_overlap (iface, srfc, boxcenter, boxhalfsize);
        if (intersec) {
          *current_map_loc = iface;
          current_map_loc++;
        }
      }
      srfc->native_tracer->cell2map[cellid + 1] = current_map_loc;  // set end of this cell, i.e. start of next cell
      if (ldebug) {
        fprintf (stderr,
                 "writing end of cell mapping for %lu (%lu,%lu,%lu) to cell2map at %p with triangle_map target %p\n",
                 cellid + 1,
                 i,
                 j,
                 k,
                 srfc->native_tracer->cell2map + cellid + 1,
                 current_map_loc);
      }
    }
  }

  if (ldebug) {  // debug output, this shows how to determine all triangles that intersect or lie inside a cell bounding box
    for (size_t icell = 0; icell < srfc->native_tracer->Ncells; ++icell) {
      size_t* start = srfc->native_tracer->cell2map[icell];
      size_t* end   = srfc->native_tracer->cell2map[icell + 1];
      fprintf (stderr, "cell %ld (start/end pointers in map %p %p) intersects with triangles", icell, start, end);
      for (size_t* ptr = start; ptr < end; ++ptr)
        fprintf (stderr, ", %lu", *ptr);
      fprintf (stderr, "\n");
    }
  }

  if (!quiet)
    fprintf (stderr, " ... finished Native-Engine-Tracer setup\n");
  return 0;
}

/* Compute ray/triangle intersection with the native tracer engine */
int native_trace_ray (                    // return true if iface is hit
  const t_triangular_surface* srfc,       // triangle data
  const size_t                iface,      // triangle index
  const double                origin[3],  // cartesian origin of ray
  const double                dir[3],     // cartesian direction of ray
  const double                range[2],   // range in which we look for intersections
  double*                     distance    // output distance to intersection
) {

  // indices of vertices of triangle face
  const size_t vA = srfc->triangles[iface][0];
  const size_t vB = srfc->triangles[iface][1];
  const size_t vC = srfc->triangles[iface][2];
  // coordinates of vertices
  const double A[3] = {srfc->vertices[vA][0], srfc->vertices[vA][1], srfc->vertices[vA][2]};
  const double B[3] = {srfc->vertices[vB][0], srfc->vertices[vB][1], srfc->vertices[vB][2]};
  const double C[3] = {srfc->vertices[vC][0], srfc->vertices[vC][1], srfc->vertices[vC][2]};

  double    hit[3];
  const int lhit = triangle_intersection_woop (origin, dir, A, B, C, range[0], range[1], distance, hit);
  return lhit;
}

/*******************************************************************/
/* Read surface data and setup triangular surface.                 */
/*******************************************************************/

int setup_triangular_surface (char* filename,  // netcdf file with grid info, see doc or example for necessary structure
                              t_triangular_surface* triangular_surface,  // output struct to hold necessary info
                              int                   quiet) {

#if HAVE_NETCDF4
  const int ldebug = 0;  // enable to see a lot of debug output
  int       ncid   = 0;
  int       idd;

  size_t n = 0;

  /* open as netCDF file */
  int status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status == NC_NOERR) {

    if (!quiet)
      fprintf (stderr, " ... reading triangular surface data from netCDF file %s\n", filename);

    RUN_NC (nc_inq_dimid (ncid, "Nvert", &idd), "reading Nvert");
    RUN_NC (nc_inq_dimlen (ncid, idd, &n), "reading Nvert");
    triangular_surface->N_vertices = n;

    RUN_NC (nc_inq_dimid (ncid, "Ntriangles", &idd), "reading Ntriangles");
    RUN_NC (nc_inq_dimlen (ncid, idd, &n), "reading Ntriangles");
    triangular_surface->N_triangles = n;

    RUN_NC (nc_inq_dimid (ncid, "N_materials", &idd), "reading N_materials");
    RUN_NC (nc_inq_dimlen (ncid, idd, &n), "reading N_materials");
    triangular_surface->N_materials = n;

    RUN_NC (nc_inq_varid (ncid, "vertices", &idd), "reading coordinates of vertices");
    if (!(triangular_surface->vertices = malloc (triangular_surface->N_vertices * sizeof (double[3]))))
      CHKERROUT (-1, "could not allocate memory for triangular_surface->vertices");
    RUN_NC (nc_get_var_double (ncid, idd, *triangular_surface->vertices), "reading coordinates of vertices");

    RUN_NC (nc_inq_varid (ncid, "triangles", &idd), "reading face indices");
    if (!(triangular_surface->triangles = malloc (triangular_surface->N_triangles * sizeof (unsigned int[3]))))
      CHKERROUT (-1, "could not allocate memory for triangular_surface->triangles");
    RUN_NC (nc_get_var_uint (ncid, idd, *triangular_surface->triangles), "reading face indices");

    RUN_NC (nc_inq_varid (ncid, "material_albedo", &idd), "material_albedo");
    if (!(triangular_surface->material_albedo = malloc (triangular_surface->N_materials * sizeof (double))))
      CHKERROUT (-1, "could not allocate memory for triangular_surface->material_albedo");
    RUN_NC (nc_get_var_double (ncid, idd, triangular_surface->material_albedo), "reading material_albedo");

    RUN_NC (nc_inq_varid (ncid, "material_of_triangle", &idd), "reading material indices of triangles");
    if (!(triangular_surface->material_of_triangle = malloc (triangular_surface->N_triangles * sizeof (unsigned int))))
      CHKERROUT (-1, "could not allocate memory for triangular_surface->material_of_triangle");
    RUN_NC (nc_get_var_uint (ncid, idd, triangular_surface->material_of_triangle), "reading material indices of triangles");

    // Optional input data, the temperature of each triangle
    triangular_surface->temperature_of_triangle = NULL;
    triangular_surface->Bplanck_of_triangle     = NULL;
    {
      int ierr = nc_inq_varid (ncid, "temperature_of_triangle", &idd);
      if (ierr == NC_NOERR) {
        if (!(triangular_surface->temperature_of_triangle = malloc (triangular_surface->N_triangles * sizeof (double))))
          CHKERROUT (-1, "could not allocate memory for triangular_surface->temperature_of_triangle");
        if (!(triangular_surface->Bplanck_of_triangle = malloc (triangular_surface->N_triangles * sizeof (double))))
          CHKERROUT (-1, "could not allocate memory for triangular_surface->Bplanck_of_triangle");
        RUN_NC (nc_get_var_double (ncid, idd, triangular_surface->temperature_of_triangle), "reading temperature of triangles");
      }
    }
    nc_close (ncid);

    if (0) {  // For testing purposes: wiggle at the vertices if they are axis aligned surfaces
      int ierr = wiggle_vertices (triangular_surface, FLT_EPSILON);
      CHKERR (ierr);
    }

    if (ldebug) {

      for (size_t i = 0; i < triangular_surface->N_vertices; ++i)
        fprintf (stderr,
                 "vertex[%5lu] coords (%8.2f %8.2f %8.2f)\n",
                 i,
                 triangular_surface->vertices[i][0],
                 triangular_surface->vertices[i][1],
                 triangular_surface->vertices[i][2]);

      for (size_t i = 0; i < triangular_surface->N_triangles; ++i)
        fprintf (stderr,
                 "triangle[%5lu] coords (%5u %5u %5u) material %5d albedo %7.6f temperature %6.2f\n",
                 i,
                 triangular_surface->triangles[i][0],
                 triangular_surface->triangles[i][1],
                 triangular_surface->triangles[i][2],
                 triangular_surface->material_of_triangle[i],
                 triangular_surface->material_albedo[triangular_surface->material_of_triangle[i]],
                 triangular_surface->temperature_of_triangle[i]);

      for (size_t i = 0; i < triangular_surface->N_materials; ++i)
        fprintf (stderr, "material_albedo[%5lu] %.6f\n", i, triangular_surface->material_albedo[i]);
    }

  } else {
    fprintf (stderr,
             "Error: unable to read triangular surface properties file %s. It is expected to be in netcdf format. \n",
             filename);
    CHKERR (-1);
  }
#else
  fprintf (stderr, "Error: libRadtran was built without NETCDF support. Unable to read netcdf format file %s\n", filename);
  CHKERR (-1);
#endif

  triangular_surface->star_engine_tracer = NULL;
#if HAVE_STAR_ENGINE
  status = init_star_engine (triangular_surface, &(triangular_surface->star_engine_tracer), quiet);
  CHKERR (status);
  CHKPOINTEROUT (triangular_surface->star_engine_tracer, "Could not initialize Star-Engine tracer");
#endif

  return 0;
}

// allocate memory for triangular surface results
int init_triangular_surface_result_struct (const size_t N_triangles, t_triangle_radiation_field** result) {
  if (*result)
    CHKERROUT (-1, "triangular_surface_result_struct already allocated?!");
  if (N_triangles == 0)
    return 0;
  *result = malloc (sizeof (t_triangle_radiation_field));
  CHKPOINTER (*result);
  (*result)->N_triangles = N_triangles;
  (*result)->ndir        = calloc (N_triangles, sizeof (size_t));
  CHKPOINTER ((*result)->ndir);
  (*result)->ndn = calloc (N_triangles, sizeof (size_t));
  CHKPOINTER ((*result)->ndn);
  (*result)->nup = calloc (N_triangles, sizeof (size_t));
  CHKPOINTER ((*result)->nup);
  (*result)->edir = calloc (N_triangles, sizeof (double));
  CHKPOINTER ((*result)->edir);
  (*result)->edn = calloc (N_triangles, sizeof (double));
  CHKPOINTER ((*result)->edn);
  (*result)->eup = calloc (N_triangles, sizeof (double));
  CHKPOINTER ((*result)->eup);
  return 0;
}

int reset_triangular_surface_result_struct (t_triangle_radiation_field* result) {
  if (!result)
    CHKERROUT (-1, "triangular_surface_result_struct not allocated?!");
  for (size_t i = 0; i < result->N_triangles; ++i) {
    result->ndir[i] = 0;
    result->ndn[i]  = 0;
    result->nup[i]  = 0;
    result->edir[i] = 0;
    result->edn[i]  = 0;
    result->eup[i]  = 0;
  }
  return 0;
}

// allocate memory for triangular surface results for spectral intervals
int init_spectral_triangular_surface_result_struct (const size_t                  N_wavelength,
                                                    const size_t                  N_triangles,
                                                    t_triangle_radiation_field*** results) {

  if (*results)
    CHKERROUT (-1, "spectral triangular_surface_result_struct already allocated?!");
  if (N_triangles == 0)
    return 0;
  *results = malloc (N_wavelength * sizeof (t_triangle_radiation_field*));
  for (size_t iv = 0; iv < N_wavelength; ++iv) {
    (*results)[iv] = NULL;
    const int ierr = init_triangular_surface_result_struct (N_triangles, (*results) + iv);
    CHKERR (ierr);
  }
  return 0;
}

int free_triangular_surface_result_struct (t_triangle_radiation_field* result) {
  if (result) {
    free (result->ndir);
    result->ndir = NULL;
    free (result->ndn);
    result->ndn = NULL;
    free (result->nup);
    result->nup = NULL;
    free (result->edir);
    result->edir = NULL;
    free (result->edn);
    result->edn = NULL;
    free (result->eup);
    result->eup = NULL;
    free (result);
    result = NULL;
  }
  return 0;
}

int add_triangular_surface_result (const double factor, t_triangle_radiation_field* add, t_triangle_radiation_field* result) {
  if (!add)
    CHKERROUT (-1, "triangle_radiation_field* add not allocated");
  if (!result)
    CHKERROUT (-1, "triangle_radiation_field* result not allocated");
  CHKERROUT (add->N_triangles - result->N_triangles, "triangle_results have different sizes");

  for (size_t id = 0; id < add->N_triangles; ++id) {
    result->ndir[id] += add->ndir[id];
    result->ndn[id] += add->ndn[id];
    result->nup[id] += add->nup[id];
  }
  int ierr;
  ierr = vec_AXPY (result->N_triangles, factor, add->edir, result->edir);
  CHKERR (ierr);
  ierr = vec_AXPY (result->N_triangles, factor, add->edn, result->edn);
  CHKERR (ierr);
  ierr = vec_AXPY (result->N_triangles, factor, add->eup, result->eup);
  CHKERR (ierr);

  const int ldebug = 0;
  if (ldebug) {
    for (size_t id = 0; id < add->N_triangles; ++id) {
      fprintf (stderr,
               "                                                                Triangle result on patch %lu (factor %f) add (%f "
               "%f %f) is now (%f %f %f) count (%lu %lu %lu)\n",
               id,
               factor,
               factor * add->edir[id],
               factor * add->edn[id],
               factor * add->eup[id],
               result->edir[id],
               result->edn[id],
               result->eup[id],
               result->ndir[id],
               result->ndn[id],
               result->nup[id]);
    }
  }
  return 0;
}

/* cross_triangular_surface

Parameter:
+ srfc - triangular_surface context
+ origin - cartesian position of ray origin
+ dir - cartesian direction of ray
+ range - min and max distance from origin at which to search for intersecting faces
+ kij - cell indices of bounding box (ignored if starengine is used)
- distance - distance to first intersection (range[1] if no intersection happened)
- primitive_id - index of face at which the intersection ocurred (-1, i.e. huge(size_t) if no intersection happened)
 */
int cross_triangular_surface (const t_triangular_surface* srfc,
                              const double                origin[3],
                              const double                dir[3],
                              const double                range[2],
                              const size_t                kij[3],
                              double*                     distance,
                              size_t*                     primitive_id) {

  const int ldebug = 0;  // enable to see a lot of debug output

  *distance     = range[1];
  *primitive_id = -1;

#if HAVE_STAR_ENGINE
  if (srfc->star_engine_tracer) {
    struct s3d_hit hit        = S3D_HIT_NULL;
    const float    forigin[3] = {origin[0], origin[1], origin[2]};
    const float    fdir[3]    = {dir[0], dir[1], dir[2]};
    const float    frange[2]  = {range[0], range[1]};

    s3d_scene_view_trace_ray (srfc->star_engine_tracer->view, forigin, fdir, frange, NULL, &hit);
    if (!S3D_HIT_NONE (&hit)) {
      *distance     = hit.distance;
      *primitive_id = hit.prim.prim_id;
    }
    return 0;
  }
#endif

  // if we dont have STAR_ENGINE, use native O(N**2) triangle intersection code....
  size_t cellid;
  size_t k = kij[0];
  if (srfc->native_tracer->layer_is_3D[k]) {
    cellid = srfc->native_tracer->kij2cellid[k][kij[1]][kij[2]];
  } else {
    cellid = srfc->native_tracer->k2cellid[k];
  }

  if (ldebug)
    fprintf (stderr,
             "loc (%8.4f %8.4f %8.4f) dir (%8.4f %8.4f %8.4f) range (%8.4f %8.4f)\n",
             origin[0],
             origin[1],
             origin[2],
             dir[0],
             dir[1],
             dir[2],
             range[0],
             range[1]);

  if (0) {  // debug mode, doing a global search
    for (size_t iface = 0; iface < srfc->N_triangles; ++iface) {
      double    this_distance;
      const int lhit = native_trace_ray (srfc, iface, origin, dir, range, &this_distance);
      if (ldebug)
        fprintf (stderr, "face %6lu hit %i in distance %g current_min_distance %g \n", iface, lhit, this_distance, *distance);
      if (lhit && (this_distance < *distance)) {
        *distance     = this_distance;
        *primitive_id = iface;
      }
    }
  } else {
    size_t* start = srfc->native_tracer->cell2map[cellid];
    size_t* end   = srfc->native_tracer->cell2map[cellid + 1];
    if (ldebug)
      fprintf (stderr,
               "kij (%4lu %4lu %4lu) 3D? %i cellid %6lu => start/end %p %p\n",
               kij[0],
               kij[1],
               kij[2],
               srfc->native_tracer->layer_is_3D[k],
               cellid,
               start,
               end);
    for (size_t* ptr = start; ptr < end; ++ptr) {
      const size_t iface = *ptr;
      double       this_distance;
      const int    lhit = native_trace_ray (srfc, iface, origin, dir, range, &this_distance);
      if (ldebug)
        fprintf (stderr, "face %6lu hit %i in distance %8.4f current_min_distance %8.4f \n", iface, lhit, this_distance, *distance);
      if (lhit && (this_distance < *distance)) {
        *distance     = this_distance;
        *primitive_id = iface;
      }
    }
  }

  if (*primitive_id == -1) {
    *distance = POS_INFINITY;
  }
  return 0;
}

int get_triangle_albedo (const t_triangular_surface* srfc, const size_t primitive_id, double* triangle_albedo) {
  const int ldebug = 0;  // enable to see a lot of debug output
  if (primitive_id >= srfc->N_triangles) {
    CHKERROUT (primitive_id, "primitive_id index out of bounds");
  }
  const size_t material_id = srfc->material_of_triangle[primitive_id];
  *triangle_albedo         = srfc->material_albedo[material_id];
  if (ldebug)
    fprintf (stderr, "face %12lu has material id %4lu and albedo %8.4f \n", primitive_id, material_id, *triangle_albedo);
  return 0;
}

int get_global_vertices_bounding_box (const t_triangular_surface* srfc, double min[3], double max[3]) {
  for (size_t k = 0; k < 3; ++k) {
    min[k] = POS_INFINITY;
    max[k] = -POS_INFINITY;
  }
  for (size_t i = 0; i < srfc->N_vertices; ++i) {
    for (size_t k = 0; k < 3; ++k) {
      min[k] = MIN (min[k], srfc->vertices[i][k]);
      max[k] = MAX (max[k], srfc->vertices[i][k]);
    }
  }
  return 0;
}

int get_triangle_bounding_box (const double A[3], const double B[3], const double C[3], double min[3], double max[3]) {
  for (size_t k = 0; k < 3; ++k) {
    min[k] = MIN (MIN (A[k], B[k]), C[k]);
    max[k] = MAX (MAX (A[k], B[k]), C[k]);
  }
  return 0;
}

int get_triangle_normal (const t_triangular_surface* srfc, const size_t primitive_id, double triangle_normal[3]) {
  const int ldebug = 0;
  if (!(primitive_id < srfc->N_triangles)) {
    CHKERROUT (primitive_id, "primitive_id index out of bounds");
  }
  const double* A = srfc->vertices[srfc->triangles[primitive_id][0]];
  const double* B = srfc->vertices[srfc->triangles[primitive_id][1]];
  const double* C = srfc->vertices[srfc->triangles[primitive_id][2]];

  const double u[3] = {B[0] - A[0], B[1] - A[1], B[2] - A[2]};

  const double v[3] = {C[0] - A[0], C[1] - A[1], C[2] - A[2]};

  v_cross_product (u, v, triangle_normal);

  const double norm_t = vec3_norm2 (triangle_normal);
  int          ierr   = vec3_scale (1. / norm_t, triangle_normal);
  CHKERR (ierr);  // normalize new vector

  if (ldebug) {
    fprintf (stderr,
             "triangle id %lu\n"
             "normal A (%7.2f %7.2f %7.2f)\n"
             "normal B (%7.2f %7.2f %7.2f)\n"
             "normal C (%7.2f %7.2f %7.2f)\n"
             "normal u (%7.4f %7.4f %7.4f)\n"
             "normal v (%7.4f %7.4f %7.4f)\n"
             "normal n (%7.4f %7.4f %7.4f) -- norm_t %7.4f\n",
             primitive_id,
             A[0],
             A[1],
             A[2],
             B[0],
             B[1],
             B[2],
             C[0],
             C[1],
             C[2],
             u[0],
             u[1],
             u[2],
             v[0],
             v[1],
             v[2],
             triangle_normal[0],
             triangle_normal[1],
             triangle_normal[2],
             norm_t);
  }
  return 0;
}

/* Watertight ray/triangle intersection
   Sven Woop, Carsten Benthin, and Ingo Wald. Journal of Computer Graphics Techniques (JCGT), 2(1):65–82, 2013.
   -> triangle intersection code from http://jcgt.org/published/0002/01/05/
   with additional remarks from: https://www.allthingsphi.com/blog/2016/12/18/watertight-ray-triangle-intersection.html
*/
int triangle_intersection_woop (const double origin[3],  // ray origin
                                const double dir[3],     // ray direction, resulting distance is measured in norm2(dir) units
                                const double vert_A[3],  // first vertex cartesian coordinates
                                const double vert_B[3],  // second vertex
                                const double vert_C[3],  // third vertex
                                const double tnear,      // will test intersection between
                                const double tfar,       //   distance [tnear, tfar]
                                double*      distance,   // distance to intersection in units of [norm2(dir)]
                                double       hit[3])           // cartesian coordinates of intersection
{
  *distance = NOT_A_NUMBER;
  hit[0]    = NOT_A_NUMBER;
  hit[1]    = NOT_A_NUMBER;
  hit[2]    = NOT_A_NUMBER;

  int ierr;

  /* calculate dimension where the ray direction is maximal */
  int kz;
  {
    double fabs_dir[3];
    vec_fabs (3, dir, fabs_dir);
    kz = max_loc (fabs_dir, 0, 2, &ierr);
  }
  int kx = kz + 1;
  if (kx == 3)
    kx = 0;
  int ky = kx + 1;
  if (ky == 3)
    ky = 0;

  /* swap kx and ky dimension to preserve winding direction of triangles */
  if (dir[kz] < 0.0) {
    int tmp = kx;
    kx      = ky;
    ky      = tmp;
  }

  /* calculate shear constants */
  const double Sx = dir[kx] / dir[kz];
  const double Sy = dir[ky] / dir[kz];
  const double Sz = 1.0 / dir[kz];

  /* calculate vertices relative to ray origin */
  double A[3];
  v_sub (vert_A, origin, A);
  double B[3];
  v_sub (vert_B, origin, B);
  double C[3];
  v_sub (vert_C, origin, C);

  /* perform shear and scaleof vertices */
  const double Ax = A[kx] - Sx * A[kz];
  const double Ay = A[ky] - Sy * A[kz];
  const double Bx = B[kx] - Sx * B[kz];
  const double By = B[ky] - Sy * B[kz];
  const double Cx = C[kx] - Sx * C[kz];
  const double Cy = C[ky] - Sy * C[kz];

  /* calculate scaled barycentric coordinates */
  const double U = Cx * By - Cy * Bx;
  const double V = Ax * Cy - Ay * Cx;
  const double W = Bx * Ay - By * Ax;

  /* Perform edge tests. Moving this test before and at the end of the previous conditional gives higher performance. */
#ifdef BACKFACE_CULLING
  if (U < 0 || V < 0 || W < 0)
    return 0;
#else
  if ((U < 0 || V < 0 || W < 0) && (U > 0 || V > 0 || W > 0))
    return 0;
#endif

  /* calculate determinant */
  const double det = U + V + W;
  if (det == 0)
    return 0;
  const double rcpDet = 1.0 / det;

  /* Calculate scale dz−coordinates of vertices and use them to calculate the hit distance. */
  const double Az = Sz * A[kz];
  const double Bz = Sz * B[kz];
  const double Cz = Sz * C[kz];
  const double T  = U * Az + V * Bz + W * Cz;

  *distance = T * rcpDet;

  if (*distance < tnear)
    return 0;
  if (*distance >= tfar)
    return 0;

  vec_AXPBYPCZ (3, 1, origin, *distance, dir, 0, hit);
  return 1;
}
