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

/************************************************************************/
/* GLmystic.h                                                           */
/*                                                                      */
/* GLmystic - 3D graphics support for MYSTIC                            */
/*                                                                      */
/* Author: Bernhard Mayer                                               */
/*                                                                      */
/************************************************************************/

#ifndef __GLmystic_h
#define __GLmystic_h

#if defined(__cplusplus)
extern "C" {
#endif

#include <sys/shm.h> /* inter process communication */
#include "mystic.h"

/* legacy OpenGL with glBegin() and glEnd() is deprecated */
#define GL_SILENCE_DEPRECATION

#define GLMYSTIC_PHOTONS 5000
#define GLMYSTIC_ELEVATION 6000
#define GLMYSTIC_CLOUD 7000
#define GLMYSTIC_RESULT 8000

#if HAVE_OPENGL
#if HAVE_WINDOWS_H && defined(_WIN32)
#include <windows.h>
#endif

#ifdef HAVE_GL_GL_H
#include <GL/gl.h>
#elif defined(HAVE_OPENGL_GL_H)
#include <OpenGL/gl.h>
#endif

#ifdef HAVE_GL_GLU_H
#include <GL/glu.h>
#elif defined(HAVE_OPENGL_GLU_H)
#include <OpenGL/glu.h>
#endif

#ifdef HAVE_GL_GLUT_H
#include <GL/glut.h>
#elif defined(HAVE_GLUT_GLUT_H)
#include <GLUT/glut.h>
#endif

/* inter process communication */
#include <sys/shm.h>
#include <sys/wait.h>
#endif

void GLmystic_init (void);

void GLmystic_free_shared (key_t key);
void GLmystic_calloc_shared (key_t key, int size);

void GLmystic_setdomain (key_t key, double xmin, double ymin, double zmin, double xmax, double ymax, double zmax);

void GLmystic_write_elevation (double** Z, int Nx, int Ny);

float* GL_access_float (key_t key, int mode);
void   GL_release_float (float* px);

int* GL_access_int (key_t key, int mode);
void GL_release_int (int* px);

void* GL_access (key_t key, int mode);
void  GL_release (void* px);

void GLmystic_start_photon_path (photon_struct* p);
void GLmystic_add_node_to_photon_path (double x, double y, double z);
void GLmystic_add_periodic_nodes_to_photon_path (photon_struct*     p,
                                                 atmosphere_struct* atmos,
                                                 double             x,
                                                 double             y,
                                                 double             z,
                                                 double             step);
void GLmystic_write_global (photon_struct* p, sample_struct* sample, result_struct* result);
void GLmystic_write_caoth (atmosphere_struct* atmos);

#endif
