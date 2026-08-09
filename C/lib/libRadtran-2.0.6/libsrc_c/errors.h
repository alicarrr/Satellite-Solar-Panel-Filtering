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

#ifndef __errors_h
#define __errors_h

#define ERROR_POSITION __LINE__, (const char*) __func__, __FILE__

#define DEPRECATE(year) deprecate (year, __LINE__, __func__, __FILE__)

#define CHKERR(int_status)                                                                                                         \
  if (int_status) {                                                                                                                \
    return chkerr (int_status, __LINE__, __func__, __FILE__);                                                                      \
  }

#define CHKERROUT(int_status, description)                                                                                         \
  if (int_status) {                                                                                                                \
    return chkerr_out (int_status, description, __LINE__, __func__, __FILE__);                                                     \
  }

#define CHKERRPTR(int_status)                                                                                                      \
  if (int_status) {                                                                                                                \
    chkerr (int_status, __LINE__, __func__, __FILE__);                                                                             \
    return NULL;                                                                                                                   \
  }

#define CHKERROUTPTR(int_status, description)                                                                                      \
  if (int_status) {                                                                                                                \
    chkerr_out (int_status, description, __LINE__, __func__, __FILE__);                                                            \
    return NULL;                                                                                                                   \
  }

#define CHKPOINTER(pointer)                                                                                                        \
  if (!pointer) {                                                                                                                  \
    return chkpointer (pointer, __LINE__, __func__, __FILE__);                                                                     \
  }

#define CHKPOINTEROUT(pointer, description)                                                                                        \
  if (!pointer) {                                                                                                                  \
    return chkpointer_out (pointer, description, __LINE__, __func__, __FILE__);                                                    \
  }

#define VALIDPTR(pointer)                                                                                                          \
  if (!pointer) {                                                                                                                  \
    return valid_ptr (pointer, __LINE__, __func__, __FILE__);                                                                      \
  }

#define VALIDPTROUT(pointer, description)                                                                                          \
  if (!pointer) {                                                                                                                  \
    return valid_ptr_out (pointer, description, __LINE__, __func__, __FILE__);                                                     \
  }

#if defined(__cplusplus)
extern "C" {
#endif

int err_out (const char* output, int status);

int mem_err_out (const char* output, int line, const char* fff, const char* file);

int fct_err_out (int status, const char* function_name, int line, const char* func, const char* filename);

int err_file_out (int status, const char* filename, int line, const char* func, const char* cfilename);

int chkerr (int status, int line, const char* func, const char* filename);
int chkerr_out (int status, const char* description, int line, const char* func, const char* filename);

int chkpointer (void* ptr, int line, const char* func, const char* filename);
int chkpointer_out (void* ptr, const char* description, int line, const char* func, const char* filename);

void* valid_ptr (void* ptr, int line, const char* func, const char* filename);
void* valid_ptr_out (void* ptr, const char* description, int line, const char* func, const char* filename);

void deprecate (int year, int line, const char* func, const char* filename);

#if defined(__cplusplus)
}
#endif

#endif /* _ERRORS_H */
