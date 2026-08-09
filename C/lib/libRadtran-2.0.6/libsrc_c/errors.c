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
#include "errors.h"

/***********************************************************************************/
/* Functions: err_out                                                     @62_30i@ */
/* Description:                                                                    */
/*   short cuts error output and return status                                     */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int err_out (const char* output, int status) {
  fprintf (stderr, output, status);
  return status;
}

/***********************************************************************************/
/* Functions: mem_err_out                                                 @62_30i@ */
/* Description:                                                                    */
/*   short cuts memory allocation error output and return status                   */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int mem_err_out (const char* output, int line, const char* fff, const char* file) {
  fprintf (stderr, "Error when trying to allocate %s in (line %d, function '%s' in '%s')\n", output, line, fff, file);
  return -1;
}

/***********************************************************************************/
/* Functions: fct_err_out                                                 @62_30i@ */
/* Description:                                                                    */
/*   short cuts error output and return status.                                    */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int fct_err_out (int status, const char* function_name, int line, const char* func, const char* filename) {
  fprintf (stderr,
           "Error %d returned by function %s in (line %d, function '%s', file '%s')\n",
           status,
           function_name,
           line,
           func,
           filename);
  return status;
}

/***********************************************************************************/
/* Function: err_file_out                                                 @62_30i@ */
/* Description:                                                                    */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: RPB                                                                     */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int err_file_out (int status, const char* filename, int line, const char* func, const char* cfilename) {
  fprintf (stderr,
           "Error %d returned when reading file %s in (line %d, function '%s', file '%s')\n",
           status,
           filename,
           line,
           func,
           cfilename);
  return status;
}

/***********************************************************************************/
/* Functions: chkerr_out                                                 @62_30i@ */
/* Description:                                                                    */
/*   short cuts error output and return status.                                    */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Fabian Jakub                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int chkerr (int status, int line, const char* func, const char* filename) {
  if (status) {
    fprintf (stderr, "Error %d in (function '%s', file '%s', line %d)\n", status, func, filename, line);
  }
  return status;
}

int chkerr_out (int status, const char* description, int line, const char* func, const char* filename) {
  if (status) {
    fprintf (stderr,
             "Error %d in (function '%s', file '%s', line %d) with reason: %s\n",
             status,
             func,
             filename,
             line,
             description);
  }
  return status;
}

int chkpointer (void* ptr, int line, const char* func, const char* filename) {
  if (!ptr) {
    fprintf (stderr, "Found NULL pointer %p in (function '%s', file '%s', line %d)\n", ptr, func, filename, line);
    return 1;
  }
  return 0;
}

int chkpointer_out (void* ptr, const char* description, int line, const char* func, const char* filename) {
  if (!ptr) {
    fprintf (stderr,
             "Found NULL pointer %p in (function '%s', file '%s', line %d) with reason: %s\n",
             ptr,
             func,
             filename,
             line,
             description);
    return 1;
  }
  return 0;
}

void* valid_ptr (void* ptr, int line, const char* func, const char* filename) {
  if (!ptr) {
    fprintf (stderr, "Found NULL pointer %p in (function '%s', file '%s', line %d)\n", ptr, func, filename, line);
    return NULL;
  }
  return ptr;
}

void* valid_ptr_out (void* ptr, const char* description, int line, const char* func, const char* filename) {
  if (!ptr) {
    fprintf (stderr,
             "Found NULL pointer %p in (function '%s', file '%s', line %d) with reason: %s\n",
             ptr,
             func,
             filename,
             line,
             description);
    return NULL;
  }
  return ptr;
}

// print a deprecation warning
void deprecate (int year, int line, const char* func, const char* filename) {
  fprintf (stdout,
           "WARNING, this part of code has been marked deprecated as of %d and will be deleted in a future release. If you want to "
           "keep it, please reach out to the maintainers. (function '%s', file '%s', line %d)\n",
           year,
           func,
           filename,
           line);
}
