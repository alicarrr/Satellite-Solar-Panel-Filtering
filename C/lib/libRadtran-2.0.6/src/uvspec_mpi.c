#ifdef HAVE_TENSTREAM
#include <mpi.h>
#include "petscsys.h"
#include "tenstream.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#include "errors.h"

#ifdef HAVE_TENSTREAM
static char help[] = "This is the C wrapper interface to the Tenstream solver environment.\n\n";
#endif

int calluvspec (char* outfilename);

int file_exist (char* filename) {
  FILE* fh = fopen (filename, "r");
  if (!fh) {
    return 0;
  } else {
    fclose (fh);
  }
  return 1;
}

int main (int argc, char** argv) {
  const int ldebug = 0;
#ifdef HAVE_TENSTREAM
  extern FILE* yyin;
  char         outfilename[FILENAME_MAX] = "";
#endif
  int  myid;
  char infilename[FILENAME_MAX] = "";
  int  ierr                     = 0;

  int have_inputfile = 0;

#ifdef HAVE_TENSTREAM
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &myid);
#else
  myid = 0;
#endif

  if (argc >= 1)
    have_inputfile = 1;

  if (!have_inputfile) {
    if (myid == 0)
      CHKERROUT (1, "We need a Libradtran Input file to work with, please specify one with -f <filename>");
  }

  sprintf (infilename, "%s", argv[1]);
  if (file_exist (infilename)) {
    if (myid == 0 && ldebug)
      fprintf (stderr, "Using File %s as libradtran input\n", infilename);
  } else {
    if (myid == 0)
      fprintf (stderr, "LibRadtran Input File %s does not exist\n", infilename);
    CHKERR (2);
  }

#ifdef HAVE_TENSTREAM
  PetscInitialize (&argc, &argv, (char*) 0, help);
  PetscInitializeFortran();

  if (myid == 0) {
    sprintf (outfilename, "std.out");

    yyin = fopen (infilename, "r");
    ierr = calluvspec (outfilename);

    if (ierr)
      fprintf (stderr, "Error, calluvspec returned %d\n", ierr);

    fclose (yyin);

    int func_index = MPI_FUNC_TENSTREAM_FINALIZE;
    MPI_Bcast (&func_index, 1, MPI_INT, 0, MPI_COMM_WORLD);

  } else {
    ierr = tenstream_slave_loop();
  }

  PetscFinalize();
  MPI_Finalize();
#else
  CHKERROUT (-1, "It seems your installation does currently not support the Tenstream solver... check your installation...");
#endif
  return ierr;
}
