#include <stdio.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "legendrePoly.h"
#include "residual.h"
#include "solver.h"
#include "burgers.h"





//-------------------------------------------------
// Evaluate BC function
//
int
boundaryCondition(const double time, double *UB)
{ 
  // Just a placeholder for now
  UB[0] =  1.0;
  UB[1] = -1.0;

  return 0;
}

//-------------------------------------------------
// Set initial state
//
int
InitialCondition(const int N, double *U)
{ 
  // Just a placeholder for now
  for( int ii=0; ii<N; ii++ ) U[ii] = 0.0;
  return 0;
}


//-------------------------------------------------
// Main routine.  Reads command line inputs
// and calls appropriate solver.
//
int
main( int argc, const char * argv[] )
{
  int ierr, nread, N, NiterMax=50;
  double nu;
  burgers B, *pB;

  pB = &B;

  // Grab command line arguments
  if( argc != 3 ){
    printf("Usage: burgers N nu\n");
    printf("  N = number of modes\n");
    printf("  nu = viscosity\n");
    printf("\n"); fflush(stdout);
    return -1;
  }

  nread = sscanf(argv[1], "%d", &N);
  if( nread != 1 ){
    printf("Usage: burgers N nu\n");
    printf("  N = number of modes\n");
    printf("  nu = viscosity\n");
    printf("\n"); fflush(stdout);
    return -1;
  }

  nread = sscanf(argv[2], "%lf", &nu);
  if( nread != 1 ){
    printf("Usage: burgers N nu\n");
    printf("  N = number of modes\n");
    printf("  nu = viscosity\n");
    printf("\n"); fflush(stdout);
    return -1;
  }

  if( (N<1) || (nu<0.0) ){
    printf("Error: N>=1 and nu>=0.0 are required!\n"); fflush(stdout);
  }

  // Allocate and initialize burgers struct
  pB->N = N;
  pB->nu = nu;

  pB->U = gsl_vector_calloc( (size_t)N );
  pB->R = gsl_vector_calloc( (size_t)N );

  // Call solver (only steady Newton currently available)
  ierr = steadyNewton(NiterMax, pB);
  if( ierr != 0 ) return ierr;

  // Write solution to screen         
  FILE *fp = fopen ("test.dat", "w");
  gsl_vector_fprintf(fp, pB->U, "%.15E");
  fclose(fp);

  // Clean up memory
  gsl_vector_free(pB->U);
  gsl_vector_free(pB->R);

  printf("\n");

  return 0;
}
