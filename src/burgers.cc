#include <stdio.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "legendrePoly.h"
#include "residual.h"

// Structure to allow for easy passing/maniuplation of some data
typedef struct
{
  int N; // number of modes
  double nu; // viscosity
  gsl_vector *U; // state vector
  gsl_vector *R; // residual vector

} burgers;

//-------------------------------------------------
// Solve steady Burgers with Newton-Raphson
//
int
steadyNewton(const int NiterMax, burgers *pB)
{
  int ierr, iiter;
  const int N = pB->N;
  double Rnorm=0.0;
  const double RTOL = 1e-12;
  const double UB[2] = {1.0, -1.0};

  gsl_vector *U = pB->U;
  gsl_vector *R = pB->R;
  gsl_vector *dU;
  gsl_matrix *R_U;

  int s;
  gsl_permutation *p;

  // Allocate storage
  dU = gsl_vector_calloc((size_t)N);
  R_U = gsl_matrix_calloc((size_t)N, (size_t)N);
  p = gsl_permutation_alloc((size_t)N);

  // Set residual to zero
  gsl_vector_set_zero(R);

  // Evaluate residual
  ierr = interiorResidual(N, pB->nu, U->data, UB, R->data, R_U->data);
  if( ierr != 0 ) return ierr;

  // Compute residual 2-norm
  Rnorm = gsl_blas_dnrm2(R);

  iiter = 0;

  printf("Iteration \t Residual Norm\n");
  printf("--------------------------\n");
  printf("%d \t %.12E\n", iiter, Rnorm);

  while( (Rnorm>RTOL) && (iiter<NiterMax) ){
    
    // Solve for update: dU = -R_U\R
    gsl_linalg_LU_decomp(R_U, p, &s);
    gsl_linalg_LU_solve(R_U, p, R, dU);
    gsl_vector_scale(dU, -1.0);

    // Update state: U += dU
    gsl_vector_add(U, dU);

    // Set residual and Jacobian to zero
    gsl_vector_set_zero(R);
    gsl_matrix_set_zero(R_U);

    // Evaluate residual
    ierr = interiorResidual(N, pB->nu, U->data, UB, R->data, R_U->data);
    if( ierr != 0 ) return ierr;

    // Compute residual 2-norm
    Rnorm = gsl_blas_dnrm2(R);

    iiter++;
    
    printf("%d \t %.12E\n", iiter, Rnorm);
  }

  // Clean up
  gsl_vector_free(dU);
  gsl_matrix_free(R_U);
  gsl_permutation_free(p);

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
