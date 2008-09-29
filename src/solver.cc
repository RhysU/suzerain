#include <stdio.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "legendrePoly.h"
#include "residual.h"
#include "burgers.h"

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
  double UB[2];

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

  // Get BCs (assumed constant)
  ierr = boundaryCondition(0.0, UB);
  if( ierr != 0 ) return ierr;

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
// March in time with 4th order Runge-Kutta
//
int
unsteadyRK4(const int Nstep, burgers *pB)
{  
  int ierr;
  const int N = pB->N;
  double time=0.0, tmptime=0.0;
  double dt = 1.0/(N+1); // hardcoded time step for now
  double UB[2], UB0[2], UB1[2];

  gsl_vector *U = pB->U;
  gsl_vector *Utmp;
  gsl_vector *R = pB->R;
  gsl_vector *K1, *K2, *K3, *K4;
  gsl_vector *dU;

  gsl_matrix *iMM;

  // Allocate storage
  dU   = gsl_vector_calloc((size_t)N);
  Utmp = gsl_vector_calloc((size_t)N);
  K1   = gsl_vector_calloc((size_t)N);
  K2   = gsl_vector_calloc((size_t)N);
  K3   = gsl_vector_calloc((size_t)N);
  K4   = gsl_vector_calloc((size_t)N);

  // Compute inverse mass matrix
  iMM = gsl_matrix_alloc((size_t)N, (size_t)N);
  ierr = legendre0InverseMassMatrix(N, iMM->data);
  if( ierr != 0 ) return ierr;
  


  for( int istep=0; istep<Nstep; istep++ ){

    
    //---------------------------------------
    // STAGE 1
    //---------------------------------------

    // Compute/store stage 1 state, time, and BCs
    gsl_vector_memcpy(Utmp, U);
    tmptime = time;

    ierr = boundaryCondition(tmptime, UB);
    if( ierr != 0 ) return ierr;

    // Set residual and Jacobian to zero
    gsl_vector_set_zero(R);

    // Evaluate residual
    ierr = interiorResidual(N, pB->nu, Utmp->data, UB, R->data, (double *)NULL);
    if( ierr != 0 ) return ierr;
    
    // K1 = -MM\R;
    gsl_blas_dgemv(CblasNoTrans, -1.0, iMM, R, 0.0, K1);


    //---------------------------------------
    // STAGE 2
    //---------------------------------------

    // Compute/store stage 2 state, time, and BCs
    gsl_vector_memcpy(Utmp, U);
    gsl_blas_daxpy(0.5*dt, K1, Utmp);

    tmptime = time + 0.5*dt;

    ierr = boundaryCondition(tmptime, UB);
    if( ierr != 0 ) return ierr;

    // Set residual and Jacobian to zero
    gsl_vector_set_zero(R);

    // Evaluate residual
    ierr = interiorResidual(N, pB->nu, Utmp->data, UB, R->data, (double *)NULL);
    if( ierr != 0 ) return ierr;
    
    // K2 = -MM\R;
    gsl_blas_dgemv(CblasNoTrans, -1.0, iMM, R, 0.0, K2);


    //---------------------------------------
    // STAGE 3
    //---------------------------------------

    // Compute/store stage 2 state, time, and BCs
    gsl_vector_memcpy(Utmp, U);
    gsl_blas_daxpy(0.5*dt, K2, Utmp);

    tmptime = time + 0.5*dt;

    ierr = boundaryCondition(tmptime, UB);
    if( ierr != 0 ) return ierr;

    // Set residual and Jacobian to zero
    gsl_vector_set_zero(R);

    // Evaluate residual
    ierr = interiorResidual(N, pB->nu, Utmp->data, UB, R->data, (double *)NULL);
    if( ierr != 0 ) return ierr;
    
    // K3 = -MM\R;
    gsl_blas_dgemv(CblasNoTrans, -1.0, iMM, R, 0.0, K3);


    //---------------------------------------
    // STAGE 4
    //---------------------------------------

    // Compute/store stage 2 state, time, and BCs
    gsl_vector_memcpy(Utmp, U);
    gsl_blas_daxpy(dt, K3, Utmp);

    tmptime = time + dt;

    ierr = boundaryCondition(tmptime, UB);
    if( ierr != 0 ) return ierr;


    // Set residual and Jacobian to zero
    gsl_vector_set_zero(R);

    // Evaluate residual
    ierr = interiorResidual(N, pB->nu, Utmp->data, UB, R->data, (double *)NULL);
    if( ierr != 0 ) return ierr;
    
    // K4 = -MM\R;
    gsl_blas_dgemv(CblasNoTrans, -1.0, iMM, R, 0.0, K4);


    //---------------------------------------
    // UNSTEADY BC TERM
    //---------------------------------------

    // Evaluate boundary conditions at endpoints of interval
    ierr = boundaryCondition(time, UB0);
    if( ierr != 0 ) return ierr;

    ierr = boundaryCondition(time+dt, UB1);
    if( ierr != 0 ) return ierr;

    // Set residual and Jacobian to zero
    gsl_vector_set_zero(R);

    ierr = unsteadyBoundaryResidual(N, UB0, UB1, R->data);
    if( ierr != 0 ) return ierr;


    //---------------------------------------
    // SUM CONTRIBUTIONS TO UPDATE
    //---------------------------------------

    // Unsteady BC piece
    gsl_blas_dgemv(CblasNoTrans, -1.0, iMM, R, 0.0, dU);

    // Standard residual piece
    gsl_blas_daxpy(dt/6.0, K1, dU);
    gsl_blas_daxpy(dt/3.0, K2, dU);
    gsl_blas_daxpy(dt/3.0, K3, dU);
    gsl_blas_daxpy(dt/6.0, K4, dU);


    //---------------------------------------
    // UPDATE STATE AND TIME
    //---------------------------------------
    gsl_blas_daxpy(1.0, dU, U);
    time += dt;
    
  }

  return 0;
}
