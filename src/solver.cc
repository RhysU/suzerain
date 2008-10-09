#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "legendrePoly.h"
#include "residual.h"
#include "tools.h"
#include "burgers.h"

#define VERBOSE 1

//-------------------------------------------------
// Solve steady Burgers with Newton-Raphson
//
int
steadyNewton(const int NiterMax, burgersSteady *pB)
{
  int ierr, iiter;
  const int Nmode = pB->Nmode;
  double Rnorm=0.0;
  const double RTOL = 1e-12;
  double UB[2] = {pB->UB[0], pB->UB[1]};

  gsl_vector *U = pB->U;
  gsl_vector *R = pB->R;
  gsl_vector *dU;
  gsl_matrix *R_U;

  quadBasis *pQB;

  int s;
  gsl_permutation *p;

  // Allocate storage
  dU = gsl_vector_calloc((size_t)Nmode);
  R_U = gsl_matrix_calloc((size_t)Nmode, (size_t)Nmode);
  p = gsl_permutation_alloc((size_t)Nmode);

  // Set residual to zero
  gsl_vector_set_zero(R);

  // evaluate quad points and weights
  ierr = evaluateQuadratureAndBasisForResidual(Nmode, &pQB);
  if( ierr != 0 ) return ierr;

  // Evaluate residual
  ierr = interiorResidual(Nmode, pB->nu, U->data, UB, pQB, pB->RABFlag, pB->kappa, R->data, R_U->data);
  if( ierr != 0 ) return ierr;

  // Compute residual 2-norm
  Rnorm = gsl_blas_dnrm2(R);

  iiter = 0;

#if VERBOSE==0
  printf("Iteration \t Residual Norm\n");
  printf("--------------------------\n");
  printf("%d \t %.12E\n", iiter, Rnorm);
#endif

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
    ierr = interiorResidual(Nmode, pB->nu, U->data, UB, pQB, pB->RABFlag, pB->kappa, R->data, R_U->data);
    if( ierr != 0 ) return ierr;

    // Compute residual 2-norm
    Rnorm = gsl_blas_dnrm2(R);

    iiter++;
    
#if VERBOSE==0    
    printf("%d \t %.12E\n", iiter, Rnorm);
#endif
  }

  if( Rnorm > RTOL ){
    printf("WARNING: Solution only converged to R=%.6E in %d iterations\n", Rnorm, iiter);
    fflush(stdout);
  }

  // Clean up
  gsl_vector_free(dU);
  gsl_matrix_free(R_U);
  gsl_permutation_free(p);  
  freeQuadratureAndBasisForResidual(pQB);
  free(pQB);


  return 0;
}

//------------------------------------------------------
// Solve steady Burgers with Newton-Raphson
// Quad points, weights, and basis must be pre-computed
//
int
steadyNewtonPreComputeBasis(const int NiterMax, burgersSteady *pB, quadBasis *pQB)
{
  int ierr, iiter;
  const int Nmode = pB->Nmode;
  double Rnorm=0.0;
  const double RTOL = 1e-12;
  double UB[2] = {pB->UB[0], pB->UB[1]};

  gsl_vector *U = pB->U;
  gsl_vector *R = pB->R;
  gsl_vector *dU;
  gsl_matrix *R_U;

  int s;
  gsl_permutation *p;

  // Allocate storage
  dU = gsl_vector_calloc((size_t)Nmode);
  R_U = gsl_matrix_calloc((size_t)Nmode, (size_t)Nmode);
  p = gsl_permutation_alloc((size_t)Nmode);

  // Set residual to zero
  gsl_vector_set_zero(R);

  // Evaluate residual
  ierr = interiorResidual(Nmode, pB->nu, U->data, UB, pQB, pB->RABFlag, pB->kappa, R->data, R_U->data);
  if( ierr != 0 ) return ierr;

  // Compute residual 2-norm
  Rnorm = gsl_blas_dnrm2(R);

  iiter = 0;

#if VERBOSE==0
  printf("Iteration \t Residual Norm\n");
  printf("--------------------------\n");
  printf("%d \t %.12E\n", iiter, Rnorm);
#endif

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
    ierr = interiorResidual(Nmode, pB->nu, U->data, UB, pQB, pB->RABFlag, pB->kappa, R->data, R_U->data);
    if( ierr != 0 ) return ierr;

    // Compute residual 2-norm
    Rnorm = gsl_blas_dnrm2(R);

    iiter++;
    
#if VERBOSE==0    
    printf("%d \t %.12E\n", iiter, Rnorm);
#endif
  }

  if( Rnorm > RTOL ){
    printf("WARNING: Solution only converged to R=%.6E in %d iterations\n", Rnorm, iiter);
    fflush(stdout);
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
unsteadyRK4(burgersUnsteady *pB)
{  
  int ierr, Nstep = pB->Nstep, Nwrite = pB->Nwrite, Nstat = pB->Nstat;
  const int Nmode = pB->Nmode;
  double time=pB->time, tmptime=time;
  double UB[2], UB0[2], UBmid[2], UB1[2], UBavg[2]={0.0, 0.0};
  double Rnorm;
  double dt = pB->dt;

  bool FirstAvg = true;
  double tstart, tend;

  gsl_vector *U = pB->U;
  gsl_vector *Uavg = pB->Uavg;
  gsl_vector *Upre, *Utmp, *Uwri, *dU;
  gsl_vector *R;
  gsl_vector *K1, *K2, *K3, *K4;

  gsl_matrix *iMM;

  quadBasis *pQB;

  FILE *fp;

  // Allocate storage
  dU   = gsl_vector_calloc((size_t)Nmode);
  Utmp = gsl_vector_calloc((size_t)Nmode);
  Upre = gsl_vector_calloc((size_t)Nmode);
  Uwri = gsl_vector_calloc((size_t)Nmode);
  R    = gsl_vector_calloc((size_t)Nmode);
  K1   = gsl_vector_calloc((size_t)Nmode);
  K2   = gsl_vector_calloc((size_t)Nmode);
  K3   = gsl_vector_calloc((size_t)Nmode);
  K4   = gsl_vector_calloc((size_t)Nmode);
  iMM  = gsl_matrix_alloc((size_t)Nmode, (size_t)Nmode);

  // evaluate quad points and weights
  ierr = evaluateQuadratureAndBasisForResidual(Nmode, &pQB);
  if( ierr != 0 ) return ierr;

  // Compute inverse mass matrix
  ierr = legendre0InverseMassMatrix(Nmode, iMM->data);
  if( ierr != 0 ) return ierr;

  // write time
  printf("Writing solution at time = %.6E\n", time); fflush(stdout);
  fp = fopen("time.dat", "w");
  fprintf(fp, "%.15E\n", time);
  fclose(fp);

  // write solution
  fp = fopen("soln.dat", "w");
  gsl_vector_fprintf(fp, pB->U, "%.15E");
  fclose(fp);


  // write bcs
  fp = fopen("bcs.dat", "w");
  fprintf(fp, "%.15E %.15E\n", pB->UB[0], pB->UB[1]);
  fclose(fp);

  // zero out Uavg in preparation for averaging
  gsl_vector_set_zero(Uavg);

  // RK4
  for( int istep=0; istep<Nstep; istep++ ){
    
    // save current U (used for averaging)
    gsl_vector_memcpy(Upre, U);


    //---------------------------------------
    // STAGE 1
    //---------------------------------------

    // Compute/store stage 1 state, time, and BCs
    gsl_vector_memcpy(Utmp, U);
    tmptime = time;

    if( istep == 0 ){
      UB0[0] = pB->UB[0];
      UB0[1] = pB->UB[1];
    } else{
      UB0[0] = UB1[0];
      UB0[1] = UB1[1];
    }

    // Set residual and Jacobian to zero
    gsl_vector_set_zero(R);

    // Evaluate residual
    ierr = interiorResidual(Nmode, pB->nu, Utmp->data, UB0, pQB, false, -1.0, R->data, (double *)NULL);
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

    ierr = boundaryCondition(pB->BCname, time, tmptime, UB0, UBmid);
    if( ierr != 0 ) return ierr;

    // Set residual and Jacobian to zero
    gsl_vector_set_zero(R);

    // Evaluate residual
    ierr = interiorResidual(Nmode, pB->nu, Utmp->data, UBmid, pQB, false, -1.0, R->data, (double *)NULL);
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

//     ierr = boundaryCondition(pB->BCname, tmptime, UBmid);
//     if( ierr != 0 ) return ierr;

    // Set residual and Jacobian to zero
    gsl_vector_set_zero(R);

    // Evaluate residual
    ierr = interiorResidual(Nmode, pB->nu, Utmp->data, UBmid, pQB, false, -1.0, R->data, (double *)NULL);
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

    ierr = boundaryCondition(pB->BCname, tmptime-0.5*dt, tmptime, UBmid, UB1);
    if( ierr != 0 ) return ierr;


    // Set residual and Jacobian to zero
    gsl_vector_set_zero(R);

    // Evaluate residual
    ierr = interiorResidual(Nmode, pB->nu, Utmp->data, UB1, pQB, false, -1.0, R->data, (double *)NULL);
    if( ierr != 0 ) return ierr;

    // K4 = -MM\R;
    gsl_blas_dgemv(CblasNoTrans, -1.0, iMM, R, 0.0, K4);


    //---------------------------------------
    // UNSTEADY BC TERM
    //---------------------------------------

//     // Evaluate boundary conditions at endpoints of interval
//     ierr = boundaryCondition(pB->BCname, time, UB0);
//     if( ierr != 0 ) return ierr;

//     ierr = boundaryCondition(pB->BCname, time+dt, UB1);
//     if( ierr != 0 ) return ierr;

    // Set residual and Jacobian to zero
    gsl_vector_set_zero(R);

    ierr = unsteadyBoundaryResidual(Nmode, UB0, UB1, pQB, R->data);
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

    // Write solution is requested
    if( !((istep+1)%Nwrite) ){
      printf("Writing solution at time = %.6E\n", time); fflush(stdout);

      // write time
      fp = fopen("time.dat", "a");
      fprintf(fp, "%.15E\n", time);
      fclose(fp);
      
      // write solution
      fp = fopen("soln.dat", "a");
      gsl_vector_fprintf(fp, U, "%.15E");
      fclose(fp);

      // write bcs
      fp = fopen("bcs.dat", "a");
      fprintf(fp, "%.15E %.15E\n", UB1[0], UB1[1]);
      fclose(fp);

    }

    // Update average
    if( istep > Nstat ){
      
      // average state
      if( FirstAvg ) tstart = time-dt;
      FirstAvg = false;
      tend = time;

      gsl_blas_daxpy(0.5, Upre, Uavg);
      gsl_blas_daxpy(0.5, U   , Uavg);

      // average bcs
      UBavg[0] += 0.5*(UB0[0] + UB1[0]);
      UBavg[1] += 0.5*(UB0[1] + UB1[1]);

      // write averaged solution
      if( !((istep+1)%Nwrite) ){
	fp = fopen("running_average_bcs.dat", "a");
	fprintf(fp, "%.15E, %.15E\n", dt*UBavg[0]/(tend-tstart), dt*UBavg[1]/(tend-tstart));
	fclose(fp);

	gsl_vector_memcpy(Uwri, Uavg);
	gsl_vector_scale(Uwri, dt/(tend-tstart));
	
	fp = fopen("running_average_soln.dat", "a");
	gsl_vector_fprintf(fp, Uwri, "%.15E");
	fclose(fp);
      }
      
    }

  }

  // scale averages
  gsl_vector_scale(Uavg, dt/(tend-tstart));
  UBavg[0] *= dt/(tend-tstart);
  UBavg[1] *= dt/(tend-tstart);

  // write averaged solution
  fp = fopen("average_soln.dat", "w");
  fprintf(fp, "%.15E\n", time);
  fprintf(fp, "%.15E %.15E\n", UBavg[0], UBavg[1]);
  fprintf(fp, "%d\n", Nmode);
  gsl_vector_fprintf(fp, Uavg, "%.15E");
  fclose(fp);

  // update for saving
  pB->time = time;
  pB->UB[0] = UB1[0];
  pB->UB[1] = UB1[1];
  pB->UBavg[0] = UBavg[0];
  pB->UBavg[1] = UBavg[1];

  // update bcs for same

  // Clean up
  gsl_vector_free(dU);
  gsl_vector_free(Utmp);
  gsl_vector_free(Upre);
  gsl_vector_free(Uwri);
  gsl_vector_free(R);
  gsl_vector_free(K1);
  gsl_vector_free(K2);
  gsl_vector_free(K3);
  gsl_vector_free(K4);
  gsl_matrix_free(iMM);
  freeQuadratureAndBasisForResidual(pQB);
  free(pQB);

  return 0;
}


// //-------------------------------------------------
// // Computes the eigen-decomposition of the
// // integrating factor matrix MM\KK
// // where MM = mass matrix and 
// // KK = stiffness matrix
// //
// // The eigenvalues are stored in the vector lambda,
// // and the eigenvectors are stored in the matrix RR.
// // The matrix iRR stores inv(RR).
// //
// int
// eigenDecompIntFactorMatrix( const int N, const double *MM, const double *KK,
// 			    double *lambda, double *RR, double *iRR)
// {


//   gsl_matrix *MMtmp, *KKtmp, *RRtmp;
//   gsl_matrix_const_view gslMM = gsl_matrix_const_view_array(MM, N, N);
//   gsl_matrix_const_view gslKK = gsl_matrix_const_view_array(KK, N, N);

//   gsl_vector_view gslLambda = gsl_vector_view_array(lambda, N);
//   gsl_matrix_view gslRR = gsl_matrix_view_array(RR, N, N);
//   gsl_matrix_view gsliRR = gsl_matrix_view_array(iRR, N, N);

//   // Initialize output
//   for( int ii=0; ii<N; ii++ ) lambda[ii] = 0.0;
//   for( int ii=0; ii<N*N; ii++ ) iRR[ii] = RR[ii] = 0.0;
  
//   // copy MM and KK to temporary storage (b/c eigen calc overwrites inputs)
//   MMtmp = gsl_matrix_alloc((size_t)N, (size_t)N);
//   KKtmp = gsl_matrix_alloc((size_t)N, (size_t)N);

//   gsl_matrix_memcpy(MMtmp, &gslMM.matrix);
//   gsl_matrix_memcpy(KKtmp, &gslKK.matrix);

//   // set up and compute eigen-decomposition
//   gsl_eigen_gensymmv_workspace *eig_workspace = gsl_eigen_gensymmv_alloc(N);
//   int ierr = gsl_eigen_gensymmv(KKtmp, MMtmp, &gslLambda.vector, &gslRR.matrix, eig_workspace);
//   if( ierr != 0 ) return ierr;

//   // invert RR
//   RRtmp = gsl_matrix_alloc((size_t)N, (size_t)N);
//   gsl_matrix_memcpy(RRtmp, &gslRR.matrix);
//   int s;
//   gsl_permutation * p = gsl_permutation_alloc(N);
//   gsl_linalg_LU_decomp(RRtmp, p, &s);

//   gsl_vector *b;
//   b = gsl_vector_alloc(N);

//   for( int ii=0; ii<N; ii++ ){
//     gsl_vector_view gsliRR = gsl_vector_view_array_with_stride(&(iRR[ii]), N, N);
//     gsl_vector_set_basis(b, ii);

//     ierr = gsl_linalg_LU_solve (RRtmp, p, b, &gsliRR.vector);
//     if( ierr != 0 ){
//       printf("Cannot solve system!\n"); fflush(stdout);
//       return -1;
//     }
//   }

//   // clean up
//   gsl_eigen_gensymmv_free(eig_workspace);
//   gsl_matrix_free(MMtmp);
//   gsl_matrix_free(KKtmp);

//   return 0;
// }


// //-------------------------------------------------
// // intFactor = RR*exp(nu*t*diag(lambda))*iRR
// //
// int
// integratingFactor( const int N, const double *RR, const double *lambda, const double *iRR, 
// 		   const double nu, const double time,
// 		   double *intFactor)
// {

//   gsl_matrix *mLambda, *tempMatrix;
//   gsl_matrix_const_view gslRR = gsl_matrix_const_view_array(RR, N, N);
//   gsl_matrix_const_view gsliRR = gsl_matrix_const_view_array(iRR, N, N);
//   gsl_matrix_view gslIntFactor = gsl_matrix_view_array(intFactor, N, N);

//   tempMatrix = gsl_matrix_calloc(N,N);

//   // allocate and initialize diagonal eigenvalue matrix
//   mLambda = gsl_matrix_calloc(N,N);
//   for( int ii=0; ii<N; ii++ ) gsl_matrix_set(mLambda, ii, ii, exp(nu*time*lambda[ii]));

//   // matrix multiply: intFactor = RR*mLambda*iRR
//   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &gslRR.matrix, mLambda, 0.0, tempMatrix);
//   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tempMatrix, &gsliRR.matrix, 0.0, &gslIntFactor.matrix);

//   // clean up
//   gsl_matrix_free(mLambda);
//   gsl_matrix_free(tempMatrix);

//   return 0;
// }

