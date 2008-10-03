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

  int s;
  gsl_permutation *p;

  // Allocate storage
  dU = gsl_vector_calloc((size_t)Nmode);
  R_U = gsl_matrix_calloc((size_t)Nmode, (size_t)Nmode);
  p = gsl_permutation_alloc((size_t)Nmode);

  // Set residual to zero
  gsl_vector_set_zero(R);

  // Evaluate residual
  ierr = interiorResidual(Nmode, pB->nu, U->data, UB, R->data, R_U->data);
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
    ierr = interiorResidual(Nmode, pB->nu, U->data, UB, R->data, R_U->data);
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
    ierr = interiorResidual(Nmode, pB->nu, Utmp->data, UB0, R->data, (double *)NULL);
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
    ierr = interiorResidual(Nmode, pB->nu, Utmp->data, UBmid, R->data, (double *)NULL);
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
    ierr = interiorResidual(Nmode, pB->nu, Utmp->data, UBmid, R->data, (double *)NULL);
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
    ierr = interiorResidual(Nmode, pB->nu, Utmp->data, UB1, R->data, (double *)NULL);
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

    ierr = unsteadyBoundaryResidual(Nmode, UB0, UB1, R->data);
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
      fprintf(fp, "%.15E, %.15E\n", UB1[0], UB1[1]);
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
  fprintf(fp, "%.15E, %.15E\n", UBavg[0], UBavg[1]);
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

  return 0;
}
