#include <stdio.h>
#include <math.h>
#include <string.h>

#include "legendrePoly.h"
#include "burgers.h"

//------------------------------------------------
// Evaluates n wave solution to Burgers eqn.
// Use for testing unsteady solver.
//
int 
nWave( const double x, const double t, const double nu, const double c, const double a, const double t0,
       double *u )
{

  // Shift x and t
  const double xx = x - c*t;
  const double tt = t + t0;

  // Some intermediate quantities
  const double rat = xx/tt;
  const double pow = -0.25*xx*xx/(nu*tt);
  const double frt = sqrt(a/tt);

  const double num = rat*frt*exp(pow);
  const double den = 1.0 + frt*exp(pow);

  // Solution
  (*u) = num/den + c;

  return 0;
}

// random number generator is global
const gsl_rng_type * global_bc_T;
gsl_rng * global_bc_r;

int
initializeRandNumGen()
{
  gsl_rng_env_setup();
  
  global_bc_T = gsl_rng_default;
  global_bc_r = gsl_rng_alloc (global_bc_T);
  return 0;
}


//-------------------------------------------------
// Evaluate BC function
//
int
boundaryCondition(const char *fname, const double time0, const double time1, const double UB0[2],
		  double UB1[2])
{ 
  int ierr;
  const double PI = 3.141592653589793e+00;
  
  // pick requested bcs
  if( strcmp(fname, "nwave")==0 ){
    ierr = nWave(-1.0, time1, 1e-2, 1.0, 16.0, 1.0, &(UB1[0]));
    ierr = nWave( 1.0, time1, 1e-2, 1.0, 16.0, 1.0, &(UB1[1]));
  } else if( strcmp(fname, "zeroFunction")==0 ){
    UB1[0] = UB1[1] = 0.0;
  } else if( strcmp(fname, "steadyShock")==0 ){
    UB1[0] =  1.0;
    UB1[1] = -1.0;  
  } else if( strcmp(fname, "steadyBL")==0 ){
    UB1[0] = 1.0;
    UB1[1] = 0.0;
  } else if( strcmp(fname, "sineOscillation")==0 ){
    double ramp = 1e3*time1;
    if( ramp > 1.0 ) ramp = 1.0;
    UB1[0] =  1.0 + ramp*0.5*sin(2.0*2.0*PI*time1);
    UB1[1] = 0.0;
  } else if( strcmp(fname, "Langevin")==0 ){
    if( UB0 == NULL ){
      UB1[0] = 1.0;
      UB1[1] = 0.0;
    } else{
      double TL = 1e-4;
      double dt = time1 - time0;
      double sig2 = 0.1;
      double dU = (UB0[0] - 1.0)*(1.0 - dt/TL) + sqrt(2.0*sig2*dt/TL)*gsl_ran_gaussian(global_bc_r, 1.0);

      UB1[0] =  1.0 + dU;
      UB1[1] = 0.0;
    }
  } else{
    printf("BC fcn %s is not known.  Exiting.\n", fname); fflush(stdout);
    return -1;
  }

  return 0;
}

//-------------------------------------------------
// Set initial state
//
int
initialCondition(const char *fname, const double UB[2], const int N, double *U)
{ 

  // L2 projection of ic function onto solution space
  int ierr;
  const int solnOrder = N-1+2;
  const int quadOrder = 4*solnOrder + 1; // exact for poly of order 3*N 
                                         // (of course, ic may not be poly, but hopefully error is small enough)
  const int Nquad = quadOrder/2 + 1;

  double xq[Nquad], wq[Nquad];
  double phi[N];
  double rhs[N], u, u0;

  gsl_matrix *iMM;

  // Get quadrature points
  ierr = legendreGaussQuad(Nquad, xq, wq);
  if( ierr != 0 ) return ierr;

  // zero rhs
  for( int ii=0; ii<N; ii++ ) rhs[ii] = 0.0;

  // Loop over quadrature points to form right hand side vector
  for( int iquad=0; iquad<Nquad; iquad++ ){
  
    // Evaluate IC
    if( strcmp(fname, "nwave")==0 ){
      ierr = nWave(xq[iquad], 0.0, 1e-2, 1.0, 16.0, 1.0, &u);
      if( ierr != 0 ) return ierr;
    } else if( strcmp(fname, "zeroFunction")==0 ){
      u = 0.0;
    } else{
      printf("Unknown IC function.  Exiting.\n"); fflush(stdout);
      return -1;
    }

    u0 = u - (UB[0]*0.5*(1.0-xq[iquad]) + UB[1]*0.5*(1.0+xq[iquad]));

    // Evaluate basis
    ierr = legendrePolyZero(N, xq[iquad], phi, (double *)NULL);
    if( ierr != 0 ) return ierr;

    // Add to right hand side
    for( int ii=0; ii<N; ii++ ) rhs[ii] += wq[iquad]*phi[ii]*u0;
  }

  // Invert system to compute initial state
  iMM = gsl_matrix_calloc((size_t)N, (size_t)N);
  ierr = legendre0InverseMassMatrix(N, iMM->data);

  gsl_vector_view gslU = gsl_vector_view_array(U  , N);
  gsl_vector_view gslR = gsl_vector_view_array(rhs, N);

  gsl_blas_dgemv(CblasNoTrans, 1.0, iMM, &gslR.vector, 0.0, &gslU.vector);

  // Clean up
  gsl_matrix_free(iMM);

  return 0;
}
