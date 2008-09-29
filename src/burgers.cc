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


//-------------------------------------------------
// Evaluate BC function
//
int
boundaryCondition(const double time, double *UB)
{ 
  int ierr;
  
//   // Steady
//   UB[0] =  1.0;
//   UB[1] = -1.0;

  // n wave
  ierr = nWave(-1.0, time, 1e-2, 1.0, 16.0, 1.0, &(UB[0]));
  ierr = nWave( 1.0, time, 1e-2, 1.0, 16.0, 1.0, &(UB[1]));

  return 0;
}

//-------------------------------------------------
// Set initial state
//
int
initialCondition(const int N, double *U)
{ 
  // Just a placeholder for now
  //for( int ii=0; ii<N; ii++ ) U[ii] = 0.0;

//   // Steady solution (33 modes, nu=0.1)
//   U[ 0] = 1.629157349371976E-15;
//   U[ 1] = -8.669500522642464E-01;
//   U[ 2] = -3.456551422036912E-15;
//   U[ 3] = 6.069143024375602E-01;
//   U[ 4] = 3.682725646668788E-15;
//   U[ 5] = -3.893348038666424E-01;
//   U[ 6] = -3.080021664329973E-15;
//   U[ 7] = 2.411008158747520E-01;
//   U[ 8] = 2.892392795861849E-15;
//   U[ 9] = -1.454095440474666E-01;
//   U[10] = -1.779750909111454E-15;
//   U[11] = 8.619460994958301E-02;
//   U[12] = 1.465217249531183E-15;
//   U[13] = -5.019382023701988E-02;
//   U[14] = -9.255112289908806E-16;
//   U[15] = 2.912267097927002E-02;
//   U[16] = 4.635076343276719E-16;
//   U[17] = -1.652185694847467E-02;
//   U[18] = -2.854236426905883E-16;
//   U[19] = 9.573571092176521E-03;
//   U[20] = 2.455572488794693E-16;
//   U[21] = -5.225493403147587E-03;
//   U[22] = -9.632359084548544E-17;
//   U[23] = 3.154854342162039E-03;
//   U[24] = 2.885938241697400E-16;
//   U[25] = -1.542924832957964E-03;
//   U[26] = 7.812360199810137E-17;
//   U[27] = 1.111393254515644E-03;
//   U[28] = 8.289866806804461E-17;
//   U[29] = -3.576972413726743E-04;
//   U[30] = 1.732423311034709E-16;
//   U[31] = 4.876094701709677E-04;
//   U[32] = -6.126450964786364E-18;

  // L2 projection of n wave solution at t=0 onto solution space
  int ierr;
  const int solnOrder = N-1+2;
  const int quadOrder = 4*solnOrder + 1; // exact for n wave order 3*N 
                                         // (of course, n wave is not poly, but hopefully error is small enough)
  const int Nquad = quadOrder/2 + 1;

  double xq[Nquad], wq[Nquad];
  double phi[N];
  double UB[2], rhs[N], u, u0;

  gsl_matrix *iMM;

  // Get quadrature points
  ierr = legendreGaussQuad(Nquad, xq, wq);
  if( ierr != 0 ) return ierr;

  // Get BCs at time 0
  ierr = boundaryCondition(0.0, UB);
  if( ierr != 0 ) return ierr;

  FILE *fp = fopen ("nWave.dat", "w");
  fprintf(fp, "%.15E, %.15E\n", UB[0], UB[1]);
  fclose(fp);

  // zero rhs
  for( int ii=0; ii<N; ii++ ) rhs[ii] = 0.0;

  // Loop over quadrature points to form right hand side vector
  for( int iquad=0; iquad<Nquad; iquad++ ){
  
    // Evaluate IC
    ierr = nWave(xq[iquad], 0.0, 1e-2, 1.0, 16.0, 1.0, &u);
    if( ierr != 0 ) return ierr;

    // Write solution to screen         
    FILE *fp = fopen ("nWave.dat", "a");
    fprintf(fp, "%.15E, %.15E\n", xq[iquad], u);
    fclose(fp);

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

  // Set initial condition
  ierr = initialCondition(N, pB->U->data);
  if( ierr != 0 ) return ierr;

  // Call solver (only steady Newton currently available)
//   ierr = steadyNewton(NiterMax, pB);
//   if( ierr != 0 ) return ierr;

  ierr = unsteadyRK4(10000, pB);
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
