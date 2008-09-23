#include <stdio.h>
#include <math.h>

#include "legendrePoly.h"

#define PI 3.14159265358979323846

/*----------------------------------------------------------------*/
/* Evaluates Lagrange polynomials using 3 term recursion.         */
int
legendrePoly(const int N, const double xi, double * L, double * L_xi)
{
  double nn;

  if( N<1 ){
    printf("ERROR (legendrePoly): Must request at least linear!\n");
    return -1;
  }

  // values
  L[0] = 1.0;
  L[1] = xi;
  for( int ii=1; ii<N-1; ii++ ){
    nn = (double)ii;
    L[ii+1] = ((2.0*nn+1.0)/(nn+1.0))*xi*L[ii] - (nn/(nn+1.0))*L[ii-1];
  }

   // derivatives
   if( L_xi != (double *)NULL ){
     L_xi[0] = 0.0;
     L_xi[1] = 1.0;
     for( int ii=1; ii<N-1; ii++ ){
       nn = (double)ii;
       L_xi[ii+1] = ((2.0*nn+1.0)/(nn+1.0))*(L[ii] + xi*L_xi[ii]) - (nn/(nn+1.0))*L_xi[ii-1];
     }
   }

  return 0;
}


/*----------------------------------------------------------------*/
/* Evaluates (1-xi^2)*Li, where Li is the ith Legendre poly       */
int
legendrePolyZero(const int N, const double xi, double * phi, double * phi_xi)
{

  const int ierr = legendrePoly(N, xi, phi, phi_xi);
  if( ierr != 0 ) return ierr;

  const double imx2 = (1.0 - xi*xi);
  const double n2xi  = -2.0*xi;
  for( int ii=0; ii<N; ii++ ){
    // Do derivatives first (so we don't have to store L separately)
    if( phi_xi != (double *)NULL ) phi_xi[ii] = n2xi*phi[ii] + imx2*phi_xi[ii];
    phi[ii] *= imx2;
  }

  return 0;
}


/*----------------------------------------------------------------*/
/* Computes the nodes and weights for Legendre-Gauss quadrature   */
int
legendreGaussQuad(const int N, double * xq, double * wq)
{

  const int ITERMAX = 200;
  const double RTOL = 1e-16;

  int ierr, iiter;
  double L[N], L_x[N], res, xtmp, dx;

  // Newton solve for roots of Legendre polynomial of order N
  for( int iquad=0; iquad<N; iquad++ ){

    // initial guess (taken from matlab implementation by Greg von Winckel)
    xtmp = -cos( (2.0*iquad+1.0)*PI/(2*N) ) - (0.27/N)*sin( PI*(-1.0+2.0*iquad/(N-1))*(N-1)/(N+1) );

    iiter = 0;

    do{
      ierr = legendrePoly(N+1, xtmp, L, L_x);
      if( ierr != 0 ) return ierr;
      
      res = L[N];
      
      dx = -res/L_x[N];
      xtmp += dx;

      iiter++;

    } while( (fabs(res)>RTOL) && (iiter < ITERMAX) );

    // warn if point is not at least close to converged
    if( fabs(res)>1e3*RTOL ){
      printf("WARNING: Quad point did not converge.\n");
      printf("         x = %.6E, residual = %.6E\n", xtmp, res );
      fflush(stdout);
    }

    xq[iquad] = xtmp;

    // compute weight
    wq[iquad] = 2.0/( (1.0-xtmp*xtmp)*L_x[N]*L_x[N] );
  }

  return 0;
}
