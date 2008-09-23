#include <stdio.h>

#include "legendrePoly.h"

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
