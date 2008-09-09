#include <stdio.h>

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


#ifdef UNIT_TEST
#include "legendrePoly_test.cc"
#endif
