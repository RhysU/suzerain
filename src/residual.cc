#include <stdio.h>
#include <math.h>

#include "legendrePoly.h"
#include "residual.h"


//------------------------------------------------------------------------
// Compute residual vector (in simplest, probably slowest, way possible).
// Should refactor to improve speed and modularity.                      
int
interiorResidual( const int N, const double nu, const double *U, const double *UB,
		  double *R, double *R_U )
{
  int ierr;
  const int solnOrder = N-1+2;
  const int quadOrder = 3*solnOrder - 1;
  const int Nquad = 3*quadOrder/2 + 1; // + 1 to be conservative

  double phi[N], phi_xi[N];
  double xq[Nquad], wq[Nquad];

  double u, u_x, F, F_u, F_u_x;

  // Get quadrature points
  ierr = legendreGaussQuad(Nquad, xq, wq);
  if( ierr != 0 ) return ierr;

  // Loop over quadrature points
  for( int iquad=0; iquad<Nquad; iquad++ ){

    // Evaluate basis
    ierr = legendrePolyZero(N, xq[iquad], phi, phi_xi);
    if( ierr != 0 ) return ierr;
  
    // Interpolate state/gradient to current point
    u = 0.0; u_x = 0.0;
    for( int ii=0; ii<N; ii++ ){
      u += U[ii]*phi[ii];
      u_x += U[ii]*phi_xi[ii];
    }

    // Add BC contribution to state/gradient
    u   += UB[0]*0.5*(1.0-xq[iquad]) + UB[1]*0.5*(1.0+xq[iquad]);
    u_x += UB[0]*0.5*(   -1.0      ) + UB[1]*0.5*(1.0          );


    // Compute flux and derivatives
    F     = 0.5*u*u - nu*u_x;
    F_u   =     u           ;
    F_u_x =         - nu    ;
    
    // Add to residual integral
    for( int ii=0; ii<N; ii++ ) R[ii] -= wq[iquad]*phi_xi[ii]*F;
    
    // Add to Jacobian integral (if requested)
    if( R_U != (double *)NULL ){
      for( int ii=0; ii<N; ii++ ){
	for( int jj=0; jj<N; jj++ ){
	  R_U[N*ii+jj] -= wq[iquad]*phi_xi[ii]*(F_u*phi[jj] + F_u_x*phi_xi[jj]);
	}
      }
    }

  } // end loop over quad points

  return 0;
}
