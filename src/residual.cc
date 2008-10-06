#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "legendrePoly.h"
#include "tools.h"
#include "residual.h"


//------------------------------------------------------------------------
// Compute residual vector (in simplest, probably slowest, way possible).
// Should refactor to improve speed and modularity.                      
int
interiorResidual( const int N, const double nu, const double *U, const double *UB, quadBasis *pQB,
		  double *R, double *R_U )
{
  int ierr, Nquad;
  double *xq, *wq, *phi, *phi_xi;

  double u, u_x, F, F_u, F_u_x;
  quadBasis *pQBtmp;

  // if necessary, compute quad points, weights and basis fcns
  if( pQB == NULL ){
    ierr = evaluateQuadratureAndBasisForResidual(N, &pQBtmp);
    if( ierr != 0 ) return ierr;
  } else{
    pQBtmp = pQB;
  }

  Nquad = pQBtmp->Nquad;
  xq = pQBtmp->xq;
  wq = pQBtmp->wq;

  // Loop over quadrature points
  for( int iquad=0; iquad<Nquad; iquad++ ){

    // Grab basis fcns evaluated at this quadrature point
    phi    = pQBtmp->phi    + iquad*N;
    phi_xi = pQBtmp->phi_xi + iquad*N;
  
    // Interpolate state/gradient to current point
    u = u_x = 0.0;
    for( int ii=0; ii<N; ii++ ){
      u += U[ii]*phi[ii];
      u_x += U[ii]*phi_xi[ii];
    }

    // Add BC contribution to state/gradient
    u   += UB[0]*0.5*(1.0-xq[iquad]) + UB[1]*0.5*(1.0+xq[iquad]);
    u_x += UB[0]*0.5*(   -1.0      ) + UB[1]*0.5*(1.0          );


    // Compute flux and derivatives
    F     = 0.5*u*u - nu*u_x;

    if( R_U != (double *)NULL ){
      F_u   =     u           ;
      F_u_x =         - nu    ;
    }
    
    // Add to residual integral
    F *= wq[iquad];
    for( int ii=0; ii<N; ii++ ) R[ii] -= phi_xi[ii]*F;
    
    // Add to Jacobian integral (if requested)
    if( R_U != (double *)NULL ){
      for( int ii=0; ii<N; ii++ ){
	for( int jj=0; jj<N; jj++ ){
	  R_U[N*ii+jj] -= wq[iquad]*phi_xi[ii]*(F_u*phi[jj] + F_u_x*phi_xi[jj]);
	}
      }
    }

  } // end loop over quad points

  // if necessary, clean up
  if( pQB==NULL ){
    freeQuadratureAndBasisForResidual(pQBtmp);
    free(pQBtmp);
  }

  return 0;
}


//------------------------------------------------------------------------
// Compute contribution of unsteady BCs to total residual
//
int
unsteadyBoundaryResidual( const int N, const double *UB0, const double *UB1, quadBasis *pQB, double *R )
{

  int ierr, Nquad;
  double ub0, ub1, dub;
  double *xq, *wq, *phi;
  quadBasis *pQBtmp;

  // if necessary, compute quad points, weights and basis fcns
  if( pQB == NULL ){
    ierr = evaluateQuadratureAndBasisForResidual(N, &pQBtmp);
    if( ierr != 0 ) return ierr;
  } else{
    pQBtmp = pQB;
  }

  Nquad = pQBtmp->Nquad;
  xq = pQBtmp->xq;
  wq = pQBtmp->wq;

  // Loop over quad points
  for( int iquad=0; iquad<Nquad; iquad++ ){

    phi = pQBtmp->phi + iquad*N;

    // Interpolate boundary function
    ub0 = UB0[0]*0.5*(1.0-xq[iquad]) + UB0[1]*0.5*(1.0+xq[iquad]);
    ub1 = UB1[0]*0.5*(1.0-xq[iquad]) + UB1[1]*0.5*(1.0+xq[iquad]);

    dub = ub1 - ub0;

    // Add to integral
    for( int ii=0; ii<N; ii++ ) R[ii] += wq[iquad]*phi[ii]*dub;

  }

  // if necessary, clean up
  if( pQB==NULL ){
    freeQuadratureAndBasisForResidual(pQBtmp);
    free(pQBtmp);
  }

  return 0;
}
