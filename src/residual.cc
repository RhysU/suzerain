#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "legendrePoly.h"
#include "tools.h"
#include "residual.h"

//------------------------------------------------------------------------
// Compute eddy viscosity for RAB simulations using
// an algebraic burgulence model.
int
eddyViscosity(const double xq, const double u, const double u_x, double *nut, double *nut_u, double *nut_u_x)
{
  double kappa = 0.3;
  double dist = (1.0-xq);
  double lmix = kappa*dist;

  double epsilon = 1.0e-6;
  double om = sqrt( u_x*u_x + epsilon*epsilon ); // smooth approx of abs(u_x)

  double om_u = 0.0;
  double om_u_x = u_x/om;

  *nut     = lmix*lmix*om;
  *nut_u   = lmix*lmix*om_u;
  *nut_u_x = lmix*lmix*om_u_x;

  return 0;
}


//------------------------------------------------------------------------
// Compute residual vector (in simplest, probably slowest, way possible).
// Should refactor to improve speed and modularity.                      
int
interiorResidual( const int N, const double nu, const double *U, const double *UB, quadBasis *pQB,
		  bool RABFlag, double *R, double *R_U )
{
  int ierr, Nquad;
  double *xq, *wq, *phi, *phi_xi;

  double u, u_x, F, F_u, F_u_x, nut, nut_u, nut_u_x;
  double nue = nu, nue_u = 0.0, nue_u_x = 0.0;
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
    

    if( RABFlag ){
      ierr = eddyViscosity(xq[iquad], u, u_x, &nut, &nut_u, &nut_u_x);
      if( ierr != 0 ) return ierr;

      // add to effective viscosity
      nue     = nu + nut;
      nue_u   =      nut_u;
      nue_u_x =      nut_u_x;
    }

    // Compute flux and derivatives
    F     = 0.5*u*u - nue*u_x;

    if( R_U != (double *)NULL ){
      F_u   =     u - nue_u  *u_x      ;
      F_u_x =       - nue_u_x*u_x - nue;
    }
    
    // Add to residual integral
//     F *= wq[iquad];
//     for( int ii=0; ii<N; ii++ ) R[ii] -= phi_xi[ii]*F;

    //for( int ii=0; ii<N; ii++ ) R[ii] -= wq[iquad]*( phi_xi[ii]*F - phi[ii]*0.5*u );
    for( int ii=0; ii<N; ii++ ) R[ii] -= wq[iquad]*( phi_xi[ii]*F - phi[ii]*(1.0 - xq[iquad])/4.0 );

    
    // Add to Jacobian integral (if requested)
    if( R_U != (double *)NULL ){
      for( int ii=0; ii<N; ii++ ){
	for( int jj=0; jj<N; jj++ ){
	  R_U[N*ii+jj] -= wq[iquad]*phi_xi[ii]*(F_u*phi[jj] + F_u_x*phi_xi[jj]);
	  //R_U[N*ii+jj] -= wq[iquad]*( phi_xi[ii]*(F_u*phi[jj] + F_u_x*phi_xi[jj]) - phi[ii]*0.5*phi[jj] );
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
