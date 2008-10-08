#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "legendrePoly.h"
#include "residual.h"
#include "solver.h"
#include "tools.h"
#include "burgers.h"

//-----------------------------------------------
// Compute soln to RAB equation (with algebraic
// burgulence model) at requested x locations 
// using burgulence model parameter kappa
//
int
solveForStateAtXLocations(const double *xx, const int Nx, const double kappa, double *u)
{

  int ierr, Nmode = 100, nIter = 20;
  double *U, phi[Nmode];
  burgersSteady Bsteady, *pBsteady = &Bsteady;

  // Set up problem data
  pBsteady->RABFlag = true;
  pBsteady->Nmode = Nmode;
  pBsteady->nu = 1e-2;
  pBsteady->kappa = kappa;
  pBsteady->UB[0] = 1.0;
  pBsteady->UB[1] = 0.0;
  pBsteady->U = gsl_vector_calloc(pBsteady->Nmode);
  pBsteady->R = gsl_vector_calloc(pBsteady->Nmode);

  ierr = initialCondition("linear", pBsteady->UB, pBsteady->Nmode, pBsteady->U->data);
  if( ierr != 0 ) return ierr;

  // Solve
  ierr = steadyNewton(nIter, pBsteady);
  if( ierr != 0 ) return ierr;

  // Interpolate solution to desired x locations
  U = pBsteady->U->data;
  for( int ixloc=0; ixloc<Nx; ixloc++ ){

    ierr = legendrePolyZero(pBsteady->Nmode, xx[ixloc], phi, (double *)NULL);
    if( ierr != 0 ) return ierr;

    u[ixloc] = 0.0;
    for( int imode=0; imode<pBsteady->Nmode; imode++ ) u[ixloc] += U[imode]*phi[imode];
  }

  // Clean up memory
  gsl_vector_free(pBsteady->U);
  gsl_vector_free(pBsteady->R);

  return 0;
}
