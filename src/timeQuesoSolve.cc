#include <stdlib.h>

#include "tools.h"
#include "quesoFunctions.h"


//-------------------------------------------------
// Main routine.  Run steady solver many many times
// to get timings and profile.
//
int
main( int argc, char * argv[] )
{
  int ierr;
  const double xx[10] = {-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
  const int Nx = 10;
  const double kappa = 0.3;
  quadBasis *pQB;
  double u[10];

  ierr = evaluateQuadratureAndBasisForResidual(100, &pQB);
  if( ierr != 0 ) return ierr;

  for( int ii=0; ii<200; ii++ ){
    ierr = solveForStateAtXLocations(xx, Nx, kappa, pQB, u);
    if( ierr != 0 ) return ierr;
  }

  freeQuadratureAndBasisForResidual(pQB);
  free(pQB);

  return 0;
}
