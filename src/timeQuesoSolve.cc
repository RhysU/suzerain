#include <stdlib.h>
#include <math.h>

#include "burgers.h"
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
  const int Nx = 10;
  const double xx[10] = {-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
  const double kappa0 = 0.3;
  double kappa;
  gsl_vector *U;
  quadBasis *pQB;
  double u[10];

  ierr = evaluateQuadratureAndBasisForResidual(100, &pQB);
  if( ierr != 0 ) return ierr;

  U = gsl_vector_calloc(100);

  for( int ii=0; ii<200; ii++ ){
    kappa = kappa0 + 0.25*sin((double)ii);
    printf("kappa = %.6E\n", kappa); fflush(stdout);
    ierr = solveForStateAtXLocations(xx, Nx, 1e-2, kappa, U, pQB, u);
    if( ierr != 0 ) return ierr;
  }

  freeQuadratureAndBasisForResidual(pQB);
  free(pQB);

  return 0;
}
