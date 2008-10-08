#include "tools.h"
int interiorResidual( const int N, const double nu, const double *U, const double *UB, quadBasis *pQB,
		      const bool RABFlag, const double kappa, double *R, double *R_U );
int unsteadyBoundaryResidual( const int N, const double *UB0, const double *UB1, quadBasis *pQB,
			      double *R );
