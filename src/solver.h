#include "burgers.h"
int steadyNewton(const int NiterMax, burgersSteady *pB);
int steadyNewtonPreComputeBasis(const int NiterMax, burgersSteady *pB, quadBasis *pQB);
int unsteadyRK4(burgersUnsteady *pB);
int eigenDecompIntFactorMatrix( const int N, const double *MM, const double *KK,
				double *lambda, double *RR, double *iRR);
int integratingFactor( const int N, const double *RR, const double *lambda, const double *iRR, 
		       const double nu, const double time, double *intFactor);
