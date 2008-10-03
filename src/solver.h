#include "burgers.h"
int steadyNewton(const int NiterMax, burgersSteady *pB);
int unsteadyRK4(burgersUnsteady *pB);
int eigenDecompIntFactorMatrix( const int N, const double *MM, const double *KK,
				double *lambda, double *RR, double *iRR);
