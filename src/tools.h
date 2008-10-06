#ifndef INCLUDE_TOOLS_H

int nWave( const double x, const double t, const double nu, const double c, const double a, const double t0,
	   double *u );
int initializeRandNumGen();
int boundaryCondition(const char *fname, const double time0, const double time1, const double UB0[2], double UB1[2]); 
int initialCondition(const char *fname, const double UB[2], const int N, double *U);


// Structure for storing quad points and basis evaled at quad points
typedef struct
{
  int Nquad; // number of quadrature points
  double *xq; // quadrature point locations
  double *wq; // quadrature weights
  double *phi; // basis evaled at quad points
  double *phi_xi; // basis gradient evaled at quad points

} quadBasis;

int allocateQuadratureAndBasisForResidual(const int Nmode, quadBasis **pQB);
int evaluateQuadratureAndBasisForResidual(const int Nmode, quadBasis **pQB);
int freeQuadratureAndBasisForResidual(quadBasis *pQB);

#define INCLUDE_TOOLS_H 1
#endif
