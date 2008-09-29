
#ifndef INCLUDE_BURGERS_H

// Structure to allow for easy passing/maniuplation of some data
typedef struct
{
  int N; // number of modes
  double nu; // viscosity
  gsl_vector *U; // state vector
  gsl_vector *R; // residual vector

} burgers;


int boundaryCondition(const double time, double *UB);
int InitialCondition(const int N, double *U);

#define INCLUDE_BURGERS_H 1
#endif
