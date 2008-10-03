
#ifndef INCLUDE_BURGERS_H

// Structures to allow for easy passing/maniuplation of some data

// For steady cases
typedef struct
{
  int Nmode; // number of modes
  double nu; // viscosity
  double UB[2]; // boundary conditions 
  gsl_vector *U; // state vector (for DNS this is running average)
  gsl_vector *R; // residual vector

} burgersSteady;


// For unsteady cases
typedef struct
{
  int Nmode; // number of modes

  int Nstep; // number of time steps to take
  int Nwrite; // write solution every Nwrite steps
  int Nstat; // start computing statistics at step Nstat

  double nu; // viscosity
  double time; // time
  double dt; // time step length

  double UB[2]; // boundary conditions at current time
  gsl_vector *U; // state vector at current time

  double UBavg[2]; // running average of boundary conditions
  gsl_vector *Uavg; // running average of state

  char BCname[80];

} burgersUnsteady;

int boundaryCondition(const char *fname, const double time0, const double time1, const double UB0[2], double UB1[2]); 
int initialCondition(const char *fname, const int N, double *U);

#define INCLUDE_BURGERS_H 1
#endif
