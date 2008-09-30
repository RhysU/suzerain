#include "burgers.h"
int steadyNewton(const int NiterMax, burgers *pB);
int unsteadyRK4(const int Nstep, const double dt, burgers *pB);
