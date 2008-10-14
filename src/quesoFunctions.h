int solveForStateAtXLocations(const double *xx, const int Nx, const double nu, const double kappa, 
			      gsl_vector *Uic, quadBasis *pQB, double *u);
int computeGradientAtOne(const double kappa, quadBasis *pQB, double *u_x1);
