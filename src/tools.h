int nWave( const double x, const double t, const double nu, const double c, const double a, const double t0,
	   double *u );
int initializeRandNumGen();
int boundaryCondition(const char *fname, const double time0, const double time1, const double UB0[2], double UB1[2]); 
int initialCondition(const char *fname, const double UB[2], const int N, double *U);
