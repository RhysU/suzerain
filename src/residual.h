int interiorResidual( const int N, const double nu, const double *U, const double *UB, double *R, double *R_U );
int unsteadyBoundaryResidual( const int N, const double *UB0, const double *UB1, double *R );
