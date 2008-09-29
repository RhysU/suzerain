#define BOOST_TEST_MODULE $Id: residual_test.cc 214 2008-09-25 oliver $
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <stdio.h>
#include "residual.h"


BOOST_AUTO_TEST_CASE( test_residual )
{
  int ierr;
  double U[2] = {0.0, 0.0};
  double UB[2] = {0.0, 0.0};
  double R[2] = {0.0, 0.0};
  const double tol=1e-12;

  // u=0 is exact solution for homogeneous Dirichlet BCs.  Thus, R=0
  ierr = interiorResidual(2, 1.0, U, UB, R, (double *)NULL);
  BOOST_REQUIRE( ierr == 0 );
  
  BOOST_CHECK_SMALL(R[0], tol);
  BOOST_CHECK_SMALL(R[1], tol);


  // u=-x, zero viscosity
  R[0] = R[1] = 0.0;
  UB[0] = 1.0; UB[1] = -1.0;
  ierr = interiorResidual(2, 0.0, U, UB, R, (double *)NULL);
  BOOST_REQUIRE( ierr == 0 );
  
  BOOST_CHECK_SMALL(R[0], tol);
  BOOST_CHECK_CLOSE(R[1], 4.0/15.0, tol);

  // u=(1-x^2), zero viscosity
  R[0] = R[1] = 0.0;
  U[0] = 1.0; U[1] = 0.0;
  UB[0] = 0.0; UB[1] = 0.0;
  ierr = interiorResidual(2, 0.0, U, UB, R, (double *)NULL);
  BOOST_REQUIRE( ierr == 0 );
  
  BOOST_CHECK_SMALL(R[0], tol);
  BOOST_CHECK_CLOSE(R[1], -32.0/105.0, tol);
}


BOOST_AUTO_TEST_CASE( ping_residual )
{

  int ierr;
  double U0[3] =  {0.5, 1.0, 0.5};
  double U[3]  =  {0.5, 1.0, 0.5};
  const double UB[2] = {1.0, -1.0};
  const double eps[2] = {1e-2, 1e-3};
  double R0[3], R_U0[9], Rp[3], R_Up[9];
  double analytic, finitedi, diff[3][2], rate;

  // Zero things out
  for( int ii=0; ii<3; ii++ ) R0[ii] = Rp[ii] = 0.0;
  for( int ii=0; ii<9; ii++ ) R_U0[ii] = R_Up[ii] = 0.0;
  
  // Baseline
  ierr = interiorResidual(3, 0.1, U, UB, R0, R_U0);
  BOOST_REQUIRE( ierr == 0 );

  for( int ii=0; ii<3; ii++ ){ // ping for each state component

    for( int jj=0; jj<2; jj++ ){
      
      // perturb the jjth state component
      for( int kk=0; kk<3; kk++ ) U[kk] = U0[kk];
      U[ii] += eps[jj];

      // zero things out
      for( int kk=0; kk<3; kk++ ) Rp[kk] = 0.0;
      for( int kk=0; kk<9; kk++ ) R_Up[kk] = 0.0;

      // Residual and Jacobian
      ierr = interiorResidual(3, 0.1, U, UB, Rp, R_Up);
      BOOST_REQUIRE( ierr == 0 );

      // compute difference
      for( int kk=0; kk<3; kk++ ){
	finitedi = (Rp[kk] - R0[kk])/eps[jj];
	analytic = 0.5*(R_Up[3*kk+ii] + R_U0[3*kk+ii]);
	diff[kk][jj] = fabs(analytic - finitedi);
	//printf("diff[%d][%d](%d) = %.6E\n", kk, ii, jj, diff[kk][jj]); fflush(stdout);
      }
    }

    // check convergence (either error is small or rate is correct)
    for( int kk=0; kk<3; kk++ ){
      if( !(diff[kk][0] < 1e-12) ){
	rate = log(diff[kk][1]/diff[kk][0])/log(eps[2]/eps[1]);
	BOOST_CHECK( rate > 1.9 );
      } else{
	BOOST_CHECK( (diff[kk][0] < 1e-12) && (diff[kk][1] < 1e-12) );
      }
    }

  }

}


BOOST_AUTO_TEST_CASE( test_unsteadyBoundaryResidual )
{
  int ierr;
  double UB0[2] = {0.0, 0.0};
  double UB1[2] = {0.0, 0.0};
  double R[2] = {0.0, 0.0};
  const double tol=1e-12;

  // BCs do not change
  ierr = unsteadyBoundaryResidual(2, UB0, UB1, R);
  BOOST_REQUIRE( ierr == 0 );
  
  BOOST_CHECK_SMALL(R[0], tol);
  BOOST_CHECK_SMALL(R[1], tol);

  // change the same on both sides
  UB0[0] = 0.0; UB0[1] = 0.0;
  UB1[0] = 1.5; UB1[1] = 1.5;

  R[0] = R[1] = 0.0;
  ierr = unsteadyBoundaryResidual(2, UB0, UB1, R);
  BOOST_REQUIRE( ierr == 0 );
  
  BOOST_CHECK_CLOSE(R[0], 2.0, tol);
  BOOST_CHECK_SMALL(R[1], tol);

  // different changes on each side
  UB0[0] = 0.0; UB0[1] = 0.0;
  UB1[0] = 1.0; UB1[1] = 2.0;

  R[0] = R[1] = 0.0;
  ierr = unsteadyBoundaryResidual(2, UB0, UB1, R);
  BOOST_REQUIRE( ierr == 0 );
  
  BOOST_CHECK_CLOSE(R[0], 2.0, tol);
  BOOST_CHECK_CLOSE(R[1], 2.0/15.0, tol);

}
