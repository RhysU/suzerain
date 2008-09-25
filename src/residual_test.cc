#define BOOST_TEST_MODULE $Id: residual_test.cc 214 2008-09-25 oliver $
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <stdio.h>
#include "residual.h"


BOOST_AUTO_TEST_CASE( test_residual_zero )
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
