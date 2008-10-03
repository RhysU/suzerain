#define BOOST_TEST_MODULE $Id: legendrePoly_test.cc 173 2008-09-09 oliver $
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <stdio.h>
#include "legendrePoly.h"
#include "solver.h"

BOOST_AUTO_TEST_CASE( test_EigenDecomp )
{
  int ierr;
  const double tol=1e-12;
  double identity[9], MM[9], KK[9];
  double lambda[3], RR[9], iRR[9];

  // Identity
  for( int ii=0; ii<9; ii++ ) identity[ii] = 0.0;
  identity[0] = identity[4] = identity[8] = 1.0;

  ierr = eigenDecompIntFactorMatrix(3, identity, identity, lambda, RR, iRR);
  BOOST_REQUIRE( ierr == 0 );
  for( int ii=0; ii<3; ii++ ) BOOST_CHECK_CLOSE(lambda[ii], 1.0, tol);
  BOOST_CHECK_CLOSE(RR[0], 1.0, tol);
  BOOST_CHECK_CLOSE(RR[4], 1.0, tol);
  BOOST_CHECK_CLOSE(RR[8], 1.0, tol);

  // MM\KK (comparing to matlab)
  ierr = legendre0MassMatrix(3, MM);
  BOOST_REQUIRE( ierr == 0 );

  ierr = legendre0StiffnessMatrix(3, KK);
  BOOST_REQUIRE( ierr == 0 );

  ierr = eigenDecompIntFactorMatrix(3, MM, KK, lambda, RR, iRR);
  BOOST_REQUIRE( ierr == 0 );
  BOOST_CHECK_CLOSE(lambda[0], 2.467437405329207e+00, tol);
  BOOST_CHECK_CLOSE(lambda[1], 2.553256259467077e+01, tol);
  BOOST_CHECK_CLOSE(lambda[2], 1.049999999999999e+01, tol);

  BOOST_CHECK_CLOSE(RR[0],  9.876163114899669e-01, tol);
  BOOST_CHECK_CLOSE(RR[1], -2.829609272384573e-01, tol);
  BOOST_CHECK_SMALL(RR[2], 1e-15);
  BOOST_CHECK_SMALL(RR[3], 1e-15);
  BOOST_CHECK_SMALL(RR[4], 1e-15);
  BOOST_CHECK_CLOSE(RR[5], -1.000000000000000e+00, tol);
  BOOST_CHECK_CLOSE(RR[6], -1.568885632509662e-01, tol);
  BOOST_CHECK_CLOSE(RR[7], -9.591314371119072e-01, tol);
  BOOST_CHECK_SMALL(RR[8], 1e-15);
}
