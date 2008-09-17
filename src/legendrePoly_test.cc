#define BOOST_TEST_MODULE $Id: legendrePoly_test.cc 173 2008-09-09 oliver $
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "legendrePoly.h"

BOOST_AUTO_TEST_CASE( test_L0 )
{
  int ierr;
  const double tol=1e-15;
  double L, L_xi;

  ierr = legendrePoly(1, -1.0, &L, &L_xi);
  BOOST_REQUIRE( ierr == 0 );
  BOOST_CHECK_CLOSE(L, 1.0, tol);
  BOOST_CHECK_CLOSE(L_xi, 0.0, tol);

  ierr = legendrePoly(1, 0.0, &L, &L_xi);
  BOOST_REQUIRE( ierr == 0 );
  BOOST_CHECK_CLOSE(L, 1.0, tol);
  BOOST_CHECK_CLOSE(L_xi, 0.0, tol);

  ierr = legendrePoly(1, 1.0, &L, &L_xi);
  BOOST_REQUIRE( ierr == 0 );
  BOOST_CHECK_CLOSE(L, 1.0, tol);
  BOOST_CHECK_CLOSE(L_xi, 0.0, tol);
}


BOOST_AUTO_TEST_CASE( test_L1 )
{
  int ierr;
  const double tol=1e-15;
  double L[2], L_xi[2];

  ierr = legendrePoly(2, -1.0, L, L_xi);
  BOOST_REQUIRE( ierr == 0 );
  BOOST_CHECK_CLOSE(L[0], 1.0, tol);
  BOOST_CHECK_CLOSE(L[1],-1.0, tol);
  BOOST_CHECK_CLOSE(L_xi[0], 0.0, tol);
  BOOST_CHECK_CLOSE(L_xi[1], 1.0, tol);

  ierr = legendrePoly(2, 0.0, L, L_xi);
  BOOST_REQUIRE( ierr == 0 );
  BOOST_CHECK_CLOSE(L[0], 1.0, tol);
  BOOST_CHECK_CLOSE(L[1], 0.0, tol);  
  BOOST_CHECK_CLOSE(L_xi[0], 0.0, tol);
  BOOST_CHECK_CLOSE(L_xi[1], 1.0, tol);

  ierr = legendrePoly(2, 1.0, L, L_xi);
  BOOST_REQUIRE( ierr == 0 );
  BOOST_CHECK_CLOSE(L[0], 1.0, tol);
  BOOST_CHECK_CLOSE(L[1], 1.0, tol);
  BOOST_CHECK_CLOSE(L_xi[0], 0.0, tol);
  BOOST_CHECK_CLOSE(L_xi[1], 1.0, tol);
}

BOOST_AUTO_TEST_CASE( test_L2 )
{
  int ierr;
  const double tol=1e-15;
  double L[3], L_xi[3];

  ierr = legendrePoly(3, -1.0, L, L_xi);
  BOOST_REQUIRE( ierr == 0 );
  BOOST_CHECK_CLOSE(L[0], 1.0, tol);
  BOOST_CHECK_CLOSE(L[1],-1.0, tol);
  BOOST_CHECK_CLOSE(L[2], 1.0, tol);  
  BOOST_CHECK_CLOSE(L_xi[0], 0.0, tol);
  BOOST_CHECK_CLOSE(L_xi[1], 1.0, tol);
  BOOST_CHECK_CLOSE(L_xi[2],-3.0, tol);

  ierr = legendrePoly(3, 0.0, L, L_xi);
  BOOST_REQUIRE( ierr == 0 );
  BOOST_CHECK_CLOSE(L[0], 1.0, tol);
  BOOST_CHECK_CLOSE(L[1], 0.0, tol);
  BOOST_CHECK_CLOSE(L[2],-0.5, tol);  
  BOOST_CHECK_CLOSE(L_xi[0], 0.0, tol);
  BOOST_CHECK_CLOSE(L_xi[1], 1.0, tol);
  BOOST_CHECK_CLOSE(L_xi[2], 0.0, tol);

  ierr = legendrePoly(3, 1.0, L, L_xi);
  BOOST_REQUIRE( ierr == 0 );
  BOOST_CHECK_CLOSE(L[0], 1.0, tol);
  BOOST_CHECK_CLOSE(L[1], 1.0, tol);
  BOOST_CHECK_CLOSE(L[2], 1.0, tol);
  BOOST_CHECK_CLOSE(L_xi[0], 0.0, tol);
  BOOST_CHECK_CLOSE(L_xi[1], 1.0, tol);
  BOOST_CHECK_CLOSE(L_xi[2], 3.0, tol);
}
