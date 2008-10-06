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

  BOOST_CHECK_CLOSE(iRR[0],  9.672103658950829e-01, tol);
  BOOST_CHECK_SMALL(iRR[1], 1e-15);
  BOOST_CHECK_CLOSE(iRR[2], -2.853443557146049e-01, tol);
  BOOST_CHECK_CLOSE(iRR[3], -1.582100625578974e-01, tol);
  BOOST_CHECK_SMALL(iRR[4], 1e-15);
  BOOST_CHECK_CLOSE(iRR[5], -9.959351732610457e-01, tol);
  BOOST_CHECK_SMALL(iRR[6], 1e-15);
  BOOST_CHECK_CLOSE(iRR[7], -1.0                  , tol);
  BOOST_CHECK_SMALL(iRR[8], tol);

}


int
integratingFactor( const int N, const double *RR, const double *lambda, const double *iRR, 
		   const double nu, const double time,
		   double *intFactor);

BOOST_AUTO_TEST_CASE( test_IntFactor )
{

  int ierr;
  
  const double tol=1e-10;
  double identity[9], lambda[3] = {0.0, 0.0, 0.0};
  double RR[9], iRR[9];
  double intFactor[9];

  // Identity eigenvectors with zero eigenvalues
  for( int ii=0; ii<9; ii++ ) identity[ii] = 0.0;
  identity[0] = identity[4] = identity[8] = 1.0;

  ierr = integratingFactor(3, identity, lambda, identity, 1.0, 1.0, intFactor);
  BOOST_REQUIRE( ierr == 0 );

  BOOST_CHECK_CLOSE(intFactor[0], 1.0, tol);
  BOOST_CHECK_SMALL(intFactor[1], 1e-15);
  BOOST_CHECK_SMALL(intFactor[2], 1e-15);
  BOOST_CHECK_SMALL(intFactor[3], 1e-15);
  BOOST_CHECK_CLOSE(intFactor[4], 1.0, tol);
  BOOST_CHECK_SMALL(intFactor[5], 1e-15);
  BOOST_CHECK_SMALL(intFactor[6], 1e-15);
  BOOST_CHECK_SMALL(intFactor[7], 1e-15);
  BOOST_CHECK_CLOSE(intFactor[8], 1.0, tol);


  // RR, lambda from real 3x3 case
  RR[0] = -9.876163114899669e-01;
  RR[1] = -2.829609272384573e-01;
  RR[2] = 0.0;
  RR[3] = 0.0;
  RR[4] = 0.0;
  RR[5] =  1.000000000000000e+00;
  RR[6] =  1.568885632509661e-01;
  RR[7] = -9.591314371119072e-01;
  RR[8] = 0.0;

  lambda[0] = 2.467437405329207e+00;
  lambda[1] = 2.553256259467077e+01;
  lambda[2] = 1.049999999999999e+01;

  iRR[0] = -9.672103658950829e-01;
  iRR[1] = 0.0;
  iRR[2] =  2.853443557146049e-01;
  iRR[3] = -1.582100625578974e-01;
  iRR[4] = 0.0;
  iRR[5] = -9.959351732610457e-01;
  iRR[6] =  0.0;
  iRR[7] =  1.000000000000000e+00;
  iRR[8] =  0.0;

  ierr = integratingFactor(3, RR, lambda, iRR, 1.0, 1.0, intFactor);
  BOOST_REQUIRE( ierr == 0 );

  BOOST_CHECK_CLOSE(intFactor[0], 5.490496472044510e+09, tol);
  BOOST_CHECK_SMALL(intFactor[1], 1e-15);
  BOOST_CHECK_CLOSE(intFactor[2], 3.456277341037881e+10, tol);
  BOOST_CHECK_SMALL(intFactor[3], 1e-15);
  BOOST_CHECK_CLOSE(intFactor[4], 3.631550267424635e+04, tol);
  BOOST_CHECK_SMALL(intFactor[5], 1e-15);
  BOOST_CHECK_CLOSE(intFactor[6], 1.861072414405010e+10, tol);
  BOOST_CHECK_SMALL(intFactor[7], 1e-15);
  BOOST_CHECK_CLOSE(intFactor[8], 1.171548413363452e+11, tol);
}
