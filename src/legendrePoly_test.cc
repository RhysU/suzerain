#define BOOST_TEST_MODULE $Id: legendrePoly_test.cc 173 2008-09-09 oliver $
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <stdio.h>
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


BOOST_AUTO_TEST_CASE( test_phi0 )
{
  int ierr;
  const double tol=1e-15;
  double phi, phi_xi;

  ierr = legendrePolyZero(1, -1.0, &phi, &phi_xi);
  BOOST_REQUIRE( ierr == 0 );
  BOOST_CHECK_CLOSE(phi, 0.0, tol);
  BOOST_CHECK_CLOSE(phi_xi, 2.0, tol);

  ierr = legendrePolyZero(1, 0.0, &phi, &phi_xi);
  BOOST_REQUIRE( ierr == 0 );
  BOOST_CHECK_CLOSE(phi, 1.0, tol);
  BOOST_CHECK_CLOSE(phi_xi, 0.0, tol);

  ierr = legendrePolyZero(1, 1.0, &phi, &phi_xi);
  BOOST_REQUIRE( ierr == 0 );
  BOOST_CHECK_CLOSE(phi, 0.0, tol);
  BOOST_CHECK_CLOSE(phi_xi, -2.0, tol);
}

BOOST_AUTO_TEST_CASE( test_phi1 )
{
  int ierr;
  const double tol=1e-15;
  double phi[2], phi_xi[2];

  ierr = legendrePolyZero(2, -1.0, phi, phi_xi);
  BOOST_REQUIRE( ierr == 0 );
  BOOST_CHECK_CLOSE(phi[0], 0.0, tol);
  BOOST_CHECK_CLOSE(phi[1], 0.0, tol);
  BOOST_CHECK_CLOSE(phi_xi[0], 2.0, tol);
  BOOST_CHECK_CLOSE(phi_xi[1],-2.0, tol);

  ierr = legendrePolyZero(2, 0.0, phi, phi_xi);
  BOOST_REQUIRE( ierr == 0 );
  BOOST_CHECK_CLOSE(phi[0], 1.0, tol);
  BOOST_CHECK_CLOSE(phi[1], 0.0, tol);  
  BOOST_CHECK_CLOSE(phi_xi[0], 0.0, tol);
  BOOST_CHECK_CLOSE(phi_xi[1], 1.0, tol);

  ierr = legendrePolyZero(2, 1.0, phi, phi_xi);
  BOOST_REQUIRE( ierr == 0 );
  BOOST_CHECK_CLOSE(phi[0], 0.0, tol);
  BOOST_CHECK_CLOSE(phi[1], 0.0, tol);
  BOOST_CHECK_CLOSE(phi_xi[0],-2.0, tol);
  BOOST_CHECK_CLOSE(phi_xi[1],-2.0, tol);
}

BOOST_AUTO_TEST_CASE( test_phi2 )
{
  int ierr;
  const double tol=1e-15;
  double phi[3], phi_xi[3];

  ierr = legendrePolyZero(3, -1.0, phi, phi_xi);
  BOOST_REQUIRE( ierr == 0 );
  BOOST_CHECK_CLOSE(phi[0], 0.0, tol);
  BOOST_CHECK_CLOSE(phi[1], 0.0, tol);
  BOOST_CHECK_CLOSE(phi[2], 0.0, tol);  
  BOOST_CHECK_CLOSE(phi_xi[0], 2.0, tol);
  BOOST_CHECK_CLOSE(phi_xi[1],-2.0, tol);
  BOOST_CHECK_CLOSE(phi_xi[2], 2.0, tol);

  ierr = legendrePolyZero(3, 0.0, phi, phi_xi);
  BOOST_REQUIRE( ierr == 0 );
  BOOST_CHECK_CLOSE(phi[0], 1.0, tol);
  BOOST_CHECK_CLOSE(phi[1], 0.0, tol);
  BOOST_CHECK_CLOSE(phi[2],-0.5, tol);  
  BOOST_CHECK_CLOSE(phi_xi[0], 0.0, tol);
  BOOST_CHECK_CLOSE(phi_xi[1], 1.0, tol);
  BOOST_CHECK_CLOSE(phi_xi[2], 0.0, tol);

  ierr = legendrePolyZero(3, 1.0, phi, phi_xi);
  BOOST_REQUIRE( ierr == 0 );
  BOOST_CHECK_CLOSE(phi[0], 0.0, tol);
  BOOST_CHECK_CLOSE(phi[1], 0.0, tol);
  BOOST_CHECK_CLOSE(phi[2], 0.0, tol);
  BOOST_CHECK_CLOSE(phi_xi[0],-2.0, tol);
  BOOST_CHECK_CLOSE(phi_xi[1],-2.0, tol);
  BOOST_CHECK_CLOSE(phi_xi[2],-2.0, tol);
}


BOOST_AUTO_TEST_CASE( test_legendreGaussQuad )
{
  int ierr;
  double tol=1e-13;
  double xq[50], wq[50];
  
  // two point
  ierr = legendreGaussQuad(2, xq, wq);
  BOOST_REQUIRE( ierr == 0 );

  BOOST_CHECK_CLOSE(xq[0], -sqrt(3.0)/3.0, tol);
  BOOST_CHECK_CLOSE(xq[1],  sqrt(3.0)/3.0, tol);

  BOOST_CHECK_CLOSE(wq[0], 1.0, tol);
  BOOST_CHECK_CLOSE(wq[1], 1.0, tol);

  // three point
  ierr = legendreGaussQuad(3, xq, wq);
  BOOST_REQUIRE( ierr == 0 );

  BOOST_CHECK_CLOSE(xq[0], -sqrt(15.0)/5.0, tol);
  BOOST_CHECK_SMALL(xq[1],                  tol); // xq[1] = 0.0
  BOOST_CHECK_CLOSE(xq[2],  sqrt(15.0)/5.0, tol);

  BOOST_CHECK_CLOSE(wq[0], 5.0/9.0, tol);
  BOOST_CHECK_CLOSE(wq[1], 8.0/9.0, tol);
  BOOST_CHECK_CLOSE(wq[2], 5.0/9.0, tol);

  // 40 point  
  ierr = legendreGaussQuad(40, xq, wq);
  BOOST_REQUIRE( ierr == 0 );

  BOOST_CHECK_CLOSE(xq[ 0], -0.998237709710559200350, tol);
  BOOST_CHECK_CLOSE(xq[ 1], -0.990726238699457006453, tol);
  BOOST_CHECK_CLOSE(xq[ 2], -0.977259949983774262663, tol);
  BOOST_CHECK_CLOSE(xq[ 3], -0.957916819213791655805, tol);
  BOOST_CHECK_CLOSE(xq[ 4], -0.932812808278676533361, tol);
  BOOST_CHECK_CLOSE(xq[ 5], -0.902098806968874296728, tol);
  BOOST_CHECK_CLOSE(xq[ 6], -0.865959503212259503821, tol);
  BOOST_CHECK_CLOSE(xq[ 7], -0.824612230833311663196, tol);
  BOOST_CHECK_CLOSE(xq[ 8], -0.778305651426519387695, tol);
  BOOST_CHECK_CLOSE(xq[ 9], -0.727318255189927103281, tol);
  BOOST_CHECK_CLOSE(xq[10], -0.671956684614179548379, tol);
  BOOST_CHECK_CLOSE(xq[11], -0.612553889667980237953, tol);
  BOOST_CHECK_CLOSE(xq[12], -0.549467125095128202076, tol);
  BOOST_CHECK_CLOSE(xq[13], -0.483075801686178712909, tol);
  BOOST_CHECK_CLOSE(xq[14], -0.413779204371605001525, tol);
  BOOST_CHECK_CLOSE(xq[15], -0.341994090825758473007, tol);
  BOOST_CHECK_CLOSE(xq[16], -0.268152185007253681141, tol);
  BOOST_CHECK_CLOSE(xq[17], -0.192697580701371099716, tol);
  BOOST_CHECK_CLOSE(xq[18], -0.116084070675255208483, tol);
  BOOST_CHECK_CLOSE(xq[19], -0.038772417506050821933, tol);
  BOOST_CHECK_CLOSE(xq[20],  0.038772417506050821933, tol);
  BOOST_CHECK_CLOSE(xq[21],  0.116084070675255208483, tol);
  BOOST_CHECK_CLOSE(xq[22],  0.192697580701371099716, tol);
  BOOST_CHECK_CLOSE(xq[23],  0.268152185007253681141, tol);
  BOOST_CHECK_CLOSE(xq[24],  0.341994090825758473007, tol);
  BOOST_CHECK_CLOSE(xq[25],  0.413779204371605001525, tol);
  BOOST_CHECK_CLOSE(xq[26],  0.483075801686178712909, tol);
  BOOST_CHECK_CLOSE(xq[27],  0.549467125095128202076, tol);
  BOOST_CHECK_CLOSE(xq[28],  0.612553889667980237953, tol);
  BOOST_CHECK_CLOSE(xq[29],  0.671956684614179548379, tol);
  BOOST_CHECK_CLOSE(xq[30],  0.727318255189927103281, tol);
  BOOST_CHECK_CLOSE(xq[31],  0.778305651426519387695, tol);
  BOOST_CHECK_CLOSE(xq[32],  0.824612230833311663196, tol);
  BOOST_CHECK_CLOSE(xq[33],  0.865959503212259503821, tol);
  BOOST_CHECK_CLOSE(xq[34],  0.902098806968874296728, tol);
  BOOST_CHECK_CLOSE(xq[35],  0.932812808278676533361, tol);
  BOOST_CHECK_CLOSE(xq[36],  0.957916819213791655805, tol);
  BOOST_CHECK_CLOSE(xq[37],  0.977259949983774262663, tol);
  BOOST_CHECK_CLOSE(xq[38],  0.990726238699457006453, tol);
  BOOST_CHECK_CLOSE(xq[39],  0.998237709710559200350, tol);


  tol = 1e-11;
  BOOST_CHECK_CLOSE(wq[ 0], 0.004521277098533191258, tol);
  BOOST_CHECK_CLOSE(wq[ 1], 0.010498284531152813615, tol);
  BOOST_CHECK_CLOSE(wq[ 2], 0.016421058381907888713, tol);
  BOOST_CHECK_CLOSE(wq[ 3], 0.022245849194166957262, tol);
  BOOST_CHECK_CLOSE(wq[ 4], 0.027937006980023401098, tol);
  BOOST_CHECK_CLOSE(wq[ 5], 0.033460195282547847393, tol);
  BOOST_CHECK_CLOSE(wq[ 6], 0.038782167974472017640, tol);
  BOOST_CHECK_CLOSE(wq[ 7], 0.043870908185673271992, tol);
  BOOST_CHECK_CLOSE(wq[ 8], 0.048695807635072232061, tol);
  BOOST_CHECK_CLOSE(wq[ 9], 0.053227846983936824355, tol);
  BOOST_CHECK_CLOSE(wq[10], 0.057439769099391551367, tol);
  BOOST_CHECK_CLOSE(wq[11], 0.061306242492928939167, tol);
  BOOST_CHECK_CLOSE(wq[12], 0.064804013456601038075, tol);
  BOOST_CHECK_CLOSE(wq[13], 0.067912045815233903826, tol);
  BOOST_CHECK_CLOSE(wq[14], 0.070611647391286779695, tol);
  BOOST_CHECK_CLOSE(wq[15], 0.072886582395804059061, tol);
  BOOST_CHECK_CLOSE(wq[16], 0.074723169057968264200, tol);
  BOOST_CHECK_CLOSE(wq[17], 0.076110361900626242372, tol);
  BOOST_CHECK_CLOSE(wq[18], 0.077039818164247965588, tol);
  BOOST_CHECK_CLOSE(wq[19], 0.077505947978424811264, tol);
  BOOST_CHECK_CLOSE(wq[20], 0.077505947978424811264, tol);
  BOOST_CHECK_CLOSE(wq[21], 0.077039818164247965588, tol);
  BOOST_CHECK_CLOSE(wq[22], 0.076110361900626242372, tol);
  BOOST_CHECK_CLOSE(wq[23], 0.074723169057968264200, tol);
  BOOST_CHECK_CLOSE(wq[24], 0.072886582395804059061, tol);
  BOOST_CHECK_CLOSE(wq[25], 0.070611647391286779695, tol);
  BOOST_CHECK_CLOSE(wq[26], 0.067912045815233903826, tol);
  BOOST_CHECK_CLOSE(wq[27], 0.064804013456601038075, tol);
  BOOST_CHECK_CLOSE(wq[28], 0.061306242492928939167, tol);
  BOOST_CHECK_CLOSE(wq[29], 0.057439769099391551367, tol);
  BOOST_CHECK_CLOSE(wq[30], 0.053227846983936824355, tol);
  BOOST_CHECK_CLOSE(wq[31], 0.048695807635072232061, tol);
  BOOST_CHECK_CLOSE(wq[32], 0.043870908185673271992, tol);
  BOOST_CHECK_CLOSE(wq[33], 0.038782167974472017640, tol);
  BOOST_CHECK_CLOSE(wq[34], 0.033460195282547847393, tol);
  BOOST_CHECK_CLOSE(wq[35], 0.027937006980023401098, tol);
  BOOST_CHECK_CLOSE(wq[36], 0.022245849194166957262, tol);
  BOOST_CHECK_CLOSE(wq[37], 0.016421058381907888713, tol);
  BOOST_CHECK_CLOSE(wq[38], 0.010498284531152813615, tol);
  BOOST_CHECK_CLOSE(wq[39], 0.004521277098533191258, tol);

}


BOOST_AUTO_TEST_CASE( test_legendre0MassMatrix )
{
  int ierr;
  double MM[9];
  const double ptol = 1e-13;
  const double stol = 1e-16;

  ierr = legendre0MassMatrix(3, MM);
  BOOST_REQUIRE( ierr==0 );

  BOOST_CHECK_CLOSE(MM[0],  16.0/15.0 , ptol);
  BOOST_CHECK_SMALL(MM[1],              stol);
  BOOST_CHECK_CLOSE(MM[2], -32.0/105.0, ptol);
  BOOST_CHECK_SMALL(MM[3],              stol);
  BOOST_CHECK_CLOSE(MM[4],  16.0/105.0, ptol);
  BOOST_CHECK_SMALL(MM[5],              stol);
  BOOST_CHECK_CLOSE(MM[6], -32.0/105.0, ptol);
  BOOST_CHECK_SMALL(MM[7],              stol);
  BOOST_CHECK_CLOSE(MM[8],  16.0/105.0, ptol); 
}


BOOST_AUTO_TEST_CASE( test_legendre0InverseMassMatrix )
{
  int ierr;
  double iMM[9];
  const double ptol = 1e-13;
  const double stol = 1e-16;

  ierr = legendre0InverseMassMatrix(3, iMM);
  BOOST_REQUIRE( ierr==0 );

  BOOST_CHECK_CLOSE(iMM[0], 2.1875 , ptol);
  BOOST_CHECK_SMALL(iMM[1],          stol);
  BOOST_CHECK_CLOSE(iMM[2], 4.3750 , ptol);
  BOOST_CHECK_SMALL(iMM[3],          stol);
  BOOST_CHECK_CLOSE(iMM[4], 6.5625 , ptol);
  BOOST_CHECK_SMALL(iMM[5],          stol);
  BOOST_CHECK_CLOSE(iMM[6], 4.3750 , ptol);
  BOOST_CHECK_SMALL(iMM[7],          stol);
  BOOST_CHECK_CLOSE(iMM[8], 15.3125, ptol); 
}
