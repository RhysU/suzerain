//  Program to test the random number streams file:    rngstream.cpp
//  From http://www.iro.umontreal.ca/~lecuyer/myftp/streams00/c++/testRngStream.cpp
//  Modified to use Boost Test for error reporting.
//  Original authors retain their copyright, of course.

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <suzerain/rngstream.hpp>
#define BOOST_TEST_MAIN
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

using namespace std;
using suzerain::rngstream;

BOOST_AUTO_TEST_CASE( main_test )
{
   double sum;
   int  i;
   rngstream g1 ("g1");
   rngstream g2 ("g2");
   rngstream g3 ("g3");

   sum = g2.RandU01 () + g3.RandU01 ();

   g1.AdvanceState (5, 3);
   sum += g1.RandU01 ();

   g1.ResetStartStream ();
   for (i = 0;  i < 35; i++)
      g1.AdvanceState (0,1);
   sum += g1.RandU01 ();

   g1.ResetStartStream ();
   long sumi = 0;
   for (i = 0;  i < 35; i++)
      sumi += g1.RandInt (1, 10);
   sum += sumi / 100.0;

   double sum3 = 0.0;
   for (i = 0;  i < 100;  i++) {
      sum3 += g3.RandU01 ();
   }
   sum += sum3 / 10.0;

   g3.ResetStartStream ();
   for (i=1; i<=5; i++)
      sum += g3.RandU01 ();

   for (i=0; i<4; i++)
      g3.ResetNextSubstream ();
   for (i=0; i<5; i++)
      sum += g3.RandU01 ();

   g3.ResetStartSubstream ();
   for (i=0; i<5; i++)
      sum += g3.RandU01 ();

   g2.ResetNextSubstream ();
   sum3 = 0.0;
   for (i=1; i<=100000; i++)
      sum3 += g2.RandU01 ();
   sum += sum3 / 10000.0;

   g3.SetAntithetic (true);
   sum3 = 0.0;
   for (i=1; i<=100000; i++)
      sum3 += g3.RandU01 ();
   sum += sum3 / 10000.0;

   unsigned long germe[6] = { 1, 1, 1, 1, 1, 1 };
   rngstream::SetPackageSeed (germe);

   rngstream gar[4] = { "Poisson", "Laplace", "Galois", "Cantor" };
   for  (i = 0; i < 4; i++)
      sum += gar[i].RandU01 ();

   gar[2].AdvanceState (-127, 0);
   sum += gar[2].RandU01 ();

   gar[2].ResetNextSubstream ();
   gar[2].IncreasedPrecis (true);
   sum3 = 0.0;
   for  (i = 0; i < 100000; i++)
      sum3 += gar[2].RandU01 ();
   sum += sum3 / 10000.0;

   gar[2].SetAntithetic (true);
   sum3 = 0.0;
   for  (i = 0; i < 100000; i++)
      sum3 += gar[2].RandU01 ();
   sum += sum3 / 10000.0;
   gar[2].SetAntithetic (false);

   gar[2].IncreasedPrecis (false);
   for  (i = 0; i < 4; i++)
      sum += gar[i].RandU01 ();

   // Original output:
   // cout.precision(14);
   // cout << "------------------------------------------\n";
   // cout << "This program should print  39.697547445251\n";
   // cout << "Actual test result =       " << sum << "\n\n";
   // Boost.Testified with precision requirement from cout.precision above.
   BOOST_REQUIRE_SMALL(std::abs(sum - 39.697547445251), 1e-13);
}

static void test_helper_RandN01(bool increasedPrecis)
{
   rngstream s;
   s.IncreasedPrecis(increasedPrecis);

   // Check Berry--Esseen theorem holds for bound on CLT estimate of mean.
   // Constants are C = 0.4784 and rho = E[|N(0,1)|^3] = 2*sqrt(2/pi) for
   // | F_n(x) - phi(x) | <= C rho sigma^-3 n^(-1/2).  The left hand side
   // is the absolute error between the observed CDF at some point x
   // and the standard normal CDF.  Choose x = 1/2 so phi(1/2) = 1/2.
   // See http://en.wikipedia.org/wiki/Berry%E2%80%93Ess%C3%A9en_theorem

   const unsigned incrsamp[] = { 10, 90, 900, 9000, 90000, 900000 };
   unsigned nneg = 0;
   for (unsigned i = 0; i < sizeof(incrsamp)/sizeof(incrsamp[0]); ++i) {
      for (unsigned j = 0; j < incrsamp[i]; ++j) {
         nneg += (s.RandN01() < 0);
      }

      double n = 0;
      for (unsigned j = 0; j <= i; ++j) n += incrsamp[j];

      const double absdiff = abs(nneg / n - 0.5);
      const double bound = 0.4784 * 2 * sqrt(M_2_PI) * 1 / sqrt(n);
      BOOST_TEST_MESSAGE("n = " << n << " gives |Fn(0.5) - 0.5| = " << absdiff
                         << " < " << bound << " ?");
      BOOST_CHECK_LE(absdiff, bound);
   }

   // TODO Test convergence rate of sample variance
   // See J. Austral. Math. Soc. 25 (Series A) (1978), 250-256 by P. Hall
}

BOOST_AUTO_TEST_CASE( RandN01_low_precision )
{
   test_helper_RandN01(false);
}

BOOST_AUTO_TEST_CASE( RandN01_high_precision )
{
   test_helper_RandN01(true);
}
