//  Program to test the random number streams file:    RngStream.cpp
//  From http://www.iro.umontreal.ca/~lecuyer/myftp/streams00/c++/testRngStream.cpp
//  Modified to use Boost Test for error reporting.
//  Original authors retain their copyright, of course.



#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <suzerain/RngStream.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/included/unit_test.hpp>

using namespace std;
using suzerain::RngStream;

BOOST_AUTO_TEST_CASE( main_test )
{
   double sum;
   int  i;
   RngStream g1 ("g1");
   RngStream g2 ("g2");
   RngStream g3 ("g3");

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
   RngStream::SetPackageSeed (germe);

   RngStream gar[4] = { "Poisson", "Laplace", "Galois", "Cantor" };
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

BOOST_AUTO_TEST_CASE( RandN01_test )
{
   RngStream s;
   s.IncreasedPrecis(false);

   using namespace boost::accumulators;
   accumulator_set<double, stats<tag::mean, tag::lazy_variance> > acc;

   for (int i = 0; i < 100; ++i) acc(s.RandN01());
   const double m2 = mean(acc), v2 = variance(acc);

   for (int i = 0; i < 900; ++i) acc(s.RandN01());
   const double m3 = mean(acc), v3 = variance(acc);

   BOOST_CHECK_GE(std::abs(m2 - 0), std::abs(m3 - 0));
   BOOST_CHECK_GE(std::abs(v2 - 1), std::abs(v3 - 1));

   for (int i = 0; i < 9000; ++i) acc(s.RandN01());
   const double m4 = mean(acc), v4 = variance(acc);

   BOOST_CHECK_GE(std::abs(m3 - 0), std::abs(m4 - 0));
   BOOST_CHECK_GE(std::abs(v3 - 1), std::abs(v4 - 1));
}
