#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <suzerain/bsmbsm.h>
#include <boost/test/included/unit_test.hpp>

#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

BOOST_AUTO_TEST_SUITE( permutation_vector )

BOOST_AUTO_TEST_CASE( q_S5n9 )
{
    // Compare against Octave-based results for
    //      q    = mod(0:44,5)*9 + floor((0:44)/5)
    // for case when S = 5, n = 9
    const int S = 5, n = 9;
    BOOST_CHECK_EQUAL( 0, suzerain_bsmbsm_q(S, n,  0));
    BOOST_CHECK_EQUAL( 9, suzerain_bsmbsm_q(S, n,  1));
    BOOST_CHECK_EQUAL(18, suzerain_bsmbsm_q(S, n,  2));
    BOOST_CHECK_EQUAL(27, suzerain_bsmbsm_q(S, n,  3));
    BOOST_CHECK_EQUAL(36, suzerain_bsmbsm_q(S, n,  4));
    BOOST_CHECK_EQUAL( 1, suzerain_bsmbsm_q(S, n,  5));
    BOOST_CHECK_EQUAL(10, suzerain_bsmbsm_q(S, n,  6));
    BOOST_CHECK_EQUAL(19, suzerain_bsmbsm_q(S, n,  7));
    BOOST_CHECK_EQUAL(28, suzerain_bsmbsm_q(S, n,  8));
    BOOST_CHECK_EQUAL(37, suzerain_bsmbsm_q(S, n,  9));
    BOOST_CHECK_EQUAL( 2, suzerain_bsmbsm_q(S, n, 10));
    BOOST_CHECK_EQUAL(11, suzerain_bsmbsm_q(S, n, 11));
    BOOST_CHECK_EQUAL(20, suzerain_bsmbsm_q(S, n, 12));
    BOOST_CHECK_EQUAL(29, suzerain_bsmbsm_q(S, n, 13));
    BOOST_CHECK_EQUAL(38, suzerain_bsmbsm_q(S, n, 14));
    BOOST_CHECK_EQUAL( 3, suzerain_bsmbsm_q(S, n, 15));
    BOOST_CHECK_EQUAL(12, suzerain_bsmbsm_q(S, n, 16));
    BOOST_CHECK_EQUAL(21, suzerain_bsmbsm_q(S, n, 17));
    BOOST_CHECK_EQUAL(30, suzerain_bsmbsm_q(S, n, 18));
    BOOST_CHECK_EQUAL(39, suzerain_bsmbsm_q(S, n, 19));
    BOOST_CHECK_EQUAL( 4, suzerain_bsmbsm_q(S, n, 20));
    BOOST_CHECK_EQUAL(13, suzerain_bsmbsm_q(S, n, 21));
    BOOST_CHECK_EQUAL(22, suzerain_bsmbsm_q(S, n, 22));
    BOOST_CHECK_EQUAL(31, suzerain_bsmbsm_q(S, n, 23));
    BOOST_CHECK_EQUAL(40, suzerain_bsmbsm_q(S, n, 24));
    BOOST_CHECK_EQUAL( 5, suzerain_bsmbsm_q(S, n, 25));
    BOOST_CHECK_EQUAL(14, suzerain_bsmbsm_q(S, n, 26));
    BOOST_CHECK_EQUAL(23, suzerain_bsmbsm_q(S, n, 27));
    BOOST_CHECK_EQUAL(32, suzerain_bsmbsm_q(S, n, 28));
    BOOST_CHECK_EQUAL(41, suzerain_bsmbsm_q(S, n, 29));
    BOOST_CHECK_EQUAL( 6, suzerain_bsmbsm_q(S, n, 30));
    BOOST_CHECK_EQUAL(15, suzerain_bsmbsm_q(S, n, 31));
    BOOST_CHECK_EQUAL(24, suzerain_bsmbsm_q(S, n, 32));
    BOOST_CHECK_EQUAL(33, suzerain_bsmbsm_q(S, n, 33));
    BOOST_CHECK_EQUAL(42, suzerain_bsmbsm_q(S, n, 34));
    BOOST_CHECK_EQUAL( 7, suzerain_bsmbsm_q(S, n, 35));
    BOOST_CHECK_EQUAL(16, suzerain_bsmbsm_q(S, n, 36));
    BOOST_CHECK_EQUAL(25, suzerain_bsmbsm_q(S, n, 37));
    BOOST_CHECK_EQUAL(34, suzerain_bsmbsm_q(S, n, 38));
    BOOST_CHECK_EQUAL(43, suzerain_bsmbsm_q(S, n, 39));
    BOOST_CHECK_EQUAL( 8, suzerain_bsmbsm_q(S, n, 40));
    BOOST_CHECK_EQUAL(17, suzerain_bsmbsm_q(S, n, 41));
    BOOST_CHECK_EQUAL(26, suzerain_bsmbsm_q(S, n, 42));
    BOOST_CHECK_EQUAL(35, suzerain_bsmbsm_q(S, n, 43));
    BOOST_CHECK_EQUAL(44, suzerain_bsmbsm_q(S, n, 44));
}

BOOST_AUTO_TEST_CASE( qinv_S5n9 )
{
    // Compare against Octave-based results for
    //      qinv = mod(0:44,9)*5 + floor((0:44)/9);
    // for case when S = 5, n = 9
    const int S = 5, n = 9;
    BOOST_CHECK_EQUAL( 0, suzerain_bsmbsm_qinv(S, n,  0));
    BOOST_CHECK_EQUAL( 5, suzerain_bsmbsm_qinv(S, n,  1));
    BOOST_CHECK_EQUAL(10, suzerain_bsmbsm_qinv(S, n,  2));
    BOOST_CHECK_EQUAL(15, suzerain_bsmbsm_qinv(S, n,  3));
    BOOST_CHECK_EQUAL(20, suzerain_bsmbsm_qinv(S, n,  4));
    BOOST_CHECK_EQUAL(25, suzerain_bsmbsm_qinv(S, n,  5));
    BOOST_CHECK_EQUAL(30, suzerain_bsmbsm_qinv(S, n,  6));
    BOOST_CHECK_EQUAL(35, suzerain_bsmbsm_qinv(S, n,  7));
    BOOST_CHECK_EQUAL(40, suzerain_bsmbsm_qinv(S, n,  8));
    BOOST_CHECK_EQUAL( 1, suzerain_bsmbsm_qinv(S, n,  9));
    BOOST_CHECK_EQUAL( 6, suzerain_bsmbsm_qinv(S, n, 10));
    BOOST_CHECK_EQUAL(11, suzerain_bsmbsm_qinv(S, n, 11));
    BOOST_CHECK_EQUAL(16, suzerain_bsmbsm_qinv(S, n, 12));
    BOOST_CHECK_EQUAL(21, suzerain_bsmbsm_qinv(S, n, 13));
    BOOST_CHECK_EQUAL(26, suzerain_bsmbsm_qinv(S, n, 14));
    BOOST_CHECK_EQUAL(31, suzerain_bsmbsm_qinv(S, n, 15));
    BOOST_CHECK_EQUAL(36, suzerain_bsmbsm_qinv(S, n, 16));
    BOOST_CHECK_EQUAL(41, suzerain_bsmbsm_qinv(S, n, 17));
    BOOST_CHECK_EQUAL( 2, suzerain_bsmbsm_qinv(S, n, 18));
    BOOST_CHECK_EQUAL( 7, suzerain_bsmbsm_qinv(S, n, 19));
    BOOST_CHECK_EQUAL(12, suzerain_bsmbsm_qinv(S, n, 20));
    BOOST_CHECK_EQUAL(17, suzerain_bsmbsm_qinv(S, n, 21));
    BOOST_CHECK_EQUAL(22, suzerain_bsmbsm_qinv(S, n, 22));
    BOOST_CHECK_EQUAL(27, suzerain_bsmbsm_qinv(S, n, 23));
    BOOST_CHECK_EQUAL(32, suzerain_bsmbsm_qinv(S, n, 24));
    BOOST_CHECK_EQUAL(37, suzerain_bsmbsm_qinv(S, n, 25));
    BOOST_CHECK_EQUAL(42, suzerain_bsmbsm_qinv(S, n, 26));
    BOOST_CHECK_EQUAL( 3, suzerain_bsmbsm_qinv(S, n, 27));
    BOOST_CHECK_EQUAL( 8, suzerain_bsmbsm_qinv(S, n, 28));
    BOOST_CHECK_EQUAL(13, suzerain_bsmbsm_qinv(S, n, 29));
    BOOST_CHECK_EQUAL(18, suzerain_bsmbsm_qinv(S, n, 30));
    BOOST_CHECK_EQUAL(23, suzerain_bsmbsm_qinv(S, n, 31));
    BOOST_CHECK_EQUAL(28, suzerain_bsmbsm_qinv(S, n, 32));
    BOOST_CHECK_EQUAL(33, suzerain_bsmbsm_qinv(S, n, 33));
    BOOST_CHECK_EQUAL(38, suzerain_bsmbsm_qinv(S, n, 34));
    BOOST_CHECK_EQUAL(43, suzerain_bsmbsm_qinv(S, n, 35));
    BOOST_CHECK_EQUAL( 4, suzerain_bsmbsm_qinv(S, n, 36));
    BOOST_CHECK_EQUAL( 9, suzerain_bsmbsm_qinv(S, n, 37));
    BOOST_CHECK_EQUAL(14, suzerain_bsmbsm_qinv(S, n, 38));
    BOOST_CHECK_EQUAL(19, suzerain_bsmbsm_qinv(S, n, 39));
    BOOST_CHECK_EQUAL(24, suzerain_bsmbsm_qinv(S, n, 40));
    BOOST_CHECK_EQUAL(29, suzerain_bsmbsm_qinv(S, n, 41));
    BOOST_CHECK_EQUAL(34, suzerain_bsmbsm_qinv(S, n, 42));
    BOOST_CHECK_EQUAL(39, suzerain_bsmbsm_qinv(S, n, 43));
    BOOST_CHECK_EQUAL(44, suzerain_bsmbsm_qinv(S, n, 44));
}

BOOST_AUTO_TEST_CASE( identity_relation )
{
   // Checks that qinv inverts q for a large operating S, n
   for (int S = 1; S < 5; ++S) {
      for (int n = 1; n < 10; ++n) {
         for (int i = 0; i < S*n; ++i) {
            BOOST_CHECK_EQUAL(i,
               suzerain_bsmbsm_qinv(S, n,
                  suzerain_bsmbsm_q(S, n, i)));
         }
      }
   }
}

BOOST_AUTO_TEST_SUITE_END()
