#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <suzerain/kahan.h>
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( correctness_double_3 )
{
    const double eps = std::numeric_limits<double>::epsilon();
    const double d[3] = { 1, eps/2, eps/2 };
    const double expected = double(1) + eps;
    BOOST_CHECK_NE(expected, d[0] + d[1] + d[2]);
    BOOST_CHECK_EQUAL(expected, suzerain_kahan(d, sizeof(d)/sizeof(d[0])));
    BOOST_CHECK_EQUAL(expected, suzerain_kahan3(d[0], d[1], d[2]));
}

BOOST_AUTO_TEST_CASE( correctness_double_4 )
{
    const double eps = std::numeric_limits<double>::epsilon();
    const double d[4] = { 1, eps/2, eps/4, eps/4 };
    const double expected = double(1) + eps;
    BOOST_CHECK_NE(expected, d[0] + d[1] + d[2] + d[3]);
    BOOST_CHECK_EQUAL(expected, suzerain_kahan(d, sizeof(d)/sizeof(d[0])));
    BOOST_CHECK_EQUAL(expected, suzerain_kahan4(d[0], d[1], d[2], d[3]));
}

BOOST_AUTO_TEST_CASE( correctness_double_5 )
{
    const double eps = std::numeric_limits<double>::epsilon();
    const double d[5] = { 1, eps/2, eps/4, eps/8, eps/8 };
    const double expected = double(1) + eps;
    BOOST_CHECK_NE(expected, d[0] + d[1] + d[2] + d[3] + d[4]);
    BOOST_CHECK_EQUAL(expected, suzerain_kahan(d, sizeof(d)/sizeof(d[0])));
    BOOST_CHECK_EQUAL(expected, suzerain_kahan5(d[0], d[1], d[2], d[3], d[4]));
}
