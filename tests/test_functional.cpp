#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <suzerain/functional.hpp>
#include <boost/test/included/unit_test.hpp>

#pragma warning(disable:383)

BOOST_AUTO_TEST_SUITE( product )

using suzerain::functional::product;

BOOST_AUTO_TEST_CASE( integer )
{
    boost::array<int,3> a = { 2, 3, 4 };
    BOOST_CHECK_EQUAL(24,     product(a.begin(), a.end()));
    BOOST_CHECK_EQUAL(24*5.0, product(a.begin(), a.end(), 5.0));

    boost::array<int,1> b = { 5 };
    BOOST_CHECK_EQUAL(5,     product(b.begin(), b.end()));
    BOOST_CHECK_EQUAL(5*6.0, product(b.begin(), b.end(), 6.0));

    boost::array<int,0> c = { };
    BOOST_CHECK_EQUAL(1,   product(c.begin(), c.end()));
    BOOST_CHECK_EQUAL(7.0, product(c.begin(), c.end(), 7.0));
}

BOOST_AUTO_TEST_CASE( floating_point )
{
    boost::array<double,3> a = { 2.0, 3.0, 4.0 };
    BOOST_CHECK_EQUAL(24.0,     product(a.begin(), a.end()));
    BOOST_CHECK_EQUAL(24.0*5.0, product(a.begin(), a.end(), 5.0));

    boost::array<double,1> b = { 5.0 };
    BOOST_CHECK_EQUAL(5.0,     product(b.begin(), b.end()));
    BOOST_CHECK_EQUAL(5.0*6.0, product(b.begin(), b.end(), 6.0));

    boost::array<double,0> c = { };
    BOOST_CHECK_EQUAL(1.0, product(c.begin(), c.end()));
    BOOST_CHECK_EQUAL(7.0, product(c.begin(), c.end(), 7.0));
}

BOOST_AUTO_TEST_SUITE_END()
