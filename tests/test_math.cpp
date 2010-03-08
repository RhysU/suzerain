#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>
#include <suzerain/math.hpp>
#include "test_tools.hpp"

#pragma warning(push,disable:383)
BOOST_AUTO_TEST_CASE( integer_power )
{
    using suzerain::math::integer_power;

    BOOST_CHECK_EQUAL( 1u, integer_power(2u, 0u));
    BOOST_CHECK_EQUAL( 2u, integer_power(2u, 1u));
    BOOST_CHECK_EQUAL( 4u, integer_power(2u, 2u));
    BOOST_CHECK_EQUAL( 8u, integer_power(2u, 3u));
    BOOST_CHECK_EQUAL(16u, integer_power(2u, 4u));

    BOOST_CHECK_EQUAL( 1ul, integer_power(2ul, 0ul));
    BOOST_CHECK_EQUAL( 2ul, integer_power(2ul, 1ul));
    BOOST_CHECK_EQUAL( 4ul, integer_power(2ul, 2ul));
    BOOST_CHECK_EQUAL( 8ul, integer_power(2ul, 3ul));
    BOOST_CHECK_EQUAL(16ul, integer_power(2ul, 4ul));

    BOOST_CHECK_EQUAL( 1.0, integer_power(2.0, 0));
    BOOST_CHECK_EQUAL( 2.0, integer_power(2.0, 1));
    BOOST_CHECK_EQUAL( 4.0, integer_power(2.0, 2));
    BOOST_CHECK_EQUAL( 8.0, integer_power(2.0, 3));
    BOOST_CHECK_EQUAL(16.0, integer_power(2.0, 4));

    BOOST_CHECK_EQUAL( 1.0f, integer_power(2.0f, 0));
    BOOST_CHECK_EQUAL( 2.0f, integer_power(2.0f, 1));
    BOOST_CHECK_EQUAL( 4.0f, integer_power(2.0f, 2));
    BOOST_CHECK_EQUAL( 8.0f, integer_power(2.0f, 3));
    BOOST_CHECK_EQUAL(16.0f, integer_power(2.0f, 4));

    BOOST_CHECK_EQUAL( 1.0,    integer_power(2.0, -0));
    BOOST_CHECK_EQUAL( 0.5,    integer_power(2.0, -1));
    BOOST_CHECK_EQUAL( 0.25,   integer_power(2.0, -2));
    BOOST_CHECK_EQUAL( 0.125,  integer_power(2.0, -3));
    BOOST_CHECK_EQUAL( 0.0625, integer_power(2.0, -4));

    BOOST_CHECK_EQUAL( 1.0f,    integer_power(2.0f, -0));
    BOOST_CHECK_EQUAL( 0.5f,    integer_power(2.0f, -1));
    BOOST_CHECK_EQUAL( 0.25f,   integer_power(2.0f, -2));
    BOOST_CHECK_EQUAL( 0.125f,  integer_power(2.0f, -3));
    BOOST_CHECK_EQUAL( 0.0625f, integer_power(2.0f, -4));

    BOOST_CHECK_EQUAL( 1.0, integer_power(0.0, 0));
    BOOST_CHECK_EQUAL( 0.0, integer_power(0.0, 1));

    BOOST_CHECK_EQUAL( 1.0f, integer_power(0.0f, 0));
    BOOST_CHECK_EQUAL( 0.0f, integer_power(0.0f, 1));
}
#pragma warning(pop)
