#include <suzerain/config.h>
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <suzerain/state.hpp>
#include <suzerain/timestepper.hpp>
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE( RealState )

BOOST_AUTO_TEST_CASE( declare_pointer )
{
    suzerain::RealState<double> *state = NULL;
}

BOOST_AUTO_TEST_CASE( constructor )
{
    const suzerain::RealState<double> foo(1, 2, 3);
    BOOST_CHECK_EQUAL(foo.variable_count, 1);
    BOOST_CHECK_EQUAL(foo.vector_length, 2);
    BOOST_CHECK_EQUAL(foo.vector_count, 3);
}

BOOST_AUTO_TEST_SUITE_END()
