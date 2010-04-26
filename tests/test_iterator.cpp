#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <suzerain/iterator.hpp>
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE( infinite_constant )

using suzerain::iterator::infinite_constant;
using suzerain::iterator::make_infinite_constant;

BOOST_AUTO_TEST_CASE( boolean )
{
    infinite_constant<bool> i_true = make_infinite_constant(true);
    BOOST_CHECK_EQUAL(true, *i_true);
    BOOST_CHECK_EQUAL(true, *i_true++);
    BOOST_CHECK_EQUAL(true, *i_true++);
    BOOST_CHECK_EQUAL(true, *i_true++);
    BOOST_CHECK_EQUAL(true, *(++i_true));
    BOOST_CHECK_EQUAL(true, *(++i_true));
    BOOST_CHECK_EQUAL(true, *(++i_true));

    infinite_constant<bool> i_false = make_infinite_constant(false);
    BOOST_CHECK_EQUAL(false, *i_false);
    BOOST_CHECK_EQUAL(false, *i_false++);
    BOOST_CHECK_EQUAL(false, *i_false++);
    BOOST_CHECK_EQUAL(false, *i_false++);
    BOOST_CHECK_EQUAL(false, *(++i_false));
    BOOST_CHECK_EQUAL(false, *(++i_false));
    BOOST_CHECK_EQUAL(false, *(++i_false));

    BOOST_CHECK_EQUAL(i_true, i_true);
    BOOST_CHECK_EQUAL(i_true, make_infinite_constant(true));
    BOOST_CHECK_EQUAL(i_false, i_false);
    BOOST_CHECK_EQUAL(i_false, make_infinite_constant(false));
    BOOST_CHECK_NE(i_true, i_false);
    BOOST_CHECK_NE(i_false, i_true);
}

BOOST_AUTO_TEST_CASE( integer )
{
    infinite_constant<int> ic = make_infinite_constant(4);
    BOOST_CHECK_EQUAL(4, *ic);
    BOOST_CHECK_EQUAL(4, *ic++);
    BOOST_CHECK_EQUAL(4, *ic++);
    BOOST_CHECK_EQUAL(4, *ic++);
    BOOST_CHECK_EQUAL(4, *(++ic));
    BOOST_CHECK_EQUAL(4, *(++ic));
    BOOST_CHECK_EQUAL(4, *(++ic));

    BOOST_CHECK_EQUAL(ic, ic);
    BOOST_CHECK_EQUAL(ic, make_infinite_constant(4));
    BOOST_CHECK_EQUAL(make_infinite_constant(4), ic);
    BOOST_CHECK_NE(ic, make_infinite_constant(5));
}

BOOST_AUTO_TEST_SUITE_END()
