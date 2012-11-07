//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <suzerain/functional.hpp>
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE( product )

using suzerain::functional::product;

BOOST_AUTO_TEST_CASE( integer )
{
    suzerain::array<int,3> a = {{ 2, 3, 4 }};
    BOOST_CHECK_EQUAL(24,     product(a.begin(), a.end()));
    BOOST_CHECK_EQUAL(24,     product(a));
    BOOST_CHECK_EQUAL(24*5.0, product(a.begin(), a.end(), 5.0));

    suzerain::array<int,1> b = {{ 5 }};
    BOOST_CHECK_EQUAL(5,     product(b.begin(), b.end()));
    BOOST_CHECK_EQUAL(5,     product(b));
    BOOST_CHECK_EQUAL(5*6.0, product(b.begin(), b.end(), 6.0));

    suzerain::array<int,0> c;
    BOOST_CHECK_EQUAL(1,   product(c.begin(), c.end()));
    BOOST_CHECK_EQUAL(7.0, product(c.begin(), c.end(), 7.0));
}

BOOST_AUTO_TEST_CASE( floating_point )
{
    suzerain::array<double,3> a = {{ 2.0, 3.0, 4.0 }};
    BOOST_CHECK_EQUAL(24.0,     product(a.begin(), a.end()));
    BOOST_CHECK_EQUAL(24.0,     product(a));
    BOOST_CHECK_EQUAL(24.0*5.0, product(a.begin(), a.end(), 5.0));

    suzerain::array<double,1> b = {{ 5.0 }};
    BOOST_CHECK_EQUAL(5.0,     product(b.begin(), b.end()));
    BOOST_CHECK_EQUAL(5.0,     product(b));
    BOOST_CHECK_EQUAL(5.0*6.0, product(b.begin(), b.end(), 6.0));

    suzerain::array<double,0> c;
    BOOST_CHECK_EQUAL(1.0, product(c.begin(), c.end()));
    BOOST_CHECK_EQUAL(1.0, product(c));
    BOOST_CHECK_EQUAL(7.0, product(c.begin(), c.end(), 7.0));
}

BOOST_AUTO_TEST_SUITE_END()
