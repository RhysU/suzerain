//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Suzerain is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------

#include <suzerain/functional.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>

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
