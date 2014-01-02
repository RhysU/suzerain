//--------------------------------------------------------------------------
//
// Copyright (C) 2012-2014 Rhys Ulerich
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

#include <suzerain/kahan.h>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>

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
