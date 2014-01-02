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

#include <suzerain/iterator.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>

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
