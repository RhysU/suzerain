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

#include <suzerain/mpl.hpp>

#define BOOST_TEST_MAIN
#include <boost/mpl/vector_c.hpp>
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>

BOOST_AUTO_TEST_SUITE( sequence_array )

using boost::mpl::vector_c;
using suzerain::mpl::sequence_array;

BOOST_AUTO_TEST_CASE( integer )
{
    sequence_array< boost::mpl::vector_c<int,2,4,6> > a;
    BOOST_CHECK_EQUAL(a.size(), 3U);

    BOOST_CHECK_EQUAL(a[0], 2);
    BOOST_CHECK_EQUAL(a[1], 4);
    BOOST_CHECK_EQUAL(a[2], 6);

    suzerain::array<int, 3> b = {{ 2, 4, 6 }};

    BOOST_CHECK( b == a );
    BOOST_CHECK( a == b );

}

BOOST_AUTO_TEST_CASE( boolean )
{
    sequence_array< boost::mpl::vector_c<bool,true,false> > a;
    BOOST_CHECK_EQUAL(a.size(), 2U);

    BOOST_CHECK_EQUAL(a[0], true);
    BOOST_CHECK_EQUAL(a[1], false);

    suzerain::array<bool, 2> b = {{ true, false }};

    BOOST_CHECK( b == a );
    BOOST_CHECK( a == b );
}

BOOST_AUTO_TEST_CASE( empty )
{
    sequence_array< boost::mpl::vector_c<double> > a;
    BOOST_CHECK_EQUAL(a.size(), 0U);

    suzerain::array<double, 0> b;

    BOOST_CHECK( b == a );
    BOOST_CHECK( a == b );
}

BOOST_AUTO_TEST_SUITE_END()
