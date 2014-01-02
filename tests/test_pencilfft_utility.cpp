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

#include <suzerain/pencilfft.hpp>

#define BOOST_TEST_MAIN
#include <boost/assign.hpp>
#include <boost/mpl/list_c.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/unit_test.hpp>
#include <fftw3.h>

#include <suzerain/common.hpp>

#include "test_tools.hpp"

typedef boost::mpl::list_c<int,0,1> int_zero_one_type;

using namespace suzerain;

BOOST_AUTO_TEST_SUITE( increment )

using pencilfft::internal::increment;

BOOST_AUTO_TEST_CASE( increment_1d_degenerate )
{
    const int n        = 1;
    int       index[n] = { 0 };
    const int shape[n] = { 1 };

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_1d_normal )
{
    const int n        = 1;
    int       index[n] = { 0 };
    const int shape[n] = { 3 };

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_1d_normal_unsigned )
{
    const unsigned int n        = 1;
    unsigned int       index[n] = { 0 };
    const unsigned int shape[n] = { 3 };

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_degenerate_first )
{
    const int           n     = 2;
    array<int,n> index_array       = {{ 0, 0 }};
    const array<int,n> shape_array = {{ 3, 1 }};
    array<int,n>::iterator       index = index_array.begin();
    array<int,n>::const_iterator shape = shape_array.begin();

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(1)(0));

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(2)(0));

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_degenerate_second )
{
    const int                            n     = 2;
    array<int, n> index_array            = {{ 0, 0 }};
    array<std::size_t,n> shape_array     = {{ 1, 3 }};
    array<int, n>::iterator index        = index_array.begin();
    array<std::size_t,n>::iterator shape = shape_array.begin();

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(1));

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(2));

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_degenerate_all )
{
    const int  n        = 2;
    int        index[n] = { 0, 0 };
    const long shape[n] = { 1, 1 };

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_normal )
{
    const int    n        = 2;
    unsigned int index[n] = { 0, 0 };
    unsigned int shape[n] = { 2, 3 };

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_normal_usualorder )
{
    const int    n        = 2;
    unsigned int index[n] = { 0, 0 };
    unsigned int shape[n] = { 2, 3 };
    unsigned int order[n] = { 0, 1 };

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_normal_reverseorder )
{
    const int    n        = 2;
    unsigned int index[n] = { 0, 0 };
    unsigned int shape[n] = { 2, 3 };
    unsigned int order[n] = { 1, 0 };

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), false);
}

BOOST_AUTO_TEST_CASE( increment_3d_degenerate_all )
{
    const int         n        = 3;
    signed int        index[n] = { 0, 0, 0 };
    const signed long shape[n] = { 1, 1, 1 };

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_3d_degenerate_middle )
{
    const int                            n           = 3;
    long                                 index[n]    =  { 0, 0, 0 };
    const array<std::size_t,n>           shape_array = {{ 3, 1, 3 }};
    array<std::size_t,n>::const_iterator shape       = shape_array.begin();

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_3d_normal )
{
    const int n = 3;
    std::vector<short> index_array(n, 0);
    std::vector<int>   shape_array(n, 2);
    std::vector<short>::iterator     index = index_array.begin();
    std::vector<int>::const_iterator shape = shape_array.begin();

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_3d_normal_outoforder )
{
    const int n = 3;
    array<short,3> index_array = {{ 0, 0, 0 }};
    array<int,3>   shape_array = {{ 2, 1, 3 }};
    array<long,3>  order_array = {{ 2, 0, 1 }};

    array<short,3>::iterator index = index_array.begin();
    array<int,3>::iterator   shape = shape_array.begin();
    array<long,3>::iterator  order = order_array.begin();

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), false);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( decrement )

using pencilfft::internal::decrement;

BOOST_AUTO_TEST_CASE( decrement_1d_degenerate )
{
    const int n        = 1;
    int       index[n] = { 0 };
    const int shape[n] = { 1 };

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_1d_normal )
{
    const int n        = 1;
    int       index[n] = { 2 };
    const int shape[n] = { 3 };

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_1d_normal_unsigned )
{
    const unsigned int n        = 1;
    unsigned int       index[n] = { 2 };
    const unsigned int shape[n] = { 3 };

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_2d_degenerate_first )
{
    const int n = 2;
    array<int,n> index_array = {{ 2, 0 }};
    array<int,n> shape_array = {{ 3, 1 }};
    array<int,n>::iterator index = index_array.begin();
    array<int,n>::iterator shape = shape_array.begin();

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_2d_degenerate_second )
{
    const int n = 2;
    array<int, n> index_array            = {{ 0, 2 }};
    array<std::size_t,n> shape_array     = {{ 1, 3 }};
    array<int, n>::iterator index        = index_array.begin();
    array<std::size_t,n>::iterator shape = shape_array.begin();

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_2d_degenerate_all )
{
    const int  n        = 2;
    int        index[n] = { 0, 0 };
    const long shape[n] = { 1, 1 };

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_2d_normal )
{
    const int    n        = 2;
    unsigned int index[n] = { 1, 2 };
    unsigned int shape[n] = { 2, 3 };

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_2d_normal_usualorder )
{
    const int    n        = 2;
    unsigned int index[n] = { 1, 2 };
    unsigned int shape[n] = { 2, 3 };
    unsigned int order[n] = { 0, 1 };

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), false);
}

BOOST_AUTO_TEST_CASE( decrement_2d_normal_reverseorder )
{
    const int    n        = 2;
    unsigned int index[n] = { 1, 2 };
    unsigned int shape[n] = { 2, 3 };
    unsigned int order[n] = { 1, 0 };

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), false);
}

BOOST_AUTO_TEST_CASE( decrement_3d_degenerate_all )
{
    const int         n        = 3;
    signed int        index[n] = { 0, 0, 0 };
    const signed long shape[n] = { 1, 1, 1 };

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_3d_degenerate_middle )
{
    const int            n               = 3;
    long                 index[n]        =  { 2, 0, 2 };
    array<std::size_t,n> shape_array     = {{ 3, 1, 3 }};
    array<std::size_t,n>::iterator shape = shape_array.begin();

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_3d_normal )
{
    const int n = 3;
    std::vector<short> index_array(n, 1);
    std::vector<int>   shape_array(n, 2);
    std::vector<short>::iterator     index = index_array.begin();
    std::vector<int>::const_iterator shape = shape_array.begin();

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_3d_normal_outoforder )
{
    const int n = 3;
    array<short,3> index_array = {{ 1, 0, 2 }};
    array<int,3>   shape_array = {{ 2, 1, 3 }};
    array<long,3>  order_array = {{ 2, 0, 1 }};

    array<short,3>::iterator index = index_array.begin();
    array<int,3>::iterator   shape = shape_array.begin();
    array<long,3>::iterator  order = order_array.begin();

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), false);
}

BOOST_AUTO_TEST_CASE( decrement_3d_spot_check_behavior_1 )
{
    const int n = 3;
    typedef array<int,n> array_type;
    array_type           index_array = {{ 3, 0, 1 }};
    array_type::iterator index       = index_array.begin();
    array_type           shape_array = {{ 4, 1, 2 }};
    array_type::iterator shape       = shape_array.begin();
    array_type           order_array = {{ 0, 1, 2 }};
    array_type::iterator order       = order_array.begin();

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(2)(0)(1));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(1)(0)(1));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(0)(1));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(3)(0)(0));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(2)(0)(0));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(1)(0)(0));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(0)(0));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), false);
}

BOOST_AUTO_TEST_CASE( decrement_3d_spot_check_behavior_2 )
{
    const int n = 3;
    typedef array<int,n> array_type;
    array_type           index_array = {{ 3, 0, 1 }};
    array_type::iterator index       = index_array.begin();
    array_type           shape_array = {{ 4, 1, 2 }};
    array_type::iterator shape       = shape_array.begin();
    array_type           order_array = {{ 2, 1, 0 }};
    array_type::iterator order       = order_array.begin();

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(3)(0)(0));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(2)(0)(1));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(2)(0)(0));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(1)(0)(1));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(1)(0)(0));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(0)(1));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(0)(0));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), false);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( crement )

using pencilfft::internal::crement;

BOOST_AUTO_TEST_CASE( increment_1d_degenerate )
{
    const int  n             = 1;
    int        index[n]      = { 0 };
    const bool increasing[n] = { true };
    const int  shape[n]      = { 1 };

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_1d_normal )
{
    const int  n             = 1;
    int        index[n]      = { 0 };
    const bool increasing[n] = { true };
    const int  shape[n]      = { 3 };

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_1d_normal_unsigned )
{
    const unsigned int n             = 1;
    unsigned int       index[n]      = { 0 };
    const bool         increasing[n] = { true };
    const unsigned int shape[n]      = { 3 };

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_degenerate_first )
{
    const int                               n     = 2;
    array<int,n> index_array                = {{ 0, 0 }};
    const array<int,n> increasing_array     = {{ 1, 1 }};
    const array<int,n> shape_array          = {{ 3, 1 }};
    array<int,n>::iterator index            = index_array.begin();
    array<int,n>::const_iterator increasing = increasing_array.begin();
    array<int,n>::const_iterator shape      = shape_array.begin();

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(1)(0));

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(2)(0));

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_degenerate_second )
{
    const int n = 2;
    array<int, n> index_array                = {{ 0, 0 }};
    const array<bool,n> increasing_array     = {{ true, true }};
    array<std::size_t,n> shape_array         = {{ 1, 3 }};
    array<int, n>::iterator index            = index_array.begin();
    array<bool,n>::const_iterator increasing = increasing_array.begin();
    array<std::size_t,n>::iterator shape     = shape_array.begin();

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(1));

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(2));

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_degenerate_all )
{
    const int  n             = 2;
    int        index[n]      = { 0, 0 };
    bool       increasing[n] = { true, true };
    const long shape[n]      = { 1, 1 };

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_normal )
{
    const int    n             = 2;
    unsigned int index[n]      = { 0, 0 };
    unsigned int increasing[n] = { 5, 57 };
    unsigned int shape[n]      = { 2, 3 };

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_normal_usualorder )
{
    const int    n             = 2;
    unsigned int index[n]      = { 0, 0 };
    unsigned int increasing[n] = { 5, 57 };
    unsigned int shape[n]      = { 2, 3 };
    unsigned int order[n]      = { 0, 1 };

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_normal_reverseorder )
{
    const int    n             = 2;
    unsigned int index[n]      = { 0, 0 };
    const bool   increasing[n] = { true, true };
    unsigned int shape[n]      = { 2, 3 };
    unsigned int order[n]      = { 1, 0 };

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), false);
}

BOOST_AUTO_TEST_CASE( increment_3d_degenerate_all )
{
    const int         n             = 3;
    signed int        index[n]      = { 0, 0, 0 };
    const signed long increasing[n] = { 1, 1, 1 };
    const signed long shape[n]      = { 1, 1, 1 };

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_3d_degenerate_middle )
{
    const int n = 3;
    long index[n] = { 0, 0, 0 };
    array<std::size_t,n> increasing_array = {{ 3, 1, 3 }};
    array<std::size_t,n>::const_iterator increasing = increasing_array.begin();
    const array<std::size_t,n> shape_array = {{ 3, 1, 3 }};
    array<std::size_t,n>::const_iterator shape = shape_array.begin();

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_3d_normal )
{
    const int n = 3;
    std::vector<short> index_array(n, 0);
    const std::vector<bool>  increasing_array(n, 1);
    const std::vector<int>   shape_array(n, 2);
    std::vector<short>::iterator      index      = index_array.begin();
    std::vector<bool>::const_iterator increasing = increasing_array.begin();
    std::vector<int>::const_iterator  shape      = shape_array.begin();

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_3d_normal_outoforder )
{
    const int n = 3;
    array<short,3> index_array            = {{ 0, 0, 0 }};
    const array<long,3>  increasing_array = {{ -1, -1, -1 }};
    const array<int,3>   shape_array      = {{ 2, 1, 3 }};
    const array<long,3>  order_array      = {{ 2, 0, 1 }};

    array<short,3>::iterator index = index_array.begin();
    array<long,3>::const_iterator increasing = increasing_array.begin();
    array<int,3>::const_iterator  shape = shape_array.begin();
    array<long,3>::const_iterator order = order_array.begin();

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), false);
}

BOOST_AUTO_TEST_CASE( decrement_1d_degenerate )
{
    const int n             = 1;
    int       index[n]      = { 0 };
    int       increasing[n] = { 0 };
    const int shape[n]      = { 1 };

    BOOST_REQUIRE_EQUAL(
        crement<n>(index,increasing,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_1d_normal )
{
    const int n             = 1;
    int       index[n]      = { 2 };
    int       increasing[n] = { 0 };
    const int shape[n]      = { 3 };

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_1d_normal_unsigned )
{
    const unsigned int n             = 1;
    unsigned int       index[n]      = { 2 };
    unsigned int       increasing[n] = { 0 };
    const unsigned int shape[n]      = { 3 };

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_2d_degenerate_first )
{
    const int n = 2;
    array<int,n> index_array = {{ 2, 0 }};
    const array<int,n> increasing_array = {{ 0, 0 }};
    const array<int,n> shape_array = {{ 3, 1 }};
    array<int,n>::iterator index = index_array.begin();
    array<int,n>::const_iterator increasing = increasing_array.begin();
    array<int,n>::const_iterator shape = shape_array.begin();

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_2d_degenerate_second )
{
    const int n = 2;
    array<int, n> index_array            = {{ 0, 2 }};
    const bool increasing[n]                    = { false, false };
    array<std::size_t,n> shape_array     = {{ 1, 3 }};
    array<int, n>::iterator index        = index_array.begin();
    array<std::size_t,n>::iterator shape = shape_array.begin();

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_2d_degenerate_all )
{
    const int  n             = 2;
    int        index[n]      = { 0, 0 };
    const int  increasing[n] = { 0, 0 };
    const long shape[n]      = { 1, 1 };

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_2d_normal )
{
    const int    n             = 2;
    unsigned int index[n]      = { 1, 2 };
    unsigned int increasing[n] = { 0, 0 };
    unsigned int shape[n]      = { 2, 3 };

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_2d_normal_usualorder )
{
    const int    n             = 2;
    unsigned int index[n]      = { 1, 2 };
    unsigned int increasing[n] = { 0, 0 };
    unsigned int shape[n]      = { 2, 3 };
    unsigned int order[n]      = { 0, 1 };

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), false);
}

BOOST_AUTO_TEST_CASE( decrement_2d_normal_reverseorder )
{
    const int          n             = 2;
    unsigned int       index[n]      = { 1, 2 };
    const unsigned int increasing[n] = { 0, 0 };
    const unsigned int shape[n]      = { 2, 3 };
    const unsigned int order[n]      = { 1, 0 };

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), false);
}

BOOST_AUTO_TEST_CASE( decrement_3d_degenerate_all )
{
    const int         n             = 3;
    signed int        index[n]      = { 0, 0, 0 };
    const signed long increasing[n] = { 0, 0, 0 };
    const signed long shape[n]      = { 1, 1, 1 };

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_3d_degenerate_middle )
{
    const int            n                = 3;
    long                 index[n]         =  { 2, 0, 2 };
    array<std::size_t,n> increasing_array = {{ 0, 0, 0 }};
    array<std::size_t,n> shape_array      = {{ 3, 1, 3 }};
    array<std::size_t,n>::const_iterator increasing = increasing_array.begin();
    array<std::size_t,n>::const_iterator shape = shape_array.begin();

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_3d_normal )
{
    const int n = 3;
    std::vector<short> index_array(n, 1);
    std::vector<short> increasing_array(n, 0);
    std::vector<int> shape_array(n, 2);
    std::vector<short>::iterator index = index_array.begin();
    std::vector<short>::const_iterator increasing = increasing_array.begin();
    std::vector<int>::const_iterator shape = shape_array.begin();

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_3d_normal_outoforder )
{
    const int n = 3;
    array<short,3> index_array = {{ 1, 0, 2 }};
    array<bool,3> increasing_array = {{ false, false, false }};
    array<int,3> shape_array = {{ 2, 1, 3 }};
    array<long,3> order_array = {{ 2, 0, 1 }};

    array<short,3>::iterator index = index_array.begin();
    array<bool,3>::iterator increasing = increasing_array.begin();
    array<int,3>::iterator shape = shape_array.begin();
    array<long,3>::iterator order = order_array.begin();

    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), false);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( crement_3d_spot_check_behavior_1, ASCENDING, int_zero_one_type )
{
    const int n = 3;
    typedef array<int,n> array_type;

    // Like shape = { 4, 3, 2 } with the second index frozen at zero
    const array_type           shape_array = {{ 4, 1, 2 }};
    array_type::const_iterator shape       = shape_array.begin();

    const array_type order_fast_dontcare_slow[] = {
        {{ 0, 1, 2 }},
        {{ 1, 0, 2 }},
        {{ 0, 2, 1 }}
    };
    const std::size_t n_order
        = sizeof(order_fast_dontcare_slow)/sizeof(order_fast_dontcare_slow[0]);
    const array_type order_slow_dontcare_fast[n_order] = {
        {{ 1, 2, 0 }},
        {{ 2, 1, 0 }},
        {{ 2, 0, 1 }}
    };

    for (std::size_t i = 0; i < n_order; ++i) {
        array_type       index_array      = {{ 0, 0, 1 }};
        const array_type increasing_array = {{ 1, ASCENDING::type::value, 0 }};

        array_type::iterator       index      = index_array.begin();
        array_type::const_iterator increasing = increasing_array.begin();

        BOOST_TEST_MESSAGE("Testing increasing " << increasing_array
                           << " with ordering " << order_fast_dontcare_slow[i]
                           << " from initial "  << index_array );
        array_type::const_iterator order = order_fast_dontcare_slow[i].begin();

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(1)(0)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(2)(0)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(3)(0)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(0)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(1)(0)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(2)(0)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(3)(0)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), false);
    }

    for (std::size_t i = 0; i < n_order; ++i) {
        array_type       index_array      = {{ 0, 0, 1 }};
        const array_type increasing_array = {{ 1, ASCENDING::type::value, 0 }};

        array_type::iterator       index      = index_array.begin();
        array_type::const_iterator increasing = increasing_array.begin();

        BOOST_TEST_MESSAGE("Testing increasing " << increasing_array
                           << " with ordering " << order_slow_dontcare_fast[i]
                           << " from initial "  << index_array );
        array_type::const_iterator order = order_slow_dontcare_fast[i].begin();

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(0)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(1)(0)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(1)(0)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(2)(0)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(2)(0)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(3)(0)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(3)(0)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), false);
    }

    for (std::size_t i = 0; i < n_order; ++i) {
        array_type       index_array      = {{ 3, 0, 0 }};
        const array_type increasing_array = {{ 0, ASCENDING::type::value, 1 }};

        array_type::iterator       index      = index_array.begin();
        array_type::const_iterator increasing = increasing_array.begin();

        BOOST_TEST_MESSAGE("Testing increasing " << increasing_array
                           << " with ordering " << order_fast_dontcare_slow[i]
                           << " from initial "  << index_array );
        array_type::const_iterator order = order_fast_dontcare_slow[i].begin();

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(2)(0)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(1)(0)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(0)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(3)(0)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(2)(0)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(1)(0)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(0)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), false);
    }

    for (std::size_t i = 0; i < n_order; ++i) {
        array_type       index_array      = {{ 3, 0, 0 }};
        const array_type increasing_array = {{ 0, ASCENDING::type::value, 1 }};

        array_type::iterator       index      = index_array.begin();
        array_type::const_iterator increasing = increasing_array.begin();

        BOOST_TEST_MESSAGE("Testing increasing " << increasing_array
                           << " with ordering " << order_slow_dontcare_fast[i]
                           << " from initial "  << index_array );
        array_type::const_iterator order = order_slow_dontcare_fast[i].begin();

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(3)(0)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(2)(0)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(2)(0)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(1)(0)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(1)(0)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(0)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(0)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), false);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( crement_3d_spot_check_behavior_2, ASCENDING, int_zero_one_type )
{
    const int n = 3;
    typedef array<int,n> array_type;

    // Like shape = { 4, 3, 2 } with the first index frozen at zero
    const array_type           shape_array = {{ 1, 3, 2 }};
    array_type::const_iterator shape       = shape_array.begin();

    const array_type order_dontcare_fast_slow[] = {
        {{ 0, 1, 2 }},
        {{ 1, 0, 2 }},
        {{ 1, 2, 0 }}
    };
    const std::size_t n_order
        = sizeof(order_dontcare_fast_slow)/sizeof(order_dontcare_fast_slow[0]);
    const array_type order_dontcare_slow_fast[n_order] = {
        {{ 0, 2, 1 }},
        {{ 2, 0, 1 }},
        {{ 2, 1, 0 }}
    };

    for (std::size_t i = 0; i < n_order; ++i) {
        array_type       index_array      = {{ 0, 0, 1 }};
        const array_type increasing_array = {{ ASCENDING::type::value, 1, 0 }};

        array_type::iterator       index      = index_array.begin();
        array_type::const_iterator increasing = increasing_array.begin();

        BOOST_TEST_MESSAGE("Testing increasing " << increasing_array
                           << " with ordering " << order_dontcare_fast_slow[i]
                           << " from initial "  << index_array );
        array_type::const_iterator order = order_dontcare_fast_slow[i].begin();

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(1)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(2)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(0)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(1)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(2)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), false);
    }

    for (std::size_t i = 0; i < n_order; ++i) {
        array_type       index_array      = {{ 0, 0, 1 }};
        const array_type increasing_array = {{ ASCENDING::type::value, 1, 0 }};

        array_type::iterator       index      = index_array.begin();
        array_type::const_iterator increasing = increasing_array.begin();

        BOOST_TEST_MESSAGE("Testing increasing " << increasing_array
                           << " with ordering " << order_dontcare_slow_fast[i]
                           << " from initial "  << index_array );
        array_type::const_iterator order = order_dontcare_slow_fast[i].begin();

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(0)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(1)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(1)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(2)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(2)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), false);
    }

    for (std::size_t i = 0; i < n_order; ++i) {
        array_type       index_array      = {{ 0, 2, 0 }};
        const array_type increasing_array = {{ ASCENDING::type::value, 0, 1 }};

        array_type::iterator       index      = index_array.begin();
        array_type::const_iterator increasing = increasing_array.begin();

        BOOST_TEST_MESSAGE("Testing increasing " << increasing_array
                           << " with ordering " << order_dontcare_fast_slow[i]
                           << " from initial "  << index_array );
        array_type::const_iterator order = order_dontcare_fast_slow[i].begin();

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(1)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(0)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(2)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(1)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(0)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), false);
    }

    for (std::size_t i = 0; i < n_order; ++i) {
        array_type       index_array      = {{ 0, 2, 0 }};
        const array_type increasing_array = {{ ASCENDING::type::value, 0, 1 }};

        array_type::iterator       index      = index_array.begin();
        array_type::const_iterator increasing = increasing_array.begin();

        BOOST_TEST_MESSAGE("Testing increasing " << increasing_array
                           << " with ordering " << order_dontcare_slow_fast[i]
                           << " from initial "  << index_array );
        array_type::const_iterator order = order_dontcare_slow_fast[i].begin();

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(2)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(1)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(1)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(0)(0));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), true);
        BOOST_CHECK_EQUAL(index_array, boost::assign::list_of(0)(0)(1));

        BOOST_REQUIRE_EQUAL(crement<n>(index,increasing,shape,order), false);
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( indexed_element_magnitude_comparator )

using pencilfft::internal::indexed_element_magnitude_comparator;
using pencilfft::internal::make_indexed_element_magnitude_comparator;

BOOST_AUTO_TEST_CASE( indexed_element_magnitude_comparator_operator_positive )
{
    const double data[3] = { 0.0, 1.0, 2.0 };
    indexed_element_magnitude_comparator<const double *> c(&data[0]);

    BOOST_CHECK_EQUAL(true, c(0, 1));
    BOOST_CHECK_EQUAL(true, c(1, 2));
    BOOST_CHECK_EQUAL(true, c(0, 2));

    BOOST_CHECK_EQUAL(false, c(1, 0));
    BOOST_CHECK_EQUAL(false, c(2, 1));
    BOOST_CHECK_EQUAL(false, c(2, 0));
}

BOOST_AUTO_TEST_CASE( indexed_element_magnitude_comparator_operator_mixed )
{
    const double data[3] = { 0.0, -1.0, -2.0 };
    indexed_element_magnitude_comparator<const double *> c(&data[0]);

    BOOST_CHECK_EQUAL(true, c(0, 1));
    BOOST_CHECK_EQUAL(true, c(1, 2));
    BOOST_CHECK_EQUAL(true, c(0, 2));

    BOOST_CHECK_EQUAL(false, c(1, 0));
    BOOST_CHECK_EQUAL(false, c(2, 1));
    BOOST_CHECK_EQUAL(false, c(2, 0));
}

BOOST_AUTO_TEST_CASE( indexed_element_magnitude_comparator_stable_sort )
{
    const int n = 5;
    const double data[n]    = { 0, 1, 2, 3, 4 };
    std::size_t  indices[n] = { 3, 1, 4, 0, 2 };

    std::stable_sort(indices, indices + n,
            make_indexed_element_magnitude_comparator(&data[0]));
    BOOST_CHECK_EQUAL(0, indices[0]);
    BOOST_CHECK_EQUAL(1, indices[1]);
    BOOST_CHECK_EQUAL(2, indices[2]);
    BOOST_CHECK_EQUAL(3, indices[3]);
    BOOST_CHECK_EQUAL(4, indices[4]);
}

BOOST_AUTO_TEST_SUITE_END()
