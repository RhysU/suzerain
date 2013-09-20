//--------------------------------------------------------------------------
//
// Copyright (C) 2013 Rhys Ulerich
// Copyright (C) 2013 The PECOS Development Team
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

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include "largo_state.hpp"

using suzerain::perfect::largo_state;

BOOST_AUTO_TEST_CASE( size )
{
    BOOST_REQUIRE_EQUAL(sizeof(largo_state), 6*sizeof(double));
}

BOOST_AUTO_TEST_CASE( default_construct_produces_trivial_state )
{
    largo_state s;
    BOOST_CHECK_EQUAL(s.rho, 0);
    BOOST_CHECK_EQUAL(s.mx,  0);
    BOOST_CHECK_EQUAL(s.mx,  0);
    BOOST_CHECK_EQUAL(s.mx,  0);
    BOOST_CHECK_EQUAL(s.e,   0);
    BOOST_CHECK_EQUAL(s.p,   0);
    BOOST_CHECK(s.trivial());
}

// Used for as_is and rescale_one test cases, just below
static void some_library_call_for_as_is(double *d)
{
    BOOST_CHECK_EQUAL(d[0], 0); d[0] += 1;
    BOOST_CHECK_EQUAL(d[1], 1); d[1] += 1;
    BOOST_CHECK_EQUAL(d[2], 2); d[2] += 1;
    BOOST_CHECK_EQUAL(d[3], 3); d[3] += 1;
    BOOST_CHECK_EQUAL(d[4], 4); d[4] += 1;
    BOOST_CHECK_EQUAL(d[5], 5); d[5] += 1;
}

BOOST_AUTO_TEST_CASE( as_is )
{
    largo_state s;
    s.rho = 0;
    s.mx  = 1;
    s.my  = 2;
    s.mz  = 3;
    s.e   = 4;
    s.p   = 5;
    BOOST_CHECK(!s.trivial());

    some_library_call_for_as_is(s.as_is());

    BOOST_CHECK_EQUAL(s.rho, 1);
    BOOST_CHECK_EQUAL(s.mx,  2);
    BOOST_CHECK_EQUAL(s.my,  3);
    BOOST_CHECK_EQUAL(s.mz,  4);
    BOOST_CHECK_EQUAL(s.e,   5);
    BOOST_CHECK_EQUAL(s.p,   6);
}

BOOST_AUTO_TEST_CASE( rescale_one )
{
    largo_state s(4, 1, 2, 3, 0, 5);
    BOOST_REQUIRE_EQUAL(s.rho, 0);
    BOOST_REQUIRE_EQUAL(s.mx,  1);
    BOOST_REQUIRE_EQUAL(s.my,  2);
    BOOST_REQUIRE_EQUAL(s.mz,  3);
    BOOST_REQUIRE_EQUAL(s.e,   4);
    BOOST_REQUIRE_EQUAL(s.p,   5);

    some_library_call_for_as_is(s.rescale(1));

    BOOST_CHECK_EQUAL(s.rho, 1);
    BOOST_CHECK_EQUAL(s.mx,  2);
    BOOST_CHECK_EQUAL(s.my,  3);
    BOOST_CHECK_EQUAL(s.mz,  4);
    BOOST_CHECK_EQUAL(s.e,   5);
    BOOST_CHECK_EQUAL(s.p,   6);
}

// Used for rescale_two test case, just below
static void some_library_call_for_rescale_two(double *d)
{
    BOOST_CHECK_EQUAL(d[0],  0); d[0] += 1;
    BOOST_CHECK_EQUAL(d[1],  1); d[1] += 1;
    BOOST_CHECK_EQUAL(d[2],  2); d[2] += 1;
    BOOST_CHECK_EQUAL(d[3],  3); d[3] += 1;
    BOOST_CHECK_EQUAL(d[4],  8); d[4] += 1; // Note factor of two
    BOOST_CHECK_EQUAL(d[5], 10); d[5] += 1; // Note factor of two
}

BOOST_AUTO_TEST_CASE( rescale_two )
{
    largo_state s(4, 1, 2, 3, 0, 5);

    BOOST_REQUIRE_EQUAL(s.rho, 0);
    BOOST_REQUIRE_EQUAL(s.mx,  1);
    BOOST_REQUIRE_EQUAL(s.my,  2);
    BOOST_REQUIRE_EQUAL(s.mz,  3);
    BOOST_REQUIRE_EQUAL(s.e,   4);
    BOOST_REQUIRE_EQUAL(s.p,   5);
    BOOST_CHECK(!s.trivial());

    some_library_call_for_rescale_two(s.rescale(2.0));

    BOOST_CHECK_EQUAL(s.rho, 1  );
    BOOST_CHECK_EQUAL(s.mx,  2  );
    BOOST_CHECK_EQUAL(s.my,  3  );
    BOOST_CHECK_EQUAL(s.mz,  4  );
    BOOST_CHECK_EQUAL(s.e,   4.5);         // Note factor of two
    BOOST_CHECK_EQUAL(s.p,   5.5);         // Note factor of two
}
