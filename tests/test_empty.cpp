//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
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

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/mpi.hpp>

#include "test_tools.hpp"

// Files performs no strictly useful test work.
// Rather it exercises a small portion of Boost.Test.
// It is intended to be used to generate valgrind suppressions.

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);       // Tear down BLAS

struct MPIFixture {
    MPIFixture()  { MPI_Init(NULL, NULL); }
    ~MPIFixture() { MPI_Finalize(); }
};
BOOST_GLOBAL_FIXTURE(MPIFixture);               // Setup and tear down MPI


BOOST_AUTO_TEST_SUITE(one_suite)

BOOST_AUTO_TEST_CASE( totally_suite )
{
    BOOST_CHECK(true);
    BOOST_TEST_PASSPOINT();
    BOOST_TEST_MESSAGE("Twiddling thumbs");
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(two_suite)

BOOST_AUTO_TEST_CASE( super_suite )
{
    BOOST_REQUIRE(true);
    BOOST_TEST_PASSPOINT();
    BOOST_TEST_MESSAGE("Still twiddling thumbs...");
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(tres_suite)

typedef boost::mpl::list< double,float > test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( super_suite, T, test_types )
{
    BOOST_REQUIRE_EQUAL(0, T());
    BOOST_TEST_PASSPOINT();
    BOOST_TEST_MESSAGE("Twiddling templated thumbs...");
}

BOOST_AUTO_TEST_SUITE_END()
