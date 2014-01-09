//-----------------------------------------------------------------------bl-
// underling 0.3.1: an FFTW MPI-based library for 3D pencil decompositions
// http://red.ices.utexas.edu/projects/underling
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
//
// This file is part of underling.
//
// underling is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// underling is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with underling.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------el-
// $Id$

#ifdef HAVE_CONFIG_H
#include <underling/config.h>
#endif
#include <underling/underling.hpp>

#include <cstdio>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <underling/error.h>
#include <mpi.h>
#include <fftw3-mpi.h>
#include "test_tools.hpp"
#include "fixtures.hpp"

// Currently this test focuses on round-trip correctness
// TODO Test that data which should be long is in fact long
// TODO Test the provided stride information
// TODO Test reuse of grids for multiple problems
// TODO Test reuse of problems on multiple data sets
// TODO Test unidirectional (i.e. down-only) transforms

// For unary function-based test case registration
struct tc
{
    int n0, n1, n2, howmany;
    bool in_place;
    unsigned flags;
};

static void test_round_trip(tc t)
{
    BoostFailErrorHandlerFixture fix1;

    // Unpack test case parameters
    MPI_Comm comm        = MPI_COMM_WORLD;
    const int n0         = t.n0;
    const int n1         = t.n1;
    const int n2         = t.n2;
    const int howmany    = t.howmany;
    const unsigned flags = t.flags;
    const bool in_place  = t.in_place;

    int procid;
    BOOST_REQUIRE_EQUAL(MPI_SUCCESS, MPI_Comm_rank(comm, &procid));

    int nproc;
    BOOST_REQUIRE_EQUAL(MPI_SUCCESS, MPI_Comm_size(comm, &nproc));

    if (!procid) {
        BOOST_TEST_MESSAGE("Testing "
                           << (in_place ? "in-place" : "out-of-place")
                           << " howmany = " << howmany
                           << " on " << n0 << "x" << n1 << "x" << n2
                           << " using " << nproc << " processor"
                           << (nproc > 1 ? "s" : "") );
    }

    // Establish the test fixture
    UnderlingFixture f(comm, n0, n1, n2, howmany, flags, in_place);
    const size_t local_memory = f.problem.local_memory();

    // Sanity check the buffer vs local extent information
    // Also gives us a reason to look up stride information
    underling::extents long_n[3];
    for (std::size_t i = 0; i < sizeof(long_n)/sizeof(long_n[0]); ++i) {
        long_n[i] = f.problem.local_extents(i);
    }
    BOOST_REQUIRE_LE(long_n[0].extent, local_memory);
    BOOST_REQUIRE_LE(long_n[1].extent, local_memory);
    BOOST_REQUIRE_LE(long_n[2].extent, local_memory);

    // Sync up processors to output trace information
    for (int i = 0; i < nproc; ++i) {
        MPI_Barrier(comm);
        if (procid == i) {
            // Dump out the process' portion of the global grid
            std::ostringstream oss;
            oss << "Rank " << procid << " has ";
            for (int j = 0; j < 3; ++j) {
                oss << "long_n" << j << ": " << long_n[j];
                if (j < 2) oss << ", ";
            }
            BOOST_TEST_MESSAGE(oss.str());
            std::cout.flush();
            std::cerr.flush();
        }
    }

    // Initialize buffers while long in n2
    std::copy(boost::make_counting_iterator(procid*10000.0),
              boost::make_counting_iterator(procid*10000.0 + local_memory),
              f.in);
    if (!f.in_place) f.fill_out_with_NaNs();

    // Transform from long in n2 to long in n1
    BOOST_REQUIRE_EQUAL(UNDERLING_SUCCESS,
            f.plan.execute_long_n2_to_long_n1(f.in, f.out));
    if (!f.in_place) f.fill_in_with_NaNs(); // Poison old buffer

    // Transform from long in n1 to long in n0
    BOOST_REQUIRE_EQUAL(UNDERLING_SUCCESS,
            f.plan.execute_long_n1_to_long_n0(f.out, f.in));
    if (!f.in_place) f.fill_out_with_NaNs(); // Poison old buffer

    // Transform from long in n0 to long in n1
    BOOST_REQUIRE_EQUAL(UNDERLING_SUCCESS,
            f.plan.execute_long_n0_to_long_n1(f.in, f.out));
    if (!f.in_place) f.fill_in_with_NaNs();  // Poison old buffer

    // Transform from long in n1 to long in n2
    BOOST_REQUIRE_EQUAL(UNDERLING_SUCCESS,
            f.plan.execute_long_n1_to_long_n2(f.out, f.in));
    if (!f.in_place) f.fill_out_with_NaNs();  // Poison old buffer

    // Ensure we successfully round-tripped back to long_n2 storage.
    // Buffer regions beyond total_extent are excluded from the check.
    BOOST_REQUIRE_EQUAL_COLLECTIONS(
            f.in, f.in + long_n[2].extent,
            boost::make_counting_iterator(procid*10000.0),
            boost::make_counting_iterator(procid*10000.0 + long_n[2].extent));
}

boost::unit_test::test_suite*
init_unit_test_suite( int argc, char* argv[] )
{
    MPI_Init(&argc, &argv);            // Initialize MPI
    atexit((void (*)()) MPI_Finalize); // Register finalize MPI
    underling_init(&argc, &argv, 0);   // Initialize underling prereqs
    atexit(&underling_cleanup);        // Register finalize underling prereqs

    boost::unit_test::framework::master_test_suite().p_name.value = __FILE__;

    // Size of global extents
    const int extents[][3] = { { 2, 3, 5 },
                               { 7, 5, 3 },
                               { 8, 8, 8 } };

    // Number of scalars to transpose
    const int howmanys[] = { 1, 2, 3, 4, 5, 7, 11, /* degenerate */0 };

    // In-place vs out-of-place
    const bool places[] = { true, false };

    // Transposed in/out-like flags
    using underling::transposed::long_n2;
    using underling::transposed::long_n0;
    const unsigned flags[] = { 0, long_n2, long_n0, long_n2 | long_n0 };

    // Create an outer product of all the cases we want to run
    const size_t ncases = sizeof(extents)/sizeof(extents[0])
                        * sizeof(howmanys)/sizeof(howmanys[0])
                        * sizeof(places)/sizeof(places[0])
                        * sizeof(flags)/sizeof(flags[0]);

    tc * const cases = new tc[ncases];
    tc *       c     = cases;

    for (size_t e = 0; e < sizeof(extents)/sizeof(extents[0]); ++e)
    for (size_t h = 0; h < sizeof(howmanys)/sizeof(howmanys[0]); ++h)
    for (size_t p = 0; p < sizeof(places)/sizeof(places[0]); ++p)
    for (size_t f = 0; f < sizeof(flags)/sizeof(flags[0]); ++f)
    {
        c->n0       = extents[e][0];
        c->n1       = extents[e][1];
        c->n2       = extents[e][2];
        c->howmany  = howmanys[h];
        c->in_place = places[p];
        c->flags    = flags[f];
        ++c;
    }

    // Register the outer product
    boost::unit_test::framework::master_test_suite().add(
            BOOST_PARAM_TEST_CASE( &test_round_trip, cases, cases + ncases),
            /* timeout in seconds */ 30 );

    delete[] cases;

    return 0;
}

static void test_extents_consistency(const bool in_place = true)
{
    UnderlingFixture f(MPI_COMM_WORLD, 2, 3, 5, 7, 0 /*flags*/, in_place);

    const underling::extents long_n[3] = {
        f.problem.local_extents(0),
        f.problem.local_extents(1),
        f.problem.local_extents(2)
    };

    for (int i = 0; i < 3; ++i) {
        BOOST_CHECK_EQUAL(long_n[i].extent, f.problem.local(i));

        int start[4];

        BOOST_CHECK_EQUAL(long_n[i].extent,
                          f.problem.local(i, start));
        BOOST_CHECK_EQUAL_COLLECTIONS(long_n[i].start, long_n[i].start + 4,
                                      start, start + 4);

        int size[4];

        BOOST_CHECK_EQUAL(long_n[i].extent,
                          f.problem.local(i, start, size));
        BOOST_CHECK_EQUAL_COLLECTIONS(long_n[i].size, long_n[i].size + 4,
                                      size, size + 4);

        int stride[4];

        BOOST_CHECK_EQUAL(long_n[i].extent,
                          f.problem.local(i, start, size, stride));
        BOOST_CHECK_EQUAL_COLLECTIONS(long_n[i].stride, long_n[i].stride + 4,
                                      stride, stride + 4);

        int order[4];

        BOOST_CHECK_EQUAL(long_n[i].extent,
                          f.problem.local(i, start, size, stride, order));
        BOOST_CHECK_EQUAL_COLLECTIONS(long_n[i].stride, long_n[i].stride + 4,
                                      stride, stride + 4);
    }
}

BOOST_AUTO_TEST_CASE( extents_consistency )
{
    test_extents_consistency(true);  // In-place
    test_extents_consistency(false); // Out-of-place
}
