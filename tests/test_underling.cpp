#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>
#include <suzerain/error.h>
#include <suzerain/mpi.hpp>
#include <suzerain/underling.hpp>
#include <fftw3-mpi.h>
#include "test_tools.hpp"

#pragma warning(disable:383)

// Contains UnderlingFixture, FFTWMPIFixture
#include "test_underling_tools.hpp"
BOOST_GLOBAL_FIXTURE(FFTWMPIFixture);

// Useful namespace import
namespace underling = suzerain::underling;

BOOST_AUTO_TEST_SUITE(RoundTrip)

// Currently this test focuses on round-trip correctness
// TODO Test that data which should be long is in fact long
// TODO Test the provided stride information
// TODO Test reuse of grids for multiple problems
// TODO Test reuse of problems on multiple data sets
// TODO Test unidirectional (i.e. down-only) transforms

static void test_round_trip(MPI_Comm comm,
                            const int n0, const int n1, const int n2,
                            const int howmany,
                            const unsigned transposed_flags)
{
    int procid;
    BOOST_REQUIRE_EQUAL(MPI_SUCCESS, MPI_Comm_rank(comm, &procid));

    int nproc;
    BOOST_REQUIRE_EQUAL(MPI_SUCCESS, MPI_Comm_size(comm, &nproc));

    if (!procid) {
        BOOST_TEST_MESSAGE("Testing howmany = " << howmany
                           << " on " << n0 << "x" << n1 << "x" << n2
                           << " using " << nproc << " processor"
                           << (nproc > 1 ? "s" : "") );
    }

    // Establish the test fixture
    UnderlingFixture f(comm, n0, n1, n2, howmany, transposed_flags);
    const size_t local_memory = f.problem.local_memory();

    // Sanity check the buffer vs local extent information
    // Also gives us a reason to look up stride information
    underling::extents long_n[3];
    for (int i = 0; i < sizeof(long_n)/sizeof(long_n[0]); ++i) {
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
        }
    }

    // Initialize entire local memory buffer while long in n2
    std::copy(boost::make_counting_iterator(procid*10000.0),
              boost::make_counting_iterator(procid*10000.0 + local_memory),
              f.data.get());

    // Transform from long in n2 to long in n0
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, f.plan.execute_long_n2_to_long_n1());
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, f.plan.execute_long_n1_to_long_n0());
    // Transform from long in n0 to long in n2
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, f.plan.execute_long_n0_to_long_n1());
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, f.plan.execute_long_n1_to_long_n2());

    // Ensure we successfully round-tripped back to long_n2 storage.
    // Buffer regions beyond total_extent are excluded from the check.
    BOOST_REQUIRE_EQUAL_COLLECTIONS(
            f.data.get(),
            f.data.get() + long_n[2].extent,
            boost::make_counting_iterator(procid*10000.0),
            boost::make_counting_iterator(procid*10000.0 + long_n[2].extent));
}

BOOST_AUTO_TEST_CASE( roundtrip8x8x8 )
{
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  1, 0);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  2, 0);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  3, 0);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  5, 0);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  7, 0);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8, 11, 0);
}

BOOST_AUTO_TEST_CASE( roundtrip8x8x8_transposed_long_n2 )
{
    using suzerain::underling::transposed::long_n2;

    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  1, long_n2);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  2, long_n2);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  3, long_n2);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  5, long_n2);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  7, long_n2);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8, 11, long_n2);
}

BOOST_AUTO_TEST_CASE( roundtrip8x8x8_transposed_long_n0 )
{
    using suzerain::underling::transposed::long_n0;

    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  1, long_n0);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  2, long_n0);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  3, long_n0);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  5, long_n0);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  7, long_n0);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8, 11, long_n0);
}

BOOST_AUTO_TEST_CASE( roundtrip8x8x8_transposed_both )
{
    using suzerain::underling::transposed::long_n2;
    using suzerain::underling::transposed::long_n0;

    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  1, long_n2 | long_n0);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  2, long_n2 | long_n0);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  3, long_n2 | long_n0);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  5, long_n2 | long_n0);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  7, long_n2 | long_n0);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8, 11, long_n2 | long_n0);
}

BOOST_AUTO_TEST_CASE( roundtrip8x8x8_degenerate_howmany )
{
    using suzerain::underling::transposed::long_n2;
    using suzerain::underling::transposed::long_n0;

    test_round_trip(MPI_COMM_WORLD, 8, 8, 8, 0, 0);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8, 0, long_n2);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8, 0, long_n0);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8, 0, long_n2 | long_n0);
}

BOOST_AUTO_TEST_CASE( roundtrip2x3x5 )
{
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  1, 0);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  2, 0);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  3, 0);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  5, 0);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  7, 0);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5, 11, 0);
}

BOOST_AUTO_TEST_CASE( roundtrip2x3x5_transposed_long_n2 )
{
    using suzerain::underling::transposed::long_n2;

    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  1, long_n2);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  2, long_n2);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  3, long_n2);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  5, long_n2);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  7, long_n2);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5, 11, long_n2);
}

BOOST_AUTO_TEST_CASE( roundtrip2x3x5_transposed_long_n0 )
{
    using suzerain::underling::transposed::long_n0;

    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  1, long_n0);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  2, long_n0);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  3, long_n0);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  5, long_n0);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  7, long_n0);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5, 11, long_n0);
}

BOOST_AUTO_TEST_CASE( roundtrip2x3x5_transposed_long_both )
{
    using suzerain::underling::transposed::long_n2;
    using suzerain::underling::transposed::long_n0;

    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  1, long_n2 | long_n0);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  2, long_n2 | long_n0);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  3, long_n2 | long_n0);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  5, long_n2 | long_n0);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  7, long_n2 | long_n0);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5, 11, long_n2 | long_n0);
}

BOOST_AUTO_TEST_CASE( roundtrip2x3x5_degenerate_howmany )
{
    using suzerain::underling::transposed::long_n2;
    using suzerain::underling::transposed::long_n0;

    test_round_trip(MPI_COMM_WORLD, 2, 3, 5, 0, 0);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5, 0, long_n2);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5, 0, long_n0);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5, 0, long_n2 | long_n0);
}

BOOST_AUTO_TEST_CASE( extents_consistency )
{
    UnderlingFixture f(MPI_COMM_WORLD, 2, 3, 5, 7);

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

BOOST_AUTO_TEST_SUITE_END()
