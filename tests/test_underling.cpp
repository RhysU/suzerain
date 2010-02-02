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
#include <suzerain/underling_fft.hpp>
#include <fftw3-mpi.h>

// A test fixture to setup and teardown MPI and FFTW MPI
struct FFTWMPIFixture {

    FFTWMPIFixture() {
        MPI_Init(NULL, NULL); // NULL valid per MPI standard section 8.7
        fftw_mpi_init();
    }

    ~FFTWMPIFixture() {
        fftw_mpi_cleanup();
        MPI_Finalize();
    }
};

BOOST_GLOBAL_FIXTURE(FFTWMPIFixture);

// A test fixture to setup and teardown an underling use case
struct UnderlingFixture {

    UnderlingFixture(MPI_Comm comm,
                     const int n0, const int n1, const int n2,
                     const int howmany)
        : grid(comm, n0, n1, n2),
          problem(grid, howmany),
          data((underling_real *) fftw_malloc(
                    problem.local_memory()*sizeof(underling_real)),
                &fftw_free),
          plan(problem, data.get(),
               suzerain::underling::transpose::all, FFTW_ESTIMATE)
    {
        BOOST_REQUIRE(grid);
        BOOST_REQUIRE(problem);
        BOOST_REQUIRE(data);
        BOOST_REQUIRE(plan);
    }

    suzerain::underling::grid grid;
    suzerain::underling::problem problem;
    boost::shared_ptr<underling_real> data;
    suzerain::underling::plan plan;

};

BOOST_AUTO_TEST_SUITE(RoundTrip)

// Currently this test focuses on round-trip correctness
// TODO Test that data which should be long is in fact long
// TODO Test the provided stride information
// TODO Test reuse of grids for multiple problems
// TODO Test reuse of problems on multiple data sets
// TODO Test unidirectional (i.e. down-only) transforms

void test_round_trip(MPI_Comm comm,
                     const int n0, const int n1, const int n2,
                     const int howmany)
{
    namespace underling = suzerain::underling;

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
    UnderlingFixture f(comm, n0, n1, n2, howmany);
    const size_t local_memory = f.problem.local_memory();

    // Sanity check the buffer vs local extent information
    // Also gives us a reason to look up stride information
    underling_extents long_n[3];
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
            for (int i = 0; i < 3; ++i) {
                oss << "long_n" << i << ": " << long_n[i];
                if (i < 2) oss << ", ";
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
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  1);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  2);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  3);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  5);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8,  7);
    test_round_trip(MPI_COMM_WORLD, 8, 8, 8, 11);
}

BOOST_AUTO_TEST_CASE( roundtrip2x3x5 )
{
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  1);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  2);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  3);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  5);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5,  7);
    test_round_trip(MPI_COMM_WORLD, 2, 3, 5, 11);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( underling_fft )

void test_c2c(MPI_Comm comm,
              const int n0, const int n1, const int n2,
              const int howmany,
              const int long_i)
{
    UnderlingFixture f(comm, n0, n1, n2, howmany);
}

BOOST_AUTO_TEST_CASE( underling_fft_c2c )
{
    test_c2c(MPI_COMM_WORLD, 2, 3, 5, 1, 0);
}

BOOST_AUTO_TEST_SUITE_END()
