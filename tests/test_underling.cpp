#include <suzerain/config.h>
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <mpi.h>
#include <fftw3-mpi.h>
#include <boost/array.hpp>
#include <boost/mpl/list_c.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <suzerain/error.h>
#include <suzerain/underling.hpp>

// Currently this test focuses on round-trip correctness
// TODO Test that data which should be long is in fact long
// TODO Test the provided stride information
// TODO Test reuse of grids for multiple problems
// TODO Test reuse of problems on multiple data sets
// TODO Test unidirectional (i.e. down-only) transforms

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

void round_trip_test(const int howmany,
                     const int NX, const int NY, const int NZ)
{
    namespace underling = suzerain::underling;

    int procid;
    BOOST_REQUIRE_EQUAL(MPI_SUCCESS, MPI_Comm_rank(MPI_COMM_WORLD, &procid));

    int nproc;
    BOOST_REQUIRE_EQUAL(MPI_SUCCESS, MPI_Comm_size(MPI_COMM_WORLD, &nproc));

    if (!procid) {
        BOOST_TEST_MESSAGE("Testing howmany = " << howmany
                           << " on " << NX << "x" << NY << "x" << NZ
                           << " using " << nproc << " processor"
                           << (nproc > 1 ? "s" : "") );
    }

    underling::grid grid(MPI_COMM_WORLD, NX, NY, NZ);
    BOOST_REQUIRE(grid);

    underling::problem problem(grid, howmany);
    BOOST_REQUIRE(problem);

    // Allocate space that will be automagically cleaned up
    const size_t local_memory = problem.local_memory();
    boost::shared_ptr<underling_real>
        data((underling_real *) fftw_malloc(
                    local_memory*sizeof(underling_real)),
             &fftw_free);
    BOOST_REQUIRE(data);

    // Sanity check the buffer vs local extent information
    // Also gives us a reason to look up stride information
    underling_extents long_n[3];
    for (int i = 0; i < sizeof(long_n)/sizeof(long_n[0]); ++i) {
        long_n[i] = problem.local_extents(i);
    }
    BOOST_REQUIRE_LE(long_n[0].extent, local_memory);
    BOOST_REQUIRE_LE(long_n[1].extent, local_memory);
    BOOST_REQUIRE_LE(long_n[2].extent, local_memory);

    // Sync up processors to output trace information
    for (int i = 0; i < nproc; ++i) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (procid == i) {
            // Dump out the process' portion of the global grid
            boost::format formatter("long_n%d: [%d,%d)x[%d,%d)x[%d,%d)");
            std::ostringstream oss;
            oss << "Rank " << procid << " has ";
            for (int i = 0; i < 3; ++i) {
                formatter % i;
                for (int j = 0; j < 3; ++j) {
                    formatter % long_n[i].start[j]
                              % (long_n[i].start[j] + long_n[i].size[j]);
                }
                oss << formatter.str();
                if (i < 2) oss << ", ";
            }
            BOOST_TEST_MESSAGE(oss.str());
        }
    }

    underling::plan plan(problem,
                         data.get(),
                         underling::transpose::all,
                         FFTW_ESTIMATE);
    BOOST_REQUIRE(plan);

    // Initialize entire local memory buffer while long in n2
    std::copy(boost::make_counting_iterator(procid*10000.0),
              boost::make_counting_iterator(procid*10000.0 + local_memory),
              data.get());

    // Transform from long in n2 to long in n0
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, plan.execute_long_n2_to_long_n1());
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, plan.execute_long_n1_to_long_n0());
    // Transform from long in n0 to long in n2
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, plan.execute_long_n0_to_long_n1());
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, plan.execute_long_n1_to_long_n2());

    // Ensure we successfully round-tripped back to long_n2 storage.
    // Buffer reqions beyond total_extent are excluded from the check.
    BOOST_REQUIRE_EQUAL_COLLECTIONS(
            data.get(),
            data.get() + long_n[2].extent,
            boost::make_counting_iterator(procid*10000.0),
            boost::make_counting_iterator(procid*10000.0 + long_n[2].extent));
}

typedef boost::mpl::list_c<int,1,2,3,5,7,11> howmany_values;

BOOST_AUTO_TEST_CASE_TEMPLATE(roundtrip8x8x8, HOWMANY, howmany_values)
{
    round_trip_test(HOWMANY::type::value, 8, 8, 8);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(roundtrip2x3x5, HOWMANY, howmany_values)
{
    round_trip_test(HOWMANY::type::value, 2, 3, 5);
}
