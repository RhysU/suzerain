#include <suzerain/config.h>
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <mpi.h>
#include <fftw3-mpi.h>
#include <suzerain/error.h>
#include <suzerain/underling.h>
#include <boost/mpl/list_c.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/test_case_template.hpp>

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

typedef boost::mpl::list_c<int,1,2,3,5,7,11> howmany_value_list;

BOOST_AUTO_TEST_CASE_TEMPLATE( two_three_five, HOWMANY, howmany_value_list )
{
    using boost::make_counting_iterator;

    const int howmany = HOWMANY::type::value;
    const int NX = 2, NY = 3, NZ = 5;

    int procid;
    BOOST_REQUIRE_EQUAL(MPI_SUCCESS, MPI_Comm_rank(MPI_COMM_WORLD, &procid));

    // Create a grid that will be automagically cleaned up
    boost::shared_ptr<boost::remove_pointer<underling_grid>::type>
        grid(underling_grid_create(MPI_COMM_WORLD, NX, NY, NZ, 0, 0),
             &underling_grid_destroy);
    BOOST_REQUIRE(grid);

    // Create a problem that will be automagically cleaned up
    boost::shared_ptr<boost::remove_pointer<underling_problem>::type>
        problem(underling_problem_create(grid.get(), howmany),
                &underling_problem_destroy);
    BOOST_REQUIRE(problem);

    // Allocate space that will be automagically cleaned up
    const size_t local_memory = underling_local_memory(problem.get());
    boost::shared_ptr<underling_real>
        data((underling_real *) fftw_malloc(
                    local_memory*sizeof(underling_real)),
             &fftw_free);
    BOOST_REQUIRE(data);

    // Sanity check the buffer vs local extent information
    // Also gives us a reason to look up stride information
    underling_extents long_n[3];
    size_t extent_n[3];
    for (int i = 0; i < sizeof(long_n)/sizeof(long_n[0]); ++i) {
        long_n[i] = underling_local_extents(problem.get(), i);
        extent_n[i]
            = long_n[i].size[0] * long_n[i].size[1] * long_n[i].size[2];
    }
    BOOST_REQUIRE_LE(extent_n[0], local_memory);
    BOOST_REQUIRE_LE(extent_n[1], local_memory);
    BOOST_REQUIRE_LE(extent_n[2], local_memory);

    // Dump out the process' portion of the global grid
    if (howmany == 1) {
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

    // Create a plan that will be automagically cleaned up
    boost::shared_ptr<boost::remove_pointer<underling_plan>::type>
        plan(underling_plan_create(
                    problem.get(), data.get(),
                    UNDERLING_TRANSPOSE_ALL, FFTW_ESTIMATE),
             &underling_plan_destroy);
    BOOST_REQUIRE(plan);

    // Initialize entire local memory buffer while long in n2
    std::copy(make_counting_iterator(procid*10000.0),
              make_counting_iterator(procid*10000.0 + local_memory),
              data.get());

    // Transform from long in n2 to long in n0
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS,
                        underling_execute_long_n2_to_long_n1(plan.get()));
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS,
                        underling_execute_long_n1_to_long_n0(plan.get()));
    // Transform from long in n0 to long in n2
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS,
                        underling_execute_long_n0_to_long_n1(plan.get()));
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS,
                        underling_execute_long_n1_to_long_n2(plan.get()));

    // Ensure we successfully round-tripped back to long_n2 storage.
    // Buffer reqions beyond total_extent are excluded from the check.
    BOOST_REQUIRE_EQUAL_COLLECTIONS(
            data.get(),
            data.get() + extent_n[2],
            make_counting_iterator(procid*10000.0),
            make_counting_iterator(procid*10000.0 + extent_n[2]));
}
