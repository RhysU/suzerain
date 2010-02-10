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
#include "test_tools.hpp"

// Useful namespace import
namespace underling = suzerain::underling;

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
                     const int howmany,
                     const unsigned transposed_flags = 0)
        : grid(comm, n0, n1, n2),
          problem(grid, howmany, transposed_flags),
          data((underling_real *) fftw_malloc(
                    problem.local_memory()*sizeof(underling_real)),
                &fftw_free),
          plan(problem, data.get(),
               underling::transpose::all, FFTW_ESTIMATE)
    {
        BOOST_REQUIRE(grid);
        BOOST_REQUIRE(problem);
        BOOST_REQUIRE(data);
        BOOST_REQUIRE(plan);
    }

    underling::grid grid;
    underling::problem problem;
    boost::shared_array<underling_real> data;
    underling::plan plan;

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

BOOST_AUTO_TEST_CASE( extents_consistency )
{
    UnderlingFixture f(MPI_COMM_WORLD, 2, 3, 5, 7);

    const underling_extents long_n[3] = {
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


BOOST_AUTO_TEST_SUITE( underling_fft )

void test_c2c(MPI_Comm comm,
              const int n0, const int n1, const int n2,
              const int howmany,
              const int long_i)
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

    const underling_real close
        = std::numeric_limits<underling_real>::epsilon()*10*n0*n1*n2;

    UnderlingFixture f(comm, n0, n1, n2, howmany);
    underling::fft::plan forward(underling::fft::plan::c2c_forward(),
                                 f.problem,
                                 long_i,
                                 f.data.get(),
                                 FFTW_ESTIMATE);
    BOOST_REQUIRE(forward);
    underling::fft::plan backward(underling::fft::plan::c2c_forward(),
                                  f.problem,
                                  long_i,
                                  f.data.get(),
                                  FFTW_ESTIMATE);
    BOOST_REQUIRE(backward);

    const underling_extents e = f.problem.local_extents(long_i);

    // Load up sample data that is constant on each pencil
    for (int i = 0; i < e.size[e.order[3]]; ++i) {
        for (int j = 0; j < e.size[e.order[2]]; ++j) {
            for (int k = 0; k < e.size[e.order[0]]; k += 2) {

                const underling_real v_re = (i+1)*(j+1)*(k+1);
                const underling_real v_im = -v_re;

                underling_real * const base = &f.data[
                      i*e.stride[e.order[3]]
                    + j*e.stride[e.order[2]]
                    + k
                ];

                for (int l = 0; l < e.size[e.order[1]]; ++l) {
                    base[ l*e.size[3]     ] = v_re;
                    base[ l*e.size[3] + 1 ] = v_im;
                }
            }
        }
    }

    // Transform from physical to wave space
    forward.execute();

    // Check the sample data transformed as expected
    for (int i = 0; i < e.size[e.order[3]]; ++i) {
        for (int j = 0; j < e.size[e.order[2]]; ++j) {
            for (int k = 0; k < e.size[e.order[0]]; k += 2) {

                // Constant magnitude mode is scaled by transform length
                const underling_real expected_re
                    = (i+1)*(j+1)*(k+1) * e.size[e.order[1]];
                const underling_real expected_im
                    = -expected_re;

                const underling_real * const base = &f.data[
                      i*e.stride[e.order[3]]
                    + j*e.stride[e.order[2]]
                    + k
                ];

                // Do we see expected zero (constant) mode magnitude?
                BOOST_CHECK_CLOSE(base[0], expected_re, close);
                BOOST_CHECK_CLOSE(base[1], expected_im, close);

                // Do we see no other modes?
                for (int l = 1; l < e.size[e.order[1]]; ++l) {
                    BOOST_CHECK_SMALL(
                            base[ l*e.stride[e.order[1]]     ], close);
                    BOOST_CHECK_SMALL(
                            base[ l*e.stride[e.order[1]] + 1 ], close);
                }
            }
        }
    }

    // Transform from wave space to physical space
    backward.execute();

    // Check that we recovered the original sample data
    for (int i = 0; i < e.size[e.order[3]]; ++i) {
        for (int j = 0; j < e.size[e.order[2]]; ++j) {
            for (int k = 0; k < e.size[e.order[0]]; k += 2) {

                // Result is scaled by transform length
                const underling_real expected_re
                    = (i+1)*(j+1)*(k+1) * e.size[e.order[1]];
                const underling_real expected_im
                    = -expected_re;

                const underling_real * const base = &f.data[
                      i*e.stride[e.order[3]]
                    + j*e.stride[e.order[2]]
                    + k
                ];

                for (int l = 0; l < e.size[e.order[1]]; ++l) {
                    BOOST_CHECK_CLOSE(
                            expected_re,
                            base[ l*e.stride[e.order[1]]     ],
                            close);
                    BOOST_CHECK_CLOSE(
                            expected_im,
                            base[ l*e.stride[e.order[1]] + 1 ],
                            close);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( underling_fft_c2c )
{
    // Cubic domain
    test_c2c(MPI_COMM_WORLD, 8, 8, 8, 2, 0); // 1 field
    test_c2c(MPI_COMM_WORLD, 8, 8, 8, 2, 1);
    test_c2c(MPI_COMM_WORLD, 8, 8, 8, 2, 2);
    test_c2c(MPI_COMM_WORLD, 8, 8, 8, 4, 0); // 2 fields
    test_c2c(MPI_COMM_WORLD, 8, 8, 8, 4, 1);
    test_c2c(MPI_COMM_WORLD, 8, 8, 8, 4, 2);
    test_c2c(MPI_COMM_WORLD, 8, 8, 8, 6, 0); // 3 fields
    test_c2c(MPI_COMM_WORLD, 8, 8, 8, 6, 1);
    test_c2c(MPI_COMM_WORLD, 8, 8, 8, 6, 2);

    // Non-cubic domain
    test_c2c(MPI_COMM_WORLD, 2, 3, 5, 2, 0); // 1 field
    test_c2c(MPI_COMM_WORLD, 2, 3, 5, 2, 1);
    test_c2c(MPI_COMM_WORLD, 2, 3, 5, 2, 2);
    test_c2c(MPI_COMM_WORLD, 2, 3, 5, 4, 0); // 2 fields
    test_c2c(MPI_COMM_WORLD, 2, 3, 5, 4, 1);
    test_c2c(MPI_COMM_WORLD, 2, 3, 5, 4, 2);
    test_c2c(MPI_COMM_WORLD, 2, 3, 5, 6, 0); // 3 fields
    test_c2c(MPI_COMM_WORLD, 2, 3, 5, 6, 1);
    test_c2c(MPI_COMM_WORLD, 2, 3, 5, 6, 2);
}

void test_c2r(MPI_Comm comm,
              const int n0, const int n1, const int n2,
              const int howmany,
              const int long_i)
{
    UnderlingFixture f(comm, n0, n1, n2, howmany);

    underling::fft::plan backward(underling::fft::plan::c2r_backward(),
                                  f.problem,
                                  long_i,
                                  f.data.get(),
                                  FFTW_ESTIMATE);
    BOOST_REQUIRE(backward);

    underling_extents e = f.problem.local_extents(long_i);
    const double close_enough
        =   std::numeric_limits<double>::epsilon()
          * 100*e.size[long_i]*e.size[long_i]*e.size[long_i];

    // Load up non-trivial sample data; note k indexing
    for (int i = 0; i < e.size[e.order[3]]; ++i) {
        for (int j = 0; j < e.size[e.order[2]]; ++j) {
            for (int k = 0; k < e.size[e.order[0]]; k += 2) {

                const periodic_function<double,int> pf(
                        2*(e.size[e.order[1]]-1), e.size[e.order[1]],
                        M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                underling_real * const base = &f.data[
                      i*e.stride[e.order[3]]
                    + j*e.stride[e.order[2]]
                    + k
                ];

                for (int l = 0; l < e.size[e.order[1]]; ++l) {
                    std::complex<double> val = pf.wave(l);
                    base[ l*e.stride[e.order[1]] ]
                        = val.real();
                    base[ l*e.stride[e.order[1]] + 1 ] = val.imag();
                }
            }
        }
    }

    backward.execute();

    // Check data transformed as expected; note k indexing
    for (int i = 0; i < e.size[e.order[3]]; ++i) {
        for (int j = 0; j < e.size[e.order[2]]; ++j) {
            for (int k = 0; k < e.size[e.order[0]] / 2; ++k) {

                const periodic_function<double,int> pf(
                        2*(e.size[e.order[1]]-1), e.size[e.order[1]],
                        M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*((k*2)+1));

                underling_real * const base = &f.data[
                      i*e.stride[e.order[3]]
                    + j*e.stride[e.order[2]]
                    + k
                ];

                for (int l = 0; l < 2*(e.size[e.order[1]]-1); ++l) {
                    const double expected = pf.physical(l);
                    const double actual   = base[ l*(e.stride[e.order[1]]/2) ];
                    BOOST_CHECK_CLOSE(expected, actual, close_enough);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( underling_fft_c2r_simple_n0 )
{
    // Long in n0, transform a single pencil
    test_c2r(MPI_COMM_SELF, 3, 1, 1, 2, 0);
    test_c2r(MPI_COMM_SELF, 3, 1, 1, 4, 0);
    test_c2r(MPI_COMM_SELF, 3, 1, 1, 6, 0);

    test_c2r(MPI_COMM_SELF, 4, 1, 1, 2, 0);
    test_c2r(MPI_COMM_SELF, 4, 1, 1, 4, 0);
    test_c2r(MPI_COMM_SELF, 4, 1, 1, 6, 0);

    test_c2r(MPI_COMM_SELF, 5, 1, 1, 2, 0);
    test_c2r(MPI_COMM_SELF, 5, 1, 1, 4, 0);
    test_c2r(MPI_COMM_SELF, 5, 1, 1, 6, 0);

    test_c2r(MPI_COMM_SELF, 6, 1, 1, 2, 0);
    test_c2r(MPI_COMM_SELF, 6, 1, 1, 4, 0);
    test_c2r(MPI_COMM_SELF, 6, 1, 1, 6, 0);

    // Long in n0, transform multiple pencils
    test_c2r(MPI_COMM_SELF, 3, 2, 2, 2, 0);
    test_c2r(MPI_COMM_SELF, 3, 2, 3, 4, 0);
    test_c2r(MPI_COMM_SELF, 3, 3, 2, 6, 0);

    test_c2r(MPI_COMM_SELF, 4, 2, 2, 2, 0);
    test_c2r(MPI_COMM_SELF, 4, 2, 3, 4, 0);
    test_c2r(MPI_COMM_SELF, 4, 3, 2, 6, 0);

    test_c2r(MPI_COMM_SELF, 5, 2, 2, 2, 0);
    test_c2r(MPI_COMM_SELF, 5, 2, 3, 4, 0);
    test_c2r(MPI_COMM_SELF, 5, 3, 2, 6, 0);

    test_c2r(MPI_COMM_SELF, 6, 2, 2, 2, 0);
    test_c2r(MPI_COMM_SELF, 6, 2, 3, 4, 0);
    test_c2r(MPI_COMM_SELF, 6, 3, 2, 6, 0);
}

BOOST_AUTO_TEST_CASE( underling_fft_c2r_simple_n1 )
{
    // Long in n1, transform a single pencil
    test_c2r(MPI_COMM_SELF, 1, 3, 1, 2, 1);
    test_c2r(MPI_COMM_SELF, 1, 3, 1, 4, 1);
    test_c2r(MPI_COMM_SELF, 1, 3, 1, 6, 1);

    test_c2r(MPI_COMM_SELF, 1, 4, 1, 2, 1);
    test_c2r(MPI_COMM_SELF, 1, 4, 1, 4, 1);
    test_c2r(MPI_COMM_SELF, 1, 4, 1, 6, 1);

    test_c2r(MPI_COMM_SELF, 1, 5, 1, 2, 1);
    test_c2r(MPI_COMM_SELF, 1, 5, 1, 4, 1);
    test_c2r(MPI_COMM_SELF, 1, 5, 1, 6, 1);

    test_c2r(MPI_COMM_SELF, 1, 6, 1, 2, 1);
    test_c2r(MPI_COMM_SELF, 1, 6, 1, 4, 1);
    test_c2r(MPI_COMM_SELF, 1, 6, 1, 6, 1);

    // Long in n1, transform multiple pencils
    test_c2r(MPI_COMM_SELF, 2, 3, 2, 2, 1);
    test_c2r(MPI_COMM_SELF, 2, 3, 3, 4, 1);
    test_c2r(MPI_COMM_SELF, 3, 3, 2, 6, 1);

    test_c2r(MPI_COMM_SELF, 2, 4, 2, 2, 1);
    test_c2r(MPI_COMM_SELF, 2, 4, 3, 4, 1);
    test_c2r(MPI_COMM_SELF, 3, 4, 2, 6, 1);

    test_c2r(MPI_COMM_SELF, 2, 5, 2, 2, 1);
    test_c2r(MPI_COMM_SELF, 2, 5, 3, 4, 1);
    test_c2r(MPI_COMM_SELF, 3, 5, 2, 6, 1);

    test_c2r(MPI_COMM_SELF, 2, 6, 2, 2, 1);
    test_c2r(MPI_COMM_SELF, 2, 6, 3, 4, 1);
    test_c2r(MPI_COMM_SELF, 3, 6, 2, 6, 1);
}

BOOST_AUTO_TEST_CASE( underling_fft_c2r_simple_n2 )
{
    // Long in n2, transform a single pencil
    test_c2r(MPI_COMM_SELF, 1, 1, 3, 2, 2);
    test_c2r(MPI_COMM_SELF, 1, 1, 3, 4, 2);
    test_c2r(MPI_COMM_SELF, 1, 1, 3, 6, 2);

    test_c2r(MPI_COMM_SELF, 1, 1, 4, 2, 2);
    test_c2r(MPI_COMM_SELF, 1, 1, 4, 4, 2);
    test_c2r(MPI_COMM_SELF, 1, 1, 4, 6, 2);

    test_c2r(MPI_COMM_SELF, 1, 1, 5, 2, 2);
    test_c2r(MPI_COMM_SELF, 1, 1, 5, 4, 2);
    test_c2r(MPI_COMM_SELF, 1, 1, 5, 6, 2);

    test_c2r(MPI_COMM_SELF, 1, 1, 6, 2, 2);
    test_c2r(MPI_COMM_SELF, 1, 1, 6, 4, 2);
    test_c2r(MPI_COMM_SELF, 1, 1, 6, 6, 2);

    // Long in n2, transform multiple pencils
    test_c2r(MPI_COMM_SELF, 2, 2, 3, 2, 2);
    test_c2r(MPI_COMM_SELF, 2, 3, 3, 4, 2);
    test_c2r(MPI_COMM_SELF, 3, 2, 3, 6, 2);

    test_c2r(MPI_COMM_SELF, 2, 2, 4, 2, 2);
    test_c2r(MPI_COMM_SELF, 2, 3, 4, 4, 2);
    test_c2r(MPI_COMM_SELF, 3, 2, 4, 6, 2);

    test_c2r(MPI_COMM_SELF, 2, 2, 5, 2, 2);
    test_c2r(MPI_COMM_SELF, 2, 3, 5, 4, 2);
    test_c2r(MPI_COMM_SELF, 3, 2, 5, 6, 2);

    test_c2r(MPI_COMM_SELF, 2, 2, 6, 2, 2);
    test_c2r(MPI_COMM_SELF, 2, 3, 6, 4, 2);
    test_c2r(MPI_COMM_SELF, 3, 2, 6, 6, 2);
}

void test_r2c(MPI_Comm comm,
              const int n0, const int n1, const int n2,
              const int howmany,
              const int long_i)
{
    UnderlingFixture f(comm, n0, n1, n2, howmany);

    underling::fft::plan forward(underling::fft::plan::r2c_forward(),
                                 f.problem,
                                 long_i,
                                 f.data.get(),
                                 FFTW_ESTIMATE);
    BOOST_REQUIRE(forward);

    underling_extents e = f.problem.local_extents(long_i);
    const double close_enough
        =   std::numeric_limits<double>::epsilon()
          * 100*e.size[long_i]*e.size[long_i]*e.size[long_i];

    // Load up non-trivial sample data; note k indexing
    for (int i = 0; i < e.size[e.order[3]]; ++i) {
        for (int j = 0; j < e.size[e.order[2]]; ++j) {
            for (int k = 0; k < e.size[e.order[0]] / 2; ++k) {

                const periodic_function<double,int> pf(
                        2*(e.size[e.order[1]]-1), e.size[e.order[1]],
                        M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                underling_real * const base = &f.data[
                      i*e.stride[e.order[3]]
                    + j*e.stride[e.order[2]]
                    + k
                ];

                for (int l = 0; l < 2*(e.size[e.order[1]]-1); ++l) {
                    base[ l*(e.stride[e.order[1]]/2) ] = pf.physical(l);
                }
            }
        }
    }

    forward.execute();

    // Check data transformed as expected; note k indexing
    for (int i = 0; i < e.size[e.order[3]]; ++i) {
        for (int j = 0; j < e.size[e.order[2]]; ++j) {
            for (int k = 0; k < e.size[e.order[0]]; k += 2) {

                const periodic_function<double,int> pf(
                        2*(e.size[e.order[1]]-1), e.size[e.order[1]],
                        M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*((k/2)+1));

                underling_real * const base = &f.data[
                      i*e.stride[e.order[3]]
                    + j*e.stride[e.order[2]]
                    + k
                ];

                for (int l = 0; l < (e.size[e.order[1]]/2+1); ++l) {
                    std::complex<double> expected
                        = pf.wave(l) * (2.0*(e.size[e.order[1]]-1));
                    BOOST_CHECK_CLOSE(
                            expected.real(),
                            base[ l*e.stride[e.order[1]]     ],
                            close_enough);
                    BOOST_CHECK_CLOSE(
                            expected.imag(),
                            base[ l*e.stride[e.order[1]] + 1 ],
                            close_enough);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( underling_fft_r2c_simple_n0 )
{
    // Long in n0, transform a single pencil
    test_r2c(MPI_COMM_SELF, 3, 1, 1, 2, 0);
    test_r2c(MPI_COMM_SELF, 3, 1, 1, 4, 0);
    test_r2c(MPI_COMM_SELF, 3, 1, 1, 6, 0);

    test_r2c(MPI_COMM_SELF, 4, 1, 1, 2, 0);
    test_r2c(MPI_COMM_SELF, 4, 1, 1, 4, 0);
    test_r2c(MPI_COMM_SELF, 4, 1, 1, 6, 0);

    test_r2c(MPI_COMM_SELF, 5, 1, 1, 2, 0);
    test_r2c(MPI_COMM_SELF, 5, 1, 1, 4, 0);
    test_r2c(MPI_COMM_SELF, 5, 1, 1, 6, 0);

    test_r2c(MPI_COMM_SELF, 6, 1, 1, 2, 0);
    test_r2c(MPI_COMM_SELF, 6, 1, 1, 4, 0);
    test_r2c(MPI_COMM_SELF, 6, 1, 1, 6, 0);

    // Long in n0, transform multiple pencils
    test_r2c(MPI_COMM_SELF, 3, 2, 2, 2, 0);
    test_r2c(MPI_COMM_SELF, 3, 2, 3, 4, 0);
    test_r2c(MPI_COMM_SELF, 3, 3, 2, 6, 0);

    test_r2c(MPI_COMM_SELF, 4, 2, 2, 2, 0);
    test_r2c(MPI_COMM_SELF, 4, 2, 3, 4, 0);
    test_r2c(MPI_COMM_SELF, 4, 3, 2, 6, 0);

    test_r2c(MPI_COMM_SELF, 5, 2, 2, 2, 0);
    test_r2c(MPI_COMM_SELF, 5, 2, 3, 4, 0);
    test_r2c(MPI_COMM_SELF, 5, 3, 2, 6, 0);

    test_r2c(MPI_COMM_SELF, 6, 2, 2, 2, 0);
    test_r2c(MPI_COMM_SELF, 6, 2, 3, 4, 0);
    test_r2c(MPI_COMM_SELF, 6, 3, 2, 6, 0);
}

BOOST_AUTO_TEST_CASE( underling_fft_r2c_simple_n1 )
{
    // Long in n1, transform a single pencil
    test_r2c(MPI_COMM_SELF, 1, 3, 1, 2, 1);
    test_r2c(MPI_COMM_SELF, 1, 3, 1, 4, 1);
    test_r2c(MPI_COMM_SELF, 1, 3, 1, 6, 1);

    test_r2c(MPI_COMM_SELF, 1, 4, 1, 2, 1);
    test_r2c(MPI_COMM_SELF, 1, 4, 1, 4, 1);
    test_r2c(MPI_COMM_SELF, 1, 4, 1, 6, 1);

    test_r2c(MPI_COMM_SELF, 1, 5, 1, 2, 1);
    test_r2c(MPI_COMM_SELF, 1, 5, 1, 4, 1);
    test_r2c(MPI_COMM_SELF, 1, 5, 1, 6, 1);

    test_r2c(MPI_COMM_SELF, 1, 6, 1, 2, 1);
    test_r2c(MPI_COMM_SELF, 1, 6, 1, 4, 1);
    test_r2c(MPI_COMM_SELF, 1, 6, 1, 6, 1);

    // Long in n1, transform multiple pencils
    test_r2c(MPI_COMM_SELF, 2, 3, 2, 2, 1);
    test_r2c(MPI_COMM_SELF, 2, 3, 3, 4, 1);
    test_r2c(MPI_COMM_SELF, 3, 3, 2, 6, 1);

    test_r2c(MPI_COMM_SELF, 2, 4, 2, 2, 1);
    test_r2c(MPI_COMM_SELF, 2, 4, 3, 4, 1);
    test_r2c(MPI_COMM_SELF, 3, 4, 2, 6, 1);

    test_r2c(MPI_COMM_SELF, 2, 5, 2, 2, 1);
    test_r2c(MPI_COMM_SELF, 2, 5, 3, 4, 1);
    test_r2c(MPI_COMM_SELF, 3, 5, 2, 6, 1);

    test_r2c(MPI_COMM_SELF, 2, 6, 2, 2, 1);
    test_r2c(MPI_COMM_SELF, 2, 6, 3, 4, 1);
    test_r2c(MPI_COMM_SELF, 3, 6, 2, 6, 1);
}

BOOST_AUTO_TEST_CASE( underling_fft_r2c_simple_n2 )
{
    // Long in n2, transform a single pencil
    test_r2c(MPI_COMM_SELF, 1, 1, 3, 2, 2);
    test_r2c(MPI_COMM_SELF, 1, 1, 3, 4, 2);
    test_r2c(MPI_COMM_SELF, 1, 1, 3, 6, 2);

    test_r2c(MPI_COMM_SELF, 1, 1, 4, 2, 2);
    test_r2c(MPI_COMM_SELF, 1, 1, 4, 4, 2);
    test_r2c(MPI_COMM_SELF, 1, 1, 4, 6, 2);

    test_r2c(MPI_COMM_SELF, 1, 1, 5, 2, 2);
    test_r2c(MPI_COMM_SELF, 1, 1, 5, 4, 2);
    test_r2c(MPI_COMM_SELF, 1, 1, 5, 6, 2);

    test_r2c(MPI_COMM_SELF, 1, 1, 6, 2, 2);
    test_r2c(MPI_COMM_SELF, 1, 1, 6, 4, 2);
    test_r2c(MPI_COMM_SELF, 1, 1, 6, 6, 2);

    // Long in n2, transform multiple pencils
    test_r2c(MPI_COMM_SELF, 2, 2, 3, 2, 2);
    test_r2c(MPI_COMM_SELF, 2, 3, 3, 4, 2);
    test_r2c(MPI_COMM_SELF, 3, 2, 3, 6, 2);

    test_r2c(MPI_COMM_SELF, 2, 2, 4, 2, 2);
    test_r2c(MPI_COMM_SELF, 2, 3, 4, 4, 2);
    test_r2c(MPI_COMM_SELF, 3, 2, 4, 6, 2);

    test_r2c(MPI_COMM_SELF, 2, 2, 5, 2, 2);
    test_r2c(MPI_COMM_SELF, 2, 3, 5, 4, 2);
    test_r2c(MPI_COMM_SELF, 3, 2, 5, 6, 2);

    test_r2c(MPI_COMM_SELF, 2, 2, 6, 2, 2);
    test_r2c(MPI_COMM_SELF, 2, 3, 6, 4, 2);
    test_r2c(MPI_COMM_SELF, 3, 2, 6, 6, 2);
}

BOOST_AUTO_TEST_SUITE_END()
