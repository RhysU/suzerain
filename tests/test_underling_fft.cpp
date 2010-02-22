#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>
#include <suzerain/error.h>
#include <suzerain/mpi.hpp>
#include <suzerain/underling_fft.hpp>
#include <fftw3-mpi.h>
#include "test_tools.hpp"

// Contains UnderlingFixture, FFTWMPIFixture
#include "test_underling_tools.hpp"
BOOST_GLOBAL_FIXTURE(FFTWMPIFixture);

// Useful namespace import
namespace underling = suzerain::underling;

// Pull out the two slow directions, excluding long_ni
void slow_non_long_directions(
        const int * const order, const int long_ni, int &dir_i, int &dir_j)
{
    BOOST_REQUIRE(order);

    if (order[2] == long_ni) {
        dir_i = order[4];
        dir_j = order[3];
    } else if (order[3] == long_ni) {
        dir_i = order[4];
        dir_j = order[2];
    } else if (order[4] == long_ni) {
        dir_i = order[3];
        dir_j = order[2];
    } else {
        BOOST_FAIL("long_ni not in {2, 3, 4}");
    }
}

// Ensure that FFTW can handle compressing directions, which will require at
// patch above FFTW 3.3alpha1 submitted to Stephen and Matteo directly.
// Suspect the FFTW installation has not been patched if this fails.
void ensureFFTWTensor7PatchInPlace() {
    double buffer[30];
    const fftw_iodim howmany_dims[] = {
        { 5,   6,   1 },
        { 2,   3,  15 }, // This direction...
        { 3,   1,   5 }, // ...and this one should be merged.
    };
    int howmany_rank = sizeof(howmany_dims)/sizeof(howmany_dims[0]);
    // If the two directions are not merged, the in-place transpose planner
    // cannot handle the case and we will get a NULL plan.
    const fftw_plan plan = fftw_plan_guru_r2r(
            0, NULL, howmany_rank, howmany_dims, buffer, buffer,
            NULL, FFTW_ESTIMATE);
    BOOST_REQUIRE_MESSAGE(
            plan, "Critical FFTW3 patch may not have been applied.");
}


BOOST_FIXTURE_TEST_SUITE( underling_fft_general,
                          BoostFailErrorHandlerFixture )

BOOST_AUTO_TEST_CASE( extents_consistency )
{
    UnderlingFixture f(MPI_COMM_WORLD, 2, 3, 5, 6);

    underling::fft::plan backward(underling::fft::plan::c2c_forward(),
                                  f.problem,
                                  0,
                                  f.data.get(),
                                  FFTW_ESTIMATE);

    const int N = 5;

    // Check input information
    {
        const underling::fft::extents input = backward.local_extents_input();

        int start[N];
        backward.local_input(start);
        BOOST_CHECK_EQUAL_COLLECTIONS(input.start, input.start + N,
                                      start, start + N);

        int size[N];
        backward.local_input(start, size);
        BOOST_CHECK_EQUAL_COLLECTIONS(input.size, input.size + N,
                                      size, size + N);

        int stride[5];
        backward.local_input(start, size, stride);
        BOOST_CHECK_EQUAL_COLLECTIONS(input.stride, input.stride + N,
                                      stride, stride + N);

        int order[5];
        backward.local_input(start, size, stride, order);
        BOOST_CHECK_EQUAL_COLLECTIONS(input.stride, input.stride + N,
                                      stride, stride + N);
    }

    // Check output information
    {
        const underling::fft::extents output = backward.local_extents_output();

        int start[N];
        backward.local_output(start);
        BOOST_CHECK_EQUAL_COLLECTIONS(output.start, output.start + N,
                                      start, start + N);

        int size[N];
        backward.local_output(start, size);
        BOOST_CHECK_EQUAL_COLLECTIONS(output.size, output.size + N,
                                      size, size + N);

        int stride[5];
        backward.local_output(start, size, stride);
        BOOST_CHECK_EQUAL_COLLECTIONS(output.stride, output.stride + N,
                                      stride, stride + N);

        int order[5];
        backward.local_output(start, size, stride, order);
        BOOST_CHECK_EQUAL_COLLECTIONS(output.stride, output.stride + N,
                                      stride, stride + N);
    }
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_FIXTURE_TEST_SUITE( underling_fft_c2c_forward,
                          BoostFailErrorHandlerFixture )

// Forward physical-to-wave followed by wave-to-physical
void test_c2c_forward(MPI_Comm comm,
                      const int n0, const int n1, const int n2,
                      const int howmany,
                      const int long_ni,
                      const unsigned transposed_flags = 0)

{
    int procid;
    BOOST_REQUIRE_EQUAL(MPI_SUCCESS, MPI_Comm_rank(comm, &procid));

    int nproc;
    BOOST_REQUIRE_EQUAL(MPI_SUCCESS, MPI_Comm_size(comm, &nproc));

    if (!procid) {
        BOOST_TEST_MESSAGE("Testing howmany = " << howmany
                           << " on " << n0 << "x" << n1 << "x" << n2
                           << " using " << nproc << " processor"
                           << (nproc > 1 ? "s" : "")
                           << " when long in " << long_ni
                           << " with flags " << transposed_flags);
    }

    const underling_real close
        = std::numeric_limits<underling_real>::epsilon()*100*n0*n1*n2;

    UnderlingFixture f(comm, n0, n1, n2, howmany, transposed_flags);
    const underling::extents extents = f.problem.local_extents(long_ni);
    underling::fft::plan forward(underling::fft::plan::c2c_forward(),
                                 f.problem,
                                 long_ni,
                                 f.data.get(),
                                 FFTW_ESTIMATE);
    BOOST_REQUIRE(forward);
    underling::fft::plan backward(forward,         // Inverse constructor!
                                  f.data.get(),
                                  FFTW_ESTIMATE);
    BOOST_REQUIRE(backward);

    // Stride information consistency check
    BOOST_REQUIRE_EQUAL(forward.local_extents_output(),
                        backward.local_extents_input());
    BOOST_REQUIRE_EQUAL(backward.local_extents_output(),
                        forward.local_extents_input());

    // Load up sample data
    {
        const underling::fft::extents e = forward.local_extents_input();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(e.size[long_ni],
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    underling_real * const pencil = &f.data[
                          i*e.stride[dir_i]
                        + j*e.stride[dir_j]
                        + k*e.stride[dir_k]
                    ];

                    for (int l = 0; l < e.size[long_ni]; ++l) {
                        underling_real * const base
                            = pencil + l*e.stride[long_ni];
                        const underling_real v_re =  pf.physical(l);
                        const underling_real v_im = -v_re;
                        base[0]                    = v_re;
                        base[e.stride[e.order[0]]] = v_im;
                    }
                }
            }
        }
    }

    // Transform from physical to wave space
    forward.execute();

    // Check the sample data transformed as expected
    {
        const underling::fft::extents e = forward.local_extents_output();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(e.size[long_ni],
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    const underling_real * const pencil = &f.data[
                          i*e.stride[dir_i]
                        + j*e.stride[dir_j]
                        + k*e.stride[dir_k]
                    ];

                    for (int l = 0; l < e.size[long_ni]; ++l) {
                        const underling_real * const base
                            = pencil + l*e.stride[long_ni];
                        const std::complex<double> mode
                            = pf.wave(l) * (double) pf.N;
                        const double expected_re = mode.real() + mode.imag();
                        const double expected_im = mode.imag() - mode.real();
                        BOOST_CHECK_CLOSE(
                                expected_re,
                                base[0],
                                close);
                        BOOST_CHECK_CLOSE(
                                expected_im,
                                base[e.stride[e.order[0]]],
                                close);
                    }
                }
            }
        }
    }

    // Transform from wave space to physical space
    backward.execute();

    // Check that we recovered the scaled sample data
    {
        const underling::fft::extents e = backward.local_extents_output();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(e.size[long_ni],
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    const underling_real * const pencil = &f.data[
                          i*e.stride[dir_i]
                        + j*e.stride[dir_j]
                        + k*e.stride[dir_k]
                    ];

                    for (int l = 0; l < e.size[long_ni]; ++l) {
                        const underling_real * const base
                            = pencil + l*e.stride[long_ni];
                        const double expected_re =   pf.physical(l) * pf.N;
                        const double expected_im = - expected_re;
                        BOOST_CHECK_CLOSE(
                                expected_re,
                                base[0],
                                close);
                        BOOST_CHECK_CLOSE(
                                expected_im,
                                base[e.stride[e.order[0]]],
                                close);
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( underling_fft_c2c_forward )
{
    // Non-cubic domain
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 2, 0); // 1 field
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 2, 1);
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 2, 2);
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 4, 0); // 2 fields
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 4, 1);
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 4, 2);
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 6, 0); // 3 fields
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 6, 1);
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 6, 2);

    // Cubic domain
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 2, 0); // 1 field
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 2, 1);
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 2, 2);
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 4, 0); // 2 fields
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 4, 1);
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 4, 2);
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 6, 0); // 3 fields
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 6, 1);
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 6, 2);
}

BOOST_AUTO_TEST_CASE( underling_fft_c2c_forward_transposed_long_n2 )
{
    ensureFFTWTensor7PatchInPlace();

    using suzerain::underling::transposed::long_n2;

    // Non-cubic domain
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 2, 0, long_n2); // 1 field
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 2, 1, long_n2);
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 2, 2, long_n2);
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 4, 0, long_n2); // 2 fields
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 4, 1, long_n2);
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 4, 2, long_n2);
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 6, 0, long_n2); // 3 fields
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 6, 1, long_n2);
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 6, 2, long_n2);

    // Cubic domain
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 2, 0, long_n2); // 1 field
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 2, 1, long_n2);
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 2, 2, long_n2);
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 4, 0, long_n2); // 2 fields
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 4, 1, long_n2);
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 4, 2, long_n2);
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 6, 0, long_n2); // 3 fields
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 6, 1, long_n2);
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 6, 2, long_n2);
}

BOOST_AUTO_TEST_CASE( underling_fft_c2c_forward_transposed_long_n0 )
{
    ensureFFTWTensor7PatchInPlace();

    using suzerain::underling::transposed::long_n0;

    // Non-cubic domain
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 2, 0, long_n0); // 1 field
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 2, 1, long_n0);
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 2, 2, long_n0);
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 4, 0, long_n0); // 2 fields
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 4, 1, long_n0);
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 4, 2, long_n0);
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 6, 0, long_n0); // 3 fields
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 6, 1, long_n0);
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 6, 2, long_n0);

    // Cubic domain
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 2, 0, long_n0); // 1 field
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 2, 1, long_n0);
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 2, 2, long_n0);
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 4, 0, long_n0); // 2 fields
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 4, 1, long_n0);
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 4, 2, long_n0);
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 6, 0, long_n0); // 3 fields
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 6, 1, long_n0);
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 6, 2, long_n0);
}

BOOST_AUTO_TEST_CASE( underling_fft_c2c_forward_transposed_both )
{
    ensureFFTWTensor7PatchInPlace();

    using suzerain::underling::transposed::long_n2;
    using suzerain::underling::transposed::long_n0;

    // Non-cubic domain
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 2, 0, long_n2 | long_n0); // 1
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 2, 1, long_n2 | long_n0);
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 2, 2, long_n2 | long_n0);
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 4, 0, long_n2 | long_n0); // 2
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 4, 1, long_n2 | long_n0);
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 4, 2, long_n2 | long_n0);
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 6, 0, long_n2 | long_n0); // 3
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 6, 1, long_n2 | long_n0);
    test_c2c_forward(MPI_COMM_WORLD, 3, 5, 7, 6, 2, long_n2 | long_n0);

    // Cubic domain
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 2, 0, long_n2 | long_n0); // 1
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 2, 1, long_n2 | long_n0);
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 2, 2, long_n2 | long_n0);
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 4, 0, long_n2 | long_n0); // 2
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 4, 1, long_n2 | long_n0);
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 4, 2, long_n2 | long_n0);
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 6, 0, long_n2 | long_n0); // 3
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 6, 1, long_n2 | long_n0);
    test_c2c_forward(MPI_COMM_WORLD, 8, 8, 8, 6, 2, long_n2 | long_n0);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_FIXTURE_TEST_SUITE( underling_fft_c2c_backward,
                          BoostFailErrorHandlerFixture )

// Backward wave-to-physical followed by physical-to-wave
void test_c2c_backward(MPI_Comm comm,
                       const int n0, const int n1, const int n2,
                       const int howmany,
                       const int long_ni,
                       const unsigned transposed_flags = 0)
{
    int procid;
    BOOST_REQUIRE_EQUAL(MPI_SUCCESS, MPI_Comm_rank(comm, &procid));

    int nproc;
    BOOST_REQUIRE_EQUAL(MPI_SUCCESS, MPI_Comm_size(comm, &nproc));

    if (!procid) {
        BOOST_TEST_MESSAGE("Testing howmany = " << howmany
                           << " on " << n0 << "x" << n1 << "x" << n2
                           << " using " << nproc << " processor"
                           << (nproc > 1 ? "s" : "")
                           << " when long in " << long_ni
                           << " with flags " << transposed_flags);
    }

    const underling_real close
        = std::numeric_limits<underling_real>::epsilon()*100*n0*n1*n2;

    UnderlingFixture f(comm, n0, n1, n2, howmany, transposed_flags);
    underling::fft::plan backward(underling::fft::plan::c2c_backward(),
                                  f.problem,
                                  long_ni,
                                  f.data.get(),
                                  FFTW_ESTIMATE);
    BOOST_REQUIRE(backward);
    const underling::extents extents = f.problem.local_extents(long_ni);
    underling::fft::plan forward(backward,         // Inverse constructor!
                                 f.data.get(),
                                 FFTW_ESTIMATE);
    BOOST_REQUIRE(forward);

    // Stride information consistency check
    BOOST_REQUIRE_EQUAL(backward.local_extents_output(),
                        forward.local_extents_input());
    BOOST_REQUIRE_EQUAL(forward.local_extents_output(),
                        backward.local_extents_input());

    // Load up sample data
    {
        const underling::fft::extents e = backward.local_extents_input();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(e.size[long_ni],
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    underling_real * const pencil = &f.data[
                          i*e.stride[dir_i]
                        + j*e.stride[dir_j]
                        + k*e.stride[dir_k]
                    ];

                    for (int l = 0; l < e.size[long_ni]; ++l) {
                        underling_real * const base
                            = pencil + l*e.stride[long_ni];
                        const std::complex<double> mode = pf.wave(l);
                        base[0]                    = mode.real() + mode.imag();
                        base[e.stride[e.order[0]]] = mode.imag() - mode.real();
                    }
                }
            }
        }
    }

    // Transform from wave space to physical space
    backward.execute();

    // Check the sample data transformed as expected
    {
        const underling::fft::extents e = backward.local_extents_output();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(e.size[long_ni],
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    const underling_real * const pencil = &f.data[
                          i*e.stride[dir_i]
                        + j*e.stride[dir_j]
                        + k*e.stride[dir_k]
                    ];

                    for (int l = 0; l < e.size[long_ni]; ++l) {
                        const underling_real * const base
                            = pencil + l*e.stride[long_ni];
                        const underling_real expected_re =  pf.physical(l);
                        const underling_real expected_im = -expected_re;
                        BOOST_CHECK_CLOSE(
                                expected_re,
                                base[0],
                                close);
                        BOOST_CHECK_CLOSE(
                                expected_im,
                                base[e.stride[e.order[0]]],
                                close);
                    }
                }
            }
        }
    }

    // Transform from physical to wave space
    forward.execute();

    // Check the sample data transformed as expected
    {
        const underling::fft::extents e = forward.local_extents_output();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(e.size[long_ni],
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    const underling_real * const pencil = &f.data[
                          i*e.stride[dir_i]
                        + j*e.stride[dir_j]
                        + k*e.stride[dir_k]
                    ];

                    for (int l = 0; l < e.size[long_ni]; ++l) {
                        const underling_real * const base
                            = pencil + l*e.stride[long_ni];
                        const std::complex<double> mode
                            = pf.wave(l) * (double) pf.N;
                        const double expected_re = mode.real() + mode.imag();
                        const double expected_im = mode.imag() - mode.real();
                        BOOST_CHECK_CLOSE(
                                expected_re,
                                base[0],
                                close);
                        BOOST_CHECK_CLOSE(
                                expected_im,
                                base[e.stride[e.order[0]]],
                                close);
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( underling_fft_c2c_backward )
{
    // Non-cubic domain
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 2, 0); // 1 field
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 2, 1);
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 2, 2);
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 4, 0); // 2 fields
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 4, 1);
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 4, 2);
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 6, 0); // 3 fields
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 6, 1);
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 6, 2);

    // Cubic domain
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 2, 0); // 1 field
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 2, 1);
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 2, 2);
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 4, 0); // 2 fields
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 4, 1);
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 4, 2);
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 6, 0); // 3 fields
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 6, 1);
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 6, 2);
}

BOOST_AUTO_TEST_CASE( underling_fft_c2c_backward_transposed_long_n2 )
{
    ensureFFTWTensor7PatchInPlace();

    using suzerain::underling::transposed::long_n2;

    // Non-cubic domain
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 2, 0, long_n2); // 1 field
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 2, 1, long_n2);
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 2, 2, long_n2);
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 4, 0, long_n2); // 2 fields
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 4, 1, long_n2);
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 4, 2, long_n2);
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 6, 0, long_n2); // 3 fields
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 6, 1, long_n2);
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 6, 2, long_n2);

    // Cubic domain
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 2, 0, long_n2); // 1 field
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 2, 1, long_n2);
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 2, 2, long_n2);
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 4, 0, long_n2); // 2 fields
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 4, 1, long_n2);
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 4, 2, long_n2);
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 6, 0, long_n2); // 3 fields
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 6, 1, long_n2);
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 6, 2, long_n2);
}

BOOST_AUTO_TEST_CASE( underling_fft_c2c_backward_transposed_long_n0 )
{
    ensureFFTWTensor7PatchInPlace();

    using suzerain::underling::transposed::long_n0;

    // Non-cubic domain
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 2, 0, long_n0); // 1 field
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 2, 1, long_n0);
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 2, 2, long_n0);
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 4, 0, long_n0); // 2 fields
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 4, 1, long_n0);
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 4, 2, long_n0);
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 6, 0, long_n0); // 3 fields
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 6, 1, long_n0);
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 6, 2, long_n0);

    // Cubic domain
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 2, 0, long_n0); // 1 field
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 2, 1, long_n0);
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 2, 2, long_n0);
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 4, 0, long_n0); // 2 fields
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 4, 1, long_n0);
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 4, 2, long_n0);
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 6, 0, long_n0); // 3 fields
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 6, 1, long_n0);
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 6, 2, long_n0);
}

BOOST_AUTO_TEST_CASE( underling_fft_c2c_backward_transposed_both )
{
    ensureFFTWTensor7PatchInPlace();

    using suzerain::underling::transposed::long_n2;
    using suzerain::underling::transposed::long_n0;

    // Non-cubic domain
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 2, 0, long_n2 | long_n0); // 1
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 2, 1, long_n2 | long_n0);
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 2, 2, long_n2 | long_n0);
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 4, 0, long_n2 | long_n0); // 2
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 4, 1, long_n2 | long_n0);
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 4, 2, long_n2 | long_n0);
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 6, 0, long_n2 | long_n0); // 3
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 6, 1, long_n2 | long_n0);
    test_c2c_backward(MPI_COMM_WORLD, 3, 5, 7, 6, 2, long_n2 | long_n0);

    // Cubic domain
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 2, 0, long_n2 | long_n0); // 1
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 2, 1, long_n2 | long_n0);
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 2, 2, long_n2 | long_n0);
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 4, 0, long_n2 | long_n0); // 2
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 4, 1, long_n2 | long_n0);
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 4, 2, long_n2 | long_n0);
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 6, 0, long_n2 | long_n0); // 3
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 6, 1, long_n2 | long_n0);
    test_c2c_backward(MPI_COMM_WORLD, 8, 8, 8, 6, 2, long_n2 | long_n0);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_FIXTURE_TEST_SUITE( underling_fft_c2r,
                          BoostFailErrorHandlerFixture )

// Test wave to physical transformation and inverse transform
void test_c2r(MPI_Comm comm,
              const int n0, const int n1, const int n2,
              const int howmany,
              const int long_ni,
              const unsigned transposed_flags = 0)
{
    int procid;
    BOOST_REQUIRE_EQUAL(MPI_SUCCESS, MPI_Comm_rank(comm, &procid));

    int nproc;
    BOOST_REQUIRE_EQUAL(MPI_SUCCESS, MPI_Comm_size(comm, &nproc));

    if (!procid) {
        BOOST_TEST_MESSAGE("Testing howmany = " << howmany
                           << " on " << n0 << "x" << n1 << "x" << n2
                           << " using " << nproc << " processor"
                           << (nproc > 1 ? "s" : "")
                           << " when long in " << long_ni
                           << " with flags " << transposed_flags);
    }
    UnderlingFixture f(comm, n0, n1, n2, howmany, transposed_flags);

    underling::fft::plan backward(underling::fft::plan::c2r_backward(),
                                  f.problem,
                                  long_ni,
                                  f.data.get(),
                                  FFTW_ESTIMATE);
    BOOST_REQUIRE(backward);
    underling::fft::plan forward(backward,       // Inverse constructor!
                                 f.data.get(),
                                 FFTW_ESTIMATE);
    BOOST_REQUIRE(forward);

    // Stride information consistency check
    BOOST_REQUIRE_EQUAL(backward.local_extents_output(),
                        forward.local_extents_input());
    BOOST_REQUIRE_EQUAL(forward.local_extents_output(),
                        backward.local_extents_input());

    const underling::extents extents = f.problem.local_extents(long_ni);
    const double close_enough
        =   std::numeric_limits<double>::epsilon()*100
          * extents.size[long_ni]*extents.size[long_ni]*extents.size[long_ni];

    // Load up sample data
    {
        const underling::fft::extents e = backward.local_extents_input();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(2*(e.size[long_ni]-1),
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    underling_real * const pencil = &f.data[
                          i*e.stride[dir_i]
                        + j*e.stride[dir_j]
                        + k*e.stride[dir_k]
                    ];

                    for (int l = 0; l < e.size[long_ni]; ++l) {
                        underling_real * const base
                            = pencil + l*e.stride[long_ni];
                        const std::complex<double> val = pf.wave(l);
                        base[0]                        = val.real();
                        base[e.stride[e.order[0]]]     = val.imag();
                    }
                }
            }
        }
    }

    // Transform from wave to physical space
    backward.execute();

    // Check data transformed as expected
    {
        const underling::fft::extents e = backward.local_extents_output();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(e.size[long_ni],
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    underling_real * const base = &f.data[
                          i*e.stride[dir_i]
                        + j*e.stride[dir_j]
                        + k*e.stride[dir_k]
                    ];

                    for (int l = 0; l < e.size[long_ni]; ++l) {
                        const double expected = pf.physical(l);
                        const double actual   = base[ l*e.stride[long_ni] ];
                        BOOST_CHECK_CLOSE(expected, actual, close_enough);
                    }
                }
            }
        }
    }

    // Transform physical to wave space
    forward.execute();

    // Check data transformed as expected
    {
        const underling::fft::extents e = backward.local_extents_input();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(2*(e.size[long_ni]-1),
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    const underling_real * const pencil = &f.data[
                          i*e.stride[dir_i]
                        + j*e.stride[dir_j]
                        + k*e.stride[dir_k]
                    ];

                    for (int l = 0; l < e.size[long_ni]; ++l) {
                        const underling_real * const base
                            = pencil + l*e.stride[long_ni];
                        const std::complex<double> expected
                            = pf.wave(l) * (double) pf.N;
                        BOOST_CHECK_CLOSE(
                                expected.real(),
                                base[0],
                                close_enough);
                        BOOST_CHECK_CLOSE(
                                expected.imag(),
                                base[e.stride[e.order[0]]],
                                close_enough);
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( underling_fft_c2r_simple_n0 )
{
    // Long in n0, transform a single pencil
    test_c2r(MPI_COMM_WORLD, 3, 1, 1, 2, 0);
    test_c2r(MPI_COMM_WORLD, 3, 1, 1, 4, 0);
    test_c2r(MPI_COMM_WORLD, 3, 1, 1, 6, 0);

    test_c2r(MPI_COMM_WORLD, 4, 1, 1, 2, 0);
    test_c2r(MPI_COMM_WORLD, 4, 1, 1, 4, 0);
    test_c2r(MPI_COMM_WORLD, 4, 1, 1, 6, 0);

    test_c2r(MPI_COMM_WORLD, 5, 1, 1, 2, 0);
    test_c2r(MPI_COMM_WORLD, 5, 1, 1, 4, 0);
    test_c2r(MPI_COMM_WORLD, 5, 1, 1, 6, 0);

    test_c2r(MPI_COMM_WORLD, 6, 1, 1, 2, 0);
    test_c2r(MPI_COMM_WORLD, 6, 1, 1, 4, 0);
    test_c2r(MPI_COMM_WORLD, 6, 1, 1, 6, 0);

    // Long in n0, transform multiple pencils
    test_c2r(MPI_COMM_WORLD, 3, 2, 2, 2, 0);
    test_c2r(MPI_COMM_WORLD, 3, 2, 3, 4, 0);
    test_c2r(MPI_COMM_WORLD, 3, 3, 2, 6, 0);

    test_c2r(MPI_COMM_WORLD, 4, 2, 2, 2, 0);
    test_c2r(MPI_COMM_WORLD, 4, 2, 3, 4, 0);
    test_c2r(MPI_COMM_WORLD, 4, 3, 2, 6, 0);

    test_c2r(MPI_COMM_WORLD, 5, 2, 2, 2, 0);
    test_c2r(MPI_COMM_WORLD, 5, 2, 3, 4, 0);
    test_c2r(MPI_COMM_WORLD, 5, 3, 2, 6, 0);

    test_c2r(MPI_COMM_WORLD, 6, 2, 2, 2, 0);
    test_c2r(MPI_COMM_WORLD, 6, 2, 3, 4, 0);
    test_c2r(MPI_COMM_WORLD, 6, 3, 2, 6, 0);
}

BOOST_AUTO_TEST_CASE( underling_fft_c2r_simple_n0_transposed_long_n2 )
{
    ensureFFTWTensor7PatchInPlace();

    using suzerain::underling::transposed::long_n2;

    // Long in n0, transform a single pencil
    test_c2r(MPI_COMM_WORLD, 3, 1, 1, 2, 0, long_n2);
    test_c2r(MPI_COMM_WORLD, 3, 1, 1, 4, 0, long_n2);
    test_c2r(MPI_COMM_WORLD, 3, 1, 1, 6, 0, long_n2);

    test_c2r(MPI_COMM_WORLD, 4, 1, 1, 2, 0, long_n2);
    test_c2r(MPI_COMM_WORLD, 4, 1, 1, 4, 0, long_n2);
    test_c2r(MPI_COMM_WORLD, 4, 1, 1, 6, 0, long_n2);

    test_c2r(MPI_COMM_WORLD, 5, 1, 1, 2, 0, long_n2);
    test_c2r(MPI_COMM_WORLD, 5, 1, 1, 4, 0, long_n2);
    test_c2r(MPI_COMM_WORLD, 5, 1, 1, 6, 0, long_n2);

    test_c2r(MPI_COMM_WORLD, 6, 1, 1, 2, 0, long_n2);
    test_c2r(MPI_COMM_WORLD, 6, 1, 1, 4, 0, long_n2);
    test_c2r(MPI_COMM_WORLD, 6, 1, 1, 6, 0, long_n2);

    // Long in n0, transform multiple pencils
    test_c2r(MPI_COMM_WORLD, 3, 2, 2, 2, 0, long_n2);
    test_c2r(MPI_COMM_WORLD, 3, 2, 3, 4, 0, long_n2);
    test_c2r(MPI_COMM_WORLD, 3, 3, 2, 6, 0, long_n2);

    test_c2r(MPI_COMM_WORLD, 4, 2, 2, 2, 0, long_n2);
    test_c2r(MPI_COMM_WORLD, 4, 2, 3, 4, 0, long_n2);
    test_c2r(MPI_COMM_WORLD, 4, 3, 2, 6, 0, long_n2);

    test_c2r(MPI_COMM_WORLD, 5, 2, 2, 2, 0, long_n2);
    test_c2r(MPI_COMM_WORLD, 5, 2, 3, 4, 0, long_n2);
    test_c2r(MPI_COMM_WORLD, 5, 3, 2, 6, 0, long_n2);

    test_c2r(MPI_COMM_WORLD, 6, 2, 2, 2, 0, long_n2);
    test_c2r(MPI_COMM_WORLD, 6, 2, 3, 4, 0, long_n2);
    test_c2r(MPI_COMM_WORLD, 6, 3, 2, 6, 0, long_n2);
}

BOOST_AUTO_TEST_CASE( underling_fft_c2r_simple_n0_transposed_long_n0 )
{
    ensureFFTWTensor7PatchInPlace();

    using suzerain::underling::transposed::long_n0;

    // Long in n0, transform a single pencil
    test_c2r(MPI_COMM_WORLD, 3, 1, 1, 2, 0, long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 1, 1, 4, 0, long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 1, 1, 6, 0, long_n0);

    test_c2r(MPI_COMM_WORLD, 4, 1, 1, 2, 0, long_n0);
    test_c2r(MPI_COMM_WORLD, 4, 1, 1, 4, 0, long_n0);
    test_c2r(MPI_COMM_WORLD, 4, 1, 1, 6, 0, long_n0);

    test_c2r(MPI_COMM_WORLD, 5, 1, 1, 2, 0, long_n0);
    test_c2r(MPI_COMM_WORLD, 5, 1, 1, 4, 0, long_n0);
    test_c2r(MPI_COMM_WORLD, 5, 1, 1, 6, 0, long_n0);

    test_c2r(MPI_COMM_WORLD, 6, 1, 1, 2, 0, long_n0);
    test_c2r(MPI_COMM_WORLD, 6, 1, 1, 4, 0, long_n0);
    test_c2r(MPI_COMM_WORLD, 6, 1, 1, 6, 0, long_n0);

    // Long in n0, transform multiple pencils
    test_c2r(MPI_COMM_WORLD, 3, 2, 2, 2, 0, long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 2, 3, 4, 0, long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 3, 2, 6, 0, long_n0);

    test_c2r(MPI_COMM_WORLD, 4, 2, 2, 2, 0, long_n0);
    test_c2r(MPI_COMM_WORLD, 4, 2, 3, 4, 0, long_n0);
    test_c2r(MPI_COMM_WORLD, 4, 3, 2, 6, 0, long_n0);

    test_c2r(MPI_COMM_WORLD, 5, 2, 2, 2, 0, long_n0);
    test_c2r(MPI_COMM_WORLD, 5, 2, 3, 4, 0, long_n0);
    test_c2r(MPI_COMM_WORLD, 5, 3, 2, 6, 0, long_n0);

    test_c2r(MPI_COMM_WORLD, 6, 2, 2, 2, 0, long_n0);
    test_c2r(MPI_COMM_WORLD, 6, 2, 3, 4, 0, long_n0);
    test_c2r(MPI_COMM_WORLD, 6, 3, 2, 6, 0, long_n0);
}

BOOST_AUTO_TEST_CASE( underling_fft_c2r_simple_n0_transposed_long_both )
{
    ensureFFTWTensor7PatchInPlace();

    using suzerain::underling::transposed::long_n2;
    using suzerain::underling::transposed::long_n0;

    // Long in n0, transform a single pencil
    test_c2r(MPI_COMM_WORLD, 3, 1, 1, 2, 0, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 1, 1, 4, 0, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 1, 1, 6, 0, long_n2 | long_n0);

    test_c2r(MPI_COMM_WORLD, 4, 1, 1, 2, 0, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 4, 1, 1, 4, 0, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 4, 1, 1, 6, 0, long_n2 | long_n0);

    test_c2r(MPI_COMM_WORLD, 5, 1, 1, 2, 0, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 5, 1, 1, 4, 0, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 5, 1, 1, 6, 0, long_n2 | long_n0);

    test_c2r(MPI_COMM_WORLD, 6, 1, 1, 2, 0, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 6, 1, 1, 4, 0, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 6, 1, 1, 6, 0, long_n2 | long_n0);

    // Long in n0, transform multiple pencils
    test_c2r(MPI_COMM_WORLD, 3, 2, 2, 2, 0, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 2, 3, 4, 0, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 3, 2, 6, 0, long_n2 | long_n0);

    test_c2r(MPI_COMM_WORLD, 4, 2, 2, 2, 0, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 4, 2, 3, 4, 0, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 4, 3, 2, 6, 0, long_n2 | long_n0);

    test_c2r(MPI_COMM_WORLD, 5, 2, 2, 2, 0, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 5, 2, 3, 4, 0, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 5, 3, 2, 6, 0, long_n2 | long_n0);

    test_c2r(MPI_COMM_WORLD, 6, 2, 2, 2, 0, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 6, 2, 3, 4, 0, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 6, 3, 2, 6, 0, long_n2 | long_n0);
}

BOOST_AUTO_TEST_CASE( underling_fft_c2r_simple_n1 )
{
    // Long in n1, transform a single pencil
    test_c2r(MPI_COMM_WORLD, 1, 3, 1, 2, 1);
    test_c2r(MPI_COMM_WORLD, 1, 3, 1, 4, 1);
    test_c2r(MPI_COMM_WORLD, 1, 3, 1, 6, 1);

    test_c2r(MPI_COMM_WORLD, 1, 4, 1, 2, 1);
    test_c2r(MPI_COMM_WORLD, 1, 4, 1, 4, 1);
    test_c2r(MPI_COMM_WORLD, 1, 4, 1, 6, 1);

    test_c2r(MPI_COMM_WORLD, 1, 5, 1, 2, 1);
    test_c2r(MPI_COMM_WORLD, 1, 5, 1, 4, 1);
    test_c2r(MPI_COMM_WORLD, 1, 5, 1, 6, 1);

    test_c2r(MPI_COMM_WORLD, 1, 6, 1, 2, 1);
    test_c2r(MPI_COMM_WORLD, 1, 6, 1, 4, 1);
    test_c2r(MPI_COMM_WORLD, 1, 6, 1, 6, 1);

    // Long in n1, transform multiple pencils
    test_c2r(MPI_COMM_WORLD, 2, 3, 2, 2, 1);
    test_c2r(MPI_COMM_WORLD, 2, 3, 3, 4, 1);
    test_c2r(MPI_COMM_WORLD, 3, 3, 2, 6, 1);

    test_c2r(MPI_COMM_WORLD, 2, 4, 2, 2, 1);
    test_c2r(MPI_COMM_WORLD, 2, 4, 3, 4, 1);
    test_c2r(MPI_COMM_WORLD, 3, 4, 2, 6, 1);

    test_c2r(MPI_COMM_WORLD, 2, 5, 2, 2, 1);
    test_c2r(MPI_COMM_WORLD, 2, 5, 3, 4, 1);
    test_c2r(MPI_COMM_WORLD, 3, 5, 2, 6, 1);

    test_c2r(MPI_COMM_WORLD, 2, 6, 2, 2, 1);
    test_c2r(MPI_COMM_WORLD, 2, 6, 3, 4, 1);
    test_c2r(MPI_COMM_WORLD, 3, 6, 2, 6, 1);
}

BOOST_AUTO_TEST_CASE( underling_fft_c2r_simple_n1_transposed_long_n2 )
{
    using suzerain::underling::transposed::long_n2;

    // Long in n1, transform a single pencil
    test_c2r(MPI_COMM_WORLD, 1, 3, 1, 2, 1, long_n2);
    test_c2r(MPI_COMM_WORLD, 1, 3, 1, 4, 1, long_n2);
    test_c2r(MPI_COMM_WORLD, 1, 3, 1, 6, 1, long_n2);

    test_c2r(MPI_COMM_WORLD, 1, 4, 1, 2, 1, long_n2);
    test_c2r(MPI_COMM_WORLD, 1, 4, 1, 4, 1, long_n2);
    test_c2r(MPI_COMM_WORLD, 1, 4, 1, 6, 1, long_n2);

    test_c2r(MPI_COMM_WORLD, 1, 5, 1, 2, 1, long_n2);
    test_c2r(MPI_COMM_WORLD, 1, 5, 1, 4, 1, long_n2);
    test_c2r(MPI_COMM_WORLD, 1, 5, 1, 6, 1, long_n2);

    test_c2r(MPI_COMM_WORLD, 1, 6, 1, 2, 1, long_n2);
    test_c2r(MPI_COMM_WORLD, 1, 6, 1, 4, 1, long_n2);
    test_c2r(MPI_COMM_WORLD, 1, 6, 1, 6, 1, long_n2);

    // Long in n1, transform multiple pencils
    test_c2r(MPI_COMM_WORLD, 2, 3, 2, 2, 1, long_n2);
    test_c2r(MPI_COMM_WORLD, 2, 3, 3, 4, 1, long_n2);
    test_c2r(MPI_COMM_WORLD, 3, 3, 2, 6, 1, long_n2);

    test_c2r(MPI_COMM_WORLD, 2, 4, 2, 2, 1, long_n2);
    test_c2r(MPI_COMM_WORLD, 2, 4, 3, 4, 1, long_n2);
    test_c2r(MPI_COMM_WORLD, 3, 4, 2, 6, 1, long_n2);

    test_c2r(MPI_COMM_WORLD, 2, 5, 2, 2, 1, long_n2);
    test_c2r(MPI_COMM_WORLD, 2, 5, 3, 4, 1, long_n2);
    test_c2r(MPI_COMM_WORLD, 3, 5, 2, 6, 1, long_n2);

    test_c2r(MPI_COMM_WORLD, 2, 6, 2, 2, 1, long_n2);
    test_c2r(MPI_COMM_WORLD, 2, 6, 3, 4, 1, long_n2);
    test_c2r(MPI_COMM_WORLD, 3, 6, 2, 6, 1, long_n2);
}

BOOST_AUTO_TEST_CASE( underling_fft_c2r_simple_n1_transposed_long_n0 )
{
    using suzerain::underling::transposed::long_n0;

    // Long in n1, transform a single pencil
    test_c2r(MPI_COMM_WORLD, 1, 3, 1, 2, 1, long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 3, 1, 4, 1, long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 3, 1, 6, 1, long_n0);

    test_c2r(MPI_COMM_WORLD, 1, 4, 1, 2, 1, long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 4, 1, 4, 1, long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 4, 1, 6, 1, long_n0);

    test_c2r(MPI_COMM_WORLD, 1, 5, 1, 2, 1, long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 5, 1, 4, 1, long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 5, 1, 6, 1, long_n0);

    test_c2r(MPI_COMM_WORLD, 1, 6, 1, 2, 1, long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 6, 1, 4, 1, long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 6, 1, 6, 1, long_n0);

    // Long in n1, transform multiple pencils
    test_c2r(MPI_COMM_WORLD, 2, 3, 2, 2, 1, long_n0);
    test_c2r(MPI_COMM_WORLD, 2, 3, 3, 4, 1, long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 3, 2, 6, 1, long_n0);

    test_c2r(MPI_COMM_WORLD, 2, 4, 2, 2, 1, long_n0);
    test_c2r(MPI_COMM_WORLD, 2, 4, 3, 4, 1, long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 4, 2, 6, 1, long_n0);

    test_c2r(MPI_COMM_WORLD, 2, 5, 2, 2, 1, long_n0);
    test_c2r(MPI_COMM_WORLD, 2, 5, 3, 4, 1, long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 5, 2, 6, 1, long_n0);

    test_c2r(MPI_COMM_WORLD, 2, 6, 2, 2, 1, long_n0);
    test_c2r(MPI_COMM_WORLD, 2, 6, 3, 4, 1, long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 6, 2, 6, 1, long_n0);
}

BOOST_AUTO_TEST_CASE( underling_fft_c2r_simple_n1_transposed_both )
{
    using suzerain::underling::transposed::long_n2;
    using suzerain::underling::transposed::long_n0;

    // Long in n1, transform a single pencil
    test_c2r(MPI_COMM_WORLD, 1, 3, 1, 2, 1, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 3, 1, 4, 1, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 3, 1, 6, 1, long_n2 | long_n0);

    test_c2r(MPI_COMM_WORLD, 1, 4, 1, 2, 1, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 4, 1, 4, 1, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 4, 1, 6, 1, long_n2 | long_n0);

    test_c2r(MPI_COMM_WORLD, 1, 5, 1, 2, 1, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 5, 1, 4, 1, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 5, 1, 6, 1, long_n2 | long_n0);

    test_c2r(MPI_COMM_WORLD, 1, 6, 1, 2, 1, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 6, 1, 4, 1, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 6, 1, 6, 1, long_n2 | long_n0);

    // Long in n1, transform multiple pencils
    test_c2r(MPI_COMM_WORLD, 2, 3, 2, 2, 1, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 2, 3, 3, 4, 1, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 3, 2, 6, 1, long_n2 | long_n0);

    test_c2r(MPI_COMM_WORLD, 2, 4, 2, 2, 1, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 2, 4, 3, 4, 1, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 4, 2, 6, 1, long_n2 | long_n0);

    test_c2r(MPI_COMM_WORLD, 2, 5, 2, 2, 1, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 2, 5, 3, 4, 1, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 5, 2, 6, 1, long_n2 | long_n0);

    test_c2r(MPI_COMM_WORLD, 2, 6, 2, 2, 1, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 2, 6, 3, 4, 1, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 6, 2, 6, 1, long_n2 | long_n0);
}

BOOST_AUTO_TEST_CASE( underling_fft_c2r_simple_n2 )
{
    // Long in n2, transform a single pencil
    test_c2r(MPI_COMM_WORLD, 1, 1, 3, 2, 2);
    test_c2r(MPI_COMM_WORLD, 1, 1, 3, 4, 2);
    test_c2r(MPI_COMM_WORLD, 1, 1, 3, 6, 2);

    test_c2r(MPI_COMM_WORLD, 1, 1, 4, 2, 2);
    test_c2r(MPI_COMM_WORLD, 1, 1, 4, 4, 2);
    test_c2r(MPI_COMM_WORLD, 1, 1, 4, 6, 2);

    test_c2r(MPI_COMM_WORLD, 1, 1, 5, 2, 2);
    test_c2r(MPI_COMM_WORLD, 1, 1, 5, 4, 2);
    test_c2r(MPI_COMM_WORLD, 1, 1, 5, 6, 2);

    test_c2r(MPI_COMM_WORLD, 1, 1, 6, 2, 2);
    test_c2r(MPI_COMM_WORLD, 1, 1, 6, 4, 2);
    test_c2r(MPI_COMM_WORLD, 1, 1, 6, 6, 2);

    // Long in n2, transform multiple pencils
    test_c2r(MPI_COMM_WORLD, 2, 2, 3, 2, 2);
    test_c2r(MPI_COMM_WORLD, 2, 3, 3, 4, 2);
    test_c2r(MPI_COMM_WORLD, 3, 2, 3, 6, 2);

    test_c2r(MPI_COMM_WORLD, 2, 2, 4, 2, 2);
    test_c2r(MPI_COMM_WORLD, 2, 3, 4, 4, 2);
    test_c2r(MPI_COMM_WORLD, 3, 2, 4, 6, 2);

    test_c2r(MPI_COMM_WORLD, 2, 2, 5, 2, 2);
    test_c2r(MPI_COMM_WORLD, 2, 3, 5, 4, 2);
    test_c2r(MPI_COMM_WORLD, 3, 2, 5, 6, 2);

    test_c2r(MPI_COMM_WORLD, 2, 2, 6, 2, 2);
    test_c2r(MPI_COMM_WORLD, 2, 3, 6, 4, 2);
    test_c2r(MPI_COMM_WORLD, 3, 2, 6, 6, 2);
}

BOOST_AUTO_TEST_CASE( underling_fft_c2r_simple_n2_transposed_long_n2 )
{
    ensureFFTWTensor7PatchInPlace();

    using suzerain::underling::transposed::long_n2;

    // Long in n2, transform a single pencil
    test_c2r(MPI_COMM_WORLD, 1, 1, 3, 2, 2, long_n2);
    test_c2r(MPI_COMM_WORLD, 1, 1, 3, 4, 2, long_n2);
    test_c2r(MPI_COMM_WORLD, 1, 1, 3, 6, 2, long_n2);

    test_c2r(MPI_COMM_WORLD, 1, 1, 4, 2, 2, long_n2);
    test_c2r(MPI_COMM_WORLD, 1, 1, 4, 4, 2, long_n2);
    test_c2r(MPI_COMM_WORLD, 1, 1, 4, 6, 2, long_n2);

    test_c2r(MPI_COMM_WORLD, 1, 1, 5, 2, 2, long_n2);
    test_c2r(MPI_COMM_WORLD, 1, 1, 5, 4, 2, long_n2);
    test_c2r(MPI_COMM_WORLD, 1, 1, 5, 6, 2, long_n2);

    test_c2r(MPI_COMM_WORLD, 1, 1, 6, 2, 2, long_n2);
    test_c2r(MPI_COMM_WORLD, 1, 1, 6, 4, 2, long_n2);
    test_c2r(MPI_COMM_WORLD, 1, 1, 6, 6, 2, long_n2);

    // Long in n2, transform multiple pencils
    test_c2r(MPI_COMM_WORLD, 2, 2, 3, 2, 2, long_n2);
    test_c2r(MPI_COMM_WORLD, 2, 3, 3, 4, 2, long_n2);
    test_c2r(MPI_COMM_WORLD, 3, 2, 3, 6, 2, long_n2);

    test_c2r(MPI_COMM_WORLD, 2, 2, 4, 2, 2, long_n2);
    test_c2r(MPI_COMM_WORLD, 2, 3, 4, 4, 2, long_n2);
    test_c2r(MPI_COMM_WORLD, 3, 2, 4, 6, 2, long_n2);

    test_c2r(MPI_COMM_WORLD, 2, 2, 5, 2, 2, long_n2);
    test_c2r(MPI_COMM_WORLD, 2, 3, 5, 4, 2, long_n2);
    test_c2r(MPI_COMM_WORLD, 3, 2, 5, 6, 2, long_n2);

    test_c2r(MPI_COMM_WORLD, 2, 2, 6, 2, 2, long_n2);
    test_c2r(MPI_COMM_WORLD, 2, 3, 6, 4, 2, long_n2);
    test_c2r(MPI_COMM_WORLD, 3, 2, 6, 6, 2, long_n2);
}

BOOST_AUTO_TEST_CASE( underling_fft_c2r_simple_n2_transposed_long_n0 )
{
    using suzerain::underling::transposed::long_n0;

    // Long in n2, transform a single pencil
    test_c2r(MPI_COMM_WORLD, 1, 1, 3, 2, 2, long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 1, 3, 4, 2, long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 1, 3, 6, 2, long_n0);

    test_c2r(MPI_COMM_WORLD, 1, 1, 4, 2, 2, long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 1, 4, 4, 2, long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 1, 4, 6, 2, long_n0);

    test_c2r(MPI_COMM_WORLD, 1, 1, 5, 2, 2, long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 1, 5, 4, 2, long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 1, 5, 6, 2, long_n0);

    test_c2r(MPI_COMM_WORLD, 1, 1, 6, 2, 2, long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 1, 6, 4, 2, long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 1, 6, 6, 2, long_n0);

    // Long in n2, transform multiple pencils
    test_c2r(MPI_COMM_WORLD, 2, 2, 3, 2, 2, long_n0);
    test_c2r(MPI_COMM_WORLD, 2, 3, 3, 4, 2, long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 2, 3, 6, 2, long_n0);

    test_c2r(MPI_COMM_WORLD, 2, 2, 4, 2, 2, long_n0);
    test_c2r(MPI_COMM_WORLD, 2, 3, 4, 4, 2, long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 2, 4, 6, 2, long_n0);

    test_c2r(MPI_COMM_WORLD, 2, 2, 5, 2, 2, long_n0);
    test_c2r(MPI_COMM_WORLD, 2, 3, 5, 4, 2, long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 2, 5, 6, 2, long_n0);

    test_c2r(MPI_COMM_WORLD, 2, 2, 6, 2, 2, long_n0);
    test_c2r(MPI_COMM_WORLD, 2, 3, 6, 4, 2, long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 2, 6, 6, 2, long_n0);
}

BOOST_AUTO_TEST_CASE( underling_fft_c2r_simple_n2_transposed_both )
{
    ensureFFTWTensor7PatchInPlace();

    using suzerain::underling::transposed::long_n2;
    using suzerain::underling::transposed::long_n0;

    // Long in n2, transform a single pencil
    test_c2r(MPI_COMM_WORLD, 1, 1, 3, 2, 2, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 1, 3, 4, 2, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 1, 3, 6, 2, long_n2 | long_n0);

    test_c2r(MPI_COMM_WORLD, 1, 1, 4, 2, 2, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 1, 4, 4, 2, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 1, 4, 6, 2, long_n2 | long_n0);

    test_c2r(MPI_COMM_WORLD, 1, 1, 5, 2, 2, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 1, 5, 4, 2, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 1, 5, 6, 2, long_n2 | long_n0);

    test_c2r(MPI_COMM_WORLD, 1, 1, 6, 2, 2, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 1, 6, 4, 2, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 1, 1, 6, 6, 2, long_n2 | long_n0);

    // Long in n2, transform multiple pencils
    test_c2r(MPI_COMM_WORLD, 2, 2, 3, 2, 2, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 2, 3, 3, 4, 2, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 2, 3, 6, 2, long_n2 | long_n0);

    test_c2r(MPI_COMM_WORLD, 2, 2, 4, 2, 2, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 2, 3, 4, 4, 2, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 2, 4, 6, 2, long_n2 | long_n0);

    test_c2r(MPI_COMM_WORLD, 2, 2, 5, 2, 2, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 2, 3, 5, 4, 2, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 2, 5, 6, 2, long_n2 | long_n0);

    test_c2r(MPI_COMM_WORLD, 2, 2, 6, 2, 2, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 2, 3, 6, 4, 2, long_n2 | long_n0);
    test_c2r(MPI_COMM_WORLD, 3, 2, 6, 6, 2, long_n2 | long_n0);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_FIXTURE_TEST_SUITE( underling_fft_r2c,
                          BoostFailErrorHandlerFixture )

// Test physical to wave transformation and inverse transform
void test_r2c(MPI_Comm comm,
              const int n0, const int n1, const int n2,
              const int howmany,
              const int long_ni,
              const unsigned transposed_flags = 0)
{
    int procid;
    BOOST_REQUIRE_EQUAL(MPI_SUCCESS, MPI_Comm_rank(comm, &procid));

    int nproc;
    BOOST_REQUIRE_EQUAL(MPI_SUCCESS, MPI_Comm_size(comm, &nproc));

    if (!procid) {
        BOOST_TEST_MESSAGE("Testing howmany = " << howmany
                           << " on " << n0 << "x" << n1 << "x" << n2
                           << " using " << nproc << " processor"
                           << (nproc > 1 ? "s" : "")
                           << " when long in " << long_ni
                           << " with flags " << transposed_flags);
    }
    UnderlingFixture f(comm, n0, n1, n2, howmany, transposed_flags);

    underling::fft::plan forward(underling::fft::plan::r2c_forward(),
                                 f.problem,
                                 long_ni,
                                 f.data.get(),
                                 FFTW_ESTIMATE);
    BOOST_REQUIRE(forward);
    underling::fft::plan backward(forward,        // Inverse constructor!
                                  f.data.get(),
                                  FFTW_ESTIMATE);
    BOOST_REQUIRE(backward);

    // Stride information consistency check
    BOOST_REQUIRE_EQUAL(forward.local_extents_output(),
                        backward.local_extents_input());
    BOOST_REQUIRE_EQUAL(backward.local_extents_output(),
                        forward.local_extents_input());

    const underling::extents extents = f.problem.local_extents(long_ni);
    const double close_enough
        =   std::numeric_limits<double>::epsilon()*100
          * extents.size[long_ni]*extents.size[long_ni]*extents.size[long_ni];

    // Load up sample data
    {
        const underling::fft::extents e = forward.local_extents_input();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(e.size[long_ni],
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    underling_real * const base = &f.data[
                          i*e.stride[dir_i]
                        + j*e.stride[dir_j]
                        + k*e.stride[dir_k]
                    ];

                    for (int l = 0; l < e.size[long_ni]; ++l) {
                        base[ l*e.stride[long_ni] ] = pf.physical(l);
                    }
                }
            }
        }
    }

    // Transform physical to wave space
    forward.execute();

    // Check data transformed as expected
    {
        const underling::fft::extents e = forward.local_extents_output();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(2*(e.size[long_ni]-1),
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    const underling_real * const pencil = &f.data[
                          i*e.stride[dir_i]
                        + j*e.stride[dir_j]
                        + k*e.stride[dir_k]
                    ];

                    for (int l = 0; l < e.size[long_ni]; ++l) {
                        const underling_real * const base
                            = pencil + l*e.stride[long_ni];
                        const std::complex<double> expected
                            = pf.wave(l) * (double) pf.N;
                        BOOST_CHECK_CLOSE(
                                expected.real(),
                                base[0],
                                close_enough);
                        BOOST_CHECK_CLOSE(
                                expected.imag(),
                                base[e.stride[e.order[0]]],
                                close_enough);
                    }
                }
            }
        }
    }

    // Transform wave to physical
    backward.execute();

    // Checking data transformed as expected
    {
        const underling::fft::extents e = forward.local_extents_input();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(e.size[long_ni],
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    underling_real * const base = &f.data[
                          i*e.stride[dir_i]
                        + j*e.stride[dir_j]
                        + k*e.stride[dir_k]
                    ];

                    for (int l = 0; l < e.size[long_ni]; ++l) {
                        const double expected
                            = pf.physical(l) * pf.N;
                        const double actual
                            = base[ l*e.stride[long_ni] ];
                        BOOST_CHECK_CLOSE(expected, actual, close_enough);
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( underling_fft_r2c_simple_n0 )
{
    // Long in n0, transform a single pencil
    test_r2c(MPI_COMM_WORLD, 3, 1, 1, 2, 0);
    test_r2c(MPI_COMM_WORLD, 3, 1, 1, 4, 0);
    test_r2c(MPI_COMM_WORLD, 3, 1, 1, 6, 0);

    test_r2c(MPI_COMM_WORLD, 4, 1, 1, 2, 0);
    test_r2c(MPI_COMM_WORLD, 4, 1, 1, 4, 0);
    test_r2c(MPI_COMM_WORLD, 4, 1, 1, 6, 0);

    test_r2c(MPI_COMM_WORLD, 5, 1, 1, 2, 0);
    test_r2c(MPI_COMM_WORLD, 5, 1, 1, 4, 0);
    test_r2c(MPI_COMM_WORLD, 5, 1, 1, 6, 0);

    test_r2c(MPI_COMM_WORLD, 6, 1, 1, 2, 0);
    test_r2c(MPI_COMM_WORLD, 6, 1, 1, 4, 0);
    test_r2c(MPI_COMM_WORLD, 6, 1, 1, 6, 0);

    // Long in n0, transform multiple pencils
    test_r2c(MPI_COMM_WORLD, 3, 2, 2, 2, 0);
    test_r2c(MPI_COMM_WORLD, 3, 2, 3, 4, 0);
    test_r2c(MPI_COMM_WORLD, 3, 3, 2, 6, 0);

    test_r2c(MPI_COMM_WORLD, 4, 2, 2, 2, 0);
    test_r2c(MPI_COMM_WORLD, 4, 2, 3, 4, 0);
    test_r2c(MPI_COMM_WORLD, 4, 3, 2, 6, 0);

    test_r2c(MPI_COMM_WORLD, 5, 2, 2, 2, 0);
    test_r2c(MPI_COMM_WORLD, 5, 2, 3, 4, 0);
    test_r2c(MPI_COMM_WORLD, 5, 3, 2, 6, 0);

    test_r2c(MPI_COMM_WORLD, 6, 2, 2, 2, 0);
    test_r2c(MPI_COMM_WORLD, 6, 2, 3, 4, 0);
    test_r2c(MPI_COMM_WORLD, 6, 3, 2, 6, 0);
}

BOOST_AUTO_TEST_CASE( underling_fft_r2c_simple_n0_transposed_long_n2 )
{
    using suzerain::underling::transposed::long_n2;

    // Long in n0, transform a single pencil
    test_r2c(MPI_COMM_WORLD, 3, 1, 1, 2, 0, long_n2);
    test_r2c(MPI_COMM_WORLD, 3, 1, 1, 4, 0, long_n2);
    test_r2c(MPI_COMM_WORLD, 3, 1, 1, 6, 0, long_n2);

    test_r2c(MPI_COMM_WORLD, 4, 1, 1, 2, 0, long_n2);
    test_r2c(MPI_COMM_WORLD, 4, 1, 1, 4, 0, long_n2);
    test_r2c(MPI_COMM_WORLD, 4, 1, 1, 6, 0, long_n2);

    test_r2c(MPI_COMM_WORLD, 5, 1, 1, 2, 0, long_n2);
    test_r2c(MPI_COMM_WORLD, 5, 1, 1, 4, 0, long_n2);
    test_r2c(MPI_COMM_WORLD, 5, 1, 1, 6, 0, long_n2);

    test_r2c(MPI_COMM_WORLD, 6, 1, 1, 2, 0, long_n2);
    test_r2c(MPI_COMM_WORLD, 6, 1, 1, 4, 0, long_n2);
    test_r2c(MPI_COMM_WORLD, 6, 1, 1, 6, 0, long_n2);

    // Long in n0, transform multiple pencils
    test_r2c(MPI_COMM_WORLD, 3, 2, 2, 2, 0, long_n2);
    test_r2c(MPI_COMM_WORLD, 3, 2, 3, 4, 0, long_n2);
    test_r2c(MPI_COMM_WORLD, 3, 3, 2, 6, 0, long_n2);

    test_r2c(MPI_COMM_WORLD, 4, 2, 2, 2, 0, long_n2);
    test_r2c(MPI_COMM_WORLD, 4, 2, 3, 4, 0, long_n2);
    test_r2c(MPI_COMM_WORLD, 4, 3, 2, 6, 0, long_n2);

    test_r2c(MPI_COMM_WORLD, 5, 2, 2, 2, 0, long_n2);
    test_r2c(MPI_COMM_WORLD, 5, 2, 3, 4, 0, long_n2);
    test_r2c(MPI_COMM_WORLD, 5, 3, 2, 6, 0, long_n2);

    test_r2c(MPI_COMM_WORLD, 6, 2, 2, 2, 0, long_n2);
    test_r2c(MPI_COMM_WORLD, 6, 2, 3, 4, 0, long_n2);
    test_r2c(MPI_COMM_WORLD, 6, 3, 2, 6, 0, long_n2);
}

BOOST_AUTO_TEST_CASE( underling_fft_r2c_simple_n0_transposed_long_n0 )
{
    ensureFFTWTensor7PatchInPlace();

    using suzerain::underling::transposed::long_n0;

    // Long in n0, transform a single pencil
    test_r2c(MPI_COMM_WORLD, 3, 1, 1, 2, 0, long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 1, 1, 4, 0, long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 1, 1, 6, 0, long_n0);

    test_r2c(MPI_COMM_WORLD, 4, 1, 1, 2, 0, long_n0);
    test_r2c(MPI_COMM_WORLD, 4, 1, 1, 4, 0, long_n0);
    test_r2c(MPI_COMM_WORLD, 4, 1, 1, 6, 0, long_n0);

    test_r2c(MPI_COMM_WORLD, 5, 1, 1, 2, 0, long_n0);
    test_r2c(MPI_COMM_WORLD, 5, 1, 1, 4, 0, long_n0);
    test_r2c(MPI_COMM_WORLD, 5, 1, 1, 6, 0, long_n0);

    test_r2c(MPI_COMM_WORLD, 6, 1, 1, 2, 0, long_n0);
    test_r2c(MPI_COMM_WORLD, 6, 1, 1, 4, 0, long_n0);
    test_r2c(MPI_COMM_WORLD, 6, 1, 1, 6, 0, long_n0);

    // Long in n0, transform multiple pencils
    test_r2c(MPI_COMM_WORLD, 3, 2, 2, 2, 0, long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 2, 3, 4, 0, long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 3, 2, 6, 0, long_n0);

    test_r2c(MPI_COMM_WORLD, 4, 2, 2, 2, 0, long_n0);
    test_r2c(MPI_COMM_WORLD, 4, 2, 3, 4, 0, long_n0);
    test_r2c(MPI_COMM_WORLD, 4, 3, 2, 6, 0, long_n0);

    test_r2c(MPI_COMM_WORLD, 5, 2, 2, 2, 0, long_n0);
    test_r2c(MPI_COMM_WORLD, 5, 2, 3, 4, 0, long_n0);
    test_r2c(MPI_COMM_WORLD, 5, 3, 2, 6, 0, long_n0);

    test_r2c(MPI_COMM_WORLD, 6, 2, 2, 2, 0, long_n0);
    test_r2c(MPI_COMM_WORLD, 6, 2, 3, 4, 0, long_n0);
    test_r2c(MPI_COMM_WORLD, 6, 3, 2, 6, 0, long_n0);
}

BOOST_AUTO_TEST_CASE( underling_fft_r2c_simple_n0_transposed_both )
{
    ensureFFTWTensor7PatchInPlace();

    using suzerain::underling::transposed::long_n2;
    using suzerain::underling::transposed::long_n0;

    // Long in n0, transform a single pencil
    test_r2c(MPI_COMM_WORLD, 3, 1, 1, 2, 0, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 1, 1, 4, 0, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 1, 1, 6, 0, long_n2 | long_n0);

    test_r2c(MPI_COMM_WORLD, 4, 1, 1, 2, 0, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 4, 1, 1, 4, 0, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 4, 1, 1, 6, 0, long_n2 | long_n0);

    test_r2c(MPI_COMM_WORLD, 5, 1, 1, 2, 0, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 5, 1, 1, 4, 0, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 5, 1, 1, 6, 0, long_n2 | long_n0);

    test_r2c(MPI_COMM_WORLD, 6, 1, 1, 2, 0, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 6, 1, 1, 4, 0, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 6, 1, 1, 6, 0, long_n2 | long_n0);

    // Long in n0, transform multiple pencils
    test_r2c(MPI_COMM_WORLD, 3, 2, 2, 2, 0, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 2, 3, 4, 0, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 3, 2, 6, 0, long_n2 | long_n0);

    test_r2c(MPI_COMM_WORLD, 4, 2, 2, 2, 0, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 4, 2, 3, 4, 0, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 4, 3, 2, 6, 0, long_n2 | long_n0);

    test_r2c(MPI_COMM_WORLD, 5, 2, 2, 2, 0, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 5, 2, 3, 4, 0, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 5, 3, 2, 6, 0, long_n2 | long_n0);

    test_r2c(MPI_COMM_WORLD, 6, 2, 2, 2, 0, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 6, 2, 3, 4, 0, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 6, 3, 2, 6, 0, long_n2 | long_n0);
}

BOOST_AUTO_TEST_CASE( underling_fft_r2c_simple_n1 )
{
    // Long in n1, transform a single pencil
    test_r2c(MPI_COMM_WORLD, 1, 3, 1, 2, 1);
    test_r2c(MPI_COMM_WORLD, 1, 3, 1, 4, 1);
    test_r2c(MPI_COMM_WORLD, 1, 3, 1, 6, 1);

    test_r2c(MPI_COMM_WORLD, 1, 4, 1, 2, 1);
    test_r2c(MPI_COMM_WORLD, 1, 4, 1, 4, 1);
    test_r2c(MPI_COMM_WORLD, 1, 4, 1, 6, 1);

    test_r2c(MPI_COMM_WORLD, 1, 5, 1, 2, 1);
    test_r2c(MPI_COMM_WORLD, 1, 5, 1, 4, 1);
    test_r2c(MPI_COMM_WORLD, 1, 5, 1, 6, 1);

    test_r2c(MPI_COMM_WORLD, 1, 6, 1, 2, 1);
    test_r2c(MPI_COMM_WORLD, 1, 6, 1, 4, 1);
    test_r2c(MPI_COMM_WORLD, 1, 6, 1, 6, 1);

    // Long in n1, transform multiple pencils
    test_r2c(MPI_COMM_WORLD, 2, 3, 2, 2, 1);
    test_r2c(MPI_COMM_WORLD, 2, 3, 3, 4, 1);
    test_r2c(MPI_COMM_WORLD, 3, 3, 2, 6, 1);

    test_r2c(MPI_COMM_WORLD, 2, 4, 2, 2, 1);
    test_r2c(MPI_COMM_WORLD, 2, 4, 3, 4, 1);
    test_r2c(MPI_COMM_WORLD, 3, 4, 2, 6, 1);

    test_r2c(MPI_COMM_WORLD, 2, 5, 2, 2, 1);
    test_r2c(MPI_COMM_WORLD, 2, 5, 3, 4, 1);
    test_r2c(MPI_COMM_WORLD, 3, 5, 2, 6, 1);

    test_r2c(MPI_COMM_WORLD, 2, 6, 2, 2, 1);
    test_r2c(MPI_COMM_WORLD, 2, 6, 3, 4, 1);
    test_r2c(MPI_COMM_WORLD, 3, 6, 2, 6, 1);
}

BOOST_AUTO_TEST_CASE( underling_fft_r2c_simple_n1_transposed_long_n2 )
{
    using suzerain::underling::transposed::long_n2;

    // Long in n1, transform a single pencil
    test_r2c(MPI_COMM_WORLD, 1, 3, 1, 2, 1, long_n2);
    test_r2c(MPI_COMM_WORLD, 1, 3, 1, 4, 1, long_n2);
    test_r2c(MPI_COMM_WORLD, 1, 3, 1, 6, 1, long_n2);

    test_r2c(MPI_COMM_WORLD, 1, 4, 1, 2, 1, long_n2);
    test_r2c(MPI_COMM_WORLD, 1, 4, 1, 4, 1, long_n2);
    test_r2c(MPI_COMM_WORLD, 1, 4, 1, 6, 1, long_n2);

    test_r2c(MPI_COMM_WORLD, 1, 5, 1, 2, 1, long_n2);
    test_r2c(MPI_COMM_WORLD, 1, 5, 1, 4, 1, long_n2);
    test_r2c(MPI_COMM_WORLD, 1, 5, 1, 6, 1, long_n2);

    test_r2c(MPI_COMM_WORLD, 1, 6, 1, 2, 1, long_n2);
    test_r2c(MPI_COMM_WORLD, 1, 6, 1, 4, 1, long_n2);
    test_r2c(MPI_COMM_WORLD, 1, 6, 1, 6, 1, long_n2);

    // Long in n1, transform multiple pencils
    test_r2c(MPI_COMM_WORLD, 2, 3, 2, 2, 1, long_n2);
    test_r2c(MPI_COMM_WORLD, 2, 3, 3, 4, 1, long_n2);
    test_r2c(MPI_COMM_WORLD, 3, 3, 2, 6, 1, long_n2);

    test_r2c(MPI_COMM_WORLD, 2, 4, 2, 2, 1, long_n2);
    test_r2c(MPI_COMM_WORLD, 2, 4, 3, 4, 1, long_n2);
    test_r2c(MPI_COMM_WORLD, 3, 4, 2, 6, 1, long_n2);

    test_r2c(MPI_COMM_WORLD, 2, 5, 2, 2, 1, long_n2);
    test_r2c(MPI_COMM_WORLD, 2, 5, 3, 4, 1, long_n2);
    test_r2c(MPI_COMM_WORLD, 3, 5, 2, 6, 1, long_n2);

    test_r2c(MPI_COMM_WORLD, 2, 6, 2, 2, 1, long_n2);
    test_r2c(MPI_COMM_WORLD, 2, 6, 3, 4, 1, long_n2);
    test_r2c(MPI_COMM_WORLD, 3, 6, 2, 6, 1, long_n2);
}

BOOST_AUTO_TEST_CASE( underling_fft_r2c_simple_n1_transposed_long_n0 )
{
    using suzerain::underling::transposed::long_n0;

    // Long in n1, transform a single pencil
    test_r2c(MPI_COMM_WORLD, 1, 3, 1, 2, 1, long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 3, 1, 4, 1, long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 3, 1, 6, 1, long_n0);

    test_r2c(MPI_COMM_WORLD, 1, 4, 1, 2, 1, long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 4, 1, 4, 1, long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 4, 1, 6, 1, long_n0);

    test_r2c(MPI_COMM_WORLD, 1, 5, 1, 2, 1, long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 5, 1, 4, 1, long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 5, 1, 6, 1, long_n0);

    test_r2c(MPI_COMM_WORLD, 1, 6, 1, 2, 1, long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 6, 1, 4, 1, long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 6, 1, 6, 1, long_n0);

    // Long in n1, transform multiple pencils
    test_r2c(MPI_COMM_WORLD, 2, 3, 2, 2, 1, long_n0);
    test_r2c(MPI_COMM_WORLD, 2, 3, 3, 4, 1, long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 3, 2, 6, 1, long_n0);

    test_r2c(MPI_COMM_WORLD, 2, 4, 2, 2, 1, long_n0);
    test_r2c(MPI_COMM_WORLD, 2, 4, 3, 4, 1, long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 4, 2, 6, 1, long_n0);

    test_r2c(MPI_COMM_WORLD, 2, 5, 2, 2, 1, long_n0);
    test_r2c(MPI_COMM_WORLD, 2, 5, 3, 4, 1, long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 5, 2, 6, 1, long_n0);

    test_r2c(MPI_COMM_WORLD, 2, 6, 2, 2, 1, long_n0);
    test_r2c(MPI_COMM_WORLD, 2, 6, 3, 4, 1, long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 6, 2, 6, 1, long_n0);
}

BOOST_AUTO_TEST_CASE( underling_fft_r2c_simple_n1_transposed_both )
{
    using suzerain::underling::transposed::long_n2;
    using suzerain::underling::transposed::long_n0;

    // Long in n1, transform a single pencil
    test_r2c(MPI_COMM_WORLD, 1, 3, 1, 2, 1, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 3, 1, 4, 1, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 3, 1, 6, 1, long_n2 | long_n0);

    test_r2c(MPI_COMM_WORLD, 1, 4, 1, 2, 1, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 4, 1, 4, 1, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 4, 1, 6, 1, long_n2 | long_n0);

    test_r2c(MPI_COMM_WORLD, 1, 5, 1, 2, 1, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 5, 1, 4, 1, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 5, 1, 6, 1, long_n2 | long_n0);

    test_r2c(MPI_COMM_WORLD, 1, 6, 1, 2, 1, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 6, 1, 4, 1, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 6, 1, 6, 1, long_n2 | long_n0);

    // Long in n1, transform multiple pencils
    test_r2c(MPI_COMM_WORLD, 2, 3, 2, 2, 1, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 2, 3, 3, 4, 1, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 3, 2, 6, 1, long_n2 | long_n0);

    test_r2c(MPI_COMM_WORLD, 2, 4, 2, 2, 1, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 2, 4, 3, 4, 1, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 4, 2, 6, 1, long_n2 | long_n0);

    test_r2c(MPI_COMM_WORLD, 2, 5, 2, 2, 1, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 2, 5, 3, 4, 1, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 5, 2, 6, 1, long_n2 | long_n0);

    test_r2c(MPI_COMM_WORLD, 2, 6, 2, 2, 1, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 2, 6, 3, 4, 1, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 6, 2, 6, 1, long_n2 | long_n0);
}

BOOST_AUTO_TEST_CASE( underling_fft_r2c_simple_n2 )
{
    // Long in n2, transform a single pencil
    test_r2c(MPI_COMM_WORLD, 1, 1, 3, 2, 2);
    test_r2c(MPI_COMM_WORLD, 1, 1, 3, 4, 2);
    test_r2c(MPI_COMM_WORLD, 1, 1, 3, 6, 2);

    test_r2c(MPI_COMM_WORLD, 1, 1, 4, 2, 2);
    test_r2c(MPI_COMM_WORLD, 1, 1, 4, 4, 2);
    test_r2c(MPI_COMM_WORLD, 1, 1, 4, 6, 2);

    test_r2c(MPI_COMM_WORLD, 1, 1, 5, 2, 2);
    test_r2c(MPI_COMM_WORLD, 1, 1, 5, 4, 2);
    test_r2c(MPI_COMM_WORLD, 1, 1, 5, 6, 2);

    test_r2c(MPI_COMM_WORLD, 1, 1, 6, 2, 2);
    test_r2c(MPI_COMM_WORLD, 1, 1, 6, 4, 2);
    test_r2c(MPI_COMM_WORLD, 1, 1, 6, 6, 2);

    // Long in n2, transform multiple pencils
    test_r2c(MPI_COMM_WORLD, 2, 2, 3, 2, 2);
    test_r2c(MPI_COMM_WORLD, 2, 3, 3, 4, 2);
    test_r2c(MPI_COMM_WORLD, 3, 2, 3, 6, 2);

    test_r2c(MPI_COMM_WORLD, 2, 2, 4, 2, 2);
    test_r2c(MPI_COMM_WORLD, 2, 3, 4, 4, 2);
    test_r2c(MPI_COMM_WORLD, 3, 2, 4, 6, 2);

    test_r2c(MPI_COMM_WORLD, 2, 2, 5, 2, 2);
    test_r2c(MPI_COMM_WORLD, 2, 3, 5, 4, 2);
    test_r2c(MPI_COMM_WORLD, 3, 2, 5, 6, 2);

    test_r2c(MPI_COMM_WORLD, 2, 2, 6, 2, 2);
    test_r2c(MPI_COMM_WORLD, 2, 3, 6, 4, 2);
    test_r2c(MPI_COMM_WORLD, 3, 2, 6, 6, 2);
}

BOOST_AUTO_TEST_CASE( underling_fft_r2c_simple_n2_transposed_long_n2 )
{
    ensureFFTWTensor7PatchInPlace();

    using suzerain::underling::transposed::long_n2;

    // Long in n2, transform a single pencil
    test_r2c(MPI_COMM_WORLD, 1, 1, 3, 2, 2, long_n2);
    test_r2c(MPI_COMM_WORLD, 1, 1, 3, 4, 2, long_n2);
    test_r2c(MPI_COMM_WORLD, 1, 1, 3, 6, 2, long_n2);

    test_r2c(MPI_COMM_WORLD, 1, 1, 4, 2, 2, long_n2);
    test_r2c(MPI_COMM_WORLD, 1, 1, 4, 4, 2, long_n2);
    test_r2c(MPI_COMM_WORLD, 1, 1, 4, 6, 2, long_n2);

    test_r2c(MPI_COMM_WORLD, 1, 1, 5, 2, 2, long_n2);
    test_r2c(MPI_COMM_WORLD, 1, 1, 5, 4, 2, long_n2);
    test_r2c(MPI_COMM_WORLD, 1, 1, 5, 6, 2, long_n2);

    test_r2c(MPI_COMM_WORLD, 1, 1, 6, 2, 2, long_n2);
    test_r2c(MPI_COMM_WORLD, 1, 1, 6, 4, 2, long_n2);
    test_r2c(MPI_COMM_WORLD, 1, 1, 6, 6, 2, long_n2);

    // Long in n2, transform multiple pencils
    test_r2c(MPI_COMM_WORLD, 2, 2, 3, 2, 2, long_n2);
    test_r2c(MPI_COMM_WORLD, 2, 3, 3, 4, 2, long_n2);
    test_r2c(MPI_COMM_WORLD, 3, 2, 3, 6, 2, long_n2);

    test_r2c(MPI_COMM_WORLD, 2, 2, 4, 2, 2, long_n2);
    test_r2c(MPI_COMM_WORLD, 2, 3, 4, 4, 2, long_n2);
    test_r2c(MPI_COMM_WORLD, 3, 2, 4, 6, 2, long_n2);

    test_r2c(MPI_COMM_WORLD, 2, 2, 5, 2, 2, long_n2);
    test_r2c(MPI_COMM_WORLD, 2, 3, 5, 4, 2, long_n2);
    test_r2c(MPI_COMM_WORLD, 3, 2, 5, 6, 2, long_n2);

    test_r2c(MPI_COMM_WORLD, 2, 2, 6, 2, 2, long_n2);
    test_r2c(MPI_COMM_WORLD, 2, 3, 6, 4, 2, long_n2);
    test_r2c(MPI_COMM_WORLD, 3, 2, 6, 6, 2, long_n2);
}

BOOST_AUTO_TEST_CASE( underling_fft_r2c_simple_n2_transposed_long_n0 )
{
    using suzerain::underling::transposed::long_n0;

    // Long in n2, transform a single pencil
    test_r2c(MPI_COMM_WORLD, 1, 1, 3, 2, 2, long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 1, 3, 4, 2, long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 1, 3, 6, 2, long_n0);

    test_r2c(MPI_COMM_WORLD, 1, 1, 4, 2, 2, long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 1, 4, 4, 2, long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 1, 4, 6, 2, long_n0);

    test_r2c(MPI_COMM_WORLD, 1, 1, 5, 2, 2, long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 1, 5, 4, 2, long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 1, 5, 6, 2, long_n0);

    test_r2c(MPI_COMM_WORLD, 1, 1, 6, 2, 2, long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 1, 6, 4, 2, long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 1, 6, 6, 2, long_n0);

    // Long in n2, transform multiple pencils
    test_r2c(MPI_COMM_WORLD, 2, 2, 3, 2, 2, long_n0);
    test_r2c(MPI_COMM_WORLD, 2, 3, 3, 4, 2, long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 2, 3, 6, 2, long_n0);

    test_r2c(MPI_COMM_WORLD, 2, 2, 4, 2, 2, long_n0);
    test_r2c(MPI_COMM_WORLD, 2, 3, 4, 4, 2, long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 2, 4, 6, 2, long_n0);

    test_r2c(MPI_COMM_WORLD, 2, 2, 5, 2, 2, long_n0);
    test_r2c(MPI_COMM_WORLD, 2, 3, 5, 4, 2, long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 2, 5, 6, 2, long_n0);

    test_r2c(MPI_COMM_WORLD, 2, 2, 6, 2, 2, long_n0);
    test_r2c(MPI_COMM_WORLD, 2, 3, 6, 4, 2, long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 2, 6, 6, 2, long_n0);
}

BOOST_AUTO_TEST_CASE( underling_fft_r2c_simple_n2_transposed_both )
{
    ensureFFTWTensor7PatchInPlace();

    using suzerain::underling::transposed::long_n2;
    using suzerain::underling::transposed::long_n0;

    // Long in n2, transform a single pencil
    test_r2c(MPI_COMM_WORLD, 1, 1, 3, 2, 2, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 1, 3, 4, 2, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 1, 3, 6, 2, long_n2 | long_n0);

    test_r2c(MPI_COMM_WORLD, 1, 1, 4, 2, 2, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 1, 4, 4, 2, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 1, 4, 6, 2, long_n2 | long_n0);

    test_r2c(MPI_COMM_WORLD, 1, 1, 5, 2, 2, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 1, 5, 4, 2, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 1, 5, 6, 2, long_n2 | long_n0);

    test_r2c(MPI_COMM_WORLD, 1, 1, 6, 2, 2, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 1, 6, 4, 2, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 1, 1, 6, 6, 2, long_n2 | long_n0);

    // Long in n2, transform multiple pencils
    test_r2c(MPI_COMM_WORLD, 2, 2, 3, 2, 2, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 2, 3, 3, 4, 2, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 2, 3, 6, 2, long_n2 | long_n0);

    test_r2c(MPI_COMM_WORLD, 2, 2, 4, 2, 2, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 2, 3, 4, 4, 2, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 2, 4, 6, 2, long_n2 | long_n0);

    test_r2c(MPI_COMM_WORLD, 2, 2, 5, 2, 2, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 2, 3, 5, 4, 2, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 2, 5, 6, 2, long_n2 | long_n0);

    test_r2c(MPI_COMM_WORLD, 2, 2, 6, 2, 2, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 2, 3, 6, 4, 2, long_n2 | long_n0);
    test_r2c(MPI_COMM_WORLD, 3, 2, 6, 6, 2, long_n2 | long_n0);
}

BOOST_AUTO_TEST_SUITE_END()
