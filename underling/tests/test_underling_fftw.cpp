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
#include <underling/underling_fftw.hpp>

#include <boost/test/included/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <boost/program_options.hpp>
#include <mpi.h>
#include <fftw3-mpi.h>
#include <underling/error.h>
#include <underling/periodic_function.hpp>
#include "test_tools.hpp"
#include "fixtures.hpp"

using underling::periodic_function;

// For unary function-based test case registration
struct tc
{
    int n0, n1, n2, howmany, long_ni;
    bool in_place;
    unsigned flags;
    unsigned packed;
};

// Pull out the two slow directions, excluding long_ni
static void slow_non_long_directions(
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
static void ensureFFTWTensor7PatchInPlace() {
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

// Forward physical-to-wave followed by wave-to-physical
static void test_c2c_forward(tc t)
{
    BoostFailErrorHandlerFixture fix1;

    // Unpack test case parameters
    MPI_Comm comm         = MPI_COMM_WORLD;
    const int n0          = t.n0;
    const int n1          = t.n1;
    const int n2          = t.n2;
    const int howmany     = t.howmany;
    const int long_ni     = t.long_ni;
    const unsigned flags  = t.flags;
    const unsigned packed = t.packed;
    const bool in_place   = t.in_place;

    if (flags) ensureFFTWTensor7PatchInPlace();

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
                           << (nproc > 1 ? "s" : "")
                           << " when long in " << long_ni
                           << " with flags " << flags
                           << " and packed " << packed);
    }
    if (in_place && long_ni == 2 && packed & underling::fftw::packed::long_n2) {
        if (!procid) BOOST_TEST_MESSAGE("Skipping invalid test");
        return;
    }
    if (in_place && long_ni == 0 && packed & underling::fftw::packed::long_n0) {
        if (!procid) BOOST_TEST_MESSAGE("Skipping invalid test");
        return;
    }
    if (nproc > 1 && (n0 == 1 || n1 == 1 || n2 == 1)) {
        if (!procid) BOOST_TEST_MESSAGE("Skipping degenerate parallel test");
        return;
    }

    const underling_real close
        = std::numeric_limits<underling_real>::epsilon()*150*n0*n1*n2;

    UnderlingFixture f(comm, n0, n1, n2, howmany, flags, in_place);
    underling::fftw::plan forward(underling::fftw::plan::c2c_forward(),
                                  f.problem,
                                  long_ni,
                                  f.in,
                                  f.out,
                                  FFTW_ESTIMATE,
                                  packed);
    BOOST_REQUIRE(forward);
    underling::fftw::plan backward(forward,         // Inverse constructor!
                                   f.out,
                                   f.in,
                                   FFTW_ESTIMATE);
    BOOST_REQUIRE(backward);

    // Stride information consistency check
    BOOST_REQUIRE_EQUAL(forward.local_extents_output(),
                        backward.local_extents_input());
    BOOST_REQUIRE_EQUAL(backward.local_extents_output(),
                        forward.local_extents_input());

    // Ensure FFT output is fastest in the long_ni direction
    BOOST_REQUIRE_EQUAL(4,       forward.local_extents_output().order[0]);
    BOOST_REQUIRE_EQUAL(3,       forward.local_extents_output().order[1]);
    BOOST_REQUIRE_EQUAL(long_ni, forward.local_extents_output().order[2]);

    // Load up sample data
    {
        const underling::fftw::extents e = forward.local_extents_input();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(e.size[long_ni],
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    underling_real * const pencil = &f.in[
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
    if (!f.in_place) f.fill_out_with_NaNs();

    // Transform from physical to wave space
    forward.execute(f.in, f.out);
    if (!f.in_place) f.fill_in_with_NaNs(); // Poison old buffer

    // Check the sample data transformed as expected
    {
        const underling::fftw::extents e = forward.local_extents_output();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(e.size[long_ni],
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    const underling_real * const pencil = &f.out[
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
    backward.execute(f.out, f.in);
    if (!f.in_place) f.fill_out_with_NaNs(); // Poison old buffer

    // Check that we recovered the scaled sample data
    {
        const underling::fftw::extents e = backward.local_extents_output();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(e.size[long_ni],
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    const underling_real * const pencil = &f.in[
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

// Backward wave-to-physical followed by physical-to-wave
static void test_c2c_backward(tc t)
{
    BoostFailErrorHandlerFixture fix1;

    // Unpack test case parameters
    MPI_Comm comm         = MPI_COMM_WORLD;
    const int n0          = t.n0;
    const int n1          = t.n1;
    const int n2          = t.n2;
    const int howmany     = t.howmany;
    const int long_ni     = t.long_ni;
    const unsigned flags  = t.flags;
    const unsigned packed = t.packed;
    const bool in_place   = t.in_place;

    if (flags) ensureFFTWTensor7PatchInPlace();

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
                           << (nproc > 1 ? "s" : "")
                           << " when long in " << long_ni
                           << " with flags " << flags
                           << " and packed " << packed);
    }
    if (in_place && long_ni == 2 && packed & underling::fftw::packed::long_n2) {
        if (!procid) BOOST_TEST_MESSAGE("Skipping invalid test");
        return;
    }
    if (in_place && long_ni == 0 && packed & underling::fftw::packed::long_n0) {
        if (!procid) BOOST_TEST_MESSAGE("Skipping invalid test");
        return;
    }
    if (nproc > 1 && (n0 == 1 || n1 == 1 || n2 == 1)) {
        if (!procid) BOOST_TEST_MESSAGE("Skipping degenerate parallel test");
        return;
    }

    const underling_real close
        = std::numeric_limits<underling_real>::epsilon()*150*n0*n1*n2;

    UnderlingFixture f(comm, n0, n1, n2, howmany, flags, in_place);
    underling::fftw::plan backward(underling::fftw::plan::c2c_backward(),
                                   f.problem,
                                   long_ni,
                                   f.in,
                                   f.out,
                                   FFTW_ESTIMATE,
                                   packed);
    BOOST_REQUIRE(backward);
    underling::fftw::plan forward(backward,         // Inverse constructor!
                                  f.out,
                                  f.in,
                                  FFTW_ESTIMATE);
    BOOST_REQUIRE(forward);

    // Stride information consistency check
    BOOST_REQUIRE_EQUAL(backward.local_extents_output(),
                        forward.local_extents_input());
    BOOST_REQUIRE_EQUAL(forward.local_extents_output(),
                        backward.local_extents_input());

    // Ensure FFT backward is fastest in the long_ni direction
    BOOST_REQUIRE_EQUAL(4,       backward.local_extents_output().order[0]);
    BOOST_REQUIRE_EQUAL(3,       backward.local_extents_output().order[1]);
    BOOST_REQUIRE_EQUAL(long_ni, backward.local_extents_output().order[2]);

    // Load up sample data
    {
        const underling::fftw::extents e = backward.local_extents_input();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(e.size[long_ni],
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    underling_real * const pencil = &f.in[
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
    if (!f.in_place) f.fill_out_with_NaNs();

    // Transform from wave space to physical space
    backward.execute(f.in, f.out);
    if (!f.in_place) f.fill_in_with_NaNs(); // Poison old buffer

    // Check the sample data transformed as expected
    {
        const underling::fftw::extents e = backward.local_extents_output();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(e.size[long_ni],
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    const underling_real * const pencil = &f.out[
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
    forward.execute(f.out, f.in);
    if (!f.in_place) f.fill_out_with_NaNs(); // Poison old buffer

    // Check the sample data transformed as expected
    {
        const underling::fftw::extents e = forward.local_extents_output();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(e.size[long_ni],
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    const underling_real * const pencil = &f.in[
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

// Test wave to physical transformation and inverse transform
static void test_c2r(tc t)
{
    BoostFailErrorHandlerFixture fix1;

    // Unpack test case parameters
    MPI_Comm comm         = MPI_COMM_WORLD;
    const int n0          = t.n0;
    const int n1          = t.n1;
    const int n2          = t.n2;
    const int howmany     = t.howmany;
    const int long_ni     = t.long_ni;
    const unsigned flags  = t.flags;
    const unsigned packed = t.packed;
    const bool in_place   = t.in_place;

    if (flags) ensureFFTWTensor7PatchInPlace();

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
                           << (nproc > 1 ? "s" : "")
                           << " when long in " << long_ni
                           << " with flags " << flags
                           << " and packed " << packed);
    }
    namespace transposed = underling::transposed;
    if (in_place && long_ni == 2 && flags & transposed::long_n2) {
        if (!procid) BOOST_TEST_MESSAGE("Skipping invalid test");
        return;
    }
    if (in_place && long_ni == 0 && flags & transposed::long_n0) {
        if (!procid) BOOST_TEST_MESSAGE("Skipping invalid test");
        return;
    }
    if (in_place && long_ni == 2 && packed & underling::fftw::packed::long_n2) {
        if (!procid) BOOST_TEST_MESSAGE("Skipping invalid test");
        return;
    }
    if (in_place && long_ni == 0 && packed & underling::fftw::packed::long_n0) {
        if (!procid) BOOST_TEST_MESSAGE("Skipping invalid test");
        return;
    }
    if (nproc > 1 && (n0 == 1 || n1 == 1 || n2 == 1)) {
        if (!procid) BOOST_TEST_MESSAGE("Skipping degenerate parallel test");
        return;
    }

    UnderlingFixture f(comm, n0, n1, n2, howmany, flags, in_place);

    underling::fftw::plan backward(underling::fftw::plan::c2r_backward(),
                                   f.problem,
                                   long_ni,
                                   f.in,
                                   f.out,
                                   FFTW_ESTIMATE,
                                   packed);
    BOOST_REQUIRE(backward);
    underling::fftw::plan forward(backward,       // Inverse constructor!
                                  f.out,
                                  f.in,
                                  FFTW_ESTIMATE);
    BOOST_REQUIRE(forward);

    // Stride information consistency check
    BOOST_REQUIRE_EQUAL(backward.local_extents_output(),
                        forward.local_extents_input());
    BOOST_REQUIRE_EQUAL(forward.local_extents_output(),
                        backward.local_extents_input());

    // Ensure FFT output is fastest in the real-valued long_ni direction
    BOOST_REQUIRE_EQUAL(4,       backward.local_extents_output().order[0]);
    BOOST_REQUIRE_EQUAL(3,       backward.local_extents_output().order[1]);
    BOOST_REQUIRE_EQUAL(long_ni, backward.local_extents_output().order[2]);

    const underling::extents extents = f.problem.local_extents(long_ni);
    const double close_enough
        =   std::numeric_limits<double>::epsilon()*150
          * extents.size[long_ni]*extents.size[long_ni]*extents.size[long_ni];

    // Load up sample data
    {
        const underling::fftw::extents e = backward.local_extents_input();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(
                            e.size[long_ni] == 1 ? 1 : 2*(e.size[long_ni]-1),
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    underling_real * const pencil = &f.in[
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
    if (!f.in_place) f.fill_out_with_NaNs();

    // Transform from wave to physical space
    backward.execute(f.in, f.out);
    if (!f.in_place) f.fill_in_with_NaNs(); // Poison old buffer

    // Check data transformed as expected
    {
        const underling::fftw::extents e = backward.local_extents_output();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(e.size[long_ni],
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    underling_real * const base = &f.out[
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
    forward.execute(f.out, f.in);
    if (!f.in_place) f.fill_out_with_NaNs(); // Poison old buffer

    // Check data transformed as expected
    {
        const underling::fftw::extents e = backward.local_extents_input();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(
                            e.size[long_ni] == 1 ? 1 : 2*(e.size[long_ni]-1),
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    const underling_real * const pencil = &f.in[
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

// Test physical to wave transformation and inverse transform
static void test_r2c(tc t)
{
    BoostFailErrorHandlerFixture fix1;

    // Unpack test case parameters
    MPI_Comm comm         = MPI_COMM_WORLD;
    const int n0          = t.n0;
    const int n1          = t.n1;
    const int n2          = t.n2;
    const int howmany     = t.howmany;
    const int long_ni     = t.long_ni;
    const unsigned flags  = t.flags;
    const unsigned packed = t.packed;
    const bool in_place   = t.in_place;

    if (flags) ensureFFTWTensor7PatchInPlace();

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
                           << (nproc > 1 ? "s" : "")
                           << " when long in " << long_ni
                           << " with flags " << flags
                           << " and packed " << packed);
    }
    namespace transposed = underling::transposed;
    if (in_place && long_ni == 2 && flags & transposed::long_n2) {
        if (!procid) BOOST_TEST_MESSAGE("Skipping invalid test");
        return;
    }
    if (in_place && long_ni == 0 && flags & transposed::long_n0) {
        if (!procid) BOOST_TEST_MESSAGE("Skipping invalid test");
        return;
    }
    if (in_place && long_ni == 2 && packed & underling::fftw::packed::long_n2) {
        if (!procid) BOOST_TEST_MESSAGE("Skipping invalid test");
        return;
    }
    if (in_place && long_ni == 0 && packed & underling::fftw::packed::long_n0) {
        if (!procid) BOOST_TEST_MESSAGE("Skipping invalid test");
        return;
    }
    if (nproc > 1 && (n0 == 1 || n1 == 1 || n2 == 1)) {
        if (!procid) BOOST_TEST_MESSAGE("Skipping degenerate parallel test");
        return;
    }

    UnderlingFixture f(comm, n0, n1, n2, howmany, flags, in_place);

    underling::fftw::plan forward(underling::fftw::plan::r2c_forward(),
                                  f.problem,
                                  long_ni,
                                  f.in,
                                  f.out,
                                  FFTW_ESTIMATE,
                                  packed);
    BOOST_REQUIRE(forward);
    underling::fftw::plan backward(forward,        // Inverse constructor!
                                   f.out,
                                   f.in,
                                   FFTW_ESTIMATE);
    BOOST_REQUIRE(backward);

    // Stride information consistency check
    BOOST_REQUIRE_EQUAL(forward.local_extents_output(),
                        backward.local_extents_input());
    BOOST_REQUIRE_EQUAL(backward.local_extents_output(),
                        forward.local_extents_input());

    // Ensure FFT input is fastest in the real-valued long_ni direction
    BOOST_REQUIRE_EQUAL(4,       forward.local_extents_input().order[0]);
    BOOST_REQUIRE_EQUAL(3,       forward.local_extents_input().order[1]);
    BOOST_REQUIRE_EQUAL(long_ni, forward.local_extents_input().order[2]);

    const underling::extents extents = f.problem.local_extents(long_ni);
    const double close_enough
        =   std::numeric_limits<double>::epsilon()*150
          * extents.size[long_ni]*extents.size[long_ni]*extents.size[long_ni];

    // Load up sample data
    {
        const underling::fftw::extents e = forward.local_extents_input();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(e.size[long_ni],
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    underling_real * const base = &f.in[
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
    if (!f.in_place) f.fill_out_with_NaNs();

    // Transform physical to wave space
    forward.execute(f.in, f.out);
    if (!f.in_place) f.fill_in_with_NaNs(); // Poison old buffer

    // Check data transformed as expected
    {
        const underling::fftw::extents e = forward.local_extents_output();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(
                            e.size[long_ni] == 1 ? 1 : 2*(e.size[long_ni]-1),
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    const underling_real * const pencil = &f.out[
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
    backward.execute(f.out, f.in);
    if (!f.in_place) f.fill_out_with_NaNs(); // Poison old buffer

    // Checking data transformed as expected
    {
        const underling::fftw::extents e = forward.local_extents_input();
        int dir_i, dir_j;
        slow_non_long_directions(e.order, long_ni, dir_i, dir_j);
        const int dir_k = e.order[1];

        for (int i = 0; i < e.size[dir_i]; ++i) {
            for (int j = 0; j < e.size[dir_j]; ++j) {
                for (int k = 0; k < e.size[dir_k]; ++k) {

                    const periodic_function<double> pf(e.size[long_ni],
                            -1, M_PI/3.0, 2.0*M_PI, (i+1)*(j+1)*(k+1));

                    underling_real * const base = &f.in[
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

boost::unit_test::test_suite*
init_unit_test_suite( int argc, char* argv[] )
{

    MPI_Init(&argc, &argv);            // Initialize MPI
    atexit((void (*)()) MPI_Finalize); // Register finalize MPI
    underling_init(&argc, &argv, 0);   // Initialize underling prereqs
    atexit(&underling_cleanup);        // Register finalize underling prereqs

    boost::unit_test::framework::master_test_suite().p_name.value = __FILE__;

    // Size of global extents
    const int extents[][3] = {  { 2, 3, 5 }
                               ,{ 7, 5, 3 }
                               ,{ 8, 6, 4 }
                               ,{ 4, 6, 8 }
                               ,{ 6, 6, 6 }
                               ,{ 8, 1, 1 }  // Reproduce bug #2059
                               ,{ 1, 8, 1 }  // Reproduce bug #2059
                               ,{ 1, 1, 8 }  // Reproduce bug #2059
    };

    // Number of real-valued scalars to transpose
    const int howmanys[] = { 2, 4, 6 };

    // Long in which direction for the tests?
    const int long_nis[] = { 0, 1, 2 };

    // In-place vs out-of-place
    const bool places[] = { true, false };

    // Transposed in/out-like flags
    using underling::transposed::long_n2;
    using underling::transposed::long_n0;
    const unsigned flags[] = { 0, long_n2, long_n0, long_n2 | long_n0 };

    // Packed storage flags
    const unsigned packed[]
        = { underling::fftw::packed::none, underling::fftw::packed::all };

    // Create an outer product of all the cases we want to run
    const size_t ncases = sizeof(extents)/sizeof(extents[0])
                        * sizeof(howmanys)/sizeof(howmanys[0])
                        * sizeof(long_nis)/sizeof(long_nis[0])
                        * sizeof(places)/sizeof(places[0])
                        * sizeof(flags)/sizeof(flags[0])
                        * sizeof(packed)/sizeof(packed[0]);

    tc * const cases = new tc[ncases];
    tc *       c     = cases;

    for (size_t e = 0; e < sizeof(extents)/sizeof(extents[0]); ++e)
    for (size_t h = 0; h < sizeof(howmanys)/sizeof(howmanys[0]); ++h)
    for (size_t l = 0; l < sizeof(long_nis)/sizeof(long_nis[0]); ++l)
    for (size_t p = 0; p < sizeof(places)/sizeof(places[0]); ++p)
    for (size_t f = 0; f < sizeof(flags)/sizeof(flags[0]); ++f)
    for (size_t k = 0; k < sizeof(packed)/sizeof(packed[0]); ++k)
    {
        c->n0       = extents[e][0];
        c->n1       = extents[e][1];
        c->n2       = extents[e][2];
        c->howmany  = howmanys[h];
        c->long_ni  = long_nis[l];
        c->in_place = places[p];
        c->flags    = flags[f];
        c->packed   = packed[k];
        ++c;
    }

    boost::unit_test::test_suite* ts1 = BOOST_TEST_SUITE( "c2c_forward" );
    ts1->add(BOOST_PARAM_TEST_CASE( &test_c2c_forward, cases, cases + ncases),
             /* timeout in seconds */ 30 );
    boost::unit_test::framework::master_test_suite().add( ts1 );

    boost::unit_test::test_suite* ts2 = BOOST_TEST_SUITE( "c2c_backward" );
    ts2->add(BOOST_PARAM_TEST_CASE( &test_c2c_backward, cases, cases + ncases),
             /* timeout in seconds */ 30 );
    boost::unit_test::framework::master_test_suite().add( ts2 );

    boost::unit_test::test_suite* ts3 = BOOST_TEST_SUITE( "c2r" );
    ts3->add(BOOST_PARAM_TEST_CASE( &test_c2r, cases, cases + ncases),
             /* timeout in seconds */ 30 );
    boost::unit_test::framework::master_test_suite().add( ts3 );

    boost::unit_test::test_suite* ts4 = BOOST_TEST_SUITE( "r2c" );
    ts4->add(BOOST_PARAM_TEST_CASE( &test_r2c, cases, cases + ncases),
             /* timeout in seconds */ 30 );
    boost::unit_test::framework::master_test_suite().add( ts4 );

    delete[] cases;

    return 0;
}

static void test_extents_consistency(const bool in_place,
                                     const unsigned packed)
{
    UnderlingFixture f(MPI_COMM_WORLD, 2, 3, 5, 6, /*flags*/0, in_place);

    underling::fftw::plan backward(underling::fftw::plan::c2c_forward(),
                                   f.problem,
                                   0,
                                   f.in,
                                   f.out,
                                   FFTW_ESTIMATE,
                                   packed);

    const int N = 5;

    // Check input information
    {
        const underling::fftw::extents input = backward.local_extents_input();

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
        const underling::fftw::extents output = backward.local_extents_output();

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

BOOST_AUTO_TEST_CASE( extents_consistency )
{
    test_extents_consistency(true,  0);
    test_extents_consistency(false, underling::fftw::packed::none);
    test_extents_consistency(false, underling::fftw::packed::all);
}
