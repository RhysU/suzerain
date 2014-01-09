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
#include <underling/periodic_function.hpp>

#define BOOST_TEST_MODULE $Id$
#include <fftw3.h>
#include <boost/test/included/unit_test.hpp>
#include "test_tools.hpp"

using underling::periodic_function;

// Not so much a test as a way to check our periodic test function definition.
// Getting it correct is key to all other FFT-related tests

static void test_c2c_forward(const int N, const int max_mode_exclusive)
{
    BOOST_TEST_MESSAGE("Testing N = " << N
                       << " with max_mode_exclusive = " << max_mode_exclusive);

    const double close = std::numeric_limits<double>::epsilon()*100*N*N*N;
    typedef std::complex<double> complex_type;
    const complex_type I(0, 1);
    using boost::scoped_array;
    using boost::shared_ptr;

    scoped_array<complex_type> buf(new complex_type[N]);

    boost::shared_ptr<boost::remove_pointer<fftw_plan>::type> forward(
        fftw_plan_dft_1d(N,
                         (fftw_complex *) buf.get(),
                         (fftw_complex *) buf.get(),
                         FFTW_FORWARD, FFTW_ESTIMATE),
        &fftw_destroy_plan);
    BOOST_REQUIRE(forward);

    periodic_function<double,int> pf(N, max_mode_exclusive);
    for (int i = 0; i < N; ++i) {
        const complex_type::value_type val = pf.physical(i);
        buf[i] = complex_type(val, -val);
    }

    fftw_execute(forward.get());

    for (int i = 0; i < N; ++i)  {
        const complex_type val = pf.wave(i);
        const complex_type expected = ((double) N) * (val - I*val);
        if (std::abs(expected.real()) < close) {
            BOOST_CHECK_SMALL(buf[i].real(), close);
        } else {
            BOOST_CHECK_CLOSE(buf[i].real(), expected.real(), close);
        }
        if (std::abs(expected.imag()) < close) {
            BOOST_CHECK_SMALL(buf[i].imag(), close);
        } else {
            BOOST_CHECK_CLOSE(buf[i].imag(), expected.imag(), close);
        }
    }

    boost::shared_ptr<boost::remove_pointer<fftw_plan>::type> backward(
        fftw_plan_dft_1d(N,
                         (fftw_complex *) buf.get(),
                         (fftw_complex *) buf.get(),
                         FFTW_BACKWARD, FFTW_ESTIMATE),
        &fftw_destroy_plan);
    BOOST_REQUIRE(backward);

    fftw_execute(backward.get());

    for (int i = 0; i < N; ++i) {
        const complex_type::value_type val = ((double) N) * pf.physical(i);
        const complex_type expected(val, -val);
        if (std::abs(expected.real()) < close) {
            BOOST_CHECK_SMALL(buf[i].real(), close);
        } else {
            BOOST_CHECK_CLOSE(buf[i].real(), expected.real(), close);
        }
        if (std::abs(expected.imag()) < close) {
            BOOST_CHECK_SMALL(buf[i].imag(), close);
        } else {
            BOOST_CHECK_CLOSE(buf[i].imag(), expected.imag(), close);
        }
    }
}

BOOST_AUTO_TEST_CASE( check_c2c_forward )
{
    for (int i = 1; i < 17; ++i) {
        for (int j = 0; j < i/2+1; ++j) {
            test_c2c_forward(i, j);
        }
    }
}

static void test_c2c_backward(const int N, const int max_mode_exclusive)
{
    BOOST_TEST_MESSAGE("Testing N = " << N
                       << " with max_mode_exclusive = " << max_mode_exclusive);

    const double close = std::numeric_limits<double>::epsilon()*150*N*N*N;
    typedef std::complex<double> complex_type;
    const complex_type I(0, 1);
    using boost::scoped_array;
    using boost::shared_ptr;

    scoped_array<complex_type> buf(new complex_type[N]);

    boost::shared_ptr<boost::remove_pointer<fftw_plan>::type> backward(
        fftw_plan_dft_1d(N,
                         (fftw_complex *) buf.get(),
                         (fftw_complex *) buf.get(),
                         FFTW_BACKWARD, FFTW_ESTIMATE),
        &fftw_destroy_plan);
    BOOST_REQUIRE(backward);

    periodic_function<double,int> pf(N, max_mode_exclusive);
    for (int i = 0; i < N; ++i) {
        const complex_type val = pf.wave(i);
        buf[i] = val - I*val;
    }

    fftw_execute(backward.get());

    for (int i = 0; i < N; ++i)  {
        const complex_type::value_type val = pf.physical(i);
        const complex_type expected(val, -val);
        if (std::abs(expected.real()) < close) {
            BOOST_CHECK_SMALL(buf[i].real(), close);
        } else {
            BOOST_CHECK_CLOSE(buf[i].real(), expected.real(), close);
        }
        if (std::abs(expected.imag()) < close) {
            BOOST_CHECK_SMALL(buf[i].imag(), close);
        } else {
            BOOST_CHECK_CLOSE(buf[i].imag(), expected.imag(), close);
        }
    }

    boost::shared_ptr<boost::remove_pointer<fftw_plan>::type> forward(
        fftw_plan_dft_1d(N,
                         (fftw_complex *) buf.get(),
                         (fftw_complex *) buf.get(),
                         FFTW_FORWARD, FFTW_ESTIMATE),
        &fftw_destroy_plan);
    BOOST_REQUIRE(forward);

    fftw_execute(forward.get());

    for (int i = 0; i < N; ++i) {
        const complex_type val = pf.wave(i);
        const complex_type expected = ((double) N) * (val - I*val);
        if (std::abs(expected.real()) < close) {
            BOOST_CHECK_SMALL(buf[i].real(), close);
        } else {
            BOOST_CHECK_CLOSE(buf[i].real(), expected.real(), close);
        }
        if (std::abs(expected.imag()) < close) {
            BOOST_CHECK_SMALL(buf[i].imag(), close);
        } else {
            BOOST_CHECK_CLOSE(buf[i].imag(), expected.imag(), close);
        }
    }
}

BOOST_AUTO_TEST_CASE( check_c2c_backward )
{
    for (int i = 1; i < 17; ++i) {
        for (int j = 0; j < i/2+1; ++j) {
            test_c2c_backward(i, j);
        }
    }
}

static void test_r2c_forward(const int N, const int max_mode_exclusive)
{
    BOOST_TEST_MESSAGE("Testing N = " << N
                       << " with max_mode_exclusive = " << max_mode_exclusive);

    const double close = std::numeric_limits<double>::epsilon()*50*N*N*N;
    typedef std::complex<double> complex_type;
    using boost::scoped_array;
    using boost::shared_ptr;

    scoped_array<double> buf(new double[2*(N/2+1)]);
    double       * const rbuf = buf.get();
    complex_type * const cbuf = (complex_type *) buf.get();

    boost::shared_ptr<boost::remove_pointer<fftw_plan>::type> forward(
        fftw_plan_dft_r2c_1d(N, rbuf, (fftw_complex *)cbuf, FFTW_ESTIMATE),
        &fftw_destroy_plan);
    BOOST_REQUIRE(forward);

    periodic_function<double,int> pf(N, max_mode_exclusive);
    for (int i = 0; i < N; ++i) {
        rbuf[i] = pf.physical(i);
    }

    fftw_execute(forward.get());

    for (int i = 0; i < (N/2+1); ++i)  {
        const complex_type expected = ((double) N) * pf.wave(i);
        if (std::abs(expected.real()) < close) {
            BOOST_CHECK_SMALL(cbuf[i].real(), close);
        } else {
            BOOST_CHECK_CLOSE(cbuf[i].real(), expected.real(), close);
        }
        if (std::abs(expected.imag()) < close) {
            BOOST_CHECK_SMALL(cbuf[i].imag(), close);
        } else {
            BOOST_CHECK_CLOSE(cbuf[i].imag(), expected.imag(), close);
        }
    }
}

BOOST_AUTO_TEST_CASE( check_r2c_forward )
{
    for (int i = 1; i < 17; ++i) {
        for (int j = 0; j < i/2+1; ++j) {
            test_r2c_forward(i, j);
        }
    }
}

static void test_c2r_backward(const int N, const int max_mode_exclusive)
{
    BOOST_TEST_MESSAGE("Testing N = " << N
                       << " with max_mode_exclusive = " << max_mode_exclusive);

    const double close = std::numeric_limits<double>::epsilon()*150*N*N*N;
    typedef std::complex<double> complex_type;
    using boost::scoped_array;
    using boost::shared_ptr;

    scoped_array<double> buf(new double[2*(N/2+1)]);
    double       * const rbuf = buf.get();
    complex_type * const cbuf = (complex_type *) buf.get();

    boost::shared_ptr<boost::remove_pointer<fftw_plan>::type> forward(
        fftw_plan_dft_c2r_1d(N, (fftw_complex *)cbuf, rbuf, FFTW_ESTIMATE),
        &fftw_destroy_plan);
    BOOST_REQUIRE(forward);

    periodic_function<double,int> pf(N, max_mode_exclusive);
    for (int i = 0; i < (N/2+1); ++i) {
        cbuf[i] = pf.wave(i);
    }

    fftw_execute(forward.get());

    for (int i = 0; i < N; ++i)  {
        const double expected = pf.physical(i);
        if (std::abs(expected) < close) {
            BOOST_CHECK_SMALL(rbuf[i], close);
        } else {
            BOOST_CHECK_CLOSE(rbuf[i], expected, close);
        }
    }
}

BOOST_AUTO_TEST_CASE( check_c2r_backward )
{
    for (int i = 1; i < 17; ++i) {
        for (int j = 0; j < i/2+1; ++j) {
            test_c2r_backward(i, j);
        }
    }
}
