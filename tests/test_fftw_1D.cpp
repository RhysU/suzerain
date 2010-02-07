#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>
#include <fftw3.h>
#include "test_tools.hpp"

// Not so much a test as a way to check our periodic test function definition.

void test_c2c_forward(const int N, const int max_mode_exclusive)
{
    const double close = std::numeric_limits<double>::epsilon()*10*N*N*N;
    typedef std::complex<double> complex_type;
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
        buf[i] = pf.physical(i);
    }

    fftw_execute(forward.get());

    for (int i = 0; i < N; ++i)  {
        const complex_type expected = ((double) N) * pf.wave(i);
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
    test_c2c_forward(4, 1);
    test_c2c_forward(4, 2);
    test_c2c_forward(4, 3);
    test_c2c_forward(5, 1);
    test_c2c_forward(5, 2);
    test_c2c_forward(5, 3);
    test_c2c_forward(6, 1);
    test_c2c_forward(6, 2);
    test_c2c_forward(6, 3);
    test_c2c_forward(6, 4);
}

void test_c2c_backward(const int N, const int max_mode_exclusive)
{
    const double close = std::numeric_limits<double>::epsilon()*10*N*N*N;
    typedef std::complex<double> complex_type;
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
        buf[i] = pf.wave(i);
    }

    fftw_execute(backward.get());

    for (int i = 0; i < N; ++i)  {
        const complex_type expected = pf.physical(i);
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
    test_c2c_backward(4, 1);
    test_c2c_backward(4, 2);
    test_c2c_backward(4, 3);
    test_c2c_backward(5, 1);
    test_c2c_backward(5, 2);
    test_c2c_backward(5, 3);
    test_c2c_backward(6, 1);
    test_c2c_backward(6, 2);
    test_c2c_backward(6, 3);
    test_c2c_backward(6, 4);
}

void test_r2c_forward(const int N, const int max_mode_exclusive)
{
    const double close = std::numeric_limits<double>::epsilon()*10*N*N*N;
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
    test_r2c_forward(4, 1);
    test_r2c_forward(4, 2);
    test_r2c_forward(4, 3);
    test_r2c_forward(5, 1);
    test_r2c_forward(5, 2);
    test_r2c_forward(5, 3);
    test_r2c_forward(6, 1);
    test_r2c_forward(6, 2);
    test_r2c_forward(6, 3);
    test_r2c_forward(6, 4);
}

void test_c2r_backward(const int N, const int max_mode_exclusive)
{
    const double close = std::numeric_limits<double>::epsilon()*10*N*N*N;
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
    test_c2r_backward(4, 1);
    test_c2r_backward(4, 2);
    test_c2r_backward(4, 3);
    test_c2r_backward(5, 1);
    test_c2r_backward(5, 2);
    test_c2r_backward(5, 3);
    test_c2r_backward(6, 1);
    test_c2r_backward(6, 2);
    test_c2r_backward(6, 3);
    test_c2r_backward(6, 4);
}
