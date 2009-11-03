#include <suzerain/config.h>
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>

typedef std::complex<double> complex;

// Ensure we can use std::complex<double> as two consecutive doubles
BOOST_AUTO_TEST_CASE( shared_c_array )
{
    // Assumption must hold true for any of this scheme to work
    BOOST_CHECK_EQUAL( sizeof(complex), 2*sizeof(double) );

    const std::size_t N = 6;
    double carray_double[N] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
    complex *carray_complex = reinterpret_cast<complex *>(carray_double);

    for (std::size_t i = 0; i < N / 2; ++i) {
        BOOST_CHECK_EQUAL(carray_double[2*i],   carray_complex[i].real());
        BOOST_CHECK_EQUAL(carray_double[2*i+1], carray_complex[i].imag());
    }

    for (std::size_t i = 0; i < N / 2; ++i) {
        carray_complex[i] += carray_complex[0];
    }

    for (std::size_t i = 0; i < N / 2; ++i) {
        BOOST_CHECK_EQUAL(carray_double[2*i],   carray_complex[i].real());
        BOOST_CHECK_EQUAL(carray_double[2*i+1], carray_complex[i].imag());
    }
}

BOOST_AUTO_TEST_CASE( real_and_imag_return_lvalues )
{
    complex x(1,2);
    x.real() += 1;
    x.imag() += 2;
    BOOST_CHECK_EQUAL(x.real(), 2);
    BOOST_CHECK_EQUAL(x.imag(), 4);
    BOOST_CHECK_EQUAL(x, complex(2,4));
}
