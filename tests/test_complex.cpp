#include <suzerain/config.h>
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <fftw3.h>
#include <suzerain/complex.hpp>
#include <boost/test/included/unit_test.hpp>

// Ensure we can use std::complex<double> as two consecutive doubles
BOOST_AUTO_TEST_CASE( shared_c_array )
{
    typedef std::complex<double> complex;

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
    // Officially, std::complex do not return lvalues,
    // but many implementations do so.
    typedef std::complex<double> complex;

    complex x(1,2);
    x.real() += 1;
    x.imag() += 2;
    BOOST_CHECK_EQUAL(x.real(), 2);
    BOOST_CHECK_EQUAL(x.imag(), 4);
    BOOST_CHECK_EQUAL(x, complex(2,4));
}

BOOST_AUTO_TEST_CASE( is_complex )
{
    using suzerain::complex::traits::is_complex;

    BOOST_CHECK(is_complex<std::complex<double> >::value == true);
    BOOST_CHECK(is_complex<std::complex<float> >::value == true);
    BOOST_CHECK(is_complex<std::complex<int> >::value == true);

    BOOST_CHECK(is_complex<double>::value == false);
    BOOST_CHECK(is_complex<float>::value == false);
    BOOST_CHECK(is_complex<int>::value == false);

    BOOST_CHECK(is_complex<fftwf_complex>::value == true);
    BOOST_CHECK(is_complex<fftw_complex>::value == true);
    BOOST_CHECK(is_complex<fftwl_complex>::value == true);
}

BOOST_AUTO_TEST_CASE( real_type )
{
    using boost::is_same;
    using suzerain::complex::traits::real;

    BOOST_CHECK((is_same<real<std::complex<float> >::type,float>::value));
    BOOST_CHECK((is_same<real<std::complex<double> >::type,double>::value));
    BOOST_CHECK((is_same<real<std::complex<int> >::type,int>::value));

    BOOST_CHECK((is_same<real<fftwf_complex>::type,float>::value));
    BOOST_CHECK((is_same<real<fftw_complex>::type,double>::value));
    BOOST_CHECK((is_same<real<fftwl_complex>::type,long double>::value));
}

BOOST_AUTO_TEST_CASE( NaN )
{
    using suzerain::complex::NaN;

    {
        std::complex<double> z = NaN<std::complex<double> >();
        BOOST_CHECK(z.real() != z.real());
        BOOST_CHECK(z.imag() != z.imag());
    }

    {
        std::complex<double> z = NaN<double>();
        BOOST_CHECK(z.real() != z.real());
        BOOST_CHECK(z.imag() != z.imag());
    }

    {
        std::complex<float> z = NaN<std::complex<float> >();
        BOOST_CHECK(z.real() != z.real());
        BOOST_CHECK(z.imag() != z.imag());
    }

    {
        std::complex<float> z = NaN<float>();
        BOOST_CHECK(z.real() != z.real());
        BOOST_CHECK(z.imag() != z.imag());
    }
}

BOOST_AUTO_TEST_CASE( assign_complex_with_complex_source )
{
    using suzerain::complex::assign_complex;

    fftw_complex a, b;
    typedef BOOST_TYPEOF(a[0]) fftw_real;
    std::complex<fftw_real> c;

    // fftw_complex from std::complex
    c.real() = 1.0;
    c.imag() = 2.0;
    assign_complex(b, c);
    BOOST_CHECK_EQUAL(b[0], 1.0);
    BOOST_CHECK_EQUAL(b[1], 2.0);

    // fftw_complex from fftw_complex
    assign_complex(a, b);
    BOOST_CHECK_EQUAL(a[0], 1.0);
    BOOST_CHECK_EQUAL(a[1], 2.0);

    // fftw_complex from components
    assign_complex(b, 3.0, 4.0);
    BOOST_CHECK_EQUAL(b[0], 3.0);
    BOOST_CHECK_EQUAL(b[1], 4.0);

    // std::complex from fftw_complex
    assign_complex(c, b);
    BOOST_CHECK_EQUAL(c.real(), 3.0);
    BOOST_CHECK_EQUAL(c.imag(), 4.0);

    // std::complex from components
    assign_complex(c, 3.0, 4.0);
    BOOST_CHECK_EQUAL(c.real(), 3.0);
    BOOST_CHECK_EQUAL(c.imag(), 4.0);
}

BOOST_AUTO_TEST_CASE( assign_complex_with_arithmetic_source )
{
    using suzerain::complex::assign_complex;

    fftw_complex a;
    a[0] = 123.0;
    a[1] = 456.0;
    assign_complex(a,789.0);
    BOOST_CHECK_EQUAL(a[0], 789.0);
    BOOST_CHECK_EQUAL(a[1], 0.0);

    std::complex<double> b;
    b.real() = 123.0;
    b.imag() = 456.0;
    assign_complex(b, 789.0);
    BOOST_CHECK_EQUAL(b.real(), 789.0);
    BOOST_CHECK_EQUAL(b.imag(), 0.0);
}

BOOST_AUTO_TEST_CASE( assign_complex_scaled )
{
    using suzerain::complex::assign_complex_scaled;

    fftw_complex a, b;
    typedef BOOST_TYPEOF(a[0]) fftw_real;
    std::complex<fftw_real> c, d;

    // fftw_complex from std::complex
    c.real() = 1.0;
    c.imag() = 2.0;
    assign_complex_scaled(b, c, 2.0);
    BOOST_CHECK_EQUAL(b[0], 2.0);
    BOOST_CHECK_EQUAL(b[1], 4.0);

    // fftw_complex from fftw_complex
    assign_complex_scaled(a, b, 2.0);
    BOOST_CHECK_EQUAL(a[0], 4.0);
    BOOST_CHECK_EQUAL(a[1], 8.0);

    // std::complex from fftw_complex
    b[0] = 6.0;
    b[1] = 8.0;
    assign_complex_scaled(c, b, 2.0);
    BOOST_CHECK_EQUAL(c.real(), 12.0);
    BOOST_CHECK_EQUAL(c.imag(), 16.0);

    // std::complex from std::complex
    assign_complex_scaled(d, c, 2.0);
    BOOST_CHECK_EQUAL(d.real(), 24.0);
    BOOST_CHECK_EQUAL(d.imag(), 32.0);
}

BOOST_AUTO_TEST_CASE( assign_complex_scaled_ipower )
{
    using suzerain::complex::assign_complex_scaled_ipower;

    fftw_complex a, b;
    typedef BOOST_TYPEOF(a[0]) fftw_real;
    std::complex<fftw_real> c;

    // fftw_complex from std::complex
    c.real() = 1.0;
    c.imag() = 2.0;
    assign_complex_scaled_ipower(b, c, 2.0, -4);
    BOOST_CHECK_EQUAL(b[0],  2.0);
    BOOST_CHECK_EQUAL(b[1],  4.0);
    assign_complex_scaled_ipower(b, c, 2.0, -3);
    BOOST_CHECK_EQUAL(b[0], -4.0);
    BOOST_CHECK_EQUAL(b[1],  2.0);
    assign_complex_scaled_ipower(b, c, 2.0, -2);
    BOOST_CHECK_EQUAL(b[0], -2.0);
    BOOST_CHECK_EQUAL(b[1], -4.0);
    assign_complex_scaled_ipower(b, c, 2.0, -1);
    BOOST_CHECK_EQUAL(b[0],  4.0);
    BOOST_CHECK_EQUAL(b[1], -2.0);
    assign_complex_scaled_ipower(b, c, 2.0, 0);
    BOOST_CHECK_EQUAL(b[0],  2.0);
    BOOST_CHECK_EQUAL(b[1],  4.0);
    assign_complex_scaled_ipower(b, c, 2.0, 1);
    BOOST_CHECK_EQUAL(b[0], -4.0);
    BOOST_CHECK_EQUAL(b[1],  2.0);
    assign_complex_scaled_ipower(b, c, 2.0, 2);
    BOOST_CHECK_EQUAL(b[0], -2.0);
    BOOST_CHECK_EQUAL(b[1], -4.0);
    assign_complex_scaled_ipower(b, c, 2.0, 3);
    BOOST_CHECK_EQUAL(b[0],  4.0);
    BOOST_CHECK_EQUAL(b[1], -2.0);
    assign_complex_scaled_ipower(b, c, 2.0, 4);
    BOOST_CHECK_EQUAL(b[0],  2.0);
    BOOST_CHECK_EQUAL(b[1],  4.0);

    // fftw_complex from fftw_complex
    b[0] = 1.0;
    b[1] = 2.0;
    assign_complex_scaled_ipower(a, b, 2.0, -4);
    BOOST_CHECK_EQUAL(a[0],  2.0);
    BOOST_CHECK_EQUAL(a[1],  4.0);
    assign_complex_scaled_ipower(a, b, 2.0, -3);
    BOOST_CHECK_EQUAL(a[0], -4.0);
    BOOST_CHECK_EQUAL(a[1],  2.0);
    assign_complex_scaled_ipower(a, b, 2.0, -2);
    BOOST_CHECK_EQUAL(a[0], -2.0);
    BOOST_CHECK_EQUAL(a[1], -4.0);
    assign_complex_scaled_ipower(a, b, 2.0, -1);
    BOOST_CHECK_EQUAL(a[0],  4.0);
    BOOST_CHECK_EQUAL(a[1], -2.0);
    assign_complex_scaled_ipower(a, b, 2.0, 0);
    BOOST_CHECK_EQUAL(a[0],  2.0);
    BOOST_CHECK_EQUAL(a[1],  4.0);
    assign_complex_scaled_ipower(a, b, 2.0, 1);
    BOOST_CHECK_EQUAL(a[0], -4.0);
    BOOST_CHECK_EQUAL(a[1],  2.0);
    assign_complex_scaled_ipower(a, b, 2.0, 2);
    BOOST_CHECK_EQUAL(a[0], -2.0);
    BOOST_CHECK_EQUAL(a[1], -4.0);
    assign_complex_scaled_ipower(a, b, 2.0, 3);
    BOOST_CHECK_EQUAL(a[0],  4.0);
    BOOST_CHECK_EQUAL(a[1], -2.0);
    assign_complex_scaled_ipower(a, b, 2.0, 4);
    BOOST_CHECK_EQUAL(a[0],  2.0);
    BOOST_CHECK_EQUAL(a[1],  4.0);

    // std::complex from fftw_complex
    b[0] = 1.0;
    b[1] = 2.0;
    assign_complex_scaled_ipower(c, b, 2.0, -4);
    BOOST_CHECK_EQUAL(c.real(),  2.0);
    BOOST_CHECK_EQUAL(c.imag(),  4.0);
    assign_complex_scaled_ipower(c, b, 2.0, -3);
    BOOST_CHECK_EQUAL(c.real(), -4.0);
    BOOST_CHECK_EQUAL(c.imag(),  2.0);
    assign_complex_scaled_ipower(c, b, 2.0, -2);
    BOOST_CHECK_EQUAL(c.real(), -2.0);
    BOOST_CHECK_EQUAL(c.imag(), -4.0);
    assign_complex_scaled_ipower(c, b, 2.0, -1);
    BOOST_CHECK_EQUAL(c.real(),  4.0);
    BOOST_CHECK_EQUAL(c.imag(), -2.0);
    assign_complex_scaled_ipower(c, b, 2.0, 0);
    BOOST_CHECK_EQUAL(c.real(),  2.0);
    BOOST_CHECK_EQUAL(c.imag(),  4.0);
    assign_complex_scaled_ipower(c, b, 2.0, 1);
    BOOST_CHECK_EQUAL(c.real(), -4.0);
    BOOST_CHECK_EQUAL(c.imag(),  2.0);
    assign_complex_scaled_ipower(c, b, 2.0, 2);
    BOOST_CHECK_EQUAL(c.real(), -2.0);
    BOOST_CHECK_EQUAL(c.imag(), -4.0);
    assign_complex_scaled_ipower(c, b, 2.0, 3);
    BOOST_CHECK_EQUAL(c.real(),  4.0);
    BOOST_CHECK_EQUAL(c.imag(), -2.0);
    assign_complex_scaled_ipower(c, b, 2.0, 4);
    BOOST_CHECK_EQUAL(c.real(),  2.0);
    BOOST_CHECK_EQUAL(c.imag(),  4.0);
}

BOOST_AUTO_TEST_CASE( assign_components )
{
    using suzerain::complex::assign_components;

    fftw_complex a;
    typedef BOOST_TYPEOF(a[0]) fftw_real;
    std::complex<fftw_real> c;
    fftw_real s1, s2;

    // from std::complex
    c.real() = 1.0;
    c.imag() = 2.0;
    assign_components(s1, s2, c);
    BOOST_CHECK_EQUAL(s1, 1.0);
    BOOST_CHECK_EQUAL(s2, 2.0);

    // from fftw_complex
    a[0] = 3.0;
    a[1] = 4.0;
    assign_components(s1, s2, a);
    BOOST_CHECK_EQUAL(s1, 3.0);
    BOOST_CHECK_EQUAL(s2, 4.0);
}
