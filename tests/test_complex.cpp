//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Suzerain is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------

#include <suzerain/complex.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <fftw3.h>

#include <suzerain/common.hpp>

#pragma warning(disable:1572)

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

BOOST_AUTO_TEST_CASE( casting_sanity_check_one )
{
    std::complex<double> z(3,4);
    double (*a)[2] = reinterpret_cast<double (*)[2]>(&z);
    BOOST_CHECK_EQUAL((*a)[0], 3.0);
    BOOST_CHECK_EQUAL((*a)[1], 4.0);
}

BOOST_AUTO_TEST_CASE( casting_sanity_check_two )
{
    std::complex<double> z(3,4);
    double a[2];
    std::memcpy(&a, &z, sizeof(std::complex<double>));
    BOOST_CHECK_EQUAL(a[0], 3.0);
    BOOST_CHECK_EQUAL(a[1], 4.0);
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

    BOOST_CHECK(is_complex<const std::complex<double> >::value == true);
    BOOST_CHECK(is_complex<const std::complex<float> >::value == true);
    BOOST_CHECK(is_complex<const std::complex<int> >::value == true);

    BOOST_CHECK(is_complex<std::complex<const double> >::value == true);
    BOOST_CHECK(is_complex<std::complex<const float> >::value == true);
    BOOST_CHECK(is_complex<std::complex<const int> >::value == true);

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

BOOST_AUTO_TEST_CASE( real_type_from_const )
{
    using boost::is_same;
    using suzerain::complex::traits::real;

    BOOST_CHECK((is_same<real<
                const std::complex<float> >::type,float>::value));
    BOOST_CHECK((is_same<real<
                const std::complex<double> >::type,double>::value));
    BOOST_CHECK((is_same<real<
                const std::complex<int> >::type,int>::value));

    BOOST_CHECK((is_same<real<const fftwf_complex>::type,float>::value));
    BOOST_CHECK((is_same<real<const fftw_complex>::type,double>::value));
    BOOST_CHECK((is_same<real<const fftwl_complex>::type,long double>::value));
}

BOOST_AUTO_TEST_CASE( is_complex_float )
{
    using suzerain::complex::traits::is_complex_float;

    BOOST_CHECK(is_complex_float<std::complex<long double> >::value == false);
    BOOST_CHECK(is_complex_float<std::complex<double> >::value == false);
    BOOST_CHECK(is_complex_float<std::complex<float> >::value == true);
    BOOST_CHECK(is_complex_float<std::complex<int> >::value == false);

    BOOST_CHECK(is_complex_float<double>::value == false);
    BOOST_CHECK(is_complex_float<float>::value == false);
    BOOST_CHECK(is_complex_float<int>::value == false);

    BOOST_CHECK(is_complex_float<fftwf_complex>::value == true);
    BOOST_CHECK(is_complex_float<fftw_complex>::value == false);
    BOOST_CHECK(is_complex_float<fftwl_complex>::value == false);
}

BOOST_AUTO_TEST_CASE( is_complex_double )
{
    using suzerain::complex::traits::is_complex_double;

    BOOST_CHECK(is_complex_double<std::complex<long double> >::value == false);
    BOOST_CHECK(is_complex_double<std::complex<double> >::value == true);
    BOOST_CHECK(is_complex_double<std::complex<float> >::value == false);
    BOOST_CHECK(is_complex_double<std::complex<int> >::value == false);

    BOOST_CHECK(is_complex_double<double>::value == false);
    BOOST_CHECK(is_complex_double<float>::value == false);
    BOOST_CHECK(is_complex_double<int>::value == false);

    BOOST_CHECK(is_complex_double<fftwf_complex>::value == false);
    BOOST_CHECK(is_complex_double<fftw_complex>::value == true);
    BOOST_CHECK(is_complex_double<fftwl_complex>::value == false);
}

BOOST_AUTO_TEST_CASE( is_complex_t )
{
    using suzerain::complex::traits::is_complex_t;

    BOOST_CHECK(is_complex_t<std::complex<suzerain::real_t> >::value == true);
    BOOST_CHECK(is_complex_t<std::complex<int> >::value == false);
}

BOOST_AUTO_TEST_CASE( is_complex_long_double )
{
    using suzerain::complex::traits::is_complex_long_double;

    BOOST_CHECK(is_complex_long_double<std::complex<long double> >::value == true);
    BOOST_CHECK(is_complex_long_double<std::complex<double> >::value == false);
    BOOST_CHECK(is_complex_long_double<std::complex<float> >::value == false);
    BOOST_CHECK(is_complex_long_double<std::complex<int> >::value == false);

    BOOST_CHECK(is_complex_long_double<double>::value == false);
    BOOST_CHECK(is_complex_long_double<float>::value == false);
    BOOST_CHECK(is_complex_long_double<int>::value == false);

    BOOST_CHECK(is_complex_long_double<fftwf_complex>::value == false);
    BOOST_CHECK(is_complex_long_double<fftw_complex>::value == false);
    BOOST_CHECK(is_complex_long_double<fftwl_complex>::value == true);
}

// When GCC's -ffast-math is enabled the usual (x != x) test for NaN fails
// See http://gcc.gnu.org/onlinedocs/gcc-4.3.5/gcc/Optimize-Options.html
#if defined __GNUC__ && defined __FAST_MATH__
BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES( NaN, 8 );
#endif

BOOST_AUTO_TEST_CASE( NaN )
{
    using suzerain::complex::NaN;

    {
        std::complex<double> z = NaN<std::complex<double> >();
        BOOST_CHECK_NE(z.real(), z.real());
        BOOST_CHECK_NE(z.imag(), z.imag());
    }

    {
        std::complex<double> z = NaN<double>();
        BOOST_CHECK_NE(z.real(), z.real());
        BOOST_CHECK_NE(z.imag(), z.imag());
    }

    {
        std::complex<float> z = NaN<std::complex<float> >();
        BOOST_CHECK_NE(z.real(), z.real());
        BOOST_CHECK_NE(z.imag(), z.imag());
    }

    {
        std::complex<float> z = NaN<float>();
        BOOST_CHECK_NE(z.real(), z.real());
        BOOST_CHECK_NE(z.imag(), z.imag());
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

BOOST_AUTO_TEST_CASE( accumulate_complex )
{
    using suzerain::complex::accumulate_complex;
    using suzerain::complex::real;
    using suzerain::complex::imag;
    typedef std::complex<double> complex_type;
    typedef double real_type;

    { // CCC
        const complex_type beta(2,3);
        complex_type       dest(5,7);
        const complex_type alpha(11,13);
        const complex_type src(17,19);

        accumulate_complex(beta, dest, alpha, src);

        BOOST_CHECK_EQUAL(- 71, real(dest));
        BOOST_CHECK_EQUAL( 459, imag(dest));
    }

    { // CCR
        const complex_type beta(2,3);
        complex_type       dest(5,7);
        const complex_type alpha(11,13);
        const real_type    src(17);

        accumulate_complex(beta, dest, alpha, src);

        BOOST_CHECK_EQUAL( 176, real(dest));
        BOOST_CHECK_EQUAL( 250, imag(dest));
    }

    { // CRC
        const complex_type beta(2,3);
        complex_type       dest(5,7);
        const real_type    alpha(11);
        const complex_type src(17,19);

        accumulate_complex(beta, dest, alpha, src);

        BOOST_CHECK_EQUAL( 176, real(dest));
        BOOST_CHECK_EQUAL( 238, imag(dest));
    }

    { // CRR
        const complex_type beta(2,3);
        complex_type       dest(5,7);
        const real_type    alpha(11);
        const real_type    src(17);

        accumulate_complex(beta, dest, alpha, src);

        BOOST_CHECK_EQUAL( 176, real(dest));
        BOOST_CHECK_EQUAL(  29, imag(dest));
    }

    { // RCC
        const real_type    beta(2);
        complex_type       dest(5,7);
        const complex_type alpha(11,13);
        const complex_type src(17,19);

        accumulate_complex(beta, dest, alpha, src);

        BOOST_CHECK_EQUAL(- 50, real(dest));
        BOOST_CHECK_EQUAL( 444, imag(dest));
    }

    { // RCR
        const real_type    beta(2);
        complex_type       dest(5,7);
        const complex_type alpha(11,13);
        const real_type    src(17);

        accumulate_complex(beta, dest, alpha, src);

        BOOST_CHECK_EQUAL( 197, real(dest));
        BOOST_CHECK_EQUAL( 235, imag(dest));
    }

    { // RRC
        const real_type    beta(2);
        complex_type       dest(5,7);
        const real_type    alpha(11);
        const complex_type src(17,19);

        accumulate_complex(beta, dest, alpha, src);

        BOOST_CHECK_EQUAL( 197, real(dest));
        BOOST_CHECK_EQUAL( 223, imag(dest));
    }

    { // RRR
        const real_type    beta(2);
        complex_type       dest(5,7);
        const real_type    alpha(11);
        const real_type    src(17);

        accumulate_complex(beta, dest, alpha, src);

        BOOST_CHECK_EQUAL( 197, real(dest));
        BOOST_CHECK_EQUAL(  14, imag(dest));
    }
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
