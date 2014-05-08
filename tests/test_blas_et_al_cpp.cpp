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

#include <suzerain/blas_et_al/blas_et_al.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>
#include <suzerain/complex.hpp>

#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

// Real-valued types to test
typedef boost::mpl::list<
    double
   ,float
> real_types;

// Complex-valued types to test
typedef boost::mpl::list<
    std::complex<double>
   ,std::complex<float>
   ,double[2]
   ,float[2]
> complex_types;

// Introduce some shorthand
namespace blas = suzerain::blas;
using suzerain::array;
using suzerain::complex::assign_complex;
using suzerain::complex::real;
using suzerain::complex::imag;


// Ensure our traits classes recognize what they should.  Lots and lots of
// error emitted if these do not compile and pass cleanly.
BOOST_AUTO_TEST_SUITE( traits_sanity )

namespace traits = suzerain::complex::traits;

BOOST_AUTO_TEST_CASE_TEMPLATE( trait_is_complex, T, complex_types )
{
   BOOST_CHECK(traits::is_complex<T>::value);
}

BOOST_AUTO_TEST_CASE( traits_is_complex_real )
{
    typedef double d_array[2];
    BOOST_CHECK(traits::is_complex_double<d_array>::value);
    BOOST_CHECK(traits::is_complex_double<double[2]>::value);

    typedef float f_array[2];
    BOOST_CHECK(traits::is_complex_float<f_array>::value);
    BOOST_CHECK(traits::is_complex_float<float[2]>::value);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( swap )

BOOST_AUTO_TEST_CASE_TEMPLATE( real_valued, T, real_types )
{
    array<T,3> x = {{ 1, 2, 3 }};
    array<T,3> y = {{ 4, 5, 6 }};
    blas::swap(x.size(), x.c_array(), 1, y.c_array(), 1);

    const array<T,3> x_expected = {{ 4, 5, 6 }};
    const array<T,3> y_expected = {{ 1, 2, 3 }};
    BOOST_CHECK_EQUAL_COLLECTIONS(
            x.begin(), x.end(), x_expected.begin(), x_expected.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(
            y.begin(), y.end(), y_expected.begin(), y_expected.end());
}

BOOST_AUTO_TEST_CASE_TEMPLATE( complex_valued, T, complex_types )
{
    const long incx = 2;
    array<T,6> x;
    assign_complex(x[0],   1, -  1);
    assign_complex(x[1], 555, -555);
    assign_complex(x[2],   2, -  2);
    assign_complex(x[3], 555, -555);
    assign_complex(x[4],   3, -  3);
    assign_complex(x[5], 555, -555);

    const long incy = 1;
    array<T,3> y;
    assign_complex(y[0],   4, -  4);
    assign_complex(y[1],   5, -  5);
    assign_complex(y[2],   6, -  6);

    BOOST_REQUIRE_EQUAL(x.size()/incx, y.size()/incy);
    blas::swap(x.size()/incx, x.c_array(), incx, y.c_array(), incy);

    BOOST_CHECK_EQUAL(real(x[0]),  4);
    BOOST_CHECK_EQUAL(imag(x[0]), -4);
    BOOST_CHECK_EQUAL(real(x[1]),  555);
    BOOST_CHECK_EQUAL(imag(x[1]), -555);
    BOOST_CHECK_EQUAL(real(x[2]),  5);
    BOOST_CHECK_EQUAL(imag(x[2]), -5);
    BOOST_CHECK_EQUAL(real(x[3]),  555);
    BOOST_CHECK_EQUAL(imag(x[3]), -555);
    BOOST_CHECK_EQUAL(real(x[4]),  6);
    BOOST_CHECK_EQUAL(imag(x[4]), -6);
    BOOST_CHECK_EQUAL(real(x[5]),  555);
    BOOST_CHECK_EQUAL(imag(x[5]), -555);

    BOOST_CHECK_EQUAL(real(y[0]),  1);
    BOOST_CHECK_EQUAL(imag(y[0]), -1);
    BOOST_CHECK_EQUAL(real(y[1]),  2);
    BOOST_CHECK_EQUAL(imag(y[1]), -2);
    BOOST_CHECK_EQUAL(real(y[2]),  3);
    BOOST_CHECK_EQUAL(imag(y[2]), -3);
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE( scal )

BOOST_AUTO_TEST_CASE_TEMPLATE( real_valued, T, real_types )
{
    { // Contiguous
        array<T,3> x = {{ 1, 2, 3 }};
        blas::scal(x.size(), 2, x.c_array(), 1);

        const array<T,3> x_expected = {{ 2, 4, 6 }};
        BOOST_CHECK_EQUAL_COLLECTIONS(
                x.begin(), x.end(), x_expected.begin(), x_expected.end());
    }

    { // Strided
        array<T,6> x = {{ 1, -1, 2, -1, 3, -1 }};
        const long inc = 2;
        blas::scal(x.size()/inc, 2, x.c_array(), inc);

        const array<T,6> x_expected = {{ 2, -1, 4, -1, 6, -1 }};
        BOOST_CHECK_EQUAL_COLLECTIONS(
                x.begin(), x.end(), x_expected.begin(), x_expected.end());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( complex_valued_complex, T, complex_types )
{
    // Complex scaling factor
    T alpha;
    assign_complex(alpha, 0.25, -.50);

    const long inc = 2;
    array<T,6> x;
    assign_complex(x[0],   1, -  1);
    assign_complex(x[1], 555, -555);
    assign_complex(x[2],   2, -  2);
    assign_complex(x[3], 555, -555);
    assign_complex(x[4],   3, -  3);
    assign_complex(x[5], 555, -555);

    blas::scal(x.size()/inc, alpha, x.c_array(), inc);
    BOOST_CHECK_EQUAL(real(x[0]), (1.0* 0.25)-(-1.0*-0.50));
    BOOST_CHECK_EQUAL(imag(x[0]), (1.0*-0.50)+(-1.0* 0.25));
    BOOST_CHECK_EQUAL(real(x[1]),  555);
    BOOST_CHECK_EQUAL(imag(x[1]), -555);
    BOOST_CHECK_EQUAL(real(x[2]), (2.0* 0.25)-(-2.0*-0.50));
    BOOST_CHECK_EQUAL(imag(x[2]), (2.0*-0.50)+(-2.0* 0.25));
    BOOST_CHECK_EQUAL(real(x[3]),  555);
    BOOST_CHECK_EQUAL(imag(x[3]), -555);
    BOOST_CHECK_EQUAL(real(x[4]), (3.0* 0.25)-(-3.0*-0.50));
    BOOST_CHECK_EQUAL(imag(x[4]), (3.0*-0.50)+(-3.0* 0.25));
    BOOST_CHECK_EQUAL(real(x[5]),  555);
    BOOST_CHECK_EQUAL(imag(x[5]), -555);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( complex_valued_real, T, complex_types )
{
    // Real scaling factor
    typename suzerain::complex::traits::real<T>::type alpha = 0.25;

    const long inc = 2;
    array<T,6> x;
    assign_complex(x[0],   1, -  1);
    assign_complex(x[1], 555, -555);
    assign_complex(x[2],   2, -  2);
    assign_complex(x[3], 555, -555);
    assign_complex(x[4],   3, -  3);
    assign_complex(x[5], 555, -555);

    blas::scal(x.size()/inc, alpha, x.c_array(), inc);
    BOOST_CHECK_EQUAL(real(x[0]), ( 1.0*0.25));
    BOOST_CHECK_EQUAL(imag(x[0]), (-1.0*0.25));
    BOOST_CHECK_EQUAL(real(x[1]),  555);
    BOOST_CHECK_EQUAL(imag(x[1]), -555);
    BOOST_CHECK_EQUAL(real(x[2]), ( 2.0*0.25));
    BOOST_CHECK_EQUAL(imag(x[2]), (-2.0*0.25));
    BOOST_CHECK_EQUAL(real(x[3]),  555);
    BOOST_CHECK_EQUAL(imag(x[3]), -555);
    BOOST_CHECK_EQUAL(real(x[4]), ( 3.0*0.25));
    BOOST_CHECK_EQUAL(imag(x[4]), (-3.0*0.25));
    BOOST_CHECK_EQUAL(real(x[5]),  555);
    BOOST_CHECK_EQUAL(imag(x[5]), -555);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( copy )

BOOST_AUTO_TEST_CASE_TEMPLATE( real_valued, T, real_types )
{
    const array<T,6> x = {{   1, -555,   2, -555,   3, -555 }};
    array<T,3> y       = {{ 555,  555, 555 }};
    const long incx = 2, incy = 1;

    BOOST_REQUIRE_EQUAL(x.size()/incx, y.size()/incy);
    blas::copy(x.size()/incx, x.data(), incx, y.c_array(), incy);

    const array<T,3> expected = {{ 1, 2, 3 }};
    BOOST_CHECK_EQUAL_COLLECTIONS(
            y.begin(), y.end(), expected.begin(), expected.end());
}

BOOST_AUTO_TEST_CASE_TEMPLATE( complex_valued, T, complex_types )
{
    const long incx = 2;
    array<T,6> x;
    assign_complex(x[0],   1, -  1);
    assign_complex(x[1], 555, -555);
    assign_complex(x[2],   2, -  2);
    assign_complex(x[3], 555, -555);
    assign_complex(x[4],   3, -  3);
    assign_complex(x[5], 555, -555);

    const long incy = 1;
    array<T,3> y;
    assign_complex(y[0], 777, -777);
    assign_complex(y[1], 777, -777);
    assign_complex(y[2], 777, -777);

    BOOST_REQUIRE_EQUAL(x.size()/incx, y.size()/incy);
    blas::copy(x.size()/incx, x.data(), incx, y.c_array(), incy);

    BOOST_CHECK_EQUAL(real(y[0]),  1);
    BOOST_CHECK_EQUAL(imag(y[0]), -1);
    BOOST_CHECK_EQUAL(real(y[1]),  2);
    BOOST_CHECK_EQUAL(imag(y[1]), -2);
    BOOST_CHECK_EQUAL(real(y[2]),  3);
    BOOST_CHECK_EQUAL(imag(y[2]), -3);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( iamax_and_iamin )

BOOST_AUTO_TEST_CASE_TEMPLATE( real_valued, T, real_types )
{
    const array<T,6> x = {{  -1, -555, 3, 0, 2, -555 }};
    const long incx = 2;
    const size_t N = x.size() / incx;

    BOOST_CHECK_EQUAL( 1, blas::iamax(N, x.data(), incx));
    BOOST_CHECK_EQUAL( 0, blas::iamin(N, x.data(), incx));

    BOOST_CHECK_EQUAL(-1, blas::iamax(0, x.data(), incx));
    BOOST_CHECK_EQUAL(-1, blas::iamin(0, x.data(), incx));
}

BOOST_AUTO_TEST_CASE_TEMPLATE( complex_valued, T, complex_types )
{
    array<T,6> x;
    assign_complex(x[0],   1, -  1);
    assign_complex(x[1], 555, -555);
    assign_complex(x[2],   3, -  3);
    assign_complex(x[3],   0,    0);
    assign_complex(x[4],   2, -  2);
    assign_complex(x[5], 555, -555);
    const long incx = 2;
    const size_t N = x.size() / incx;

    BOOST_CHECK_EQUAL( 1, blas::iamax(N, x.data(), incx));
    BOOST_CHECK_EQUAL( 0, blas::iamin(N, x.data(), incx));

    BOOST_CHECK_EQUAL(-1, blas::iamax(0, x.data(), incx));
    BOOST_CHECK_EQUAL(-1, blas::iamin(0, x.data(), incx));
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( axpy )

BOOST_AUTO_TEST_CASE_TEMPLATE( real_valued, T, real_types )
{
    const array<T,6> x = {{ 1, -1, 2, -1, 3, -1 }};
    const long incx = 2;
    array<T,3> y = {{ 4, 5, 6 }};
    const int incy = 1;

    BOOST_REQUIRE_EQUAL(x.size()/incx, y.size()/incy);
    blas::axpy(x.size()/incx, 3, x.data(), incx, y.c_array(), incy);

    const array<T,3> expected = {{ 1*3+4, 2*3+5, 3*3+6 }};
    BOOST_CHECK_EQUAL_COLLECTIONS(
            y.begin(), y.end(), expected.begin(), expected.end());
}

BOOST_AUTO_TEST_CASE_TEMPLATE( complex_valued_complex, T, complex_types )
{
    // Complex scaling factor
    T alpha;
    assign_complex(alpha, 0.25, -.50);

    const long incx = 2;
    array<T,6> x;
    assign_complex(x[0],   1, -  1);
    assign_complex(x[1], 555, -555);
    assign_complex(x[2],   2, -  2);
    assign_complex(x[3], 555, -555);
    assign_complex(x[4],   3, -  3);
    assign_complex(x[5], 555, -555);

    const long incy = 1;
    array<T,3> y;
    assign_complex(y[0],   4, -  4);
    assign_complex(y[1],   5, -  5);
    assign_complex(y[2],   6, -  6);

    BOOST_REQUIRE_EQUAL(x.size()/incx, y.size()/incy);
    blas::axpy(x.size()/incx, alpha, x.data(), incx, y.c_array(), incy);

    BOOST_CHECK_EQUAL(real(y[0]), (1.0* 0.25)-(-1.0*-0.50) + 4);
    BOOST_CHECK_EQUAL(imag(y[0]), (1.0*-0.50)+(-1.0* 0.25) - 4);
    BOOST_CHECK_EQUAL(real(y[1]), (2.0* 0.25)-(-2.0*-0.50) + 5);
    BOOST_CHECK_EQUAL(imag(y[1]), (2.0*-0.50)+(-2.0* 0.25) - 5);
    BOOST_CHECK_EQUAL(real(y[2]), (3.0* 0.25)-(-3.0*-0.50) + 6);
    BOOST_CHECK_EQUAL(imag(y[2]), (3.0*-0.50)+(-3.0* 0.25) - 6);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( complex_valued_real, T, complex_types )
{
    // Real scaling factor
    const typename suzerain::complex::traits::real<T>::type alpha = 0.25;

    const long incx = 2;
    array<T,6> x;
    assign_complex(x[0],   1, -  1);
    assign_complex(x[1], 555, -555);
    assign_complex(x[2],   2, -  2);
    assign_complex(x[3], 555, -555);
    assign_complex(x[4],   3, -  3);
    assign_complex(x[5], 555, -555);

    const long incy = 1;
    array<T,3> y;
    assign_complex(y[0],   4, -  4);
    assign_complex(y[1],   5, -  5);
    assign_complex(y[2],   6, -  6);

    BOOST_REQUIRE_EQUAL(x.size()/incx, y.size()/incy);
    blas::axpy(x.size()/incx, alpha, x.data(), incx, y.c_array(), incy);

    BOOST_CHECK_EQUAL(real(y[0]), ( 1.0*0.25) + 4);
    BOOST_CHECK_EQUAL(imag(y[0]), (-1.0*0.25) - 4);
    BOOST_CHECK_EQUAL(real(y[1]), ( 2.0*0.25) + 5);
    BOOST_CHECK_EQUAL(imag(y[1]), (-2.0*0.25) - 5);
    BOOST_CHECK_EQUAL(real(y[2]), ( 3.0*0.25) + 6);
    BOOST_CHECK_EQUAL(imag(y[2]), (-3.0*0.25) - 6);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( dot )

BOOST_AUTO_TEST_CASE_TEMPLATE( real_valued, T, real_types )
{
    const array<T,6> x = {{ 1, -1, 2, -1, 3, -1 }};
    const int incx = 2;
    const array<T,3> y = {{ 4, 5, 6 }};
    const long incy = 1;

    BOOST_REQUIRE_EQUAL(x.size()/incx, y.size()/incy);
    const T result = blas::dot(x.size()/incx, x.data(), incx, y.data(), incy);

    const T expected = (1*4) + (2*5) + (3*6);
    BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( complex_valued, T, complex_types )
{
    array<T,6> x;
    assign_complex(x[0],   1,    2);
    assign_complex(x[1], 555, -555);
    assign_complex(x[2],   2,    3);
    assign_complex(x[3], 555, -555);
    assign_complex(x[4],   3,    4);
    assign_complex(x[5], 555, -555);
    const int incx = 2;

    array<T,3> y;
    assign_complex(y[0],   2,    3);
    assign_complex(y[1],   3,    4);
    assign_complex(y[2],   4,    5);
    const long incy = 1;

    BOOST_REQUIRE_EQUAL(x.size()/incx, y.size()/incy);

    typedef std::complex<
            typename suzerain::complex::traits::real<T>::type
         > std_complex_t;
    const std_complex_t result
       = blas::dot(x.size()/incx, x.data(), incx, y.data(), incy);

    const std_complex_t expected(58, -3);
    BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( nrm2 )

BOOST_AUTO_TEST_CASE_TEMPLATE( real_valued, T, real_types )
{
    const array<T,6> x = {{ 3, -555, 4, -555 }};
    const long inc = 2;
    const T result = blas::nrm2(x.size()/inc, x.data(), inc);

    const T expected = 5;
    BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( complex_valued, T, complex_types )
{
    array<T,4> x;
    assign_complex(x[0],   3,    4);
    assign_complex(x[1], 555, -555);
    assign_complex(x[2],   5,    6);
    assign_complex(x[3], 555, -555);
    const long inc = 2;

    typedef typename suzerain::complex::traits::real<T>::type real_t;
    const real_t result   = blas::nrm2(x.size()/inc, x.data(), inc);
    const real_t expected = std::sqrt(86.0);
    const real_t epsilon  = std::numeric_limits<real_t>::epsilon();
    BOOST_CHECK_CLOSE(result, expected, epsilon);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( asum )

BOOST_AUTO_TEST_CASE_TEMPLATE( real_valued, T, real_types )
{
    const array<T,6> x = {{ -1, -555, 2, -555, -3, -555 }};
    const long inc = 2;
    const T result = blas::asum(x.size()/inc, x.data(), inc);

    const T expected = 1 + 2 + 3;
    BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( axpby )

BOOST_AUTO_TEST_CASE_TEMPLATE( real_valued, T, real_types )
{
    const array<T,6> x = {{ 1, -1, 2, -1, 3, -1 }};
    const long incx = 2;
    array<T,3> y = {{ 4, 5, 6 }};
    const int incy = 1;

    BOOST_REQUIRE_EQUAL(x.size()/incx, y.size()/incy);
    blas::axpby(x.size()/incx, 3, x.data(), incx, 7, y.c_array(), incy);

    const array<T,3> expected = {{ 1*3+7*4, 2*3+7*5, 3*3+7*6 }};
    BOOST_CHECK_EQUAL_COLLECTIONS(
            y.begin(), y.end(), expected.begin(), expected.end());
}

BOOST_AUTO_TEST_SUITE_END()
