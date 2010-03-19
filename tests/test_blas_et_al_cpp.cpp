#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <suzerain/complex.hpp>
#include <suzerain/blas_et_al.hpp>
#include <boost/test/included/unit_test.hpp>

#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

typedef boost::mpl::list<double,float> real_types;
typedef boost::mpl::list< std::complex<double>,
                          std::complex<float>  > complex_types;

// Introduce some shorthand
namespace blas = suzerain::blas;
using boost::array;

BOOST_AUTO_TEST_SUITE( swap )

BOOST_AUTO_TEST_CASE_TEMPLATE( real_valued, T, real_types )
{
    array<T,3> x = { 1, 2, 3 };
    array<T,3> y = { 4, 5, 6 };
    blas::swap(x.size(), x.c_array(), 1, y.c_array(), 1);

    const array<T,3> x_expected = { 4, 5, 6 };
    const array<T,3> y_expected = { 1, 2, 3 };
    BOOST_CHECK_EQUAL_COLLECTIONS(
            x.begin(), x.end(), x_expected.begin(), x_expected.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(
            y.begin(), y.end(), y_expected.begin(), y_expected.end());
}

// TODO Add complex_valued swap test

BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE( scal )

BOOST_AUTO_TEST_CASE_TEMPLATE( real_valued, T, real_types )
{
    { // Contiguous
        array<T,3> x = { 1, 2, 3 };
        blas::scal(x.size(), 2, x.c_array(), 1);

        const array<T,3> x_expected = { 2, 4, 6 };
        BOOST_CHECK_EQUAL_COLLECTIONS(
                x.begin(), x.end(), x_expected.begin(), x_expected.end());
    }

    { // Strided
        array<T,6> x = { 1, -1, 2, -1, 3, -1 };
        const long inc = 2;
        blas::scal(x.size()/inc, 2, x.c_array(), inc);

        const array<T,6> x_expected = { 2, -1, 4, -1, 6, -1 };
        BOOST_CHECK_EQUAL_COLLECTIONS(
                x.begin(), x.end(), x_expected.begin(), x_expected.end());
    }
}

// TODO Add complex_valued scal test with complex scaling factor
// TODO Add complex_valued scal test with real scaling factor

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( copy )

BOOST_AUTO_TEST_CASE_TEMPLATE( real_valued, T, real_types )
{
    const array<T,6> x = {   1, -555,   2, -555,   3, -555 };
    array<T,3> y       = { 555,  555, 555 };
    const long incx = 2, incy = 1;
    blas::copy(x.size()/incx, x.data(), incx, y.c_array(), incy);

    const array<T,3> expected = { 1, 2, 3 };
    BOOST_CHECK_EQUAL_COLLECTIONS(
            y.begin(), y.end(), expected.begin(), expected.end());
}

// TODO Add complex_valued copy test

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( axpy )

BOOST_AUTO_TEST_CASE_TEMPLATE( real_valued, T, real_types )
{
    const array<T,6> x = { 1, -1, 2, -1, 3, -1 };
    const long incx = 2;
    array<T,3> y = { 4, 5, 6 };
    const int incy = 1;
    blas::axpy(x.size()/incx, 3, x.data(), incx, y.c_array(), incy);

    const array<T,3> expected = { 1*3+4, 2*3+5, 3*3+6 };
    BOOST_CHECK_EQUAL_COLLECTIONS(
            y.begin(), y.end(), expected.begin(), expected.end());
}

// TODO Add complex_valued axpy test with complex-valued alpha
// TODO Add complex_valued axpy test with real-valued alpha

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( dot )

BOOST_AUTO_TEST_CASE_TEMPLATE( real_valued, T, real_types )
{
    const array<T,6> x = { 1, -1, 2, -1, 3, -1 };
    const int incx = 2;
    const array<T,3> y = { 4, 5, 6 };
    const long incy = 1;
    const T result = blas::dot(x.size()/incx, x.data(), incx, y.data(), incy);

    const T expected = (1*4) + (2*5) + (3*6);
    BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( asum )

BOOST_AUTO_TEST_CASE_TEMPLATE( real_valued, T, real_types )
{
    const array<T,6> x = { -1, -555, 2, -555, -3, -555 };
    const long inc = 2;
    const T result = blas::asum(x.size()/inc, x.data(), inc);

    const T expected = 1 + 2 + 3;
    BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( axpby )

BOOST_AUTO_TEST_CASE_TEMPLATE( real_valued, T, real_types )
{
    const array<T,6> x = { 1, -1, 2, -1, 3, -1 };
    const long incx = 2;
    array<T,3> y = { 4, 5, 6 };
    const int incy = 1;
    blas::axpby(x.size()/incx, 3, x.data(), incx, 7, y.c_array(), incy);

    const array<T,3> expected = { 1*3+7*4, 2*3+7*5, 3*3+7*6 };
    BOOST_CHECK_EQUAL_COLLECTIONS(
            y.begin(), y.end(), expected.begin(), expected.end());
}

BOOST_AUTO_TEST_SUITE_END()
