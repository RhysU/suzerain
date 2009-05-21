#define BOOST_TEST_MODULE $Id$

#include <suzerain/config.h>

#include <algorithm>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/foreach.hpp>
#include <complex>
#include <typeinfo>

namespace ublas = boost::numeric::ublas;

typedef std::complex<double>                            complex;
typedef ublas::shallow_array_adaptor<double>            shallow_adaptor_double;
typedef ublas::vector<double, shallow_adaptor_double>   shallow_vector_double;
typedef ublas::shallow_array_adaptor<complex>           shallow_adaptor_complex;
typedef ublas::vector<complex, shallow_adaptor_complex> shallow_vector_complex;

// Ensure we can use std::complex<double> as two consecutive doubles
BOOST_AUTO_TEST_CASE( shared_c_array )
{
    // Assumption must hold true for any of this scheme to work
    BOOST_REQUIRE_EQUAL( sizeof(complex), 2*sizeof(double) );

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

// Ensure we can use a C-array of doubles as a ublas::vector.
BOOST_AUTO_TEST_CASE( shallow_array_adaptor_double )
{
    const std::size_t N = 3;
    double raw[N];

    shallow_adaptor_double adaptor(N, raw);
    shallow_vector_double  vec(N, adaptor);

    BOOST_CHECK_EQUAL( vec.size(), N );

    std::fill(&raw[0], &raw[N], 1.0);
    BOOST_CHECK_EQUAL_COLLECTIONS(&raw[0], &raw[N], vec.begin(), vec.end());

    std::fill(vec.begin(), vec.end(), 2.0);
    BOOST_CHECK_EQUAL_COLLECTIONS(&raw[0], &raw[N], vec.begin(), vec.end());
}

// Ensure we can use a C-array of std::complex<double> as a ublas::vector.
BOOST_AUTO_TEST_CASE( shallow_array_adaptor_complex )
{
    const std::size_t N = 3;
    complex raw[N];

    shallow_adaptor_complex adaptor(N, raw);
    shallow_vector_complex  vec(N, adaptor);

    BOOST_CHECK_EQUAL( vec.size(), N );

    std::fill(&raw[0], &raw[N], complex(1.0, -1.0));
    BOOST_CHECK_EQUAL_COLLECTIONS(&raw[0], &raw[N], vec.begin(), vec.end());

    std::fill(vec.begin(), vec.end(), complex(2.0, -2.0));
    BOOST_CHECK_EQUAL_COLLECTIONS(&raw[0], &raw[N], vec.begin(), vec.end());
}

// Ensure we can use the same data as either double or complex ublas::vector
BOOST_AUTO_TEST_CASE( shallow_array_adaptor_shared )
{
    const std::size_t N = 6;
    double carray_double[N] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
    complex *carray_complex = reinterpret_cast<complex *>(carray_double);

    shallow_adaptor_double adaptor_double(N, carray_double);
    shallow_vector_double  vec_double(N, adaptor_double);

    shallow_adaptor_complex adaptor_complex(N / 2, carray_complex);
    shallow_vector_complex  vec_complex(N / 2, adaptor_complex);

    for (std::size_t i = 0; i < N / 2; ++i) {
        BOOST_CHECK_EQUAL(vec_double[2*i],   vec_complex[i].real());
        BOOST_CHECK_EQUAL(vec_double[2*i+1], vec_complex[i].imag());
    }

    for (std::size_t i = 0; i < N / 2; ++i) {
        vec_complex[i] += vec_complex[0];
    }

    for (std::size_t i = 0; i < N / 2; ++i) {
        BOOST_CHECK_EQUAL(vec_double[2*i],   vec_complex[i].real());
        BOOST_CHECK_EQUAL(vec_double[2*i+1], vec_complex[i].imag());
    }

    for (std::size_t i = 0; i < N; ++i) {
        vec_double[i] += vec_double[N-1];
    }

    for (std::size_t i = 0; i < N / 2; ++i) {
        BOOST_CHECK_EQUAL(vec_double[2*i],   vec_complex[i].real());
        BOOST_CHECK_EQUAL(vec_double[2*i+1], vec_complex[i].imag());
    }
}
