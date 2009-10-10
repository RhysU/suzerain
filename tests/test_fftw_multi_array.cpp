#define BOOST_TEST_MODULE $Id$
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <numeric>
#include <vector>
#include <boost/array.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/multi_array.hpp>
#include <suzerain/fftw_multi_array.hpp>

#include "test_tools.hpp"

/* DEBUG */
#include <iostream>
#include <iterator>

using namespace pecos::suzerain;

BOOST_AUTO_TEST_CASE( increment_1d_degenerate )
{
    const int n        = 1;
    int       index[n] = { 0 };
    const int shape[n] = { 1 };

    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_1d_normal )
{
    const int n        = 1;
    int       index[n] = { 0 };
    const int shape[n] = { 3 };

    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_degenerate_first )
{
    const int           n     = 2;
    boost::array<int,n> index = {{ 0, 0 }};
    boost::array<int,n> shape = {{ 3, 1 }};

    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_degenerate_second )
{
    const int                   n     = 2;
    boost::array<int, n>        index = {{ 0, 0 }};
    boost::array<std::size_t,n> shape = {{ 1, 3 }};

    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 2);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_degenerate_all )
{
    const int  n        = 2;
    int        index[n] = { 0, 0 };
    const long shape[n] = { 1, 1 };

    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_normal )
{
    const int    n        = 2;
    unsigned int index[n] = { 0, 0 };
    unsigned int shape[n] = { 2, 2 };

    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_3d_degenerate_all )
{
    const int         n        = 3;
    signed int        index[n] = { 0, 0, 0 };
    const signed long shape[n] = { 1, 1, 1 };

    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_3d_degenerate_middle )
{
    const int                   n        = 3;
    long                        index[n] =  { 0, 0, 0 };
    boost::array<std::size_t,n> shape    = {{ 3, 1, 3 }};

    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 2);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 2);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 2);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_3d_normal )
{
    const int n = 3;
    std::vector<short> index(n, 0);
    std::vector<int> shape(n, 2);

    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

template<class ComplexMultiArray>
void check_1D_conjugate_symmetry(ComplexMultiArray &in, ComplexMultiArray &out)
{
    BOOST_CHECK_EQUAL(in.num_dimensions(), 1); // Test sanity
    const double close_enough = std::numeric_limits<double>::epsilon();
    const int N = in.shape()[0];

    // Load a real-valued function into the input array
    for (int i = 0; i < N; ++i) {
        in[i].real() = 1;
        in[i].imag() = 0;
        const double xi = i*2*M_PI/N;
        for (int j = 0; j < N/2; ++j) {
            in[i].real() += i*sin(i*xi) - i*cos(i*xi);
        }
    }

    fftw_multi_array::c2c_transform(0, in, out, FFTW_FORWARD);

    // Real input should exhibit conjugate symmetry in wave space
    BOOST_CHECK_SMALL(out[0].imag(), close_enough);
    if (N%2 == 0) {
        BOOST_CHECK_SMALL(out[N/2].imag(), close_enough);
    }
    for (int i = 1; i < N/2; ++i) {
        BOOST_CHECK_CLOSE(out[i].real(),  out[N-i].real(), close_enough);
        BOOST_CHECK_CLOSE(out[i].imag(), -out[N-i].imag(), close_enough);
    }

    // Load an imaginary-valued function into the input array
    for (int i = 0; i < N; ++i) {
        in[i].real() = 0;
        in[i].imag() = 1;
        const double xi = i*2*M_PI/N;
        for (int j = 0; j < N/2; ++j) {
            in[i].imag() += i*sin(i*xi) - i*cos(i*xi);
        }
    }

    fftw_multi_array::c2c_transform(0, in, out, FFTW_FORWARD);

    // Imaginary input should exhibit a similar symmetry in wave space:
    // Re(X_k) = - Re(X_{N-k}), Im(X_k) = Im(X_{N-k})
    BOOST_CHECK_SMALL(out[0].real(), close_enough);
    if (N%2 == 0) {
        BOOST_CHECK_SMALL(out[N/2].real(), close_enough);
    }
    for (int i = 1; i < N/2; ++i) {
        BOOST_CHECK_CLOSE(out[i].real(), -out[N-i].real(), close_enough);
        BOOST_CHECK_CLOSE(out[i].imag(),  out[N-i].imag(), close_enough);
    }
}

BOOST_AUTO_TEST_CASE( c2c_transform_1d_real_input )
{
    typedef boost::multi_array<std::complex<double>,1> array_type;
    array_type in, out;

    const int sizes[] = {
        4,     8,  16,   32,   64,        // Easy: powers of 2
        2*3, 2*5, 2*7, 2*11, 2*13, 2*17   // Moderate: Even lengths
// DISABLED 5, 9, 17, 33, 65  // Harder
    };
    for (int i = 0; i < sizeof(sizes)/sizeof(sizes[0]); ++i) {
        in.resize(boost::extents[sizes[i]]);
        out.resize(boost::extents[sizes[i]]);
        check_1D_conjugate_symmetry(in, out);
    }
}
