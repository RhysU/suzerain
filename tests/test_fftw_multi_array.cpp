#define BOOST_TEST_MODULE $Id$
#include <cmath>
#include <complex>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <string>
#include <vector>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/static_assert.hpp>
#include <boost/test/included/unit_test.hpp>
#include <suzerain/fftw_multi_array.hpp>
#include <fftw3.h>

#include "test_tools.hpp"

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

template<class MultiArray>
void debug_dump(const std::string &prefix, const MultiArray &x)
{
    BOOST_STATIC_ASSERT(MultiArray::dimensionality == 1);
    using namespace std;

    cout << prefix;
    copy(x.begin(), x.end(),
         ostream_iterator<typename MultiArray::element>(std::cout, " "));
    cout << endl;
}

void debug_dump(const std::string &prefix,
               const fftw_complex * begin,
               const fftw_complex * end) {
    using namespace std;

    cout << prefix;
    while (begin != end) {
        cout << "("
             << (*begin)[0]
             << ","
             << (*begin)[1]
             << ")";
        begin++;
        if (begin != end) cout << ", ";
    }
    cout << endl;
}

// Helper function that kicks the tires of a 1D c2c transform
template<class ComplexMultiArray>
void check_1D_complex(ComplexMultiArray &in, ComplexMultiArray &out)
{
    BOOST_STATIC_ASSERT(ComplexMultiArray::dimensionality == 1);
    const int N = in.shape()[0];
    const double close_enough
        = std::numeric_limits<double>::epsilon()*10*N*log(N);

    // Load a real-valued function into the input array
    for (int i = 0; i < N; ++i) {
        in[i].real() = 1;
        in[i].imag() = 0;
        const double xi = i*2*M_PI/N;
        for (int j = 0; j <= (N+1)/2; ++j) {
            in[i].real() += j*sin(j*xi) - j*cos(j*xi);
        }
    }

    fftw_multi_array::c2c_transform(0, in, out, FFTW_FORWARD);

    // Real input should exhibit conjugate symmetry in wave space
    BOOST_CHECK_SMALL(out[0].imag(), close_enough);
    for (int i = 1; i <= (N+1)/2; ++i) {
        BOOST_CHECK_CLOSE(out[i].real(),  out[N-i].real(), close_enough);
        BOOST_CHECK_CLOSE(out[i].imag(), -out[N-i].imag(), close_enough);
    }

    // Load an imaginary-valued function into the input array
    for (int i = 0; i < N; ++i) {
        in[i].real() = 0;
        in[i].imag() = 1;
        const double xi = i*2*M_PI/N;
        for (int j = 0; j <= (N+1)/2; ++j) {
            in[i].imag() += j*sin(j*xi) - j*cos(j*xi);
        }
    }

    fftw_multi_array::c2c_transform(0, in, out, FFTW_FORWARD);

    // Imaginary input should exhibit a similar symmetry in wave space:
    // Re(X_k) = - Re(X_{N-k}), Im(X_k) = Im(X_{N-k})
    BOOST_CHECK_SMALL(out[0].real(), close_enough);
    for (int i = 1; i <= (N+1)/2; ++i) {
        BOOST_CHECK_CLOSE(out[i].real(), -out[N-i].real(), close_enough);
        BOOST_CHECK_CLOSE(out[i].imag(),  out[N-i].imag(), close_enough);
    }

    // Compare our raw double results with FFTW's directly computed result
    // Plan before loading in the data since planning overwrites in
    boost::shared_array<fftw_complex> buffer(
        static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex)*N)),
        std::ptr_fun(fftw_free));
    BOOST_CHECK(buffer.get() != NULL);
    boost::shared_ptr<boost::remove_pointer<fftw_plan>::type> plan(
            fftw_plan_dft_1d(N,
                             reinterpret_cast<fftw_complex*>(in.data()),
                             buffer.get(),
                             FFTW_FORWARD,
                             FFTW_PRESERVE_INPUT),
            std::ptr_fun(fftw_destroy_plan));
    BOOST_CHECK(plan.get() != NULL);

    // Load a complex-valued function into the input array
    for (int i = 0; i < N; ++i) {
        in[i].real() = N;
        in[i].imag() = N;
        const double xi = i*2*M_PI/N;
        for (int j = 0; j <= (N+1)/2; ++j) {
            in[i].real() += j*cos(j*xi + M_PI/3);
            in[i].imag() -= j*cos(j*xi + M_PI/5);
        }
    }

    // Transform input FFTW directly and also our wrapper
    fftw_execute(plan.get());  // Important to be first for in == out
    fftw_multi_array::c2c_transform(0, in, out, FFTW_FORWARD);

    // Ensure we got exactly the same result
    for (int i = 0; i < N; ++i) {
        BOOST_CHECK_EQUAL(out[i].real(), buffer[i][0]);
        BOOST_CHECK_EQUAL(out[i].imag(), buffer[i][1]);
    }
}

BOOST_AUTO_TEST_CASE( c2c_transform_1d )
{
    typedef boost::multi_array<std::complex<double>,1> array_type;

    const int sizes[] = {
        2,  4,  8, 16, 32, 64                 // Easy: powers of 2
        , 2*3, 2*5, 2*7, 2*11, 2*13, 2*17     // Moderate: Even lengths
        , 3, 5, 7, 9, 11, 13, 17, 19, 23, 29  // Hard: Primes
    };

    // Out of place transformations
    {
        array_type in, out;
        for (int i = 0; i < sizeof(sizes)/sizeof(sizes[0]); ++i) {
            in.resize(boost::extents[sizes[i]]);
            out.resize(boost::extents[sizes[i]]);
            check_1D_complex(in, out);
        }
    }

    // In place transformations
    {
        array_type both;
        for (int i = 0; i < sizeof(sizes)/sizeof(sizes[0]); ++i) {
            both.resize(boost::extents[sizes[i]]);
            check_1D_complex(both, both);
        }
    }
}
