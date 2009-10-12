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
#include <boost/scoped_array.hpp>
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
template<class ComplexMultiArray1, class ComplexMultiArray2>
void check_1D_complex_forward(ComplexMultiArray1 &in,
                              ComplexMultiArray2 &out,
                              const double dealias_by = 1.0)
{
    BOOST_STATIC_ASSERT(ComplexMultiArray1::dimensionality == 1);
    BOOST_STATIC_ASSERT(ComplexMultiArray2::dimensionality == 1);
    const int N = in.shape()[0];
    const double close_enough
        = std::numeric_limits<double>::epsilon()*10*N*N;

    // Load a real-valued function into the input array
    // Function includes only resolvable modes
    for (int i = 0; i < N; ++i) {
        double real_part = 1, imag_part = 0;
        const double xi = i*2*M_PI/N;
        for (int j = 0; j <= (N+1)/2; ++j) {
            real_part += j*sin(j*xi) - j*cos(j*xi);
        }
        fftw_multi_array::detail::assign_complex(in[i], real_part, imag_part);
    }

    fftw_multi_array::c2c_transform(0, in, out, FFTW_FORWARD, dealias_by);

    // Real input should exhibit conjugate symmetry in wave space
    BOOST_CHECK_SMALL(out[0].imag(), close_enough);
    for (int i = 1; i <= (N+1)/2; ++i) {
        double a_real, a_imag, b_real, b_imag;
        fftw_multi_array::detail::assign_components(a_real, a_imag, out[i]);
        fftw_multi_array::detail::assign_components(b_real, b_imag, out[N-i]);
        BOOST_CHECK_CLOSE(a_real,  b_real, close_enough);
        BOOST_CHECK_CLOSE(a_imag, -b_imag, close_enough);
    }

    // Load an imaginary-valued function into the input array
    // Function includes only resolvable modes
    for (int i = 0; i < N; ++i) {
        double real_part = 0, imag_part = 1;
        const double xi = i*2*M_PI/N;
        for (int j = 0; j <= (N+1)/2; ++j) {
            imag_part += j*sin(j*xi) - j*cos(j*xi);
        }
        fftw_multi_array::detail::assign_complex(in[i], real_part, imag_part);
    }

    fftw_multi_array::c2c_transform(0, in, out, FFTW_FORWARD, dealias_by);

    // Imaginary input should exhibit a similar symmetry in wave space:
    // Re(X_k) = - Re(X_{N-k}), Im(X_k) = Im(X_{N-k})
    BOOST_CHECK_SMALL(out[0].real(), close_enough);
    for (int i = 1; i <= (N+1)/2; ++i) {
        double a_real, a_imag, b_real, b_imag;
        fftw_multi_array::detail::assign_components(a_real, a_imag, out[i]);
        fftw_multi_array::detail::assign_components(b_real, b_imag, out[N-i]);
        BOOST_CHECK_CLOSE(a_real, -b_real, close_enough);
        BOOST_CHECK_CLOSE(a_imag,  b_imag, close_enough);
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
        double real_part = N, imag_part = N;
        const double xi = i*2*M_PI/N;
        for (int j = 0; j <= (N+1)/2; ++j) {
            real_part += j*cos(j*xi + M_PI/3);
            imag_part -= j*cos(j*xi + M_PI/5);
        }
        fftw_multi_array::detail::assign_complex(in[i], real_part, imag_part);
    }

    // Transform input FFTW directly and also our wrapper
    fftw_execute(plan.get());  // Important to be first for in == out
    fftw_multi_array::c2c_transform(0, in, out, FFTW_FORWARD, dealias_by);

    // Ensure we got exactly the same result
    for (int i = 0; i < N; ++i) {
        double out_real, out_imag;
        fftw_multi_array::detail::assign_components(
                out_real, out_imag, out[i]);
        BOOST_CHECK_EQUAL(out_real, buffer[i][0]);
        BOOST_CHECK_EQUAL(out_imag, buffer[i][1]);
    }
}

// Helper function that kicks the tires of a 1D c2c transform
template<class ComplexMultiArray1, class ComplexMultiArray2>
void check_1D_complex_backward(ComplexMultiArray1 &in,
                               ComplexMultiArray2 &out)
{
    BOOST_STATIC_ASSERT(ComplexMultiArray1::dimensionality == 1);
    BOOST_STATIC_ASSERT(ComplexMultiArray2::dimensionality == 1);
    const int N = in.shape()[0];
    const double close_enough
        = std::numeric_limits<double>::epsilon()*10*N*N;

    // Load a conjugate-symmetric function into the input array
    fftw_multi_array::detail::assign_complex(in[0], N, 0);
    for (int i = 1; i < (N+1)/2; ++i) {
        fftw_multi_array::detail::assign_complex(in[i],   i,  i);
        fftw_multi_array::detail::assign_complex(in[N-i], i, -i);
    }
    if (N%2 == 0) fftw_multi_array::detail::assign_complex(in[N/2], N/2, 0);

    fftw_multi_array::c2c_transform(0, in, out, FFTW_BACKWARD);

    // Output should be a real-valued function in physical space
    for (int i = 0; i < N; ++i) {
        double a_real, a_imag;
        fftw_multi_array::detail::assign_components(a_real, a_imag, out[i]);
        BOOST_CHECK_SMALL(a_imag, close_enough);
    }

    // Load a real-symmetric function into the input array
    fftw_multi_array::detail::assign_complex(in[0], 0, N);
    for (int i = 1; i <= (N+1)/2; ++i) {
        fftw_multi_array::detail::assign_complex(in[i],    i, i);
        fftw_multi_array::detail::assign_complex(in[N-i], -i, i);
    }
    if (N%2 == 0) fftw_multi_array::detail::assign_complex(in[N/2], 0, N/2);

    fftw_multi_array::c2c_transform(0, in, out, FFTW_BACKWARD);

    // Output should be an imaginary-valued function in physical space
    for (int i = 0; i < N; ++i) {
        double a_real, a_imag;
        fftw_multi_array::detail::assign_components(a_real, a_imag, out[i]);
        BOOST_CHECK_SMALL(a_real, close_enough);
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
                             FFTW_BACKWARD,
                             FFTW_PRESERVE_INPUT),
            std::ptr_fun(fftw_destroy_plan));
    BOOST_CHECK(plan.get() != NULL);

    // Load a complex-valued function into the input array
    for (int i = 0; i < N; ++i) {
        fftw_multi_array::detail::assign_complex(in[i], i, i);
    }

    // Transform input FFTW directly and also our wrapper
    fftw_execute(plan.get());  // Important to be first for in == out
    fftw_multi_array::c2c_transform(0, in, out, FFTW_BACKWARD);

    // Ensure we got exactly the same result after normalization
    for (int i = 0; i < N; ++i) {
        double out_real, out_imag;
        fftw_multi_array::detail::assign_components(
                out_real, out_imag, out[i]);
        BOOST_CHECK_EQUAL(out_real, buffer[i][0] / N);
        BOOST_CHECK_EQUAL(out_imag, buffer[i][1] / N);
    }
}

BOOST_AUTO_TEST_CASE( c2c_transform_1d )
{
    const int sizes[] = {
        2,  4,  8, 16, 32, 64                 // Easy: powers of 2
        , 2*3, 2*5, 2*7, 2*11, 2*13, 2*17     // Moderate: Even lengths
        , 3, 5, 7, 9, 11, 13, 17, 19, 23, 29  // Hard: Primes
    };

    // Out of place transformations on std::complex
    {
        typedef boost::multi_array<std::complex<double>,1> array_type;
        array_type in, out;
        for (int i = 0; i < sizeof(sizes)/sizeof(sizes[0]); ++i) {
            in.resize(boost::extents[sizes[i]]);
            out.resize(boost::extents[sizes[i]]);
            check_1D_complex_forward(in, out);
            check_1D_complex_backward(in, out);
        }
    }

    // In place transformations on std::complex
    {
        typedef boost::multi_array<std::complex<double>,1> array_type;
        array_type both;
        for (int i = 0; i < sizeof(sizes)/sizeof(sizes[0]); ++i) {
            both.resize(boost::extents[sizes[i]]);
            check_1D_complex_forward(both, both);
            check_1D_complex_backward(both, both);
        }
    }

    // TODO Out of place transformations on fftw_complex
    // TODO In place transformations on fftw_complex
    // Broken: boost::multi_array cannot handle fftw_complex elements
    // Awaiting response from boost-users mailing list
    // http://article.gmane.org/gmane.comp.lib.boost.user/52327
}
