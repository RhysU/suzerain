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
#include <boost/preprocessor/comparison/greater.hpp>
#include <boost/preprocessor/comparison/less.hpp>
#include <boost/preprocessor/control/iif.hpp>
#include <boost/preprocessor/facilities/empty.hpp>
#include <boost/preprocessor/seq/elem.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/for_each_product.hpp>
#include <boost/preprocessor/seq/size.hpp>
#include <boost/preprocessor/seq/to_tuple.hpp>
#include <boost/scoped_array.hpp>
#include <boost/static_assert.hpp>
#include <boost/test/included/unit_test.hpp>
#include <suzerain/fftw_multi_array.hpp>
#include <fftw3.h>

#include "test_tools.hpp"

using namespace pecos::suzerain;

BOOST_AUTO_TEST_SUITE( increment )

BOOST_AUTO_TEST_CASE( increment_1d_degenerate )
{
    const int n        = 1;
    int       index[n] = { 0 };
    const int shape[n] = { 1 };

    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_1d_normal )
{
    const int n        = 1;
    int       index[n] = { 0 };
    const int shape[n] = { 3 };

    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 2);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_degenerate_first )
{
    const int           n     = 2;
    boost::array<int,n> index = {{ 0, 0 }};
    boost::array<int,n> shape = {{ 3, 1 }};

    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 2);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_degenerate_second )
{
    const int                   n     = 2;
    boost::array<int, n>        index = {{ 0, 0 }};
    boost::array<std::size_t,n> shape = {{ 1, 3 }};

    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_degenerate_all )
{
    const int  n        = 2;
    int        index[n] = { 0, 0 };
    const long shape[n] = { 1, 1 };

    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_normal )
{
    const int    n        = 2;
    unsigned int index[n] = { 0, 0 };
    unsigned int shape[n] = { 2, 2 };

    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_3d_degenerate_all )
{
    const int         n        = 3;
    signed int        index[n] = { 0, 0, 0 };
    const signed long shape[n] = { 1, 1, 1 };

    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_3d_degenerate_middle )
{
    const int                   n        = 3;
    long                        index[n] =  { 0, 0, 0 };
    boost::array<std::size_t,n> shape    = {{ 3, 1, 3 }};

    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 2);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 2);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 2);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_3d_normal )
{
    const int n = 3;
    std::vector<short> index(n, 0);
    std::vector<int> shape(n, 2);

    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_SUITE_END()

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

// Produce a real signal with known frequency content
template<typename FPT, typename Integer>
FPT real_test_function(const Integer NR,
                       const Integer max_mode,
                       const Integer i,
                       const FPT shift = M_PI/3.0) {
    const FPT xi = i*2*M_PI/NR;
    FPT retval = (max_mode >= 0) ? NR : 0;
    for (Integer i = 1; i <= max_mode; ++i) {
        retval += i*sin(i*xi + shift);
    }
    return retval;
}

template<class ComplexMultiArray>
void fill_with_complex_NaN(ComplexMultiArray &x)
{
    typedef typename ComplexMultiArray::element element;
    const typename element::value_type quiet_NaN
        = std::numeric_limits<typename element::value_type>::quiet_NaN();
    const element nan_value(quiet_NaN, quiet_NaN);
    std::fill(x.data(), x.data() + x.num_elements(), nan_value);
}

// Helper function that kicks the tires of a 1D c2c transform
template<class ComplexMultiArray1, class ComplexMultiArray2>
void check_1D_complex_forward(ComplexMultiArray1 &in, ComplexMultiArray2 &out)
{
    BOOST_STATIC_ASSERT(ComplexMultiArray1::dimensionality == 1);
    BOOST_STATIC_ASSERT(ComplexMultiArray2::dimensionality == 1);
    const int NR = in.shape()[0];
    const int NC = out.shape()[0];
    const double close_enough
        = std::numeric_limits<double>::epsilon()*10*NR*NR;

    // Load a real-valued function into the input array and transform it
    fill_with_complex_NaN(in);
    fill_with_complex_NaN(out);
    for (int i = 0; i < NR; ++i) {
        fftw_multi_array::detail::assign_complex(
                in[i], real_test_function<double>(NR, (NR+1)/2, i), 0.0);
    }
    fftw_multi_array::c2c_transform(0, in, out, FFTW_FORWARD);

    // Real input should exhibit conjugate symmetry in wave space...
    BOOST_REQUIRE_SMALL(out[0].imag(), close_enough);
    for (int i = 1; i < (std::min(NC,NR)+1)/2; ++i) { // ...up to grid modes
        double a_real, a_imag, b_real, b_imag;
        fftw_multi_array::detail::assign_components(a_real, a_imag, out[i]);
        fftw_multi_array::detail::assign_components(b_real, b_imag, out[NC-i]);
        BOOST_REQUIRE_CLOSE(a_real,  b_real, close_enough);
        BOOST_REQUIRE_CLOSE(a_imag, -b_imag, close_enough);
    }

    // Load an imaginary-valued function into the input array and transform it
    fill_with_complex_NaN(in);
    fill_with_complex_NaN(out);
    for (int i = 0; i < NR; ++i) {
        fftw_multi_array::detail::assign_complex(
                in[i], 0.0, real_test_function<double>(NR, (NR+1)/2, i));
    }
    fftw_multi_array::c2c_transform(0, in, out, FFTW_FORWARD);

    // Imaginary input should exhibit a similar symmetry in wave space...
    // Re(X_k) = - Re(X_{N-k}), Im(X_k) = Im(X_{N-k})
    BOOST_REQUIRE_SMALL(out[0].real(), close_enough);
    for (int i = 1; i < (std::min(NC,NR)+1)/2; ++i) { // ...up to grid modes
        double a_real, a_imag, b_real, b_imag;
        fftw_multi_array::detail::assign_components(a_real, a_imag, out[i]);
        fftw_multi_array::detail::assign_components(b_real, b_imag, out[NC-i]);
        BOOST_REQUIRE_CLOSE(a_real, -b_real, close_enough);
        BOOST_REQUIRE_CLOSE(a_imag,  b_imag, close_enough);
    }
}

// Compare our results with FFTW's directly computed result
template<class ComplexMultiArray1, class ComplexMultiArray2>
void compare_1D_complex_forward(ComplexMultiArray1 &in,
                                ComplexMultiArray2 &out)
{
    BOOST_STATIC_ASSERT(ComplexMultiArray1::dimensionality == 1);
    BOOST_STATIC_ASSERT(ComplexMultiArray2::dimensionality == 1);
    const int NC = in.shape()[0];
    const int NR = out.shape()[0];
    const double close_enough
        = std::numeric_limits<double>::epsilon()*10*NC*NC;

    // Plan before loading in the data since planning overwrites in
    boost::shared_array<fftw_complex> buffer(
        static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex)*NR)),
        std::ptr_fun(fftw_free));
    BOOST_REQUIRE(buffer.get() != NULL);
    typename ComplexMultiArray1::element * const in_data
        = in.origin() + std::inner_product(
            in.index_bases(), in.index_bases()+1, in.strides(), 0);
    boost::shared_ptr<boost::remove_pointer<fftw_plan>::type> plan(
            fftw_plan_many_dft(1,
                               &NR,
                               1,
                               reinterpret_cast<fftw_complex*>(in_data),
                               NULL,
                               in.strides()[0],
                               NR,
                               buffer.get(),
                               NULL,
                               1,
                               NR,
                               FFTW_FORWARD,
                               FFTW_PRESERVE_INPUT),
            std::ptr_fun(fftw_destroy_plan));
    BOOST_REQUIRE(plan.get() != NULL);

    // Load a complex-valued function into the input array
    fill_with_complex_NaN(in);
    fill_with_complex_NaN(out);
    for (int i = 0; i < NR; ++i) {
        fftw_multi_array::detail::assign_complex(
                in[i],
                real_test_function<double>(NR, (NR+1)/2, i, M_PI/3.0),
                real_test_function<double>(NR, (NR+1)/2, i, M_PI/5.0));
    }

    // Transform input FFTW directly and also our wrapper
    fftw_execute(plan.get());  // Important to be first for in == out
    fftw_multi_array::c2c_transform(0, in, out, FFTW_FORWARD);

    // Ensure we got the same result
    for (int i = 0; i < NR; ++i) {
        double out_real, out_imag;
        fftw_multi_array::detail::assign_components(
                out_real, out_imag, out[i]);
        // FFTW with a stride gives a different result than with stride 1
        // BOOST_REQUIRE_EQUAL would be nice, but it fails here
        BOOST_REQUIRE_CLOSE(out_real, buffer[i][0], close_enough/(NR*NR));
        BOOST_REQUIRE_CLOSE(out_imag, buffer[i][1], close_enough/(NR*NR));
    }
}

// Helper function that kicks the tires of a 1D c2c transform
template<class ComplexMultiArray1, class ComplexMultiArray2>
void check_1D_complex_backward(ComplexMultiArray1 &in, ComplexMultiArray2 &out)
{
    BOOST_STATIC_ASSERT(ComplexMultiArray1::dimensionality == 1);
    BOOST_STATIC_ASSERT(ComplexMultiArray2::dimensionality == 1);
    const int NC = in.shape()[0];
    const int NR = out.shape()[0];
    const double close_enough
        = std::numeric_limits<double>::epsilon()*10*NC*NC;

    // Load a conjugate-symmetric function into the input array...
    fill_with_complex_NaN(in);
    fill_with_complex_NaN(out);
    fftw_multi_array::detail::assign_complex(in[0], NC, 0);
    for (int i = 1; i < (NC+1)/2; ++i) {
        if (i < (NR+1)/2) { // ...up to grid modes
            fftw_multi_array::detail::assign_complex(in[i],    i,  i);
            fftw_multi_array::detail::assign_complex(in[NC-i], i, -i);
        } else {
            fftw_multi_array::detail::assign_complex(in[i],    0, 0);
            fftw_multi_array::detail::assign_complex(in[NC-i], 0, 0);
        }
    }
    if (NC%2 == 0) {
        fftw_multi_array::detail::assign_complex(in[NC/2],
                (NC >= NR) ? NC/2 : 0, 0);
    }

    fftw_multi_array::c2c_transform(0, in, out, FFTW_BACKWARD);

    // Output should be a real-valued function in physical space
    for (int i = 0; i < NR; ++i) {
        double a_real, a_imag;
        fftw_multi_array::detail::assign_components(a_real, a_imag, out[i]);
        BOOST_REQUIRE_SMALL(a_imag, close_enough);
    }

    // Load a real-symmetric function into the input array...
    fill_with_complex_NaN(in);
    fill_with_complex_NaN(out);
    fftw_multi_array::detail::assign_complex(in[0], 0, NC);
    for (int i = 1; i < (NC+1)/2; ++i) {
        if (i < (NR+1)/2) { // ...up to grid modes
            fftw_multi_array::detail::assign_complex(in[i],     i, i);
            fftw_multi_array::detail::assign_complex(in[NC-i], -i, i);
        } else {
            fftw_multi_array::detail::assign_complex(in[i],    0, 0);
            fftw_multi_array::detail::assign_complex(in[NC-i], 0, 0);
        }
    }
    if (NC%2 == 0) {
        fftw_multi_array::detail::assign_complex(in[NC/2],
                0, (NC >= NR) ? NC/2 : 0);
    }

    fftw_multi_array::c2c_transform(0, in, out, FFTW_BACKWARD);

    // Output should be an imaginary-valued function in physical space
    for (int i = 0; i < NR; ++i) {
        double a_real, a_imag;
        fftw_multi_array::detail::assign_components(a_real, a_imag, out[i]);
        BOOST_REQUIRE_SMALL(a_real, close_enough);
    }
}

// Compare our results with FFTW's directly computed result
template<class ComplexMultiArray1, class ComplexMultiArray2>
void compare_1D_complex_backward(ComplexMultiArray1 &in,
                                 ComplexMultiArray2 &out)
{
    BOOST_STATIC_ASSERT(ComplexMultiArray1::dimensionality == 1);
    BOOST_STATIC_ASSERT(ComplexMultiArray2::dimensionality == 1);
    const int NC = in.shape()[0];
    const int NR = out.shape()[0];
    const double close_enough
        = std::numeric_limits<double>::epsilon()*10*NC*NC;

    // Plan before loading in the data since planning overwrites in
    boost::shared_array<fftw_complex> buffer(
        static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex)*NC)),
        std::ptr_fun(fftw_free));
    BOOST_REQUIRE(buffer.get() != NULL);
    typename ComplexMultiArray1::element * const in_data
        = in.origin() + std::inner_product(
            in.index_bases(), in.index_bases()+1, in.strides(), 0);
    boost::shared_ptr<boost::remove_pointer<fftw_plan>::type> plan(
            fftw_plan_many_dft(1,
                               &NC,
                               1,
                               reinterpret_cast<fftw_complex*>(in_data),
                               NULL,
                               in.strides()[0],
                               NC,
                               buffer.get(),
                               NULL,
                               1,
                               NC,
                               FFTW_BACKWARD,
                               FFTW_PRESERVE_INPUT),
            std::ptr_fun(fftw_destroy_plan));
    BOOST_REQUIRE(plan.get() != NULL);

    // Load a complex-valued function into the input array
    fill_with_complex_NaN(in);
    fill_with_complex_NaN(out);
    for (int i = 0; i < NC; ++i) {
        fftw_multi_array::detail::assign_complex(in[i], i, i);
    }

    // Transform input FFTW directly and also our wrapper
    fftw_execute(plan.get());  // Important to be first for in == out
    fftw_multi_array::c2c_transform(0, in, out, FFTW_BACKWARD);

    // Ensure we got exactly the same result after normalization
    for (int i = 0; i < NR; ++i) {
        double out_real, out_imag;
        fftw_multi_array::detail::assign_components(
                out_real, out_imag, out[i]);
        // FFTW with a stride gives a different result than with stride 1
        // BOOST_REQUIRE_EQUAL would be nice, but it fails here
        BOOST_REQUIRE_CLOSE(out_real, buffer[i][0] / NR, close_enough/10);
        BOOST_REQUIRE_CLOSE(out_imag, buffer[i][1] / NR, close_enough/10);
    }
}


/* powers of 2; simple composites of form 2*prime; prime numbers; 1 */
#define TRANSFORM_1D_SIZE_SEQ \
        (2)(4)(8)(16)(32)(64) \
        (6)(10)(14)(22)(26)(34) \
        (3)(5)(7)(9)(11)(13)(17)(19)(23)(29) \
        (1)

BOOST_AUTO_TEST_SUITE( c2c_1d_out_of_place );
#define TEST_C2C_1D_OUT_OF_PLACE(r, data, elem) \
        BOOST_AUTO_TEST_CASE( BOOST_PP_CAT(c2c_1d_out_of_place_,elem) ) \
        { c2c_1d_out_of_place(elem); }
void c2c_1d_out_of_place(const int N)
{
    typedef boost::multi_array<std::complex<double>,1> array_type;
    array_type in(boost::extents[N]), out(boost::extents[N]);

    // No dealiasing in effect: in == out
    check_1D_complex_forward(in, out);
    compare_1D_complex_forward(in, out);
    check_1D_complex_backward(in, out);
    compare_1D_complex_backward(in, out);
}
BOOST_PP_SEQ_FOR_EACH(TEST_C2C_1D_OUT_OF_PLACE,_,TRANSFORM_1D_SIZE_SEQ);
BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE( c2c_1d_out_of_place_one_reversed );
#define TEST_C2C_1D_OUT_OF_PLACE_ONE_REVERSED(r, data, elem) \
        BOOST_AUTO_TEST_CASE( \
          BOOST_PP_CAT(c2c_1d_out_of_place_one_reversed_,elem) ) \
        { c2c_1d_out_of_place_one_reversed(elem); }
void c2c_1d_out_of_place_one_reversed(const int N)
{
    typedef boost::multi_array<std::complex<double>,1> array_type;
    typedef boost::general_storage_order<array_type::dimensionality> storage;
    array_type::size_type ordering[array_type::dimensionality] = { 0 };
    const bool ascending[array_type::dimensionality] = { false };

    array_type in(boost::extents[N], storage(ordering, ascending));
    array_type out(boost::extents[N]);

    check_1D_complex_forward(in, out);
    compare_1D_complex_forward(in, out);
    check_1D_complex_backward(in, out);
    compare_1D_complex_backward(in, out);
}
BOOST_PP_SEQ_FOR_EACH(TEST_C2C_1D_OUT_OF_PLACE_ONE_REVERSED,\
    _,TRANSFORM_1D_SIZE_SEQ);
BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE( c2c_1d_out_of_place_two_reversed );
#define TEST_C2C_1D_OUT_OF_PLACE_TWO_REVERSED(r, data, elem) \
        BOOST_AUTO_TEST_CASE( \
          BOOST_PP_CAT(c2c_1d_out_of_place_two_reversed_,elem) ) \
        { c2c_1d_out_of_place_two_reversed(elem); }
void c2c_1d_out_of_place_two_reversed(const int N)
{
    typedef boost::multi_array<std::complex<double>,1> array_type;
    typedef boost::general_storage_order<array_type::dimensionality> storage;
    array_type::size_type ordering[array_type::dimensionality] = { 0 };
    const bool ascending[array_type::dimensionality] = { false };

    array_type in(boost::extents[N], storage(ordering, ascending));
    array_type out(boost::extents[N], storage(ordering, ascending));

    check_1D_complex_forward(in, out);
    compare_1D_complex_forward(in, out);
    check_1D_complex_backward(in, out);
    compare_1D_complex_backward(in, out);
}
BOOST_PP_SEQ_FOR_EACH(TEST_C2C_1D_OUT_OF_PLACE_TWO_REVERSED,\
    _,TRANSFORM_1D_SIZE_SEQ);
BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE( c2c_1d_in_place );
#define TEST_C2C_1D_IN_PLACE(r, data, elem) \
        BOOST_AUTO_TEST_CASE( BOOST_PP_CAT(c2c_1d_in_place_,elem) ) \
        { c2c_1d_in_place(elem); }
void c2c_1d_in_place(const int N)
{
    typedef boost::multi_array<std::complex<double>,1> array_type;
    array_type both(boost::extents[N]);

    // No dealiasing in effect for in place transform: NR == NC
    check_1D_complex_forward(both, both);
    compare_1D_complex_forward(both, both);
    check_1D_complex_backward(both, both);
    compare_1D_complex_forward(both, both);
}
BOOST_PP_SEQ_FOR_EACH(TEST_C2C_1D_IN_PLACE,_,TRANSFORM_1D_SIZE_SEQ);
BOOST_AUTO_TEST_SUITE_END();

// TODO Out of place transformations on fftw_complex
// TODO In place transformations on fftw_complex
// Broken: boost::multi_array cannot handle fftw_complex elements
// Awaiting response from boost-users mailing list
// http://article.gmane.org/gmane.comp.lib.boost.user/52327
// Helper function that kicks the tires of a 1D c2c transform

BOOST_AUTO_TEST_SUITE( c2c_1d_out_of_place_dealiased );
#define TEST_C2C_1D_OUT_OF_PLACE_DEALIASED(r, product) \
        BOOST_AUTO_TEST_CASE( \
          BOOST_PP_CAT(c2c_1d_out_of_place_dealiased_, \
            BOOST_PP_CAT(BOOST_PP_SEQ_ELEM(0,product), \
              BOOST_PP_CAT(_, \
                BOOST_PP_CAT(BOOST_PP_SEQ_ELEM(1,product), \
                  BOOST_PP_CAT(_,BOOST_PP_SEQ_ELEM(2,product))))))) \
        {\
            BOOST_PP_CAT(c2c_1d_out_of_place_dealiased_, \
                         BOOST_PP_SEQ_ELEM(0,product))( \
                BOOST_PP_SEQ_ELEM(1,product), \
                BOOST_PP_SEQ_ELEM(2,product) ); \
        }
#define TEST_C2C_1D_OUT_OF_PLACE_DEALIASED_LESS_ONLY(r,product) \
    BOOST_PP_IIF(BOOST_PP_LESS(BOOST_PP_SEQ_ELEM(1,product),\
                BOOST_PP_SEQ_ELEM(2,product)), \
                TEST_C2C_1D_OUT_OF_PLACE_DEALIASED(r,product), \
                BOOST_PP_EMPTY());
#define TEST_C2C_1D_OUT_OF_PLACE_DEALIASED_GREATER_ONLY(r,product) \
    BOOST_PP_IIF(BOOST_PP_GREATER(BOOST_PP_SEQ_ELEM(1,product),\
                BOOST_PP_SEQ_ELEM(2,product)), \
                TEST_C2C_1D_OUT_OF_PLACE_DEALIASED(r,product), \
                BOOST_PP_EMPTY());
void c2c_1d_out_of_place_dealiased_forward(const int NR, const int NC)
{
    typedef boost::multi_array<std::complex<double>,1> array_type;
    array_type in(boost::extents[NR]), out(boost::extents[NC]);
    check_1D_complex_forward(in, out); // Dealiasing in effect
}
void c2c_1d_out_of_place_dealiased_backward(const int NC, const int NR)
{
    typedef boost::multi_array<std::complex<double>,1> array_type;
    array_type in(boost::extents[NC]), out(boost::extents[NR]);
    check_1D_complex_backward(in, out); // Dealiasing in effect
}
BOOST_PP_SEQ_FOR_EACH_PRODUCT(\
        TEST_C2C_1D_OUT_OF_PLACE_DEALIASED_LESS_ONLY, \
        ((forward))(TRANSFORM_1D_SIZE_SEQ)(TRANSFORM_1D_SIZE_SEQ) );
BOOST_PP_SEQ_FOR_EACH_PRODUCT(\
        TEST_C2C_1D_OUT_OF_PLACE_DEALIASED_GREATER_ONLY, \
        ((forward))(TRANSFORM_1D_SIZE_SEQ)(TRANSFORM_1D_SIZE_SEQ) );
BOOST_PP_SEQ_FOR_EACH_PRODUCT(\
        TEST_C2C_1D_OUT_OF_PLACE_DEALIASED_LESS_ONLY, \
        ((backward))(TRANSFORM_1D_SIZE_SEQ)(TRANSFORM_1D_SIZE_SEQ) );
BOOST_PP_SEQ_FOR_EACH_PRODUCT(\
        TEST_C2C_1D_OUT_OF_PLACE_DEALIASED_GREATER_ONLY, \
        ((backward))(TRANSFORM_1D_SIZE_SEQ)(TRANSFORM_1D_SIZE_SEQ) );
BOOST_AUTO_TEST_SUITE_END();
