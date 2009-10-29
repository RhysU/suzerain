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

BOOST_AUTO_TEST_SUITE( c2c_2d_forward_simple )

// Helper function testing directional transforms for a small 2D grid
template<class ComplexMultiArray1, class ComplexMultiArray2>
void c2c_2d_forward_4_by_3(ComplexMultiArray1 &in, ComplexMultiArray2 &out)
{
    using fftw_multi_array::detail::assign_complex;
    using fftw_multi_array::detail::assign_components;

    BOOST_STATIC_ASSERT(ComplexMultiArray1::dimensionality == 2);
    BOOST_STATIC_ASSERT(ComplexMultiArray2::dimensionality == 2);
    const int M = 4, N = 3;
    BOOST_REQUIRE_EQUAL(in.shape()[0], M);
    BOOST_REQUIRE_EQUAL(in.shape()[1], N);
    BOOST_REQUIRE_EQUAL(out.shape()[0], M);
    BOOST_REQUIRE_EQUAL(out.shape()[1], N);

    const double data[M][N] = {1,2,3, 4,5,6, 7,8,9, 10,11,12};
    const double close_enough
        = std::numeric_limits<double>::epsilon()*10*M*M*N*N;


    { // Transform zeroth dimension and test against expected
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < N; ++j)
                assign_complex(in[i][j], data[i][j], 0);

        typedef std::complex<double> z;
        const z expected0[M][N] = {
            z( 5.5, 0.),  z( 6.5, 0.),  z( 7.5, 0.),
            z(-1.5, 1.5), z(-1.5, 1.5), z(-1.5, 1.5),
            z(-1.5, 0.),  z(-1.5, 0.),  z(-1.5, 0.),
            z(-1.5,-1.5), z(-1.5,-1.5), z(-1.5,-1.5)
        };
        fftw_multi_array::forward_c2c(0, in, out);
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                double e_real, e_imag;
                assign_components(e_real, e_imag, expected0[i][j]);
                if (fabs(e_real) < close_enough) {
                    BOOST_CHECK_SMALL(real(out[i][j]), close_enough);
                } else {
                    BOOST_CHECK_CLOSE(e_real, real(out[i][j]), close_enough);
                }
                if (fabs(e_imag) < close_enough) {
                    BOOST_CHECK_SMALL(imag(out[i][j]), close_enough);
                } else {
                    BOOST_CHECK_CLOSE(e_imag, imag(out[i][j]), close_enough);
                }
            }
        }
    }

    { // Transform first dimension and test against expected
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < N; ++j)
                assign_complex(in[i][j], data[i][j], 0);

        typedef std::complex<double> z;
        const z expected1[M][N] = {
            z( 2.,0.), z(-0.5,0.288675134594813), z(-0.5,-0.288675134594813),
            z( 5.,0.), z(-0.5,0.288675134594813), z(-0.5,-0.288675134594813),
            z( 8.,0.), z(-0.5,0.288675134594813), z(-0.5,-0.288675134594813),
            z(11.,0.), z(-0.5,0.288675134594813), z(-0.5,-0.288675134594813)
        };
        fftw_multi_array::forward_c2c(1, in, out);
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                double e_real, e_imag;
                assign_components(e_real, e_imag, expected1[i][j]);
                if (fabs(e_real) < close_enough) {
                    BOOST_CHECK_SMALL(real(out[i][j]), close_enough);
                } else {
                    BOOST_CHECK_CLOSE(e_real, real(out[i][j]), close_enough);
                }
                if (fabs(e_imag) < close_enough) {
                    BOOST_CHECK_SMALL(imag(out[i][j]), close_enough);
                } else {
                    BOOST_CHECK_CLOSE(e_imag, imag(out[i][j]), close_enough);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( c2c_2d_complex_forward_out_of_place_c_storage )
{
    typedef boost::multi_array<std::complex<double>,2> array_type;

    array_type in(boost::extents[4][3],  boost::c_storage_order());
    array_type out(boost::extents[4][3], boost::c_storage_order());
    c2c_2d_forward_4_by_3(in, out);
}

BOOST_AUTO_TEST_CASE( c2c_2d_complex_forward_out_of_place_fortran_storage )
{
    typedef boost::multi_array<std::complex<double>,2> array_type;

    array_type in(boost::extents[4][3],  boost::fortran_storage_order());
    array_type out(boost::extents[4][3], boost::fortran_storage_order());
    c2c_2d_forward_4_by_3(in, out);
}

BOOST_AUTO_TEST_CASE( c2c_2d_complex_forward_out_of_place_c2fortran_storage )
{
    typedef boost::multi_array<std::complex<double>,2> array_type;

    array_type in(boost::extents[4][3],  boost::c_storage_order());
    array_type out(boost::extents[4][3], boost::fortran_storage_order());
    c2c_2d_forward_4_by_3(in, out);
}

BOOST_AUTO_TEST_CASE( c2c_2d_complex_forward_out_of_place_fortran2c_storage )
{
    typedef boost::multi_array<std::complex<double>,2> array_type;

    array_type in(boost::extents[4][3],  boost::fortran_storage_order());
    array_type out(boost::extents[4][3], boost::c_storage_order());
    c2c_2d_forward_4_by_3(in, out);
}

BOOST_AUTO_TEST_CASE( c2c_2d_complex_forward_in_place_c_storage )
{
    typedef boost::multi_array<std::complex<double>,2> array_type;

    array_type both(boost::extents[4][3], boost::c_storage_order());
    c2c_2d_forward_4_by_3(both, both);
}

BOOST_AUTO_TEST_CASE( c2c_2d_complex_forward_in_place_fortran_storage )
{
    typedef boost::multi_array<std::complex<double>,2> array_type;

    array_type both(boost::extents[4][3], boost::fortran_storage_order());
    c2c_2d_forward_4_by_3(both, both);
}

BOOST_AUTO_TEST_SUITE_END()

