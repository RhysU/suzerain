#include <suzerain/config.h>
#include <suzerain/common.hpp>
#include <suzerain/fftw_multi_array.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <fftw3.h>
#include "test_tools.hpp"

using namespace suzerain;

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

    const double data[M][N] = {{1,2,3}, {4,5,6}, {7,8,9}, {10,11,12}};
    const double close = std::numeric_limits<double>::epsilon()*10*M*M*N*N;

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
                if (fabs(e_real) < close) {
                    BOOST_CHECK_SMALL(real(out[i][j]), close);
                } else {
                    BOOST_CHECK_CLOSE(e_real, real(out[i][j]), close);
                }
                if (fabs(e_imag) < close) {
                    BOOST_CHECK_SMALL(imag(out[i][j]), close);
                } else {
                    BOOST_CHECK_CLOSE(e_imag, imag(out[i][j]), close);
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
                if (fabs(e_real) < close) {
                    BOOST_CHECK_SMALL(real(out[i][j]), close);
                } else {
                    BOOST_CHECK_CLOSE(e_real, real(out[i][j]), close);
                }
                if (fabs(e_imag) < close) {
                    BOOST_CHECK_SMALL(imag(out[i][j]), close);
                } else {
                    BOOST_CHECK_CLOSE(e_imag, imag(out[i][j]), close);
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

BOOST_AUTO_TEST_SUITE( r2c_2d_forward_simple )

// Helper function testing directional transforms for a small 2D grid
template<class RealArray, class ComplexArray>
void r2c_2d_forward_4_by_3(RealArray &in, ComplexArray &out)
{
    using fftw_multi_array::detail::assign_complex;
    using fftw_multi_array::detail::assign_components;

    BOOST_STATIC_ASSERT(RealArray::dimensionality == 2);
    BOOST_STATIC_ASSERT(ComplexArray::dimensionality == 2);
    const int M = 4, N = 3;
    BOOST_REQUIRE_EQUAL(in.shape()[0], M);
    BOOST_REQUIRE_EQUAL(in.shape()[1], N);
    BOOST_REQUIRE_GE(out.shape()[0], std::max(M,M/2+1));
    BOOST_REQUIRE_GE(out.shape()[1], std::max(N,N/2+1));

    const double data[M][N] = {1,2,3, 4,5,6, 7,8,9, 10,11,12};
    const double close = std::numeric_limits<double>::epsilon()*10*M*M*N*N;

    // Need to create views of out to match expected dimensions
    typedef typename ComplexArray::index_range range;
    typedef typename ComplexArray::template array_view<2>::type complex_view;

    { // Transform zeroth dimension and test against expected
        {
            // Deliberately convoluted to track down some valgrind problems
            boost::array<typename RealArray::index,2> idx;
            for (idx[0] = 0; idx[0] < M; ++idx[0]) {
                for (idx[1] = 0; idx[1] < N; ++idx[1]) {
                    typename RealArray::element &val = in(idx);
                    val = data[idx[0]][idx[1]];
                }
            }
        }

        typedef std::complex<double> z;
        const z expected0[M/2+1][N] = {
            z( 5.5, 0.),  z( 6.5, 0.),  z( 7.5, 0.),
            z(-1.5, 1.5), z(-1.5, 1.5), z(-1.5, 1.5),
            z(-1.5, 0.),  z(-1.5, 0.),  z(-1.5, 0.)
        };

        complex_view out_view
            = out[boost::indices[range().finish(M/2+1)][range()]];
        fftw_multi_array::forward_r2c(0, in, out_view);

        for (int i = 0; i < out_view.shape()[0]; ++i) {
            for (int j = 0; j < out_view.shape()[1]; ++j) {
                double e_real, e_imag;
                assign_components(e_real, e_imag, expected0[i][j]);
                if (fabs(e_real) < close) {
                    BOOST_CHECK_SMALL(real(out[i][j]), close);
                } else {
                    BOOST_CHECK_CLOSE(e_real, real(out[i][j]), close);
                }
                if (fabs(e_imag) < close) {
                    BOOST_CHECK_SMALL(imag(out[i][j]), close);
                } else {
                    BOOST_CHECK_CLOSE(e_imag, imag(out[i][j]), close);
                }
            }
        }
    }

    { // Transform first dimension and test against expected
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < N; ++j)
                in[i][j] = data[i][j];

        typedef std::complex<double> z;
        const z expected1[M][N/2+1] = {
            z( 2.,0.), z(-0.5,0.288675134594813),
            z( 5.,0.), z(-0.5,0.288675134594813),
            z( 8.,0.), z(-0.5,0.288675134594813),
            z(11.,0.), z(-0.5,0.288675134594813)
        };

        complex_view out_view
            = out[boost::indices[range()][range().finish(N/2+1)]];
        fftw_multi_array::forward_r2c(1, in, out_view);

        for (int i = 0; i < out_view.shape()[0]; ++i) {
            for (int j = 0; j < out_view.shape()[1]; ++j) {
                double e_real, e_imag;
                assign_components(e_real, e_imag, expected1[i][j]);
                if (fabs(e_real) < close) {
                    BOOST_CHECK_SMALL(real(out_view[i][j]), close);
                } else {
                    BOOST_CHECK_CLOSE(e_real, real(out_view[i][j]), close);
                }
                if (fabs(e_imag) < close) {
                    BOOST_CHECK_SMALL(imag(out_view[i][j]), close);
                } else {
                    BOOST_CHECK_CLOSE(e_imag, imag(out_view[i][j]), close);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( r2c_2d_complex_forward_out_of_place_c_storage )
{
    typedef boost::multi_array<std::complex<double>,2> complex_array;
    typedef boost::multi_array<complex_array::element::value_type,2> real_array;

    real_array    in( boost::extents[4][3], boost::c_storage_order());
    complex_array out(boost::extents[4][3], boost::c_storage_order());
    r2c_2d_forward_4_by_3(in, out);
}

BOOST_AUTO_TEST_CASE( r2c_2d_complex_forward_out_of_place_fortran_storage )
{
    typedef boost::multi_array<std::complex<double>,2> complex_array;
    typedef boost::multi_array<complex_array::element::value_type,2> real_array;

    real_array    in( boost::extents[4][3], boost::fortran_storage_order());
    complex_array out(boost::extents[4][3], boost::fortran_storage_order());
    r2c_2d_forward_4_by_3(in, out);
}

BOOST_AUTO_TEST_CASE( r2c_2d_complex_forward_out_of_place_c2fortran_storage )
{
    typedef boost::multi_array<std::complex<double>,2> complex_array;
    typedef boost::multi_array<complex_array::element::value_type,2> real_array;

    real_array    in( boost::extents[4][3], boost::c_storage_order());
    complex_array out(boost::extents[4][3], boost::fortran_storage_order());
    r2c_2d_forward_4_by_3(in, out);
}

BOOST_AUTO_TEST_CASE( r2c_2d_complex_forward_out_of_place_fortran2c_storage )
{
    typedef boost::multi_array<std::complex<double>,2> complex_array;
    typedef boost::multi_array<complex_array::element::value_type,2> real_array;

    real_array    in( boost::extents[4][3], boost::fortran_storage_order());
    complex_array out(boost::extents[4][3], boost::c_storage_order());
    r2c_2d_forward_4_by_3(in, out);
}

void test_r2c_2d_complex_forward_in_place(const std::size_t (&ordering)[2],
                                          const bool        (&ascending)[2])
{
    typedef std::complex<double> complex;
    const std::size_t size_ratio = sizeof(complex)/sizeof(complex::value_type);

    // Logical two dimensional array size in real space
    const std::size_t M = 4, N = 3;
    // Find maximum required complex extents for transforms in each direction
    const std::size_t complexM = std::max(M, M/2+1);
    const std::size_t complexN = std::max(N, N/2+1);
    // Convert to maximum required real extents sizes
    const std::size_t realM = complexM * size_ratio;
    const std::size_t realN = complexN * size_ratio;
    // Determine amount of complex storage to allocate
    const std::size_t storage_amount
        = std::max(complexM*complexN, realM*realN/size_ratio);

    // Create appropriately typed views of the same raw storage
    boost::scoped_array<complex> raw(new complex[storage_amount]);
    typedef boost::multi_array_ref<complex::value_type,2> real_array;
    typedef boost::multi_array_ref<complex,2>             complex_array;
    typedef boost::general_storage_order<real_array::dimensionality> storage;
    real_array    in( reinterpret_cast<complex::value_type *>(raw.get()),
                      boost::extents[realM][realN],
                      storage(ordering, ascending));
    complex_array out(raw.get(),
                      boost::extents[complexM][complexN],
                      storage(ordering, ascending));

    typedef real_array::index_range range;
    typedef boost::array_view_gen<real_array,2>::type real_view;
    real_view in_view = in[boost::indices[range(0,M)][range(0,N)]];

    r2c_2d_forward_4_by_3(in_view, out);
}

BOOST_AUTO_TEST_CASE( r2c_2d_complex_forward_in_place_c_storage )
{
    const std::size_t ordering[2]  = { 1, 0 };
    const bool        ascending[2] = { true, true };
    test_r2c_2d_complex_forward_in_place(ordering, ascending);
}

BOOST_AUTO_TEST_CASE( r2c_2d_complex_forward_in_place_fortran_storage )
{
    const std::size_t ordering[2]  = { 0, 1 };
    const bool        ascending[2] = { true, true };
    test_r2c_2d_complex_forward_in_place(ordering, ascending);
}

BOOST_AUTO_TEST_CASE( r2c_2d_complex_forward_in_place_c_storage_reversed )
{
    const std::size_t ordering[2]  = { 1, 0 }; // C ordering
    {
        const bool ascending[2] = { true, false };
        test_r2c_2d_complex_forward_in_place(ordering, ascending);
    }
    {
        const bool ascending[2] = { false, true };
        test_r2c_2d_complex_forward_in_place(ordering, ascending);
    }
    {
        const bool ascending[2] = { false, false };
        test_r2c_2d_complex_forward_in_place(ordering, ascending);
    }
}

BOOST_AUTO_TEST_CASE( r2c_2d_complex_forward_in_place_fortran_storage_reversed )
{
    const std::size_t ordering[2]  = { 0, 1 }; // Fortran ordering
    {
        const bool ascending[2] = { true, false };
        test_r2c_2d_complex_forward_in_place(ordering, ascending);
    }
    {
        const bool ascending[2] = { false, true };
        test_r2c_2d_complex_forward_in_place(ordering, ascending);
    }
    {
        const bool ascending[2] = { false, false };
        test_r2c_2d_complex_forward_in_place(ordering, ascending);
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( c2r_2d_backward_simple )

// Helper function testing directional transforms for a small 2D grid
template<class ComplexArray, class RealArray>
void c2r_2d_backward_4_by_3(ComplexArray &in, RealArray &out)
{
    BOOST_STATIC_ASSERT(ComplexArray::dimensionality == 2);
    BOOST_STATIC_ASSERT(RealArray::dimensionality == 2);
    const int M = 4, N = 3;
    BOOST_REQUIRE_GE(in.shape()[0], std::max(M,M/2+1));
    BOOST_REQUIRE_GE(in.shape()[1], std::max(N,N/2+1));
    BOOST_REQUIRE_EQUAL(out.shape()[0], M);
    BOOST_REQUIRE_EQUAL(out.shape()[1], N);

    const double expected[M][N] = {1,2,3, 4,5,6, 7,8,9, 10,11,12};
    const double close = std::numeric_limits<double>::epsilon()*10*M*M*N*N;

    // Need to create views of in to match expected dimensions
    typedef typename ComplexArray::index index;
    typedef typename ComplexArray::index_range range;
    typedef typename ComplexArray::template array_view<2>::type complex_view;

    { // Transform zeroth dimension and test against expected
        typedef std::complex<double> z;
        const z data[M/2+1][N] = {
            z( 5.5, 0.),  z( 6.5, 0.),  z( 7.5, 0.),
            z(-1.5, 1.5), z(-1.5, 1.5), z(-1.5, 1.5),
            z(-1.5, 0.),  z(-1.5, 0.),  z(-1.5, 0.)
        };
        complex_view in_view
            = in[boost::indices[range().finish(M/2+1)][range()]];
        for (index i = 0; i < in_view.shape()[0]; ++i)
            for (index j = 0; j < in_view.shape()[1]; ++j)
                in_view[i][j] = data[i][j];

        fftw_multi_array::backward_c2r(0, in_view, out);

        for (int i = 0; i < out.shape()[0]; ++i) {
            for (int j = 0; j < out.shape()[1]; ++j) {
                if (fabs(expected[i][j]) < close) {
                    BOOST_CHECK_SMALL(out[i][j], close);
                } else {
                    BOOST_CHECK_CLOSE(expected[i][j], out[i][j], close);
                }
            }
        }
    }

    { // Transform first dimension and test against expected
        typedef std::complex<double> z;
        const z data[M][N/2+1] = {
            z( 2.,0.), z(-0.5,0.288675134594813),
            z( 5.,0.), z(-0.5,0.288675134594813),
            z( 8.,0.), z(-0.5,0.288675134594813),
            z(11.,0.), z(-0.5,0.288675134594813)
        };
        complex_view in_view
            = in[boost::indices[range()][range().finish(N/2+1)]];
        for (index i = 0; i < in_view.shape()[0]; ++i)
            for (index j = 0; j < in_view.shape()[1]; ++j)
                in_view[i][j] = data[i][j];

        fftw_multi_array::backward_c2r(1, in_view, out);

        for (int i = 0; i < out.shape()[0]; ++i) {
            for (int j = 0; j < out.shape()[1]; ++j) {
                if (fabs(expected[i][j]) < close) {
                    BOOST_CHECK_SMALL(out[i][j], close);
                } else {
                    BOOST_CHECK_CLOSE(expected[i][j], out[i][j], close);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( c2r_2d_complex_backward_out_of_place_c_storage )
{
    typedef boost::multi_array<std::complex<double>,2> complex_array;
    typedef boost::multi_array<complex_array::element::value_type,2> real_array;

    complex_array in( boost::extents[4][3], boost::c_storage_order());
    real_array    out(boost::extents[4][3], boost::c_storage_order());
    c2r_2d_backward_4_by_3(in, out);
}

BOOST_AUTO_TEST_CASE( c2r_2d_complex_backward_out_of_place_fortran_storage )
{
    typedef boost::multi_array<std::complex<double>,2> complex_array;
    typedef boost::multi_array<complex_array::element::value_type,2> real_array;

    complex_array in( boost::extents[4][3], boost::fortran_storage_order());
    real_array    out(boost::extents[4][3], boost::fortran_storage_order());
    c2r_2d_backward_4_by_3(in, out);
}

BOOST_AUTO_TEST_CASE( c2r_2d_complex_backward_out_of_place_c2fortran_storage )
{
    typedef boost::multi_array<std::complex<double>,2> complex_array;
    typedef boost::multi_array<complex_array::element::value_type,2> real_array;

    complex_array in( boost::extents[4][3], boost::fortran_storage_order());
    real_array    out(boost::extents[4][3], boost::c_storage_order());
    c2r_2d_backward_4_by_3(in, out);
}

BOOST_AUTO_TEST_CASE( c2r_2d_complex_backward_out_of_place_fortran2c_storage )
{
    typedef boost::multi_array<std::complex<double>,2> complex_array;
    typedef boost::multi_array<complex_array::element::value_type,2> real_array;

    complex_array in( boost::extents[4][3], boost::c_storage_order());
    real_array    out(boost::extents[4][3], boost::fortran_storage_order());
    c2r_2d_backward_4_by_3(in, out);
}

void test_c2r_2d_complex_backward_in_place(const std::size_t (&ordering)[2],
                                           const bool        (&ascending)[2])
{
    typedef std::complex<double> complex;
    const std::size_t size_ratio = sizeof(complex)/sizeof(complex::value_type);

    // Logical two dimensional array size in real space
    const std::size_t M = 4, N = 3;
    // Find maximum required complex extents for transforms in each direction
    const std::size_t complexM = std::max(M, M/2+1);
    const std::size_t complexN = std::max(N, N/2+1);
    // Convert to maximum required real extents sizes
    const std::size_t realM = complexM * size_ratio;
    const std::size_t realN = complexN * size_ratio;
    // Determine amount of complex storage to allocate
    const std::size_t storage_amount
        = std::max(complexM*complexN, realM*realN/size_ratio);

    // Create appropriately typed views of the same raw storage
    boost::scoped_array<complex> raw(new complex[storage_amount]);
    typedef boost::multi_array_ref<complex,2>             complex_array;
    typedef boost::multi_array_ref<complex::value_type,2> real_array;
    typedef boost::general_storage_order<real_array::dimensionality> storage;
    complex_array in( raw.get(),
                      boost::extents[complexM][complexN],
                      storage(ordering, ascending));
    real_array    out(reinterpret_cast<complex::value_type *>(raw.get()),
                      boost::extents[realM][realN],
                      storage(ordering, ascending));

    typedef real_array::index_range range;
    typedef boost::array_view_gen<real_array,2>::type real_view;
    real_view out_view = out[boost::indices[range(0,M)][range(0,N)]];

    c2r_2d_backward_4_by_3(in, out_view);
}

BOOST_AUTO_TEST_CASE( c2r_2d_complex_backward_in_place_c_storage )
{
    const std::size_t ordering[2]  = { 1, 0 };
    const bool        ascending[2] = { true, true };
    test_c2r_2d_complex_backward_in_place(ordering, ascending);
}

BOOST_AUTO_TEST_CASE( c2r_2d_complex_backward_in_place_fortran_storage )
{
    const std::size_t ordering[2]  = { 0, 1 };
    const bool        ascending[2] = { true, true };
    test_c2r_2d_complex_backward_in_place(ordering, ascending);
}

BOOST_AUTO_TEST_CASE( c2r_2d_complex_backward_in_place_c_storage_reversed )
{
    const std::size_t ordering[2]  = { 1, 0 }; // C ordering
    {
        const bool ascending[2] = { true, false };
        test_c2r_2d_complex_backward_in_place(ordering, ascending);
    }
    {
        const bool ascending[2] = { false, true };
        test_c2r_2d_complex_backward_in_place(ordering, ascending);
    }
    {
        const bool ascending[2] = { false, false };
        test_c2r_2d_complex_backward_in_place(ordering, ascending);
    }
}

BOOST_AUTO_TEST_CASE( c2r_2d_complex_backward_in_place_fortran_storage_reversed )
{
    const std::size_t ordering[2]  = { 0, 1 }; // Fortran ordering
    {
        const bool ascending[2] = { true, false };
        test_c2r_2d_complex_backward_in_place(ordering, ascending);
    }
    {
        const bool ascending[2] = { false, true };
        test_c2r_2d_complex_backward_in_place(ordering, ascending);
    }
    {
        const bool ascending[2] = { false, false };
        test_c2r_2d_complex_backward_in_place(ordering, ascending);
    }
}

BOOST_AUTO_TEST_SUITE_END()
