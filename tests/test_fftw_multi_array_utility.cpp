#include <suzerain/config.h>
#include <suzerain/common.hpp>
#include <suzerain/fftw_multi_array.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/assign.hpp>
#include <boost/test/included/unit_test.hpp>
#include <fftw3.h>
#include "test_tools.hpp"

using namespace suzerain;

BOOST_AUTO_TEST_SUITE( integer_power )
BOOST_AUTO_TEST_CASE( detail_integer_power )
{
    using fftw_multi_array::detail::integer_power;

    BOOST_CHECK_EQUAL( 1u, integer_power(2u, 0u));
    BOOST_CHECK_EQUAL( 2u, integer_power(2u, 1u));
    BOOST_CHECK_EQUAL( 4u, integer_power(2u, 2u));
    BOOST_CHECK_EQUAL( 8u, integer_power(2u, 3u));
    BOOST_CHECK_EQUAL(16u, integer_power(2u, 4u));

    BOOST_CHECK_EQUAL( 1ul, integer_power(2ul, 0ul));
    BOOST_CHECK_EQUAL( 2ul, integer_power(2ul, 1ul));
    BOOST_CHECK_EQUAL( 4ul, integer_power(2ul, 2ul));
    BOOST_CHECK_EQUAL( 8ul, integer_power(2ul, 3ul));
    BOOST_CHECK_EQUAL(16ul, integer_power(2ul, 4ul));

    BOOST_CHECK_EQUAL( 1.0, integer_power(2.0, 0));
    BOOST_CHECK_EQUAL( 2.0, integer_power(2.0, 1));
    BOOST_CHECK_EQUAL( 4.0, integer_power(2.0, 2));
    BOOST_CHECK_EQUAL( 8.0, integer_power(2.0, 3));
    BOOST_CHECK_EQUAL(16.0, integer_power(2.0, 4));

    BOOST_CHECK_EQUAL( 1.0,    integer_power(2.0, -0));
    BOOST_CHECK_EQUAL( 0.5,    integer_power(2.0, -1));
    BOOST_CHECK_EQUAL( 0.25,   integer_power(2.0, -2));
    BOOST_CHECK_EQUAL( 0.125,  integer_power(2.0, -3));
    BOOST_CHECK_EQUAL( 0.0625, integer_power(2.0, -4));

    BOOST_CHECK_EQUAL( 1.0, integer_power(0.0, 0));
    BOOST_CHECK_EQUAL( 0.0, integer_power(0.0, 1));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( complex_helpers )

BOOST_AUTO_TEST_CASE( assign_complex )
{
    using fftw_multi_array::detail::assign_complex;

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

BOOST_AUTO_TEST_CASE( assign_complex_scaled )
{
    using fftw_multi_array::detail::assign_complex_scaled;

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
    using fftw_multi_array::detail::assign_complex_scaled_ipower;

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
    using fftw_multi_array::detail::assign_components;

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

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( increment )

BOOST_AUTO_TEST_CASE( increment_1d_degenerate )
{
    using fftw_multi_array::detail::increment;

    const int n        = 1;
    int       index[n] = { 0 };
    const int shape[n] = { 1 };

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_1d_normal )
{
    using fftw_multi_array::detail::increment;

    const int n        = 1;
    int       index[n] = { 0 };
    const int shape[n] = { 3 };

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_1d_normal_unsigned )
{
    using fftw_multi_array::detail::increment;

    const unsigned int n        = 1;
    unsigned int       index[n] = { 0 };
    const unsigned int shape[n] = { 3 };

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_degenerate_first )
{
    using fftw_multi_array::detail::increment;

    const int           n     = 2;
    boost::array<int,n> index_array = {{ 0, 0 }};
    boost::array<int,n> shape_array = {{ 3, 1 }};
    boost::array<int,n>::iterator index = index_array.begin();
    boost::array<int,n>::iterator shape = shape_array.begin();

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index_array, boost::assign::list_of(1)(0));

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index_array, boost::assign::list_of(2)(0));

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_degenerate_second )
{
    using fftw_multi_array::detail::increment;

    const int                   n     = 2;
    boost::array<int, n> index_array            = {{ 0, 0 }};
    boost::array<std::size_t,n> shape_array     = {{ 1, 3 }};
    boost::array<int, n>::iterator index        = index_array.begin();
    boost::array<std::size_t,n>::iterator shape = shape_array.begin();

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index_array, boost::assign::list_of(0)(1));

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index_array, boost::assign::list_of(0)(2));

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_degenerate_all )
{
    using fftw_multi_array::detail::increment;

    const int  n        = 2;
    int        index[n] = { 0, 0 };
    const long shape[n] = { 1, 1 };

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_normal )
{
    using fftw_multi_array::detail::increment;

    const int    n        = 2;
    unsigned int index[n] = { 0, 0 };
    unsigned int shape[n] = { 2, 3 };

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_normal_usualorder )
{
    using fftw_multi_array::detail::increment;

    const int    n        = 2;
    unsigned int index[n] = { 0, 0 };
    unsigned int shape[n] = { 2, 3 };
    unsigned int order[n] = { 0, 1 };

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_normal_reverseorder )
{
    using fftw_multi_array::detail::increment;

    const int    n        = 2;
    unsigned int index[n] = { 0, 0 };
    unsigned int shape[n] = { 2, 3 };
    unsigned int order[n] = { 1, 0 };

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), false);
}

BOOST_AUTO_TEST_CASE( increment_3d_degenerate_all )
{
    using fftw_multi_array::detail::increment;

    const int         n        = 3;
    signed int        index[n] = { 0, 0, 0 };
    const signed long shape[n] = { 1, 1, 1 };

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_3d_degenerate_middle )
{
    using fftw_multi_array::detail::increment;

    const int                   n               = 3;
    long                        index[n]        =  { 0, 0, 0 };
    boost::array<std::size_t,n> shape_array     = {{ 3, 1, 3 }};
    boost::array<std::size_t,n>::iterator shape = shape_array.begin();

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 2);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 2);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 2);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_3d_normal )
{
    using fftw_multi_array::detail::increment;

    const int n = 3;
    std::vector<short> index_array(n, 0);
    std::vector<int>   shape_array(n, 2);
    std::vector<short>::iterator     index = index_array.begin();
    std::vector<int>::const_iterator shape = shape_array.begin();

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_3d_normal_outoforder )
{
    using fftw_multi_array::detail::increment;

    const int n = 3;
    boost::array<short,3> index_array = {{ 0, 0, 0 }};
    boost::array<int,3>   shape_array = {{ 2, 1, 3 }};
    boost::array<long,3>  order_array = {{ 2, 0, 1 }};

    boost::array<short,3>::iterator index = index_array.begin();
    boost::array<int,3>::iterator   shape = shape_array.begin();
    boost::array<long,3>::iterator  order = order_array.begin();

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape,order), false);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( decrement )

BOOST_AUTO_TEST_CASE( decrement_1d_degenerate )
{
    using fftw_multi_array::detail::decrement;

    const int n        = 1;
    int       index[n] = { 0 };
    const int shape[n] = { 1 };

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_1d_normal )
{
    using fftw_multi_array::detail::decrement;

    const int n        = 1;
    int       index[n] = { 2 };
    const int shape[n] = { 3 };

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_1d_normal_unsigned )
{
    using fftw_multi_array::detail::decrement;

    const unsigned int n        = 1;
    unsigned int       index[n] = { 2 };
    const unsigned int shape[n] = { 3 };

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_2d_degenerate_first )
{
    using fftw_multi_array::detail::decrement;

    const int n = 2;
    boost::array<int,n> index_array = {{ 2, 0 }};
    boost::array<int,n> shape_array = {{ 3, 1 }};
    boost::array<int,n>::iterator index = index_array.begin();
    boost::array<int,n>::iterator shape = shape_array.begin();

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_2d_degenerate_second )
{
    using fftw_multi_array::detail::decrement;

    const int n = 2;
    boost::array<int, n> index_array            = {{ 0, 2 }};
    boost::array<std::size_t,n> shape_array     = {{ 1, 3 }};
    boost::array<int, n>::iterator index        = index_array.begin();
    boost::array<std::size_t,n>::iterator shape = shape_array.begin();

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_2d_degenerate_all )
{
    using fftw_multi_array::detail::decrement;

    const int  n        = 2;
    int        index[n] = { 0, 0 };
    const long shape[n] = { 1, 1 };

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_2d_normal )
{
    using fftw_multi_array::detail::decrement;

    const int    n        = 2;
    unsigned int index[n] = { 1, 2 };
    unsigned int shape[n] = { 2, 3 };

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_2d_normal_usualorder )
{
    using fftw_multi_array::detail::decrement;

    const int    n        = 2;
    unsigned int index[n] = { 1, 2 };
    unsigned int shape[n] = { 2, 3 };
    unsigned int order[n] = { 0, 1 };

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), false);
}

BOOST_AUTO_TEST_CASE( decrement_2d_normal_reverseorder )
{
    using fftw_multi_array::detail::decrement;

    const int    n        = 2;
    unsigned int index[n] = { 1, 2 };
    unsigned int shape[n] = { 2, 3 };
    unsigned int order[n] = { 1, 0 };

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 2);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), false);
}

BOOST_AUTO_TEST_CASE( decrement_3d_degenerate_all )
{
    using fftw_multi_array::detail::decrement;

    const int         n        = 3;
    signed int        index[n] = { 0, 0, 0 };
    const signed long shape[n] = { 1, 1, 1 };

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_3d_degenerate_middle )
{
    using fftw_multi_array::detail::decrement;

    const int                   n               = 3;
    long                        index[n]        =  { 2, 0, 2 };
    boost::array<std::size_t,n> shape_array     = {{ 3, 1, 3 }};
    boost::array<std::size_t,n>::iterator shape = shape_array.begin();

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 2);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 2);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_3d_normal )
{
    using fftw_multi_array::detail::decrement;

    const int n = 3;
    std::vector<short> index_array(n, 1);
    std::vector<int>   shape_array(n, 2);
    std::vector<short>::iterator     index = index_array.begin();
    std::vector<int>::const_iterator shape = shape_array.begin();

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( decrement_3d_normal_outoforder )
{
    using fftw_multi_array::detail::decrement;

    const int n = 3;
    boost::array<short,3> index_array = {{ 1, 0, 2 }};
    boost::array<int,3>   shape_array = {{ 2, 1, 3 }};
    boost::array<long,3>  order_array = {{ 2, 0, 1 }};

    boost::array<short,3>::iterator index = index_array.begin();
    boost::array<int,3>::iterator   shape = shape_array.begin();
    boost::array<long,3>::iterator  order = order_array.begin();

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 2);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 1);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(index[2], 0);
    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), false);
}

BOOST_AUTO_TEST_CASE( decrement_3d_spot_check_behavior_1 )
{
    using fftw_multi_array::detail::decrement;

    const int n = 3;
    typedef boost::array<int,n> array_type;
    array_type           index_array = {{ 3, 0, 1 }};
    array_type::iterator index       = index_array.begin();
    array_type           shape_array = {{ 4, 1, 2 }};
    array_type::iterator shape       = shape_array.begin();
    array_type           order_array = {{ 0, 1, 2 }};
    array_type::iterator order       = order_array.begin();

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index_array, boost::assign::list_of(2)(0)(1));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index_array, boost::assign::list_of(1)(0)(1));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index_array, boost::assign::list_of(0)(0)(1));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index_array, boost::assign::list_of(3)(0)(0));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index_array, boost::assign::list_of(2)(0)(0));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index_array, boost::assign::list_of(1)(0)(0));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index_array, boost::assign::list_of(0)(0)(0));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), false);
}

BOOST_AUTO_TEST_CASE( decrement_3d_spot_check_behavior_2 )
{
    using fftw_multi_array::detail::decrement;

    const int n = 3;
    typedef boost::array<int,n> array_type;
    array_type           index_array = {{ 3, 0, 1 }};
    array_type::iterator index       = index_array.begin();
    array_type           shape_array = {{ 4, 1, 2 }};
    array_type::iterator shape       = shape_array.begin();
    array_type           order_array = {{ 2, 1, 0 }};
    array_type::iterator order       = order_array.begin();

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index_array, boost::assign::list_of(3)(0)(0));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index_array, boost::assign::list_of(2)(0)(1));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index_array, boost::assign::list_of(2)(0)(0));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index_array, boost::assign::list_of(1)(0)(1));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index_array, boost::assign::list_of(1)(0)(0));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index_array, boost::assign::list_of(0)(0)(1));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), true);
    BOOST_REQUIRE_EQUAL(index_array, boost::assign::list_of(0)(0)(0));

    BOOST_REQUIRE_EQUAL(decrement<n>(index,shape,order), false);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( indexed_element_magnitude_comparator )

BOOST_AUTO_TEST_CASE( indexed_element_magnitude_comparator_operator_positive )
{
    using fftw_multi_array::detail::indexed_element_magnitude_comparator;

    const double data[3] = { 0.0, 1.0, 2.0 };
    indexed_element_magnitude_comparator<const double *> c(&data[0]);

    BOOST_CHECK_EQUAL(true, c(0, 1));
    BOOST_CHECK_EQUAL(true, c(1, 2));
    BOOST_CHECK_EQUAL(true, c(0, 2));

    BOOST_CHECK_EQUAL(false, c(1, 0));
    BOOST_CHECK_EQUAL(false, c(2, 1));
    BOOST_CHECK_EQUAL(false, c(2, 0));
}

BOOST_AUTO_TEST_CASE( indexed_element_magnitude_comparator_operator_mixed )
{
    using fftw_multi_array::detail::indexed_element_magnitude_comparator;

    const double data[3] = { 0.0, -1.0, -2.0 };
    indexed_element_magnitude_comparator<const double *> c(&data[0]);

    BOOST_CHECK_EQUAL(true, c(0, 1));
    BOOST_CHECK_EQUAL(true, c(1, 2));
    BOOST_CHECK_EQUAL(true, c(0, 2));

    BOOST_CHECK_EQUAL(false, c(1, 0));
    BOOST_CHECK_EQUAL(false, c(2, 1));
    BOOST_CHECK_EQUAL(false, c(2, 0));
}

BOOST_AUTO_TEST_CASE( indexed_element_magnitude_comparator_stable_sort )
{
    using fftw_multi_array::detail::make_indexed_element_magnitude_comparator;

    const int n = 5;
    const double data[n]    = { 0, 1, 2, 3, 4 };
    std::size_t  indices[n] = { 3, 1, 4, 0, 2 };

    std::stable_sort(indices, indices + n,
            make_indexed_element_magnitude_comparator(&data[0]));
    BOOST_CHECK_EQUAL(0, indices[0]);
    BOOST_CHECK_EQUAL(1, indices[1]);
    BOOST_CHECK_EQUAL(2, indices[2]);
    BOOST_CHECK_EQUAL(3, indices[3]);
    BOOST_CHECK_EQUAL(4, indices[4]);
}

BOOST_AUTO_TEST_SUITE_END()
