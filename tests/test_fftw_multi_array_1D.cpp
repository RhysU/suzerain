#include <suzerain/config.h>
#include <suzerain/common.hpp>
#include <suzerain/fftw_multi_array.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>
#include <fftw3.h>
#include "test_tools.hpp"

using namespace pecos::suzerain;

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

BOOST_AUTO_TEST_CASE( increment_2d_degenerate_first )
{
    using fftw_multi_array::detail::increment;

    const int           n     = 2;
    boost::array<int,n> index_array = {{ 0, 0 }};
    boost::array<int,n> shape_array = {{ 3, 1 }};
    boost::array<int,n>::iterator index = index_array.begin();
    boost::array<int,n>::iterator shape = shape_array.begin();

    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 1);
    BOOST_REQUIRE_EQUAL(index[1], 0);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 2);
    BOOST_REQUIRE_EQUAL(index[1], 0);
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
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 1);
    BOOST_REQUIRE_EQUAL(increment<n>(index,shape), true);
    BOOST_REQUIRE_EQUAL(index[0], 0);
    BOOST_REQUIRE_EQUAL(index[1], 2);
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

BOOST_AUTO_TEST_SUITE( indexed_element_comparator )

BOOST_AUTO_TEST_CASE( indexed_element_comparator_operator )
{
    using fftw_multi_array::detail::indexed_element_comparator;

    const double data[3] = { 0.0, 1.0, 2.0 };
    indexed_element_comparator<const double *> c(&data[0]);

    BOOST_CHECK_EQUAL(true, c(0, 1));
    BOOST_CHECK_EQUAL(true, c(1, 2));
    BOOST_CHECK_EQUAL(true, c(0, 2));

    BOOST_CHECK_EQUAL(false, c(1, 0));
    BOOST_CHECK_EQUAL(false, c(2, 1));
    BOOST_CHECK_EQUAL(false, c(2, 0));
}

BOOST_AUTO_TEST_CASE( indexed_element_comparator_stable_sort )
{
    using fftw_multi_array::detail::make_indexed_element_comparator;

    const int n = 5;
    const double data[n]    = { 0, 1, 2, 3, 4 };
    std::size_t  indices[n] = { 3, 1, 4, 0, 2 };

    std::stable_sort(indices, indices + n,
            make_indexed_element_comparator(&data[0]));
    BOOST_CHECK_EQUAL(0, indices[0]);
    BOOST_CHECK_EQUAL(1, indices[1]);
    BOOST_CHECK_EQUAL(2, indices[2]);
    BOOST_CHECK_EQUAL(3, indices[3]);
    BOOST_CHECK_EQUAL(4, indices[4]);
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

// Produce a periodic real signal with known frequency content
// on domain of supplied length
template<typename FPT, typename Integer>
FPT real_test_function(const Integer NR,
                       const Integer max_mode_exclusive,
                       const Integer i,
                       const FPT shift = M_PI/3.0,
                       const FPT length = 2.0*M_PI,
                       const Integer derivative = 0) {
    const FPT xi = i*length/NR;
    FPT retval = (max_mode_exclusive > 0 && derivative == 0 )
        ? 17.0 : 0; // Zero mode fixed
    for (Integer i = 1; i < max_mode_exclusive; ++i) {
        switch (derivative % 4) {
            case 0:
                retval +=   i * pow(i*(2.0*M_PI/length), derivative)
                          * sin(i*(2.0*M_PI/length)*xi + shift);
                break;
            case 1:
                retval +=   i * pow(i*(2.0*M_PI/length), derivative)
                          * cos(i*(2.0*M_PI/length)*xi + shift);
                break;
            case 2:
                retval -=   i * pow(i*(2.0*M_PI/length), derivative)
                          * sin(i*(2.0*M_PI/length)*xi + shift);
                break;
            case 3:
                retval -=   i * pow(i*(2.0*M_PI/length), derivative)
                          * cos(i*(2.0*M_PI/length)*xi + shift);
                break;
            default:
                BOOST_ERROR("Unexpected derivative % 4 result");
        }
    }
    return retval;
}

template<class ComplexMultiArray>
void fill_with_complex_NaN(ComplexMultiArray &x)
{
    typedef typename ComplexMultiArray::element element;
    element NaN_value;
    typedef typename fftw_multi_array::detail::transform_traits<
            element>::real_type real_type;
    const real_type quiet_NaN = std::numeric_limits<real_type>::quiet_NaN();

    const element * const end = x.data() + x.num_elements();
    element *it = x.data();
    while (it != end) {
        fftw_multi_array::detail::assign_complex(*it++, quiet_NaN, quiet_NaN);
    }
}

template<class ComplexMultiArray>
void fill_with_complex_zero(ComplexMultiArray &x)
{
    typedef typename ComplexMultiArray::element element;

    const element * const end = x.data() + x.num_elements();
    element *it = x.data();
    while (it != end) {
        fftw_multi_array::detail::assign_complex(*it++, 0, 0);
    }
}

// Helper function that kicks the tires of a 1D c2c transform
template<class ComplexMultiArray1, class ComplexMultiArray2>
void symmetry_1D_complex_forward(ComplexMultiArray1 &in, ComplexMultiArray2 &out)
{
    BOOST_STATIC_ASSERT(ComplexMultiArray1::dimensionality == 1);
    BOOST_STATIC_ASSERT(ComplexMultiArray2::dimensionality == 1);
    typedef typename fftw_multi_array::detail::transform_traits<
        typename ComplexMultiArray1::element>::real_type real_type;
    const int NR = in.shape()[0];
    const int NC = out.shape()[0];
    const real_type close_enough
        = std::numeric_limits<real_type>::epsilon()*10*NR*NR;
    const real_type shift = M_PI/3.0;

    // Load a real-valued function into the input array and transform it
    fill_with_complex_NaN(in);
    fill_with_complex_NaN(out);
    for (int i = 0; i < NR; ++i) {
        fftw_multi_array::detail::assign_complex(in[i],
                real_test_function<real_type>(NR, (NR+1)/2, i, shift), 0.0);
    }
    fftw_multi_array::forward_c2c(0, in, out);

    // Real input should exhibit conjugate symmetry in wave space...
    for (int i = 1; i < (std::min(NC,NR)+1)/2; ++i) { // ...up to grid modes
        real_type a_real, a_imag, b_real, b_imag;
        fftw_multi_array::detail::assign_components(a_real, a_imag, out[i]);
        fftw_multi_array::detail::assign_components(b_real, b_imag, out[NC-i]);
        BOOST_REQUIRE_CLOSE(a_real,  b_real, close_enough);
        BOOST_REQUIRE_CLOSE(a_imag, -b_imag, close_enough);

        // We should also see the expected frequency content
        const real_type real_expected = i * sin(shift) * 1.0/2.0;
        const real_type imag_expected = i * cos(shift) * 1.0/2.0;
        BOOST_REQUIRE_CLOSE(b_real, real_expected, sqrt(close_enough));
        BOOST_REQUIRE_CLOSE(b_imag, imag_expected, sqrt(close_enough));
    }
    // Ensure we see the expected zero mode magnitude
    {
        real_type z_real, z_imag;
        fftw_multi_array::detail::assign_components(z_real, z_imag, out[0]);
        BOOST_REQUIRE_SMALL(z_imag, close_enough);
        BOOST_REQUIRE_CLOSE(z_real, 17.0, close_enough);
    }

    // Load an imaginary-valued function into the input array and transform it
    fill_with_complex_NaN(in);
    fill_with_complex_NaN(out);
    for (int i = 0; i < NR; ++i) {
        fftw_multi_array::detail::assign_complex(in[i], 0.0,
                real_test_function<real_type>(NR, (NR+1)/2, i, shift));
    }
    fftw_multi_array::forward_c2c(0, in, out);

    // Imaginary input should exhibit a similar symmetry in wave space...
    // Re(X_k) = - Re(X_{N-k}), Im(X_k) = Im(X_{N-k})
    for (int i = 1; i < (std::min(NC,NR)+1)/2; ++i) { // ...up to grid modes
        real_type a_real, a_imag, b_real, b_imag;
        fftw_multi_array::detail::assign_components(a_real, a_imag, out[i]);
        fftw_multi_array::detail::assign_components(b_real, b_imag, out[NC-i]);
        BOOST_REQUIRE_CLOSE(a_real, -b_real, close_enough);
        BOOST_REQUIRE_CLOSE(a_imag,  b_imag, close_enough);

        // We should also see the expected frequency content
        const real_type real_expected = - i * cos(shift) * 1.0/2.0;
        const real_type imag_expected =   i * sin(shift) * 1.0/2.0;
        BOOST_REQUIRE_CLOSE(b_real, real_expected, sqrt(close_enough));
        BOOST_REQUIRE_CLOSE(b_imag, imag_expected, sqrt(close_enough));
    }
    // Ensure we see the expected zero mode magnitude
    {
        real_type z_real, z_imag;
        fftw_multi_array::detail::assign_components(z_real, z_imag, out[0]);
        BOOST_REQUIRE_SMALL(z_real, close_enough);
        BOOST_REQUIRE_CLOSE(z_imag, 17, close_enough);
    }
}

// Helper function that kicks the tires of a 1D r2c transform
template<class RealMultiArray, class ComplexMultiArray>
void check_1D_half_forward(RealMultiArray &in, ComplexMultiArray &out)
{
    BOOST_STATIC_ASSERT(RealMultiArray::dimensionality == 1);
    BOOST_STATIC_ASSERT(ComplexMultiArray::dimensionality == 1);
    typedef typename fftw_multi_array::detail::transform_traits<
        typename ComplexMultiArray::element>::real_type real_type;
    const int NR = in.shape()[0];
    const int NC = out.shape()[0];
    const real_type close_enough
        = std::numeric_limits<real_type>::epsilon()*10*NR*NR;
    const real_type shift = M_PI/3.0;

    // Load a real-valued function into the input array and transform it
    fill_with_complex_NaN(out);
    for (int i = 0; i < NR; ++i) {
        in[i] = real_test_function<real_type>(NR, (NR+1)/2, i, shift);
    }
    fftw_multi_array::forward_r2c(0, in, out);

    // We should see the expected frequency content up to grid modes
    for (int i = 1; i < (std::min(NC,NR)+1)/2; ++i) {
        real_type z_real, z_imag;
        fftw_multi_array::detail::assign_components(z_real, z_imag, out[i]);
        const real_type real_expected =  i * sin(shift) * 1.0/2.0;
        const real_type imag_expected = -i * cos(shift) * 1.0/2.0;
        BOOST_REQUIRE_CLOSE(z_real, real_expected, sqrt(close_enough));
        BOOST_REQUIRE_CLOSE(z_imag, imag_expected, sqrt(close_enough));
    }
    // Ensure we see the expected zero mode magnitude
    {
        real_type z_real, z_imag;
        fftw_multi_array::detail::assign_components(z_real, z_imag, out[0]);
        BOOST_REQUIRE_SMALL(z_imag, close_enough);
        BOOST_REQUIRE_CLOSE(z_real, 17.0, close_enough);
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
        = std::numeric_limits<double>::epsilon()*1e2*NC*NC;

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
    fftw_multi_array::forward_c2c(0, in, out);

    // Ensure we got the same result, up to scaling differences
    for (int i = 0; i < NR; ++i) {
        double out_real, out_imag;
        fftw_multi_array::detail::assign_components(
                out_real, out_imag, out[i]);
        // FFTW with a stride gives a different result than with stride 1
        // BOOST_REQUIRE_EQUAL would be nice, but it fails here
        if (fabs(buffer[i][0]) < close_enough) {
            BOOST_REQUIRE_SMALL(out_real*NC, close_enough);
        } else {
            BOOST_REQUIRE_CLOSE(out_real*NC, buffer[i][0], close_enough);
        }
        if (fabs(buffer[i][1]) < close_enough) {
            BOOST_REQUIRE_SMALL(out_imag*NC, close_enough);
        } else {
            BOOST_REQUIRE_CLOSE(out_imag*NC, buffer[i][1], close_enough);
        }
    }
}

// Helper function that kicks the tires of a 1D c2c transform
template<class ComplexMultiArray1, class ComplexMultiArray2>
void symmetry_1D_complex_backward(ComplexMultiArray1 &in, ComplexMultiArray2 &out)
{
    BOOST_STATIC_ASSERT(ComplexMultiArray1::dimensionality == 1);
    BOOST_STATIC_ASSERT(ComplexMultiArray2::dimensionality == 1);
    typedef typename fftw_multi_array::detail::transform_traits<
        typename ComplexMultiArray2::element>::real_type real_type;
    const int NC = in.shape()[0];
    const int NR = out.shape()[0];
    const real_type close_enough
        = std::numeric_limits<real_type>::epsilon()*10*NC*NC;

    // Load a conjugate-symmetric function into the input array...
    // ...with known frequency content and constant offset 17
    fill_with_complex_NaN(out);
    fill_with_complex_zero(in);
    fftw_multi_array::detail::assign_complex(in[0], 17.0, 0.0);
    for (int i = 1; i < (std::min(NC,NR)+1)/2; ++i) {
        const real_type val = i/2.0;
        fftw_multi_array::detail::assign_complex(in[i],    val, -val);
        fftw_multi_array::detail::assign_complex(in[NC-i], val,  val);
    }
    // TODO Test behavior of highest even half mode

    fftw_multi_array::backward_c2c(0, in, out);

    // Output should be a real-valued function in physical space...
    for (int i = 0; i < NR; ++i) {
        real_type a_real, a_imag;
        fftw_multi_array::detail::assign_components(a_real, a_imag, out[i]);
        BOOST_REQUIRE_SMALL(a_imag, close_enough);

        // ...with known physical space values
        const real_type xi = i*2*M_PI/NR;
        real_type expected_real = 17;
        for (int j = 1; j < (std::min(NC,NR)+1)/2; ++j) {
            expected_real += j*sin(j*xi) + j*cos(j*xi);
        }
        BOOST_REQUIRE_CLOSE(a_real, expected_real, sqrt(close_enough));
    }

    // Load a real-symmetric function into the input array...
    // ...with known frequency content and constant offset 17
    fill_with_complex_NaN(out);
    fill_with_complex_zero(in);
    fftw_multi_array::detail::assign_complex(in[0], 0.0, 17.0);
    for (int i = 1; i < (std::min(NC,NR)+1)/2; ++i) {
        const real_type val = i/2.0;
        fftw_multi_array::detail::assign_complex(in[i],     val, val);
        fftw_multi_array::detail::assign_complex(in[NC-i], -val, val);
    }
    // TODO Test behavior of highest even half mode

    fftw_multi_array::backward_c2c(0, in, out);

    // Output should be an imaginary-valued function in physical space...
    for (int i = 0; i < NR; ++i) {
        real_type a_real, a_imag;
        fftw_multi_array::detail::assign_components(a_real, a_imag, out[i]);
        BOOST_REQUIRE_SMALL(a_real, close_enough);

        // ...with known physical space values
        const real_type xi = i*2*M_PI/NR;
        real_type expected_imag = 17;
        for (int j = 1; j < (std::min(NC,NR)+1)/2; ++j) {
            expected_imag += j*sin(j*xi) + j*cos(j*xi);
        }
        BOOST_REQUIRE_CLOSE(a_imag, expected_imag, sqrt(close_enough));
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
        = std::numeric_limits<double>::epsilon()*1e2*NC*NC;

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
    // Important for FFTW to be first to handle in === out

    // Ramp up amplitudes for FFTW's backwards transform
    for (int i = 0; i < NC; ++i) {
        fftw_multi_array::detail::assign_complex_scaled(in[i], in[i], NC);
    }
    fftw_execute(plan.get());
    // Ramp down amplitudes for our backwards transform
    for (int i = 0; i < NC; ++i) {
        fftw_multi_array::detail::assign_complex_scaled(in[i], in[i], 1.0/NC);
    }
    fftw_multi_array::backward_c2c(0, in, out);

    // Ensure we got exactly the same result after normalization
    for (int i = 0; i < NR; ++i) {
        double out_real, out_imag;
        fftw_multi_array::detail::assign_components(
                out_real, out_imag, out[i]);
        // FFTW with a stride gives a different result than with stride 1
        // BOOST_REQUIRE_EQUAL would be nice, but it fails here
        BOOST_REQUIRE_CLOSE(out_real, buffer[i][0] / NR, close_enough);
        BOOST_REQUIRE_CLOSE(out_imag, buffer[i][1] / NR, close_enough);
    }
}

template<class ComplexMultiArray1, class ComplexMultiArray2>
void differentiate_on_forward_1D_complex(ComplexMultiArray1 &in,
                                         ComplexMultiArray2 &out)
{
    BOOST_STATIC_ASSERT(ComplexMultiArray1::dimensionality == 1);
    BOOST_STATIC_ASSERT(ComplexMultiArray2::dimensionality == 1);
    const int NR = in.shape()[0];
    const int NC = out.shape()[0];
    const double close_enough
        = std::numeric_limits<double>::epsilon()*10*NR*NR;
    const double shift = M_PI/3.0;

    const double length[2] = { 2.0 * M_PI, 10.0 };
    for (int l = 0; l < sizeof(length)/sizeof(length[0]); ++l) {
        for (int derivative = 0; derivative < 8; ++derivative) {
            // Load a complex function into the input array
            fill_with_complex_NaN(in);
            fill_with_complex_NaN(out);
            for (int i = 0; i < NR; ++i) {
                const double val = real_test_function<double>(
                        NR, (std::min(NR,NC)+1)/2, i, shift, length[l], 0);
                fftw_multi_array::detail::assign_complex(in[i], val, -val);
            }

            // Forward transform and differentiate
            fftw_multi_array::forward_c2c(0, in, out, length[l], derivative);

            // Backwards transform without differentiating
            fftw_multi_array::backward_c2c(0, out, in, length[l], 0);

            // Ensure we see what we expect
            for (int i = 0; i < NR; ++i) {
                const double expected_val = real_test_function<double>(
                        NR, (std::min(NR,NC)+1)/2, i, shift, length[l],
                        derivative);
                double real, imag;
                fftw_multi_array::detail::assign_components(real, imag, in[i]);
                if (fabs(expected_val) < close_enough) {
                    BOOST_REQUIRE_SMALL(real, close_enough);
                    BOOST_REQUIRE_SMALL(imag, close_enough);
                } else {
                    BOOST_REQUIRE_CLOSE(
                            expected_val, real, sqrt(close_enough));
                    BOOST_REQUIRE_CLOSE(
                           -expected_val, imag, sqrt(close_enough));
                }
            }
        }
    }
}

template<class ComplexMultiArray1, class ComplexMultiArray2>
void differentiate_on_backward_1D_complex(ComplexMultiArray1 &in,
                                          ComplexMultiArray2 &out)
{
    BOOST_STATIC_ASSERT(ComplexMultiArray1::dimensionality == 1);
    BOOST_STATIC_ASSERT(ComplexMultiArray2::dimensionality == 1);
    const int NR = in.shape()[0];
    const int NC = out.shape()[0];
    const double close_enough
        = std::numeric_limits<double>::epsilon()*10*NR*NR;
    const double shift = M_PI/3.0;

    const double length[2] = { 2.0 * M_PI, 10.0 };
    for (int l = 0; l < sizeof(length)/sizeof(length[0]); ++l) {
        for (int derivative = 0; derivative < 8; ++derivative) {
            // Load a complex function into the input array
            fill_with_complex_NaN(in);
            fill_with_complex_NaN(out);
            for (int i = 0; i < NR; ++i) {
                const double val = real_test_function<double>(
                        NR, (std::min(NR,NC)+1)/2, i, shift, length[l], 0);
                fftw_multi_array::detail::assign_complex(in[i], val, -val);
            }

            // Forward transform without differentiating
            fftw_multi_array::forward_c2c(0, in, out, length[l], 0);

            // Backwards transform and differentiate
            fftw_multi_array::backward_c2c(0, out, in, length[l], derivative);

            // Ensure we see what we expect as the derivative
            for (int i = 0; i < NR; ++i) {
                const double expected_val = real_test_function<double>(
                        NR, (std::min(NR,NC)+1)/2, i, shift, length[l],
                        derivative);
                double real, imag;
                fftw_multi_array::detail::assign_components(real, imag, in[i]);
                if (fabs(expected_val) < close_enough) {
                    BOOST_REQUIRE_SMALL(real, close_enough);
                    BOOST_REQUIRE_SMALL(imag, close_enough);
                } else {
                    BOOST_REQUIRE_CLOSE(
                             expected_val, real, sqrt(close_enough));
                    BOOST_REQUIRE_CLOSE(
                            -expected_val, imag, sqrt(close_enough));
                }
            }
        }
    }
}


/* powers of 2; simple composites of form 2*prime; prime numbers; 1 */
#define TRANSFORM_1D_SIZE_SEQ \
        (2)(4)(8)(16)(32)(64) \
        (6)(10)(14)(22)(26)(34) \
        (3)(5)(7)(9)(11)(13)(17)(19)(23)(29) \
        (1)

BOOST_AUTO_TEST_SUITE( test_1d_out_of_place );
#define TEST_1D_OUT_OF_PLACE(r, data, elem) \
        BOOST_AUTO_TEST_CASE( BOOST_PP_CAT(test_1d_out_of_place_,elem) ) \
        { test_1d_out_of_place(elem); }
void test_1d_out_of_place(const int N)
{
    // C2C: Test multi_array using std::complex
    {
        typedef boost::multi_array<std::complex<double>,1> array_type;
        array_type in(boost::extents[N]), out(boost::extents[N]);

        // No dealiasing in effect: size in == size out
        symmetry_1D_complex_forward(in, out);
        compare_1D_complex_forward(in, out);
        symmetry_1D_complex_backward(in, out);
        compare_1D_complex_backward(in, out);
        differentiate_on_forward_1D_complex(in, out);
        differentiate_on_backward_1D_complex(in, out);
    }

    // C2C: Test multi_array_ref using fftw_complex
    {
        boost::scoped_array<fftw_complex> in_data(new fftw_complex[N]);
        boost::scoped_array<fftw_complex> out_data(new fftw_complex[N]);
        typedef boost::multi_array_ref<fftw_complex, 1> array_ref_type;
        array_ref_type in(in_data.get(), boost::extents[N]);
        array_ref_type out(out_data.get(), boost::extents[N]);

        // No dealiasing in effect: size in == size out
        symmetry_1D_complex_forward(in, out);
        compare_1D_complex_forward(in, out);
        symmetry_1D_complex_backward(in, out);
        compare_1D_complex_backward(in, out);
        differentiate_on_forward_1D_complex(in, out);
        differentiate_on_backward_1D_complex(in, out);
    }

    // R2C: Test multi_array using std::complex
    {
        boost::multi_array<double,1>               in(boost::extents[N]);
        boost::multi_array<std::complex<double>,1> out(boost::extents[N/2+1]);

        // No dealiasing in effect
        check_1D_half_forward(in, out);
    }

    // R2C: Test multi_array using fftw_complex
    {
        boost::multi_array<double,1>      in(boost::extents[N]);
        boost::scoped_array<fftw_complex> out_data(new fftw_complex[N/2+1]);
        boost::multi_array_ref<fftw_complex, 1> out(
                out_data.get(), boost::extents[N/2+1]);

        // No dealiasing in effect
        check_1D_half_forward(in, out);
    }
}
BOOST_PP_SEQ_FOR_EACH(TEST_1D_OUT_OF_PLACE,_,TRANSFORM_1D_SIZE_SEQ);
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

    symmetry_1D_complex_forward(in, out);
    compare_1D_complex_forward(in, out);
    symmetry_1D_complex_backward(in, out);
    compare_1D_complex_backward(in, out);
    differentiate_on_forward_1D_complex(in, out);
    differentiate_on_backward_1D_complex(in, out);
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

    symmetry_1D_complex_forward(in, out);
    compare_1D_complex_forward(in, out);
    symmetry_1D_complex_backward(in, out);
    compare_1D_complex_backward(in, out);
    differentiate_on_forward_1D_complex(in, out);
    differentiate_on_backward_1D_complex(in, out);
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
    symmetry_1D_complex_forward(both, both);
    compare_1D_complex_forward(both, both);
    symmetry_1D_complex_backward(both, both);
    compare_1D_complex_forward(both, both);
    differentiate_on_forward_1D_complex(both, both);
    differentiate_on_backward_1D_complex(both, both);
}
BOOST_PP_SEQ_FOR_EACH(TEST_C2C_1D_IN_PLACE,_,TRANSFORM_1D_SIZE_SEQ);
BOOST_AUTO_TEST_SUITE_END();

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
    symmetry_1D_complex_forward(in, out); // Dealiasing in effect
    differentiate_on_forward_1D_complex(in, out);
}
void c2c_1d_out_of_place_dealiased_backward(const int NC, const int NR)
{
    typedef boost::multi_array<std::complex<double>,1> array_type;
    array_type in(boost::extents[NC]), out(boost::extents[NR]);
    symmetry_1D_complex_backward(in, out); // Dealiasing in effect
    differentiate_on_backward_1D_complex(out, in); // NB reversed
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
