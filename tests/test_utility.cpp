#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/assign.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <suzerain/utility.hpp>

BOOST_AUTO_TEST_CASE( to_yxz )
{
    using suzerain::to_yxz;
    {
        boost::array<char,3> xyz = {{ 'x', 'y', 'z' }};
        BOOST_CHECK_EQUAL(to_yxz(xyz),
                          boost::assign::list_of('y')('x')('z'));
    }
    {
        boost::array<char,3> xyz = {{ 'x', 'y', 'z' }};
        BOOST_CHECK_EQUAL(to_yxz('w', xyz),
                          boost::assign::list_of('w')('y')('x')('z'));
    }
}

BOOST_AUTO_TEST_CASE( to_xzy )
{
    using suzerain::to_xzy;
    {
        boost::array<char,3> xyz = {{ 'x', 'y', 'z' }};
        BOOST_CHECK_EQUAL(to_xzy(xyz),
                          boost::assign::list_of('x')('z')('y'));
    }
    {
        boost::array<char,3> xyz = {{ 'x', 'y', 'z' }};
        BOOST_CHECK_EQUAL(to_xzy('w', xyz),
                          boost::assign::list_of('w')('x')('z')('y'));
    }
}

BOOST_AUTO_TEST_CASE( prepend )
{
    using suzerain::prepend;
    {
        boost::array<int,0> a;
        BOOST_CHECK_EQUAL(prepend(0, a),
                          boost::assign::list_of(0));
    }

    {
        boost::array<int,1> a = {{ 2 }};
        BOOST_CHECK_EQUAL(prepend(0, a),
                          boost::assign::list_of(0)(2));
    }

    {
        boost::array<int,2> a = {{ 2, 3 }};
        BOOST_CHECK_EQUAL(prepend(0, a),
                          boost::assign::list_of(0)(2)(3));
    }

    {
        boost::array<int,3> a = {{ 2, 3, 4 }};
        BOOST_CHECK_EQUAL(prepend(0, a),
                          boost::assign::list_of(0)(2)(3)(4));
    }
}

BOOST_AUTO_TEST_CASE( strides_cm )
{
    using suzerain::strides_cm;

    {
        boost::array<int,1> extents = {{ 2 }};
        BOOST_CHECK_EQUAL(strides_cm(extents),
                          boost::assign::list_of(1));
    }

    {
        boost::array<int,2> extents = {{ 2, 3 }};
        BOOST_CHECK_EQUAL(strides_cm(extents),
                          boost::assign::list_of(1)(2));
    }

    {
        boost::array<int,3> extents = {{ 2, 3, 4 }};
        BOOST_CHECK_EQUAL(strides_cm(extents),
                          boost::assign::list_of(1)(2)(6));
    }

    {
        boost::array<int,3> extents = {{ 13, 16, 24 }};
        BOOST_CHECK_EQUAL(strides_cm(extents),
                          boost::assign::list_of(1)(13)(208));
    }

    {
        boost::array<int,4> extents = {{ 2, 3, 4, 5 }};
        BOOST_CHECK_EQUAL(strides_cm(extents),
                          boost::assign::list_of(1)(2)(6)(24));
    }
}

BOOST_AUTO_TEST_CASE( strides_rm )
{
    using suzerain::strides_rm;

    {
        boost::array<int,1> extents = {{ 2 }};
        BOOST_CHECK_EQUAL(strides_rm(extents),
                          boost::assign::list_of(1));
    }

    {
        boost::array<int,2> extents = {{ 2, 3 }};
        BOOST_CHECK_EQUAL(strides_rm(extents),
                          boost::assign::list_of(3)(1));
    }

    {
        boost::array<int,3> extents = {{ 2, 3, 4 }};
        BOOST_CHECK_EQUAL(strides_rm(extents),
                          boost::assign::list_of(12)(4)(1));
    }

    {
        boost::array<int,4> extents = {{ 2, 3, 4, 5 }};
        BOOST_CHECK_EQUAL(strides_rm(extents),
                          boost::assign::list_of(60)(20)(5)(1));
    }
}

BOOST_AUTO_TEST_CASE( is_nonnegative )
{
   BOOST_CHECK(!suzerain::is_nonnegative<int>(-1) );
   BOOST_CHECK( suzerain::is_nonnegative<int>( 0) );
   BOOST_CHECK( suzerain::is_nonnegative<int>( 1) );

   BOOST_CHECK(!suzerain::is_nonnegative<long>(-1) );
   BOOST_CHECK( suzerain::is_nonnegative<long>( 0) );
   BOOST_CHECK( suzerain::is_nonnegative<long>( 1) );

   BOOST_CHECK( suzerain::is_nonnegative<unsigned int>(0U) );
   BOOST_CHECK( suzerain::is_nonnegative<unsigned int>(1U) );

   BOOST_CHECK( suzerain::is_nonnegative<unsigned long>(0U) );
   BOOST_CHECK( suzerain::is_nonnegative<unsigned long>(1U) );

   BOOST_CHECK(!suzerain::is_nonnegative<float>(-1.0) );
   BOOST_CHECK( suzerain::is_nonnegative<float>( 0.0) );
   BOOST_CHECK( suzerain::is_nonnegative<float>( 1.0) );

   BOOST_CHECK(!suzerain::is_nonnegative<double>(-1.0) );
   BOOST_CHECK( suzerain::is_nonnegative<double>( 0.0) );
   BOOST_CHECK( suzerain::is_nonnegative<double>( 1.0) );
}
