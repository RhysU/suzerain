//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <suzerain/utility.hpp>
#define BOOST_TEST_MAIN
#include <boost/assign.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( to_yxz )
{
    using suzerain::to_yxz;
    {
        boost::array<char,3> xyz = {{ 'x', 'y', 'z' }};
        BOOST_CHECK_EQUAL(to_yxz(xyz),
                          boost::assign::list_of('y')('x')('z'));
    }
    {
        Eigen::Vector3i xyz(1,2,3);
        BOOST_CHECK_EQUAL(to_yxz(xyz),
                          boost::assign::list_of(2)(1)(3));
    }
    {
        boost::array<char,3> xyz = {{ 'x', 'y', 'z' }};
        BOOST_CHECK_EQUAL(to_yxz('w', xyz),
                          boost::assign::list_of('w')('y')('x')('z'));
    }
    {
        Eigen::Vector3i xyz(1,2,3);
        BOOST_CHECK_EQUAL(to_yxz(4,xyz),
                          boost::assign::list_of(4)(2)(1)(3));
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
        Eigen::Vector3i xyz(1,2,3);
        BOOST_CHECK_EQUAL(to_xzy(xyz),
                          boost::assign::list_of(1)(3)(2));
    }
    {
        boost::array<char,3> xyz = {{ 'x', 'y', 'z' }};
        BOOST_CHECK_EQUAL(to_xzy('w', xyz),
                          boost::assign::list_of('w')('x')('z')('y'));
    }
    {
        Eigen::Vector3i xyz(1,2,3);
        BOOST_CHECK_EQUAL(to_xzy(4, xyz),
                          boost::assign::list_of(4)(1)(3)(2));
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

BOOST_AUTO_TEST_CASE( any )
{
   boost::array<bool,0> empty;
   boost::array<bool,1> f     = {{ false }};
   boost::array<bool,1> t     = {{ true  }};
   boost::array<bool,2> ff    = {{ false, false }};
   boost::array<bool,2> ft    = {{ false, true  }};
   boost::array<bool,2> tf    = {{ true,  false }};
   boost::array<bool,2> tt    = {{ true,  true  }};

   BOOST_CHECK(!suzerain::any(empty.begin(), empty.end()));
   BOOST_CHECK(!suzerain::any(f.begin(),     f.end()    ));
   BOOST_CHECK( suzerain::any(t.begin(),     t.end()    ));
   BOOST_CHECK(!suzerain::any(ff.begin(),    ff.end()   ));
   BOOST_CHECK( suzerain::any(ft.begin(),    ft.end()   ));
   BOOST_CHECK( suzerain::any(tf.begin(),    tf.end()   ));
   BOOST_CHECK( suzerain::any(tt.begin(),    tt.end()   ));
}

BOOST_AUTO_TEST_CASE( all )
{
   boost::array<bool,0> empty;
   boost::array<bool,1> f     = {{ false }};
   boost::array<bool,1> t     = {{ true  }};
   boost::array<bool,2> ff    = {{ false, false }};
   boost::array<bool,2> ft    = {{ false, true  }};
   boost::array<bool,2> tf    = {{ true,  false }};
   boost::array<bool,2> tt    = {{ true,  true  }};

   BOOST_CHECK( suzerain::all(empty.begin(), empty.end()));
   BOOST_CHECK(!suzerain::all(f.begin(),     f.end()    ));
   BOOST_CHECK( suzerain::all(t.begin(),     t.end()    ));
   BOOST_CHECK(!suzerain::all(ff.begin(),    ff.end()   ));
   BOOST_CHECK(!suzerain::all(ft.begin(),    ft.end()   ));
   BOOST_CHECK(!suzerain::all(tf.begin(),    tf.end()   ));
   BOOST_CHECK( suzerain::all(tt.begin(),    tt.end()   ));
}
