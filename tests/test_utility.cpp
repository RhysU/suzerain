//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Suzerain is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------

#include <suzerain/utility.hpp>

#define BOOST_TEST_MAIN
#include <boost/assign.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>

BOOST_AUTO_TEST_CASE( to_yxz )
{
    using suzerain::to_yxz;
    {
        suzerain::array<char,3> xyz = {{ 'x', 'y', 'z' }};
        BOOST_CHECK_EQUAL(to_yxz(xyz),
                          boost::assign::list_of('y')('x')('z'));
    }
    {
        suzerain::Vector3i xyz(1,2,3);
        BOOST_CHECK_EQUAL(to_yxz(xyz),
                          boost::assign::list_of(2)(1)(3));
    }
    {
        suzerain::array<char,3> xyz = {{ 'x', 'y', 'z' }};
        BOOST_CHECK_EQUAL(to_yxz('w', xyz),
                          boost::assign::list_of('w')('y')('x')('z'));
    }
    {
        suzerain::Vector3i xyz(1,2,3);
        BOOST_CHECK_EQUAL(to_yxz(4,xyz),
                          boost::assign::list_of(4)(2)(1)(3));
    }
}

BOOST_AUTO_TEST_CASE( to_xzy )
{
    using suzerain::to_xzy;
    {
        suzerain::array<char,3> xyz = {{ 'x', 'y', 'z' }};
        BOOST_CHECK_EQUAL(to_xzy(xyz),
                          boost::assign::list_of('x')('z')('y'));
    }
    {
        suzerain::Vector3i xyz(1,2,3);
        BOOST_CHECK_EQUAL(to_xzy(xyz),
                          boost::assign::list_of(1)(3)(2));
    }
    {
        suzerain::array<char,3> xyz = {{ 'x', 'y', 'z' }};
        BOOST_CHECK_EQUAL(to_xzy('w', xyz),
                          boost::assign::list_of('w')('x')('z')('y'));
    }
    {
        suzerain::Vector3i xyz(1,2,3);
        BOOST_CHECK_EQUAL(to_xzy(4, xyz),
                          boost::assign::list_of(4)(1)(3)(2));
    }
}

BOOST_AUTO_TEST_CASE( prepend )
{
    using suzerain::prepend;
    {
        suzerain::array<int,0> a;
        BOOST_CHECK_EQUAL(prepend(0, a),
                          boost::assign::list_of(0));
    }

    {
        suzerain::array<int,1> a = {{ 2 }};
        BOOST_CHECK_EQUAL(prepend(0, a),
                          boost::assign::list_of(0)(2));
    }

    {
        suzerain::array<int,2> a = {{ 2, 3 }};
        BOOST_CHECK_EQUAL(prepend(0, a),
                          boost::assign::list_of(0)(2)(3));
    }

    {
        suzerain::array<int,3> a = {{ 2, 3, 4 }};
        BOOST_CHECK_EQUAL(prepend(0, a),
                          boost::assign::list_of(0)(2)(3)(4));
    }
}

BOOST_AUTO_TEST_CASE( strides_cm )
{
    using suzerain::strides_cm;

    {
        suzerain::array<int,1> extents = {{ 2 }};
        BOOST_CHECK_EQUAL(strides_cm(extents),
                          boost::assign::list_of(1));
    }

    {
        suzerain::array<int,2> extents = {{ 2, 3 }};
        BOOST_CHECK_EQUAL(strides_cm(extents),
                          boost::assign::list_of(1)(2));
    }

    {
        suzerain::array<int,3> extents = {{ 2, 3, 4 }};
        BOOST_CHECK_EQUAL(strides_cm(extents),
                          boost::assign::list_of(1)(2)(6));
    }

    {
        suzerain::array<int,3> extents = {{ 13, 16, 24 }};
        BOOST_CHECK_EQUAL(strides_cm(extents),
                          boost::assign::list_of(1)(13)(208));
    }

    {
        suzerain::array<int,4> extents = {{ 2, 3, 4, 5 }};
        BOOST_CHECK_EQUAL(strides_cm(extents),
                          boost::assign::list_of(1)(2)(6)(24));
    }
}

BOOST_AUTO_TEST_CASE( strides_rm )
{
    using suzerain::strides_rm;

    {
        suzerain::array<int,1> extents = {{ 2 }};
        BOOST_CHECK_EQUAL(strides_rm(extents),
                          boost::assign::list_of(1));
    }

    {
        suzerain::array<int,2> extents = {{ 2, 3 }};
        BOOST_CHECK_EQUAL(strides_rm(extents),
                          boost::assign::list_of(3)(1));
    }

    {
        suzerain::array<int,3> extents = {{ 2, 3, 4 }};
        BOOST_CHECK_EQUAL(strides_rm(extents),
                          boost::assign::list_of(12)(4)(1));
    }

    {
        suzerain::array<int,4> extents = {{ 2, 3, 4, 5 }};
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
   suzerain::array<bool,0> empty;
   suzerain::array<bool,1> f     = {{ false }};
   suzerain::array<bool,1> t     = {{ true  }};
   suzerain::array<bool,2> ff    = {{ false, false }};
   suzerain::array<bool,2> ft    = {{ false, true  }};
   suzerain::array<bool,2> tf    = {{ true,  false }};
   suzerain::array<bool,2> tt    = {{ true,  true  }};

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
   suzerain::array<bool,0> empty;
   suzerain::array<bool,1> f     = {{ false }};
   suzerain::array<bool,1> t     = {{ true  }};
   suzerain::array<bool,2> ff    = {{ false, false }};
   suzerain::array<bool,2> ft    = {{ false, true  }};
   suzerain::array<bool,2> tf    = {{ true,  false }};
   suzerain::array<bool,2> tt    = {{ true,  true  }};

   BOOST_CHECK( suzerain::all(empty.begin(), empty.end()));
   BOOST_CHECK(!suzerain::all(f.begin(),     f.end()    ));
   BOOST_CHECK( suzerain::all(t.begin(),     t.end()    ));
   BOOST_CHECK(!suzerain::all(ff.begin(),    ff.end()   ));
   BOOST_CHECK(!suzerain::all(ft.begin(),    ft.end()   ));
   BOOST_CHECK(!suzerain::all(tf.begin(),    tf.end()   ));
   BOOST_CHECK( suzerain::all(tt.begin(),    tt.end()   ));
}
