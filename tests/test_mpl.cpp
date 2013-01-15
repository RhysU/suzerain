//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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

#include <suzerain/mpl.hpp>

#define BOOST_TEST_MAIN
#include <boost/mpl/vector_c.hpp>
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>

BOOST_AUTO_TEST_SUITE( sequence_array )

using boost::mpl::vector_c;
using suzerain::mpl::sequence_array;

BOOST_AUTO_TEST_CASE( integer )
{
    sequence_array< boost::mpl::vector_c<int,2,4,6> > a;
    BOOST_CHECK_EQUAL(a.size(), 3U);

    BOOST_CHECK_EQUAL(a[0], 2);
    BOOST_CHECK_EQUAL(a[1], 4);
    BOOST_CHECK_EQUAL(a[2], 6);

    suzerain::array<int, 3> b = {{ 2, 4, 6 }};

    BOOST_CHECK( b == a );
    BOOST_CHECK( a == b );

}

BOOST_AUTO_TEST_CASE( boolean )
{
    sequence_array< boost::mpl::vector_c<bool,true,false> > a;
    BOOST_CHECK_EQUAL(a.size(), 2U);

    BOOST_CHECK_EQUAL(a[0], true);
    BOOST_CHECK_EQUAL(a[1], false);

    suzerain::array<bool, 2> b = {{ true, false }};

    BOOST_CHECK( b == a );
    BOOST_CHECK( a == b );
}

BOOST_AUTO_TEST_CASE( empty )
{
    sequence_array< boost::mpl::vector_c<double> > a;
    BOOST_CHECK_EQUAL(a.size(), 0U);

    suzerain::array<double, 0> b;

    BOOST_CHECK( b == a );
    BOOST_CHECK( a == b );
}

BOOST_AUTO_TEST_SUITE_END()
