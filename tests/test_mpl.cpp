#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <suzerain/mpl.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE( sequence_array )

using boost::mpl::vector_c;
using suzerain::mpl::sequence_array;

BOOST_AUTO_TEST_CASE( integer )
{
    sequence_array< boost::mpl::vector_c<int,2,4,6> > a;
    BOOST_CHECK_EQUAL(a.size(), 3);

    BOOST_CHECK_EQUAL(a[0], 2);
    BOOST_CHECK_EQUAL(a[1], 4);
    BOOST_CHECK_EQUAL(a[2], 6);

    boost::array<int, 3> b = {{ 2, 4, 6 }};

    BOOST_CHECK( b == a );
    BOOST_CHECK( a == b );

}

BOOST_AUTO_TEST_CASE( boolean )
{
    sequence_array< boost::mpl::vector_c<bool,true,false> > a;
    BOOST_CHECK_EQUAL(a.size(), 2);

    BOOST_CHECK_EQUAL(a[0], true);
    BOOST_CHECK_EQUAL(a[1], false);

    boost::array<bool, 2> b = {{ true, false }};

    BOOST_CHECK( b == a );
    BOOST_CHECK( a == b );
}

BOOST_AUTO_TEST_CASE( empty )
{
    sequence_array< boost::mpl::vector_c<double> > a;
    BOOST_CHECK_EQUAL(a.size(), 0);

    boost::array<double, 0> b;

    BOOST_CHECK( b == a );
    BOOST_CHECK( a == b );
}

BOOST_AUTO_TEST_SUITE_END()
