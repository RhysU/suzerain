#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <suzerain/traits.hpp>
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#pragma warning(disable:1572)

BOOST_AUTO_TEST_CASE( component_type_real )
{
    using boost::is_same;
    using suzerain::traits::component;

    BOOST_CHECK((is_same<component<float >::type, float>::value));
    BOOST_CHECK((is_same<component<double>::type, double>::value));
    BOOST_CHECK((is_same<component<int   >::type, int>::value));
}

BOOST_AUTO_TEST_CASE( component_type_complex )
{
    using boost::is_same;
    using std::complex;
    using suzerain::traits::component;

    BOOST_CHECK((is_same<component<complex<float>  >::type, float>::value));
    BOOST_CHECK((is_same<component<complex<double> >::type, double>::value));
    BOOST_CHECK((is_same<component<complex<int>    >::type, int>::value));
}
