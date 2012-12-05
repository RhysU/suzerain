//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
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

#include <suzerain/traits.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>

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
