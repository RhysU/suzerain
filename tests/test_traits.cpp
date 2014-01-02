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
