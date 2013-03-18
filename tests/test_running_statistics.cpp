//--------------------------------------------------------------------------
//
// Copyright (C) 2012, 2013 Rhys Ulerich
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

#include <suzerain/running_statistics.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>

// Explicit instantiation to flush out compilation errors
template class suzerain::running_statistics<double, 5>;
// template class suzerain::running_statistics<double, 1>;
// template class suzerain::running_statistics<double, 0>;
// template class suzerain::running_statistics<float,  3>;

// TODO Implement
BOOST_AUTO_TEST_SUITE(one_suite)

BOOST_AUTO_TEST_CASE( totally_suite )
{
    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()
