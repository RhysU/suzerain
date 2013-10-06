/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2013 Rhys Ulerich
 * Copyright (C) 2013 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Suzerain is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 */

#include <suzerain/bl.h>

#define BOOST_TEST_MAIN
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>
#include <suzerain/bspline.hpp>

#include "test_tools.hpp"

using suzerain::bspline;
using suzerain::bsplineop;
using suzerain::bsplineop_lu;

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

#pragma warning(disable:1572 2014 2015)

BOOST_AUTO_TEST_SUITE(bl_compute_viscous)
// FIXME Implement
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(bl_compute_thick)

// FIXME Remove
BOOST_AUTO_TEST_CASE( placeholder ) { BOOST_CHECK(true); }

// FIXME Test suzerain_bl_find_edge

// FIXME Test suzerain_bl_compute_deltastar

// FIXME Test suzerain_bl_compute_theta

// FIXME Test suzerain_bl_compute_thick

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(bl_compute_qoi)
// FIXME Implement
BOOST_AUTO_TEST_SUITE_END()
