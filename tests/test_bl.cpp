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
#include <suzerain/countof.h>
#include <suzerain/blasius.h>
#include <suzerain/bspline.hpp>

#include "test_tools.hpp"

using suzerain::bspline;
using suzerain::bsplineop;
using suzerain::bsplineop_lu;
using suzerain::shared_array;
using suzerain::shared_ptr;

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

#pragma warning(disable:1572 2014 2015)

BOOST_AUTO_TEST_SUITE(bl_compute_viscous)
// FIXME Implement
BOOST_AUTO_TEST_SUITE_END()

// A test fixture making test profile(s) and B-splines available
struct ProfileFixture {

    static const double breakpts[10]; // Init just below
    gsl_spline * const blasius_u_vs_eta;
    gsl_interp_accel * accel;
    bspline      b;
    bsplineop    op;
    bsplineop_lu lu;

    ProfileFixture()
        : blasius_u_vs_eta(suzerain_blasius_u_vs_eta())
        , accel(gsl_interp_accel_alloc())
        , b(8, bspline::from_breakpoints(),
            SUZERAIN_COUNTOF(breakpts), breakpts)
        , op(b, 2, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE)
        , lu(op)
    {
        lu.factor_mass(op);
    }

    ~ProfileFixture()
    {
        gsl_spline_free(blasius_u_vs_eta);
        gsl_interp_accel_free(accel);
    }

};

// As the data for Ganapol's Blasius profile runs up to 8.8 instead of 5.0, and
// the routines compute edge quantities from the profiles, our basis runs up to
// 8.8 as well.
const double ProfileFixture::breakpts[10] = {
    0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.8
};

BOOST_FIXTURE_TEST_SUITE(bl_compute_thick, ProfileFixture)

// FIXME Test suzerain_bl_find_edge

BOOST_AUTO_TEST_CASE( blasius_deltastar )
{
    // Prepare the Blasius velocity profile as coefficients on basis
    shared_array<double> rho_u(new double[b.n()]);
    for (int i = 0; i < b.n(); ++i) {
        rho_u[i] = gsl_spline_eval(blasius_u_vs_eta,
                                   b.collocation_point(i), accel);
    }
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.solve(1, rho_u.get(), 1, b.n()));

    // Prepare integration working storage
    shared_ptr<gsl_matrix> dB(
            gsl_matrix_alloc(b.k(), 1),
            gsl_matrix_free);
    BOOST_REQUIRE(dB);
    shared_ptr<gsl_integration_workspace> iw(
            gsl_integration_workspace_alloc(256),
            gsl_integration_workspace_free);
    BOOST_REQUIRE(iw);

    // Integrate for deltastar
    double deltastar = GSL_NAN;
    double abserr    = GSL_NAN;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_compute_deltastar(
        b.collocation_point(b.n()-1), rho_u.get(), &deltastar, dB.get(), b.bw,
        b.dbw, iw.get(), GSL_SQRT_DBL_EPSILON, GSL_SQRT_DBL_EPSILON, &abserr));

    // Check against good value
    // Good value taken from White, Fluid Mechanics, 4th Edition eqn (7.31).
    BOOST_CHECK_CLOSE(1.721, deltastar, 0.013);
    BOOST_CHECK_LE(abserr, GSL_SQRT_DBL_EPSILON);
}

// FIXME Test suzerain_bl_compute_theta

// FIXME Test suzerain_bl_compute_thick

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(bl_compute_qoi)
// FIXME Implement
BOOST_AUTO_TEST_SUITE_END()
