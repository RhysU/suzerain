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

// A fixture exposing Ganapol's Blasius profile using eta breakpoints
template <int k>
struct UniformGanapolFixture {

    bspline      b;
    bsplineop    op;
    bsplineop_lu lu;

    UniformGanapolFixture()
        : b(k, bspline::from_breakpoints(),
            SUZERAIN_COUNTOF(suzerain_blasius_ganapol_eta),
            suzerain_blasius_ganapol_eta)
        , op(b, 0, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE)
        , lu(op)
    {
        lu.factor_mass(op);
    }

};

// A fixture exposing Ganapol's Blasius profile using non-uniform breakpoints
template <int k>
struct NonuniformGanapolFixture {

    static const double breakpts[10]; // Init just below
    gsl_spline * const blasius_u;
    gsl_interp_accel * accel;
    bspline      b;
    bsplineop    op;
    bsplineop_lu lu;

    NonuniformGanapolFixture()
        : blasius_u(suzerain_blasius_u())
        , accel(gsl_interp_accel_alloc())
        , b(k, bspline::from_breakpoints(),
            SUZERAIN_COUNTOF(breakpts), breakpts)
        , op(b, 0, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE)
        , lu(op)
    {
        lu.factor_mass(op);
    }

    ~NonuniformGanapolFixture()
    {
        gsl_spline_free(blasius_u);
        gsl_interp_accel_free(accel);
    }

};

// As the data for Ganapol's Blasius profile runs up to 8.8 instead of 5.0, and
// the routines compute edge quantities from the profiles, our basis runs up to
// 8.8 as well.
template <int k>
const double NonuniformGanapolFixture<k>::breakpts[10] = {
    0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.8
};

BOOST_FIXTURE_TEST_SUITE(bl_compute_thick_linear, UniformGanapolFixture<2>)

BOOST_AUTO_TEST_CASE( blasius_deltastar )
{
    // Prepare the Blasius velocity profile as coefficients on basis
    // For a linear B-spline basis, the collocation points are the breakpoints
    shared_array<double> u(new double[b.n()]);
    std::memcpy(u.get(),
                suzerain_blasius_ganapol_fp,
                sizeof(suzerain_blasius_ganapol_fp));
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.solve(1, u.get(), 1, b.n()));

    // Prepare integration working storage
    shared_ptr<gsl_matrix> dB(
            gsl_matrix_alloc(b.k(), 1),
            gsl_matrix_free);
    BOOST_REQUIRE(dB);

    // Integrate for deltastar
    double deltastar = GSL_NAN;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_compute_deltastar(
        b.collocation_point(b.n()-1), u.get(),
        &deltastar, dB.get(), b.bw, b.dbw));

    // Check against good value found using Octave's trapz on Ganapol data
    BOOST_CHECK_SMALL((deltastar - 1.72189445179000), 5e-6);
}

BOOST_AUTO_TEST_CASE( blasius_theta )
{
    // Prepare the Blasius velocity profile as coefficients on basis
    // For a linear B-spline basis, the collocation points are the breakpoints
    shared_array<double> u(new double[b.n()]);
    std::memcpy(u.get(),
                suzerain_blasius_ganapol_fp,
                sizeof(suzerain_blasius_ganapol_fp));
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.solve(1, u.get(), 1, b.n()));

    // Prepare integration working storage
    shared_ptr<gsl_matrix> dB(
            gsl_matrix_alloc(b.k(), 1),
            gsl_matrix_free);
    BOOST_REQUIRE(dB);

    // Integrate for theta
    // Pretend density is nondimensionally one so rho_u == u
    // The absolute error behavior on this integral is unsatisfying
    // though there's no reason adaptive results should match Octave's trapz.
    double theta  = GSL_NAN;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_compute_theta(
        b.collocation_point(b.n()-1), u.get(), u.get(), &theta, dB.get(),
        b.bw, b.dbw));

    // Check against good value found using Octave's trapz on Ganapol data
    BOOST_CHECK_SMALL((theta - 0.663007750711612), 0.0020);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_FIXTURE_TEST_SUITE(bl_compute_thick_quadratic, UniformGanapolFixture<4>)

BOOST_AUTO_TEST_CASE( blasius_find_edge )
{
    const double Re_x = 1e3;

    // Prepare B-spline coefficients for Blasius profile kinetic energy
    shared_array<double> ke(new double[b.n()]);
    {
        shared_ptr<gsl_spline> fit(suzerain_blasius_ke(Re_x),
                                   gsl_spline_free);
        shared_ptr<gsl_interp_accel> a(gsl_interp_accel_alloc(),
                                       gsl_interp_accel_free);
        for (int i = 0; i < b.n(); ++i) {
            ke[i] = gsl_spline_eval(fit.get(), b.collocation_point(i), a.get());
        }
    }
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.solve(1, ke.get(), 1, b.n()));

    // Prepare working storage
    shared_ptr<gsl_matrix> dB(gsl_matrix_alloc(b.k(), 3), gsl_matrix_free);
    BOOST_REQUIRE(dB);

    // Find edge using kinetic energy profile
    double location = GSL_NAN;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_find_edge(
        ke.get(), &location, dB.get(), b.bw, b.dbw));

    // Tolerance from by eyeballing results computed in Octave:
    //   source writeups/notebooks/blasius.m
    //   Re=1000; plot(eta, blasius_kepp(Re), eta, zeros(size(eta)));
    // There's no reason Octave plots should produce exactly this value.
    BOOST_REQUIRE_SMALL((6.485 - location), 0.01);

    // Finally, as a sanity check, be sure the second derivative of KE is small
    // The tests both our zero-crossing logic as well as the ke__yy profile.
    shared_ptr<gsl_spline> fit(suzerain_blasius_ke__yy(Re_x),
                               gsl_spline_free);
    shared_ptr<gsl_interp_accel> a(gsl_interp_accel_alloc(),
                                   gsl_interp_accel_free);
    const double ke__yy = gsl_spline_eval(fit.get(), location, a.get());
    BOOST_REQUIRE_SMALL(ke__yy, 0.01);
}

BOOST_AUTO_TEST_CASE( blasius_compute_thick )
{
    const double Re_x = 1e5;

    shared_ptr<gsl_interp_accel> a(gsl_interp_accel_alloc(),
                                   gsl_interp_accel_free);

    // Prepare B-spline coefficients for kinetic energy
    shared_array<double> ke(new double[b.n()]);
    {
        shared_ptr<gsl_spline> fit(suzerain_blasius_ke(Re_x), gsl_spline_free);
        gsl_interp_accel_reset(a.get());
        for (int i = 0; i < b.n(); ++i)
            ke[i] = gsl_spline_eval(fit.get(), b.collocation_point(i), a.get());
    }
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.solve(1, ke.get(), 1, b.n()));

    // Prepare B-spline coefficients for streamwise velocity
    shared_array<double> u(new double[b.n()]);
    {
        shared_ptr<gsl_spline> fit(suzerain_blasius_u(), gsl_spline_free);
        gsl_interp_accel_reset(a.get());
        for (int i = 0; i < b.n(); ++i)
            u[i] = gsl_spline_eval(fit.get(), b.collocation_point(i), a.get());
    }
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.solve(1, u.get(), 1, b.n()));

    // Assume uniform density of 0.5 and scale to get streamwise momentum
    shared_array<double> rho_u(new double[b.n()]);
    for (int i = 0; i < b.n(); ++i) {
        rho_u[i] = u[i] / 2;
    }

    // Compute a bunch of thickness-related quantities
    // Known good values computed by Octave using trapz from Ganapol data
    size_t cnt = 0;
    suzerain_bl_thick thick;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_compute_thick(
        ke.get() /* \approx H_0 */, rho_u.get(), u.get(), &thick, b.bw, b.dbw));
    BOOST_CHECK_CLOSE(thick.delta,     8.22,              0.25); ++cnt;
    BOOST_CHECK_CLOSE(thick.deltastar, 1.72189445179000,  0.10); ++cnt;
    BOOST_CHECK_CLOSE(thick.theta,     0.663007750711612, 0.25); ++cnt;
    BOOST_CHECK_EQUAL(cnt, sizeof(thick)/sizeof(thick.delta));
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_FIXTURE_TEST_SUITE(bl_compute_thick_splined, NonuniformGanapolFixture<8>)

BOOST_AUTO_TEST_CASE( blasius_deltastar )
{
    // Prepare the Blasius velocity profile as coefficients on basis
    shared_array<double> rho_u(new double[b.n()]);
    for (int i = 0; i < b.n(); ++i) {
        rho_u[i] = gsl_spline_eval(blasius_u,
                                   b.collocation_point(i), accel);
    }
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.solve(1, rho_u.get(), 1, b.n()));

    // Prepare integration working storage
    shared_ptr<gsl_matrix> dB(
            gsl_matrix_alloc(b.k(), 1),
            gsl_matrix_free);
    BOOST_REQUIRE(dB);

    // Integrate for deltastar
    double deltastar = GSL_NAN;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_compute_deltastar(
        b.collocation_point(b.n()-1), rho_u.get(),
        &deltastar, dB.get(), b.bw, b.dbw));

    // Check against good value
    // Good value taken from White, Fluid Mechanics, 4th Edition eqn (7.31).
    // This tolerance is admittedly larger than I would like.
    BOOST_CHECK_CLOSE(1.721, deltastar, 0.015);
}

BOOST_AUTO_TEST_CASE( blasius_theta )
{
    // Prepare the Blasius velocity profile as coefficients on basis
    // Pretend that density is uniformly two throughout profile.  Yes, two.
    shared_array<double> rho_u(new double[b.n()]);
    shared_array<double> u    (new double[b.n()]);
    for (int i = 0; i < b.n(); ++i) {
        u[i] = gsl_spline_eval(blasius_u, b.collocation_point(i), accel);
        rho_u[i] = 2*u[i];
    }
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.solve(1, rho_u.get(), 1, b.n()));
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.solve(1,     u.get(), 1, b.n()));

    // Prepare integration working storage
    shared_ptr<gsl_matrix> dB(
            gsl_matrix_alloc(b.k(), 1),
            gsl_matrix_free);
    BOOST_REQUIRE(dB);

    // Integrate for theta
    double theta  = GSL_NAN;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_compute_theta(
        b.collocation_point(b.n()-1), rho_u.get(), u.get(), &theta, dB.get(),
        b.bw, b.dbw));

    // Check against good value
    // Good value taken from White, Fluid Mechanics, 4th Edition eqn (7.31).
    // This tolerance is admittedly larger than I would like.
    BOOST_CHECK_CLOSE(0.664, theta, 0.018);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( qoi )



// Test data taken from
// https://svn.ices.utexas.edu/repos/pecos/turbulence/heatshield_bl/trunk/laminar/wall.dat
// for dstat ~ 2.5007933304570438.
struct wall_type : public suzerain_bl_local
{
    wall_type()
    {
        std::fill_n(reinterpret_cast<double *>(this),
                    sizeof(*this)/sizeof(double),
                    std::numeric_limits<double>::quiet_NaN());
        a      =  843.98854040857316;
        gamma  =  1.3527314891502726;
        mu     =  4.8813735289922836e-05;
        Pr     =  0.65543907074081864;
        p__x   =  -1231.0557214607243;
        rho    =  0.018546113877544138;
        T      =  1391.8731000995472;
        u      =  0.0023004630235243829;
        u__x   =  3117.168135000341;
        u__y   =  333239.70652878482;
        v      =  0.27873944160103337;
    }
};
static const wall_type wall;

// Test data taken from
// https://svn.ices.utexas.edu/repos/pecos/turbulence/heatshield_bl/trunk/laminar/edge.dat
// for dstat ~ 2.5007933304570438.
struct edge_type : public suzerain_bl_local
{
    edge_type()
    {
        std::fill_n(reinterpret_cast<double *>(this),
                    sizeof(*this)/sizeof(double),
                    std::numeric_limits<double>::quiet_NaN());
        a      =  1956.3958663603282;
        gamma  =  1.4083595370046604;
        mu     =  0.00016225807140364439;
        Pr     =  0.80552596752550176;
        p__x   =  -1860.4745416352641;
        rho    =  0.0037307784953988427;
        T      =  5840.4009311559321;
        u      =  1396.7581826189837;
        u__x   =  499.26639968207024;
        u__y   =  4634.7550551015656;
        v      =  -41.964917478845166;
    }
};
static const edge_type edge;


// Test data taken from
// https://svn.ices.utexas.edu/repos/pecos/turbulence/heatshield_bl/trunk/laminar/scenario.dat
// for dstat ~ 2.5007933304570438.
BOOST_AUTO_TEST_CASE( compute_viscous )
{
    suzerain_bl_viscous viscous;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS,
                        suzerain_bl_compute_viscous(&wall, &viscous));

    size_t cnt = 0; // Tracks if all quantities were tested
    const double tol = GSL_SQRT_DBL_EPSILON;
    BOOST_CHECK_CLOSE(viscous.tau_w,    16.266674822587674,     tol); ++cnt;
    BOOST_CHECK_CLOSE(viscous.u_tau,    29.615763373028074,     tol); ++cnt;
    BOOST_CHECK_CLOSE(viscous.delta_nu, 8.8872252594154481e-05, tol); ++cnt;
    BOOST_CHECK_EQUAL(cnt, sizeof(viscous)/sizeof(viscous.tau_w));
}

// Test data taken from
// https://svn.ices.utexas.edu/repos/pecos/turbulence/heatshield_bl/trunk/laminar/scenario.dat
// for dstat ~ 2.5007933304570438.
BOOST_AUTO_TEST_CASE( compute_qoi )
{
    suzerain_bl_viscous viscous;            // Answers from just above
    std::fill_n(reinterpret_cast<double *>(&viscous),
                sizeof(viscous)/sizeof(double),
                std::numeric_limits<double>::quiet_NaN());
    viscous.tau_w    = 16.266674822587674;
    viscous.u_tau    = 29.615763373028074;
    viscous.delta_nu = 8.8872252594154481e-05;

    suzerain_bl_thick thick;                // Answers from scenario.dat
    std::fill_n(reinterpret_cast<double *>(&thick),
                sizeof(thick)/sizeof(double),
                std::numeric_limits<double>::quiet_NaN());
    thick.delta     = 0.046525678647201738;
    thick.deltastar = 0.0044202837563584669;
    thick.theta     = 0.0059005327804110153;

    const double code_Ma = 1;               // Data from a dimensional code
    const double code_Re = 1;               // Ditto
    suzerain_bl_qoi qoi;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_compute_qoi(
            code_Ma, code_Re, &wall, &viscous, &edge, &thick, &qoi));

    size_t cnt = 0; // Tracks if all quantities were tested
    const double tol = GSL_SQRT_DBL_EPSILON;
    BOOST_CHECK_CLOSE(qoi.beta,         -0.50556278312573966,   tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.Cf,           0.0044697874046917899,  tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.gamma_e,      1.4083595370046604,     tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.K_e,          1.1130040269123832e-05, tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.K_s,          24.849095784497852,     tol*1e4); ++cnt; // FIXME Why?
    BOOST_CHECK_CLOSE(qoi.K_w,          3.3483624867674195e-06, tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.Lambda_n,     5.3212990115980237,     tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.Ma_e,         0.71394455827465408,    tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.p_ex,         -0.011892537649319856,  tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.Pr_w,         0.65543907074081864,    tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.Re_delta,     1494.1943713234461,     tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.Re_deltastar, 141.95952214875473,     tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.Re_theta,     189.49842591559681,     tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.ratio_rho,    edge.rho / wall.rho,    tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.ratio_nu,       (edge.mu / edge.rho)
                                        / (wall.mu / wall.rho), tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.ratio_T,      4.1960728537236802,     tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.shapefactor,  0.7491329886401481,     tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.v_wallplus,   0.0094118607746200931,  tol); ++cnt;
    BOOST_CHECK_EQUAL(cnt, sizeof(qoi)/sizeof(qoi.beta));
}

BOOST_AUTO_TEST_SUITE_END()
