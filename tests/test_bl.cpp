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

// A fixture exposing Blasius profile information at local Re_x
// B-spline workspaces are prepared using Blasius-related breakpoints
template <int k, int Re_x>
struct BlasiusFixture {

    gsl_spline * const blasius_u;
    gsl_spline * const blasius_v;
    gsl_interp_accel * const accel;
    bspline      b;
    bsplineop    op;
    bsplineop_lu lu;

    BlasiusFixture()
        : blasius_u(suzerain_blasius_u(Re_x))
        , blasius_v(suzerain_blasius_v(Re_x))
        , accel(gsl_interp_accel_alloc())
          // "size - 1" occurs because suzerain_blasius_extended_*
          // contains extrapolations at "infinity" which perturb
          // the computation of thick.deltaH far below.
        , b(k, bspline::from_breakpoints(), blasius_u->size - 1, blasius_u->x)
        , op(b, 0, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE)
        , lu(op)
    {
        lu.factor_mass(op);
    }

    ~BlasiusFixture()
    {
        gsl_interp_accel_free(accel);
        gsl_spline_free(blasius_v);
        gsl_spline_free(blasius_u);
    }

};

typedef BlasiusFixture<2,1> linear_fixture;
BOOST_FIXTURE_TEST_SUITE(bl_compute_thick_linear, linear_fixture)

BOOST_AUTO_TEST_CASE( blasius_delta1 )
{
    // Prepare the Blasius velocity profile as coefficients on basis
    // For a linear B-spline basis, the collocation points are the breakpoints
    shared_array<double> u(new double[b.n()]);
    std::memcpy(u.get(), blasius_u->y, b.n()*sizeof(double));
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.solve(1, u.get(), 1, b.n()));

    // Prepare integration working storage
    shared_ptr<gsl_matrix> dB(gsl_matrix_alloc(b.k(), 1), gsl_matrix_free);
    BOOST_REQUIRE(dB);

    // Integrate for delta1
    double delta1 = GSL_NAN;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_displacement_thickness(
        u.get(), &delta1, dB.get(), b.bw, b.dbw));

    // Good value found using Octave's trapz on extended data up to eta = 13.6
    BOOST_CHECK_SMALL((delta1 - 1.72189451530768), GSL_SQRT_DBL_EPSILON);
}

BOOST_AUTO_TEST_CASE( blasius_delta2 )
{
    // Prepare the Blasius velocity profile as coefficients on basis
    // For a linear B-spline basis, the collocation points are the breakpoints
    shared_array<double> u(new double[b.n()]);
    std::memcpy(u.get(), blasius_u->y, b.n()*sizeof(double));
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.solve(1, u.get(), 1, b.n()));

    // Prepare integration working storage
    shared_ptr<gsl_matrix> dB(gsl_matrix_alloc(b.k(), 1), gsl_matrix_free);
    BOOST_REQUIRE(dB);

    // Integrate for delta2
    // Pretend density is nondimensionally one so rho_u == u
    double delta2  = GSL_NAN;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_momentum_thickness(
        u.get(), u.get(), &delta2, dB.get(), b.bw, b.dbw));

    // "Good" value found using Octave's trapz on extended data generated by
    // running 'blasius 38.5 0.025' against https://github.com/RhysU/blasius
    const double tol      = 0.11;
    const double expected = 0.664045493818590;
    BOOST_REQUIRE_CLOSE(delta2, expected, tol);

    // The enthalpy thickness is mathematically identical to the
    // momentum thickness when one substitutes H0 for u.  So sanity
    // check check that routine using the same test machinery.
    double deltaH  = GSL_NAN;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_enthalpy_thickness(
        u.get(), u.get(), &deltaH, dB.get(), b.bw, b.dbw));
    BOOST_CHECK_CLOSE(deltaH, expected, tol);
}

BOOST_AUTO_TEST_CASE( blasius_delta3 )
{
    // Prepare the Blasius velocity profile as coefficients on basis
    // For a linear B-spline basis, the collocation points are the breakpoints
    shared_array<double> u(new double[b.n()]);
    std::memcpy(u.get(), blasius_u->y, b.n()*sizeof(double));
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.solve(1, u.get(), 1, b.n()));

    // Prepare integration working storage
    shared_ptr<gsl_matrix> dB(gsl_matrix_alloc(b.k(), 1), gsl_matrix_free);
    BOOST_REQUIRE(dB);

    // Integrate for delta3
    // Pretend density is nondimensionally one so rho_u == u
    double delta3  = GSL_NAN;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_energy_thickness(
        u.get(), u.get(), &delta3, dB.get(), b.bw, b.dbw));

    // "Good" value found using Octave's trapz on extended data generated by
    // running 'blasius 38.5 0.025' against https://github.com/RhysU/blasius
    BOOST_CHECK_CLOSE(delta3, 1.04430629471857, 0.11);
}

BOOST_AUTO_TEST_SUITE_END()


typedef BlasiusFixture<4,1000> fixture_four_thousand;
BOOST_FIXTURE_TEST_CASE( blasius_find_edge, fixture_four_thousand )
{
    // Prepare B-spline coefficients for Blasius profile kinetic energy
    // Kinetic energy should be evaluated consistently with velocity fits
    shared_array<double> ke(new double[b.n()]);
    for (int i = 0; i < b.n(); ++i) {
        const double u = gsl_spline_eval(
                blasius_u, b.collocation_point(i), accel);
        const double v = gsl_spline_eval(
                blasius_v, b.collocation_point(i), accel);
        ke[i] = (u*u + v*v) / 2;
    }
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.solve(1, ke.get(), 1, b.n()));

    // Prepare working storage
    shared_ptr<gsl_matrix> dB(gsl_matrix_alloc(b.k(), 3), gsl_matrix_free);
    BOOST_REQUIRE(dB);

    // Find edge using kinetic energy profile
    double location = GSL_NAN;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_find_edge(
        ke.get() /* \approx H_0 */, &location, dB.get(), b.bw, b.dbw));

    // Thickness from eyeballing results computed in Octave
    //   source writeups/notebooks/blasius.m
    //   Re=1000; plot(blasius_y(Re), blasius_kepp(Re),
    //                 blasius_y(Re), zeros(size(blasius_eta)))
    //   ylim([-eps eps])
    // where machine epsilon arises from knowing the internals of
    // suzerain_bl_find_edge's invocation of suzerain_bspline_crossing.
    // This admittedly is a lousy test highly dependent on input data.
    BOOST_REQUIRE_GT(location, 0.40);
    BOOST_REQUIRE_LT(location, 0.80);
}

typedef BlasiusFixture<4,10000> fixture_four_ten_thousand;
BOOST_FIXTURE_TEST_CASE( blasius_compute_thicknesses, fixture_four_ten_thousand)
{
    // Prepare Blasius profile information on the collocation points
    shared_array<double> u    (new double[b.n()]);
    shared_array<double> v    (new double[b.n()]);
    shared_array<double> ke   (new double[b.n()]);
    shared_array<double> rho_u(new double[b.n()]);
    for (int i = 0; i < b.n(); ++i) {
        u[i]     = gsl_spline_eval(blasius_u, b.collocation_point(i), accel);
        v[i]     = gsl_spline_eval(blasius_v, b.collocation_point(i), accel);
        ke[i]    = (u[i]*u[i] + v[i]*v[i]) / 2;
        rho_u[i] = u[i] / 2;                     // Density uniformly 0.5

        // std::cerr << b.collocation_point(i) << '\t'  // If you need it...
        //           << u[i]                   << '\t'
        //           << v[i]                   << '\t'
        //           << ke[i]                  << '\t'
        //           << rho_u[i]               << '\n';
    }

    // Convert to B-spline coefficients
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.solve(1, u    .get(), 1, b.n()));
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.solve(1, v    .get(), 1, b.n()));
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.solve(1, ke   .get(), 1, b.n()));
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.solve(1, rho_u.get(), 1, b.n()));

    // Compute a bunch of thickness-related quantities
    // Thickness from eyeballing results computed in Octave:
    //   source writeups/notebooks/blasius.m
    //   Re=10000; plot(blasius_y(Re), blasius_kepp(Re),
    //                  blasius_y(Re), zeros(size(blasius_eta)))
    // Integrals found using Octave's trapz against this data.
    // Again, this is a very poor test for thick.delta.
    // Notice on this dataset that delta3 and deltaH are the same thing!
    size_t cnt = 0;
    suzerain_bl_thicknesses thick;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_compute_thicknesses(
        ke.get() /* \approx H_0 */, rho_u.get(), u.get(), &thick, b.bw, b.dbw));
    BOOST_CHECK_GT   (thick.delta,  0.12);                        ++cnt;
    BOOST_CHECK_CLOSE(thick.delta1, 0.0172085683613221,  0.01  ); ++cnt;
    BOOST_CHECK_CLOSE(thick.delta2, 0.00664045493818580, 0.0105); ++cnt;
    BOOST_CHECK_CLOSE(thick.delta3, 0.0104430629471855,  0.01  ); ++cnt;
    BOOST_CHECK_CLOSE(thick.deltaH, 0.0104435139441593,  0.01  ); ++cnt;
    BOOST_CHECK_EQUAL(cnt, sizeof(thick)/sizeof(thick.delta));
}


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
        rho    =  0.018546113877544138;
        T      =  1391.8731000995472;
        u      =  0.0023004630235243829;
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
        rho    =  0.0037307784953988427;
        T      =  5840.4009311559321;
        u      =  1396.7581826189837;
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
BOOST_AUTO_TEST_CASE( compute_qoi_and_pg )
{
    suzerain_bl_viscous viscous;            // Answers from just above
    std::fill_n(reinterpret_cast<double *>(&viscous),
                sizeof(viscous)/sizeof(double),
                std::numeric_limits<double>::quiet_NaN());
    viscous.tau_w    = 16.266674822587674;
    viscous.u_tau    = 29.615763373028074;
    viscous.delta_nu = 8.8872252594154481e-05;

    suzerain_bl_thicknesses thick;                // Answers from scenario.dat
    std::fill_n(reinterpret_cast<double *>(&thick),
                sizeof(thick)/sizeof(double),
                std::numeric_limits<double>::quiet_NaN());
    thick.delta  = 0.046525678647201738;
    thick.delta1 = 0.0044202837563584669;
    thick.delta2 = 0.0059005327804110153;

    const double code_Ma = 1;               // Data from a dimensional code
    const double code_Re = 1;               // Ditto
    suzerain_bl_qoi qoi;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_compute_qoi(
            code_Ma, code_Re, &wall, &viscous, &edge, &thick, &qoi));

    size_t cnt = 0; // Tracks if all quantities were tested
    const double tol = GSL_SQRT_DBL_EPSILON;
    const double meh = sqrt(tol);
    BOOST_CHECK_CLOSE(qoi.cf,          0.0044697874046917899,  tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.gamma_e,     1.4083595370046604,     tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.Ma_e,        0.71394455827465408,    tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.Pr_w,        0.65543907074081864,    tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.ratio_rho,   edge.rho / wall.rho,    tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.ratio_nu,      (edge.mu / edge.rho)
                                       / (wall.mu / wall.rho), tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.ratio_T,     4.1960728537236802,     tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.Re_delta,    1494.1943713234461,     tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.Re_delta1,   141.95952214875473,     tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.Re_delta2,   189.49842591559681,     tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.shapefactor, 0.7491329886401481,     tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.v_wallplus,  0.0094118607746200931,  tol); ++cnt;
    BOOST_CHECK_EQUAL(cnt, sizeof(qoi)/sizeof(qoi.cf));

    cnt = 0; // Tracks if all quantities were tested
    const double edge_p__x   =  -1860.4745416352641;  // Corresponds to edge
    const double edge_u__x   =  499.26639968207024;
    suzerain_bl_pg pg;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_compute_pg(
            code_Ma, code_Re, &wall, &viscous,
            &edge, edge_p__x, edge_u__x, &thick, &pg));
    BOOST_CHECK_CLOSE(pg.Clauser,      -0.50556278312573966,   tol); ++cnt;
    BOOST_CHECK_CLOSE(pg.Lambda_n,     5.3212990115980237,     tol); ++cnt;
    BOOST_CHECK_CLOSE(pg.Launder_e,    1.1130040269123832e-05, tol); ++cnt;
    BOOST_CHECK_CLOSE(pg.Launder_w,    3.3483624867674195e-06, tol); ++cnt;
    BOOST_CHECK_CLOSE(pg.p_ex,         -0.011892537649319856,  tol); ++cnt;
    BOOST_CHECK_CLOSE(pg.Pohlhausen,   24.849095784497852,     meh); ++cnt;
    BOOST_CHECK_EQUAL(cnt, sizeof(pg)/sizeof(pg.Clauser));
}

BOOST_AUTO_TEST_SUITE_END()
