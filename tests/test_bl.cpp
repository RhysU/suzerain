/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2013-2014 Rhys Ulerich
 * Copyright (C) 2013-2014 The PECOS Development Team
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

    double code_Ma;
    gsl_spline * const blasius_u;
    gsl_spline * const blasius_v;
    gsl_interp_accel * const accel;
    bspline      b;
    bsplineop    op;
    bsplineop_lu lu;
    shared_array<double> u;
    shared_array<double> v;
    shared_array<double> ke;
    shared_array<double> H0;

    BlasiusFixture()
        : code_Ma(1.5)
        , blasius_u(suzerain_blasius_u(Re_x))
        , blasius_v(suzerain_blasius_v(Re_x))
        , accel(gsl_interp_accel_alloc())
        , b(k,
            bspline::from_breakpoints(),
            suzerain_blasius_extended_size,
            blasius_u->x + /*UGLY HACK uses Nnegative from blasius.c */10)
        , op(b, 0, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE)
        , lu(op)
        , u (new double[b.n()])
        , v (new double[b.n()])
        , ke(new double[b.n()])
        , H0(new double[b.n()])
    {
        lu.factor_mass(op);

        // Prepare Blasius profile information as coefficients
        for (int i = 0; i < b.n(); ++i) {
            u[i]  = gsl_spline_eval(blasius_u, b.collocation_point(i), accel);
            v[i]  = gsl_spline_eval(blasius_v, b.collocation_point(i), accel);
            ke[i] = (u[i]*u[i] + v[i]*v[i]) / 2;
            H0[i] = u[i] + code_Ma*code_Ma*ke[i]; // h = u => delta2 = deltaH0
        }

        BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.solve(1, u .get(), 1, b.n()));
        BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.solve(1, v .get(), 1, b.n()));
        BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.solve(1, ke.get(), 1, b.n()));
        BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.solve(1, H0.get(), 1, b.n()));
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

BOOST_AUTO_TEST_CASE( blasius_find_edge99 )
{
    // Prepare working storage
    shared_ptr<gsl_matrix> dB(gsl_matrix_alloc(b.k(), 1), gsl_matrix_free);

    // Find edge using kinetic energy profile
    double location = GSL_NAN;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_find_edge99(
        u.get(),
        gsl_bspline_breakpoint(0,              b.bw),
        gsl_bspline_breakpoint(b.bw->nbreak-1, b.bw),
        &location, dB.get(), b.bw, b.dbw));

    // For the classical Blasius layer truncated in the usual fashion,
    // we should see something around 4.91 when x = Re_x = 1 per
    // https://en.wikipedia.org/wiki/Boundary-layer_thickness.
    BOOST_REQUIRE_GT(location, 4.90);
    BOOST_REQUIRE_LT(location, 4.92);
}

BOOST_AUTO_TEST_CASE( blasius_delta1 )
{
    // Prepare integration working storage
    shared_ptr<gsl_vector> Bk(gsl_vector_alloc(b.k()), gsl_vector_free);
    shared_ptr<gsl_integration_glfixed_table> tbl(
            gsl_integration_glfixed_table_alloc((b.k() + 1)/2),
            gsl_integration_glfixed_table_free);

    // Integrate for delta1
    double delta1 = GSL_NAN;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_displacement_thickness(
        u.get(), &delta1, Bk.get(), b.bw, tbl.get()));

    // Good value found using Octave's trapz on extended data up to eta = 13.6
    // Compare Schlichting 8th Edition of 5.0*0.34 on p162 per Blasius
    BOOST_CHECK_SMALL((delta1 - 1.72189451530768), GSL_SQRT_DBL_EPSILON);
}

BOOST_AUTO_TEST_CASE( blasius_delta2 )
{
    // Prepare integration working storage
    shared_ptr<gsl_vector> Bk(gsl_vector_alloc(b.k()), gsl_vector_free);
    shared_ptr<gsl_integration_glfixed_table> tbl(
            gsl_integration_glfixed_table_alloc(b.k()),
            gsl_integration_glfixed_table_free);

    // Integrate for delta2
    // Pretend density is nondimensionally one so rhou == u
    double delta2  = GSL_NAN;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_momentum_thickness(
        u.get(), u.get(), &delta2, Bk.get(), b.bw, tbl.get()));

    // "Good" value found using Octave's trapz on extended data generated by
    // running 'blasius 38.5 0.025' against https://github.com/RhysU/blasius
    // Compare Schlichting 8th Edition of 5.0*0.13 on p162 per Blasius
    const double tol      = 0.11;
    const double expected = 0.664045493818590;
    BOOST_REQUIRE_CLOSE(delta2, expected, tol);
}

BOOST_AUTO_TEST_CASE( blasius_delta3 )
{
    // Prepare integration working storage
    shared_ptr<gsl_vector> Bk(gsl_vector_alloc(b.k()), gsl_vector_free);
    shared_ptr<gsl_integration_glfixed_table> tbl(
            gsl_integration_glfixed_table_alloc(b.k()),
            gsl_integration_glfixed_table_free);

    // Integrate for delta3
    // Pretend density is nondimensionally one so rhou == u
    double delta3  = GSL_NAN;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_energy_thickness(
        ke.get(), u.get(), &delta3, Bk.get(), b.bw, tbl.get()));

    // "Good" value found using Octave's trapz on extended data generated by
    // running 'blasius 38.5 0.05' against https://github.com/RhysU/blasius.
    //
    // Compare Schlichting 8th Edition of 5.0*0.20 on p162 per Blasius and
    // notice that this result is considerably thicker.  That is because we
    // included the wall-normal velocity in the kinetic energy computation.
    BOOST_CHECK_CLOSE(delta3, 1.30347780585580, 0.11);

    // The enthalpy thickness is mathematically identical to the energy
    // thickness when one substitutes u^2 for H0.  This relies on the
    // fact that u^2 goes to zero at the wall so that the reference
    // "viscous H0" value is zero there.
    double deltaH0 = GSL_NAN;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_enthalpy_thickness(
        ke.get(), u.get(), &deltaH0, Bk.get(), b.bw, tbl.get()));
    BOOST_CHECK_CLOSE(deltaH0, delta3, GSL_SQRT_DBL_EPSILON);

    // To convince you that the wall-normal velocity makes a difference,
    // lets remove it from the kinetic energy and try again...
    for (int i = 0; i < b.n(); ++i) {
        double u = gsl_spline_eval(blasius_u, b.collocation_point(i), accel);
        ke[i] = u*u/2;                   // Neglect v^2/2
    }
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.solve(1, u .get(), 1, b.n()));
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_energy_thickness(
        ke.get(), u.get(), &delta3, Bk.get(), b.bw, tbl.get()));
    BOOST_CHECK_CLOSE(delta3, 1.0, 4.5); // Convinced?
}

BOOST_AUTO_TEST_SUITE_END()


typedef BlasiusFixture<4,1000> fixture_four_thousand;
BOOST_FIXTURE_TEST_CASE( blasius_find_edge, fixture_four_thousand )
{
    // Prepare working storage
    shared_ptr<gsl_matrix> dB(gsl_matrix_alloc(b.k(), 3), gsl_matrix_free);

    // Find edge using kinetic energy profile
    double location = GSL_NAN;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_find_edge(
        ke.get() /* \approx H_0 */,
        gsl_bspline_breakpoint(0,              b.bw),
        gsl_bspline_breakpoint(b.bw->nbreak-1, b.bw),
        &location, dB.get(), b.bw, b.dbw));

    // Thickness from eyeballing results computed in Octave
    //   source notebooks/blasius.m
    //   Re=1000; plot(blasius_y(Re), blasius_kepp(Re),
    //                 blasius_y(Re), zeros(size(blasius_eta)))
    //   ylim([-eps eps])
    // where machine epsilon arises from knowing the internals of
    // suzerain_bl_find_edge's invocation of suzerain_bspline_crossing.
    // This admittedly is a lousy test highly dependent on input data.
    BOOST_REQUIRE_GT(location, 0.40);
    BOOST_REQUIRE_LT(location, 0.80);
}

BOOST_AUTO_TEST_SUITE( challenging_edge )

// Observed from temporal BL as described in Redmine ticket #2980.
struct ChallengingFixture {

    enum { k = 8, ndof = 64, nbreak = ndof + 2 - k };
    static const double breakpoints[nbreak];
    static const double coeffs_H0[ndof];
    bspline b;

    ChallengingFixture()
        : b(k, bspline::from_breakpoints(), nbreak, breakpoints)
    {}
};

const double ChallengingFixture::breakpoints[ChallengingFixture::nbreak] = {
    0, 0.0011002998920519591, 0.0023220284341718411, 0.003678499490802567,
    0.0051844664196041279, 0.0068562729740557771, 0.008712018914787345,
    0.010771741498605492, 0.013057614040625332, 0.015594162752043728,
    0.018408503037206314, 0.021530596380496814, 0.02499352885616668,
    0.028833812140274917, 0.033091707678748472, 0.037811574351668931,
    0.043042239550508121, 0.048837393028061227, 0.055256002162267004,
    0.062362746363366028, 0.070228467213387491, 0.078930629519176065,
    0.08855378674389569, 0.099190042215263841, 0.11093949505120637,
    0.12391065685977676, 0.1382208219345924, 0.15399636987091303,
    0.17137297528737894, 0.19049569470740324, 0.21151889573520699,
    0.23460598862326432, 0.25992891542249441, 0.28766734748754019,
    0.31800753864890252, 0.35114077946464239, 0.38726139836483631,
    0.42656425906481599, 0.46924171130479353, 0.51547996476424474,
    0.56545487482890033, 0.61932715449011133, 0.6772370593931627,
    0.73929863270616791, 0.80559364203647155, 0.87616539005293781,
    0.95101263060831287, 1.0300838686543241, 1.1132723597480392,
    1.2004121474761016, 1.2912754786903955, 1.3855719119255943,
    1.48294938046729, 1.5829973877764068, 1.6852524024193267,
    1.789205389231395, 1.8943112736435577, 2
};

const double ChallengingFixture::coeffs_H0[ChallengingFixture::ndof] = {
    2.5000000000000004, 2.5015268558121302, 2.5047511184390201,
    2.509865814549654, 2.5170901003392805, 2.5266738233951562,
    2.5389032165290208, 2.5541080243959753, 2.5710934224380222,
    2.5900825791456299, 2.6113297926985801, 2.6351252186340504,
    2.6617999958327667, 2.691731967235095, 2.7253513848974569,
    2.7631467337657138, 2.8056693015066374, 2.853536455740735,
    2.9074308003817704, 2.9680949365740878, 3.0363165842710358,
    3.1129040399612911, 3.1986435945198108, 3.2942416211251215,
    3.4002408767356096, 3.5169239687324505, 3.6441973128644234,
    3.7814966125704705, 3.9277263930628576, 4.0813327357629516,
    4.2405460597019839, 4.4039157109946876, 4.5710229588406426,
    4.7432084909816705, 4.9237502816225893, 5.1172299686449563,
    5.3280738148537354, 5.5592065199489209, 5.8114409875075408,
    6.0838872200362246, 6.37447894808966, 6.6802250259972453,
    6.9969191398022135, 7.318828612602398, 7.6383865570797553,
    7.9460052984357548, 8.2301150859323808, 8.4770726599191288,
    8.6717623945827498, 8.7976982750921344, 8.838355344139142,
    8.7775120658681285, 8.6017626830522964, 8.3011520648003838,
    7.8715826377277249, 7.3179825383914974, 6.6535754596501695,
    5.9045881732060925, 5.2264162816800335, 4.6618845833331015,
    4.2266019762251634, 3.9171887070027975, 3.7206427930339285,
    3.6249768899238313
};

BOOST_FIXTURE_TEST_CASE( find_edge, ChallengingFixture )
{
    // Sanity check the provided data amount
    BOOST_REQUIRE_EQUAL(b.n(), static_cast<int>(ndof));

    // Prepare working storage
    shared_ptr<gsl_matrix> dB(gsl_matrix_alloc(b.k(), 3), gsl_matrix_free);

    // Find edge using H0 already expressed as coefficients
    double location = GSL_NAN;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_find_edge(
        coeffs_H0,
        gsl_bspline_breakpoint(0,              b.bw),
        gsl_bspline_breakpoint(b.bw->nbreak-1, b.bw),
        &location, dB.get(), b.bw, b.dbw));

    // Ensure we found an actual edge in this profile
    BOOST_CHECK(!(boost::math::isnan)(location));
}

BOOST_AUTO_TEST_SUITE_END();

typedef BlasiusFixture<4,10000> fixture_four_ten_thousand;

BOOST_FIXTURE_TEST_CASE( blasius_thicknesses_reynolds,
                         fixture_four_ten_thousand )
{
    // Prepare, beyond the fixture, uniform density of 0.5 by scaling u coeffs
    shared_array<double> rhou(new double[b.n()]);
    for (int i = 0; i < b.n(); ++i) {
        rhou[i] = u[i] / 2;
    }

    // Compute a bunch of thickness-related quantities
    // Thickness from eyeballing results computed in Octave:
    //   source notebooks/blasius.m
    //   Re=10000; plot(blasius_y(Re), blasius_kepp(Re),
    //                  blasius_y(Re), zeros(size(blasius_eta)))
    // Integrals found using Octave's trapz against this data.
    // Again, this is a very poor test for thick.delta.
    // Notice on this dataset that delta2 and deltaH0 are the same thing!
    size_t cnt = 0; // Tracks if all quantities were tested
    suzerain_bl_thicknesses thick;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_compute_thicknesses(
        H0.get(), ke.get(), rhou.get(), u.get(),
        &thick, b.bw, b.dbw));
    BOOST_CHECK_GT   (thick.delta,   0.12);                        ++cnt;
    BOOST_CHECK_CLOSE(thick.delta1,  0.0172085683613221,  0.01  ); ++cnt;
    BOOST_CHECK_CLOSE(thick.delta2,  0.00664045493818580, 0.0105); ++cnt;
    BOOST_CHECK_CLOSE(thick.delta3,  0.0104435139441593,  0.01  ); ++cnt;
    BOOST_CHECK_GT   (thick.deltaH0, thick.delta2);                ++cnt;
    /* Quantity delta99 not tested */                              ++cnt;
    BOOST_CHECK_EQUAL(cnt, sizeof(thick)/sizeof(thick.delta));

    // Compute the same values but now with a baseflow-friendly routine using a
    // stride trick to compute with a constant baseflow.  Here we test for
    // consistency against simpler routine rather than against golden value.
    const double tol = GSL_SQRT_DBL_EPSILON;
    suzerain_bl_thicknesses baseflow;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS,
                        suzerain_bl_compute_thicknesses_baseflow(
                 H0.get(),    ke.get(), rhou.get(), u.get(),
        0,       &H0[b.n()-1], &rhou[b.n()-1], &u[b.n()-1], &v[b.n()-1],
        &baseflow, b.bw, b.dbw));
    cnt = 0;
    BOOST_CHECK_EQUAL(thick.delta,   baseflow.delta);        ++cnt;
    BOOST_CHECK_CLOSE(thick.delta1,  baseflow.delta1,  tol); ++cnt;
    BOOST_CHECK_CLOSE(thick.delta2,  baseflow.delta2,  tol); ++cnt;
    BOOST_CHECK_CLOSE(thick.delta3,  baseflow.delta3,  tol); ++cnt;
    BOOST_CHECK_CLOSE(thick.deltaH0, baseflow.deltaH0, tol); ++cnt;
    BOOST_CHECK_CLOSE(thick.delta99, baseflow.delta99, tol); ++cnt;
    BOOST_CHECK_EQUAL(cnt, sizeof(thick)/sizeof(thick.delta));

    // Prepare expected Reynolds number results
    const double code_Re = 777;
    suzerain_bl_local edge;
    {
        edge.mu   = 5;
        edge.rho  = rhou[b.n()-1] / u[b.n()-1];
        edge.u    = u[b.n()-1];
        edge.y    = thick.delta;
    }
    suzerain_bl_local edge99 = edge; // Pretend edge and edge99
    thick.delta99 = thick.delta;     // are the same location.
    suzerain_bl_reynolds expected_Re;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_compute_reynolds(
                code_Re, &edge, &edge99, &thick, &expected_Re));

    // Test the baseflow-based Reynolds number computations using
    // the simpler, uniform inviscid flow as a sanity check.
    cnt = 0;
    suzerain_bl_reynolds reynolds;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS,
                        suzerain_bl_compute_reynolds_baseflow(
           code_Re,
           H0.get(), ke.get(), rhou.get(), u.get(),
        0, &H0[b.n()-1], &rhou[b.n()-1], &u[b.n()-1], &v[b.n()-1],
        &edge, &edge99, &reynolds, b.bw));
    BOOST_CHECK_CLOSE(expected_Re.delta,   reynolds.delta,   tol); ++cnt;
    BOOST_CHECK_CLOSE(expected_Re.delta1,  reynolds.delta1,  tol); ++cnt;
    BOOST_CHECK_CLOSE(expected_Re.delta2,  reynolds.delta2,  tol); ++cnt;
    BOOST_CHECK_CLOSE(expected_Re.delta3,  reynolds.delta3,  tol); ++cnt;
    BOOST_CHECK_CLOSE(expected_Re.deltaH0, reynolds.deltaH0, tol); ++cnt;
    BOOST_CHECK_CLOSE(expected_Re.delta99, reynolds.delta99, tol); ++cnt;
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
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_compute_viscous(
                /*code_Re*/1, &wall, &viscous));

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
    thick.delta   = const_cast<double&>(edge.y) = 0.046525678647201738;
    thick.delta1  = 0.0044202837563584669;
    thick.delta2  = 0.0059005327804110153;

    const double code_Ma = 1;               // Data from a dimensional code
    const double code_Re = 1;               // Ditto

    // Pretend edge and edge99 are the same location for testing convenience
    suzerain_bl_local edge99 = edge;
    thick.delta99 = thick.delta;

    suzerain_bl_reynolds reynolds;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_compute_reynolds(
            code_Re, &edge, &edge99, &thick, &reynolds));

    size_t cnt = 0; // Tracks if all quantities were tested
    const double tol = GSL_SQRT_DBL_EPSILON;
    const double meh = sqrt(tol);
    BOOST_CHECK_CLOSE(reynolds.delta,    1494.1943713234461,   tol); ++cnt;
    BOOST_CHECK_CLOSE(reynolds.delta1,   141.95952214875473,   tol); ++cnt;
    BOOST_CHECK_CLOSE(reynolds.delta2,   189.49842591559681,   tol); ++cnt;
    /* Quantity reynolds.delta3 not tested */                        ++cnt;
    /* Quantity reynolds.deltaH0 not tested */                       ++cnt;
    BOOST_CHECK_CLOSE(reynolds.delta99,  1494.1943713234461,   tol); ++cnt;
    BOOST_CHECK_EQUAL(cnt, sizeof(reynolds)/sizeof(reynolds.delta));

    suzerain_bl_qoi qoi;
    BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bl_compute_qoi(
            code_Ma, code_Re, &wall, &viscous, &edge, &thick, &qoi));

    cnt = 0; // Tracks if all quantities were tested
    BOOST_CHECK_CLOSE(qoi.cf,          0.0044697874046917899,  tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.gamma_e,     1.4083595370046604,     tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.Ma_e,        0.71394455827465408,    tol); ++cnt;
    /* Quantity Ma_tau  not tested */                                ++cnt;
    BOOST_CHECK_CLOSE(qoi.Pr_w,        0.65543907074081864,    tol); ++cnt;
    /* Quantity Bq not tested */                                     ++cnt;
    BOOST_CHECK_CLOSE(qoi.ratio_rho,   edge.rho / wall.rho,    tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.ratio_nu,      (edge.mu / edge.rho)
                                       / (wall.mu / wall.rho), tol); ++cnt;
    BOOST_CHECK_CLOSE(qoi.ratio_T,     4.1960728537236802,     tol); ++cnt;
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
