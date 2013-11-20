/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2013 Rhys Ulerich
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

#include <suzerain/radial_nozzle.h>

#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_vector.h>

#include <suzerain/common.h>
#include <suzerain/error.h>

// Helper seiing if the solution satisfies the governing equation.
// If the governing equation is not satisfied, then all bets are off...
static
void check_radial_nozzle_residual(
    const char * who,
    const suzerain_radial_nozzle_solution * const s,
    const double tol)
{
    const double Ma0  = s->Ma0;
    const double gam0 = s->gam0;
    for (size_t i = 0; i < s->size; ++i) {
        const double R   = s->state[i].R;
        const double u   = s->state[i].u,   up   = s->state[i].up;
        const double p   = s->state[i].p,   pp   = s->state[i].pp;
        const double rho = s->state[i].rho, rhop = s->state[i].rhop;
        const double a2  = s->state[i].a2;
        (void) p;

        // Relative test between u' and RHS as (u' - RHS) may be large
        double up_rhs = - (u / R)
                      * (2 + Ma0*Ma0*(gam0-1) - Ma0*Ma0*(gam0-1)*u*u)
                      / (2 + Ma0*Ma0*(gam0-1) - Ma0*Ma0*(gam0+1)*u*u);
        gsl_test_rel(up, up_rhs, tol, "%s: up res[%d] ", who, i);

        // Relative test between p' and RHS (momentum equation)
        gsl_test_rel(pp, -Ma0*Ma0*rho*u*up,
                     tol, "%s: pp res[%d] ", who, i);

        // Relative test between (log rho)' and RHS
        gsl_test_rel(rhop / rho, -Ma0*Ma0*u*up/a2,
                     tol, "%s: rhop res[%d] ", who, i);

        // Energy residual from baseflow.tex writeup, nozzle.m
        gsl_test_rel(a2, 1 + Ma0*Ma0*(gam0-1)/2*(1-u*u),
                     tol, "%s: res_energy[%d] ", who, i);

        // Continuity polar residual from baseflow.tex writeup, nozzle.m
        gsl_test_abs(0, rho*up + u*rhop + rho*u/R,
                     tol, "%s: res_mass[%d] ", who, i);
    }
}

// Helper seeing if the solution satisfies the ideal gas equation of state.
// This is only approximately satisfied because of constant stagnation energy
// assumption used in the derivation.  Larger radii should provide better
// agreement.
static
void check_ideal_gas_approximation(
    const char * who,
    const suzerain_radial_nozzle_solution * const s,
    const double tol)
{
    for (size_t i = 0; i < s->size; ++i) {
        gsl_test_rel(s->state[i].rho * s->state[i].a2,
                     s->gam0 * s->state[i].p,
                     tol, "%s: ideal_EOS[%d] at %g", who, i, s->state[i].R);

    }
}

// Test conversion to Cartesian primitive state by checking the residual
// of the Euler spatial operator computed in primitive form.  Residual
// should be small as the radial nozzle problem produces a stationary flow.
//
// Technically rho*a2 and not gam0*p should appear in A and B because the
// nozzle problem employed a constant stagnation energy approximation rather
// than requiring the ideal gas EOS to hold.  However, for large enough radii
// the approximation should be sane.  but the ideal gas EOS should
// approximately hold.  Test that with check_ideal_gas_approximation() first to
// discover tolerances.
static
void check_euler_primitive(
    const char * who,
    const suzerain_radial_nozzle_solution * const s,
    const double Ma,
    const double tol)
{
    for (size_t i = 0; i < s->size; ++i) {

        // Compute primitive state and derivatives for Ma, s->gam0
        double rho, u, v, p, rho_xi, u_xi, v_xi, p_xi, rho_y , u_y , v_y , p_y;
        suzerain_radial_nozzle_cartesian_primitive(s, i, Ma, &rho, &u, &v,
                &p, &rho_xi, &u_xi, &v_xi, &p_xi, &rho_y, &u_y, &v_y, &p_y);

        // In nondimensional primitive variables with Ma dependence, that is
        const double U  [4] = { rho,    u,    v,    p    }; (void) U;
        const double U_x[4] = { rho_xi, u_xi, v_xi, p_xi };
        const double U_y[4] = { rho_y,  u_y,  v_y,  p_y  };
        // when a_0 != u_0, the 2D Euler equations take the
        // form \partial_t U + A \partial_x U + B \partial_y U = 0 with
        const double A[4][4] = { { u, rho,          0, 0           },
                                 { 0, u,            0, 1/rho/Ma/Ma },
                                 { 0, 0,            u, 0           },
                                 { 0, s->gam0*p, 0, u              } };
        // and
        const double B[4][4] = { { v, 0, rho,        0           },
                                 { 0, v, 0,          0           },
                                 { 0, 0, v,          1/rho/Ma/Ma },
                                 { 0, 0, s->gam0*p,  v           } };
        // Computing the two matrix-vector products,
        double AU_x[4] = { 0, 0, 0, 0 };
        double BU_y[4] = { 0, 0, 0, 0 };
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                AU_x[i] += A[i][j] * U_x[j];
                BU_y[i] += B[i][j] * U_y[j];
            }
        }
        // we've now got the spatial residual of the Euler equations.
        // As we should observe a steady solution, check these balance.
        const double R = s->state[i].R;
        gsl_test_rel(AU_x[0], -BU_y[0], tol, "%s: rho_t Ma=%g, R=%g", who,Ma,R);
        gsl_test_rel(AU_x[1], -BU_y[1], tol, "%s: u_t   Ma=%g, R=%g", who,Ma,R);
        gsl_test_rel(AU_x[2], -BU_y[2], tol, "%s: v_t   Ma=%g, R=%g", who,Ma,R);
        gsl_test_rel(AU_x[3], -BU_y[3], tol, "%s: p_t   Ma=%g, R=%g", who,Ma,R);
    }
}

// Test conversion to Cartesian conserved state by checking the residual of the
// Euler spatial operator computed in primitive form.  Residual should be small
// as the radial nozzle problem produces a stationary flow.  The same ideal gas
// EOS caveats from check_euler_primitive apply here, except magnified since
// the EOS is used to form total energy here.
static
void check_euler_conserved(
    const char * who,
    const suzerain_radial_nozzle_solution * const s,
    const double Ma,
    const double tol)
{
    for (size_t i = 0; i < s->size; ++i) {

        // Compute conserved state, pressure, and derivatives for Ma, s->gam0
        double r,    ru,    rv,    rE,    p;
        double r_xi, ru_xi, rv_xi, rE_xi, p_xi;
        double r_y,  ru_y,  rv_y,  rE_y,  p_y;
        suzerain_radial_nozzle_cartesian_conserved(s, i, Ma, &r, &ru, &rv,
                &rE, &p, &r_xi, &ru_xi, &rv_xi, &rE_xi, &p_xi, &r_y, &ru_y,
                &rv_y, &rE_y, &p_y);

        // Compute fluxes using conserved quantities and pressure
        const double f_xi[4] = {   ru_xi
                               ,   p_xi/Ma/Ma + ru/r*(2*ru_xi - ru*r_xi/r)
                               ,   (ru_xi*rv + ru*rv_xi - ru*rv*r_xi/r)/r
                               ,   (ru_xi*rE + ru*rE_xi - ru*rE*r_xi/r)/r
                                 + (ru_xi*p + ru*p_xi - ru*p*r_xi/r)/r/Ma/Ma };
        const double f_y[4] = {   rv_y
                              ,   (ru_y*rv + ru*rv_y - ru*rv*r_y/r)/r
                              ,   p_y/Ma/Ma  + rv/r*(2*rv_y  - rv*r_y /r)
                              ,   (rv_y *rE + rv*rE_y  - rv*rE*r_y /r)/r
                                + (rv_y *p + rv*p_y  - rv*p*r_y /r)/r/Ma/Ma };

        // Check absolute balance of the streamwise and wall-normal fluxes
        const double R = s->state[i].R;
        gsl_test_abs(f_xi[0], -f_y[0], tol, "%s: rho_t  Ma=%g, R=%g", who,Ma,R);
        gsl_test_abs(f_xi[1], -f_y[1], tol, "%s: rhou_t Ma=%g, R=%g", who,Ma,R);
        gsl_test_abs(f_xi[2], -f_y[2], tol, "%s: rhov_t Ma=%g, R=%g", who,Ma,R);
        gsl_test_abs(f_xi[3], -f_y[3], tol, "%s: rhoE_t Ma=%g, R=%g", who,Ma,R);
    }
}

// Subsonic verification test from notebooks/nozzle.m
// Beware the slightly different argument order relative to that code
static
void test_subsonic()
{
    const double R[]  = {10.0, 10.1, 10.2, 10.3, 10.4, 10.5};
    const size_t N    = sizeof(R)/sizeof(R[0]);
    const double Ma0  = 1.5;
    const double gam0 = 1.4;
    const double rho1 =  9./10;
    const double u1   = -2./ 7;
    const double p1   = rho1/gam0 *(1+(gam0-1)/2*Ma0*Ma0*(1-u1*u1));
    suzerain_radial_nozzle_solution * s = suzerain_radial_nozzle_solver(
            Ma0, gam0, rho1, u1, p1, R, N);

    // Check that the scenario parameters were stored correctly
    gsl_test_abs(s->Ma0,  Ma0,  GSL_DBL_EPSILON, "%s Ma0 ", __func__);
    gsl_test_abs(s->gam0, gam0, GSL_DBL_EPSILON, "%s gam0", __func__);

    // Check that the initial conditions were stored correctly
    suzerain_radial_nozzle_state ini = s->state[0];
    gsl_test_abs(ini.R,   R[0], GSL_DBL_EPSILON, "%s init R   ", __func__);
    gsl_test_abs(ini.u,   u1,   GSL_DBL_EPSILON, "%s init u   ", __func__);
    gsl_test_abs(ini.rho, rho1, GSL_DBL_EPSILON, "%s init rho ", __func__);
    gsl_test_abs(ini.p,   p1,   GSL_DBL_EPSILON, "%s init p   ", __func__);

    // Expected results computed by notebooks/nozzle.m using Octave
    double tol = GSL_SQRT_DBL_EPSILON;
    suzerain_radial_nozzle_state fin = s->state[N-1];
    gsl_test_rel(fin.R,    R[N-1],             tol, "%s final R   ", __func__);
    gsl_test_rel(fin.u,   -0.270256147749861,  tol, "%s final u   ", __func__);
    gsl_test_rel(fin.a2,   1.41713272657153,   tol, "%s final a2  ", __func__);
    gsl_test_rel(fin.up,   0.0291149687163293, tol, "%s final up  ", __func__);
    // Octave and GSL RKF45 adaptive control differs, hence lower tolerances...
    tol = sqrt(tol);
    gsl_test_rel(fin.rho,  0.906169799365092,  tol, "%s final rho ", __func__);
    gsl_test_rel(fin.p,    0.917259198936451,  tol, "%s final p   ", __func__);
    gsl_test_rel(fin.rhop, 0.0113207052839328, tol, "%s final rhop", __func__);
    gsl_test_rel(fin.pp,   0.0160429419457325, tol, "%s final pp  ", __func__);

    // Does the pointwise solution satisfy the governing equations?
    check_radial_nozzle_residual (__func__, s, 100*GSL_DBL_EPSILON);

    // Test edge Mach and pressure gradient parameter computations
    // Expected from notebooks/nozzle_qoi.m for delta = sqrt(10.5**2 - 10**2)
    const double Mae  = suzerain_radial_nozzle_qoi_Mae (s, s->size-1);
    const double pexi = suzerain_radial_nozzle_qoi_pexi(s, s->size-1);
    gsl_test_rel(Mae,   0.324318914847395, tol, "%s qoi_Mae ", __func__);
    gsl_test_rel(pexi, -0.362152908606146, tol, "%s qoi_pexi", __func__);

    // Do results approximately satisfy Cartesian Euler at various Ma?
    // Radii are sufficiently large that ideal gas EOS holds nicely.
    // The implicit conversion to conserved state incurs more error.
    check_ideal_gas_approximation(__func__, s,      100*GSL_DBL_EPSILON);
    check_euler_primitive        (__func__, s, Ma0, 100*GSL_DBL_EPSILON);
    check_euler_primitive        (__func__, s, 5.0, 100*GSL_DBL_EPSILON);
    check_euler_conserved        (__func__, s, Ma0, 0.025);
    check_euler_conserved        (__func__, s, 5.0, 0.025);

    free(s);
}

// Supersonic test from notebooks/nozzle.m
// Beware the slightly different argument order relative to that code
static
void test_supersonic1()
{
    const double R[]  = {1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.};
    const size_t N    = sizeof(R)/sizeof(R[0]);
    const double Ma0  = 1.0;
    const double gam0 = 1.4;
    const double rho1 = 1.0;
    const double u1   = 1/Ma0 + GSL_SQRT_DBL_EPSILON;
    const double p1   = rho1/gam0 *(1+(gam0-1)/2*Ma0*Ma0*(1-u1*u1));
    suzerain_radial_nozzle_solution * s = suzerain_radial_nozzle_solver(
            Ma0, gam0, rho1, u1, p1, R, N);

    // Check that the scenario parameters were stored correctly
    gsl_test_abs(s->Ma0,  Ma0,  GSL_DBL_EPSILON, "%s Ma0 ", __func__);
    gsl_test_abs(s->gam0, gam0, GSL_DBL_EPSILON, "%s gam0", __func__);

    // Check that the initial conditions were stored correctly
    suzerain_radial_nozzle_state ini = s->state[0];
    gsl_test_abs(ini.R,   R[0], GSL_DBL_EPSILON, "%s init R   ", __func__);
    gsl_test_abs(ini.u,   u1,   GSL_DBL_EPSILON, "%s init u   ", __func__);
    gsl_test_abs(ini.rho, rho1, GSL_DBL_EPSILON, "%s init rho ", __func__);
    gsl_test_abs(ini.p,   p1,   GSL_DBL_EPSILON, "%s init p   ", __func__);

    // Expected results computed by notebooks/nozzle.m using Octave
    double tol = GSL_SQRT_DBL_EPSILON;
    suzerain_radial_nozzle_state fin = s->state[N-1];
    gsl_test_rel(fin.R,     R[N-1],            tol, "%s final R   ", __func__);
    gsl_test_rel(fin.u,     1.71679859932442,  tol, "%s final u   ", __func__);
    gsl_test_rel(fin.a2,    0.610520513871545, tol, "%s final a2  ", __func__);
    gsl_test_rel(fin.up,    0.224261011684323, tol, "%s final up  ", __func__);
    // Octave and GSL RKF45 adaptive control differs, hence lower tolerances...
    tol = sqrt(tol);
    gsl_test_rel(fin.rho,   0.291239756609692, tol, "%s final rho ", __func__);
    gsl_test_rel(fin.p,     0.127005603875867, tol, "%s final p   ", __func__);
    gsl_test_rel(fin.rhop, -0.183663783096991, tol, "%s final rhop", __func__);
    gsl_test_rel(fin.pp,   -0.112130507235967, tol, "%s final pp  ", __func__);

    // Does the pointwise solution satisfy the governing equations?
    check_radial_nozzle_residual (__func__, s, 100*GSL_DBL_EPSILON);

    // Test edge Mach and pressure gradient parameter computations
    // Expected results by notebooks/nozzle_qoi.m for delta = sqrt(3)
    const double Mae  = suzerain_radial_nozzle_qoi_Mae (s, s->size-1);
    const double pexi = suzerain_radial_nozzle_qoi_pexi(s, s->size-1);
    gsl_test_rel(Mae,   1.09859906253134,  tol, "%s qoi_Mae ", __func__);
    gsl_test_rel(pexi, -0.452506737297551, tol, "%s qoi_pexi", __func__);

    // Do results approximately satisfy Cartesian Euler at various Ma?
    // Small radii case the ideal gas EOS to not be quite-so-satisfied.
    check_ideal_gas_approximation(__func__, s,          GSL_SQRT_DBL_EPSILON);
    check_euler_primitive        (__func__, s, Ma0, 100*GSL_SQRT_DBL_EPSILON);
    check_euler_primitive        (__func__, s, 5.0, 100*GSL_SQRT_DBL_EPSILON);
    // check_euler_conserved() is a disaster for these radii so do not bother.

    free(s);
}

// Supersonic with large radius so ideal gas EOS is better approximated
// relative to test_supersonic1().
static
void test_supersonic2()
{
    const double R[]  = {1000.0, 1000.1, 1000.2, 1000.3, 1000.4, 1000.5};
    const size_t N    = sizeof(R)/sizeof(R[0]);
    const double Ma0  = 1.1;
    const double gam0 = 1.4;
    const double rho1 = 1.0;
    const double u1   = 1/Ma0 + GSL_SQRT_DBL_EPSILON;
    const double p1   = rho1/gam0 *(1+(gam0-1)/2*Ma0*Ma0*(1-u1*u1));
    suzerain_radial_nozzle_solution * s = suzerain_radial_nozzle_solver(
            Ma0, gam0, rho1, u1, p1, R, N);

    // Does the pointwise solution satisfy the governing equations?
    check_radial_nozzle_residual (__func__, s, 100*GSL_DBL_EPSILON);

    // Do results approximately satisfy Cartesian Euler at various Ma?
    // Small radii case the ideal gas EOS to not be quite-so-satisfied.
    check_ideal_gas_approximation(__func__, s,      1e4*GSL_DBL_EPSILON);
    check_euler_primitive        (__func__, s, Ma0, GSL_SQRT_DBL_EPSILON);
    check_euler_primitive        (__func__, s, 5.0, GSL_SQRT_DBL_EPSILON);
    check_euler_conserved        (__func__, s, Ma0, 0.0015);
////check_euler_conserved        (__func__, s, 5.0, 0.025);

    free(s);
}

int main(int argc, char **argv)
{
    SUZERAIN_UNUSED(argc);
    SUZERAIN_UNUSED(argv);
    gsl_ieee_env_setup();

    if (argc > 1) {
        gsl_test_verbose(1);
    }

    test_subsonic();
    test_supersonic1();
    test_supersonic2();

    exit(gsl_test_summary());
}
