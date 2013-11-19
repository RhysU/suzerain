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

// Helper method testing if the solution satisfies the governing equation.
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
        double u      = s->state[i].u;
        double up_rhs = - (u / s->state[i].R)
                      * (2 + Ma0*Ma0*(gam0-1) - Ma0*Ma0*(gam0-1)*u*u)
                      / (2 + Ma0*Ma0*(gam0-1) - Ma0*Ma0*(gam0+1)*u*u);
        double pp_rhs = -Ma0*Ma0*s->state[i].rho*u*s->state[i].up;
        double rp_rhs = -Ma0*Ma0*u*s->state[i].up/s->state[i].a2;

        // Relative test between u' and RHS as (u' - RHS) may be large
        gsl_test_rel(s->state[i].up, up_rhs, tol, "%s: up res[%d] ", who, i);

        // Relative test between p' and RHS
        gsl_test_rel(s->state[i].pp, pp_rhs, tol, "%s: pp res[%d] ", who, i);

        // Relative test between (log rho)' and RHS
        gsl_test_rel(s->state[i].rhop / s->state[i].rho, rp_rhs, tol,
                     "%s: rhop res[%d] ", who, i);

    }
}

// Helper method if the solution satisfies Euler in polar coordinates.
// If not, then the computed solutions aren't particularly useful to us...
static
void check_radial_euler_residual(
    const char * who,
    const suzerain_radial_nozzle_solution * const s,
    const double tol)
{
    for (size_t i = 0; i < s->size; ++i) {
        const double R   = s->state[i].R;
        const double u   = s->state[i].u,   up   = s->state[i].up;
        const double p   = s->state[i].p,   pp   = s->state[i].pp;
        const double rho = s->state[i].rho, rhop = s->state[i].rhop;
        (void) p;

        // Continuity residual from baseflow.tex writeup, nozzle.m
        gsl_test_abs(0, rho*up + u*rhop + rho*u/R,
                     tol, "%s: up mass[%d] ", who, i);

        // Momentum residual from baseflow.tex writeup, nozzle.m
        gsl_test_rel(pp, -rho*u*up,
                     tol, "%s: up momentum[%d] ", who, i);

        // Energy residual from baseflow.tex writeup, nozzle.m
        gsl_test_rel(s->state[i].a2,
                     1 + s->Ma0*s->Ma0*(s->gam0-1)/2*(1-u*u),
                     tol, "%s: up energy[%d] ", who, i);
    }
}

// Test conversion to Cartesian primitive state by checking the residual
// of the Euler spatial operator computed in primitive form.  Residual
// should be small as the radial nozzle problem produces a stationary flow.
static
void check_cartesian_primitive(
    const char * who,
    const suzerain_radial_nozzle_solution * const s,
    const double Ma,
    const double tol)
{
    // Compute primitive state at the outermost radius for Ma, s->gam0
    double rho, u, v, p, rho_xi, u_xi, v_xi, p_xi, rho_y , u_y , v_y , p_y;
    suzerain_radial_nozzle_cartesian_primitive(s, s->size-1, Ma, &rho, &u, &v,
            &p, &rho_xi, &u_xi, &v_xi, &p_xi, &rho_y, &u_y, &v_y, &p_y);

    // In nondimensional primitive variables with Ma dependence, that is
    const double U  [4] = { rho,    u,    v,    p    }; (void) U;
    const double U_x[4] = { rho_xi, u_xi, v_xi, p_xi };
    const double U_y[4] = { rho_y,  u_y,  v_y,  p_y  };
    // when a_0 != u_0, the 2D Euler equations take the
    // form \partial_t U + A \partial_x U + B \partial_y U = 0 with
    const double A[4][4] = { { u, rho,          0, 0     },
                             { 0, u,            0, 1/rho },
                             { 0, 0,            u, 0     },
                             { 0, s->gam0/Ma/Ma*p, 0, u     } };
    // and
    const double B[4][4] = { { v, 0, rho,              0     },
                             { 0, v, 0,                0     },
                             { 0, 0, v,                1/rho },
                             { 0, 0, s->gam0/Ma/Ma*p,  v     } };
    // which may be seen writeups/notebooks/Giles_BC_Nondimensional.nb
    // under "Sanity check the linearized evolution equation".  All together,
    double U_t[4] = { 0, 0, 0, 0};
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            U_t[i] += A[i][j] * U_x[j] + B[i][j] * U_y[j];
        }
    }
    // where we've now got the spatial residual of the Euler equations in U_t.

////// FIXME #2503
////// As we should observe a steady solution, now check against zero:
////const double R = s->state[s->size-1].R;
////gsl_test_abs(U_t[0], 0.0, tol, "%s: rho_t at R = %g", who, R);
////gsl_test_abs(U_t[1], 0.0, tol, "%s: u_t   at R = %g", who, R);
////gsl_test_abs(U_t[2], 0.0, tol, "%s: v_t   at R = %g", who, R);
////gsl_test_abs(U_t[3], 0.0, tol, "%s: p_t   at R = %g", who, R);
}

// Second subsonic test from writeups/notebooks/nozzle.m
// Beware the slightly different argument order relative to that code
static
void test_subsonic()
{
    const double R[]  = {1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.};
    const size_t N    = sizeof(R)/sizeof(R[0]);
    const double Ma0  = 1.0;
    const double gam0 = 1.4;
    const double rho1 = 1.0;
    const double u1   = -1/Ma0 + GSL_SQRT_DBL_EPSILON;
    const double p1   = 1.0;
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

    // Expected results computed by writeups/notebooks/nozzle.m using Octave
    double tol = GSL_SQRT_DBL_EPSILON;
    suzerain_radial_nozzle_state fin = s->state[N-1];
    gsl_test_rel(fin.R,    R[N-1],             tol, "%s final R   ", __func__);
    gsl_test_rel(fin.u,   -0.332008421365410,  tol, "%s final u   ", __func__);
    gsl_test_rel(fin.a2,   1.17795408162849,   tol, "%s final a2  ", __func__);
    gsl_test_rel(fin.up,   0.183142130216718,  tol, "%s final up  ", __func__);
    // Octave and GSL RKF45 adaptive control differs, hence lower tolerances...
    tol = sqrt(tol);
    gsl_test_rel(fin.rho,  1.50598587203636,   tol, "%s final rho ", __func__);
    gsl_test_rel(fin.p,    1.55284442659135,   tol, "%s final p   ", __func__);
    gsl_test_rel(fin.rhop, 0.0777373796467802, tol, "%s final rhop", __func__);
    gsl_test_rel(fin.pp,   0.0915710636498338, tol, "%s final pp  ", __func__);

    // Does the pointwise solution satisfy the radial governing equations?
    check_radial_nozzle_residual(__func__, s, GSL_SQRT_DBL_EPSILON);
    check_radial_euler_residual (__func__, s, GSL_SQRT_DBL_EPSILON);

    // Test edge Mach and pressure gradient parameter computations
    // Expected results by writeups/notebooks/nozzle_qoi.m for delta = sqrt(3)
    const double Mae  = suzerain_radial_nozzle_qoi_Mae (s, s->size-1);
    const double pexi = suzerain_radial_nozzle_qoi_pexi(s, s->size-1);
    gsl_test_rel(Mae,   0.152951916589060, tol, "%s qoi_Mae ", __func__);
    gsl_test_rel(pexi, -1.91086402711018,  tol, "%s qoi_pexi", __func__);

    // Can the results be converted to a Cartesian frame correctly?
    check_cartesian_primitive(__func__, s, 1.5, GSL_SQRT_DBL_EPSILON);
    // TODO check_cartesian_conserved

    free(s);
}

// Supersonic test from writeups/notebooks/nozzle.m
// Beware the slightly different argument order relative to that code
static
void test_supersonic()
{
    const double R[]  = {1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.};
    const size_t N    = sizeof(R)/sizeof(R[0]);
    const double Ma0  = 1.0;
    const double gam0 = 1.4;
    const double rho1 = 1.0;
    const double u1   = 1/Ma0 + GSL_SQRT_DBL_EPSILON;
    const double p1   = 1.0;
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

    // Expected results computed by writeups/notebooks/nozzle.m using Octave
    double tol = GSL_SQRT_DBL_EPSILON;
    suzerain_radial_nozzle_state fin = s->state[N-1];
    gsl_test_rel(fin.R,     R[N-1],             tol, "%s final R   ", __func__);
    gsl_test_rel(fin.u,     1.71679859932923,   tol, "%s final u   ", __func__);
    gsl_test_rel(fin.a2,    0.610520513868242,  tol, "%s final a2  ", __func__);
    gsl_test_rel(fin.up,    0.224261011681836,  tol, "%s final up  ", __func__);
    // Octave and GSL RKF45 adaptive control differs, hence lower tolerances...
    tol = sqrt(tol);
    gsl_test_rel(fin.rho,   0.291239756946575,  tol, "%s final rho ", __func__);
    gsl_test_rel(fin.p,     0.412719892947868,  tol, "%s final p   ", __func__);
    gsl_test_rel(fin.rhop, -0.183663783311630,  tol, "%s final rhop", __func__);
    gsl_test_rel(fin.pp,   -0.112130507369524,  tol, "%s final pp  ", __func__);

    // Does the pointwise solution satisfy the radial governing equations?
    check_radial_nozzle_residual(__func__, s, GSL_SQRT_DBL_EPSILON);
    check_radial_euler_residual (__func__, s, GSL_SQRT_DBL_EPSILON);

    // Test edge Mach and pressure gradient parameter computations
    // Expected results by writeups/notebooks/nozzle_qoi.m for delta = sqrt(3)
    const double Mae  = suzerain_radial_nozzle_qoi_Mae (s, s->size-1);
    const double pexi = suzerain_radial_nozzle_qoi_pexi(s, s->size-1);
    gsl_test_rel(Mae,   1.09859906253134,  tol, "%s qoi_Mae ", __func__);
    gsl_test_rel(pexi, -0.452506737297551, tol, "%s qoi_pexi", __func__);

    // Can the results be converted to a Cartesian frame correctly?
    check_cartesian_primitive(__func__, s, 1.5, GSL_SQRT_DBL_EPSILON);
    // TODO check_cartesian_conserved

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
    test_supersonic();

    exit(gsl_test_summary());
}
