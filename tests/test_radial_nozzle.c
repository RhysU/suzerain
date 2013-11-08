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
#include <gsl/gsl_test.h>

#include <suzerain/common.h>
#include <suzerain/error.h>

// Second subsonic test from writeups/notebooks/nozzle1.m
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
    gsl_test_rel(fin.rho,  1.61808760864002,   tol, "%s final rho ", __func__);
    gsl_test_rel(fin.p,    1.68133559416742,   tol, "%s final p   ", __func__);
    gsl_test_rel(fin.rhop, 0.0835239513576762, tol, "%s final rhop", __func__);
    gsl_test_rel(fin.pp,   0.0983873794154783, tol, "%s final pp  ", __func__);

    // Test edge Mach and pressure gradient parameter computations
    // Expected results by writeups/notebooks/nozzle_qoi.m for delta = sqrt(3)
    const double Mae  = suzerain_radial_nozzle_qoi_Mae (s, s->size-1);
    const double pexi = suzerain_radial_nozzle_qoi_pexi(s, s->size-1);
    gsl_test_rel(Mae,   0.152951916589060, tol, "%s qoi_Mae ", __func__);
    gsl_test_rel(pexi, -1.91086402711018,  tol, "%s qoi_pexi", __func__);

    free(s);
}

// Supersonic test from writeups/notebooks/nozzle1.m
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
    gsl_test_rel(fin.rho,   0.188891739106648,  tol, "%s final rho ", __func__);
    gsl_test_rel(fin.p,     0.336677203999661,  tol, "%s final p   ", __func__);
    gsl_test_rel(fin.rhop, -0.119120314492264,  tol, "%s final rhop", __func__);
    gsl_test_rel(fin.pp,   -0.0727253956159634, tol, "%s final pp  ", __func__);

    // Test edge Mach and pressure gradient parameter computations
    // Expected results by writeups/notebooks/nozzle_qoi.m for delta = sqrt(3)
    const double Mae  = suzerain_radial_nozzle_qoi_Mae (s, s->size-1);
    const double pexi = suzerain_radial_nozzle_qoi_pexi(s, s->size-1);
    gsl_test_rel(Mae,   1.09859906253134,  tol, "%s qoi_Mae ", __func__);
    gsl_test_rel(pexi, -0.452506737297551, tol, "%s qoi_pexi", __func__);

    free(s);
}

// TODO suzerain_radial_nozzle_cartesian_primitive

// TODO suzerain_radial_nozzle_cartesian_conserved

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
