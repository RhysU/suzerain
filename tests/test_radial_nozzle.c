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
    const double R[] = {1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};
    const size_t N   = sizeof(R)/sizeof(R[0]);
    suzerain_radial_nozzle_solution * s = suzerain_radial_nozzle_solver(
            1.0, 1.4, 1, -1 + GSL_SQRT_DBL_EPSILON, 1, R, N);

    // Expected values computed by writeups/notebooks/nozzle1.m via Octave
    const double tol = GSL_SQRT_DBL_EPSILON;
    suzerain_radial_nozzle_state fin = s->state[N-1];
    gsl_test_abs(fin.R,   R[N-1],           tol, "%s final R",   __func__);
    gsl_test_abs(fin.u,  -0.33200842136541, tol, "%s final u",   __func__);
    gsl_test_abs(fin.a2,  1.17795408162849, tol, "%s final a2",  __func__);
    gsl_test_abs(fin.rho, 1.61933117083039, tol, "%s final rho", __func__);
    gsl_test_abs(fin.p,   1.68254700209343, tol, "%s final p",   __func__);

    free(s);
}

// From writeups/notebooks/nozzle1.m test cases
// Beware the slightly different argument order relative to that code
static
void test_supersonic()
{
    const double R[] = {1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};
    const size_t N   = sizeof(R)/sizeof(R[0]);
    suzerain_radial_nozzle_solution * s = suzerain_radial_nozzle_solver(
            1.0, 1.4, 1, 1 + GSL_SQRT_DBL_EPSILON, 1, R, N);

    // TODO
    gsl_test_abs(1.0, 1.0, GSL_DBL_EPSILON, "NOP");

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

    // FIXME test_subsonic();
    // FIXME test_supersonic();

    exit(gsl_test_summary());
}
