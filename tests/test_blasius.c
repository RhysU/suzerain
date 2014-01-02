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

#include <suzerain/blasius.h>

#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_test.h>

#include <suzerain/common.h>
#include <suzerain/countof.h>

static void test_extended();
static void test_wall_derivative();

int
main(int argc, char **argv)
{
    SUZERAIN_UNUSED(argc);
    SUZERAIN_UNUSED(argv);

    gsl_ieee_env_setup();

    test_extended();
    test_wall_derivative();

    exit(gsl_test_summary());
}

/* Check that extended Blasius data matches Ganapol well in overlap */
static void test_extended()
{
    for (size_t i = 0;
         i < SUZERAIN_COUNTOF(suzerain_blasius_ganapol_eta);
         ++i) {

        gsl_test_rel(suzerain_blasius_extended_eta[i],
                     suzerain_blasius_ganapol_eta[i],
                     GSL_SQRT_DBL_EPSILON/10,
                     "eta[%d]", i);

        gsl_test_rel(suzerain_blasius_extended_f[i],
                     suzerain_blasius_ganapol_f[i],
                     GSL_SQRT_DBL_EPSILON/10, "f  [%d] at eta = %g",
                     i, suzerain_blasius_ganapol_eta[i]);

        gsl_test_rel(suzerain_blasius_extended_fp[i],
                     suzerain_blasius_ganapol_fp[i],
                     GSL_SQRT_DBL_EPSILON/10, "fp [%d] at eta = %g",
                     i, suzerain_blasius_ganapol_eta[i]);

        gsl_test_rel(suzerain_blasius_extended_fpp[i],
                     suzerain_blasius_ganapol_fpp[i],
                     GSL_SQRT_DBL_EPSILON/10, "fpp[%d] at eta = %g",
                     i, suzerain_blasius_ganapol_eta[i]);

    }
}

/* Check the derivative of u at the wall */
static void test_wall_derivative()
{
    gsl_spline * const blasius_u = suzerain_blasius_u(/*Re_x = */1);
    gsl_interp_accel * const accel = gsl_interp_accel_alloc();

    const double u  = gsl_spline_eval      (blasius_u, 0, accel);
    const double du = gsl_spline_eval_deriv(blasius_u, 0, accel);

    gsl_test_abs(u,  0.0,                 GSL_DBL_EPSILON,      "u[0.0]");
    gsl_test_abs(du, 0.33205733621519630, GSL_SQRT_DBL_EPSILON, "du[0.0]");

    gsl_interp_accel_free(accel);
    gsl_spline_free(blasius_u);
}
