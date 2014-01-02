/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2014 Rhys Ulerich
 * Copyright (C) 2012-2014 The PECOS Development Team
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

#include <suzerain/htstretch.h>

#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_test.h>

#include <suzerain/common.h>
#include <suzerain/error.h>

#define CASE(f, x, y, z, e) \
    gsl_test_abs(suzerain_ ## f((x),(y),(z)), (e), GSL_DBL_EPSILON, \
        #f "(%g, %g, %g)", (x), (y), (z) );

static
void test_htstretch1()
{
    /* Check boundaries for varying delta values */
    for (int i = 1; i < 10; i += 2) {
        CASE(htstretch1, (double) i, 2., 0., 0.);
        CASE(htstretch1, (double) i, 2., 2., 1.);
    }

    /* Spot check evaluation and derivatives from Mathematica N[expr, 16] */
    CASE(htstretch1,        7., 3., 1./ 9., +1.130109424675245e-6);
    CASE(htstretch1_ddelta, 7., 3., 1./ 9., -2.053316509625048e-6);
    CASE(htstretch1_dL,     7., 3., 1./ 9., -4.827691460278840e-7);
    CASE(htstretch1_dx,     7., 3., 1./ 9., +0.00001303476694275287);
    CASE(htstretch1,        3., 2., 1./11., +0.001553191024452532);
    CASE(htstretch1_ddelta, 3., 2., 1./11., -0.002512942959289942);
    CASE(htstretch1_dL,     3., 2., 1./11., -0.0008866151849355177);
    CASE(htstretch1_dx,     3., 2., 1./11., +0.01950553406858139);

    /* Spot check evaluation near wall from Mathematica N[expr, 16] */
    CASE(htstretch1,        4., 1., 1.0e-03, +5.388915112675453e-6);
    CASE(htstretch1_ddelta, 4., 1., 1.0e-03, -9.425182395066187e-6);
    CASE(htstretch1_dL,     4., 1., 1.0e-03, -5.410484940505296e-6);
    CASE(htstretch1_dx,     4., 1., 1.0e-03, +0.005410484940505296);
    CASE(htstretch1,        4., 1., 1.0e-08, +5.367402865013899e-11);
    CASE(htstretch1_ddelta, 4., 1., 1.0e-08, -9.392957375935146e-11);
    CASE(htstretch1_dL,     4., 1., 1.0e-08, -5.367403079566020e-11);
    CASE(htstretch1_dx,     4., 1., 1.0e-08, +0.005367403079566020);

    /* Spot check evaluation against analytic delta = 0 degeneration */
    CASE(htstretch1,        0., 2., 1./ 11., +1./22.);
    CASE(htstretch1_ddelta, 0., 2., 1./ 11., +0.);
    CASE(htstretch1_dL,     0., 2., 1./ 11., -1./44.);
    CASE(htstretch1_dx,     0., 2., 1./ 11., +1./2.);
}

static
void test_htstretch2()
{
    /* Check boundaries for varying delta values */
    for (int i = 1; i < 10; i += 2) {
        CASE(htstretch2, (double) i, 2., 0., 0.);
        CASE(htstretch2, (double) i, 2., 2., 1.);
    }

    /* Spot check evaluation and derivatives from Mathematica N[expr, 16] */
    CASE(htstretch2,        7., 3., 1./ 9., +0.0006192752099952831);
    CASE(htstretch2_ddelta, 7., 3., 1./ 9., -0.0005055862703645726);
    CASE(htstretch2_dL,     7., 3., 1./ 9., -0.0002643833885364380);
    CASE(htstretch2_dx,     7., 3., 1./ 9., +0.007138351490483826);
    CASE(htstretch2,        3., 2., 1./11., +0.01541983224264333);
    CASE(htstretch2_ddelta, 3., 2., 1./11., -0.009494629614989877);
    CASE(htstretch2_dL,     3., 2., 1./11., -0.008679927388636600);
    CASE(htstretch2_dx,     3., 2., 1./11., +0.1909584025500052);

    /* Spot check evaluation near wall from Mathematica N[expr, 16] */
    CASE(htstretch2,        4., 1., 1.0e-03, +0.0001471408880498442);
    CASE(htstretch2_ddelta, 4., 1., 1.0e-03, -0.0001102915457626845);
    CASE(htstretch2_dL,     4., 1., 1.0e-03, -0.0001477088973243801);
    CASE(htstretch2_dx,     4., 1., 1.0e-03, +0.1477088973243801);
    CASE(htstretch2,        4., 1., 1.0e-08, +1.465742869555286e-9);
    CASE(htstretch2_ddelta, 4., 1., 1.0e-08, -1.100290869880845e-9);
    CASE(htstretch2_dL,     4., 1., 1.0e-08, -1.465742926075948e-9);
    CASE(htstretch2_dx,     4., 1., 1.0e-08, +0.1465742926075948);

    /* Spot check evaluation against analytic delta = 0 degeneration */
    CASE(htstretch2,        0., 3., 1./9., +1./27.);
    CASE(htstretch2_ddelta, 0., 3., 1./9., +0);
    CASE(htstretch2_dL,     0., 3., 1./9., -1./81.);
    CASE(htstretch2_dx,     0., 3., 1./9., 1./3.);
}

static
void test_htstretch1_find_delta()
{
    /* Test case found by direct evaluation in Mathematica */
    const double L        = 2.0;
    const double crit_x   = 0.5;
    const double crit_val = 0.04012773950620386L;

    const double expected = 2.345678901234567;
    const int maxiter     = 50;

    double delta;

    /* Attempt to recover expected */
    gsl_test_int(suzerain_htstretch1_find_delta(
                    L, crit_x, crit_val, GSL_DBL_EPSILON, maxiter, &delta),
                 SUZERAIN_SUCCESS,
            "htstretch1_find_delta returns SUZERAIN_SUCCESS");
    gsl_test_abs(delta, expected, 100*GSL_DBL_EPSILON,
            "htstretch1_find_delta(delta, %g, %g, %g)", L, crit_x, crit_val);
}

static
void test_htstretch2_find_delta()
{
    /* Test case found by direct evaluation in Mathematica */
    const double L        = 2.0;
    const double crit_x   = 0.5;
    const double crit_val = 0.03141345703210717;

    const double expected = 6.7890123456789012;
    const int maxiter     = 50;

    double delta;

    /* Attempt to recover expected */
    gsl_test_int(suzerain_htstretch2_find_delta(
                    L, crit_x, crit_val, GSL_DBL_EPSILON, maxiter, &delta),
                 SUZERAIN_SUCCESS,
            "htstretch2_find_delta returns SUZERAIN_SUCCESS");
    gsl_test_abs(delta, expected, 100*GSL_DBL_EPSILON,
            "htstretch2_find_delta(delta, %g, %g, %g)", L, crit_x, crit_val);
}

int main(int argc, char **argv)
{
    SUZERAIN_UNUSED(argc);
    SUZERAIN_UNUSED(argv);
    gsl_ieee_env_setup();

    if (argc > 1) {
        gsl_test_verbose(1);
    }

    test_htstretch1();
    test_htstretch2();

    test_htstretch1_find_delta();
    test_htstretch2_find_delta();

    exit(gsl_test_summary());
}
