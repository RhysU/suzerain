/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2012 The PECOS Development Team
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

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.h>
#pragma hdrstop
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_test.h>
#include <suzerain/filterop.h>

static void test_filterop();

int
main(int argc, char **argv)
{
    SUZERAIN_UNUSED(argc);
    SUZERAIN_UNUSED(argv);

    gsl_ieee_env_setup();

    test_filterop();

    exit(gsl_test_summary());
}

// static void alloc_workspaces(
//         size_t k, size_t nbreak, const double * breakpoints,
//         gsl_bspline_workspace **w, gsl_bspline_deriv_workspace **dw,
//         gsl_matrix **scratch)
// {
//         gsl_vector_const_view bv
//             = gsl_vector_const_view_array(breakpoints, nbreak);
//         *w = gsl_bspline_alloc(k, nbreak);
//         gsl_bspline_knots(&bv.vector, *w);
//         *dw = gsl_bspline_deriv_alloc(k);
//         *scratch = gsl_matrix_alloc(k, k);
// }

// static void free_workspaces(
//         gsl_bspline_workspace **w, gsl_bspline_deriv_workspace **dw,
//         gsl_matrix **scratch)
// {
//     gsl_matrix_free(*scratch);
//     *scratch = NULL;
//     gsl_bspline_deriv_free(*dw);
//     *dw = NULL;
//     gsl_bspline_free(*w);
//     *w = NULL;
// }

static void test_filterop()
{
    // FIXME Test it
    // alloc_workspaces(3, sizeof(b)/sizeof(b[0]), b, &w, &dw, &scratch);
    // free_workspaces(&w, &dw, &scratch);
}
