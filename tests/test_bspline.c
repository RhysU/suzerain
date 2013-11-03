/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2011, 2012, 2013 Rhys Ulerich
 * Copyright (C) 2012, 2013 The PECOS Development Team
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

#include <suzerain/bspline.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_test.h>

#include <suzerain/common.h>

static void alloc_workspaces(
        size_t k, size_t nbreak, const double * breakpoints,
        gsl_bspline_workspace **w, gsl_bspline_deriv_workspace **dw,
        gsl_matrix **scratch);
static void free_workspaces(
        gsl_bspline_workspace **w, gsl_bspline_deriv_workspace **dw,
        gsl_matrix **scratch);

static void test_integration_coefficients();
static void test_linear_combination();
static void test_linear_combination_complex();
static void test_crossing();
static void test_spacing_greville_abscissae();
static void test_spacing_breakpoints();
static void test_bspline_htstretch2_evdeltascale_greville_abscissae();
static void test_bspline_htstretch1_evdeltascale_greville_abscissae();
static void test_bspline_htstretch2_evdeltascale_breakpoints();
static void test_bspline_htstretch1_evdeltascale_breakpoints();

int
main(int argc, char **argv)
{
    SUZERAIN_UNUSED(argc);
    SUZERAIN_UNUSED(argv);

    gsl_ieee_env_setup();

    test_integration_coefficients();
    test_linear_combination();
    test_linear_combination_complex();
    test_crossing();
    test_spacing_greville_abscissae();
    test_spacing_breakpoints();
    test_bspline_htstretch2_evdeltascale_greville_abscissae();
    test_bspline_htstretch1_evdeltascale_greville_abscissae();
    test_bspline_htstretch2_evdeltascale_breakpoints();
    test_bspline_htstretch1_evdeltascale_breakpoints();

    exit(gsl_test_summary());
}

static void alloc_workspaces(
        size_t k, size_t nbreak, const double * breakpoints,
        gsl_bspline_workspace **w, gsl_bspline_deriv_workspace **dw,
        gsl_matrix **scratch)
{
        gsl_vector_const_view bv
            = gsl_vector_const_view_array(breakpoints, nbreak);
        *w = gsl_bspline_alloc(k, nbreak);
        gsl_bspline_knots(&bv.vector, *w);
        *dw = gsl_bspline_deriv_alloc(k);
        *scratch = gsl_matrix_alloc(k, k);
}

static void free_workspaces(
        gsl_bspline_workspace **w, gsl_bspline_deriv_workspace **dw,
        gsl_matrix **scratch)
{
    gsl_matrix_free(*scratch);
    *scratch = NULL;
    gsl_bspline_deriv_free(*dw);
    *dw = NULL;
    gsl_bspline_free(*w);
    *w = NULL;
}

/************************************************************************
 * Exact integration coefficients found using Mathematica 7 as follows:
 *  k = 2;
 *  b = {0, 1, 2, 3};
 *  n = Length[b] + k - 2;
 *  knots = Join[ ConstantArray[First[b], k - 1],
 *                b, ConstantArray[Last[b], k - 1]];
 *  B[i_, x_] := BSplineBasis[{k - 1, knots}, i, x]
 *  Table[Integrate[B[i, x], {x, Min[b], Max[b]}], {i, 0, n - 1}]
 ***********************************************************************/
static void test_integration_coefficients()
{
    const double b[] = { 0.0, 1.0, 2.0, 3.0 };
    gsl_bspline_workspace *w;
    gsl_bspline_deriv_workspace *dw;
    gsl_matrix *scratch;

    /* Piecewise linears with strided coefficient storage */
    {
        alloc_workspaces(2, sizeof(b)/sizeof(b[0]), b, &w, &dw, &scratch);
        const size_t inc = 2;
        const double e0[] = { .5, -555, 1, -555, 1, -555, .5, -555};
        const double e1[] = { -1, -555, 0, -555, 0, -555,  1, -555};
        gsl_vector *coeffs = gsl_vector_alloc(sizeof(e0)/sizeof(e0[0]));

        /* Zeroth derivative, integrating from low to high */
        gsl_vector_set_all(coeffs, -555.0);
        gsl_test(suzerain_bspline_integration_coefficients(
                0 , gsl_vector_ptr(coeffs,0), inc,
                -GSL_DBL_MAX, GSL_DBL_MAX, scratch, w, dw),
                "integration_coefficients 0 linear");
        for (size_t i = 0; i < w->n; ++i) {
            gsl_test_rel(gsl_vector_get(coeffs, 2*i), e0[2*i],
                    GSL_DBL_EPSILON*1000,
                    "integration_coefficients 0 linear %d value", 2*i);
            gsl_test(gsl_vector_get(coeffs, 2*i+1) != e0[2*i+1],
                    "integration_coefficients 0 linear %d stride", 2*i + 1);
        }

        /* First derivative, integrating from low to high */
        gsl_vector_set_all(coeffs, -555.0);
        gsl_test(suzerain_bspline_integration_coefficients(
                1 , gsl_vector_ptr(coeffs,0), inc,
                -GSL_DBL_MAX, GSL_DBL_MAX, scratch, w, dw),
                "integration_coefficients 1 linear");
        for (size_t i = 0; i < w->n; ++i) {
            gsl_test_rel(gsl_vector_get(coeffs, 2*i), e1[2*i],
                    GSL_DBL_EPSILON*1000,
                    "integration_coefficients 1 linear %d value", 2*i);
            gsl_test(gsl_vector_get(coeffs, 2*i+1) != e1[2*i+1],
                    "integration_coefficients 1 linear %d stride", 2*i + 1);
        }

        /* Zeroth derivative, integrating from high to low */
        gsl_vector_set_all(coeffs, -555.0);
        gsl_test(suzerain_bspline_integration_coefficients(
                0 , gsl_vector_ptr(coeffs,0), inc,
                GSL_DBL_MAX, -GSL_DBL_MAX, scratch, w, dw),
                "integration_coefficients 0 linear");
        for (size_t i = 0; i < w->n; ++i) {
            gsl_test_rel(gsl_vector_get(coeffs, 2*i), -1*e0[2*i],
                    GSL_DBL_EPSILON*1000,
                    "integration_coefficients 0 linear %d value", 2*i);
            gsl_test(gsl_vector_get(coeffs, 2*i+1) != e0[2*i+1],
                    "integration_coefficients 0 linear %d stride", 2*i + 1);
        }

        /* First derivative, integrating from high to low */
        gsl_vector_set_all(coeffs, -555.0);
        gsl_test(suzerain_bspline_integration_coefficients(
                1 , gsl_vector_ptr(coeffs,0), inc,
                GSL_DBL_MAX, -GSL_DBL_MAX, scratch, w, dw),
                "integration_coefficients 1 linear");
        for (size_t i = 0; i < w->n; ++i) {
            gsl_test_rel(gsl_vector_get(coeffs, 2*i), -1*e1[2*i],
                    GSL_DBL_EPSILON*1000,
                    "integration_coefficients 1 linear %d value", 2*i);
            gsl_test(gsl_vector_get(coeffs, 2*i+1) != e1[2*i+1],
                    "integration_coefficients 1 linear %d stride", 2*i + 1);
        }

        gsl_vector_free(coeffs);
        free_workspaces(&w, &dw, &scratch);
    }

    /* Piecewise quadratics */
    {
        alloc_workspaces(3, sizeof(b)/sizeof(b[0]), b, &w, &dw, &scratch);
        const size_t inc = 1;
        const double e0[] = { 1./3, 2./3, 1, 2./3, 1./3 };
        const double e1[] = {   -1,    0, 0,    0,    1 };
        const double e2[] = {    2,   -2, 0,   -2,    2 };
        gsl_vector *coeffs = gsl_vector_alloc(sizeof(e0)/sizeof(e0[0]));

        gsl_vector_set_all(coeffs, -555.0);
        gsl_test(suzerain_bspline_integration_coefficients(
                0 , gsl_vector_ptr(coeffs,0), inc,
                -GSL_DBL_MAX, GSL_DBL_MAX, scratch, w, dw),
                "integration_coefficients 0 quadratic");
        for (size_t i = 0; i < w->n; ++i) {
            gsl_test_rel(gsl_vector_get(coeffs, i), e0[i],
                    GSL_DBL_EPSILON*1000,
                    "integration_coefficients 0 quadratic %d value", i);
        }

        gsl_vector_set_all(coeffs, -555.0);
        gsl_test(suzerain_bspline_integration_coefficients(
                1 , gsl_vector_ptr(coeffs,0), inc,
                -GSL_DBL_MAX, GSL_DBL_MAX, scratch, w, dw),
                "integration_coefficients 1 quadratic");
        for (size_t i = 0; i < w->n; ++i) {
            gsl_test_rel(gsl_vector_get(coeffs, i), e1[i],
                    GSL_DBL_EPSILON*1000,
                    "integration_coefficients 1 quadratic %d value", i);
        }

        gsl_vector_set_all(coeffs, -555.0);
        gsl_test(suzerain_bspline_integration_coefficients(
                2 , gsl_vector_ptr(coeffs,0), inc,
                -GSL_DBL_MAX, GSL_DBL_MAX, scratch, w, dw),
                "integration_coefficients 2 quadratic");
        for (size_t i = 0; i < w->n; ++i) {
            gsl_test_rel(gsl_vector_get(coeffs, i), e2[i],
                    GSL_DBL_EPSILON*1000,
                    "integration_coefficients 2 quadratic %d value", i);
        }

        gsl_vector_free(coeffs);
        free_workspaces(&w, &dw, &scratch);
    }

    /* Piecewise cubics */
    {
        alloc_workspaces(4, sizeof(b)/sizeof(b[0]), b, &w, &dw, &scratch);
        const size_t inc = 1;
        const double e0[] = { 1./4, 1./2, 3./4, 3./4, 1./2, 1./4 };
        const double e1[] = {   -1,    0,    0,    0,    0,    1 };
        const double e2[] = {    3,   -3,    0,    0,   -3,    3 };
        const double e3[] = {   -6,    9,   -3,    3,   -9,    6 };

        gsl_vector *coeffs = gsl_vector_alloc(sizeof(e0)/sizeof(e0[0]));

        gsl_vector_set_all(coeffs, -555.0);
        gsl_test(suzerain_bspline_integration_coefficients(
                0 , gsl_vector_ptr(coeffs,0), inc,
                -GSL_DBL_MAX, GSL_DBL_MAX, scratch, w, dw),
                "integration_coefficients 0 cubic");
        for (size_t i = 0; i < w->n; ++i) {
            gsl_test_rel(gsl_vector_get(coeffs, i), e0[i],
                    GSL_DBL_EPSILON*1000,
                    "integration_coefficients 0 cubic %d value", i);
        }

        gsl_vector_set_all(coeffs, -555.0);
        gsl_test(suzerain_bspline_integration_coefficients(
                1 , gsl_vector_ptr(coeffs,0), inc,
                -GSL_DBL_MAX, GSL_DBL_MAX, scratch, w, dw),
                "integration_coefficients 1 cubic");
        for (size_t i = 0; i < w->n; ++i) {
            gsl_test_rel(gsl_vector_get(coeffs, i), e1[i],
                    GSL_DBL_EPSILON*1000,
                    "integration_coefficients 1 cubic %d value", i);
        }

        gsl_vector_set_all(coeffs, -555.0);
        gsl_test(suzerain_bspline_integration_coefficients(
                2 , gsl_vector_ptr(coeffs,0), inc,
                -GSL_DBL_MAX, GSL_DBL_MAX, scratch, w, dw),
                "integration_coefficients 2 cubic");
        for (size_t i = 0; i < w->n; ++i) {
            gsl_test_rel(gsl_vector_get(coeffs, i), e2[i],
                    GSL_DBL_EPSILON*1000,
                    "integration_coefficients 2 cubic %d value", i);
        }

        gsl_vector_set_all(coeffs, -555.0);
        gsl_test(suzerain_bspline_integration_coefficients(
                3 , gsl_vector_ptr(coeffs,0), inc,
                -GSL_DBL_MAX, GSL_DBL_MAX, scratch, w, dw),
                "integration_coefficients 3 cubic");
        for (size_t i = 0; i < w->n; ++i) {
            gsl_test_rel(gsl_vector_get(coeffs, i), e3[i],
                    GSL_DBL_EPSILON*1000,
                    "integration_coefficients 3 cubic %d value", i);
        }

        gsl_vector_free(coeffs);
        free_workspaces(&w, &dw, &scratch);
    }
}

static void test_linear_combination()
{
    const double b[] = { 0.0, 1.0, 2.0, 3.0 };
    gsl_bspline_workspace *w;
    gsl_bspline_deriv_workspace *dw;
    gsl_matrix *scratch;
    alloc_workspaces(4, sizeof(b)/sizeof(b[0]), b, &w, &dw, &scratch);
    gsl_vector *coeffs = gsl_vector_alloc(w->n);

    /* Check that we have a partition of unity at the breakpoints */
    /* All derivatives should be zero as well */
    {
        gsl_vector_set_all(coeffs, 1.0);

        const size_t ldvalues = sizeof(b)/sizeof(b[0]);
        gsl_vector *values = gsl_vector_alloc(w->k * ldvalues);
        gsl_vector_set_all(values, -555.0);

        /* Evaluate all derivatives */
        gsl_test(suzerain_bspline_linear_combination(
                    w->k - 1, gsl_vector_const_ptr(coeffs, 0),
                    sizeof(b)/sizeof(b[0]), b, gsl_vector_ptr(values, 0),
                    ldvalues, scratch, w, dw),
                "linear_combination unity multiple");
        for (size_t i = 0; i < ldvalues; ++i) {
            gsl_test_rel(gsl_vector_get(values, i), 1.0, GSL_DBL_EPSILON*1000,
                    "linear_combination unity value %i", i);
        }
        for (size_t i = ldvalues; i < values->size; ++i) {
            gsl_test_abs(gsl_vector_get(values, i), 0.0, GSL_DBL_EPSILON,
                    "linear_combination unity derivs value %i", i);

        }

        /* Evaluate only the (k-1)th derivative */
        gsl_vector_set_all(values, -555.0);
        gsl_test(suzerain_bspline_linear_combination(
                    w->k - 1, gsl_vector_const_ptr(coeffs, 0),
                    sizeof(b)/sizeof(b[0]), b, gsl_vector_ptr(values, 0),
                    0, scratch, w, dw),
                "linear_combination unity single");
        for (size_t i = 0; i < ldvalues; ++i) {
            gsl_test_abs(gsl_vector_get(values, i), 0.0, GSL_DBL_EPSILON,
                    "linear_combination unity derivs value %i", i);
        }
        for (size_t i = ldvalues; i < values->size; ++i) {
            gsl_test_rel(gsl_vector_get(values, i), -555.0, GSL_DBL_EPSILON,
                    "linear_combination unity untouched %i", i);

        }

        gsl_vector_free(values);
    }

    /* Spot checks on basis function evaluation and derivatives */
    {
        const double points[] = { 1.0 };
        const size_t ldvalues = sizeof(points)/sizeof(points[0]);
        gsl_vector *values = gsl_vector_alloc(w->k * ldvalues);

        /* Investigate second basis function */
        gsl_vector_set_basis(coeffs, 1);
        gsl_vector_set_all(values, -555.0);
        gsl_test(suzerain_bspline_linear_combination(
                    w->k - 1, gsl_vector_const_ptr(coeffs, 0),
                    ldvalues,  points, gsl_vector_ptr(values, 0),
                    ldvalues, scratch, w, dw),
                "linear_combination basis 2");
        gsl_test_rel(gsl_vector_get(values, 0), 1./4.,  GSL_DBL_EPSILON*1000,
                "linear_combination basis 2 deriv %d", 0);
        gsl_test_rel(gsl_vector_get(values, 1), -3./4., GSL_DBL_EPSILON*1000,
                "linear_combination basis 2 deriv %d", 1);
        gsl_test_rel(gsl_vector_get(values, 2), 3./2., GSL_DBL_EPSILON*1000,
                "linear_combination basis 2 deriv %d", 2);

        /* Investigate third basis function */
        gsl_vector_set_basis(coeffs, 2);
        gsl_vector_set_all(values, -555.0);
        gsl_test(suzerain_bspline_linear_combination(
                w->k - 1, gsl_vector_const_ptr(coeffs, 0), ldvalues,  points,
                gsl_vector_ptr(values, 0), ldvalues, scratch, w, dw),
                "linear_combination basis 3");
        gsl_test_rel(gsl_vector_get(values, 0), 7./12., GSL_DBL_EPSILON*1000,
                "linear_combination basis 3 deriv %d", 0);
        gsl_test_rel(gsl_vector_get(values, 1), 1./4.,  GSL_DBL_EPSILON*1000,
                "linear_combination basis 3 deriv %d", 1);
        gsl_test_rel(gsl_vector_get(values, 2), -5./2., GSL_DBL_EPSILON*1000,
                "linear_combination basis 3 deriv %d", 2);

        /* Investigate fourth basis function */
        gsl_vector_set_basis(coeffs, 3);
        gsl_vector_set_all(values, -555.0);
        gsl_test(suzerain_bspline_linear_combination(
                w->k - 1, gsl_vector_const_ptr(coeffs, 0), ldvalues,  points,
                gsl_vector_ptr(values, 0), ldvalues, scratch, w, dw),
                "linear_combination basis 4");
        gsl_test_rel(gsl_vector_get(values, 0), 1./6., GSL_DBL_EPSILON*1000,
                "linear_combination basis 4 deriv %d", 0);
        gsl_test_rel(gsl_vector_get(values, 1), 1./2., GSL_DBL_EPSILON*1000,
                "linear_combination basis 4 deriv %d", 1);
        gsl_test_rel(gsl_vector_get(values, 2), 1., GSL_DBL_EPSILON*1000,
                "linear_combination basis 4 deriv %d", 2);

        /* Investigate fifth basis function */
        /* All zero due to influence of endpoints/repeated knots */
        gsl_vector_set_basis(coeffs, 4);
        gsl_vector_set_all(values, -555.0);
        gsl_test(suzerain_bspline_linear_combination(
                w->k - 1, gsl_vector_const_ptr(coeffs, 0), ldvalues,  points,
                gsl_vector_ptr(values, 0), ldvalues, scratch, w, dw),
                "linear_combination basis 5");
        gsl_test_abs(gsl_vector_get(values, 0), 0., GSL_DBL_EPSILON*1000,
                "linear_combination basis 5 deriv %d", 0);
        gsl_test_abs(gsl_vector_get(values, 1), 0., GSL_DBL_EPSILON*1000,
                "linear_combination basis 5 deriv %d", 1);
        gsl_test_abs(gsl_vector_get(values, 2), 0., GSL_DBL_EPSILON*1000,
                "linear_combination basis 5 deriv %d", 2);

        /* Check behavior of a linear combination of basis functions */
        /* Include some poison values to check bounds access correct */
        gsl_vector_set(coeffs, 0, GSL_NAN);
        gsl_vector_set(coeffs, 1, 1.0);
        gsl_vector_set(coeffs, 2, 2.0);
        gsl_vector_set(coeffs, 3, 3.0);
        gsl_vector_set(coeffs, 4, 4.0);
        gsl_vector_set(coeffs, 5, GSL_NAN);
        gsl_test(suzerain_bspline_linear_combination(
                w->k - 1, gsl_vector_const_ptr(coeffs, 0), ldvalues,  points,
                gsl_vector_ptr(values, 0), ldvalues, scratch, w, dw),
                "linear_combination linear combination multiple");
        gsl_test_abs(gsl_vector_get(values, 0), 23./12., GSL_DBL_EPSILON*1000,
                "linear_combination linear combination deriv %d", 0);
        gsl_test_abs(gsl_vector_get(values, 1), 5./4., GSL_DBL_EPSILON*1000,
                "linear_combination linear combination deriv %d", 1);
        gsl_test_abs(gsl_vector_get(values, 2), -1./2., GSL_DBL_EPSILON*1000,
                "linear_combination linear combination deriv %d", 2);

        /* Linear combination of basis functions for second deriv */
        gsl_test(suzerain_bspline_linear_combination(
                2, gsl_vector_const_ptr(coeffs, 0), ldvalues,  points,
                gsl_vector_ptr(values, 0), 0, scratch, w, dw),
                "linear_combination linear combination single");
        gsl_test_abs(gsl_vector_get(values, 0), -1./2., GSL_DBL_EPSILON*1000,
                "linear_combination linear combination single value");

        gsl_vector_free(values);
    }

    gsl_vector_free(coeffs);
    free_workspaces(&w, &dw, &scratch);
}

static void test_linear_combination_complex()
{
    const double b[] = { 0.0, 1.0, 2.0, 3.0 };
    gsl_bspline_workspace *w;
    gsl_bspline_deriv_workspace *dw;
    gsl_matrix *scratch;
    alloc_workspaces(4, sizeof(b)/sizeof(b[0]), b, &w, &dw, &scratch);
    gsl_vector_complex *coeffs = gsl_vector_complex_alloc(w->n);

    /* Check that we have a partition of unity at the breakpoints */
    /* separately on the real- and imaginary- parts of the coefficients */
    {
        const size_t ldvalues = sizeof(b)/sizeof(b[0]);
        gsl_vector_complex *values = gsl_vector_complex_alloc(w->k * ldvalues);
        gsl_vector_complex_set_all(values, gsl_complex_rect(-555.0, -555.0));

        /* Set real-only components to all ones */
        gsl_vector_complex_set_all(coeffs, gsl_complex_rect(1,0));

        /* Evaluate all derivatives */
        gsl_test(suzerain_bspline_linear_combination_complex(
                w->k - 1,
                (complex_double *) gsl_vector_complex_const_ptr(coeffs, 0),
                sizeof(b)/sizeof(b[0]), b,
                (complex_double *) gsl_vector_complex_ptr(values, 0),
                ldvalues, scratch, w, dw),
                "linear_combination_complex real unity multiple");
        for (size_t i = 0; i < ldvalues; ++i) {
            gsl_test_rel(GSL_REAL(gsl_vector_complex_get(values, i)), 1.0,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex real unity real value %i", i);
            gsl_test_abs(GSL_IMAG(gsl_vector_complex_get(values, i)), 0.0,
                GSL_DBL_EPSILON,
                "linear_combination_complex real unity imag value %i", i);
        }
        for (size_t i = ldvalues; i < values->size; ++i) {
            gsl_test_abs(GSL_REAL(gsl_vector_complex_get(values, i)), 0.0,
                GSL_DBL_EPSILON,
                "linear_combination_complex real unity real value %i", i);
            gsl_test_abs(GSL_IMAG(gsl_vector_complex_get(values, i)), 0.0,
                GSL_DBL_EPSILON,
                "linear_combination_complex real unity imag value %i", i);
        }

        /* Evaluate only the (k-1)th derivative */
        gsl_test(suzerain_bspline_linear_combination_complex(
                w->k - 1,
                (complex_double *) gsl_vector_complex_const_ptr(coeffs, 0),
                sizeof(b)/sizeof(b[0]), b,
                (complex_double *) gsl_vector_complex_ptr(values, 0),
                0, scratch, w, dw),
                "linear_combination_complex real unity single");
        for (size_t i = 0; i < ldvalues; ++i) {
            gsl_test_abs(GSL_REAL(gsl_vector_complex_get(values, i)), 0.0,
                GSL_DBL_EPSILON,
                "linear_combination_complex real unity real value %i", i);
            gsl_test_abs(GSL_IMAG(gsl_vector_complex_get(values, i)), 0.0,
                GSL_DBL_EPSILON,
                "linear_combination_complex real unity imag value %i", i);
        }

        /* Set imag-only components to all ones */
        gsl_vector_complex_set_all(coeffs, gsl_complex_rect(0,1));

        /* Evaluate all derivatives */
        gsl_test(suzerain_bspline_linear_combination_complex(
                w->k - 1,
                (complex_double *) gsl_vector_complex_const_ptr(coeffs, 0),
                sizeof(b)/sizeof(b[0]), b,
                (complex_double *) gsl_vector_complex_ptr(values, 0),
                ldvalues, scratch, w, dw),
                "linear_combination_complex imag unity multiple");
        for (size_t i = 0; i < ldvalues; ++i) {
            gsl_test_abs(GSL_REAL(gsl_vector_complex_get(values, i)), 0.0,
                GSL_DBL_EPSILON,
                "linear_combination_complex imag unity real value %i", i);
            gsl_test_rel(GSL_IMAG(gsl_vector_complex_get(values, i)), 1.0,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex imag unity imag value %i", i);
        }
        for (size_t i = ldvalues; i < values->size; ++i) {
            gsl_test_abs(GSL_REAL(gsl_vector_complex_get(values, i)), 0.0,
                GSL_DBL_EPSILON,
                "linear_combination_complex imag unity real value %i", i);
            gsl_test_abs(GSL_IMAG(gsl_vector_complex_get(values, i)), 0.0,
                GSL_DBL_EPSILON,
                "linear_combination_complex imag unity imag value %i", i);
        }

        /* Evaluate only the (k-1)th derivative */
        gsl_test(suzerain_bspline_linear_combination_complex(
                w->k - 1,
                (complex_double *) gsl_vector_complex_const_ptr(coeffs, 0),
                sizeof(b)/sizeof(b[0]), b,
                (complex_double *) gsl_vector_complex_ptr(values, 0),
                0, scratch, w, dw),
                "linear_combination_complex imag unity single");
        for (size_t i = 0; i < ldvalues; ++i) {
            gsl_test_abs(GSL_REAL(gsl_vector_complex_get(values, i)), 0.0,
                GSL_DBL_EPSILON,
                "linear_combination_complex imag unity real value %i", i);
            gsl_test_abs(GSL_IMAG(gsl_vector_complex_get(values, i)), 0.0,
                GSL_DBL_EPSILON,
                "linear_combination_complex imag unity imag value %i", i);
        }

        gsl_vector_complex_free(values);
    }

    /* Spot checks on basis function evaluation and derivatives */
    {
        const double points[] = { 1.0 };
        const size_t ldvalues = sizeof(points)/sizeof(points[0]);
        gsl_vector_complex *values = gsl_vector_complex_alloc(w->k * ldvalues);

        /* Investigate second basis function */
        gsl_vector_complex_set_zero(coeffs);
        gsl_vector_complex_set(coeffs, 1, gsl_complex_rect(1.0, -0.5));
        gsl_vector_complex_set_all(values, gsl_complex_rect(-555.0,-555.0));
        gsl_test(suzerain_bspline_linear_combination_complex(
                    w->k - 1,
                    (complex_double *) gsl_vector_complex_const_ptr(coeffs, 0),
                    ldvalues, points,
                    (complex_double *) gsl_vector_complex_ptr(values, 0),
                    ldvalues, scratch, w, dw),
                "linear_combination_complex basis 2");
        gsl_test_rel(GSL_REAL(gsl_vector_complex_get(values, 0)), 1./4.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex basis 2 deriv %d real", 0);
        gsl_test_rel(GSL_IMAG(gsl_vector_complex_get(values, 0)), -1./8.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex basis 2 deriv %d imag", 0);
        gsl_test_rel(GSL_REAL(gsl_vector_complex_get(values, 1)), -3./4.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex basis 2 deriv %d real", 1);
        gsl_test_rel(GSL_IMAG(gsl_vector_complex_get(values, 1)), 3./8.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex basis 2 deriv %d imag", 1);
        gsl_test_rel(GSL_REAL(gsl_vector_complex_get(values, 2)), 3./2.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex basis 2 deriv %d real", 2);
        gsl_test_rel(GSL_IMAG(gsl_vector_complex_get(values, 2)), -3./4.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex basis 2 deriv %d imag", 2);

        /* Investigate third basis function */
        gsl_vector_complex_set_zero(coeffs);
        gsl_vector_complex_set(coeffs, 2, gsl_complex_rect(1.0, -0.5));
        gsl_vector_complex_set_all(values, gsl_complex_rect(-555.0,-555.0));
        gsl_test(suzerain_bspline_linear_combination_complex(
                    w->k - 1,
                    (complex_double *) gsl_vector_complex_const_ptr(coeffs, 0),
                    ldvalues, points,
                    (complex_double *) gsl_vector_complex_ptr(values, 0),
                    ldvalues, scratch, w, dw),
                "linear_combination_complex basis 3");
        gsl_test_rel(GSL_REAL(gsl_vector_complex_get(values, 0)), 7./12.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex basis 3 deriv %d real", 0);
        gsl_test_rel(GSL_IMAG(gsl_vector_complex_get(values, 0)), -7./24.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex basis 3 deriv %d imag", 0);
        gsl_test_rel(GSL_REAL(gsl_vector_complex_get(values, 1)), 1./4.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex basis 3 deriv %d real", 1);
        gsl_test_rel(GSL_IMAG(gsl_vector_complex_get(values, 1)), -1./8.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex basis 3 deriv %d imag", 1);
        gsl_test_rel(GSL_REAL(gsl_vector_complex_get(values, 2)), -5./2.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex basis 3 deriv %d real", 2);
        gsl_test_rel(GSL_IMAG(gsl_vector_complex_get(values, 2)),  5./4.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex basis 3 deriv %d imag", 2);

        /* Investigate fourth basis function */
        gsl_vector_complex_set_zero(coeffs);
        gsl_vector_complex_set(coeffs, 3, gsl_complex_rect(1.0, -0.5));
        gsl_vector_complex_set_all(values, gsl_complex_rect(-555.0,-555.0));
        gsl_test(suzerain_bspline_linear_combination_complex(
                    w->k - 1,
                    (complex_double *) gsl_vector_complex_const_ptr(coeffs, 0),
                    ldvalues, points,
                    (complex_double *) gsl_vector_complex_ptr(values, 0),
                    ldvalues, scratch, w, dw),
                "linear_combination_complex basis 4");
        gsl_test_rel(GSL_REAL(gsl_vector_complex_get(values, 0)), 1./6.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex basis 4 deriv %d real", 0);
        gsl_test_rel(GSL_IMAG(gsl_vector_complex_get(values, 0)), -1./12.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex basis 4 deriv %d imag", 0);
        gsl_test_rel(GSL_REAL(gsl_vector_complex_get(values, 1)), 1./2.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex basis 4 deriv %d real", 1);
        gsl_test_rel(GSL_IMAG(gsl_vector_complex_get(values, 1)), -1./4.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex basis 4 deriv %d imag", 1);
        gsl_test_rel(GSL_REAL(gsl_vector_complex_get(values, 2)), 1.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex basis 4 deriv %d real", 2);
        gsl_test_rel(GSL_IMAG(gsl_vector_complex_get(values, 2)),  -0.5,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex basis 4 deriv %d imag", 2);

        /* Investigate fifth basis function */
        /* All zero due to influence of endpoints/repeated knots */
        gsl_vector_complex_set_zero(coeffs);
        gsl_vector_complex_set(coeffs, 4, gsl_complex_rect(1.0, -0.5));
        gsl_vector_complex_set_all(values, gsl_complex_rect(-555.0,-555.0));
        gsl_test(suzerain_bspline_linear_combination_complex(
                    w->k - 1,
                    (complex_double *) gsl_vector_complex_const_ptr(coeffs, 0),
                    ldvalues, points,
                    (complex_double *) gsl_vector_complex_ptr(values, 0),
                    ldvalues, scratch, w, dw),
                "linear_combination_complex basis 5");
        gsl_test_abs(GSL_REAL(gsl_vector_complex_get(values, 0)), 0.,
                GSL_DBL_EPSILON,
                "linear_combination_complex basis 5 deriv %d real", 0);
        gsl_test_abs(GSL_IMAG(gsl_vector_complex_get(values, 0)), 0.,
                GSL_DBL_EPSILON,
                "linear_combination_complex basis 5 deriv %d imag", 0);
        gsl_test_abs(GSL_REAL(gsl_vector_complex_get(values, 1)), 0.,
                GSL_DBL_EPSILON,
                "linear_combination_complex basis 5 deriv %d real", 1);
        gsl_test_abs(GSL_IMAG(gsl_vector_complex_get(values, 1)), 0.,
                GSL_DBL_EPSILON,
                "linear_combination_complex basis 5 deriv %d imag", 1);
        gsl_test_abs(GSL_REAL(gsl_vector_complex_get(values, 2)), 0.,
                GSL_DBL_EPSILON,
                "linear_combination_complex basis 5 deriv %d real", 2);
        gsl_test_abs(GSL_IMAG(gsl_vector_complex_get(values, 2)), 0.,
                GSL_DBL_EPSILON,
                "linear_combination_complex basis 5 deriv %d imag", 2);

        /* Check behavior of a linear combination of basis functions */
        /* Include some poison values to check bounds access correct */
        gsl_vector_complex_set(coeffs, 0, gsl_complex_rect(GSL_NAN, GSL_NAN));
        gsl_vector_complex_set(coeffs, 1, gsl_complex_rect(1.0, -0.5));
        gsl_vector_complex_set(coeffs, 2, gsl_complex_rect(2.0, -1.0));
        gsl_vector_complex_set(coeffs, 3, gsl_complex_rect(3.0, -1.5));
        gsl_vector_complex_set(coeffs, 4, gsl_complex_rect(4.0, -2.0));
        gsl_vector_complex_set(coeffs, 5, gsl_complex_rect(GSL_NAN, GSL_NAN));
        gsl_test(suzerain_bspline_linear_combination_complex(
                    w->k - 1,
                    (complex_double *) gsl_vector_complex_const_ptr(coeffs, 0),
                    ldvalues, points,
                    (complex_double *) gsl_vector_complex_ptr(values, 0),
                    ldvalues, scratch, w, dw),
                "linear_combination_complex linear combination");
        gsl_test_rel(GSL_REAL(gsl_vector_complex_get(values, 0)), 23./12.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex linear combination deriv %d real", 0);
        gsl_test_rel(GSL_IMAG(gsl_vector_complex_get(values, 0)), -23./24.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex linear combination deriv %d imag", 0);
        gsl_test_rel(GSL_REAL(gsl_vector_complex_get(values, 1)), 5./4.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex linear combination deriv %d real", 1);
        gsl_test_rel(GSL_IMAG(gsl_vector_complex_get(values, 1)), -5./8.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex linear combination deriv %d imag", 1);
        gsl_test_rel(GSL_REAL(gsl_vector_complex_get(values, 2)),  -1./2.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex linear combination deriv %d real", 2);
        gsl_test_rel(GSL_IMAG(gsl_vector_complex_get(values, 2)),  1./4.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex linear combination deriv %d imag", 2);

        /* Linear combination of basis functions for second deriv */
        gsl_test(suzerain_bspline_linear_combination_complex(
                    2,
                    (complex_double *) gsl_vector_complex_const_ptr(coeffs, 0),
                    ldvalues, points,
                    (complex_double *) gsl_vector_complex_ptr(values, 0),
                    0, scratch, w, dw),
                "linear_combination_complex linear combination single");
        gsl_test_rel(GSL_REAL(gsl_vector_complex_get(values, 0)),  -1./2.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex linear combination deriv %d real", 2);
        gsl_test_rel(GSL_IMAG(gsl_vector_complex_get(values, 0)),  1./4.,
                GSL_DBL_EPSILON*1000,
                "linear_combination_complex linear combination deriv %d imag", 2);

        gsl_vector_complex_free(values);
    }

    gsl_vector_complex_free(coeffs);
    free_workspaces(&w, &dw, &scratch);
}

static void test_crossing()
{
    /* Test case uses B-spline expansion of x**2 on [-2.0, 2.5] */
    /* Coefficients computed analytically via Mathematica */
    const size_t k        = 5;
    const double b[]      = { -2.0, -0.5, 0.0, 1.5, 2.5 };
    const double coeffs[] = { 4., 5./2, 1, -11./24, 7./24, 55./24, 5., 25./4 };

    gsl_bspline_workspace *w;
    gsl_bspline_deriv_workspace *dw;
    gsl_matrix *scratch;
    alloc_workspaces(k, sizeof(b)/sizeof(b[0]), b, &w, &dw, &scratch);

    /* Check, via linear_combination, that we can recover x**2 at abscissae. */
    /* This is a zeroth-order test that the coefficients are sane. */
    /* If they're not, nothing else in within this function will work. */
    {
        double g[] = {-2., -13./8, -9./8, -1./4, 7./8, 13./8, 9./4, 5./2};
        double v[] = { 0.,     0.,    0.,    0.,   0.,    0.,   0.,   0.};
        suzerain_bspline_linear_combination(0, coeffs, sizeof(g)/sizeof(g[0]),
                                            g, v, 0, scratch, w, dw);
        for (size_t i = 0; i < sizeof(g)/sizeof(g[0]); ++i) {
            gsl_test_rel(v[i], g[i]*g[i], 100*GSL_DBL_EPSILON,
                         "%s: reproduce x**2 at collocation point %d",
                          __func__, i);
        }
    }

    /* Successful test for 0th derivative checking range */
    {
        double lower = -0.5, upper = 2.5, location = GSL_NAN;
        const int status = suzerain_bspline_crossing(
                0, coeffs, 1.0, &lower, &upper, 100,
                GSL_DBL_EPSILON, GSL_DBL_EPSILON, &location, scratch, w, dw);
        gsl_test(status, "GSL_SUCCESS at %s:%d", __FILE__, __LINE__);
        gsl_test_rel(location, 1.0, GSL_SQRT_DBL_EPSILON,
                     "Test for 0th derivative 1-crossing of x**2 at x = 1");
    }

    /* Successful test for 0th derivative checking range */
    {
        double lower = -2.0, upper = 0.5, location = GSL_NAN;
        const int status = suzerain_bspline_crossing(
                0, coeffs, 1.0, &lower, &upper, 100,
                GSL_DBL_EPSILON, GSL_DBL_EPSILON, &location, scratch, w, dw);
        gsl_test(status, "GSL_SUCCESS at %s:%d", __FILE__, __LINE__);
        gsl_test_rel(location, -1.0, GSL_SQRT_DBL_EPSILON,
                     "Test for 0th derivative 1-crossing of x**2 at x = -1");
    }

    /* Successful test for 1st derivative across whole range */
    {
        double lower = -2.0, upper = 2.5, location = GSL_NAN;
        const int status = suzerain_bspline_crossing(
                1, coeffs, 1.0, &lower, &upper, 100,
                GSL_DBL_EPSILON, GSL_DBL_EPSILON, &location, scratch, w, dw);
        gsl_test(status, "GSL_SUCCESS at %s:%d", __FILE__, __LINE__);
        gsl_test_rel(location, 0.5, GSL_SQRT_DBL_EPSILON,
                     "Test for 1st derivative 1-crossing of x**2 at x = 1/2");
    }

    /* Ensure looking for a nonexistent root will not bring down a binary */
    {
        double lower = -0.5, upper = 2.5, location = 555;
        const int status = suzerain_bspline_crossing(
                0, coeffs, 9.0, &lower, &upper, 100,
                GSL_DBL_EPSILON, GSL_DBL_EPSILON, &location, scratch, w, dw);
        gsl_test(!status, "Failure expected at %s:%d", __FILE__, __LINE__);
        gsl_test(!gsl_isnan(location), "Failure produces NaN at %s:%d",
                 __FILE__, __LINE__);
    }

    free_workspaces(&w, &dw, &scratch);
}

static void test_spacing_greville_abscissae()
{
    const char msg[] = "spacing_greville_abscissae for k=%d, i=%d";

    const double b[] = { -3, 2, 9, 12, 21 };
    gsl_bspline_workspace *w;
    gsl_bspline_deriv_workspace *dw;
    gsl_matrix *scratch;

    {
        const int k = 2;
        alloc_workspaces(k, sizeof(b)/sizeof(b[0]), b, &w, &dw, &scratch);
        gsl_test_rel(suzerain_bspline_spacing_greville_abscissae(0, w),
                     5, GSL_DBL_EPSILON, msg, k, 0);
        gsl_test_rel(suzerain_bspline_spacing_greville_abscissae(1, w),
                     5, GSL_DBL_EPSILON, msg, k, 1);
        gsl_test_rel(suzerain_bspline_spacing_greville_abscissae(2, w),
                     3, GSL_DBL_EPSILON, msg, k, 2);
        gsl_test_rel(suzerain_bspline_spacing_greville_abscissae(3, w),
                     3, GSL_DBL_EPSILON, msg, k, 3);
        gsl_test_rel(suzerain_bspline_spacing_greville_abscissae(4, w),
                     9, GSL_DBL_EPSILON, msg, k, 4);
        free_workspaces(&w, &dw, &scratch);
    }

    {
        const int k = 3;
        alloc_workspaces(k, sizeof(b)/sizeof(b[0]), b, &w, &dw, &scratch);
        gsl_test_rel(suzerain_bspline_spacing_greville_abscissae(0, w),
                     2.5, GSL_DBL_EPSILON, msg, k, 0);
        gsl_test_rel(suzerain_bspline_spacing_greville_abscissae(1, w),
                     2.5, GSL_DBL_EPSILON, msg, k, 1);
        gsl_test_rel(suzerain_bspline_spacing_greville_abscissae(2, w),
                     5.0, GSL_DBL_EPSILON, msg, k, 2);
        gsl_test_rel(suzerain_bspline_spacing_greville_abscissae(3, w),
                     5.0, GSL_DBL_EPSILON, msg, k, 3);
        gsl_test_rel(suzerain_bspline_spacing_greville_abscissae(4, w),
                     4.5, GSL_DBL_EPSILON, msg, k, 4);
        gsl_test_rel(suzerain_bspline_spacing_greville_abscissae(5, w),
                     4.5, GSL_DBL_EPSILON, msg, k, 5);
        free_workspaces(&w, &dw, &scratch);
    }
}

static void test_spacing_breakpoints()
{
    const char msg[] = "spacing_breakpoints for k=%d, i=%d";

    const double b[] = { -3, 2, 9, 12, 21 };
    gsl_bspline_workspace *w;
    gsl_bspline_deriv_workspace *dw;
    gsl_matrix *scratch;

    {
        const int k = 2;
        alloc_workspaces(k, sizeof(b)/sizeof(b[0]), b, &w, &dw, &scratch);
        gsl_test_rel(suzerain_bspline_spacing_breakpoints(0, w),
                     5, GSL_DBL_EPSILON, msg, k, 0);
        gsl_test_rel(suzerain_bspline_spacing_breakpoints(1, w),
                     5, GSL_DBL_EPSILON, msg, k, 1);
        gsl_test_rel(suzerain_bspline_spacing_breakpoints(2, w),
                     3, GSL_DBL_EPSILON, msg, k, 2);
        gsl_test_rel(suzerain_bspline_spacing_breakpoints(3, w),
                     3, GSL_DBL_EPSILON, msg, k, 3);
        gsl_test_rel(suzerain_bspline_spacing_breakpoints(4, w),
                     9, GSL_DBL_EPSILON, msg, k, 4);
        free_workspaces(&w, &dw, &scratch);
    }

    {
        const int k = 3;
        alloc_workspaces(k, sizeof(b)/sizeof(b[0]), b, &w, &dw, &scratch);
        gsl_test_rel(suzerain_bspline_spacing_breakpoints(0, w),
                     5.0, GSL_DBL_EPSILON, msg, k, 0);
        gsl_test_rel(suzerain_bspline_spacing_breakpoints(1, w),
                     5.0, GSL_DBL_EPSILON, msg, k, 1);
        gsl_test_rel(suzerain_bspline_spacing_breakpoints(2, w),
                     3.0, GSL_DBL_EPSILON, msg, k, 2);
        gsl_test_rel(suzerain_bspline_spacing_breakpoints(3, w),
                     3.0, GSL_DBL_EPSILON, msg, k, 3);
        gsl_test_rel(suzerain_bspline_spacing_breakpoints(4, w),
                     3.0, GSL_DBL_EPSILON, msg, k, 4);
        gsl_test_rel(suzerain_bspline_spacing_breakpoints(5, w),
                     9.0, GSL_DBL_EPSILON, msg, k, 5);
        free_workspaces(&w, &dw, &scratch);
    }
}

// Test that curve fits reasonably reproduce the source data.
// Mainly meant as a sanity check on the coded coefficients.
// Test cases chosen pseudo-randomly from the original fit data.
// Choice required that we could nail at least one output per k.
static void test_bspline_htstretch2_evdeltascale_greville_abscissae()
{
    static const double data[][5] = {
        /* 0*/ {4   , 1.25  , 37    , 6.73374  , 3.07899},
        /* 1*/ {4   , 1.5   , 379   , 5.76858  , 2.81976},
        /* 2*/ {4   , 2.25  , 307   , 6.03399  , 2.89108},
        /* 3*/ {4   , 2.75  , 419   , 6.02432  , 2.88836},
        /* 4*/ {4   , 3     , 509   , 5.99355  , 2.88014},
        /* 5*/ {5   , 0     , 307   , 6.52402  , 4.05568},
        /* 6*/ {5   , 1.5   , 1024  , 6.73901  , 4.10416},
        /* 7*/ {5   , 1.75  , 139   , 7.55957  , 4.30002},
        /* 8*/ {5   , 2.5   , 859   , 6.94602  , 4.15221},
        /* 9*/ {6   , 0     , 163   , 7.56162  , 4.96626},
        /*10*/ {6   , 0     , 773   , 7.55908  , 4.96645},
        /*11*/ {6   , 1.25  , 191   , 8.23777  , 5.14309},
        /*12*/ {6   , 1.75  , 919   , 7.90075  , 5.04966},
        /*13*/ {7   , 0.5   , 853   , 8.6917   , 6.03185},
        /*14*/ {7   , 1     , 727   , 8.84793  , 6.07643},
        /*15*/ {7   , 1.75  , 193   , 9.85756  , 6.40897},
        /*16*/ {8   , 0.25  , 661   , 9.68683  , 7.00981},
        /*17*/ {8   , 0.5   , 331   , 9.83226  , 7.0566 },
        /*18*/ {8   , 1.5   , 643   , 10.1602  , 7.16693},
        /*19*/ {8   , 1.75  , 197   , 11.1179  , 7.53374},
        /*20*/ {8   , 2     , 277   , 10.9846  , 7.4787 },
        /*21*/ {8   , 2.25  , 359   , 10.9075  , 7.44735},
        /*22*/ {8   , 2.25  , 839   , 10.3245  , 7.22522},
        /*23*/ {8   , 2.25  , 1021  , 10.2353  , 7.19308},
        /*24*/ {9   , 0     , 89    , 10.6955  , 8.00178},
        /*25*/ {9   , 0.25  , 743   , 10.7245  , 8.01779},
        /*26*/ {9   , 1     , 1021  , 10.9371  , 8.08999},
        /*27*/ {9   , 1.75  , 229   , 12.2281  , 8.61095},
        /*28*/ {9   , 2.5   , 137   , 14.2547  , 9.63656},
        /*29*/ {10  , 0     , 359   , 11.707   , 8.9994 },
        /*30    10  , 1.5   , 32    , 20.9344  , 13.1556}, Spectral! */
        /*31*/ {10  , 2     , 353   , 13.2322  , 9.62829},
        /*32*/ {11  , 0.25  , 768   , 12.7973  , 10.0216},
        /*33*/ {11  , 0.5   , 419   , 12.9793  , 10.0897},
        /*34*/ {11  , 1.25  , 256   , 13.9588  , 10.5142},
        /*35*/ {11  , 1.25  , 640   , 13.3646  , 10.2449},
        /*36*/ {11  , 1.75  , 359   , 14.2029  , 10.6295},
        /*37*/ {11  , 2.25  , 97    , 18.5825  , 13.5221},
        /*38*/ {11  , 2.5   , 419   , 14.7094  , 10.886 },
        /*39*/ {11  , 2.75  , 271   , 15.8173  , 11.5164},
        /*40   {12  , 0     , 32    , 14.7039  , 10.9998}, Spectral! */
        /*41*/ {12  , 0     , 673   , 13.763   , 10.9999},
        /*42*/ {12  , 0.5   , 367   , 14.0612  , 11.1148},
        /*43*/ {12  , 0.5   , 557   , 13.9874  , 11.0843},
        /*44*/ {12  , 0.75  , 727   , 14.0898  , 11.1249},
        /*45*/ {12  , 1     , 883   , 14.1829  , 11.1631},
        /*46*/ {12  , 2     , 251   , 16.3478  , 12.2724},
        /*47*/ {12  , 2.5   , 397   , 16.0987  , 12.126 },
        /*48*/ {13  , 1.25  , 96    , 18.1713  , 13.8094},
        /*49*/ {13  , 1.5   , 251   , 16.8028  , 12.9783},
        /*50*/ {13  , 1.5   , 743   , 15.6825  , 12.385 },
        /*51*/ {13  , 2.5   , 101   , 23.3395  , 20.412 },
        /*52*/ {14  , 0.25  , 191   , 16.0141  , 13.0845},
        /*53*/ {14  , 0.5   , 677   , 16.0528  , 13.0951},
        /*54*/ {14  , 1     , 227   , 17.1979  , 13.655 },
        /*55*/ {14  , 2     , 347   , 18.333   , 14.3153},
        /*56*/ {14  , 2     , 739   , 17.2246  , 13.6636},
        /*57*/ {14  , 2.5   , 389   , 18.8008  , 14.6226},
        /*58*/ {14  , 3     , 83    , 30.7556  , 28.1666},
        /*59*/ {15  , 1     , 937   , 17.3792  , 14.2334},
        /*60*/ {15  , 2     , 479   , 19.0074  , 15.1276},
        /*61*/ {15  , 2.5   , 929   , 18.5206  , 14.8345},
        /*62*/ {16  , 0     , 1019  , 17.8576  , 15     },
        /*63*/ {16  , 0.25  , 941   , 17.944   , 15.0338},
        /*64*/ {16  , 0.75  , 367   , 18.6397  , 15.3596},
        /*65*/ {16  , 1.25  , 757   , 18.7903  , 15.4339},
        /*66*/ {16  , 2     , 571   , 19.9442  , 16.0989},
        /*67*/ {16  , 2.5   , 383   , 21.5722  , 17.2482},
        /*68*/ {17  , 0     , 523   , 18.8795  , 16     },
        /*69*/ {17  , 0.25  , 859   , 18.9776  , 16.0398},
        /*70*/ {17  , 0.75  , 181   , 20.31    , 16.7308},
        /*71*/ {17  , 1.25  , 293   , 20.9575  , 17.1175},
        /*72*/ {17  , 2.25  , 239   , 24.2356  , 19.7779},
        /*73*/ {17  , 3     , 479   , 23.0293  , 18.6634}
    };
    static const size_t ncases = sizeof(data)/sizeof(data[0]);

    for (size_t i = 0; i < ncases; ++i) {
        const int k          = (int) data[i][0];
        const double htdelta =       data[i][1];
        const int N          = (int) data[i][2];
        const double cc1     =       data[i][3];
        const double cc2     =       data[i][4];

        double C, Clow, Chigh;

        // For C^{(1)} we want...
        gsl_test(suzerain_bspline_htstretch2_evdeltascale_greville_abscissae(
                        1, k, htdelta, N, &C, &Clow, &Chigh),
                 "%s(%d, %d, %g, %d, ...)", __func__, 1, k, htdelta, N);

        gsl_test(!(Clow <= C && C <= Chigh),
                 "%s:%d self-consistency1 i=%d %g <= %g <= %g",
                 __FILE__, __LINE__, i, Clow, C, Chigh);

        gsl_test(!(Clow <= cc1 && cc1 <= Chigh),
                "%s:%d accuracy1 i=%d %g <= %g <= %g",
                 __FILE__, __LINE__, i, Clow, cc1, Chigh);

        // For C^{(2)} we want...
        gsl_test(suzerain_bspline_htstretch2_evdeltascale_greville_abscissae(
                        2, k, htdelta, N, &C, &Clow, &Chigh),
                 "%s(%d, %d, %g, %d, ...)", __func__, 2, k, htdelta, N);

        gsl_test(!(Clow <= C && C <= Chigh),
                 "%s:%d self-consistency2 i=%d %g <= %g <= %g",
                 __FILE__, __LINE__, i, Clow, C, Chigh);

        gsl_test(!(Clow <= cc2 && cc2 <= Chigh),
                 "%s:%d accuracy2 i=%d %g <= %g <= %g",
                 __FILE__, __LINE__, i, Clow, cc2, Chigh);

    }
}

// Test that curve fits reasonably reproduce the source data.
// Mainly meant as a sanity check on the coded coefficients.
// Test cases chosen pseudo-randomly from the original fit data.
// Choice required that we could nail at least one output per k.
static void test_bspline_htstretch1_evdeltascale_greville_abscissae()
{
    static const double data[][5] = {
        /*00*/{ 4, 0.00,  523, 5.44153, 2.72074},
        /*01*/{ 4, 0.50,  571, 5.53726, 2.75301},
        /*02*/{ 4, 0.75,  631, 5.58746, 2.76814},
        /*03*/{ 4, 2.50,  631, 5.88828, 2.85192},
        /*04*/{ 4, 2.75,  839, 5.83512, 2.83757},
        /*05*/{ 5, 2.00,  433, 7.16452, 4.20355},
        /*06*/{ 5, 2.25,  313, 7.41171, 4.26241},
        /*07*/{ 5, 2.50,  103, 8.75905, 4.5958 },
        /*08*/{ 5, 2.75,  191, 8.01886, 4.40985},
        /*09*/{ 5, 3.00,  283, 7.71521, 4.33556},
        /*10*/{ 6, 0.50,  929, 7.66702, 4.9902 },
        /*11*/{ 6, 1.50,  953, 7.91727, 5.0538 },
        /*12*/{ 6, 2.75,  443, 8.56929, 5.23537},
        /*13*/{ 6, 2.75,   67, 12.0328, 6.49388},
        /*14*/{ 7, 0.50,   64, 9.44346, 6.26963},
        /*15*/{ 7, 1.00,  619, 8.99354, 6.11986},
        /*16*/{ 7, 1.75, 1013, 9.07768, 6.14568},
        /*17*/{ 7, 1.75,  317, 9.6947 , 6.34944},
        /*18*/{ 7, 2.00, 1019, 9.13086, 6.16234},
        /*19*/{ 7, 2.00,  577, 9.39406, 6.24746},
        /*20*/{ 7, 2.75,  631, 9.5542 , 6.30105},
        /*21*/{ 7, 3.00,  853, 9.41929, 6.25575},
        /*22*/{ 8, 0.25,  449, 9.73924, 7.02597},
        /*23*/{ 8, 0.25,  800, 9.70722, 7.01586},
        /*24*/{ 8, 0.50,  557, 9.8547 , 7.06287},
        /*25*/{ 8, 1.25,  773, 10.1364, 7.158  },
        /*26*/{ 8, 1.50,  199, 11.2443, 7.58351},
        /*27*/{ 8, 1.50,  211, 11.174 , 7.5544 },
        /*28*/{ 8, 1.50,   83, 12.8824, 8.34102},
        /*29*/{ 8, 2.50,   48, 18.3515, 13.5572},
        /*30*/{ 8, 2.50,  569, 10.7747, 7.394  },
        /*31*/{ 9, 0.50,  631, 10.9025, 8.07774},
        /*32*/{ 9, 1.75,  311, 12.2285, 8.60864},
        /*33*/{ 9, 2.00,  463, 11.9688, 8.49518},
        /*34*/{ 9, 3.00,  967, 11.7089, 8.38593},
        /*35*/{10, 0.50,  487, 12.0177, 9.1099 },
        /*36*/{10, 0.75,  149, 12.9502, 9.49965},
        /*37*/{10, 1.00,  941, 12.1611, 9.16431},
        /*38*/{10, 1.50,  797, 12.4811, 9.29331},
        /*39*/{10, 2.00,  179, 14.8647, 10.4886},
        /*40*/{10, 2.00,  307, 13.7581, 9.88353},
        /*41*/{10, 2.00,  739, 12.7627, 9.41349},
        /*42*/{11, 0.50,  929, 12.9597, 10.0807},
        /*43*/{11, 0.75,  911, 13.1135, 10.1406},
        /*44*/{11, 1.25,  277, 14.3534, 10.7028},
        /*45*/{11, 1.50,  139, 16.1456, 11.723 },
        /*46*/{11, 2.75,  229, 16.6575, 12.0705},
        /*47*/{12, 0.00,  251, 13.7662, 10.9999},
        /*48*/{12, 0.25, 1019, 13.856 , 11.0331},
        /*49*/{12, 0.50,  557, 14.1221, 11.1377},
        /*50*/{12, 1.00,  883, 14.369 , 11.2421},
        /*51*/{12, 1.25,  739, 14.6393, 11.3639},
        /*52*/{12, 2.00,  673, 15.2226, 11.6479},
        /*53*/{12, 2.00,  811, 15.0321, 11.552 },
        /*54*/{12, 2.75,  461, 16.3028, 12.246 },
        /*55*/{12, 3.00,  349, 17.1661, 12.7956},
        /*56*/{12, 3.00,  433, 16.6248, 12.4431},
        /*57*/{13, 0.00,  571, 14.789 , 12     },
        /*58*/{13, 0.50,  751, 15.1119, 12.1265},
        /*59*/{13, 0.50,  919, 15.069 , 12.1085},
        /*60*/{13, 0.75,  571, 15.4487, 12.2747},
        /*61*/{13, 1.25,  491, 16.1169, 12.6012},
        /*62*/{13, 2.00,  256, 18.296 , 13.933 },
        /*63*/{13, 2.25,  613, 16.7307, 12.9345},
        /*64*/{14, 0.00,  512, 15.8136, 13     },
        /*65*/{14, 0.00,  853, 15.8128, 13     },
        /*66*/{14, 0.50,  797, 16.1548, 13.1376},
        /*67*/{14, 0.50,  977, 16.109 , 13.1178},
        /*68*/{14, 0.75,  641, 16.4847, 13.2874},
        /*69*/{14, 1.25,  401, 17.5379, 13.837 },
        /*70*/{14, 1.50,  919, 16.9218, 13.5028},
        /*71*/{14, 2.00,  461, 18.2572, 14.2677},
        /*72*/{14, 2.50,  953, 17.4952, 13.8121},
        /*73*/{15, 0.00,  349, 16.8382, 14     },
        /*74*/{15, 2.25,  571, 19.3718, 15.3602},
        /*75*/{15, 3.00,  101, 31.6167, 29.2606},
        /*76*/{15, 3.00, 1021, 18.8744, 15.0447},
        /*77*/{15, 3.00,  419, 20.9666, 16.5416},
        /*78*/{16, 0.25,  181, 18.3123, 15.1995},
        /*79*/{16, 0.25,  691, 18.034 , 15.0708},
        /*80*/{16, 0.25,  768, 18.0215, 15.0655},
        /*81*/{16, 0.75,  631, 18.6763, 15.3741},
        /*82*/{16, 1.00,  192, 20.7962, 16.6665},
        /*83*/{16, 1.00,  787, 18.8393, 15.4583},
        /*84*/{16, 2.00,  499, 20.6654, 16.5773},
        /*85*/{16, 2.50,  103, 31.7243, 29.2853},
        /*86*/{16, 3.00, 1021, 20.1004, 16.1978},
        /*87*/{17, 0.75, 1024, 19.504 , 16.282 },
        /*88*/{17, 1.25,  853, 20.167 , 16.6365},
        /*89*/{17, 1.50,  127, 26.4158, 23.6473},
        /*90*/{17, 1.75,  443, 21.8753, 17.7455},
        /*91*/{17, 2.00, 1009, 20.6553, 16.925 },
        /*92*/{17, 2.00,  439, 22.2945, 18.0627},
        /*93*/{17, 2.00,  977, 20.6996, 16.9524},
        /*94*/{17, 2.75,  457, 23.2321, 18.8513},
        /*95*/{17, 3.00, 1024, 21.3283, 17.36  }
    };
    static const size_t ncases = sizeof(data)/sizeof(data[0]);

    for (size_t i = 0; i < ncases; ++i) {
        const int k          = (int) data[i][0];
        const double htdelta =       data[i][1];
        const int N          = (int) data[i][2];
        const double cc1     =       data[i][3];
        const double cc2     =       data[i][4];

        double C, Clow, Chigh;

        // For C^{(1)} we want...
        gsl_test(suzerain_bspline_htstretch1_evdeltascale_greville_abscissae(
                        1, k, htdelta, N, &C, &Clow, &Chigh),
                 "%s(%d, %d, %g, %d, ...)", __func__, 1, k, htdelta, N);

        gsl_test(!(Clow <= C && C <= Chigh),
                 "%s:%d self-consistency1 i=%d %g <= %g <= %g",
                 __FILE__, __LINE__, i, Clow, C, Chigh);

        gsl_test(!(Clow <= cc1 && cc1 <= Chigh),
                "%s:%d accuracy1 i=%d %g <= %g <= %g",
                 __FILE__, __LINE__, i, Clow, cc1, Chigh);

        // For C^{(2)} we want...
        gsl_test(suzerain_bspline_htstretch1_evdeltascale_greville_abscissae(
                        2, k, htdelta, N, &C, &Clow, &Chigh),
                 "%s(%d, %d, %g, %d, ...)", __func__, 2, k, htdelta, N);

        gsl_test(!(Clow <= C && C <= Chigh),
                 "%s:%d self-consistency2 i=%d %g <= %g <= %g",
                 __FILE__, __LINE__, i, Clow, C, Chigh);

        gsl_test(!(Clow <= cc2 && cc2 <= Chigh),
                 "%s:%d accuracy2 i=%d %g <= %g <= %g",
                 __FILE__, __LINE__, i, Clow, cc2, Chigh);

    }
}

// Test that curve fits reasonably reproduce the source data.
// Mainly meant as a sanity check on the coded coefficients.
// Test cases chosen pseudo-randomly from the original fit data.
// Choice required that we could nail at least one output per k.
static void test_bspline_htstretch2_evdeltascale_breakpoints()
{
    static const double data[][5] = {
        /*  0*/ {  4,   0.25,  613,   1.82154,  0.909792  },
        /*  1*/ {  4,   0.5,   600,   1.8337,   0.913887  },
        /*  2*/ {  4,   1.25,  73,    2.08151,  0.983002  },
        /*  3*/ {  4,   1,     487,   1.87,     0.924979  },
        /*  4*/ {  5,   0,     421,   1.63094,  1.01392   },
        /*  5*/ {  5,   0.75,  311,   1.68115,  1.02535   },
        /*  6*/ {  5,   2.25,  421,   1.78693,  1.04998   },
        /*  7*/ {  5,   2.25,  761,   1.73416,  1.03752   },
        /*  8*/ {  5,   2.5,   864,   1.73608,  1.03795   },
        /*  9*/ {  5,   2.75,  439,   1.81642,  1.05695   },
        /* 10*/ {  5,   2.75,  521,   1.79523,  1.05191   },
        /* 11*/ {  5,   3,     709,   1.77445,  1.04697   },
        /* 12*/ {  6,   0.25,  223,   1.52525,  0.996263  },
        /* 13*/ {  6,   0.5,   653,   1.52933,  0.997101  },
        /* 14*/ {  6,   0.75,  593,   1.5439,   1.00067   },
        /* 15*/ {  6,   1.25,  401,   1.5928,   1.01339   },
        /* 16*/ {  6,   1.25,  751,   1.56446,  1.00587   },
        /* 17*/ {  6,   1.5,   61,    1.91163,  1.10975   },
        /* 18*/ {  6,   1.75,  73,    1.93902,  1.11826   },
        /* 19*/ {  6,   2,     641,   1.61368,  1.01892   },
        /* 20*/ {  6,   2.75,  192,   1.85428,  1.09      },
        /* 21*/ {  7,   0.25,  281,   1.4456,   1.0046    },
        /* 22*/ {  7,   0.25,  881,   1.43961,  1.00292   },
        /* 23*/ {  7,   0.5,   821,   1.44899,  1.00541   },
        /* 24*/ {  7,   0.75,  113,   1.53257,  1.0315    },
        /* 25*/ {  7,   2.75,  173,   1.81471,  1.13312   },
        /* 26*/ {  7,   2.75,  433,   1.62428,  1.06126   },
        /* 27*/ {  7,   2,     863,   1.51724,  1.02567   },
        /* 28*/ {  7,   3,     829,   1.56328,  1.04048   },
        /* 29*/ {  8,   0,     167,   1.37772,  0.999359  },
        /* 30*/ {  8,   0.25,  751,   1.3833,   1.00123   },
        /* 31*/ {  8,   0.5,   157,   1.42208,  1.01424   },
        /* 32*/ {  8,   0.75,  47,    1.56092,  1.06841   },
        /* 33*/ {  8,   0.75,  83,    1.5011,   1.04353   },
        /* 34*/ {  8,   1,     43,    1.68923,  1.1204    },
        /* 35*/ {  8,   1.5,   709,   1.44652,  1.02211   },
        /* 36*/ {  8,   1.75,  239,   1.55991,  1.06476   },
        /* 37*/ {  8,   2.5,   719,   1.49907,  1.04105   },
        /* 38*/ {  9,   0,     523,   1.33438,  1.0002    },
        /* 39*/ {  9,   1.25,  683,   1.39209,  1.02032   },
        /* 40*/ {  9,   1.5,   239,   1.48916,  1.05932   },
        /* 41*/ {  9,   1.5,   727,   1.40369,  1.02465   },
        /* 42*/ {  9,   1.75,  859,   1.40839,  1.02642   },
        /* 43*/ {  9,   2,     331,   1.50726,  1.0669    },
        /* 44*/ {  9,   2,     503,   1.4613,   1.04743   },
        /* 45*/ {  9,   2,     593,   1.44696,  1.04158   },
        /* 46*/ {  9,   2.75,  107,   1.95121,  1.31153   },
        /* 47*/ {  9,   3,     311,   1.61623,  1.11697   },
        /* 48*/ {  9,   3,     541,   1.51911,  1.07193   },
        /* 49*/ {  10,  0,     197,   1.30113,  0.999936  },
        /* 50*/ {  10,  0.25,  149,   1.31768,  1.00628   },
        /* 51*/ {  10,  0.25,  443,   1.30943,  1.00294   },
        /* 52*/ {  10,  0.25,  83,    1.3221,   1.00904   },
        /* 53*/ {  10,  0.25,  911,   1.30615,  1.00176   },
        /* 54*/ {  10,  1.25,  512,   1.37288,  1.02702   },
        /* 55*/ {  10,  1.25,  541,   1.3701,   1.02588   },
        /* 56*/ {  10,  1.5,   64,    1.76244,  1.22675   },
        /* 57*/ {  10,  1.5,   881,   1.36238,  1.02271   },
        /* 58*/ {  10,  1.5,   937,   1.35975,  1.02165   },
        /* 59*/ {  10,  1.75,  599,   1.39868,  1.0377    },
        /* 60*/ {  10,  1.75,  613,   1.39705,  1.03701   },
        /* 61*/ {  10,  2.25,  293,   1.52473,  1.09629   },
        /* 62*/ {  10,  2.25,  859,   1.40105,  1.03865   },
        /* 63*/ {  10,  2.75,  431,   1.50785,  1.08781   },
        /* 64*/ {  10,  2,     919,   1.38463,  1.03174   },
        /* 65*/ {  10,  3,     433,   1.52617,  1.09693   },
        /* 66*/ {  10,  3,     773,   1.44623,  1.05855   },
        /* 67*/ {  11,  0.25,  827,   1.27944,  1.00206   },
        /* 68*/ {  11,  1.5,   139,   1.52376,  1.11748   },
        /* 69*/ {  11,  1.5,   631,   1.35362,  1.03188   },
        /* 70*/ {  11,  1,     593,   1.32311,  1.01893   },
        /* 71*/ {  11,  2.25,  1021,  1.36419,  1.03648   },
        /* 72*/ {  11,  2.25,  739,   1.38841,  1.04757   },
        /* 73*/ {  11,  2.5,   647,   1.41556,  1.0605    },
        /* 74*/ {  11,  2,     719,   1.37627,  1.04198   },
        /* 75*/ {  11,  2.75,  281,   1.57288,  1.14626   },
        /* 76*/ {  11,  2.75,  983,   1.38838,  1.04751   },
        /* 77*/ {  11,  3,     631,   1.44789,  1.0766    },
        /* 78*/ {  11,  3,     941,   1.40272,  1.05428   },
        /* 79*/ {  12,  0.5,   193,   1.2931,   1.01693   },
        /* 80*/ {  12,  0,     563,   1.2512,   0.999992  },
        /* 81*/ {  12,  0,     64,    1.25967,  0.99999   },
        /* 82*/ {  12,  0.75,  757,   1.28007,  1.01102   },
        /* 83*/ {  12,  1.5,   181,   1.46117,  1.10112   },
        /* 84*/ {  12,  2.5,   859,   1.36925,  1.05157   },
        /* 85*/ {  13,  0.25,  911,   1.23804,  1.00207   },
        /* 86*/ {  13,  0.5,   960,   1.24669,  1.00541   },
        /* 87*/ {  13,  0.75,  593,   1.26729,  1.01403   },
        /* 88*/ {  13,  1.5,   761,   1.30559,  1.03146   },
        /* 89*/ {  13,  2.5,   283,   1.51981,  1.15724   },
        /* 90*/ {  13,  2.75,  317,   1.52373,  1.16012   },
        /* 91*/ {  13,  2.75,  521,   1.42867,  1.09809   },
        /* 92*/ {  13,  2.75,  743,   1.38177,  1.0709    },
        /* 93*/ {  13,  3,     191,   1.72314,  1.33138   },
        /* 94*/ {  14,  0.75,  521,   1.25536,  1.01636   },
        /* 95*/ {  14,  0.75,  593,   1.25195,  1.01479   },
        /* 96*/ {  14,  0.75,  89,    1.35687,  1.07011   },
        /* 97*/ {  14,  1.25,  503,   1.29658,  1.03626   },
        /* 98*/ {  14,  1.25,  727,   1.27767,  1.02679   },
        /* 99*/ {  14,  1.5,   853,   1.28529,  1.03049   },
        /*100*/ {  14,  2,     307,   1.42998,  1.11389   },
        /*101*/ {  14,  2.5,   929,   1.33367,  1.05576   },
        /*102*/ {  14,  2,     631,   1.33867,  1.05858   },
        /*103*/ {  14,  3,     313,   1.55098,  1.20354   },
        /*104*/ {  14,  3,     379,   1.50224,  1.165     },
        /*105*/ {  14,  3,     607,   1.4127,   1.10278   },
        /*106*/ {  15,  0.25,  83,    1.22584,  1.01084   },
        /*107*/ {  15,  0.5,   864,   1.21843,  1.00638   },
        /*108*/ {  15,  0,     739,   1.20257,  1         },
        /*109*/ {  15,  1.5,   384,   1.3306,   1.06425   },
        /*110*/ {  15,  1.75,  401,   1.35358,  1.07803   },
        /*111*/ {  15,  1.75,  953,   1.28115,  1.03652   },
        /*112*/ {  15,  2.25,  229,   1.52696,  1.20454   },
        /*113*/ {  15,  2.25,  256,   1.49796,  1.18037   },
        /*114*/ {  15,  2.25,  337,   1.43829,  1.13469   },
        /*115*/ {  15,  2.25,  349,   1.43172,  1.12997   },
        /*116*/ {  15,  2.25,  439,   1.39342,  1.10354   },
        /*117*/ {  15,  2.25,  859,   1.31634,  1.05584   },
        /*118*/ {  15,  2.5,   107,   1.94637,  1.75001   },
        /*119*/ {  15,  2,     829,   1.30471,  1.04929   },
        /*120*/ {  16,  0.25,  283,   1.20307,  1.00533   },
        /*121*/ {  16,  0.25,  541,   1.1988,   1.00335   },
        /*122*/ {  16,  0.5,   101,   1.26117,  1.03519   },
        /*123*/ {  16,  0.5,   941,   1.20571,  1.00623   },
        /*124*/ {  16,  0.75,  113,   1.31482,  1.06498   },
        /*125*/ {  16,  0,     757,   1.19052,  1         },
        /*126*/ {  16,  1.75,  384,   1.35109,  1.08707   },
        /*127*/ {  16,  2,     317,   1.41188,  1.12947   },
        /*128*/ {  16,  2.75,  569,   1.39056,  1.1142    },
        /*129*/ {  16,  3,     541,   1.4187,   1.135     },
        /*130*/ {  16,  3,     911,   1.34191,  1.08107   },
        /*131*/ {  17,  0.25,  149,   1.19859,  1.00869   },
        /*132*/ {  17,  0.5,   107,   1.24942,  1.03538   },
        /*133*/ {  17,  0.5,   173,   1.22993,  1.02413   },
        /*134*/ {  17,  0.5,   743,   1.1981,   1.00773   },
        /*135*/ {  17,  0.75,  773,   1.21112,  1.01385   },
        /*136*/ {  17,  1,     137,   1.35151,  1.09707   },
        /*137*/ {  17,  1.25,  277,   1.31565,  1.07357   },
        /*138*/ {  17,  1.25,  709,   1.24657,  1.03217   },
        /*139*/ {  17,  1.25,  919,   1.23501,  1.02592   },
        /*140*/ {  17,  1.5,   967,   1.24707,  1.03239   },
        /*141*/ {  17,  2.25,  103,   1.923,    1.75354   },
        /*142*/ {  17,  2.25,  907,   1.29462,  1.06012   },
        /*143*/ {  17,  2.5,   449,   1.40409,  1.13703   },
        /*144*/ {  17,  2,     823,   1.28773,  1.0559    }
    };
    static const size_t ncases = sizeof(data)/sizeof(data[0]);

    for (size_t i = 0; i < ncases; ++i) {
        const int k          = (int) data[i][0];
        const double htdelta =       data[i][1];
        const int N          = (int) data[i][2];
        const double bc1     =       data[i][3];
        const double bc2     =       data[i][4];

        double C, Clow, Chigh;

        // For C^{(1)} we want...
        gsl_test(suzerain_bspline_htstretch2_evdeltascale_breakpoints(
                        1, k, htdelta, N, &C, &Clow, &Chigh),
                 "%s(%d, %d, %g, %d, ...)", __func__, 1, k, htdelta, N);

        gsl_test(!(Clow <= C && C <= Chigh),
                 "%s:%d self-consistency1 i=%d %g <= %g <= %g",
                 __FILE__, __LINE__, i, Clow, C, Chigh);

        gsl_test(!(Clow <= bc1 && bc1 <= Chigh),
                "%s:%d accuracy1 i=%d %g <= %g <= %g",
                 __FILE__, __LINE__, i, Clow, bc1, Chigh);

        // For C^{(2)} we want...
        gsl_test(suzerain_bspline_htstretch2_evdeltascale_breakpoints(
                        2, k, htdelta, N, &C, &Clow, &Chigh),
                 "%s(%d, %d, %g, %d, ...)", __func__, 2, k, htdelta, N);

        gsl_test(!(Clow <= C && C <= Chigh),
                 "%s:%d self-consistency2 i=%d %g <= %g <= %g",
                 __FILE__, __LINE__, i, Clow, C, Chigh);

        gsl_test(!(Clow <= bc2 && bc2 <= Chigh),
                 "%s:%d accuracy2 i=%d %g <= %g <= %g",
                 __FILE__, __LINE__, i, Clow, bc2, Chigh);

    }
}

// Test that curve fits reasonably reproduce the source data.
// Mainly meant as a sanity check on the coded coefficients.
// Test cases chosen pseudo-randomly from the original fit data.
// Choice required that we could nail at least one output per k.
static void test_bspline_htstretch1_evdeltascale_breakpoints()
{
    static const double data[][5] = {
        /*  0*/  {  4,   0.25,  751,   1.82471,  0.910874  },
        /*  1*/  {  4,   0.50,  401,   1.85424,  0.920272  },
        /*  2*/  {  4,   0.75,  479,   1.87238,  0.925625  },
        /*  3*/  {  4,   0.75,  864,   1.85325,  0.919945  },
        /*  4*/  {  4,   1.00,  283,   1.92899,  0.941544  },
        /*  5*/  {  4,   1.00,  977,   1.86363,  0.923041  },
        /*  6*/  {  4,   1.50,  467,   1.93541,  0.943258  },
        /*  7*/  {  4,   1.75,  163,   2.10351,  0.987787  },
        /*  8*/  {  4,   1.75,  491,   1.94797,  0.946662  },
        /*  9*/  {  4,   1.75,  600,   1.93076,  0.941979  },
        /* 10*/  {  4,   2.00,  96,    2.29091,  1.03649   },
        /* 11*/  {  4,   2.75,  479,   2.00695,  0.962418  },
        /* 12*/  {  5,   0.00,  359,   1.63097,  1.01392   },
        /* 13*/  {  5,   0.25,  929,   1.63987,  1.01583   },
        /* 14*/  {  5,   0.50,  431,   1.66788,  1.02218   },
        /* 15*/  {  5,   0.75,  211,   1.7302,   1.03668   },
        /* 16*/  {  5,   0.75,  937,   1.66672,  1.02189   },
        /* 17*/  {  5,   1.25,  107,   1.91772,  1.08142   },
        /* 18*/  {  5,   1.50,  127,   1.9339,   1.08526   },
        /* 19*/  {  5,   2.00,  883,   1.7282,   1.03608   },
        /* 20*/  {  5,   2.00,  89,    2.15004,  1.13887   },
        /* 21*/  {  5,   2.25,  359,   1.83205,  1.06061   },
        /* 22*/  {  5,   2.25,  563,   1.77699,  1.04753   },
        /* 23*/  {  5,   2.25,  953,   1.73203,  1.03697   },
        /* 24*/  {  5,   2.50,  960,   1.7398,   1.03879   },
        /* 25*/  {  5,   2.75,  139,   2.10771,  1.1281    },
        /* 26*/  {  5,   2.75,  941,   1.74931,  1.04102   },
        /* 27*/  {  6,   0.25,  233,   1.53386,  0.998225  },
        /* 28*/  {  6,   0.75,  103,   1.67688,  1.03684   },
        /* 29*/  {  6,   2.00,  41,    2.55764,  1.37052   },
        /* 30*/  {  6,   3.00,  131,   2.05985,  1.15899   },
        /* 31*/  {  7,   0.25,  151,   1.46346,  1.00962   },
        /* 32*/  {  7,   0.50,  48,    1.6077,   1.05689   },
        /* 33*/  {  7,   1.50,  53,    2.08559,  1.25466   },
        /* 34*/  {  7,   2.25,  919,   1.53771,  1.03212   },
        /* 35*/  {  7,   2.50,  787,   1.55963,  1.03924   },
        /* 36*/  {  7,   2.75,  283,   1.72132,  1.09641   },
        /* 37*/  {  7,   3.00,  353,   1.6939,   1.08616   },
        /* 38*/  {  7,   3.00,  907,   1.56402,  1.04068   },
        /* 39*/  {  8,   0.00,  271,   1.37734,  0.999357  },
        /* 40*/  {  8,   0.00,  96,    1.37911,  0.999363  },
        /* 41*/  {  8,   0.75,  733,   1.41911,  1.01269   },
        /* 42*/  {  8,   1.00,  1013,  1.42363,  1.01419   },
        /* 43*/  {  8,   1.25,  1009,  1.43591,  1.01836   },
        /* 44*/  {  8,   1.25,  67,    1.83324,  1.18772   },
        /* 45*/  {  8,   1.50,  241,   1.57523,  1.07062   },
        /* 46*/  {  8,   1.50,  479,   1.49632,  1.03999   },
        /* 47*/  {  8,   1.75,  673,   1.48426,  1.03553   },
        /* 48*/  {  8,   2.50,  349,   1.61075,  1.08512   },
        /* 49*/  {  8,   3.00,  640,   1.54801,  1.05971   },
        /* 50*/  {  9,   0.00,  643,   1.33436,  1.0002    },
        /* 51*/  {  9,   0.25,  64,    1.38863,  1.01994   },
        /* 52*/  {  9,   0.75,  653,   1.38049,  1.016     },
        /* 53*/  {  9,   1.00,  47,    1.8395,   1.23653   },
        /* 54*/  {  9,   1.00,  947,   1.38379,  1.01718   },
        /* 55*/  {  9,   1.25,  521,   1.4302,   1.03483   },
        /* 56*/  {  9,   1.75,  384,   1.49997,  1.06357   },
        /* 57*/  {  9,   2.00,  587,   1.46998,  1.05089   },
        /* 58*/  {  9,   2.00,  600,   1.46781,  1.04999   },
        /* 59*/  {  9,   2.50,  137,   1.84635,  1.24477   },
        /* 60*/  {  9,   2.50,  499,   1.51787,  1.07133   },
        /* 61*/  {  9,   2.50,  683,   1.47959,  1.05488   },
        /* 62*/  {  10,  0.25,  79,    1.34864,  1.01801   },
        /* 63*/  {  10,  0.50,  64,    1.4547,   1.06348   },
        /* 64*/  {  10,  0.75,  227,   1.40139,  1.03893   },
        /* 65*/  {  10,  1.00,  577,   1.37238,  1.0267    },
        /* 66*/  {  10,  1.00,  739,   1.36072,  1.022     },
        /* 67*/  {  10,  1.75,  509,   1.43855,  1.05505   },
        /* 68*/  {  10,  2.25,  647,   1.44344,  1.05724   },
        /* 69*/  {  10,  2.75,  167,   1.79141,  1.25622   },
        /* 70*/  {  10,  3.00,  384,   1.5681,   1.11862   },
        /* 71*/  {  11,  0.00,  283,   1.27376,  1.00002   },
        /* 72*/  {  11,  0.00,  499,   1.27359,  1.00002   },
        /* 73*/  {  11,  0.00,  733,   1.27355,  1.00002   },
        /* 74*/  {  11,  0.00,  881,   1.27354,  1.00002   },
        /* 75*/  {  11,  0.25,  257,   1.29506,  1.00784   },
        /* 76*/  {  11,  0.75,  223,   1.37834,  1.04305   },
        /* 77*/  {  11,  0.75,  233,   1.37495,  1.04149   },
        /* 78*/  {  11,  1.00,  197,   1.43739,  1.07139   },
        /* 79*/  {  11,  1.00,  911,   1.32631,  1.02013   },
        /* 80*/  {  11,  1.50,  839,   1.35844,  1.03389   },
        /* 81*/  {  11,  1.75,  199,   1.5688,   1.14387   },
        /* 82*/  {  11,  1.75,  256,   1.51458,  1.1123    },
        /* 83*/  {  11,  1.75,  383,   1.4495,   1.0774    },
        /* 84*/  {  11,  2.00,  251,   1.55025,  1.13282   },
        /* 85*/  {  12,  0.25,  853,   1.26074,  1.00342   },
        /* 86*/  {  12,  0.75,  661,   1.29961,  1.01914   },
        /* 87*/  {  12,  1.00,  167,   1.44398,  1.09127   },
        /* 88*/  {  12,  1.00,  307,   1.37119,  1.05259   },
        /* 89*/  {  12,  1.00,  480,   1.33708,  1.03603   },
        /* 90*/  {  12,  1.00,  827,   1.30893,  1.02318   },
        /* 91*/  {  12,  1.75,  1019,  1.33798,  1.03639   },
        /* 92*/  {  12,  1.75,  191,   1.5682,   1.16826   },
        /* 93*/  {  12,  1.75,  677,   1.3687,   1.05125   },
        /* 94*/  {  12,  1.75,  853,   1.35014,  1.04217   },
        /* 95*/  {  12,  1.75,  997,   1.33939,  1.03705   },
        /* 96*/  {  12,  2.25,  911,   1.36783,  1.05081   },
        /* 97*/  {  12,  2.75,  907,   1.38842,  1.06122   },
        /* 98*/  {  13,  0.00,  809,   1.23238,  1         },
        /* 99*/  {  13,  0.25,  619,   1.24449,  1.00453   },
        /*100*/  {  13,  1.00,  523,   1.31501,  1.03597   },
        /*101*/  {  13,  1.75,  349,   1.4336,   1.10107   },
        /*102*/  {  13,  2.00,  829,   1.34889,  1.05309   },
        /*103*/  {  13,  2.25,  647,   1.38765,  1.07415   },
        /*104*/  {  13,  2.25,  971,   1.34654,  1.05186   },
        /*105*/  {  13,  2.75,  241,   1.63847,  1.25194   },
        /*106*/  {  13,  3.00,  991,   1.37351,  1.06629   },
        /*107*/  {  14,  0.25,  631,   1.2285,   1.00466   },
        /*108*/  {  14,  1.00,  257,   1.36213,  1.07194   },
        /*109*/  {  14,  1.00,  739,   1.28192,  1.02877   },
        /*110*/  {  14,  1.25,  577,   1.31695,  1.04669   },
        /*111*/  {  14,  1.25,  919,   1.28751,  1.03151   },
        /*112*/  {  14,  1.25,  96,    1.6555,   1.29816   },
        /*113*/  {  14,  2.25,  179,   1.67462,  1.32575   },
        /*114*/  {  14,  2.75,  439,   1.47244,  1.14333   },
        /*115*/  {  15,  0.50,  229,   1.26906,  1.03043   },
        /*116*/  {  15,  0.75,  337,   1.28774,  1.04008   },
        /*117*/  {  15,  0.75,  479,   1.26801,  1.02972   },
        /*118*/  {  15,  0.75,  881,   1.24446,  1.01807   },
        /*119*/  {  15,  0.75,  907,   1.24359,  1.01766   },
        /*120*/  {  15,  1.50,  128,   1.63417,  1.31465   },
        /*121*/  {  15,  1.50,  541,   1.33312,  1.06558   },
        /*122*/  {  15,  2.00,  541,   1.37325,  1.09036   },
        /*123*/  {  15,  2.00,  857,   1.32222,  1.05918   },
        /*124*/  {  15,  2.50,  521,   1.41535,  1.11864   },
        /*125*/  {  15,  3.00,  463,   1.47444,  1.16269   },
        /*126*/  {  16,  0.00,  359,   1.19068,  1         },
        /*127*/  {  16,  0.00,  619,   1.19055,  1         },
        /*128*/  {  16,  0.00,  997,   1.19051,  1         },
        /*129*/  {  16,  0.25,  467,   1.20596,  1.00634   },
        /*130*/  {  16,  0.75,  331,   1.27903,  1.04321   },
        /*131*/  {  16,  1.50,  233,   1.45492,  1.16369   },
        /*132*/  {  16,  2.00,  619,   1.34841,  1.08531   },
        /*133*/  {  16,  2.00,  800,   1.31992,  1.06717   },
        /*134*/  {  16,  2.25,  947,   1.31618,  1.06487   },
        /*135*/  {  17,  0.25,  431,   1.19653,  1.00701   },
        /*136*/  {  17,  0.50,  960,   1.20423,  1.01048   },
        /*137*/  {  17,  1.00,  127,   1.4695,   1.1915    },
        /*138*/  {  17,  1.00,  277,   1.32945,  1.08258   },
        /*139*/  {  17,  1.75,  271,   1.46076,  1.18585   },
        /*140*/  {  17,  1.75,  600,   1.32714,  1.08103   },
        /*141*/  {  17,  2.00,  317,   1.45926,  1.18465   },
        /*142*/  {  17,  2.25,  347,   1.46794,  1.19285   },
        /*143*/  {  17,  2.25,  577,   1.36995,  1.11114   },
        /*144*/  {  17,  2.75,  491,   1.43639,  1.16436   }
    };
    static const size_t ncases = sizeof(data)/sizeof(data[0]);

    for (size_t i = 0; i < ncases; ++i) {
        const int k          = (int) data[i][0];
        const double htdelta =       data[i][1];
        const int N          = (int) data[i][2];
        const double bc1     =       data[i][3];
        const double bc2     =       data[i][4];

        double C, Clow, Chigh;

        // For C^{(1)} we want...
        gsl_test(suzerain_bspline_htstretch1_evdeltascale_breakpoints(
                        1, k, htdelta, N, &C, &Clow, &Chigh),
                 "%s(%d, %d, %g, %d, ...)", __func__, 1, k, htdelta, N);

        gsl_test(!(Clow <= C && C <= Chigh),
                 "%s:%d self-consistency1 i=%d %g <= %g <= %g",
                 __FILE__, __LINE__, i, Clow, C, Chigh);

        gsl_test(!(Clow <= bc1 && bc1 <= Chigh),
                "%s:%d accuracy1 i=%d %g <= %g <= %g",
                 __FILE__, __LINE__, i, Clow, bc1, Chigh);

        // For C^{(2)} we want...
        gsl_test(suzerain_bspline_htstretch1_evdeltascale_breakpoints(
                        2, k, htdelta, N, &C, &Clow, &Chigh),
                 "%s(%d, %d, %g, %d, ...)", __func__, 2, k, htdelta, N);

        gsl_test(!(Clow <= C && C <= Chigh),
                 "%s:%d self-consistency2 i=%d %g <= %g <= %g",
                 __FILE__, __LINE__, i, Clow, C, Chigh);

        gsl_test(!(Clow <= bc2 && bc2 <= Chigh),
                 "%s:%d accuracy2 i=%d %g <= %g <= %g",
                 __FILE__, __LINE__, i, Clow, bc2, Chigh);

    }
}
