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

        gsl_vector_set_all(coeffs, -555.0);
        gsl_test(suzerain_bspline_integration_coefficients(
                0 , gsl_vector_ptr(coeffs,0), inc, scratch, w, dw),
                "integration_coefficients 0 linear");
        for (size_t i = 0; i < w->n; ++i) {
            gsl_test_rel(gsl_vector_get(coeffs, 2*i), e0[2*i],
                    GSL_DBL_EPSILON*1000,
                    "integration_coefficients 0 linear %d value", 2*i);
            gsl_test(gsl_vector_get(coeffs, 2*i+1) != e0[2*i+1],
                    "integration_coefficients 0 linear %d stride", 2*i + 1);
        }

        gsl_vector_set_all(coeffs, -555.0);
        gsl_test(suzerain_bspline_integration_coefficients(
                1 , gsl_vector_ptr(coeffs,0), inc, scratch, w, dw),
                "integration_coefficients 1 linear");
        for (size_t i = 0; i < w->n; ++i) {
            gsl_test_rel(gsl_vector_get(coeffs, 2*i), e1[2*i],
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
                0 , gsl_vector_ptr(coeffs,0), inc, scratch, w, dw),
                "integration_coefficients 0 quadratic");
        for (size_t i = 0; i < w->n; ++i) {
            gsl_test_rel(gsl_vector_get(coeffs, i), e0[i],
                    GSL_DBL_EPSILON*1000,
                    "integration_coefficients 0 quadratic %d value", i);
        }

        gsl_vector_set_all(coeffs, -555.0);
        gsl_test(suzerain_bspline_integration_coefficients(
                1 , gsl_vector_ptr(coeffs,0), inc, scratch, w, dw),
                "integration_coefficients 1 quadratic");
        for (size_t i = 0; i < w->n; ++i) {
            gsl_test_rel(gsl_vector_get(coeffs, i), e1[i],
                    GSL_DBL_EPSILON*1000,
                    "integration_coefficients 1 quadratic %d value", i);
        }

        gsl_vector_set_all(coeffs, -555.0);
        gsl_test(suzerain_bspline_integration_coefficients(
                2 , gsl_vector_ptr(coeffs,0), inc, scratch, w, dw),
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
                0 , gsl_vector_ptr(coeffs,0), inc, scratch, w, dw),
                "integration_coefficients 0 cubic");
        for (size_t i = 0; i < w->n; ++i) {
            gsl_test_rel(gsl_vector_get(coeffs, i), e0[i],
                    GSL_DBL_EPSILON*1000,
                    "integration_coefficients 0 cubic %d value", i);
        }

        gsl_vector_set_all(coeffs, -555.0);
        gsl_test(suzerain_bspline_integration_coefficients(
                1 , gsl_vector_ptr(coeffs,0), inc, scratch, w, dw),
                "integration_coefficients 1 cubic");
        for (size_t i = 0; i < w->n; ++i) {
            gsl_test_rel(gsl_vector_get(coeffs, i), e1[i],
                    GSL_DBL_EPSILON*1000,
                    "integration_coefficients 1 cubic %d value", i);
        }

        gsl_vector_set_all(coeffs, -555.0);
        gsl_test(suzerain_bspline_integration_coefficients(
                2 , gsl_vector_ptr(coeffs,0), inc, scratch, w, dw),
                "integration_coefficients 2 cubic");
        for (size_t i = 0; i < w->n; ++i) {
            gsl_test_rel(gsl_vector_get(coeffs, i), e2[i],
                    GSL_DBL_EPSILON*1000,
                    "integration_coefficients 2 cubic %d value", i);
        }

        gsl_vector_set_all(coeffs, -555.0);
        gsl_test(suzerain_bspline_integration_coefficients(
                3 , gsl_vector_ptr(coeffs,0), inc, scratch, w, dw),
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
/////////* 0*/ {4   , 1.25  , 37    , 6.73374  , 3.07899},
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
        gsl_test(suzerain_bspline_htstretch2_evdeltascale_breakpoints(
                        1, k, htdelta, N, &C, &Clow, &Chigh),
                 "%s(%d, %d, %g, %d, ...)", __func__, 1, k, htdelta, N);

        gsl_test(!(Clow <= C && C <= Chigh),
                 "%s:%d self-consistency1 i=%d %g <= %g <= %g",
                 __FILE__, __LINE__, i, Clow, C, Chigh);

        gsl_test(!(Clow <= cc1 && cc1 <= Chigh),
                "%s:%d accuracy1 i=%d %g <= %g <= %g",
                 __FILE__, __LINE__, i, Clow, cc1, Chigh);

        // For C^{(2)} we want...
        gsl_test(suzerain_bspline_htstretch2_evdeltascale_breakpoints(
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
static void test_bspline_htstretch1_evdeltascale_breakpoints()
{
    static const double data[][5] = {
/////////*00*/{ 4, 0.00,  523, 5.44153, 2.72074},
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
        gsl_test(suzerain_bspline_htstretch1_evdeltascale_breakpoints(
                        1, k, htdelta, N, &C, &Clow, &Chigh),
                 "%s(%d, %d, %g, %d, ...)", __func__, 1, k, htdelta, N);

        gsl_test(!(Clow <= C && C <= Chigh),
                 "%s:%d self-consistency1 i=%d %g <= %g <= %g",
                 __FILE__, __LINE__, i, Clow, C, Chigh);

        gsl_test(!(Clow <= cc1 && cc1 <= Chigh),
                "%s:%d accuracy1 i=%d %g <= %g <= %g",
                 __FILE__, __LINE__, i, Clow, cc1, Chigh);

        // For C^{(2)} we want...
        gsl_test(suzerain_bspline_htstretch1_evdeltascale_breakpoints(
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
