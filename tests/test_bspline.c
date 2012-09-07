/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2011, 2012 Rhys Ulerich
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
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_test.h>
#include <suzerain/bspline.h>

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
static void test_bspline_htstretch2_evdeltascale();

int
main(int argc, char **argv)
{
    SUZERAIN_UNUSED(argc);
    SUZERAIN_UNUSED(argv);

    gsl_ieee_env_setup();

    test_integration_coefficients();
    test_linear_combination();
    test_linear_combination_complex();
    test_bspline_htstretch2_evdeltascale();

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

// Test that curve fits reasonably reproduce the source data
// Mainly meant as a sanity check on the coded coefficients
// Test cases chosen pseudo-randomly from the original fit data
static void test_bspline_htstretch2_evdeltascale()
{
    static const double data[14][5] = {
        {  4, 0.25, 197,  9.4995 , 16284.1 },
        {  5, 1.25, 1021, 13.9698, 25.5979 },
        {  6, 2.00, 757,  19.4491, 11.1824 },
        {  7, 1.75, 491,  25.4483, 8.72098 },
        {  8, 0.00, 727,  29.7909, 7.72303 },
        {  9, 2.25, 263,  41.7296, 7.12156 },
        { 10, 2.00, 809,  45.4335, 6.7925  },
        { 11, 1.25, 797,  51.6752, 6.54841 },
        { 12, 2.25, 839,  61.5748, 6.36001 },
        { 13, 1.25, 181,  72.9038, 6.21279 },
        { 14, 1.00, 373,  76.9372, 6.11104 },
        { 15, 1.00,  48,  104.342, 5.98301 },
        { 16, 1.25, 571,  95.1626, 5.94511 },
        { 17, 2.75, 421,  119.228, 5.86762 },
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
        gsl_test(suzerain_bspline_htstretch2_evdeltascale(
                        1, k, htdelta, N, &C, &Clow, &Chigh),
                 "%s(%d, %d, %g, %d, ...)", __func__, 1, k, htdelta, N);

        gsl_test(!(Clow <= C && C <= Chigh),
                 "%s:%d self-consistency1 k=%d %g <= %g <= %g",
                 __FILE__, __LINE__, k, Clow, C, Chigh);

        gsl_test(!(Clow <= cc1 && cc1 <= Chigh),
                "%s:%d accuracy1 k=%d %g <= %g <= %g",
                 __FILE__, __LINE__, k, Clow, cc1, Chigh);

        // For C^{(2)} we want...
        gsl_test(suzerain_bspline_htstretch2_evdeltascale(
                        2, k, htdelta, N, &C, &Clow, &Chigh),
                 "%s(%d, %d, %g, %d, ...)", __func__, 2, k, htdelta, N);

        gsl_test(!(Clow <= C && C <= Chigh),
                 "%s:%d self-consistency2 k=%d %g <= %g <= %g",
                 __FILE__, __LINE__, k, Clow, C, Chigh);

        gsl_test(!(Clow <= cc2 && cc2 <= Chigh),
                 "%s:%d accuracy2 k=%d %g <= %g <= %g",
                 __FILE__, __LINE__, k, Clow, cc2, Chigh);

    }

}
