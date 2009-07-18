#define BOOST_TEST_MODULE $Id$

#include <suzerain/config.h>

#include <string.h>
#include <boost/format.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <gsl/gsl_poly.h>
#include <suzerain/blas_et_al.h>
#include <suzerain/bspline.h>
#include <suzerain/error.h>
#include <suzerain/function.h>
#include <log4cxx/logger.h>

log4cxx::LoggerPtr logger = log4cxx::Logger::getRootLogger();

#include "test_tools.hpp"

// TODO Check suzerain_bspline_lu_form_general for nonuniform ku/kl
// TODO compute_derivatives_of_a_general_polynomial for more orders

BOOST_AUTO_TEST_CASE( allocation_okay )
{
    const double breakpoints[] = { 0.0, 1.0, 2.0, 3.0 };
    const int order  = 4;
    const int nderiv = 2;
    const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);

    suzerain_bspline_workspace *w
        = suzerain_bspline_alloc(order, nderiv, nbreak, breakpoints,
            SUZERAIN_BSPLINE_COLLOCATION_GREVILLE);

    suzerain_bspline_lu_workspace *luw
        = suzerain_bspline_lu_alloc(w);
    suzerain_bspline_lu_free(luw);

    suzerain_bspline_free(w);
}

// Check a simple piecewise linear case's general banded storage
BOOST_AUTO_TEST_CASE( piecewise_linear_memory_application_solution )
{
    const double breakpoints[] = { 0.0, 1.0, 2.0, 3.0 };
    const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);
    const int order  = 2;
    const int nderiv = 2;

    suzerain_bspline_workspace *w
        = suzerain_bspline_alloc(order, nderiv, nbreak, breakpoints,
            SUZERAIN_BSPLINE_COLLOCATION_GREVILLE);

    {
        /* Check w->D[0], the mass matrix, against known good solution:
         *   1   0   0   0
         *   0   1   0   0
         *   0   0   1   0
         *   0   0   0   1
         * Known good is in general banded matrix column-major order.
         */
        const double good_D0[] = { 1,
                                   1,
                                   1,
                                   1 };
        BOOST_CHECK_EQUAL(0, w->ku[0]);
        BOOST_CHECK_EQUAL(0, w->kl[0]);
        BOOST_CHECK_EQUAL_COLLECTIONS(
            good_D0, good_D0 + sizeof(good_D0)/sizeof(good_D0[0]),
            w->D[0] + w->ku[0], w->D[0] + w->storagesize[0] - w->kl[0]);

        /* Check w->D[0] application against multiple vectors */
        const int nrhs = 2;
        double vector[] = { 1, 2, 3, 4,
                            5, 6, 7, 8 };
        const double good_result[] = { 1, 2, 3, 4,
                                       5, 6, 7, 8 };
        const int ldb = sizeof(vector)/(sizeof(vector[0]))/nrhs;
        suzerain_bspline_apply_operator(0, nrhs, vector, ldb, w);
        BOOST_CHECK_EQUAL_COLLECTIONS(
            good_result, good_result + sizeof(good_result)/sizeof(good_result[0]),
            vector, vector + sizeof(vector)/sizeof(vector[0]));
    }

    {
        /* Check w->D[1], the first derivative matrix, against known good:
         *  -1   1   0   0
         *   0  -1   1   0
         *   0   0  -1   1
         *   0   0  -1   1
         * Known good is in general banded matrix column-major order.
         */
        const double good_D1[] = { /*DK*/  -1,     0,
                                       1,  -1,     0,
                                       1,  -1,    -1,
                                       1,   1  /*DK*/ };
        BOOST_CHECK_EQUAL(1, w->ku[1]);
        BOOST_CHECK_EQUAL(1, w->kl[1]);
        BOOST_CHECK_EQUAL_COLLECTIONS(
            good_D1, good_D1 + sizeof(good_D1)/sizeof(good_D1[0]),
            w->D[1] + w->ku[1], w->D[1] + w->storagesize[1] - w->kl[1]);

        /* Check w->D[0] application against multiple vectors */
        const int nrhs = 2;
        double vector[] = { 1, 3, 2, 4,
                            7, 6, 5, 8 };
        const double good_result[] = {  2, -1, 2, 2,
                                       -1, -1, 3, 3 };
        const int ldb = sizeof(vector)/(sizeof(vector[0]))/nrhs;
        suzerain_bspline_apply_operator(1, nrhs, vector, ldb, w);
        BOOST_CHECK_EQUAL_COLLECTIONS(
            good_result, good_result + sizeof(good_result)/sizeof(good_result[0]),
            vector, vector + sizeof(vector)/sizeof(vector[0]));
    }

    {
        /* Check w->D[2], the second derivative matrix, against zero result.
         */
        const double good_D2[] = { 0,
                                   0,
                                   0,
                                   0 };
        BOOST_CHECK_EQUAL(0, w->ku[2]);
        BOOST_CHECK_EQUAL(0, w->kl[2]);
        BOOST_CHECK_EQUAL_COLLECTIONS(
            good_D2, good_D2 + sizeof(good_D2)/sizeof(good_D2[0]),
            w->D[2] + w->ku[2], w->D[2] + w->storagesize[2] - w->kl[2]);
    }

    suzerain_bspline_lu_workspace *luw
        = suzerain_bspline_lu_alloc(w);

    /* Form 2*D[0] - 3*D[1] operator in LU-ready banded storage.  Answer is
     *   5   -3    0     0
     *   0    5   -3     0
     *   0    0    5    -3
     *   0    0    3    -1
     * which, in LU-form where L has ones on the main diagonal, is
     *   5   -3    0     0
     *   0    5   -3     0
     *   0    0    5    -3
     *   0    0    0.6   0.8
     * The pivot matrix is eye(4).  Check it in octave using [l,u,p] = lu(A).
     * Known good is in general banded matrix column-major order with
     * additional superdiagonal to allow for LU factorization fill-in.
     */
    const double good_A0[] = { /*DK*/   /*DK*/   5,     0,
                               /*DK*/0,    -3,   5,     0,
                                     0,    -3,   5,     0.6,
                                     0,    -3,   0.8  /*DK*/ };
    const double coeff[] = { 2.0, -3.0 };
    suzerain_bspline_lu_form_general(
        sizeof(coeff)/sizeof(coeff[0]), coeff, w, luw);
    check_close_collections(
        good_A0, good_A0 + sizeof(good_A0)/sizeof(good_A0[0]),
        luw->A + luw->ku, luw->A + luw->storagesize - luw->kl,
        1e-12);

    /* Check that multiple rhs solution works for operator found just above */
    {
        const int nrhs = 2;
        double vector[] = { 1,  2, 3, 4,
                           -4, -1, 1, 3};
        const double good_result[] = {  1.25, 1.75, 2.25, 2.75,
                                       -0.2,  1,    2,    3 };
        const int ldb = sizeof(vector)/(sizeof(vector[0]))/nrhs;
        suzerain_bspline_lu_solve(nrhs, vector, ldb, luw);
        check_close_collections(
            good_result, good_result + sizeof(good_result)/sizeof(good_result[0]),
            vector, vector + sizeof(vector)/sizeof(vector[0]),
            1.0e-12);
    }

    suzerain_bspline_lu_free(luw);
    suzerain_bspline_free(w);
}

// Check a simple piecewise quadratic case's general banded storage
BOOST_AUTO_TEST_CASE( piecewise_quadratic_memory_application_solution )
{
    const double breakpoints[] = { 0.0, 1.0, 2.0, 3.0 };
    const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);
    const int order  = 3;
    const int nderiv = 1;

    suzerain_bspline_workspace *w
        = suzerain_bspline_alloc(order, nderiv, nbreak, breakpoints,
            SUZERAIN_BSPLINE_COLLOCATION_GREVILLE);

    {
        /* Check w->D[0], the mass matrix, against known good solution:
         *  1.    0.    0.    0.    0.
         *  1./4. 5./8. 1./8. 0.    0.
         *  0.    1./8. 3./4. 1./8. 0.
         *  0.    0.    1./8. 5./8. 1./4.
         *  0.    0.    0.    0.    1.
         * Known good is in general banded matrix column-major order.
         */
        const double good_D0[] = { /*DK*/      1.,       1./4.,
                                      0.,      5./8.,    1./8.,
                                      1./8.,   3./4.,    1./8.,
                                      1./8.,   5./8.,    0.,
                                      1./4.,       1.     /*DK*/ };
        BOOST_CHECK_EQUAL_COLLECTIONS(
            good_D0, good_D0 + sizeof(good_D0)/sizeof(good_D0[0]),
            w->D[0] + w->ku[0], w->D[0] + w->storagesize[0] - w->kl[0]);

        /* Check w->D[0] application against multiple vectors */
        const int nrhs = 2;
        double vector[] = { 1, 2, 3, 4, 5,
                            5, 6, 7, 8, 9 };
        const double good_result[] = { 1., 15./8., 3., 33./8., 5.,
                                       5., 47./8., 7., 65./8., 9. };
        const int ldb = sizeof(vector)/(sizeof(vector[0]))/nrhs;
        suzerain_bspline_apply_operator(0, nrhs, vector, ldb, w);
        BOOST_CHECK_EQUAL_COLLECTIONS(
            good_result, good_result + sizeof(good_result)/sizeof(good_result[0]),
            vector, vector + sizeof(vector)/sizeof(vector[0]));
    }

    suzerain_bspline_free(w);
}


// Check a piecewise cubic case's general banded storage
BOOST_AUTO_TEST_CASE( piecewise_cubic_memory_application_solution )
{
    const double breakpoints[] = { 0.0, 1.0, 2.0, 3.0 };
    const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);
    const int order  = 4;
    const int nderiv = 2;

    suzerain_bspline_workspace *w
        = suzerain_bspline_alloc(order, nderiv, nbreak, breakpoints,
            SUZERAIN_BSPLINE_COLLOCATION_GREVILLE);

    {
        /* Check w->D[0], the mass matrix, against known good solution:
         * 1.       0.        0.       0.        0.       0.
         * 8./27.  61./108.  43./324.  1./162.   0.       0.
         * 0.       1./4.     7./12.   1./6.     0.       0.
         * 0.       0.        1./6.    7./12.    1./4.    0.
         * 0.       0.        1./162. 43./324.  61./108.  8./27.
         * 0.       0.        0.       0.        0.       1.
         * Known good is in general banded matrix column-major order.
         */
        const double good_D0[] = {
            /*DK*/      /*DK*/          1.,       8./27.,         0,
            /*DK*/0,        0,         61./108.,  1./4.,          0,
                  0,       43./324.,    7./12.,   1./6.,          1./162.,
                  1./162.,  1./6.,      7./12.,   43./324.,       0,
                  0,        1./4.,     61./108.,  0,        /*DK*/0,
                  0,        8./27.,     1.        /*DK*/    /*DK*/
        };
        check_close_collections(
            good_D0, good_D0 + sizeof(good_D0)/sizeof(good_D0[0]),
            w->D[0] + w->ku[0], w->D[0] + w->storagesize[0] - w->kl[0],
            1.0e-12);

        /* Check w->D[0] application against multiple vectors */
        {
            const int nrhs = 2;
            double vector[] = { 1, 2, 3, 4, 5, 6,
                                4, 5, 6, 7, 8, 9  };
            const double good_result[] = {
                1., 599./324., 35./12., 49./12., 1669./324., 6.,
                4., 1571./324., 71./12., 85./12., 2641./324., 9.
            };
            const int ldb = sizeof(vector)/(sizeof(vector[0]))/nrhs;
            suzerain_bspline_apply_operator(0, nrhs, vector, ldb, w);
            check_close_collections(
                good_result, good_result + sizeof(good_result)/sizeof(good_result[0]),
                vector, vector + sizeof(vector)/sizeof(vector[0]),
                1.0e-12);
        }
    }

    suzerain_bspline_free(w);
}

// Polynomial test helpers
typedef struct { int n; double c[]; } poly_params; // Flexible array
double poly_f(double x, void *params)
{
    poly_params *p = (poly_params *) params;
    return gsl_poly_eval(p->c, p->n, x);
}
void poly_params_differentiate(poly_params *params)
{
    for (int i = 1; i < params->n; ++i) params->c[i-1] = params->c[i] *i;
    params->c[params->n-1] = 0;
}

// Sanity check the polynomial test helpers
BOOST_AUTO_TEST_CASE( gsl_poly_eval_and_deriv )
{
    poly_params *p = (poly_params *)
                      malloc(sizeof(poly_params) + 3*sizeof(double));
    p->n    = 3;
    p->c[0] = 1.1; // Constant
    p->c[1] = 2.2; // Linear
    p->c[2] = 3.3; // Quadratic
    suzerain_function f = {poly_f, p};

    const double   dc[] = { 2.2, 6.6, 0.0 }; // (d/dx)(p->c)
    const double  ddc[] = { 6.6, 0.0, 0.0 }; // (d^2/dx^2)(p->c)
    const double dddc[] = { 0.0, 0.0, 0.0 }; // (d^3/dx^3)(p->c)

    BOOST_CHECK_CLOSE(SUZERAIN_FN_EVAL(&f,1.0),  6.6, 0.001);
    BOOST_CHECK_CLOSE(SUZERAIN_FN_EVAL(&f,2.0), 18.7, 0.001);

    poly_params_differentiate(p);
    BOOST_CHECK_EQUAL_COLLECTIONS(p->c, p->c + p->n, dc, dc + p->n);

    poly_params_differentiate(p);
    BOOST_CHECK_EQUAL_COLLECTIONS(p->c, p->c + p->n, ddc, ddc + p->n);

    BOOST_CHECK_CLOSE(SUZERAIN_FN_EVAL(&f,1.0), 6.6, 0.001);
    BOOST_CHECK_CLOSE(SUZERAIN_FN_EVAL(&f,2.0), 6.6, 0.001);

    // Differentiate a constant
    poly_params_differentiate(p);
    BOOST_CHECK_EQUAL_COLLECTIONS(p->c, p->c + p->n, dddc, dddc + p->n);

    BOOST_CHECK_EQUAL(SUZERAIN_FN_EVAL(&f,1.0), 0.0);
    BOOST_CHECK_EQUAL(SUZERAIN_FN_EVAL(&f,2.0), 0.0);

    // Differentiate zero
    poly_params_differentiate(p);
    BOOST_CHECK_EQUAL_COLLECTIONS(p->c, p->c + p->n, dddc, dddc + p->n);

    BOOST_CHECK_EQUAL(SUZERAIN_FN_EVAL(&f,1.0), 0.0);
    BOOST_CHECK_EQUAL(SUZERAIN_FN_EVAL(&f,2.0), 0.0);

    free(p);
}

BOOST_AUTO_TEST_CASE( compute_derivatives_of_a_general_polynomial )
{
    // Test parameters
    const double breakpoints[] = { 0.0, 1.0, 2.0, 3.0 };
    const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);
    const int order  = 7; /* Comparatively higher order than above tests */
    const int nderiv = 7;

    // Initialize workspaces
    suzerain_bspline_workspace *w
        = suzerain_bspline_alloc(order, nderiv, nbreak, breakpoints,
            SUZERAIN_BSPLINE_COLLOCATION_GREVILLE);
    const int ncoeff = suzerain_bspline_ncoefficients(w);

    // Initialize mass matrix in factored form
    suzerain_bspline_lu_workspace *mass
        = suzerain_bspline_lu_alloc(w);
    suzerain_bspline_lu_form_mass(w, mass);

    // Initialize test function
    poly_params *p = (poly_params *)
                      malloc(sizeof(poly_params) + order*sizeof(double));
    p->n = order;
    p->c[0] = 1.9;
    for (int i = 1; i < order; ++i) {
        p->c[i] = p->c[i-1] + 0.9;
    }
    suzerain_function f = {poly_f, p};

    // Compute expected coefficients for derivatives [0...nderiv]
    // by directly differentiating the polynomial test function.
    double * const expected
        = (double *) malloc((nderiv+1) * ncoeff * sizeof(double));
    for (int i = 0; i <= nderiv; ++i) {
        suzerain_bspline_find_coefficient_rhs(
                &f, expected + i*ncoeff, w);
        poly_params_differentiate(p);  // Drop the polynomial order by one
    }
    suzerain_bspline_lu_solve(nderiv+1, expected, ncoeff, mass);

    // Make copies of the zeroth derivative coefficients
    double * const actual
        = (double *) malloc((nderiv+1) * ncoeff * sizeof(double));
    for (int i = 0; i <= nderiv; ++i) {
        memcpy(actual + i*ncoeff, expected, ncoeff * sizeof(actual[0]));
    }

    // Solve M*x' = D*x ...
    // ...starting by applying the derivative operators
    for (int i = 0; i <= nderiv; ++i) {
        suzerain_bspline_apply_operator(i, 1, actual + i*ncoeff, ncoeff, w);
    }
    // ...finish by solving with the mass matrix
    suzerain_bspline_lu_solve(nderiv+1, actual, ncoeff, mass);

    // See if we got anywhere close
    for (int i = 0; i <= nderiv; ++i) {
        check_close_collections(
                expected + i*ncoeff, expected + (i+1)*ncoeff,
                actual + i*ncoeff, actual + (i+1)*ncoeff, 1.0e-09);
    }

    free(actual);
    free(expected);
    suzerain_bspline_lu_free(mass);
    suzerain_bspline_free(w);
    free(p);
}

BOOST_AUTO_TEST_CASE( derivatives_of_a_piecewise_cubic_representation )
{
    const double breakpoints[] = { 0.0, 1.0, 2.0, 3.0 };
    const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);
    const int order  = 4;
    const int nderiv = 3;

    poly_params *p = (poly_params *)
                      malloc(sizeof(poly_params) + 4*sizeof(double));
    p->n = 4;
    suzerain_function f = {poly_f, p};

    suzerain_bspline_workspace *w
        = suzerain_bspline_alloc(order, nderiv, nbreak, breakpoints,
            SUZERAIN_BSPLINE_COLLOCATION_GREVILLE);
    const int ncoeff = suzerain_bspline_ncoefficients(w);

    // Form the mass matrix M
    suzerain_bspline_lu_workspace * const mass
        = suzerain_bspline_lu_alloc(w);
    suzerain_bspline_lu_form_mass(w, mass);

    {
        const int derivative = 1;

        p->c[0] = 1.2; // Constant
        p->c[1] = 3.4; // Linear
        p->c[2] = 0.0; // Quadratic
        p->c[3] = 0.0; // Cubic

        // Compute the right hand side coefficients for M x = b
        double * coefficient = (double *) malloc(ncoeff * sizeof(double));
        suzerain_bspline_find_coefficient_rhs(&f, coefficient, w);

        // Solve for function coefficients using the mass matrix
        suzerain_bspline_lu_solve(1, coefficient, ncoeff, mass);

        // Take the n-th derivative of the coefficients using M x' = D x
        suzerain_bspline_apply_operator(derivative, 1, coefficient, ncoeff, w);
        suzerain_bspline_lu_solve(1, coefficient, ncoeff, mass);

        // Ensure we recover the leading order, scaled monomial coefficients
        for (int i = 0; i < ncoeff; ++i) {
            BOOST_CHECK_CLOSE(1.0 * p->c[1], coefficient[i], 1e-12);
        }

        free(coefficient);
    }

    {
        const int derivative = 2;

        p->c[0] = 1.2; // Constant
        p->c[1] = 3.4; // Linear
        p->c[2] = 5.6; // Quadratic
        p->c[3] = 0.0; // Cubic

        // Compute the right hand side coefficients for M x = b
        double * coefficient = (double *) malloc(ncoeff * sizeof(double));
        suzerain_bspline_find_coefficient_rhs(&f, coefficient, w);

        // Solve for function coefficients using the mass matrix
        suzerain_bspline_lu_solve(1, coefficient, ncoeff, mass);

        // Take the n-th derivative of the coefficients using M x' = D x
        suzerain_bspline_apply_operator(derivative, 1, coefficient, ncoeff, w);
        suzerain_bspline_lu_solve(1, coefficient, ncoeff, mass);

        // Ensure we recover the leading order, scaled monomial coefficients
        for (int i = 0; i < ncoeff; ++i) {
            BOOST_CHECK_CLOSE(2.0 * p->c[2], coefficient[i], 1e-11);
        }

        free(coefficient);
    }

    {
        const int derivative = 3;

        p->c[0] = 1.2; // Constant
        p->c[1] = 3.4; // Linear
        p->c[2] = 5.6; // Quadratic
        p->c[3] = 7.8; // Cubic

        // Compute the right hand side coefficients for M x = b
        double * coefficient = (double *) malloc(ncoeff * sizeof(double));
        suzerain_bspline_find_coefficient_rhs(&f, coefficient, w);

        // Solve for function coefficients using the mass matrix
        suzerain_bspline_lu_solve(1, coefficient, ncoeff, mass);

        // Take the n-th derivative of the coefficients using M x' = D x
        suzerain_bspline_apply_operator(derivative, 1, coefficient, ncoeff, w);
        suzerain_bspline_lu_solve(1, coefficient, ncoeff, mass);

        // Ensure we recover the leading order, scaled monomial coefficients
        for (int i = 0; i < ncoeff; ++i) {
            BOOST_CHECK_CLOSE(6.0 * p->c[3], coefficient[i], 1e-11);
        }

        free(coefficient);
    }

    suzerain_bspline_lu_free(mass);
    suzerain_bspline_free(w);
    free(p);
}

void
log4cxx_error_handler(const char *reason, const char *file,
                      int line, int err)
{
    LOG4CXX_ERROR(logger,
                  boost::format("%s caught [%s:%d: %s (%d)]")
                  % __func__ % file % line % reason % err);
}

// Intended to ensure our bandwidth routines are okay for
// high order bsplines and high order derivatives
BOOST_AUTO_TEST_CASE( ensure_create_operation_in_alloc_succeeds )
{
    const double breakpoints[] = { 0.0, 1.0, 2.0, 3.0 };
    const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);

    suzerain_error_handler_t * previous_handler
        = suzerain_set_error_handler(&log4cxx_error_handler);

    const int maxorder = 21;
    for (int order = 1; order <= maxorder; ++order) {
        const int maxnderiv = maxorder;
        for (int nderiv = 0; nderiv <= maxnderiv; ++nderiv) {
            suzerain_bspline_workspace *w
                = suzerain_bspline_alloc(
                    order, nderiv, nbreak, breakpoints,
                    SUZERAIN_BSPLINE_COLLOCATION_GREVILLE);
            BOOST_CHECK_MESSAGE(w != NULL, boost::format(
                "Error allocating operator for order %d, nderiv %d")
                % order % nderiv);

            if (order == nderiv && w != NULL) {
                for (int k = 0; k <= w->nderivatives; ++k) {
                LOG4CXX_TRACE(logger, boost::format(
                    "Bandwidth details: order=%2d, deriv=%2d, kl=%2d, ku=%2d")
                    % order % k % w->kl[k] % w->ku[k]);
                }
            }

            suzerain_bspline_free(w);  // Should accept w == NULL
        }
    }

    suzerain_set_error_handler(previous_handler);
}

BOOST_AUTO_TEST_CASE( bspline_evaluation_routine )
{
    const double breakpoints[] = { 0.0, 1.0, 2.0, 3.0 };
    const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);
    const int order  = 4;
    const int nderiv = 2;
    const int ncoeff = 6;

    suzerain_bspline_workspace *w
        = suzerain_bspline_alloc(order, nderiv, nbreak, breakpoints,
            SUZERAIN_BSPLINE_COLLOCATION_GREVILLE);

    // Sanity check on fixed storage assumption for test case
    BOOST_REQUIRE_EQUAL(ncoeff, suzerain_bspline_ncoefficients(w));

    // Check that we have a partition of unity
    {
        double coefficients[ncoeff];
        for (int i = 0; i < ncoeff; ++i) coefficients[i] = 1.0;

        const double * points = breakpoints;
        const int npoints = sizeof(breakpoints)/sizeof(breakpoints[0]);

        const int ldvalues = npoints;
        const int nvalues = (nderiv+1) * ldvalues;
        double values[nvalues];

        suzerain_bspline_evaluate(
                nderiv, coefficients, npoints, points, values, ldvalues, w);

        for (int i=0; i < nvalues; ++i) {
// TODO FIXME
//            BOOST_CHECK_EQUAL(values[i], 1.0);
        }
    }

    suzerain_bspline_free(w);
}
