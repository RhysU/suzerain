#define BOOST_TEST_MODULE $Id$

#include <config.h>

#include <boost/format.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <mkl_blas.h>
#include <log4cxx/logger.h>
#include <suzerain/bspline_operator.h>

using boost::format;
log4cxx::LoggerPtr logger = log4cxx::Logger::getRootLogger();

BOOST_AUTO_TEST_CASE( allocation_okay )
{
    const int order  = 4, nbreak = 10, nderiv = 2;

    suzerain_bspline_operator_workspace * w
        = suzerain_bspline_operator_alloc(order, nbreak, nderiv,
            SUZERAIN_BSPLINE_OPERATOR_COLLOCATION_GREVILLE);

    suzerain_bspline_operator_lu_workspace *luw
        = suzerain_bspline_operator_lu_alloc(w);
    suzerain_bspline_operator_lu_free(luw);

    suzerain_bspline_operator_free(w);
}

// Check a simple piecewise linear case's general banded storage
BOOST_AUTO_TEST_CASE( memory_layout )
{
    const double breakpoints[] = { 0.0, 1.0, 2.0, 3.0 };
    const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);
    const int order  = 2;
    const int nderiv = 1;

    suzerain_bspline_operator_workspace * w
        = suzerain_bspline_operator_alloc(order, nbreak, nderiv,
            SUZERAIN_BSPLINE_OPERATOR_COLLOCATION_GREVILLE);

    suzerain_bspline_operator_create(breakpoints, w);

    /* Check w->D[0], the mass matrix, against known good solution:
     *   1              0              0              0
     *   0              1              0              0
     *   0              0              1              0
     *   0              0              0              1
     * Known good is in general banded matrix column-major order.
     */
    const double good_D0[] = { /*DK*/   1,     0,
                                   0,   1,     0,
                                   0,   1,     0,
                                   0,   1  /*DK*/ };
    BOOST_CHECK_EQUAL_COLLECTIONS(
        good_D0, good_D0 + sizeof(good_D0)/sizeof(good_D0[0]),
        w->D[0] + 1, w->D[0] + w->storagesize - 1);

    /* Check w->D[1], the first derivative matrix, against known good:
     *       -1              1              0              0
     *        0             -1              1              0
     *        0              0             -1              1
     *        0              0             -1              1
     * Known good is in general banded matrix column-major order.
     */
    const double good_D1[] = { /*DK*/  -1,     0,
                                   1,  -1,     0,
                                   1,  -1,    -1,
                                   1,   1  /*DK*/ };
    BOOST_CHECK_EQUAL_COLLECTIONS(
        good_D1, good_D1 + sizeof(good_D1)/sizeof(good_D1[0]),
        w->D[1] + 1, w->D[1] + w->storagesize - 1);

    suzerain_bspline_operator_free(w);
}

// Polynomial test helpers
typedef struct { int n; double c[]; } poly_params; // Flexible array
double poly_f(double x, void *params)
{
    poly_params * p = (poly_params *) params;
    return gsl_poly_eval(p->c, p->n, x);
}
void poly_params_differentiate(poly_params *params)
{
    for (int i = 1; i < params->n; ++i) params->c[i-1] = params->c[i] * i;
    params->c[params->n-1] = 0;
}

// Sanity check the polynomial test helpers
BOOST_AUTO_TEST_CASE( gsl_poly_eval_and_deriv )
{
    poly_params * p = (poly_params *)
                      malloc(sizeof(poly_params) + 3*sizeof(double));
    p->n    = 3;
    p->c[0] = 1.1; // Constant
    p->c[1] = 2.2; // Linear
    p->c[2] = 3.3; // Quadratic
    gsl_function f = {poly_f, p};

    const double  dc[] = { 2.2, 6.6, 0.0 }; // (d/dx)(p->c)
    const double ddc[] = { 6.6, 0.0, 0.0 }; // (d^2/dx^2)(p->c)

    BOOST_REQUIRE_CLOSE(GSL_FN_EVAL(&f,1.0),  6.6, 0.001);
    BOOST_REQUIRE_CLOSE(GSL_FN_EVAL(&f,2.0), 18.7, 0.001);

    poly_params_differentiate(p);
    BOOST_CHECK_EQUAL_COLLECTIONS(p->c, p->c + p->n, dc, dc + p->n);

    poly_params_differentiate(p);
    BOOST_CHECK_EQUAL_COLLECTIONS(p->c, p->c + p->n, ddc, ddc + p->n);

    BOOST_REQUIRE_CLOSE(GSL_FN_EVAL(&f,1.0), 6.6, 0.001);
    BOOST_REQUIRE_CLOSE(GSL_FN_EVAL(&f,2.0), 6.6, 0.001);

    free(p);
}

// BOOST_AUTO_TEST_CASE( differentiate_representable_polynomial )
// {
//     const double breakpoints[] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };
//     const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);
//     const int order  = 3;
//     const int nderiv = 4;
// 
//     suzerain_bspline_operator_workspace * w
//         = suzerain_bspline_operator_alloc(order, nbreak, nderiv,
//             SUZERAIN_BSPLINE_OPERATOR_COLLOCATION_GREVILLE);
// 
//     suzerain_bspline_operator_create(breakpoints, w);
// 
// 
// 
//     /* Check w->D[0], the mass matrix, against known good solution:
//      *   1              0              0              0
//      *   0              1              0              0
//      *   0              0              1              0
//      *   0              0              0              1
//      * Known good is in general banded matrix column-major order.
//      */
//     const double good_D0[] = { /*DK*/   1,     0,
//                                    0,   1,     0,
//                                    0,   1,     0,
//                                    0,   1  /*DK*/ };
//     BOOST_CHECK_EQUAL_COLLECTIONS(
//         good_D0, good_D0 + sizeof(good_D0)/sizeof(good_D0[0]),
//         w->D[0] + 1, w->D[0] + w->storagesize - 1);
// 
//     /* Check w->D[1], the first derivative matrix, against known good:
//      *       -1              1              0              0
//      *        0             -1              1              0
//      *        0              0             -1              1
//      *        0              0             -1              1
//      * Known good is in general banded matrix column-major order.
//      */
//     const double good_D1[] = { /*DK*/  -1,     0,
//                                    1,  -1,     0,
//                                    1,  -1,    -1,
//                                    1,   1  /*DK*/ };
//     BOOST_CHECK_EQUAL_COLLECTIONS(
//         good_D1, good_D1 + sizeof(good_D1)/sizeof(good_D1[0]),
//         w->D[1] + 1, w->D[1] + w->storagesize - 1);
// 
//     suzerain_bspline_operator_free(w);
// }
