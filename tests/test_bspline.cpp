#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_poly.h>
#include <suzerain/blas_et_al.h>
#include <suzerain/bspline.h>
#include <suzerain/error.h>
#include <suzerain/function.h>
#include <suzerain/math.hpp>
#include <log4cxx/logger.h>

// We test C API but include C++ API to ensure it compiles
// TODO Exercise C++ API, including bspline, bspline_lu, bspline_luz
#include <suzerain/bspline.hpp>

log4cxx::LoggerPtr logger = log4cxx::Logger::getRootLogger();

#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

#pragma warning(disable:1572 2014 2015)

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
    BOOST_REQUIRE(w != NULL);

    suzerain_bspline_lu_workspace *luw
        = suzerain_bspline_lu_alloc(w);
    BOOST_REQUIRE(luw != NULL);

    suzerain_bspline_luz_workspace *luzw
        = suzerain_bspline_luz_alloc(w);
    BOOST_REQUIRE(luzw != NULL);

    suzerain_bspline_luz_free(luzw);
    suzerain_bspline_lu_free(luw);
    suzerain_bspline_free(w);
}

// Check a simple piecewise linear case's general banded storage
BOOST_AUTO_TEST_CASE( collocation_piecewise_linear_memory_application_soln )
{
    const double breakpoints[] = { 0.0, 1.0, 2.0, 3.0 };
    const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);
    const int order  = 2;
    const int nderiv = 2;

    suzerain_bspline_workspace *w
        = suzerain_bspline_alloc(order, nderiv, nbreak, breakpoints,
            SUZERAIN_BSPLINE_COLLOCATION_GREVILLE);
    BOOST_CHECK_EQUAL(suzerain_bspline_ndof(w), w->ndof);

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
        CHECK_GBMATRIX_CLOSE(
                  4,       4,        0,        0, good_D0,     1,
            w->ndof, w->ndof, w->kl[0], w->ku[0], w->D[0], w->ld,
            1e-12);

        {
            /* Check w->D[0] application against multiple real vectors */
            const int nrhs = 2;
            double b[] = { 1, 2, 3, 4,
                           5, 6, 7, 8 };
            const int ldb = sizeof(b)/(sizeof(b[0]))/nrhs;
            const int incb = 1;
            suzerain_bspline_apply_operator(0, nrhs, 3.0, b, incb, ldb, w);
            const double b_good[] = { 3*1, 3*2, 3*3, 3*4,
                                      3*5, 3*6, 3*7, 3*8 };
            BOOST_CHECK_EQUAL_COLLECTIONS(
                b_good, b_good + sizeof(b_good)/sizeof(b_good[0]),
                b, b + sizeof(b)/sizeof(b[0]));
        }

        {
            /* Check w->D[0] application against multiple complex vectors */
            const int nrhs = 2;
            double b[][2] = {
                { 1,  2}, { 3,  4}, { 5,  6}, { 7,  8},
                {11, 12}, {13, 14}, {15, 16}, {17, 18}
            };
            const int ldb = sizeof(b)/(sizeof(b[0]))/nrhs;
            const int incb = 1;
            suzerain_bspline_zapply_operator(0, nrhs, 1.0, b, incb, ldb, w);
            const double b_good[][2] = {
                { 1,  2}, { 3,  4}, { 5,  6}, { 7,  8},
                {11, 12}, {13, 14}, {15, 16}, {17, 18}
            };
            BOOST_CHECK_EQUAL_COLLECTIONS(
                (double *) b_good,
                (double *)(b_good + sizeof(b_good)/sizeof(b_good[0])),
                (double *) b,
                (double *)(b + sizeof(b)/sizeof(b[0])));
        }

        {
            /* Check w->D[0] accumulation against multiple real vectors */
            const int nrhs = 2;
            const double x[] = { 1, 2, 3, 4,
                                 5, 6, 7, 8 };
            const int ldx  = sizeof(x)/sizeof(x[0])/nrhs;
            const int incx = 1;
            double y[] = {  1,  1,  1,  1,
                           -1, -1, -1, -1 };
            const int ldy  = sizeof(y)/sizeof(y[0])/nrhs;
            const int incy = 1;
            suzerain_bspline_accumulate_operator(
                0, nrhs, 1.0, x, incx, ldx, 1.0, y, incy, ldy, w);
            const double y_good[] = { 2, 3, 4, 5,
                                      4, 5, 6, 7 };
            BOOST_CHECK_EQUAL_COLLECTIONS(
                y_good,
                y_good + sizeof(y_good)/sizeof(y_good[0]),
                y,
                y + sizeof(y)/sizeof(y[0]));
        }

        {
            /* Check w->D[0] accumulation against multiple complex vectors */
            /* using complex-valued coefficients                           */
            const int nrhs = 2;
            const double x[][2] = {
                { 1,  2}, { 3,  4}, { 5,  6}, { 7,  8},
                {11, 12}, {13, 14}, {15, 16}, {17, 18}
            };
            const int ldx  = sizeof(x)/sizeof(x[0])/nrhs;
            const int incx = 1;
            double y[][2] = {
                { 1,  1}, { 1,  1}, {-1, -1}, {-1, -1},
                {-1, -1}, {-1, -1}, { 1,  1}, { 1,  1}
            };
            const int ldy  = sizeof(y)/sizeof(y[0])/nrhs;
            const int incy = 1;
            const double alpha[2] = { 2, 3 }; /* NB complex-valued */
            const double beta[2]  = { 5, 7 }; /* NB complex-valued */
            suzerain_bspline_zaccumulate_operator(
                0, nrhs, alpha, x, incx, ldx, beta, y, incy, ldy, w);
            const double y_good[][2] = {
                 { -6, 19},  { -8, 29},  { -6, 15},  { -8, 25},
                 {-12, 45},  {-14, 55},  {-20, 89},  {-22, 99}
            };
            BOOST_CHECK_EQUAL_COLLECTIONS(
                (double *) y_good,
                (double *)(y_good + sizeof(y_good)/sizeof(y_good[0])),
                (double *) y,
                (double *)(y + sizeof(y)/sizeof(y[0])));
        }

        {
            /* Check w->D[0] accumulation against multiple complex vectors */
            /* using real-valued coefficients                              */
            const int nrhs = 2;
            const double x[][2] = {
                { 1,  2}, { 3,  4}, { 5,  6}, { 7,  8},
                {11, 12}, {13, 14}, {15, 16}, {17, 18}
            };
            const int ldx  = sizeof(x)/sizeof(x[0])/nrhs;
            const int incx = 1;
            double y[][2] = {
                { 1,  1}, { 1,  1}, {-1, -1}, {-1, -1},
                {-1, -1}, {-1, -1}, { 1,  1}, { 1,  1}
            };
            const int ldy  = sizeof(y)/sizeof(y[0])/nrhs;
            const int incy = 1;
            const double alpha[2] = { 2, 0 }; /* NB real-valued */
            const double beta[2]  = { 5, 0 }; /* NB real-valued */
            suzerain_bspline_zaccumulate_operator(
                0, nrhs, alpha, x, incx, ldx, beta, y, incy, ldy, w);
            const double y_good[][2] = {
                { 7,  9}, {11, 13}, { 5,  7}, { 9, 11},
                {17, 19}, {21, 23}, {35, 37}, {39, 41}
            };
            BOOST_CHECK_EQUAL_COLLECTIONS(
                (double *) y_good,
                (double *)(y_good + sizeof(y_good)/sizeof(y_good[0])),
                (double *) y,
                (double *)(y + sizeof(y)/sizeof(y[0])));
        }
    }

    {
        /* Check w->D[1], the first derivative matrix, against known good:
         *  -1   1   0   0
         *   0  -1   1   0
         *   0   0  -1   1
         *   0   0  -1   1
         * Known good is in general banded matrix column-major order.
         */
        const double good_D1[] = { /*DK*/0,   -1,       0,
                                         1,   -1,       0,
                                         1,   -1,      -1,
                                         1,    1, /*DK*/0 };
        CHECK_GBMATRIX_CLOSE(
                  4,       4,        1,        1, good_D1,     3,
            w->ndof, w->ndof, w->kl[1], w->ku[1], w->D[1], w->ld,
            1e-12);

        /* Check w->D[1] application against multiple vectors */
        /* Includes b having non-unit stride */
        const int nrhs = 2;
        double vapply[] = {
            1, /*DK*/555, 3, /*DK*/555, 2, /*DK*/555, 4,
            7, /*DK*/555, 6, /*DK*/555, 5, /*DK*/555, 8
        };
        const double vapply_good[] = {
            2, /*DK*/555, -1, /*DK*/555, 2, /*DK*/555, 2,
           -1, /*DK*/555, -1, /*DK*/555, 3, /*DK*/555, 3
        };
        const int ldb = sizeof(vapply)/(sizeof(vapply[0]))/nrhs;
        const int incb = 2;
        suzerain_bspline_apply_operator(1, nrhs, 1.0, vapply, incb, ldb, w);
        BOOST_CHECK_EQUAL_COLLECTIONS(
            vapply_good,
            vapply_good + sizeof(vapply_good)/sizeof(vapply_good[0]),
            vapply,
            vapply + sizeof(vapply)/sizeof(vapply[0]));
    }

    {
        /* Check w->D[2], the second derivative matrix, against zero result.
         */
        const double good_D2[] = { 0,
                                   0,
                                   0,
                                   0 };
        CHECK_GBMATRIX_CLOSE(
                  4,       4,        0,        0, good_D2,     1,
            w->ndof, w->ndof, w->kl[2], w->ku[2], w->D[2], w->ld,
            1e-12);
    }

    /************************************************************************
     * Integration functionality
     * Exact integration solutions found using Mathematica 7 as follows:
     *  k = 2;
     *  b = {0, 1, 2, 3};
     *  n = Length[b] + k - 2;
     *  knots = Join[ ConstantArray[First[b], k - 1],
     *                b, ConstantArray[Last[b], k - 1]];
     *  B[i_, x_] := BSplineBasis[{k - 1, knots}, i, x]
     *  Table[Integrate[B[i, x], {x, Min[b], Max[b]}], {i, 0, n - 1}]
     ***********************************************************************/
    {
        const double expected[] = { .5, 1, 1, .5 };
        check_close_collections(
            expected, expected + sizeof(expected)/sizeof(expected[0]),
            w->I, w->I + w->ndof, 1e-12);

        const double coeffs[] = { 1, 777, 888,
                                  2, 777, 888,
                                  3, 777, 888,
                                  4, 777, 888 }; // Padded, incx = 3
        double value;
        suzerain_bspline_integrate(coeffs, /*incx*/3, &value, w);
        BOOST_CHECK_CLOSE(value, 15./2, 1e-12);

        const double zcoeffs[][2] = { {1,5}, {777,888},
                                      {2,6}, {777,888},
                                      {3,7}, {777,888},
                                      {4,8}, {777,888} }; // Padded, incx = 2
        double zvalue[2];
        suzerain_bspline_zintegrate(zcoeffs, /*incx*/2, &zvalue, w);
        BOOST_CHECK_CLOSE(zvalue[0], 15./2, 1e-12);
        BOOST_CHECK_CLOSE(zvalue[1], 39./2, 1e-12);
    }

    /********************************/
    /* Real-valued LU functionality */
    /********************************/

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
    {
        const double good_A0[] = { /*DK*/0, /*DK*/0,   5,         0,
                                   /*DK*/0,      -3,   5,         0,
                                         0,      -3,   5,         0.6,
                                         0,      -3,   0.8, /*DK*/0    };
        const double coeff[] = { 2.0, -3.0, 999.0 };
        suzerain_bspline_lu_form_general(
            sizeof(coeff)/sizeof(coeff[0]), coeff, w, luw);
        CHECK_GBMATRIX_CLOSE(
                    4,         4,       1,       2, good_A0,       4,
            luw->ndof, luw->ndof, luw->kl, luw->ku,  luw->A, luw->ld,
            1e-12);
    }

    /* Check that multiple rhs solution works for operator found just above */
    /* Also ensures that non-unit stride works as expected */
    {
        const int nrhs = 2;
        double b[] = {
             1, /*DK*/555,  2, /*DK*/555, 3, /*DK*/555, 4,
            -4, /*DK*/555, -1, /*DK*/555, 1, /*DK*/555, 3
        };
        const double b_good[] = {
            1.25, /*DK*/555, 1.75, /*DK*/555, 2.25, /*DK*/555, 2.75,
           -0.20, /*DK*/555, 1.00, /*DK*/555, 2.00, /*DK*/555, 3.00
        };
        const int ldb = sizeof(b)/(sizeof(b[0]))/nrhs;
        const int incb = 2;
        suzerain_bspline_lu_solve(nrhs, b, incb, ldb, luw);
        check_close_collections(
            b_good,
            b_good + sizeof(b_good)/sizeof(b_good[0]),
            b,
            b + sizeof(b)/sizeof(b[0]),
            1.0e-12);
    }

    suzerain_bspline_lu_free(luw);

    /***********************************/
    /* Complex-valued LU functionality */
    /***********************************/

    suzerain_bspline_luz_workspace *luzw
        = suzerain_bspline_luz_alloc(w);

    /* Form (2-3*i)*D[0] + (7-5*i)*D[1] operator in LU-ready banded storage.
     * Answer is
     *   -5+2i   7-5i   0-0i   0-0i
     *    0-0i  -5+2i   7-5i   0-0i
     *    0-0i   0-0i  -5+2i   7-5i
     *    0-0i   0-0i  -7+5i   9-8i
     * which, in LU-form where L has ones on the main diagonal, is
     *   -5+2i   7-5i       0+    0i       0+     0i
     *    0+0i  -5+2i       7-    5i       0+     0i
     *    0+0i   0+0i      -7+    5i       9-     8i
     *    0+0i   0+0i   45/74+11/74i   25/74-109/74i
     * The pivot matrix is
     *    1   0   0   0
     *    0   1   0   0
     *    0   0   0   1
     *    0   0   1   0
     * Check it in octave using [l,u,p] = lu(A).
     * Known good is in general banded matrix column-major order with
     * additional superdiagonal to allow for LU factorization fill-in.
     */
    {
        BOOST_TEST_MESSAGE("suzerain_bspline_luz_form_general");
        const double good_A0[][2] = {
            /*DK*/{0,0}, /*DK*/{0, 0}, {     -5,2},               {0,0},
            /*DK*/{0,0},       {7,-5}, {     -5,2},               {0,0},
                  {0,0},       {7,-5}, {     -7,5},         {45./74.,11./74.},
                  {0,0},       {9,-8}, {25./74.,-109./74.}, /*DK*/{0,0}
        };
        const double coeff[][2] = { {2.0, -3.0}, {7.0, -5.0}, {999.0, -999.0} };
        suzerain_bspline_luz_form_general(
            sizeof(coeff)/sizeof(coeff[0]), coeff, w, luzw);
        CHECK_GBMATRIX_CLOSE(
                     4,          4,        1,        2,  good_A0,        4,
            luzw->ndof, luzw->ndof, luzw->kl, luzw->ku,  luzw->A, luzw->ld,
            1e-12);
    }

    /* Check that multiple rhs solution works for operator found just above */
    {
        BOOST_TEST_MESSAGE("suzerain_bspline_luz_solve contiguous");
        const int nrhs = 2;
        double b[][2] = {
            { 1,7}, { 2,5}, {3,-2}, {4,-6},
            {-4,6}, {-1,3}, {1,-1}, {3, 6}
        };
        const double b_good[][2] = {
            // Vector 1
            {-3889./1033.,    88./1869.},
            {-4186./1565.,   377./1897.},
            {- 305./ 169.,    56./ 169.},
            {- 123./ 169.,  -  9./ 169.},
            // Vector 2
            { 2157./ 232., -1235./ 102.},
            { 4869./ 730., -7076./1245.},
            {  778./ 169., - 380./ 169.},
            {  557./ 169., - 120./ 169.}
        };
        const int incb = 1;
        const int ldb = sizeof(b)/(sizeof(b[0]))/nrhs;
        suzerain_bspline_luz_solve(nrhs, b, incb, ldb, luzw);
        /* Tolerance requirement adequate? condest(A) ~= 122.7 */
        /* Also, using approximate rationals via 'format rat'  */
        check_close_complex_collections(
            b_good, b_good + sizeof(b_good)/sizeof(b_good[0]),
            b, b + sizeof(b)/sizeof(b[0]),
            1.0e-5);
    }

    /* Check that multiple rhs solution works for operator found just above */
    {
        BOOST_TEST_MESSAGE("suzerain_bspline_luz_solve noncontiguous");
        const int nrhs = 2;
        double b[][2] = {
            { 1,7}, {55,-55}, { 2,5}, {66,-66}, {3,-2}, {77,-77}, {4,-6},
            {-4,6}, {55,-55}, {-1,3}, {66,-66}, {1,-1}, {77,-77}, {3, 6}
        };
        const double b_good[][2] = {
            // Vector 1
            {-3889./1033.,    88./1869.},
            {55,-55},
            {-4186./1565.,   377./1897.},
            {66,-66},
            {- 305./ 169.,    56./ 169.},
            {77,-77},
            {- 123./ 169.,  -  9./ 169.},
            // Vector 2
            { 2157./ 232., -1235./ 102.},
            {55,-55},
            { 4869./ 730., -7076./1245.},
            {66,-66},
            {  778./ 169., - 380./ 169.},
            {77,-77},
            {  557./ 169., - 120./ 169.}
        };
        const int incb = 2;
        const int ldb = sizeof(b)/(sizeof(b[0]))/nrhs;
        suzerain_bspline_luz_solve(nrhs, b, incb, ldb, luzw);
        /* Tolerance requirement adequate? condest(A) ~= 122.7 */
        /* Also, using approximate rationals via 'format rat'  */
        check_close_complex_collections(
            b_good, b_good + sizeof(b_good)/sizeof(b_good[0]),
            b, b + sizeof(b)/sizeof(b[0]),
            1.0e-5);
    }

    suzerain_bspline_luz_free(luzw);

    suzerain_bspline_free(w);
}

// Check a simple piecewise quadratic case's general banded storage
BOOST_AUTO_TEST_CASE( collocation_piecewise_quadratic_memory_application_soln )
{
    const double breakpoints[] = { 0.0, 1.0, 2.0, 3.0 };
    const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);
    const int order  = 3;
    const int nderiv = 1;

    suzerain_bspline_workspace *w
        = suzerain_bspline_alloc(order, nderiv, nbreak, breakpoints,
            SUZERAIN_BSPLINE_COLLOCATION_GREVILLE);
    BOOST_CHECK_EQUAL(suzerain_bspline_ndof(w), w->ndof);

    {
        /* Check w->D[0], the mass matrix, against known good solution:
         *  1.    0.    0.    0.    0.
         *  1./4. 5./8. 1./8. 0.    0.
         *  0.    1./8. 3./4. 1./8. 0.
         *  0.    0.    1./8. 5./8. 1./4.
         *  0.    0.    0.    0.    1.
         * Known good is in general banded matrix column-major order.
         */
        const double good_D0[] = { /*DK*/0.,   1.,      1./4.,
                                         0.,   5./8.,   1./8.,
                                      1./8.,   3./4.,   1./8.,
                                      1./8.,   5./8.,   0.,
                                      1./4.,      1.,   /*DK*/0. };
        CHECK_GBMATRIX_CLOSE(
                  5,       5,        1,        1, good_D0,     3,
            w->ndof, w->ndof, w->kl[0], w->ku[0], w->D[0], w->ld,
            1e-12);

        {
            /* Check w->D[0] application against multiple real vectors */
            const int nrhs = 2;
            double b[] = { 1, 2, 3, 4, 5,
                           5, 6, 7, 8, 9 };
            const double b_good[] = { 2*1., 2*15./8., 2*3., 2*33./8., 2*5.,
                                      2*5., 2*47./8., 2*7., 2*65./8., 2*9. };
            const int ldb = sizeof(b)/(sizeof(b[0]))/nrhs;
            suzerain_bspline_apply_operator(0, nrhs, 2.0, b, 1, ldb, w);
            BOOST_CHECK_EQUAL_COLLECTIONS(
                b_good, b_good + sizeof(b_good)/sizeof(b_good[0]),
                b, b + sizeof(b)/sizeof(b[0]));
        }

        {
            /* Check w->D[0] application against multiple complex vectors */
            const int nrhs = 2;
            double b[][2] = {
                { 1,  2}, { 3,  4}, { 5,  6}, { 7,  8}, { 9, 10},
                {11, 12}, {13, 14}, {15, 16}, {17, 18}, {19, 20}
            };
            const double b_good[][2] = {
                { 1., 2.}, { 2.75, 3.75}, { 5., 6.}, { 7.25, 8.25}, { 9.,10.},
                {11.,12.}, {12.75,13.75}, {15.,16.}, {17.25,18.25}, {19.,20.}
            };
            const int ldb = sizeof(b)/(sizeof(b[0]))/nrhs;
            suzerain_bspline_zapply_operator(0, nrhs, 1.0, b, 1, ldb, w);
            BOOST_CHECK_EQUAL_COLLECTIONS(
                (double *)  b_good,
                (double *) (b_good + sizeof(b_good)/sizeof(b_good[0])),
                (double *)  b,
                (double *) (b + sizeof(b)/sizeof(b[0])));
        }

        {
            /* Check w->D[0] accumulation against multiple real vectors */
            const int nrhs = 2;
            const double x[] = { 1, 2, 3, 4, 5,
                                 5, 6, 7, 8, 9 };
            const int ldx  = sizeof(x)/sizeof(x[0])/nrhs;
            const int incx = 1;
            double y[] = {  3, -4,  5, -6,  7,
                           -3,  4, -5,  6, -7 };
            const int ldy  = sizeof(y)/sizeof(y[0])/nrhs;
            const int incy = 1;
            const double y_good[] = {
                1.+3, 15./8.-4, 3.+5, 33./8.-6, 5.+7,
                5.-3, 47./8.+4, 7.-5, 65./8.+6, 9.-7
            };
            suzerain_bspline_accumulate_operator(
                0, nrhs, 1.0, x, incx, ldx, 1.0, y, incy, ldy, w);
            BOOST_CHECK_EQUAL_COLLECTIONS(
                y_good, y_good + sizeof(y_good)/sizeof(y_good[0]),
                y, y + sizeof(y)/sizeof(y[0]));
        }

        {
            /* Check w->D[0] accumulation against multiple complex vectors */
            /* using complex-valued coefficients                           */
            const int nrhs = 2;
            const double x[][2] = {
                { 1,  2}, { 3,  4}, { 5,  6}, { 7,  8}, { 9, 10},
                {11, 12}, {13, 14}, {15, 16}, {17, 18}, {19, 20}
            };
            const int ldx  = sizeof(x)/sizeof(x[0])/nrhs;
            const int incx = 1;
            double y[][2] = {
                { 1,  1}, { 1,  1}, { 1, -1}, {-1, -1}, {-1, -1},
                {-1, -1}, {-1, -1}, {-1,  1}, { 1,  1}, { 1,  1}
            };
            const int ldy  = sizeof(y)/sizeof(y[0])/nrhs;
            const int incy = 1;
            const double alpha[2] = { 2, 3 }; /* NB complex-valued */
            const double beta[2]  = { 5, 7 }; /* NB complex-valued */
            suzerain_bspline_zaccumulate_operator(
                0, nrhs, alpha, x, incx, ldx, beta, y, incy, ldy, w);
            const double y_good[][2] = {
                {- 6.00, 19.00}, {- 7.75,  27.75}, {  4.00,  29.00},
                                 {- 8.25,  26.25}, {-10.00,  35.00},
                {-12.00, 45.00}, {-13.75,  53.75}, {-30.00,  75.00},
                                 {-22.25, 100.25}, {-24.00, 109.00}
            };
            BOOST_CHECK_EQUAL_COLLECTIONS(
                (double *) y_good,
                (double *)(y_good + sizeof(y_good)/sizeof(y_good[0])),
                (double *) y,
                (double *)(y + sizeof(y)/sizeof(y[0])));
        }

        {
            /* Check w->D[0] accumulation against multiple complex vectors */
            /* using real-valued coefficients                              */
            const int nrhs = 2;
            const double x[][2] = {
                { 1,  2}, { 3,  4}, { 5,  6}, { 7,  8}, { 9, 10},
                {11, 12}, {13, 14}, {15, 16}, {17, 18}, {19, 20}
            };
            const int ldx  = sizeof(x)/sizeof(x[0])/nrhs;
            const int incx = 1;
            double y[][2] = {
                { 1,  1}, { 1,  1}, { 1, -1}, {-1, -1}, {-1, -1},
                {-1, -1}, {-1, -1}, {-1,  1}, { 1,  1}, { 1,  1}
            };
            const int ldy  = sizeof(y)/sizeof(y[0])/nrhs;
            const int incy = 1;
            const double alpha[2] = { 2, 0 }; /* NB real-valued */
            const double beta[2]  = { 5, 0 }; /* NB real-valued */
            suzerain_bspline_zaccumulate_operator(
                0, nrhs, alpha, x, incx, ldx, beta, y, incy, ldy, w);
            const double y_good[][2] = {
                { 7.,  9.}, {10.5, 12.5}, {15.,  7.}, { 9.5, 11.5}, {13., 15.},
                {17., 19.}, {20.5, 22.5}, {25., 37.}, {39.5, 41.5}, {43., 45.}
            };
            BOOST_CHECK_EQUAL_COLLECTIONS(
                (double *) y_good,
                (double *)(y_good + sizeof(y_good)/sizeof(y_good[0])),
                (double *) y,
                (double *)(y + sizeof(y)/sizeof(y[0])));
        }
    }

    /************************************************************************
     * Integration functionality
     * Exact integration solutions found using Mathematica 7 as follows:
     *  k = 3;
     *  b = {0, 1, 2, 3};
     *  n = Length[b] + k - 2;
     *  knots = Join[ ConstantArray[First[b], k - 1],
     *                b, ConstantArray[Last[b], k - 1]];
     *  B[i_, x_] := BSplineBasis[{k - 1, knots}, i, x]
     *  Table[Integrate[B[i, x], {x, Min[b], Max[b]}], {i, 0, n - 1}]
     ***********************************************************************/
    {
        const double expected[] = { 1./3, 2./3, 1, 2./3, 1./3 };
        check_close_collections(
            expected, expected + sizeof(expected)/sizeof(expected[0]),
            w->I, w->I + w->ndof, 1e-12);

        const double coeffs[] = { 1, 2, 3, 4, 5 };
        double value;
        suzerain_bspline_integrate(coeffs, /*incx*/1, &value, w);
        BOOST_CHECK_CLOSE(value, 9.0, 1e-12);

        const double zcoeffs[][2] = { {1,6}, {2,7}, {3,8}, {4,9}, {5,10} };
        double zvalue[2];
        suzerain_bspline_zintegrate(zcoeffs, /*incx*/1, &zvalue, w);
        BOOST_CHECK_CLOSE(zvalue[0],  9.0, 1e-12);
        BOOST_CHECK_CLOSE(zvalue[1], 24.0, 1e-12);
    }

    suzerain_bspline_free(w);
}

// Check a piecewise cubic case's general banded storage
BOOST_AUTO_TEST_CASE( collocation_piecewise_cubic_memory_application_soln )
{
    const double breakpoints[] = { 0.0, 1.0, 2.0, 3.0 };
    const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);
    const int order  = 4;
    const int nderiv = 2;

    suzerain_bspline_workspace *w
        = suzerain_bspline_alloc(order, nderiv, nbreak, breakpoints,
            SUZERAIN_BSPLINE_COLLOCATION_GREVILLE);
    BOOST_CHECK_EQUAL(suzerain_bspline_ndof(w), w->ndof);

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
            /*DK*/0,  /*DK*/0,          1.,       8./27.,         0,
            /*DK*/0,        0,         61./108.,  1./4.,          0,
                  0,       43./324.,    7./12.,   1./6.,          1./162.,
                  1./162.,  1./6.,      7./12.,   43./324.,       0,
                  0,        1./4.,     61./108.,  0,        /*DK*/0,
                  0,        8./27.,     1.,       /*DK*/0,  /*DK*/0
        };
        CHECK_GBMATRIX_CLOSE(
                  6,       6,        2,        2, good_D0,     5,
            w->ndof, w->ndof, w->kl[0], w->ku[0], w->D[0], w->ld,
            1e-12);

        /* Check w->D[0] application against multiple vectors */
        {
            const int nrhs = 2;
            double vapply[] = { 1, 2, 3, 4, 5, 6,
                                4, 5, 6, 7, 8, 9  };
            const double vapply_good[] = {
                1., 599./324., 35./12., 49./12., 1669./324., 6.,
                4., 1571./324., 71./12., 85./12., 2641./324., 9.
            };
            const int ldb = sizeof(vapply)/(sizeof(vapply[0]))/nrhs;
            suzerain_bspline_apply_operator(0, nrhs, 1.0, vapply, 1, ldb, w);
            check_close_collections(
                vapply_good,
                vapply_good + sizeof(vapply_good)/sizeof(vapply_good[0]),
                vapply,
                vapply + sizeof(vapply)/sizeof(vapply[0]),
                1.0e-12);
        }
    }

    /************************************************************************
     * Integration functionality
     * Exact integration solutions found using Mathematica 7 as follows:
     *  k = 4;
     *  b = {0, 1, 2, 3};
     *  n = Length[b] + k - 2;
     *  knots = Join[ ConstantArray[First[b], k - 1],
     *                b, ConstantArray[Last[b], k - 1]];
     *  B[i_, x_] := BSplineBasis[{k - 1, knots}, i, x]
     *  Table[Integrate[B[i, x], {x, Min[b], Max[b]}], {i, 0, n - 1}]
     ***********************************************************************/
    {
        const double expected[] = { 1./4, 1./2, 3./4, 3./4, 1./2, 1./4 };
        check_close_collections(
            expected, expected + sizeof(expected)/sizeof(expected[0]),
            w->I, w->I + w->ndof, 1e-12);

        const double coeffs[] = { 1, 2, 3, 4, 5, 6 };
        double value;
        suzerain_bspline_integrate(coeffs, /*incx*/1, &value, w);
        BOOST_CHECK_CLOSE(value, 21./2, 1e-12);

        const double zcoeffs[][2] = {
            {1,7}, {2,8}, {3,9}, {4,10}, {5,11}, {6,12}
        };
        double zvalue[2];
        suzerain_bspline_zintegrate(zcoeffs, /*incx*/1, &zvalue, w);
        BOOST_CHECK_CLOSE(zvalue[0], 21./2, 1e-12);
        BOOST_CHECK_CLOSE(zvalue[1], 57./2, 1e-12);
    }

    suzerain_bspline_free(w);
}

// Check a simple piecewise constant case's general banded storage
// See http://www.scribd.com/doc/52035371/Finding-Galerkin-L-2-based-Operators-for-B-spline-discretizations for the details.
BOOST_AUTO_TEST_CASE( galerkin_piecewise_constant_memory )
{
    const double breakpoints[] = { 0.0, 1.0, 9.0/8.0, 3.0/2.0, 2.0, 3.0 };
    const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);
    const int order  = 1;
    const int nderiv = 0;

    suzerain_bspline_workspace *w
        = suzerain_bspline_alloc(order, nderiv, nbreak, breakpoints,
            SUZERAIN_BSPLINE_GALERKIN_L2);
    BOOST_CHECK_EQUAL(suzerain_bspline_ndof(w), w->ndof);

    /* Check w->D[0], the mass matrix, against known good solution:
        *   1   0   0   0   0
        *   0  1/8  0   0   0
        *   0   0  3/8  0   0
        *   0   0   0  1/2  0
        *   0   0   0   0   1
        * Known good is in general banded matrix column-major order.
        */
    const double good_D0[] = { 1., 1./8., 3./8., 1./2., 1. };
    CHECK_GBMATRIX_CLOSE(
                5,       5,        0,        0, good_D0,     1,
        w->ndof, w->ndof, w->kl[0], w->ku[0], w->D[0], w->ld,
        1e-12);

    suzerain_bspline_free(w);
}

// Check a simple piecewise linear case's general banded storage
// See http://www.scribd.com/doc/52035371/Finding-Galerkin-L-2-based-Operators-for-B-spline-discretizations for the details.
BOOST_AUTO_TEST_CASE( galerkin_piecewise_linear_memory )
{
    const double breakpoints[] = { 0.0, 1.0, 9.0/8.0, 3.0/2.0, 2.0, 3.0 };
    const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);
    const int order  = 2;
    const int nderiv = 1;

    suzerain_bspline_workspace *w
        = suzerain_bspline_alloc(order, nderiv, nbreak, breakpoints,
            SUZERAIN_BSPLINE_GALERKIN_L2);
    BOOST_CHECK_EQUAL(suzerain_bspline_ndof(w), w->ndof);

    /* Check w->D[0], the mass matrix, against known good:
     * 1/3  1/6    0    0    0    0
     * 1/6  3/8   1/48  0    0    0
     *  0   1/48  1/6  1/16  0    0
     *  0    0    1/16 7/24 1/12  0
     *  0    0    0    1/12 1/2  1/6
     *  0    0    0    0    1/6  1/3
     * Known good is in general banded matrix column-major order.
     */
    const double good_D0[] = {0.    , 1./3. , 1./6. ,
                              1./6. , 3./8. , 1./48.,
                              1./48., 1./6. , 1./16.,
                              1./16., 7./24., 1./12.,
                              1./12., 1./2. , 1./6. ,
                              1./6. , 1./3. , 0.    };
    CHECK_GBMATRIX_CLOSE(
                6,       6,        1,        1, good_D0,   3,
        w->ndof, w->ndof, w->kl[0], w->ku[0], w->D[0], w->ld,
        1e-12);


    /* Check w->D[0], the first derivative matrix, against known good:
     * -1/2  1/2   0    0    0   0
     * -1/2   0   1/2   0    0   0
     *   0  -1/2   0   1/2   0   0
     *   0    0  -1/2   0   1/2  0
     *   0    0    0  -1/2   0  1/2
     *   0    0    0   0   -1/2 1/2
     * Known good is in general banded matrix column-major order.
     */
    const double good_D1[] = { 0.   , -1./2., -1./2.,
                               1./2.,  0.   , -1./2.,
                               1./2.,  0.   , -1./2.,
                               1./2.,  0.   , -1./2.,
                               1./2.,  0.   , -1./2.,
                               1./2.,  1./2.,  0.  };
    CHECK_GBMATRIX_CLOSE(
                6,       6,        1,        1, good_D1,   3,
        w->ndof, w->ndof, w->kl[1], w->ku[1], w->D[1], w->ld,
        1e-12);

    suzerain_bspline_free(w);
}

// Polynomial test helpers
typedef struct { int n; double c[]; } poly_params; // Flexible array

static
double poly_f(double x, void *params)
{
    poly_params *p = (poly_params *) params;
    return gsl_poly_eval(p->c, p->n, x);
}

static
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
        = (double *) malloc((nderiv+1) * w->ndof * sizeof(double));
    for (int i = 0; i <= nderiv; ++i) {
        suzerain_bspline_find_interpolation_problem_rhs(
                &f, expected + i*w->ndof, w);
        poly_params_differentiate(p);  // Drop the polynomial order by one
    }
    suzerain_bspline_lu_solve(nderiv+1, expected, 1, w->ndof, mass);

    // Make copies of the zeroth derivative coefficients
    double * const actual
        = (double *) malloc((nderiv+1) * w->ndof * sizeof(double));
    for (int i = 0; i <= nderiv; ++i) {
        memcpy(actual + i*w->ndof, expected, w->ndof * sizeof(actual[0]));
    }

    // Solve M*x' = D*x ...
    // ...starting by applying the derivative operators
    for (int i = 0; i <= nderiv; ++i) {
        suzerain_bspline_apply_operator(
                i, 1, 1.0, actual + i*w->ndof, 1, w->ndof, w);
    }
    // ...finish by solving with the mass matrix
    suzerain_bspline_lu_solve(nderiv+1, actual, 1, w->ndof, mass);

    // See if we got anywhere close
    for (int i = 0; i <= nderiv; ++i) {
        check_close_collections(
                expected + i*w->ndof, expected + (i+1)*w->ndof,
                actual + i*w->ndof, actual + (i+1)*w->ndof, 1.0e-09);
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
        double * coefficient = (double *) malloc(w->ndof * sizeof(double));
        suzerain_bspline_find_interpolation_problem_rhs(&f, coefficient, w);

        // Solve for function coefficients using the mass matrix
        suzerain_bspline_lu_solve(1, coefficient, 1, w->ndof, mass);

        // Take the n-th derivative of the coefficients using M x' = D x
        suzerain_bspline_apply_operator(
                derivative, 1, 1.0, coefficient, 1, w->ndof, w);
        suzerain_bspline_lu_solve(1, coefficient, 1, w->ndof, mass);

        // Ensure we recover the leading order, scaled monomial coefficients
        for (int i = 0; i < w->ndof; ++i) {
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
        double * coefficient = (double *) malloc(w->ndof * sizeof(double));
        suzerain_bspline_find_interpolation_problem_rhs(&f, coefficient, w);

        // Solve for function coefficients using the mass matrix
        suzerain_bspline_lu_solve(1, coefficient, 1, w->ndof, mass);

        // Take the n-th derivative of the coefficients using M x' = D x
        suzerain_bspline_apply_operator(
                derivative, 1, 1.0, coefficient, 1, w->ndof, w);
        suzerain_bspline_lu_solve(1, coefficient, 1, w->ndof, mass);

        // Ensure we recover the leading order, scaled monomial coefficients
        for (int i = 0; i < w->ndof; ++i) {
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
        double * coefficient = (double *) malloc(w->ndof * sizeof(double));
        suzerain_bspline_find_interpolation_problem_rhs(&f, coefficient, w);

        // Solve for function coefficients using the mass matrix
        suzerain_bspline_lu_solve(1, coefficient, 1, w->ndof, mass);

        // Take the n-th derivative of the coefficients using M x' = D x
        suzerain_bspline_apply_operator(
                derivative, 1, 1.0, coefficient, 1, w->ndof, w);
        suzerain_bspline_lu_solve(1, coefficient, 1, w->ndof, mass);

        // Ensure we recover the leading order, scaled monomial coefficients
        for (int i = 0; i < w->ndof; ++i) {
            BOOST_CHECK_CLOSE(6.0 * p->c[3], coefficient[i], 1e-11);
        }

        free(coefficient);
    }

    suzerain_bspline_lu_free(mass);
    suzerain_bspline_free(w);
    free(p);
}

static
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

BOOST_AUTO_TEST_CASE( bspline_evaluation_routine_for_real_coefficients )
{
    const double breakpoints[] = { 0.0, 1.0, 2.0, 3.0 };
    const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);
    const int order  = 4;
    const int nderiv = 2;
    const int ndof = 6;

    suzerain_bspline_workspace *w
        = suzerain_bspline_alloc(order, nderiv, nbreak, breakpoints,
            SUZERAIN_BSPLINE_COLLOCATION_GREVILLE);

    // Sanity check on fixed storage assumption for test case
    BOOST_REQUIRE_EQUAL(ndof, w->ndof);

    // Check that we have a partition of unity
    {
        double coefficients[ndof];
        for (int i = 0; i < ndof; ++i) coefficients[i] = 1.0;

        const double * points = breakpoints;
        const int npoints = sizeof(breakpoints)/sizeof(breakpoints[0]);

        const int ldvalues = npoints;
        const int nvalues = (nderiv+1) * ldvalues;
        double values[nvalues];

        // Compute 0...nderiv derivatives
        suzerain_bspline_evaluate(
                nderiv, coefficients, npoints, points, values, ldvalues, w);
        for (int i=0; i < ldvalues; ++i) {
            BOOST_CHECK_CLOSE(values[i], 1.0, 1.0e-13); // Partition of unity
        }
        for (int i=ldvalues; i < nvalues; ++i) {
            BOOST_CHECK_EQUAL(values[i], 0.0); // Higher derivatives all zero
        }

        // Compute only nderiv-th derivative
        suzerain_bspline_evaluate(
                nderiv, coefficients, npoints, points, values, 0, w);
        for (int i=0; i < ldvalues; ++i) {
            BOOST_CHECK_EQUAL(values[i], 0.0); // Partition of unity
        }
    }

    // Check basis function evaluation, derivatives, and ddot behavior
    // Note we assume that GSL bspline functionality is sound and
    // only do some minor spot checks here
    {
        BOOST_REQUIRE_EQUAL(6, w->ndof); // Sanity
        BOOST_REQUIRE_EQUAL(2, nderiv);  // Sanity

        double       coefficients[ndof];
        const double points[]             = { 1.0 };
        const int    npoints              = sizeof(points)/sizeof(points[0]);
        const int    ldvalues             = npoints;
        const int    nvalues              = (nderiv+1) *ldvalues;
        double       values[nvalues];

        // Investigate second basis function
        for (int i = 0; i < ndof; ++i) coefficients[i] = 0.0;
        coefficients[1] = 1.0;
        suzerain_bspline_evaluate(
                nderiv, coefficients, npoints, points, values, ldvalues, w);
        BOOST_CHECK_EQUAL( 1./4., values[0]); // 0th derivative
        BOOST_CHECK_EQUAL(-3./4., values[1]); // 1st derivative
        BOOST_CHECK_EQUAL( 3./2., values[2]); // 2nd derivative

        // Investigate third basis function
        for (int i = 0; i < ndof; ++i) coefficients[i] = 0.0;
        coefficients[2] = 1.0;
        suzerain_bspline_evaluate(
                nderiv, coefficients, npoints, points, values, ldvalues, w);
        BOOST_CHECK_CLOSE( 7./12., values[0], 1.0e-13); // 0th derivative
        BOOST_CHECK_EQUAL( 1./ 4., values[1]);          // 1st derivative
        BOOST_CHECK_EQUAL(-5./ 2., values[2]);          // 2nd derivative

        // Investigate fourth basis function
        for (int i = 0; i < ndof; ++i) coefficients[i] = 0.0;
        coefficients[3] = 1.0;
        suzerain_bspline_evaluate(
                nderiv, coefficients, npoints, points, values, ldvalues, w);
        BOOST_CHECK_EQUAL( 1./ 6., values[0]); // 0th derivative
        BOOST_CHECK_EQUAL( 1./ 2., values[1]); // 1st derivative
        BOOST_CHECK_EQUAL(     1., values[2]); // 2nd derivative

        // Investigate fifth basis function
        for (int i = 0; i < ndof; ++i) coefficients[i] = 0.0;
        coefficients[4] = 1.0;
        suzerain_bspline_evaluate(
                nderiv, coefficients, npoints, points, values, ldvalues, w);
        // All zero due to influence of endpoints/repeated knots
        BOOST_CHECK_EQUAL( 0., values[0]); // 0th derivative
        BOOST_CHECK_EQUAL( 0., values[1]); // 1st derivative
        BOOST_CHECK_EQUAL( 0., values[2]); // 2nd derivative

        // Check behavior of linear combinations of basis functions
        coefficients[0] = GSL_NAN; // Poison value checks ddot bounds
        coefficients[1] = 1.0;
        coefficients[2] = 2.0;
        coefficients[3] = 3.0;
        coefficients[4] = 4.0;
        coefficients[5] = GSL_NAN; // Poison value checks ddot bounds
        suzerain_bspline_evaluate(
                nderiv, coefficients, npoints, points, values, ldvalues, w);
        BOOST_CHECK_CLOSE(23./12., values[0], 1.0e-13); // 0th derivative
        BOOST_CHECK_CLOSE( 5./ 4., values[1], 1.0e-13); // 1st derivative
        BOOST_CHECK_CLOSE(-1./ 2., values[2], 1.0e-13); // 2nd derivative

        // Check behavior of linear combinations of basis functions
        // when ldvalues == 0 so only highest derivative is computed
        suzerain_bspline_evaluate(
                nderiv, coefficients, npoints, points, values, 0, w);
        BOOST_CHECK_CLOSE(-1./ 2., values[0], 1.0e-13); // 2nd derivative
    }

    suzerain_bspline_free(w);
}

BOOST_AUTO_TEST_CASE( bspline_evaluation_routine_for_complex_coefficients )
{
    const double breakpoints[] = { 0.0, 1.0, 2.0, 3.0 };
    const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);
    const int order  = 4;
    const int nderiv = 2;
    const int ndof = 6;

    suzerain_bspline_workspace *w
        = suzerain_bspline_alloc(order, nderiv, nbreak, breakpoints,
            SUZERAIN_BSPLINE_COLLOCATION_GREVILLE);

    // Sanity check on fixed storage assumption for test case
    BOOST_REQUIRE_EQUAL(ndof, w->ndof);

    // Check that we have a partition of unity for real-only coefficients
    {
        double coefficients[ndof][2];
        for (int i = 0; i < ndof; ++i) {
            coefficients[i][0] = 1.0;
            coefficients[i][1] = 0.0;
        }

        const double * points = breakpoints;
        const int npoints = sizeof(breakpoints)/sizeof(breakpoints[0]);

        const int ldvalues = npoints;
        const int nvalues = (nderiv+1) * ldvalues;
        double values[nvalues][2];

        // Compute 0...nderiv derivatives
        suzerain_bspline_zevaluate(
                nderiv, coefficients, npoints, points, values, ldvalues, w);
        for (int i=0; i < ldvalues; ++i) {
            BOOST_CHECK_CLOSE(values[i][0], 1.0, 1.0e-13); // Partition of 1
            BOOST_CHECK_EQUAL(values[i][1], 0.0);          // No imag content
        }
        for (int i=ldvalues; i < nvalues; ++i) {
            BOOST_CHECK_EQUAL(values[i][0], 0.0); // Derivatives all zero
            BOOST_CHECK_EQUAL(values[i][1], 0.0); // Derivatives all zero
        }

        // Compute only nderiv-th derivative
        suzerain_bspline_zevaluate(
                nderiv, coefficients, npoints, points, values, 0, w);
        for (int i=0; i < ldvalues; ++i) {
            BOOST_CHECK_EQUAL(values[i][0], 0.0); // Partition of unity
            BOOST_CHECK_EQUAL(values[i][1], 0.0); // Partition of unity
        }
    }

    // Check that we have a partition of unity for imag-only coefficients
    {
        double coefficients[ndof][2];
        for (int i = 0; i < ndof; ++i) {
            coefficients[i][0] = 0.0;
            coefficients[i][1] = 1.0;
        }

        const double * points = breakpoints;
        const int npoints = sizeof(breakpoints)/sizeof(breakpoints[0]);

        const int ldvalues = npoints;
        const int nvalues = (nderiv+1) * ldvalues;
        double values[nvalues][2];

        // Compute 0...nderiv derivatives
        suzerain_bspline_zevaluate(
                nderiv, coefficients, npoints, points, values, ldvalues, w);
        for (int i=0; i < ldvalues; ++i) {
            BOOST_CHECK_EQUAL(values[i][0], 0.0);          // No real content
            BOOST_CHECK_CLOSE(values[i][1], 1.0, 1.0e-13); // Partition of 1
        }
        for (int i=ldvalues; i < nvalues; ++i) {
            BOOST_CHECK_EQUAL(values[i][0], 0.0); // Derivatives all zero
            BOOST_CHECK_EQUAL(values[i][1], 0.0); // Derivatives all zero
        }

        // Compute only nderiv-th derivative
        suzerain_bspline_zevaluate(
                nderiv, coefficients, npoints, points, values, 0, w);
        for (int i=0; i < ldvalues; ++i) {
            BOOST_CHECK_EQUAL(values[i][0], 0.0); // Partition of unity
            BOOST_CHECK_EQUAL(values[i][1], 0.0); // Partition of unity
        }
    }

    // Check basis function evaluation, derivatives, and ddot behavior
    // Note we assume that GSL bspline functionality is sound and
    // only do some minor spot checks here
    {
        BOOST_REQUIRE_EQUAL(6, w->ndof); // Sanity
        BOOST_REQUIRE_EQUAL(2, nderiv);  // Sanity

        double       coefficients[ndof][2];
        const double points[]               = { 1.0 };
        const int    npoints                = sizeof(points)/sizeof(points[0]);
        const int    ldvalues               = npoints;
        const int    nvalues                = (nderiv+1) *ldvalues;
        double       values[nvalues][2];

        // Investigate second basis function
        for (int i = 0; i < ndof; ++i) {
            coefficients[i][0] = 0.0;
            coefficients[i][1] = 0.0;
        }
        coefficients[1][0] =  1.0;
        coefficients[1][1] = -0.5;

        suzerain_bspline_zevaluate(
                nderiv, coefficients, npoints, points, values, ldvalues, w);
        BOOST_CHECK_EQUAL( 1./4., values[0][0]); // Re(0th deriv)
        BOOST_CHECK_EQUAL(-1./8., values[0][1]); // Im(0th deriv)
        BOOST_CHECK_EQUAL(-3./4., values[1][0]); // Re(1st deriv)
        BOOST_CHECK_EQUAL( 3./8., values[1][1]); // Im(1st deriv)
        BOOST_CHECK_EQUAL( 3./2., values[2][0]); // Re(2nd deriv)
        BOOST_CHECK_EQUAL(-3./4., values[2][1]); // Im(2nd deriv)

        // Investigate third basis function
        for (int i = 0; i < ndof; ++i) {
            coefficients[i][0] = 0.0;
            coefficients[i][1] = 0.0;
        }
        coefficients[2][0] =  1.0;
        coefficients[2][1] = -0.5;
        suzerain_bspline_zevaluate(
                nderiv, coefficients, npoints, points, values, ldvalues, w);
        BOOST_CHECK_CLOSE( 7./12., values[0][0], 1.0e-13); // Re(0th deriv)
        BOOST_CHECK_CLOSE(-7./24., values[0][1], 1.0e-13); // Im(0th deriv)
        BOOST_CHECK_EQUAL( 1./ 4., values[1][0]);          // Re(1st deriv)
        BOOST_CHECK_EQUAL(-1./ 8., values[1][1]);          // Im(1st deriv)
        BOOST_CHECK_EQUAL(-5./ 2., values[2][0]);          // Re(2nd deriv)
        BOOST_CHECK_EQUAL( 5./ 4., values[2][1]);          // Im(2nd deriv)

        // Investigate fourth basis function
        for (int i = 0; i < ndof; ++i) {
            coefficients[i][0] = 0.0;
            coefficients[i][1] = 0.0;
        }
        coefficients[3][0] =  1.0;
        coefficients[3][1] = -0.5;
        suzerain_bspline_zevaluate(
                nderiv, coefficients, npoints, points, values, ldvalues, w);
        BOOST_CHECK_EQUAL( 1./ 6., values[0][0]); // Re(0th deriv)
        BOOST_CHECK_EQUAL(-1./12., values[0][1]); // Im(0th deriv)
        BOOST_CHECK_EQUAL( 1./ 2., values[1][0]); // Re(1st deriv)
        BOOST_CHECK_EQUAL(-1./ 4., values[1][1]); // Im(1st deriv)
        BOOST_CHECK_EQUAL(     1., values[2][0]); // Re(2nd deriv)
        BOOST_CHECK_EQUAL(  -0.5 , values[2][1]); // Im(2nd deriv)

        // Investigate fifth basis function
        for (int i = 0; i < ndof; ++i) {
            coefficients[i][0] = 0.0;
            coefficients[i][1] = 0.0;
        }
        coefficients[4][0] =  1.0;
        coefficients[4][1] = -0.5;
        suzerain_bspline_zevaluate(
                nderiv, coefficients, npoints, points, values, ldvalues, w);
        // All zero due to influence of endpoints/repeated knots
        BOOST_CHECK_EQUAL( 0., values[0][0]); // Re(0th deriv)
        BOOST_CHECK_EQUAL( 0., values[0][1]); // Im(0th deriv)
        BOOST_CHECK_EQUAL( 0., values[1][0]); // Re(1st deriv)
        BOOST_CHECK_EQUAL( 0., values[1][1]); // Im(1st deriv)
        BOOST_CHECK_EQUAL( 0., values[2][0]); // Re(2nd deriv)
        BOOST_CHECK_EQUAL( 0., values[2][1]); // Im(2nd deriv)

        // Check behavior of linear combinations of basis functions
        coefficients[0][0] = GSL_NAN; // Poison value checks ddot bounds
        coefficients[0][1] = GSL_NAN; // Poison value checks ddot bounds
        coefficients[1][0] =  1.0;
        coefficients[1][1] = -0.5;
        coefficients[2][0] =  2.0;
        coefficients[2][1] = -1.0;
        coefficients[3][0] =  3.0;
        coefficients[3][1] = -1.5;
        coefficients[4][0] =  4.0;
        coefficients[4][1] = -2.0;
        coefficients[5][0] = GSL_NAN; // Poison value checks ddot bounds
        coefficients[5][1] = GSL_NAN; // Poison value checks ddot bounds
        suzerain_bspline_zevaluate(
                nderiv, coefficients, npoints, points, values, ldvalues, w);
        BOOST_CHECK_CLOSE( 23./12., values[0][0], 1.0e-13); // Re(0th deriv)
        BOOST_CHECK_CLOSE(-23./24., values[0][1], 1.0e-13); // Im(0th deriv)
        BOOST_CHECK_CLOSE(  5./ 4., values[1][0], 1.0e-13); // Re(1st deriv)
        BOOST_CHECK_CLOSE(- 5./ 8., values[1][1], 1.0e-13); // Im(1st deriv)
        BOOST_CHECK_CLOSE(- 1./ 2., values[2][0], 1.0e-13); // Re(2nd deriv)
        BOOST_CHECK_CLOSE(  1./ 4., values[2][1], 1.0e-13); // Im(2nd deriv)

        // Check behavior of linear combinations of basis functions
        // when ldvalues == 0 so only highest deriv is computed
        suzerain_bspline_zevaluate(
                nderiv, coefficients, npoints, points, values, 0, w);
        BOOST_CHECK_CLOSE(-1./ 2., values[0][0], 1.0e-13); // Re(2nd deriv)
        BOOST_CHECK_CLOSE( 1./ 4., values[0][1], 1.0e-13); // Im(2nd deriv)
    }

    suzerain_bspline_free(w);
}

BOOST_AUTO_TEST_CASE( collocation_point_evaluation_is_operator_application )
{
    double breakpoints[35];
    const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);
    const int order  = 7;
    const int ndof   = nbreak + order - 2;
    const int nderiv = order - 1;
    double points[ndof];
    double coefficients[ndof];
    double values_eval[ndof];
    double values_apply[ndof];

    // Establish workspace and get collocation point information
    suzerain::math::logspace(0.1, 3.0, nbreak, breakpoints);
    suzerain_bspline_workspace *w
        = suzerain_bspline_alloc(order, nderiv, nbreak, breakpoints,
            SUZERAIN_BSPLINE_COLLOCATION_GREVILLE);
    BOOST_REQUIRE_EQUAL(ndof, w->ndof);
    suzerain_bspline_collocation_points(points, 1, w);

    // Generate "random" coefficients
    srand(634949092u);
    for (int i = 0; i < ndof; ++i) {
        coefficients[i] = 100.0 * (rand() / (RAND_MAX + 1.0)) - 50.0;
    }

    for (int k = 0; k <= nderiv; ++k) {
        BOOST_TEST_MESSAGE("Testing nderiv = " << k);

        // Evaluate coefficients times derivatives of basis functions
        BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bspline_evaluate(
                k, coefficients, ndof, points, values_eval, 0, w));

        // Apply derivative matrices directly
        std::memcpy(values_apply, coefficients, ndof*sizeof(double));
        BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, suzerain_bspline_apply_operator(
                k, 1, 1.0, values_apply, 1, ndof, w));

        // Check that both mechanisms give the same result
        check_close_collections(values_eval, values_eval + ndof,
                                values_apply, values_apply + ndof,
                                1.0e-12 * pow(10.0, k));
    }

    suzerain_bspline_free(w);
}
