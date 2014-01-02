/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2011-2014 Rhys Ulerich
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

// FIXME Missing BOOST_AUTO_TEST_SUITE( test_suite_name ) invocations ?

#include <suzerain/bspline.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <gsl/gsl_poly.h>

#include <suzerain/common.hpp>
#include <suzerain/countof.h>
#include <suzerain/function.h>
#include <suzerain/math.hpp>

#include "test_tools.hpp"

using suzerain::bspline;
using suzerain::bsplineop;
using suzerain::bsplineop_lu;
using suzerain::bsplineop_luz;

#define COUNTOF(x) SUZERAIN_COUNTOF(x)

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

#pragma warning(disable:1572 2014 2015)

BOOST_AUTO_TEST_CASE( allocation_okay )
{
    const double breakpts[] = { 0.0, 1.0, 2.0, 3.0 };
    bspline b(4, bspline::from_breakpoints(), COUNTOF(breakpts), breakpts);

    {
        bsplineop op(b, 2, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE);
        bsplineop_lu lu(op);
        bsplineop_luz luz(op);
    }

    {
        bsplineop op(b, 2, SUZERAIN_BSPLINEOP_GALERKIN_L2);
        bsplineop_lu lu(op);
        bsplineop_luz luz(op);
    }
}

BOOST_AUTO_TEST_CASE( collocation_piecewise_linear )
{
    const double breakpts[] = { 0.0, 1.0, 2.0, 3.0 };
    bspline b(2, bspline::from_breakpoints(), COUNTOF(breakpts), breakpts);
    bsplineop op(b, 2, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE);
    BOOST_REQUIRE_EQUAL(op.n(), b.n());

    {
        /* Check D[0], the mass matrix, against known good solution:
         *   1   0   0   0
         *   0   1   0   0
         *   0   0   1   0
         *   0   0   0   1
         * Known good is transposed in general banded column-major storage.
         */
        const double good_D_T0[] = { 1, 1, 1, 1 };
        CHECK_GBMATRIX_CLOSE(
                 4,      4,        0,       0,  good_D_T0,       1,
            op.n(), op.n(), op.kl(0), op.ku(0), op.D_T(0), op.ld(),
            std::numeric_limits<double>::epsilon()*1000);

        {
            /* Check D[0] application against multiple real vectors */
            const int nrhs = 2;
            double b[] = { 1, 2, 3, 4,
                           5, 6, 7, 8 };
            const int ldb = sizeof(b)/(sizeof(b[0]))/nrhs;
            const int incb = 1;
            op.apply(0, nrhs, 3.0, b, incb, ldb);
            const double b_good[] = { 3*1, 3*2, 3*3, 3*4,
                                      3*5, 3*6, 3*7, 3*8 };
            BOOST_CHECK_EQUAL_COLLECTIONS(
                b_good, b_good + sizeof(b_good)/sizeof(b_good[0]),
                b, b + sizeof(b)/sizeof(b[0]));
        }

        {
            /* Check D[0] application against multiple complex vectors */
            const int nrhs = 2;
            double b[][2] = {
                { 1,  2}, { 3,  4}, { 5,  6}, { 7,  8},
                {11, 12}, {13, 14}, {15, 16}, {17, 18}
            };
            const int ldb = sizeof(b)/(sizeof(b[0]))/nrhs;
            const int incb = 1;
            op.apply(0, nrhs, 1.0, b, incb, ldb);
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
            /* Check D[0] accumulation against multiple real vectors */
            const int nrhs = 2;
            const double x[] = { 1, 2, 3, 4,
                                 5, 6, 7, 8 };
            const int ldx  = sizeof(x)/sizeof(x[0])/nrhs;
            const int incx = 1;
            double y[] = {  1,  1,  1,  1,
                           -1, -1, -1, -1 };
            const int ldy  = sizeof(y)/sizeof(y[0])/nrhs;
            const int incy = 1;
            op.accumulate(0, nrhs, 1.0, x, incx, ldx, 1.0, y, incy, ldy);
            const double y_good[] = { 2, 3, 4, 5,
                                      4, 5, 6, 7 };
            BOOST_CHECK_EQUAL_COLLECTIONS(
                y_good,
                y_good + sizeof(y_good)/sizeof(y_good[0]),
                y,
                y + sizeof(y)/sizeof(y[0]));
        }

        {
            /* Check D[0] accumulation against multiple complex vectors */
            /* using complex-valued coefficients                        */
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
            op.accumulate(0, nrhs, alpha, x, incx, ldx, beta, y, incy, ldy);
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
            /* Check D[0] accumulation against multiple complex vectors */
            /* using real-valued coefficients                           */
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
            op.accumulate(0, nrhs, alpha, x, incx, ldx, beta, y, incy, ldy);
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
        /* Check D[1], the first derivative matrix, against known good:
         *  -1   1   0   0
         *   0  -1   1   0
         *   0   0  -1   1
         *   0   0  -1   1
         * Known good is transposed in general banded column-major storage.
         */
        const double good_D_T1[] = { /*DK*/0,   -1,       1,
                                           0,   -1,       1,
                                           0,   -1,       1,
                                          -1,    1, /*DK*/0 };
        CHECK_GBMATRIX_CLOSE(
                 4,      4,        1,        1, good_D_T1,       3,
            op.n(), op.n(), op.kl(1), op.ku(1), op.D_T(1), op.ld(),
            std::numeric_limits<double>::epsilon()*1000);

        /* Check D[1] application against multiple vectors */
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
        op.apply(1, nrhs, 1.0, vapply, incb, ldb);
        BOOST_CHECK_EQUAL_COLLECTIONS(
            vapply_good,
            vapply_good + sizeof(vapply_good)/sizeof(vapply_good[0]),
            vapply,
            vapply + sizeof(vapply)/sizeof(vapply[0]));
    }

    {
        /* Check D[2], the second derivative matrix, against zero result. */
        /* Known good is transposed in general banded column-major storage. */
        const double good_D_T2[] = { 0, 0, 0, 0 };
        CHECK_GBMATRIX_CLOSE(
                 4,      4,        0,        0, good_D_T2,       1,
            op.n(), op.n(), op.kl(2), op.ku(2), op.D_T(2), op.ld(),
            std::numeric_limits<double>::epsilon()*1000);
    }

    /********************************/
    /* Real-valued LU functionality */
    /********************************/

    bsplineop_lu lu(op);
    BOOST_REQUIRE_EQUAL(lu.n(), op.n());

    /* Form 2*D[0] - 3*D[1] operator which (not including any transpose) is:
     *   5  -3   0   0
     *   0   5  -3   0
     *   0   0   5  -3
     *   0   0   3  -1
     * Skip known good tests as they are tedious and error prone to code.
     * Instead we will sanity check factorization and sanity check the
     * action of the inverse on an extended collection of basis vectors.
     */
    {
        const double coeff[] = { 2.0, -3.0, 999.0 };
        lu.opform(sizeof(coeff)/sizeof(coeff[0]), coeff, op);

        // Factorize, checking the norm and conditioning are sane
        double norm = -123, rcond = -456;
        BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.opnorm(norm));
        BOOST_CHECK_GT(norm, 0.0);
        BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.factor());
        BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, lu.rcond(norm, rcond));
        BOOST_CHECK_GT(rcond, 0.0);
        BOOST_CHECK_LE(rcond, 1.0);

        const int nrhs = 4;
        double b[] = { 1, 0, 0, 0,
                       0, 1, 0, 0,
                       0, 0, 1, 0,
                       0, 0, 0, 1 };
        const double b_good[] = {
                         1./5.,       0,   0,        0,
                        3./25.,   1./5.,   0,        0,
                      -9./100., -3./20., -1./4., -3./4.,
                      27./100.,  9./20.,  3./4.,  5./4.
        };
        const int incb = 1;
        const int ldb = sizeof(b)/(sizeof(b[0]))/nrhs;
        lu.solve(nrhs, b, incb, ldb);
        /* Tolerance requirement adequate? condest(A) ~= 29.9 */
        /* Also, using approximate rationals via 'format rat' */
        check_close_collections(
            b_good, b_good + sizeof(b_good)/sizeof(b_good[0]),
            b, b + sizeof(b)/sizeof(b[0]),
            1.0e-10);
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

        lu.solve(nrhs, b, incb, ldb);
        check_close_collections(
            b_good,
            b_good + sizeof(b_good)/sizeof(b_good[0]),
            b,
            b + sizeof(b)/sizeof(b[0]),
            std::numeric_limits<double>::epsilon()*1000);
    }

    /***********************************/
    /* Complex-valued LU functionality */
    /***********************************/

    bsplineop_luz luz(op);
    BOOST_REQUIRE_EQUAL(luz.n(), op.n());

    /* Form (2-3*i)*D[0] + (7-5*i)*D[1]) in LU-ready banded storage.
     * The answer (not including any transpose used for storage) is:
     *   -5+2*i  7-5*i  0-0*i  0-0*i
     *    0-0*i -5+2*i  7-5*i  0-0*i
     *    0-0*i  0-0*i -5+2*i  7-5*i
     *    0-0*i  0-0*i -7+5*i  9-8*i
     * Skip known good tests as they are tedious and error prone to code.
     * Instead we will sanity check factorization and sanity check the
     * action of the inverse on an extended collection of basis vectors.
     */
    {
        BOOST_TEST_MESSAGE("suzerain_bspline_luz_form");
        const double coeff[][2] = {{2.0, -3.0}, {7.0, -5.0}, {999.0, -999.0}};
        luz.opform(sizeof(coeff)/sizeof(coeff[0]), coeff, op);

        // Factorize, checking the norm and conditioning are sane
        double norm = -123, rcond = -456;
        BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, luz.opnorm(norm));
        BOOST_CHECK_GT(norm, 0.0);
        BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, luz.factor());
        BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, luz.rcond(norm, rcond));
        BOOST_CHECK_GT(rcond, 0.0);
        BOOST_CHECK_LE(rcond, 1.0);

        const int nrhs = 8;
        double b[][2] = {
            {1,0}, {0,0}, {0,0}, {0,0},
            {0,1}, {0,0}, {0,0}, {0,0},
            {0,0}, {1,0}, {0,0}, {0,0},
            {0,0}, {0,1}, {0,0}, {0,0},
            {0,0}, {0,0}, {1,0}, {0,0},
            {0,0}, {0,0}, {0,1}, {0,0},
            {0,0}, {0,0}, {0,0}, {1,0},
            {0,0}, {0,0}, {0,0}, {0,1}
        };
        const double b_good[][2] = {
            {-5./29,-2./29}, {0,0}, {0,0}, {0,0},
            { 2./29,-5./29}, {0,0}, {0,0}, {0,0},
            {-247./841,-35./841}, {-5./29,-2./29}, {0,0}, {0,0},
            {35./841,-247./841}, {2./29,-5./29}, {0,0}, {0,0},
            {4059./2368,760./467}, {3923./4901,6099./4901}, {51./169,148./169}, {25./169,109./169},
            {-760./467,4059./2368}, {-56./45,3923./4901}, {-148./169,51./169}, {-109./169,25./169},
            {-430./393,-4824./3751}, {-643./1356,-4630./4901}, {-25./169,-109./169}, {1./169,-70./169},
            {4824./3751,-3963./3622}, {4630./4901,-2324./4901}, {109./169,-25./169}, {70./169,1./169}
        };
        const int incb = 1;
        const int ldb = sizeof(b)/(sizeof(b[0]))/nrhs;
        luz.solve(nrhs, b, incb, ldb);
        /* Tolerance requirement adequate? condest(A) ~= 122.7 */
        /* Also, using approximate rationals via 'format rat'  */
        check_close_collections(
            b_good, b_good + sizeof(b_good)/sizeof(b_good[0]),
            b, b + sizeof(b)/sizeof(b[0]),
            1.0e-5);
    }

    /* Check that multiple rhs solution works for operator found just above */
    /* Technically not necessary given basis results above, but can't hurt */
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
        luz.solve(nrhs, b, incb, ldb);
        /* Tolerance requirement adequate? condest(A) ~= 122.7 */
        /* Also, using approximate rationals via 'format rat'  */
        check_close_collections(
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
        luz.solve(nrhs, b, incb, ldb);
        /* Tolerance requirement adequate? condest(A) ~= 122.7 */
        /* Also, using approximate rationals via 'format rat'  */
        check_close_collections(
            b_good, b_good + sizeof(b_good)/sizeof(b_good[0]),
            b, b + sizeof(b)/sizeof(b[0]),
            1.0e-5);
    }
}

BOOST_AUTO_TEST_CASE( collocation_piecewise_quadratic)
{
    const double breakpts[] = { 0.0, 1.0, 2.0, 3.0 };
    bspline b(3, bspline::from_breakpoints(), COUNTOF(breakpts), breakpts);
    bsplineop op(b, 1, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE);
    BOOST_REQUIRE_EQUAL(op.n(), b.n());

    {
        /* Check w->D[0], the mass matrix, against known good solution:
         *  1.    0.    0.    0.    0.
         *  1./4. 5./8. 1./8. 0.    0.
         *  0.    1./8. 3./4. 1./8. 0.
         *  0.    0.    1./8. 5./8. 1./4.
         *  0.    0.    0.    0.    1.
         * Known good is transposed in general banded column-major storage.
         */
        const double good_D_T0[] = { /*DK*/0.,     1.,     0.,
                                        1./4.,  5./8.,  1./8.,
                                        1./8.,  3./4.,  1./8.,
                                        1./8.,  5./8.,  1./4.,
                                           0.,     1.,  /*DK*/0. };
        CHECK_GBMATRIX_CLOSE(
                 5,      5,       1,        1,  good_D_T0,       3,
            op.n(), op.n(), op.kl(0), op.ku(0), op.D_T(0), op.ld(),
            std::numeric_limits<double>::epsilon()*1000);

        {
            /* Check D[0] application against multiple real vectors */
            const int nrhs = 2;
            double b[] = { 1, 2, 3, 4, 5,
                           5, 6, 7, 8, 9 };
            const double b_good[] = { 2*1., 2*15./8., 2*3., 2*33./8., 2*5.,
                                      2*5., 2*47./8., 2*7., 2*65./8., 2*9. };
            const int ldb = sizeof(b)/(sizeof(b[0]))/nrhs;
            op.apply(0, nrhs, 2.0, b, 1, ldb);
            BOOST_CHECK_EQUAL_COLLECTIONS(
                b_good, b_good + sizeof(b_good)/sizeof(b_good[0]),
                b, b + sizeof(b)/sizeof(b[0]));
        }

        {
            /* Check D[0] application against multiple complex vectors */
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
            op.apply(0, nrhs, 1.0, b, 1, ldb);
            BOOST_CHECK_EQUAL_COLLECTIONS(
                (double *)  b_good,
                (double *) (b_good + sizeof(b_good)/sizeof(b_good[0])),
                (double *)  b,
                (double *) (b + sizeof(b)/sizeof(b[0])));
        }

        {
            /* Check D[0] accumulation against multiple real vectors */
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
            op.accumulate(0, nrhs, 1.0, x, incx, ldx, 1.0, y, incy, ldy);
            BOOST_CHECK_EQUAL_COLLECTIONS(
                y_good, y_good + sizeof(y_good)/sizeof(y_good[0]),
                y, y + sizeof(y)/sizeof(y[0]));
        }

        {
            /* Check D[0] accumulation against multiple complex vectors */
            /* using complex-valued coefficients                        */
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
            op.accumulate(0, nrhs, alpha, x, incx, ldx, beta, y, incy, ldy);
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
            /* Check D[0] accumulation against multiple complex vectors */
            /* using real-valued coefficients                           */
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
            op.accumulate(0, nrhs, alpha, x, incx, ldx, beta, y, incy, ldy);
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
}

BOOST_AUTO_TEST_CASE( collocation_piecewise_cubic )
{
    const double breakpts[] = { 0.0, 1.0, 2.0, 3.0 };
    bspline b(4, bspline::from_breakpoints(), COUNTOF(breakpts), breakpts);
    bsplineop op(b, 2, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE);
    BOOST_REQUIRE_EQUAL(op.n(), b.n());

    {
        /* Check D[0], the mass matrix, against known good solution:
         * 1.       0.        0.       0.        0.       0.
         * 8./27.  61./108.  43./324.  1./162.   0.       0.
         * 0.       1./4.     7./12.   1./6.     0.       0.
         * 0.       0.        1./6.    7./12.    1./4.    0.
         * 0.       0.        1./162. 43./324.  61./108.  8./27.
         * 0.       0.        0.       0.        0.       1.
         * Known good is transposed in general banded column-major storage.
         */
        const double good_D_T0[] = {
            /*DK*/0,  /*DK*/0,          1.,       0.,             0,
            /*DK*/0,   8./27.,    61./108., 43./324.,       1./162.,
                  0,    1./4.,      7./12.,    1./6.,             0,
                  0,    1./6.,      7./12.,    1./4.,             0,
            1./162., 43./324.,    61./108.,   8./27.,       /*DK*/0,
                  0,        0,          1.,  /*DK*/0,       /*DK*/0
        };
        CHECK_GBMATRIX_CLOSE(
                 6,      6,        2,        2, good_D_T0,       5,
            op.n(), op.n(), op.kl(0), op.ku(0), op.D_T(0), op.ld(),
            std::numeric_limits<double>::epsilon()*1000);

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
            op.apply(0, nrhs, 1.0, vapply, 1, ldb);
            check_close_collections(
                vapply_good,
                vapply_good + sizeof(vapply_good)/sizeof(vapply_good[0]),
                vapply,
                vapply + sizeof(vapply)/sizeof(vapply[0]),
                1.0e-12);
        }
    }
}

/* Check correctness when band storage is inferior to dense storage. */
BOOST_AUTO_TEST_CASE( collocation_piecewise_cubic_almost_dense )
{
    const double breakpts[] = { 0.0, 1.0 };
    bspline b(4, bspline::from_breakpoints(), COUNTOF(breakpts), breakpts);
    bsplineop op(b, 2, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE);
    BOOST_REQUIRE_EQUAL(op.n(), b.n());

    /* Check D[0], the mass matrix, against known good solution */
    /* Known good is transposed in general banded column-major storage. */
    const double good_D_T0[] = {
        0., 0.,         0.,     1.,     0.,     0.,     0.,
        0., 0.,     8./27.,  4./9.,  2./9., 1./27.,     0.,
        0., 1./27.,  2./9.,  4./9., 8./27.,     0.,     0.,
        0., 0.,         0.,     1.,     0.,     0.,     0.
    };
    CHECK_GBMATRIX_CLOSE(
                4,      4,        3,       3,  good_D_T0,       7,
        op.n(), op.n(), op.kl(0), op.ku(0), op.D_T(0), op.ld(),
        std::numeric_limits<double>::epsilon()*1000);

    /* Check D[1], the derivative matrix, against known good solution */
    /* Known good is transposed in general banded column-major storage. */
    const double good_D_T1[] = {
        0.,     0.,     0., -3.,    3.,    0., 0.,
        0.,     0., -4./3.,  0.,    1., 1./3., 0.,
        0., -1./3.,    -1.,  0., 4./3.,    0., 0.,
        0.,     0.,    -3.,  3.,    0.,    0., 0.
    };
    CHECK_GBMATRIX_CLOSE(
                4,      4,        3,       3,  good_D_T1,       7,
        op.n(), op.n(), op.kl(1), op.ku(1), op.D_T(1), op.ld(),
        std::numeric_limits<double>::epsilon()*1000);
}

/* Check correctness when band storage is inferior to dense storage. */
BOOST_AUTO_TEST_CASE( collocation_piecewise_septic_dense_and_then_some )
{
    const double breakpts[] = { 0.0, 1.0 };
    bspline b(8, bspline::from_breakpoints(), COUNTOF(breakpts), breakpts);
    bsplineop op(b, 2, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE);
    BOOST_REQUIRE_EQUAL(op.n(), b.n());

    /* Matrices are a mess to input so we check only their actions */

    {
        /* Check w->D[0] application against multiple real vectors */
        const int nrhs = 2;
        double b[] = { 2, 3, 5, 7, 11, 13, 17, 19,
                       3, 4, 6, 8, 12, 14, 18, 20 };
        const int ldb = sizeof(b)/(sizeof(b[0]))/nrhs;
        const int incb = 1;
        op.apply(0, nrhs, 823543, b, incb, ldb);
        const double b_good[] = {
            1647086,2767369,4378912,6397757,8694586,11121601,13549484,15647317,
            2470629,3590912,5202455,7221300,9518129,11945144,14373027,16470860
        };
        check_close_collections(
            b_good, b_good + sizeof(b_good)/sizeof(b_good[0]),
            b, b + sizeof(b)/sizeof(b[0]),
            std::numeric_limits<double>::epsilon()*1000);
    }

    {
        /* Check w->D[1] application against multiple real vectors */
        const int nrhs = 2;
        double b[] = { 2, 3, 5, 7, 11, 13, 17, 19,
                       19, 17, 13, 11, 7, 5, 3, 2 };
        const int ldb = sizeof(b)/(sizeof(b[0]))/nrhs;
        const int incb = 1;
        op.apply(1, nrhs, 16807, b, incb, ldb);
        const double b_good[] = {
            117649,197354,261593,311986,340553,350234,337249,235298,
           -235298,-337249,-350234,-340553,-311986,-261593,-197354,-117649
        };
        check_close_collections(
            b_good, b_good + sizeof(b_good)/sizeof(b_good[0]),
            b, b + sizeof(b)/sizeof(b[0]),
            std::numeric_limits<double>::epsilon()*1000);
    }
}

// See http://www.scribd.com/doc/52035371/Finding-Galerkin-L-2-based-Operators-for-B-spline-discretizations for the details.
BOOST_AUTO_TEST_CASE( galerkin_piecewise_constant )
{
    const double breakpts[] = { 0.0, 1.0, 9.0/8.0, 3.0/2.0, 2.0, 3.0 };
    bspline b(1, bspline::from_breakpoints(), COUNTOF(breakpts), breakpts);
    bsplineop op(b, 0, SUZERAIN_BSPLINEOP_GALERKIN_L2);
    BOOST_REQUIRE_EQUAL(op.n(), b.n());

    /* Check w->D[0], the mass matrix, against known good solution:
        *   1   0   0   0   0
        *   0  1/8  0   0   0
        *   0   0  3/8  0   0
        *   0   0   0  1/2  0
        *   0   0   0   0   1
        * Known good is in general banded matrix column-major order.
        */
    const double good_D_T0[] = { 1., 1./8., 3./8., 1./2., 1. };
    CHECK_GBMATRIX_CLOSE(
                 5,      5,        0,        0, good_D_T0,     1,
            op.n(), op.n(), op.kl(0), op.ku(0), op.D_T(0), op.ld(),
            std::numeric_limits<double>::epsilon()*1000);
    CHECK_GBMATRIX_SYMMETRIC( // Mass matrix is analytically symmetric
            op.n(), op.n(), op.kl(0), op.ku(0), op.D_T(0), op.ld());
}

// See http://www.scribd.com/doc/52035371/Finding-Galerkin-L-2-based-Operators-for-B-spline-discretizations for the details.
BOOST_AUTO_TEST_CASE( galerkin_piecewise_linear )
{
    const double breakpts[] = { 0.0, 1.0, 9.0/8.0, 3.0/2.0, 2.0, 3.0 };
    bspline b(2, bspline::from_breakpoints(), COUNTOF(breakpts), breakpts);
    bsplineop op(b, 1, SUZERAIN_BSPLINEOP_GALERKIN_L2);
    BOOST_REQUIRE_EQUAL(op.n(), b.n());

    /* Check D[0], the mass matrix, against known good:
     * 1/3  1/6    0    0    0    0
     * 1/6  3/8   1/48  0    0    0
     *  0   1/48  1/6  1/16  0    0
     *  0    0    1/16 7/24 1/12  0
     *  0    0    0    1/12 1/2  1/6
     *  0    0    0    0    1/6  1/3
     * Known good is transposed in general banded column-major storage.
     */
    const double good_D_T0[] = {0.    , 1./3. , 1./6. ,
                                1./6. , 3./8. , 1./48.,
                                1./48., 1./6. , 1./16.,
                                1./16., 7./24., 1./12.,
                                1./12., 1./2. , 1./6. ,
                                1./6. , 1./3. , 0.    };
    CHECK_GBMATRIX_CLOSE(
                 6,       6,       1,        1, good_D_T0,       3,
            op.n(), op.n(), op.kl(0), op.ku(0), op.D_T(0), op.ld(),
            std::numeric_limits<double>::epsilon()*1000);
    CHECK_GBMATRIX_SYMMETRIC( // Mass matrix is analytically symmetric
            op.n(), op.n(), op.kl(0), op.ku(0), op.D_T(0), op.ld());

    // FIXME What did this test block accomplish exactly?
    for (int i = 0; i < 5; ++i) {
        const int offset = (i * op.ld()) + (op.ku(0) + 1);
        BOOST_CHECK_EQUAL(op.D_T(0)[offset], op.D_T(0)[offset + op.kl(0)]);
    }

    /* Check D[1], the first derivative matrix, against known good:
     * -1/2  1/2   0    0    0   0
     * -1/2   0   1/2   0    0   0
     *   0  -1/2   0   1/2   0   0
     *   0    0  -1/2   0   1/2  0
     *   0    0    0  -1/2   0  1/2
     *   0    0    0   0   -1/2 1/2
     * Known good is in general banded matrix column-major order.
     */
    const double good_D_T1[] = {  0.   , -1./2., 1./2.,
                                 -1./2.,  0.   , 1./2.,
                                 -1./2.,  0.   , 1./2.,
                                 -1./2.,  0.   , 1./2.,
                                 -1./2.,  0.   , 1./2.,
                                 -1./2.,  1./2.,  0.  };
    CHECK_GBMATRIX_CLOSE(
                6,       6,        1,        1, good_D_T1,       3,
            op.n(), op.n(), op.kl(1), op.ku(1), op.D_T(1), op.ld(),
            std::numeric_limits<double>::epsilon()*1000);
}

// See http://www.scribd.com/doc/52035371/Finding-Galerkin-L-2-based-Operators-for-B-spline-discretizations for the details.
BOOST_AUTO_TEST_CASE( galerkin_piecewise_quadratic )
{
    const double breakpts[] = { 0.0, 1.0, 9.0/8.0, 3.0/2.0, 2.0, 3.0 };
    bspline b(3, bspline::from_breakpoints(), COUNTOF(breakpts), breakpts);
    bsplineop op(b, 2, SUZERAIN_BSPLINEOP_GALERKIN_L2);
    BOOST_REQUIRE_EQUAL(op.n(), b.n());

    /* Check D[0], the mass matrix, against known good
     * transposed in general banded matrix column-major order.
     */
    const double good_D_T0[] = {
            0.      , 0.          , 1./5.    ,  14./135.    ,  4./135. ,
            0.      , 14./135.    , 19./120. ,  65./576.    ,  1./8640.,
            4./135. , 65./576.    , 107./360.,  851./15120. ,  9./2240.,
            1./8640., 851./15120. , 5./28.   ,  1919./20160.,  1./315. ,
            9./2240., 1919./20160., 59./168. ,  16./105.    ,  1./45.  ,
            1./315. , 16./105.    , 7./30.   ,  1./9.       ,  0.      ,
            1./45.  , 1./9.       , 1./5.    ,  0.          ,  0.
    };
    CHECK_GBMATRIX_CLOSE(
                 7,      7,        2,        2, good_D_T0,       5,
            op.n(), op.n(), op.kl(0), op.ku(0), op.D_T(0), op.ld(),
            std::numeric_limits<double>::epsilon()*1000);
    CHECK_GBMATRIX_SYMMETRIC( // Mass matrix is analytically symmetric
            op.n(), op.n(), op.kl(0), op.ku(0), op.D_T(0), op.ld());

    {
        /* Check D[0] application against multiple real vectors */
        const int nrhs = 2;
        double b[] = {  1,  2,  3,  4,  5,  6,  7,
                       11, 12, 13, 14, 15, 16, 17 };
        const int ldb = sizeof(b)/(sizeof(b[0]))/nrhs;
        const int incb = 1;
        op.apply(0, nrhs, 60480.0, b, incb, ldb);
        const double b_good[] = {
             30016,  45927,  84201,  83363, 194661, 178560, 131712,
            231616, 272727, 386601, 284963, 572661, 480960, 333312
        };
        check_close_collections(
            b_good, b_good + sizeof(b_good)/sizeof(b_good[0]),
            b, b + sizeof(b)/sizeof(b[0]),
            std::numeric_limits<double>::epsilon()*1000);
    }

    {
        /* Check D[0] application against multiple complex vectors */
        const int nrhs = 2;
        double b[][2] = {
            {  1, 2}, { 2, 3}, { 3, 4}, { 4, 5}, { 5, 6}, { 6, 7}, {7,  8},
            { 11,12}, {12,13}, {13,14}, {14,15}, {15,16}, {16,17}, {17,18}
        };
        const int ldb = sizeof(b)/(sizeof(b[0]))/nrhs;
        const int incb = 1;
        op.apply(0, nrhs, 60480.0, b, incb, ldb);
        const double b_good[][2] = {
            {30016,50176},{45927,68607},{84201,114441},{83363,103523},{194661,232461},{178560,208800},{131712,151872},
            {231616,251776},{272727,295407},{386601,416841},{284963,305123},{572661,610461},{480960,511200},{333312,353472}
        };
        check_close_collections(
            (double *) b_good,
            (double *)(b_good + sizeof(b_good)/sizeof(b_good[0])),
            (double *) b,
            (double *)(b + sizeof(b)/sizeof(b[0])),
            std::numeric_limits<double>::epsilon()*1000);
    }

    {
        /* Check D[0] accumulation against multiple real vectors */
        const int nrhs = 2;
        double x[] = {  1,  2,  3,  4,  5,  6,  7,
                       11, 12, 13, 14, 15, 16, 17 };
        const int ldx  = sizeof(x)/sizeof(x[0])/nrhs;
        const int incx = 1;
        double y[] = {  1,  2,  3,  4,  5,  6,  7,
                       -1, -2, -3, -4, -5, -6, -7 };
        const int ldy  = sizeof(y)/sizeof(y[0])/nrhs;
        const int incy = 1;
        op.accumulate(0, nrhs, 60480.0, x, incx, ldx, 1.0, y, incy, ldy);
        const double y_good[] = {
                30017,45929,84204,83367,194666,178566,131719,
                231615,272725,386598,284959,572656,480954,333305
        };
        check_close_collections(
            y_good, y_good + sizeof(y_good)/sizeof(y_good[0]),
            y, y + sizeof(y)/sizeof(y[0]),
            std::numeric_limits<double>::epsilon()*1000);
    }

    {
        /* Check D[0] accumulation against multiple complex vectors */
        /* using complex-valued coefficients                        */
        const int nrhs = 2;
        const double x[][2] = {
            { 1, 2},{ 2, 3},{ 3, 4},{ 4, 5},{ 5, 6},{ 6, 7},{ 7, 8},
            {11,12},{12,13},{13,14},{14,15},{15,16},{16,17},{17,18}
        };
        const int ldx  = sizeof(x)/sizeof(x[0])/nrhs;
        const int incx = 1;
        double y[][2] = {
            { 1, 2},{ 2, 3},{ 3, 4},{ 4, 5},{ 5, 6},{ 6, 7},{ 7, 8},
            {-1,-2},{-2,-3},{-3,-4},{-4,-5},{-5,-6},{-6,-7},{-7,-8},
        };
        const int ldy  = sizeof(y)/sizeof(y[0])/nrhs;
        const int incy = 1;
        const double alpha[2] = { 60480, -2*60480 }; /* NB complex-valued */
        const double beta[2]  = {     6,        7 }; /* NB complex-valued */
        op.accumulate(0, nrhs, alpha, x, incx, ldx, beta, y, incy, ldy);
        const double y_good[][2] = {
            {130360,-9837},{183132,-23215},{313073,-53916},{290398,-63145},{659571,-156790},{596147,-148236},{435442,-111455},
            {735176,-211475},{863550,-250079},{1220293,-356406},{895220,-264861},{1793595,-534932},{1503373,-450804},{1040270,-313249}
        };
        check_close_collections(
            (double *) y_good,
            (double *)(y_good + sizeof(y_good)/sizeof(y_good[0])),
            (double *) y,
            (double *)(y + sizeof(y)/sizeof(y[0])),
            std::numeric_limits<double>::epsilon()*1000);
    }

    {
        /* Check D[0] accumulation against multiple complex vectors */
        /* using real-valued coefficients                           */
        const int nrhs = 2;
        const double x[][2] = {
            { 1, 2},{ 2, 3},{ 3, 4},{ 4, 5},{ 5, 6},{ 6, 7},{ 7, 8},
            {11,12},{12,13},{13,14},{14,15},{15,16},{16,17},{17,18}
        };
        const int ldx  = sizeof(x)/sizeof(x[0])/nrhs;
        const int incx = 1;
        double y[][2] = {
            { 1, 2},{ 2, 3},{ 3, 4},{ 4, 5},{ 5, 6},{ 6, 7},{ 7, 8},
            {-1,-2},{-2,-3},{-3,-4},{-4,-5},{-5,-6},{-6,-7},{-7,-8},
        };
        const int ldy  = sizeof(y)/sizeof(y[0])/nrhs;
        const int incy = 1;
        const double alpha = 60480;
        const double beta  =     6;
        op.accumulate(0, nrhs, alpha, x, incx, ldx, beta, y, incy, ldy);
        const double y_good[][2] = {
            {30022,50188},{45939,68625},{84219,114465},{83387,103553},{194691,232497},{178596,208842},{131754,151920},
            {231610,251764},{272715,295389},{386583,416817},{284939,305093},{572631,610425},{480924,511158},{333270,353424}
        };
        check_close_collections(
            (double *) y_good,
            (double *)(y_good + sizeof(y_good)/sizeof(y_good[0])),
            (double *) y,
            (double *)(y + sizeof(y)/sizeof(y[0])),
            std::numeric_limits<double>::epsilon()*1000);
    }

    /* Check D[1], the first derivative matrix, against known good
     * transposed in general banded matrix column-major order.
     */
    const double good_D_T1[] = {
            0.      ,         0.  ,  -1./2., 19./54.  , 4./27. ,
            0.      ,   -19./54.  ,   0.   , 25./72.  , 1./216.,
            -4./27. ,   -25./72.  ,   0.   , 167./378., 3./56. ,
            -1./216.,   -167./378.,   0.   , 209./504., 2./63. ,
            -3./56. ,   -209./504.,   0.   , 5./14.   , 1./9.  ,
            -2./63. ,   -5./14.   ,   0.   , 7./18.   , 0.     ,
            -1./9.  ,   -7./18.   ,   1./2., 0.       , 0.
    };
    CHECK_GBMATRIX_CLOSE(
                7,       7,        2,        2, good_D_T1,       5,
            op.n(), op.n(), op.kl(1), op.ku(1), op.D_T(1), op.ld(),
            std::numeric_limits<double>::epsilon()*1000);

    {
        /* Check D[1] application against multiple real vectors */
        const int nrhs = 2;
        double b[] = {  2,  3,  5,  7,  6,  4,  1,
                        1,  4,  6,  7,  5,  3,  2 };
        const int ldb = sizeof(b)/(sizeof(b[0]))/nrhs;
        const int incb = 1;
        op.apply(1, nrhs, 1512.0, b, incb, ldb);
        const double b_good[] = {
            1204,1610,3139,593,-2466,-2988,-2604,
            2716,2667,2757,-757,-2919,-1860,-1092
        };
        check_close_collections(
            b_good, b_good + sizeof(b_good)/sizeof(b_good[0]),
            b, b + sizeof(b)/sizeof(b[0]),
            std::numeric_limits<double>::epsilon()*10000); /* tol drops */
    }

    {
        /* Check D[1] application against multiple complex vectors */
        const int nrhs = 2;
        double b[][2] = {
            { 2, 1},{ 3, 4},{ 5, 6},{ 7, 7},{ 6, 5},{ 4, 3},{ 1, 2},
            { 1, 2},{ 4, 3},{ 6, 5},{ 7, 7},{ 5, 6},{ 3, 4},{ 2, 1}
        };
        const int ldb = sizeof(b)/(sizeof(b[0]))/nrhs;
        const int incb = 1;
        op.apply(1, nrhs, 1512.0, b, incb, ldb);
        const double b_good[][2] = {
            {1204,+2716},{1610,+2667},{3139,+2757},{593,-757},{-2466,-2919},{-2988,-1860},{-2604,-1092},
            {2716,+1204},{2667,+1610},{2757,+3139},{-757,+593},{-2919,-2466},{-1860,-2988},{-1092,-2604}
        };
        check_close_collections(
            (double *) b_good,
            (double *)(b_good + sizeof(b_good)/sizeof(b_good[0])),
            (double *) b,
            (double *)(b + sizeof(b)/sizeof(b[0])),
            std::numeric_limits<double>::epsilon()*10000); // drop tolerance
    }

    {
        /* Check D[1] accumulation against multiple real vectors */
        const int nrhs = 2;
        double x[] = {  2,  3,  5,  7,  6,  4,  1,
                        1,  4,  6,  7,  5,  3,  2 };
        const int ldx  = sizeof(x)/sizeof(x[0])/nrhs;
        const int incx = 1;
        double y[] = {  1,  2,  3,  4,  5,  6,  7,
                       -1, -2, -3, -4, -5, -6, -7 };
        const int ldy  = sizeof(y)/sizeof(y[0])/nrhs;
        const int incy = 1;
        op.accumulate(1, nrhs, 1512.0, x, incx, ldx, 1.0, y, incy, ldy);
        const double y_good[] = {
            1205,1612,3142,597,-2461,-2982,-2597,
            2715,2665,2754,-761,-2924,-1866,-1099
        };
        check_close_collections(
            y_good, y_good + sizeof(y_good)/sizeof(y_good[0]),
            y, y + sizeof(y)/sizeof(y[0]),
            std::numeric_limits<double>::epsilon()*10000); // drop tolerance
    }

    {
        /* Check D[1] accumulation against multiple complex vectors */
        /* using complex-valued coefficients                        */
        const int nrhs = 2;
        const double x[][2] = {
            { 2, 1},{ 3, 4},{ 5, 6},{ 7, 7},{ 6, 5},{ 4, 3},{ 1, 2},
            { 1, 2},{ 4, 3},{ 6, 5},{ 7, 7},{ 5, 6},{ 3, 4},{ 2, 1}
        };
        const int ldx  = sizeof(x)/sizeof(x[0])/nrhs;
        const int incx = 1;
        double y[][2] = {
            { 1, 2},{ 2, 3},{ 3, 4},{ 4, 5},{ 5, 6},{ 6, 7},{ 7, 8},
            {-1,-2},{-2,-3},{-3,-4},{-4,-5},{-5,-6},{-6,-7},{-7,-8},
        };
        const int ldy  = sizeof(y)/sizeof(y[0])/nrhs;
        const int incy = 1;
        const double alpha[2] = { 1512, -2*1512 }; /* NB complex-valued */
        const double beta[2]  = {    6,       7 }; /* NB complex-valued */
        op.accumulate(1, nrhs, alpha, x, incx, ldx, beta, y, incy, ldy);
        const double y_good[][2] = {
            {6628,+327},{6935,-521},{8643,-3476},{-932,-1885},{-8316,+2084},{-6721,+4200},{-4802,+4213},
            {5132,-4247},{5896,-3756},{9045,-2420},{440,+2049},{-7839,+3301},{-7823,+648},{-6286,-517}
        };
        check_close_collections(
            (double *) y_good,
            (double *)(y_good + sizeof(y_good)/sizeof(y_good[0])),
            (double *) y,
            (double *)(y + sizeof(y)/sizeof(y[0])),
            std::numeric_limits<double>::epsilon()*10000);
    }

    {
        /* Check D[1] accumulation against multiple complex vectors */
        /* using real-valued coefficients                           */
        const int nrhs = 2;
        const double x[][2] = {
            { 2, 1},{ 3, 4},{ 5, 6},{ 7, 7},{ 6, 5},{ 4, 3},{ 1, 2},
            { 1, 2},{ 4, 3},{ 6, 5},{ 7, 7},{ 5, 6},{ 3, 4},{ 2, 1}
        };
        const int ldx  = sizeof(x)/sizeof(x[0])/nrhs;
        const int incx = 1;
        double y[][2] = {
            { 1, 2},{ 2, 3},{ 3, 4},{ 4, 5},{ 5, 6},{ 6, 7},{ 7, 8},
            {-1,-2},{-2,-3},{-3,-4},{-4,-5},{-5,-6},{-6,-7},{-7,-8},
        };
        const int ldy  = sizeof(y)/sizeof(y[0])/nrhs;
        const int incy = 1;
        const double alpha = 1512;
        const double beta  =    6;
        op.accumulate(1, nrhs, alpha, x, incx, ldx, beta, y, incy, ldy);
        const double y_good[][2] = {
            {1210,+2728},{1622,+2685},{3157,+2781},{617,-727},{-2436,-2883},{-2952,-1818},{-2562,-1044},
            {2710,+1192},{2655,+1592},{2739,+3115},{-781,+563},{-2949,-2502},{-1896,-3030},{-1134,-2652}
        };
        check_close_collections(
            (double *) y_good,
            (double *)(y_good + sizeof(y_good)/sizeof(y_good[0])),
            (double *) y,
            (double *)(y + sizeof(y)/sizeof(y[0])),
            std::numeric_limits<double>::epsilon()*10000);
    }

    /* Check D[2], the second derivative matrix, against known good
     * transposed in general banded matrix column-major order.
     */
    const double good_D_T2[] = {
                 0.,         0.,     2./3.,   -34./27.,  16./27.,
                 0.,    20./27.,    -4./3.,      4./9.,   4./27.,
            16./27.,      4./9.,   -32./9.,  368./189.,    4./7.,
             4./27.,  368./189.,  -64./21.,    44./63.,  16./63.,
              4./7.,    44./63.,  -40./21.,     4./21.,    4./9.,
            16./63.,     4./21.,    -4./3.,      8./9.,       0.,
              4./9.,    -10./9.,     2./3.,         0.,       0.
    };
    CHECK_GBMATRIX_CLOSE(
                7,       7,        2,        2, good_D_T2,       5,
            op.n(), op.n(), op.kl(2), op.ku(2), op.D_T(2), op.ld(),
            std::numeric_limits<double>::epsilon()*1000);
}

// See http://www.scribd.com/doc/52035371/Finding-Galerkin-L-2-based-Operators-for-B-spline-discretizations for the details.
BOOST_AUTO_TEST_CASE( galerkin_piecewise_cubic )
{
    const double breakpts[] = { 0.0, 1.0, 9.0/8.0, 3.0/2.0, 2.0, 3.0 };
    bspline b(4, bspline::from_breakpoints(), COUNTOF(breakpts), breakpts);
    bsplineop op(b, 3, SUZERAIN_BSPLINEOP_GALERKIN_L2);
    BOOST_REQUIRE_EQUAL(op.n(), b.n());

    /* Check D[0], the mass matrix, against known good
     * transposed in general banded matrix column-major order.
     */
    const double good_D_T0[] = {
                 0.,                 0.,                 0.,            1./7.,         121./1620.,           16./567.,       4./945.,
                 0.,                 0.,         121./1620.,       257./2520.,       7901./96768.,    66553./2903040.,   1./2903040.,
                 0.,           16./567.,       7901./96768.,     9263./60480.,   111205./1016064.,   63109./25401600.,   27./627200.,
            4./945.,    66553./2903040.,   111205./1016064.,   56569./211680.,     18281./211680.,        347./35840.,      1./8820.,
        1./2903040.,   63109./25401600.,     18281./211680.,   85759./352800.,   507307./3763200.,         271./8820.,      4./1575.,
        27./627200.,        347./35840.,   507307./3763200.,         75./392.,         647./5880.,            4./175.,            0.,
           1./8820.,         271./8820.,         647./5880.,         16./105.,         103./1260.,                 0.,            0.,
           4./1575.,            4./175.,         103./1260.,            1./7.,                 0.,                 0.,            0.
    };
    CHECK_GBMATRIX_CLOSE(
                8,       8,        3,        3, good_D_T0,       7,
            op.n(), op.n(), op.kl(0), op.ku(0), op.D_T(0), op.ld(),
            std::numeric_limits<double>::epsilon()*1000);
    CHECK_GBMATRIX_SYMMETRIC( // Mass matrix is analytically symmetric
            op.n(), op.n(), op.kl(0), op.ku(0), op.D_T(0), op.ld());

    /* Check D[1], the first derivative matrix, against known good
     * transposed in general banded matrix column-major order.
     */
    const double good_D_T1[] = {
               0.,              0.,              0.,    -1./2.,       257./810.,         62./405.,    4./135.,
               0.,              0.,      -257./810.,        0.,      353./1728.,     5857./51840.,  1./51840.,
               0.,       -62./405.,     -353./1728.,        0.,   29489./90720.,   14293./453600.,  9./11200.,
         -4./135.,   -5857./51840.,  -29489./90720.,        0.,     1861./4725.,     4853./67200.,    1./630.,
       -1./51840., -14293./453600.,    -1861./4725.,        0.,   19093./67200.,       389./3150.,    4./225.,
       -9./11200.,   -4853./67200.,  -19093./67200.,        0.,       121./525.,         19./150.,         0.,
         -1./630.,     -389./3150.,      -121./525.,        0.,         16./45.,               0.,         0.,
         -4./225.,       -19./150.,        -16./45.,     1./2.,              0.,               0.,         0.
    };
    CHECK_GBMATRIX_CLOSE(
                8,       8,        3,        3, good_D_T1,       7,
            op.n(), op.n(), op.kl(1), op.ku(1), op.D_T(1), op.ld(),
            std::numeric_limits<double>::epsilon()*1000);

    /* Check D[2], the second derivative matrix, against known good
     * transposed in general banded matrix column-major order.
     */
    const double good_D_T2[] = {
               0.,           0.,           0.,         6./5.,   -274./135.,      88./135.,      8./45.,
               0.,           0.,    131./135.,      -19./15.,    -23./180.,    457./1080.,    1./1080.,
               0.,     88./135.,    -23./180.,       -10./9.,   473./1890.,   3061./9450.,     9./700.,
           8./45.,   457./1080.,   473./1890.,    -668./315.,      52./63.,       17./40.,     2./105.,
         1./1080.,  3061./9450.,      52./63.,    -836./525.,   -37./1400.,      38./105.,      8./75.,
          9./700.,      17./40.,   -37./1400.,        -6./7.,      -4./35.,       14./25.,          0.,
          2./105.,     38./105.,      -4./35.,        -7./5.,      17./15.,            0.,          0.,
           8./75.,      14./25.,     -28./15.,         6./5.,           0.,            0.,          0.
    };
    // Reduced precision as -37/1400 won't pass at the usual one
    CHECK_GBMATRIX_CLOSE(
                8,       8,        3,        3, good_D_T2,       7,
            op.n(), op.n(), op.kl(2), op.ku(2), op.D_T(2), op.ld(),
            std::numeric_limits<double>::epsilon()*10000);

    /* Check D[3], the third derivative matrix, against known good
     * transposed in general banded matrix column-major order.
     */
    const double good_D_T3[] = {
           0.,           0.,           0.,    -3./2.,      217./54.,     -92./27.,     8./9.,
           0.,           0.,     -91./54.,     9./2.,       -34./9.,      25./27.,    1./27.,
           0.,     -52./27.,       34./9.,        0.,    -844./189.,   2308./945.,    6./35.,
       -8./9.,     -25./27.,    844./189.,        0.,   -1408./315.,      57./35.,    4./21.,
      -1./27.,  -2308./945.,   1408./315.,        0.,     -103./35.,     44./105.,    8./15.,
      -6./35.,     -57./35.,     103./35.,        0.,     -103./35.,        9./5.,        0.,
      -4./21.,    -44./105.,     103./35.,    -9./2.,        13./6.,           0.,        0.,
      -8./15.,       11./5.,      -19./6.,     3./2.,            0.,           0.,        0.
    };
    CHECK_GBMATRIX_CLOSE(
                8,       8,        3,        3, good_D_T3,       7,
            op.n(), op.n(), op.kl(3), op.ku(3), op.D_T(3), op.ld(),
            std::numeric_limits<double>::epsilon()*1000);
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

BOOST_AUTO_TEST_CASE( derivatives_of_a_piecewise_cubic_representation )
{
    const double breakpts[] = { 0.0, 1.0, 2.0, 3.0 };
    bspline b(4, bspline::from_breakpoints(), COUNTOF(breakpts), breakpts);
    bsplineop op(b, 3, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE);

    poly_params *p = (poly_params *)
                      malloc(sizeof(poly_params) + 4*sizeof(double));
    p->n = 4;
    suzerain_function f = {poly_f, p};

    // Form and factorize the mass matrix M
    bsplineop_lu mass(op);
    mass.factor_mass(op);

    {
        const int derivative = 1;

        p->c[0] = 1.2; // Constant
        p->c[1] = 3.4; // Linear
        p->c[2] = 0.0; // Quadratic
        p->c[3] = 0.0; // Cubic

        // Compute the right hand side coefficients for M x = b
        double * coefficient = (double *) malloc(b.n() * sizeof(double));
        op.interpolation_rhs(&f, coefficient, b);

        // Solve for function coefficients using the mass matrix
        mass.solve(1, coefficient, 1, b.n());

        // Take the n-th derivative of the coefficients using M x' = D x
        op.apply(derivative, 1, 1.0, coefficient, 1, b.n());
        mass.solve(1, coefficient, 1, b.n());

        // Ensure we recover the leading order, scaled monomial coefficients
        for (int i = 0; i < b.n(); ++i) {
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
        double * coefficient = (double *) malloc(b.n() * sizeof(double));
        op.interpolation_rhs(&f, coefficient, b);

        // Solve for function coefficients using the mass matrix
        mass.solve(1, coefficient, 1, b.n());

        // Take the n-th derivative of the coefficients using M x' = D x
        op.apply(derivative, 1, 1.0, coefficient, 1, b.n());
        mass.solve(1, coefficient, 1, b.n());

        // Ensure we recover the leading order, scaled monomial coefficients
        for (int i = 0; i < b.n(); ++i) {
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
        double * coefficient = (double *) malloc(b.n() * sizeof(double));
        op.interpolation_rhs(&f, coefficient, b);

        // Solve for function coefficients using the mass matrix
        mass.solve(1, coefficient, 1, b.n());

        // Take the n-th derivative of the coefficients using M x' = D x
        op.apply(derivative, 1, 1.0, coefficient, 1, b.n());
        mass.solve(1, coefficient, 1, b.n());

        // Ensure we recover the leading order, scaled monomial coefficients
        for (int i = 0; i < b.n(); ++i) {
            BOOST_CHECK_CLOSE(6.0 * p->c[3], coefficient[i], 1e-11);
        }

        free(coefficient);
    }

    free(p);
}

// Intended to ensure our bandwidth routines are okay for
// high order bsplines and high order derivatives
BOOST_AUTO_TEST_CASE( ensure_creation_succeeds_at_high_order )
{
    const double breakpts[] = { 0.0, 1.0, 2.0, 3.0 };
    const int nbreak = COUNTOF(breakpts);

    const int maxk = 21;
    for (int k = 1; k <= maxk; ++k) {
        bspline b(k, bspline::from_breakpoints(), nbreak, breakpts);
        BOOST_TEST_MESSAGE("b: k = " << k);
        const int maxnderiv = maxk;
        for (int nderiv = 0; nderiv <= maxnderiv; ++nderiv) {
            bsplineop op1(b, nderiv, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE);
            BOOST_TEST_MESSAGE(
                    "\tcollocation: k = " << k << ", nderiv = " << nderiv);
            bsplineop op2(b, nderiv, SUZERAIN_BSPLINEOP_GALERKIN_L2);
            BOOST_TEST_MESSAGE(
                    "\tgalerkin:    k = " << k << ", nderiv = " << nderiv);
        }
    }
}

BOOST_AUTO_TEST_CASE( collocation_point_evaluation_is_operator_application )
{
    const int nbreak = 35;
    const int order  = 7;
    const int ndof   = nbreak + order - 2;
    const int nderiv = order - 1;
    double breakpts[nbreak];
    double points[ndof];
    double coefficients[ndof];
    double values_eval[ndof];
    double values_apply[ndof];

    // Establish workspace and get collocation point information
    suzerain::math::logspace(0.1, 3.0, nbreak, breakpts);
    bspline b(order, bspline::from_breakpoints(), nbreak, breakpts);
    bsplineop op(b, nderiv, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE);
    BOOST_REQUIRE_EQUAL(ndof, b.n());
    for (int i = 0; i < ndof; ++i) points[i] = b.collocation_point(i);

    // Generate "random" coefficients
    srand(634949092u);
    for (int i = 0; i < ndof; ++i) {
        coefficients[i] = 100.0 * (rand() / (RAND_MAX + 1.0)) - 50.0;
    }

    for (int k = 0; k <= nderiv; ++k) {
        BOOST_TEST_MESSAGE("Testing nderiv = " << k);

        // Evaluate coefficients times derivatives of basis functions
        BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, b.linear_combination(
                k, coefficients, ndof, points, values_eval, 0));

        // Apply derivative matrices directly
        std::memcpy(values_apply, coefficients, ndof*sizeof(double));
        BOOST_REQUIRE_EQUAL(SUZERAIN_SUCCESS, op.apply(
                k, 1, 1.0, values_apply, 1, ndof));

        // Check that both mechanisms give the same result
        check_close_collections(values_eval, values_eval + ndof,
                                values_apply, values_apply + ndof,
                                1.0e-12 * pow(10.0, k));
    }
}

static
void real_polynomial_interpolation(const int k,
                                   const int nderiv,
                                   const int nbreak,
                                   const double *breakpts,
                                   suzerain_bsplineop_method method,
                                   const double tol)
{
    // Prepare a descriptive header to bracket results
    std::ostringstream msg;
    msg << __FUNCTION__ << "(k=" << k << ",b=[";
    std::copy(breakpts, breakpts + nbreak,
                std::ostream_iterator<double>(msg, " "));
    msg << "],m=" << method << ",tol=" << tol << ")";

    BOOST_TEST_MESSAGE("Entering parameters " << msg.str());

    // Initialize B-spline basis, operators, and factored mass matrix
    bspline b(k, bspline::from_breakpoints(), nbreak, breakpts);
    bsplineop op(b, k, method);
    bsplineop_lu mass(op);
    mass.factor_mass(op);
    const int ndof = b.n();

    // Initialize polynomial test function which we should recapture exactly
    // Polynomial concocted to be reasonably behaved on domains of interest.
    suzerain::shared_ptr<poly_params> p(
        (poly_params *) malloc(sizeof(poly_params) + b.k()*sizeof(double)),
        free);
    p->n = b.k();
    for (int i = 0; i < p->n; ++i) {
        p->c[i] = pow(-1.0, (double) i) * (i + 10) / pow(4.0, (double) i + 1);
    }
    suzerain_function f = {poly_f, p.get()};

    // Compute expected coefficients for derivatives [0...nderiv] inclusive
    // by directly differentiating the polynomial test function.
    const std::size_t nstorage = (nderiv + 1) * ndof;
    suzerain::scoped_array<double> expected(new double[nstorage]);
    for (int i = 0; i <= nderiv; ++i) {
        op.interpolation_rhs(&f, &expected[i*ndof], b);
        poly_params_differentiate(p.get());
    }
    mass.solve(nderiv+1, expected.get(), 1, mass.n());

    // Test in-place application and solution
    {
        suzerain::scoped_array<double> result(new double[nstorage]);

        // Make multiple copies of the zeroth derivative coefficients
        for (int i = 0; i <= nderiv; ++i) {
            std::copy(&expected[0], &expected[ndof], &result[i*ndof]);
        }

        // Solve M*x' = D*x ...
        // ...starting by applying the derivative operators
        for (int i = 0; i <= nderiv; ++i) {
            op.apply(i, 1, 1.0, &result[i*ndof], 1, ndof);
        }
        // ...finishing by solving with the mass matrix
        mass.solve(nderiv+1, result.get(), 1, ndof);

        // See if we got anywhere close
        for (int i = 0; i <= nderiv; ++i) {
            BOOST_TEST_MESSAGE("\tApply-and-solve nderiv=" << i);
            check_close_collections(
                    &expected[i*ndof], &expected[(i+1)*ndof],
                    &result[i*ndof],   &result[(i+1)*ndof], tol);
        }
    }

    // Test out-of-place application and solution
    {
        suzerain::scoped_array<double> result(new double[nstorage]);

        // Solve M*x' = D*x ...
        // ...starting by accumulating the derivative operators
        for (int i = 0; i <= nderiv; ++i) {
            op.accumulate(i, 1, 1.0, &expected[0],    1, ndof,
                                0.0, &result[i*ndof], 1, ndof);
        }
        // ...finishing by solving with the mass matrix
        mass.solve(nderiv+1, result.get(), 1, ndof);

        // See if we got anywhere close
        for (int i = 0; i <= nderiv; ++i) {
            BOOST_TEST_MESSAGE("\tAccumulate-and-solve nderiv=" << i);
            check_close_collections(
                    &expected[i*ndof], &expected[(i+1)*ndof],
                    &result[i*ndof],   &result[(i+1)*ndof], tol);
        }
    }

    BOOST_TEST_MESSAGE("Leaving parameters " << msg.str());
}

BOOST_AUTO_TEST_CASE( compute_derivatives_of_a_general_polynomial )
{
    // Comparatively loose tolerance required on these tests
    // Indicative of a less-than-optimal implementation
    const double tol = std::sqrt(std::numeric_limits<double>::epsilon());

    BOOST_TEST_MESSAGE("Uniform breakpoints");
    {
        const double breakpts[] = { 0.0, 1.0, 2.0, 3.0, 4.0 };
        const int maxnderiv = 3;

        for (int k = 5; k < 9; ++k) {
            real_polynomial_interpolation(k, std::min(k-1, maxnderiv),
                    COUNTOF(breakpts), breakpts,
                    SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE, tol * k);
        }

        for (int k = 5; k < 9; ++k) {
            real_polynomial_interpolation(k, std::min(k-1, maxnderiv),
                    COUNTOF(breakpts), breakpts,
                    SUZERAIN_BSPLINEOP_GALERKIN_L2, 3 * tol * k * k);
        }
    }

    BOOST_TEST_MESSAGE("Spectral-like single interval with high order");
    {
        const double breakpts[] = { 3.0, 4.0 };
        const int maxnderiv = 2;

        for (int k = 11; k < 13; ++k) {
            real_polynomial_interpolation(k, std::min(k-1, maxnderiv),
                    COUNTOF(breakpts), breakpts,
                    SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE, tol * k);
        }

        for (int k = 11; k < 13; ++k) {
            real_polynomial_interpolation(k, std::min(k-1, maxnderiv),
                    COUNTOF(breakpts), breakpts,
                    SUZERAIN_BSPLINEOP_GALERKIN_L2, 7 * tol * k * k);
        }
    }
}

struct RealMassFixture {

    static const double breakpts[10]; // Init just below
    bspline      b;
    bsplineop    op;
    bsplineop_lu lu;
    bsplineop_lu t;

    RealMassFixture()
        : b(7, bspline::from_breakpoints(), COUNTOF(breakpts), breakpts),
          op(b, 2, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE),
          lu(op),
          t(op)
    { lu.factor_mass(op); }

};

const double RealMassFixture::breakpts[10] = {
    0.0, 0.5, 1.0, 3.0, 3.3, 7.0, 12.0, 12.5, 14.0, 15.0
};

// Tests for consistency among ways we might acquire just a mass matrix
BOOST_FIXTURE_TEST_SUITE( opaccumulate_real, RealMassFixture )

static const double zero = 0.0;
static const double half = 0.5;
static const double one  = 1.0;

BOOST_AUTO_TEST_CASE( form_mass )
{
    t.opform_mass(op);
    t.factor();

    CHECK_GBMATRIX_CLOSE(
        lu.n(), lu.n(), lu.kl(), lu.ku()+lu.kl(), lu.A_T(), lu.ld(),
         t.n(),  t.n(),  t.kl(),  t.ku()+ t.kl(),  t.A_T(),  t.ld(),
        std::numeric_limits<double>::epsilon()*1000);
}

BOOST_AUTO_TEST_CASE( form )
{
    t.opform(1, &one, op);
    t.factor();

    CHECK_GBMATRIX_CLOSE(
        lu.n(), lu.n(), lu.kl(), lu.ku()+lu.kl(), lu.A_T(), lu.ld(),
         t.n(),  t.n(),  t.kl(),  t.ku()+ t.kl(),  t.A_T(),  t.ld(),
        std::numeric_limits<double>::epsilon()*1000);
}

BOOST_AUTO_TEST_CASE( opaccumulate_simple )
{
    t.opaccumulate(1, &one, op, zero);
    t.factor();

    CHECK_GBMATRIX_CLOSE(
        lu.n(), lu.n(), lu.kl(), lu.ku()+lu.kl(), lu.A_T(), lu.ld(),
         t.n(),  t.n(),  t.kl(),  t.ku()+ t.kl(),  t.A_T(),  t.ld(),
        std::numeric_limits<double>::epsilon()*1000);
}

BOOST_AUTO_TEST_CASE( opaccumulate_simple_twice )
{
    t.opaccumulate(1, &half, op, zero);
    t.opaccumulate(1, &half, op, one);
    t.factor();

    CHECK_GBMATRIX_CLOSE(
        lu.n(), lu.n(), lu.kl(), lu.ku()+lu.kl(), lu.A_T(), lu.ld(),
         t.n(),  t.n(),  t.kl(),  t.ku()+ t.kl(),  t.A_T(),  t.ld(),
        std::numeric_limits<double>::epsilon()*1000);
}

BOOST_AUTO_TEST_CASE( opaccumulate_scale )
{
    t.opaccumulate(1, &half, op, zero);
    t.opaccumulate(0, NULL,  op, 2.0); // Scaling operation
    t.factor();

    CHECK_GBMATRIX_CLOSE(
        lu.n(), lu.n(), lu.kl(), lu.ku()+lu.kl(), lu.A_T(), lu.ld(),
         t.n(),  t.n(),  t.kl(),  t.ku()+ t.kl(),  t.A_T(),  t.ld(),
        std::numeric_limits<double>::epsilon()*1000);
}

BOOST_AUTO_TEST_CASE( opaccumulate_repeated )
{
    // Accumulate our way to a mass matrix
    { double c[3] = { 5,  6,  4 }; t.opaccumulate(3, c, op, one ); }
    { double c[3] = { 0, -2, 14 }; t.opaccumulate(3, c, op, zero); } // 0!
    { double c[3] = { 2,  0,  0 }; t.opaccumulate(3, c, op, one ); }
    { double c[3] = { 0,  1, -7 }; t.opaccumulate(3, c, op, half); }
    t.factor();

    CHECK_GBMATRIX_CLOSE(
        lu.n(), lu.n(), lu.kl(), lu.ku()+lu.kl(), lu.A_T(), lu.ld(),
         t.n(),  t.n(),  t.kl(),  t.ku()+ t.kl(),  t.A_T(),  t.ld(),
        std::numeric_limits<double>::epsilon()*1000000);

    // Accumulate our way to a mass matrix one more time
    { double c[3] = { 0,  4, -8 }; t.opaccumulate(3, c, op, zero); } // 0!
    { double c[3] = { 2,  0,  0 }; t.opaccumulate(3, c, op, one ); }
    { double c[3] = { 0, -2,  4 }; t.opaccumulate(3, c, op, half); }
    t.factor();

    CHECK_GBMATRIX_CLOSE(
        lu.n(), lu.n(), lu.kl(), lu.ku()+lu.kl(), lu.A_T(), lu.ld(),
         t.n(),  t.n(),  t.kl(),  t.ku()+ t.kl(),  t.A_T(),  t.ld(),
        std::numeric_limits<double>::epsilon()*1000000);
}

BOOST_AUTO_TEST_SUITE_END()

struct ComplexMassFixture {

    static const double breakpts[10]; // Init just below
    bspline       b;
    bsplineop     op;
    bsplineop_luz luz;
    bsplineop_luz t;

    ComplexMassFixture()
        : b(7, bspline::from_breakpoints(), COUNTOF(breakpts), breakpts),
          op(b, 2, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE),
          luz(op),
          t(op)
    { luz.factor_mass(op); }

};

const double ComplexMassFixture::breakpts[10] = {
    0.0, 0.5, 1.0, 3.0, 3.3, 7.0, 12.0, 12.5, 14.0, 15.0
};

// Tests for consistency among ways we might acquire just a mass matrix
BOOST_FIXTURE_TEST_SUITE( opaccumulate_complex, ComplexMassFixture )

static const double zero[2] = { 0.0, 0.0 };
static const double half[2] = { 0.5, 0.0 };
static const double one[2]  = { 1.0, 0.0 };

BOOST_AUTO_TEST_CASE( form_mass )
{
    t.opform_mass(op);
    t.factor();

    CHECK_GBMATRIX_CLOSE(
        luz.n(), luz.n(), luz.kl(), luz.ku()+luz.kl(), luz.A_T(), luz.ld(),
         t.n(),    t.n(),   t.kl(),   t.ku()+  t.kl(),   t.A_T(),   t.ld(),
        std::numeric_limits<double>::epsilon()*1000);
}

BOOST_AUTO_TEST_CASE( form )
{
    t.opform(1, &one, op);
    t.factor();

    CHECK_GBMATRIX_CLOSE(
        luz.n(), luz.n(), luz.kl(), luz.ku()+luz.kl(), luz.A_T(), luz.ld(),
         t.n(),    t.n(),   t.kl(),   t.ku()+  t.kl(),   t.A_T(),   t.ld(),
        std::numeric_limits<double>::epsilon()*1000);
}

BOOST_AUTO_TEST_CASE( opaccumulate_simple )
{
    t.opaccumulate(1, &one, op, zero);
    t.factor();

    CHECK_GBMATRIX_CLOSE(
        luz.n(), luz.n(), luz.kl(), luz.ku()+luz.kl(), luz.A_T(), luz.ld(),
         t.n(),    t.n(),   t.kl(),   t.ku()+  t.kl(),   t.A_T(),   t.ld(),
        std::numeric_limits<double>::epsilon()*1000);
}

BOOST_AUTO_TEST_CASE( opaccumulate_simple_twice )
{
    t.opaccumulate(1, &half, op, zero);
    t.opaccumulate(1, &half, op, one);
    t.factor();

    CHECK_GBMATRIX_CLOSE(
        luz.n(), luz.n(), luz.kl(), luz.ku()+luz.kl(), luz.A_T(), luz.ld(),
         t.n(),    t.n(),   t.kl(),   t.ku()+  t.kl(),   t.A_T(),   t.ld(),
        std::numeric_limits<double>::epsilon()*1000);
}

BOOST_AUTO_TEST_CASE( opaccumulate_scale )
{
    const double coeff[2] = { 0.08, 0.04 };
    const double scale[2] = { 10.0, -5.0 };
    const double (*my_null)[2] = NULL; // Proper type required

    t.opaccumulate(1, &coeff,  op, zero);
    t.opaccumulate(0, my_null, op, scale); // Scaling operation
    t.factor();

    CHECK_GBMATRIX_CLOSE(
        luz.n(), luz.n(), luz.kl(), luz.ku()+luz.kl(), luz.A_T(), luz.ld(),
         t.n(),    t.n(),   t.kl(),   t.ku()+  t.kl(),   t.A_T(),   t.ld(),
        std::numeric_limits<double>::epsilon()*1000);
}

BOOST_AUTO_TEST_CASE( opaccumulate_repeated )
{
    // Accumulate our way to a mass matrix
    { double c[3][2] = { {5,0}, { 6,0}, { 4,0} }; t.opaccumulate(3, c, op, one ); }
    { double c[3][2] = { {0,0}, {-2,0}, {14,0} }; t.opaccumulate(3, c, op, zero); } // 0!
    { double c[3][2] = { {2,0}, { 0,0}, { 0,0} }; t.opaccumulate(3, c, op, one ); }
    { double c[3][2] = { {0,0}, { 1,0}, {-7,0} }; t.opaccumulate(3, c, op, half); }
    t.factor();

    CHECK_GBMATRIX_CLOSE(
        luz.n(), luz.n(), luz.kl(), luz.ku()+luz.kl(), luz.A_T(), luz.ld(),
         t.n(),    t.n(),   t.kl(),   t.ku()+  t.kl(),   t.A_T(),   t.ld(),
        std::numeric_limits<double>::epsilon()*1000000);

    // Accumulate our way to a mass matrix a second time using imag higher derivs
    { double c[3][2] = { {0,0}, {0,-2}, {0,14} }; t.opaccumulate(3, c, op, zero); } // 0!
    { double c[3][2] = { {2,0}, {0, 0}, {0, 0} }; t.opaccumulate(3, c, op, one ); }
    { double c[3][2] = { {0,0}, {0, 1}, {0,-7} }; t.opaccumulate(3, c, op, half); }
    t.factor();

    CHECK_GBMATRIX_CLOSE(
        luz.n(), luz.n(), luz.kl(), luz.ku()+luz.kl(), luz.A_T(), luz.ld(),
         t.n(),    t.n(),   t.kl(),   t.ku()+  t.kl(),   t.A_T(),   t.ld(),
        std::numeric_limits<double>::epsilon()*1000000);

}

BOOST_AUTO_TEST_SUITE_END()
