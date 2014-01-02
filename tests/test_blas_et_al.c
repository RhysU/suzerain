//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Suzerain is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------

#include <suzerain/blas_et_al.h>

#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_test.h>

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#ifdef HAVE_MKL
#include <mkl.h>
#endif

#ifdef SUZERAIN_HAVE_MKL
#include <mkl_service.h>
#endif

#include <suzerain/common.h>

static
void
test_daxpby()
{
    int i;

    const double alpha      = 2.0;
    const double x[]        = {1.0, 2.0, 3.0};
    const int    incx       = 1;
    const double beta       = 3.0;
    double       y[]        = {4.0, -1, 5.0, -2, 6.0, -3};
    const int    incy       = 2;
    const int    nx         = sizeof(x)/sizeof(x[0]);
    const int    ny         = sizeof(y)/sizeof(y[0]);
    const double expected[] = {
        alpha*x[0] + beta*y[0],
        y[1],
        alpha*x[1] + beta*y[2],
        y[3],
        alpha*x[2] + beta*y[4],
        y[5]
    };
    const int nexpected = sizeof(expected)/sizeof(expected[0]);

    gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
    gsl_test_int(ny, nexpected, "Expected results' length");

    suzerain_blas_daxpby(nx/incx, alpha, x, incx, beta, y, incy);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(y[i], expected[i], GSL_DBL_EPSILON,
                "%s:%d index %d", __func__, __LINE__, i);
    }
}

static
void
test_saxpby()
{
    int i;

    const float  alpha      = 2.0;
    const float  x[]        = {1.0, 2.0, 3.0};
    const int    incx       = 1;
    const float  beta       = 3.0;
    float        y[]        = {4.0, -1, 5.0, -2, 6.0, -3};
    const int    incy       = 2;
    const int    nx         = sizeof(x)/sizeof(x[0]);
    const int    ny         = sizeof(y)/sizeof(y[0]);
    const float expected[] = {
        alpha*x[0] + beta*y[0],
        y[1],
        alpha*x[1] + beta*y[2],
        y[3],
        alpha*x[2] + beta*y[4],
        y[5]
    };
    const int nexpected = sizeof(expected)/sizeof(expected[0]);

    gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
    gsl_test_int(ny, nexpected, "Expected results' length");

    suzerain_blas_saxpby(nx/incx, alpha, x, incx, beta, y, incy);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(y[i], expected[i], GSL_FLT_EPSILON,
                "%s:%d index %d", __func__, __LINE__, i);
    }
}

static
void
test_caxpby()
{
    int i;

    const float alpha_re = 2.;
    const float alpha_im = -0.25;
    const float x[][2]   = {{1., -1.}, {2., -2.}, {3., -3.}};
    const int   incx     = 1;
    const float beta_re  = 3.;
    const float beta_im  = -0.50;
    float       y[][2]   = {{4.,-4.}, {-1,1},{5.,-5},{-2,2},{6.,-6},{-3,3}};
    const int   incy     = 2;
    const int   nx       = sizeof(x)/sizeof(x[0]);
    const int   ny       = sizeof(y)/sizeof(y[0]);
    const float  expected[][2] = {
        {
            alpha_re*x[0][0]-alpha_im*x[0][1]+beta_re*y[0][0]-beta_im*y[0][1],
            alpha_re*x[0][1]+alpha_im*x[0][0]+beta_re*y[0][1]+beta_im*y[0][0]
        },
        {
            y[1][0],
            y[1][1]
        },
        {
            alpha_re*x[1][0]-alpha_im*x[1][1]+beta_re*y[2][0]-beta_im*y[2][1],
            alpha_re*x[1][1]+alpha_im*x[1][0]+beta_re*y[2][1]+beta_im*y[2][0]
        },
        {
            y[3][0],
            y[3][1]
        },
        {
            alpha_re*x[2][0]-alpha_im*x[2][1]+beta_re*y[4][0]-beta_im*y[4][1],
            alpha_re*x[2][1]+alpha_im*x[2][0]+beta_re*y[4][1]+beta_im*y[4][0]
        },
        {
            y[5][0],
            y[5][1]
        }
    };
    const int nexpected = sizeof(expected)/sizeof(expected[0]);

    gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
    gsl_test_int(ny, nexpected, "Expected results' length");

    suzerain_blas_caxpby(nx/incx,
            complex_float(alpha_re, alpha_im), (complex_float *) x, incx,
            complex_float(beta_re,  beta_im),  (complex_float *) y, incy);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(y[i][0], expected[i][0], GSL_FLT_EPSILON,
                "%s:%d real index %d", __func__, __LINE__, i);
        gsl_test_abs(y[i][1], expected[i][1], GSL_FLT_EPSILON,
                "%s:%d imag index %d", __func__, __LINE__, i);
    }
}

static
void
test_zaxpby()
{
    int i;

    const double alpha_re = 2.;
    const double alpha_im = -0.25;
    const double x[][2]   = {{1., -1.}, {2., -2.}, {3., -3.}};
    const int    incx     = 1;
    const double beta_re  = 3.;
    const double beta_im  = -0.50;
    double       y[][2]   = {{4.,-4.},{-1,1},{5.,-5},{-2,2},{6.,-6},{-3,3}};
    const int    incy     = 2;
    const int    nx       = sizeof(x)/sizeof(x[0]);
    const int    ny       = sizeof(y)/sizeof(y[0]);
    const double expected[][2] = {
        {
            alpha_re*x[0][0]-alpha_im*x[0][1]+beta_re*y[0][0]-beta_im*y[0][1],
            alpha_re*x[0][1]+alpha_im*x[0][0]+beta_re*y[0][1]+beta_im*y[0][0]
        },
        {
            y[1][0],
            y[1][1]
        },
        {
            alpha_re*x[1][0]-alpha_im*x[1][1]+beta_re*y[2][0]-beta_im*y[2][1],
            alpha_re*x[1][1]+alpha_im*x[1][0]+beta_re*y[2][1]+beta_im*y[2][0]
        },
        {
            y[3][0],
            y[3][1]
        },
        {
            alpha_re*x[2][0]-alpha_im*x[2][1]+beta_re*y[4][0]-beta_im*y[4][1],
            alpha_re*x[2][1]+alpha_im*x[2][0]+beta_re*y[4][1]+beta_im*y[4][0]
        },
        {
            y[5][0],
            y[5][1]
        }
    };
    const int nexpected = sizeof(expected)/sizeof(expected[0]);

    gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
    gsl_test_int(ny, nexpected, "Expected results' length");

    suzerain_blas_zaxpby(nx/incx,
            complex_double(alpha_re, alpha_im), (complex_double *) x, incx,
            complex_double(beta_re,  beta_im),  (complex_double *) y, incy);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(y[i][0], expected[i][0], GSL_DBL_EPSILON,
                "%s:%d real index %d", __func__, __LINE__, i);
        gsl_test_abs(y[i][1], expected[i][1], GSL_DBL_EPSILON,
                "%s:%d imag index %d", __func__, __LINE__, i);
    }
}

static
void
test_dwaxpby()
{
    int i;

    const double alpha      = 2.0;
    const double x[]        = {1.0, 2.0, 3.0};
    const int    incx       = 1;
    const int    nx         = sizeof(x)/sizeof(x[0]);
    const double beta       = 3.0;
    const double y[]        = {4.0, -1, 5.0, -2, 6.0, -3};
    const int    incy       = 2;
    const int    ny         = sizeof(y)/sizeof(y[0]);
    double       w[nx/incx];
    const int    nw         = sizeof(w)/sizeof(w[0]);
    const int    incw       = 1;
    const double expected[] = {
        alpha*x[0] + beta*y[0],
        alpha*x[1] + beta*y[2],
        alpha*x[2] + beta*y[4],
    };
    const int nexpected = sizeof(expected)/sizeof(expected[0]);

    gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
    gsl_test_int(nw, nexpected, "Expected results' length");

    suzerain_blas_dwaxpby(nx/incx, alpha, x, incx, beta, y, incy, w, incw);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(w[i], expected[i], GSL_DBL_EPSILON,
                "%s:%d index %d", __func__, __LINE__, i);
    }
}

static
void
test_swaxpby()
{
    int i;

    const float  alpha      = 2.0;
    const float  x[]        = {1.0, 2.0, 3.0};
    const int    incx       = 1;
    const int    nx         = sizeof(x)/sizeof(x[0]);
    const float  beta       = 3.0;
    const float  y[]        = {4.0, -1, 5.0, -2, 6.0, -3};
    const int    incy       = 2;
    const int    ny         = sizeof(y)/sizeof(y[0]);
    float        w[nx/incx];
    const int    nw         = sizeof(w)/sizeof(w[0]);
    const int    incw       = 1;
    const float  expected[] = {
        alpha*x[0] + beta*y[0],
        alpha*x[1] + beta*y[2],
        alpha*x[2] + beta*y[4],
    };
    const int nexpected = sizeof(expected)/sizeof(expected[0]);

    gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
    gsl_test_int(nw, nexpected, "Expected results' length");

    suzerain_blas_swaxpby(nx/incx, alpha, x, incx, beta, y, incy, w, incw);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(w[i], expected[i], GSL_FLT_EPSILON,
                "%s:%d index %d", __func__, __LINE__, i);
    }
}

static
void
test_daxpby_nop()
{
    int i;

    const double alpha     = 0.0;  /* NOP zero */
    const double x[]       = {1.0, 2.0, 3.0};
    const int    incx      = 1;
    const double beta      = 1.0;  /* NOP one */
    double       y[]       = {4.0, -1, 5.0, -2, 6.0, -3};
    const int    incy      = 2;
    const int    nx        = sizeof(x)/sizeof(x[0]);
    const int    ny        = sizeof(y)/sizeof(y[0]);
    const double *expected = y;
    const int    nexpected = ny;

    gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
    gsl_test_int(ny, nexpected, "Expected results' length");

    suzerain_blas_daxpby(nx/incx, alpha, x, incx, beta, y, incy);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(y[i], expected[i], GSL_DBL_EPSILON,
                "%s:%d index %d", __func__, __LINE__, i);
    }
}

static
void
test_saxpby_nop()
{
    int i;

    const float  alpha     = 0.0;  /* NOP zero */
    const float  x[]       = {1.0, 2.0, 3.0};
    const int    incx      = 1;
    const float  beta      = 1.0;  /* NOP one */
    float        y[]       = {4.0, -1, 5.0, -2, 6.0, -3};
    const int    incy      = 2;
    const int    nx        = sizeof(x)/sizeof(x[0]);
    const int    ny        = sizeof(y)/sizeof(y[0]);
    const float *expected = y;
    const int    nexpected = ny;

    gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
    gsl_test_int(ny, nexpected, "Expected results' length");

    suzerain_blas_saxpby(nx/incx, alpha, x, incx, beta, y, incy);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(y[i], expected[i], GSL_FLT_EPSILON,
                "%s:%d index %d", __func__, __LINE__, i);
    }
}

static
void
test_dswap()
{
    int i;

    double       x[]    = {1.0, 2.0, 3.0};
    double       y[]    = {4.0, 5.0, 6.0};
    const int    incx   = 1;
    const int    incy   = 1;
    const int    nx     = sizeof(x)/sizeof(x[0]);
    const int    ny     = sizeof(y)/sizeof(y[0]);
    const double x_expected[] = { 4.0, 5.0, 6.0 };
    const double y_expected[] = { 1.0, 2.0, 3.0 };
    const int nx_expected = sizeof(x_expected)/sizeof(x_expected[0]);
    const int ny_expected = sizeof(y_expected)/sizeof(y_expected[0]);

    gsl_test_int(nx/incx, ny/incy, "Consistent vector lengths");
    gsl_test_int(nx, nx_expected, "Expected x results' length");
    gsl_test_int(ny, ny_expected, "Expected y results' length");

    suzerain_blas_dswap(nx/incx, x, incx, y, incy);
    for (i = 0; i < nx_expected; ++i) {
        gsl_test_abs(x[i], x_expected[i], GSL_DBL_EPSILON,
                "%s:%d x index %d", __func__, __LINE__, i);
    }
    for (i = 0; i < ny_expected; ++i) {
        gsl_test_abs(y[i], y_expected[i], GSL_DBL_EPSILON,
                "%s:%d y index %d", __func__, __LINE__, i);
    }
}

static
void
test_sswap()
{
    int i;

    float       x[]          = {1.0, 2.0, 3.0};
    float       y[]          = {4.0, 5.0, 6.0};
    const int   incx         = 1;
    const int   incy         = 1;
    const int   nx           = sizeof(x)/sizeof(x[0]);
    const int   ny           = sizeof(y)/sizeof(y[0]);
    const float x_expected[] = { 4.0, 5.0, 6.0 };
    const float y_expected[] = { 1.0, 2.0, 3.0 };
    const int nx_expected = sizeof(x_expected)/sizeof(x_expected[0]);
    const int ny_expected = sizeof(y_expected)/sizeof(y_expected[0]);

    gsl_test_int(nx/incx, ny/incy, "Consistent vector lengths");
    gsl_test_int(nx, nx_expected, "Expected x results' length");
    gsl_test_int(ny, ny_expected, "Expected y results' length");

    suzerain_blas_sswap(nx/incx, x, incx, y, incy);
    for (i = 0; i < nx_expected; ++i) {
        gsl_test_abs(x[i], x_expected[i], GSL_DBL_EPSILON,
                "%s:%d x index %d", __func__, __LINE__, i);
    }
    for (i = 0; i < ny_expected; ++i) {
        gsl_test_abs(y[i], y_expected[i], GSL_DBL_EPSILON,
                "%s:%d y index %d", __func__, __LINE__, i);
    }
}

static
void
test_cswap()
{
    int i;

    float        x[][2]    = {{1.0, -1.0}, {2.0, -2.0}, {3.0, -3.0}};
    float        y[][2]    = {{4.0, -4.0}, {5.0, -5.0}, {6.0, -6.0}};
    const int    incx   = 1;
    const int    incy   = 1;
    const int    nx     = sizeof(x)/sizeof(x[0]);
    const int    ny     = sizeof(y)/sizeof(y[0]);
    const float  x_expected[][2] = { {4.0, -4.0}, {5.0, -5.0}, {6.0, -6.0} };
    const float  y_expected[][2] = { {1.0, -1.0}, {2.0, -2.0}, {3.0, -3.0} };
    const int nx_expected = sizeof(x_expected)/sizeof(x_expected[0]);
    const int ny_expected = sizeof(y_expected)/sizeof(y_expected[0]);

    gsl_test_int(nx/incx, ny/incy, "Consistent vector lengths");
    gsl_test_int(nx, nx_expected, "Expected x results' length");
    gsl_test_int(ny, ny_expected, "Expected y results' length");

    suzerain_blas_cswap(nx/incx, (complex_float *) x, incx,
                                 (complex_float *) y, incy);
    for (i = 0; i < nx_expected; ++i) {
        gsl_test_abs(x[i][0], x_expected[i][0], GSL_FLT_EPSILON,
                "%s:%d x real index %d", __func__, __LINE__, i);
        gsl_test_abs(x[i][1], x_expected[i][1], GSL_FLT_EPSILON,
                "%s:%d x imag index %d", __func__, __LINE__, i);
    }
    for (i = 0; i < ny_expected; ++i) {
        gsl_test_abs(y[i][0], y_expected[i][0], GSL_FLT_EPSILON,
                "%s:%d y real index %d", __func__, __LINE__, i);
        gsl_test_abs(y[i][1], y_expected[i][1], GSL_FLT_EPSILON,
                "%s:%d y imag index %d", __func__, __LINE__, i);
    }
}

static
void
test_zswap()
{
    int i;

    double       x[][2]    = {{1.0, -1.0}, {2.0, -2.0}, {3.0, -3.0}};
    double       y[][2]    = {{4.0, -4.0}, {5.0, -5.0}, {6.0, -6.0}};
    const int    incx   = 1;
    const int    incy   = 1;
    const int    nx     = sizeof(x)/sizeof(x[0]);
    const int    ny     = sizeof(y)/sizeof(y[0]);
    const double x_expected[][2] = { {4.0, -4.0}, {5.0, -5.0}, {6.0, -6.0} };
    const double y_expected[][2] = { {1.0, -1.0}, {2.0, -2.0}, {3.0, -3.0} };
    const int nx_expected = sizeof(x_expected)/sizeof(x_expected[0]);
    const int ny_expected = sizeof(y_expected)/sizeof(y_expected[0]);

    gsl_test_int(nx/incx, ny/incy, "Consistent vector lengths");
    gsl_test_int(nx, nx_expected, "Expected x results' length");
    gsl_test_int(ny, ny_expected, "Expected y results' length");

    suzerain_blas_zswap(nx/incx, (complex_double*) x, incx,
                                 (complex_double*) y, incy);
    for (i = 0; i < nx_expected; ++i) {
        gsl_test_abs(x[i][0], x_expected[i][0], GSL_DBL_EPSILON,
                "%s:%d x real index %d", __func__, __LINE__, i);
        gsl_test_abs(x[i][1], x_expected[i][1], GSL_DBL_EPSILON,
                "%s:%d x imag index %d", __func__, __LINE__, i);
    }
    for (i = 0; i < ny_expected; ++i) {
        gsl_test_abs(y[i][0], y_expected[i][0], GSL_DBL_EPSILON,
                "%s:%d y real index %d", __func__, __LINE__, i);
        gsl_test_abs(y[i][1], y_expected[i][1], GSL_DBL_EPSILON,
                "%s:%d y imag index %d", __func__, __LINE__, i);
    }
}

static
void
test_dscal()
{
    int i;

    {
        const double alpha  = 2.0;
        double       x[]    = {1.0, 2.0, 3.0};
        const int    incx   = 1;
        const int    nx     = sizeof(x)/sizeof(x[0]);
        const double x_expected[] = { alpha*x[0], alpha*x[1], alpha*x[2] };
        const int nx_expected = sizeof(x_expected)/sizeof(x_expected[0]);

        gsl_test_int(nx, nx_expected, "Expected results' length");
        suzerain_blas_dscal(nx/incx, alpha, x, incx);
        for (i = 0; i < nx_expected; ++i) {
            gsl_test_abs(x[i], x_expected[i], GSL_DBL_EPSILON,
                    "%s:%d index %d", __func__, __LINE__, i);
        }
    }

    {
        const double beta   = 3.0;
        double       y[]    = {4.0, -1, 5.0, -2, 6.0, -3};
        const int    incy   = 2;
        const int    ny     = sizeof(y)/sizeof(y[0]);
        const double y_expected[] = {
            beta*y[0], y[1], beta*y[2], y[3], beta*y[4], y[5] };
        const int ny_expected = sizeof(y_expected)/sizeof(y_expected[0]);

        gsl_test_int(ny, ny_expected, "Expected results' length");
        suzerain_blas_dscal(ny/incy, beta, y, incy);
        for (i = 0; i < ny_expected; ++i) {
            gsl_test_abs(y[i], y_expected[i], GSL_DBL_EPSILON,
                    "%s:%d index %d", __func__, __LINE__, i);
        }
    }
}

static
void
test_sscal()
{
    int i;

    {
        const float  alpha  = 2.0;
        float        x[]    = {1.0, 2.0, 3.0};
        const int    incx   = 1;
        const int    nx     = sizeof(x)/sizeof(x[0]);
        const float  x_expected[] = { alpha*x[0], alpha*x[1], alpha*x[2] };
        const int nx_expected = sizeof(x_expected)/sizeof(x_expected[0]);

        gsl_test_int(nx, nx_expected, "Expected results' length");
        suzerain_blas_sscal(nx/incx, alpha, x, incx);
        for (i = 0; i < nx_expected; ++i) {
            gsl_test_abs(x[i], x_expected[i], GSL_FLT_EPSILON,
                    "%s:%d index %d", __func__, __LINE__, i);
        }
    }

    {
        const float  beta   = 3.0;
        float        y[]    = {4.0, -1, 5.0, -2, 6.0, -3};
        const int    incy   = 2;
        const int    ny     = sizeof(y)/sizeof(y[0]);
        const float  y_expected[] = {
            beta*y[0], y[1], beta*y[2], y[3], beta*y[4], y[5] };
        const int ny_expected = sizeof(y_expected)/sizeof(y_expected[0]);

        gsl_test_int(ny, ny_expected, "Expected results' length");
        suzerain_blas_sscal(ny/incy, beta, y, incy);
        for (i = 0; i < ny_expected; ++i) {
            gsl_test_abs(y[i], y_expected[i], GSL_FLT_EPSILON,
                    "%s:%d index %d", __func__, __LINE__, i);
        }
    }
}

static
void
test_zscal()
{
    int i;

    const double alpha[2]      = {2.0, -0.25};
    double       x[][2]        = {{4.0,-4.0}, {1,-1},
                                  {5.0,-5.0}, {2,-2},
                                  {6.0,-6.0}, {3,-3}};
    const int    incx          = 2;
    const int    nx            = sizeof(x)/sizeof(x[0]);
    const double expected[][2] = {
        {
            alpha[0]*x[0][0] - alpha[1]*x[0][1],
            alpha[0]*x[0][1] + alpha[1]*x[0][0]
        },
        {
            x[1][0],
            x[1][1]
        },
        {
            alpha[0]*x[2][0] - alpha[1]*x[2][1],
            alpha[0]*x[2][1] + alpha[1]*x[2][0]
        },
        {
            x[3][0],
            x[3][1]
        },
        {
            alpha[0]*x[4][0] - alpha[1]*x[4][1],
            alpha[0]*x[4][1] + alpha[1]*x[4][0]
        },
        {
            x[5][0],
            x[5][1]
        }
    };
    const int nexpected = sizeof(expected)/sizeof(expected[0]);

    gsl_test_int(nx, nexpected, "Expected results' length");

    suzerain_blas_zscal(nx/incx, complex_double(alpha[0], alpha[1]),
                        (complex_double *) x, incx);

    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(x[i][0], expected[i][0], GSL_DBL_EPSILON,
                "%s:%d real index %d", __func__, __LINE__, i);
        gsl_test_abs(x[i][1], expected[i][1], GSL_DBL_EPSILON,
                "%s:%d imag index %d", __func__, __LINE__, i);
    }
}

static
void
test_cscal()
{
    int i;

    const float  alpha[2]      = {2.0, -0.25};
    float        x[][2]        = {{4.0,-4.0}, {1,-1},
                                  {5.0,-5.0}, {2,-2},
                                  {6.0,-6.0}, {3,-3}};
    const int    incx          = 2;
    const int    nx            = sizeof(x)/sizeof(x[0]);
    const float  expected[][2] = {
        {
            alpha[0]*x[0][0] - alpha[1]*x[0][1],
            alpha[0]*x[0][1] + alpha[1]*x[0][0]
        },
        {
            x[1][0],
            x[1][1]
        },
        {
            alpha[0]*x[2][0] - alpha[1]*x[2][1],
            alpha[0]*x[2][1] + alpha[1]*x[2][0]
        },
        {
            x[3][0],
            x[3][1]
        },
        {
            alpha[0]*x[4][0] - alpha[1]*x[4][1],
            alpha[0]*x[4][1] + alpha[1]*x[4][0]
        },
        {
            x[5][0],
            x[5][1]
        }
    };
    const int nexpected = sizeof(expected)/sizeof(expected[0]);

    gsl_test_int(nx, nexpected, "Expected results' length");

    suzerain_blas_cscal(nx/incx, complex_float(alpha[0], alpha[1]),
                        (complex_float *) x, incx);

    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(x[i][0], expected[i][0], GSL_FLT_EPSILON,
                "%s:%d real index %d", __func__, __LINE__, i);
        gsl_test_abs(x[i][1], expected[i][1], GSL_FLT_EPSILON,
                "%s:%d imag index %d", __func__, __LINE__, i);
    }
}

static
void
test_dcopy()
{
    int i;

    const double x[]    = {1.0, 2.0, 3.0};
    double       y[]    = {555.0, 555.0, 555.0};
    const int    incx   = 1;
    const int    incy   = 1;
    const int    nx     = sizeof(x)/sizeof(x[0]);
    const int    ny     = sizeof(y)/sizeof(y[0]);
    const double y_expected[] = { 1.0, 2.0, 3.0 };
    const int ny_expected = sizeof(y_expected)/sizeof(y_expected[0]);

    gsl_test_int(nx/incx, ny/incy, "Consistent vector lengths");
    gsl_test_int(ny, ny_expected, "Expected y results' length");

    suzerain_blas_dcopy(nx/incx, x, incx, y, incy);
    for (i = 0; i < ny_expected; ++i) {
        gsl_test_abs(y[i], y_expected[i], GSL_DBL_EPSILON,
                "%s:%d y index %d", __func__, __LINE__, i);
    }
}

static
void
test_dcopy_fill()
{
    int i;

    const double x      = 3;
    double       y[]    = {555.0, 555.0, 555.0};
    const int    incy   = 1;
    const int    ny     = sizeof(y)/sizeof(y[0]);
    const double y_expected[] = { 3.0, 3.0, 3.0 };
    const int ny_expected = sizeof(y_expected)/sizeof(y_expected[0]);

    gsl_test_int(ny, ny_expected, "Expected y results' length");

    suzerain_blas_dcopy(ny/incy, &x, 0 /* zero! */, y, incy);
    for (i = 0; i < ny_expected; ++i) {
        gsl_test_abs(y[i], y_expected[i], GSL_DBL_EPSILON,
                "%s:%d y index %d", __func__, __LINE__, i);
    }
}

static
void
test_scopy()
{
    int i;

    const float  x[]    = {1.0, 2.0, 3.0};
    float        y[]    = {555.0, 555.0, 555.0};
    const int    incx   = 1;
    const int    incy   = 1;
    const int    nx     = sizeof(x)/sizeof(x[0]);
    const int    ny     = sizeof(y)/sizeof(y[0]);
    const float  y_expected[] = { 1.0, 2.0, 3.0 };
    const int ny_expected = sizeof(y_expected)/sizeof(y_expected[0]);

    gsl_test_int(nx/incx, ny/incy, "Consistent vector lengths");
    gsl_test_int(ny, ny_expected, "Expected y results' length");

    suzerain_blas_scopy(nx/incx, x, incx, y, incy);
    for (i = 0; i < ny_expected; ++i) {
        gsl_test_abs(y[i], y_expected[i], GSL_DBL_EPSILON,
                "%s:%d y index %d", __func__, __LINE__, i);
    }
}

static
void
test_zcopy()
{
    int i;

    const double x[][2]        = {{1.0,-1.0}, {2.0,-2.0}, {3.0,-3.0}};
    const int    incx          = 1;
    double       y[][2]        = {{555.0,555.0}, {777.0,777.0},
                                  {555.0,555.0}, {777.0,777.0},
                                  {555.0,555.0}, {777.0,777.0}};
    const int incy             = 2;
    const int nx               = sizeof(x)/sizeof(x[0]);
    const int ny               = sizeof(y)/sizeof(y[0]);
    const double expected[][2] = {{  1.0, -1.0}, {777.0,777.0},
                                  {  2.0, -2.0}, {777.0,777.0},
                                  {  3.0, -3.0}, {777.0,777.0}};
    const int nexpected = sizeof(expected)/sizeof(expected[0]);

    gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
    gsl_test_int(ny, nexpected, "Expected results' length");

    suzerain_blas_zcopy(nx/incx, (complex_double *) x, incx,
                                 (complex_double *) y, incy);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(y[i][0], expected[i][0], GSL_DBL_EPSILON,
                "%s:%d real index %d", __func__, __LINE__, i);
        gsl_test_abs(y[i][1], expected[i][1], GSL_DBL_EPSILON,
                "%s:%d imag index %d", __func__, __LINE__, i);
    }
}

static
void
test_ccopy()
{
    int i;

    const float  x[][2]        = {{1.0,-1.0}, {2.0,-2.0}, {3.0,-3.0}};
    const int    incx          = 1;
    float        y[][2]        = {{555.0,555.0}, {777.0,777.0},
                                  {555.0,555.0}, {777.0,777.0},
                                  {555.0,555.0}, {777.0,777.0}};
    const int incy             = 2;
    const int nx               = sizeof(x)/sizeof(x[0]);
    const int ny               = sizeof(y)/sizeof(y[0]);
    const float  expected[][2] = {{  1.0, -1.0}, {777.0,777.0},
                                  {  2.0, -2.0}, {777.0,777.0},
                                  {  3.0, -3.0}, {777.0,777.0}};
    const int nexpected = sizeof(expected)/sizeof(expected[0]);

    gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
    gsl_test_int(ny, nexpected, "Expected results' length");

    suzerain_blas_ccopy(nx/incx, (complex_float *) x, incx,
                                 (complex_float *) y, incy);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(y[i][0], expected[i][0], GSL_FLT_EPSILON,
                "%s:%d real index %d", __func__, __LINE__, i);
        gsl_test_abs(y[i][1], expected[i][1], GSL_FLT_EPSILON,
                "%s:%d imag index %d", __func__, __LINE__, i);
    }
}

static
void
test_snrm2()
{
    const float  x[]      = {3.0, 555.0, 4.0, 555};
    const int    incx     = 2;
    const int    nx       = sizeof(x)/sizeof(x[0])/incx;
    const float  expected = 5;

    gsl_test_abs(suzerain_blas_snrm2(nx, x, incx), expected, GSL_FLT_EPSILON,
            "%s:%d result", __func__);
}

static
void
test_dnrm2()
{
    const double x[]      = {3.0, 555.0, 4.0, 555};
    const int    incx     = 2;
    const int    nx       = sizeof(x)/sizeof(x[0])/incx;
    const double expected = 5;

    gsl_test_abs(suzerain_blas_dnrm2(nx, x, incx), expected, GSL_DBL_EPSILON,
            "%s:%d result", __func__);
}

static
void
test_scnrm2()
{
    const float  x[][2]   = {{3.0,4.0},{555,555},{5.0,6.0},{555,555}};
    const int    incx     = 2;
    const int    nx       = sizeof(x)/sizeof(x[0])/incx;
    const float  expected = 9.27361849549570375251641607399; /* Sqrt 86 */

    gsl_test_abs(suzerain_blas_scnrm2(nx, (complex_float *) x, incx),
                 expected, GSL_FLT_EPSILON, "%s:%d result", __func__);
}

static
void
test_dznrm2()
{
    const double x[][2]   = {{3.0,4.0},{555,555},{5.0,6.0},{555,555}};
    const int    incx     = 2;
    const int    nx       = sizeof(x)/sizeof(x[0])/incx;
    const double expected = 9.27361849549570375251641607399; /* Sqrt 86 */

    gsl_test_abs(suzerain_blas_dznrm2(nx, (complex_double *) x, incx),
                 expected, GSL_FLT_EPSILON, "%s:%d result", __func__);
}

static
void
test_daxpy()
{
    int i;

    const double alpha      = 2.0;
    const double x[]        = {1.0, 2.0, 3.0};
    const int    incx       = 1;
    double       y[]        = {4.0, -1, 5.0, -2, 6.0, -3};
    const int    incy       = 2;
    const int    nx         = sizeof(x)/sizeof(x[0]);
    const int    ny         = sizeof(y)/sizeof(y[0]);
    const double expected[] = {
        alpha*x[0] + y[0],
        y[1],
        alpha*x[1] + y[2],
        y[3],
        alpha*x[2] + y[4],
        y[5]
    };
    const int nexpected = sizeof(expected)/sizeof(expected[0]);

    gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
    gsl_test_int(ny, nexpected, "Expected results' length");

    suzerain_blas_daxpy(nx/incx, alpha, x, incx, y, incy);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(y[i], expected[i], GSL_DBL_EPSILON,
                "%s:%d index %d", __func__, __LINE__, i);
    }
}

static
void
test_daxpy_const()
{
    int i;

    const double alpha      = 1.0;
    const double x          = 1.0;
    double       y[]        = {4.0, -1, 5.0, -2, 6.0, -3};
    const int    incy       = 2;
    const int    ny         = sizeof(y)/sizeof(y[0]);
    const double expected[] = {
        alpha*x + y[0],
        y[1],
        alpha*x + y[2],
        y[3],
        alpha*x + y[4],
        y[5]
    };
    const int nexpected = sizeof(expected)/sizeof(expected[0]);

    gsl_test_int(ny, nexpected, "Expected results' length");

    suzerain_blas_daxpy(ny/incy, alpha, &x, 0 /* zero! */, y, incy);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(y[i], expected[i], GSL_DBL_EPSILON,
                "%s:%d index %d", __func__, __LINE__, i);
    }
}

static
void
test_saxpy()
{
    int i;

    const float  alpha      = 2.0;
    const float  x[]        = {1.0, 2.0, 3.0};
    const int    incx       = 1;
    float        y[]        = {4.0, -1, 5.0, -2, 6.0, -3};
    const int    incy       = 2;
    const int    nx         = sizeof(x)/sizeof(x[0]);
    const int    ny         = sizeof(y)/sizeof(y[0]);
    const float expected[] = {
        alpha*x[0] + y[0],
        y[1],
        alpha*x[1] + y[2],
        y[3],
        alpha*x[2] + y[4],
        y[5]
    };
    const int nexpected = sizeof(expected)/sizeof(expected[0]);

    gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
    gsl_test_int(ny, nexpected, "Expected results' length");

    suzerain_blas_saxpy(nx/incx, alpha, x, incx, y, incy);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(y[i], expected[i], GSL_FLT_EPSILON,
                "%s:%d index %d", __func__, __LINE__, i);
    }
}

static
void
test_zaxpy()
{
    int i;

    const double alpha[2] = {2.0, -0.25};
    const double x[][2]   = {{1.0,-1.0}, {2.0,-2.0}, {3.0,-3.0}};
    const int    incx     = 1;
    double       y[][2]
        = {{4.0,-4.0}, {-1,-1}, {5.0,-5.0}, {-2,-2}, {6.0,-6.0}, {-3,-3}};
    const int incy        = 2;
    const int nx          = sizeof(x)/sizeof(x[0]);
    const int ny          = sizeof(y)/sizeof(y[0]);
    const double expected[][2] = {
        {
            alpha[0]*x[0][0] - alpha[1]*x[0][1] + y[0][0],
            alpha[0]*x[0][1] + alpha[1]*x[0][0] + y[0][1]
        },
        {
            y[1][0],
            y[1][1]
        },
        {
            alpha[0]*x[1][0] - alpha[1]*x[1][1] + y[2][0],
            alpha[0]*x[1][1] + alpha[1]*x[1][0] + y[2][1]
        },
        {
            y[3][0],
            y[3][1]
        },
        {
            alpha[0]*x[2][0] - alpha[1]*x[2][1] + y[4][0],
            alpha[0]*x[2][1] + alpha[1]*x[2][0] + y[4][1]
        },
        {
            y[5][0],
            y[5][1]
        }
    };
    const int nexpected = sizeof(expected)/sizeof(expected[0]);

    gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
    gsl_test_int(ny, nexpected, "Expected results' length");

    suzerain_blas_zaxpy(nx/incx, complex_double(alpha[0],alpha[1]),
                        (complex_double *) x, incx,
                        (complex_double *) y, incy);

    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(y[i][0], expected[i][0], GSL_DBL_EPSILON,
                "%s:%d real index %d", __func__, __LINE__, i);
        gsl_test_abs(y[i][1], expected[i][1], GSL_DBL_EPSILON,
                "%s:%d imag index %d", __func__, __LINE__, i);
    }
}

static
void
test_caxpy()
{
    int i;

    const float  alpha[2] = {2.0, -0.25};
    const float  x[][2]   = {{1.0,-1.0}, {2.0,-2.0}, {3.0,-3.0}};
    const int    incx     = 1;
    float        y[][2]
        = {{4.0,-4.0}, {-1,-1}, {5.0,-5.0}, {-2,-2}, {6.0,-6.0}, {-3,-3}};
    const int incy        = 2;
    const int nx          = sizeof(x)/sizeof(x[0]);
    const int ny          = sizeof(y)/sizeof(y[0]);
    const float  expected[][2] = {
        {
            alpha[0]*x[0][0] - alpha[1]*x[0][1] + y[0][0],
            alpha[0]*x[0][1] + alpha[1]*x[0][0] + y[0][1]
        },
        {
            y[1][0],
            y[1][1]
        },
        {
            alpha[0]*x[1][0] - alpha[1]*x[1][1] + y[2][0],
            alpha[0]*x[1][1] + alpha[1]*x[1][0] + y[2][1]
        },
        {
            y[3][0],
            y[3][1]
        },
        {
            alpha[0]*x[2][0] - alpha[1]*x[2][1] + y[4][0],
            alpha[0]*x[2][1] + alpha[1]*x[2][0] + y[4][1]
        },
        {
            y[5][0],
            y[5][1]
        }
    };
    const int nexpected = sizeof(expected)/sizeof(expected[0]);

    gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
    gsl_test_int(ny, nexpected, "Expected results' length");

    suzerain_blas_caxpy(nx/incx, complex_float(alpha[0],alpha[1]),
                        (complex_float *) x, incx,
                        (complex_float *) y, incy);

    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(y[i][0], expected[i][0], GSL_FLT_EPSILON,
                "%s:%d real index %d", __func__, __LINE__, i);
        gsl_test_abs(y[i][1], expected[i][1], GSL_FLT_EPSILON,
                "%s:%d imag index %d", __func__, __LINE__, i);
    }
}

static
void
test_ddot()
{
    const double x[]        = {1.0, 2.0, 3.0};
    const int    incx       = 1;
    double       y[]        = {4.0, -1, 5.0, -2, 6.0, -3};
    const int    incy       = 2;
    const int    nx         = sizeof(x)/sizeof(x[0]);
    const int    ny         = sizeof(y)/sizeof(y[0]);
    const double expected   = (x[0]*y[0]) + (x[1]*y[2]) + (x[2]*y[4]);

    gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");

    const double result = suzerain_blas_ddot(nx/incx, x, incx, y, incy);
    gsl_test_abs(result, expected, GSL_DBL_EPSILON, "%s:%d result", __func__);
}

static
void
test_sdot()
{
    const float  x[]        = {1.0, 2.0, 3.0};
    const int    incx       = 1;
    float        y[]        = {4.0, -1, 5.0, -2, 6.0, -3};
    const int    incy       = 2;
    const int    nx         = sizeof(x)/sizeof(x[0]);
    const int    ny         = sizeof(y)/sizeof(y[0]);
    const float  expected   = (x[0]*y[0]) + (x[1]*y[2]) + (x[2]*y[4]);

    gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");

    const float  result = suzerain_blas_sdot(nx/incx, x, incx, y, incy);
    gsl_test_abs(result, expected, GSL_DBL_EPSILON, "%s:%d result", __func__);
}

static
void
test_zdotc()
{
    const double x[][2]      = {{1.0, 2.0}, {2.0, 3.0}, {3.0, 4.0}};
    const int    incx        = 1;
    const double y[][2]      = { {2.0, 3.0}, {-5,-5},
                                 {3.0, 4.0}, {-5,-5},
                                 {4.0, 5.0}, {-5,-5} };
    const int    incy        = 2;
    const int    nx          = sizeof(x)/sizeof(x[0]);
    const int    ny          = sizeof(y)/sizeof(y[0]);
    const double expected[2] = { 58.0, -3.0 };

    gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");

    double result[2] = {555, 555};
    suzerain_blas_zdotc(nx/incx, (complex_double *) x, incx,
                                 (complex_double *) y, incy,
                                 (complex_double *) result);
    gsl_test_rel(result[0], expected[0], GSL_DBL_EPSILON,
            "%s:%d real result", __func__);
    gsl_test_rel(result[1], expected[1], GSL_DBL_EPSILON,
            "%s:%d imag result", __func__);
}

static
void
test_cdotc()
{
    const float  x[][2]      = {{1.0, 2.0}, {2.0, 3.0}, {3.0, 4.0}};
    const int    incx        = 1;
    const float  y[][2]      = { {2.0, 3.0}, {-5,-5},
                                 {3.0, 4.0}, {-5,-5},
                                 {4.0, 5.0}, {-5,-5} };
    const int    incy        = 2;
    const int    nx          = sizeof(x)/sizeof(x[0]);
    const int    ny          = sizeof(y)/sizeof(y[0]);
    const float  expected[2] = { 58.0, -3.0 };

    gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");

    float  result[2] = {555, 555};
    suzerain_blas_cdotc(nx/incx, (complex_float *) x, incx,
                                 (complex_float *) y, incy,
                                 (complex_float *) result);
    gsl_test_rel(result[0], expected[0], GSL_FLT_EPSILON,
            "%s:%d real result", __func__);
    gsl_test_rel(result[1], expected[1], GSL_FLT_EPSILON,
            "%s:%d imag result", __func__);
}

static
void
test_dgb_acc1()
{
    int i;

    const int    m     = 7;
    const int    n     = 7;
    const int    ku    = 2;
    const int    kl    = 1;
    const int    lda   = 5;
    const int    ldb   = 6;
    const double alpha = 2.0;
    const double beta  = 3.0;

    /* Negative one represents values outside of the band */
    /* Decimal parts flag regions outside the matrix entirely */
    const double a_data[]  = {
        /*lda buffer*/   /*ku2*/  /*ku1*/ /*diag*/   /*kl1*/
                 -1.1,     -1.0,    -1.0,     1.0,      2.0,
                 -1.1,     -1.0,     3.0,     4.0,      5.0,
                 -1.1,      6.0,     7.0,     8.0,      9.0,
                 -1.1,     10.0,    11.0,    12.0,      1.0,
                 -1.1,      2.0,     3.0,     4.0,      5.0,
                 -1.1,      6.0,     7.0,     8.0,      9.0,
                 -1.1,     10.0,    11.0,    12.0,     -1.0
    };
    double b_data[]  = {
        /*ldb buffer*/   /*ku2*/  /*ku1*/ /*diag*/   /*kl1*/
         -2.2,   -1.1,     -1.0,    -1.0,     1.0,      2.0,
         -2.2,   -1.1,     -1.0,     3.0,     4.0,      5.0,
         -2.2,   -1.1,      6.0,     7.0,     8.0,      9.0,
         -2.2,   -1.1,     10.0,    11.0,    12.0,      1.0,
         -2.2,   -1.1,      2.0,     3.0,     4.0,      5.0,
         -2.2,   -1.1,      6.0,     7.0,     8.0,      9.0,
         -2.2,   -1.1,     10.0,    11.0,    12.0,     -1.0
    };
    const double expected_data[]  = {
        /*ldb buffer*/   /*ku2*/  /*ku1*/ /*diag*/   /*kl1*/
         -2.2,   -1.1,     -1.0,    -1.0,     5.0,     10.0,
         -2.2,   -1.1,     -1.0,    15.0,    20.0,     25.0,
         -2.2,   -1.1,     30.0,    35.0,    40.0,     45.0,
         -2.2,   -1.1,     50.0,    55.0,    60.0,      5.0,
         -2.2,   -1.1,     10.0,    15.0,    20.0,     25.0,
         -2.2,   -1.1,     30.0,    35.0,    40.0,     45.0,
         -2.2,   -1.1,     50.0,    55.0,    60.0,     -1.0
    };
    const double *a = a_data + lda-(ku+1+kl);
    double       *b = b_data + ldb-(ku+1+kl);

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = sizeof(expected_data)/sizeof(expected_data[0]);
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blas_dgb_acc( m, n, kl, ku, alpha, a, lda, beta, b, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i], expected_data[i], GSL_DBL_EPSILON,
                "%s:%d index %d", __func__, __LINE__, i);
    }
}

static
void
test_dgb_acc2()
{
    int i;

    const int    m     = 2;
    const int    n     = 2;
    const int    ku    = 2;
    const int    kl    = 1;
    const int    lda   = ku + 1 + kl;
    const int    ldb   = ku + 1 + kl;
    const double alpha = 2.0;
    const double beta  = 3.0;

    /* Negative one represents values outside of the band */
    const double a_data[]  = {
        /*ku2*/  /*ku1*/ /*diag*/   /*kl1*/
          -1.0,    -1.0,     1.0,      2.0,
          -1.0,     3.0,     4.0,     -5.0
    };
    double b_data[]  = {
        /*ku2*/  /*ku1*/ /*diag*/   /*kl1*/
          -1.0,    -1.0,     1.0,      2.0,
          -1.0,     3.0,     4.0,     -5.0
    };
    const double expected_data[]  = {
        /*ku2*/  /*ku1*/ /*diag*/   /*kl1*/
          -1.0,    -1.0,     5.0,     10.0,
          -1.0,    15.0,    20.0,     -5.0
    };
    const double *a = a_data + lda-(ku+1+kl);
    double       *b = b_data + ldb-(ku+1+kl);

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = sizeof(expected_data)/sizeof(expected_data[0]);
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blas_dgb_acc( m, n, kl, ku, alpha, a, lda, beta, b, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i], expected_data[i], GSL_DBL_EPSILON,
                "%s:%d index %d", __func__, __LINE__, i);
    }
}

static
void
test_sgb_acc1()
{
    int i;

    const int    m     = 4;
    const int    n     = 4;
    const int    ku    = 2;
    const int    kl    = 1;
    const int    lda   = 5;
    const int    ldb   = 6;
    const float  alpha = 2.0;
    const float  beta  = 3.0;

    /* Negative one represents values outside of the band */
    /* Decimal parts flag regions outside the matrix entirely */
    const float  a_data[]  = {
        /*lda buffer*/   /*ku2*/  /*ku1*/ /*diag*/   /*kl1*/
                 -1.1,     -1.0,    -1.0,     1.0,      2.0,
                 -1.1,     -1.0,     3.0,     4.0,      5.0,
                 -1.1,      6.0,     7.0,     8.0,      9.0,
                 -1.1,     10.0,    11.0,    12.0,     -1.0
    };
    float  b_data[]  = {
        /*ldb buffer*/   /*ku2*/  /*ku1*/ /*diag*/   /*kl1*/
         -2.2,   -1.1,     -1.0,    -1.0,     1.0,      2.0,
         -2.2,   -1.1,     -1.0,     3.0,     4.0,      5.0,
         -2.2,   -1.1,      6.0,     7.0,     8.0,      9.0,
         -2.2,   -1.1,     10.0,    11.0,    12.0,     -1.0
    };
    const float  expected_data[]  = {
        /*ldb buffer*/   /*ku2*/  /*ku1*/ /*diag*/   /*kl1*/
         -2.2,   -1.1,     -1.0,    -1.0,     5.0,     10.0,
         -2.2,   -1.1,     -1.0,    15.0,    20.0,     25.0,
         -2.2,   -1.1,     30.0,    35.0,    40.0,     45.0,
         -2.2,   -1.1,     50.0,    55.0,    60.0,     -1.0
    };
    const float  *a = a_data + lda-(ku+1+kl);
    float        *b = b_data + ldb-(ku+1+kl);

    const float  nb        = sizeof(b_data)/sizeof(b_data[0]);
    const float  nexpected = sizeof(expected_data)/sizeof(expected_data[0]);
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blas_sgb_acc( m, n, kl, ku, alpha, a, lda, beta, b, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i], expected_data[i], GSL_FLT_EPSILON,
                "%s:%d index %d", __func__, __LINE__, i);
    }
}

static
void
test_sgb_acc2()
{
    int i;

    const int    m     = 5;
    const int    n     = 5;
    const int    ku    = 2;
    const int    kl    = 0;
    const int    lda   = ku + 1 + kl;
    const int    ldb   = ku + 1 + kl;
    const float  alpha = 2.0;
    const float  beta  = 3.0;

    /* Negative one represents values outside of the band */
    const float  a_data[] = {
        /*ku2*/  /*ku1*/ /*diag*/
          -1.0,    -1.0,     1.0,
          -1.0,     3.0,     4.0,
           6.0,     7.0,     8.0,
          10.0,    11.0,    12.0,
          10.0,    11.0,    12.0
    };
    float  b_data[] = {
        /*ku2*/  /*ku1*/ /*diag*/
          -1.0,    -1.0,     1.0,
          -1.0,     3.0,     4.0,
           6.0,     7.0,     8.0,
          10.0,    11.0,    12.0,
          10.0,    11.0,    12.0
    };
    const float  expected_data[] = {
        /*ku2*/  /*ku1*/ /*diag*/
          -1.0,    -1.0,     5.0,
          -1.0,    15.0,    20.0,
          30.0,    35.0,    40.0,
          50.0,    55.0,    60.0,
          50.0,    55.0,    60.0
    };
    const float  *a = a_data + lda-(ku+1+kl);
    float        *b = b_data + ldb-(ku+1+kl);

    const float  nb        = sizeof(b_data)/sizeof(b_data[0]);
    const float  nexpected = sizeof(expected_data)/sizeof(expected_data[0]);
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blas_sgb_acc( m, n, kl, ku, alpha, a, lda, beta, b, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i], expected_data[i], GSL_FLT_EPSILON,
                "%s:%d index %d", __func__, __LINE__, i);
    }
}

static
void
test_cgb_acc1()
{
    int i;

    const int    m        = 4;
    const int    n        = 4;
    const int    ku       = 2;
    const int    kl       = 1;
    const int    lda      = 5;
    const int    ldb      = 6;
    const complex_float alpha = 2.0; /* TODO Nontrivial imaginary part */
    const complex_float beta  = 3.0; /* TODO Nontrivial imaginary part */

    /* Negative one represents values outside of the band */
    /* Decimal parts flag regions outside the matrix entirely */
    const float  a_data[][2]  = {
        /*lda buffer*/ /*ku2*/  /*ku1*/  /*diag*/ /*kl1*/
        {-1.1,-1.1},   {-1,-1}, {-1,-1}, { 1, 1}, { 2, 2},
        {-1.1,-1.1},   {-1,-1}, { 3, 3}, { 4, 4}, { 5, 5},
        {-1.1,-1.1},   { 6, 6}, { 7, 7}, { 8, 8}, { 9, 9},
        {-1.1,-1.1},   {10,10}, {11,11}, {12,12}, {-1,-1}
    };
    float  b_data[][2]  = {
        /*ldb buffer*/            /*ku2*/  /*ku1*/  /*diag*/ /*kl1*/
        {-2.2,-2.2}, {-1.1,-1.1}, {-1,-1}, {-1,-1}, { 1, 1}, { 2, 2},
        {-2.2,-2.2}, {-1.1,-1.1}, {-1,-1}, { 3, 3}, { 4, 4}, { 5, 5},
        {-2.2,-2.2}, {-1.1,-1.1}, { 6, 6}, { 7, 7}, { 8, 8}, { 9, 9},
        {-2.2,-2.2}, {-1.1,-1.1}, {10,10}, {11,11}, {12,12}, {-1,-1}
    };
    const float  expected_data[][2]  = {
        /*ldb buffer*/            /*ku2*/  /*ku1*/  /*diag*/ /*kl1*/
        {-2.2,-2.2}, {-1.1,-1.1}, {-1,-1}, {-1,-1}, { 5, 5}, {10,10},
        {-2.2,-2.2}, {-1.1,-1.1}, {-1,-1}, {15,15}, {20,20}, {25,25},
        {-2.2,-2.2}, {-1.1,-1.1}, {30,30}, {35,35}, {40,40}, {45,45},
        {-2.2,-2.2}, {-1.1,-1.1}, {50,50}, {55,55}, {60,60}, {-1,-1}
    };
    const float  (*a)[2] = a_data + lda-(ku+1+kl);
    float        (*b)[2] = b_data + ldb-(ku+1+kl);

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = sizeof(expected_data)/sizeof(expected_data[0]);
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blas_cgb_acc( m, n, kl, ku, alpha, (complex_float *) a, lda,
                                         beta,  (complex_float *) b, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i][0], expected_data[i][0], GSL_FLT_EPSILON,
                "%s:%d real index %d", __func__, __LINE__, i);
        gsl_test_abs(b_data[i][1], expected_data[i][1], GSL_FLT_EPSILON,
                "%s:%d imag index %d", __func__, __LINE__, i);
    }
}

static
void
test_cgb_acc2()
{
    int i;

    const int    m        = 7;
    const int    n        = 7;
    const int    ku       = 1;
    const int    kl       = 1;
    const int    lda      = ku + 1 + kl;
    const int    ldb      = ku + 1 + kl;
    const complex_float alpha = 2.0; /* TODO Nontrivial imaginary part */
    const complex_float beta  = 3.0; /* TODO Nontrivial imaginary part */

    /* Negative one represents values outside of the band */
    const float  a_data[][2]  = {
        /*ku1*/  /*diag*/ /*kl1*/
        {-1,-1}, { 1, 1}, { 2, 2},
        { 3, 3}, { 4, 4}, { 5, 5},
        { 7, 7}, { 8, 8}, { 9, 9},
        {11,11}, {12,12}, { 1, 1},
        { 3, 3}, { 4, 4}, { 5, 5},
        { 7, 7}, { 8, 8}, { 9, 9},
        {11,11}, {12,12}, {-1,-1}
    };
    float  b_data[][2]  = {
        /*ku1*/  /*diag*/ /*kl1*/
        {-1,-1}, { 1, 1}, { 2, 2},
        { 3, 3}, { 4, 4}, { 5, 5},
        { 7, 7}, { 8, 8}, { 9, 9},
        {11,11}, {12,12}, { 1, 1},
        { 3, 3}, { 4, 4}, { 5, 5},
        { 7, 7}, { 8, 8}, { 9, 9},
        {11,11}, {12,12}, {-1,-1}
    };
    const float  expected_data[][2]  = {
        /*ku1*/  /*diag*/ /*kl1*/
        {-1,-1}, { 5, 5}, {10,10},
        {15,15}, {20,20}, {25,25},
        {35,35}, {40,40}, {45,45},
        {55,55}, {60,60}, { 5, 5},
        {15,15}, {20,20}, {25,25},
        {35,35}, {40,40}, {45,45},
        {55,55}, {60,60}, {-1,-1}
    };
    const float  (*a)[2] = a_data + lda-(ku+1+kl);
    float        (*b)[2] = b_data + ldb-(ku+1+kl);

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = sizeof(expected_data)/sizeof(expected_data[0]);
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blas_cgb_acc( m, n, kl, ku, alpha, (complex_float *) a, lda,
                                         beta,  (complex_float *) b, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i][0], expected_data[i][0], GSL_FLT_EPSILON,
                "%s:%d real index %d", __func__, __LINE__, i);
        gsl_test_abs(b_data[i][1], expected_data[i][1], GSL_FLT_EPSILON,
                "%s:%d imag index %d", __func__, __LINE__, i);
    }
}

static
void
test_zgb_acc1()
{
    int i;

    const int    m        = 4;
    const int    n        = 4;
    const int    ku       = 2;
    const int    kl       = 1;
    const int    lda      = 5;
    const int    ldb      = 6;
    const complex_double alpha = 2.0; /* TODO Nontrivial imaginary part */
    const complex_double beta  = 3.0; /* TODO Nontrivial imaginary part */

    /* Negative one represents values outside of the band */
    /* Decimal parts flag regions outside the matrix entirely */
    const double a_data[][2]  = {
        /*lda buffer*/ /*ku2*/  /*ku1*/  /*diag*/ /*kl1*/
        {-1.1,-1.1},   {-1,-1}, {-1,-1}, { 1, 1}, { 2, 2},
        {-1.1,-1.1},   {-1,-1}, { 3, 3}, { 4, 4}, { 5, 5},
        {-1.1,-1.1},   { 6, 6}, { 7, 7}, { 8, 8}, { 9, 9},
        {-1.1,-1.1},   {10,10}, {11,11}, {12,12}, {-1,-1}
    };
    double b_data[][2]  = {
        /*ldb buffer*/            /*ku2*/  /*ku1*/  /*diag*/ /*kl1*/
        {-2.2,-2.2}, {-1.1,-1.1}, {-1,-1}, {-1,-1}, { 1, 1}, { 2, 2},
        {-2.2,-2.2}, {-1.1,-1.1}, {-1,-1}, { 3, 3}, { 4, 4}, { 5, 5},
        {-2.2,-2.2}, {-1.1,-1.1}, { 6, 6}, { 7, 7}, { 8, 8}, { 9, 9},
        {-2.2,-2.2}, {-1.1,-1.1}, {10,10}, {11,11}, {12,12}, {-1,-1}
    };
    const double expected_data[][2]  = {
        /*ldb buffer*/            /*ku2*/  /*ku1*/  /*diag*/ /*kl1*/
        {-2.2,-2.2}, {-1.1,-1.1}, {-1,-1}, {-1,-1}, { 5, 5}, {10,10},
        {-2.2,-2.2}, {-1.1,-1.1}, {-1,-1}, {15,15}, {20,20}, {25,25},
        {-2.2,-2.2}, {-1.1,-1.1}, {30,30}, {35,35}, {40,40}, {45,45},
        {-2.2,-2.2}, {-1.1,-1.1}, {50,50}, {55,55}, {60,60}, {-1,-1}
    };
    const double (*a)[2] = a_data + lda-(ku+1+kl);
    double       (*b)[2] = b_data + ldb-(ku+1+kl);

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = sizeof(expected_data)/sizeof(expected_data[0]);
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blas_zgb_acc( m, n, kl, ku, alpha, (complex_double *) a, lda,
                                         beta,  (complex_double *) b, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i][0], expected_data[i][0], GSL_FLT_EPSILON,
                "%s:%d real index %d", __func__, __LINE__, i);
        gsl_test_abs(b_data[i][1], expected_data[i][1], GSL_FLT_EPSILON,
                "%s:%d imag index %d", __func__, __LINE__, i);
    }
}

static
void
test_zgb_acc2()
{
    int i;

    const int    m        = 3;
    const int    n        = 3;
    const int    ku       = 2;
    const int    kl       = 1;
    const int    lda      = ku + 1 + kl;
    const int    ldb      = ku + 1 + kl;
    const complex_double alpha = 2.0; /* TODO Nontrivial imaginary part */
    const complex_double beta  = 3.0; /* TODO Nontrivial imaginary part */

    /* Negative one represents values outside of the band */
    /* Decimal parts flag regions outside the matrix entirely */
    const double a_data[][2]  = {
        /*ku2*/  /*ku1*/  /*diag*/ /*kl1*/
        {-1,-1}, {-1,-1}, { 1, 1}, { 2, 2},
        {-1,-1}, { 3, 3}, { 4, 4}, { 5, 5},
        { 6, 6}, { 7, 7}, { 8, 8}, {-9,-9}
    };
    double b_data[][2]  = {
        /*ku2*/  /*ku1*/  /*diag*/ /*kl1*/
        {-1,-1}, {-1,-1}, { 1, 1}, { 2, 2},
        {-1,-1}, { 3, 3}, { 4, 4}, { 5, 5},
        { 6, 6}, { 7, 7}, { 8, 8}, {-9,-9}
    };
    const double expected_data[][2]  = {
        /*ku2*/  /*ku1*/  /*diag*/ /*kl1*/
        {-1,-1}, {-1,-1}, { 5, 5}, {10,10},
        {-1,-1}, {15,15}, {20,20}, {25,25},
        {30,30}, {35,35}, {40,40}, {-9,-9}
    };
    const double (*a)[2] = a_data + lda-(ku+1+kl);
    double       (*b)[2] = b_data + ldb-(ku+1+kl);

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = sizeof(expected_data)/sizeof(expected_data[0]);
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blas_zgb_acc( m, n, kl, ku, alpha, (complex_double *) a, lda,
                                         beta,  (complex_double *) b, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i][0], expected_data[i][0], GSL_FLT_EPSILON,
                "%s:%d real index %d", __func__, __LINE__, i);
        gsl_test_abs(b_data[i][1], expected_data[i][1], GSL_FLT_EPSILON,
                "%s:%d imag index %d", __func__, __LINE__, i);
    }
}

static
void
test_dgb_acc_nop()
{
    int i;

    const int    m     = 4;
    const int    n     = m;
    const int    ku    = 2;
    const int    kl    = 1;
    const int    lda   = 5;
    const int    ldb   = 6;
    const double alpha = 0.0; /* NOP zero */
    const double beta  = 1.0; /* NOP one */

    /* Negative one represents values outside of the band */
    /* Decimal parts flag regions outside the matrix entirely */
    const double a_data[]  = {
        /*lda buffer*/   /*ku2*/  /*ku1*/ /*diag*/   /*kl1*/
                 -1.1,     -1.0,    -1.0,     1.0,      2.0,
                 -1.1,     -1.0,     3.0,     4.0,      5.0,
                 -1.1,      6.0,     7.0,     8.0,      9.0,
                 -1.1,     10.0,    11.0,    12.0,     -1.0
    };
    double b_data[]  = {
        /*ldb buffer*/   /*ku2*/  /*ku1*/ /*diag*/   /*kl1*/
         -2.2,   -1.1,     -1.0,    -1.0,     1.0,      2.0,
         -2.2,   -1.1,     -1.0,     3.0,     4.0,      5.0,
         -2.2,   -1.1,      6.0,     7.0,     8.0,      9.0,
         -2.2,   -1.1,     10.0,    11.0,    12.0,     -1.0
    };
    const double *a = a_data + lda-(ku+1+kl);
    double       *b = b_data + ldb-(ku+1+kl);
    double       *expected_data;

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = nb;
    gsl_test_int(nb, nexpected, "Expected results' length");

    /* Make a clean copy of the original data */
    expected_data = malloc(nb*sizeof(b_data[0]));
    memcpy(expected_data, b_data, nb*sizeof(b_data[0]));

    suzerain_blas_dgb_acc( m, n, kl, ku, alpha, a, lda, beta, b, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i], expected_data[i], GSL_DBL_EPSILON,
                "%s:%d index %d", __func__, __LINE__, i);
    }

    free(expected_data);
}

static
void
test_sgb_acc_nop()
{
    int i;

    const int    m     = 4;
    const int    n     = m;
    const int    ku    = 2;
    const int    kl    = 1;
    const int    lda   = 5;
    const int    ldb   = 6;
    const float alpha  = 0.0; /* NOP zero */
    const float beta   = 1.0; /* NOP one */

    /* Negative one represents values outside of the band */
    /* Decimal parts flag regions outside the matrix entirely */
    const float a_data[]  = {
        /*lda buffer*/   /*ku2*/  /*ku1*/ /*diag*/   /*kl1*/
                 -1.1,     -1.0,    -1.0,     1.0,      2.0,
                 -1.1,     -1.0,     3.0,     4.0,      5.0,
                 -1.1,      6.0,     7.0,     8.0,      9.0,
                 -1.1,     10.0,    11.0,    12.0,     -1.0
    };
    float b_data[]  = {
        /*ldb buffer*/   /*ku2*/  /*ku1*/ /*diag*/   /*kl1*/
         -2.2,   -1.1,     -1.0,    -1.0,     1.0,      2.0,
         -2.2,   -1.1,     -1.0,     3.0,     4.0,      5.0,
         -2.2,   -1.1,      6.0,     7.0,     8.0,      9.0,
         -2.2,   -1.1,     10.0,    11.0,    12.0,     -1.0
    };
    const float *a = a_data + lda-(ku+1+kl);
    float       *b = b_data + ldb-(ku+1+kl);
    float       *expected_data;

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = nb;
    gsl_test_int(nb, nexpected, "Expected results' length");

    /* Make a clean copy of the original data */
    expected_data = malloc(nb*sizeof(b_data[0]));
    memcpy(expected_data, b_data, nb*sizeof(b_data[0]));

    suzerain_blas_sgb_acc( m, n, kl, ku, alpha, a, lda, beta, b, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i], expected_data[i], GSL_DBL_EPSILON,
                "%s:%d index %d", __func__, __LINE__, i);
    }

    free(expected_data);
}

static
void
test_cgb_acc_nop()
{
    int i;

    const int    m       = 4;
    const int    n       = m;
    const int    ku      = 2;
    const int    kl      = 1;
    const int    lda     = 5;
    const int    ldb     = 6;
    const complex_float alpha = 0.0; /* NOP zero */
    const complex_float beta  = 1.0; /* NOP one */

    /* Negative one represents values outside of the band */
    /* Decimal parts flag regions outside the matrix entirely */
    const float  a_data[][2]  = {
        /*lda buffer*/ /*ku2*/  /*ku1*/  /*diag*/ /*kl1*/
        {-1.1,-1.1},   {-1,-1}, {-1,-1}, { 1, 1}, { 2, 2},
        {-1.1,-1.1},   {-1,-1}, { 3, 3}, { 4, 4}, { 5, 5},
        {-1.1,-1.1},   { 6, 6}, { 7, 7}, { 8, 8}, { 9, 9},
        {-1.1,-1.1},   {10,10}, {11,11}, {12,12}, {-1,-1}
    };
    float  b_data[][2]  = {
        /*ldb buffer*/            /*ku2*/  /*ku1*/  /*diag*/ /*kl1*/
        {-2.2,-2.2}, {-1.1,-1.1}, {-1,-1}, {-1,-1}, { 1, 1}, { 2, 2},
        {-2.2,-2.2}, {-1.1,-1.1}, {-1,-1}, { 3, 3}, { 4, 4}, { 5, 5},
        {-2.2,-2.2}, {-1.1,-1.1}, { 6, 6}, { 7, 7}, { 8, 8}, { 9, 9},
        {-2.2,-2.2}, {-1.1,-1.1}, {10,10}, {11,11}, {12,12}, {-1,-1}
    };
    const float  (*a)[2] = a_data + lda-(ku+1+kl);
    float        (*b)[2] = b_data + ldb-(ku+1+kl);
    float        (*expected_data)[2];

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = nb;
    gsl_test_int(nb, nexpected, "Expected results' length");

    /* Make a clean copy of the original data */
    expected_data = malloc(nb*sizeof(b_data[0]));
    memcpy(expected_data, b_data, nb*sizeof(b_data[0]));

    suzerain_blas_cgb_acc( m, n, kl, ku, alpha, (complex_float *) a, lda,
                                         beta,  (complex_float *) b, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i][0], expected_data[i][0], GSL_FLT_EPSILON,
                "%s:%d real index %d", __func__, __LINE__, i);
        gsl_test_abs(b_data[i][1], expected_data[i][1], GSL_FLT_EPSILON,
                "%s:%d imag index %d", __func__, __LINE__, i);
    }

    free(expected_data);
}

static
void
test_zgb_acc_nop()
{
    int i;

    const int    m        = 4;
    const int    n        = m;
    const int    ku       = 2;
    const int    kl       = 1;
    const int    lda      = 5;
    const int    ldb      = 6;
    const complex_double alpha = 0.0; /* NOP zero */
    const complex_double beta  = 1.0; /* NOP one */

    /* Negative one represents values outside of the band */
    /* Decimal parts flag regions outside the matrix entirely */
    const double a_data[][2]  = {
        /*lda buffer*/ /*ku2*/  /*ku1*/  /*diag*/ /*kl1*/
        {-1.1,-1.1},   {-1,-1}, {-1,-1}, { 1, 1}, { 2, 2},
        {-1.1,-1.1},   {-1,-1}, { 3, 3}, { 4, 4}, { 5, 5},
        {-1.1,-1.1},   { 6, 6}, { 7, 7}, { 8, 8}, { 9, 9},
        {-1.1,-1.1},   {10,10}, {11,11}, {12,12}, {-1,-1}
    };
    double b_data[][2]  = {
        /*ldb buffer*/            /*ku2*/  /*ku1*/  /*diag*/ /*kl1*/
        {-2.2,-2.2}, {-1.1,-1.1}, {-1,-1}, {-1,-1}, { 1, 1}, { 2, 2},
        {-2.2,-2.2}, {-1.1,-1.1}, {-1,-1}, { 3, 3}, { 4, 4}, { 5, 5},
        {-2.2,-2.2}, {-1.1,-1.1}, { 6, 6}, { 7, 7}, { 8, 8}, { 9, 9},
        {-2.2,-2.2}, {-1.1,-1.1}, {10,10}, {11,11}, {12,12}, {-1,-1}
    };
    const double (*a)[2] = a_data + lda-(ku+1+kl);
    double       (*b)[2] = b_data + ldb-(ku+1+kl);
    double       (*expected_data)[2];

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = nb;
    gsl_test_int(nb, nexpected, "Expected results' length");

    /* Make a clean copy of the original data */
    expected_data = malloc(nb*sizeof(b_data[0]));
    memcpy(expected_data, b_data, nb*sizeof(b_data[0]));

    suzerain_blas_zgb_acc( m, n, kl, ku, alpha, (complex_double *) a, lda,
                                         beta,  (complex_double *) b, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i][0], expected_data[i][0], GSL_FLT_EPSILON,
                "%s:%d real index %d", __func__, __LINE__, i);
        gsl_test_abs(b_data[i][1], expected_data[i][1], GSL_FLT_EPSILON,
                "%s:%d imag index %d", __func__, __LINE__, i);
    }

    free(expected_data);
}

static
void
test_dgb_diag_scale_acc1(const double beta)
{
    int i;

    const int    m     = 10;
    const int    n     = 10;
    const int    ku    = 1;
    const int    kl    = 2;
    const int    lda   = 5;
    const int    ldb   = 6;
    const double alpha = 2.0;

    /* Negative one represents values outside of the band */
    /* Decimal parts flag regions outside the matrix entirely */
    const double a_data[] = {
        /*lda buffer*/ /*ku1*/ /*diag*/ /*kl1*/ /*kl2*/
                 -1.1,     -1,       0,    10,      20,
                 -1.1,      1,      11,    21,      31,
                 -1.1,     12,      22,    32,      42,
                 -1.1,     23,      33,    43,      53,
                 -1.1,     34,      44,    54,      64,
                 -1.1,     45,      55,    65,      75,
                 -1.1,     56,      66,    76,      86,
                 -1.1,     67,      77,    87,      97,
                 -1.1,     78,      88,    98,      -1,
                 -1.1,     89,      99,    -1,      -1
    };

    double b_data[] = {
        /*lda buffer*/ /*ku1*/ /*diag*/ /*kl1*/ /*kl2*/
          -2.2,  -1.1,     -1,    0.00,  0.10,    0.20,
          -2.2,  -1.1,   0.01,    0.11,  0.21,    0.31,
          -2.2,  -1.1,   0.12,    0.22,  0.32,    0.42,
          -2.2,  -1.1,   0.23,    0.33,  0.43,    0.53,
          -2.2,  -1.1,   0.34,    0.44,  0.54,    0.64,
          -2.2,  -1.1,   0.45,    0.55,  0.65,    0.75,
          -2.2,  -1.1,   0.56,    0.66,  0.76,    0.86,
          -2.2,  -1.1,   0.67,    0.77,  0.87,    0.97,
          -2.2,  -1.1,   0.78,    0.88,  0.98,      -1,
          -2.2,  -1.1,   0.89,    0.99,    -1,      -1
    };

    const double d_data[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 };

#define R(d,a,b) (alpha*d_data[d]*a + beta*b)
    const double expected_data[] = {
        /*lda buffer*/       /*ku1*/      /*diag*/       /*kl1*/       /*kl2*/
          -2.2,  -1.1,           -1, R(0, 0,0.00), R(0,10,0.10), R(0,20,0.20),
          -2.2,  -1.1, R(1, 1,0.01), R(1,11,0.11), R(1,21,0.21), R(1,31,0.31),
          -2.2,  -1.1, R(2,12,0.12), R(2,22,0.22), R(2,32,0.32), R(2,42,0.42),
          -2.2,  -1.1, R(3,23,0.23), R(3,33,0.33), R(3,43,0.43), R(3,53,0.53),
          -2.2,  -1.1, R(4,34,0.34), R(4,44,0.44), R(4,54,0.54), R(4,64,0.64),
          -2.2,  -1.1, R(5,45,0.45), R(5,55,0.55), R(5,65,0.65), R(5,75,0.75),
          -2.2,  -1.1, R(6,56,0.56), R(6,66,0.66), R(6,76,0.76), R(6,86,0.86),
          -2.2,  -1.1, R(7,67,0.67), R(7,77,0.77), R(7,87,0.87), R(7,97,0.97),
          -2.2,  -1.1, R(8,78,0.78), R(8,88,0.88), R(8,98,0.98),           -1,
          -2.2,  -1.1, R(9,89,0.89), R(9,99,0.99),         -1,             -1
    };
#undef R

    const double *a = a_data + lda-(ku+1+kl);
    double       *b = b_data + ldb-(ku+1+kl);

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = sizeof(expected_data)/sizeof(expected_data[0]);
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blasext_dgb_diag_scale_acc('R', m, n, kl, ku,
                                        alpha, d_data, 1, a, 1, lda,
                                        beta,             b, 1, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i], expected_data[i], GSL_DBL_EPSILON,
                "%s:%d index %d", __func__, __LINE__, i);
    }
}

static
void
test_dgb_diag_scale_acc2(const double beta)
{
    int i;

    const int    m     = 10;
    const int    n     = 10;
    const int    ku    = 1;
    const int    kl    = 2;
    const int    lda   = 5;
    const int    ldb   = 6;
    const double alpha = 2.0;

    /* Negative one represents values outside of the band */
    /* Decimal parts flag regions outside the matrix entirely */
    const double a_data[] = {
        /*lda buffer*/ /*ku1*/ /*diag*/ /*kl1*/ /*kl2*/
                 -1.1,     -1,       0,    10,      20,
                 -1.1,      1,      11,    21,      31,
                 -1.1,     12,      22,    32,      42,
                 -1.1,     23,      33,    43,      53,
                 -1.1,     34,      44,    54,      64,
                 -1.1,     45,      55,    65,      75,
                 -1.1,     56,      66,    76,      86,
                 -1.1,     67,      77,    87,      97,
                 -1.1,     78,      88,    98,      -1,
                 -1.1,     89,      99,    -1,      -1
    };

    double b_data[] = {
        /*lda buffer*/ /*ku1*/ /*diag*/ /*kl1*/ /*kl2*/
          -2.2,  -1.1,     -1,    0.00,  0.10,    0.20,
          -2.2,  -1.1,   0.01,    0.11,  0.21,    0.31,
          -2.2,  -1.1,   0.12,    0.22,  0.32,    0.42,
          -2.2,  -1.1,   0.23,    0.33,  0.43,    0.53,
          -2.2,  -1.1,   0.34,    0.44,  0.54,    0.64,
          -2.2,  -1.1,   0.45,    0.55,  0.65,    0.75,
          -2.2,  -1.1,   0.56,    0.66,  0.76,    0.86,
          -2.2,  -1.1,   0.67,    0.77,  0.87,    0.97,
          -2.2,  -1.1,   0.78,    0.88,  0.98,      -1,
          -2.2,  -1.1,   0.89,    0.99,    -1,      -1
    };

    const double d_data[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 };

#define R(d,a,b) (alpha*d_data[d]*a + beta*b)
    const double expected_data[] = {
        /*lda buffer*/       /*ku1*/      /*diag*/       /*kl1*/       /*kl2*/
          -2.2,  -1.1,           -1, R(0, 0,0.00), R(1,10,0.10), R(2,20,0.20),
          -2.2,  -1.1, R(0, 1,0.01), R(1,11,0.11), R(2,21,0.21), R(3,31,0.31),
          -2.2,  -1.1, R(1,12,0.12), R(2,22,0.22), R(3,32,0.32), R(4,42,0.42),
          -2.2,  -1.1, R(2,23,0.23), R(3,33,0.33), R(4,43,0.43), R(5,53,0.53),
          -2.2,  -1.1, R(3,34,0.34), R(4,44,0.44), R(5,54,0.54), R(6,64,0.64),
          -2.2,  -1.1, R(4,45,0.45), R(5,55,0.55), R(6,65,0.65), R(7,75,0.75),
          -2.2,  -1.1, R(5,56,0.56), R(6,66,0.66), R(7,76,0.76), R(8,86,0.86),
          -2.2,  -1.1, R(6,67,0.67), R(7,77,0.77), R(8,87,0.87), R(9,97,0.97),
          -2.2,  -1.1, R(7,78,0.78), R(8,88,0.88), R(9,98,0.98),           -1,
          -2.2,  -1.1, R(8,89,0.89), R(9,99,0.99),         -1,             -1
    };
#undef R

    const double *a = a_data + lda-(ku+1+kl);
    double       *b = b_data + ldb-(ku+1+kl);

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = sizeof(expected_data)/sizeof(expected_data[0]);
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blasext_dgb_diag_scale_acc('L', m, n, kl, ku,
                                        alpha, d_data, 1, a, 1, lda,
                                        beta,             b, 1, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i], expected_data[i], GSL_DBL_EPSILON,
                "%s:%d index %d", __func__, __LINE__, i);
    }
}

static
void
test_sgb_diag_scale_acc1(const float  beta)
{
    int i;

    const int    m     = 10;
    const int    n     = 10;
    const int    ku    = 1;
    const int    kl    = 2;
    const int    lda   = 5;
    const int    ldb   = 6;
    const float  alpha = 2.0;

    /* Negative one represents values outside of the band */
    /* Decimal parts flag regions outside the matrix entirely */
    const float  a_data[] = {
        /*lda buffer*/ /*ku1*/ /*diag*/ /*kl1*/ /*kl2*/
                 -1.1,     -1,       0,    10,      20,
                 -1.1,      1,      11,    21,      31,
                 -1.1,     12,      22,    32,      42,
                 -1.1,     23,      33,    43,      53,
                 -1.1,     34,      44,    54,      64,
                 -1.1,     45,      55,    65,      75,
                 -1.1,     56,      66,    76,      86,
                 -1.1,     67,      77,    87,      97,
                 -1.1,     78,      88,    98,      -1,
                 -1.1,     89,      99,    -1,      -1
    };

    float  b_data[] = {
        /*lda buffer*/ /*ku1*/ /*diag*/ /*kl1*/ /*kl2*/
          -2.2,  -1.1,     -1,    0.00,  0.10,    0.20,
          -2.2,  -1.1,   0.01,    0.11,  0.21,    0.31,
          -2.2,  -1.1,   0.12,    0.22,  0.32,    0.42,
          -2.2,  -1.1,   0.23,    0.33,  0.43,    0.53,
          -2.2,  -1.1,   0.34,    0.44,  0.54,    0.64,
          -2.2,  -1.1,   0.45,    0.55,  0.65,    0.75,
          -2.2,  -1.1,   0.56,    0.66,  0.76,    0.86,
          -2.2,  -1.1,   0.67,    0.77,  0.87,    0.97,
          -2.2,  -1.1,   0.78,    0.88,  0.98,      -1,
          -2.2,  -1.1,   0.89,    0.99,    -1,      -1
    };

    const float  d_data[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 };

#define R(d,a,b) (alpha*d_data[d]*a + beta*b)
    const float  expected_data[] = {
        /*lda buffer*/       /*ku1*/      /*diag*/       /*kl1*/       /*kl2*/
          -2.2,  -1.1,           -1, R(0, 0,0.00), R(0,10,0.10), R(0,20,0.20),
          -2.2,  -1.1, R(1, 1,0.01), R(1,11,0.11), R(1,21,0.21), R(1,31,0.31),
          -2.2,  -1.1, R(2,12,0.12), R(2,22,0.22), R(2,32,0.32), R(2,42,0.42),
          -2.2,  -1.1, R(3,23,0.23), R(3,33,0.33), R(3,43,0.43), R(3,53,0.53),
          -2.2,  -1.1, R(4,34,0.34), R(4,44,0.44), R(4,54,0.54), R(4,64,0.64),
          -2.2,  -1.1, R(5,45,0.45), R(5,55,0.55), R(5,65,0.65), R(5,75,0.75),
          -2.2,  -1.1, R(6,56,0.56), R(6,66,0.66), R(6,76,0.76), R(6,86,0.86),
          -2.2,  -1.1, R(7,67,0.67), R(7,77,0.77), R(7,87,0.87), R(7,97,0.97),
          -2.2,  -1.1, R(8,78,0.78), R(8,88,0.88), R(8,98,0.98),           -1,
          -2.2,  -1.1, R(9,89,0.89), R(9,99,0.99),         -1,             -1
    };
#undef R

    const float  *a = a_data + lda-(ku+1+kl);
    float        *b = b_data + ldb-(ku+1+kl);

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = sizeof(expected_data)/sizeof(expected_data[0]);
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blasext_sgb_diag_scale_acc('R', m, n, kl, ku,
                                        alpha, d_data, 1, a, 1, lda,
                                        beta,             b, 1, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i], expected_data[i], GSL_FLT_EPSILON,
                "%s:%d index %d", __func__, __LINE__, i);
    }
}

static
void
test_sgb_diag_scale_acc2(const float  beta)
{
    int i;

    const int    m     = 10;
    const int    n     = 10;
    const int    ku    = 1;
    const int    kl    = 2;
    const int    lda   = 5;
    const int    ldb   = 6;
    const float  alpha = 2.0;

    /* Negative one represents values outside of the band */
    /* Decimal parts flag regions outside the matrix entirely */
    const float  a_data[] = {
        /*lda buffer*/ /*ku1*/ /*diag*/ /*kl1*/ /*kl2*/
                 -1.1,     -1,       0,    10,      20,
                 -1.1,      1,      11,    21,      31,
                 -1.1,     12,      22,    32,      42,
                 -1.1,     23,      33,    43,      53,
                 -1.1,     34,      44,    54,      64,
                 -1.1,     45,      55,    65,      75,
                 -1.1,     56,      66,    76,      86,
                 -1.1,     67,      77,    87,      97,
                 -1.1,     78,      88,    98,      -1,
                 -1.1,     89,      99,    -1,      -1
    };

    float  b_data[] = {
        /*lda buffer*/ /*ku1*/ /*diag*/ /*kl1*/ /*kl2*/
          -2.2,  -1.1,     -1,    0.00,  0.10,    0.20,
          -2.2,  -1.1,   0.01,    0.11,  0.21,    0.31,
          -2.2,  -1.1,   0.12,    0.22,  0.32,    0.42,
          -2.2,  -1.1,   0.23,    0.33,  0.43,    0.53,
          -2.2,  -1.1,   0.34,    0.44,  0.54,    0.64,
          -2.2,  -1.1,   0.45,    0.55,  0.65,    0.75,
          -2.2,  -1.1,   0.56,    0.66,  0.76,    0.86,
          -2.2,  -1.1,   0.67,    0.77,  0.87,    0.97,
          -2.2,  -1.1,   0.78,    0.88,  0.98,      -1,
          -2.2,  -1.1,   0.89,    0.99,    -1,      -1
    };

    const float  d_data[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 };

#define R(d,a,b) (alpha*d_data[d]*a + beta*b)
    const float  expected_data[] = {
        /*lda buffer*/       /*ku1*/      /*diag*/       /*kl1*/       /*kl2*/
          -2.2,  -1.1,           -1, R(0, 0,0.00), R(1,10,0.10), R(2,20,0.20),
          -2.2,  -1.1, R(0, 1,0.01), R(1,11,0.11), R(2,21,0.21), R(3,31,0.31),
          -2.2,  -1.1, R(1,12,0.12), R(2,22,0.22), R(3,32,0.32), R(4,42,0.42),
          -2.2,  -1.1, R(2,23,0.23), R(3,33,0.33), R(4,43,0.43), R(5,53,0.53),
          -2.2,  -1.1, R(3,34,0.34), R(4,44,0.44), R(5,54,0.54), R(6,64,0.64),
          -2.2,  -1.1, R(4,45,0.45), R(5,55,0.55), R(6,65,0.65), R(7,75,0.75),
          -2.2,  -1.1, R(5,56,0.56), R(6,66,0.66), R(7,76,0.76), R(8,86,0.86),
          -2.2,  -1.1, R(6,67,0.67), R(7,77,0.77), R(8,87,0.87), R(9,97,0.97),
          -2.2,  -1.1, R(7,78,0.78), R(8,88,0.88), R(9,98,0.98),           -1,
          -2.2,  -1.1, R(8,89,0.89), R(9,99,0.99),         -1,             -1
    };
#undef R

    const float  *a = a_data + lda-(ku+1+kl);
    float        *b = b_data + ldb-(ku+1+kl);

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = sizeof(expected_data)/sizeof(expected_data[0]);
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blasext_sgb_diag_scale_acc('L', m, n, kl, ku,
                                        alpha, d_data, 1, a, 1, lda,
                                        beta,             b, 1, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i], expected_data[i], GSL_FLT_EPSILON,
                "%s:%d index %d", __func__, __LINE__, i);
    }
}

static
void
test_zaxpy_d()
{
    int i;
    const complex_double alpha = complex_double(2.0, 3.0);

    {
        /* Both x and y stride 1 */
        const double x[]      = {1.0, 2.0, 3.0};
        const int    incx     = 1;
        double       y[][2]   = {{4,-4}, {5,-5}, {6,-6}};
        const int    incy     = 1;
        const int    nx       = sizeof(x)/sizeof(x[0]);
        const int    ny       = sizeof(y)/sizeof(y[0]);
        const double expected[][2] = {
            {2*1+4,3*1-4}, {2*2+5,3*2-5}, {2*3+6,3*3-6}
        };
        const int nexpected = sizeof(expected)/sizeof(expected[0]);

        gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
        gsl_test_int(ny, nexpected, "Expected results' length");

        suzerain_blas_zaxpy_d(nx/incx, alpha, x, incx,
                                (complex_double *) y, incy);
        for (i = 0; i < nexpected; ++i) {
            gsl_test_abs(y[i][0], expected[i][0], GSL_DBL_EPSILON,
                    "%s:%d real index %d", __func__, __LINE__, i);
            gsl_test_abs(y[i][1], expected[i][1], GSL_DBL_EPSILON,
                    "%s:%d imag index %d", __func__, __LINE__, i);
        }
    }

    {
        /* x stride 2 */
        const double x[]      = {1.0, -1, 2.0, -2, 3.0, -3};
        const int    incx     = 2;
        double       y[][2]   = {{4,-4}, {5,-5}, {6,-6}};
        const int    incy     = 1;
        const int    nx       = sizeof(x)/sizeof(x[0]);
        const int    ny       = sizeof(y)/sizeof(y[0]);
        const double expected[][2] = {
            {2*1+4,3*1-4}, {2*2+5,3*2-5}, {2*3+6,3*3-6}
        };
        const int nexpected = sizeof(expected)/sizeof(expected[0]);

        gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
        gsl_test_int(ny, nexpected, "Expected results' length");

        suzerain_blas_zaxpy_d(nx/incx, alpha, x, incx,
                                (complex_double *) y, incy);
        for (i = 0; i < nexpected; ++i) {
            gsl_test_abs(y[i][0], expected[i][0], GSL_DBL_EPSILON,
                    "%s:%d real index %d", __func__, __LINE__, i);
            gsl_test_abs(y[i][1], expected[i][1], GSL_DBL_EPSILON,
                    "%s:%d imag index %d", __func__, __LINE__, i);
        }
    }

    {
        /* y stride 2 */
        const double x[]      = {1.0, 2.0, 3.0};
        const int    incx     = 1;
        double       y[][2]   = {{4,-4},{1,-1},{5,-5},{2,-2},{6,-6},{3,-3}};
        const int    incy     = 2;
        const int    nx       = sizeof(x)/sizeof(x[0]);
        const int    ny       = sizeof(y)/sizeof(y[0]);
        const double expected[][2] = {
            {2*1+4,3*1-4}, {1,-1}, {2*2+5,3*2-5}, {2,-2}, {2*3+6,3*3-6}, {3,-3}
        };
        const int nexpected = sizeof(expected)/sizeof(expected[0]);

        gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
        gsl_test_int(ny, nexpected, "Expected results' length");

        suzerain_blas_zaxpy_d(nx/incx, alpha, x, incx,
                                (complex_double *) y, incy);
        for (i = 0; i < nexpected; ++i) {
            gsl_test_abs(y[i][0], expected[i][0], GSL_DBL_EPSILON,
                    "%s:%d real index %d", __func__, __LINE__, i);
            gsl_test_abs(y[i][1], expected[i][1], GSL_DBL_EPSILON,
                    "%s:%d imag index %d", __func__, __LINE__, i);
        }
    }

    {
        /* Both x and y stride 2 */
        const double x[]      = {1.0, -555, 2.0, -666, 3.0, -777};
        const int    incx     = 2;
        double       y[][2]   = {{4,-4},{1,-1},{5,-5},{2,-2},{6,-6},{3,-3}};
        const int    incy     = 2;
        const int    nx       = sizeof(x)/sizeof(x[0]);
        const int    ny       = sizeof(y)/sizeof(y[0]);
        const double expected[][2] = {
            {2*1+4,3*1-4}, {1,-1}, {2*2+5,3*2-5}, {2,-2}, {2*3+6,3*3-6}, {3,-3}
        };
        const int nexpected = sizeof(expected)/sizeof(expected[0]);

        gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
        gsl_test_int(ny, nexpected, "Expected results' length");

        suzerain_blas_zaxpy_d(nx/incx, alpha, x, incx,
                                (complex_double *) y, incy);
        for (i = 0; i < nexpected; ++i) {
            gsl_test_abs(y[i][0], expected[i][0], GSL_DBL_EPSILON,
                    "%s:%d real index %d", __func__, __LINE__, i);
            gsl_test_abs(y[i][1], expected[i][1], GSL_DBL_EPSILON,
                    "%s:%d imag index %d", __func__, __LINE__, i);
        }
    }
}

static
void
test_zaxpby_d()
{
    int i;
    const complex_double alpha = complex_double(2.0, 3.0);
    const complex_double beta  = complex_double(5.0, 7.0);

    {
        /* Both x and y stride 1 */
        const double x[]      = {1.0, 2.0, 3.0};
        const int    incx     = 1;
        double       y[][2]   = {{4,-4}, {5,-5}, {6,-6}};
        const int    incy     = 1;
        const int    nx       = sizeof(x)/sizeof(x[0]);
        const int    ny       = sizeof(y)/sizeof(y[0]);
        const double expected[][2] = {
            {50, 11}, {64, 16}, {78, 21}
        };
        const int nexpected = sizeof(expected)/sizeof(expected[0]);

        gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
        gsl_test_int(ny, nexpected, "Expected results' length");

        suzerain_blas_zaxpby_d(nx/incx, alpha, x, incx, beta,
                                  (complex_double *) y, incy);
        for (i = 0; i < nexpected; ++i) {
            gsl_test_abs(y[i][0], expected[i][0], GSL_DBL_EPSILON,
                    "%s:%d real index %d", __func__, __LINE__, i);
            gsl_test_abs(y[i][1], expected[i][1], GSL_DBL_EPSILON,
                    "%s:%d imag index %d", __func__, __LINE__, i);
        }
    }

    {
        /* x stride 2 */
        const double x[]      = {1.0, -1, 2.0, -2, 3.0, -3};
        const int    incx     = 2;
        double       y[][2]   = {{4,-4}, {5,-5}, {6,-6}};
        const int    incy     = 1;
        const int    nx       = sizeof(x)/sizeof(x[0]);
        const int    ny       = sizeof(y)/sizeof(y[0]);
        const double expected[][2] = {
            {50, 11}, {64, 16}, {78, 21}
        };
        const int nexpected = sizeof(expected)/sizeof(expected[0]);

        gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
        gsl_test_int(ny, nexpected, "Expected results' length");

        suzerain_blas_zaxpby_d(nx/incx, alpha, x, incx, beta,
                                  (complex_double *) y, incy);
        for (i = 0; i < nexpected; ++i) {
            gsl_test_abs(y[i][0], expected[i][0], GSL_DBL_EPSILON,
                    "%s:%d real index %d", __func__, __LINE__, i);
            gsl_test_abs(y[i][1], expected[i][1], GSL_DBL_EPSILON,
                    "%s:%d imag index %d", __func__, __LINE__, i);
        }
    }

    {
        /* y stride 2 */
        const double x[]      = {1.0, 2.0, 3.0};
        const int    incx     = 1;
        double       y[][2]   = {{4,-4},{1,-1},{5,-5},{2,-2},{6,-6},{3,-3}};
        const int    incy     = 2;
        const int    nx       = sizeof(x)/sizeof(x[0]);
        const int    ny       = sizeof(y)/sizeof(y[0]);
        const double expected[][2] = {
            {50, 11}, {1,-1}, {64, 16}, {2,-2}, {78, 21}, {3,-3}
        };
        const int nexpected = sizeof(expected)/sizeof(expected[0]);

        gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
        gsl_test_int(ny, nexpected, "Expected results' length");

        suzerain_blas_zaxpby_d(nx/incx, alpha, x, incx, beta,
                                  (complex_double *) y, incy);
        for (i = 0; i < nexpected; ++i) {
            gsl_test_abs(y[i][0], expected[i][0], GSL_DBL_EPSILON,
                    "%s:%d real index %d", __func__, __LINE__, i);
            gsl_test_abs(y[i][1], expected[i][1], GSL_DBL_EPSILON,
                    "%s:%d imag index %d", __func__, __LINE__, i);
        }
    }

    {
        /* Both x and y stride 2 */
        const double x[]      = {1.0, -555, 2.0, -666, 3.0, -777};
        const int    incx     = 2;
        double       y[][2]   = {{4,-4},{1,-1},{5,-5},{2,-2},{6,-6},{3,-3}};
        const int    incy     = 2;
        const int    nx       = sizeof(x)/sizeof(x[0]);
        const int    ny       = sizeof(y)/sizeof(y[0]);
        const double expected[][2] = {
            {50, 11}, {1,-1}, {64, 16}, {2,-2}, {78, 21}, {3,-3}
        };
        const int nexpected = sizeof(expected)/sizeof(expected[0]);

        gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
        gsl_test_int(ny, nexpected, "Expected results' length");

        suzerain_blas_zaxpby_d(nx/incx, alpha, x, incx, beta,
                                  (complex_double *) y, incy);
        for (i = 0; i < nexpected; ++i) {
            gsl_test_abs(y[i][0], expected[i][0], GSL_DBL_EPSILON,
                    "%s:%d real index %d", __func__, __LINE__, i);
            gsl_test_abs(y[i][1], expected[i][1], GSL_DBL_EPSILON,
                    "%s:%d imag index %d", __func__, __LINE__, i);
        }
    }
}

static
void
test_blasext_i2s_zaxpby2()
{
    { /* Test pure copy for single, stride one vectors */

        const int N = 4;
        double z[N*2];
        double z_re[N];
        double z_im[N];

        for (int i = 0; i < N; ++i) {
            z[2*i]   = i+1;   /* Real part */
            z[2*i+1] = i+1.5; /* Imag part */
        }

        suzerain_blasext_i2s_zaxpby2(
                1, N, NULL, z, 1, N,
                NULL, z_re, 1, N, z_im, 1, N);

        for (int i = 0; i < N; ++i) {
            gsl_test_abs(z_re[i], i+1, GSL_DBL_EPSILON,
                        "%s:%d real part of index %d",
                        __func__, __LINE__, i);
            gsl_test_abs(z_im[i], i+1.5, GSL_DBL_EPSILON,
                        "%s:%d imag part of index %d",
                        __func__, __LINE__, i);
        }
    }

    { /* Test pure copy for multiple, strided vectors */
        const int M       = 3;
        const int N       = 4;
        const int incz    = 2;
        const int incz_re = 3;
        const int incz_im = 4;
        const int ldz     = N*incz + 2;
        const int ldz_re  = N*incz_re + 3;
        const int ldz_im  = N*incz_im + 5;
        double z[2*(N*incz)*(M*ldz)];
        double z_re[(N*incz_re)*(M*ldz_re)];
        double z_im[(N*incz_im)*(M*ldz_im)];

        memset(z_re, 'r', sizeof(z_re)/sizeof(z_re[0])); // Avoid NaNs
        memset(z_im, 'r', sizeof(z_im)/sizeof(z_im[0])); // Avoid NaNs

        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < N; ++i) {
                z[2*(i*incz + j*ldz)]   = i+1;   /* Real part */
                z[2*(i*incz + j*ldz)+1] = i+1.5; /* Imag part */
            }
        }

        suzerain_blasext_i2s_zaxpby2(
                M, N, NULL,
                z, incz, ldz,
                NULL,
                z_re, incz_re, ldz_re,
                z_im, incz_im, ldz_im);

        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < N; ++i) {
                gsl_test_abs(z_re[i*incz_re+j*ldz_re], i+1, GSL_DBL_EPSILON,
                        "%s:%d real part of index %d",
                        __func__, __LINE__, i);
                gsl_test_abs(z_im[i*incz_im+j*ldz_im], i+1.5, GSL_DBL_EPSILON,
                        "%s:%d imag part of index %d",
                        __func__, __LINE__, i);
            }
        }
    }

    { /* Test real scaling for single, stride one vectors */
        const int N = 4;
        const double alpha[2] = {2.0, 0.0};
        double z[2*N];
        double z_re[N];
        double z_im[N];

        for (int i = 0; i < N; ++i) {
            z[2*i]   = i+1;   /* Real part */
            z[2*i+1] = i+1.5; /* Imag part */
        }

        suzerain_blasext_i2s_zaxpby2(
                1, N, alpha, z, 1, N,
                NULL,
                z_re, 1, N, z_im, 1, N);

        for (int i = 0; i < N; ++i) {
            gsl_test_abs(z_re[i], 2*(i+1), GSL_DBL_EPSILON,
                        "%s:%d real part of index",
                        __func__, __LINE__, i);
            gsl_test_abs(z_im[i], 2*(i+1.5), GSL_DBL_EPSILON,
                        "%s:%d imag part of index %d",
                        __func__, __LINE__, i);
        }
    }

    { /* Test real scaling for multiple, strided vectors */
        const int M       = 3;
        const int N       = 4;
        const double alpha[2] = {2.0, 0.0};
        const int incz    = 2;
        const int incz_re = 3;
        const int incz_im = 4;
        const int ldz     = N*incz + 2;
        const int ldz_re  = N*incz_re + 3;
        const int ldz_im  = N*incz_im + 5;
        double z[2*(N*incz)*(M*ldz)];
        double z_re[(N*incz_re)*(M*ldz_re)];
        double z_im[(N*incz_im)*(M*ldz_im)];

        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < N; ++i) {
                z[2*(i*incz + j*ldz)]   = i+1;   /* Real part */
                z[2*(i*incz + j*ldz)+1] = i+1.5; /* Imag part */
            }
        }

        suzerain_blasext_i2s_zaxpby2(
                M, N, alpha,
                z, incz, ldz,
                NULL,
                z_re, incz_re, ldz_re,
                z_im, incz_im, ldz_im);

        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < N; ++i) {
                gsl_test_abs(z_re[i*incz_re+j*ldz_re],
                        2*(i+1), GSL_DBL_EPSILON,
                        "%s:%d real part of index %d",
                        __func__, __LINE__, i);
                gsl_test_abs(z_im[i*incz_im+j*ldz_im],
                        2*(i+1.5), GSL_DBL_EPSILON,
                        "%s:%d imag part of index %d",
                        __func__, __LINE__, i);
            }
        }
    }

    { /* Test complex scaling for contiguous vectors */
        const int N = 4;
        const double alpha[2] = {3.0, 2.0};
        const double beta[2]  = {4.0, 5.0};
        double z[2*N];
        double z_re[N];
        double z_im[N];

        for (int i = 0; i < N; ++i) {
            z[2*(i)] = (i+1)*(i+1);
            z[2*(i)+1] = i+1.5;
            z_re[i] = i+7;
            z_im[i] = i+11;
        }

        suzerain_blasext_i2s_zaxpby2(
                1, N, alpha, z, 1, N,
                beta, z_re, 1, N, z_im, 1, N);

        for (int i = 0; i < N; ++i) {
            gsl_test_abs(z_re[i], 3*i*i + 3*i - 27,
                        GSL_DBL_EPSILON,
                        "%s:%d real part of index %d",
                        __func__, __LINE__, i);
            gsl_test_abs(z_im[i], 2*i*i + 16*i + 85.5,
                        GSL_DBL_EPSILON,
                        "%s:%d imag part of index %d",
                        __func__, __LINE__, i);
        }
    }

    { /* Test complex scaling for multiple, strided vectors */
        const int M       = 3;
        const int N       = 4;
        const double alpha[2] = {3.0, 2.0};
        const double beta[2]  = {4.0, 5.0};
        const int incz    = 2;
        const int incz_re = 3;
        const int incz_im = 4;
        const int ldz     = N*incz + 2;
        const int ldz_re  = N*incz_re + 3;
        const int ldz_im  = N*incz_im + 5;
        double z[2*(N*incz)*(M*ldz)];
        double z_re[(N*incz_re)*(M*ldz_re)];
        double z_im[(N*incz_im)*(M*ldz_im)];

        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < N; ++i) {
                z[2*(i*incz + j*ldz)] = (i+1)*(i+1);
                z[2*(i*incz + j*ldz)+1] = i+1.5;
                z_re[i*incz_re + j*ldz_re] = i+7;
                z_im[i*incz_im + j*ldz_im] = i+11;
            }
        }

        suzerain_blasext_i2s_zaxpby2(
                M, N, alpha,
                z, incz, ldz,
                beta,
                z_re, incz_re, ldz_re,
                z_im, incz_im, ldz_im);

        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < N; ++i) {
                gsl_test_abs(z_re[i*incz_re+j*ldz_re], 3*i*i + 3*i - 27,
                        GSL_DBL_EPSILON,
                        "%s:%d real part of index %d",
                        __func__, __LINE__, i);
                gsl_test_abs(z_im[i*incz_im+j*ldz_im], 2*i*i + 16*i + 85.5,
                        GSL_DBL_EPSILON,
                        "%s:%d imag part of index %d",
                        __func__, __LINE__, i);
            }
        }
    }
}

static
void
test_zgb_acc_d1()
{
    int i;

    const int    m        = 5;
    const int    n        = 5;
    const int    ku       = 2;
    const int    kl       = 1;
    const int    lda      = 5;
    const int    ldb      = 6;
    const complex_double alpha = complex_double(2.0, -3.0);
    const complex_double beta  = complex_double(7.0, -5.0);

    /* Negative one represents values outside of the band */
    /* Decimal parts flag regions outside the matrix entirely */
    const double a_data[] = {
        /*lda buffer*/ /*ku2*/ /*ku1*/ /*diag*/ /*kl1*/
        -1.1,          -1,     -1,      1,       2,
        -1.1,          -1,      3,      4,       5,
        -1.1,           6,      7,      8,       9,
        -1.1,           6,      7,      8,       9,
        -1.1,          10,     11,     12,      -1
    };
    double b_data[][2]  = {
        /*ldb buffer*/            /*ku2*/  /*ku1*/  /*diag*/ /*kl1*/
        {-2.2,-2.2}, {-1.1,-1.1}, {-1,-1}, {-1,-1}, { 1, 1}, { 2, 2},
        {-2.2,-2.2}, {-1.1,-1.1}, {-1,-1}, { 3, 3}, { 4, 4}, { 5, 5},
        {-2.2,-2.2}, {-1.1,-1.1}, { 6, 6}, { 7, 7}, { 8, 8}, { 9, 9},
        {-2.2,-2.2}, {-1.1,-1.1}, { 6, 6}, { 7, 7}, { 8, 8}, { 9, 9},
        {-2.2,-2.2}, {-1.1,-1.1}, {10,10}, {11,11}, {12,12}, {-1,-1}
    };
    const double expected_data[][2]  = {
        /*ldb buffer*/            /*ku2*/  /*ku1*/  /*diag*/ /*kl1*/
        {-2.2,-2.2}, {-1.1,-1.1}, { -1, -1}, { -1,- 1}, { 14,- 1}, { 28,-2},
        {-2.2,-2.2}, {-1.1,-1.1}, { -1, -1}, { 42,- 3}, { 56,- 4}, { 70,-5},
        {-2.2,-2.2}, {-1.1,-1.1}, { 84, -6}, { 98,- 7}, {112,- 8}, {126,-9},
        {-2.2,-2.2}, {-1.1,-1.1}, { 84, -6}, { 98,- 7}, {112,- 8}, {126,-9},
        {-2.2,-2.2}, {-1.1,-1.1}, {140,-10}, {154,-11}, {168,-12}, { -1,-1}
    };
    const double  *a     = a_data + lda-(ku+1+kl);
    double       (*b)[2] = b_data + ldb-(ku+1+kl);

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = sizeof(expected_data)/sizeof(expected_data[0]);
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blas_zgb_acc_d( m, n, kl, ku, alpha, a, lda,
                              beta, (complex_double *) b, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i][0], expected_data[i][0], GSL_FLT_EPSILON,
                "%s:%d real index %d", __func__, __LINE__, i);
        gsl_test_abs(b_data[i][1], expected_data[i][1], GSL_FLT_EPSILON,
                "%s:%d imag index %d", __func__, __LINE__, i);
    }
}

static
void
test_zgb_acc_d2()
{
    int i;

    const int    m        = 5;
    const int    n        = 5;
    const int    ku       = 2;
    const int    kl       = 1;
    const int    lda      = ku + 1 + kl;
    const int    ldb      = ku + 1 + kl;
    const complex_double alpha = complex_double(2.0, -3.0);
    const complex_double beta  = complex_double(7.0, -5.0);

    /* Negative one represents values outside of the band */
    /* Decimal parts flag regions outside the matrix entirely */
    const double a_data[] = {
        /*ku2*/ /*ku1*/ /*diag*/ /*kl1*/
        -1,     -1,      1,       2,
        -1,      3,      4,       5,
         6,      7,      8,       9,
         6,      7,      8,       9,
        10,     11,     12,      -1
    };
    double b_data[][2]  = {
        /*ku2*/  /*ku1*/  /*diag*/ /*kl1*/
        {-1,-1}, {-1,-1}, { 1, 1}, { 2, 2},
        {-1,-1}, { 3, 3}, { 4, 4}, { 5, 5},
        { 6, 6}, { 7, 7}, { 8, 8}, { 9, 9},
        { 6, 6}, { 7, 7}, { 8, 8}, { 9, 9},
        {10,10}, {11,11}, {12,12}, {-1,-1}
    };
    const double expected_data[][2]  = {
        /*ku2*/  /*ku1*/  /*diag*/ /*kl1*/
        { -1,- 1}, { -1,- 1}, { 14,- 1}, { 28,-2},
        { -1,- 1}, { 42,- 3}, { 56,- 4}, { 70,-5},
        { 84,- 6}, { 98,- 7}, {112,- 8}, {126,-9},
        { 84,- 6}, { 98,- 7}, {112,- 8}, {126,-9},
        {140,-10}, {154,-11}, {168,-12}, {- 1,-1}
    };
    const double  *a     = a_data + lda-(ku+1+kl);
    double       (*b)[2] = b_data + ldb-(ku+1+kl);

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = sizeof(expected_data)/sizeof(expected_data[0]);
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blas_zgb_acc_d(m, n, kl, ku, alpha, a, lda,
                              beta, (complex_double *) b, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i][0], expected_data[i][0], GSL_FLT_EPSILON,
                "%s:%d real index %d", __func__, __LINE__, i);
        gsl_test_abs(b_data[i][1], expected_data[i][1], GSL_FLT_EPSILON,
                "%s:%d imag index %d", __func__, __LINE__, i);
    }
}

static
void
test_zgb_acc_d_nop()
{
    int i;

    const int    m        = 25;
    const int    n        = 5;
    const int    ku       = 2;
    const int    kl       = 1;
    const int    lda      = ku + 1 + kl;
    const int    ldb      = ku + 1 + kl;
    const complex_double alpha = 0; /* NOP zero */
    const complex_double beta  = 1; /* NOP one */

    /* Negative one represents values outside of the band */
    /* Decimal parts flag regions outside the matrix entirely */
    const double a_data[] = {
        /*ku2*/ /*ku1*/ /*diag*/ /*kl1*/
        -1,     -1,      1,       2,
        -1,      3,      4,       5,
         6,      7,      8,       9,
         6,      7,      8,       9,
        10,     11,     12,      -1
    };
    double b_data[][2]  = {
        /*ku2*/  /*ku1*/  /*diag*/ /*kl1*/
        {-1,-1}, {-1,-1}, { 1, 1}, { 2, 2},
        {-1,-1}, { 3, 3}, { 4, 4}, { 5, 5},
        { 6, 6}, { 7, 7}, { 8, 8}, { 9, 9},
        { 6, 6}, { 7, 7}, { 8, 8}, { 9, 9},
        {10,10}, {11,11}, {12,12}, {-1,-1}
    };
    const double  *a     = a_data + lda-(ku+1+kl);
    double       (*b)[2] = b_data + ldb-(ku+1+kl);
    double       (*expected_data)[2];

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = nb;
    gsl_test_int(nb, nexpected, "Expected results' length");

    /* Make a clean copy of the original data */
    expected_data = malloc(nb*sizeof(b_data[0]));
    memcpy(expected_data, b_data, nb*sizeof(b_data[0]));

    suzerain_blas_zgb_acc_d(m, n, kl, ku, alpha, a, lda,
                              beta, (complex_double *) b, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i][0], expected_data[i][0], GSL_FLT_EPSILON,
                "%s:%d real index %d", __func__, __LINE__, i);
        gsl_test_abs(b_data[i][1], expected_data[i][1], GSL_FLT_EPSILON,
                "%s:%d imag index %d", __func__, __LINE__, i);
    }

    free(expected_data);
}

static
void
test_blasext_sgbnorm1()
{
    const int m = 4, n = 5, kl = 3, ku = 4, lda = ku + 1 + kl;
    const float a[] = { 555, 555, 555, 555, 6,   -1,  -2,  -4,
                        555, 555, 555, -4,  10,  1,   -5,  555,
                        555, 555, 10,  10,  -7,  5,   555, 555,
                        555, -10, -1,  -9,  -2,  555, 555, 555,
                        -7,  3,   -5,  3,   555, 555, 555, 555 };

    float norm1;
    const int status = suzerain_blasext_sgbnorm1(m, n, kl, ku, a, lda, &norm1);
    gsl_test_int(status, 0, "%s:%d call success", __func__, __LINE__);
    gsl_test_abs(norm1, 32.0, GSL_FLT_EPSILON, "%s:%d norm result %d",
                 __func__, __LINE__, norm1);
}

static
void
test_blasext_dgbnorm1()
{
    const int m = 4, n = 5, kl = 3, ku = 4, lda = ku + 1 + kl;
    const double a[] = { 555, 555, 555, 555, 6,   -1,  -2,  -4,
                         555, 555, 555, -4,  10,  1,   -5,  555,
                         555, 555, 10,  10,  -7,  5,   555, 555,
                         555, -10, -1,  -9,  -2,  555, 555, 555,
                         -7,  3,   -5,  3,   555, 555, 555, 555 };

    double norm1;
    const int status = suzerain_blasext_dgbnorm1(m, n, kl, ku, a, lda, &norm1);
    gsl_test_int(status, 0, "%s:%d call success", __func__, __LINE__);
    gsl_test_abs(norm1, 32.0, GSL_DBL_EPSILON, "%s:%d norm result %d",
                 __func__, __LINE__, norm1);
}

static
void
test_blasext_cgbnorm1()
{
    const int m = 4, n = 5, kl = 3, ku = 4, lda = ku + 1 + kl;
    const float a[][2] = {
        {555,0},{555,0}, {555,0},{555,0}, {6,3}, {-1,5},  {-2,5},  {-4,-5},
        {555,0},{555,0}, {555,0},{-4,-4},{10,-7},{1,-2},  {-5,-7}, {555,0},
        {555,0},{555,0}, {10,-2},{10,-4},{-7,7}, {5,6},   {555,0}, {555,0},
        {555,0},{-10,-4},{-1,2}, {-9,0}, {-2,-7},{555,0}, {555,0}, {555,0},
        {-7,7}, {3,6},   {-5,-7},{3,-3}, {555,0},{555,0}, {555,0}, {555,0}
    };
    const float expected = 7.*sqrt(2.)+2.*sqrt(26.)+2.*sqrt(29.)+sqrt(61.);

    float norm1;
    const int status = suzerain_blasext_cgbnorm1(
            m, n, kl, ku, (complex_float *) a, lda, &norm1);
    gsl_test_int(status, 0, "%s:%d call success", __func__, __LINE__);
    gsl_test_abs(norm1, expected, GSL_FLT_EPSILON, "%s:%d norm result %d",
                 __func__, __LINE__, norm1);
}

static
void
test_blasext_zgbnorm1()
{
    const int m = 4, n = 5, kl = 3, ku = 4, lda = ku + 1 + kl;
    const double a[][2] = {
        {555,0},{555,0}, {555,0},{555,0}, {6,3}, {-1,5},  {-2,5},  {-4,-5},
        {555,0},{555,0}, {555,0},{-4,-4},{10,-7},{1,-2},  {-5,-7}, {555,0},
        {555,0},{555,0}, {10,-2},{10,-4},{-7,7}, {5,6},   {555,0}, {555,0},
        {555,0},{-10,-4},{-1,2}, {-9,0}, {-2,-7},{555,0}, {555,0}, {555,0},
        {-7,7}, {3,6},   {-5,-7},{3,-3}, {555,0},{555,0}, {555,0}, {555,0}
    };
    const double expected = 7.*sqrt(2.)+2.*sqrt(26.)+2.*sqrt(29.)+sqrt(61.);

    double norm1;
    const int status = suzerain_blasext_zgbnorm1(
            m, n, kl, ku, (complex_double *) a, lda, &norm1);
    gsl_test_int(status, 0, "%s:%d call success", __func__, __LINE__);
    gsl_test_abs(norm1, expected, GSL_DBL_EPSILON*100, "%s:%d norm result %d",
                 __func__, __LINE__, norm1);
}

static
void
test_lapack_dgbcon()
{
    /*
     * Test matrix is
     *      5  -3   0  0
     *      0   5  -3  0
     *      0   0   5 -3
     *      0   0   3 -1
     * which has a 1-norm of 11.  It has one subdiagonal but we'll use kl = 2
     * to test storage-padding details.  Extra kl = 2 rows of padding present
     * to store factorization multipliers.
     */

    const int kl = 2;
    const int ku = 1;
    const int n  = 4;
    double ab[] = { 0, 0,  0,  5, 0, 0,
                    0, 0, -3,  5, 0, 0,
                    0, 0, -3,  5, 3, 0,
                    0, 0, -3, -1, 0, 0 };
    int ipiv[n];
    const double norm1 = 11;

    // Prepare factorization
    const int f = suzerain_lapack_dgbtrf(n, n, kl, ku, ab, 2*kl + ku + 1, ipiv);
    gsl_test_int(f, 0, "%s:%d factorization success", __func__, __LINE__);

    // Compute reciprocal of condition number
    double rcond = -555;
    double work[3*n];
    int    iwork[n];
    const int g = suzerain_lapack_dgbcon('1', n, kl, ku, ab, 2*kl + ku + 1,
                                         ipiv, norm1, &rcond, work, iwork);
    gsl_test_int(g, 0, "%s:%d condition number estimation success",
                 __func__, __LINE__);

    // Check result against expected
    gsl_test_abs(rcond, 25.0/748.0, GSL_DBL_EPSILON,
            "%s:%d condition number estimation result %d",
            __func__, __LINE__, rcond);
}

static
SUZERAIN_DONTINLINE
void
test_blasext_ddemote()
{
    int info;
    double x[]   = {1, 2, 3, 4, 5, 6, 7};
    float  e[]   = {1, 2, 3, 4, 5, 6, 7};
    const int n = sizeof(x)/sizeof(x[0]);

    info = suzerain_blasext_ddemote(n, x);
    gsl_test_int(info, 0, "%s.%d report success", __func__, __LINE__);

    info = memcmp(x, e, sizeof(e));
    gsl_test_int(info, 0, "%s.%d matches expected", __func__, __LINE__);
}

static
SUZERAIN_DONTINLINE
void
test_blasext_dpromote()
{
    int info;
    double  e[]   = {1, 2, 3, 4, 5, 6, 7};
    float   x[]   = {1, 2, 3, 4, 5, 6, 7, 0, 0, 0, 0, 0, 0, 0};
    const int n = sizeof(e)/sizeof(e[0]);
    gsl_test_int(2*n, sizeof(x)/sizeof(x[0]),
            "%s.%d consistency", __func__, __LINE__);

    info = suzerain_blasext_dpromote(n, (void *) x);
    gsl_test_int(info, 0, "%s.%d report success", __func__, __LINE__);

    info = memcmp(x, e, sizeof(e));
    gsl_test_int(info, 0, "%s.%d matches expected", __func__, __LINE__);
}

static
SUZERAIN_DONTINLINE
void
test_blasext_zdemote()
{
    int info;
    double x[][2]   = {{1,-1}, {2,-2}, {3,-3}, {4,-4}, {5,-5}, {6,-6}, {7,-7}};
    float  e[][2]   = {{1,-1}, {2,-2}, {3,-3}, {4,-4}, {5,-5}, {6,-6}, {7,-7}};
    const int n = sizeof(x)/sizeof(x[0]);

    info = suzerain_blasext_zdemote(n, (void *) x);
    gsl_test_int(info, 0, "%s.%d report success", __func__, __LINE__);

    info = memcmp(x, e, sizeof(e));
    gsl_test_int(info, 0, "%s.%d matches expected", __func__, __LINE__);
}

static
SUZERAIN_DONTINLINE
void
test_blasext_zpromote()
{
    int info;
    float  x[][2]   = {{1,-1}, {2,-2}, {3,-3}, {4,-4}, {5,-5}, {6,-6}, {7,-7},
                       {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}};
    double e[][2]   = {{1,-1}, {2,-2}, {3,-3}, {4,-4}, {5,-5}, {6,-6}, {7,-7}};
    const int n = sizeof(e)/sizeof(e[0]);
    gsl_test_int(2*n, sizeof(x)/sizeof(x[0]),
            "%s.%d consistency", __func__, __LINE__);

    info = suzerain_blasext_zpromote(n, (void *) x);
    gsl_test_int(info, 0, "%s.%d report success", __func__, __LINE__);

    info = memcmp(x, e, sizeof(e));
    gsl_test_int(info, 0, "%s.%d matches expected", __func__, __LINE__);
}

// Form nastily conditioned test problems per Mark Lotkin, A Set of Test
// Matrices (1955) available at http://www.jstor.org/stable/2002051.
//
// Specifically, form the N-th test matrix as a general banded matrix A in ab,
// x = [1, ..., 1], and b = A*x computed in long double precision.
static void lotkin1955(int N,
                       double *ab,
                       double *b,
                       double *x)
{
    // The Lotkin matrices are square, but the banded routines don't care
    // provided we store them using the correct banded layout
    const int kl = N - 1, ku = kl, ldab = kl + 1 + ku;

    // First row of A where i == 0 contains only ones
    for (int j = 0; j < N; ++j) {
        ab[j*ldab+(ku+0-j)] = 1;
    }
    x[0] = 1;
    b[0] = N;

    // Second and subsequent rows of A
    // Form b_i using long doubles for extra precision
    // Reverse iteration on j sums from smallest to largest for b_i
    for (int i = 1; i < N; ++i) {
        x[i] = 1;
        long double b_i = 0;
        for (int j = N; j --> 0 ;) {
            const long double a_ij = 1.0L / (i+j+1);
            b_i += a_ij;
            ab[j*ldab+(ku+i-j)] = (double) a_ij;
        }
        b[i] = (double) b_i;
    }
}

void test_lapackext_dsgbsvx()
{
#define MAX_N (20)

    // Working storage for solving Lotkin problems
    double ab [MAX_N * (  (MAX_N - 1) + 1 + (MAX_N - 1))];
    double afb[MAX_N * (2*(MAX_N - 1) + 1 + (MAX_N - 1))];
    double b  [MAX_N], x[MAX_N], r[MAX_N];
    int    piv[MAX_N];

    enum fact_type { exact, approx };
    enum iter_type { sd, s, d };

    // Indexing is like fact[exact_type][iter_type][Lotkin matrix number]
    char   fact [2][3][MAX_N] = {};
    int    apprx[2][3][MAX_N] = {};
    int    aiter[2][3][MAX_N] = {};
    double afrob[2][3][MAX_N] = {};
    int    siter[2][3][MAX_N] = {};
    int    diter[2][3][MAX_N] = {};
    double tolsc[2][3][MAX_N] = {};
    double res  [2][3][MAX_N] = {};
    double ores [2][3][MAX_N] = {};  // Observed residual

    // Form Lotkin-based test problems and then solve them using DSGBSVX
    for (int n = 0; n < MAX_N; ++n) {
        int info;
        int i, j;  // Used to track fact_type, iter_type

        // Form the Lotkin test problem
        lotkin1955(n, ab, b, x);

        // -----------------------TEST1--------------------------
        // Compute solution forcing a factorization along the way
        // ------------------------------------------------------
        i = exact; j = sd;
        fact [i][j][n] = 'N';
        apprx[i][j][n] =   0;
        aiter[i][j][n] =   5;
        afrob[i][j][n] = - 1;
        siter[i][j][n] =  25;
        diter[i][j][n] =  10;
        tolsc[i][j][n] =   1;
        res  [i][j][n] = - 1;
        if (!n) continue;

        info = suzerain_lapackext_dsgbsvx(
                fact[i][j]+n, apprx[i][j]+n, aiter[i][j][n], 'n', n, n-1, n-1,
                ab, afrob[i][j]+n, afb, piv, b, x, siter[i][j]+n,
                diter[i][j]+n, tolsc[i][j]+n, r, res[i][j]+n);
        gsl_test(info,
                "%s:%d returns zero on Lotkin matrix %d",
                 __func__, __LINE__, n);
        gsl_test(tolsc[i][j][n] > 1,
                "%s:%d Lotkin matrix %d gives tolsc <= 1: %g",
                 __func__, __LINE__, n, tolsc[i][j][n]);

        // Check solution bests 10 times the reported tolerance using
        //     r + Ax - b = b - Ax + Ax - b = 0
        // Beware that b was found in higher precision by lotkin1955(...)
        suzerain_blas_dgbmv('N', n, n, n-1, n-1,
                            1.0, ab, (n-1)+1+(n-1), x, 1, 1, r, 1);
        suzerain_blas_daxpy(n, -1.0, b, 1, r, 1);
        ores[i][j][n] = suzerain_blas_dnrm2(n, r, 1);
        gsl_test(ores[i][j][n] > /*Empirical*/ 10*res[i][j][n],
                 "%s:%d Lotkin matrix %d residual satisfied: %g vs reported %g",
                 __func__, __LINE__, n, ores[i][j][n], res[i][j][n]);

        // ---------------------------TEST2------------------------------
        // Perturb the factorization by swapping first two pivot entries
        // This turns the factorization into mush regardless of precision
        // --------------------------------------------------------------
        { int t = piv[0]; piv[0] = piv[1 % n]; piv[1 % n] = t; }
        i = approx; j = sd;
        fact [i][j][n] = fact [exact][sd][n];
        apprx[i][j][n] =   1;
        aiter[i][j][n] =  75;
        afrob[i][j][n] = afrob[exact][sd][n];
        siter[i][j][n] =  75;
        diter[i][j][n] =  75;
        tolsc[i][j][n] =   1;
        res  [i][j][n] = - 1;

        info = suzerain_lapackext_dsgbsvx(
                fact[i][j]+n, apprx[i][j]+n, aiter[i][j][n], 'n', n, n-1, n-1,
                ab, afrob[i][j]+n, afb, piv, b, x, siter[i][j]+n,
                diter[i][j]+n, tolsc[i][j]+n, r, res[i][j]+n);
        gsl_test(info,
                "%s:%d returns zero on Lotkin matrix %d",
                 __func__, __LINE__, n);
        gsl_test(tolsc[i][j][n] > 1,
                "%s:%d Lotkin matrix %d gives tolsc <= 1: %g",
                 __func__, __LINE__, n, tolsc[i][j][n]);
        if (n != 2) { // Flipping pivots for n == 2 is too much it seems
            gsl_test(apprx[i][j][n] == 0,
                    "%s:%d matrix %d solved with perturbed factorization",
                    __func__, __LINE__, n);
        }

        // Check solution bests 10 times the reported tolerance using
        //     r + Ax - b = b - Ax + Ax - b = 0
        // Beware that b was found in higher precision by lotkin1955(...)
        suzerain_blas_dgbmv('N', n, n, n-1, n-1,
                            1.0, ab, (n-1)+1+(n-1), x, 1, 1, r, 1);
        suzerain_blas_daxpy(n, -1.0, b, 1, r, 1);
        ores[i][j][n] = suzerain_blas_dnrm2(n, r, 1);
        gsl_test(ores[i][j][n] > /*Empirical*/ 5*res[i][j][n],
                 "%s:%d Lotkin matrix %d residual satisfied: %g vs reported %g",
                 __func__, __LINE__, n, ores[i][j][n], res[i][j][n]);

        // ----------------TEST3------------------
        // Permit only double precision iterations
        // ---------------------------------------
        i = exact; j = d;
        fact [i][j][n] = 'N';
        apprx[i][j][n] =   0;
        aiter[i][j][n] =   5;
        afrob[i][j][n] = - 1;
        siter[i][j][n] = - 1;
        diter[i][j][n] =  10;
        tolsc[i][j][n] =   1;
        res  [i][j][n] = - 1;

        info = suzerain_lapackext_dsgbsvx(
                fact[i][j]+n, apprx[i][j]+n, aiter[i][j][n], 'n', n, n-1, n-1,
                ab, afrob[i][j]+n, afb, piv, b, x, siter[i][j]+n,
                diter[i][j]+n, tolsc[i][j]+n, r, res[i][j]+n);
        gsl_test(info,
                "%s:%d returns zero on Lotkin matrix %d",
                 __func__, __LINE__, n);
        gsl_test(tolsc[i][j][n] > 1,
                "%s:%d Lotkin matrix %d gives tolsc <= 1: %g",
                 __func__, __LINE__, n, tolsc[i][j][n]);
        gsl_test(fact[i][j][n] == 'S',
                "%s:%d Lotkin matrix %d reports no mixed precision in fact",
                 __func__, __LINE__, n);
        gsl_test(siter[i][j][n] != 0,
                "%s:%d Lotkin matrix %d reports no mixed precision in siter",
                 __func__, __LINE__, n);

        // Check solution bests 10 times the reported tolerance using
        //     r + Ax - b = b - Ax + Ax - b = 0
        // Beware that b was found in higher precision by lotkin1955(...)
        suzerain_blas_dgbmv('N', n, n, n-1, n-1,
                            1.0, ab, (n-1)+1+(n-1), x, 1, 1, r, 1);
        suzerain_blas_daxpy(n, -1.0, b, 1, r, 1);
        ores[i][j][n] = suzerain_blas_dnrm2(n, r, 1);
        gsl_test(ores[i][j][n] > /*Empirical*/ 10*res[i][j][n],
                 "%s:%d Lotkin matrix %d residual satisfied: %g vs reported %g",
                 __func__, __LINE__, n, ores[i][j][n], res[i][j][n]);

        // ----------------TEST4-----------------
        // Permit only mixed precision iterations
        // --------------------------------------
        i = exact; j = s;
        fact [i][j][n] = 'N';
        apprx[i][j][n] =   0;
        aiter[i][j][n] =  20;
        afrob[i][j][n] = - 1;
        siter[i][j][n] =  30;
        diter[i][j][n] = - 1;
        tolsc[i][j][n] =   1;
        res  [i][j][n] = - 1;

        info = suzerain_lapackext_dsgbsvx(
                fact[i][j]+n, apprx[i][j]+n, aiter[i][j][n], 'n', n, n-1, n-1,
                ab, afrob[i][j]+n, afb, piv, b, x, siter[i][j]+n,
                diter[i][j]+n, tolsc[i][j]+n, r, res[i][j]+n);
        gsl_test(info,
                "%s:%d returns zero on Lotkin matrix %d",
                 __func__, __LINE__, n);
        // Conditioning is just too nasty beyond n == 7 so don't require tolsc
        // However, the residual and everything else must still hold though
        if (n < 8) {
            gsl_test(tolsc[i][j][n] > 1,
                    "%s:%d Lotkin matrix %d gives tolsc <= 1: %g",
                    __func__, __LINE__, n, tolsc[i][j][n]);
        }
        gsl_test(fact[i][j][n] == 'D',
                "%s:%d Lotkin matrix %d reports no double precision in fact",
                 __func__, __LINE__, n);
        gsl_test(diter[i][j][n] != 0,
                "%s:%d Lotkin matrix %d reports no double precision in diter",
                 __func__, __LINE__, n);

        // Check solution bests 10 times the reported tolerance using
        //     r + Ax - b = b - Ax + Ax - b = 0
        // Beware that b was found in higher precision by lotkin1955(...)
        suzerain_blas_dgbmv('N', n, n, n-1, n-1,
                            1.0, ab, (n-1)+1+(n-1), x, 1, 1, r, 1);
        suzerain_blas_daxpy(n, -1.0, b, 1, r, 1);
        ores[i][j][n] = suzerain_blas_dnrm2(n, r, 1);
        gsl_test(ores[i][j][n] > /*Empirical*/ 10*res[i][j][n],
                 "%s:%d Lotkin matrix %d residual satisfied: %g vs reported %g",
                 __func__, __LINE__, n, ores[i][j][n], res[i][j][n]);

        // ----------------TEST5-----------------------------------
        // Permit only double precision iterations with approximate
        // --------------------------------------------------------
        // Unnecessary from a coverage perspective
        i = approx; j = d;

        // ----------------TEST6----------------------------------
        // Permit only mixed precision iterations with approximate
        // -------------------------------------------------------
        // Unnecessary from a coverage perspective
        i = approx; j = d;

    }

    const int BREAK_HERE_WITH_A_DEBUGGER_IF_YOU_LIKE = 123;
    (void) BREAK_HERE_WITH_A_DEBUGGER_IF_YOU_LIKE;

#undef MAX_N
}

int
main(int argc, char **argv)
{
    SUZERAIN_UNUSED(argc);
    SUZERAIN_UNUSED(argv);

    gsl_ieee_env_setup();

    test_daxpby();
    test_daxpby_nop();

    test_saxpby();
    test_saxpby_nop();

    test_caxpby();
    /* TODO test_caxpby_nop */

    test_zaxpby();
    /* TODO test_zaxpby_nop */

    test_dwaxpby();
    test_swaxpby();

    test_dswap();
    test_sswap();
    test_zswap();
    test_cswap();

    test_dscal();
    test_sscal();
    test_zscal();
    test_cscal();

    test_dcopy();
    test_dcopy_fill();
    test_scopy();
    test_zcopy();
    test_ccopy();

    test_snrm2();
    test_dnrm2();
    test_scnrm2();
    test_dznrm2();

    test_daxpy();
    test_daxpy_const();
    test_saxpy();
    test_caxpy();
    test_zaxpy();

    test_ddot();
    test_sdot();
    test_cdotc();
    test_zdotc();

    test_dgb_acc1();
    test_dgb_acc2();
    test_dgb_acc_nop();

    test_sgb_acc1();
    test_sgb_acc2();
    test_sgb_acc_nop();

    test_zgb_acc1();
    test_zgb_acc2();
    test_zgb_acc_nop();

    test_cgb_acc1();
    test_cgb_acc2();
    test_cgb_acc_nop();

    test_dgb_diag_scale_acc1(/* beta = */ 6.0);
    test_dgb_diag_scale_acc1(/* beta = */ 1.0);
    test_dgb_diag_scale_acc2(/* beta = */ 6.0);
    test_dgb_diag_scale_acc2(/* beta = */ 1.0);

    test_sgb_diag_scale_acc1(/* beta = */ 6.0);
    test_sgb_diag_scale_acc1(/* beta = */ 1.0);
    test_sgb_diag_scale_acc2(/* beta = */ 6.0);
    test_sgb_diag_scale_acc2(/* beta = */ 1.0);

    test_zaxpy_d();
    test_zaxpby_d();

    /* TODO Add suzerain_blas_{d,s}{s,g}bmv test cases                  */
    /* Exercised fairly well in test_bsplineop via                      */
    /* accumulate{,_complex} and apply{,_complex}                       */

    /* TODO Add suzerain_blasext_d{s,g}bmzv test cases                  */
    /* suzerain_blasext_d{s,g}bmzv exercised somewhat in test_bsplineop */
    /* via accumulate_complex and apply_complex                         */

    test_zgb_acc_d1();
    test_zgb_acc_d2();
    test_zgb_acc_d_nop();

    test_blasext_i2s_zaxpby2();
    test_blasext_sgbnorm1();
    test_blasext_dgbnorm1();
    test_blasext_cgbnorm1();
    test_blasext_zgbnorm1();

    test_blasext_ddemote();
    test_blasext_dpromote();
    test_blasext_zdemote();
    test_blasext_zpromote();

    /* TODO Add test_lapack_{c,z}gbtr{f,s} */
    /* Already exercised to some extent in test_bsplineop */

    test_lapack_dgbcon();

    /* TODO Add test_lapack_{s,c,z}gbcon */
    /* zgbcon already exercised to some extent in test_bsplineop */

    test_lapackext_dsgbsvx();

#ifdef SUZERAIN_HAVE_MKL
#if INTEL_MKL_VERSION < 110002
    MKL_FreeBuffers();
#else
    mkl_free_buffers();
#endif
#endif

    exit(gsl_test_summary());
}
