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

#include <suzerain/blas_et_al/blas.h>

#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_test.h>

#ifdef SUZERAIN_HAVE_MKL
#include <mkl.h>
#include <mkl_service.h>
#endif

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

    test_zaxpy_d();
    test_zaxpby_d();

    /* TODO Add suzerain_blas_{d,s}{s,g}bmv test cases                  */
    /* Exercised fairly well in test_bsplineop via                      */
    /* accumulate{,_complex} and apply{,_complex}                       */

    test_zgb_acc_d1();
    test_zgb_acc_d2();
    test_zgb_acc_d_nop();

#ifdef SUZERAIN_HAVE_MKL
#if INTEL_MKL_VERSION < 110002
    MKL_FreeBuffers();
#else
    mkl_free_buffers();
#endif
#endif

    exit(gsl_test_summary());
}
