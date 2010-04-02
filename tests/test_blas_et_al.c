#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.h>
#pragma hdrstop
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_test.h>
#include <suzerain/blas_et_al.h>

#ifdef SUZERAIN_HAVE_MKL
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
        gsl_test_abs(y[i], expected[i], GSL_DBL_EPSILON, "daxpby index %d", i);
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
        gsl_test_abs(y[i], expected[i], GSL_FLT_EPSILON, "saxpby index %d", i);
    }
}

static
void
test_caxpby()
{
    int i;

    const float alpha[2] = {2., -0.25};
    const float x[][2]   = {{1., -1.}, {2., -2.}, {3., -3.}};
    const int   incx     = 1;
    const float beta[2]  = {3., -0.50};
    float       y[][2]   = {{4.,-4.}, {-1,1},{5.,-5},{-2,2},{6.,-6},{-3,3}};
    const int   incy     = 2;
    const int   nx       = sizeof(x)/sizeof(x[0]);
    const int   ny       = sizeof(y)/sizeof(y[0]);
    const float  expected[][2] = {
        {
            alpha[0]*x[0][0]-alpha[1]*x[0][1]+beta[0]*y[0][0]-beta[1]*y[0][1],
            alpha[0]*x[0][1]+alpha[1]*x[0][0]+beta[0]*y[0][1]+beta[1]*y[0][0]
        },
        {
            y[1][0],
            y[1][1]
        },
        {
            alpha[0]*x[1][0]-alpha[1]*x[1][1]+beta[0]*y[2][0]-beta[1]*y[2][1],
            alpha[0]*x[1][1]+alpha[1]*x[1][0]+beta[0]*y[2][1]+beta[1]*y[2][0]
        },
        {
            y[3][0],
            y[3][1]
        },
        {
            alpha[0]*x[2][0]-alpha[1]*x[2][1]+beta[0]*y[4][0]-beta[1]*y[4][1],
            alpha[0]*x[2][1]+alpha[1]*x[2][0]+beta[0]*y[4][1]+beta[1]*y[4][0]
        },
        {
            y[5][0],
            y[5][1]
        }
    };
    const int nexpected = sizeof(expected)/sizeof(expected[0]);

    gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
    gsl_test_int(ny, nexpected, "Expected results' length");

    suzerain_blas_caxpby(nx/incx, alpha, x, incx, beta, y, incy);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(y[i][0], expected[i][0], GSL_FLT_EPSILON,
                "caxpby real index %d", i);
        gsl_test_abs(y[i][1], expected[i][1], GSL_FLT_EPSILON,
                "caxpby imag index %d", i);
    }
}

static
void
test_zaxpby()
{
    int i;

    const double alpha[2] = {2., -0.25};
    const double x[][2]   = {{1., -1.}, {2., -2.}, {3., -3.}};
    const int    incx     = 1;
    const double beta[2]  = {3., -0.50};
    double       y[][2]   = {{4.,-4.},{-1,1},{5.,-5},{-2,2},{6.,-6},{-3,3}};
    const int    incy     = 2;
    const int    nx       = sizeof(x)/sizeof(x[0]);
    const int    ny       = sizeof(y)/sizeof(y[0]);
    const double expected[][2] = {
        {
            alpha[0]*x[0][0]-alpha[1]*x[0][1]+beta[0]*y[0][0]-beta[1]*y[0][1],
            alpha[0]*x[0][1]+alpha[1]*x[0][0]+beta[0]*y[0][1]+beta[1]*y[0][0]
        },
        {
            y[1][0],
            y[1][1]
        },
        {
            alpha[0]*x[1][0]-alpha[1]*x[1][1]+beta[0]*y[2][0]-beta[1]*y[2][1],
            alpha[0]*x[1][1]+alpha[1]*x[1][0]+beta[0]*y[2][1]+beta[1]*y[2][0]
        },
        {
            y[3][0],
            y[3][1]
        },
        {
            alpha[0]*x[2][0]-alpha[1]*x[2][1]+beta[0]*y[4][0]-beta[1]*y[4][1],
            alpha[0]*x[2][1]+alpha[1]*x[2][0]+beta[0]*y[4][1]+beta[1]*y[4][0]
        },
        {
            y[5][0],
            y[5][1]
        }
    };
    const int nexpected = sizeof(expected)/sizeof(expected[0]);

    gsl_test_int(nx/incx, ny/incy, "Vectors of equivalent lengths");
    gsl_test_int(ny, nexpected, "Expected results' length");

    suzerain_blas_zaxpby(nx/incx, alpha, x, incx, beta, y, incy);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(y[i][0], expected[i][0], GSL_DBL_EPSILON,
                "zaxpby real index %d", i);
        gsl_test_abs(y[i][1], expected[i][1], GSL_DBL_EPSILON,
                "zaxpby imag index %d", i);
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
        gsl_test_abs(w[i], expected[i], GSL_DBL_EPSILON, "dwaxpby index %d", i);
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
        gsl_test_abs(w[i], expected[i], GSL_FLT_EPSILON, "swaxpby index %d", i);
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
        gsl_test_abs(y[i], expected[i], GSL_DBL_EPSILON, "daxpby index %d", i);
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
        gsl_test_abs(y[i], expected[i], GSL_FLT_EPSILON, "saxpby index %d", i);
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
                "dswap x index %d", i);
    }
    for (i = 0; i < ny_expected; ++i) {
        gsl_test_abs(y[i], y_expected[i], GSL_DBL_EPSILON,
                "dswap y index %d", i);
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
                "sswap x index %d", i);
    }
    for (i = 0; i < ny_expected; ++i) {
        gsl_test_abs(y[i], y_expected[i], GSL_DBL_EPSILON,
                "sswap y index %d", i);
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

    suzerain_blas_cswap(nx/incx, x, incx, y, incy);
    for (i = 0; i < nx_expected; ++i) {
        gsl_test_abs(x[i][0], x_expected[i][0], GSL_FLT_EPSILON,
                "cswap x real index %d", i);
        gsl_test_abs(x[i][1], x_expected[i][1], GSL_FLT_EPSILON,
                "cswap x imag index %d", i);
    }
    for (i = 0; i < ny_expected; ++i) {
        gsl_test_abs(y[i][0], y_expected[i][0], GSL_FLT_EPSILON,
                "cswap y real index %d", i);
        gsl_test_abs(y[i][1], y_expected[i][1], GSL_FLT_EPSILON,
                "cswap y imag index %d", i);
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

    suzerain_blas_zswap(nx/incx, x, incx, y, incy);
    for (i = 0; i < nx_expected; ++i) {
        gsl_test_abs(x[i][0], x_expected[i][0], GSL_DBL_EPSILON,
                "zswap x real index %d", i);
        gsl_test_abs(x[i][1], x_expected[i][1], GSL_DBL_EPSILON,
                "zswap x imag index %d", i);
    }
    for (i = 0; i < ny_expected; ++i) {
        gsl_test_abs(y[i][0], y_expected[i][0], GSL_DBL_EPSILON,
                "zswap y real index %d", i);
        gsl_test_abs(y[i][1], y_expected[i][1], GSL_DBL_EPSILON,
                "zswap y imag index %d", i);
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
                    "dscal index %d", i);
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
                    "dscal index %d", i);
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
                    "sscal index %d", i);
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
                    "sscal index %d", i);
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

    suzerain_blas_zscal(nx/incx, alpha, x, incx);

    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(x[i][0], expected[i][0], GSL_DBL_EPSILON,
                "zscal real index %d", i);
        gsl_test_abs(x[i][1], expected[i][1], GSL_DBL_EPSILON,
                "zscal imag index %d", i);
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

    suzerain_blas_cscal(nx/incx, alpha, x, incx);

    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(x[i][0], expected[i][0], GSL_FLT_EPSILON,
                "zscal real index %d", i);
        gsl_test_abs(x[i][1], expected[i][1], GSL_FLT_EPSILON,
                "zscal imag index %d", i);
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
                "dcopy y index %d", i);
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
                "scopy y index %d", i);
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

    suzerain_blas_zcopy(nx/incx, x, incx, y, incy);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(y[i][0], expected[i][0], GSL_DBL_EPSILON,
                "zcopy real index %d", i);
        gsl_test_abs(y[i][1], expected[i][1], GSL_DBL_EPSILON,
                "zcopy imag index %d", i);
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

    suzerain_blas_ccopy(nx/incx, x, incx, y, incy);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(y[i][0], expected[i][0], GSL_FLT_EPSILON,
                "ccopy real index %d", i);
        gsl_test_abs(y[i][1], expected[i][1], GSL_FLT_EPSILON,
                "ccopy imag index %d", i);
    }
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
        gsl_test_abs(y[i], expected[i], GSL_DBL_EPSILON, "daxpy index %d", i);
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
        gsl_test_abs(y[i], expected[i], GSL_FLT_EPSILON, "saxpy index %d", i);
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

    suzerain_blas_zaxpy(nx/incx, alpha, x, incx, y, incy);

    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(y[i][0], expected[i][0], GSL_DBL_EPSILON,
                "zaxpy real index %d", i);
        gsl_test_abs(y[i][1], expected[i][1], GSL_DBL_EPSILON,
                "zaxpy imag index %d", i);
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

    suzerain_blas_caxpy(nx/incx, alpha, x, incx, y, incy);

    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(y[i][0], expected[i][0], GSL_FLT_EPSILON,
                "caxpy real index %d", i);
        gsl_test_abs(y[i][1], expected[i][1], GSL_FLT_EPSILON,
                "caxpy imag index %d", i);
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
    gsl_test_abs(result, expected, GSL_DBL_EPSILON, "ddot result");
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
    gsl_test_abs(result, expected, GSL_DBL_EPSILON, "ddot result");
}

static
void
test_dgb_acc()
{
    int i;

    const int    m     = 4;
    const int    n     = m;
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
                 -1.1,     10.0,    11.0,    12.0,     -1.0
    };
    double b_data[]  = {
        /*ldb buffer*/   /*ku2*/  /*ku1*/ /*diag*/   /*kl1*/
         -2.2,   -1.1,     -1.0,    -1.0,     1.0,      2.0,
         -2.2,   -1.1,     -1.0,     3.0,     4.0,      5.0,
         -2.2,   -1.1,      6.0,     7.0,     8.0,      9.0,
         -2.2,   -1.1,     10.0,    11.0,    12.0,     -1.0
    };
    const double expected_data[]  = {
        /*ldb buffer*/   /*ku2*/  /*ku1*/ /*diag*/   /*kl1*/
         -2.2,   -1.1,     -5.0,    -5.0,     5.0,     10.0,
         -2.2,   -1.1,     -5.0,    15.0,    20.0,     25.0,
         -2.2,   -1.1,     30.0,    35.0,    40.0,     45.0,
         -2.2,   -1.1,     50.0,    55.0,    60.0,     -5.0
    };
    const double *a = a_data + lda-(ku+1+kl);
    double       *b = b_data + ldb-(ku+1+kl);

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = sizeof(expected_data)/sizeof(expected_data[0]);
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blas_dgb_acc( m, n, kl, ku, alpha, a, lda, beta, b, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i], expected_data[i], GSL_DBL_EPSILON,
                "dgb_acc index %d", i);
    }
}

static
void
test_sgb_acc()
{
    int i;

    const int    m     = 4;
    const int    n     = m;
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
         -2.2,   -1.1,     -5.0,    -5.0,     5.0,     10.0,
         -2.2,   -1.1,     -5.0,    15.0,    20.0,     25.0,
         -2.2,   -1.1,     30.0,    35.0,    40.0,     45.0,
         -2.2,   -1.1,     50.0,    55.0,    60.0,     -5.0
    };
    const float  *a = a_data + lda-(ku+1+kl);
    float        *b = b_data + ldb-(ku+1+kl);

    const float  nb        = sizeof(b_data)/sizeof(b_data[0]);
    const float  nexpected = sizeof(expected_data)/sizeof(expected_data[0]);
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blas_sgb_acc( m, n, kl, ku, alpha, a, lda, beta, b, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i], expected_data[i], GSL_FLT_EPSILON,
                "dgb_acc index %d", i);
    }
}

static
void
test_cgb_acc()
{
    int i;

    const int    m        = 4;
    const int    n        = m;
    const int    ku       = 2;
    const int    kl       = 1;
    const int    lda      = 5;
    const int    ldb      = 6;
    const float  alpha[2] = {2.0, 0.0}; /* TODO Nontrivial imaginary part */
    const float  beta[2]  = {3.0, 0.0}; /* TODO Nontrivial imaginary part */

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
        {-2.2,-2.2}, {-1.1,-1.1}, {-5,-5}, {-5,-5}, { 5, 5}, {10,10},
        {-2.2,-2.2}, {-1.1,-1.1}, {-5,-5}, {15,15}, {20,20}, {25,25},
        {-2.2,-2.2}, {-1.1,-1.1}, {30,30}, {35,35}, {40,40}, {45,45},
        {-2.2,-2.2}, {-1.1,-1.1}, {50,50}, {55,55}, {60,60}, {-5,-5}
    };
    const float  (*a)[2] = a_data + lda-(ku+1+kl);
    float        (*b)[2] = b_data + ldb-(ku+1+kl);

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = sizeof(expected_data)/sizeof(expected_data[0]);
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blas_cgb_acc( m, n, kl, ku, alpha, a, lda, beta, b, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i][0], expected_data[i][0], GSL_FLT_EPSILON,
                "cgb_acc real index %d", i);
        gsl_test_abs(b_data[i][1], expected_data[i][1], GSL_FLT_EPSILON,
                "cgb_acc imag index %d", i);
    }
}

static
void
test_zgb_acc()
{
    int i;

    const int    m        = 4;
    const int    n        = m;
    const int    ku       = 2;
    const int    kl       = 1;
    const int    lda      = 5;
    const int    ldb      = 6;
    const double alpha[2] = {2.0, 0.0}; /* TODO Nontrivial imaginary part */
    const double beta[2]  = {3.0, 0.0}; /* TODO Nontrivial imaginary part */

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
        {-2.2,-2.2}, {-1.1,-1.1}, {-5,-5}, {-5,-5}, { 5, 5}, {10,10},
        {-2.2,-2.2}, {-1.1,-1.1}, {-5,-5}, {15,15}, {20,20}, {25,25},
        {-2.2,-2.2}, {-1.1,-1.1}, {30,30}, {35,35}, {40,40}, {45,45},
        {-2.2,-2.2}, {-1.1,-1.1}, {50,50}, {55,55}, {60,60}, {-5,-5}
    };
    const double (*a)[2] = a_data + lda-(ku+1+kl);
    double       (*b)[2] = b_data + ldb-(ku+1+kl);

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = sizeof(expected_data)/sizeof(expected_data[0]);
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blas_zgb_acc( m, n, kl, ku, alpha, a, lda, beta, b, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i][0], expected_data[i][0], GSL_FLT_EPSILON,
                "zgb_acc real index %d", i);
        gsl_test_abs(b_data[i][1], expected_data[i][1], GSL_FLT_EPSILON,
                "zgb_acc imag index %d", i);
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
                "dgb_acc index %d", i);
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
                "sgb_acc index %d", i);
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
    const float alpha[2] = {0.0, 0.0}; /* NOP zero */
    const float beta[2]  = {1.0, 0.0}; /* NOP one */

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

    suzerain_blas_cgb_acc( m, n, kl, ku, alpha, a, lda, beta, b, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i][0], expected_data[i][0], GSL_FLT_EPSILON,
                "cgb_acc real index %d", i);
        gsl_test_abs(b_data[i][1], expected_data[i][1], GSL_FLT_EPSILON,
                "cgb_acc imag index %d", i);
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
    const double alpha[2] = {0.0, 0.0}; /* NOP zero */
    const double beta[2]  = {1.0, 0.0}; /* NOP one */

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

    suzerain_blas_zgb_acc( m, n, kl, ku, alpha, a, lda, beta, b, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i][0], expected_data[i][0], GSL_FLT_EPSILON,
                "zgb_acc real index %d", i);
        gsl_test_abs(b_data[i][1], expected_data[i][1], GSL_FLT_EPSILON,
                "zgb_acc imag index %d", i);
    }

    free(expected_data);
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
                        "i2s_zaxpby2 real part of index %d at line %d",
                        i, __LINE__);
            gsl_test_abs(z_im[i], i+1.5, GSL_DBL_EPSILON,
                        "i2s_zaxpby2 imag part of index %d at line %d",
                        i, __LINE__);
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
                        "i2s_zaxpby2 real part of index %d at line %d",
                        i, __LINE__);
                gsl_test_abs(z_im[i*incz_im+j*ldz_im], i+1.5, GSL_DBL_EPSILON,
                        "i2s_zaxpby2 imag part of index %d at line %d",
                        i, __LINE__);
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
                        "i2s_zaxpby2 real part of index %d at line %d",
                        i, __LINE__);
            gsl_test_abs(z_im[i], 2*(i+1.5), GSL_DBL_EPSILON,
                        "i2s_zaxpby2 imag part of index %d at line %d",
                        i, __LINE__);
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
                        "i2s_zaxpby2 real part of index %d at line %d",
                        i, __LINE__);
                gsl_test_abs(z_im[i*incz_im+j*ldz_im],
                        2*(i+1.5), GSL_DBL_EPSILON,
                        "i2s_zaxpby2 imag part of index %d at line %d",
                        i, __LINE__);
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
                        "i2s_zaxpby2 real part of index %d at line %d",
                        i, __LINE__);
            gsl_test_abs(z_im[i], 2*i*i + 16*i + 85.5,
                        GSL_DBL_EPSILON,
                        "i2s_zaxpby2 imag part of index %d at line %d",
                        i, __LINE__);
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
                        "i2s_zaxpby2 real part of index %d at line %d",
                        i, __LINE__);
                gsl_test_abs(z_im[i*incz_im+j*ldz_im], 2*i*i + 16*i + 85.5,
                        GSL_DBL_EPSILON,
                        "i2s_zaxpby2 imag part of index %d at line %d",
                        i, __LINE__);
            }
        }
    }
}

int
main(int argc, char **argv)
{
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
    test_scopy();
    test_zcopy();
    test_ccopy();

    test_daxpy();
    test_saxpy();
    test_caxpy();
    test_zaxpy();

    test_ddot();
    test_sdot();

    test_dgb_acc();
    test_dgb_acc_nop();

    test_sgb_acc();
    test_sgb_acc_nop();

    test_zgb_acc();
    test_zgb_acc_nop();

    test_cgb_acc();
    test_cgb_acc_nop();

    /* TODO Add suzerain_blasext_dgbmzv test cases                    */
    /* suzerain_blasext_dgbmzv exercised somewhat in test_bspline via */
    /* zaccumulate and zapply                                         */

    test_blasext_i2s_zaxpby2();

    /* TODO Add test_{c,z}gbtr{f,s} */
    /* Already exercised to some extent in test_bspline */

#ifdef SUZERAIN_HAVE_MKL
    MKL_FreeBuffers();
#endif

    exit(gsl_test_summary());
}
