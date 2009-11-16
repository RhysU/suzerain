#include <suzerain/config.h>
#include <suzerain/common.h>
#pragma hdrstop
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_test.h>
#include <suzerain/blas_et_al.h>

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

    const double nb        = sizeof(b_data)/sizeof(b_data[0]);
    const double nexpected = sizeof(expected_data)/sizeof(expected_data[0]);
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blas_dgb_acc( m, n, kl, ku, alpha, a, lda, beta, b, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i], expected_data[i], GSL_DBL_EPSILON,
                "dgb_acc index %d", i);
    }
}

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
    const double *expected_data = b_data;
    const double *a = a_data + lda-(ku+1+kl);
    double       *b = b_data + ldb-(ku+1+kl);

    const double nb        = sizeof(b_data)/sizeof(b_data[0]);
    const double nexpected = nb;
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blas_dgb_acc( m, n, kl, ku, alpha, a, lda, beta, b, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i], expected_data[i], GSL_DBL_EPSILON,
                "dgb_acc index %d", i);
    }
}

int
main(int argc, char **argv)
{
    gsl_ieee_env_setup();

    test_daxpby();
    test_daxpby_nop();

    test_dwaxpby();

    test_dgb_acc();
    test_dgb_acc_nop();

    exit(gsl_test_summary());
}
