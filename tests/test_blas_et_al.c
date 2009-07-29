#include <suzerain/config.h>

#include <stdlib.h>
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

int
main(int argc, char **argv)
{
    gsl_ieee_env_setup();

    test_daxpby();

    exit(gsl_test_summary());
}
