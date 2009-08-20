#include <suzerain/config.h>

#include <stdlib.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_test.h>
#include <suzerain/timestepper.h>


void
check_scalar_substeps()
{
    const int           n       = 1;
    const int           kl      = 0;
    const int           ku      = 0;
    const double        M[1]    = { 1.0 };
    const int           ldM     = 1;
    const int           nD      = 3;
    const double        xi[3]   = { 2.0, 3.0, 5.0 };
    const double        D0[1]   = { 7.0 };
    const double        D1[1]   = { 11.0 };
    const double        D2[1]   = { 13.0 };
    const double *const D[3]    = { D0, D1, D2 };
    const int           ldD     = 1;
    double              delta_t = 17.0;
    const int           nrhs    = 1;
    /* a declared below */
    const int           inca    = 1;
    const int           lda     = 1;
    const double b[1]    = { 23.0 };
    const int           incb    = 1;
    const int           ldb     = 1;
    /* c declared below */
    const int           incc    = 1;
    const int           ldc     = 1;

    {
        double    a[1]    = { 19.0 };
        double    c[1]    = { 29.0 };
        const int substep = 0;
        suzerain_smr91_substep(
                n, kl, ku,
                M, ldM,
                nD, xi, D, ldD,
                delta_t, nrhs,
                a, inca, lda,
                b, incb, ldb,
                c, incc, ldc,
                substep);
        gsl_test_abs(*a, -111557.0/4393.0, 1.0e2*GSL_DBL_EPSILON,
                "%s substep %d result", __func__, substep);
    }

    {
        double    a[1]    = { 19.0 };
        double    c[1]    = { 29.0 };
        const int substep = 1;
        suzerain_smr91_substep(
                n, kl, ku,
                M, ldM,
                nD, xi, D, ldD,
                delta_t, nrhs,
                a, inca, lda,
                b, incb, ldb,
                c, incc, ldc,
                substep);
        gsl_test_abs(*a, 80129.0/11870.0, 1.0e2*GSL_DBL_EPSILON,
                "%s substep %d result", __func__, substep);
    }

    {
        double    a[1]    = { 19.0 };
        double    c[1]    = { 29.0 };
        const int substep = 2;
        suzerain_smr91_substep(
                n, kl, ku,
                M, ldM,
                nD, xi, D, ldD,
                delta_t, nrhs,
                a, inca, lda,
                b, incb, ldb,
                c, incc, ldc,
                substep);
        gsl_test_abs(*a, -36817.0/1898.0, 1.0e2*GSL_DBL_EPSILON,
                "%s substep %d result", __func__, substep);
    }
}

int
main(int argc, char **argv)
{
    gsl_ieee_env_setup();

    check_scalar_substeps();

    exit(gsl_test_summary());
}
