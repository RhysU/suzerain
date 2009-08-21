#include <suzerain/config.h>

#include <stdlib.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_test.h>
#include <suzerain/timestepper.h>

void
check_smr91_constants()
{
    const suzerain_lsrk_method m = suzerain_lsrk_smr91;

    for (int i = 0; i < m.substeps; ++i) {
        const double res = m.alpha[i] + m.beta[i] - m.gamma[i] - m.zeta[i];
        gsl_test_abs(res, 0.0, GSL_DBL_EPSILON,
                "Coefficient residual for %s substep %d of %d",
                m.name, i, m.substeps);
    }
}

void
check_smr91_scalareqn_substeps()
{
    /* Checks that operations are as expected */

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
    const double b[1]           = { 23.0 };
    const int           incb    = 1;
    const int           ldb     = 1;
    /* c declared below */
    const int           incc    = 1;
    const int           ldc     = 1;

    {
        double    a[1]    = { 19.0 };
        double    c[1]    = { 29.0 };
        const int substep = 0;
        suzerain_lsrk_substep(
                suzerain_lsrk_smr91,
                n, kl, ku,
                M, ldM,
                nD, xi, D, ldD,
                delta_t, nrhs,
                a, inca, lda,
                b, incb, ldb,
                c, incc, ldc,
                substep);
        gsl_test_abs(*a, -111557.0/4393.0, 1.0e3*GSL_DBL_EPSILON,
                "%s substep %d result", __func__, substep);
    }

    {
        double    a[1]    = { 19.0 };
        double    c[1]    = { 29.0 };
        const int substep = 1;
        suzerain_lsrk_substep(
                suzerain_lsrk_smr91,
                n, kl, ku,
                M, ldM,
                nD, xi, D, ldD,
                delta_t, nrhs,
                a, inca, lda,
                b, incb, ldb,
                c, incc, ldc,
                substep);
        gsl_test_abs(*a, 80129.0/11870.0, 1.0e3*GSL_DBL_EPSILON,
                "%s substep %d result", __func__, substep);
    }

    {
        double    a[1]    = { 19.0 };
        double    c[1]    = { 29.0 };
        const int substep = 2;
        suzerain_lsrk_substep(
                suzerain_lsrk_smr91,
                n, kl, ku,
                M, ldM,
                nD, xi, D, ldD,
                delta_t, nrhs,
                a, inca, lda,
                b, incb, ldb,
                c, incc, ldc,
                substep);
        gsl_test_abs(*a, -36817.0/1898.0, 1.0e3*GSL_DBL_EPSILON,
                "%s substep %d result", __func__, substep);
    }
}

void
check_smr91_matrixeqn_substeps()
{
    /* Checks band storage, increments, leading dimensions, etc. */

    const int           n       = 2;
    const int           kl      = 1;
    const int           ku      = 0;
    const double        M[4]    = { 2.0, 3.0, 5.0, /*DK*/555 };
    const int           ldM     = kl + 1 + ku;
    const int           nD      = 1;
    const double        xi[1]   = { 7.0 };
    const double        D0[4]   = { 11.0, 13.0, 17.0, /*DK*/555 };
    const double *const D[1]    = { D0 };
    const int           ldD     = kl + 1 + ku;
    double              delta_t = 19.0;
    const int           nrhs    = 2;
    /* a declared below */
    const double b[]            = {
        23.0, /*DK*/555, 29.0, /*DK*/555, /*DK*/555,
        31.0, /*DK*/555, 37.0, /*DK*/555, /*DK*/555
    };
    const int           incb    = 2;
    const int           ldb     = 5;
    /* c declared below */

    {
        double    a[]     = {
            41.0, 43.0, /*DK*/555,
            47.0, 53.0, /*DK*/555
        };
        const int inca    = 1;
        const int lda     = 3;
        double    c[]     = {
            59.0, 61.0,
            67.0, 71.0
        };
        const int incc    = 1;
        const int ldc     = 2;
        const int substep = 0;
        suzerain_lsrk_substep(
                suzerain_lsrk_smr91,
                n, kl, ku,
                M, ldM,
                nD, xi, D, ldD,
                delta_t, nrhs,
                a, inca, lda,
                b, incb, ldb,
                c, incc, ldc,
                substep);

        double    a_expected[]     = {
            -8960639.0/161433.0, -807513949909.0/13375854081.0, /*DK*/555,
            -1146337.0/ 17937.0, -110806751027.0/ 1486206009.0, /*DK*/555
        };

        for (int i = 0; i < sizeof(a_expected)/sizeof(a_expected[0]); ++i) {
            gsl_test_abs(a[i], a_expected[i], 1.0e3*GSL_DBL_EPSILON,
                    "%s substep %d result a[%d]", __func__, substep, i);
        }
    }

    {
        double    a[]     = {
            41.0, 43.0, /*DK*/555,
            47.0, 53.0, /*DK*/555
        };
        const int inca    = 1;
        const int lda     = 3;
        double    c[]     = {
            59.0, 61.0,
            67.0, 71.0
        };
        const int incc    = 1;
        const int ldc     = 2;
        const int substep = 1;
        suzerain_lsrk_substep(
                suzerain_lsrk_smr91,
                n, kl, ku,
                M, ldM,
                nD, xi, D, ldD,
                delta_t, nrhs,
                a, inca, lda,
                b, incb, ldb,
                c, incc, ldc,
                substep);

        double    a_expected[]     = {
            112507.0/ 7267.0, 1332605257.0/ 81281395.0, /*DK*/555,
            635233.0/36335.0, 8043298231.0/406406975.0, /*DK*/555
        };

        for (int i = 0; i < sizeof(a_expected)/sizeof(a_expected[0]); ++i) {
            gsl_test_abs(a[i], a_expected[i], 1.0e3*GSL_DBL_EPSILON,
                    "%s substep %d result a[%d]", __func__, substep, i);
        }
    }

    {
        double    a[]     = {
            41.0, 43.0, /*DK*/555,
            47.0, 53.0, /*DK*/555
        };
        const int inca    = 1;
        const int lda     = 3;
        double    c[]     = {
            59.0, 61.0,
            67.0, 71.0
        };
        const int incc    = 1;
        const int ldc     = 2;
        const int substep = 2;
        suzerain_lsrk_substep(
                suzerain_lsrk_smr91,
                n, kl, ku,
                M, ldM,
                nD, xi, D, ldD,
                delta_t, nrhs,
                a, inca, lda,
                b, incb, ldb,
                c, incc, ldc,
                substep);

        double    a_expected[]     = {
            -58803.0/1451.0, -139589509.0/3237181.0, /*DK*/555,
            -68261.0/1451.0, -174698270.0/3237181.0, /*DK*/555
        };

        for (int i = 0; i < sizeof(a_expected)/sizeof(a_expected[0]); ++i) {
            gsl_test_abs(a[i], a_expected[i], 1.0e3*GSL_DBL_EPSILON,
                    "%s substep %d result a[%d]", __func__, substep, i);
        }
    }
}

int
main(int argc, char **argv)
{
    gsl_ieee_env_setup();

    check_smr91_constants();
    check_smr91_scalareqn_substeps();
    check_smr91_matrixeqn_substeps();

    exit(gsl_test_summary());
}
