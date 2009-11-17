#include <suzerain/config.h>
#include <suzerain/common.h>
#pragma hdrstop
#include <suzerain/blas_et_al.h>
#include <suzerain/timestepper.h>
#include <suzerain/richardson.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_vector.h>

void
check_smr91_constants()
{
    const suzerain_lsrk_method * const m = suzerain_lsrk_smr91;

    for (int i = 0; i < m->substeps; ++i) {
        const double res = m->alpha[i] + m->beta[i] - m->gamma[i] - m->zeta[i];
        gsl_test_abs(res, 0.0, GSL_DBL_EPSILON,
                "Coefficient residual for %s substep %d of %d",
                m->name, i, m->substeps);
    }
}

void
check_smr91_imex_scalareqn_substeps()
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
        suzerain_lsrk_imex_substep(
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
        suzerain_lsrk_imex_substep(
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
        suzerain_lsrk_imex_substep(
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
check_smr91_ex_scalareqn_substeps()
{
    /* Checks that operations are as expected */

    const int           n       = 1;
    const int           kl      = 0;
    const int           ku      = 0;
    double              M[1]    = { 1.0 };
    const int           ldM     = 1;
    int                 ipivM[1];
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

    /* Factor M prior to invoking any substeps */
    suzerain_lapack_dgbtrf(n, n, kl, ku, M, ldM, ipivM);

    {
        double    a[1]    = { 19.0 };
        double    c[1]    = { 29.0 };
        const int substep = 0;
        suzerain_lsrk_ex_substep(
                suzerain_lsrk_smr91,
                n, kl, ku,
                M, ldM, ipivM,
                delta_t, nrhs,
                a, inca, lda,
                b, incb, ldb,
                c, incc, ldc,
                substep);
        gsl_test_abs(*a, 3413.0/15.0, 1.0e3*GSL_DBL_EPSILON,
                "%s substep %d result", __func__, substep);
    }

    {
        double    a[1]    = { 19.0 };
        double    c[1]    = { 29.0 };
        const int substep = 1;
        suzerain_lsrk_ex_substep(
                suzerain_lsrk_smr91,
                n, kl, ku,
                M, ldM, ipivM,
                delta_t, nrhs,
                a, inca, lda,
                b, incb, ldb,
                c, incc, ldc,
                substep);
        gsl_test_abs(*a, 1267.0/30.0, 1.0e3*GSL_DBL_EPSILON,
                "%s substep %d result", __func__, substep);
    }

    {
        double    a[1]    = { 19.0 };
        double    c[1]    = { 29.0 };
        const int substep = 2;
        suzerain_lsrk_ex_substep(
                suzerain_lsrk_smr91,
                n, kl, ku,
                M, ldM, ipivM,
                delta_t, nrhs,
                a, inca, lda,
                b, incb, ldb,
                c, incc, ldc,
                substep);
        gsl_test_abs(*a, 641.0/6.0, 1.0e3*GSL_DBL_EPSILON,
                "%s substep %d result", __func__, substep);
    }
}

void
check_smr91_imex_matrixeqn_substeps()
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
        suzerain_lsrk_imex_substep(
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
        suzerain_lsrk_imex_substep(
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
        suzerain_lsrk_imex_substep(
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

void
check_smr91_ex_matrixeqn_substeps()
{
    /* Checks band storage, increments, leading dimensions, etc. */

    const int           n       = 2;
    const int           kl      = 1;
    const int           ku      = 0;
    /* Additional space required to perform LU using GBTRF */
    double              M[6]
        = { /*DK*/555, 2.0, 3.0, /*DK*/555, 5.0, /*DK*/555};
    const int           ldM     = 2*kl + 1 + ku;
    int                 ipivM[2];
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

    /* Factor M prior to invoking any substeps */
    suzerain_lapack_dgbtrf(n, n, kl, ku, M, ldM, ipivM);

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
        suzerain_lsrk_ex_substep(
                suzerain_lsrk_smr91,
                n, kl, ku,
                M, ldM, ipivM,
                delta_t, nrhs,
                a, inca, lda,
                b, incb, ldb,
                c, incc, ldc,
                substep);

        double    a_expected[]     = {
            2363.0/15.0, 2389.0/75.0, /*DK*/555,
            3061.0/15.0, 2531.0/75.0, /*DK*/555
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
        suzerain_lsrk_ex_substep(
                suzerain_lsrk_smr91,
                n, kl, ku,
                M, ldM, ipivM,
                delta_t, nrhs,
                a, inca, lda,
                b, incb, ldb,
                c, incc, ldc,
                substep);

        double    a_expected[]     = {
            -803.0/30.0,  639.0/10.0, /*DK*/555,
            -319.0/30.0, 1743.0/25.0, /*DK*/555
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
        suzerain_lsrk_ex_substep(
                suzerain_lsrk_smr91,
                n, kl, ku,
                M, ldM, ipivM,
                delta_t, nrhs,
                a, inca, lda,
                b, incb, ldb,
                c, incc, ldc,
                substep);

        double    a_expected[]     = {
            -86.0/3.0, 1063.0/15.0, /*DK*/555,
              8.0/3.0, 2179.0/30.0, /*DK*/555
        };

        for (int i = 0; i < sizeof(a_expected)/sizeof(a_expected[0]); ++i) {
            gsl_test_abs(a[i], a_expected[i], 1.0e3*GSL_DBL_EPSILON,
                    "%s substep %d result a[%d]", __func__, substep, i);
        }
    }
}

void
check_smr91_ex_matrixeqn_substeps_contiguous()
{
    /* Checks correct computation under contiguous memory optimizations */

    const int           n       = 2;
    const int           kl      = 1;
    const int           ku      = 0;
    /* Additional space required to perform LU using GBTRF */
    double              M[6]
        = { /*DK*/555, 2.0, 3.0, /*DK*/555, 5.0, /*DK*/555};
    const int           ldM     = 2*kl + 1 + ku;
    int                 ipivM[2];
    double              delta_t = 19.0;
    const int           nrhs    = 2;
    /* a declared below */
    const double b[]            = {
        23.0,  29.0,
        31.0,  37.0
    };
    const int           incb    = 1;
    const int           ldb     = 2;
    /* c declared below */

    /* Factor M prior to invoking any substeps */
    suzerain_lapack_dgbtrf(n, n, kl, ku, M, ldM, ipivM);

    {
        double    a[]     = {
            41.0, 43.0,
            47.0, 53.0
        };
        const int inca    = 1;
        const int lda     = 2;
        double    c[]     = {
            59.0, 61.0,
            67.0, 71.0
        };
        const int incc    = 1;
        const int ldc     = 2;
        const int substep = 0;
        suzerain_lsrk_ex_substep(
                suzerain_lsrk_smr91,
                n, kl, ku,
                M, ldM, ipivM,
                delta_t, nrhs,
                a, inca, lda,
                b, incb, ldb,
                c, incc, ldc,
                substep);

        double    a_expected[]     = {
            2363.0/15.0, 2389.0/75.0,
            3061.0/15.0, 2531.0/75.0
        };

        for (int i = 0; i < sizeof(a_expected)/sizeof(a_expected[0]); ++i) {
            gsl_test_abs(a[i], a_expected[i], 1.0e3*GSL_DBL_EPSILON,
                    "%s substep %d result a[%d]", __func__, substep, i);
        }
    }

    {
        double    a[]     = {
            41.0, 43.0,
            47.0, 53.0
        };
        const int inca    = 1;
        const int lda     = 2;
        double    c[]     = {
            59.0, 61.0,
            67.0, 71.0
        };
        const int incc    = 1;
        const int ldc     = 2;
        const int substep = 1;
        suzerain_lsrk_ex_substep(
                suzerain_lsrk_smr91,
                n, kl, ku,
                M, ldM, ipivM,
                delta_t, nrhs,
                a, inca, lda,
                b, incb, ldb,
                c, incc, ldc,
                substep);

        double    a_expected[]     = {
            -803.0/30.0,  639.0/10.0,
            -319.0/30.0, 1743.0/25.0
        };

        for (int i = 0; i < sizeof(a_expected)/sizeof(a_expected[0]); ++i) {
            gsl_test_abs(a[i], a_expected[i], 1.0e3*GSL_DBL_EPSILON,
                    "%s substep %d result a[%d]", __func__, substep, i);
        }
    }

    {
        double    a[]     = {
            41.0, 43.0,
            47.0, 53.0
        };
        const int inca    = 1;
        const int lda     = 2;
        double    c[]     = {
            59.0, 61.0,
            67.0, 71.0
        };
        const int incc    = 1;
        const int ldc     = 2;
        const int substep = 2;
        suzerain_lsrk_ex_substep(
                suzerain_lsrk_smr91,
                n, kl, ku,
                M, ldM, ipivM,
                delta_t, nrhs,
                a, inca, lda,
                b, incb, ldb,
                c, incc, ldc,
                substep);

        double    a_expected[]     = {
            -86.0/3.0, 1063.0/15.0,
              8.0/3.0, 2179.0/30.0
        };

        for (int i = 0; i < sizeof(a_expected)/sizeof(a_expected[0]); ++i) {
            gsl_test_abs(a[i], a_expected[i], 1.0e3*GSL_DBL_EPSILON,
                    "%s substep %d result a[%d]", __func__, substep, i);
        }
    }
}

double
riccati_equation_solution(double a, double b, double c, double t)
{
    return a + 1.0/(-1.0/(2*a+b) + c*exp(-(2*a+b)*t));
}

double
riccati_equation_nonlinear_operator(double y, double a, double b, double c)
{
    return pow(y, 2.0) - a*a - a*b;
}

void
check_smr91_imex_convergence_rate_riccati_equation()
{
    /* Solves (d/dt) y = y^2 + b y - a^2 -a b and compares against the
     * solution y(t) = a + (-(2*a+b)^(-1) + c*exp(-(2*a+b)*t))^(-1). */

    /* Time parameters; small duration puts us in asymptotic regime */
    const double t_initial = 0.140;
    const double t_final   = 0.145;

    /* Problem coefficients */
    const double coeff_a =   2.0;
    const double coeff_b =   2.0;
    const double coeff_c = -50.0;

    /* Problem exact solution */
    const double exact_initial
        = riccati_equation_solution(coeff_a, coeff_b, coeff_c, t_initial);
    const double exact_final
        = riccati_equation_solution(coeff_a, coeff_b, coeff_c, t_final);

    /* Chosen grid sizes */
    const int coarse_nsteps = 16;
    const int finer_nsteps  = 2*coarse_nsteps;

    /* Linear operator details so L = b = M^-1 xi[1] D */
    const int           n       = 1;
    const int           kl      = 0;
    const int           ku      = 0;
    const double        M[1]    = { 2.0 }; /* Compare xi */
    const int           ldM     = 1;
    const int           nD      = 1;
    const double        xi[1]   = { 2.0 }; /* Compare M */
    const double        D0[1]   = { coeff_b };
    const double *const D[1]    = { D0 };
    const int           ldD     = 1;
    const int           nrhs    = 1;

    /* Working array and storage parameters */
    double a[1] = { 0.0 }, b[1] = { 0.0 }, c[1] = { 0.0 };
    const int    inca = 1, incb = 1, incc = 1;
    const int    lda  = 1, ldb  = 1, ldc  = 1;

    /* Coarser grid calculation */
    double coarse_final;
    {
        const double delta_t = (t_final - t_initial)/coarse_nsteps;
        a[0] = exact_initial;

        for (int i = 0; i < coarse_nsteps; ++i) {
            for (int substep = 0;
                 substep < suzerain_lsrk_smr91->substeps;
                 ++substep) {
                b[0] = riccati_equation_nonlinear_operator(
                        a[0], coeff_a, coeff_b, coeff_c);
                gsl_test(suzerain_lsrk_imex_substep(
                            suzerain_lsrk_smr91,
                            n, kl, ku,
                            M, ldM,
                            nD, xi, D, ldD,
                            delta_t, nrhs,
                            a, inca, lda,
                            b, incb, ldb,
                            c, incc, ldc,
                            substep),
                        "Unexpected error reported in %s", __func__);
                c[0] = b[0];
            }
        }

        /* Next tolerance found using Octave code */
        gsl_test_abs(a[0], exact_final, 1.0e-10,
                "%s solution at t=%f using %d steps",
                __func__, t_final, coarse_nsteps);
    }
    coarse_final = a[0];
    const double coarse_error = fabs(exact_final - coarse_final);

    /* Finer grid calculation */
    double finer_final;
    {
        const double delta_t = (t_final - t_initial)/finer_nsteps;
        a[0] = exact_initial;

        for (int i = 0; i < finer_nsteps; ++i) {
            for (int substep = 0;
                 substep < suzerain_lsrk_smr91->substeps;
                 ++substep) {
                b[0] = riccati_equation_nonlinear_operator(
                        a[0], coeff_a, coeff_b, coeff_c);
                gsl_test(suzerain_lsrk_imex_substep(
                            suzerain_lsrk_smr91,
                            n, kl, ku,
                            M, ldM,
                            nD, xi, D, ldD,
                            delta_t, nrhs,
                            a, inca, lda,
                            b, incb, ldb,
                            c, incc, ldc,
                            substep),
                        "Unexpected error reported in %s", __func__);
                c[0] = b[0];
            }
        }

        /* Next tolerance found using Octave code */
        gsl_test_abs(a[0], exact_final, 1.0e-11,
                "%s solution at t=%f using %d steps",
                __func__, t_final, finer_nsteps);
    }
    finer_final = a[0];
    const double finer_error = fabs(exact_final - finer_final);

    /* SMR91 is second order against the Riccati problem */
    const double expected_order = 1.98; /* allows for floating point losses */
    const double observed_order
        = log(coarse_error/finer_error)/log(finer_nsteps/coarse_nsteps);
    gsl_test(gsl_isnan(observed_order),
            "%s observed convergence not NAN: %g", __func__, observed_order);
    gsl_test(observed_order < expected_order,
            "%s observed convergence order of %f not lower than %f",
            __func__, observed_order, expected_order);

    /* Richardson extrapolation should show h^2 term elimination
     * gives a better result than h^1, h^3, etc... */
    {
        gsl_matrix * data = gsl_matrix_alloc(1,2);
        gsl_vector * k = gsl_vector_alloc(1);
        gsl_matrix * normtable = gsl_matrix_alloc(data->size2, data->size2);
        gsl_vector * exact = gsl_vector_alloc(1);
        gsl_vector_set(exact, 0, exact_final);

        double richardson_h_error[4];
        for (int i = 0;
             i < sizeof(richardson_h_error)/sizeof(richardson_h_error[0]);
             ++i) {
            gsl_matrix_set(data, 0, 0, coarse_final);
            gsl_matrix_set(data, 0, 1, finer_final);
            gsl_vector_set(k,0,i+1);
            gsl_test(suzerain_richardson_extrapolation(
                        data, finer_nsteps/coarse_nsteps, k, normtable, exact),
                    "Unexpected error reported in %s");
            richardson_h_error[i] = gsl_matrix_get(
                    normtable, normtable->size2-1, normtable->size2-1);
        }

        gsl_matrix_free(normtable);
        gsl_vector_free(exact);
        gsl_vector_free(k);
        gsl_matrix_free(data);

        gsl_test(richardson_h_error[1] > richardson_h_error[0],
                "%s extrapolation should show h^2 leading error but %g > %g",
                __func__, richardson_h_error[1], richardson_h_error[0]);
        gsl_test(richardson_h_error[1] > richardson_h_error[2],
                "%s extrapolation should show h^2 leading error but %g > %g",
                __func__, richardson_h_error[1], richardson_h_error[2]);
        gsl_test(richardson_h_error[1] > richardson_h_error[3],
                "%s extrapolation should show h^2 leading error but %g > %g",
                __func__, richardson_h_error[1], richardson_h_error[3]);
    }
}

void
check_smr91_imex_convergence_rate_only_explicit_operator()
{
    /* Solves (d/dt) y = a*y and compares against y0*exp(a*t)
     * using only an explicit operator */

    /* Time parameters; small duration puts us in asymptotic regime */
    const double t_initial = 0.140;
    const double t_final   = 0.145;

    /* Problem coefficients */
    const double coeff_a   = 2.0;
    const double coeff_y0  = 1.0;

    /* Problem exact solution */
    const double exact_initial = coeff_y0*exp(coeff_a*t_initial);
    const double exact_final   = coeff_y0*exp(coeff_a*t_final);

    /* Chosen grid sizes */
    const int coarse_nsteps = 16;
    const int finer_nsteps  = 2*coarse_nsteps;

    /* Linear operator details so L = 0 = M^-1 xi[1] D */
    const int           n       = 1;
    const int           kl      = 0;
    const int           ku      = 0;
    const double        M[1]    = { 1.0 };
    const int           ldM     = 1;
    const int           nD      = 1;
    const double        xi[1]   = { 0.0 }; /* NOP L operator */
    const double        D0[1]   = { 1.0 };
    const double *const D[1]    = { D0 };
    const int           ldD     = 1;
    const int           nrhs    = 1;

    /* Working array and storage parameters */
    double a[1] = { 0.0 }, b[1] = { 0.0 }, c[1] = { 0.0 };
    const int    inca = 1, incb = 1, incc = 1;
    const int    lda  = 1, ldb  = 1, ldc  = 1;

    /* Coarser grid calculation */
    double coarse_final;
    {
        const double delta_t = (t_final - t_initial)/coarse_nsteps;
        a[0] = exact_initial;

        for (int i = 0; i < coarse_nsteps; ++i) {
            for (int substep = 0;
                 substep < suzerain_lsrk_smr91->substeps;
                 ++substep) {
                b[0] = coeff_a * a[0]; /* Explicit operator */
                gsl_test(suzerain_lsrk_imex_substep(
                            suzerain_lsrk_smr91,
                            n, kl, ku,
                            M, ldM,
                            nD, xi, D, ldD,
                            delta_t, nrhs,
                            a, inca, lda,
                            b, incb, ldb,
                            c, incc, ldc,
                            substep),
                        "Unexpected error reported in %s", __func__);
                c[0] = b[0];
            }
        }

        /* Next tolerance found using Octave code */
        gsl_test_abs(a[0], exact_final, 1.0e-12,
                "%s solution at t=%f using %d steps",
                __func__, t_final, coarse_nsteps);
    }
    coarse_final = a[0];
    const double coarse_error = fabs(exact_final - coarse_final);

    /* Finer grid calculation */
    double finer_final;
    {
        const double delta_t = (t_final - t_initial)/finer_nsteps;
        a[0] = exact_initial;

        for (int i = 0; i < finer_nsteps; ++i) {
            for (int substep = 0;
                 substep < suzerain_lsrk_smr91->substeps;
                 ++substep) {
                b[0] = coeff_a * a[0]; /* Explicit operator */
                gsl_test(suzerain_lsrk_imex_substep(
                            suzerain_lsrk_smr91,
                            n, kl, ku,
                            M, ldM,
                            nD, xi, D, ldD,
                            delta_t, nrhs,
                            a, inca, lda,
                            b, incb, ldb,
                            c, incc, ldc,
                            substep),
                        "Unexpected error reported in %s", __func__);
                c[0] = b[0];
            }
        }

        /* Next tolerance found using Octave code */
        gsl_test_abs(a[0], exact_final, 1.0e-13,
                "%s solution at t=%f using %d steps",
                __func__, t_final, finer_nsteps);
    }
    finer_final = a[0];
    const double finer_error = fabs(exact_final - finer_final);

    /* SMR91 is third order against the explicit-only problem */
    const double expected_order = 2.98; /* allows for floating point losses */
    const double observed_order
        = log(coarse_error/finer_error)/log(finer_nsteps/coarse_nsteps);
    gsl_test(gsl_isnan(observed_order),
            "%s observed convergence not NAN: %g", __func__, observed_order);
    gsl_test(observed_order < expected_order,
            "%s observed convergence order of %f not lower than %f",
            __func__, observed_order, expected_order);

    /* Richardson extrapolation should show h^3 term elimination
     * gives a better result than h^1, h^2, h^4, etc... */
    {
        gsl_matrix * data = gsl_matrix_alloc(1,2);
        gsl_vector * k = gsl_vector_alloc(1);
        gsl_matrix * normtable = gsl_matrix_alloc(data->size2, data->size2);
        gsl_vector * exact = gsl_vector_alloc(1);
        gsl_vector_set(exact, 0, exact_final);

        double richardson_h_error[4];
        for (int i = 0;
             i < sizeof(richardson_h_error)/sizeof(richardson_h_error[0]);
             ++i) {
            gsl_matrix_set(data, 0, 0, coarse_final);
            gsl_matrix_set(data, 0, 1, finer_final);
            gsl_vector_set(k,0,i+1);
            gsl_test(suzerain_richardson_extrapolation(
                        data, finer_nsteps/coarse_nsteps, k, normtable, exact),
                    "Unexpected error reported in %s");
            richardson_h_error[i] = gsl_matrix_get(
                    normtable, normtable->size2-1, normtable->size2-1);
        }

        gsl_matrix_free(normtable);
        gsl_vector_free(exact);
        gsl_vector_free(k);
        gsl_matrix_free(data);

        gsl_test(richardson_h_error[2] > richardson_h_error[0],
                "%s extrapolation should show h^2 leading error but %g > %g",
                __func__, richardson_h_error[2], richardson_h_error[0]);
        gsl_test(richardson_h_error[2] > richardson_h_error[1],
                "%s extrapolation should show h^2 leading error but %g > %g",
                __func__, richardson_h_error[2], richardson_h_error[1]);
        gsl_test(richardson_h_error[2] > richardson_h_error[3],
                "%s extrapolation should show h^2 leading error but %g > %g",
                __func__, richardson_h_error[2], richardson_h_error[3]);
    }
}

void
check_smr91_ex_convergence_rate_only_explicit_operator()
{
    /* Solves (d/dt) y = a*y and compares against y0*exp(a*t)
     * using only an explicit operator */

    /* Time parameters; small duration puts us in asymptotic regime */
    const double t_initial = 0.140;
    const double t_final   = 0.145;

    /* Problem coefficients */
    const double coeff_a   = 2.0;
    const double coeff_y0  = 1.0;

    /* Problem exact solution */
    const double exact_initial = coeff_y0*exp(coeff_a*t_initial);
    const double exact_final   = coeff_y0*exp(coeff_a*t_final);

    /* Chosen grid sizes */
    const int coarse_nsteps = 16;
    const int finer_nsteps  = 2*coarse_nsteps;

    /* Linear operator details so L = 0 = M^-1 xi[1] D */
    const int           n       = 1;
    const int           kl      = 0;
    const int           ku      = 0;
    double              M[1]    = { 1.0 };
    const int           ldM     = 1;
    int                 ipivM[1];
    const int           nrhs    = 1;

    /* Working array and storage parameters */
    double a[1] = { 0.0 }, b[1] = { 0.0 }, c[1] = { 0.0 };
    const int    inca = 1, incb = 1, incc = 1;
    const int    lda  = 1, ldb  = 1, ldc  = 1;

    /* Coarser grid calculation */
    double coarse_final;
    {
        const double delta_t = (t_final - t_initial)/coarse_nsteps;
        a[0] = exact_initial;

        for (int i = 0; i < coarse_nsteps; ++i) {
            for (int substep = 0;
                 substep < suzerain_lsrk_smr91->substeps;
                 ++substep) {
                b[0] = coeff_a * a[0]; /* Explicit operator */
                gsl_test(suzerain_lsrk_ex_substep(
                            suzerain_lsrk_smr91,
                            n, kl, ku,
                            M, ldM, ipivM,
                            delta_t, nrhs,
                            a, inca, lda,
                            b, incb, ldb,
                            c, incc, ldc,
                            substep),
                        "Unexpected error reported in %s", __func__);
                c[0] = b[0];
            }
        }

        /* Next tolerance found using Octave code */
        gsl_test_abs(a[0], exact_final, 1.0e-12,
                "%s solution at t=%f using %d steps",
                __func__, t_final, coarse_nsteps);
    }
    coarse_final = a[0];
    const double coarse_error = fabs(exact_final - coarse_final);

    /* Finer grid calculation */
    double finer_final;
    {
        const double delta_t = (t_final - t_initial)/finer_nsteps;
        a[0] = exact_initial;

        for (int i = 0; i < finer_nsteps; ++i) {
            for (int substep = 0;
                 substep < suzerain_lsrk_smr91->substeps;
                 ++substep) {
                b[0] = coeff_a * a[0]; /* Explicit operator */
                gsl_test(suzerain_lsrk_ex_substep(
                            suzerain_lsrk_smr91,
                            n, kl, ku,
                            M, ldM, ipivM,
                            delta_t, nrhs,
                            a, inca, lda,
                            b, incb, ldb,
                            c, incc, ldc,
                            substep),
                        "Unexpected error reported in %s", __func__);
                c[0] = b[0];
            }
        }

        /* Next tolerance found using Octave code */
        gsl_test_abs(a[0], exact_final, 1.0e-13,
                "%s solution at t=%f using %d steps",
                __func__, t_final, finer_nsteps);
    }
    finer_final = a[0];
    const double finer_error = fabs(exact_final - finer_final);

    /* SMR91 is third order against the explicit-only problem */
    const double expected_order = 2.98; /* allows for floating point losses */
    const double observed_order
        = log(coarse_error/finer_error)/log(finer_nsteps/coarse_nsteps);
    gsl_test(gsl_isnan(observed_order),
            "%s observed convergence not NAN: %g", __func__, observed_order);
    gsl_test(observed_order < expected_order,
            "%s observed convergence order of %f not lower than %f",
            __func__, observed_order, expected_order);

    /* Richardson extrapolation should show h^3 term elimination
     * gives a better result than h^1, h^2, h^4, etc... */
    {
        gsl_matrix * data = gsl_matrix_alloc(1,2);
        gsl_vector * k = gsl_vector_alloc(1);
        gsl_matrix * normtable = gsl_matrix_alloc(data->size2, data->size2);
        gsl_vector * exact = gsl_vector_alloc(1);
        gsl_vector_set(exact, 0, exact_final);

        double richardson_h_error[4];
        for (int i = 0;
             i < sizeof(richardson_h_error)/sizeof(richardson_h_error[0]);
             ++i) {
            gsl_matrix_set(data, 0, 0, coarse_final);
            gsl_matrix_set(data, 0, 1, finer_final);
            gsl_vector_set(k,0,i+1);
            gsl_test(suzerain_richardson_extrapolation(
                        data, finer_nsteps/coarse_nsteps, k, normtable, exact),
                    "Unexpected error reported in %s");
            richardson_h_error[i] = gsl_matrix_get(
                    normtable, normtable->size2-1, normtable->size2-1);
        }

        gsl_matrix_free(normtable);
        gsl_vector_free(exact);
        gsl_vector_free(k);
        gsl_matrix_free(data);

        gsl_test(richardson_h_error[2] > richardson_h_error[0],
                "%s extrapolation should show h^2 leading error but %g > %g",
                __func__, richardson_h_error[2], richardson_h_error[0]);
        gsl_test(richardson_h_error[2] > richardson_h_error[1],
                "%s extrapolation should show h^2 leading error but %g > %g",
                __func__, richardson_h_error[2], richardson_h_error[1]);
        gsl_test(richardson_h_error[2] > richardson_h_error[3],
                "%s extrapolation should show h^2 leading error but %g > %g",
                __func__, richardson_h_error[2], richardson_h_error[3]);
    }
}

int
main(int argc, char **argv)
{
    gsl_ieee_env_setup();

    check_smr91_constants();

    // Using SMR91 in hybrid implicit/explicit mode
    check_smr91_imex_scalareqn_substeps();
    check_smr91_imex_matrixeqn_substeps();
    check_smr91_imex_convergence_rate_riccati_equation();
    check_smr91_imex_convergence_rate_only_explicit_operator();

    // Using SMR91 in explicit-only mode
    check_smr91_ex_scalareqn_substeps();
    check_smr91_ex_matrixeqn_substeps();
    check_smr91_ex_matrixeqn_substeps_contiguous();
    check_smr91_ex_convergence_rate_only_explicit_operator();

    exit(gsl_test_summary());
}
