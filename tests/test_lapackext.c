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

#include <suzerain/blas_et_al/lapackext.h>

#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_test.h>

#include <suzerain/blas_et_al/blas.h>

#ifdef SUZERAIN_HAVE_MKL
#include <mkl.h>
#include <mkl_service.h>
#endif

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
