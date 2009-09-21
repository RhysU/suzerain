/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2009 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
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
 *
 * timestepper.c: SMR91 timestepping implementation
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <config.h>

#include <stdlib.h>
#include <stdio.h>
#include <suzerain/blas_et_al.h>
#include <suzerain/error.h>
#include <suzerain/timestepper.h>

static const double _smr91_alpha[3] = { 29.0/ 96.0,  -3.0/40.0,  1.0/ 6.0 };
static const double _smr91_beta[3]  = { 37.0/160.0,   5.0/24.0,  1.0/ 6.0 };
static const double _smr91_gamma[3] = { 8.0/  15.0,   5.0/12.0,  3.0/ 4.0 };
static const double _smr91_zeta[3]  = {        0.0, -17.0/60.0, -5.0/12.0 };
static const suzerain_lsrk_method _suzerain_lsrk_smr91 = {
    "SMR91", 3, _smr91_alpha, _smr91_beta, _smr91_gamma, _smr91_zeta
};
const suzerain_lsrk_method * const suzerain_lsrk_smr91 = &_suzerain_lsrk_smr91;

int
suzerain_lsrk_substep(
        const suzerain_lsrk_method * const method,
        const int n,
        const int kl,
        const int ku,
        const double * const M,
        const int ldM,
        const int nD,
        const double * const xi,
        const double * const * const D,
        const int ldD,
        double delta_t,
        const int nrhs,
        double       * const a, const int inca, const int lda,
        const double * const b, const int incb, const int ldb,
        double       * const c, const int incc, const int ldc,
        const int substep)
{
    if (substep < 0 || substep >= method->substeps) {
        SUZERAIN_ERROR("requested substep out of range", SUZERAIN_EINVAL);
    }

    /* Allocate and clear working space for matrix \hat{D} */
    const int ld_hatD = kl + 1 + ku;
    double * const hatD = suzerain_blas_calloc(n*ld_hatD, sizeof(hatD[0]));
    if (hatD == NULL) {
        SUZERAIN_ERROR("failed to allocate space for hatD", SUZERAIN_ENOMEM);
    }
    /* Accumulate $\hat{D} = \sum_j (\gamma_i + \zeta_{i-1}) \xi_j D_j$ */
    for (int j = 0; j < nD; ++j) {
        const double alpha
            = (method->gamma[substep] + method->zeta[substep]) * xi[j];
        suzerain_blas_dgb_acc(
                n, n, kl, ku, alpha, D[j], ldD, 1.0, hatD, ld_hatD);
    }

    /* Allocate and clear space for matrix \hat{M} with LU-ready padding */
    const int ld_hatM = 2*kl + 1 + ku;
    double * const hatM = suzerain_blas_calloc(n * ld_hatM, sizeof(hatM[0]));
    if (hatM == NULL) {
        free(hatD);
        SUZERAIN_ERROR("failed to allocate space for hatM", SUZERAIN_ENOMEM);
    }
    /* Accumulate $\hat{M}  = M - \sum_j \Delta{}t \beta_i \xi_j D_j */
    suzerain_blas_dgb_acc(n, n, kl, ku, 1.0, M, ldM, 0.0, hatM + kl, ld_hatM);
    for (int j = 0; j < nD; ++j) {
        const double alpha = -1.0 * delta_t * method->beta[substep] * xi[j];
        suzerain_blas_dgb_acc(
            n, n, kl, ku, alpha, D[j], ldD, 1.0, hatM + kl, ld_hatM);
    }

    /* Allocate space for the LU-factorization pivot information */
    int * const ipiv = suzerain_blas_malloc(n * sizeof(ipiv[0]));
    if (ipiv == NULL) {
        free(hatM);
        free(hatD);
        SUZERAIN_ERROR("failed to allocate space for ipiv", SUZERAIN_ENOMEM);
    }

    /* Find LU factorization of \hat{M} */
    const int error_dgbtrf = suzerain_lapack_dgbtrf(
            n, n, kl, ku, hatM, ld_hatM, ipiv);
    if (error_dgbtrf) {
        free(ipiv);
        free(hatM);
        free(hatD);
        char buffer[72];
        snprintf(buffer, sizeof(buffer)/sizeof(buffer[0]),
                 "suzerain_lapack_dgbtrf reported error: %d", error_dgbtrf);
        SUZERAIN_ERROR(buffer, SUZERAIN_ESANITY);
    }

    /* Allocate working space for contiguous vector d */
    const int incd = 1;
    double * const d = suzerain_blas_malloc(n * incd * sizeof(d[0]));
    if (d == NULL) {
        free(ipiv);
        free(hatM);
        free(hatD);
        SUZERAIN_ERROR("failed to allocate space for d", SUZERAIN_ENOMEM);
    }

    for (int j = 0; j < nrhs; ++j) {
        double       * const  a_j = a + j*lda;
        const double * const  b_j = b + j*ldb;
        double       * const  c_j = c + j*ldc;

        suzerain_blas_daxpby(n,
                method->gamma[substep], b_j, incb,
                method->zeta[substep],  c_j, incc);

        suzerain_blas_dgbmv('N', n, n, kl, ku,
                1.0, M, ldM, c_j, incc,
                0.0, d, incd);

        suzerain_blas_dgbmv('N', n, n, kl, ku,
                1.0, hatD, ld_hatD, a_j, inca,
                1.0, d, incd);

        const int error_dgbtrs = suzerain_lapack_dgbtrs(
                'N', n, kl, ku, 1, hatM, ld_hatM, ipiv, d, n*incd);
        if (error_dgbtrs) {
            free(d);
            free(ipiv);
            free(hatM);
            free(hatD);
            char buffer[72];
            snprintf(buffer, sizeof(buffer)/sizeof(buffer[0]),
                    "suzerain_lapack_dgbtrs reported error: %d", error_dgbtrs);
            SUZERAIN_ERROR(buffer, SUZERAIN_ESANITY);
        }

        suzerain_blas_daxpy(n,
                delta_t, d, incd,
                a_j, inca);
    }

    free(d);
    free(ipiv);
    free(hatM);
    free(hatD);

    return SUZERAIN_SUCCESS;
}
