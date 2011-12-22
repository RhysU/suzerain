/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2011 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * gbmv.c: provides custom, BLAS-like general band matrix-vector operations
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.h>
#pragma hdrstop

#pragma warning(disable:1418 1572 2259)

// Produce general bandwidth versions

#define GBMV_STATIC
#define GBMV_FUNCTION  suzerain_gbmv_s
#define GBMV_COMPONENT float
#define GBMV_SCALAR    float
#define GBMV_KL        const int kl,
#define GBMV_KU        const int ku,
#define GBMV_LDA       const int lda,
#include "gbmv.def"

#define GBMV_STATIC
#define GBMV_FUNCTION  suzerain_gbmv_d
#define GBMV_COMPONENT double
#define GBMV_SCALAR    double
#define GBMV_KL        const int kl,
#define GBMV_KU        const int ku,
#define GBMV_LDA       const int lda,
#include "gbmv.def"

#define GBMV_STATIC    static
#define GBMV_FUNCTION  suzerain_gbmv_internal_sc
#define GBMV_COMPONENT float
#define GBMV_SCALAR    float _Complex
#define GBMV_KL        const int kl,
#define GBMV_KU        const int ku,
#define GBMV_LDA       const int lda,
#include "gbmv.def"

void
suzerain_gbmv_sc(
        const char trans,
        const int m,
        const int n,
        const int kl,
        const int ku,
        const float alpha[2],
        const float *a,
        const int lda,
        const float (*x)[2],
        const int incx,
        const float beta[2],
        float (*y)[2],
        const int incy)
{
    float _Complex alpha_c, beta_c;
    memcpy(&alpha_c, alpha, sizeof(float _Complex));
    memcpy(&beta_c,  beta,  sizeof(float _Complex));
    suzerain_gbmv_internal_sc(trans, m, n, kl, ku,
                              alpha_c, (void *) a, lda,
                                       (void *) x, incx,
                              beta_c,  (void *) y, incy);
}

#define GBMV_STATIC    static
#define GBMV_FUNCTION  suzerain_gbmv_internal_dz
#define GBMV_COMPONENT double
#define GBMV_SCALAR    double _Complex
#define GBMV_KL        const int kl,
#define GBMV_KU        const int ku,
#define GBMV_LDA       const int lda,
#include "gbmv.def"

void
suzerain_gbmv_dz(
        const char trans,
        const int m,
        const int n,
        const int kl,
        const int ku,
        const double alpha[2],
        const double *a,
        const int lda,
        const double (*x)[2],
        const int incx,
        const double beta[2],
        double (*y)[2],
        const int incy)
{
    double _Complex alpha_c, beta_c;
    memcpy(&alpha_c, alpha, sizeof(double _Complex));
    memcpy(&beta_c,  beta,  sizeof(double _Complex));
    suzerain_gbmv_internal_dz(trans, m, n, kl, ku,
                              alpha_c, (void *) a, lda,
                                       (void *) x, incx,
                              beta_c,  (void *) y, incy);
}
