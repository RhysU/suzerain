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

// File iteration is used to generate bandwidth-specific routines.  See
// http://www.boost.org/libs/preprocessor/doc/topics/file_iteration.html for an
// overview.  We iterate first to generate internal, static routines and then
// continue on to the non-iterated logic.

#if !defined(BOOST_PP_IS_ITERATING) || !BOOST_PP_IS_ITERATING

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.h>
#pragma hdrstop
#include <suzerain/gbmv.h>

#pragma warning(disable:1418 1572 2259)

// ---------------------------------------------------------
// Iterate on this file to generate fixed bandwidth routines
// ---------------------------------------------------------
#define BOOST_PP_ITERATION_LIMITS (0, 10)
#define BOOST_PP_FILENAME_1 <suzerain/gbmv.c>
#include BOOST_PP_ITERATE()

// -----------------------------------
// Generate general bandwidth routines
// -----------------------------------

#define GBMV_STATIC    static
#define GBMV_FUNCTION  suzerain_gbmv_internal_s
#define GBMV_COMPONENT float
#define GBMV_SCALAR    float
#define GBMV_KL        const int kl,
#define GBMV_KU        const int ku,
#define GBMV_LDA       const int lda,
#include "gbmv.def"

#define GBMV_STATIC    static
#define GBMV_FUNCTION  suzerain_gbmv_internal_d
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

#define GBMV_STATIC    static
#define GBMV_FUNCTION  suzerain_gbmv_internal_dz
#define GBMV_COMPONENT double
#define GBMV_SCALAR    double _Complex
#define GBMV_KL        const int kl,
#define GBMV_KU        const int ku,
#define GBMV_LDA       const int lda,
#include "gbmv.def"

// ------------------------------------------------------------------
// Provide externally callable logic dispatching to internal routines
// ------------------------------------------------------------------

int
suzerain_gbmv_s(
        const char trans,
        const int m,
        const int n,
        const int kl,
        const int ku,
        const float alpha,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    // TODO Dispatch to fixed bandwidth versions if appropriate
    return suzerain_gbmv_internal_s(trans, m, n, kl, ku,
                                    alpha, a, lda, x, incx,
                                    beta,          y, incy);
}

int
suzerain_gbmv_d(
        const char trans,
        const int m,
        const int n,
        const int kl,
        const int ku,
        const double alpha,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    // TODO Dispatch to fixed bandwidth versions if appropriate
    return suzerain_gbmv_internal_d(trans, m, n, kl, ku,
                                    alpha, a, lda, x, incx,
                                    beta,          y, incy);
}

int
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
    // TODO Dispatch to fixed bandwidth versions if appropriate
    return suzerain_gbmv_internal_sc(trans, m, n, kl, ku,
                                     alpha_c, (void *) a, lda,
                                              (void *) x, incx,
                                     beta_c,  (void *) y, incy);
}

int
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
    // TODO Dispatch to fixed bandwidth versions if appropriate
    return suzerain_gbmv_internal_dz(trans, m, n, kl, ku,
                                     alpha_c, (void *) a, lda,
                                              (void *) x, incx,
                                     beta_c,  (void *) y, incy);
}

#else

// ----------------------------------------------------------
// Generate fixed bandwidth routines using BOOST_PP_ITERATE()
// ----------------------------------------------------------

// TODO Add logic here

#endif
