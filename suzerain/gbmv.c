/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2011, 2012 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *--------------------------------------------------------------------------
 * gbmv.c: provides BLAS-like general band matrix-vector operations
 * $Id$
 */

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
#define FIXEDBW_LOWER 0
#define FIXEDBW_UPPER 15
#define BOOST_PP_ITERATION_LIMITS (FIXEDBW_LOWER, FIXEDBW_UPPER)
#define BOOST_PP_FILENAME_1 <suzerain/gbmv.c>
#include BOOST_PP_ITERATE()

// -----------------------------------
// Generate general bandwidth routines
// -----------------------------------

#define STATIC   static
#define FUNCTION suzerain_gbmv_internal_s
#define TYPE_A   float
#define TYPE_X   float
#define TYPE_Y   float
#define KL       const int kl,
#define KU       const int ku,
#include "gbmv.def"

#define STATIC   static
#define FUNCTION suzerain_gbmv_internal_ssc
#define TYPE_A   float
#define TYPE_X   float
#define TYPE_Y   complex_float
#define KL       const int kl,
#define KU       const int ku,
#include "gbmv.def"

#define STATIC   static
#define FUNCTION suzerain_gbmv_internal_scc
#define TYPE_A   float
#define TYPE_X   complex_float
#define TYPE_Y   complex_float
#define KL       const int kl,
#define KU       const int ku,
#include "gbmv.def"

#define STATIC   static
#define FUNCTION suzerain_gbmv_internal_d
#define TYPE_A   double
#define TYPE_X   double
#define TYPE_Y   double
#define KL       const int kl,
#define KU       const int ku,
#include "gbmv.def"

#define STATIC   static
#define FUNCTION suzerain_gbmv_internal_ddz
#define TYPE_A   double
#define TYPE_X   double
#define TYPE_Y   complex_double
#define KL       const int kl,
#define KU       const int ku,
#include "gbmv.def"

#define STATIC   static
#define FUNCTION suzerain_gbmv_internal_dzz
#define TYPE_A   double
#define TYPE_X   complex_double
#define TYPE_Y   complex_double
#define KL       const int kl,
#define KU       const int ku,
#include "gbmv.def"

// ------------------------------------------------------------------
// Provide externally callable logic dispatching to internal routines
// ------------------------------------------------------------------

#define FIXEDBW_CASE(z,num,prefix)                           \
    case num: return BOOST_PP_CAT(prefix,num)(               \
        trans, m, n, alpha, a, lda, x, incx, beta, y, incy);

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
    // Dispatch to fixed bandwidth specialization for small bandwidth...
    if (kl == ku) {
        switch (kl) {
            BOOST_PP_REPEAT_FROM_TO(FIXEDBW_LOWER, BOOST_PP_INC(FIXEDBW_UPPER),
                                    FIXEDBW_CASE, suzerain_gbmv_internal_s)
        }
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_gbmv_internal_s(trans, m, n, kl, ku,
                                    alpha, a, lda, x, incx,
                                    beta,          y, incy);
}

int
suzerain_gbmv_ssc(
        const char trans,
        const int m,
        const int n,
        const int kl,
        const int ku,
        const complex_float alpha,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy)
{
    // Dispatch to fixed bandwidth specialization for small bandwidth...
    if (kl == ku) {
        switch (kl) {
            BOOST_PP_REPEAT_FROM_TO(FIXEDBW_LOWER, BOOST_PP_INC(FIXEDBW_UPPER),
                                    FIXEDBW_CASE, suzerain_gbmv_internal_ssc)
        }
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_gbmv_internal_ssc(trans, m, n, kl, ku,
                                      alpha, a, lda, x, incx,
                                      beta,          y, incy);
}

int
suzerain_gbmv_scc(
        const char trans,
        const int m,
        const int n,
        const int kl,
        const int ku,
        const complex_float alpha,
        const float *a,
        const int lda,
        const complex_float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy)
{
    // Dispatch to fixed bandwidth specialization for small bandwidth...
    if (kl == ku) {
        switch (kl) {
            BOOST_PP_REPEAT_FROM_TO(FIXEDBW_LOWER, BOOST_PP_INC(FIXEDBW_UPPER),
                                    FIXEDBW_CASE, suzerain_gbmv_internal_scc)
        }
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_gbmv_internal_scc(trans, m, n, kl, ku,
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
    // Dispatch to fixed bandwidth specialization for small bandwidth...
    if (kl == ku) {
        switch (kl) {
            BOOST_PP_REPEAT_FROM_TO(FIXEDBW_LOWER, BOOST_PP_INC(FIXEDBW_UPPER),
                                    FIXEDBW_CASE, suzerain_gbmv_internal_d)
        }
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_gbmv_internal_d(trans, m, n, kl, ku,
                                    alpha, a, lda, x, incx,
                                    beta,          y, incy);
}

int
suzerain_gbmv_ddz(
        const char trans,
        const int m,
        const int n,
        const int kl,
        const int ku,
        const complex_double alpha,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy)
{
    // Dispatch to fixed bandwidth specialization for small bandwidth...
    if (kl == ku) {
        switch (kl) {
            BOOST_PP_REPEAT_FROM_TO(FIXEDBW_LOWER, BOOST_PP_INC(FIXEDBW_UPPER),
                                    FIXEDBW_CASE, suzerain_gbmv_internal_ddz)
        }
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_gbmv_internal_ddz(trans, m, n, kl, ku,
                                      alpha, a, lda, x, incx,
                                      beta,          y, incy);
}

int
suzerain_gbmv_dzz(
        const char trans,
        const int m,
        const int n,
        const int kl,
        const int ku,
        const complex_double alpha,
        const double *a,
        const int lda,
        const complex_double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy)
{
    // Dispatch to fixed bandwidth specialization for small bandwidth...
    if (kl == ku) {
        switch (kl) {
            BOOST_PP_REPEAT_FROM_TO(FIXEDBW_LOWER, BOOST_PP_INC(FIXEDBW_UPPER),
                                    FIXEDBW_CASE, suzerain_gbmv_internal_dzz)
        }
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_gbmv_internal_dzz(trans, m, n, kl, ku,
                                      alpha, a, lda, x, incx,
                                      beta,          y, incy);
}

#else

// ------------------------------------------------------------
// Generate fixed bandwidth routines using BOOST_PP_ITERATION()
// ------------------------------------------------------------
#define kl BOOST_PP_ITERATION()
#define ku BOOST_PP_ITERATION()

#define STATIC   static
#define FUNCTION BOOST_PP_CAT(suzerain_gbmv_internal_s, BOOST_PP_ITERATION())
#define TYPE_A   float
#define TYPE_X   float
#define TYPE_Y   float
#define KL
#define KU
#include "gbmv.def"

#define STATIC   static
#define FUNCTION BOOST_PP_CAT(suzerain_gbmv_internal_ssc, BOOST_PP_ITERATION())
#define TYPE_A   float
#define TYPE_X   float
#define TYPE_Y   complex_float
#define KL
#define KU
#include "gbmv.def"

#define STATIC   static
#define FUNCTION BOOST_PP_CAT(suzerain_gbmv_internal_scc, BOOST_PP_ITERATION())
#define TYPE_A   float
#define TYPE_X   complex_float
#define TYPE_Y   complex_float
#define KL
#define KU
#include "gbmv.def"

#define STATIC   static
#define FUNCTION BOOST_PP_CAT(suzerain_gbmv_internal_d, BOOST_PP_ITERATION())
#define TYPE_A   double
#define TYPE_X   double
#define TYPE_Y   double
#define KL
#define KU
#include "gbmv.def"

#define STATIC   static
#define FUNCTION BOOST_PP_CAT(suzerain_gbmv_internal_ddz, BOOST_PP_ITERATION())
#define TYPE_A   double
#define TYPE_X   double
#define TYPE_Y   complex_double
#define KL
#define KU
#include "gbmv.def"

#define STATIC   static
#define FUNCTION BOOST_PP_CAT(suzerain_gbmv_internal_dzz, BOOST_PP_ITERATION())
#define TYPE_A   double
#define TYPE_X   complex_double
#define TYPE_Y   complex_double
#define KL
#define KU
#include "gbmv.def"

#undef kl
#undef ku

#endif
