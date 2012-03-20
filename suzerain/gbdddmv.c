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
 * gbdddmv.c: provides BLAS-like general band matrix-vector operations
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
#include <suzerain/gbdddmv.h>

#pragma warning(disable:1418 1572 2259)

// ---------------------------------------------------------
// Iterate on this file to generate fixed bandwidth routines
// ---------------------------------------------------------
#define FIXEDBW_LOWER 0
#define FIXEDBW_UPPER 15
#define BOOST_PP_ITERATION_LIMITS (FIXEDBW_LOWER, FIXEDBW_UPPER)
#define BOOST_PP_FILENAME_1 <suzerain/gbdddmv.c>
#include BOOST_PP_ITERATE()

// -----------------------------------
// Generate general bandwidth routines
// -----------------------------------

#define GBDDDMV_STATIC    static
#define GBDDDMV_FUNCTION  suzerain_gbdddmv_internal_s
#define GBDDDMV_COMPONENT float
#define GBDDDMV_SCALAR    float
#define GBDDDMV_KL        const int kl,
#define GBDDDMV_KU        const int ku,
#include "gbdddmv.def"

#define GBDDDMV_STATIC    static
#define GBDDDMV_FUNCTION  suzerain_gbdddmv_internal_d
#define GBDDDMV_COMPONENT double
#define GBDDDMV_SCALAR    double
#define GBDDDMV_KL        const int kl,
#define GBDDDMV_KU        const int ku,
#include "gbdddmv.def"

#define GBDDDMV_STATIC    static
#define GBDDDMV_FUNCTION  suzerain_gbdddmv_internal_sc
#define GBDDDMV_COMPONENT float
#define GBDDDMV_SCALAR    complex_float
#define GBDDDMV_KL        const int kl,
#define GBDDDMV_KU        const int ku,
#include "gbdddmv.def"

#define GBDDDMV_STATIC    static
#define GBDDDMV_FUNCTION  suzerain_gbdddmv_internal_dz
#define GBDDDMV_COMPONENT double
#define GBDDDMV_SCALAR    complex_double
#define GBDDDMV_KL        const int kl,
#define GBDDDMV_KU        const int ku,
#include "gbdddmv.def"

// ------------------------------------------------------------------
// Provide externally callable logic dispatching to internal routines
// ------------------------------------------------------------------

#define FIXEDBW_CASE(z,num,prefix)                                      \
    case num: return BOOST_PP_CAT(prefix,num)(                          \
        trans, n, alpha0, d0, ldd0, alpha1, d1, ldd1, alpha2, d2, ldd2, \
        a, lda, x, incx, beta, y, incy);

int
suzerain_gbdddmv_s(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha0,
        const float *d0,
        const int ldd0,
        const float alpha1,
        const float *d1,
        const int ldd1,
        const float alpha2,
        const float *d2,
        const int ldd2,
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
                                    FIXEDBW_CASE, suzerain_gbdddmv_internal_s)
        }
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_gbdddmv_internal_s(trans, n, kl, ku,
                                       alpha0, d0, ldd0,
                                       alpha1, d1, ldd1,
                                       alpha2, d2, ldd2,
                                       a, lda, x, incx, beta, y, incy);
}

int
suzerain_gbdddmv_d(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha0,
        const double *d0,
        const int ldd0,
        const double alpha1,
        const double *d1,
        const int ldd1,
        const double alpha2,
        const double *d2,
        const int ldd2,
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
                                    FIXEDBW_CASE, suzerain_gbdddmv_internal_d)
        }
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_gbdddmv_internal_d(trans, n, kl, ku,
                                       alpha0, d0, ldd0,
                                       alpha1, d1, ldd1,
                                       alpha2, d2, ldd2,
                                       a, lda, x, incx, beta, y, incy);
}

int
suzerain_gbdddmv_sc(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const complex_float alpha0,
        const float *d0,
        const int ldd0,
        const complex_float alpha1,
        const float *d1,
        const int ldd1,
        const complex_float alpha2,
        const float *d2,
        const int ldd2,
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
                                    FIXEDBW_CASE, suzerain_gbdddmv_internal_sc)
        }
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_gbdddmv_internal_sc(trans, n, kl, ku,
                                        alpha0, d0, ldd0,
                                        alpha1, d1, ldd1,
                                        alpha2, d2, ldd2,
                                        a, lda, x, incx, beta, y, incy);
}

int
suzerain_gbdddmv_dz(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const complex_double alpha0,
        const double *d0,
        const int ldd0,
        const complex_double alpha1,
        const double *d1,
        const int ldd1,
        const complex_double alpha2,
        const double *d2,
        const int ldd2,
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
                                    FIXEDBW_CASE, suzerain_gbdddmv_internal_dz)
        }
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_gbdddmv_internal_dz(trans, n, kl, ku,
                                        alpha0, d0, ldd0,
                                        alpha1, d1, ldd1,
                                        alpha2, d2, ldd2,
                                        a, lda, x, incx, beta, y, incy);
}

#else

// ------------------------------------------------------------
// Generate fixed bandwidth routines using BOOST_PP_ITERATION()
// ------------------------------------------------------------
#define kl BOOST_PP_ITERATION()
#define ku BOOST_PP_ITERATION()

#define GBDDDMV_STATIC    static
#define GBDDDMV_FUNCTION  BOOST_PP_CAT(suzerain_gbdddmv_internal_s, \
                                       BOOST_PP_ITERATION())
#define GBDDDMV_COMPONENT float
#define GBDDDMV_SCALAR    float
#define GBDDDMV_KL
#define GBDDDMV_KU
#include "gbdddmv.def"

#define GBDDDMV_STATIC    static
#define GBDDDMV_FUNCTION  BOOST_PP_CAT(suzerain_gbdddmv_internal_d, \
                                       BOOST_PP_ITERATION())
#define GBDDDMV_COMPONENT double
#define GBDDDMV_SCALAR    double
#define GBDDDMV_KL
#define GBDDDMV_KU
#include "gbdddmv.def"

#define GBDDDMV_STATIC    static
#define GBDDDMV_FUNCTION  BOOST_PP_CAT(suzerain_gbdddmv_internal_sc, \
                                       BOOST_PP_ITERATION())
#define GBDDDMV_COMPONENT float
#define GBDDDMV_SCALAR    complex_float
#define GBDDDMV_KL
#define GBDDDMV_KU
#include "gbdddmv.def"

#define GBDDDMV_STATIC    static
#define GBDDDMV_FUNCTION  BOOST_PP_CAT(suzerain_gbdddmv_internal_dz, \
                                       BOOST_PP_ITERATION())
#define GBDDDMV_COMPONENT double
#define GBDDDMV_SCALAR    complex_double
#define GBDDDMV_KL
#define GBDDDMV_KU
#include "gbdddmv.def"

#undef kl
#undef ku

#endif
