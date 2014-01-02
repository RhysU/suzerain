/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2011-2014 Rhys Ulerich
 * Copyright (C) 2012-2014 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
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
 */

/** @file
 * @copydoc gbdddmv.h
 */

// File iteration is used to generate bandwidth-specific routines.  See
// http://www.boost.org/libs/preprocessor/doc/topics/file_iteration.html for an
// overview.  We iterate first to generate internal, static routines and then
// continue on to the non-iterated logic.

#if !defined(BOOST_PP_IS_ITERATING) || !BOOST_PP_IS_ITERATING

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/gbdddmv.h>

#include <suzerain/common.h>

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

#define STATIC   static
#define FUNCTION suzerain_gbdddmv_internal_s
#define TYPE_D   float
#define TYPE_A   float
#define TYPE_X   float
#define TYPE_Y   float
#define KL       const int kl,
#define KU       const int ku,
#include "gbdddmv.def"

#define STATIC   static
#define FUNCTION suzerain_gbdddmv_internal_ssc
#define TYPE_D   float
#define TYPE_A   float
#define TYPE_X   float
#define TYPE_Y   complex_float
#define KL       const int kl,
#define KU       const int ku,
#include "gbdddmv.def"

#define STATIC   static
#define FUNCTION suzerain_gbdddmv_internal_scc
#define TYPE_D   float
#define TYPE_A   float
#define TYPE_X   complex_float
#define TYPE_Y   complex_float
#define KL       const int kl,
#define KU       const int ku,
#include "gbdddmv.def"

#define STATIC   static
#define FUNCTION suzerain_gbdddmv_internal_d
#define TYPE_D   double
#define TYPE_A   double
#define TYPE_X   double
#define TYPE_Y   double
#define KL       const int kl,
#define KU       const int ku,
#include "gbdddmv.def"

#define STATIC   static
#define FUNCTION suzerain_gbdddmv_internal_ddz
#define TYPE_D   double
#define TYPE_A   double
#define TYPE_X   double
#define TYPE_Y   complex_double
#define KL       const int kl,
#define KU       const int ku,
#include "gbdddmv.def"

#define STATIC   static
#define FUNCTION suzerain_gbdddmv_internal_dzz
#define TYPE_D   double
#define TYPE_A   double
#define TYPE_X   complex_double
#define TYPE_Y   complex_double
#define KL       const int kl,
#define KU       const int ku,
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
suzerain_gbdddmv_ssc(
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
                                    FIXEDBW_CASE, suzerain_gbdddmv_internal_ssc)
        }
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_gbdddmv_internal_ssc(trans, n, kl, ku,
                                         alpha0, d0, ldd0,
                                         alpha1, d1, ldd1,
                                         alpha2, d2, ldd2,
                                         a, lda, x, incx, beta, y, incy);
}

int
suzerain_gbdddmv_scc(
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
                                    FIXEDBW_CASE, suzerain_gbdddmv_internal_scc)
        }
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_gbdddmv_internal_scc(trans, n, kl, ku,
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
suzerain_gbdddmv_ddz(
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
                                    FIXEDBW_CASE, suzerain_gbdddmv_internal_ddz)
        }
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_gbdddmv_internal_ddz(trans, n, kl, ku,
                                         alpha0, d0, ldd0,
                                         alpha1, d1, ldd1,
                                         alpha2, d2, ldd2,
                                         a, lda, x, incx, beta, y, incy);
}

int
suzerain_gbdddmv_dzz(
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
                                    FIXEDBW_CASE, suzerain_gbdddmv_internal_dzz)
        }
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_gbdddmv_internal_dzz(trans, n, kl, ku,
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

#define STATIC   static
#define FUNCTION BOOST_PP_CAT(suzerain_gbdddmv_internal_s, \
                              BOOST_PP_ITERATION())
#define TYPE_D   float
#define TYPE_A   float
#define TYPE_X   float
#define TYPE_Y   float
#define KL
#define KU
#include "gbdddmv.def"

#define STATIC   static
#define FUNCTION BOOST_PP_CAT(suzerain_gbdddmv_internal_ssc, \
                              BOOST_PP_ITERATION())
#define TYPE_D   float
#define TYPE_A   float
#define TYPE_X   float
#define TYPE_Y   complex_float
#define KL
#define KU
#include "gbdddmv.def"

#define STATIC   static
#define FUNCTION BOOST_PP_CAT(suzerain_gbdddmv_internal_scc, \
                              BOOST_PP_ITERATION())
#define TYPE_D   float
#define TYPE_A   float
#define TYPE_X   complex_float
#define TYPE_Y   complex_float
#define KL
#define KU
#include "gbdddmv.def"

#define STATIC   static
#define FUNCTION BOOST_PP_CAT(suzerain_gbdddmv_internal_d, \
                              BOOST_PP_ITERATION())
#define TYPE_D   double
#define TYPE_A   double
#define TYPE_X   double
#define TYPE_Y   double
#define KL
#define KU
#include "gbdddmv.def"

#define STATIC   static
#define FUNCTION BOOST_PP_CAT(suzerain_gbdddmv_internal_ddz, \
                              BOOST_PP_ITERATION())
#define TYPE_D   double
#define TYPE_A   double
#define TYPE_X   double
#define TYPE_Y   complex_double
#define KL
#define KU
#include "gbdddmv.def"

#define STATIC   static
#define FUNCTION BOOST_PP_CAT(suzerain_gbdddmv_internal_dzz, \
                              BOOST_PP_ITERATION())
#define TYPE_D   double
#define TYPE_A   double
#define TYPE_X   complex_double
#define TYPE_Y   complex_double
#define KL
#define KU
#include "gbdddmv.def"

#undef kl
#undef ku

#endif
