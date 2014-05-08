/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2014 Rhys Ulerich
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
 * @copydoc blasext.h
 */

#include <suzerain/blas_et_al/blasext.h>

#include <suzerain/common.h>
#include <suzerain/blas_et_al/blas.h>
#include <suzerain/blas_et_al/gbdddddmv.h>
#include <suzerain/blas_et_al/gbddddmv.h>
#include <suzerain/blas_et_al/gbdddmv.h>
#include <suzerain/blas_et_al/gbddmv.h>
#include <suzerain/blas_et_al/gbdmv.h>
#include <suzerain/blas_et_al/lapack.h>

// Shorthand
static inline int imin(int a, int b) { return a < b ? a : b; }
static inline int imax(int a, int b) { return a > b ? a : b; }
#define UNLIKELY(expr) SUZERAIN_UNLIKELY(expr)

// Many of the short methods have "inline" though their declarations do not.
// to allow inlining them later within this particular translation unit.

// Thank you captain obvious...
#pragma warning(disable:981)

int
suzerain_blasext_sgbdmv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha,
        const float *d,
        const int ldd,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    float * const z = suzerain_blas_malloc(n*sizeof(float));
    if (UNLIKELY(!z)) return suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_sgbmv_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha, d, ldd, z, 1, beta, y, incy);
    suzerain_blas_free(z);
    return 0;
}

int
suzerain_blasext_sgbdmv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha,
        const float *d,
        const int ldd,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    const int info = suzerain_gbdmv_s(trans, n, kl, ku,
                                      alpha, d, ldd,
                                      a, lda, x, incx,
                                      beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

int
suzerain_blasext_dgbdmv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha,
        const double *d,
        const int ldd,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    double * const z = suzerain_blas_malloc(n*sizeof(double));
    if (UNLIKELY(!z)) return suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_dgbmv_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha, d, ldd, z, 1, beta, y, incy);
    suzerain_blas_free(z);
    return 0;
}

inline int
suzerain_blasext_dgbdmv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha,
        const double *d,
        const int ldd,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    const int info = suzerain_gbdmv_d(trans, n, kl, ku,
                                      alpha, d, ldd,
                                      a, lda, x, incx,
                                      beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

int
suzerain_blasext_cgbdmv_s_c_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const complex_float alpha,
        const float *d,
        const int ldd,
        const float *a,
        const int lda,
        const complex_float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy)
{
    complex_float * const z = suzerain_blas_malloc(n*sizeof(complex_float));
    if (UNLIKELY(!z)) return suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_cgbmv_s_c_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_csbmv_s_c_external(
            'U', n, 0, alpha, d, ldd, z, 1, beta, y, incy);
    suzerain_blas_free(z);
    return 0;
}

inline int
suzerain_blasext_cgbdmv_s_c(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const complex_float alpha,
        const float *d,
        const int ldd,
        const float *a,
        const int lda,
        const complex_float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy)
{
    const int info = suzerain_gbdmv_scc(trans, n, kl, ku,
                                        alpha, d, ldd,
                                        a, lda, x, incx,
                                        beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return 0;
}

inline int
suzerain_blasext_cgbdmv_s_s(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const complex_float alpha,
        const float *d,
        const int ldd,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy)
{
    const int info = suzerain_gbdmv_ssc(trans, n, kl, ku,
                                        alpha, d, ldd,
                                        a, lda, x, incx,
                                        beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return 0;
}

int
suzerain_blasext_zgbdmv_d_z_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const complex_double alpha,
        const double *d,
        const int ldd,
        const double *a,
        const int lda,
        const complex_double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy)
{
    complex_double * const z = suzerain_blas_malloc(n*sizeof(complex_double));
    if (UNLIKELY(!z)) return suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_zgbmv_d_z_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_zsbmv_d_z_external(
            'U', n, 0, alpha, d, ldd, z, 1, beta, y, incy);
    suzerain_blas_free(z);
    return 0;
}

inline int
suzerain_blasext_zgbdmv_d_z(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const complex_double alpha,
        const double *d,
        const int ldd,
        const double *a,
        const int lda,
        const complex_double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy)
{
    const int info = suzerain_gbdmv_dzz(trans, n, kl, ku,
                                        alpha, d, ldd,
                                        a, lda, x, incx,
                                        beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

inline int
suzerain_blasext_zgbdmv_d_d(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const complex_double alpha,
        const double *d,
        const int ldd,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy)
{
    const int info = suzerain_gbdmv_ddz(trans, n, kl, ku,
                                        alpha, d, ldd,
                                        a, lda, x, incx,
                                        beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

int
suzerain_blasext_sgbddmv_external(
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
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    float * const z = suzerain_blas_malloc(n*sizeof(float));
    if (UNLIKELY(!z)) return suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_sgbmv_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha0, d0, ldd0, z, 1, beta, y, incy);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha1, d1, ldd1, z, 1, 1,    y, incy);
    suzerain_blas_free(z);
    return 0;
}

inline int
suzerain_blasext_sgbddmv(
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
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    const int info = suzerain_gbddmv_s(trans, n, kl, ku,
                                       alpha0, d0, ldd0,
                                       alpha1, d1, ldd1,
                                       a, lda, x, incx,
                                       beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

int
suzerain_blasext_dgbddmv_external(
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
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    double * const z = suzerain_blas_malloc(n*sizeof(double));
    if (UNLIKELY(!z)) return suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_dgbmv_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha0, d0, ldd0, z, 1, beta, y, incy);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha1, d1, ldd1, z, 1, 1,    y, incy);
    suzerain_blas_free(z);
    return 0;
}

inline int
suzerain_blasext_dgbddmv(
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
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    const int info = suzerain_gbddmv_d(trans, n, kl, ku,
                                       alpha0, d0, ldd0,
                                       alpha1, d1, ldd1,
                                       a, lda, x, incx,
                                       beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

int
suzerain_blasext_cgbddmv_s_c_external(
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
        const float *a,
        const int lda,
        const complex_float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy)
{
    complex_float * const z = suzerain_blas_malloc(n*sizeof(complex_float));
    if (UNLIKELY(!z)) return suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_cgbmv_s_c_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_csbmv_s_c_external(
            'U', n, 0, alpha0, d0, ldd0, z, 1, beta, y, incy);
    suzerain_blas_csbmv_s_c_external(
            'U', n, 0, alpha1, d1, ldd1, z, 1, 1,    y, incy);
    suzerain_blas_free(z);
    return 0;
}

inline int
suzerain_blasext_cgbddmv_s_c(
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
        const float *a,
        const int lda,
        const complex_float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy)
{
    const int info = suzerain_gbddmv_scc(trans, n, kl, ku,
                                         alpha0, d0, ldd0,
                                         alpha1, d1, ldd1,
                                         a, lda, x, incx,
                                         beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

inline int
suzerain_blasext_cgbddmv_s_s(
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
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy)
{
    const int info = suzerain_gbddmv_ssc(trans, n, kl, ku,
                                         alpha0, d0, ldd0,
                                         alpha1, d1, ldd1,
                                         a, lda, x, incx,
                                         beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

int
suzerain_blasext_zgbddmv_d_z_external(
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
        const double *a,
        const int lda,
        const complex_double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy)
{
    complex_double * const z = suzerain_blas_malloc(n*sizeof(complex_double));
    if (UNLIKELY(!z)) return suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_zgbmv_d_z_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_zsbmv_d_z_external(
            'U', n, 0, alpha0, d0, ldd0, z, 1, beta, y, incy);
    suzerain_blas_zsbmv_d_z_external(
            'U', n, 0, alpha1, d1, ldd1, z, 1, 1,    y, incy);
    suzerain_blas_free(z);
    return 0;
}

inline int
suzerain_blasext_zgbddmv_d_z(
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
        const double *a,
        const int lda,
        const complex_double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy)
{
    const int info = suzerain_gbddmv_dzz(trans, n, kl, ku,
                                         alpha0, d0, ldd0,
                                         alpha1, d1, ldd1,
                                         a, lda, x, incx,
                                         beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return 0;
}

inline int
suzerain_blasext_zgbddmv_d_d(
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
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy)
{
    const int info = suzerain_gbddmv_ddz(trans, n, kl, ku,
                                         alpha0, d0, ldd0,
                                         alpha1, d1, ldd1,
                                         a, lda, x, incx,
                                         beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return 0;
}

int
suzerain_blasext_sgbdddmv_external(
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
    float * const z = suzerain_blas_malloc(n*sizeof(float));
    if (UNLIKELY(!z)) return suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_sgbmv_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha0, d0, ldd0, z, 1, beta, y, incy);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha1, d1, ldd1, z, 1, 1,    y, incy);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha2, d2, ldd2, z, 1, 1,    y, incy);
    suzerain_blas_free(z);
    return 0;
}

inline int
suzerain_blasext_sgbdddmv(
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
    const int info = suzerain_gbdddmv_s(trans, n, kl, ku,
                                        alpha0, d0, ldd0,
                                        alpha1, d1, ldd1,
                                        alpha2, d2, ldd2,
                                        a, lda, x, incx,
                                        beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

int
suzerain_blasext_dgbdddmv_external(
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
    double * const z = suzerain_blas_malloc(n*sizeof(double));
    if (UNLIKELY(!z)) return suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_dgbmv_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha0, d0, ldd0, z, 1, beta, y, incy);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha1, d1, ldd1, z, 1, 1,    y, incy);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha2, d2, ldd2, z, 1, 1,    y, incy);
    suzerain_blas_free(z);
    return 0;
}

inline int
suzerain_blasext_dgbdddmv(
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
    const int info = suzerain_gbdddmv_d(trans, n, kl, ku,
                                        alpha0, d0, ldd0,
                                        alpha1, d1, ldd1,
                                        alpha2, d2, ldd2,
                                        a, lda, x, incx,
                                        beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

int
suzerain_blasext_cgbdddmv_s_c_external(
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
    complex_float * const z = suzerain_blas_malloc(n*sizeof(complex_float));
    if (UNLIKELY(!z)) return suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_cgbmv_s_c_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_csbmv_s_c_external(
            'U', n, 0, alpha0, d0, ldd0, z, 1, beta, y, incy);
    suzerain_blas_csbmv_s_c_external(
            'U', n, 0, alpha1, d1, ldd1, z, 1, 1,    y, incy);
    suzerain_blas_csbmv_s_c_external(
            'U', n, 0, alpha2, d2, ldd2, z, 1, 1,    y, incy);
    suzerain_blas_free(z);
    return 0;
}

inline int
suzerain_blasext_cgbdddmv_s_c(
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
    const int info = suzerain_gbdddmv_scc(trans, n, kl, ku,
                                          alpha0, d0, ldd0,
                                          alpha1, d1, ldd1,
                                          alpha2, d2, ldd2,
                                          a, lda, x, incx,
                                          beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

inline int
suzerain_blasext_cgbdddmv_s_s(
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
    const int info = suzerain_gbdddmv_ssc(trans, n, kl, ku,
                                          alpha0, d0, ldd0,
                                          alpha1, d1, ldd1,
                                          alpha2, d2, ldd2,
                                          a, lda, x, incx,
                                          beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

int
suzerain_blasext_zgbdddmv_d_z_external(
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
    complex_double * const z = suzerain_blas_malloc(n*sizeof(complex_double));
    if (UNLIKELY(!z)) return suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_zgbmv_d_z_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_zsbmv_d_z_external(
            'U', n, 0, alpha0, d0, ldd0, z, 1, beta, y, incy);
    suzerain_blas_zsbmv_d_z_external(
            'U', n, 0, alpha1, d1, ldd1, z, 1, 1,    y, incy);
    suzerain_blas_zsbmv_d_z_external(
            'U', n, 0, alpha2, d2, ldd2, z, 1, 1,    y, incy);
    suzerain_blas_free(z);
    return 0;
}

inline int
suzerain_blasext_zgbdddmv_d_z(
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
    const int info = suzerain_gbdddmv_dzz(trans, n, kl, ku,
                                          alpha0, d0, ldd0,
                                          alpha1, d1, ldd1,
                                          alpha2, d2, ldd2,
                                          a, lda, x, incx,
                                          beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

inline int
suzerain_blasext_zgbdddmv_d_d(
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
    const int info = suzerain_gbdddmv_ddz(trans, n, kl, ku,
                                          alpha0, d0, ldd0,
                                          alpha1, d1, ldd1,
                                          alpha2, d2, ldd2,
                                          a, lda, x, incx,
                                          beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

int
suzerain_blasext_sgbddddmv_external(
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
        const float alpha3,
        const float *d3,
        const int ldd3,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    float * const z = suzerain_blas_malloc(n*sizeof(float));
    if (UNLIKELY(!z)) return suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_sgbmv_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha0, d0, ldd0, z, 1, beta, y, incy);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha1, d1, ldd1, z, 1, 1,    y, incy);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha2, d2, ldd2, z, 1, 1,    y, incy);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha3, d3, ldd3, z, 1, 1,    y, incy);
    suzerain_blas_free(z);
    return 0;
}

inline int
suzerain_blasext_sgbddddmv(
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
        const float alpha3,
        const float *d3,
        const int ldd3,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    const int info = suzerain_gbddddmv_s(trans, n, kl, ku,
                                         alpha0, d0, ldd0,
                                         alpha1, d1, ldd1,
                                         alpha2, d2, ldd2,
                                         alpha3, d3, ldd3,
                                         a, lda, x, incx,
                                         beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

int
suzerain_blasext_dgbddddmv_external(
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
        const double alpha3,
        const double *d3,
        const int ldd3,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    double * const z = suzerain_blas_malloc(n*sizeof(double));
    if (UNLIKELY(!z)) return suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_dgbmv_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha0, d0, ldd0, z, 1, beta, y, incy);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha1, d1, ldd1, z, 1, 1,    y, incy);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha2, d2, ldd2, z, 1, 1,    y, incy);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha3, d3, ldd3, z, 1, 1,    y, incy);
    suzerain_blas_free(z);
    return 0;
}

inline int
suzerain_blasext_dgbddddmv(
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
        const double alpha3,
        const double *d3,
        const int ldd3,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    const int info = suzerain_gbddddmv_d(trans, n, kl, ku,
                                        alpha0, d0, ldd0,
                                        alpha1, d1, ldd1,
                                        alpha2, d2, ldd2,
                                        alpha3, d3, ldd3,
                                        a, lda, x, incx,
                                        beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

int
suzerain_blasext_cgbddddmv_s_c_external(
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
        const complex_float alpha3,
        const float *d3,
        const int ldd3,
        const float *a,
        const int lda,
        const complex_float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy)
{
    complex_float * const z = suzerain_blas_malloc(n*sizeof(complex_float));
    if (UNLIKELY(!z)) return suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_cgbmv_s_c_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_csbmv_s_c_external(
            'U', n, 0, alpha0, d0, ldd0, z, 1, beta, y, incy);
    suzerain_blas_csbmv_s_c_external(
            'U', n, 0, alpha1, d1, ldd1, z, 1, 1,    y, incy);
    suzerain_blas_csbmv_s_c_external(
            'U', n, 0, alpha2, d2, ldd2, z, 1, 1,    y, incy);
    suzerain_blas_csbmv_s_c_external(
            'U', n, 0, alpha3, d3, ldd3, z, 1, 1,    y, incy);
    suzerain_blas_free(z);
    return 0;
}

inline int
suzerain_blasext_cgbddddmv_s_c(
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
        const complex_float alpha3,
        const float *d3,
        const int ldd3,
        const float *a,
        const int lda,
        const complex_float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy)
{
    const int info = suzerain_gbddddmv_scc(trans, n, kl, ku,
                                           alpha0, d0, ldd0,
                                           alpha1, d1, ldd1,
                                           alpha2, d2, ldd2,
                                           alpha3, d3, ldd3,
                                           a, lda, x, incx,
                                           beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

inline int
suzerain_blasext_cgbddddmv_s_s(
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
        const complex_float alpha3,
        const float *d3,
        const int ldd3,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy)
{
    const int info = suzerain_gbddddmv_ssc(trans, n, kl, ku,
                                           alpha0, d0, ldd0,
                                           alpha1, d1, ldd1,
                                           alpha2, d2, ldd2,
                                           alpha3, d3, ldd3,
                                           a, lda, x, incx,
                                           beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

int
suzerain_blasext_zgbddddmv_d_z_external(
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
        const complex_double alpha3,
        const double *d3,
        const int ldd3,
        const double *a,
        const int lda,
        const complex_double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy)
{
    complex_double * const z = suzerain_blas_malloc(n*sizeof(complex_double));
    if (UNLIKELY(!z)) return suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_zgbmv_d_z_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_zsbmv_d_z_external(
            'U', n, 0, alpha0, d0, ldd0, z, 1, beta, y, incy);
    suzerain_blas_zsbmv_d_z_external(
            'U', n, 0, alpha1, d1, ldd1, z, 1, 1,    y, incy);
    suzerain_blas_zsbmv_d_z_external(
            'U', n, 0, alpha2, d2, ldd2, z, 1, 1,    y, incy);
    suzerain_blas_zsbmv_d_z_external(
            'U', n, 0, alpha3, d3, ldd3, z, 1, 1,    y, incy);
    suzerain_blas_free(z);
    return 0;
}

inline int
suzerain_blasext_zgbddddmv_d_z(
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
        const complex_double alpha3,
        const double *d3,
        const int ldd3,
        const double *a,
        const int lda,
        const complex_double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy)
{
    const int info = suzerain_gbddddmv_dzz(trans, n, kl, ku,
                                           alpha0, d0, ldd0,
                                           alpha1, d1, ldd1,
                                           alpha2, d2, ldd2,
                                           alpha3, d3, ldd3,
                                           a, lda, x, incx,
                                           beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

inline int
suzerain_blasext_zgbddddmv_d_d(
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
        const complex_double alpha3,
        const double *d3,
        const int ldd3,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy)
{
    const int info = suzerain_gbddddmv_ddz(trans, n, kl, ku,
                                           alpha0, d0, ldd0,
                                           alpha1, d1, ldd1,
                                           alpha2, d2, ldd2,
                                           alpha3, d3, ldd3,
                                           a, lda, x, incx,
                                           beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

inline int
suzerain_blasext_sgbdddddmv_external(
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
        const float alpha3,
        const float *d3,
        const int ldd3,
        const float alpha4,
        const float *d4,
        const int ldd4,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    float * const z = suzerain_blas_malloc(n*sizeof(float));
    if (UNLIKELY(!z)) return suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_sgbmv_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha0, d0, ldd0, z, 1, beta, y, incy);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha1, d1, ldd1, z, 1, 1,    y, incy);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha2, d2, ldd2, z, 1, 1,    y, incy);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha3, d3, ldd3, z, 1, 1,    y, incy);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha4, d4, ldd4, z, 1, 1,    y, incy);
    suzerain_blas_free(z);
    return 0;
}

inline int
suzerain_blasext_sgbdddddmv(
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
        const float alpha3,
        const float *d3,
        const int ldd3,
        const float alpha4,
        const float *d4,
        const int ldd4,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    const int info = suzerain_gbdddddmv_s(trans, n, kl, ku,
                                          alpha0, d0, ldd0,
                                          alpha1, d1, ldd1,
                                          alpha2, d2, ldd2,
                                          alpha3, d3, ldd3,
                                          alpha4, d4, ldd4,
                                          a, lda, x, incx,
                                          beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

int
suzerain_blasext_dgbdddddmv_external(
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
        const double alpha3,
        const double *d3,
        const int ldd3,
        const double alpha4,
        const double *d4,
        const int ldd4,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    double * const z = suzerain_blas_malloc(n*sizeof(double));
    if (UNLIKELY(!z)) return suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_dgbmv_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha0, d0, ldd0, z, 1, beta, y, incy);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha1, d1, ldd1, z, 1, 1,    y, incy);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha2, d2, ldd2, z, 1, 1,    y, incy);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha3, d3, ldd3, z, 1, 1,    y, incy);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha4, d4, ldd4, z, 1, 1,    y, incy);
    suzerain_blas_free(z);
    return 0;
}

inline int
suzerain_blasext_dgbdddddmv(
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
        const double alpha3,
        const double *d3,
        const int ldd3,
        const double alpha4,
        const double *d4,
        const int ldd4,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    const int info = suzerain_gbdddddmv_d(trans, n, kl, ku,
                                          alpha0, d0, ldd0,
                                          alpha1, d1, ldd1,
                                          alpha2, d2, ldd2,
                                          alpha3, d3, ldd3,
                                          alpha4, d4, ldd4,
                                          a, lda, x, incx,
                                          beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

int
suzerain_blasext_cgbdddddmv_s_c_external(
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
        const complex_float alpha3,
        const float *d3,
        const int ldd3,
        const complex_float alpha4,
        const float *d4,
        const int ldd4,
        const float *a,
        const int lda,
        const complex_float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy)
{
    complex_float * const z = suzerain_blas_malloc(n*sizeof(complex_float));
    if (UNLIKELY(!z)) return suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_cgbmv_s_c_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_csbmv_s_c_external(
            'U', n, 0, alpha0, d0, ldd0, z, 1, beta, y, incy);
    suzerain_blas_csbmv_s_c_external(
            'U', n, 0, alpha1, d1, ldd1, z, 1, 1,    y, incy);
    suzerain_blas_csbmv_s_c_external(
            'U', n, 0, alpha2, d2, ldd2, z, 1, 1,    y, incy);
    suzerain_blas_csbmv_s_c_external(
            'U', n, 0, alpha3, d3, ldd3, z, 1, 1,    y, incy);
    suzerain_blas_csbmv_s_c_external(
            'U', n, 0, alpha4, d4, ldd4, z, 1, 1,    y, incy);
    suzerain_blas_free(z);
    return 0;
}

inline int
suzerain_blasext_cgbdddddmv_s_c(
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
        const complex_float alpha3,
        const float *d3,
        const int ldd3,
        const complex_float alpha4,
        const float *d4,
        const int ldd4,
        const float *a,
        const int lda,
        const complex_float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy)
{
    const int info = suzerain_gbdddddmv_scc(trans, n, kl, ku,
                                            alpha0, d0, ldd0,
                                            alpha1, d1, ldd1,
                                            alpha2, d2, ldd2,
                                            alpha3, d3, ldd3,
                                            alpha4, d4, ldd4,
                                            a, lda, x, incx,
                                            beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

inline int
suzerain_blasext_cgbdddddmv_s_s(
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
        const complex_float alpha3,
        const float *d3,
        const int ldd3,
        const complex_float alpha4,
        const float *d4,
        const int ldd4,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy)
{
    const int info = suzerain_gbdddddmv_ssc(trans, n, kl, ku,
                                            alpha0, d0, ldd0,
                                            alpha1, d1, ldd1,
                                            alpha2, d2, ldd2,
                                            alpha3, d3, ldd3,
                                            alpha4, d4, ldd4,
                                            a, lda, x, incx,
                                            beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

int
suzerain_blasext_zgbdddddmv_d_z_external(
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
        const complex_double alpha3,
        const double *d3,
        const int ldd3,
        const complex_double alpha4,
        const double *d4,
        const int ldd4,
        const double *a,
        const int lda,
        const complex_double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy)
{
    complex_double * const z = suzerain_blas_malloc(n*sizeof(complex_double));
    if (UNLIKELY(!z)) return suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_zgbmv_d_z_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_zsbmv_d_z_external(
            'U', n, 0, alpha0, d0, ldd0, z, 1, beta, y, incy);
    suzerain_blas_zsbmv_d_z_external(
            'U', n, 0, alpha1, d1, ldd1, z, 1, 1,    y, incy);
    suzerain_blas_zsbmv_d_z_external(
            'U', n, 0, alpha2, d2, ldd2, z, 1, 1,    y, incy);
    suzerain_blas_zsbmv_d_z_external(
            'U', n, 0, alpha3, d3, ldd3, z, 1, 1,    y, incy);
    suzerain_blas_zsbmv_d_z_external(
            'U', n, 0, alpha4, d4, ldd4, z, 1, 1,    y, incy);
    suzerain_blas_free(z);
    return 0;
}

inline int
suzerain_blasext_zgbdddddmv_d_z(
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
        const complex_double alpha3,
        const double *d3,
        const int ldd3,
        const complex_double alpha4,
        const double *d4,
        const int ldd4,
        const double *a,
        const int lda,
        const complex_double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy)
{
    const int info = suzerain_gbdddddmv_dzz(trans, n, kl, ku,
                                            alpha0, d0, ldd0,
                                            alpha1, d1, ldd1,
                                            alpha2, d2, ldd2,
                                            alpha3, d3, ldd3,
                                            alpha4, d4, ldd4,
                                            a, lda, x, incx,
                                            beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

inline int
suzerain_blasext_zgbdddddmv_d_d(
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
        const complex_double alpha3,
        const double *d3,
        const int ldd3,
        const complex_double alpha4,
        const double *d4,
        const int ldd4,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy)
{
    const int info = suzerain_gbdddddmv_ddz(trans, n, kl, ku,
                                            alpha0, d0, ldd0,
                                            alpha1, d1, ldd1,
                                            alpha2, d2, ldd2,
                                            alpha3, d3, ldd3,
                                            alpha4, d4, ldd4,
                                            a, lda, x, incx,
                                            beta,   y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

int
suzerain_blasext_sgb_diag_scale_acc(
        char side,
        int m,
        int n,
        int kl,
        int ku,
        const float alpha,
        const float *d,
        int ldd,
        const float *a,
        int inca,
        int lda,
        const float beta,
        float *b,
        int incb,
        int ldb)
{
    if (UNLIKELY(kl   < 0       )) return suzerain_blas_xerbla(__func__,  3);
    if (UNLIKELY(ku   < 0       )) return suzerain_blas_xerbla(__func__,  4);
    if (UNLIKELY(inca < 1       )) return suzerain_blas_xerbla(__func__,  7);
    if (UNLIKELY(lda  <= kl + ku)) return suzerain_blas_xerbla(__func__,  8);
    if (UNLIKELY(ldd  < 0       )) return suzerain_blas_xerbla(__func__, 10);
    if (UNLIKELY(incb < 1       )) return suzerain_blas_xerbla(__func__, 13);
    if (UNLIKELY(ldb  <= kl + ku)) return suzerain_blas_xerbla(__func__, 14);

#pragma warning(push,disable:1572)
    const _Bool alpha_is_zero = (alpha == 0.0f);
    const _Bool beta_is_one   = (beta  == 1.0f);
#pragma warning(pop)

    // If necessary, recast side == 'L' details into a side == 'R' traversal
    switch (toupper(side)) {
    case 'R': break;
    case 'L': inca = lda - inca; a += ku - kl*inca;  // Traverse A by rows
              incb = ldb - incb; b += ku - kl*incb;  // Traverse B by rows
              kl ^= ku; ku ^= kl; kl ^= ku;          // Swap kl and ku
              break;
    default:  return suzerain_blas_xerbla(__func__, 1);
    }

    // Quick return if possible
    if (UNLIKELY((alpha_is_zero && beta_is_one) || m <= 0 || n <= 0))
        return 0;

    // Banded matrix dereference has form a[(ku + i)*inca + j*(lda - inca)].
    // Incorporate the ku offset and decrement ldX to speed indexing in loops.
    // Further, increment kl anticipating expressions like imin(m, j + kl + 1).
    a += ku*inca; lda -= inca;
    b += ku*incb; ldb -= incb;
    ++kl;

    if (beta_is_one) {
        for (int j = 0; j < n; a += lda, d += ldd, b += ldb, ++j) {
            const int il = imax(0, j - ku);
            const int iu = imin(m, j + kl);
            suzerain_blas_saxpy (iu - il, alpha*(*d), a + il*inca, inca,
                                                      b + il*incb, incb);
        }
    } else {
        for (int j = 0; j < n; a += lda, d += ldd, b += ldb, ++j) {
            const int il = imax(0, j - ku);
            const int iu = imin(m, j + kl);
            suzerain_blas_saxpby(iu - il, alpha*(*d), a + il*inca, inca,
                                          beta,       b + il*incb, incb);
        }
    }

    return 0;
}

int
suzerain_blasext_dgb_diag_scale_acc(
        char side,
        int m,
        int n,
        int kl,
        int ku,
        const double alpha,
        const double *d,
        int ldd,
        const double *a,
        int inca,
        int lda,
        const double beta,
        double *b,
        int incb,
        int ldb)
{
    if (UNLIKELY(kl   < 0       )) return suzerain_blas_xerbla(__func__,  4);
    if (UNLIKELY(ku   < 0       )) return suzerain_blas_xerbla(__func__,  5);
    if (UNLIKELY(ldd  < 0       )) return suzerain_blas_xerbla(__func__,  8);
    if (UNLIKELY(inca < 1       )) return suzerain_blas_xerbla(__func__, 10);
    if (UNLIKELY(lda  <= kl + ku)) return suzerain_blas_xerbla(__func__, 11);
    if (UNLIKELY(incb < 1       )) return suzerain_blas_xerbla(__func__, 14);
    if (UNLIKELY(ldb  <= kl + ku)) return suzerain_blas_xerbla(__func__, 15);

#pragma warning(push,disable:1572)
    const _Bool alpha_is_zero = (alpha == 0.0);
    const _Bool beta_is_one   = (beta  == 1.0);
#pragma warning(pop)

    // If necessary, recast side == 'L' details into a side == 'R' traversal
    switch (toupper(side)) {
    case 'R': break;
    case 'L': inca = lda - inca; a += ku - kl*inca;  // Traverse A by rows
              incb = ldb - incb; b += ku - kl*incb;  // Traverse B by rows
              kl ^= ku; ku ^= kl; kl ^= ku;          // Swap kl and ku
              break;
    default:  return suzerain_blas_xerbla(__func__, 1);
    }

    // Quick return if possible
    if (UNLIKELY((alpha_is_zero && beta_is_one) || m <= 0 || n <= 0))
        return 0;

    // Banded matrix dereference has form a[(ku + i)*inca + j*(lda - inca)].
    // Incorporate the ku offset and decrement ldX to speed indexing in loops.
    // Further, increment kl anticipating expressions like imin(m, j + kl + 1).
    a += ku*inca; lda -= inca;
    b += ku*incb; ldb -= incb;
    ++kl;

    if (beta_is_one) {
        for (int j = 0; j < n; a += lda, d += ldd, b += ldb, ++j) {
            const int il = imax(0, j - ku);
            const int iu = imin(m, j + kl);
            suzerain_blas_daxpy (iu - il, alpha*(*d), a + il*inca, inca,
                                                      b + il*incb, incb);
        }
    } else {
        for (int j = 0; j < n; a += lda, d += ldd, b += ldb, ++j) {
            const int il = imax(0, j - ku);
            const int iu = imin(m, j + kl);
            suzerain_blas_daxpby(iu - il, alpha*(*d), a + il*inca, inca,
                                          beta,       b + il*incb, incb);
        }
    }

    return 0;
}

int
suzerain_blasext_cgb_diag_scale_acc(
        char side,
        int m,
        int n,
        int kl,
        int ku,
        const complex_float alpha,
        const complex_float *d,
        int ldd,
        const complex_float *a,
        int inca,
        int lda,
        const complex_float beta,
        complex_float *b,
        int incb,
        int ldb)
{
    if (UNLIKELY(kl   < 0       )) return suzerain_blas_xerbla(__func__,  4);
    if (UNLIKELY(ku   < 0       )) return suzerain_blas_xerbla(__func__,  5);
    if (UNLIKELY(ldd  < 0       )) return suzerain_blas_xerbla(__func__,  8);
    if (UNLIKELY(inca < 1       )) return suzerain_blas_xerbla(__func__, 10);
    if (UNLIKELY(lda  <= kl + ku)) return suzerain_blas_xerbla(__func__, 11);
    if (UNLIKELY(incb < 1       )) return suzerain_blas_xerbla(__func__, 14);
    if (UNLIKELY(ldb  <= kl + ku)) return suzerain_blas_xerbla(__func__, 15);

#pragma warning(push,disable:1572)
    const _Bool alpha_is_zero = (alpha == 0.0f);
    const _Bool beta_is_one   = (beta  == 1.0f);
#pragma warning(pop)

    // If necessary, recast side == 'L' details into a side == 'R' traversal
    switch (toupper(side)) {
    case 'R': break;
    case 'L': inca = lda - inca; a += ku - kl*inca;  // Traverse A by rows
              incb = ldb - incb; b += ku - kl*incb;  // Traverse B by rows
              kl ^= ku; ku ^= kl; kl ^= ku;          // Swap kl and ku
              break;
    default:  return suzerain_blas_xerbla(__func__, 1);
    }

    // Quick return if possible
    if (UNLIKELY((alpha_is_zero && beta_is_one) || m <= 0 || n <= 0))
        return 0;

    // Banded matrix dereference has form a[(ku + i)*inca + j*(lda - inca)].
    // Incorporate the ku offset and decrement ldX to speed indexing in loops.
    // Further, increment kl anticipating expressions like imin(m, j + kl + 1).
    a += ku*inca; lda -= inca;
    b += ku*incb; ldb -= incb;
    ++kl;

    if (beta_is_one) {
        for (int j = 0; j < n; a += lda, d += ldd, b += ldb, ++j) {
            const int il = imax(0, j - ku);
            const int iu = imin(m, j + kl);
            const complex_float tmp = alpha*(*d);
            suzerain_blas_caxpy (iu - il, tmp,  a + il*inca, inca,
                                                b + il*incb, incb);
        }
    } else {
        for (int j = 0; j < n; a += lda, d += ldd, b += ldb, ++j) {
            const int il = imax(0, j - ku);
            const int iu = imin(m, j + kl);
            const complex_float tmp = alpha*(*d);
            suzerain_blas_caxpby(iu - il, tmp,  a + il*inca, inca,
                                          beta, b + il*incb, incb);
        }
    }

    return 0;
}

int
suzerain_blasext_zgb_diag_scale_acc(
        char side,
        int m,
        int n,
        int kl,
        int ku,
        const complex_double alpha,
        const complex_double *d,
        int ldd,
        const complex_double *a,
        int inca,
        int lda,
        const complex_double beta,
        complex_double *b,
        int incb,
        int ldb)
{
    if (UNLIKELY(kl   < 0       )) return suzerain_blas_xerbla(__func__,  4);
    if (UNLIKELY(ku   < 0       )) return suzerain_blas_xerbla(__func__,  5);
    if (UNLIKELY(ldd  < 0       )) return suzerain_blas_xerbla(__func__,  8);
    if (UNLIKELY(inca < 1       )) return suzerain_blas_xerbla(__func__, 10);
    if (UNLIKELY(lda  <= kl + ku)) return suzerain_blas_xerbla(__func__, 11);
    if (UNLIKELY(incb < 1       )) return suzerain_blas_xerbla(__func__, 14);
    if (UNLIKELY(ldb  <= kl + ku)) return suzerain_blas_xerbla(__func__, 15);

#pragma warning(push,disable:1572)
    const _Bool alpha_is_zero = (alpha == 0.0);
    const _Bool beta_is_one   = (beta  == 1.0);
#pragma warning(pop)

    // If necessary, recast side == 'L' details into a side == 'R' traversal
    switch (toupper(side)) {
    case 'R': break;
    case 'L': inca = lda - inca; a += ku - kl*inca;  // Traverse A by rows
              incb = ldb - incb; b += ku - kl*incb;  // Traverse B by rows
              kl ^= ku; ku ^= kl; kl ^= ku;          // Swap kl and ku
              break;
    default:  return suzerain_blas_xerbla(__func__, 1);
    }

    // Quick return if possible
    if (UNLIKELY((alpha_is_zero && beta_is_one) || m <= 0 || n <= 0))
        return 0;

    // Banded matrix dereference has form a[(ku + i)*inca + j*(lda - inca)].
    // Incorporate the ku offset and decrement ldX to speed indexing in loops.
    // Further, increment kl anticipating expressions like imin(m, j + kl + 1).
    a += ku*inca; lda -= inca;
    b += ku*incb; ldb -= incb;
    ++kl;

    if (beta_is_one) {
        for (int j = 0; j < n; a += lda, d += ldd, b += ldb, ++j) {
            const int il = imax(0, j - ku);
            const int iu = imin(m, j + kl);
            const complex_double tmp = alpha*(*d);
            suzerain_blas_zaxpy (iu - il, tmp,  a + il*inca, inca,
                                                b + il*incb, incb);
        }
    } else {
        for (int j = 0; j < n; a += lda, d += ldd, b += ldb, ++j) {
            const int il = imax(0, j - ku);
            const int iu = imin(m, j + kl);
            const complex_double tmp = alpha*(*d);
            suzerain_blas_zaxpby(iu - il, tmp,  a + il*inca, inca,
                                          beta, b + il*incb, incb);
        }
    }

    return 0;
}

int
suzerain_blasext_zgb_diag_scale_acc_d(
        char side,
        int m,
        int n,
        int kl,
        int ku,
        const complex_double alpha,
        const double *d,
        int ldd,
        const double *a,
        int inca,
        int lda,
        const complex_double beta,
        complex_double *b,
        int incb,
        int ldb)
{
    if (UNLIKELY(kl   < 0       )) return suzerain_blas_xerbla(__func__,  4);
    if (UNLIKELY(ku   < 0       )) return suzerain_blas_xerbla(__func__,  5);
    if (UNLIKELY(ldd  < 0       )) return suzerain_blas_xerbla(__func__,  8);
    if (UNLIKELY(inca < 1       )) return suzerain_blas_xerbla(__func__, 10);
    if (UNLIKELY(lda  <= kl + ku)) return suzerain_blas_xerbla(__func__, 11);
    if (UNLIKELY(incb < 1       )) return suzerain_blas_xerbla(__func__, 14);
    if (UNLIKELY(ldb  <= kl + ku)) return suzerain_blas_xerbla(__func__, 15);

#pragma warning(push,disable:1572)
    const _Bool alpha_is_zero = (alpha == 0.0);
    const _Bool beta_is_one   = (beta  == 1.0);
#pragma warning(pop)

    // If necessary, recast side == 'L' details into a side == 'R' traversal
    switch (toupper(side)) {
    case 'R': break;
    case 'L': inca = lda - inca; a += ku - kl*inca;  // Traverse A by rows
              incb = ldb - incb; b += ku - kl*incb;  // Traverse B by rows
              kl ^= ku; ku ^= kl; kl ^= ku;          // Swap kl and ku
              break;
    default:  return suzerain_blas_xerbla(__func__, 1);
    }

    // Quick return if possible
    if (UNLIKELY((alpha_is_zero && beta_is_one) || m <= 0 || n <= 0))
        return 0;

    // Banded matrix dereference has form a[(ku + i)*inca + j*(lda - inca)].
    // Incorporate the ku offset and decrement ldX to speed indexing in loops.
    // Further, increment kl anticipating expressions like imin(m, j + kl + 1).
    a += ku*inca; lda -= inca;
    b += ku*incb; ldb -= incb;
    ++kl;

    if (beta_is_one) {
        for (int j = 0; j < n; a += lda, d += ldd, b += ldb, ++j) {
            const int il = imax(0, j - ku);
            const int iu = imin(m, j + kl);
            const complex_double tmp = alpha*(*d);
            suzerain_blas_zaxpy_d (iu - il, tmp,  a + il*inca, inca,
                                                  b + il*incb, incb);
        }
    } else {
        for (int j = 0; j < n; a += lda, d += ldd, b += ldb, ++j) {
            const int il = imax(0, j - ku);
            const int iu = imin(m, j + kl);
            const complex_double tmp = alpha*(*d);
            suzerain_blas_zaxpby_d(iu - il, tmp,  a + il*inca, inca,
                                            beta, b + il*incb, incb);
        }
    }

    return 0;
}

int
suzerain_blasext_zgb_ddiag_scale_acc_d(
        char side,
        int m,
        int n,
        int kl,
        int ku,
        const complex_double alpha0,
        const double *d0,
        int ldd0,
        const complex_double alpha1,
        const double *d1,
        int ldd1,
        const double *a,
        int inca,
        int lda,
        const complex_double beta,
        complex_double *b,
        int incb,
        int ldb)
{
    if (UNLIKELY(kl   < 0       )) return suzerain_blas_xerbla(__func__,  4);
    if (UNLIKELY(ku   < 0       )) return suzerain_blas_xerbla(__func__,  5);
    if (UNLIKELY(ldd0 < 0       )) return suzerain_blas_xerbla(__func__,  8);
    if (UNLIKELY(ldd1 < 0       )) return suzerain_blas_xerbla(__func__, 11);
    if (UNLIKELY(inca < 1       )) return suzerain_blas_xerbla(__func__, 13);
    if (UNLIKELY(lda  <= kl + ku)) return suzerain_blas_xerbla(__func__, 14);
    if (UNLIKELY(incb < 1       )) return suzerain_blas_xerbla(__func__, 17);
    if (UNLIKELY(lda  <= kl + ku)) return suzerain_blas_xerbla(__func__, 18);

#pragma warning(push,disable:1572)
    const _Bool alpha0_is_zero = (alpha0 == 0.0);
    const _Bool alpha1_is_zero = (alpha1 == 0.0);
    const _Bool beta_is_one    = (beta   == 1.0);
#pragma warning(pop)

    // If necessary, recast side == 'L' details into a side == 'R' traversal
    switch (toupper(side)) {
    case 'R': break;
    case 'L': inca = lda - inca; a += ku - kl*inca;  // Traverse A by rows
              incb = ldb - incb; b += ku - kl*incb;  // Traverse B by rows
              kl ^= ku; ku ^= kl; kl ^= ku;          // Swap kl and ku
              break;
    default:  return suzerain_blas_xerbla(__func__, 1);
    }

    // Quick return if possible
    if (UNLIKELY(    (alpha0_is_zero && alpha1_is_zero && beta_is_one)
                  || m <= 0 || n <= 0))
        return 0;

    // Banded matrix dereference has form a[(ku + i)*inca + j*(lda - inca)].
    // Incorporate the ku offset and decrement ldX to speed indexing in loops.
    // Further, increment kl anticipating expressions like imin(m, j + kl + 1).
    a += ku*inca; lda -= inca;
    b += ku*incb; ldb -= incb;
    ++kl;

    if (beta_is_one) {
        for (int j = 0; j < n;
             a += lda, d0 += ldd0, d1 += ldd1, b += ldb, ++j) {
            const int il = imax(0, j - ku);
            const int iu = imin(m, j + kl);
            const complex_double tmp = alpha0*(*d0) + alpha1*(*d1);
            suzerain_blas_zaxpy_d (iu - il, tmp,  a + il*inca, inca,
                                                  b + il*incb, incb);
        }
    } else {
        for (int j = 0; j < n;
             a += lda, d0 += ldd0, d1 += ldd1, b += ldb, ++j) {
            const int il = imax(0, j - ku);
            const int iu = imin(m, j + kl);
            const complex_double tmp = alpha0*(*d0) + alpha1*(*d1);
            suzerain_blas_zaxpby_d(iu - il, tmp,  a + il*inca, inca,
                                            beta, b + il*incb, incb);
        }
    }

    return 0;
}

int
suzerain_blasext_zgb_dddiag_scale_acc_d(
        char side,
        int m,
        int n,
        int kl,
        int ku,
        const complex_double alpha0,
        const double *d0,
        int ldd0,
        const complex_double alpha1,
        const double *d1,
        int ldd1,
        const complex_double alpha2,
        const double *d2,
        int ldd2,
        const double *a,
        int inca,
        int lda,
        const complex_double beta,
        complex_double *b,
        int incb,
        int ldb)
{
    if (UNLIKELY(kl   < 0       )) return suzerain_blas_xerbla(__func__,  4);
    if (UNLIKELY(ku   < 0       )) return suzerain_blas_xerbla(__func__,  5);
    if (UNLIKELY(ldd0 < 0       )) return suzerain_blas_xerbla(__func__,  8);
    if (UNLIKELY(ldd1 < 0       )) return suzerain_blas_xerbla(__func__, 11);
    if (UNLIKELY(ldd2 < 0       )) return suzerain_blas_xerbla(__func__, 14);
    if (UNLIKELY(inca < 1       )) return suzerain_blas_xerbla(__func__, 16);
    if (UNLIKELY(lda  <= kl + ku)) return suzerain_blas_xerbla(__func__, 17);
    if (UNLIKELY(incb < 1       )) return suzerain_blas_xerbla(__func__, 20);
    if (UNLIKELY(lda  <= kl + ku)) return suzerain_blas_xerbla(__func__, 21);

#pragma warning(push,disable:1572)
    const _Bool alpha0_is_zero = (alpha0 == 0.0);
    const _Bool alpha1_is_zero = (alpha1 == 0.0);
    const _Bool alpha2_is_zero = (alpha2 == 0.0);
    const _Bool beta_is_one    = (beta   == 1.0);
#pragma warning(pop)

    // If necessary, recast side == 'L' details into a side == 'R' traversal
    switch (toupper(side)) {
    case 'R': break;
    case 'L': inca = lda - inca; a += ku - kl*inca;  // Traverse A by rows
              incb = ldb - incb; b += ku - kl*incb;  // Traverse B by rows
              kl ^= ku; ku ^= kl; kl ^= ku;          // Swap kl and ku
              break;
    default:  return suzerain_blas_xerbla(__func__, 1);
    }

    // Quick return if possible
    if (UNLIKELY(    (alpha0_is_zero && alpha1_is_zero && alpha2_is_zero
                                     && beta_is_one)
                  || m <= 0 || n <= 0))
        return 0;

    // Banded matrix dereference has form a[(ku + i)*inca + j*(lda - inca)].
    // Incorporate the ku offset and decrement ldX to speed indexing in loops.
    // Further, increment kl anticipating expressions like imin(m, j + kl + 1).
    a += ku*inca; lda -= inca;
    b += ku*incb; ldb -= incb;
    ++kl;

    if (beta_is_one) {
        for (int j = 0; j < n;
             a += lda, d0 += ldd0, d1 += ldd1, d2 += ldd2, b += ldb, ++j) {
            const int il = imax(0, j - ku);
            const int iu = imin(m, j + kl);
            const complex_double tmp = alpha0*(*d0)+alpha1*(*d1)+alpha2*(*d2);
            suzerain_blas_zaxpy_d (iu - il, tmp,  a + il*inca, inca,
                                                  b + il*incb, incb);
        }
    } else {
        for (int j = 0; j < n;
             a += lda, d0 += ldd0, d1 += ldd1, d2 += ldd2, b += ldb, ++j) {
            const int il = imax(0, j - ku);
            const int iu = imin(m, j + kl);
            const complex_double tmp = alpha0*(*d0)+alpha1*(*d1)+alpha2*(*d2);
            suzerain_blas_zaxpby_d(iu - il, tmp,  a + il*inca, inca,
                                            beta, b + il*incb, incb);
        }
    }

    return 0;
}

int
suzerain_blasext_zgb_ddddiag_scale_acc_d(
        char side,
        int m,
        int n,
        int kl,
        int ku,
        const complex_double alpha0,
        const double *d0,
        int ldd0,
        const complex_double alpha1,
        const double *d1,
        int ldd1,
        const complex_double alpha2,
        const double *d2,
        int ldd2,
        const complex_double alpha3,
        const double *d3,
        int ldd3,
        const double *a,
        int inca,
        int lda,
        const complex_double beta,
        complex_double *b,
        int incb,
        int ldb)
{
    if (UNLIKELY(kl   < 0       )) return suzerain_blas_xerbla(__func__,  4);
    if (UNLIKELY(ku   < 0       )) return suzerain_blas_xerbla(__func__,  5);
    if (UNLIKELY(ldd0 < 0       )) return suzerain_blas_xerbla(__func__,  8);
    if (UNLIKELY(ldd1 < 0       )) return suzerain_blas_xerbla(__func__, 11);
    if (UNLIKELY(ldd2 < 0       )) return suzerain_blas_xerbla(__func__, 14);
    if (UNLIKELY(ldd3 < 0       )) return suzerain_blas_xerbla(__func__, 17);
    if (UNLIKELY(inca < 1       )) return suzerain_blas_xerbla(__func__, 19);
    if (UNLIKELY(lda  <= kl + ku)) return suzerain_blas_xerbla(__func__, 20);
    if (UNLIKELY(incb < 1       )) return suzerain_blas_xerbla(__func__, 23);
    if (UNLIKELY(lda  <= kl + ku)) return suzerain_blas_xerbla(__func__, 24);

#pragma warning(push,disable:1572)
    const _Bool alpha0_is_zero = (alpha0 == 0.0);
    const _Bool alpha1_is_zero = (alpha1 == 0.0);
    const _Bool alpha2_is_zero = (alpha2 == 0.0);
    const _Bool alpha3_is_zero = (alpha3 == 0.0);
    const _Bool beta_is_one    = (beta   == 1.0);
#pragma warning(pop)

    // If necessary, recast side == 'L' details into a side == 'R' traversal
    switch (toupper(side)) {
    case 'R': break;
    case 'L': inca = lda - inca; a += ku - kl*inca;  // Traverse A by rows
              incb = ldb - incb; b += ku - kl*incb;  // Traverse B by rows
              kl ^= ku; ku ^= kl; kl ^= ku;          // Swap kl and ku
              break;
    default:  return suzerain_blas_xerbla(__func__, 1);
    }

    // Quick return if possible
    if (UNLIKELY(    (alpha0_is_zero && alpha1_is_zero && alpha2_is_zero
                                     && alpha3_is_zero && beta_is_one)
                  || m <= 0 || n <= 0))
        return 0;

    // Banded matrix dereference has form a[(ku + i)*inca + j*(lda - inca)].
    // Incorporate the ku offset and decrement ldX to speed indexing in loops.
    // Further, increment kl anticipating expressions like imin(m, j + kl + 1).
    a += ku*inca; lda -= inca;
    b += ku*incb; ldb -= incb;
    ++kl;

    if (beta_is_one) {
        for (int j = 0; j < n; ++j) {
            const int il = imax(0, j - ku);
            const int iu = imin(m, j + kl);
            const complex_double tmp = alpha0*(*d0)
                                     + alpha1*(*d1)
                                     + alpha2*(*d2)
                                     + alpha3*(*d3);
            suzerain_blas_zaxpy_d (iu - il, tmp,  a + il*inca, inca,
                                                  b + il*incb, incb);
            a += lda; d0 += ldd0; d1 += ldd1; d2 += ldd2; d3 += ldd3; b += ldb;
        }
    } else {
        for (int j = 0; j < n; ++j) {
            const int il = imax(0, j - ku);
            const int iu = imin(m, j + kl);
            const complex_double tmp = alpha0*(*d0)
                                     + alpha1*(*d1)
                                     + alpha2*(*d2)
                                     + alpha3*(*d3);
            suzerain_blas_zaxpby_d(iu - il, tmp,  a + il*inca, inca,
                                            beta, b + il*incb, incb);
            a += lda; d0 += ldd0; d1 += ldd1; d2 += ldd2; d3 += ldd3; b += ldb;
        }
    }

    return 0;
}

int
suzerain_blasext_zgb_dddddiag_scale_acc_d(
        char side,
        int m,
        int n,
        int kl,
        int ku,
        const complex_double alpha0,
        const double *d0,
        int ldd0,
        const complex_double alpha1,
        const double *d1,
        int ldd1,
        const complex_double alpha2,
        const double *d2,
        int ldd2,
        const complex_double alpha3,
        const double *d3,
        int ldd3,
        const complex_double alpha4,
        const double *d4,
        int ldd4,
        const double *a,
        int inca,
        int lda,
        const complex_double beta,
        complex_double *b,
        int incb,
        int ldb)
{
    if (UNLIKELY(kl   < 0       )) return suzerain_blas_xerbla(__func__,  4);
    if (UNLIKELY(ku   < 0       )) return suzerain_blas_xerbla(__func__,  5);
    if (UNLIKELY(ldd0 < 0       )) return suzerain_blas_xerbla(__func__,  8);
    if (UNLIKELY(ldd1 < 0       )) return suzerain_blas_xerbla(__func__, 11);
    if (UNLIKELY(ldd2 < 0       )) return suzerain_blas_xerbla(__func__, 14);
    if (UNLIKELY(ldd3 < 0       )) return suzerain_blas_xerbla(__func__, 17);
    if (UNLIKELY(ldd4 < 0       )) return suzerain_blas_xerbla(__func__, 20);
    if (UNLIKELY(inca < 1       )) return suzerain_blas_xerbla(__func__, 22);
    if (UNLIKELY(lda  <= kl + ku)) return suzerain_blas_xerbla(__func__, 23);
    if (UNLIKELY(incb < 1       )) return suzerain_blas_xerbla(__func__, 26);
    if (UNLIKELY(lda  <= kl + ku)) return suzerain_blas_xerbla(__func__, 27);

#pragma warning(push,disable:1572)
    const _Bool alpha0_is_zero = (alpha0 == 0.0);
    const _Bool alpha1_is_zero = (alpha1 == 0.0);
    const _Bool alpha2_is_zero = (alpha2 == 0.0);
    const _Bool alpha3_is_zero = (alpha3 == 0.0);
    const _Bool alpha4_is_zero = (alpha4 == 0.0);
    const _Bool beta_is_one    = (beta   == 1.0);
#pragma warning(pop)

    // If necessary, recast side == 'L' details into a side == 'R' traversal
    switch (toupper(side)) {
    case 'R': break;
    case 'L': inca = lda - inca; a += ku - kl*inca;  // Traverse A by rows
              incb = ldb - incb; b += ku - kl*incb;  // Traverse B by rows
              kl ^= ku; ku ^= kl; kl ^= ku;          // Swap kl and ku
              break;
    default:  return suzerain_blas_xerbla(__func__, 1);
    }

    // Quick return if possible
    if (UNLIKELY(    (alpha0_is_zero && alpha1_is_zero && alpha2_is_zero
                                     && alpha3_is_zero && alpha4_is_zero
                                     && beta_is_one)
                  || m <= 0 || n <= 0))
        return 0;

    // Banded matrix dereference has form a[(ku + i)*inca + j*(lda - inca)].
    // Incorporate the ku offset and decrement ldX to speed indexing in loops.
    // Further, increment kl anticipating expressions like imin(m, j + kl + 1).
    a += ku*inca; lda -= inca;
    b += ku*incb; ldb -= incb;
    ++kl;

    if (beta_is_one) {
        for (int j = 0; j < n; ++j) {
            const int il = imax(0, j - ku);
            const int iu = imin(m, j + kl);
            const complex_double tmp = alpha0*(*d0)
                                     + alpha1*(*d1)
                                     + alpha2*(*d2)
                                     + alpha3*(*d3)
                                     + alpha4*(*d4);
            suzerain_blas_zaxpy_d (iu - il, tmp,  a + il*inca, inca,
                                                  b + il*incb, incb);
            a += lda;
            d0 += ldd0; d1 += ldd1; d2 += ldd2; d3 += ldd3; d4 += ldd4;
            b += ldb;
        }
    } else {
        for (int j = 0; j < n; ++j) {
            const int il = imax(0, j - ku);
            const int iu = imin(m, j + kl);
            const complex_double tmp = alpha0*(*d0)
                                     + alpha1*(*d1)
                                     + alpha2*(*d2)
                                     + alpha3*(*d3)
                                     + alpha4*(*d4);
            suzerain_blas_zaxpby_d(iu - il, tmp,  a + il*inca, inca,
                                            beta, b + il*incb, incb);
            a += lda;
            d0 += ldd0; d1 += ldd1; d2 += ldd2; d3 += ldd3; d4 += ldd4;
            b += ldb;
        }
    }

    return 0;
}

void
suzerain_blasext_i2s_zaxpby2(
        const int m,
        const int n,
        const double * const alpha,
        const double * const x,
        const int incx,
        const int ldx,
        const double * const beta,
        double * const y_re,
        const int incy_re,
        const int ldy_re,
        double * const y_im,
        const int incy_im,
        const int ldy_im)
{
    assert(incx    >= 0); // TODO Handle negative incx
    assert(incy_re >= 0); // TODO Handle negative incy_re
    assert(incy_im >= 0); // TODO Handle negative incy_im

    /* alpha == NULL defaults to 1 */
    const double alpha_re = (alpha == NULL) ? 1.0 : alpha[0];
    const double alpha_im = (alpha == NULL) ? 0.0 : alpha[1];
    /* beta == NULL defaults to 0 */
    const double beta_re  = (beta == NULL)  ? 0.0 : beta[0];
    const double beta_im  = (beta == NULL)  ? 0.0 : beta[1];

#pragma warning(push,disable:1572)
    if (alpha_im == 0.0 && beta_im == 0.0) {
#pragma warning(pop)

        // Avoid complex arithmetic for real-valued scaling
        for (int j = 0; j < m; ++j) {

            const double * p_x_re = x + 2*j*ldx;
            const double * p_x_im = x + 2*j*ldx + 1;
            double * p_y_re = y_re + j*ldy_re;
            double * p_y_im = y_im + j*ldy_im;

#pragma unroll
            for (int i = 0; i < n; ++i) {

                /* Compute complex-valued alpha*x */
                const double alpha_x_re = (*p_x_re)*alpha_re;
                const double alpha_x_im = (*p_x_im)*alpha_re;

                /* Compute complex-valued beta*y */
                const double beta_y_re =  (*p_y_re)*beta_re;
                const double beta_y_im =  (*p_y_im)*beta_re;

                /* Store result alpha*x + beta*y in y */
                *p_y_re = alpha_x_re + beta_y_re;
                *p_y_im = alpha_x_im + beta_y_im;

                /* Increment pointers for next iteration */
                p_x_re += 2*incx;
                p_x_im += 2*incx;
                p_y_re += incy_re;
                p_y_im += incy_im;
            }
        }

    } else {

        // Complex arithmetic required
        for (int j = 0; j < m; ++j) {

            const double * p_x_re = x + 2*j*ldx;
            const double * p_x_im = x + 2*j*ldx + 1;
            double * p_y_re = y_re + j*ldy_re;
            double * p_y_im = y_im + j*ldy_im;

#pragma unroll
            for (int i = 0; i < n; ++i) {

                /* Compute complex-valued alpha*x */
                const double alpha_x_re
                    = (*p_x_re)*alpha_re - (*p_x_im)*alpha_im;
                const double alpha_x_im
                    = (*p_x_re)*alpha_im + (*p_x_im)*alpha_re;

                /* Compute complex-valued beta*y */
                const double beta_y_re
                    =  (*p_y_re)*beta_re - (*p_y_im)*beta_im;
                const double beta_y_im
                    =  (*p_y_re)*beta_im + (*p_y_im)*beta_re;

                /* Store result alpha*x + beta*y in y */
                *p_y_re = alpha_x_re + beta_y_re;
                *p_y_im = alpha_x_im + beta_y_im;

                /* Increment pointers for next iteration */
                p_x_re += 2*incx;
                p_x_im += 2*incx;
                p_y_re += incy_re;
                p_y_im += incy_im;
            }
        }

    }
}

int
suzerain_blasext_sgbnorm1(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const float *a,
        const int lda,
        float *norm1)
{
    // Defer to LAPACK for the square case
    if (m == n) {
        *norm1 = suzerain_lapack_slangb('1', n, kl, ku, a, lda, NULL);
        return 0;
    }

    if (UNLIKELY(m  < 0))        return suzerain_blas_xerbla(__func__, -1);
    if (UNLIKELY(n  < 0))        return suzerain_blas_xerbla(__func__, -2);
    if (UNLIKELY(kl < 0))        return suzerain_blas_xerbla(__func__, -3);
    if (UNLIKELY(ku < 0))        return suzerain_blas_xerbla(__func__, -4);
    if (UNLIKELY(!a))            return suzerain_blas_xerbla(__func__, -5);
    if (UNLIKELY(lda < kl+ku+1)) return suzerain_blas_xerbla(__func__, -6);
    if (UNLIKELY(!norm1))        return suzerain_blas_xerbla(__func__, -7);

    *norm1 = 0;
    for (int j = 0; j < n; ++j) {
        const float *a_j  = a + j *lda;

        float     s    = 0;
        const int ibgn = imax(0, ku - j);
        const int iend = ku + imin(kl + 1, m - j);
        for (int i = ibgn; i < iend; ++i) {
            s += fabsf(a_j[i]);
        }
        *norm1 = fmaxf(s, *norm1);
    }

    return 0;
}

int
suzerain_blasext_dgbnorm1(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const double *a,
        const int lda,
        double *norm1)
{
    // Defer to LAPACK for the square case
    if (m == n) {
        *norm1 = suzerain_lapack_dlangb('1', n, kl, ku, a, lda, NULL);
        return 0;
    }

    if (UNLIKELY(m  < 0))        return suzerain_blas_xerbla(__func__, -1);
    if (UNLIKELY(n  < 0))        return suzerain_blas_xerbla(__func__, -2);
    if (UNLIKELY(kl < 0))        return suzerain_blas_xerbla(__func__, -3);
    if (UNLIKELY(ku < 0))        return suzerain_blas_xerbla(__func__, -4);
    if (UNLIKELY(!a))            return suzerain_blas_xerbla(__func__, -5);
    if (UNLIKELY(lda < kl+ku+1)) return suzerain_blas_xerbla(__func__, -6);
    if (UNLIKELY(!norm1))        return suzerain_blas_xerbla(__func__, -7);

    *norm1 = 0;
    for (int j = 0; j < n; ++j) {
        const double *a_j  = a + j *lda;

        double    s    = 0;
        const int ibgn = imax(0, ku - j);
        const int iend = ku + imin(kl + 1, m - j);
        for (int i = ibgn; i < iend; ++i) {
            s += fabs(a_j[i]);
        }
        *norm1 = fmax(s, *norm1);
    }

    return 0;
}

int
suzerain_blasext_cgbnorm1(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const complex_float *a,
        const int lda,
        float *norm1)
{
    // Defer to LAPACK for the square case
    if (m == n) {
        *norm1 = suzerain_lapack_clangb('1', n, kl, ku, a, lda, NULL);
        return 0;
    }

    if (UNLIKELY(m  < 0))        return suzerain_blas_xerbla(__func__, -1);
    if (UNLIKELY(n  < 0))        return suzerain_blas_xerbla(__func__, -2);
    if (UNLIKELY(kl < 0))        return suzerain_blas_xerbla(__func__, -3);
    if (UNLIKELY(ku < 0))        return suzerain_blas_xerbla(__func__, -4);
    if (UNLIKELY(!a))            return suzerain_blas_xerbla(__func__, -5);
    if (UNLIKELY(lda < kl+ku+1)) return suzerain_blas_xerbla(__func__, -6);
    if (UNLIKELY(!norm1))        return suzerain_blas_xerbla(__func__, -7);

    *norm1 = 0;
    for (int j = 0; j < n; ++j) {
        const complex_float *a_j  = a + j *lda;

        float     s    = 0;
        const int ibgn = imax(0, ku - j);
        const int iend = ku + imin(kl + 1, m - j);
        for (int i = ibgn; i < iend; ++i) {
            s += cabsf(a_j[i]);
        }
        *norm1 = fmaxf(s, *norm1);
    }

    return 0;
}

int
suzerain_blasext_zgbnorm1(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const complex_double *a,
        const int lda,
        double *norm1)
{
    // Defer to LAPACK for the square case
    if (m == n) {
        *norm1 = suzerain_lapack_zlangb('1', n, kl, ku, a, lda, NULL);
        return 0;
    }

    if (UNLIKELY(m  < 0))        return suzerain_blas_xerbla(__func__, -1);
    if (UNLIKELY(n  < 0))        return suzerain_blas_xerbla(__func__, -2);
    if (UNLIKELY(kl < 0))        return suzerain_blas_xerbla(__func__, -3);
    if (UNLIKELY(ku < 0))        return suzerain_blas_xerbla(__func__, -4);
    if (UNLIKELY(!a))            return suzerain_blas_xerbla(__func__, -5);
    if (UNLIKELY(lda < kl+ku+1)) return suzerain_blas_xerbla(__func__, -6);
    if (UNLIKELY(!norm1))        return suzerain_blas_xerbla(__func__, -7);

    *norm1 = 0;
    for (int j = 0; j < n; ++j) {
        const complex_double *a_j  = a + j *lda;

        double    s    = 0;
        const int ibgn = imax(0, ku - j);
        const int iend = ku + imin(kl + 1, m - j);
        for (int i = ibgn; i < iend; ++i) {
            s += cabs(a_j[i]);
        }
        *norm1 = fmax(s, *norm1);
    }

    return 0;
}

inline int
suzerain_blasext_ddemote(
        const int n,
        void *x)
{
    if (UNLIKELY(n < 0)) return -1;
    if (UNLIKELY(!x))    return -2;

    // Written to be -ansi-alias, -Wstrict-aliasing friendly
    char * const xchar = (char *) x;
    for (int i = 0; i < n; ++i) {
        double in;
        memcpy(&in, xchar + i*sizeof(in), sizeof(in));
#pragma warning(push,disable:2259)
        float out = in;
#pragma warning(pop)
        memcpy(xchar + i*sizeof(out), &out, sizeof(out));
    }

    return 0;
}

inline int
suzerain_blasext_dpromote(
        const int n,
        void *x)
{
    if (UNLIKELY(n < 0)) return -1;
    if (UNLIKELY(!x))    return -2;

    // Written to be -ansi-alias, -Wstrict-aliasing friendly
    char * const xchar = (char *) x;
    for (int i = n; i --> 0;) {
        float in;
        memcpy(&in, xchar + i*sizeof(in), sizeof(in));
        double out = in;
        memcpy(xchar + i*sizeof(out), &out, sizeof(out));
    }

    return 0;
}

inline int
suzerain_blasext_zdemote(
        const int n,
        void *x)
{
    if (UNLIKELY(n < 0)) return -1;
    if (UNLIKELY(!x))    return -2;

    // Written to be -ansi-alias, -Wstrict-aliasing friendly
    char * const xchar = (char *) x;
    for (int i = 0; i < n; ++i) {
        complex_double in;
        memcpy(&in, xchar + i*sizeof(in), sizeof(in));
#pragma warning(push,disable:2259)
        complex_float out = in;
#pragma warning(pop)
        memcpy(xchar + i*sizeof(out), &out, sizeof(out));
    }

    return 0;
}

inline int
suzerain_blasext_zpromote(
        const int n,
        void *x)
{
    if (UNLIKELY(n < 0)) return -1;
    if (UNLIKELY(!x))    return -2;

    // Written to be -ansi-alias, -Wstrict-aliasing friendly
    char * const xchar = (char *) x;
    for (int i = n; i --> 0;) {
        complex_float in;
        memcpy(&in, xchar + i*sizeof(in), sizeof(in));
        complex_double out = in;
        memcpy(xchar + i*sizeof(out), &out, sizeof(out));
    }

    return 0;
}
