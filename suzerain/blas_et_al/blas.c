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
 * @copydoc blas.h
 */

#include <suzerain/blas_et_al/blas.h>

#include <suzerain/common.h>
#include <suzerain/blas_et_al/blasext.h>
#include <suzerain/blas_et_al/gbmv.h>
#include <suzerain/blas_et_al/sbmv.h>

// Odd looking, unused, anonymous enumerations are globally-scoped versions
// of "Compile Time Assertions" by Ralf Holly (http://drdobbs.com/184401873)
enum { // Required for strict aliasing workarounds
    assert_floatp  = 1/(sizeof(float*)  == sizeof(complex_float*) ),
    assert_doublep = 1/(sizeof(double*) == sizeof(complex_double*))
};

#ifdef SUZERAIN_HAVE_MKL
# include <mkl_types.h>
# include <mkl_blas.h>
# define BLAS_FUNC(name,NAME)   name             /* Think AC_F77_WRAPPERS */
enum { // Required for binary interoperability with MKL
    assert_mkl_i = 1 /(sizeof(MKL_INT)       == sizeof(int)           ),
    assert_mkl_c = 1 /(sizeof(MKL_Complex8)  == sizeof(complex_float) ),
    assert_mkl_z = 1 /(sizeof(MKL_Complex16) == sizeof(complex_double))
};
#else
# error "No suitable BLAS and/or LAPACK library found during configuration"
#endif

// Shorthand
#define CAST_VOIDP(l)       SUZERAIN_CAST_VOIDP(l)
#define CONST_CAST_VOIDP(l) SUZERAIN_CONST_CAST_VOIDP(l)
#define UNLIKELY(expr)      SUZERAIN_UNLIKELY(expr)

// Many of the short methods have "inline" though their declarations do not.
// to allow inlining them later within this particular translation unit.

// Thank you captain obvious...
#pragma warning(disable:981)

int
suzerain_blas_xerbla(const char *srname, const int info)
{
    const int lsrname = srname ? strlen(srname) : 0;
    BLAS_FUNC(xerbla,XERBLA)(srname, &info, lsrname);
    return info;
}

void *
suzerain_blas_malloc(size_t size)
{
    void * p = NULL;

    /* We do not use MKL_malloc to avoid later needing MKL_free calls. */
    const int status = posix_memalign(&p, SUZERAIN_BLAS_ALIGNMENT, size);

    switch (status) {
    case 0:
#ifndef NDEBUG  /* Debug?  Fill likely double-valued buffers with NaNs. */
    {
        const double NaN = nan("NAN");
        const size_t quo = size / sizeof(double);
        for (size_t i = 0; i < quo; ++i) ((double *)p)[i] = NaN;
    }
#endif
        return p;
    case ENOMEM:
        return NULL;
    case EINVAL:
    default:
        fprintf(stderr, "Fatal posix_memalign error at %s:%d-- %s\n",
                __FILE__, __LINE__, strerror(status));
        abort();
    }
}

void *
suzerain_blas_calloc(size_t nmemb, size_t size)
{
    const size_t total_bytes = nmemb * size;
    void * p = suzerain_blas_malloc(total_bytes);
    if (p != NULL) {
        memset(p, 0, total_bytes);
    }
    return p;
}

void
suzerain_blas_free(void *ptr)
{
    // Must match suzerain_blas_malloc's malloc-like routine!
    if (ptr) free(ptr);
}

inline void
suzerain_blas_sscal(
        const int n,
        const float alpha,
        float *x,
        const int incx)
{
    return BLAS_FUNC(sscal,SSCAL)(&n, &alpha, x, &incx);
}

inline void
suzerain_blas_dscal(
        const int n,
        const double alpha,
        double *x,
        const int incx)
{
    return BLAS_FUNC(dscal,DSCAL)(&n, &alpha, x, &incx);
}

inline void
suzerain_blas_cscal(
        const int n,
        const complex_float alpha,
        complex_float *x,
        const int incx)
{
    return BLAS_FUNC(cscal,CSCAL)(
            &n, CONST_CAST_VOIDP(&alpha), CAST_VOIDP(x), &incx);
}

inline void
suzerain_blas_zscal(
        const int n,
        const complex_double alpha,
        complex_double *x,
        const int incx)
{
    return BLAS_FUNC(zscal,ZSCAL)(
            &n, CONST_CAST_VOIDP(&alpha), CAST_VOIDP(x), &incx);
}

inline void
suzerain_blas_sswap(
        const int n,
        float *x,
        const int incx,
        float *y,
        const int incy)
{
    return BLAS_FUNC(sswap,SSWAP)(&n, x, &incx, y, &incy);
}

inline void
suzerain_blas_dswap(
        const int n,
        double *x,
        const int incx,
        double *y,
        const int incy)
{
    return BLAS_FUNC(dswap,DSWAP)(&n, x, &incx, y, &incy);
}

inline void
suzerain_blas_cswap(
        const int n,
        complex_float *x,
        const int incx,
        complex_float *y,
        const int incy)
{
    return BLAS_FUNC(cswap,CSWAP)(
                &n, CAST_VOIDP(x), &incx, CAST_VOIDP(y), &incy);
}

inline void
suzerain_blas_zswap(
        const int n,
        complex_double *x,
        const int incx,
        complex_double *y,
        const int incy)
{
    return BLAS_FUNC(zswap,ZSWAP)(
                &n, CAST_VOIDP(x), &incx, CAST_VOIDP(y), &incy);
}

inline void
suzerain_blas_scopy(
        const int n,
        const float *x,
        const int incx,
        float *y,
        const int incy)
{
    return BLAS_FUNC(scopy,SCOPY)(&n, x, &incx, y, &incy);
}

inline void
suzerain_blas_dcopy(
        const int n,
        const double *x,
        const int incx,
        double *y,
        const int incy)
{
    return BLAS_FUNC(dcopy,DCOPY)(&n, x, &incx, y, &incy);
}

inline void
suzerain_blas_ccopy(
        const int n,
        const complex_float *x,
        const int incx,
        complex_float *y,
        const int incy)
{
    return BLAS_FUNC(ccopy,CCOPY)(
            &n, CONST_CAST_VOIDP(x), &incx, CAST_VOIDP(y), &incy);
}

inline void
suzerain_blas_zcopy(
        const int n,
        const complex_double *x,
        const int incx,
        complex_double *y,
        const int incy)
{
    return BLAS_FUNC(zcopy,ZCOPY)(
            &n, CONST_CAST_VOIDP(x), &incx, CAST_VOIDP(y), &incy);
}

inline float
suzerain_blas_sdot(
        const int n,
        const float *x,
        const int incx,
        const float *y,
        const int incy)
{
    return BLAS_FUNC(sdot,SDOT)(&n, x, &incx, y, &incy);
}

inline double
suzerain_blas_ddot(
        const int n,
        const double *x,
        const int incx,
        const double *y,
        const int incy)
{
    return BLAS_FUNC(ddot,DDOT)(&n, x, &incx, y, &incy);
}

inline void
suzerain_blas_cdotc(
        const int n,
        const complex_float *x,
        const int incx,
        const complex_float *y,
        const int incy,
        complex_float *dotc)
{
    return BLAS_FUNC(cdotc,CDOTC)(      CAST_VOIDP(dotc), &n,
                                  CONST_CAST_VOIDP(x   ), &incx,
                                  CONST_CAST_VOIDP(y   ), &incy);
}

inline void
suzerain_blas_zdotc(
        const int n,
        const complex_double *x,
        const int incx,
        const complex_double *y,
        const int incy,
        complex_double *dotc)
{
    return BLAS_FUNC(zdotc,ZDOTC)(      CAST_VOIDP(dotc), &n,
                                  CONST_CAST_VOIDP(x   ), &incx,
                                  CONST_CAST_VOIDP(y   ), &incy);
}

inline float
suzerain_blas_snrm2(
        const int n,
        const float *x,
        const int incx)
{
    return BLAS_FUNC(snrm2,SNRM2)(&n, x, &incx);
}

inline double
suzerain_blas_dnrm2(
        const int n,
        const double *x,
        const int incx)
{
    return BLAS_FUNC(dnrm2,DNRM2)(&n, x, &incx);
}

inline float
suzerain_blas_scnrm2(
        const int n,
        const complex_float *x,
        const int incx)
{
    return BLAS_FUNC(scnrm2,SCNRM2)(&n, CONST_CAST_VOIDP(x), &incx);
}

inline double
suzerain_blas_dznrm2(
        const int n,
        const complex_double *x,
        const int incx)
{
    return BLAS_FUNC(dznrm2,DZNRM2)(&n, CONST_CAST_VOIDP(x), &incx);
}

inline float
suzerain_blas_sasum(
        const int n,
        const float *x,
        const int incx)
{
    return BLAS_FUNC(sasum,SASUM)(&n, x, &incx);
}

inline double
suzerain_blas_dasum(
        const int n,
        const double *x,
        const int incx)
{
    return BLAS_FUNC(dasum,DASUM)(&n, x, &incx);
}

inline float
suzerain_blas_scasum(
        const int n,
        const complex_float *x,
        const int incx)
{
    return BLAS_FUNC(scasum,SCASUM)(&n, CAST_VOIDP(x), &incx);
}

inline double
suzerain_blas_dzasum(
        const int n,
        const complex_double *x,
        const int incx)
{
    return BLAS_FUNC(dzasum,DZASUM)(&n, CAST_VOIDP(x), &incx);
}

inline int
suzerain_blas_isamax(
        const int n,
        const float *x,
        const int incx)
{
#ifdef SUZERAIN_HAVE_MKL
    return isamax(&n, x, &incx) - 1 /* zero-indexed */;
#else
#error "Sanity: suzerain_blas_isamax not implemented for BLAS"
#endif
}

inline int
suzerain_blas_idamax(
        const int n,
        const double *x,
        const int incx)
{
#ifdef SUZERAIN_HAVE_MKL
    return idamax(&n, x, &incx) - 1 /* zero-indexed */;
#else
#error "Sanity: suzerain_blas_idamax not implemented for BLAS"
#endif
}

inline int
suzerain_blas_icamax(
        const int n,
        const complex_float *x,
        const int incx)
{
#ifdef SUZERAIN_HAVE_MKL
    return icamax(&n, CONST_CAST_VOIDP(x), &incx) - 1 /* zero-indexed */;
#else
#error "Sanity: suzerain_blas_icamax not implemented for BLAS"
#endif
}

inline int
suzerain_blas_izamax(
        const int n,
        const complex_double *x,
        const int incx)
{
#ifdef SUZERAIN_HAVE_MKL
    return izamax(&n, CONST_CAST_VOIDP(x), &incx) - 1 /* zero-indexed */;
#else
#error "Sanity: suzerain_blas_izamax not implemented for BLAS"
#endif
}

inline int
suzerain_blas_isamin(
        const int n,
        const float *x,
        const int incx)
{
#ifdef SUZERAIN_HAVE_MKL
    return isamin(&n, x, &incx) - 1 /* zero-indexed */;
#else
#error "Sanity: suzerain_blas_isamin not implemented for BLAS"
#endif
}

inline int
suzerain_blas_idamin(
        const int n,
        const double *x,
        const int incx)
{
#ifdef SUZERAIN_HAVE_MKL
    return idamin(&n, x, &incx) - 1 /* zero-indexed */;
#else
#error "Sanity: suzerain_blas_idamin not implemented for BLAS"
#endif
}

inline int
suzerain_blas_icamin(
        const int n,
        const complex_float *x,
        const int incx)
{
#ifdef SUZERAIN_HAVE_MKL
    return icamin(&n, CONST_CAST_VOIDP(x), &incx) - 1 /* zero-indexed */;
#else
#error "Sanity: suzerain_blas_icamin not implemented for BLAS"
#endif
}

inline int
suzerain_blas_izamin(
        const int n,
        const complex_double *x,
        const int incx)
{
#ifdef SUZERAIN_HAVE_MKL
    return izamin(&n, CONST_CAST_VOIDP(x), &incx) - 1 /* zero-indexed */;
#else
#error "Sanity: suzerain_blas_izamin not implemented for BLAS"
#endif
}

inline void
suzerain_blas_saxpy(
        const int n,
        const float alpha,
        const float *x,
        const int incx,
        float *y,
        const int incy)
{
    return BLAS_FUNC(saxpy,SAXPY)(&n, &alpha, x, &incx, y, &incy);
}

inline void
suzerain_blas_daxpy(
        const int n,
        const double alpha,
        const double *x,
        const int incx,
        double *y,
        const int incy)
{
    return BLAS_FUNC(daxpy,DAXPY)(&n, &alpha, x, &incx, y, &incy);
}

inline void
suzerain_blas_caxpy(
        const int n,
        const complex_float alpha,
        const complex_float *x,
        const int incx,
        complex_float *y,
        const int incy)
{
    return BLAS_FUNC(caxpy,CAXPY)(&n, CONST_CAST_VOIDP(&alpha),
                                      CONST_CAST_VOIDP(x     ), &incx,
                                            CAST_VOIDP(y     ), &incy);
}

inline void
suzerain_blas_zaxpy(
        const int n,
        const complex_double alpha,
        const complex_double *x,
        const int incx,
        complex_double *y,
        const int incy)
{
    return BLAS_FUNC(zaxpy,ZAXPY)(&n, CONST_CAST_VOIDP(&alpha),
                                      CONST_CAST_VOIDP(x     ), &incx,
                                            CAST_VOIDP(y     ), &incy);
}

inline void
suzerain_blas_zaxpy_d(
        const int n,
        const complex_double alpha,
        const double *x,
        const int incx,
        complex_double *y,
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    if (incx == 1 && incy == 1) {  // Unit strides

#pragma unroll
        for (int i = 0; i < n; ++i) {
            y[i] += alpha*x[i];
        }

    } else {                       // General strides

        // Adjust for possibly negative incx and incy
        int ix = (incx < 0) ? (1 - n)*incx : 0;
        int iy = (incy < 0) ? (1 - n)*incy : 0;
#pragma unroll
        for (int i = 0; i < n; ++i, ix += incx, iy += incy) {
            y[iy] += alpha*x[ix];
        }

    }
#else
#error "Sanity: suzerain_blas_zaxpy_d not implemented for BLAS"
#endif
}

inline void
suzerain_blas_saxpby(
        const int n,
        const float alpha,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    /* Simulate saxpby since MKL lacks the routine. */
#pragma warning(push,disable:1572)
    if (beta != 1.0f)
#pragma warning(pop)
        sscal(&n, &beta, y, &incy);
    saxpy(&n, &alpha, x, &incx, y, &incy);
#else
    BLAS_FUNC(saxpby,SAXPBY)(&n, &alpha, x, &incx, &beta, y, &incy);
#endif
}

inline void
suzerain_blas_daxpby(
        const int n,
        const double alpha,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    /* Simulate daxpby since MKL lacks the routine. */
#pragma warning(push,disable:1572)
    if (beta != 1.0)
#pragma warning(pop)
        dscal(&n, &beta, y, &incy);
    daxpy(&n, &alpha, x, &incx, y, &incy);
#else
    BLAS_FUNC(daxpby,DAXPBY)(&n, &alpha, x, &incx, &beta, y, &incy);
#endif
}

inline void
suzerain_blas_caxpby(
        const int n,
        const complex_float alpha,
        const complex_float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    /* Simulate caxpby since MKL lacks the routine. */
#pragma warning(push,disable:1572)
    if (beta != 1.0f)
#pragma warning(pop)
        suzerain_blas_cscal(n, beta, y, incy);
    suzerain_blas_caxpy(n, alpha, x, incx, y, incy);
#else
    BLAS_FUNC(caxpby,CAXPBY)(
            &n,
            CONST_CAST_VOIDP(&alpha), CONST_CAST_VOIDP(x), &incx,
            CONST_CAST_VOIDP(&beta ),       CAST_VOIDP(y), &incy);
#endif
}

inline void
suzerain_blas_zaxpby(
        const int n,
        const complex_double alpha,
        const complex_double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    /* Simulate caxpby since MKL lacks the routine. */
#pragma warning(push,disable:1572)
    if (beta != 1.0)
#pragma warning(pop)
        suzerain_blas_zscal(n, beta, y, incy);
    suzerain_blas_zaxpy(n, alpha, x, incx, y, incy);
#else
    BLAS_FUNC(zaxpby,ZAXPBY)(
            &n,
            CONST_CAST_VOIDP(&alpha), CONST_CAST_VOIDP(x), &incx,
            CONST_CAST_VOIDP(&beta ),       CAST_VOIDP(y), &incy);
#endif
}

inline void
suzerain_blas_zaxpby_d(
        const int n,
        const complex_double alpha,
        const double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL

#pragma warning(push,disable:1572)
    const _Bool beta_is_zero = (beta == 0);
#pragma warning(pop)

    if (incx == 1 && incy == 1) {  // Unit strides

        if (UNLIKELY(beta_is_zero)) {
#pragma unroll
            for (int i = 0; i < n; ++i) {
                y[i]  = alpha*x[i];
            }
        } else {
#pragma unroll
            for (int i = 0; i < n; ++i) {
                y[i] *= beta;
                y[i] += alpha*x[i];
            }
        }

    } else {                       // General strides

        // Adjust for possibly negative incx and incy
        int ix = (incx < 0) ? (1 - n)*incx : 0;
        int iy = (incy < 0) ? (1 - n)*incy : 0;
        if (UNLIKELY(beta_is_zero)) {
#pragma unroll
            for (int i = 0; i < n; ++i, ix += incx, iy += incy) {
                y[iy]  = alpha*x[ix];
            }
        } else {
#pragma unroll
            for (int i = 0; i < n; ++i, ix += incx, iy += incy) {
                y[iy] *= beta;
                y[iy] += alpha*x[ix];
            }
        }

    }
#else
#error "Sanity: suzerain_blas_zaxpby_d not implemented for BLAS"
#endif
}

inline void
suzerain_blas_swaxpby(
        const int n,
        const float alpha,
        const float *x,
        const int incx,
        const float beta,
        const float *y,
        const int incy,
        float *w,
        const int incw)
{
#ifdef SUZERAIN_HAVE_MKL
    /* Simulate swaxpby since MKL lacks the routine. */
    scopy(&n, y, &incy, w, &incw);
#pragma warning(push,disable:1572)
    if (beta != 1.0f)
#pragma warning(pop)
        sscal(&n, &beta, w, &incw);
    saxpy(&n, &alpha, x, &incx, w, &incw);
#else
    BLAS_FUNC(swaxpby,SWAXPBY)(&n,
                               &alpha, x, &incx,
                               &beta,  y, &incy,
                                       w, &incw);
#endif
}

inline void
suzerain_blas_dwaxpby(
        const int n,
        const double alpha,
        const double *x,
        const int incx,
        const double beta,
        const double *y,
        const int incy,
        double *w,
        const int incw)
{
#ifdef SUZERAIN_HAVE_MKL
    /* Simulate dwaxpby since MKL lacks the routine. */
    dcopy(&n, y, &incy, w, &incw);
#pragma warning(push,disable:1572)
    if (beta != 1.0)
#pragma warning(pop)
        dscal(&n, &beta, w, &incw);
    daxpy(&n, &alpha, x, &incx, w, &incw);
#else
    BLAS_FUNC(dwaxpby,DWAXPBY)(&n,
                               &alpha, x, &incx,
                               &beta,  y, &incy,
                                       w, &incw);
#endif
}

inline void
suzerain_blas_sgbmv_external(
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
    return BLAS_FUNC(sgbmv,SGBMV)(&trans, &m, &n, &kl, &ku,
                                  &alpha, a, &lda, x, &incx,
                                  &beta,           y, &incy);
}

inline void
suzerain_blas_dgbmv_external(
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
    return BLAS_FUNC(dgbmv,DGBMV)(&trans, &m, &n, &kl, &ku,
                                  &alpha, a, &lda, x, &incx,
                                  &beta,           y, &incy);
}

inline void
suzerain_blas_cgbmv_external(
        const char trans,
        const int m,
        const int n,
        const int kl,
        const int ku,
        const complex_float alpha,
        const complex_float *a,
        const int lda,
        const complex_float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy)
{
    return BLAS_FUNC(cgbmv,CGBMV)(
            &trans, &m, &n, &kl, &ku,
            CAST_VOIDP(&alpha), CAST_VOIDP(a), &lda, CAST_VOIDP(x), &incx,
            CAST_VOIDP(&beta ),                      CAST_VOIDP(y), &incy);
}

inline void
suzerain_blas_zgbmv_external(
        const char trans,
        const int m,
        const int n,
        const int kl,
        const int ku,
        const complex_double alpha,
        const complex_double *a,
        const int lda,
        const complex_double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy)
{
    return BLAS_FUNC(zgbmv,ZGBMV)(
            &trans, &m, &n, &kl, &ku,
            CAST_VOIDP(&alpha), CAST_VOIDP(a), &lda, CAST_VOIDP(x), &incx,
            CAST_VOIDP(&beta ),                      CAST_VOIDP(y), &incy);
}

void
suzerain_blas_cgbmv_s_c_external(
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
    float *x_re, *y_re;
    memcpy(&x_re, &x, sizeof(x));
    memcpy(&y_re, &y, sizeof(y));

#pragma warning(push,disable:1572)
    if (cimag(alpha) == 0.0 && cimag(beta) == 0.0) {
#pragma warning(pop)
        /* Real-valued alpha and beta: scale y as we go */
        suzerain_blas_sgbmv_external(trans, m, n, kl, ku,
                                     creal(alpha), a, lda, x_re, 2*incx,
                                     creal(beta), y_re, 2*incy);
        suzerain_blas_sgbmv_external(trans, m, n, kl, ku,
                                     -cimag(alpha), a, lda, x_re+1, 2*incx,
                                     1.0, y_re, 2*incy);
        suzerain_blas_sgbmv_external(trans, m, n, kl, ku,
                                     creal(alpha), a, lda, x_re+1, 2*incx,
                                     creal(beta), y_re+1, 2*incy);
        suzerain_blas_sgbmv_external(trans, m, n, kl, ku,
                                     cimag(alpha), a, lda, x_re, 2*incx,
                                     1.0, y_re+1, 2*incy);
    } else {
        /* Complex-valued alpha and/or beta: scale y and then accumulate */
        suzerain_blas_cscal((toupper(trans) == 'N' ? m : n), beta, y, incy);
        suzerain_blas_sgbmv_external(trans, m, n, kl, ku,
                                     creal(alpha), a, lda, x_re, 2*incx,
                                     1.0, y_re, 2*incy);
        suzerain_blas_sgbmv_external(trans, m, n, kl, ku,
                                     cimag(alpha), a, lda, x_re, 2*incx,
                                     1.0, y_re+1, 2*incy);
        suzerain_blas_sgbmv_external(trans, m, n, kl, ku,
                                     creal(alpha), a, lda, x_re+1, 2*incx,
                                     1.0, y_re+1, 2*incy);
        suzerain_blas_sgbmv_external(trans, m, n, kl, ku,
                                     -cimag(alpha), a, lda, x_re+1, 2*incx,
                                     1.0, y_re, 2*incy);
    }
}

void
suzerain_blas_zgbmv_d_z_external(
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
    double *x_re, *y_re;
    memcpy(&x_re, &x, sizeof(x));
    memcpy(&y_re, &y, sizeof(y));

#pragma warning(push,disable:1572)
    if (cimag(alpha) == 0.0 && cimag(beta) == 0.0) {
#pragma warning(pop)
        /* Real-valued alpha and beta: scale y as we go */
        suzerain_blas_dgbmv_external(trans, m, n, kl, ku,
                                     creal(alpha), a, lda, x_re, 2*incx,
                                     creal(beta), y_re, 2*incy);
        suzerain_blas_dgbmv_external(trans, m, n, kl, ku,
                                     -cimag(alpha), a, lda, x_re+1, 2*incx,
                                     1.0, y_re, 2*incy);
        suzerain_blas_dgbmv_external(trans, m, n, kl, ku,
                                     creal(alpha), a, lda, x_re+1, 2*incx,
                                     creal(beta), y_re+1, 2*incy);
        suzerain_blas_dgbmv_external(trans, m, n, kl, ku,
                                     cimag(alpha), a, lda, x_re, 2*incx,
                                     1.0, y_re+1, 2*incy);
    } else {
        /* Complex-valued alpha and/or beta: scale y and then accumulate */
        suzerain_blas_zscal((toupper(trans) == 'N' ? m : n), beta, y, incy);
        suzerain_blas_dgbmv_external(trans, m, n, kl, ku,
                                     creal(alpha), a, lda, x_re, 2*incx,
                                     1.0, y_re, 2*incy);
        suzerain_blas_dgbmv_external(trans, m, n, kl, ku,
                                     cimag(alpha), a, lda, x_re, 2*incx,
                                     1.0, y_re+1, 2*incy);
        suzerain_blas_dgbmv_external(trans, m, n, kl, ku,
                                     creal(alpha), a, lda, x_re+1, 2*incx,
                                     1.0, y_re+1, 2*incy);
        suzerain_blas_dgbmv_external(trans, m, n, kl, ku,
                                     -cimag(alpha), a, lda, x_re+1, 2*incx,
                                     1.0, y_re, 2*incy);
    }
}

inline int
suzerain_blas_sgbmv(
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
    const int info = suzerain_gbmv_s(trans, m, n, kl, ku,
                                     alpha, a, lda, x, incx,
                                     beta,          y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

inline int
suzerain_blas_dgbmv(
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
    const int info = suzerain_gbmv_d(trans, m, n, kl, ku,
                                     alpha, a, lda, x, incx,
                                     beta,          y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

inline int
suzerain_blas_cgbmv(
        const char trans,
        const int m,
        const int n,
        const int kl,
        const int ku,
        const complex_float alpha,
        const complex_float *a,
        const int lda,
        const complex_float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy)
{
    const int info = suzerain_gbmv_c(trans, m, n, kl, ku,
                                     alpha, a, lda, x, incx,
                                     beta,          y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

inline int
suzerain_blas_zgbmv(
        const char trans,
        const int m,
        const int n,
        const int kl,
        const int ku,
        const complex_double alpha,
        const complex_double *a,
        const int lda,
        const complex_double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy)
{
    const int info = suzerain_gbmv_z(trans, m, n, kl, ku,
                                     alpha, a, lda, x, incx,
                                     beta,          y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

inline int
suzerain_blas_cgbmv_s_c(
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
    const int info = suzerain_gbmv_scc(trans, m, n, kl, ku,
                                       alpha, a, lda, x, incx,
                                       beta,          y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

inline int
suzerain_blas_zgbmv_d_z(
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
    const int info = suzerain_gbmv_dzz(trans, m, n, kl, ku,
                                       alpha, a, lda, x, incx,
                                       beta,          y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

inline int
suzerain_blas_cgbmv_s_s(
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
    const int info = suzerain_gbmv_ssc(trans, m, n, kl, ku,
                                       alpha, a, lda, x, incx,
                                       beta,          y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

inline int
suzerain_blas_zgbmv_d_d(
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
    const int info = suzerain_gbmv_ddz(trans, m, n, kl, ku,
                                       alpha, a, lda, x, incx,
                                       beta,          y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

inline void
suzerain_blas_ssbmv_external(
        const char uplo,
        const int n,
        const int k,
        const float alpha,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    return BLAS_FUNC(ssbmv,SSBMV)(&uplo, &n, &k, &alpha, a, &lda,
                                  x, &incx, &beta, y, &incy);
}

inline void
suzerain_blas_dsbmv_external(
        const char uplo,
        const int n,
        const int k,
        const double alpha,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    return BLAS_FUNC(dsbmv,DSBMV)(&uplo, &n, &k, &alpha, a, &lda,
                                  x, &incx, &beta, y, &incy);
}

inline int
suzerain_blas_ssbmv(
        const char uplo,
        const int n,
        const int k,
        const float alpha,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    const int info = suzerain_sbmv_s(uplo, n, k,
                                     alpha, a, lda, x, incx,
                                     beta,          y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

inline int
suzerain_blas_dsbmv(
        const char uplo,
        const int n,
        const int k,
        const double alpha,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    const int info = suzerain_sbmv_d(uplo, n, k,
                                     alpha, a, lda, x, incx,
                                     beta,          y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

inline int
suzerain_blas_csbmv_s_c(
        const char uplo,
        const int n,
        const int k,
        const complex_float alpha,
        const float *a,
        const int lda,
        const complex_float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy)
{
    const int info = suzerain_sbmv_scc(uplo, n, k,
                                      alpha, a, lda, x, incx,
                                      beta,          y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

inline int
suzerain_blas_zsbmv_d_z(
        const char uplo,
        const int n,
        const int k,
        const complex_double alpha,
        const double *a,
        const int lda,
        const complex_double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy)
{
    const int info = suzerain_sbmv_dzz(uplo, n, k,
                                       alpha, a, lda, x, incx,
                                       beta,          y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

inline int
suzerain_blas_csbmv_s_s(
        const char uplo,
        const int n,
        const int k,
        const complex_float alpha,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy)
{
    const int info = suzerain_sbmv_ssc(uplo, n, k,
                                      alpha, a, lda, x, incx,
                                      beta,          y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

inline int
suzerain_blas_zsbmv_d_d(
        const char uplo,
        const int n,
        const int k,
        const complex_double alpha,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy)
{
    const int info = suzerain_sbmv_ddz(uplo, n, k,
                                       alpha, a, lda, x, incx,
                                       beta,          y, incy);
    if (UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
    return info;
}

inline int
suzerain_blas_sgb_acc(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const float alpha,
        const float *a,
        const int lda,
        const float beta,
        float *b,
        const int ldb)
{
#ifdef SUZERAIN_HAVE_MKL
    /* Internal sgb_acc because MKL lacks the routine. */
    const float one = 1.0f;
    return suzerain_blasext_sgb_diag_scale_acc('R', m, n, kl, ku,
                                               alpha, &one, 0, a, 1, lda,
                                               beta,           b, 1, ldb);
#else
#error "Sanity: suzerain_blas_sgb_acc not implemented for BLAS"
#endif
}

inline int
suzerain_blas_dgb_acc(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const double alpha,
        const double *a,
        const int lda,
        const double beta,
        double *b,
        const int ldb)
{
#ifdef SUZERAIN_HAVE_MKL
    /* Internal dgb_acc because MKL lacks the routine. */
    const double one = 1.0;
    return suzerain_blasext_dgb_diag_scale_acc('R', m, n, kl, ku,
                                               alpha, &one, 0, a, 1, lda,
                                               beta,           b, 1, ldb);
#else
#error "Sanity: suzerain_blas_dgb_acc not implemented for BLAS"
#endif
}

inline int
suzerain_blas_cgb_acc(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const complex_float alpha,
        const complex_float *a,
        const int lda,
        const complex_float beta,
        complex_float *b,
        const int ldb)
{
#ifdef SUZERAIN_HAVE_MKL
    /* Internal cgb_acc because MKL lacks the routine. */
    const complex_float one = 1.0f;
    return suzerain_blasext_cgb_diag_scale_acc('R', m, n, kl, ku,
                                               alpha, &one, 0, a, 1, lda,
                                               beta,           b, 1, ldb);
#else
#error "Sanity: suzerain_blas_cgb_acc not implemented for BLAS"
#endif
}

inline int
suzerain_blas_zgb_acc(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const complex_double alpha,
        const complex_double *a,
        const int lda,
        const complex_double beta,
        complex_double *b,
        const int ldb)
{
#ifdef SUZERAIN_HAVE_MKL
    /* Internal zgb_acc because MKL lacks the routine. */
    const complex_double one = 1.0;
    return suzerain_blasext_zgb_diag_scale_acc('R', m, n, kl, ku,
                                               alpha, &one, 0, a, 1, lda,
                                               beta,           b, 1, ldb);
#else
#error "Sanity: suzerain_blas_zgb_acc not implemented for BLAS"
#endif
}

inline int
suzerain_blas_zgb_acc_d(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const complex_double alpha,
        const double *a,
        const int lda,
        const complex_double beta,
        complex_double *b,
        const int ldb)
{
#ifdef SUZERAIN_HAVE_MKL
    const double one = 1.0;
    return suzerain_blasext_zgb_diag_scale_acc_d('R', m, n, kl, ku,
                                                alpha, &one, 0, a, 1, lda,
                                                beta,           b, 1, ldb);
#else
#error "Sanity: suzerain_blas_zgb_acc_d not implemented for BLAS"
#endif
}

void
suzerain_blas_csbmv_s_c_external(
        const char uplo,
        const int n,
        const int k,
        const complex_float alpha,
        const float *a,
        const int lda,
        const complex_float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy)
{
    float *x_re, *y_re;
    memcpy(&x_re, &x, sizeof(x));
    memcpy(&y_re, &y, sizeof(y));

#pragma warning(push,disable:1572)
    if (cimag(alpha) == 0.0 && cimag(beta) == 0.0) {
#pragma warning(pop)
        /* Real-valued alpha and beta: scale y as we go */
        suzerain_blas_ssbmv_external(uplo, n, k,
                                     creal(alpha), a, lda, x_re, 2*incx,
                                     creal(beta), y_re, 2*incy);
        suzerain_blas_ssbmv_external(uplo, n, k,
                                     -cimag(alpha), a, lda, x_re+1, 2*incx,
                                     1.0, y_re, 2*incy);
        suzerain_blas_ssbmv_external(uplo, n, k,
                                     creal(alpha), a, lda, x_re+1, 2*incx,
                                     creal(beta), y_re+1, 2*incy);
        suzerain_blas_ssbmv_external(uplo, n, k,
                                     cimag(alpha), a, lda, x_re, 2*incx,
                                     1.0, y_re+1, 2*incy);
    } else {
        /* Complex-valued alpha and/or beta: scale y and then accumulate */
        suzerain_blas_cscal(n, beta, y, incy); /* NB cscal */
        suzerain_blas_ssbmv_external(uplo, n, k,
                                     creal(alpha), a, lda, x_re, 2*incx,
                                     1.0, y_re, 2*incy);
        suzerain_blas_ssbmv_external(uplo, n, k,
                                     cimag(alpha), a, lda, x_re, 2*incx,
                                     1.0, y_re+1, 2*incy);
        suzerain_blas_ssbmv_external(uplo, n, k,
                                     creal(alpha), a, lda, x_re+1, 2*incx,
                                     1.0, y_re+1, 2*incy);
        suzerain_blas_ssbmv_external(uplo, n, k,
                                     -cimag(alpha), a, lda, x_re+1, 2*incx,
                                     1.0, y_re, 2*incy);
    }
}

void
suzerain_blas_zsbmv_d_z_external(
        const char uplo,
        const int n,
        const int k,
        const complex_double alpha,
        const double *a,
        const int lda,
        const complex_double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy)
{
    double *x_re, *y_re;
    memcpy(&x_re, &x, sizeof(x));
    memcpy(&y_re, &y, sizeof(y));

#pragma warning(push,disable:1572)
    if (cimag(alpha) == 0.0 && cimag(beta) == 0.0) {
#pragma warning(pop)
        /* Real-valued alpha and beta: scale y as we go */
        suzerain_blas_dsbmv_external(uplo, n, k,
                                     creal(alpha), a, lda, x_re, 2*incx,
                                     creal(beta), y_re, 2*incy);
        suzerain_blas_dsbmv_external(uplo, n, k,
                                     -cimag(alpha), a, lda, x_re+1, 2*incx,
                                     1.0, y_re, 2*incy);
        suzerain_blas_dsbmv_external(uplo, n, k,
                                     creal(alpha), a, lda, x_re+1, 2*incx,
                                     creal(beta), y_re+1, 2*incy);
        suzerain_blas_dsbmv_external(uplo, n, k,
                                     cimag(alpha), a, lda, x_re, 2*incx,
                                     1.0, y_re+1, 2*incy);
    } else {
        /* Complex-valued alpha and/or beta: scale y and then accumulate */
        suzerain_blas_zscal(n, beta, y, incy); /* NB zscal */
        suzerain_blas_dsbmv_external(uplo, n, k,
                                     creal(alpha), a, lda, x_re, 2*incx,
                                     1.0, y_re, 2*incy);
        suzerain_blas_dsbmv_external(uplo, n, k,
                                     cimag(alpha), a, lda, x_re, 2*incx,
                                     1.0, y_re+1, 2*incy);
        suzerain_blas_dsbmv_external(uplo, n, k,
                                     creal(alpha), a, lda, x_re+1, 2*incx,
                                     1.0, y_re+1, 2*incy);
        suzerain_blas_dsbmv_external(uplo, n, k,
                                     -cimag(alpha), a, lda, x_re+1, 2*incx,
                                     1.0, y_re, 2*incy);
    }
}
