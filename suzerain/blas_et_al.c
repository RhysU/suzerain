/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2010, 2011, 2012 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *--------------------------------------------------------------------------
 * blas_et_al.c: wraps external implementations of BLAS, LAPACK, et al.
 * $Id$
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.h>
#pragma hdrstop
#include <suzerain/blas_et_al.h>
#include <suzerain/gbdddddmv.h>
#include <suzerain/gbddddmv.h>
#include <suzerain/gbdddmv.h>
#include <suzerain/gbddmv.h>
#include <suzerain/gbdmv.h>
#include <suzerain/gbmv.h>
#include <suzerain/sbmv.h>

// Odd looking, unused, anonymous enumerations are globally-scoped versions
// of "Compile Time Assertions" by Ralf Holly (http://drdobbs.com/184401873)
#define assert_static(e) do{ enum { assert_static__## line = 1/(e) }; }while(0)

enum { // Required for strict aliasing workarounds
    assert_floatp  = 1/(sizeof(float*)  == sizeof(complex_float*) ),
    assert_doublep = 1/(sizeof(double*) == sizeof(complex_double*))
};

#ifdef SUZERAIN_HAVE_MKL
# include <mkl_types.h>
# include <mkl_blas.h>
# include <mkl_lapack.h>
# define BLAS_FUNC(name,NAME)   name             /* Think AC_F77_WRAPPERS */
# define LAPACK_FUNC(name,NAME) name             /* Think AC_F77_WRAPPERS */
enum { // Required for binary interoperability with MKL
    assert_mkl_i = 1 /(sizeof(MKL_INT)       == sizeof(int)           ),
    assert_mkl_c = 1 /(sizeof(MKL_Complex8)  == sizeof(complex_float) ),
    assert_mkl_z = 1 /(sizeof(MKL_Complex16) == sizeof(complex_double))
};
#else
# error "No suitable BLAS and/or LAPACK library found during configuration"
#endif

// Shorthand
static inline int imin(int a, int b) { return a < b ? a : b; }
static inline int imax(int a, int b) { return a > b ? a : b; }
#define UNLIKELY(expr) SUZERAIN_UNLIKELY(expr)

// Many of the short methods have "inline" though their declarations do not.
// This allows inlining them later within this particular translation unit so
// we don't pay for needless function call overhead when our BLAS wrapper
// routines call other BLAS wrapper routines (e.g. within suzerain_blasext_*).

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
            &n, (const void *) &alpha, (void *) x, &incx);
}

inline void
suzerain_blas_zscal(
        const int n,
        const complex_double alpha,
        complex_double *x,
        const int incx)
{
    return BLAS_FUNC(zscal,ZSCAL)(
            &n, (const void *) &alpha, (void *) x, &incx);
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
    return BLAS_FUNC(cswap,CSWAP)(&n, (void *) x, &incx, (void *) y, &incy);
}

inline void
suzerain_blas_zswap(
        const int n,
        complex_double *x,
        const int incx,
        complex_double *y,
        const int incy)
{
    return BLAS_FUNC(zswap,ZSWAP)(&n, (void *) x, &incx, (void *) y, &incy);
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
            &n, (const void *) x, &incx, (void *) y, &incy);
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
            &n, (const void *) x, &incx, (void *) y, &incy);
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
    return BLAS_FUNC(cdotc,CDOTC)((      void *)dotc, &n,
                                  (const void *)x,    &incx,
                                  (const void *)y,    &incy);
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
    return BLAS_FUNC(zdotc,ZDOTC)((      void *)dotc, &n,
                                  (const void *)x,    &incx,
                                  (const void *)y,    &incy);
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
    return BLAS_FUNC(scnrm2,SCNRM2)(&n, (const void *) x, &incx);
}

inline double
suzerain_blas_dznrm2(
        const int n,
        const complex_double *x,
        const int incx)
{
    return BLAS_FUNC(dznrm2,DZNRM2)(&n, (const void *) x, &incx);
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
    return BLAS_FUNC(scasum,SCASUM)(&n, (void *) x, &incx);
}

inline double
suzerain_blas_dzasum(
        const int n,
        const complex_double *x,
        const int incx)
{
    return BLAS_FUNC(dzasum,DZASUM)(&n, (void *) x, &incx);
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
    return icamax(&n, (MKL_Complex8 *) x, &incx) - 1 /* zero-indexed */;
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
    return izamax(&n, (MKL_Complex16 *) x, &incx) - 1 /* zero-indexed */;
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
    return icamin(&n, (MKL_Complex8 *) x, &incx) - 1 /* zero-indexed */;
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
    return izamin(&n, (MKL_Complex16 *) x, &incx) - 1 /* zero-indexed */;
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
    return BLAS_FUNC(caxpy,CAXPY)(&n, (const void *) &alpha,
                                      (const void *) x, &incx,
                                      (      void *) y, &incy);
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
    return BLAS_FUNC(zaxpy,ZAXPY)(&n, (const void *) &alpha,
                                      (const void *) x, &incx,
                                      (      void *) y, &incy);
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
    BLAS_FUNC(caxpby,CAXPBY)(&n,
                             (const void *) &alpha, (const void *) x, &incx,
                             (const void *) &beta,  (      void *) y, &incy);
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
    BLAS_FUNC(zaxpby,ZAXPBY)(&n,
                             (const void *) &alpha, (const void *) x, &incx,
                             (const void *) &beta,  (      void *) y, &incy);
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

inline int
suzerain_lapack_sgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        float *ab,
        const int ldab,
        int *ipiv)
{
    int info = 0;
    LAPACK_FUNC(sgbtrf,SGBTRF)((int*)&m, (int*)&n, (int*)&kl, (int*)&ku,
                               ab, (int*)&ldab, ipiv, &info);
    return info;
}

inline int
suzerain_lapack_dgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        double *ab,
        const int ldab,
        int *ipiv)
{
    int info = 0;
    LAPACK_FUNC(dgbtrf,DGBTRF)((int*)&m, (int*)&n, (int*)&kl, (int*)&ku,
                               ab, (int*)&ldab, ipiv, &info);
    return info;
}

inline int
suzerain_lapack_cgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        complex_float *ab,
        const int ldab,
        int *ipiv)
{
    int info = 0;
    LAPACK_FUNC(cgbtrf,CGBTRF)((int*)&m, (int*)&n, (int*)&kl, (int*)&ku,
                               (void*)ab, (int *)&ldab, ipiv, &info);
    return info;
}

inline int
suzerain_lapack_zgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        complex_double *ab,
        const int ldab,
        int *ipiv)
{
    int info = 0;
    LAPACK_FUNC(zgbtrf,ZGBTRF)((int*)&m, (int*)&n, (int*)&kl, (int*)&ku,
                               (void*)ab, (int *)&ldab, ipiv, &info);
    return info;
}

inline int
suzerain_lapack_sgbtrs(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        const float *ab,
        const int ldab,
        const int *ipiv,
        float *b,
        const int ldb)
{
    int info = 0;
    LAPACK_FUNC(sgbtrs,SGBTRS)((char*)&trans, (int*)&n, (int*)&kl, (int*)&ku,
                               (int*)&nrhs, (float *)ab, (int*)&ldab,
                               (int *)ipiv, b, (int*)&ldb, &info);
    return info;
}

inline int
suzerain_lapack_dgbtrs(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        const double *ab,
        const int ldab,
        const int *ipiv,
        double *b,
        const int ldb)
{
    int info = 0;
    LAPACK_FUNC(dgbtrs,DGBTRS)((char*)&trans, (int*)&n, (int*)&kl, (int*)&ku,
                               (int*)&nrhs, (double *)ab, (int*)&ldab,
                               (int *)ipiv, b, (int*)&ldb, &info);
    return info;
}

inline int
suzerain_lapack_cgbtrs(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        const complex_float *ab,
        const int ldab,
        const int *ipiv,
        complex_float *b,
        const int ldb)
{
    int info = 0;
    LAPACK_FUNC(cgbtrs,CGBTRS)((char*)&trans, (int*)&n, (int*)&kl, (int*)&ku,
                               (int*)&nrhs, (void*)ab, (int*)&ldab,
                               (int *)ipiv, (void*)b,  (int*)&ldb, &info);
    return info;
}

inline int
suzerain_lapack_zgbtrs(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        const complex_double *ab,
        const int ldab,
        const int *ipiv,
        complex_double *b,
        const int ldb)
{
    int info = 0;
    LAPACK_FUNC(zgbtrs,ZGBTRS)((char*)&trans, (int*)&n, (int*)&kl, (int*)&ku,
                               (int*)&nrhs, (void*)ab, (int*)&ldab,
                               (int *)ipiv, (void*)b,  (int*)&ldb, &info);
    return info;
}

inline int
suzerain_lapack_sgbcon(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const float *ab,
        const int ldab,
        const int *ipiv,
        const float anorm,
        float *rcond,
        float *work,
        int *iwork)
{
    int info = 0;
    LAPACK_FUNC(sgbcon,SGBCON)((char*)&norm, (int*)&n, (int*)&kl, (int*)&ku,
                               (float*)ab, (int*)&ldab, (int*) ipiv,
                               (float*)&anorm, rcond, work, iwork, &info);
    return info;
}

inline int
suzerain_lapack_dgbcon(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const double *ab,
        const int ldab,
        const int *ipiv,
        const double anorm,
        double *rcond,
        double *work,
        int *iwork)
{
    int info = 0;
    LAPACK_FUNC(dgbcon,DGBCON)((char*)&norm, (int*)&n, (int*)&kl, (int*)&ku,
                               (double*)ab, (int*)&ldab, (int*) ipiv,
                               (double*)&anorm, rcond, work, iwork, &info);
    return info;
}

inline int
suzerain_lapack_cgbcon(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const complex_float *ab,
        const int ldab,
        const int *ipiv,
        const float anorm,
        float *rcond,
        complex_float *work,
        float  *rwork)
{
    int info = 0;
    LAPACK_FUNC(cgbcon,CGBCON)((char*)&norm, (int*)&n, (int*)&kl, (int*)&ku,
                               (void*)ab, (int*)&ldab, (int*)ipiv,
                               (float*)&anorm, rcond, (void*)work, rwork,
                               &info);
    return info;
}

inline int
suzerain_lapack_zgbcon(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const complex_double *ab,
        const int ldab,
        const int *ipiv,
        const double anorm,
        double *rcond,
        complex_double *work,
        double  *rwork)
{
    int info = 0;
    LAPACK_FUNC(zgbcon,ZGBCON)((char*)&norm, (int*)&n, (int*)&kl, (int*)&ku,
                               (void*)ab, (int*)&ldab, (int*) ipiv,
                               (double*)&anorm, rcond, (void*)work, rwork,
                               &info);
    return info;
}

inline int
suzerain_lapack_sgbsv(
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        float *ab,
        const int ldab,
        int *ipiv,
        float *b,
        const int ldb)
{
    int info;
    LAPACK_FUNC(sgbsv,SGBSV)((int*)&n, (int*)&kl, (int*)&ku, (int*)&nrhs,
                             ab, (int*)&ldab, ipiv, b, (int*)&ldb, &info);
    return info;
}

inline int
suzerain_lapack_dgbsv(
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        double *ab,
        const int ldab,
        int *ipiv,
        double *b,
        const int ldb)
{
    int info;
    LAPACK_FUNC(dgbsv,DGBSV)((int*)&n, (int*)&kl, (int*)&ku, (int*)&nrhs,
                             ab, (int*)&ldab, ipiv, b, (int*)&ldb, &info);
    return info;
}

inline int
suzerain_lapack_cgbsv(
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        complex_float *ab,
        const int ldab,
        int *ipiv,
        complex_float *b,
        const int ldb)
{
    int info;
    LAPACK_FUNC(cgbsv,CGBSV)((int*)&n, (int*)&kl, (int*)&ku, (int*)&nrhs,
                             (void*)ab, (int*)&ldab, ipiv,
                             (void*)b, (int*)&ldb, &info);
    return info;
}

inline int
suzerain_lapack_zgbsv(
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        complex_double *ab,
        const int ldab,
        int *ipiv,
        complex_double *b,
        const int ldb)
{
    int info;
    zgbsv((int*)&n, (int*)&kl, (int*)&ku, (int*)&nrhs,
          (void*)ab, (int*)&ldab, ipiv,
          (void*)b, (int*)&ldb, &info);
    return info;
}

inline int
suzerain_lapack_sgbsvx(
        const char fact,
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        float *ab,
        const int ldab,
        float *afb,
        const int ldafb,
        int *ipiv,
        char *equed,
        float *r,
        float *c,
        float *b,
        const int ldb,
        float *x,
        const int ldx,
        float *rcond,
        float *ferr,
        float *berr,
        float *work,
        int *iwork)
{
    int info;
    LAPACK_FUNC(sgbsvx,SGBSVX)((char*)&fact, (char*)&trans,
                               (int*)&n, (int*)&kl, (int*)&ku, (int*)&nrhs, ab,
                               (int*)&ldab, afb, (int*)&ldafb, ipiv, equed, r,
                               c, b, (int*)&ldb, x, (int*)&ldx, rcond, ferr,
                               berr, work, iwork, &info);
    return info;
}

inline int
suzerain_lapack_dgbsvx(
        const char fact,
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        double *ab,
        const int ldab,
        double *afb,
        const int ldafb,
        int *ipiv,
        char *equed,
        double *r,
        double *c,
        double *b,
        const int ldb,
        double *x,
        const int ldx,
        double *rcond,
        double *ferr,
        double *berr,
        double *work,
        int *iwork)
{
    int info;
    LAPACK_FUNC(dgbsvx,DGBSVX)((char*)&fact, (char*)&trans, (int*)&n,
                               (int*)&kl, (int*)&ku, (int*)&nrhs, ab,
                               (int*)&ldab, afb, (int*)&ldafb, ipiv, equed, r,
                               c, b, (int*)&ldb, x, (int*)&ldx, rcond, ferr,
                               berr, work, iwork, &info);
    return info;
}

inline int
suzerain_lapack_cgbsvx(
        const char fact,
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        complex_float *ab,
        const int ldab,
        complex_float *afb,
        const int ldafb,
        int *ipiv,
        char *equed,
        float *r,
        float *c,
        complex_float *b,
        const int ldb,
        complex_float *x,
        const int ldx,
        float *rcond,
        float *ferr,
        float *berr,
        complex_float *work,
        float *rwork)
{
    int info;
    LAPACK_FUNC(cgbsvx,CGBSVX)((char*)&fact, (char*)&trans, (int*)&n,
                               (int*)&kl, (int*)&ku, (int*)&nrhs, (void*)ab,
                               (int*)&ldab, (void*)afb, (int*)&ldafb, ipiv,
                               equed, r, c, (void*)b, (int*)&ldb, (void*)x,
                               (int*)&ldx, rcond, ferr, berr, (void*)work,
                               rwork, &info);
    return info;
}

inline int
suzerain_lapack_zgbsvx(
        const char fact,
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        complex_double *ab,
        const int ldab,
        complex_double *afb,
        const int ldafb,
        int *ipiv,
        char *equed,
        double *r,
        double *c,
        complex_double *b,
        const int ldb,
        complex_double *x,
        const int ldx,
        double *rcond,
        double *ferr,
        double *berr,
        complex_double *work,
        double *rwork)
{
    int info;
    LAPACK_FUNC(zgbsvx,ZGBSVX)((char*)&fact, (char*)&trans,
                               (int*)&n, (int*)&kl, (int*)&ku,
                               (int*)&nrhs, (void*)ab, (int*)&ldab,
                               (void*)afb, (int*)&ldafb, ipiv, equed,
                               r, c, (void*)b, (int*)&ldb, (void*)x,
                               (int*)&ldx, rcond, ferr, berr, (void*)work,
                               rwork, &info);
    return info;
}

inline int
suzerain_lapack_sgbrfs(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        float *ab,
        const int ldab,
        float *afb,
        const int ldafb,
        int *ipiv,
        float *b,
        const int ldb,
        float *x,
        const int ldx,
        float *ferr,
        float *berr,
        float *work,
        int *iwork)
{
    int info;
    LAPACK_FUNC(sgbrfs,SGBRFS)((char*)&trans,
                               (int*)&n, (int*)&kl, (int*)&ku,
                               (int*)&nrhs, ab, (int*)&ldab, afb,
                               (int*)&ldafb, ipiv, b, (int*)&ldb, x,
                               (int*)&ldx, ferr, berr, work, iwork,
                               &info);
    return info;
}

inline int
suzerain_lapack_dgbrfs(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        double *ab,
        const int ldab,
        double *afb,
        const int ldafb,
        int *ipiv,
        double *b,
        const int ldb,
        double *x,
        const int ldx,
        double *ferr,
        double *berr,
        double *work,
        int *iwork)
{
    int info;
    LAPACK_FUNC(dgbrfs,DGBRFS)((char*)&trans, (int*)&n,
                               (int*)&kl, (int*)&ku, (int*)&nrhs, ab,
                               (int*)&ldab, afb, (int*)&ldafb, ipiv, b,
                               (int*)&ldb, x, (int*)&ldx, ferr, berr,
                               work, iwork, &info);
    return info;
}

inline int
suzerain_lapack_cgbrfs(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        complex_float *ab,
        const int ldab,
        complex_float *afb,
        const int ldafb,
        int *ipiv,
        complex_float *b,
        const int ldb,
        complex_float *x,
        const int ldx,
        float *ferr,
        float *berr,
        complex_float *work,
        float *rwork)
{
    int info;
    LAPACK_FUNC(cgbrfs,CGBRFS)((char*)&trans, (int*)&n,
                               (int*)&kl, (int*)&ku, (int*)&nrhs,
                               (void*)ab, (int*)&ldab, (void*)afb,
                               (int*)&ldafb, ipiv, (void*)b, (int*)&ldb,
                               (void*)x, (int*)&ldx, ferr, berr,
                               (void*)work, rwork, &info);
    return info;
}

inline int
suzerain_lapack_zgbrfs(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        complex_double *ab,
        const int ldab,
        complex_double *afb,
        const int ldafb,
        int *ipiv,
        complex_double *b,
        const int ldb,
        complex_double *x,
        const int ldx,
        double *ferr,
        double *berr,
        complex_double *work,
        double *rwork)
{
    int info;
    LAPACK_FUNC(zgbrfs,ZGBRFS)((char*)&trans,
                               (int*)&n, (int*)&kl, (int*)&ku,
                               (int*)&nrhs, (void*)ab, (int*)&ldab,
                               (void*)afb, (int*)&ldafb, ipiv, (void*)b,
                               (int*)&ldb, (void*)x, (int*)&ldx, ferr,
                               berr, (void*)work, rwork, &info);
    return info;
}

inline float
suzerain_lapack_slangb(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const float *ab,
        const int ldab,
        float *work)
{
    return LAPACK_FUNC(slangb,SLANGB)((char*)&norm, (int*)&n, (int*)&kl,
                                      (int*)&ku, (float*)ab, (int*)&ldab,
                                      work);
}

inline double
suzerain_lapack_dlangb(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const double *ab,
        const int ldab,
        double *work)
{
    return LAPACK_FUNC(dlangb,DLANGB)((char*)&norm, (int*)&n, (int*)&kl,
                                      (int*)&ku, (double*)ab, (int*)&ldab,
                                      work);
}

inline float
suzerain_lapack_clangb(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const complex_float *ab,
        const int ldab,
        float *work)
{
    return LAPACK_FUNC(clangb,CLANGB)((char*)&norm, (int*)&n, (int*)&kl,
                                      (int*)&ku, (void*)ab, (int*)&ldab, work);
}

inline double
suzerain_lapack_zlangb(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const complex_double *ab,
        const int ldab,
        double *work)
{
    return LAPACK_FUNC(zlangb,ZLANGB)((char*)&norm, (int*)&n, (int*)&kl,
                                      (int*)&ku, (void*)ab, (int*)&ldab, work);
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
