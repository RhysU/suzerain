/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010, 2011, 2012 The PECOS Development Team
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
 * blas_et_al.c: wraps external implementations of BLAS, LAPACK, et al.
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.h>
#pragma hdrstop

#ifdef SUZERAIN_HAVE_MKL
#include <mkl_types.h>
#include <mkl_blas.h>
#include <mkl_lapack.h>
#else
#error "No suitable BLAS and/or LAPACK library found during configuration"
#endif

#include <suzerain/blas_et_al.h>
#include <suzerain/gbdddmv.h>
#include <suzerain/gbiddmv.h>
#include <suzerain/gbddmv.h>
#include <suzerain/gbidmv.h>
#include <suzerain/gbdmv.h>
#include <suzerain/gbmv.h>
#include <suzerain/sbmv.h>

static inline int imin(int a, int b) { return a < b ? a : b; }
static inline int imax(int a, int b) { return a > b ? a : b; }

void
suzerain_blas_xerbla(const char *srname, const int info)
{
#ifdef SUZERAIN_HAVE_MKL
    const int lsrname = srname ? strlen(srname) : 0;
    xerbla(srname, &info, lsrname);
#else
#error "Sanity failure"
#endif
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

void
suzerain_blas_sswap(
        const int n,
        float *x,
        const int incx,
        float *y,
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    sswap(&n, x, &incx, y, &incy);
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_dswap(
        const int n,
        double *x,
        const int incx,
        double *y,
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    dswap(&n, x, &incx, y, &incy);
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_cswap(
        const int n,
        float  (*x)[2],
        const int incx,
        float  (*y)[2],
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    cswap(&n, (MKL_Complex8 *) x, &incx, (MKL_Complex8 *) y, &incy);
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_zswap(
        const int n,
        double (*x)[2],
        const int incx,
        double (*y)[2],
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    zswap(&n, (MKL_Complex16 *) x, &incx, (MKL_Complex16 *) y, &incy);
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_scopy(
        const int n,
        const float *x,
        const int incx,
        float *y,
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    scopy(&n, x, &incx, y, &incy);
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_dcopy(
        const int n,
        const double *x,
        const int incx,
        double *y,
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    dcopy(&n, x, &incx, y, &incy);
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_ccopy(
        const int n,
        const float (*x)[2],
        const int incx,
        float (*y)[2],
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    ccopy(&n,
          (const MKL_Complex8 *) x, &incx,
          (      MKL_Complex8 *) y, &incy);
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_zcopy(
        const int n,
        const double (*x)[2],
        const int incx,
        double (*y)[2],
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    zcopy(&n,
         (const MKL_Complex16 *) x, &incx,
         (      MKL_Complex16 *) y, &incy);
#else
#error "Sanity failure"
#endif
}

float
suzerain_blas_sdot(
        const int n,
        const float *x,
        const int incx,
        const float *y,
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    return sdot(&n, x, &incx, y, &incy);
#else
#error "Sanity failure"
#endif
}

double
suzerain_blas_ddot(
        const int n,
        const double *x,
        const int incx,
        const double *y,
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    return ddot(&n, x, &incx, y, &incy);
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_cdotc(
        const int n,
        const float (*x)[2],
        const int incx,
        const float (*y)[2],
        const int incy,
        float dotc[2])
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    return cdotc((      MKL_Complex8 *)dotc, &n,
                 (const MKL_Complex8 *)x,    &incx,
                 (const MKL_Complex8 *)y,    &incy);
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_zdotc(
        const int n,
        const double (*x)[2],
        const int incx,
        const double (*y)[2],
        const int incy,
        double dotc[2])
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    return zdotc((      MKL_Complex16 *)dotc, &n,
                 (const MKL_Complex16 *)x,    &incx,
                 (const MKL_Complex16 *)y,    &incy);
#else
#error "Sanity failure"
#endif
}

float
suzerain_blas_snrm2(
        const int n,
        const float *x,
        const int incx)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    return snrm2(&n, x, &incx);
#else
#error "Sanity failure"
#endif
}

double
suzerain_blas_dnrm2(
        const int n,
        const double *x,
        const int incx)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    return dnrm2(&n, x, &incx);
#else
#error "Sanity failure"
#endif
}

float
suzerain_blas_scnrm2(
        const int n,
        const float (*x)[2],
        const int incx)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    return scnrm2(&n, (const MKL_Complex8 *) x, &incx);
#else
#error "Sanity failure"
#endif
}

double
suzerain_blas_dznrm2(
        const int n,
        const double (*x)[2],
        const int incx)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    return dznrm2(&n, (const MKL_Complex16 *) x, &incx);
#else
#error "Sanity failure"
#endif
}

float
suzerain_blas_sasum(
        const int n,
        const float *x,
        const int incx)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    return sasum(&n, x, &incx);
#else
#error "Sanity failure"
#endif
}

double
suzerain_blas_dasum(
        const int n,
        const double *x,
        const int incx)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    return dasum(&n, x, &incx);
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_saxpy(
        const int n,
        const float alpha,
        const float *x,
        const int incx,
        float *y,
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    saxpy(&n, &alpha, x, &incx, y, &incy);
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_daxpy(
        const int n,
        const double alpha,
        const double *x,
        const int incx,
        double *y,
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    daxpy(&n, &alpha, x, &incx, y, &incy);
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_caxpy(
        const int n,
        const float alpha[2],
        const float (*x)[2],
        const int incx,
        float (*y)[2],
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    caxpy(&n,
          (const MKL_Complex8 *) alpha,
          (const MKL_Complex8 *) x, &incx,
          (      MKL_Complex8 *) y, &incy);
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_zaxpy(
        const int n,
        const double alpha[2],
        const double (*x)[2],
        const int incx,
        double (*y)[2],
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    zaxpy(&n,
        (const MKL_Complex16 *) alpha,
        (const MKL_Complex16 *) x, &incx,
        (      MKL_Complex16 *) y, &incy);
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_saxpby(
        const int n,
        const float alpha,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
#pragma warning(push,disable:1572)
#ifdef SUZERAIN_HAVE_MKL
    /* Simulate saxpby since MKL lacks the routine. */
    if (SUZERAIN_UNLIKELY((alpha == 0.0f && beta == 1.0f) || n <= 0)) return;

    assert(sizeof(MKL_INT) == sizeof(int));
    if (beta != 1.0f)
        sscal(&n, &beta, y, &incy);
    saxpy(&n, &alpha, x, &incx, y, &incy);
#else
#error "Sanity failure"
#endif
#pragma warning(pop)
}

void
suzerain_blas_daxpby(
        const int n,
        const double alpha,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
#pragma warning(push,disable:1572)
#ifdef SUZERAIN_HAVE_MKL
    /* Simulate daxpby since MKL lacks the routine. */
    if (SUZERAIN_UNLIKELY((alpha == 0.0 && beta == 1.0) || n <= 0)) return;

    assert(sizeof(MKL_INT) == sizeof(int));
    if (beta != 1.0)
        dscal(&n, &beta, y, &incy);
    daxpy(&n, &alpha, x, &incx, y, &incy);
#else
#error "Sanity failure"
#endif
#pragma warning(pop)
}

void
suzerain_blas_caxpby(
        const int n,
        const float alpha[2],
        const float (*x)[2],
        const int incx,
        const float beta[2],
        float (*y)[2],
        const int incy)
{
#pragma warning(push,disable:1572)
#ifdef SUZERAIN_HAVE_MKL
    /* Simulate caxpby since MKL lacks the routine. */
    const int beta_is_one = (beta[0] == 1.0f && beta[1] == 0.0f);

    if (SUZERAIN_UNLIKELY((   alpha[0] == 0.0f && alpha[1] == 0.0f
                           && beta_is_one) || n <= 0)) {
        return;
    }

    if (!beta_is_one)
        suzerain_blas_cscal(n, beta, y, incy);
    suzerain_blas_caxpy(n, alpha, x, incx, y, incy);
#else
#error "Sanity failure"
#endif
#pragma warning(pop)
}

void
suzerain_blas_zaxpby(
        const int n,
        const double alpha[2],
        const double (*x)[2],
        const int incx,
        const double beta[2],
        double (*y)[2],
        const int incy)
{
#pragma warning(push,disable:1572)
#ifdef SUZERAIN_HAVE_MKL
    /* Simulate caxpby since MKL lacks the routine. */
    const int beta_is_one = (beta[0] == 1.0 && beta[1] == 0.0);

    if (SUZERAIN_UNLIKELY((   alpha[0] == 0.0 && alpha[1] == 0.0
                           && beta_is_one) || n <= 0)) {
        return;
    }

    if (!beta_is_one)
        suzerain_blas_zscal(n, beta, y, incy);
    suzerain_blas_zaxpy(n, alpha, x, incx, y, incy);
#else
#error "Sanity failure"
#endif
#pragma warning(pop)
}

void
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
#pragma warning(push,disable:1572)
#ifdef SUZERAIN_HAVE_MKL
    /* Simulate swaxpby since MKL lacks the routine. */
    assert(sizeof(MKL_INT) == sizeof(int));
    scopy(&n, y, &incy, w, &incw);
    if (beta != 1.0f)
        sscal(&n, &beta, w, &incw);
    saxpy(&n, &alpha, x, &incx, w, &incw);
#else
#error "Sanity failure"
#endif
#pragma warning(pop)
}

void
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
#pragma warning(push,disable:1572)
#ifdef SUZERAIN_HAVE_MKL
    /* Simulate dwaxpby since MKL lacks the routine. */
    assert(sizeof(MKL_INT) == sizeof(int));
    dcopy(&n, y, &incy, w, &incw);
    if (beta != 1.0)
        dscal(&n, &beta, w, &incw);
    daxpy(&n, &alpha, x, &incx, w, &incw);
#else
#error "Sanity failure"
#endif
#pragma warning(pop)
}

void
suzerain_blas_sscal(
        const int n,
        const float alpha,
        float *x,
        const int incx)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    sscal(&n, &alpha, x, &incx);
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_dscal(
        const int n,
        const double alpha,
        double *x,
        const int incx)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    dscal(&n, &alpha, x, &incx);
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_cscal(
        const int n,
        const float alpha[2],
        float (*x)[2],
        const int incx)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    cscal(&n,
         (const MKL_Complex8 *) alpha,
         (      MKL_Complex8 *) x, &incx);
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_zscal(
        const int n,
        const double alpha[2],
        double (*x)[2],
        const int incx)
{
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    zscal(&n,
         (const MKL_Complex16 *) alpha,
         (      MKL_Complex16 *) x, &incx);
#else
#error "Sanity failure"
#endif
}

void
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
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    sgbmv(&trans, &m, &n, &kl, &ku, &alpha, a, &lda,
          x, &incx, &beta, y, &incy);
#else
#error "Sanity failure"
#endif
}

void
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
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    dgbmv(&trans, &m, &n, &kl, &ku, &alpha, a, &lda,
          x, &incx, &beta, y, &incy);
#else
#error "Sanity failure"
#endif
}

void
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
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
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
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
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
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    ssbmv(&uplo, &n, &k, &alpha, a, &lda,
          x, &incx, &beta, y, &incy);
#else
#error "Sanity failure"
#endif
}

void
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
#ifdef SUZERAIN_HAVE_MKL
    assert(sizeof(MKL_INT) == sizeof(int));
    dsbmv(&uplo, &n, &k, &alpha, a, &lda,
          x, &incx, &beta, y, &incy);
#else
#error "Sanity failure"
#endif
}

void
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
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
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
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
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
    suzerain_blasext_sgb_diag_scale_acc(m, n, kl, ku,
                                        alpha, a, 1, lda, &one, 0,
                                        beta,  b, 1, ldb);
#else
#error "Sanity failure"
#endif
}

void
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
    suzerain_blasext_dgb_diag_scale_acc(m, n, kl, ku,
                                        alpha, a, 1, lda, &one, 0,
                                        beta,  b, 1, ldb);
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_cgb_acc(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const float alpha[2],
        const float (*a)[2],
        const int lda,
        const float beta[2],
        float (*b)[2],
        const int ldb)
{
#ifdef SUZERAIN_HAVE_MKL
    /* Internal cgb_acc because MKL lacks the routine. */
    const float one[2] = { 1.0f, 0.0f };
    suzerain_blasext_cgb_diag_scale_acc(m, n, kl, ku,
                                        alpha, a, 1, lda, &one, 0,
                                        beta,  b, 1, ldb);
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_zgb_acc(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const double alpha[2],
        const double (*a)[2],
        const int lda,
        const double beta[2],
        double (*b)[2],
        const int ldb)
{
#ifdef SUZERAIN_HAVE_MKL
    /* Internal zgb_acc because MKL lacks the routine. */
    const double one[2] = { 1.0, 0.0 };
    suzerain_blasext_zgb_diag_scale_acc(m, n, kl, ku,
                                        alpha, a, 1, lda, &one, 0,
                                        beta,  b, 1, ldb);
#else
#error "Sanity failure"
#endif
}

int
suzerain_lapack_sgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        float *ab,
        const int ldab,
        int *ipiv)
{
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    assert(sizeof(MKL_INT) == sizeof(int));
    int info = 0;
    sgbtrf((int*)&m, (int*)&n, (int*)&kl, (int*)&ku,
           ab, (int*)&ldab, ipiv, &info);
    return info;
#else
#error "Sanity failure"
#endif
}

int
suzerain_lapack_dgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        double *ab,
        const int ldab,
        int *ipiv)
{
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    assert(sizeof(MKL_INT) == sizeof(int));
    int info = 0;
    dgbtrf((int*)&m, (int*)&n, (int*)&kl, (int*)&ku,
           ab, (int*)&ldab, ipiv, &info);
    return info;
#else
#error "Sanity failure"
#endif
}

int
suzerain_lapack_cgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        float (*ab)[2],
        const int ldab,
        int *ipiv)
{
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    assert(sizeof(MKL_INT) == sizeof(int));
    int info = 0;
    cgbtrf((int*)&m, (int*)&n, (int*)&kl, (int*)&ku,
           (MKL_Complex8*)ab, (int *)&ldab, ipiv, &info);
    return info;
#else
#error "Sanity failure"
#endif
}

int
suzerain_lapack_zgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        double (*ab)[2],
        const int ldab,
        int *ipiv)
{
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    assert(sizeof(MKL_INT) == sizeof(int));
    int info = 0;
    zgbtrf((int*)&m, (int*)&n, (int*)&kl, (int*)&ku,
           (MKL_Complex16*)ab, (int *)&ldab, ipiv, &info);
    return info;
#else
#error "Sanity failure"
#endif
}

int
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
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    assert(sizeof(MKL_INT) == sizeof(int));
    int info = 0;
    sgbtrs((char*)&trans, (int*)&n, (int*)&kl, (int*)&ku, (int*)&nrhs,
        (float *)ab, (int*)&ldab, (int *)ipiv, b, (int*)&ldb, &info);
    return info;
#else
#error "Sanity failure"
#endif
}

int
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
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    assert(sizeof(MKL_INT) == sizeof(int));
    int info = 0;
    dgbtrs((char*)&trans, (int*)&n, (int*)&kl, (int*)&ku, (int*)&nrhs,
          (double *)ab, (int*)&ldab, (int *)ipiv, b, (int*)&ldb, &info);
    return info;
#else
#error "Sanity failure"
#endif
}

int
suzerain_lapack_cgbtrs(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        const float (*ab)[2],
        const int ldab,
        const int *ipiv,
        float (*b)[2],
        const int ldb)
{
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    assert(sizeof(MKL_INT) == sizeof(int));
    int info = 0;
    cgbtrs((char*)&trans, (int*)&n, (int*)&kl, (int*)&ku, (int*)&nrhs,
        (MKL_Complex8*)ab, (int*)&ldab, (int *)ipiv,
        (MKL_Complex8*)b,  (int*)&ldb, &info);
    return info;
#else
#error "Sanity failure"
#endif
}

int
suzerain_lapack_zgbtrs(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        const double (*ab)[2],
        const int ldab,
        const int *ipiv,
        double (*b)[2],
        const int ldb)
{
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    assert(sizeof(MKL_INT) == sizeof(int));
    int info = 0;
    zgbtrs((char*)&trans, (int*)&n, (int*)&kl, (int*)&ku, (int*)&nrhs,
        (MKL_Complex16*)ab, (int*)&ldab, (int *)ipiv,
        (MKL_Complex16*)b,  (int*)&ldb, &info);
    return info;
#else
#error "Sanity failure"
#endif
}

int
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
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    assert(sizeof(MKL_INT) == sizeof(int));
    int info = 0;
    sgbcon((char*)&norm, (int*)&n, (int*)&kl, (int*)&ku,
            (float*)ab, (int*)&ldab, (int*) ipiv, (float*)&anorm,
            rcond, work, iwork, &info);
    return info;
#else
#error "Sanity failure"
#endif
}

int
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
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    assert(sizeof(MKL_INT) == sizeof(int));
    int info = 0;
    dgbcon((char*)&norm, (int*)&n, (int*)&kl, (int*)&ku,
            (double*)ab, (int*)&ldab, (int*) ipiv, (double*)&anorm,
            rcond, work, iwork, &info);
    return info;
#else
#error "Sanity failure"
#endif
}

int
suzerain_lapack_cgbcon(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const float (*ab)[2],
        const int ldab,
        const int *ipiv,
        const float anorm,
        float *rcond,
        float (*work)[2],
        float  *rwork)
{
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    assert(sizeof(MKL_INT) == sizeof(int));
    int info = 0;
    cgbcon((char*)&norm, (int*)&n, (int*)&kl, (int*)&ku,
            (MKL_Complex8*)ab, (int*)&ldab, (int*)ipiv, (float*)&anorm,
            rcond, (MKL_Complex8*)work, rwork, &info);
    return info;
#else
#error "Sanity failure"
#endif
}

int
suzerain_lapack_zgbcon(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const double (*ab)[2],
        const int ldab,
        const int *ipiv,
        const double anorm,
        double *rcond,
        double (*work)[2],
        double  *rwork)
{
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    assert(sizeof(MKL_INT) == sizeof(int));
    int info = 0;
    zgbcon((char*)&norm, (int*)&n, (int*)&kl, (int*)&ku,
            (MKL_Complex16*)ab, (int*)&ldab, (int*) ipiv, (double*)&anorm,
            rcond, (MKL_Complex16*)work, rwork, &info);
    return info;
#else
#error "Sanity failure"
#endif
}

int
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
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    assert(sizeof(MKL_INT) == sizeof(int));
    int info;
    sgbsvx((char*)&fact, (char*)&trans,
           (int*)&n, (int*)&kl, (int*)&ku, (int*)&nrhs, ab, (int*)&ldab,
           afb, (int*)&ldafb, ipiv, equed, r, c, b, (int*)&ldb, x,
           (int*)&ldx, rcond, ferr, berr, work, iwork, &info);
    return info;
#else
#error "Sanity failure"
#endif
}

int
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
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    assert(sizeof(MKL_INT) == sizeof(int));
    int info;
    dgbsvx((char*)&fact, (char*)&trans,
           (int*)&n, (int*)&kl, (int*)&ku, (int*)&nrhs, ab, (int*)&ldab,
           afb, (int*)&ldafb, ipiv, equed, r, c, b, (int*)&ldb, x,
           (int*)&ldx, rcond, ferr, berr, work, iwork, &info);
    return info;
#else
#error "Sanity failure"
#endif
}

int
suzerain_lapack_cgbsvx(
        const char fact,
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        float (*ab)[2],
        const int ldab,
        float (*afb)[2],
        const int ldafb,
        int *ipiv,
        char *equed,
        float *r,
        float *c,
        float (*b)[2],
        const int ldb,
        float (*x)[2],
        const int ldx,
        float *rcond,
        float *ferr,
        float *berr,
        float (*work)[2],
        float *rwork)
{
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    assert(sizeof(MKL_INT) == sizeof(int));
    int info;
    cgbsvx((char*)&fact, (char*)&trans,
           (int*)&n, (int*)&kl, (int*)&ku, (int*)&nrhs, (MKL_Complex8*)ab,
           (int*)&ldab, (MKL_Complex8*)afb, (int*)&ldafb, ipiv, equed, r,
           c, (MKL_Complex8*)b, (int*)&ldb, (MKL_Complex8*)x, (int*)&ldx,
           rcond, ferr, berr, (MKL_Complex8*)work, rwork, &info);
    return info;
#else
#error "Sanity failure"
#endif
}

int
suzerain_lapack_zgbsvx(
        const char fact,
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        double (*ab)[2],
        const int ldab,
        double (*afb)[2],
        const int ldafb,
        int *ipiv,
        char *equed,
        double *r,
        double *c,
        double (*b)[2],
        const int ldb,
        double (*x)[2],
        const int ldx,
        double *rcond,
        double *ferr,
        double *berr,
        double (*work)[2],
        double *rwork)
{
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    assert(sizeof(MKL_INT) == sizeof(int));
    int info;
    zgbsvx((char*)&fact, (char*)&trans,
           (int*)&n, (int*)&kl, (int*)&ku, (int*)&nrhs, (MKL_Complex16*)ab,
           (int*)&ldab, (MKL_Complex16*)afb, (int*)&ldafb, ipiv, equed, r,
           c, (MKL_Complex16*)b, (int*)&ldb, (MKL_Complex16*)x, (int*)&ldx,
           rcond, ferr, berr, (MKL_Complex16*)work, rwork, &info);
    return info;
#else
#error "Sanity failure"
#endif
}

float
suzerain_lapack_slangb(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const float *ab,
        const int ldab,
        float *work)
{
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    assert(sizeof(MKL_INT) == sizeof(int));
    return slangb((char*)&norm, (int*)&n, (int*)&kl, (int*)&ku,
                  (float*)ab, (int*)&ldab, work);
#else
#error "Sanity failure"
#endif
}

double
suzerain_lapack_dlangb(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const double *ab,
        const int ldab,
        double *work)
{
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    assert(sizeof(MKL_INT) == sizeof(int));
    return dlangb((char*)&norm, (int*)&n, (int*)&kl, (int*)&ku,
                  (double*)ab, (int*)&ldab, work);
#else
#error "Sanity failure"
#endif
}

float
suzerain_lapack_clangb(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const float (*ab)[2],
        const int ldab,
        float *work)
{
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    assert(sizeof(MKL_INT) == sizeof(int));
    return clangb((char*)&norm, (int*)&n, (int*)&kl, (int*)&ku,
                  (MKL_Complex8*)ab, (int*)&ldab, work);
#else
#error "Sanity failure"
#endif
}

double
suzerain_lapack_zlangb(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const double (*ab)[2],
        const int ldab,
        double *work)
{
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    assert(sizeof(MKL_INT) == sizeof(int));
    return zlangb((char*)&norm, (int*)&n, (int*)&kl, (int*)&ku,
                  (MKL_Complex16*)ab, (int*)&ldab, work);
#else
#error "Sanity failure"
#endif
}

void
suzerain_blasext_daxpzy(
        const int n,
        const double alpha[2],
        const double * restrict x,
        const int incx,
        double (* restrict y)[2],
        const int incy)
{
    assert(incx >= 0); // TODO Handle negative incx
    assert(incy >= 0); // TODO Handle negative incy

    const double alpha_re = alpha[0];
    const double alpha_im = alpha[1];

    if (SUZERAIN_UNLIKELY(incx != 1 || incy != 1)) {
        /* General stride case */
#pragma unroll
        for (int i = 0; i < n; ++i) {
            const double      xi    = x[i * incx];
            double * restrict yi    = y[i * incy];

            yi[0] += alpha_re*xi;
            yi[1] += alpha_im*xi;
        }
    } else {
        /* Unit stride case */
#pragma unroll
        for (int i = 0; i < n; ++i) {
            const double      xi    = x[i];
            double * restrict yi    = y[i];

            yi[0] += alpha_re*xi;
            yi[1] += alpha_im*xi;
        }
    }
}

void
suzerain_blasext_daxpzby(
        const int n,
        const double alpha[2],
        const double * restrict x,
        const int incx,
        const double beta[2],
        double (* restrict y)[2],
        const int incy)
{
    assert(incx >= 0); // TODO Handle negative incx
    assert(incy >= 0); // TODO Handle negative incy

#pragma warning(push,disable:1572)
    if (SUZERAIN_UNLIKELY((beta[0] == 1.0 && beta[1] == 0.0))) {
#pragma warning(pop)
        return suzerain_blasext_daxpzy(n, alpha, x, incx, y, incy);
    }

    const double alpha_re = alpha[0];
    const double alpha_im = alpha[1];
    const double beta_re  = beta[0];
    const double beta_im  = beta[1];

    if (SUZERAIN_UNLIKELY(incx != 1 || incy != 1)) {
        /* General stride case */
#pragma unroll
        for (int i = 0; i < n; ++i) {
            const double      xi    = x[i * incx];
            double * restrict yi    = y[i * incy];
            double            yi_re = yi[0];
            double            yi_im = yi[1];

            yi[0] = beta_re*yi_re - beta_im*yi_im + alpha_re*xi;
            yi[1] = beta_re*yi_im + beta_im*yi_re + alpha_im*xi;
        }
    } else {
        /* Unit stride case */
#pragma unroll
        for (int i = 0; i < n; ++i) {
            const double      xi    = x[i];
            double * restrict yi    = y[i];
            double            yi_re = yi[0];
            double            yi_im = yi[1];

            yi[0] = beta_re*yi_re - beta_im*yi_im + alpha_re*xi;
            yi[1] = beta_re*yi_im + beta_im*yi_re + alpha_im*xi;
        }
    }
}

void
suzerain_blasext_sgbmzv_external(
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
#pragma warning(push,disable:1572)
    if (alpha[1] == 0.0 && beta[1] == 0.0) {
#pragma warning(pop)
        /* Real-valued alpha and beta: scale y as we go */
        suzerain_blas_sgbmv_external(trans, m, n, kl, ku,
                                     alpha[0], a, lda, &(x[0][0]), 2*incx,
                                     beta[0], &(y[0][0]), 2*incy);
        suzerain_blas_sgbmv_external(trans, m, n, kl, ku,
                                     -alpha[1], a, lda, &(x[0][1]), 2*incx,
                                     1.0, &(y[0][0]), 2*incy);
        suzerain_blas_sgbmv_external(trans, m, n, kl, ku,
                                     alpha[0], a, lda, &(x[0][1]), 2*incx,
                                     beta[0], &(y[0][1]), 2*incy);
        suzerain_blas_sgbmv_external(trans, m, n, kl, ku,
                                     alpha[1], a, lda, &(x[0][0]), 2*incx,
                                     1.0, &(y[0][1]), 2*incy);
    } else {
        /* Complex-valued alpha and/or beta: scale y and then accumulate */
        suzerain_blas_cscal((toupper(trans) == 'N' ? m : n), beta, y, incy);
        suzerain_blas_sgbmv_external(trans, m, n, kl, ku,
                                     alpha[0], a, lda, &(x[0][0]), 2*incx,
                                     1.0, &(y[0][0]), 2*incy);
        suzerain_blas_sgbmv_external(trans, m, n, kl, ku,
                                     alpha[1], a, lda, &(x[0][0]), 2*incx,
                                     1.0, &(y[0][1]), 2*incy);
        suzerain_blas_sgbmv_external(trans, m, n, kl, ku,
                                     alpha[0], a, lda, &(x[0][1]), 2*incx,
                                     1.0, &(y[0][1]), 2*incy);
        suzerain_blas_sgbmv_external(trans, m, n, kl, ku,
                                     -alpha[1], a, lda, &(x[0][1]), 2*incx,
                                     1.0, &(y[0][0]), 2*incy);
    }
}

void
suzerain_blasext_dgbmzv_external(
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
#pragma warning(push,disable:1572)
    if (alpha[1] == 0.0 && beta[1] == 0.0) {
#pragma warning(pop)
        /* Real-valued alpha and beta: scale y as we go */
        suzerain_blas_dgbmv_external(trans, m, n, kl, ku,
                                     alpha[0], a, lda, &(x[0][0]), 2*incx,
                                     beta[0], &(y[0][0]), 2*incy);
        suzerain_blas_dgbmv_external(trans, m, n, kl, ku,
                                     -alpha[1], a, lda, &(x[0][1]), 2*incx,
                                     1.0, &(y[0][0]), 2*incy);
        suzerain_blas_dgbmv_external(trans, m, n, kl, ku,
                                     alpha[0], a, lda, &(x[0][1]), 2*incx,
                                     beta[0], &(y[0][1]), 2*incy);
        suzerain_blas_dgbmv_external(trans, m, n, kl, ku,
                                     alpha[1], a, lda, &(x[0][0]), 2*incx,
                                     1.0, &(y[0][1]), 2*incy);
    } else {
        /* Complex-valued alpha and/or beta: scale y and then accumulate */
        suzerain_blas_zscal((toupper(trans) == 'N' ? m : n), beta, y, incy);
        suzerain_blas_dgbmv_external(trans, m, n, kl, ku,
                                     alpha[0], a, lda, &(x[0][0]), 2*incx,
                                     1.0, &(y[0][0]), 2*incy);
        suzerain_blas_dgbmv_external(trans, m, n, kl, ku,
                                     alpha[1], a, lda, &(x[0][0]), 2*incx,
                                     1.0, &(y[0][1]), 2*incy);
        suzerain_blas_dgbmv_external(trans, m, n, kl, ku,
                                     alpha[0], a, lda, &(x[0][1]), 2*incx,
                                     1.0, &(y[0][1]), 2*incy);
        suzerain_blas_dgbmv_external(trans, m, n, kl, ku,
                                     -alpha[1], a, lda, &(x[0][1]), 2*incx,
                                     1.0, &(y[0][0]), 2*incy);
    }
}

void
suzerain_blasext_sgbmzv(
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
    const int info = suzerain_gbmv_sc(trans, m, n, kl, ku,
                                      alpha, a, lda, x, incx,
                                      beta,          y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_dgbmzv(
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
    const int info = suzerain_gbmv_dz(trans, m, n, kl, ku,
                                      alpha, a, lda, x, incx,
                                      beta,          y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_ssbmzv_external(
        const char uplo,
        const int n,
        const int k,
        const float alpha[2],
        const float *a,
        const int lda,
        const float (*x)[2],
        const int incx,
        const float beta[2],
        float (*y)[2],
        const int incy)
{
#pragma warning(push,disable:1572)
    if (alpha[1] == 0.0 && beta[1] == 0.0) {
#pragma warning(pop)
        /* Real-valued alpha and beta: scale y as we go */
        suzerain_blas_ssbmv_external(uplo, n, k,
                                     alpha[0], a, lda, &(x[0][0]), 2*incx,
                                     beta[0], &(y[0][0]), 2*incy);
        suzerain_blas_ssbmv_external(uplo, n, k,
                                     -alpha[1], a, lda, &(x[0][1]), 2*incx,
                                     1.0, &(y[0][0]), 2*incy);
        suzerain_blas_ssbmv_external(uplo, n, k,
                                     alpha[0], a, lda, &(x[0][1]), 2*incx,
                                     beta[0], &(y[0][1]), 2*incy);
        suzerain_blas_ssbmv_external(uplo, n, k,
                                     alpha[1], a, lda, &(x[0][0]), 2*incx,
                                     1.0, &(y[0][1]), 2*incy);
    } else {
        /* Complex-valued alpha and/or beta: scale y and then accumulate */
        suzerain_blas_cscal(n, beta, y, incy); /* NB cscal */
        suzerain_blas_ssbmv_external(uplo, n, k,
                                     alpha[0], a, lda, &(x[0][0]), 2*incx,
                                     1.0, &(y[0][0]), 2*incy);
        suzerain_blas_ssbmv_external(uplo, n, k,
                                     alpha[1], a, lda, &(x[0][0]), 2*incx,
                                     1.0, &(y[0][1]), 2*incy);
        suzerain_blas_ssbmv_external(uplo, n, k,
                                     alpha[0], a, lda, &(x[0][1]), 2*incx,
                                     1.0, &(y[0][1]), 2*incy);
        suzerain_blas_ssbmv_external(uplo, n, k,
                                     -alpha[1], a, lda, &(x[0][1]), 2*incx,
                                     1.0, &(y[0][0]), 2*incy);
    }
}

void
suzerain_blasext_dsbmzv_external(
        const char uplo,
        const int n,
        const int k,
        const double alpha[2],
        const double *a,
        const int lda,
        const double (*x)[2],
        const int incx,
        const double beta[2],
        double (*y)[2],
        const int incy)
{
#pragma warning(push,disable:1572)
    if (alpha[1] == 0.0 && beta[1] == 0.0) {
#pragma warning(pop)
        /* Real-valued alpha and beta: scale y as we go */
        suzerain_blas_dsbmv_external(uplo, n, k,
                                     alpha[0], a, lda, &(x[0][0]), 2*incx,
                                     beta[0], &(y[0][0]), 2*incy);
        suzerain_blas_dsbmv_external(uplo, n, k,
                                     -alpha[1], a, lda, &(x[0][1]), 2*incx,
                                     1.0, &(y[0][0]), 2*incy);
        suzerain_blas_dsbmv_external(uplo, n, k,
                                     alpha[0], a, lda, &(x[0][1]), 2*incx,
                                     beta[0], &(y[0][1]), 2*incy);
        suzerain_blas_dsbmv_external(uplo, n, k,
                                     alpha[1], a, lda, &(x[0][0]), 2*incx,
                                     1.0, &(y[0][1]), 2*incy);
    } else {
        /* Complex-valued alpha and/or beta: scale y and then accumulate */
        suzerain_blas_zscal(n, beta, y, incy); /* NB zscal */
        suzerain_blas_dsbmv_external(uplo, n, k,
                                     alpha[0], a, lda, &(x[0][0]), 2*incx,
                                     1.0, &(y[0][0]), 2*incy);
        suzerain_blas_dsbmv_external(uplo, n, k,
                                     alpha[1], a, lda, &(x[0][0]), 2*incx,
                                     1.0, &(y[0][1]), 2*incy);
        suzerain_blas_dsbmv_external(uplo, n, k,
                                     alpha[0], a, lda, &(x[0][1]), 2*incx,
                                     1.0, &(y[0][1]), 2*incy);
        suzerain_blas_dsbmv_external(uplo, n, k,
                                     -alpha[1], a, lda, &(x[0][1]), 2*incx,
                                     1.0, &(y[0][0]), 2*incy);
    }
}

void
suzerain_blasext_sgbdmv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha,
        const float *d,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    float * const z = suzerain_blas_malloc(n*sizeof(float));
    if (SUZERAIN_UNLIKELY(!z)) suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_sgbmv_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha, d, 1, z, 1, beta, y, incy);
    suzerain_blas_free(z);
}

void
suzerain_blasext_sgbdmv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha,
        const float *d,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    const int info = suzerain_gbdmv_s(trans, n, kl, ku,
                                      alpha, d,
                                      a, lda, x, incx,
                                      beta,   y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_dgbdmv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha,
        const double *d,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    double * const z = suzerain_blas_malloc(n*sizeof(double));
    if (SUZERAIN_UNLIKELY(!z)) suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_dgbmv_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha, d, 1, z, 1, beta, y, incy);
    suzerain_blas_free(z);
}

void
suzerain_blasext_dgbdmv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha,
        const double *d,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    const int info = suzerain_gbdmv_d(trans, n, kl, ku,
                                      alpha, d,
                                      a, lda, x, incx,
                                      beta,   y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_sgbdmzv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha[2],
        const float *d,
        const float *a,
        const int lda,
        const float (*x)[2],
        const int incx,
        const float beta[2],
        float (*y)[2],
        const int incy)
{
    static const float one[2]  = { 1, 0 };
    static const float zero[2] = { 0, 0 };
    float (* const z)[2] = suzerain_blas_malloc(2*n*sizeof(float));
    if (SUZERAIN_UNLIKELY(!z)) suzerain_blas_xerbla(__func__, -1);
    suzerain_blasext_sgbmzv_external(
            trans, n, n, kl, ku, one, a, lda, x, incx, zero, z, 1);
    suzerain_blasext_ssbmzv_external(
            'U', n, 0, alpha, d, 1, (void *) z, 1, beta, y, incy);
    suzerain_blas_free(z);
}

void
suzerain_blasext_sgbdmzv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha[2],
        const float *d,
        const float *a,
        const int lda,
        const float (*x)[2],
        const int incx,
        const float beta[2],
        float (*y)[2],
        const int incy)
{
    const int info = suzerain_gbdmv_sc(trans, n, kl, ku,
                                       alpha, d,
                                       a, lda, x, incx,
                                       beta,   y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_dgbdmzv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha[2],
        const double *d,
        const double *a,
        const int lda,
        const double (*x)[2],
        const int incx,
        const double beta[2],
        double (*y)[2],
        const int incy)
{
    static const double one[2]  = { 1, 0 };
    static const double zero[2] = { 0, 0 };
    double (* const z)[2] = suzerain_blas_malloc(2*n*sizeof(double));
    if (SUZERAIN_UNLIKELY(!z)) suzerain_blas_xerbla(__func__, -1);
    suzerain_blasext_dgbmzv_external(
            trans, n, n, kl, ku, one, a, lda, x, incx, zero, z, 1);
    suzerain_blasext_dsbmzv_external(
            'U', n, 0, alpha, d, 1, (void *) z, 1, beta, y, incy);
    suzerain_blas_free(z);
}

void
suzerain_blasext_dgbdmzv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha[2],
        const double *d,
        const double *a,
        const int lda,
        const double (*x)[2],
        const int incx,
        const double beta[2],
        double (*y)[2],
        const int incy)
{
    const int info = suzerain_gbdmv_dz(trans, n, kl, ku,
                                       alpha, d,
                                       a, lda, x, incx,
                                       beta,   y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_sgbddmv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha0,
        const float *d0,
        const float alpha1,
        const float *d1,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    float * const z = suzerain_blas_malloc(n*sizeof(float));
    if (SUZERAIN_UNLIKELY(!z)) suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_sgbmv_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha0, d0, 1, z, 1, beta, y, incy);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha1, d1, 1, z, 1, 1,    y, incy);
    suzerain_blas_free(z);
}

void
suzerain_blasext_sgbddmv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha0,
        const float *d0,
        const float alpha1,
        const float *d1,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    const int info = suzerain_gbddmv_s(trans, n, kl, ku,
                                       alpha0, d0, alpha1, d1,
                                       a, lda, x, incx,
                                       beta,   y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_dgbddmv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha0,
        const double *d0,
        const double alpha1,
        const double *d1,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    double * const z = suzerain_blas_malloc(n*sizeof(double));
    if (SUZERAIN_UNLIKELY(!z)) suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_dgbmv_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha0, d0, 1, z, 1, beta, y, incy);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha1, d1, 1, z, 1, 1,    y, incy);
    suzerain_blas_free(z);
}

void
suzerain_blasext_dgbddmv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha0,
        const double *d0,
        const double alpha1,
        const double *d1,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    const int info = suzerain_gbddmv_d(trans, n, kl, ku,
                                       alpha0, d0, alpha1, d1,
                                       a, lda, x, incx,
                                       beta,   y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_sgbddmzv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha0[2],
        const float *d0,
        const float alpha1[2],
        const float *d1,
        const float *a,
        const int lda,
        const float (*x)[2],
        const int incx,
        const float beta[2],
        float (*y)[2],
        const int incy)
{
    static const float one[2]  = { 1, 0 };
    static const float zero[2] = { 0, 0 };
    float (* const z)[2] = suzerain_blas_malloc(2*n*sizeof(float));
    if (SUZERAIN_UNLIKELY(!z)) suzerain_blas_xerbla(__func__, -1);
    suzerain_blasext_sgbmzv_external(
            trans, n, n, kl, ku, one, a, lda, x, incx, zero, z, 1);
    suzerain_blasext_ssbmzv_external(
            'U', n, 0, alpha0, d0, 1, (void *) z, 1, beta, y, incy);
    suzerain_blasext_ssbmzv_external(
            'U', n, 0, alpha1, d1, 1, (void *) z, 1, one,  y, incy);
    suzerain_blas_free(z);
}

void
suzerain_blasext_sgbddmzv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha0[2],
        const float *d0,
        const float alpha1[2],
        const float *d1,
        const float *a,
        const int lda,
        const float (*x)[2],
        const int incx,
        const float beta[2],
        float (*y)[2],
        const int incy)
{
    const int info = suzerain_gbddmv_sc(trans, n, kl, ku,
                                        alpha0, d0, alpha1, d1,
                                        a, lda, x, incx,
                                        beta,   y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_dgbddmzv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha0[2],
        const double *d0,
        const double alpha1[2],
        const double *d1,
        const double *a,
        const int lda,
        const double (*x)[2],
        const int incx,
        const double beta[2],
        double (*y)[2],
        const int incy)
{
    static const double one[2]  = { 1, 0 };
    static const double zero[2] = { 0, 0 };
    double (* const z)[2] = suzerain_blas_malloc(2*n*sizeof(double));
    if (SUZERAIN_UNLIKELY(!z)) suzerain_blas_xerbla(__func__, -1);
    suzerain_blasext_dgbmzv_external(
            trans, n, n, kl, ku, one, a, lda, x, incx, zero, z, 1);
    suzerain_blasext_dsbmzv_external(
            'U', n, 0, alpha0, d0, 1, (void *) z, 1, beta, y, incy);
    suzerain_blasext_dsbmzv_external(
            'U', n, 0, alpha1, d1, 1, (void *) z, 1, one,  y, incy);
    suzerain_blas_free(z);
}

void
suzerain_blasext_dgbddmzv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha0[2],
        const double *d0,
        const double alpha1[2],
        const double *d1,
        const double *a,
        const int lda,
        const double (*x)[2],
        const int incx,
        const double beta[2],
        double (*y)[2],
        const int incy)
{
    const int info = suzerain_gbddmv_dz(trans, n, kl, ku,
                                        alpha0, d0, alpha1, d1,
                                        a, lda, x, incx,
                                        beta,   y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_sgbidmv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha0,
        const float alpha1,
        const float *d1,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    float * const z = suzerain_blas_malloc(n*sizeof(float));
    if (SUZERAIN_UNLIKELY(!z)) suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_sgbmv_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_saxpby(
            n, alpha0, z, 1, beta, y, incy);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha1, d1, 1, z, 1, 1, y, incy);
    suzerain_blas_free(z);
}

void
suzerain_blasext_sgbidmv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha0,
        const float alpha1,
        const float *d1,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    const int info = suzerain_gbidmv_s(trans, n, kl, ku,
                                       alpha0, alpha1, d1,
                                       a, lda, x, incx,
                                       beta,   y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_dgbidmv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha0,
        const double alpha1,
        const double *d1,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    double * const z = suzerain_blas_malloc(n*sizeof(double));
    if (SUZERAIN_UNLIKELY(!z)) suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_dgbmv_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_daxpby(
            n, alpha0, z, 1, beta, y, incy);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha1, d1, 1, z, 1, 1, y, incy);
    suzerain_blas_free(z);
}

void
suzerain_blasext_dgbidmv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha0,
        const double alpha1,
        const double *d1,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    const int info = suzerain_gbidmv_d(trans, n, kl, ku,
                                       alpha0, alpha1, d1,
                                       a, lda, x, incx,
                                       beta,   y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_sgbidmzv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha0[2],
        const float alpha1[2],
        const float *d1,
        const float *a,
        const int lda,
        const float (*x)[2],
        const int incx,
        const float beta[2],
        float (*y)[2],
        const int incy)
{
    static const float one[2]  = { 1, 0 };
    static const float zero[2] = { 0, 0 };
    float (* const z)[2] = suzerain_blas_malloc(2*n*sizeof(float));
    if (SUZERAIN_UNLIKELY(!z)) suzerain_blas_xerbla(__func__, -1);
    suzerain_blasext_sgbmzv_external(
            trans, n, n, kl, ku, one, a, lda, x, incx, zero, z, 1);
    suzerain_blas_caxpby(
            n, alpha0, (void *) z, 1, beta, y, incy);
    suzerain_blasext_ssbmzv_external(
            'U', n, 0, alpha1, d1, 1, (void *) z, 1, one, y, incy);
    suzerain_blas_free(z);
}

void
suzerain_blasext_sgbidmzv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha0[2],
        const float alpha1[2],
        const float *d1,
        const float *a,
        const int lda,
        const float (*x)[2],
        const int incx,
        const float beta[2],
        float (*y)[2],
        const int incy)
{
    const int info = suzerain_gbidmv_sc(trans, n, kl, ku,
                                        alpha0, alpha1, d1,
                                        a, lda, x, incx,
                                        beta,   y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_dgbidmzv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha0[2],
        const double alpha1[2],
        const double *d1,
        const double *a,
        const int lda,
        const double (*x)[2],
        const int incx,
        const double beta[2],
        double (*y)[2],
        const int incy)
{
    static const double one[2]  = { 1, 0 };
    static const double zero[2] = { 0, 0 };
    double (* const z)[2] = suzerain_blas_malloc(2*n*sizeof(double));
    if (SUZERAIN_UNLIKELY(!z)) suzerain_blas_xerbla(__func__, -1);
    suzerain_blasext_dgbmzv_external(
            trans, n, n, kl, ku, one, a, lda, x, incx, zero, z, 1);
    suzerain_blas_zaxpby(
            n, alpha0, (void *) z, 1, beta, y, incy);
    suzerain_blasext_dsbmzv_external(
            'U', n, 0, alpha1, d1, 1, (void *) z, 1, one, y, incy);
    suzerain_blas_free(z);
}

void
suzerain_blasext_dgbidmzv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha0[2],
        const double alpha1[2],
        const double *d1,
        const double *a,
        const int lda,
        const double (*x)[2],
        const int incx,
        const double beta[2],
        double (*y)[2],
        const int incy)
{
    const int info = suzerain_gbidmv_dz(trans, n, kl, ku,
                                        alpha0, alpha1, d1,
                                        a, lda, x, incx,
                                        beta,   y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_sgbdddmv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha0,
        const float *d0,
        const float alpha1,
        const float *d1,
        const float alpha2,
        const float *d2,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    float * const z = suzerain_blas_malloc(n*sizeof(float));
    if (SUZERAIN_UNLIKELY(!z)) suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_sgbmv_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha0, d0, 1, z, 1, beta, y, incy);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha1, d1, 1, z, 1, 1,    y, incy);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha2, d2, 1, z, 1, 1,    y, incy);
    suzerain_blas_free(z);
}

void
suzerain_blasext_sgbdddmv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha0,
        const float *d0,
        const float alpha1,
        const float *d1,
        const float alpha2,
        const float *d2,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    const int info = suzerain_gbdddmv_s(trans, n, kl, ku,
                                        alpha0, d0, alpha1, d1, alpha2, d2,
                                        a, lda, x, incx,
                                        beta,   y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_dgbdddmv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha0,
        const double *d0,
        const double alpha1,
        const double *d1,
        const double alpha2,
        const double *d2,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    double * const z = suzerain_blas_malloc(n*sizeof(double));
    if (SUZERAIN_UNLIKELY(!z)) suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_dgbmv_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha0, d0, 1, z, 1, beta, y, incy);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha1, d1, 1, z, 1, 1,    y, incy);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha2, d2, 1, z, 1, 1,    y, incy);
    suzerain_blas_free(z);
}

void
suzerain_blasext_dgbdddmv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha0,
        const double *d0,
        const double alpha1,
        const double *d1,
        const double alpha2,
        const double *d2,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    const int info = suzerain_gbdddmv_d(trans, n, kl, ku,
                                        alpha0, d0, alpha1, d1, alpha2, d2,
                                        a, lda, x, incx,
                                        beta,   y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_sgbdddmzv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha0[2],
        const float *d0,
        const float alpha1[2],
        const float *d1,
        const float alpha2[2],
        const float *d2,
        const float *a,
        const int lda,
        const float (*x)[2],
        const int incx,
        const float beta[2],
        float (*y)[2],
        const int incy)
{
    static const float one[2]  = { 1, 0 };
    static const float zero[2] = { 0, 0 };
    float (* const z)[2] = suzerain_blas_malloc(2*n*sizeof(float));
    if (SUZERAIN_UNLIKELY(!z)) suzerain_blas_xerbla(__func__, -1);
    suzerain_blasext_sgbmzv_external(
            trans, n, n, kl, ku, one, a, lda, x, incx, zero, z, 1);
    suzerain_blasext_ssbmzv_external(
            'U', n, 0, alpha0, d0, 1, (void *) z, 1, beta, y, incy);
    suzerain_blasext_ssbmzv_external(
            'U', n, 0, alpha1, d1, 1, (void *) z, 1, one,  y, incy);
    suzerain_blasext_ssbmzv_external(
            'U', n, 0, alpha2, d2, 1, (void *) z, 1, one,  y, incy);
    suzerain_blas_free(z);
}

void
suzerain_blasext_sgbdddmzv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha0[2],
        const float *d0,
        const float alpha1[2],
        const float *d1,
        const float alpha2[2],
        const float *d2,
        const float *a,
        const int lda,
        const float (*x)[2],
        const int incx,
        const float beta[2],
        float (*y)[2],
        const int incy)
{
    const int info = suzerain_gbdddmv_sc(trans, n, kl, ku,
                                         alpha0, d0, alpha1, d1, alpha2, d2,
                                         a, lda, x, incx,
                                         beta,   y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_dgbdddmzv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha0[2],
        const double *d0,
        const double alpha1[2],
        const double *d1,
        const double alpha2[2],
        const double *d2,
        const double *a,
        const int lda,
        const double (*x)[2],
        const int incx,
        const double beta[2],
        double (*y)[2],
        const int incy)
{
    static const double one[2]  = { 1, 0 };
    static const double zero[2] = { 0, 0 };
    double (* const z)[2] = suzerain_blas_malloc(2*n*sizeof(double));
    if (SUZERAIN_UNLIKELY(!z)) suzerain_blas_xerbla(__func__, -1);
    suzerain_blasext_dgbmzv_external(
            trans, n, n, kl, ku, one, a, lda, x, incx, zero, z, 1);
    suzerain_blasext_dsbmzv_external(
            'U', n, 0, alpha0, d0, 1, (void *) z, 1, beta, y, incy);
    suzerain_blasext_dsbmzv_external(
            'U', n, 0, alpha1, d1, 1, (void *) z, 1, one,  y, incy);
    suzerain_blasext_dsbmzv_external(
            'U', n, 0, alpha2, d2, 1, (void *) z, 1, one,  y, incy);
    suzerain_blas_free(z);
}

void
suzerain_blasext_dgbdddmzv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha0[2],
        const double *d0,
        const double alpha1[2],
        const double *d1,
        const double alpha2[2],
        const double *d2,
        const double *a,
        const int lda,
        const double (*x)[2],
        const int incx,
        const double beta[2],
        double (*y)[2],
        const int incy)
{
    const int info = suzerain_gbdddmv_dz(trans, n, kl, ku,
                                         alpha0, d0, alpha1, d1, alpha2, d2,
                                         a, lda, x, incx,
                                         beta,   y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_sgbiddmv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha0,
        const float alpha1,
        const float *d1,
        const float alpha2,
        const float *d2,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    float * const z = suzerain_blas_malloc(n*sizeof(float));
    if (SUZERAIN_UNLIKELY(!z)) suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_sgbmv_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_saxpby(
            n, alpha0, z, 1, beta, y, incy);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha1, d1, 1, z, 1, 1, y, incy);
    suzerain_blas_ssbmv_external(
            'U', n, 0, alpha2, d2, 1, z, 1, 1, y, incy);
    suzerain_blas_free(z);
}

void
suzerain_blasext_sgbiddmv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha0,
        const float alpha1,
        const float *d1,
        const float alpha2,
        const float *d2,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    const int info = suzerain_gbiddmv_s(trans, n, kl, ku,
                                        alpha0, alpha1, d1, alpha2, d2,
                                        a, lda, x, incx,
                                        beta,   y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_dgbiddmv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha0,
        const double alpha1,
        const double *d1,
        const double alpha2,
        const double *d2,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    double * const z = suzerain_blas_malloc(n*sizeof(double));
    if (SUZERAIN_UNLIKELY(!z)) suzerain_blas_xerbla(__func__, -1);
    suzerain_blas_dgbmv_external(
            trans, n, n, kl, ku, 1, a, lda, x, incx, 0, z, 1);
    suzerain_blas_daxpby(
            n, alpha0, z, 1, beta, y, incy);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha1, d1, 1, z, 1, 1, y, incy);
    suzerain_blas_dsbmv_external(
            'U', n, 0, alpha2, d2, 1, z, 1, 1, y, incy);
    suzerain_blas_free(z);
}

void
suzerain_blasext_dgbiddmv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha0,
        const double alpha1,
        const double *d1,
        const double alpha2,
        const double *d2,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    const int info = suzerain_gbiddmv_d(trans, n, kl, ku,
                                        alpha0, alpha1, d1, alpha2, d2,
                                        a, lda, x, incx,
                                        beta,   y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_sgbiddmzv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha0[2],
        const float alpha1[2],
        const float *d1,
        const float alpha2[2],
        const float *d2,
        const float *a,
        const int lda,
        const float (*x)[2],
        const int incx,
        const float beta[2],
        float (*y)[2],
        const int incy)
{
    static const float one[2]  = { 1, 0 };
    static const float zero[2] = { 0, 0 };
    float (* const z)[2] = suzerain_blas_malloc(2*n*sizeof(float));
    if (SUZERAIN_UNLIKELY(!z)) suzerain_blas_xerbla(__func__, -1);
    suzerain_blasext_sgbmzv_external(
            trans, n, n, kl, ku, one, a, lda, x, incx, zero, z, 1);
    suzerain_blas_caxpby(
            n, alpha0, (void *) z, 1, beta, y, incy);
    suzerain_blasext_ssbmzv_external(
            'U', n, 0, alpha1, d1, 1, (void *) z, 1, one, y, incy);
    suzerain_blasext_ssbmzv_external(
            'U', n, 0, alpha2, d2, 1, (void *) z, 1, one, y, incy);
    suzerain_blas_free(z);
}

void
suzerain_blasext_sgbiddmzv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha0[2],
        const float alpha1[2],
        const float *d1,
        const float alpha2[2],
        const float *d2,
        const float *a,
        const int lda,
        const float (*x)[2],
        const int incx,
        const float beta[2],
        float (*y)[2],
        const int incy)
{
    const int info = suzerain_gbiddmv_sc(trans, n, kl, ku,
                                         alpha0, alpha1, d1, alpha2, d2,
                                         a, lda, x, incx,
                                         beta,   y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_dgbiddmzv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha0[2],
        const double alpha1[2],
        const double *d1,
        const double alpha2[2],
        const double *d2,
        const double *a,
        const int lda,
        const double (*x)[2],
        const int incx,
        const double beta[2],
        double (*y)[2],
        const int incy)
{
    static const double one[2]  = { 1, 0 };
    static const double zero[2] = { 0, 0 };
    double (* const z)[2] = suzerain_blas_malloc(2*n*sizeof(double));
    if (SUZERAIN_UNLIKELY(!z)) suzerain_blas_xerbla(__func__, -1);
    suzerain_blasext_dgbmzv_external(
            trans, n, n, kl, ku, one, a, lda, x, incx, zero, z, 1);
    suzerain_blas_zaxpby(
            n, alpha0, (void *) z, 1, beta, y, incy);
    suzerain_blasext_dsbmzv_external(
            'U', n, 0, alpha1, d1, 1, (void *) z, 1, one, y, incy);
    suzerain_blasext_dsbmzv_external(
            'U', n, 0, alpha2, d2, 1, (void *) z, 1, one, y, incy);
    suzerain_blas_free(z);
}

void
suzerain_blasext_dgbiddmzv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha0[2],
        const double alpha1[2],
        const double *d1,
        const double alpha2[2],
        const double *d2,
        const double *a,
        const int lda,
        const double (*x)[2],
        const int incx,
        const double beta[2],
        double (*y)[2],
        const int incy)
{
    const int info = suzerain_gbiddmv_dz(trans, n, kl, ku,
                                         alpha0, alpha1, d1, alpha2, d2,
                                         a, lda, x, incx,
                                         beta,   y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_ssbmzv(
        const char uplo,
        const int n,
        const int k,
        const float alpha[2],
        const float *a,
        const int lda,
        const float (*x)[2],
        const int incx,
        const float beta[2],
        float (*y)[2],
        const int incy)
{
    const int info = suzerain_sbmv_sc(uplo, n, k,
                                      alpha, a, lda, x, incx,
                                      beta,          y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_dsbmzv(
        const char uplo,
        const int n,
        const int k,
        const double alpha[2],
        const double *a,
        const int lda,
        const double (*x)[2],
        const int incx,
        const double beta[2],
        double (*y)[2],
        const int incy)
{
    const int info = suzerain_sbmv_dz(uplo, n, k,
                                      alpha, a, lda, x, incx,
                                      beta,          y, incy);
    if (SUZERAIN_UNLIKELY(info)) suzerain_blas_xerbla(__func__, info);
}

void
suzerain_blasext_sge_diag_scale_acc(
        const int m,
        const int n,
        const float alpha,
        const float *a,
        const int inca,
        const int lda,
        const float *d,
        const int ldd,
        const float beta,
        float *b,
        const int incb,
        const int ldb)
{
    if (SUZERAIN_UNLIKELY(inca < 1)) suzerain_blas_xerbla(__func__,  5);
    if (SUZERAIN_UNLIKELY(ldd  < 0)) suzerain_blas_xerbla(__func__,  8);
    if (SUZERAIN_UNLIKELY(incb < 1)) suzerain_blas_xerbla(__func__, 11);

#pragma warning(push,disable:1572)
    const _Bool alpha_is_zero = (alpha == 0.0f);
    const _Bool beta_is_one   = (beta  == 1.0f);
#pragma warning(pop)

    if (SUZERAIN_UNLIKELY((alpha_is_zero && beta_is_one) || m <= 0 || n <= 0))
        return;

#define COMPUTE_SCALED_DIAGONAL_AS(x) \
    const float tmp = alpha*(*d)

    const float * const b_end = b + n*ldb;
    if (beta_is_one) {
        for (/* NOP */; b < b_end; a += lda, d += ldd, b += ldb) {
            COMPUTE_SCALED_DIAGONAL_AS(tmp);
            suzerain_blas_saxpy( m, tmp, a, inca,       b, incb);
        }
    } else {
        for (/* NOP */; b < b_end; a += lda, d += ldd, b += ldb) {
            COMPUTE_SCALED_DIAGONAL_AS(tmp);
            suzerain_blas_saxpby(m, tmp, a, inca, beta, b, incb);
        }
    }

#undef COMPUTE_SCALED_DIAGONAL_AS
}

void
suzerain_blasext_dge_diag_scale_acc(
        const int m,
        const int n,
        const double alpha,
        const double *a,
        const int inca,
        const int lda,
        const double *d,
        const int ldd,
        const double beta,
        double *b,
        const int incb,
        const int ldb)
{
    if (SUZERAIN_UNLIKELY(inca < 1)) suzerain_blas_xerbla(__func__,  5);
    if (SUZERAIN_UNLIKELY(ldd  < 0)) suzerain_blas_xerbla(__func__,  8);
    if (SUZERAIN_UNLIKELY(incb < 1)) suzerain_blas_xerbla(__func__, 11);

#pragma warning(push,disable:1572)
    const _Bool alpha_is_zero = (alpha == 0.0);
    const _Bool beta_is_one   = (beta  == 1.0);
#pragma warning(pop)

    if (SUZERAIN_UNLIKELY((alpha_is_zero && beta_is_one) || m <= 0 || n <= 0))
        return;

#define COMPUTE_SCALED_DIAGONAL_AS(x) \
    const double tmp = alpha*(*d)

    const double * const b_end = b + n*ldb;
    if (beta_is_one) {
        for (/* NOP */; b < b_end; a += lda, d += ldd, b += ldb) {
            COMPUTE_SCALED_DIAGONAL_AS(tmp);
            suzerain_blas_daxpy( m, tmp, a, inca,       b, incb);
        }
    } else {
        for (/* NOP */; b < b_end; a += lda, d += ldd, b += ldb) {
            COMPUTE_SCALED_DIAGONAL_AS(tmp);
            suzerain_blas_daxpby(m, tmp, a, inca, beta, b, incb);
        }
    }

#undef COMPUTE_SCALED_DIAGONAL_AS
}

void
suzerain_blasext_cge_diag_scale_acc(
        const int m,
        const int n,
        const float alpha[2],
        const float (*a)[2],
        const int inca,
        const int lda,
        const float (*d)[2],
        const int ldd,
        const float beta[2],
        float (*b)[2],
        const int incb,
        const int ldb)
{
    if (SUZERAIN_UNLIKELY(inca < 1)) suzerain_blas_xerbla(__func__,  5);
    if (SUZERAIN_UNLIKELY(ldd  < 0)) suzerain_blas_xerbla(__func__,  8);
    if (SUZERAIN_UNLIKELY(incb < 1)) suzerain_blas_xerbla(__func__, 11);

#pragma warning(push,disable:1572)
    const _Bool alpha_is_zero = (alpha[0] == 0.0f && alpha[1] == 0.0f);
    const _Bool beta_is_one   = (beta[0]  == 1.0f && beta[1]  == 0.0f);
#pragma warning(pop)

    if (SUZERAIN_UNLIKELY((alpha_is_zero && beta_is_one) || m <= 0 || n <= 0))
        return;

#define COMPUTE_SCALED_DIAGONAL_AS(x)                               \
    const float tmp[2] = {   alpha[0]*(*d)[0] - alpha[1]*(*d)[1],   \
                             alpha[0]*(*d)[1] + alpha[1]*(*d)[0] };

    /* const */ float (* const b_end)[2] = b + n*ldb;
    if (beta_is_one) {
        for (/* NOP */; b < b_end; a += lda, d += ldd, b += ldb) {
            COMPUTE_SCALED_DIAGONAL_AS(tmp);
            suzerain_blas_caxpy( m, tmp, a, inca,       b, incb);
        }
    } else {
        for (/* NOP */; b < b_end; a += lda, d += ldd, b += ldb) {
            COMPUTE_SCALED_DIAGONAL_AS(tmp);
            suzerain_blas_caxpby(m, tmp, a, inca, beta, b, incb);
        }
    }

#undef COMPUTE_SCALED_DIAGONAL_AS
}

void
suzerain_blasext_zge_diag_scale_acc(
        const int m,
        const int n,
        const double alpha[2],
        const double (*a)[2],
        const int inca,
        const int lda,
        const double (*d)[2],
        const int ldd,
        const double beta[2],
        double (*b)[2],
        const int incb,
        const int ldb)
{
    if (SUZERAIN_UNLIKELY(inca < 1)) suzerain_blas_xerbla(__func__,  5);
    if (SUZERAIN_UNLIKELY(ldd  < 0)) suzerain_blas_xerbla(__func__,  8);
    if (SUZERAIN_UNLIKELY(incb < 1)) suzerain_blas_xerbla(__func__, 11);

#pragma warning(push,disable:1572)
    const _Bool alpha_is_zero = (alpha[0] == 0.0 && alpha[1] == 0.0);
    const _Bool beta_is_one   = (beta[0]  == 1.0 && beta[1]  == 0.0);
#pragma warning(pop)

    if (SUZERAIN_UNLIKELY((alpha_is_zero && beta_is_one) || m <= 0 || n <= 0))
        return;

#define COMPUTE_SCALED_DIAGONAL_AS(x)                                \
    const double tmp[2] = {   alpha[0]*(*d)[0] - alpha[1]*(*d)[1],   \
                              alpha[0]*(*d)[1] + alpha[1]*(*d)[0] };

    /* const */ double (* const b_end)[2] = b + n*ldb;
    if (beta_is_one) {
        for (/* NOP */; b < b_end; a += lda, d += ldd, b += ldb) {
            COMPUTE_SCALED_DIAGONAL_AS(tmp);
            suzerain_blas_zaxpy( m, tmp, a, inca,       b, incb);
        }
    } else {
        for (/* NOP */; b < b_end; a += lda, d += ldd, b += ldb) {
            COMPUTE_SCALED_DIAGONAL_AS(tmp);
            suzerain_blas_zaxpby(m, tmp, a, inca, beta, b, incb);
        }
    }

#undef COMPUTE_SCALED_DIAGONAL_AS
}

void
suzerain_blasext_zge_diag_scale_dacc(
        const int m,
        const int n,
        const double alpha[2],
        const double *a,
        const int inca,
        const int lda,
        const double *d,
        const int ldd,
        const double beta[2],
        double (*b)[2],
        const int incb,
        const int ldb)
{
    if (SUZERAIN_UNLIKELY(inca < 1)) suzerain_blas_xerbla(__func__,  5);
    if (SUZERAIN_UNLIKELY(ldd  < 0)) suzerain_blas_xerbla(__func__,  8);
    if (SUZERAIN_UNLIKELY(incb < 1)) suzerain_blas_xerbla(__func__, 11);

#pragma warning(push,disable:1572)
    const _Bool alpha_is_zero = (alpha[0] == 0.0 && alpha[1] == 0.0);
    const _Bool beta_is_one   = (beta[0]  == 1.0 && beta[1]  == 0.0);
#pragma warning(pop)

    if (SUZERAIN_UNLIKELY((alpha_is_zero && beta_is_one) || m <= 0 || n <= 0))
        return;

#define COMPUTE_SCALED_DIAGONAL_AS(x)                         \
    const double tmp[2] = {   alpha[0]*(*d), alpha[1]*(*d) };

    /* const */ double (* const b_end)[2] = b + n*ldb;
    if (beta_is_one) {
        for (/* NOP */; b < b_end; a += lda, d += ldd, b += ldb) {
            COMPUTE_SCALED_DIAGONAL_AS(tmp);
            suzerain_blasext_daxpzy( m, tmp, a, inca,       b, incb);
        }
    } else {
        for (/* NOP */; b < b_end; a += lda, d += ldd, b += ldb) {
            COMPUTE_SCALED_DIAGONAL_AS(tmp);
            suzerain_blasext_daxpzby(m, tmp, a, inca, beta, b, incb);
        }
    }

#undef COMPUTE_SCALED_DIAGONAL_AS
}

void
suzerain_blasext_zge_ddiag_scale_dacc(
        const int m,
        const int n,
        const double *a,
        const int inca,
        const int lda,
        const double alpha0[2],
        const double *d0,
        const int ldd0,
        const double alpha1[2],
        const double *d1,
        const int ldd1,
        const double beta[2],
        double (*b)[2],
        const int incb,
        const int ldb)
{
    if (SUZERAIN_UNLIKELY(inca < 1)) suzerain_blas_xerbla(__func__,  4);
    if (SUZERAIN_UNLIKELY(ldd0 < 0)) suzerain_blas_xerbla(__func__,  8);
    if (SUZERAIN_UNLIKELY(ldd1 < 0)) suzerain_blas_xerbla(__func__, 11);
    if (SUZERAIN_UNLIKELY(incb < 1)) suzerain_blas_xerbla(__func__, 14);

#pragma warning(push,disable:1572)
    const _Bool alpha0_is_zero = (alpha0[0] == 0.0 && alpha0[1] == 0.0);
    const _Bool alpha1_is_zero = (alpha1[0] == 0.0 && alpha1[1] == 0.0);
    const _Bool beta_is_one    = (beta[0]   == 1.0 && beta[1]   == 0.0);
#pragma warning(pop)

    if (SUZERAIN_UNLIKELY((   alpha0_is_zero && alpha1_is_zero && beta_is_one)
                           || m <= 0
                           || n <= 0))
        return;

#define COMPUTE_SCALED_DIAGONAL_AS(x)                              \
    const double tmp[2] = {   alpha0[0]*(*d0) + alpha1[0]*(*d1),   \
                              alpha0[1]*(*d0) + alpha1[1]*(*d1) };

    /* const */ double (* const b_end)[2] = b + n*ldb;
    if (beta_is_one) {
        for (/* NOP */; b < b_end; a += lda, d0 += ldd0, d1 += ldd1, b += ldb) {
            COMPUTE_SCALED_DIAGONAL_AS(tmp);
            suzerain_blasext_daxpzy( m, tmp, a, inca,       b, incb);
        }
    } else {
        for (/* NOP */; b < b_end; a += lda, d0 += ldd0, d1 += ldd1, b += ldb) {
            COMPUTE_SCALED_DIAGONAL_AS(tmp);
            suzerain_blasext_daxpzby(m, tmp, a, inca, beta, b, incb);
        }
    }

#undef COMPUTE_SCALED_DIAGONAL_AS
}

void
suzerain_blasext_zge_dddiag_scale_dacc(
        const int m,
        const int n,
        const double *a,
        const int inca,
        const int lda,
        const double alpha0[2],
        const double *d0,
        const int ldd0,
        const double alpha1[2],
        const double *d1,
        const int ldd1,
        const double alpha2[2],
        const double *d2,
        const int ldd2,
        const double beta[2],
        double (*b)[2],
        const int incb,
        const int ldb)
{
    if (SUZERAIN_UNLIKELY(inca < 1)) suzerain_blas_xerbla(__func__,  4);
    if (SUZERAIN_UNLIKELY(ldd0 < 0)) suzerain_blas_xerbla(__func__,  8);
    if (SUZERAIN_UNLIKELY(ldd1 < 0)) suzerain_blas_xerbla(__func__, 11);
    if (SUZERAIN_UNLIKELY(ldd2 < 0)) suzerain_blas_xerbla(__func__, 14);
    if (SUZERAIN_UNLIKELY(incb < 1)) suzerain_blas_xerbla(__func__, 17);

#pragma warning(push,disable:1572)
    const _Bool alpha0_is_zero = (alpha0[0] == 0.0 && alpha0[1] == 0.0);
    const _Bool alpha1_is_zero = (alpha1[0] == 0.0 && alpha1[1] == 0.0);
    const _Bool alpha2_is_zero = (alpha2[0] == 0.0 && alpha2[1] == 0.0);
    const _Bool beta_is_one    = (beta[0]   == 1.0 && beta[1]   == 0.0);
#pragma warning(pop)

    if (SUZERAIN_UNLIKELY((   alpha0_is_zero && alpha1_is_zero
                                             && alpha2_is_zero && beta_is_one)
                           || m <= 0
                           || n <= 0))
        return;

#define COMPUTE_SCALED_DIAGONAL_AS(x)                            \
    const double tmp[2] = {                                      \
            alpha0[0]*(*d0) + alpha1[0]*(*d1) + alpha2[0]*(*d2), \
            alpha0[1]*(*d0) + alpha1[1]*(*d1) + alpha2[1]*(*d2)  \
    };

    /* const */ double (* const b_end)[2] = b + n*ldb;
    if (beta_is_one) {
        for (/* NOP */; b < b_end;
             a += lda, d0 += ldd0, d1 += ldd1, d2 += ldd2, b += ldb) {
            COMPUTE_SCALED_DIAGONAL_AS(tmp);
            suzerain_blasext_daxpzy( m, tmp, a, inca,       b, incb);
        }
    } else {
        for (/* NOP */; b < b_end;
             a += lda, d0 += ldd0, d1 += ldd1, d2 += ldd2, b += ldb) {
            COMPUTE_SCALED_DIAGONAL_AS(tmp);
            suzerain_blasext_daxpzby(m, tmp, a, inca, beta, b, incb);
        }
    }

#undef COMPUTE_SCALED_DIAGONAL_AS
}

void
suzerain_blasext_sgb_diag_scale_acc(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const float alpha,
        const float *a,
        const int inca,
        const int lda,
        const float *d,
        const int ldd,
        const float beta,
        float *b,
        const int incb,
        const int ldb)
{
    SUZERAIN_UNUSED(m);
    // Only partial checks as next routine covers remaining requirements
    if (SUZERAIN_UNLIKELY(kl < 0)) suzerain_blas_xerbla(__func__, 3);
    if (SUZERAIN_UNLIKELY(ku < 0)) suzerain_blas_xerbla(__func__, 4);
    return suzerain_blasext_sge_diag_scale_acc(ku + 1 + kl, n,
                                               alpha, a, inca, lda, d, ldd,
                                               beta,  b, incb, ldb);
}

void
suzerain_blasext_dgb_diag_scale_acc(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const double alpha,
        const double *a,
        const int inca,
        const int lda,
        const double *d,
        const int ldd,
        const double beta,
        double *b,
        const int incb,
        const int ldb)
{
    SUZERAIN_UNUSED(m);
    // Only partial checks as next routine covers remaining requirements
    if (SUZERAIN_UNLIKELY(kl < 0)) suzerain_blas_xerbla(__func__, 3);
    if (SUZERAIN_UNLIKELY(ku < 0)) suzerain_blas_xerbla(__func__, 4);
    return suzerain_blasext_dge_diag_scale_acc(ku + 1 + kl, n,
                                               alpha, a, inca, lda, d, ldd,
                                               beta,  b, incb, ldb);
}

void
suzerain_blasext_cgb_diag_scale_acc(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const float alpha[2],
        const float (*a)[2],
        const int inca,
        const int lda,
        const float (*d)[2],
        const int ldd,
        const float beta[2],
        float (*b)[2],
        const int incb,
        const int ldb)
{
    SUZERAIN_UNUSED(m);
    // Only partial checks as next routine covers remaining requirements
    if (SUZERAIN_UNLIKELY(kl < 0)) suzerain_blas_xerbla(__func__, 3);
    if (SUZERAIN_UNLIKELY(ku < 0)) suzerain_blas_xerbla(__func__, 4);
    return suzerain_blasext_cge_diag_scale_acc(ku + 1 + kl, n,
                                               alpha, a, inca, lda, d, ldd,
                                               beta,  b, incb, ldb);
}

void
suzerain_blasext_zgb_diag_scale_acc(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const double alpha[2],
        const double (*a)[2],
        const int inca,
        const int lda,
        const double (*d)[2],
        const int ldd,
        const double beta[2],
        double (*b)[2],
        const int incb,
        const int ldb)
{
    SUZERAIN_UNUSED(m);
    // Only partial checks as next routine covers remaining requirements
    if (SUZERAIN_UNLIKELY(kl < 0)) suzerain_blas_xerbla(__func__, 3);
    if (SUZERAIN_UNLIKELY(ku < 0)) suzerain_blas_xerbla(__func__, 4);
    return suzerain_blasext_zge_diag_scale_acc(ku + 1 + kl, n,
                                               alpha, a, inca, lda, d, ldd,
                                               beta,  b, incb, ldb);
}

void
suzerain_blasext_zgb_diag_scale_dacc(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const double alpha[2],
        const double *a,
        const int inca,
        const int lda,
        const double *d,
        const int ldd,
        const double beta[2],
        double (*b)[2],
        const int incb,
        const int ldb)
{
    SUZERAIN_UNUSED(m);
    // Only partial checks as next routine covers remaining requirements
    if (SUZERAIN_UNLIKELY(kl < 0)) suzerain_blas_xerbla(__func__, 3);
    if (SUZERAIN_UNLIKELY(ku < 0)) suzerain_blas_xerbla(__func__, 4);
    return suzerain_blasext_zge_diag_scale_dacc(ku + 1 + kl, n,
                                                alpha, a, inca, lda, d, ldd,
                                                beta,  b, incb, ldb);
}

void
suzerain_blasext_zgb_ddiag_scale_dacc(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const double *a,
        const int inca,
        const int lda,
        const double alpha0[2],
        const double *d0,
        const int ldd0,
        const double alpha1[2],
        const double *d1,
        const int ldd1,
        const double beta[2],
        double (*b)[2],
        const int incb,
        const int ldb)
{
    SUZERAIN_UNUSED(m);
    // Only partial checks as next routine covers remaining requirements
    if (SUZERAIN_UNLIKELY(kl < 0)) suzerain_blas_xerbla(__func__, 3);
    if (SUZERAIN_UNLIKELY(ku < 0)) suzerain_blas_xerbla(__func__, 4);
    return suzerain_blasext_zge_ddiag_scale_dacc(ku + 1 + kl, n,
                                                 a, inca, lda,
                                                 alpha0, d0, ldd0,
                                                 alpha1, d1, ldd1,
                                                 beta,    b, incb, ldb);
}

void
suzerain_blasext_zgb_dddiag_scale_dacc(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const double *a,
        const int inca,
        const int lda,
        const double alpha0[2],
        const double *d0,
        const int ldd0,
        const double alpha1[2],
        const double *d1,
        const int ldd1,
        const double alpha2[2],
        const double *d2,
        const int ldd2,
        const double beta[2],
        double (*b)[2],
        const int incb,
        const int ldb)
{
    SUZERAIN_UNUSED(m);
    // Only partial checks as next routine covers remaining requirements
    if (SUZERAIN_UNLIKELY(kl < 0)) suzerain_blas_xerbla(__func__, 3);
    if (SUZERAIN_UNLIKELY(ku < 0)) suzerain_blas_xerbla(__func__, 4);
    return suzerain_blasext_zge_dddiag_scale_dacc(ku + 1 + kl, n,
                                                  a, inca, lda,
                                                  alpha0, d0, ldd0,
                                                  alpha1, d1, ldd1,
                                                  alpha2, d2, ldd2,
                                                  beta,    b, incb, ldb);
}

void
suzerain_blasext_zgb_dacc(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const double alpha[2],
        const double *a,
        const int lda,
        const double beta[2],
        double (*b)[2],
        const int ldb)
{
    const double one = 1.0;
    suzerain_blasext_zgb_diag_scale_dacc(m, n, kl, ku,
                                         alpha, a, 1, lda, &one, 0,
                                         beta,  b, 1, ldb);
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

    if (m  < 0)        return -1;
    if (n  < 0)        return -2;
    if (kl < 0)        return -3;
    if (ku < 0)        return -4;
    if (!a)            return -5;
    if (lda < kl+ku+1) return -6;
    if (!norm1)        return -7;

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

    if (m  < 0)        return -1;
    if (n  < 0)        return -2;
    if (kl < 0)        return -3;
    if (ku < 0)        return -4;
    if (!a)            return -5;
    if (lda < kl+ku+1) return -6;
    if (!norm1)        return -7;

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
        const float (*a)[2],
        const int lda,
        float *norm1)
{
    // Defer to LAPACK for the square case
    if (m == n) {
        *norm1 = suzerain_lapack_clangb('1', n, kl, ku, a, lda, NULL);
        return 0;
    }

    if (m  < 0)        return -1;
    if (n  < 0)        return -2;
    if (kl < 0)        return -3;
    if (ku < 0)        return -4;
    if (!a)            return -5;
    if (lda < kl+ku+1) return -6;
    if (!norm1)        return -7;

    *norm1 = 0;
    for (int j = 0; j < n; ++j) {
        const float (*a_j)[2]  = a + j *lda;

        float     s    = 0;
        const int ibgn = imax(0, ku - j);
        const int iend = ku + imin(kl + 1, m - j);
        for (int i = ibgn; i < iend; ++i) {
            s += sqrtf(a_j[i][0]*a_j[i][0] + a_j[i][1]*a_j[i][1]);
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
        const double (*a)[2],
        const int lda,
        double *norm1)
{
    // Defer to LAPACK for the square case
    if (m == n) {
        *norm1 = suzerain_lapack_zlangb('1', n, kl, ku, a, lda, NULL);
        return 0;
    }

    if (m  < 0)        return -1;
    if (n  < 0)        return -2;
    if (kl < 0)        return -3;
    if (ku < 0)        return -4;
    if (!a)            return -5;
    if (lda < kl+ku+1) return -6;
    if (!norm1)        return -7;

    *norm1 = 0;
    for (int j = 0; j < n; ++j) {
        const double (*a_j)[2]  = a + j *lda;

        double    s    = 0;
        const int ibgn = imax(0, ku - j);
        const int iend = ku + imin(kl + 1, m - j);
        for (int i = ibgn; i < iend; ++i) {
            s += sqrt(a_j[i][0]*a_j[i][0] + a_j[i][1]*a_j[i][1]);
        }
        *norm1 = fmax(s, *norm1);
    }

    return 0;
}
