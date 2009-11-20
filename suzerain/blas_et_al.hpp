/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2009 The PECOS Development Team
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
 * blas_et_al.hpp: C++ wrappers around blas_et_al.h's functionality
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_BLAS_ET_AL_HPP__
#define __SUZERAIN_BLAS_ET_AL_HPP__

#include <suzerain/common.h>
#include <suzerain/blas_et_al.h>

/** @file
 * C++ wrappers for external BLAS and LAPACK routines necessary for Suzerain.
 * Provides function overloading atop the pure C interfaces in blas_et_al.h.
 */
namespace suzerain {

/**
 * Provides function templates for BLAS level 1, 2, and 3 functionality.  Also
 * includes functionality from the <a
 * href="http://www.netlib.org/blas/blast-forum/">BLAS Technical Forum
 * Standard</a>.
 */
namespace blas {

/*! \name Memory allocation
 * @{
 */

/*! @copydoc suzerain_blas_malloc(size_t) */
inline void * malloc(size_t size)
{
    return suzerain_blas_malloc(size);
}

/*! @copydoc suzerain_blas_calloc(size_t,size_t) */
inline void * calloc(size_t nmemb, size_t size)
{
    return suzerain_blas_calloc(nmemb, size);
}

/*! @} */

/*! \name BLAS level 1 operations
 * @{
 */

/*! @copydoc suzerain_blas_scopy */
template< typename FPT > void copy(
        const int n,
        const FPT *x,
        const int incx,
        FPT *y,
        const int incy);

/*! @copydoc suzerain_blas_scopy */
template<> inline void copy<float>(
        const int n,
        const float *x,
        const int incx,
        float *y,
        const int incy)
{
    return suzerain_blas_scopy(n, x, incx, y, incy);
}

/*! @copydoc suzerain_blas_scopy */
template<> inline void copy<double>(
        const int n,
        const double *x,
        const int incx,
        double *y,
        const int incy)
{
    return suzerain_blas_dcopy(n, x, incx, y, incy);
}

/*! @copydoc suzerain_blas_sdot */
template< typename FPT > FPT dot(
        const int n,
        const FPT *x,
        const int incx,
        const FPT *y,
        const int incy);

/*! @copydoc suzerain_blas_sdot */
template<> inline float dot<float>(
        const int n,
        const float *x,
        const int incx,
        const float *y,
        const int incy)
{
    return suzerain_blas_sdot(n, x, incx, y, incy);
}

/*! @copydoc suzerain_blas_sdot */
template<> inline double dot<double>(
        const int n,
        const double *x,
        const int incx,
        const double *y,
        const int incy)
{
    return suzerain_blas_ddot(n, x, incx, y, incy);
}

/*! @copydoc suzerain_blas_sasum */
template< typename FPT > FPT asum(
        const int n,
        const FPT *x,
        const int incx);

/*! @copydoc suzerain_blas_sasum */
template<> inline float asum<float>(
        const int n,
        const float *x,
        const int incx)
{
    return suzerain_blas_sasum(n, x, incx);
}

/*! @copydoc suzerain_blas_sasum */
template<> inline double asum<double>(
        const int n,
        const double *x,
        const int incx)
{
    return suzerain_blas_dasum(n, x, incx);
}

/*! @copydoc suzerain_blas_saxpy */
template< typename FPT > void axpy(
        const int n,
        const FPT alpha,
        const FPT *x,
        const int incx,
        FPT *y,
        const int incy);

/*! @copydoc suzerain_blas_saxpy */
template<> inline void axpy<float>(
        const int n,
        const float alpha,
        const float *x,
        const int incx,
        float *y,
        const int incy)
{
    return suzerain_blas_saxpy(n, alpha, x, incx, y, incy);
}

/*! @copydoc suzerain_blas_saxpy */
template<> inline void axpy<double>(
        const int n,
        const double alpha,
        const double *x,
        const int incx,
        double *y,
        const int incy)
{
    return suzerain_blas_daxpy(n, alpha, x, incx, y, incy);
}

/*! @copydoc suzerain_blas_saxpby */
template< typename FPT > void axpby(
        const int n,
        const FPT alpha,
        const FPT *x,
        const int incx,
        const FPT beta,
        FPT *y,
        const int incy);

/*! @copydoc suzerain_blas_saxpby */
template<> inline void axpby<float>(
        const int n,
        const float alpha,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    return suzerain_blas_saxpby(n, alpha, x, incx, beta, y, incy);
}

/*! @copydoc suzerain_blas_saxpby */
template<> inline void axpby<double>(
        const int n,
        const double alpha,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    return suzerain_blas_daxpby(n, alpha, x, incx, beta, y, incy);
}

/*! @copydoc suzerain_blas_swaxpby */
template< typename FPT > void waxpby(
        const int n,
        const FPT alpha,
        const FPT *x,
        const int incx,
        const FPT beta,
        const FPT *y,
        const int incy,
        FPT *w,
        const int incw);

/*! @copydoc suzerain_blas_swaxpby */
template<> inline void waxpby<float>(
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
    return suzerain_blas_swaxpby(n, alpha, x, incx, beta, y, incy, w, incw);
}

/*! @copydoc suzerain_blas_swaxpby */
template<> inline void waxpby<double>(
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
    return suzerain_blas_dwaxpby(n, alpha, x, incx, beta, y, incy, w, incw);
}

/*! @} */

/*! \name BLAS level 2 operations
 * @{
 */

/*! @copydoc suzerain_blas_sgbmv */
template< typename FPT > void gbmv(
        const char trans,
        const int m,
        const int n,
        const int kl,
        const int ku,
        const FPT alpha,
        const FPT *a,
        const int lda,
        const FPT *x,
        const int incx,
        const FPT beta,
        FPT *y,
        const int incy);

/*! @copydoc suzerain_blas_sgbmv */
template<> inline void gbmv<float>(
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
    return suzerain_blas_sgbmv(
        trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
}

/*! @copydoc suzerain_blas_sgbmv */
template<> inline void gbmv<double>(
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
    return suzerain_blas_dgbmv(
        trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
}

/*! @copydoc suzerain_blas_sgb_acc */
template< typename FPT > void gb_acc(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const FPT alpha,
        const FPT *a,
        const int lda,
        const FPT beta,
        FPT *b,
        const int ldb);

/*! @copydoc suzerain_blas_sgb_acc */
template<> inline void gb_acc<float>(
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
    return suzerain_blas_sgb_acc(
        m, n, kl, ku, alpha, a, lda, beta, b, ldb);
}

/*! @copydoc suzerain_blas_sgb_acc */
template<> inline void gb_acc<double>(
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
    return suzerain_blas_dgb_acc(
        m, n, kl, ku, alpha, a, lda, beta, b, ldb);
}

/*! @} */

/*! \name BLAS level 3 operations
 * @{
 */

/*! @} */

} // namespace blas

/**
 * Provides function templates for LAPACK functionality.
 */
namespace lapack {

/*! @copydoc suzerain_lapack_sgbtrf */
template< typename FPT > int gbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        FPT *ab,
        const int ldab,
        int *ipiv);

/*! @copydoc suzerain_lapack_sgbtrf */
template<> inline int gbtrf<float>(
        const int m,
        const int n,
        const int kl,
        const int ku,
        float *ab,
        const int ldab,
        int *ipiv)
{
    return suzerain_lapack_sgbtrf(m, n, kl, ku, ab, ldab, ipiv);
}

/*! @copydoc suzerain_lapack_sgbtrf */
template<> inline int gbtrf<double>(
        const int m,
        const int n,
        const int kl,
        const int ku,
        double *ab,
        const int ldab,
        int *ipiv)
{
    return suzerain_lapack_dgbtrf(m, n, kl, ku, ab, ldab, ipiv);
}

/*! @copydoc suzerain_lapack_sgbtrs */
template< typename FPT > int gbtrs(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        const FPT *ab,
        const int ldab,
        const int *ipiv,
        FPT *b,
        const int ldb);

/*! @copydoc suzerain_lapack_sgbtrs */
template<> inline int gbtrs<float>(
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
    return suzerain_lapack_sgbtrs(
        trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb);
}

/*! @copydoc suzerain_lapack_sgbtrs */
template<> inline int gbtrs<double>(
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
    return suzerain_lapack_dgbtrs(
        trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb);
}

} // namespace lapack

} // namespace suzerain

#endif /* __SUZERAIN_BLAS_ET_AL_HPP__ */
