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
 * blas_et_al.h: wraps external implementations of BLAS, LAPACK, et al.
 * $Id$
 */

#ifndef __SUZERAIN_BLAS_ET_AL_H
#define __SUZERAIN_BLAS_ET_AL_H

#include <suzerain/common.h>
#include <suzerain/complex.h>

/** \file
 * Wraps external BLAS and LAPACK routines necessary for Suzerain.
 * Provided to insulate the library from potential variations in
 * type signatures as well as to consolidate all Fortran-from-C
 * parameter differences.
 */

#ifdef __cplusplus
extern "C" {
#endif

/*! \name Miscellaneous utilities
 * @{
 */

/*!
 * \brief Invokes the BLAS implementation's \c XERBLA routine.
 *
 * \param srname The name of the calling routine.
 * \param info   The position of the invalid parameter in the
 *               parameter list of the calling routine.
 *
 * \return The error code \c info.
 */
int
suzerain_blas_xerbla(const char *srname, const int info);

/*! @} */

/*! \name Aligned memory allocation and deallocation
 * @{
 */

/**
 * The alignment of memory per the underlying BLAS' recommendations
 * for performance and numerical stability.
 *
 * Align at 16-byte boundaries per MKL user guide section 8.
 */
#define SUZERAIN_BLAS_ALIGNMENT 16

/*!
 * \brief Allocates memory aligned according to SUZERAIN_BLAS_ALIGNMENT.
 *
 * The memory is not cleared.  Any returned memory must later be freed using \c
 * suzerain_blas_free().
 *
 * \param size Number of bytes to allocate.
 *
 * \return On success, return a pointer to the allocated memory.  On failure,
 * return a NULL pointer.
 */
void *
suzerain_blas_malloc(size_t size);

/*!
 * \brief Allocate and clear memory aligned according to
 * SUZERAIN_BLAS_ALIGNMENT.
 *
 * The memory is set to zero.  Any returned memory must later be freed using \c
 * suzerain_blas_free.
 *
 * \param nmemb Number of elements to allocate
 * \param size Size of each member in bytes.
 *
 * \return On success, return a pointer to the allocated memory.  On failure,
 * return a NULL pointer.
 */
void *
suzerain_blas_calloc(size_t nmemb, size_t size);

/*!
 * \brief Free memory previously allocated by suzerain_blas_malloc() or
 * suzerain_blas_calloc().
 *
 * \param ptr Pointer to the previously allocated memory, which may be NULL.
 */
void
suzerain_blas_free(void *ptr);

/*! @} */

/*! \name BLAS level 1 operations
 * @{
 */

/*!
 * \brief Perform \f$ x \leftrightarrow{} y \f$ using BLAS's swap.
 *
 * \param n Number of elements in \c x and \c y.
 * \param x Source vector.
 * \param incx Source vector stride.
 * \param y Target vector.
 * \param incy Target vector stride.
 *
 * \see A BLAS reference for more details.
 */
void
suzerain_blas_sswap(
        const int n,
        float *x,
        const int incx,
        float *y,
        const int incy);

/*! \copydoc suzerain_blas_sswap */
void
suzerain_blas_dswap(
        const int n,
        double *x,
        const int incx,
        double *y,
        const int incy);

/*! \copydoc suzerain_blas_sswap */
void
suzerain_blas_cswap(
        const int n,
        complex_float *x,
        const int incx,
        complex_float *y,
        const int incy);

/*! \copydoc suzerain_blas_sswap */
void
suzerain_blas_zswap(
        const int n,
        complex_double *x,
        const int incx,
        complex_double *y,
        const int incy);

/*!
 * \brief Perform \f$ y \leftarrow{} x \f$ using BLAS's copy.
 *
 * \param n Number of elements in \c x and \c y.
 * \param x Source vector.
 * \param incx Source vector stride.
 * \param y Target vector.
 * \param incy Target vector stride.
 *
 * \see A BLAS reference for more details.
 */
void
suzerain_blas_scopy(
        const int n,
        const float *x,
        const int incx,
        float *y,
        const int incy);

/*! \copydoc suzerain_blas_scopy */
void
suzerain_blas_dcopy(
        const int n,
        const double *x,
        const int incx,
        double *y,
        const int incy);

/*! \copydoc suzerain_blas_scopy */
void
suzerain_blas_ccopy(
        const int n,
        const complex_float *x,
        const int incx,
        complex_float *y,
        const int incy);

/*! \copydoc suzerain_blas_scopy */
void
suzerain_blas_zcopy(
        const int n,
        const complex_double *x,
        const int incx,
        complex_double *y,
        const int incy);

/*!
 * \brief Compute \f$ x \cdot{} y \f$ using BLAS's dot.
 *
 * \param n Number of elements in \c x and \c y.
 * \param x First source vector.
 * \param incx First source vector stride.
 * \param y Second source vector.
 * \param incy Second source vector stride.
 *
 * \return \f$ x \cdot{} y \f$.
 *
 * \see A BLAS reference for more details.
 */
float
suzerain_blas_sdot(
        const int n,
        const float *x,
        const int incx,
        const float *y,
        const int incy);

/*! \copydoc suzerain_blas_sdot */
double
suzerain_blas_ddot(
        const int n,
        const double *x,
        const int incx,
        const double *y,
        const int incy);

/*!
 * \brief Compute \f$ x \cdot{} y \f$ using BLAS's dotc.
 *
 * \param[in]  n Number of elements in \c x and \c y.
 * \param[in]  x First source vector.
 * \param[in]  incx First source vector stride.
 * \param[in]  y Second source vector.
 * \param[in]  incy Second source vector stride.
 * \param[out] dotc The complex-valued result \f$ x \cdot{} y \f$.
 *
 * \see A BLAS reference for more details.
 */
void
suzerain_blas_cdotc(
        const int n,
        const complex_float *x,
        const int incx,
        const complex_float *y,
        const int incy,
        complex_float *dotc);

/*! \copydoc suzerain_blas_cdotc */
void
suzerain_blas_zdotc(
        const int n,
        const complex_double *x,
        const int incx,
        const complex_double *y,
        const int incy,
        complex_double *dotc);

/*!
 * \brief Compute \f$ \left|\left| x \right|\right|_{2} \f$ using BLAS's nrm2.
 *
 * \param n Number of elements in \c x.
 * \param x Source vector.
 * \param incx Source vector stride.
 *
 * \return \f$ \left|\left| x \right|\right|_{2} \f$
 *
 * \see A BLAS reference for more details.
 */
float
suzerain_blas_snrm2(
        const int n,
        const float *x,
        const int incx);

/*! \copydoc suzerain_blas_snrm2 */
double
suzerain_blas_dnrm2(
        const int n,
        const double *x,
        const int incx);

/*! \copydoc suzerain_blas_snrm2 */
float
suzerain_blas_scnrm2(
        const int n,
        const complex_float *x,
        const int incx);

/*! \copydoc suzerain_blas_snrm2 */
double
suzerain_blas_dznrm2(
        const int n,
        const complex_double *x,
        const int incx);

/*!
 * \brief Compute \f$ \left|\left| x \right|\right|_{1} \f$ using BLAS's asum.
 *
 * \param n Number of elements in \c x.
 * \param x Source vector.
 * \param incx Source vector stride.
 *
 * \return \f$ \left|\left| x \right|\right|_{1} \f$
 *
 * \see A BLAS reference for more details.
 */
float
suzerain_blas_sasum(
        const int n,
        const float *x,
        const int incx);

/*! \copydoc suzerain_blas_sasum */
double
suzerain_blas_dasum(
        const int n,
        const double *x,
        const int incx);

/*! \copydoc suzerain_blas_sasum */
float
suzerain_blas_scasum(
        const int n,
        const complex_float *x,
        const int incx);

/*! \copydoc suzerain_blas_sasum */
double
suzerain_blas_dzasum(
        const int n,
        const complex_double *x,
        const int incx);

/*!
 * \brief Find <em>zero-indexed</em> \f$i\f$ such that
 * \f$i = \operatorname{arg\,max}_{0\leq{}i<n} \left|x_i\right|\f$ using
 * BLAS's amax.
 *
 * \param n Number of elements in \c x.
 * \param x Source vector.
 * \param incx Source vector stride.
 *
 * \return Zero-indexed \f$i\f$ on success.  On nonpositive \c n returns -1.
 *
 * \see A BLAS reference for more details.
 */
int
suzerain_blas_isamax(
        const int n,
        const float *x,
        const int incx);

/*! \copydoc suzerain_blas_isamax */
int
suzerain_blas_idamax(
        const int n,
        const double *x,
        const int incx);

/*!
 * \brief Find <em>zero-indexed</em> \f$i\f$ such that \f$i =
 * \operatorname{arg\,max}_{0\leq{}i<n} \left|\mbox{Re}\left(x_i\right)\right|
 * + \left|\mbox{Im}\left(x_i\right)\right|\f$ using BLAS's amax.
 *
 * \copydetails suzerain_blas_isamax
 */
int
suzerain_blas_icamax(
        const int n,
        const complex_float *x,
        const int incx);

/*! \copydoc suzerain_blas_icamax */
int
suzerain_blas_izamax(
        const int n,
        const complex_double *x,
        const int incx);

/*!
 * \brief Find <em>zero-indexed</em> \f$i\f$ such that
 * \f$i = \operatorname{arg\,min}_{0\leq{}i<n} \left|x_i\right|\f$ using
 * BLAS's amin.
 *
 * \param n Number of elements in \c x.
 * \param x Source vector.
 * \param incx Source vector stride.
 *
 * \return Zero-indexed \f$i\f$ on success.  On nonpositive \c n returns -1.
 *
 * \see A BLAS reference for more details.
 */
int
suzerain_blas_isamin(
        const int n,
        const float *x,
        const int incx);

/*! \copydoc suzerain_blas_isamin */
int
suzerain_blas_idamin(
        const int n,
        const double *x,
        const int incx);

/*!
 * \brief Find <em>zero-indexed</em> \f$i\f$ such that \f$i =
 * \operatorname{arg\,min}_{0\leq{}i<n} \left|\mbox{Re}\left(x_i\right)\right|
 * + \left|\mbox{Im}\left(x_i\right)\right|\f$ using BLAS's amin.
 *
 * \copydetails suzerain_blas_isamin
 */
int
suzerain_blas_icamin(
        const int n,
        const complex_float *x,
        const int incx);

/*! \copydoc suzerain_blas_icamin */
int
suzerain_blas_izamin(
        const int n,
        const complex_double *x,
        const int incx);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{}x + y \f$ using BLAS's axpy.
 *
 * \param n Number of elements in \c x and \c y.
 * \param alpha Multiplicative scalar \f$ \alpha \f$
 * \param x First source vector.
 * \param incx First source vector stride.
 * \param y Second source vector and target vector.
 * \param incy Second source vector and target vector stride.
 *
 * \see A BLAS reference for more details.
 */
void
suzerain_blas_saxpy(
        const int n,
        const float alpha,
        const float *x,
        const int incx,
        float *y,
        const int incy);

/*! \copydoc suzerain_blas_saxpy */
void
suzerain_blas_daxpy(
        const int n,
        const double alpha,
        const double *x,
        const int incx,
        double *y,
        const int incy);

/*! \copydoc suzerain_blas_saxpy */
void
suzerain_blas_caxpy(
        const int n,
        const complex_float alpha,
        const complex_float *x,
        const int incx,
        complex_float *y,
        const int incy);

/*! \copydoc suzerain_blas_saxpy */
void
suzerain_blas_zaxpy(
        const int n,
        const complex_double alpha,
        const complex_double *x,
        const int incx,
        complex_double *y,
        const int incy);

/*! \copydoc suzerain_blas_saxpy */
void
suzerain_blas_zaxpy_d(
        const int n,
        const complex_double alpha,
        const double *x,
        const int incx,
        complex_double *y,
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{}x + \beta{}y \f$
 * using BLAS's axpby.
 *
 * If axpby is not available, simulate it using scal and axpy.
 *
 * \param n Number of elements in \c x and \c y.
 * \param alpha Multiplicative scalar \f$ \alpha \f$
 * \param x First source vector.
 * \param incx First source vector stride.
 * \param beta Multiplicative scalar \f$ \beta \f$
 * \param y Second source vector and target vector.
 * \param incy Second source vector and target vector stride.
 *
 * \see A BLAS reference for more details.
 */
void
suzerain_blas_saxpby(
        const int n,
        const float alpha,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy);

/*! \copydoc suzerain_blas_saxpby */
void
suzerain_blas_daxpby(
        const int n,
        const double alpha,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy);

/*! \copydoc suzerain_blas_saxpby */
void
suzerain_blas_caxpby(
        const int n,
        const complex_float alpha,
        const complex_float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy);

/*! \copydoc suzerain_blas_saxpby */
void
suzerain_blas_zaxpby(
        const int n,
        const complex_double alpha,
        const complex_double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy);

/*! \copydoc suzerain_blas_saxpby */
void
suzerain_blas_zaxpby_d(
        const int n,
        const complex_double alpha,
        const double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy);

/*!
 * \brief Compute \f$ w \leftarrow{} \alpha{}x + \beta{}y \f$
 * using BLAS's waxpby.
 *
 * If waxpby is not available, simulate it using copy, scal, and axpy.
 *
 * \param n Number of elements in \c x and \c y.
 * \param alpha Multiplicative scalar \f$ \alpha \f$
 * \param x First source vector.
 * \param incx First source vector stride.
 * \param beta Multiplicative scalar \f$ \beta \f$
 * \param y Second source vector.
 * \param incy Second source vector.
 * \param w Target vector.
 * \param incw Target vector stride.
 *
 * \see A BLAS reference for more details.
 */
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
        const int incw);

/*! \copydoc suzerain_blas_swaxpby */
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
        const int incw);

/*!
 * \brief Compute \f$ x \leftarrow{} \alpha{}x \f$
 * using BLAS's scal.
 *
 * \param n Number of elements in \c x.
 * \param alpha Multiplicative scalar \f$ \alpha \f$
 * \param x Source and target vector.
 * \param incx Source and target vector stride.
 *
 * \see A BLAS reference for more details.
 */
void
suzerain_blas_sscal(
        const int n,
        const float alpha,
        float *x,
        const int incx);

/*! \copydoc suzerain_blas_sscal */
void
suzerain_blas_dscal(
        const int n,
        const double alpha,
        double *x,
        const int incx);

/*! \copydoc suzerain_blas_sscal */
void
suzerain_blas_cscal(
        const int n,
        const complex_float alpha,
        complex_float *x,
        const int incx);

/*! \copydoc suzerain_blas_sscal */
void
suzerain_blas_zscal(
        const int n,
        const complex_double alpha,
        complex_double *x,
        const int incx);

/*! @} */

/*! \name BLAS level 2 operations
 * @{
 */

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$
 * using an external BLAS's gbmv.
 *
 * Transposes of \f$ A \f$ can be taken using the \c trans parameter.
 *
 * \param trans One of 'N', 'T', or 'C' for no transpose, a transpose,
 *      or a conjugate transpose, respectively.
 * \param m Number of rows in matrix \c a.
 * \param n Number of columns in matrix \c a.
 * \param kl Number of subdiagonals in band storage of \c a.
 * \param ku Number of superdiagonals in band storage of \c a.
 * \param alpha Multiplicative scalar \f$ \alpha \f$.
 * \param a General band storage for matrix \f$ A \f$.
 * \param lda Leading dimension of \c a.
 * \param x Vector to be multiplied.
 * \param incx Stride of vector \c x.
 * \param beta Multiplicative scalar \f$ \beta \f$.
 * \param y Vector to be added to product and to contain result.
 * \param incy Stride of vector \c y.
 *
 * \see A BLAS reference for more details, especially for general
 *      band matrix storage requirements.
 */
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
        const int incy);

/*! \copydoc suzerain_blas_sgbmv_external */
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
        const int incy);

/*! \copydoc suzerain_blas_sgbmv_external */
void
suzerain_blas_cgbmv_s_external(
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
        const int incy);

/*! \copydoc suzerain_blas_sgbmv_external */
void
suzerain_blas_zgbmv_d_external(
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$.
 * The computation may use internal, optimized kernels for some fixed bandwidth
 * cases while deferring to the BLAS for general use.
 *
 * \copydetails suzerain_blas_sgbmv_external
 */
int
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
        const int incy);

/*! \copydoc suzerain_blas_sgbmv */
int
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
        const int incy);

/*! \copydoc suzerain_blas_sgbmv */
int
suzerain_blas_cgbmv_s(
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
        const int incy);

/*! \copydoc suzerain_blas_sgbmv */
int
suzerain_blas_zgbmv_d(
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$
 * for symmetric \f$ A \f$ using an external BLAS's sbmv.
 *
 * \param uplo Either 'U'/'u' or 'L'/'l' if the upper or lower triangular
 *      part of \f$A\f$ is supplied in \c a, respectively.
 * \param n Number of rows and columns in matrix \c a.
 * \param k Number of superdiagonals in band storage of \c a.
 * \param alpha Multiplicative scalar \f$ \alpha \f$.
 * \param a Symmetric band storage for matrix \f$ A \f$.
 * \param lda Leading dimension of \c a.
 * \param x Vector to be multiplied.
 * \param incx Stride of vector \c x.
 * \param beta Multiplicative scalar \f$ \beta \f$.
 * \param y Vector to be added to product and to contain result.
 * \param incy Stride of vector \c y.
 *
 * \see A BLAS reference for more details, especially for general
 *      band matrix storage requirements.
 */
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
        const int incy);

/*! \copydoc suzerain_blas_ssbmv_external */
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
        const int incy);

/*! \copydoc suzerain_blas_ssbmv_external */
void
suzerain_blas_csbmv_s_external(
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
        const int incy);

/*! \copydoc suzerain_blas_ssbmv_external */
void
suzerain_blas_zsbmv_d_external(
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$ for
 * banded, symmetric \f$ A \f$.  The computation may use internal, optimized
 * kernels for some fixed bandwidth cases while deferring to the BLAS for
 * general use.
 *
 * \copydetails suzerain_blas_ssbmv_external
 */
int
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
        const int incy);

/*! \copydoc suzerain_blas_ssbmv */
int
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
        const int incy);

/*! \copydoc suzerain_blas_ssbmv */
int
suzerain_blas_csbmv_s(
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
        const int incy);

/*! \copydoc suzerain_blas_csbmv_s */
int
suzerain_blas_zsbmv_d(
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
        const int incy);


/*!
 * \brief Compute \f$ B \leftarrow{} \alpha{}A + \beta{}B \f$ using
 * BLAS' gb_acc.
 *
 * Matrices \f$ A \f$ and \f$ B \f$ both have
 * band storage and must have the same shape and same number
 * of super- and subdiagonals.
 *
 * \param m Number of rows in matrices \f$ A \f$ and \f$ B \f$.
 * \param n Number of columns in matrices \f$ A \f$ and \f$ B \f$.
 * \param kl Number of subdiagonals in band storage of \c ab.
 * \param ku Number of superdiagonals in band storage of \c ab.
 * \param alpha Multiplicative scalar \f$ \alpha \f$
 * \param a General band storage of the matrix \f$ A \f$.
 * \param lda Leading dimension of \c a.
 * \param beta Multiplicative scalar \f$ \beta \f$
 * \param b General band storage of the matrix \f$ B \f$.
 * \param ldb Leading dimension of \c b.
 *
 * \see A BLAS reference for more details, especially for general
 *      band matrix storage requirements.
 */
int
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
        const int ldb);

/*! \copydoc suzerain_blas_sgb_acc */
int
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
        const int ldb);

/*! \copydoc suzerain_blas_sgb_acc */
int
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
        const int ldb);

/*! \copydoc suzerain_blas_sgb_acc */
int
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
        const int ldb);

/*! \copydoc suzerain_blas_sgb_acc */
int
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
        const int ldb);


/*! @} */

/*! \name BLAS level 3 operations
 * @{
 */

/*! @} */

/*! \name LAPACK operations
 * @{
 */

/*!
 * \brief Compute the LUP decomposition of a general banded matrix using
 * LAPACK's gbtrf.
 *
 * Stores the results back into the same matrix.  Note that the banded matrix
 * storage must have \c kl extra superdiagonals available to handle the
 * factorization fill in.
 *
 * \param m Number of rows in matrix \c ab.
 * \param n Number of columns in matrix \c ab.
 * \param kl Number of subdiagonals in band storage of \c ab.
 * \param ku Number of superdiagonals in band storage of \c ab.
 * \param ab General band storage of the matrix to factor.
 * \param ldab Leading dimension of \c ab.
 * \param ipiv Pivot matrix computed in the decomposition.
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 *
 * \see gbtrs for how to solve a linear system
 *      once you have decomposed the matrix.
 * \see A LAPACK reference for more details, especially for the
 *      <tt>lda >= 2*kl + ku + 1</tt> storage requirement and the resulting
 *      factored storage format.
 */
int
suzerain_lapack_sgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        float *ab,
        const int ldab,
        int *ipiv);

/*! \copydoc suzerain_lapack_sgbtrf */
int
suzerain_lapack_dgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        double *ab,
        const int ldab,
        int *ipiv);

/*! \copydoc suzerain_lapack_sgbtrf */
int
suzerain_lapack_cgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        complex_float *ab,
        const int ldab,
        int *ipiv);

/*! \copydoc suzerain_lapack_sgbtrf */
int
suzerain_lapack_zgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        complex_double *ab,
        const int ldab,
        int *ipiv);

/*!
 * \brief Solve \f$ AX = B \f$ using the previously LUP decomposed general band
 * matrix \f$ A \f$ and LAPACK's gbtrs.
 *
 * Transposes of \f$ A \f$ can be taken using the \c trans parameter.
 *
 * \param trans One of 'N', 'T', or 'C' for no transpose, a transpose,
 *      or a conjugate transpose, respectively.
 * \param n Number of rows and columns in matrix \c ab.
 * \param kl Number of subdiagonals in band storage of \c ab.
 * \param ku Number of superdiagonals in nonfactored matrix \c ab.
 *      Note this is \e not the number of superdiagonals in the storage
 *      format of \c ab, but rather the number of superdiagonals required
 *      to store the non-factored matrix \c ab.  This is odd.
 * \param nrhs Number of right hand sides, or columns, in \c b.
 * \param ab General band storage of the matrix to factor.
 * \param ldab Leading dimension of \c ab.
 * \param ipiv Pivot matrix already computed in the decomposition.
 * \param b Matrix \f$ B \f$ containing right hand sides on invocation and
 *      solutions on return.
 * \param ldb Leading dimension of matrix \c b.
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 *
 * \see gbtrf for how to decompose the matrix \f$ A \f$.
 * \see A LAPACK reference for more details.
 */
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
        const int ldb);

/*! \copydoc suzerain_lapack_sgbtrs */
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
        const int ldb);

/*! \copydoc suzerain_lapack_sgbtrs */
int
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
        const int ldb);

/*! \copydoc suzerain_lapack_sgbtrs */
int
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
        const int ldb);

/*!
 * \brief Estimate the reciprocal of the condition number of a
 * real-valued general band matrix \f$ A \f$ using LAPACK's gbcon.
 *
 * \param norm One of '0', '1', or 'I' for the 1-norm, the 1-norm,
 *      or the infinity-norm, respectively.
 * \param n Number of rows and columns in matrix \c ab.
 * \param kl Number of subdiagonals in band storage of \c ab.
 * \param ku Number of superdiagonals in nonfactored matrix \c ab.
 *      Note this is \e not the number of superdiagonals in the storage
 *      format of \c ab, but rather the number of superdiagonals required
 *      to store the non-factored matrix \c ab.  This is odd.
 * \param ab General band storage of the matrix to factor.
 * \param ldab Leading dimension of \c ab.
 * \param ipiv Pivot matrix already computed in the decomposition.
 * \param anorm The norm of the matrix according to the norm choice
 *      made in \c norm.
 * \param rcond The reciprocal of the condition number of the matrix
 *      \f$ A \f$ computed as <tt>1/(norm(A)*norm(inv(A)))</tt>.
 * \param work Work array of dimension <tt>3*n</tt>.
 * \param iwork Work array of dimension <tt>n</tt>.
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 *
 * \see gbtrf for how to decompose the matrix \f$ A \f$.
 * \see A LAPACK reference for more details.
 */
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
        int *iwork);

/*! \copydoc suzerain_lapack_sgbcon */
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
        int *iwork);

/*!
 * \brief Estimate the reciprocal of the condition number of a
 * complex-valued general band matrix \f$ A \f$ using LAPACK's gbcon.
 *
 * \param norm One of '0', '1', or 'I' for the 1-norm, the 1-norm,
 *      or the infinity-norm, respectively.
 * \param n Number of rows and columns in matrix \c ab.
 * \param kl Number of subdiagonals in band storage of \c ab.
 * \param ku Number of superdiagonals in nonfactored matrix \c ab.
 *      Note this is \e not the number of superdiagonals in the storage
 *      format of \c ab, but rather the number of superdiagonals required
 *      to store the non-factored matrix \c ab.  This is odd.
 * \param ab General band storage of the matrix to factor.
 * \param ldab Leading dimension of \c ab.
 * \param ipiv Pivot matrix already computed in the decomposition.
 * \param anorm The norm of the matrix according to the norm choice
 *      made in \c norm.
 * \param rcond The reciprocal of the condition number of the matrix
 *      \f$ A \f$ computed as <tt>1/(norm(A)*norm(inv(A)))</tt>.
 * \param work Work array of dimension <tt>2*n</tt>.
 * \param rwork Work array of dimension <tt>n</tt>.
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 *
 * \see gbtrf for how to decompose the matrix \f$ A \f$.
 * \see A LAPACK reference for more details.
 */
int
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
        float  *rwork);

/*! \copydoc suzerain_lapack_cgbcon */
int
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
        double *rwork);

/*!
 * \brief Compute the solution to a banded system of linear equations \f$ AX=B
 * \f$ using an in-place LU factorization.
 *
 * \param[in]     n
 * \param[in]     kl
 * \param[in]     ku
 * \param[in]     nrhs
 * \param[in,out] ab    Dimension (<tt>ldab</tt>,<tt>n</tt>)
 * \param[in]     ldab  Minimum <tt>kl + ku + 1</tt>
 * \param[in,out] ipiv  Dimension \c n
 * \param[in,out] b     Dimension (<tt>ldb</tt>,<tt>nrhs</tt>)
 * \param[in]     ldb   Minimum \c n
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 *
 * \see A LAPACK reference for the details.
 */
int
suzerain_lapack_sgbsv(
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        float *ab,
        const int ldab,
        int *ipiv,
        float *b,
        const int ldb);

/*! \copydoc suzerain_lapack_sgbsv */
int
suzerain_lapack_dgbsv(
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        double *ab,
        const int ldab,
        int *ipiv,
        double *b,
        const int ldb);

/*! \copydoc suzerain_lapack_sgbsv */
int
suzerain_lapack_cgbsv(
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        complex_float *ab,
        const int ldab,
        int *ipiv,
        complex_float *b,
        const int ldb);

/*! \copydoc suzerain_lapack_sgbsv */
int
suzerain_lapack_zgbsv(
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        complex_double *ab,
        const int ldab,
        int *ipiv,
        complex_double *b,
        const int ldb);

/*!
 * \brief Compute the solution to a banded system of linear equations \f$ AX=B
 * \f$ using an out-of-place LU factorization.  Error bounds on the solution
 * and a condition estimate are also provided.
 *
 * \param[in]     fact
 * \param[in]     trans
 * \param[in]     n
 * \param[in]     kl
 * \param[in]     ku
 * \param[in]     nrhs
 * \param[in,out] ab    Dimension (<tt>ldab</tt>,<tt>n</tt>)
 * \param[in]     ldab  Minimum <tt>kl + ku + 1</tt>
 * \param[in,out] afb   Dimension (<tt>ldfab</tt>,<tt>n</tt>)
 * \param[in]     ldafb Minimum <tt>2*kl + ku + 1</tt>
 * \param[in,out] ipiv  Dimension \c n
 * \param[in,out] equed
 * \param[in,out] r     Dimension \c n when fact == 'R' or 'B'
 * \param[in,out] c     Dimension \c n when fact == 'C' or 'B'
 * \param[in,out] b     Dimension (<tt>ldb</tt>,<tt>nrhs</tt>)
 * \param[in]     ldb   Minimum \c n
 * \param[in,out] x     Dimension (<tt>ldx</tt>,<tt>nrhs</tt>)
 * \param[in]     ldx   Minimum \c n
 * \param[out]    rcond
 * \param[out]    ferr  Dimension \c nrhs
 * \param[out]    berr  Dimension \c nrhs
 * \param[out]    work  Dimension <tt>3*n</tt>
 * \param[out]    iwork Dimension \c n
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 *
 * \see A LAPACK reference for the details.
 */
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
        int *iwork);

/*! \copydoc suzerain_lapack_sgbsvx */
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
        int *iwork);

/*!
 * \brief Compute the solution to a banded system of linear equations \f$ AX=B
 * \f$ using an out-of-place LU factorization.  Error bounds on the solution
 * and a condition estimate are also provided.
 *
 * \param[in]     fact
 * \param[in]     trans
 * \param[in]     n
 * \param[in]     kl
 * \param[in]     ku
 * \param[in]     nrhs
 * \param[in,out] ab    Dimension (<tt>ldab</tt>,<tt>n</tt>)
 * \param[in]     ldab  Minimum <tt>kl + ku + 1</tt>
 * \param[in,out] afb   Dimension (<tt>ldfab</tt>,<tt>n</tt>)
 * \param[in]     ldafb Minimum <tt>2*kl + ku + 1</tt>
 * \param[in,out] ipiv  Dimension \c n
 * \param[in,out] equed
 * \param[in,out] r     Dimension \c n when fact == 'R' or 'B'
 * \param[in,out] c     Dimension \c n when fact == 'C' or 'B'
 * \param[in,out] b     Dimension (<tt>ldb</tt>,<tt>nrhs</tt>)
 * \param[in]     ldb   Minimum \c n
 * \param[in,out] x     Dimension (<tt>ldx</tt>,<tt>nrhs</tt>)
 * \param[in]     ldx   Minimum \c n
 * \param[out]    rcond
 * \param[out]    ferr  Dimension \c nrhs
 * \param[out]    berr  Dimension \c nrhs
 * \param[out]    work  Dimension <tt>2*n</tt>
 * \param[out]    rwork Dimension \c n
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 *
 * \see A LAPACK reference for the details.
 */
int
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
        float *rwork);

/*! \copydoc suzerain_lapack_cgbsvx */
int
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
        double *rwork);

/*!
 * \brief Compute the one norm, or the Frobenius norm, or the infinity
 * norm, or the element of largest absolute value of square band matrix.
 *
 * \see A LAPACK reference for more details.
 */
float
suzerain_lapack_slangb(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const float *ab,
        const int ldab,
        float *work);

/*! \copydoc suzerain_lapack_slangb */
double
suzerain_lapack_dlangb(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const double *ab,
        const int ldab,
        double *work);

/*! \copydoc suzerain_lapack_slangb */
float
suzerain_lapack_clangb(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const complex_float *ab,
        const int ldab,
        float *work);

/*! \copydoc suzerain_lapack_slangb */
double
suzerain_lapack_zlangb(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const complex_double *ab,
        const int ldab,
        double *work);

/*! @} */

/*! \name BLAS-like extensions
 * These are extensions built atop the BLAS, on vendor-specific BLAS-like routines,
 * and/or on custom coded loops.
 *
 * Some of these extensions refer to interleaved versus split complex storage.
 * Interleaved complex storage stores the imaginary component immediately after
 * the real component in memory, for example \c fftw_complex or
 * <tt>std::complex<T></tt>.  Split complex storage stores the real and
 * imaginary components in wholly separate locations.
 *
 * \see The FFTW manual for further discussion on
 * <a href="http://www.fftw.org/fftw3_doc/Interleaved-and-split-arrays.html">
 * interleaved versus split complex storage</a>.
 * @{
 */

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} D A x + \beta{} y \f$
 * using an external BLAS.
 *
 * \copydetails suzerain_gbdmv_s
 * \return On error calls suzerain_blas_xerbla().
 */
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
        const int incy);

/*! \copydoc suzerain_blasext_sgbdmv_external */
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} D A x + \beta{} y \f$ for complex
 * \f$\alpha{}\f$, \f$x\f$, \f$\beta\f$, and \f$y\f$ but real-valued \f$D\f$
 * and \f$A\f$ using an external BLAS.
 *
 * \copydetails suzerain_blasext_sgbdmv_external
 */
int
suzerain_blasext_cgbdmv_s_external(
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbdmv_s_external */
int
suzerain_blasext_zgbdmv_d_external(
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} D A x + \beta{} y \f$.
 *
 * \copydetails suzerain_gbdmv_s
 * \return On error calls suzerain_blas_xerbla().
 */
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
        const int incy);

/*! \copydoc suzerain_blasext_sgbdmv */
int
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} D A x + \beta{} y \f$ for complex
 * \f$\alpha{}\f$, \f$x\f$, \f$\beta\f$, and \f$y\f$ but real-valued \f$D\f$
 * and \f$A\f$.
 *
 * \copydetails suzerain_blasext_sgbdmv
 */
int
suzerain_blasext_cgbdmv_s(
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbdmv_s */
int
suzerain_blasext_zgbdmv_d(
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \left(\alpha_{0} D_0 + \alpha_{1}
 * D_1\right) A x + \beta{} y \f$ using an external BLAS.
 *
 * \copydetails suzerain_gbddmv_s
 * \return On error calls suzerain_blas_xerbla().
 */
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
        const int incy);

/*! \copydoc suzerain_blasext_sgbddmv_external */
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \left(\alpha{0} D_0 + \alpha_{1}
 * D_1\right) A x + \beta{} y \f$ for complex \f$\alpha{}\f$, \f$x\f$,
 * \f$\beta\f$, and \f$y\f$ but real-valued \f$D_0\f$, \f$D_1\f$, and \f$A\f$
 * using an external BLAS.
 *
 * \copydetails suzerain_blasext_sgbddmv_external
 */
int
suzerain_blasext_cgbddmv_s_external(
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbddmv_s_external */
int
suzerain_blasext_zgbddmv_d_external(
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \left(\alpha_{0} D_0 + \alpha_{1}
 * D_1\right) A x + \beta{} y \f$.
 *
 * \copydetails suzerain_blasext_sgbddmv_external
 */
int
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
        const int incy);

/*! \copydoc suzerain_blasext_sgbddmv */
int
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \left(\alpha_{0} D_0 + \alpha_{1}
 * D_1\right) A x + \beta{} y \f$ for complex \f$\alpha{}\f$, \f$x\f$,
 * \f$\beta\f$, and \f$y\f$ but real-valued \f$D_0\f$, \f$D_1\f$, and \f$A\f$.
 *
 * \copydetails suzerain_blasext_sgbddmv
 */
int
suzerain_blasext_cgbddmv_s(
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbddmv_s */
int
suzerain_blasext_zgbddmv_d(
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \left(\alpha_{0} D_0 + \alpha_{1} D_1 +
 * \alpha_{2} D_2\right) A x + \beta{} y \f$ using an external BLAS.
 *
 * \copydetails suzerain_gbdddmv_s
 * \return On error calls suzerain_blas_xerbla().
 */
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
        const int incy);

/*! \copydoc suzerain_blasext_sgbdddmv_external */
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \left(\alpha{0} D_0 + \alpha_{1} D_1 +
 * \alpha_{2} D_2\right) A x + \beta{} y \f$ for complex \f$\alpha{}\f$,
 * \f$x\f$, \f$\beta\f$, and \f$y\f$ but real-valued \f$D_0\f$, \f$D_1\f$,
 * \f$D_2\f$, and \f$A\f$ using an external BLAS.
 *
 * \copydetails suzerain_blasext_sgbdddmv_external
 */
int
suzerain_blasext_cgbdddmv_s_external(
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbdddmv_s_external */
int
suzerain_blasext_zgbdddmv_d_external(
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \left(\alpha_{0} D_0 + \alpha_{1}
 * D_1 + \alpha_{2} D_2\right) A x + \beta{} y \f$.
 *
 * \copydetails suzerain_blasext_sgbdddmv_external
 */
int
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
        const int incy);

/*! \copydoc suzerain_blasext_sgbdddmv */
int
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \left(\alpha_{0} D_0 + \alpha_{1} D_1 +
 * \alpha_{2} D_2\right) A x + \beta{} y \f$ for complex \f$\alpha{}\f$,
 * \f$x\f$, \f$\beta\f$, and \f$y\f$ but real-valued \f$D_0\f$, \f$D_1\f$,
 * \f$D_2\f$, and \f$A\f$.
 *
 * \copydetails suzerain_blasext_sgbdddmv
 */
int
suzerain_blasext_cgbdddmv_s(
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbdddmv_s */
int
suzerain_blasext_zgbdddmv_d(
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \left(\alpha_{0} D_0 + \alpha_{1} D_1 +
 * \alpha_{2} D_2 + \alpha_{3} D_3 \right) A x + \beta{} y \f$ using an
 * external BLAS.
 *
 * \copydetails suzerain_gbddddmv_s
 * \return On error calls suzerain_blas_xerbla().
 */
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
        const int incy);

/*! \copydoc suzerain_blasext_sgbddddmv_external */
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \left(\alpha{0} D_0 + \alpha_{1} D_1 +
 * \alpha_{2} D_2 + \alpha_{3} D_{3} \right) A x + \beta{} y \f$ for complex
 * \f$\alpha{}\f$, \f$x\f$, \f$\beta\f$, and \f$y\f$ but real-valued \f$D_0\f$,
 * \f$D_1\f$, \f$D_2\f$, \f$D_3\f$, and \f$A\f$ using an external BLAS.
 *
 * \copydetails suzerain_blasext_sgbddddmv_external
 */
int
suzerain_blasext_cgbddddmv_s_external(
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbddddmv_s_external */
int
suzerain_blasext_zgbddddmv_d_external(
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \left(\alpha_{0} D_0 + \alpha_{1}
 * D_1 + \alpha_{2} D_2 + \alpha_{3} D_3 \right) A x + \beta{} y \f$.
 *
 * \copydetails suzerain_blasext_sgbddddmv_external
 */
int
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
        const int incy);

/*! \copydoc suzerain_blasext_sgbddddmv */
int
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \left(\alpha_{0} D_0 + \alpha_{1} D_1 +
 * \alpha_{2} D_2 + \alpha_{3} D_3 \right) A x + \beta{} y \f$ for complex
 * \f$\alpha{}\f$, \f$x\f$, \f$\beta\f$, and \f$y\f$ but real-valued \f$D_0\f$,
 * \f$D_1\f$, \f$D_2\f$, \f$D_3\f$, and \f$A\f$.
 *
 * \copydetails suzerain_blasext_sgbddddmv
 */
int
suzerain_blasext_cgbddddmv_s(
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbddddmv_s */
int
suzerain_blasext_zgbddddmv_d(
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \left(\alpha_{0} D_0 + \alpha_{1} D_1 +
 * \alpha_{2} D_2 + \alpha_{3} D_3 + \alpha_{4} D_4 \right) A x + \beta{} y \f$
 * using an external BLAS.
 *
 * \copydetails suzerain_gbdddddmv_s
 * \return On error calls suzerain_blas_xerbla().
 */
int
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
        const int incy);

/*! \copydoc suzerain_blasext_sgbdddddmv_external */
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \left(\alpha{0} D_0 + \alpha_{1} D_1 +
 * \alpha_{2} D_2 + \alpha_{3} D_{3} + \alpha_{4} D_{4} \right) A x + \beta{} y
 * \f$ for complex \f$\alpha{}\f$, \f$x\f$, \f$\beta\f$, and \f$y\f$ but
 * real-valued \f$D_0\f$, \f$D_1\f$, \f$D_2\f$, \f$D_3\f$, and \f$A\f$ using an
 * external BLAS.
 *
 * \copydetails suzerain_blasext_sgbdddddmv_external
 */
int
suzerain_blasext_cgbdddddmv_s_external(
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbdddddmv_s_external */
int
suzerain_blasext_zgbdddddmv_d_external(
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \left(\alpha_{0} D_0 + \alpha_{1} D_1 +
 * \alpha_{2} D_2 + \alpha_{3} D_3 + \alpha_{4} D_4 \right) A x + \beta{} y
 * \f$.
 *
 * \copydetails suzerain_blasext_sgbdddddmv_external
 */
int
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
        const int incy);

/*! \copydoc suzerain_blasext_sgbdddddmv */
int
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \left(\alpha_{0} D_0 + \alpha_{1} D_1 +
 * \alpha_{2} D_2 + \alpha_{3} D_3 + \alpha_{4} D_4 \right) A x + \beta{} y \f$
 * for complex \f$\alpha{}\f$, \f$x\f$, \f$\beta\f$, and \f$y\f$ but
 * real-valued \f$D_0\f$, \f$D_1\f$, \f$D_2\f$, \f$D_3\f$, \f$D_4\f$, and
 * \f$A\f$.
 *
 * \copydetails suzerain_blasext_sgbdddddmv
 */
int
suzerain_blasext_cgbdddddmv_s(
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbdddddmv_s */
int
suzerain_blasext_zgbdddddmv_d(
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
        const int incy);

/*!
 * \brief Compute either \f$ B \leftarrow{} \alpha{} D A + \beta{}B \f$ or \f$
 * B \leftarrow{} \alpha{} A D + \beta{}B\f$ for banded \f$A\f$, diagonal
 * \f$D\f$, and banded \f$B\f$.
 *
 * Matrices \f$ A \f$ and \f$ B \f$ both have band storage and must have the
 * same shape and same number of super- and subdiagonals.  All three matrices
 * may be generally strided.  The operation and interface differs from the
 * BLAS' gb_diag_scale_acc.
 *
 * \param side One of 'L' or 'R' indicating whether \f$D\f$ should be applied
 *        to the left or right side of \f$A\f$, respectively.
 * \param m Number of rows in matrices \f$ A \f$ and \f$ B \f$
 *          and size of matrix \f$ D\f$ when <tt>side == 'L'</tt>.
 * \param n Number of columns in matrices \f$ A \f$ and \f$ B \f$
 *          and size of matrix \f$ D\f$ when <tt>side == 'R'</tt>.
 * \param kl Number of subdiagonals in band storage of \c a and \c b.
 * \param ku Number of superdiagonals in band storage of \c a and \c b.
 * \param alpha Multiplicative scalar \f$ \alpha \f$
 * \param d Diagonal storage of matrix \f$ D \f$.
 * \param ldd Nonnegative stride between diagonal entries in \c d.
 * \param a General band storage of the matrix \f$ A \f$.
 * \param inca Strictly positive stride between values in \c a.
 * \param lda Leading dimension of \c a.
 * \param beta Multiplicative scalar \f$ \beta \f$
 * \param b General band storage of the matrix \f$ B \f$.
 * \param incb Strictly positive stride between values in \c b.
 * \param ldb Leading dimension of \c b.
 *
 * \see A BLAS reference for general band matrix storage requirements.
 */
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
        int ldb);

/*! \copydoc suzerain_blasext_sgb_diag_scale_acc */
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
        int ldb);

/*! \copydoc suzerain_blasext_sgb_diag_scale_acc */
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
        int ldb);

/*! \copydoc suzerain_blasext_sgb_diag_scale_acc */
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
        int ldb);

/*!
 * \brief Compute either \f$ B \leftarrow{} \alpha{} D A + \beta{}B \f$ or \f$
 * B \leftarrow{} \alpha{} A D + \beta{}B \f$ for real banded \f$A\f$, real
 * diagonal \f$D\f$, and complex banded \f$B\f$.
 *
 * \copydetails suzerain_blasext_zgb_diag_scale_acc
 */
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
        int ldb);

/*!
 * \brief Compute either \f$ B \leftarrow{} (\alpha_0 D_0 + \alpha_1 D_1) A +
 * \beta{}B \f$ or \f$ B \leftarrow{} A (\alpha_0 D_0 + \alpha_1 D_1) +
 * \beta{}B \f$ for real banded \f$A\f$, real diagonal \f$D_0\f$, real diagonal
 * \f$D_1\f$, and complex banded \f$B\f$.
 *
 * Matrices \f$ A \f$ and \f$ B \f$ both have banded storage.  Any matrix may
 * be generally strided.  Matrices \f$ D_0 \f$ and \f$ D_1 \f$ may be aliased.
 *
 * \param side One of 'L' or 'R' indicating whether the diagonal matrices
 *        should be applied to the left or right side of \f$A\f$, respectively.
 * \param m Number of rows in matrices \f$ A \f$ and \f$ B \f$
 *          and size of diagonal matrices when <tt>side == 'L'</tt>.
 * \param n Number of columns in matrices \f$ A \f$ and \f$ B \f$
 *          and size of diagonal matrices when <tt>side == 'R'</tt>.
 * \param kl Number of subdiagonals in band storage of \c a and \c b.
 * \param ku Number of superdiagonals in band storage of \c a and \c b.
 * \param alpha0 Multiplicative scalar \f$ \alpha_0 \f$
 * \param d0 Diagonal storage of matrix \f$ D_0 \f$.
 * \param ldd0 Nonnegative stride between diagonal entries in \c d0.
 * \param alpha1 Multiplicative scalar \f$ \alpha_1 \f$
 * \param d1 Diagonal storage of matrix \f$ D_1 \f$.
 * \param ldd1 Nonnegative stride between diagonal entries in \c d1.
 * \param a General band storage of the matrix \f$ A \f$.
 * \param inca Strictly positive stride between values in \c a.
 * \param lda Leading dimension of \c a.
 * \param beta Multiplicative scalar \f$ \beta \f$
 * \param b General band storage of the matrix \f$ B \f$.
 * \param incb Strictly positive stride between values in \c b.
 * \param ldb Leading dimension of \c b.
 *
 * \see A BLAS reference for banded matrix storage requirements.
 */
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
        int ldb);

/*!
 * \brief Compute either \f$ B \leftarrow{} (\alpha_0 D_0 + \alpha_1 D_1 +
 * \alpha_2 D_2) A + \beta{}B \f$ or \f$ B \leftarrow{} A (\alpha_0 D_0 +
 * \alpha_1 D_1 + \alpha_2 D_2) + \beta{}B \f$ for real banded \f$A\f$, real
 * diagonal \f$D_0\f$, real diagonal \f$D_1\f$, real diagonal \f$D_2\f$,  and
 * complex banded \f$B\f$.
 *
 * Matrices \f$ A \f$ and \f$ B \f$ both have banded storage.  Any matrix may
 * be generally strided.  Matrices \f$ D_0 \f$, \f$ D_1 \f$, and \f$ D_2 \f$
 * may be aliased.
 *
 * \param side One of 'L' or 'R' indicating whether the diagonal matrices
 *        should be applied to the left or right side of \f$A\f$, respectively.
 * \param m Number of rows in matrices \f$ A \f$ and \f$ B \f$
 *          and size of diagonal matrices when <tt>side == 'L'</tt>.
 * \param n Number of columns in matrices \f$ A \f$ and \f$ B \f$
 *          and size of diagonal matrices when <tt>side == 'R'</tt>.
 * \param kl Number of subdiagonals in band storage of \c a and \c b.
 * \param ku Number of superdiagonals in band storage of \c a and \c b.
 * \param alpha0 Multiplicative scalar \f$ \alpha_0 \f$
 * \param d0 Diagonal storage of matrix \f$ D_0 \f$.
 * \param ldd0 Nonnegative stride between diagonal entries in \c d0.
 * \param alpha1 Multiplicative scalar \f$ \alpha_1 \f$
 * \param d1 Diagonal storage of matrix \f$ D_1 \f$.
 * \param ldd1 Nonnegative stride between diagonal entries in \c d1.
 * \param alpha2 Multiplicative scalar \f$ \alpha_2 \f$
 * \param d2 Diagonal storage of matrix \f$ D_2 \f$.
 * \param ldd2 Nonnegative stride between diagonal entries in \c d2.
 * \param a General band storage of the matrix \f$ A \f$.
 * \param inca Strictly positive stride between values in \c a.
 * \param lda Leading dimension of \c a.
 * \param beta Multiplicative scalar \f$ \beta \f$
 * \param b General band storage of the matrix \f$ B \f$.
 * \param incb Strictly positive stride between values in \c b.
 * \param ldb Leading dimension of \c b.
 *
 * \see A BLAS reference for banded matrix storage requirements.
 */
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
        int ldb);

/*!
 * \brief Compute either \f$ B \leftarrow{} (\alpha_0 D_0 + \alpha_1 D_1 +
 * \alpha_2 D_2 + \alpha_3 D_3) A + \beta{}B \f$ or \f$ B \leftarrow{} A
 * (\alpha_0 D_0 + \alpha_1 D_1 + \alpha_2 D_2 + \alpha_3 D_3) + \beta{}B \f$
 * for real banded \f$A\f$, real diagonal \f$D_0\f$, real diagonal \f$D_1\f$,
 * real diagonal \f$D_2\f$, real diagonal \f$D_3\f$ and complex banded \f$B\f$.
 *
 * Matrices \f$ A \f$ and \f$ B \f$ both have banded storage.  Any matrix may
 * be generally strided.  Matrices \f$ D_0 \f$, \f$ D_1 \f$, \f$ D_2 \f$, and
 * \f$ D_3 \f$ may be aliased.
 *
 * \param side One of 'L' or 'R' indicating whether the diagonal matrices
 *        should be applied to the left or right side of \f$A\f$, respectively.
 * \param m Number of rows in matrices \f$ A \f$ and \f$ B \f$
 *          and size of diagonal matrices when <tt>side == 'L'</tt>.
 * \param n Number of columns in matrices \f$ A \f$ and \f$ B \f$
 *          and size of diagonal matrices when <tt>side == 'R'</tt>.
 * \param kl Number of subdiagonals in band storage of \c a and \c b.
 * \param ku Number of superdiagonals in band storage of \c a and \c b.
 * \param alpha0 Multiplicative scalar \f$ \alpha_0 \f$
 * \param d0 Diagonal storage of matrix \f$ D_0 \f$.
 * \param ldd0 Nonnegative stride between diagonal entries in \c d0.
 * \param alpha1 Multiplicative scalar \f$ \alpha_1 \f$
 * \param d1 Diagonal storage of matrix \f$ D_1 \f$.
 * \param ldd1 Nonnegative stride between diagonal entries in \c d1.
 * \param alpha2 Multiplicative scalar \f$ \alpha_2 \f$
 * \param d2 Diagonal storage of matrix \f$ D_2 \f$.
 * \param ldd2 Nonnegative stride between diagonal entries in \c d2.
 * \param alpha3 Multiplicative scalar \f$ \alpha_3 \f$
 * \param d3 Diagonal storage of matrix \f$ D_3 \f$.
 * \param ldd3 Nonnegative stride between diagonal entries in \c d3.
 * \param a General band storage of the matrix \f$ A \f$.
 * \param inca Strictly positive stride between values in \c a.
 * \param lda Leading dimension of \c a.
 * \param beta Multiplicative scalar \f$ \beta \f$
 * \param b General band storage of the matrix \f$ B \f$.
 * \param incb Strictly positive stride between values in \c b.
 * \param ldb Leading dimension of \c b.
 *
 * \see A BLAS reference for banded matrix storage requirements.
 */
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
        int ldb);

/*!
 * \brief Compute two-strided \f$ Y \leftarrow{} \alpha{}X+\beta{}Y \f$
 *        where \f$ X \f$ and \f$ Y \f$ have split and interleaved
 *        complex storage, respectively.
 *
 * Accumulated interleaved complex column-major matrices into
 * two split complex column-major matrices, one each for the real and
 * imaginary parts.  The input and output arguments must not coincide.
 * Note that \c incx and \c ldx are specified in terms of complex strides,
 * and not real-valued strides.
 *
 * \param m Number of rows in column-major matrices \f$ X \f$ and \f$ Y \f$.
 * \param n Number of columns in column-major matrices \f$ X \f$ and \f$ Y \f$.
 * \param alpha Scale factor \f$ \alpha \f$ to apply to matrix \f$ X \f$.
 * \param x Beginning of interleaved complex storage for matrix \f$ X \f$.
 *          Real and imaginary components must be in adjacent locations.
 * \param incx Increment between values in a column of \f$ X \f$,
 *             specified in terms of complex elements like <tt>double[2]</tt>.
 * \param ldx  Leading distance between values in a row of \f$ X \f$,
 *             specified in terms of complex elements like <tt>double[2]</tt>.
 * \param beta Scale factor \f$ \beta \f$ to apply to matrix \f$ Y \f$
 * \param y_re Beginning of the storage for \f$\operatorname{real} Y \f$.
 * \param incy_re Increment between values in a column of
 *                \f$ \operatorname{real} Y \f$, specified in terms of
 *                <tt>sizeof(double)</tt>.
 * \param ldy_re Leading distance between values in a row of
 *                \f$ \operatorname{real} Y \f$, specified in terms of
 *                <tt>sizeof(double)</tt>.
 * \param y_im Beginning of the storage for \f$\operatorname{imag} Y \f$.
 * \param incy_im Increment between values in a column of
 *                \f$ \operatorname{imag} Y \f$, specified in terms of
 *                <tt>sizeof(double)</tt>.
 * \param ldy_im Leading distance between values in a row of
 *                \f$ \operatorname{imag} Y \f$, specified in terms of
 *                <tt>sizeof(double)</tt>.
 */
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
        const int ldy_im);

/*!
 * \brief Compute \f$ \left|\left|A\right|\right| \f$ for a general
 *        banded matrix.
 *
 * \param m Number of rows in matrix \c a.
 * \param n Number of columns in matrix \c a.
 * \param kl Number of subdiagonals in band storage of \c a.
 * \param ku Number of superdiagonals in band storage of \c a.
 * \param a General band storage for matrix \f$ A \f$.
 * \param lda Leading dimension of \c a.
 * \param norm1 The one norm of the matrix.
 *
 * \see A BLAS reference for more details, especially for general
 *      band matrix storage requirements.
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 */
int
suzerain_blasext_sgbnorm1(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const float *a,
        const int lda,
        float *norm1);

/*! \copydoc suzerain_blasext_sgbnorm1 */
int
suzerain_blasext_dgbnorm1(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const double *a,
        const int lda,
        double *norm1);

/*! \copydoc suzerain_blasext_sgbnorm1 */
int
suzerain_blasext_cgbnorm1(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const complex_float *a,
        const int lda,
        float *norm1);

/*! \copydoc suzerain_blasext_sgbnorm1 */
int
suzerain_blasext_zgbnorm1(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const complex_double *a,
        const int lda,
        double *norm1);

/*! @} */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __SUZERAIN_BLAS_ET_AL_H */
