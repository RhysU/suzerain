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

#ifndef SUZERAIN_BLAS_ET_AL_H
#define SUZERAIN_BLAS_ET_AL_H

/** @file
 * Wraps and extends external implementations of BLAS, LAPACK, et al.
 */

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
#define SUZERAIN_BLAS_ALIGNMENT (16)

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
 *
 * \see A BLAS reference for more details.
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
 *
 * \see A BLAS reference for more details.
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
 * \see A BLAS reference for general band matrix storage requirements.
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
        const int incy);

/*! \copydoc suzerain_blas_sgbmv_external */
void
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
        const int incy);

/*! \copydoc suzerain_blas_sgbmv_external */
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
        const int incy);

/*! \copydoc suzerain_blas_sgbmv_external */
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
        const int incy);

/*! \copydoc suzerain_blas_sgbmv */
int
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
        const int incy);

/*! \copydoc suzerain_blas_sgbmv */
int
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
        const int incy);

/*! \copydoc suzerain_blas_sgbmv */
int
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
        const int incy);

/*! \copydoc suzerain_blas_sgbmv */
int
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
        const int incy);

/*! \copydoc suzerain_blas_sgbmv */
int
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
 * \see A BLAS reference for general band matrix storage requirements.
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
        const int incy);

/*! \copydoc suzerain_blas_ssbmv_external */
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
        const int incy);

/*! \copydoc suzerain_blas_csbmv_s_c */
int
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
        const int incy);

/*! \copydoc suzerain_blas_ssbmv */
int
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
        const int incy);

/*! \copydoc suzerain_blas_csbmv_s_s */
int
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
 * \see A BLAS reference for general band matrix storage requirements.
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
 *
 * \see A BLAS reference for more details.
 * @{
 */

/*! @} */

/*! \name LAPACK operations
 *
 * \see A LAPACK reference for more details.
 * @{
 */

/**
 * Copies all or part of a two-dimensional matrix A to another matrix B using
 * LAPACK's LACPY.
 *
 * @param uplo Specifies the part of the matrix A to be copied to matrix B.
 *             If 'U' or 'u' the upper triangular part is copied.
 *             If 'L' or 'l' the lower triangular part is copied.
 *             Otherwise, the entire matrix is copied.
 * @param m    Number of rows of the matrix A.
 * @param n    Number of columns of the matrix A.
 * @param a    The \c m by \c n matrix A.
 * @param lda  The leading dimension of \c a.
 * @param b    The \c m by \c n matrix B.
 * @param ldb  The leading dimension of \c b.
 */
void suzerain_lapack_slacpy(
                char uplo,
                const int m,
                const int n,
                const float * a,
                const int lda,
                float *b,
                const int ldb);

/*! \copydoc suzerain_lapack_slacpy */
void suzerain_lapack_dlacpy(
                char uplo,
                const int m,
                const int n,
                const double * a,
                const int lda,
                double *b,
                const int ldb);

/*! \copydoc suzerain_lapack_slacpy */
void suzerain_lapack_clacpy(
                char uplo,
                const int m,
                const int n,
                const complex_float * a,
                const int lda,
                complex_float *b,
                const int ldb);

/*! \copydoc suzerain_lapack_slacpy */
void suzerain_lapack_zlacpy(
                char uplo,
                const int m,
                const int n,
                const complex_double * a,
                const int lda,
                complex_double *b,
                const int ldb);

/**
 * Determines machine parameters for floating-point arithmetic.
 *
 * The input parameter \c cmach should specify which of the
 * following values is desired:
 * <dl>
 * <dt><tt>E</tt>, <tt>e</tt></dt>
 * <dd><tt>eps</tt>, relative machine precision</dd>
 * <dt><tt>S</tt>, <tt>s</tt></dt>
 * <dd><tt>sfmin</tt>, safe minimum,
 *     such that <tt>1/sfmin</tt> does not overflow</dd>
 * <dt><tt>B</tt>, <tt>b</tt></dt>
 * <dd><tt>base</tt>, base of the machine</dd>
 * <dt><tt>P</tt>, <tt>p</tt></dt>
 * <dd><tt>eps*base</tt></dd>
 * <dt><tt>N</tt>, <tt>n</tt></dt>
 * <dd><tt>t</tt>, number of (base) digits in the mantissa</dd>
 * <dt><tt>R</tt>, <tt>r</tt></dt>
 * <dd><tt>rnd</tt>, 1.0 when rounding occurs in addition, 0.0 otherwise</dd>
 * <dt><tt>M</tt>, <tt>m</tt></dt>
 * <dd><tt>emin</tt>, minimum exponent before (gradual) underflow</dd>
 * <dt><tt>U</tt>, <tt>u</tt></dt>
 * <dd><tt>rmin</tt>, underflow threshold - <tt>base**(emin-1)</tt></dd>
 * <dt><tt>L</tt>, <tt>l</tt></dt>
 * <dd><tt>emax</tt>, largest exponent before overflow</dd>
 * <dt><tt>O</tt>, <tt>o</tt></dt>
 * <dd><tt>rmax</tt>, overflow threshold  - <tt>(base**emax)*(1-eps)</tt></dd>
 * </dl>
 *
 * @param cmach The desired machine-specific quantity.
 *
 * @return the value requested per parameter \c cmach.
 */
float suzerain_lapack_slamch(char cmach);

/*! \copydoc suzerain_lapack_slamch */
double suzerain_lapack_dlamch(char cmach);

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
 * \see A LAPACK reference for the <tt>lda >= 2*kl + ku + 1</tt> storage
 *      requirement and the resulting factored storage format.
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
 * \param[in,out] afb   Dimension (<tt>ldafb</tt>,<tt>n</tt>)
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
 * \param[in,out] afb   Dimension (<tt>ldafb</tt>,<tt>n</tt>)
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
 * \brief Iteratively refine the solution to a banded system of linear
 * equations \f$AX=B\f$ and obtain error estimates.
 *
 * \param[in]     trans
 * \param[in]     n
 * \param[in]     kl
 * \param[in]     ku
 * \param[in]     nrhs
 * \param[in,out] ab    Dimension (<tt>ldab</tt>,<tt>n</tt>)
 * \param[in]     ldab  Minimum <tt>kl + ku + 1</tt>
 * \param[in,out] afb   Dimension (<tt>ldafb</tt>,<tt>n</tt>)
 * \param[in]     ldafb Minimum <tt>2*kl + ku + 1</tt>
 * \param[in,out] ipiv  Dimension \c n
 * \param[in,out] b     Dimension (<tt>ldb</tt>,<tt>nrhs</tt>)
 * \param[in]     ldb   Minimum \c n
 * \param[in,out] x     Dimension (<tt>ldx</tt>,<tt>nrhs</tt>)
 * \param[in]     ldx   Minimum \c n
 * \param[out]    ferr  Dimension \c nrhs
 * \param[out]    berr  Dimension \c nrhs
 * \param[out]    work  Dimension <tt>3*n</tt>
 * \param[out]    iwork Dimension \c n
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 */
int
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
        int *iwork);

/*! \copydoc suzerain_lapack_sgbrfs */
int
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
        int *iwork);

/*!
 * \brief Iteratively refine the solution to a banded system of linear
 * equations \f$AX=B\f$ and obtain error estimates.
 *
 * \param[in]     trans
 * \param[in]     n
 * \param[in]     kl
 * \param[in]     ku
 * \param[in]     nrhs
 * \param[in,out] ab    Dimension (<tt>ldab</tt>,<tt>n</tt>)
 * \param[in]     ldab  Minimum <tt>kl + ku + 1</tt>
 * \param[in,out] afb   Dimension (<tt>ldafb</tt>,<tt>n</tt>)
 * \param[in]     ldafb Minimum <tt>2*kl + ku + 1</tt>
 * \param[in,out] ipiv  Dimension \c n
 * \param[in,out] b     Dimension (<tt>ldb</tt>,<tt>nrhs</tt>)
 * \param[in]     ldb   Minimum \c n
 * \param[in,out] x     Dimension (<tt>ldx</tt>,<tt>nrhs</tt>)
 * \param[in]     ldx   Minimum \c n
 * \param[out]    ferr  Dimension \c nrhs
 * \param[out]    berr  Dimension \c nrhs
 * \param[out]    work  Dimension <tt>2*n</tt>
 * \param[out]    rwork Dimension \c n
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 */
int
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
        float *rwork);

/*! \copydoc suzerain_lapack_cgbrfs */
int
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
        double *rwork);

/*!
 * \brief Compute the one norm, or the Frobenius norm, or the infinity
 * norm, or the element of largest absolute value of square band matrix.
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
 * These are extensions built atop the BLAS, on vendor-specific BLAS-like
 * routines, and/or on custom coded loops.  Some of these extensions resemble
 * those found in the <a href="http://www.netlib.org/xblas/">XBLAS</a>.
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbdmv_s_c_external */
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbdmv_s_c */
int
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} D A x + \beta{} y \f$ for complex
 * \f$\alpha{}\f$, \f$\beta\f$, and \f$y\f$ but real-valued \f$D\f$, \f$A\f$,
 * and \f$x\f$.
 *
 * \copydetails suzerain_blasext_sgbdmv
 */
int
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbdmv_s_s */
int
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbddmv_s_c_external */
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbddmv_s_c */
int
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \left(\alpha_{0} D_0 + \alpha_{1}
 * D_1\right) A x + \beta{} y \f$ for complex \f$\alpha{}\f$, \f$\beta\f$, and
 * \f$y\f$ but real-valued \f$D_0\f$, \f$D_1\f$, \f$A\f$, and \f$x\f$.
 *
 * \copydetails suzerain_blasext_sgbddmv
 */
int
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbddmv_s_s */
int
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbdddmv_s_c_external */
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbdddmv_s_c */
int
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \left(\alpha_{0} D_0 + \alpha_{1} D_1 +
 * \alpha_{2} D_2\right) A x + \beta{} y \f$ for complex \f$\alpha{}\f$,
 * \f$\beta\f$, and \f$y\f$ but real-valued \f$D_0\f$, \f$D_1\f$, \f$D_2\f$,
 * \f$A\f$, and \f$x\f$.
 *
 * \copydetails suzerain_blasext_sgbdddmv
 */
int
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbdddmv_s_s */
int
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbddddmv_s_c_external */
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbddddmv_s_c */
int
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \left(\alpha_{0} D_0 + \alpha_{1} D_1 +
 * \alpha_{2} D_2 + \alpha_{3} D_3 \right) A x + \beta{} y \f$ for complex
 * \f$\alpha{}\f$, \f$\beta\f$, and \f$y\f$ but real-valued \f$D_0\f$,
 * \f$D_1\f$, \f$D_2\f$, \f$D_3\f$, \f$A\f$, and \f$x\f$.
 *
 * \copydetails suzerain_blasext_sgbddddmv
 */
int
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbddddmv_s_s */
int
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbdddddmv_s_c_external */
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbdddddmv_s_c */
int
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \left(\alpha_{0} D_0 + \alpha_{1} D_1 +
 * \alpha_{2} D_2 + \alpha_{3} D_3 + \alpha_{4} D_4 \right) A x + \beta{} y \f$
 * for complex \f$\alpha{}\f$, \f$\beta\f$, and \f$y\f$ but real-valued
 * \f$D_0\f$, \f$D_1\f$, \f$D_2\f$, \f$D_3\f$, \f$D_4\f$, \f$A\f$, and \f$x\f$.
 *
 * \copydetails suzerain_blasext_sgbdddddmv
 */
int
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
        const int incy);

/*! \copydoc suzerain_blasext_cgbdddddmv_s_s */
int
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
 * \see A BLAS reference for general banded matrix storage requirements.
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
 * \see A BLAS reference for general banded matrix storage requirements.
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
 * real diagonal \f$D_2\f$, real diagonal \f$D_3\f$,and complex banded \f$B\f$.
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
 * \see A BLAS reference for general banded matrix storage requirements.
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
 * \brief Compute either \f$ B \leftarrow{} (\alpha_0 D_0 + \alpha_1 D_1 +
 * \alpha_2 D_2 + \alpha_3 D_3 + \alpha_4 D_4) A + \beta{}B \f$ or \f$ B
 * \leftarrow{} A (\alpha_0 D_0 + \alpha_1 D_1 + \alpha_2 D_2 + \alpha_3 D_3 +
 * \alpha_4 D_4) + \beta{}B \f$ for real banded \f$A\f$, real diagonal
 * \f$D_0\f$, real diagonal \f$D_1\f$, real diagonal \f$D_2\f$, real diagonal
 * \f$D_3\f$, real diagonal \f$D_4\f$, and complex banded \f$B\f$.
 *
 * Matrices \f$ A \f$ and \f$ B \f$ both have banded storage.  Any matrix may
 * be generally strided.  Matrices \f$ D_0 \f$, \f$ D_1 \f$, \f$ D_2 \f$, \f$
 * D_3 \f$, and \f$ D_4 \f$ may be aliased.
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
 * \param alpha4 Multiplicative scalar \f$ \alpha_4 \f$
 * \param d4 Diagonal storage of matrix \f$ D_4 \f$.
 * \param ldd4 Nonnegative stride between diagonal entries in \c d4.
 * \param a General band storage of the matrix \f$ A \f$.
 * \param inca Strictly positive stride between values in \c a.
 * \param lda Leading dimension of \c a.
 * \param beta Multiplicative scalar \f$ \beta \f$
 * \param b General band storage of the matrix \f$ B \f$.
 * \param incb Strictly positive stride between values in \c b.
 * \param ldb Leading dimension of \c b.
 *
 * \see A BLAS reference for general banded matrix storage requirements.
 */
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
 * \see A BLAS reference for general band matrix storage requirements.
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

/**
 * Demote a contiguous double precision vector to single
 * precision as an in-place operation.
 *
 * \param[in]     n Number of elements in \c x.
 * \param[in,out] x Source and target vector.
 *                On input, \c x contains double precision data.
 *                On output, \c x contains single precision data.
 *
 * @return On success, return zero.  A negative return
 *         value indicates an error in the corresponding
 *         argument.
 */
int
suzerain_blasext_ddemote(
        const int n,
        void *x);

/**
 * Promote a contiguous single precision vector to
 * double precision as an in-place operation.
 *
 * \param[in]     n Number of elements in \c x.
 * \param[in,out] x Source and target vector.
 *                On input, \c x contains single precision data.
 *                On output, \c x contains double precision data.
 *
 * @return On success, return zero.  A negative return
 *         value indicates an error in the corresponding
 *         argument.
 */
int
suzerain_blasext_dpromote(
        const int n,
        void *x);

/**
 * Demote a contiguous complex double precision vector to
 * complex single precision as an in-place operation.
 *
 * \param[in]     n Number of elements in \c x.
 * \param[in,out] x Source and target vector.
 *                On input, \c x contains double precision data.
 *                On output, \c x contains single precision data.
 *
 * @return On success, return zero.  A negative return
 *         value indicates an error in the corresponding
 *         argument.
 */
int
suzerain_blasext_zdemote(
        const int n,
        void *x);

/**
 * Promote a contiguous complex single precision vector to
 * complex double precision as an in-place operation.
 *
 * \param[in]     n Number of elements in \c x.
 * \param[in,out] x Source and target vector.
 *                On input, \c x contains single precision data.
 *                On output, \c x contains double precision data.
 *
 * @return On success, return zero.  A negative return
 *         value indicates an error in the corresponding
 *         argument.
 */
int
suzerain_blasext_zpromote(
        const int n,
        void *x);

/**
 * The goofiest kernel you'll ever come across:
 * \verbatim
 *
 *     /w0\            /a0 a5 a10 a15 a20\   /x0\
 *     |w1|            |a1 a6 a11 a16 a21|   |x1|
 *     |w2| += alpha * |a2 a7 a12 a17 a22| * |x2|
 *     |w3|            |a3 a8 a13 a18 a23|   |x3|
 *     \w4/            \a4 a9 a14 a19 a24/   \x4/
 *
 *                     /b0 b5 b10 b15 b20\   /x0\
 *                     |b1 b6 b11 b16 b21|   |x1|
 *           + beta  * |b2 b7 b12 b17 b22| * |x2|
 *                     |b3 b8 b13 b18 b23|   |x3|
 *                     \b4 b9 b14 b19 b24/   \x4/
 *
 *                     /c0 c5 c10 c15 c20\   /y0\
 *                     |c1 c6 c11 c16 c21|   |y1|
 *           + gamma * |c2 c7 c12 c17 c22| * |y2|
 *                     |c3 c8 c13 c18 c23|   |y3|
 *                     \c4 c9 c14 c19 c24/   \y4/
 *
 * \endverbatim
 *
 * \param [in ] alpha  Scaling factor on first matrix
 * \param [in ] a      Five-by-five matrix, column-major contiguous
 * \param [in ] beta   Scaling factor on second matrix
 * \param [in ] b      Five-by-five matrix, column-major contiguous
 * \param [in ] x      First five vector to be hit by \e two matrices
 * \param [in ] gamma  Scaling factor on third matrix
 * \param [in ] c      Five-by-five matrix, column-major contiguous
 * \param [in ] y      Second five vector to be hit by \e one matrix
 * \param [out] w0     First component of w output vector
 * \param [out] w1     Second component of w output vector
 * \param [out] w2     Third component of w output vector
 * \param [out] w3     Fourth component of w output vector
 * \param [out] w4     Fifth component of w output vector
 *
 * @return Zero, always.
 */
int
suzerain_blasext_zgedsummv55(
        const complex_double                     alpha,
        const         double * SUZERAIN_RESTRICT a,
        const complex_double                     beta,
        const         double * SUZERAIN_RESTRICT b,
        const complex_double * SUZERAIN_RESTRICT x,
        const complex_double                     gamma,
        const         double * SUZERAIN_RESTRICT c,
        const complex_double * SUZERAIN_RESTRICT y,
              complex_double * SUZERAIN_RESTRICT w0,
              complex_double * SUZERAIN_RESTRICT w1,
              complex_double * SUZERAIN_RESTRICT w2,
              complex_double * SUZERAIN_RESTRICT w3,
              complex_double * SUZERAIN_RESTRICT w4);

/*! @} */

/*! \name LAPACK-like extensions
 * These are extensions built atop LAPACK, the BLAS, and/or on custom coded
 * loops.  Some of these extensions resemble newer LAPACK features which
 * have not trickled down to banded matrices from general matrices.
 * @{
 */

/*!
 * \brief Compute the solution to a banded system of linear equations \f$ AX=B
 * \f$ using an in-place LU factorization.  This is implemented as a
 * possible call to GBTRF followed by a call to GBTRS.
 *
 * Transposes of \f$ A \f$ can be taken using the \c trans parameter.
 *
 * \param[in,out] fact  One of 'N' or 'F'.  When 'N' on entry, \f$ A \f$ is
 *                      factorized and \c fact is set to 'F' on exit.
 *                      Otherwise, \f$ A \f$ must have been factorized in-place
 *                      by a previous call to this method employing 'N'.
 * \param[in]     trans One of 'N', 'T', or 'C' for no transpose,
 *                      a transpose, or a conjugate transpose, respectively.
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
 */
int
suzerain_lapackext_sgbsv(
        char * const fact,
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        float *ab,
const int ldab,
        int *ipiv,
        float *b,
        const int ldb);

/*! \copydoc suzerain_lapackext_sgbsv */
int
suzerain_lapackext_dgbsv(
        char * const fact,
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        double *ab,
        const int ldab,
        int *ipiv,
        double *b,
        const int ldb);

/*! \copydoc suzerain_lapackext_sgbsv */
int
suzerain_lapackext_cgbsv(
        char * const fact,
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        complex_float *ab,
        const int ldab,
        int *ipiv,
        complex_float *b,
        const int ldb);

/*! \copydoc suzerain_lapackext_sgbsv */
int
suzerain_lapackext_zgbsv(
        char * const fact,
        const char trans,
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
 * \brief Compute the solution to a banded system of linear equations \f$A x =
 * b \f$ using single precision LU factorization followed by double precision
 * iterative refinement.
 *
 * The overall functionality is an extension of the capability of LAPACK's <a
 * href="http://www.netlib.org/lapack/double/dsgesv.f">DSGESV</a> but written
 * for general banded matrices.  The approach follows DSGESV which is based on
 * the DGESIRSV algorithm presented in June 2006 as <a
 * href="http://www.netlib.org/lapack/lawnspdf/lawn175.pdf">LAWN 175:
 * Exploiting the Performance of 32 bit Floating Point Arithmetic in Obtaining
 * 64 bit Accuracy (Revisiting Iterative Refinement for Linear Systems)</a> by
 * Langou et al.  Unlike their algorithm, this routine
 * <ol>
 *   <li>permits the user to provide a preexisting single
 *       or double precision factorization,</li>
 *   <li>returns the single or double precision factorization
 *       for subsequent reuse,</li>
 *   <li>permits the caller to provide the Frobenius norm of \f$A\f$,</li>
 *   <li>assumes only one right hand side is of interest
 *       (above changes, however, permit efficient repeated invocation),</li>
 *   <li>requires that matrix storage is contiguous (i.e.
 *       <tt>ldab == kl + ku + 1</tt>, <tt>ldafb == 2*kl + ku + 1</tt>),</li>
 *       for reasons of cache-friendliness,</li>
 *   <li>performs precision demotion and promotion in-place
 *       for reasons of cache-friendliness,</li>
 *   <li>performs iterative refinement after the double-precision
 *       fallback LU factorization if the direct solution fails
 *       to meet the solution tolerance,</li>
 *   <li>permits the user to scale the required tolerance by some
 *       amount to accommodate situations where lower precision
 *       results are acceptable,</li>
 *   <li>requires the user to specify the maximum number of refinements
 *       to attempt (Langou et al. suggested thirty, <tt>?gbrfsx</tt> suggests
 *       up to one hundred in aggressive circumstances like providing
 *       approximate factorizations),</li>
 *   <li>makes accessible the number iterative refinements attempted,</li>
 *   <li>and makes accessible the final residual vector and residual.</li>
 * </ol>
 * All in all, these changes are designed to permit combining single precision
 * factorization along with possibly using approximate factorizations and
 * refinement as a means to Newton iteration.  The consistent error bound
 * testing in either single precision factorization or double precision
 * factorization is intended to permit the results to be \e indistinguishable
 * regardless of solution method.
 *
 * \param[in,out] fact  On input, 'N' if factorization must be performed;
 *                      'S' if a single precision factorization is
 *                      provided in \c afb and \c ipiv; 'D' if a double
 *                      precision factorization is provided in \c afb
 *                      and \c ipiv.  On output, either 'S' or 'D'
 *                      indicating which factorization precision has
 *                      been returned in \c afb and \c ipiv.
 * \param[in,out] apprx On input, if \c fact is 'S' or 'D' and \c apprx
 *                      is nonzero, assume the factorization provided in
 *                      \c afb and \c ipiv is \f$ LUP = A + \Delta{}A\f$
 *                      for some perturbation \f$ \Delta{}A \f$.
 *                      If the provided factorization fails to deliver
 *                      adequate convergence as described under \c aiter,
 *                      \c ab will be factorized and \c apprx set to
 *                      zero on return.  If \c fact is 'N' on input,
 *                      \c apprx is ignored and set to zero on return.
 * \param[in]     aiter At most \c aiter refinement steps will be attempted
 *                      before it can be confirmed that the residual
 *                      is decaying by a factor of two at each step.
 *                      Parameter \c aiter is used during mixed-precision
 *                      refinements and during double-precision refinements
 *                      when \c apprx is nonzero.
 * \param[in]     trans If 'N' solve \f$A      x = b\f$.
 *                      If 'T' solve \f$A^\top x = b\f$.
 *                      If 'C' solve \f$A^H    x = b\f$.
 * \param[in]     n     Number of rows and columns in \f$ A \f$.
 * \param[in]     kl    Number of subdiagonals in \f$ A \f$.
 * \param[in]     ku    Number of superdiagonals in \f$ A \f$.
 * \param[in]     ab    Double precision matrix in banded storage of
 *                      dimension (<tt>ldab == kl + ku + 1</tt>, <tt>n</tt>).
 * \param[in,out] afrob The Frobenius norm of \f$ A \f$.
 *                      If negative on entry and \c tolsc is strictly
 *                      positive, the norm is computed and returned to
 *                      permit caching the result.  Otherwise, \c afrob
 *                      is not modified.
 * \param[in,out] afb   Single or double precision LU factorization of \f$A\f$
 *                      of dimension (<tt>ldafb == 2*kl + ku + 1</tt>,
 *                      <tt>n</tt>).  Whether a single or double
 *                      precision result is returned is communicated by
 *                      \c fact on return.
 * \param[in,out] ipiv  LU pivoting information for the factorization
 *                      in \c afb of dimension \c n.
 * \param[in]     b     Right hand side \f$ b \f$ of length \c n.
 * \param[out]    x     Computed solution \f$ x \f$ of length \c n.
 * \param[in,out] siter On input, the maximum number of single precision
 *                      iterative refinements to attempt.  Zero specifies
 *                      no single precision refinements are to occur.
 *                      If negative, no single precision factorization
 *                      or solve will be formed and double precision
 *                      processing occurs.  On output, the number of
 *                      such refinements performed.
 * \param[in,out] diter On input, the maximum number of double precision
 *                      iterative refinements to attempt.  Zero specifies
 *                      no double precision refinements are to occur.
 *                      If negative, no double precision factorization or
 *                      solve will be performed.  On output, the number
 *                      of such refinements performed.
 * \param[in,out] tolsc On input, a nonnegative multiplicative factor
 *                      used to scale the Langou et al. tolerance \f$
 *                      \sqrt{n} \text{eps} \|A\|_\text{Fro} \|x\|_2
 *                      \f$ against which the residual is compared as a
 *                      stopping criterion.  The recommended value from
 *                      Langou et al is \c 1.0 to regain full accuracy
 *                      as measured per backward stability.  The special
 *                      value \c 0.0 can be provided to specify that
 *                      double precision machine epsilon should be used
 *                      which, in conjunction with <tt>siter < 0</tt>,
 *                      makes the refinement process act like LAPACK's
 *                      <tt>?gbrfs</tt>.  On output, the fraction of the
 *                      tolerance represented by the returned solution.
 *                      Values greater than one indicate the desired
 *                      tolerance could not be met.
 * \param[out]    r     Solution residual \f$ b - A x \f$.
 * \param[out]    res   2-norm of the residual.
 *                      That is, \f$\|b - A x\|_2\f$.
 *
 * One or both of \c siter or \c diter must be nonnegative on entry.
 *
 * If \c tolsc is identically zero, \c siter should be strictly negative.
 * This requirement is not enforced, but it is madness to think that
 * mixed-precision iterative refinement will provide residuals better
 * than those returned by a full-precision procedure.
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 *         Errors related to the <tt>i</tt>-th argument are indicated
 *         by a return value of <tt>-i</tt>.
 */
int
suzerain_lapackext_dsgbsvx(
        char * const fact,
        int * const apprx,
        const int aiter,
        char trans,
        const int n,
        const int kl,
        const int ku,
        double * const ab,
        double * const afrob,
        double * const afb,
        int * const ipiv,
        double * const b,
        double * const x,
        int * const siter,
        int * const diter,
        double * const tolsc,
        double * const r,
        double * const res);

/*! \copydoc suzerain_lapackext_dsgbsvx */
int
suzerain_lapackext_zcgbsvx(
        char * const fact,
        int * const apprx,
        const int aiter,
        char trans,
        const int n,
        const int kl,
        const int ku,
        complex_double * const ab,
        double * const afrob,
        complex_double * const afb,
        int * const ipiv,
        complex_double * const b,
        complex_double * const x,
        int * const siter,
        int * const diter,
        double * const tolsc,
        complex_double * const r,
        double * const res);

/*! @} */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_BLAS_ET_AL_H */
