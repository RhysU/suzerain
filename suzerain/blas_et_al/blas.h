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

#ifndef SUZERAIN_BLAS_ET_AL_BLAS_H
#define SUZERAIN_BLAS_ET_AL_BLAS_H

/*!\file
 * Wraps external BLAS routines necessary for Suzerain.
 * Provided to insulate Suzerain from potential variations in
 * type signatures as well as to consolidate all Fortran-from-C
 * parameter differences.
 */

#include <suzerain/common.h>
#include <suzerain/complex.h>

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

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_BLAS_ET_AL_BLAS_H */
