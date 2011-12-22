/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2011 The PECOS Development Team
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
 * gbmv.h: provides custom, BLAS-like general band matrix-vector operations
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_GBMV_H__
#define __SUZERAIN_GBMV_H__

#ifdef __cplusplus
extern "C" {
#endif

/** \file
 * Provides custom, BLAS-like general band matrix-vector operations.
 * Includes mixed real/complex operations and fixed bandwidth kernels.
 */

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$ using an
 * internal, general bandwidth gbmv implementation.  Using an implementation
 * from BLAS is should be preferred to this implementation directly.
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
 * \return Zero on success and a BLAS-like error code otherwise.
 * \see A BLAS reference for more details, especially for general
 *      band storage matrix requirements.
 */
int
suzerain_gbmv_s(
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

/*! \copydoc suzerain_gbmv_s */
int
suzerain_gbmv_d(
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

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$ for complex
 * \f$\alpha{}\f$, \f$x\f$, \f$\beta\f$, and \f$y\f$ but real-valued \f$A\f$.
 * Real-valued strides are in units of <tt>float</tt> while complex-valued
 * strides are in units of <tt>float[2]</tt>.
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
 * \param lda Leading dimension of \c a in units of <tt>float</tt>.
 * \param x Vector to be multiplied.
 * \param incx Stride of vector \c x in units of <tt>float[2]</tt>.
 * \param beta Multiplicative scalar \f$ \beta \f$.
 * \param y Vector to be added to product and to contain result.
 * \param incy Stride of vector \c y in units of <tt>float[2]</tt>.
 *
 * \return Zero on success and a BLAS-like error code otherwise.
 * \see A BLAS reference for for general band storage matrix requirements.
 */
int
suzerain_gbmv_sc(
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$ for complex
 * \f$\alpha{}\f$, \f$x\f$, \f$\beta\f$, and \f$y\f$ but real-valued \f$A\f$.
 * Real-valued strides are in units of <tt>double</tt> while complex-valued
 * strides are in units of <tt>double[2]</tt>.
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
 * \param lda Leading dimension of \c a in units of <tt>double</tt>.
 * \param x Vector to be multiplied.
 * \param incx Stride of vector \c x in units of <tt>double[2]</tt>.
 * \param beta Multiplicative scalar \f$ \beta \f$.
 * \param y Vector to be added to product and to contain result.
 * \param incy Stride of vector \c y in units of <tt>double[2]</tt>.
 *
 * \return Zero on success and a BLAS-like error code otherwise.
 * \see A BLAS reference for for general band storage matrix requirements.
 */
int
suzerain_gbmv_dz(
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
        const int incy);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __SUZERAIN_GBMV_H__ */
