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

#ifndef SUZERAIN_GBMV_H
#define SUZERAIN_GBMV_H

#include <suzerain/complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @file
 * Provides custom, BLAS-like general band matrix-vector operations.
 * Includes mixed real/complex operations and fixed bandwidth kernels.
 */

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} \begin{bmatrix} I & 0 \\ 0 & 1 +
 * \gamma \end{bmatrix} A x + \beta{} y \f$.
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
 * \param gamma Controls factor \f$ 1+\gamma \f$ scaling only the final row.
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
        const int incy,
        const float gamma);

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
        const int incy,
        const double gamma);

/*! \copydoc suzerain_gbmv_s */
int
suzerain_gbmv_c(
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
        const int incy,
        const float gamma);

/*! \copydoc suzerain_gbmv_s */
int
suzerain_gbmv_z(
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
        const int incy,
        const double gamma);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} \begin{bmatrix} I & 0 \\ 0 & 1 +
 * \gamma \end{bmatrix} A x + \beta{} y \f$ for complex \f$\alpha{}\f$,
 * \f$\beta\f$, and \f$y\f$ but real-valued \f$\gamma\f$, \f$A\f$, and \f$x\f$.
 *
 * \copydetails suzerain_gbmv_s
 */
int
suzerain_gbmv_ssc(
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
        const int incy,
        const float gamma);

/*! \copydoc suzerain_gbmv_ssc */
int
suzerain_gbmv_ddz(
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
        const int incy,
        const double gamma);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} \begin{bmatrix} I & 0 \\ 0 & 1 +
 * \gamma \end{bmatrix} A x + \beta{} y \f$ for complex \f$\alpha{}\f$, \f$x\f$,
 * \f$\beta\f$, and \f$y\f$ but real-valued \f$\gamma\f$ and \f$A\f$.
 *
 * \copydetails suzerain_gbmv_s
 */
int
suzerain_gbmv_scc(
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
        const int incy,
        const float gamma);

/*! \copydoc suzerain_gbmv_scc */
int
suzerain_gbmv_dzz(
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
        const int incy,
        const double gamma);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_GBMV_H */
