/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2011, 2012 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *--------------------------------------------------------------------------
 * gbddmv.h: provides custom, BLAS-like general band matrix-vector operations
 * $Id$
 */

#ifndef __SUZERAIN_GBDDMV_H
#define __SUZERAIN_GBDDMV_H

#ifdef __cplusplus
extern "C" {
#endif

#include <suzerain/complex.h>

/** \file
 * Provides custom, BLAS-like diagonal matrix times general band matrix-vector
 * operations.  Includes mixed real/complex operations and fixed bandwidth
 * kernels.
 */

/*!
 * \brief Compute \f$ y \leftarrow{} \left(\alpha_{0} D_0 + \alpha_{1} D_1
 * \right) A x + \beta{} y \f$.
 *
 * Transposes of \f$ A \f$ can be taken using the \c trans parameter.
 *
 * \param trans One of 'N', 'T', or 'C' for no transpose, a transpose,
 *      or a conjugate transpose, respectively.
 * \param n Dimension of all matrices and vectors.
 * \param kl Number of subdiagonals in band storage of \c a.
 * \param ku Number of superdiagonals in band storage of \c a.
 * \param alpha0 Multiplicative scalar \f$ \alpha_0 \f$.
 * \param d0 Contiguous storage for diagonal matrix \f$ D_0 \f$.
 * \param alpha1 Multiplicative scalar \f$ \alpha_1 \f$.
 * \param d1 Contiguous storage for diagonal matrix \f$ D_1 \f$.
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
suzerain_gbddmv_s(
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
        const int incy);

/*! \copydoc suzerain_gbddmv_s */
int
suzerain_gbddmv_d(
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \left(\alpha_{0} D_0 + \alpha_{1}
 * D_1\right) A x + \beta{} y \f$ for complex \f$\alpha{}\f$, \f$x\f$,
 * \f$\beta\f$, and \f$y\f$ but real-valued \f$D_0\f$, \f$D_1\f$, and \f$A\f$.
 *
 * Transposes of \f$ A \f$ can be taken using the \c trans parameter.
 *
 * \param trans One of 'N', 'T', or 'C' for no transpose, a transpose,
 *      or a conjugate transpose, respectively.
 * \param n Dimension of all matrices and vectors.
 * \param kl Number of subdiagonals in band storage of \c a.
 * \param ku Number of superdiagonals in band storage of \c a.
 * \param alpha0 Multiplicative scalar \f$ \alpha_0 \f$.
 * \param d0 Contiguous storage for diagonal matrix \f$ D_0 \f$.
 * \param alpha1 Multiplicative scalar \f$ \alpha_1 \f$.
 * \param d1 Contiguous storage for diagonal matrix \f$ D_1 \f$.
 * \param a General band storage for matrix \f$ A \f$.
 * \param lda Leading dimension of \c a.
 * \param x Vector to be multiplied.
 * \param incx Stride of vector \c x.
 * \param beta Multiplicative scalar \f$ \beta \f$.
 * \param y Vector to be added to product and to contain result.
 * \param incy Stride of vector \c y.
 *
 * \return Zero on success and a BLAS-like error code otherwise.
 * \see A BLAS reference for for general band storage matrix requirements.
 */
int
suzerain_gbddmv_sc(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const complex_float alpha0,
        const float *d0,
        const complex_float alpha1,
        const float *d1,
        const float *a,
        const int lda,
        const complex_float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy);

/*! \copydoc suzerain_gbddmv_sc */
int
suzerain_gbddmv_dz(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const complex_double alpha0,
        const double *d0,
        const complex_double alpha1,
        const double *d1,
        const double *a,
        const int lda,
        const complex_double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __SUZERAIN_GBDDMV_H */
