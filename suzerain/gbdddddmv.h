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
 * gbdddddmv.h: provides BLAS-like general band matrix-vector operations
 * $Id$
 */

#ifndef __SUZERAIN_GBDDDDDMV_H
#define __SUZERAIN_GBDDDDDMV_H

#ifdef __cplusplus
extern "C" {
#endif

/** \file
 * Provides custom, BLAS-like diagonal matrix times general band matrix-vector
 * operations.  Includes mixed real/complex operations and fixed bandwidth
 * kernels.
 */

#include <suzerain/complex.h>

/*!
 * \brief Compute \f$ y \leftarrow{} \left(\alpha_{0} D_0 + \alpha_{1} D_1 +
 * \alpha_{2} D_2 + \alpha_{3} D_3 + \alpha{4} D_4 \right) A x + \beta{} y \f$.
 *
 * Transposes of \f$ A \f$ can be taken using the \c trans parameter.
 *
 * \param trans One of 'N', 'T', or 'C' for no transpose, a transpose,
 *      or a conjugate transpose, respectively.
 * \param n Dimension of all matrices and vectors.
 * \param kl Number of subdiagonals in band storage of \c a.
 * \param ku Number of superdiagonals in band storage of \c a.
 * \param alpha0 Multiplicative scalar \f$ \alpha_0 \f$.
 * \param d0 Storage for diagonal matrix \f$ D_0 \f$.
 * \param ldd0 Leading dimension of \c d0.
 * \param alpha1 Multiplicative scalar \f$ \alpha_1 \f$.
 * \param d1 Contiguous storage for diagonal matrix \f$ D_1 \f$.
 * \param ldd1 Leading dimension of \c d1.
 * \param alpha2 Multiplicative scalar \f$ \alpha_2 \f$.
 * \param d2 Contiguous storage for diagonal matrix \f$ D_2 \f$.
 * \param ldd2 Leading dimension of \c d2.
 * \param alpha3 Multiplicative scalar \f$ \alpha_3 \f$.
 * \param d3 Contiguous storage for diagonal matrix \f$ D_3 \f$.
 * \param ldd3 Leading dimension of \c d3.
 * \param alpha4 Multiplicative scalar \f$ \alpha_4 \f$.
 * \param d4 Contiguous storage for diagonal matrix \f$ D_4 \f$.
 * \param ldd4 Leading dimension of \c d4.
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
suzerain_gbdddddmv_s(
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

/*! \copydoc suzerain_gbdddddmv_s */
int
suzerain_gbdddddmv_d(
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
 * \alpha_{2} D_2 + \alpha_{3} D_{3} + \alpha_{4} D_{4} \right) A x + \beta{} y
 * \f$ for complex \f$\alpha{}\f$, \f$\beta\f$, and \f$y\f$ but real-valued
 * \f$D_0\f$, \f$D_1\f$, \f$D_2\f$, \f$D_3\f$, \f$D_4\f$, \f$A\f$, and \f$x\f$.
 *
 * \copydetails suzerain_gbdddddmv_s
 */
int
suzerain_gbdddddmv_ssc(
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

/*! \copydoc suzerain_gbdddddmv_ssc */
int
suzerain_gbdddddmv_ddz(
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
 * \brief Compute \f$ y \leftarrow{} \left(\alpha_{0} D_0 + \alpha_{1} D_1 +
 * \alpha_{2} D_2 + \alpha_{3} D_{3} + \alpha_{4} D_{4} \right) A x + \beta{} y
 * \f$ for complex \f$\alpha{}\f$, \f$x\f$, \f$\beta\f$, and \f$y\f$ but
 * real-valued \f$D_0\f$, \f$D_1\f$, \f$D_2\f$, \f$D_3\f$, \f$D_4\f$, and
 * \f$A\f$.
 *
 * \copydetails suzerain_gbdddddmv_s
 */
int
suzerain_gbdddddmv_scc(
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

/*! \copydoc suzerain_gbdddddmv_scc */
int
suzerain_gbdddddmv_dzz(
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

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __SUZERAIN_GBDDDDDMV_H */
