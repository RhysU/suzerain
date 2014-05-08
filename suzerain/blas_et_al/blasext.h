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

#ifndef SUZERAIN_BLAS_ET_AL_BLASEXT_H
#define SUZERAIN_BLAS_ET_AL_BLASEXT_H

/*!\file
 * BLAS-like extensions built atop the BLAS, on vendor-specific BLAS-like
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
 */

#include <suzerain/common.h>
#include <suzerain/complex.h>

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_BLAS_ET_AL_BLASEXT_H */
