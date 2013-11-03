/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2011, 2012, 2013 Rhys Ulerich
 * Copyright (C) 2012, 2013 The PECOS Development Team
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

#ifndef SUZERAIN_BSPLINE_H
#define SUZERAIN_BSPLINE_H

/** @file
 * Higher-level logic build atop the GSL B-spline routines
 */

#include <gsl/gsl_bspline.h>

#include <suzerain/complex.h>
#include <suzerain/pre_gsl.h>

/** @file
 * Provides higher-level logic built atop the <a
 * href="http://www.gnu.org/software/gsl/">GNU Scientific Library</a> (GSL) <a
 * href="http://www.gnu.org/software/gsl/manual/html_node/Basis-Splines.html">
 * B-spline routines</a>.  Many of these routines can (and should) be
 * refactored and submitted upstream to the GSL.
 *
 * Multiple threads may call these routines simultaneously provided that each
 * thread has its own workspace instance.  The behavior is undefined if two
 * threads simultaneously share a workspace.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Evaluate a function and its derivatives based upon a linear combination of
 * the basis functions.  At each \c point in \c points, evaluate the function
 * and its derivatives determined from a linear combination of the supplied
 * real-valued \c coeffs times the B-spline basis functions at \c point.
 * Derivatives \c 0 through \c nderiv (inclusive) are computed and stored
 * as columns within \c values.
 * That is, \f[
 *  \mbox{values}\left[i + k\,\mbox{ldvalues}\right] = \frac{d^k}{dx^k}
 *  \sum_{j=0}^{\mbox{ndof}} \mbox{coeffs}_{j}\,B_{j}(\mbox{points}_{i})
 * \f]
 * for \f$ i\in\left\{0,\ldots,\mbox{npoints}\right\} \f$,
 *\f$ k\in\left\{0,\ldots,\mbox{nderiv}\right\} \f$,
 * \f$\mbox{ndof} = \f$ <code>w->n</code>, and B-spline basis
 * functions \f$ B_{j} \f$.  If only a single derivative is desired, passing \c
 * 0 to \c ldvalues will cause only that single derivative to be written in the
 * first column of \c values.
 *
 * \note It is more efficient to compute a function and its derivatives
 * simultaneously than to request each derivative separately.  This is due to
 * the recurrence relationship used to compute B-spline derivatives.
 *
 * @param[in] nderiv Maximum requested derivative.  This may be higher
 *      than the number of derivatives requested in suzerain_bspline_alloc().
 * @param[in] coeffs Real-valued expansion coefficients for a function in
 *      terms of the B-spline basis.  Must be of length <code>w->n</code>.
 * @param[in] npoints Number of evaluation points.
 * @param[in] points Points at which to evaluate the function.
 * @param[out] values Matrix of real values resulting from evaluating the
 *      function and its derivatives.  Matrix dimensions are <code>w->n</code>
 *      by \c nderiv if \c ldvalues >= <code>w->n</code>.  If \c ldvalues
 *      is zero, only a single column is returned in \c values.
 * @param[in] ldvalues Leading dimension of the output matrix \c values
 *            measured as real-valued strides.
 * @param[in] dB Temporary storage to use of size <tt>w->k</tt> by
 *            no less than <tt>nderiv + 1</tt>.
 * @param[in] w Workspace to use.
 * @param[in] dw Workspace to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 * @see suzerain_bspline_linear_combination_complex() for a way to evaluate a
 *      function and its derivatives when the coefficients are complex-valued.
 */
int
suzerain_bspline_linear_combination(
    const size_t nderiv,
    const double * coeffs,
    const size_t npoints,
    const double * points,
    double * values,
    const size_t ldvalues,
    gsl_matrix *dB,
    gsl_bspline_workspace *w,
    gsl_bspline_deriv_workspace *dw);

/**
 * Evaluate a function and its derivatives based upon a linear combination of
 * the basis functions.  At each \c point in \c points, evaluate the function
 * and its derivatives determined from a linear combination of the supplied
 * complex-valued \c coeffs times the B-spline basis functions at \c
 * point.  Derivatives \c 0 through \c nderiv (inclusive) are computed and
 * stored as columns within \c values.
 * That is, \f[
 *  \mbox{values}\left[i + k\,\mbox{ldvalues}\right] = \frac{d^k}{dx^k}
 *  \sum_{j=0}^{\mbox{ndof}} \mbox{coeffs}_{j}\,B_{j}(\mbox{points}_{i})
 * \f]
 * for \f$ i\in\left\{0,\ldots,\mbox{npoints}\right\} \f$,
 *\f$ k\in\left\{0,\ldots,\mbox{nderiv}\right\} \f$,
 * \f$\mbox{ndof} = \f$ <code>w->n</code>, and B-spline basis
 * functions \f$ B_{j} \f$.  If only a single derivative is desired, passing \c
 * 0 to \c ldvalues will cause only that single derivative to be written in the
 * first column of \c values.
 *
 * \note It is more efficient to compute a function and its derivatives
 * simultaneously than to request each derivative separately.  This is due to
 * the recurrence relationship used to compute B-spline derivatives.
 *
 * @param[in] nderiv Maximum requested derivative.  This may be higher
 *      than the number of derivatives requested in suzerain_bspline_alloc().
 * @param[in] coeffs Complex-valued expansion coefficients for a function
 *      in terms of the B-spline basis.  Must be of length <code>w->n</code>.
 * @param[in] npoints Number of evaluation points.
 * @param[in] points Points at which to evaluate the function.
 * @param[out] values Matrix of complex values resulting from evaluating the
 *      function and its derivatives.  Matrix dimensions are <code>w->n</code>
 *      by \c nderiv if \c ldvalues >= <code>code</code>.  If \c ldvalues
 *      is zero, only a single column is returned in \c values.
 * @param[in] ldvalues Leading dimension of the output matrix \c values
 *            measured as complex-valued strides.
 * @param[in] dB Temporary storage to use of size <tt>w->k</tt> by
 *            no less than <tt>nderiv + 1</tt>.
 * @param[in] w Workspace to use.
 * @param[in] dw Workspace to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 * @see suzerain_bspline_linear_combination() for a way to evaluate a function
 *      and its derivatives when the coeffs are real-valued.
 */
int
suzerain_bspline_linear_combination_complex(
    const size_t nderiv,
    const complex_double *coeffs,
    const size_t npoints,
    const double * points,
    complex_double *values,
    const size_t ldvalues,
    gsl_matrix *dB,
    gsl_bspline_workspace *w,
    gsl_bspline_deriv_workspace *dw);

/**
 * Using the function specified by B-spline coefficients \c coeffs
 * for basis \c w and \c dw, find the \c location at which the
 * function's \c nderiv derivative crosses \c value within the region
 * <tt>[lower, upper]</tt>.  The crossing must be known <i>a priori</i>
 * to exist by the intermediate value theorem.  If multiple crossings
 * exist, it is undefined which is returned.  Convergence is reached
 * whenever <tt>[lower, upper]</tt> satisfy the criteria set forth in <a
 * href="http://www.gnu.org/software/gsl/manual/html_node/Search-Stopping-Parameters.html">gsl_root_test_interval</a>
 * using \c epsabs and \c epsrel.
 *
 * @param[in] nderiv Derivative of interest  This may be higher than the number
 *     of derivatives requested in suzerain_bspline_alloc().
 * @param[in] coeffs Real-valued expansion coefficients for the function with
 *     respect to B-spline basis.  Must be of length <code>w->n</code>.
 * @param[in] value  The function value at which a crossing is sought.
 * @param[in,out] lower On entry, the lower edge of the range in which
 *     a crossing is sought.  On exit, the active lower bound in the search.
 * @param[in,out] upper On entry, the upper edge of the range in which
 *     a crossing is sought.  On exit, the active upper bound in the search.
 * @param[in] maxiter Maximum number of iterations to attempt.
 * @param[in] epsabs Absolute error bound used as a stopping criterion.
 * @param[in] epsrel Relative error bound used as a stopping criterion.
 * @param[out] location Location at which a crossing is found.
 * @param[in] dB Temporary storage to use of size <tt>w->k</tt> by
 *            no less than <tt>nderiv + 1</tt>.
 * @param[in] w Workspace to use.
 * @param[in] dw Workspace to use.
 *
 * @return ::SUZERAIN_SUCCESS on success and returns the answer in
 * <code>*location</code>.  If convergence is slow, ::SUZERAIN_CONTINUE
 * will be returned and the method may be re-invoked with the now updated
 * \c lower and \c upper range details to proceed.  On recoverable
 * error sets <code>*location</code> to be <tt>NaN</tt> and returns
 * one of #suzerain_error_status.  On unrecoverable error, additionally
 * calls suzerain_error().
 */
int
suzerain_bspline_crossing(
    const size_t nderiv,
    const double * coeffs,
    const double value,
    double * lower,
    double * upper,
    const size_t maxiter,
    const double epsabs,
    const double epsrel,
    double * location,
    gsl_matrix *dB,
    gsl_bspline_workspace *w,
    gsl_bspline_deriv_workspace *dw);

/**
 * Compute the coefficients \f$ \gamma_{i} \f$ for <code>0 <= i < w->n</code>
 * such that \f$ \vec{\gamma}\cdot\vec{\beta} =
 * \int_{x_\mbox{lo}}^{x_\mbox{hi}} \sum_{i} \beta_{i}
 * B_{i}^{(\mbox{nderiv})}(x) \, dx\f$.  Integration bounds outside the support
 * of the basis as well as a reversed bounds are handled correctly.  Ranges
 * like <tt>[-DBL_MAX,DBL_MAX]</tt> are shorthand effectively indicating
 * integration over the support of the B-spline basis.
 *
 * @param[in]  nderiv The derivative to integrate.
 * @param[out] coeffs Real-valued coefficients \f$ \gamma_{i} \f$.
 * @param[in]  inc Stride between elements of \c coeffs.
 * @param[in]  lo The lower integration bound \f$x_\mbox{lo}$.
 * @param[in]  hi The upper integration bound \f$x_\mbox{hi}$.
 * @param[in]  dB Temporary storage to use of size <tt>w->k</tt> by
 *             no less than <tt>nderiv + 1</tt>.
 * @param[in]  w Workspace to use (which sets the integration bounds).
 * @param[in]  dw Workspace to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_bspline_integration_coefficients(
    const size_t nderiv,
    double * coeffs,
    size_t inc,
    double lo,
    double hi,
    gsl_matrix *dB,
    gsl_bspline_workspace *w,
    gsl_bspline_deriv_workspace *dw);

/**
 * Compute the "distance" between two B-spline workspaces.  Distance is "huge"
 * if any of the order, number of degrees of freedom, or number of knots
 * differ.  When all those criteria match the distance becomes the maximum
 * absolute difference between the knot vectors.
 *
 * @param a First B-spline workspace
 * @param b Second B-spline workspace
 *
 * @return A measure of the distance between the two workspaces
 */
double
suzerain_bspline_distance(
    const gsl_bspline_workspace *a,
    const gsl_bspline_workspace *b);

/**
 * Constant used to define distinct B-spline workspaces per
 * \ref suzerain_bspline_distance in the presence of floating point error.
 */
extern const double suzerain_bspline_distance_distinct;

// FIXME Implement per Redmine ticket #1591
#ifndef SUZERAIN_PARSED_BY_DOXYGEN
/**
 * Construct a general banded matrix \f$A\f$ that computes the vector
 * \f$\frac{d^k}{dx^k} \alpha_{i} f(x_{i})\f$ (no summation implied) when it
 * multiplies the coefficient vector \f$\beta_{j}\f$ where
 * \f$f(x)=\sum_{j}\beta_{j}B_{j}(x)\f$.  Such banded matrices are useful, for
 * example, when repeatedly evaluating B-splines on a particular set of points.
 *
 * When \c a is \c NULL, the routine computes the minimum acceptable \c kl, \c
 * ku, and \c lda necessary for storing the coefficients.  When \c a is
 * <tt>non-NULL</tt>, the banded matrix coefficients are computed and stored
 * according to whatever \c kl, \c ku, and \c lda is provided.  All other
 * values within the banded storage will be zeroed.
 *
 * The number of rows and columns is stored in \c m and \c n whenever these
 * variables are <tt>non-NULL</tt> on entry.  Of course, \c m is always equal
 * to \c npoints while \c n is always equal to the number of degrees of freedom
 * within B-spline workspace \c w.  These arguments are provided for ease of
 * interoperation with the BLAS.
 *
 * @param[in]    nderiv Derivative to evaluate.  May be zero.
 * @param[in]    npoints  Number of points on which to evaluate basis.
 * @param[in]    points   Points \f$x_{i}\f$ at which to evaluate basis.
 * @param[in]    weights  Weights \f$\alpha_{i}\f$ used to weight \f$x_{i}\f$.
 * @param[out]   m Number of rows in resulting banded matrix \f$A\f$.
 * @param[out]   n Number of columns in resulting banded matrix \f$A\f$.
 * @param[inout] kl Number of subdiagonals in the matrix \f$A\f$
 * @param[inout] ku Number of superdiagonals in the matrix \f$A\f$
 * @param[inout] a  Coefficients of \f$A\f$ stored in BLAS-compatible
 *                  general band matrix format.
 * @param[inout] lda Leading dimension between columns in \c a.
 * @param[in] dB Temporary storage to use of size <tt>w->k</tt> by
 *            no less than <tt>nderiv + 1</tt>.
 * @param[in]    w  Workspace to use for B-spline evaluation.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_bspline_linear_combination_matrix(
    int nderiv,
    int npoints,
    const double *points,
    const double *weights,
    double *m,
    double *n,
    double *kl,
    double *ku,
    double *a,
    double *lda,
    gsl_matrix *dB,
    gsl_bspline_workspace *w,
    gsl_bspline_deriv_workspace dw);
#endif /* SUZERAIN_PARSED_BY_DOXYGEN */

/**
 * Find the minimum distance to Greville abscissae <tt>i-1</tt> and
 * <tt>i+1</tt> from abscissa <tt>i</tt>.  Accounts for the edge cases
 * when <tt>i == 0</tt> or <tt>i == N-1</tt> for a basis with <tt>N</tt>
 * degrees of freedom.
 *
 * @param[in] i Abscissa of interest.
 * @param[in] w Workspace to use.
 *
 * @return Minimum distance to the next nearest abscissae.
 */
double
suzerain_bspline_spacing_greville_abscissae(
    size_t i,
    gsl_bspline_workspace *w);

/**
 * Find the minimum distance between the breakpoints comprising
 * the support of the <tt>i</tt>th B-spline basis function.
 *
 * @param[in] i Degree of freedom of interest.
 * @param[in] w Workspace to use.
 *
 * @return Minimum distance between breakpoints for basis function <tt>i</tt>.
 *         Equivalently, this is the minimum distance between knots
 *         once repeated knots are eliminated.
 */
double
suzerain_bspline_spacing_breakpoints(
    size_t i,
    gsl_bspline_workspace *w);

/**
 * Compute the grid spacing scaling factor \f$C^{(i)}\f$ necessary to make
 * \f$\left(\frac{\pi}{C\Delta{}x}\right)^i\f$ a good estimate of the maximum
 * eigenvalue for the \f$i\f$th Greville abscissae-based, collocation
 * derivative operator when using two-sided hyperbolic tangent stretching.
 * \f$\Delta{}x\f$ must measure the distance between adjacent <em>collocation
 * points</em>.  Beware that spectral-like cases, those where <tt>N <=
 * 5*k</tt>, will likely give untrustworthy results.  See the Suzerain model
 * document for more details.
 *
 * @param[in]  nderiv  The desired derivative order.
 *                     Currently only the values \c 1 and \c 2 are acceptable.
 * @param[in]  k       The B-spline order following GSL conventions.
 *                     For example, <tt>k=4</tt> denotes piecewise cubics.
 * @param[in]  htdelta The nonnegative, two-sided grid stretching parameter.
 * @param[in]  N       The number of degrees of freedom in the basis.
 * @param[out] C       The computed grid spacing scaling factor \f$C^{(i)}\f$.
 * @param[out] Clow    If not NULL, a lower error bound on \c C.  That is,
 *                     a reasonable lower bound on what \f$C^{(i)}\f$ might
 *                     truly be.
 * @param[out] Chigh   If not NULL, an upper error bound on \c C.  That is,
 *                     a reasonable upper bound on what \f$C^{(i)}\f$ might
 *                     truly be.
 *
 * @see suzerain_htstretch2() for more information on the grid stretching.
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_bspline_htstretch2_evdeltascale_greville_abscissae(
    const int nderiv,
    const int k,
    const double htdelta,
    const int N,
    double * const C,
    double * const Clow,
    double * const Chigh);

/**
 * Compute the grid spacing scaling factor \f$C^{(i)}\f$ necessary to make
 * \f$\left(\frac{\pi}{C\Delta{}x}\right)^i\f$ a good estimate of the maximum
 * eigenvalue for the \f$i\f$th Greville abscissae-based, collocation
 * derivative operator when using one-sided hyperbolic tangent stretching.
 * \f$\Delta{}x\f$ must measure the distance between adjacent <em>collocation
 * points</em>.  Beware that spectral-like cases, those where <tt>N <=
 * 5*k</tt>, will likely give untrustworthy results.  See the Suzerain model
 * document for more details.
 *
 * @param[in]  nderiv  The desired derivative order.
 *                     Currently only the values \c 1 and \c 2 are acceptable.
 * @param[in]  k       The B-spline order following GSL conventions.
 *                     For example, <tt>k=4</tt> denotes piecewise cubics.
 * @param[in]  htdelta The nonnegative, one-sided grid stretching parameter.
 * @param[in]  N       The number of degrees of freedom in the basis.
 * @param[out] C       The computed grid spacing scaling factor \f$C^{(i)}\f$.
 * @param[out] Clow    If not NULL, a lower error bound on \c C.  That is,
 *                     a reasonable lower bound on what \f$C^{(i)}\f$ might
 *                     truly be.
 * @param[out] Chigh   If not NULL, an upper error bound on \c C.  That is,
 *                     a reasonable upper bound on what \f$C^{(i)}\f$ might
 *                     truly be.
 *
 * @see suzerain_htstretch1() for more information on the grid stretching.
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_bspline_htstretch1_evdeltascale_greville_abscissae(
    const int nderiv,
    const int k,
    const double htdelta,
    const int N,
    double * const C,
    double * const Clow,
    double * const Chigh);

/**
 * Compute the grid spacing scaling factor \f$C^{(i)}\f$ necessary to make
 * \f$\left(\frac{\pi}{C\Delta{}x}\right)^i\f$ a good estimate of the maximum
 * eigenvalue for the \f$i\f$th Greville abscissae-based, collocation
 * derivative operator when using two-sided hyperbolic tangent stretching.
 * \f$\Delta{}x\f$ must measure the distance between adjacent
 * <em>breakpoints</em>.  Beware that spectral-like cases, those where <tt>N <=
 * 5*k</tt>, will likely give untrustworthy results.  See the Suzerain model
 * document for more details.
 *
 * @param[in]  nderiv  The desired derivative order.
 *                     Currently only the values \c 1 and \c 2 are acceptable.
 * @param[in]  k       The B-spline order following GSL conventions.
 *                     For example, <tt>k=4</tt> denotes piecewise cubics.
 * @param[in]  htdelta The nonnegative, two-sided grid stretching parameter.
 * @param[in]  N       The number of degrees of freedom in the basis.
 * @param[out] C       The computed grid spacing scaling factor \f$C^{(i)}\f$.
 * @param[out] Clow    If not NULL, a lower error bound on \c C.  That is,
 *                     a reasonable lower bound on what \f$C^{(i)}\f$ might
 *                     truly be.
 * @param[out] Chigh   If not NULL, an upper error bound on \c C.  That is,
 *                     a reasonable upper bound on what \f$C^{(i)}\f$ might
 *                     truly be.
 *
 * @see suzerain_htstretch2() for more information on the grid stretching.
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_bspline_htstretch2_evdeltascale_breakpoints(
    const int nderiv,
    const int k,
    const double htdelta,
    const int N,
    double * const C,
    double * const Clow,
    double * const Chigh);

/**
 * Compute the grid spacing scaling factor \f$C^{(i)}\f$ necessary to make
 * \f$\left(\frac{\pi}{C\Delta{}x}\right)^i\f$ a good estimate of the maximum
 * eigenvalue for the \f$i\f$th Greville abscissae-based, collocation
 * derivative operator when using one-sided hyperbolic tangent stretching.
 * \f$\Delta{}x\f$ must measure the distance between adjacent
 * <em>breakpoints</em>.  Beware that spectral-like cases, those where <tt>N <=
 * 5*k</tt>, will likely give untrustworthy results.  See the Suzerain model
 * document for more details.
 *
 * @param[in]  nderiv  The desired derivative order.
 *                     Currently only the values \c 1 and \c 2 are acceptable.
 * @param[in]  k       The B-spline order following GSL conventions.
 *                     For example, <tt>k=4</tt> denotes piecewise cubics.
 * @param[in]  htdelta The nonnegative, one-sided grid stretching parameter.
 * @param[in]  N       The number of degrees of freedom in the basis.
 * @param[out] C       The computed grid spacing scaling factor \f$C^{(i)}\f$.
 * @param[out] Clow    If not NULL, a lower error bound on \c C.  That is,
 *                     a reasonable lower bound on what \f$C^{(i)}\f$ might
 *                     truly be.
 * @param[out] Chigh   If not NULL, an upper error bound on \c C.  That is,
 *                     a reasonable upper bound on what \f$C^{(i)}\f$ might
 *                     truly be.
 *
 * @see suzerain_htstretch1() for more information on the grid stretching.
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_bspline_htstretch1_evdeltascale_breakpoints(
    const int nderiv,
    const int k,
    const double htdelta,
    const int N,
    double * const C,
    double * const Clow,
    double * const Chigh);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // SUZERAIN_BSPLINE_H
