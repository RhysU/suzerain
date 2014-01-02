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

#ifndef SUZERAIN_HTSTRETCH_H
#define SUZERAIN_HTSTRETCH_H

/** @file
 * Computes quantities related to the common hyperbolic tangent-based grid
 * stretching scheme associated with Marcel Vinokur's 1983 Journal of
 * Computational Physics paper.
 *
 * @see The <a href="http://dx.doi.org/10.1016/0021-9991(83)90065-7">original
 * paper</a> for more details.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Compute the one-sided hyperbolic tangent stretching function.
 * The function has domain \f$\left[0,L\right]\f$, range \f$\left[0,1\right]\f$
 * and is given by
 * \f[
 *      s\left(x\right) = 1 + \frac{
 *          \mbox{tanh}\left[\delta\left(x/L-1\right)\right]
 *      }{
 *          \mbox{tanh}\left[\delta\right]
 *      }
 * \f]
 * with the limiting \f$\delta\to{}0\f$ result being simply \f$x/L\f$.
 *
 *
 * @param delta The stretching factor \f$\delta\in\left[0,\infty\right)\f$
 * @param L     The domain endpoint \f$L\f$
 * @param x     The desired evaluation point \f$x\in\left[0,L\right]\f$
 *
 * @return The one-sided hyperbolic tangent stretching function evaluated
 *         for the given parameters.
 *
 * @see CFD Online's <a
 * href="http://www.cfd-online.com/Wiki/Structured_mesh_generation">structured
 * mesh generation</a> topic for more information.
 */
double
suzerain_htstretch1(const double delta,
                    const double L,
                    const double x);

/**
 * Compute the partial derivative of the one-sided hyperbolic tangent
 * stretching function with respect to the stretching factor \f$\delta\f$.
 *
 * @param delta The stretching factor \f$\delta\in\left[0,\infty\right)\f$
 * @param L     The domain endpoint \f$L\f$
 * @param x     The desired evaluation point \f$x\in\left[0,L\right]\f$
 *
 * @return The one-sided hyperbolic tangent stretching function partial
 *         derivative with respect to the stretching factor evaluated
 *         for the given parameters.
 *
 * @see suzerain_htstretch1() for the form of the function.
 */
double
suzerain_htstretch1_ddelta(const double delta,
                           const double L,
                           const double x);

/**
 * Compute the partial derivative of the one-sided hyperbolic tangent
 * stretching function with respect to the domain size \f$L\f$.
 *
 * @param delta The stretching factor \f$\delta\in\left[0,\infty\right)\f$
 * @param L     The domain endpoint \f$L\f$
 * @param x     The desired evaluation point \f$x\in\left[0,L\right]\f$
 *
 * @return The one-sided hyperbolic tangent stretching function partial
 *         derivative with respect to the domain size evaluated
 *         for the given parameters.
 *
 * @see suzerain_htstretch1() for the form of the function.
 */
double
suzerain_htstretch1_dL(const double delta,
                       const double L,
                       const double x);

/**
 * Compute the partial derivative of the one-sided hyperbolic tangent
 * stretching function with respect to the evaluation point \f$x\f$.
 *
 * @param delta The stretching factor \f$\delta\in\left[0,\infty\right)\f$
 * @param L     The domain endpoint \f$L\f$
 * @param x     The desired evaluation point \f$x\in\left[0,L\right]\f$
 *
 * @return The one-sided hyperbolic tangent stretching function partial
 *         derivative with respect to the evaluation point
 *         for the given parameters.
 *
 * @see suzerain_htstretch1() for the form of the function.
 */
double
suzerain_htstretch1_dx(const double delta,
                       const double L,
                       const double x);

/**
 * Find a stretching factor \f$\delta\f$ so that
 * \f$s\left(x_{\mbox{crit}}\right) \approx{} s_{\mbox{crit}}\f$.
 *
 * @param[in]  L       The domain endpoint \f$L\f$
 * @param[in]  x_crit  The desired evaluation point
 *                     \f$x_{\mbox{crit}}\in\left[0,L\right]\f$
 * @param[in]  s_crit  The desired evaluation value
 *                     \f$s_{\mbox{crit}}\in\left[0,1\right]\f$
 * @param[in]  epsabs  The absolute tolerance criteria for stopping the
 *                     iterative process per <tt>gsl_root_test_residual()</tt>
 * @param[in]  maxiter The maximum number of iterations to attempt
 * @param[out] delta   An approximation to the necessary stretching factor
 *                     \f$\delta\f$ to solve the problem.
 *                     When ::SUZERAIN_SUCCESS is returned, this meets the
 *                     given tolerance.  When ::SUZERAIN_EMAXITER is returned,
 *                     this is a best guess but does not meet the tolerance.
 *                     For any other return value, \c delta is undefined.
 *
 * @return ::SUZERAIN_SUCCESS (zero) on success.  Returns ::SUZERAIN_EMAXITER
 *         if the maximum number of iterations was reached without meeting the
 *         given tolerance requirement.  On any other error, calls
 *         suzerain_error() and returns one of #suzerain_error_status.
 *
 * @see suzerain_htstretch1() for the form of \f$s\left(x\right)\f$.
 * @see GSL's <a href="http://www.gnu.org/software/gsl/manual/html_node/Search-Stopping-Parameters.html">Search Stopping Parameters</a>
 * documentation on <tt>gsl_root_test_residual()</tt> for more information
 * on the convergence criterion.
 */
int
suzerain_htstretch1_find_delta(const double L,
                               const double x_crit,
                               const double s_crit,
                               const double epsabs,
                               const int maxiter,
                               double *delta);

/**
 * Compute the two-sided hyperbolic tangent stretching function.
 * The function has domain \f$\left[0,L\right]\f$, range \f$\left[0,1\right]\f$
 * and is given by
 * \f[
 *      u\left(x\right) = \frac{1}{2}\left(1 + \frac{
 *          \mbox{tanh}\left[\delta\left(x/L-1/2\right)\right]
 *      }{
 *          \mbox{tanh}\left[\delta/2\right]
 *      }\right)
 * \f]
 * with the limiting \f$\delta\to{}0\f$ result being simply \f$x/L\f$.
 *
 * @param delta The stretching factor \f$\delta\in\left[0,\infty\right)\f$
 * @param L     The domain endpoint \f$L\f$
 * @param x     The desired evaluation point \f$x\in\left[0,L\right]\f$
 *
 * @return The two-sided hyperbolic tangent stretching function evaluated
 *         for the given parameters.
 *
 * @see CFD Online's <a
 * href="http://www.cfd-online.com/Wiki/Structured_mesh_generation">structured
 * mesh generation</a> topic for more information.
 */
double
suzerain_htstretch2(const double delta,
                    const double L,
                    const double x);

/**
 * Compute the partial derivative of the two-sided hyperbolic tangent
 * stretching function with respect to the stretching factor \f$\delta\f$.
 *
 * @param delta The stretching factor \f$\delta\in\left[0,\infty\right)\f$
 * @param L     The domain endpoint \f$L\f$
 * @param x     The desired evaluation point \f$x\in\left[0,L\right]\f$
 *
 * @return The two-sided hyperbolic tangent stretching function partial
 *         derivative with respect to the stretching factor evaluated
 *         for the given parameters.
 *
 * @see suzerain_htstretch2() for the form of the function.
 */
double
suzerain_htstretch2_ddelta(const double delta,
                           const double L,
                           const double x);

/**
 * Compute the partial derivative of the two-sided hyperbolic tangent
 * stretching function with respect to the domain size \f$L\f$.
 *
 * @param delta The stretching factor \f$\delta\in\left[0,\infty\right)\f$
 * @param L     The domain endpoint \f$L\f$
 * @param x     The desired evaluation point \f$x\in\left[0,L\right]\f$
 *
 * @return The two-sided hyperbolic tangent stretching function partial
 *         derivative with respect to the domain size evaluated
 *         for the given parameters.
 *
 * @see suzerain_htstretch2() for the form of the function.
 */
double
suzerain_htstretch2_dL(const double delta,
                       const double L,
                       const double x);

/**
 * Compute the partial derivative of the two-sided hyperbolic tangent
 * stretching function with respect to the evaluation point \f$x\f$.
 *
 * @param delta The stretching factor \f$\delta\in\left[0,\infty\right)\f$
 * @param L     The domain endpoint \f$L\f$
 * @param x     The desired evaluation point \f$x\in\left[0,L\right]\f$
 *
 * @return The two-sided hyperbolic tangent stretching function partial
 *         derivative with respect to the evaluation point
 *         for the given parameters.
 *
 * @see suzerain_htstretch2() for the form of the function.
 */
double
suzerain_htstretch2_dx(const double delta,
                       const double L,
                       const double x);

/**
 * Find a stretching factor \f$\delta\f$ so that
 * \f$u\left(x_{\mbox{crit}}\right) \approx{} u_{\mbox{crit}}\f$.
 *
 * @param[in]  L       The domain endpoint \f$L\f$
 * @param[in]  x_crit  The desired evaluation point
 *                     \f$x_{\mbox{crit}}\in\left[0,L\right]\f$
 * @param[in]  u_crit  The desired evaluation value
 *                     \f$u_{\mbox{crit}}\in\left[0,1\right]\f$
 * @param[in]  epsabs  The absolute tolerance criteria for stopping the
 *                     iterative process per <tt>gsl_root_test_residual()</tt>
 * @param[in]  maxiter The maximum number of iterations to attempt
 * @param[out] delta   An approximation to the necessary stretching factor
 *                     \f$\delta\f$ to solve the problem.
 *                     When ::SUZERAIN_SUCCESS is returned, this meets the
 *                     given tolerance.  When ::SUZERAIN_EMAXITER is returned,
 *                     this is a best guess but does not meet the tolerance.
 *                     For any other return value, \c delta is undefined.
 *
 * @return ::SUZERAIN_SUCCESS (zero) on success.  Returns ::SUZERAIN_EMAXITER
 *         if the maximum number of iterations was reached without meeting the
 *         given tolerance requirement.  On any other error, calls
 *         suzerain_error() and returns one of #suzerain_error_status.
 *
 * @see suzerain_htstretch2() for the form of \f$u\left(x\right)\f$.
 * @see GSL's <a href="http://www.gnu.org/software/gsl/manual/html_node/Search-Stopping-Parameters.html">Search Stopping Parameters</a>
 * documentation on <tt>gsl_root_test_residual()</tt> for more information
 * on the convergence criterion.
 */
int
suzerain_htstretch2_find_delta(const double L,
                               const double x_crit,
                               const double u_crit,
                               const double epsabs,
                               const int maxiter,
                               double *delta);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_HTSTRETCH_H */
