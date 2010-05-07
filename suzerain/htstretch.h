/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
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
 * htstretch.c: Routines providing hyperbolic tangent grid stretching
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_HTSTRETCH_H__
#define __SUZERAIN_HTSTRETCH_H__

/* @file
 * Computes quantities related to the common hyperbolic tangent-based grid
 * stretching scheme associated with Marcel Vinokur.
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
 *
 * @param delta The stretching factor \f$\delta\in\left(0,\infty\right)\f$
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
 * @param delta The stretching factor \f$\delta\in\left(0,\infty\right)\f$
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
 * @param delta The stretching factor \f$\delta\in\left(0,\infty\right)\f$
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
 * @param delta The stretching factor \f$\delta\in\left(0,\infty\right)\f$
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

// TODO Document
double
suzerain_htstretch1_pick_delta(const double L,
                               const double crit_x,
                               const double crit_val);

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
 *
 * @param delta The stretching factor \f$\delta\in\left(0,\infty\right)\f$
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
 * @param delta The stretching factor \f$\delta\in\left(0,\infty\right)\f$
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
 * @param delta The stretching factor \f$\delta\in\left(0,\infty\right)\f$
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
 * @param delta The stretching factor \f$\delta\in\left(0,\infty\right)\f$
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

double
suzerain_htstretch_twosided_delta(const double L,
                                  const double xi_crit,
                                  const double u_crit);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __SUZERAIN_HTSTRETCH_H__ */
