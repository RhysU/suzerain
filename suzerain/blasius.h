/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2013 Rhys Ulerich
 * Copyright (C) 2013 The PECOS Development Team
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

#ifndef SUZERAIN_BLASIUS_H
#define SUZERAIN_BLASIUS_H

/** @file
 * Presents a Blasius laminar flow profile curve fit.
 */

#include <gsl/gsl_spline.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Obtain a Blasius profile spline fit producing \f$u / u_\infty\f$ given
 * \f$\eta = y \sqrt{\frac{u_\infty}{\nu x}}\f$.  The data is from Table 3b of
 * <a href="http://arxiv.org/format/1006.3888v1"> Highly Accurate Solutions of
 * the Blasius and {Falkner-Skan} Boundary Layer Equations via Convergence
 * Acceleration" </a> by B. D. Ganapol (2010).  The returned
 * <tt>gsl_spline*</tt> can be interrogated using <a
 * href="http://www.gnu.org/software/gsl/manual/html_node/Higher_002dlevel-Interface.html"
 * the usual routines</a>.  The return value must be subsequently cleaned up
 * using <tt>gsl_spline_free()</tt>.
 *
 * @return On success, a <tt>gsl_spline *</tt> suitable for evaluation using,
 *         for example, <tt>gsl_spline_eval()</tt>.  The return value must be
 *         subsequently cleaned up using <tt>gsl_spline_free()</tt>.  On failure
 *         \c NULL is returned.
 */
gsl_spline * suzerain_blasis_spline_fit();

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_BLASIUS_H */
