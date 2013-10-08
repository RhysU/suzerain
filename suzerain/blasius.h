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

// Obtain sqrt(3) in the right fashion depending on language
#ifdef __cplusplus
#include <cmath>
#else
#include <math.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Compute the Blasius coordinate \f$\eta\f$ given \f$y, \mbox{Re}_x\f$.
 *
 * @param y    Streamwise location \f$y\f$.
 * @param Re_x Local Reynolds number \f$\mbox{Re}_x = \frac{u_\infty}{\nu x}\f$.
 *             Generally, values greater than \f$10^6\f$ are physically
 *             invalid as the flow should be turbulent in that regime.
 *
 * @return \f$\eta = y \sqrt{\mbox{Re}_x}\f$.
 */
static inline
double suzerain_blasius_eta(double y, double Re_x)
{
#ifdef __cplusplus
    using namespace std;
#endif
    return y * sqrt(Re_x);
}

/**
 * Tabulated \f$\eta\f$ data from Table 3b of <a
 * href="http://arxiv.org/format/1006.3888v1"> Highly Accurate Solutions of the
 * Blasius and {Falkner-Skan} Boundary Layer Equations via Convergence
 * Acceleration" </a> by B. D. Ganapol (2010).
 */
extern const double suzerain_blasius_ganapol_eta[45];

/**
 * Tabulated \f$f\left(\eta\right)\f$ data from Table 3b of <a
 * href="http://arxiv.org/format/1006.3888v1"> Highly Accurate Solutions of the
 * Blasius and {Falkner-Skan} Boundary Layer Equations via Convergence
 * Acceleration" </a> by B. D. Ganapol (2010).
 *
 * @see \ref suzerain_blasius_ganapol_eta for the matching \f$\eta\f$ values.
 */
extern const double suzerain_blasius_ganapol_f[45];

/**
 * Tabulated \f$f^{\prime}\left(\eta\right)\f$ data from Table 3b of <a
 * href="http://arxiv.org/format/1006.3888v1"> Highly Accurate Solutions of the
 * Blasius and {Falkner-Skan} Boundary Layer Equations via Convergence
 * Acceleration" </a> by B. D. Ganapol (2010).
 *
 * @see \ref suzerain_blasius_ganapol_eta for the matching \f$\eta\f$ values.
 */
extern const double suzerain_blasius_ganapol_fp[45];

/**
 * Tabulated \f$f^{\prime\prime}\left(\eta\right)\f$ data from Table 3b of <a
 * href="http://arxiv.org/format/1006.3888v1"> Highly Accurate Solutions of the
 * Blasius and {Falkner-Skan} Boundary Layer Equations via Convergence
 * Acceleration" </a> by B. D. Ganapol (2010).
 *
 * @see \ref suzerain_blasius_ganapol_eta for the matching \f$\eta\f$ values.
 */
extern const double suzerain_blasius_ganapol_fpp[45];

/**
 * Obtain a Blasius profile spline fit producing nondimensional \f$u /
 * u_\infty\f$ given \f$\eta = y \sqrt{\mbox{Re}_x}\f$.  The data is
 * from \ref suzerain_blasius_ganapol_fp and \ref suzerain_blasius_ganapol_eta.
 * The returned <tt>gsl_spline*</tt> can be interrogated using <a
 * href="http://www.gnu.org/software/gsl/manual/html_node/Higher_002dlevel-Interface.html"
 * the usual routines</a>.  The return value must be subsequently cleaned up
 * using <tt>gsl_spline_free()</tt>.
 *
 * @return On success, a <tt>gsl_spline *</tt> suitable for evaluation using,
 *         for example, <tt>gsl_spline_eval()</tt>.  The return value must be
 *         subsequently cleaned up using <tt>gsl_spline_free()</tt>.
 *         On failure \c NULL is returned.
 */
gsl_spline * suzerain_blasius_u();

/**
 * Obtain a Blasius profile spline fit producing nondimensional \f$v /
 * u_\infty\f$ given \f$\eta = y \sqrt{\mbox{Re}_x}\f$.  The data is computed
 * from \ref suzerain_blasius_ganapol_f, \ref suzerain_blasius_ganapol_fp, and
 * \ref suzerain_blasius_ganapol_eta.  The returned <tt>gsl_spline*</tt> can be
 * interrogated using <a
 * href="http://www.gnu.org/software/gsl/manual/html_node/Higher_002dlevel-Interface.html"
 * the usual routines</a>.  The return value must be subsequently cleaned up
 * using <tt>gsl_spline_free()</tt>.
 *
 * @param Re_x Local Reynolds number \f$\mbox{Re}_x = \frac{u_\infty}{\nu x}\f$.
 *             Generally, values greater than \f$10^6\f$ are physically
 *             invalid as the flow should be turbulent in that regime.
 *
 * @return On success, a <tt>gsl_spline *</tt> suitable for evaluation using,
 *         for example, <tt>gsl_spline_eval()</tt>.  The return value must be
 *         subsequently cleaned up using <tt>gsl_spline_free()</tt>.
 *         On failure \c NULL is returned.
 */
gsl_spline * suzerain_blasius_v(const double Re_x);

/**
 * Obtain a Blasius profile spline fit producing nondimensional kinetic energy
 * \f$\frac{u^2 + v^2}{2 u_\infty^2}\f$ given \f$\eta = y
 * \sqrt{\mbox{Re}_x}\f$.
 *
 * @copydetails suzerain_blasius_v
 */
gsl_spline * suzerain_blasius_ke(const double Re_x);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_BLASIUS_H */
