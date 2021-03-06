/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2013-2014 Rhys Ulerich
 * Copyright (C) 2013-2014 The PECOS Development Team
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
 * Presents a Blasius laminar flow profile curve fit.  A nice discussion of the
 * solution and its properties appears in section 9.3 of <a
 * href="http://www.worldcat.org/oclc/844779335">Fluid Mechanics</a> by Kundu,
 * Cohen, and Dowling (ISBN 9780123821003).
 *
 * Throughout these routines \f$x\f$ is the streamwise location downstream from
 * a flat plate edge and \f$y\f$ is the wall-normal distance from the plate.
 * The local Reynolds number \f$\mbox{Re}_x = \frac{u_\infty}{\nu x}\f$ scales
 * the similarity coordinate \f$\eta = y \sqrt{\mbox{Re}_x}\f$.  Generally,
 * values greater than \f$10^6\f$ are physically invalid as the flow should be
 * turbulent in that regime.
 */

#ifdef __cplusplus
#include <cstddef>
#else
#include <stddef.h>
#endif

#include <gsl/gsl_spline.h>

#ifdef __cplusplus
extern "C" {
#endif

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
 * The number of contiguous data points accessible via \ref
 * suzerain_blasius_extended_eta, \ref suzerain_blasius_extended_f, \ref
 * suzerain_blasius_extended_fp, and \ref suzerain_blasius_extended_fpp.
 */
#ifdef __cplusplus
extern const std::size_t suzerain_blasius_extended_size;
#else
extern const      size_t suzerain_blasius_extended_size;
#endif

/**
 * Generated \f$\eta\f$ data suitable for use outside well beyond the classical
 * \f$\delta_{99}\f$ boundary layer thickness which matches \ref
 * suzerain_blasius_ganapol_eta to eight digits.
 *
 * @see \ref suzerain_blasius_extended_size for the number of available points.
 */
extern const double * const suzerain_blasius_extended_eta;

/**
 * Generated \f$f\left(\eta\right)\f$ data suitable for use outside well beyond
 * the classical \f$\delta_{99}\f$ boundary layer thickness which matches \ref
 * suzerain_blasius_ganapol_f to eight digits.
 *
 * @see \ref suzerain_blasius_extended_size for the number of available points.
 * @see \ref suzerain_blasius_extended_eta for the matching \f$\eta\f$ values.
 */
extern const double * const suzerain_blasius_extended_f;

/**
 * Generated \f$f^{\prime}\left(\eta\right)\f$ data suitable for use outside
 * well beyond the classical \f$\delta_{99}\f$ boundary layer thickness which
 * matches \ref suzerain_blasius_ganapol_fp to eight digits.
 *
 * @see \ref suzerain_blasius_extended_size for the number of available points.
 * @see \ref suzerain_blasius_extended_eta for the matching \f$\eta\f$ values.
 */
extern const double * const suzerain_blasius_extended_fp;

/**
 * Generated \f$f^{\prime\prime}\left(\eta\right)\f$ data suitable for use
 * outside well beyond the classical \f$\delta_{99}\f$ boundary layer thickness
 * which matches \ref suzerain_blasius_ganapol_fpp to eight digits.
 *
 * @see \ref suzerain_blasius_extended_size for the number of available points.
 * @see \ref suzerain_blasius_extended_eta for the matching \f$\eta\f$ values.
 */
extern const double * const suzerain_blasius_extended_fpp;

/**
 * Obtain a Blasius profile fit producing nondimensional \f$u / u_\infty\f$
 * given \f$y = \eta / \sqrt{\mbox{Re}_x}\f$.  The data is from \ref
 * suzerain_blasius_extended_fp and \ref suzerain_blasius_extended_eta.
 *
 * The returned <tt>gsl_spline*</tt> can be interrogated using <a
 * href="http://www.gnu.org/software/gsl/manual/html_node/Higher_002dlevel-Interface.html">
 * the usual routines</a>.  The return value must be subsequently cleaned up
 * using <tt>gsl_spline_free()</tt>.
 *
 * @param Re_x Local Reynolds number \f$\mbox{Re}_x\f$.
 *
 * @return On success, a <tt>gsl_spline *</tt> for use with
 *         <tt>gsl_spline_eval()</tt>.  On failure \c NULL is returned.
 */
gsl_spline * suzerain_blasius_u(const double Re_x);

/**
 * Obtain a Blasius profile fit producing nondimensional \f$v /
 * u_\infty\f$ given \f$y = \eta / \sqrt{\mbox{Re}_x}\f$.
 *
 * @copydetails suzerain_blasius_v
 */
gsl_spline * suzerain_blasius_v(const double Re_x);

/**
 * Obtain a Blasius profile fit producing nondimensional \f$\frac{T -
 * T_w}{T_\infty - T_w}\f$ given \f$y = \eta / \sqrt{\mbox{Re}_x \mbox{Pr}
 * }\f$.  The data used is identical to \ref suzerain_blasius_u except that it
 * is rescaled to produce a thermal boundary layer given similarity assumptions
 * employing a constant Prandtl number.
 *
 * The returned <tt>gsl_spline*</tt> can be interrogated using <a
 * href="http://www.gnu.org/software/gsl/manual/html_node/Higher_002dlevel-Interface.html">
 * the usual routines</a>.  The return value must be subsequently cleaned up
 * using <tt>gsl_spline_free()</tt>.
 *
 * @param Re_x Local Reynolds number \f$\mbox{Re}_x\f$.
 * @param Pr   Constant Prandtl number \f$\mbox{Pr}\f$.
 *
 * @return On success, a <tt>gsl_spline *</tt> for use with
 *         <tt>gsl_spline_eval()</tt>.  On failure \c NULL is returned.
 */
gsl_spline * suzerain_blasius_T(const double Re_x, const double Pr);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_BLASIUS_H */
