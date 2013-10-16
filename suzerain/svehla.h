/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
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

#ifndef SUZERAIN_SVEHLA_H
#define SUZERAIN_SVEHLA_H

/** @file
 * Curve fits from
 * <a href="http://ntrs.nasa.gov/details.jsp?R=1020958">
 * NASA TR R-132</a> by Svehla Roger (1962).
 */

#include <gsl/gsl_spline.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Obtain a lookup table of air's viscosity in Pascal-seconds versus
 * temperature in degrees Kelvin according to Svehla 1962 table IV on page 117.
 *
 * The returned <tt>gsl_spline*</tt> can be interrogated using <a
 * href="http://www.gnu.org/software/gsl/manual/html_node/Higher_002dlevel-Interface.html"
 * the usual routines</a>.  The return value must be subsequently cleaned up
 * using <tt>gsl_spline_free()</tt>.
 *
 * @return On success, a <tt>gsl_spline *</tt> suitable for evaluation using,
 *         for example, <tt>gsl_spline_eval()</tt>.  On failure \c NULL is
 *         returned.
 */
gsl_spline * suzerain_svehla_air_mu_vs_T();

/** Molar mass of dry air according to Table IV in grams per mole */
#define SUZERAIN_SVEHLA_AIR_M (28.97)

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_SVEHLA_H */
