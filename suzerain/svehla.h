/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2010, 2011, 2012 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *--------------------------------------------------------------------------
 * svehla.h: Provides information from NASA TR R-132 by Svehla Roger (1962)
 * $Id$
 */

#ifndef __SUZERAIN_SVEHLA_H
#define __SUZERAIN_SVEHLA_H

#include <gsl/gsl_spline.h>

/** @file
 * Provides information from
 * <a href="http://ntrs.nasa.gov/details.jsp?R=1020958">
 * NASA TR R-132</a> by Svehla Roger (1962).
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Obtain a lookup table of air's viscosity in Pascal-seconds versus
 * temperature in degrees Kelvin according to Svehla 1962 table IV on page 117.
 * The returned <tt>gsl_spline*</tt> can be interrogated using <a
 * href="http://www.gnu.org/software/gsl/manual/html_node/Higher_002dlevel-Interface.html"
 * the usual routines</a>.  The return value must be subsequently cleaned up
 * using <tt>gsl_spline_free()</tt>.
 *
 * @return On success, a <tt>gsl_spline *</tt> suitable for evaluation using,
 *         for example, <tt>gsl_spline_eval()</tt>.  The return value must be 
 *         subsequently cleaned up using <tt>gsl_spline_free()</tt>.  On failure
 *         \c NULL is returned.
 */
gsl_spline * suzerain_svehla_air_mu_vs_T();

/** Molar mass of dry air according to Table IV in grams per mole */
#define SUZERAIN_SVEHLA_AIR_M (28.97)

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __SUZERAIN_SVEHLA_H */
