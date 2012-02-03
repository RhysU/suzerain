/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2011, 2012 The PECOS Development Team
 *
 * This file is part of Suzerain.
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This file is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This file is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this file.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 * pre_gsl.h: functionality not present in public GSL releases
 * $Id$
 */

#ifndef __SUZERAIN_PRE_GSL_H
#define __SUZERAIN_PRE_GSL_H

#ifdef __cplusplus
extern "C" {
#endif

// gsl_integration_glfixed in GSL 1.14 but
// gsl_integration_glfixed_point is 1.14+ so build it atop 1.14's public API
// FIXME: Remove this logic once GSL 1.15 becomes widespread
#if    (!defined GSL_MAJOR_VERSION                     ) \
    || (GSL_MAJOR_VERSION < 2 && GSL_MINOR_VERSION < 15)
#include <gsl/gsl_integration.h>

int
gsl_integration_glfixed_point (
        double a,
        double b,
        size_t i,
        double *xi,
        double *wi,
        const gsl_integration_glfixed_table * t);
#endif

// gsl_bspline_knots_greville is 1.15+ so build it atop 1.15's public API
// FIXME: Remove this logic once GSL 1.16 becomes widespread
#if    (!defined GSL_MAJOR_VERSION                     ) \
    || (GSL_MAJOR_VERSION < 2 && GSL_MINOR_VERSION < 16)
#include <gsl/gsl_bspline.h>

int
gsl_bspline_knots_greville(const gsl_vector *abscissae,
                           gsl_bspline_workspace *w,
                           double *abserr);
#endif

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // __SUZERAIN_PRE_GSL_H
