/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2011 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
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
 *
 * pre_gsl.h: functionality not present in public GSL releases
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_PRE_GSL_H
#define __SUZERAIN_PRE_GSL_H

#include <gsl/gsl_integration.h>

// gsl_integration_glfixed in GSL 1.14 but
// gsl_integration_glfixed_point is 1.14+ so build it atop 1.14's public API
// Source code lifted from GSL which is cool since it's my copyright :)
// FIXME: Remove this logic once GSL 1.15 becomes widespread
#if    (!defined GSL_MAJOR_VERSION                     ) \
    || (GSL_MAJOR_VERSION < 2 && GSL_MINOR_VERSION < 15)
int
gsl_integration_glfixed_point (
        double a,
        double b,
        size_t i,
        double *xi,
        double *wi,
        const gsl_integration_glfixed_table * t);
#endif

#endif // __SUZERAIN_PRE_GSL_H
