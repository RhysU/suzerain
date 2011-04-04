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

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.h>
#pragma hdrstop
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>

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
        const gsl_integration_glfixed_table * t)
{
    const double A = (b - a) / 2;  /* Length of [a,b] */
    const double B = (a + b) / 2;  /* Midpoint of [a,b] */

    if (i >= t->n) {
        GSL_ERROR ("i must be less than t->n", GSL_EINVAL);
    }

    /* See comments above gsl_integration_glfixed for struct's x, w layout. */
    /* Simply unpack that layout into a sorted set of points, weights. */
    if (GSL_IS_ODD(t->n)) {
        const int k = ((int) i) - ((int) t->n) / 2;
        const int s = k < 0 ? -1 : 1;

        *xi = B + s*A*t->x[s*k];
        *wi =       A*t->w[s*k];
    } else if (/* GSL_IS_EVEN(t->n) && */ i < t->n / 2) {
        i = (t->n / 2) - 1 - i;
        *xi = B - A*t->x[i];
        *wi =     A*t->w[i];
    } else /* GSL_IS_EVEN(t->n) && i >= n / 2 */ {
        i  -= t->n / 2;
        *xi = B + A*t->x[i];
        *wi =     A*t->w[i];
    }

    return GSL_SUCCESS;
}
#endif
