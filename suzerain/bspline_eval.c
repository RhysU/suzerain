/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2009 The PECOS Development Team
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
 * bspline_eval.c: B-spline series evaluation routines
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include "config.h"

#include <alloca.h>
#include <stdlib.h>
#include <suzerain/bspline_eval.h>

void
suzerain_bspline_series_evaluate(
        const int n,
        const int k,
        const double * const t,
        const int d,
        const double x,
        const double * const C,
        double * const Q) {

    /* Find rmm such that x \in [t_{rmm + m}, t_{rmm + m + 1}] */
    const int rmm = suzerain_bspline_series_span_search(k, n+k-1, t, x);

    /* Compute the B-spline series on the span */
    suzerain_bspline_series_span_evaluate(k, t + rmm, d, x, C + rmm, Q);
}

void
suzerain_bspline_series_span_evaluate(
        const int k,
        const double * const t,
        const int d,
        const double x,
        const double * const A,
        double * const Q)
{
/* The computation uses E. T. Y. Lee's routine published in "A
 * Simplified B-Spline Computation Routine", Computing volume 29
 * pages 365--371 (1982).
 */
    const int m = k - 1; /* piecewise degree */

    if (m == 0) {
        Q[0] = A[0];
       return;
    }

    const size_t scratch_size = k * sizeof(x);
    double * const a = alloca(scratch_size);
    double * const b = alloca(scratch_size);

    for (int i = 1; i <= m; ++i) {
        a[i] = x - t[i];
        b[i] = t[m+i] - x;
    }

    const int fromright = b[1] > a[m];
    if (fromright) {
        for (int i = 0; i <= m; ++i) {
            Q[i] = A[i];
        }

        for (int j = 1; j <= m; ++j) {
            const int mmj = m - j;
            for (int i = 0; i <= mmj; ++i) {
                const int ipj = i + j;
                const int ip1 = i + 1;
                Q[i] = (a[ipj]*Q[ip1]+b[ip1]*Q[i])/(a[ipj]+b[ip1]);
            }
        }

        for (int j = 1; j <= d; ++j) {
            for (int i = d; i >= j; --i) {
                Q[i] = (Q[i]-Q[i-1])/(b[i-j+1]/(m-j+1));
            }
        }
    } else { /* !fromright */
        for (int i = 0; i <= m; ++i) {
            Q[i] = A[m-i];
        }

        for (int j = 1; j <= m; ++j) {
            const int mmj = m - j;
            for (int i = 0; i <= mmj; ++i) {
                const int mmi     = m - i;
                const int mmjp1mi = mmi - j + 1;
                const int ip1     = i + 1;
                Q[i] = (b[mmjp1mi]*Q[ip1]+a[mmi]*Q[i])/(a[mmi]+b[mmjp1mi]);
            }
        }

        for (int j = 1; j <= d; ++j) {
            for (int i = d; i >= j; --i) {
                Q[i] = (Q[i-1]-Q[i])/(a[m-i+j]/(m-j+1));
            }
        }
    }

    return;
}

int
suzerain_bspline_series_span_search(
        const int k,
        const int nt,
        const double * t,
        const double x)
{
    const int m = k - 1; /* piecewise degree */

    /* Perform a binary search over interior knots */
    const double * left  = t + m;
    const double * right = t + nt - m;
    int separation = right - left;
    while (separation > 1) {
        const double * const mid = left + separation/2;
        if (x < *mid) {
            right = mid;
        } else {
            left = mid;
        }
        separation = right - left;
    }

    /* Do not land between repeated knots */
    while (*left == *right) {
        /* Push to lower knots to enforce continuity: left <= x < right */
        /* At rightmost interior interval allows left <= x <= right */
        --left;
        --right;
    }

    return left - m - t;
}


