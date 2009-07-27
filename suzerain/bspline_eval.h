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
 * bspline_eval.h: B-spline series evaluation routines
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef PECOS_SUZERAIN_BSPLINE_EVAL_H
#define PECOS_SUZERAIN_BSPLINE_EVAL_H

/** @file
 * Provides B-spline series evaluation, including derivatives.
 */

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
/** Marks beginning of public declarations using C linkage for C++ compiler */
# define __BEGIN_DECLS extern "C" {
/** Marks ending of public declarations using C linkage for C++ compiler */
# define __END_DECLS }
#else
/** Marks beginning of public declarations for C compiler */
# define __BEGIN_DECLS /* empty */
/** Marks ending of public declarations for C compiler */
# define __END_DECLS /* empty */
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
__BEGIN_DECLS
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 * Evaluate a B-spline series and its derivatives at a point.
 *
 * Computes \f$ Q_j = D^j \sum_{i=0}^{n-1} C_i B_{i,k}(x) \f$ for \f$
 * j\in\left\{0,\ldots,d\right\}\f$ where \f$m=k-1\f$ is the piecewise degree
 * of the spline.
 *
 * \warning As this is a computational kernel, no error checking is performed.
 *      The caller must ensure all preconditions are satisfied.
 *
 * @param[in] n Number of coefficients in B-spline series
 * @param[in] k Spline order
 * @param[in] t Nondecreasing knot vector containing
 *      \f$ \left( t_i \right)_{i=0}^{n+m} \f$
 * @param[in] d Maximum derivative to compute where \f$d\leq{}m\f$
 * @param[in] C B-coefficients \f$ \left( C_i \right)_{i=0}^{n-1}
 * @param[in] x Point \f$x\in\left[t_{m},t_{n}\right]\f$ at which to evaluate
 *      the B-spline series.
 * @param[out] Q Result vector containing derivatives 0 through d.
 */
void
suzerain_bspline_series_evaluate(
        const int n,
        const int k,
        const double * const t,
        const int d,
        const double x,
        const double * const C,
        double * const Q);

/**
 * Evaluate a B-spline series and its derivatives on a single knot span.
 *
 * Computes \f$ Q_j = D^j \sum_{i=0}^m A_i B_{i,k}(x) \f$ for \f$
 * j\in\left\{0,\ldots,d\right\}\f$ where \f$m=k-1\f$ is the piecewise degree
 * of the spline.
 *
 * The computation uses E. T. Y. Lee's routine published in "A Simplified
 * B-Spline Computation Routine", \e Computing volume 29 pages 365&ndash;371
 * (1982).  Please note that this algorithm has stability problems compared
 * with the de Boor algorithm as discussed in Lee's "Comments on some
 * B-spline algorithms", Computing volume 36 pages 229&ndash;238 (1986).
 *
 * \warning As this is a computational kernel, no error checking is performed.
 *      The caller must ensure all preconditions are satisfied.
 *
 * @param[in] k Spline order
 * @param[in] t Knot vector containing \f$\left\{t_0,\ldots,t_{2m+1}\right\}\f$
 *      where \f$ t_m < t_{m+1} \f$.
 * @param[in] d Maximum derivative to compute where \f$d\leq{}m\f$
 * @param[in] A B-coefficients \f$\left\{A_0,\ldots,A_m\right\}\f$
 * @param[in] x Point \f$x\in\left[t_m,t_{m+1}\right]\f$ at which to evaluate
 *      the B-spline series.
 * @param[out] Q Result vector containing derivatives 0 through d.
 */
void
suzerain_bspline_series_span_evaluate(
        const int k,
        const double * const t,
        const int d,
        const double x,
        const double * const A,
        double * const Q);

/**
 * Find the offset of knot \f$ t_0 \f$ in vector \c t such that \f$ x\in\left[t_m,
 * t_{m+1}\right] \f$ for \f$m=k-1\f$.
 *
 * \warning As this is a computational kernel, no error checking is performed.
 *
 * @param[in] k Spline order
 * @param[in] nt Number of knots present in \c t.  This count must include the
 *      \f$ m \f$ repeated knots at each end of the knot vector.
 * @param[in] t Nondecreasing vector of knot locations
 * @param[in] x Value to find.
 *
 * @return The offset of knot \f$ t_0\f$ in vector \c t.
 */
int
suzerain_bspline_series_span_search(
        const int k,
        const int nt,
        const double * t,
        const double x);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
__END_DECLS
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif // PECOS_SUZERAIN_BSPLINE_EVAL_H
