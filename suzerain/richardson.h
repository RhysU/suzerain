/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2014 Rhys Ulerich
 * Copyright (C) 2012-2014 The PECOS Development Team
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

#ifndef SUZERAIN_RICHARDSON_H
#define SUZERAIN_RICHARDSON_H

/** @file
 * Generalized Richardson extrapolation built atop the <a
 * href="http://www.gnu.org/software/gsl/">GNU Scientific Library</a> (GSL).
 */

#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector_double.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Perform Richardson extrapolation on \f$A_{i}(h)\f$ and \f$A_{i}(h/t)\f$ to
 * get \f$A_{i+1}(h)\f$.
 *
 * Given \f$A_{i}(h)\f$, an approximation of \f$A\f$ such that \f$A - A(h) =
 * a_{i} h^k_{i} + a_{i+1} h^k_{i+1} + \dots\f$ use Richardson extrapolation to
 * find an approximation to leading order \f$k_{i+1}\f$ from evaluations at
 * \f$A_{i}(h)\f$ and \f$A_{i}(h/t)\f$ for \f$t > 0\f$.  Specifically,
 * \f$A_{i+1}(h) = (t^{k_{i}} A_{i}(h/t) - A_{i}(h)) / (t^{k_{i}} - 1)\f$.
 *
 * @param[in,out] Ah On entry, the coarser approximation \f$A_{i}(h)\f$.
 *      On exit, the extrapolated approximation \f$A_{i+1}(h)\f$.
 * @param[in] Aht The finer approximation \f$A_{i}(h/t)\f$
 * @param[in] ki The leading term order in the Taylor expansion \f$k_{i}\f$
 * @param[in] t The refinement factor between the two approximations
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.  The parameter \c Ah will be
 *      in an undefined state after any error.
 */
int
suzerain_richardson_extrapolation_step(
        gsl_vector * Ah,
        const gsl_vector * Aht,
        const double t,
        const double ki);

/**
 * Perform Richardson extrapolation on a sequence of refined approximations.
 * The routine can perform multiple step extrapolation, e.g. using
 * \f$ A_{0}(h) \f$, \f$ A_{0}(h/2) \f$, and \f$ A_{0}(h/4) \f$ to compute
 * \f$ A_{2}(h) \f$.
 *
 * @see suzerain_richardson_extrapolation_step() for details on each
 *      extrapolation step and the terminology used.
 *
 * @param[in,out] A On entry, column \c i of \c A contains
 *              \f$ A_{0}(h/t^{i})\f$.  On exit, column \c 0 of \c A
 *              contains \f$ A_{A->size2}(h) \f$.  Other columns will be
 *              overwritten.
 * @param[in] t The refinement factor between the two approximations.
 * @param[in] k Leading error term orders.  Assumed to start at 1 if \c NULL
 *              provided.  Assumed to increment by 1 if no second element
 *              provided.  If the vector is shorter than <tt>A->size2</tt>, any
 *              unspecified higher index entries are assumed to grow by
 *              <tt>k(k->size-1) - k(k->size-2)</tt>.
 * @param[out] normtable If non-<tt>NULL</tt>, the routine outputs a table
 *              showing the \f$ \ell_2 \f$ error at each step in the
 *              extrapolation process calculated against the provided
 *              \c exact parameter.
 * @param[in] exact Exact solution used to calculate \c normtable.  If \c NULL,
 *              \c exact is treated as if it contained all zeros.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.  The parameters \c A and
 *      \c normtable will be in an undefined state after any error.
 */
int
suzerain_richardson_extrapolation(
        gsl_matrix * const A,
        const double t,
        const gsl_vector * k,
        gsl_matrix * normtable,
        const gsl_vector * const exact);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_RICHARDSON_H */
