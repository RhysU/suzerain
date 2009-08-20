/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2009 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * richardson.h: generalized Richardson extrapolation routines
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_RICHARDSON_H__
#define __SUZERAIN_RICHARDSON_H__

#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector_double.h>

/** @file
 * Provides generalized Richardson extrapolation built atop the <a
 * href="http://www.gnu.org/software/gsl/">GNU Scientific Library</a> (GSL).
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

/**
 * Perform Richardson extrapolation on \f$A_{i}(h)\f$ and \f$A_{i}(h/t)\f$ to
 * get \f$A_{i+1}(h)\f$.
 *
 * Given \f$A_{i}(h)\f$, an approximation of \f$A\f$ such that \f$A - A(h) =
 * a_{i}*h^k_{i} + a_{i+1}*h^k_{i+1} + \dots\f$ use Richardson extrapolation to
 * find an approximation to leading order \f$k_{i+1}\f$ from evaluations at
 * \f$A_{i}(h)\f$ and \f$A_{i}(h/t)\f$ for \f$t > 0\f$.  Specifically,
 * \f$A_{i+1}(h) = (t**k_{i}*A_{i}(h/t) - A_{i}(h)) / (t**k_{i} - 1)\f$.
 *
 * @param Aih the coarse approximation \f$A_{i}(h)\f$
 * @param Aiht the finer approximation \f$A_{i}(h/t)\f$
 * @param ki the leading error term in the Taylor expansion \f$k_{i}\f$
 * @param t the refinement factor between the two approximations
 * @param Aip1h the extrapolated result \f$A_{i+1}(h)\f$
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_richardson_step(
    const gsl_vector * Aih,
    const gsl_vector * Aiht,
    const int ki,
    const double t,
    gsl_vector * Aip1h);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
__END_DECLS
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif /* __SUZERAIN_RICHARDSON_H__ */
