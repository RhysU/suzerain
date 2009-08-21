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
 * timestepper.h: SMR91 timestepping implementation
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef PECOS_SUZERAIN_TIMESTEPPER_H
#define PECOS_SUZERAIN_TIMESTEPPER_H

/** @file
 * FIXME
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

typedef struct suzerain_lsrk_method {
    const char   * const name;
    const int            substeps;
    const double * const alpha;
    const double * const beta;
    const double * const gamma;
    const double * const zeta;
} suzerain_lsrk_method;

extern const suzerain_lsrk_method suzerain_lsrk_smr91;

int
suzerain_lsrk_substep(
        const suzerain_lsrk_method method,
        const int n,
        const int kl,
        const int ku,
        const double * const M,
        const int ldM,
        const int nD,
        const double * const xi,
        const double * const * const D,
        const int ldD,
        double delta_t,
        const int nrhs,
        double       * const a, const int inca, const int lda,
        const double * const b, const int incb, const int ldb,
        double       * const c, const int incc, const int ldc,
        const int substep);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
__END_DECLS
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif // PECOS_SUZERAIN_TIMESTEPPER_H
