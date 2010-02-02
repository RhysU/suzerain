/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
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
 * underling_fft.h: Convenience wrappers around FFTW-like planning
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_UNDERLING_FFT_H
#define __SUZERAIN_UNDERLING_FFT_H

#include <suzerain/underling.h>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif
__BEGIN_DECLS
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** @file
 * Provides FFTW-like planning routines atop Underling's pencil decomposition
 * information via the underling_extents struct.
 */

/**
 * A type encapsulating FFTW-like planning information.
 */
typedef struct underling_fftplan_s *underling_fftplan;

/**
 * Create an FFT plan to perform a complex-to-complex transpose on the given
 * data when long in the <tt>i</tt>th direction.  The problem must have had
 * <tt>howmany</tt> specified as a multiple of two at creation time, and it
 * will be treated as <tt>howmany/2</tt> complex fields.  Note that the
 * transform is not normalized.
 *
 * @param problem Problem to use for layout and stride information.
 * @param data Pointer to the start of the memory allocated to execute
 *             the problem.
 * @param i    Direction across which to perform the FFT, which is assumed
 *             to be long and therefore stride one whenever the returned
 *             plan is executed.
 * @param fftw_sign Either FFTW_FORWARD or FFTW_BACKWARD.
 * @param fftw_rigor_flags One of FFTW's rigor planning flags, e.g.
 *                         FFTW_MEASURE.  Note that \c data is overwritten
 *                         during the planning process for any value other
 *                         than FFTW_ESTIMATE.
 *
 * @return On success, return a valid \c underling_fftplan.  On failure, calls
 *         suzerain_error and returns NULL.
 * @see The method underling_fftplan_destroy for how to destroy an instance.
 */
underling_fftplan
underling_fftplan_create_c2c(
        const underling_problem problem,
        underling_real * data,
        int i,
        unsigned fftw_sign,
        unsigned fftw_rigor_flags);

// FIXME Implement
// FIXME Document
underling_fftplan
underling_fftplan_create_r2c(
        const underling_problem problem,
        underling_real * data,
        int i,
        unsigned fftw_rigor_flags);

// FIXME Implement
// FIXME Document
underling_fftplan
underling_fftplan_create_c2r(
        const underling_problem problem,
        underling_real * data,
        int i,
        unsigned fftw_rigor_flags);

/**
 * Perform a previously planned FFT.  Appropriate calls to the underlying FFT
 * implementation will occur.
 *
 * @param fftplan Plan to be executed.
 *
 * @return SUZERAIN_SUCCESS (zero) on success and non-zero on failure.
 */
int
underling_fftplan_execute(
        const underling_fftplan fftplan);
/**
 * Destroy all resources associated with the given fftplan.
 *
 * @param fftplan Plan to be destroyed.
 */
void
underling_fftplan_destroy(
        underling_fftplan fftplan);

/**
 * Dump an instance's internals in a debugging-friendly format.
 *
 * @param fftplan Plan to be dumped.
 * @param output_file Desired output handle,
 *                    which may be \c stdout or \c stderr.
 */
void
underling_fprint_fftplan(
        const underling_fftplan fftplan,
        FILE *output_file);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
__END_DECLS
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif // __SUZERAIN_UNDERLING_FFT_H
