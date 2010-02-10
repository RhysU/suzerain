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

// FIXME Document underling_fft_extents

typedef struct underling_fft_extents {
    int start[5];

    int size[5];

    int stride[5];

    int order[5];
} underling_fft_extents;

/** A static instance used to communicate wholly invalid extents */
extern const underling_fft_extents UNDERLING_FFT_EXTENTS_INVALID;

/**
 * A type encapsulating FFTW-like planning information.
 */
typedef struct underling_fft_plan_s *underling_fft_plan;

/**
 * Create a plan to perform a forward complex-to-complex FFT on the given data
 * when long in the <tt>long_ni</tt>th direction.  The problem must have had
 * <tt>howmany</tt> specified as a multiple of two at creation time, and it
 * will be treated as <tt>howmany/2</tt> complex fields.  Note that the
 * transform is not normalized.
 *
 * @param problem Problem to use for layout and stride information.
 * @param long_ni Direction across which to perform the FFT, which is assumed
 *                to be long whenever the returned plan is executed.
 * @param data Pointer to the start of the memory allocated to execute
 *             the problem.
 * @param fftw_rigor_flags One of FFTW's rigor planning flags, e.g.
 *                         FFTW_MEASURE.  Note that \c data is overwritten
 *                         during the planning process for any value other
 *                         than FFTW_ESTIMATE.
 *
 * @return On success, return a valid \c underling_fft_plan.  On failure, calls
 *         suzerain_error and returns NULL.
 * @see The method underling_fft_plan_destroy for how to destroy an instance.
 */
underling_fft_plan
underling_fft_plan_create_c2c_forward(
        const underling_problem problem,
        int long_ni,
        underling_real * data,
        unsigned fftw_rigor_flags);

/**
 * Create a plan to perform a backward complex-to-complex FFT on the given data
 * when long in the <tt>long_ni</tt>th direction.  The problem must have had
 * <tt>howmany</tt> specified as a multiple of two at creation time, and it
 * will be treated as <tt>howmany/2</tt> complex fields.  Note that the
 * transform is not normalized.
 *
 * @param problem Problem to use for layout and stride information.
 * @param long_ni Direction across which to perform the FFT, which is assumed
 *                to be long whenever the returned plan is executed.
 * @param data Pointer to the start of the memory allocated to execute
 *             the problem.
 * @param fftw_rigor_flags One of FFTW's rigor planning flags, e.g.
 *                         FFTW_MEASURE.  Note that \c data is overwritten
 *                         during the planning process for any value other
 *                         than FFTW_ESTIMATE.
 *
 * @return On success, return a valid \c underling_fft_plan.  On failure, calls
 *         suzerain_error and returns NULL.
 * @see The method underling_fft_plan_destroy for how to destroy an instance.
 */
underling_fft_plan
underling_fft_plan_create_c2c_backward(
        const underling_problem problem,
        int long_ni,
        underling_real * data,
        unsigned fftw_rigor_flags);

// FIXME Add expected input stride information

// FIXME Implement
// FIXME Document
underling_fft_plan
underling_fft_plan_create_r2c_forward(
        const underling_problem problem,
        int long_ni,
        underling_real * data,
        unsigned fftw_rigor_flags);

// FIXME Add relevant transform size information to documentation

/**
 * Create a plan to perform a backward complex-to-real FFT on the given data
 * when long in the <tt>long_ni</tt>th direction.  The problem must have had
 * <tt>howmany</tt> specified as a multiple of two at creation time, and it
 * will be treated as <tt>howmany/2</tt> complex fields.  Note that the
 * transform is not normalized.
 *
 * @param problem Problem to use for layout and stride information.
 * @param long_ni Direction across which to perform the FFT, which is assumed
 *                to be long whenever the returned plan is executed.
 * @param data Pointer to the start of the memory allocated to execute
 *             the problem.
 * @param fftw_rigor_flags One of FFTW's rigor planning flags, e.g.
 *                         FFTW_MEASURE.  Note that \c data is overwritten
 *                         during the planning process for any value other
 *                         than FFTW_ESTIMATE.
 *
 * @return On success, return a valid \c underling_fft_plan.  On failure, calls
 *         suzerain_error and returns NULL.
 * @see The method underling_fft_plan_destroy for how to destroy an instance.
 */
underling_fft_plan
underling_fft_plan_create_c2r_backward(
        const underling_problem problem,
        int long_ni,
        underling_real * data,
        unsigned fftw_rigor_flags);

underling_fft_extents
underling_fft_local_extents_input(
        const underling_fft_plan plan);

underling_fft_extents
underling_fft_local_extents_output(
        const underling_fft_plan plan);

void
underling_fft_local_input(
        const underling_fft_plan plan,
        int *start,
        int *size,
        int *stride,
        int *order);

void
underling_fft_local_output(
        const underling_fft_plan plan,
        int *start,
        int *size,
        int *stride,
        int *order);

/**
 * Perform a previously planned FFT.  Appropriate calls to the underlying FFT
 * implementation will occur.
 *
 * @param plan Plan to be executed.
 *
 * @return SUZERAIN_SUCCESS (zero) on success and non-zero on failure.
 */
int
underling_fft_plan_execute(
        const underling_fft_plan plan);
/**
 * Destroy all resources associated with the given plan.
 *
 * @param plan Plan to be destroyed.
 */
void
underling_fft_plan_destroy(
        underling_fft_plan plan);

/**
 * Dump an instance's internals in a debugging-friendly format.
 *
 * @param plan Plan to be dumped.
 * @param output_file Desired output handle,
 *                    which may be \c stdout or \c stderr.
 */
void
underling_fft_fprint_plan(
        const underling_fft_plan plan,
        FILE *output_file);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
__END_DECLS
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif // __SUZERAIN_UNDERLING_FFT_H
