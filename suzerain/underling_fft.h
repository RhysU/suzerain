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

#ifdef __cplusplus
extern "C" {
#endif

/** @file
 * Provides FFTW-like planning routines atop Underling's pencil decomposition
 * information via the underling_extents struct.
 */

/**
 * A transparent type storing the local sizes, strides, and storage when the
 * data is long in a particular direction \c n0, \c n1, or \c n2.
 *
 * The global data stored locally in direction <tt>i</tt> in {0,1,2,3,4} is
 * <tt>[start[i],start[i] + size[i])</tt> where the lower index is inclusive
 * and the upper index exclusive.  If both indices are the same then no data is
 * stored locally in the given direction.
 *
 * Indices {0,1,2} correspond to information about directions n{0,1,2},
 * respectively.  Index <tt>3</tt> gives information regarding the interleaved
 * data fields whose number is determined by the <tt>howmany</tt> parameter to
 * underling_problem_create. Index <tt>4</tt> gives information about the
 * scalar components that comprise index 3.  It always has size two for
 * a complex-valued field and size one for a real-valued field.
 */
typedef struct underling_fft_extents {
    /**
     * The inclusive global data starting offset in directions
     * <tt>n{0,1,2}</tt>.  Indices <tt>3</tt> and <tt>4</tt> give information
     * on the interleaved data fields and is always equal to zero.
     **/
    int start[5];

    /**
     * The amount of global data stored locally in directions
     * <tt>n{0,1,2}</tt>.  Index <tt>3</tt> gives information on the
     * interleaved data fields and is always equal to the <tt>howmany/2</tt>
     * parameter provided to underling_problem_create.  Index <tt>4</tt> is
     * <tt>2</tt> for a complex-valued field and <tt>1</tt> for a real-valued
     * field.
     */
    int size[5];

    /**
     * The stride between adjacent elements in directions <tt>n{0,1,2}</tt>.
     * Indices <tt>3</tt> and <tt>4</tt> give information on the interleaved
     * data fields.
     */
    int stride[5];

    /**
     * The storage order from the fastest index to the slowest.  That is,
     * <tt>order[0]</tt> gives the index in {0,1,2,3,4} of the fastest
     * direction, <tt>order[1]</tt> gives the index of the next fastest
     * direction, etc.  Useful in generic algorithms which are independent of
     * which direction is long but in which you need to walk memory optimally.
     */
    int order[5];
} underling_fft_extents;

/** A static instance used to communicate wholly invalid extents */
extern const underling_fft_extents UNDERLING_FFT_EXTENTS_INVALID;

/**
 * Compare two underling_fft_extents instances using lexicographic ordering.
 *
 * @param e1 First instance to compare.
 * @param e2 Second instance to compare.
 *
 * @return Returns an integer less than, equal to, or greater than zero if
 *         <tt>*e1</tt> is found, respectively, to be less than, to match,
 *         or be greater than <tt>*e2</tt>.
 */
int
underling_fft_extents_cmp(const underling_fft_extents * const e1,
                          const underling_fft_extents * const e2);

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
 * @see The method underling_fft_plan_create_inverse for how to create
 *      the corresponding inverse complex-to-complex backward FFT.  It is
 *      <em>incorrect</em> to use underling_fft_create_c2c_backward for
 *      that purpose.
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
 * @see The method underling_fft_plan_create_inverse for how to create
 *      the corresponding inverse complex-to-complex forward FFT.  It is
 *      <em>incorrect</em> to use underling_fft_create_c2c_forward for
 *      that purpose.
 */
underling_fft_plan
underling_fft_plan_create_c2c_backward(
        const underling_problem problem,
        int long_ni,
        underling_real * data,
        unsigned fftw_rigor_flags);

/**
 * Create a plan to perform a forward real-to-complex FFT on the given data
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
 * @see The method underling_fft_plan_create_inverse for how to create
 *      the corresponding inverse complex-to-real backward FFT.  It is
 *      <em>incorrect</em> to use underling_fft_create_c2r_backward for
 *      that purpose.
 */
underling_fft_plan
underling_fft_plan_create_r2c_forward(
        const underling_problem problem,
        int long_ni,
        underling_real * data,
        unsigned fftw_rigor_flags);

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
 * @see The method underling_fft_plan_create_inverse for how to create
 *      the corresponding inverse real-to-complex forward FFT.  It is
 *      <em>incorrect</em> to use underling_fft_create_r2c_forward for
 *      that purpose.
 */
underling_fft_plan
underling_fft_plan_create_c2r_backward(
        const underling_problem problem,
        int long_ni,
        underling_real * data,
        unsigned fftw_rigor_flags);

/**
 * Create a plan to invert another underling_fft_plan.  Inverse plans
 * appropriately account for all input ordering issues stemming from use of
 * flags like UNDERLING_TRANSPOSED_LONG_N2 or UNDERLING_TRANSPOSED_LONG_N0.
 * Note that the inverse transform is not normalized.  Plan pairs created using
 * this method will have compatible input and output underling_fft_extents
 * information.
 *
 * @param plan_to_invert Prior plan to use for layout and stride information.
 * @param data Pointer to the start of the memory allocated to execute
 *             the problem.
 * @param fftw_rigor_flags One of FFTW's rigor planning flags, e.g.
 *                         FFTW_MEASURE.  Note that \c data is overwritten
 *                         during the planning process for any value other
 *                         than FFTW_ESTIMATE.
 *
 * @return On success, return a valid \c underling_fft_plan which inverts
 *         plan_to_invert up to normalization.  On failure, calls
 *         suzerain_error and returns NULL.
 * @see The method underling_fft_plan_destroy for how to destroy an instance.
 */
underling_fft_plan
underling_fft_plan_create_inverse(
        const underling_fft_plan plan_to_invert,
        underling_real * data,
        unsigned fftw_rigor_flags);

/**
 * Obtain local size, stride, and storage information for the input data
 * to a given plan.  This information should be used to "load" the data
 * prior to executing the plan with underling_fft_plan_execute.
 *
 * @param plan Plan for which to retrieve information.
 *
 * @return a valid underling_fft_extents structure on success.  On failure,
 *         calls suzerain_error and returns UNDERLING_EXTENTS_INVALID.
 * @see The method underling_fft_local_input for a way to obtain only a subset
 *      of this information, or for a more Fortran-ready interface.
 */
underling_fft_extents
underling_fft_local_extents_input(
        const underling_fft_plan plan);

/**
 * Obtain local size, stride, and storage information for the input data
 * to a given plan.  This information should be used to process the data
 * after executing the plan with underling_fft_plan_execute.
 *
 * @param plan Plan for which to retrieve information.
 *
 * @return a valid underling_fft_extents structure on success.  On failure,
 *         calls suzerain_error and returns UNDERLING_EXTENTS_INVALID.
 * @see The method underling_fft_local_output for a way to obtain only a subset
 *      of this information, or for a more Fortran-ready interface.
 */
underling_fft_extents
underling_fft_local_extents_output(
        const underling_fft_plan plan);

/**
 * Obtain the processor-local sizes, storage details, and global starting
 * offsets for the given plan's input data when long in the direction for
 * which the plan was created.  This is identical to the data
 * obtainable via underling_fft_local_extents_input but is provided in a more
 * Fortran-ready interface.  All strides and sizes are given in units of
 * underling_real.
 *
 * @param[in]     plan Plan for which to retrieve information.
 * @param[in,out] start If non-NULL on entry, contains
 *                underling_extents.start on successful return.
 * @param[in,out] size If non-NULL on entry, contains
 *                underling_extents.size on successful return.
 * @param[in,out] stride If non-NULL on entry, contains
 *                underling_extents.stride on successful return.
 * @param[in,out] order If non-NULL on entry, contains
 *                underling_extents.order on successful return.
 *
 * @return SUZERAIN_SUCCESS (zero) on success and non-zero on failure.
 *
 * @see The method underling_fft_local_extents_input for a more C-friendly and
 *      const-correct capable way to obtain all of this information.
 */
int
underling_fft_local_input(
        const underling_fft_plan plan,
        int *start,
        int *size,
        int *stride,
        int *order);

/**
 * Obtain the processor-local sizes, storage details, and global starting
 * offsets for the given plan's output data when long in the direction for
 * which the plan was created.  This is identical to the data
 * obtainable via underling_fft_local_extents_output but is provided in a more
 * Fortran-ready interface.  All strides and sizes are given in units of
 * underling_real.
 *
 * @param[in]     plan Plan for which to retrieve information.
 * @param[in,out] start If non-NULL on entry, contains
 *                underling_extents.start on successful return.
 * @param[in,out] size If non-NULL on entry, contains
 *                underling_extents.size on successful return.
 * @param[in,out] stride If non-NULL on entry, contains
 *                underling_extents.stride on successful return.
 * @param[in,out] order If non-NULL on entry, contains
 *                underling_extents.order on successful return.
 *
 * @return SUZERAIN_SUCCESS (zero) on success and non-zero on failure.
 *
 * @see The method underling_fft_local_extents_output for a more C-friendly and
 *      const-correct capable way to obtain all of this information.
 */
int
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
 * @param extents Extents to dump.
 * @param output_file Desired output handle,
 *                    which may be \c stdout or \c stderr.
 */
void
underling_fft_fprint_extents(
        const underling_fft_extents *extents,
        FILE *output_file);

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

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // __SUZERAIN_UNDERLING_FFT_H
