/***************************************************************************
//-----------------------------------------------------------------------bl-
// underling 0.3.1: an FFTW MPI-based library for 3D pencil decompositions
// http://red.ices.utexas.edu/projects/underling
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
//
// This file is part of underling.
//
// underling is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// underling is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with underling.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------el-
// $Id$
***************************************************************************/

#ifndef UNDERLING_FFTW_H
#define UNDERLING_FFTW_H

#include <underling/underling.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @file
 * Provides FFTW-like planning routines atop Underling's pencil decomposition
 * information via the underling_extents struct.
 */

/** @name Creation, execution, and destruction of underling_fftw_plans
 * @{
 */

/**
 * A type encapsulating FFTW-like planning information.
 */
typedef struct underling_fftw_plan_s *underling_fftw_plan;

/**
 * Flag indicating the FFT transform output and input for the "long in n2"
 * direction is packed contiguously in memory.  Transform input must likewise
 * be packed.  This flag may only be used for out-of-place transforms and
 * likely will incur additional memory access cost with each transform.  The
 * flag is provided for situations where other compute kernels benefit greatly
 * from having packed, contiguous storage and/or for compatibility with third
 * party libraries.
 *
 * When combined with UNDERLING_TRANSPOSED_LONG_N2, the FFT transform output
 * and input <em>is not</em> stored row-major <tt>n2 x (n0/pB x n1/pA)</tt> but
 * rather is stored packed contiguously as <tt>(n0/pB x n1/pA) x n2</tt>.
 * In-memory reshuffling is minimized in this circumstance.
 *
 * @see The documentation for underling_fftw_plan_create_c2c_forward,
 * underling_fftw_plan_create_c2c_backward,
 * underling_fftw_plan_create_r2c_forward, or
 * underling_fftw_plan_create_c2r_backward for more details on creating plans.
 */
#define UNDERLING_FFTW_PACKED_LONG_N2 (1U << 6)

/**
 * Flag indicating the FFT transform output and input for the "long in n0"
 * direction is packed contiguously in memory.  Transform input must likewise
 * be packed.  This flag may only be used for out-of-place transforms and
 * likely will incur additional memory access cost with each transform.  The
 * flag is provided for situations where other compute kernels benefit greatly
 * from having packed, contiguous storage and/or for compatibility with third
 * party libraries.
 *
 * When combined with UNDERLING_TRANSPOSED_LONG_N0, the FFT transform output
 * and input <em>is not</em> stored row-major <tt>n0 x (n1/pB x * n2/pA)</tt>
 * but rather is stored packed contiguously as <tt>(n1/pB x n2/pA) x n0</tt>.
 * In-memory reshuffling is minimized in this circumstance.
 *
 * @see The documentation for underling_fftw_plan_create_c2c_forward,
 * underling_fftw_plan_create_c2c_backward,
 * underling_fftw_plan_create_r2c_forward, or
 * underling_fftw_plan_create_c2r_backward for more details on creating plans.
 */
#define UNDERLING_FFTW_PACKED_LONG_N0 (1U << 7)

/** Convenience flag indicating packed transform output whenever possible. */
#define UNDERLING_FFTW_PACKED_ALL                                       \
        (UNDERLING_FFTW_PACKED_LONG_N2 | UNDERLING_FFTW_PACKED_LONG_N0)

/**
 * Flag indicating the FFT transform output for no direction is necessarily
 * packed contiguously in memory.  Under some circumstances, transform output
 * and input may be packed inadvertently (e.g. complex-to-complex transforms on
 * regular grids).
 *
 * @see UNDERLING_FFTW_PACKED_LONG_N2, UNDERLING_FFTW_PACKED_LONG_N0,
 *      and UNDERLING_FFTW_PACKED_ALL for alternatives.
 */
#define UNDERLING_FFTW_PACKED_NONE (1U << 8)

/**
 * Create a plan to perform a forward complex-to-complex FFT on the given data
 * when long in the <tt>long_ni</tt>th direction.  The problem must have had
 * <tt>howmany</tt> specified as a multiple of two at creation time, and it
 * will be treated as <tt>howmany/2</tt> complex fields.  Note that the
 * transform is not normalized.
 *
 * Out-of-place plans are created by specifying input and output buffers such
 * that <tt>in != out</tt>.  Executing an out-of-place plan will always
 * destroy the contents of the input buffer \c in.  In-place plans can be
 * created by specifying <tt>in == out</tt>.  In-place plans always use less
 * memory but may run more slowly than out-of-place plans.
 *
 * @param problem Problem to use for layout and stride information.
 * @param long_ni Direction across which to perform the FFT, which is assumed
 *                to be long whenever the returned plan is executed.
 * @param in  Input buffer containing the source data to be transformed.
 * @param out Output buffer to contain the data after transformation.
 * @param fftw_rigor_flags One of FFTW's rigor planning flags, e.g.
 *        FFTW_ESTIMATE.  Specifying zero is equivalent to providing
 *        FFTW_MEASURE.  Note that the buffers \c in and \c out are overwritten
 *        during the planning process for any value other than FFTW_ESTIMATE.
 * @param packed_flags One of UNDERLING_FFTW_PACKED_LONG_N2,
 *        UNDERLING_FFTW_PACKED_LONG_N0, UNDERLING_FFTW_PACKED_ALL, or
 *        UNDERLING_FFTW_PACKED_NONE.  Specifying zero is equivalent to
 *        providing UNDERLING_FFTW_PACKED_NONE.  Note only
 *        UNDERLING_FFTW_PACKED_NONE is valid when creating an in-place plan.
 *
 * @return On success, return a valid \c underling_fftw_plan.  On failure,
 *         calls underling_error and returns NULL.
 * @see The method underling_fftw_plan_destroy for how to destroy an instance.
 * @see The method underling_fftw_plan_create_inverse for how to create
 *      the corresponding inverse FFT.  It is <em>incorrect</em> to use
 *      any other way to invert the return FFT.
 */
underling_fftw_plan
underling_fftw_plan_create_c2c_forward(
        const underling_problem problem,
        int long_ni,
        underling_real * in,
        underling_real * out,
        unsigned fftw_rigor_flags,
        unsigned packed_flags) UNDERLING_API;

/**
 * Create a plan to perform a backward complex-to-complex FFT on the given data
 * when long in the <tt>long_ni</tt>th direction.
 * @copydetails underling_fftw_plan_create_c2c_forward
 */
underling_fftw_plan
underling_fftw_plan_create_c2c_backward(
        const underling_problem problem,
        int long_ni,
        underling_real * in,
        underling_real * out,
        unsigned fftw_rigor_flags,
        unsigned packed_flags) UNDERLING_API;

/**
 * Create a plan to perform a forward real-to-complex FFT on the given data
 * when long in the <tt>long_ni</tt>th direction.
 * @copydetails underling_fftw_plan_create_c2c_forward
 */
underling_fftw_plan
underling_fftw_plan_create_r2c_forward(
        const underling_problem problem,
        int long_ni,
        underling_real * in,
        underling_real * out,
        unsigned fftw_rigor_flags,
        unsigned packed_flags) UNDERLING_API;

/**
 * Create a plan to perform a backward complex-to-real FFT on the given data
 * when long in the <tt>long_ni</tt>th direction.
 * @copydetails underling_fftw_plan_create_c2c_forward
 */
underling_fftw_plan
underling_fftw_plan_create_c2r_backward(
        const underling_problem problem,
        int long_ni,
        underling_real * in,
        underling_real * out,
        unsigned fftw_rigor_flags,
        unsigned packed_flags) UNDERLING_API;

/**
 * Create a plan to invert another underling_fftw_plan.  Inverse plans
 * appropriately account for all input ordering issues stemming from use of
 * flags like UNDERLING_TRANSPOSED_LONG_N2, UNDERLING_TRANSPOSED_LONG_N0,
 * UNDERLING_FFTW_PACKED_LONG_N2, or UNDERLING_FFTW_PACKED_LONG_N0.  Note that
 * the inverse transform is not normalized.  Plan pairs created using this
 * method will have compatible input and output underling_fftw_extents
 * information.
 *
 * Out-of-place plans are created by specifying input and output buffers such
 * that <tt>in != out</tt>.  Executing an out-of-place plan will always
 * destroy the contents of the input buffer \c in.  In-place plans can be
 * created by specifying <tt>in == out</tt>.  In-place plans always use less
 * memory but will often run more slowly than out-of-place plans.
 *
 * @param plan_to_invert Prior plan to use for layout and stride information.
 * @param in  Input buffer containing the source data to be transformed.
 * @param out Output buffer to contain the data after transformation.
 * @param fftw_rigor_flags One of FFTW's rigor planning flags, e.g.
 *        FFTW_ESTIMATE.  Specifying zero is equivalent to providing
 *        FFTW_MEASURE.  Note that the buffers \c in and \c out are overwritten
 *        during the planning process for any value other than FFTW_ESTIMATE.
 *
 * @return On success, return a valid \c underling_fftw_plan which inverts
 *         plan_to_invert up to normalization.  On failure, calls
 *         underling_error and returns NULL.
 * @see The method underling_fftw_plan_destroy for how to destroy an instance.
 */
underling_fftw_plan
underling_fftw_plan_create_inverse(
        const underling_fftw_plan plan_to_invert,
        underling_real * in,
        underling_real * out,
        unsigned fftw_rigor_flags) UNDERLING_API;

/**
 * Perform a previously planned FFT.  Appropriate calls to the underlying FFT
 * implementation will occur.  The input and output buffers <em>must</em> be
 * aligned identically to the input and output buffers provided during
 * planning.
 *
 * @param plan Plan to be executed.
 * @param in  Input buffer on which to execute the plan.  For out-of-place
 *            transforms, this buffer's contents will be destroyed.
 * @param out Output buffer on which to execute the plan.  For in-place
 *            transforms, one must specify <tt>out == in</tt>.
 *
 * @return UNDERLING_SUCCESS (zero) on success and non-zero on failure.
 */
int
underling_fftw_plan_execute(
        const underling_fftw_plan plan,
        underling_real * in,
        underling_real * out) UNDERLING_API;

/**
 * Destroy all resources associated with the given plan.
 *
 * @param plan Plan to be destroyed.
 */
void
underling_fftw_plan_destroy(
        underling_fftw_plan plan) UNDERLING_API;

/**
 * Dump an instance's internals in a debugging-friendly format.
 *
 * @param plan Plan to be dumped.
 * @param output_file Desired output handle,
 *                    which may be \c stdout or \c stderr.
 */
void
underling_fftw_fprint_plan(
        const underling_fftw_plan plan,
        FILE *output_file) UNDERLING_API;

/**@}*/

/** @name Obtaining storage details for an underling_fftw_plan
 * @{
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
typedef struct underling_fftw_extents {
    /**
     * The inclusive global data starting offset in directions
     * <tt>n{0,1,2}</tt>.  Indices <tt>3</tt> and <tt>4</tt> give information
     * on the interleaved data fields and are always equal to zero.
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
} underling_fftw_extents;

/** A static instance used to communicate wholly invalid extents */
extern const underling_fftw_extents UNDERLING_FFTW_EXTENTS_INVALID;

/**
 * Compare two underling_fftw_extents instances using lexicographic ordering.
 *
 * @param e1 First instance to compare.
 * @param e2 Second instance to compare.
 *
 * @return Returns an integer less than, equal to, or greater than zero if
 *         <tt>*e1</tt> is found, respectively, to be less than, to match,
 *         or be greater than <tt>*e2</tt>.
 */
int
underling_fftw_extents_cmp(
        const underling_fftw_extents * const e1,
        const underling_fftw_extents * const e2) UNDERLING_API;

/**
 * Dump an instance's internals in a debugging-friendly format.
 *
 * @param extents Extents to dump.
 * @param output_file Desired output handle,
 *                    which may be \c stdout or \c stderr.
 */
void
underling_fftw_fprint_extents(
        const underling_fftw_extents *extents,
        FILE *output_file) UNDERLING_API;

/**
 * Obtain local size, stride, and storage information for the input data
 * to a given plan.  This information should be used to "load" the data
 * prior to executing the plan with underling_fftw_plan_execute.
 *
 * @param plan Plan for which to retrieve information.
 *
 * @return a valid underling_fftw_extents structure on success.  On failure,
 *         calls underling_error and returns UNDERLING_EXTENTS_INVALID.
 * @see The method underling_fftw_local_input for a way to obtain only a subset
 *      of this information, or for a more Fortran-ready interface.
 */
underling_fftw_extents
underling_fftw_local_extents_input(
        const underling_fftw_plan plan) UNDERLING_API;

/**
 * Obtain local size, stride, and storage information for the input data
 * to a given plan.  This information should be used to process the data
 * after executing the plan with underling_fftw_plan_execute.
 *
 * @param plan Plan for which to retrieve information.
 *
 * @return a valid underling_fftw_extents structure on success.  On failure,
 *         calls underling_error and returns UNDERLING_EXTENTS_INVALID.
 * @see The method underling_fftw_local_output for a way to obtain only a subset
 *      of this information, or for a more Fortran-ready interface.
 */
underling_fftw_extents
underling_fftw_local_extents_output(
        const underling_fftw_plan plan) UNDERLING_API;

/**
 * Obtain the processor-local sizes, storage details, and global starting
 * offsets for the given plan's input data when long in the direction for
 * which the plan was created.  This is identical to the data
 * obtainable via underling_fftw_local_extents_input but is provided in a more
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
 * @return UNDERLING_SUCCESS (zero) on success and non-zero on failure.
 *
 * @see The method underling_fftw_local_extents_input for a more C-friendly and
 *      const-correct capable way to obtain all of this information.
 */
int
underling_fftw_local_input(
        const underling_fftw_plan plan,
        int *start,
        int *size,
        int *stride,
        int *order) UNDERLING_API;

/**
 * Obtain the processor-local sizes, storage details, and global starting
 * offsets for the given plan's output data when long in the direction for
 * which the plan was created.  This is identical to the data
 * obtainable via underling_fftw_local_extents_output but is provided in a more
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
 * @return UNDERLING_SUCCESS (zero) on success and non-zero on failure.
 *
 * @see The method underling_fftw_local_extents_output for a more C-friendly and
 *      const-correct capable way to obtain all of this information.
 */
int
underling_fftw_local_output(
        const underling_fftw_plan plan,
        int *start,
        int *size,
        int *stride,
        int *order) UNDERLING_API;

/**@}*/

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // UNDERLING_FFTW_H
