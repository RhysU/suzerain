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
 * underling.h: A parallel, three dimensional FFT library atop MPI
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_UNDERLING_H
#define __SUZERAIN_UNDERLING_H

#include <suzerain/mpi.h>

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

// TODO Group Doxygen documentation by class-like concepts

/** @file
 * Provides a parallel, pencil-based 3D Cartesian domain decomposition atop
 * FFTW MPI.  The decomposition splits the 3D domain \c n0 by \c n1 by \c n2
 * across the 2D Cartesian topology given by \c pA by \c pB.  Parallel
 * transposes and local reordering are used to make the domain "long" and
 * stride one in each direction.  Arbitrary-length interleaved data sets may be
 * transposed.  Methods which must be called collectively in the MPI sense are
 * marked in their descriptions.
 */

/**
 * The real-valued scalar type used throughout the library.
 * All stride and size information is given in terms of this type.
 */
typedef double underling_real;

/** A type encapsulating a reusable domain-to-processor mapping. */
typedef struct underling_grid_s *underling_grid;

/**
 * A type encapsulating all parallel decomposition information,
 * including required local storage and stride details.
 */
typedef struct underling_problem_s *underling_problem;

/**
 * A type encapsulating the FFTW MPI invocations necessary
 * to transition from being long in one direction to long in another.
 * Similar in nature to an \c fftw_plan.
 */
typedef struct underling_plan_s *underling_plan;

/**
 * A transparent type storing the local sizes, strides, and storage when the
 * data is long in a particular direction \c n0, \c n1, or \c n2.
 *
 * The global data stored locally in direction <tt>i</tt> in {0,1,2,3} is
 * <tt>[start[i],start[i] + size[i])</tt> where the lower index is inclusive
 * and the upper index exclusive.  If both indices are the same then no data is
 * stored locally in the given direction.
 *
 * Indices {0,1,2} correspond to information about directions n{0,1,2},
 * respectively.  Index <tt>3</tt> gives information regarding the interleaved
 * data fields whose number is determined by the <tt>howmany</tt> parameter to
 * underling_problem_create.
 */
typedef struct underling_extents {
    /**
     * The inclusive global data starting offset in directions
     * <tt>n{0,1,2}</tt>.  Index <tt>3</tt> gives information on the
     * interleaved data fields and is always equal to zero.
     **/
    int start[4];

    /**
     * The amount of global data stored locally in directions
     * <tt>n{0,1,2}</tt>.  Index <tt>3</tt> gives information on the
     * interleaved data fields and is always equal to the <tt>howmany</tt>
     * parameter provided to underling_problem_create.
     */
    int size[4];

    /**
     * The stride between adjacent elements in directions <tt>n{0,1,2}</tt>.
     * Index <tt>3</tt> gives information on the interleaved data fields.
     */
    int stride[4];

    /**
     * The storage order from the fastest index to the slowest.  That is,
     * <tt>order[0]</tt> gives the index in {0,1,2,3} of the fastest
     * direction, <tt>order[1]</tt> gives the index of the next fastest
     * direction, etc.  Useful in generic algorithms which are independent of
     * which direction is long but in which you need to walk memory optimally.
     */
    int order[4];

    /**
     * The total extent of all local information, excluding communication
     * buffers, in units of underling_real.  Computation over all data when
     * long in a direction, e.g. uniform scaling, should use this size.
     *
     * It is technically redundant to store this information separately since
     * <tt>extent == size[0] * size[1] * size[2] * size[3]</tt>. However, is
     * very convenient to have the result handy.
     */
    size_t extent;
} underling_extents;

/** A static instance used to communicate wholly invalid extents */
extern const underling_extents UNDERLING_EXTENTS_INVALID;

/**
 * Compare two underling_extents instances using lexicographic ordering.
 *
 * @param e1 First instance to compare.
 * @param e2 Second instance to compare.
 *
 * @return Returns an integer less than, equal to, or greater than zero if
 *         <tt>*e1</tt> is found, respectively, to be less than, to match,
 *         or be greater than <tt>*e2</tt>.
 */
int
underling_extents_cmp(const underling_extents * const e1,
                      const underling_extents * const e2);

/**
 * Flag indicating a transform from long in \c n2 to long in \c n1.
 *
 * @see The documentation for underling_grid_create for more
 *      details on the associated storage order.
 **/
#define UNDERLING_TRANSPOSE_LONG_N2_TO_LONG_N1 (1U << 0)

/**
 * Flag indicating a transform from long in \c n1 to long in \c n0.
 *
 * @see The documentation for underling_grid_create for more
 *      details on the associated storage order.
 */
#define UNDERLING_TRANSPOSE_LONG_N1_TO_LONG_N0 (1U << 1)

/**
 * Flag indicating a transform from long in \c n0 to long in \c n1.
 *
 * @see The documentation for underling_grid_create for more
 *      details on the associated storage order.
 */
#define UNDERLING_TRANSPOSE_LONG_N0_TO_LONG_N1 (1U << 2)

/**
 * Flag indicating a transform from long in \c n1 to long in \c n2.
 *
 * @see The documentation for underling_grid_create for more
 *      details on the associated storage order.
 */
#define UNDERLING_TRANSPOSE_LONG_N1_TO_LONG_N2 (1U << 3)

/** Convenience flag indicating all transform directions */
#define UNDERLING_TRANSPOSE_ALL \
         (   UNDERLING_TRANSPOSE_LONG_N2_TO_LONG_N1   \
           | UNDERLING_TRANSPOSE_LONG_N1_TO_LONG_N0   \
           | UNDERLING_TRANSPOSE_LONG_N0_TO_LONG_N1   \
           | UNDERLING_TRANSPOSE_LONG_N1_TO_LONG_N2 )

/**
 * Flag indicating that "long in n2" storage is stored as <tt>n2 x (n0/pB x
 * n1/pA)</tt> in row-major order.  This differs from the normal storage
 * documented at underling_grid_create.  Using this flag may speed up
 * underling_execute_long_n2_to_long_n1 and
 * underling_execute_long_n1_to_long_n2.  Both underling_local_extents and
 * underling_local will automatically return information reflecting this
 * storage choice.
 *
 * @see The documentation for underling_grid_create for more
 *      details on storage ordering.
 */
#define UNDERLING_TRANSPOSED_LONG_N2 (1U << 4)

/**
 * Flag indicating that "long in n0" storage is stored as <tt>n0 x (n1/pB x
 * n2/pA)</tt> in row-major order.  This differs from the normal storage
 * documented at underling_grid_create.  Using this flag may speed up
 * underling_execute_long_n1_to_long_n0 and
 * underling_execute_long_n0_to_long_n1.  Both underling_local_extents and
 * underling_local will automatically return information reflecting this
 * storage choice.
 *
 * @see The documentation for underling_grid_create for more
 *      details on storage ordering.
 */
#define UNDERLING_TRANSPOSED_LONG_N0 (1U << 5)

/**
 * Collectively create a reusable domain-to-processor mapping across the given
 * MPI communicator.
 *
 * The 3D domain \c n0 by \c n1 by \c n2 will be split across the 2D Cartesian
 * topology given by \c pA by \c pB as follows:
 *  - When "long in n2" storage is <tt>(n0/pB x n1/pA) x n2</tt>.
 *  - When "long in n1" storage is <tt>(n2/pA x n0/pB) x n1</tt>.
 *  - When "long in n0" storage is <tt>(n1/pB x n2/pA) x n0</tt>.
 *
 * All orders given are row-major; the rightmost index is fastest.  Expressions
 * like <tt>n{0,1,2}/p{A,B}</tt> indicate that the <tt>n{0,1,2}</tt> direction
 * is decomposed across a grid of size <tt>p{A,B}</tt>.  The decomposition is
 * balanced given communication overhead expectations.
 *
 * \note It is much, much easier to obtain the appropriate storage order and
 * stride information from underling_local_extents or underling_local rather
 * than to compute it yourself.
 *
 * Specifying zero for either or both of \c pA and \c pB results in an
 * automatic decomposition of the communicator into a 2D Cartesian grid using
 * <tt>MPI_Dims_create</tt>.
 *
 * @param comm MPI communicator indicating the processes to be used
 *             for the parallel domain decomposition.  The communicator
 *             is cloned.
 * @param n0 Global size of the domain in the \c n0 direction.
 * @param n1 Global size of the domain in the \c n1 direction.
 * @param n2 Global size of the domain in the \c n2 direction.
 * @param pA Processor grid size in the \c pA direction.  Providing
 *           zero causes automatic size selection.
 * @param pB Processor grid size in the \c pB direction.  Providing
 *           zero causes automatic size selection.
 *
 * @return A valid, non-NULL underling_grid instance on success.  On failure,
 *         suzerain_error is invoked and NULL is returned.
 *
 * @see The method underling_grid_destroy for how to destroy an instance.
 */
underling_grid
underling_grid_create(
        MPI_Comm comm,
        int n0,
        int n1,
        int n2,
        int pA,
        int pB);

/**
 * Obtain the size of the processor grid in the \c pA direction.
 *
 * @param grid Grid for which to retrieve information.
 *
 * @return On success, the runtime size of the \c pA direction.
 *         On failure, calls suzerain_error and returns null.
 * @see The method underling_grid_create for how and when data
 *      is split across the \c pA processor grid direction.
 */
int
underling_grid_pA_size(
        const underling_grid grid);

/**
 * Obtain the size of the processor grid in the \c pB direction.
 *
 * @param grid Grid for which to retrieve information.
 *
 * @return On success, the runtime size of the \c pB direction.
 *         On failure, calls suzerain_error and returns null.
 * @see The method underling_grid_create for how and when data
 *      is split across the \c pB processor grid direction.
 */
int
underling_grid_pB_size(
        const underling_grid grid);

/**
 * Destroy all resources associated with the given grid.
 *
 * @param grid Grid to be destroyed.
 */
void
underling_grid_destroy(
        underling_grid grid);

/**
 * Collectively create an instance encapsulating the parallel decomposition
 * details for a decomposition over \c grid containing \c howmany
 * scalar fields of type underling_real.
 *
 * @param grid Domain-to-processor mapping to use.  One grid may be
 *        used to create multiple underling_problem instances.
 * @param howmany Number of interleaved scalar fields of type underling_real to
 *        simultaneously transpose.
 * @param transposed_flags Either zero or some combination of
 *        UNDERLING_TRANSPOSED_LONG_N2 and UNDERLING_TRANSPOSED_LONG_N0.
 *
 * @return A valid, non-NULL underling_problem instance on success.  On failure,
 *         suzerain_error is invoked and NULL is returned.
 *
 * @see UNDERLING_TRANSPOSED_LONG_N2 and/or UNDERLING_TRANSPOSED_LONG_N0
 *      for information on how their usage effects storage ordering.
 * @see The method underling_problem_destroy for how to destroy an instance.
 */
underling_problem
underling_problem_create(
        underling_grid grid,
        int howmany,
        unsigned transposed_flags);

/**
 * Destroy all resources associated with the given problem.
 *
 * @param problem Problem to be destroyed.
 */
void
underling_problem_destroy(
        underling_problem problem);

/**
 * Obtain the processor-local sizes, storage details, and global starting
 * offsets when the data is long in the n<tt>i</tt> direction.
 * All strides and sizes are given in units of underling_real.
 *
 * @param problem Problem for which to retrieve information.
 * @param i Retrieve information when long in n<tt>i</tt> for \c i in
 *        <tt>{0, 1, 2}</tt>.
 *
 * @return a valid underling_extents structure on success.  On failure,
 *         calls suzerain_error and returns UNDERLING_EXTENTS_INVALID.
 *
 * @see The method underling_local for a way to obtain only a subset of
 *      this information, or for a more Fortran-ready interface.
 */
underling_extents
underling_local_extents(
        const underling_problem problem,
        int i);

/**
 * Obtain the processor-local sizes, storage details, and global starting
 * offsets when the data is long in the n<tt>i</tt> direction.  This is
 * identical to the data obtainable via underling_local_extents but is provided
 * in a more Fortran-ready interface.  All strides and sizes are given in units
 * of underling_real.
 *
 * @param[in] problem Problem for which to retrieve information.
 * @param[in] i Retrieve information when long in n<tt>i</tt> for \c i in
 *            <tt>{0, 1, 2}</tt>.
 * @param[in,out] start If non-NULL on entry, contains
 *                underling_extents.start on successful return.
 * @param[in,out] size If non-NULL on entry, contains
 *                underling_extents.size on successful return.
 * @param[in,out] stride If non-NULL on entry, contains
 *                underling_extents.stride on successful return.
 * @param[in,out] order If non-NULL on entry, contains
 *                underling_extents.order on successful return.
 *
 * @return The value of underling_extents.extent on successful return.
 *         On error calls suzerain_error and returns 0.  Note that a
 *         0 return value may not indicate an error as the processor
 *         may have zero storage in some unbalanced decompositions.
 *
 * @see The method underling_local_extents for a more C-friendly and
 *      const-correct capable way to obtain all of this information.
 */
size_t
underling_local(
        const underling_problem problem,
        int i,
        int *start,
        int *size,
        int *stride,
        int *order);

/**
 * Find the amount of local storage necessary to plan and execute a problem.
 * This includes the storage to hold data in the long in <tt>n{0,1,2}</tt>
 * directions as well as any additional communication buffer space required.
 *
 * @param problem Problem for which to retrieve information.
 *
 * @return On success, the amount of storage required in units of
 *         underling_real.  On failure, calls suzerain_error and returns zero.
 */
size_t
underling_local_memory(
        const underling_problem problem);

/**
 * Find the amount of local storage necessary to merely hold a problem in
 * memory.  This includes holding the data long in the <tt>n{0,1,2}</tt>
 * directions but not any additional communication buffer overhead.
 *
 * @param problem Problem for which to retrieve information.
 *
 * @return On success, the amount of storage required in units of
 *         underling_real.  On failure, calls suzerain_error and returns zero.
 */
size_t
underling_local_memory_optimum(
        const underling_problem problem);

/**
 * Collectively find the maximum amount of per-processor memory required
 * to handle the problem across all processors in the grid.
 *
 * @param grid Grid for which to retrieve information.
 * @param problem Problem for which to retrieve information.
 *
 * @return On success, the maximal memory per-processor required to
 *         handle the problem.  On failure, calls suzerain_error
 *         and returns zero.
 */
size_t
underling_local_memory_maximum(
        const underling_grid    grid,
        const underling_problem problem);

/**
 * Collectively find the minimal amount of per-processor memory required
 * to handle the problem across all processors in the grid.
 *
 * @param grid Grid for which to retrieve information.
 * @param problem Problem for which to retrieve information.
 *
 * @return On success, the minimal memory per-processor required to
 *         handle the problem.  On failure, calls suzerain_error
 *         and returns zero.
 */
size_t
underling_local_memory_minimum(
        const underling_grid    grid,
        const underling_problem problem);

/**
 * Collectively find the global amount of memory required to
 * handle the problem across all processors in the grid.
 *
 * @param grid Grid for which to retrieve information.
 * @param problem Problem for which to retrieve information.
 *
 * @return On success, the global memory across all processors
 *         required to handle the problem.  On failure, calls
 *         suzerain_error and returns zero.
 */
size_t
underling_global_memory(
        const underling_grid    grid,
        const underling_problem problem);

/**
 * Find the theoretical optimum (minimum) amount of memory required to
 * store the given problem on all processors on the given grid assuming no
 * communications buffer overhead.
 *
 * @param grid Grid for which to retrieve information.
 * @param problem Problem for which to retrieve information.
 *
 * @return On success, the minimal global memory required to handle
 *         the problem.  On failure, calls suzerain_error and returns zero.
 */
size_t
underling_global_memory_optimum(
        const underling_grid    grid,
        const underling_problem problem);

/**
 * Collectively create an execution plan to solve the given decomposition
 * problem on the given data.  Creating a plan may have significant cost.  Once
 * a plan is created, it may be repeatedly used without incurring this one time
 * overhead.
 *
 * Planning cost can be reduced by only requesting the transform capabilities
 * you require using \c transform_flags.  It should contain the bitwise OR of
 * one or more of the following:
 *  - UNDERLING_TRANSPOSE_LONG_N2_TO_LONG_N1
 *  - UNDERLING_TRANSPOSE_LONG_N1_TO_LONG_N0
 *  - UNDERLING_TRANSPOSE_LONG_N0_TO_LONG_N1
 *  - UNDERLING_TRANSPOSE_LONG_N1_TO_LONG_N2
 *  - UNDERLING_TRANSPOSE_ALL
 *
 * It is an error to execute a transpose when the corresponding flag was
 * not provided to this method.
 *
 * Planning cost may be modified using \c fftw_rigor_flags.  It must
 * be one of the following:
 *  - FFTW_ESTIMATE
 *  - FFTW_MEASURE
 *  - FFTW_PATIENT
 *  - FFTW_EXHAUSTIVE
 *
 * Longer planning will usually result in shorter execution time and higher
 * performance.  See the <a
 * href="http://www.fftw.org/fftw3.3alpha_doc/Planner-Flags.html"> FFTW manual
 * regarding planner flags</a> for more details.  Note that using any value
 * other than FFTW_ESTIMATE will cause \c data to be overwritten during the
 * invocation.
 *
 * @param problem Problem for which a plan is sought.
 * @param data Storage of size at least <tt>underling_local_memory(problem)</tt>
 *             <tt>underling_real</tt>s in which all data movement occurs.
 * @param transform_flags Desired transforms to plan.  Specifying zero
 *                        is equivalent to specifying UNDERLING_TRANSPOSE_ALL.
 * @param fftw_rigor_flags Desired FFTW planning rigor.  Specifying zero
 *                         is equivalent to specifying FFTW_MEASURE.
 *
 * @return On success, returns a valid underling_plan suitable for execution.
 *         On failure, calls suzerain_error and returns NULL.
 *
 * @see The method underling_plan_destroy for how to destroy an instance.
 */
underling_plan
underling_plan_create(
        const underling_problem problem,
        underling_real * data,
        unsigned transform_flags,
        unsigned fftw_rigor_flags);

/**
 * Destroy all resources associated with the given plan.
 *
 * @param plan Plan to be destroyed.
 */
void
underling_plan_destroy(
        underling_plan plan);

/**
 * Collectively transform the data provided at plan creation time from being
 * long in \c n2 to long in \c n1.  Appropriate MPI calls and data reordering
 * will occur.
 *
 * @param plan Plan to be executed.
 *
 * @return SUZERAIN_SUCCESS (zero) on success and non-zero on failure.
 *
 * @see The documentation for underling_grid_create for more details on the
 *      associated storage orders.
 * @see The methods underling_local_extents or underling_local for how to
 *      obtain local storage details in either the input or output layout.
 */
int
underling_execute_long_n2_to_long_n1(
        const underling_plan plan);

/**
 * Collectively transform the data provided at plan creation time from being
 * long in \c n1 to long in \c n0.  Appropriate MPI calls and data reordering
 * will occur.
 *
 * @param plan Plan to be executed.
 *
 * @return SUZERAIN_SUCCESS (zero) on success and non-zero on failure.
 *
 * @see The documentation for underling_grid_create for more details on the
 *      associated storage orders.
 * @see The methods underling_local_extents or underling_local for how to
 *      obtain local storage details in either the input or output layout.
 */
int
underling_execute_long_n1_to_long_n0(
        const underling_plan plan);

/**
 * Collectively transform the data provided at plan creation time from being
 * long in \c n0 to long in \c n1.  Appropriate MPI calls and data reordering
 * will occur.
 *
 * @param plan Plan to be executed.
 *
 * @return SUZERAIN_SUCCESS (zero) on success and non-zero on failure.
 *
 * @see The documentation for underling_grid_create for more details on the
 *      associated storage orders.
 * @see The methods underling_local_extents or underling_local for how to
 *      obtain local storage details in either the input or output layout.
 */
int
underling_execute_long_n0_to_long_n1(
        const underling_plan plan);

/**
 * Collectively transform the data provided at plan creation time from being
 * long in \c n1 to long in \c n2.  Appropriate MPI calls and data reordering
 * will occur.
 *
 * @param plan Plan to be executed.
 *
 * @return SUZERAIN_SUCCESS (zero) on success and non-zero on failure.
 *
 * @see The documentation for underling_grid_create for more details on the
 *      associated storage orders.
 * @see The methods underling_local_extents or underling_local for how to
 *      obtain local storage details in either the input or output layout.
 */
int
underling_execute_long_n1_to_long_n2(
        const underling_plan plan);

/**
 * Dump an instance's internals in a debugging-friendly format.
 *
 * @param grid Grid to dump.
 * @param output_file Desired output handle,
 *                    which may be \c stdout or \c stderr.
 */
void
underling_fprint_grid(
        const underling_grid grid,
        FILE *output_file);

/**
 * Dump an instance's internals in a debugging-friendly format.
 *
 * @param problem Problem to dump.
 * @param output_file Desired output handle,
 *                    which may be \c stdout or \c stderr.
 */
void
underling_fprint_problem(
        const underling_problem problem,
        FILE *output_file);

/**
 * Dump an instance's internals in a debugging-friendly format.
 *
 * @param extents Extents to dump.
 * @param output_file Desired output handle,
 *                    which may be \c stdout or \c stderr.
 */
void
underling_fprint_extents(
        const underling_extents *extents,
        FILE *output_file);

/**
 * Dump an instance's internals in a debugging-friendly format.
 *
 * @param plan Plan to dump.
 * @param output_file Desired output handle,
 *                    which may be \c stdout or \c stderr.
 */
void
underling_fprint_plan(
        const underling_plan plan,
        FILE *output_file);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
__END_DECLS
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif // __SUZERAIN_UNDERLING_H
