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
 * underling.h: A parallel, three dimensional FFT library atop MPI
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_UNDERLING_H
#define __SUZERAIN_UNDERLING_H

#include <suzerain/common.h>
#include <suzerain/error.h>
#include <mpi.h>
#include <fftw3-mpi.h>

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

/* TODO Add appropriate grid/problem/plan const-ness */

typedef double underling_real;     /**< Real-valued scalar */
typedef struct underling_grid_s    *underling_grid;
typedef struct underling_problem_s *underling_problem;
typedef struct underling_plan_s    *underling_plan;

typedef struct underling_extents {
    int    start[3];
    int    size[3];
    int    stride[3];
    size_t total_extent;
} underling_extents;

extern const underling_extents UNDERLING_EXTENTS_INVALID;

/** Flag indicating a transform from long in \c n0 to long in \c n2. */
#define UNDERLING_DIRECTION_FORWARD  (1U << 0)
/** Flag indicating a transform from long in \c n2 to long in \c n0. */
#define UNDERLING_DIRECTION_BACKWARD (1U << 1)
/** Convenience flag indicating transforming in both directions */
#define UNDERLING_DIRECTION_BOTH \
        (UNDERLING_DIRECTION_FORWARD | UNDERLING_DIRECTION_BACKWARD)

underling_grid
underling_grid_create(
        MPI_Comm comm,
        int n0,
        int n1,
        int n2,
        int pA,
        int pB);

void
underling_grid_destroy(
        underling_grid grid);

underling_problem
underling_problem_create(
        underling_grid grid,
        int howmany);

void
underling_problem_destroy(
        underling_problem problem);

size_t
underling_local(
        const underling_problem problem,
        int n,
        int *start,
        int *size,
        int *stride);

underling_extents
underling_local_extents(
        const underling_problem problem,
        int n);

size_t
underling_local_memory(
        const underling_problem problem);

size_t
underling_local_memory_optimum(
        const underling_grid    grid,
        const underling_problem problem);

size_t
underling_local_memory_maximum(
        const underling_grid    grid,
        const underling_problem problem);

size_t
underling_local_memory_minimum(
        const underling_grid    grid,
        const underling_problem problem);

size_t
underling_global_memory(
        const underling_grid    grid,
        const underling_problem problem);

size_t
underling_global_memory_optimum(
        const underling_grid    grid,
        const underling_problem problem);

underling_plan
underling_plan_create(
        const underling_problem problem,
        underling_real * data,
        unsigned direction_flags,
        unsigned fftw_rigor_flags);

void
underling_plan_destroy(
        underling_plan plan);

int
underling_execute_backward(
        const underling_plan plan);

int
underling_execute_forward(
        const underling_plan plan);

void
underling_fprint_grid(
        const underling_grid grid,
        FILE *output_file);

void
underling_fprint_problem(
        const underling_problem problem,
        FILE *output_file);

void
underling_fprint_extents(
        const underling_extents *extents,
        FILE *output_file);

void
underling_fprint_plan(
        const underling_plan plan,
        FILE *output_file);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
__END_DECLS
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif // __SUZERAIN_UNDERLING_H
