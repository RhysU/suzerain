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

typedef double         underling_real;
typedef underling_real underling_complex[2];

typedef struct underling_grid {
    int np0;
    int nw0;
    int n1;
    int n2;
    int p0;
    int p1;
    int block_a;
    MPI_Comm g_comm;
    MPI_Comm p0_comm;
    MPI_Comm p1_comm;
} underling_grid;

underling_grid *
underling_grid_create(
        MPI_Comm comm,
        int n0,
        int n1,
        int n2,
        int p0,
        int p1);

void
underling_grid_destroy(
        underling_grid * grid);

typedef struct underling_transpose_details {
    ptrdiff_t n[2];
    ptrdiff_t howmany;                        // # of real-valued components
    MPI_Comm comm;                            // underling_grid owns resource
    unsigned flags;
    ptrdiff_t local_n0;
    ptrdiff_t local_n0_start;
    ptrdiff_t local_n1;
    ptrdiff_t local_n1_start;
    ptrdiff_t local_size;
} underling_transpose_details;

typedef struct underling_problem {
    int nfields;                              // # of complex-valued state
    underling_transpose_details tophysical_A; // n2 long to n1 long
    underling_transpose_details tophysical_B; // n1 long to n0 long
    underling_transpose_details towave_B;     // n0 long to n1 long
    underling_transpose_details towave_A;     // n1 long to n2 long
    ptrdiff_t local_size;                     // Max of all local sizes
} underling_problem;

underling_problem *
underling_problem_create(
        underling_grid *grid,
        int nfields);

void
underling_problem_destroy(
        underling_problem * problem);

size_t
underling_local_size(
        underling_problem * problem);

typedef struct underling_plan {
    underling_problem * p;             // underling_problem owns resources
    fftw_plan transpose_tophysical_A;
    fftw_plan c2c_tophysical_n1;
    fftw_plan transpose_tophysical_B;
    fftw_plan c2r_tophysical_n0;
    fftw_plan r2c_towave_n0;
    fftw_plan transpose_towave_A;
    fftw_plan c2c_towave_n1;
    fftw_plan transpose_towave_B;
    underling_real *data;              // API end-user owns resources
} underling_plan;

underling_plan *
underling_plan_create(
        underling_problem * problem,
        underling_real * data,
        int will_perform_c2r,
        int will_perform_r2c,
        unsigned flags);

void
underling_plan_destroy(
        underling_plan * plan);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
__END_DECLS
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif // __SUZERAIN_UNDERLING_H
