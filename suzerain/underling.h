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

typedef struct underling_problem {
    int howmany;
} underling_problem;

underling_problem *
underling_problem_create(
        underling_grid *grid,
        int howmany);

void
underling_problem_destroy(
        underling_problem * problem);

size_t
underling_size_local(
        underling_grid * grid,
        underling_problem * problem);

typedef struct underling_plan {
    underling_problem * p;
    int howmany;
    fftw_plan plan12;
    fftw_plan plan23;
    fftw_plan plan32;
    fftw_plan plan21;
    underling_real *data;
} underling_plan;

underling_plan *
underling_plan_create(
        underling_problem * problem,
        int howmany,
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
