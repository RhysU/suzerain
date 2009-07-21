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
#ifndef PECOS_SUZERAIN_UNDERLING_H
#define PECOS_SUZERAIN_UNDERLING_H

#include <stdlib.h>
#include <suzerain/error.h>

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

__BEGIN_DECLS

typedef enum {
    underling_state_uninitialized  = 0,
    underling_state_physical       = 1,
    underling_state_wave           = 2,
    underling_state_nottransformed = 4
} underling_state;

typedef struct underling_dimension underling_dimension;
struct underling_dimension {
    int size;
    int stride;
    int global_size;
    int global_start;
    double dealias_by;
    underling_dimension *transformed;
};

typedef struct {
    int                  ndim;
    underling_state     *state;
    underling_dimension *dim_p;
    underling_dimension *dim_w;
} underling_workspace;

typedef double         underling_real;
typedef underling_real underling_complex[2];

underling_workspace *underling_workspace_alloc(int ndim);
void                 underling_workspace_free(underling_workspace *w);

int underling_prepare_physical_size(underling_workspace *w,
                                    const int *physical_size);
int underling_prepare_link(underling_workspace *w,
                           int dim_physical,
                           int dim_wave);
int underling_prepare_state(underling_workspace *w,
                            const underling_state *state);

__END_DECLS

#endif // PECOS_SUZERAIN_UNDERLING_H
