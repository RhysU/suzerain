/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
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
#include <suzerain/underling_error.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS   /* empty */
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
    int global_offset;
    underling_dimension * transformed;
};

typedef struct {
    int                   ndim;
    int                 * global_size;
    double              * dealias_by;
    underling_state     * state;
    underling_dimension * dim_r;
    underling_dimension * dim_c;
} underling_workspace;

typedef double         underling_real;
typedef underling_real underling_complex[2];

underling_workspace * underling_workspace_alloc(int ndim);
void                  underling_workspace_free(underling_workspace * w);

__END_DECLS

#endif // PECOS_SUZERAIN_UNDERLING_H
