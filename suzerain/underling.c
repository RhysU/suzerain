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
 * underling.c: A parallel, three dimensional FFT library atop MPI
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <stdlib.h>
#include <suzerain/underling.h>

underling_workspace *
underling_workspace_alloc(int ndim) {
    int i;
    underling_workspace * w;

    if (ndim < 1) {
        UNDERLING_ERROR_NULL("ndim must be at least 1", UNDERLING_EINVAL);
    }

    w = malloc(sizeof(underling_workspace));
    if (w == NULL) {
        UNDERLING_ERROR_NULL("failed to allocate space for workspace",
                             UNDERLING_ENOMEM);
    }

    w->ndim = ndim;

    w->global_size = malloc(w->ndim * sizeof(w->global_size[0]));
    if (w->global_size == NULL) {
        UNDERLING_ERROR_NULL("failed to allocate space for global_size",
                             UNDERLING_ENOMEM);
        free(w);
    }

    w->dealias_by = malloc(w->ndim * sizeof(w->dealias_by[0]));
    if (w->dealias_by == NULL) {
        UNDERLING_ERROR_NULL("failed to allocate space for dealias_by",
                             UNDERLING_ENOMEM);
        free(w->global_size);
        free(w);
    }

    w->state = malloc(w->ndim * sizeof(w->state[0]));
    if (w->state == NULL) {
        UNDERLING_ERROR_NULL("failed to allocate space for state",
                             UNDERLING_ENOMEM);
        free(w->global_size);
        free(w->dealias_by);
        free(w);
    }

    w->dim_r = malloc(w->ndim * sizeof(w->dim_r[0]));
    if (w->dim_r == NULL) {
        UNDERLING_ERROR_NULL("failed to allocate space for dim_r",
                             UNDERLING_ENOMEM);
        free(w->global_size);
        free(w->dealias_by);
        free(w->state);
        free(w);
    }

    w->dim_c = malloc(w->ndim * sizeof(w->dim_c[0]));
    if (w->dim_c == NULL) {
        UNDERLING_ERROR_NULL("failed to allocate space for dim_c",
                             UNDERLING_ENOMEM);
        free(w->global_size);
        free(w->dealias_by);
        free(w->state);
        free(w->dim_r);
        free(w);
    }

    for (i = 0; i < w->ndim; ++i) {
        w->global_size[i]         = 0;
        w->dealias_by[i]          = 3.0/2.0;
        w->state[i]               = underling_state_uninitialized;
        w->dim_r[i].size          = 0;
        w->dim_r[i].stride        = 0;
        w->dim_r[i].global_offset = 0;
        w->dim_r[i].transformed   = NULL;
        w->dim_c[i].size          = 0;
        w->dim_c[i].stride        = 0;
        w->dim_c[i].global_offset = 0;
        w->dim_c[i].transformed   = NULL;
    }

    return w;
}

void
underling_workspace_free(underling_workspace * w) {
    free(w->global_size);
    free(w->dealias_by);
    free(w->state);
    free(w->dim_r);
    free(w->dim_c);
    free(w);
}
