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

#include <assert.h>
#include <stdlib.h>
#include <suzerain/underling.h>

underling_workspace *
underling_workspace_alloc(int ndim)
{
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

    w->state = malloc(w->ndim * sizeof(w->state[0]));
    if (w->state == NULL) {
        UNDERLING_ERROR_NULL("failed to allocate space for state",
                             UNDERLING_ENOMEM);
        free(w);
    }

    w->dim_p = malloc(w->ndim * sizeof(w->dim_p[0]));
    if (w->dim_p == NULL) {
        UNDERLING_ERROR_NULL("failed to allocate space for dim_p",
                             UNDERLING_ENOMEM);
        free(w->state);
        free(w);
    }

    w->dim_w = malloc(w->ndim * sizeof(w->dim_w[0]));
    if (w->dim_w == NULL) {
        UNDERLING_ERROR_NULL("failed to allocate space for dim_w",
                             UNDERLING_ENOMEM);
        free(w->state);
        free(w->dim_p);
        free(w);
    }

    for (i = 0; i < w->ndim; ++i) {
        w->state[i]               = underling_state_uninitialized;

        w->dim_p[i].size          = 0;
        w->dim_p[i].stride        = 0;
        w->dim_p[i].global_size   = 0;
        w->dim_p[i].global_start  = 0;
        w->dim_p[i].dealias_by    = 3.0/2.0;
        w->dim_p[i].transformed   = NULL;

        w->dim_w[i].size          = 0;
        w->dim_w[i].stride        = 0;
        w->dim_w[i].global_size   = 0;
        w->dim_w[i].global_start  = 0;
        w->dim_p[i].dealias_by    = 1.0;
        w->dim_w[i].transformed   = NULL;
    }

    return w;
}

void
underling_workspace_free(underling_workspace * w)
{
    free(w->state);
    free(w->dim_p);
    free(w->dim_w);
    free(w);
}

int
underling_prepare_physical_size(underling_workspace *w,
                                const int *physical_size)
{
    int i;

    for (i = 0; i < w->ndim; ++i) {
        if (physical_size[i] < 1) {
            UNDERLING_ERROR("physical_size < 1", UNDERLING_EINVAL);
        }
    }

    for (i = 0; i < w->ndim; ++i) {
        w->dim_p[i].global_size = physical_size[i];
    }

    return UNDERLING_SUCCESS;
}


int
underling_prepare_link(underling_workspace *w,
                       int dim_physical,
                       int dim_wave)
{
    if (dim_physical < 0 || dim_physical >= w->ndim) {
        UNDERLING_ERROR("dim_physical index out of range", UNDERLING_EINVAL);
    }
    if (w->dim_p[dim_physical].transformed != NULL) {
        UNDERLING_ERROR("dim_physical already linked", UNDERLING_EINVAL);
    }

    if (dim_wave < 0 || dim_wave >= w->ndim) {
        UNDERLING_ERROR("dim_wave index out of range", UNDERLING_EINVAL);
    }
    if (w->dim_w[dim_wave].transformed != NULL) {
        UNDERLING_ERROR("dim_wave already linked", UNDERLING_EINVAL);
    }

    w->dim_p[dim_physical].transformed = &(w->dim_w[dim_wave]);
    w->dim_w[dim_wave].transformed     = &(w->dim_p[dim_physical]);

    return UNDERLING_SUCCESS;
}

int
underling_prepare_state(underling_workspace *w,
                        const underling_state *state)
{
    int i;

    for (i = 0; i < w->ndim; ++i) {
        w->state[i] = state[i];
    }

    return UNDERLING_SUCCESS;
}
