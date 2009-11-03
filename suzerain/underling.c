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
 * underling.c: A parallel, three dimensional FFT library atop MPI
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <config.h>
#include <suzerain/common.h>
#pragma hdrstop
#include "gl_array_list.h"
#include <gsl/gsl_combination.h>
#include <suzerain/error.h>
#include <suzerain/underling.h>

bool
scalar_to_physical_gl_listelement_equals_fn(const void *elt1,
                                            const void *elt2)
{
    const underling_scalar_to_physical * const stp1
        = (underling_scalar_to_physical *) elt1;
    const underling_scalar_to_physical * const stp2
        = (underling_scalar_to_physical *) elt2;

    const char * const name1
        = !(stp1) ? NULL : !(stp1->name) ? NULL : stp1->name[0];
    const char * const name2
        = !(stp2) ? NULL : !(stp2->name) ? NULL : stp2->name[0];

    return 0 == strcmp(name1, name2);
}

void
scalar_to_physical_free(underling_scalar_to_physical * stp)
{
    /* FIXME: Implement */
}

underling_scalar_to_physical *
scalar_to_physical_alloc(const underling_workspace * const w,
                         const char * const field_name,
                         const int nderivative)
{
    underling_scalar_to_physical * const retval =
        calloc(1, sizeof(underling_scalar_to_physical));
    if (retval == NULL) {
        scalar_to_physical_free(retval);
        SUZERAIN_ERROR_NULL("failed to allocate space for scalar_to_physical",
                             SUZERAIN_ENOMEM);
    }
    retval->stage       = w->nstage-1; /* Start in all wave space */
    retval->nderivative = nderivative;
    retval->nfield      = 0;

    retval->index = calloc(nderivative+1, sizeof(retval->index[0]));
    if (retval == NULL) {
        scalar_to_physical_free(retval);
        SUZERAIN_ERROR_NULL("failed to allocate space for scalar_to_physical",
                             SUZERAIN_ENOMEM);
    }
    for (int i = 0; i <= nderivative; ++i) {
        /* FIXME STARTHERE */
        retval->index[i] = gsl_combination_calloc(w->ndim, i+1);
    }

    return retval;
}

underling_workspace *
underling_workspace_alloc(const int ndim, const int nstage)
{
    if (ndim < 1) {
        SUZERAIN_ERROR_NULL("ndim must be at least 1", SUZERAIN_EINVAL);
    }
    if (nstage < 1) {
        SUZERAIN_ERROR_NULL("nstage must be at least 1", SUZERAIN_EINVAL);
    }
    if (nstage > ndim+1) {
        SUZERAIN_ERROR_NULL("nstage must be less than ndim+1", SUZERAIN_EINVAL);
    }

    underling_workspace * w = calloc(1, sizeof(underling_workspace));
    if (w == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate space for workspace",
                             SUZERAIN_ENOMEM);
    }
    w->ndim   = ndim;
    w->nstage = nstage;

    w->stage = calloc(nstage, sizeof(w->stage[0]));
    if (w->stage == NULL) {
        underling_workspace_free(w);
        SUZERAIN_ERROR_NULL("failed to allocate space for stage",
                             SUZERAIN_ENOMEM);
    }

    for (int i = 0; i < nstage; ++i) {
        w->stage[i].dim = calloc(ndim, sizeof(w->stage[0].dim[0]));
        if (w->stage[i].dim == NULL) {
            underling_workspace_free(w);
            SUZERAIN_ERROR_NULL("failed to allocate space for stage's dim",
                                SUZERAIN_ENOMEM);
        }
    }

    /* In stage 0, ndim - (nstage-1) dimensions are in physical space. */
    for (int j = 0; j < ndim; ++j) {
        if (j < nstage-1) {
            w->stage[0].dim[j].state = UNDERLING_STATE_PHYSICAL;
        } else {
            w->stage[0].dim[j].state = UNDERLING_STATE_NOTTRANSFORMED;
        }
    }

    for (int i = 0; i < nstage-1; ++i) {
        for (int j = 0; j < ndim; ++j) {
            /* Establish pointers to this dimension in next and previous stage */
            const int j_next_r2c = (j + ndim - 1) % ndim;
            const int j_next_c2r = (j + 1) % ndim;
            w->stage[i].dim[j].next_r2c = w->stage[i+1].dim + j_next_r2c;
            w->stage[i+1].dim[j].next_c2r = w->stage[i].dim + j_next_c2r;

            /* Assume the dimension stays in the same state... */
            w->stage[i].dim[j].next_r2c->state = w->stage[i].dim[j].state;
        }

        /* ...unless it was the zeroth dimension in this stage */
        switch (w->stage[i].dim[0].state) {
            case UNDERLING_STATE_PHYSICAL:
                w->stage[i].dim[0].next_r2c->state = UNDERLING_STATE_WAVE;
                break;
            case UNDERLING_STATE_NOTTRANSFORMED:
                /* NOP */
                break;
            default:
                SUZERAIN_ERROR_NULL("unexpected state in allocation logic",
                                    SUZERAIN_ESANITY);
        }
    }

    return w;
}

void
underling_workspace_free(underling_workspace * w)
{
    if (w) {
        if (w->stage) {
            for (int i = 0; i < w->nstage; ++i) {
                if (w->stage[i].dim) {
                    /* Free items allocated per-dimension  */
                    for (int j = 0; j < w->ndim; ++j) {
                        if (w->stage[i].dim[j].name) {
                            free(w->stage[i].dim[j].name);
                            w->stage[i].dim[j].name = NULL;
                        }
                    }
                    /* Free dimensions previously allocated in a block */
                    free(w->stage[i].dim);
                    w->stage[i].dim = NULL;
                }
            }
            free(w->stage);
            w->stage = NULL;
        }
        free(w);
    }
}

int
underling_name_dimension(underling_workspace * const w,
                         const int ndim,
                         const int nstage,
                         const char * const name)
{
    if (ndim < 0 || ndim > w->ndim+1) {
        SUZERAIN_ERROR("ndim out of range", SUZERAIN_EINVAL);
    }
    if (nstage < 0 || nstage > w->nstage+1) {
        SUZERAIN_ERROR("nstage out of range", SUZERAIN_EINVAL);
    }

    underling_dimension * dim;

    dim = &(w->stage[nstage].dim[ndim]);
    while (dim) {
        dim->name = strdup(name);
        if (! dim->name) {
            SUZERAIN_ERROR("Unable to allocate space for name",
                           SUZERAIN_ENOMEM);
        }

        dim = dim->next_c2r;
    }

    dim = w->stage[nstage].dim[ndim].next_r2c;
    while (dim) {
        dim->name = strdup(name);
        if (! dim->name) {
            SUZERAIN_ERROR("Unable to allocate space for name",
                           SUZERAIN_ENOMEM);
        }

        dim = dim->next_r2c;
    }

    return SUZERAIN_SUCCESS;
}


int
underling_scalar_to_physical_add(underling_workspace * const w,
                                 const char * const name,
                                 const int max_derivative)
{
    if (max_derivative < 0) {
        SUZERAIN_ERROR("max_derivative must be nonnegative", SUZERAIN_EINVAL);
    }
    if (max_derivative > 2) {
        SUZERAIN_ERROR("Unable to handle max_derivative > 2", SUZERAIN_ESANITY);
    }
    if (!name || !strlen(name)) {
        SUZERAIN_ERROR("Unable to add a nameless scalar_to_physical field",
                       SUZERAIN_EINVAL);
    }

    /* Go to end of scalar_to_physical list, looking for duplicate names */

    return SUZERAIN_SUCCESS;
}

int
underling_prepare_physical_size(underling_workspace *w,
                                const int *physical_size)
{
    for (int i = 0; i < w->ndim; ++i) {
        if (physical_size[i] < 1) {
            SUZERAIN_ERROR("physical_size < 1", SUZERAIN_EINVAL);
        }
    }

    // FIXME
/*    for (int i = 0; i < w->ndim; ++i) {*/
/*        w->dim_p[i].global_size = physical_size[i];*/
/*    }*/

    return SUZERAIN_SUCCESS;
}
