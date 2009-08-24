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
 * richardson.c: generalized Richardson extrapolation routines
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include "config.h"

#include <stdlib.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_nan.h>
#include <suzerain/error.h>
#include <suzerain/richardson.h>

int
suzerain_richardson_extrapolation_step(
    gsl_vector * Ah,
    const gsl_vector * Aht,
    const double t,
    const int ki)
{
    const double tki       = gsl_pow_int(t, ki);
    const double inv_tkim1 = 1.0/(tki-1.0);

    gsl_blas_dscal(-inv_tkim1, Ah);

    int error = gsl_blas_daxpy(tki*inv_tkim1, Aht, Ah);
    if (error) return error;

    return SUZERAIN_SUCCESS;
}

int
suzerain_richardson_extrapolation(
        gsl_matrix * const A,
        const double t,
        const gsl_vector_int * k,
        gsl_matrix * normtable,
        const gsl_vector * const exact)
{
    gsl_vector * scratch = NULL;
    const size_t k_size = (k) ? k->size : 0;

    if (exact) {
        if (!normtable) {
            SUZERAIN_ERROR("exact provided but normtable was not",
                    SUZERAIN_EINVAL);
        }
        if (exact->size != A->size1) {
            SUZERAIN_ERROR("exact->size does not match A->size1",
                    SUZERAIN_EINVAL);
        }
    }

    if (normtable) {
        if (normtable->size1 != A->size2) {
            SUZERAIN_ERROR("normtable->size1 does not match A->size2",
                    SUZERAIN_EINVAL);
        }
        if (normtable->size2 != A->size2) {
            SUZERAIN_ERROR("normtable->size2 does not match A->size2",
                    SUZERAIN_EINVAL);
        }

        gsl_matrix_set_all(normtable, GSL_NAN);

        /* Compute the first column in the norm table */
        if (exact) {
            scratch = gsl_vector_alloc(A->size1);
        }
        for (int i = 0; i < normtable->size2; ++i) {
            double norm;
            if (exact) {
                gsl_matrix_get_col(scratch, A, i);
                gsl_blas_daxpy(-1.0, exact, scratch);
                norm = gsl_blas_dnrm2(scratch);
            } else {
                gsl_vector_view Ai = gsl_matrix_column(A, i);
                norm = gsl_blas_dnrm2(&Ai.vector);
            }
            gsl_matrix_set(normtable, i, 0, norm);
        }
    }

    for (int i = 0; i < A->size2; ++i) {
        for (int j = 0; j < A->size2 - i - 1; ++j) {
            gsl_vector_view Aih  = gsl_matrix_column(A, j);
            gsl_vector_view Aiht = gsl_matrix_column(A, j+1);

            /* Provide automagic around the k parameter */
            int ki;
            if (k_size) {
                if (i < k_size) {
                    ki = gsl_vector_int_get(k, i);
                } else if (k_size == 1) {
                    ki = gsl_vector_int_get(k, 0) + i;
                } else {
                    const int increment = gsl_vector_int_get(k, k_size-1)
                        - gsl_vector_int_get(k, k_size-2);
                    ki = (i - (k_size-1))*increment;
                }
            } else {
                ki = i+1;
            }

            const int error = suzerain_richardson_extrapolation_step(
                    &Aih.vector, &Aiht.vector, t, ki);
            if (error) {
                if (scratch) {
                    gsl_vector_free(scratch);
                }
                return error;
            }

            if (normtable) {
                double norm;
                if (exact) {
                    gsl_vector_memcpy(scratch, &Aih.vector);
                    gsl_blas_daxpy(-1.0, exact, scratch);
                    norm = gsl_blas_dnrm2(scratch);
                } else {
                    norm = gsl_blas_dnrm2(&Aih.vector);
                }
                gsl_matrix_set(normtable, i+j+1, i+1, norm);
            }
        }
    }

    if (scratch) {
        gsl_vector_free(scratch);
    }

    return SUZERAIN_SUCCESS;
}
