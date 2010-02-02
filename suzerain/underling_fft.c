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
 * underling_fft.c: Convenience wrappers around FFTW-like planning
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.h>
#pragma hdrstop
#include <suzerain/error.h>
#include <suzerain/underling_fft.h>
#include <fftw3.h>

// *******************************************************************
// INTERNAL STRUCTS INTERNAL STRUCTS INTERNAL STRUCTS INTERNAL STRUCTS
// *******************************************************************

struct underling_fftplan_s {
    fftw_plan fftwp;
};

// **************************************************************************
// IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION
// **************************************************************************

underling_fftplan
underling_fftw_plan_c2c(
        const underling_problem problem,
        underling_real * data,
        int i,
        unsigned fftw_sign,
        unsigned fftw_rigor_flags)
{
    if (SUZERAIN_UNLIKELY(problem == NULL)) {
        SUZERAIN_ERROR_NULL("problem == NULL", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(data == NULL)) {
        SUZERAIN_ERROR_NULL("data == NULL", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(i < 0 || i > 2)) {
        SUZERAIN_ERROR_NULL("i < 0 or i > 2", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(   fftw_sign != FFTW_FORWARD
                          || fftw_sign != FFTW_BACKWARD)) {
        SUZERAIN_ERROR_NULL(
                "fftw_sign not one of FFTW_{FORWARD,BACKWARD}",
                SUZERAIN_EINVAL);
    }
    const unsigned non_rigor_mask =   ~FFTW_ESTIMATE
                                    & ~FFTW_MEASURE
                                    & ~FFTW_PATIENT
                                    & ~FFTW_EXHAUSTIVE
                                    & ~FFTW_WISDOM_ONLY;
    if (SUZERAIN_UNLIKELY(fftw_rigor_flags & non_rigor_mask)) {
        SUZERAIN_ERROR_NULL("FFTW non-rigor bits disallowed", SUZERAIN_EINVAL);
    }

    fftw_plan fftwp = NULL; // TODO Implement and sanity check

    // Create and initialize the grid workspace
    underling_fftplan f = calloc(1, sizeof(struct underling_fftplan_s));
    if (SUZERAIN_UNLIKELY(f == NULL)) {
        SUZERAIN_ERROR_NULL("failed to allocate space for fftplan",
                             SUZERAIN_ENOMEM);
    }
    // Copy the fftplan parameters to the fftplan workspace
    f->fftwp = fftwp;

    return f;
}

int
underling_fftplan_execute(
        const underling_fftplan fftplan)
{
    if (SUZERAIN_UNLIKELY(fftplan == NULL)) {
        SUZERAIN_ERROR("fftplan == NULL", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(fftplan->fftwp == NULL)) {
        SUZERAIN_ERROR("fftplan->fftwp == NULL", SUZERAIN_EINVAL);
    }

    fftw_execute(fftplan->fftwp);

    return SUZERAIN_SUCCESS;
}

void
underling_fftplan_destroy(
        underling_fftplan fftplan)
{
    if (fftplan) {
        if (fftplan->fftwp) {
            fftw_destroy_plan(fftplan->fftwp);
            fftplan->fftwp = NULL;
        }
        free(fftplan);
    }
}

void
underling_fprint_fftplan(
        const underling_fftplan fftplan,
        FILE *output_file)
{
    fprintf(output_file, "{underling_fftplan:");
    if (!fftplan) {
        fprintf(output_file, "NULL");
    } else {
        if (fftplan->fftwp) {
            fprintf(output_file, "{fftwp:");
            fftw_fprint_plan(fftplan->fftwp, output_file);
            fprintf(output_file, "}");
        }
    }
}
