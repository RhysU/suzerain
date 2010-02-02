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
        int i,
        underling_real * data,
        unsigned fftw_sign,
        unsigned fftw_rigor_flags)
{
    // Sanity check input arguments
    if (SUZERAIN_UNLIKELY(problem == NULL)) {
        SUZERAIN_ERROR_VAL("problem == NULL", SUZERAIN_EINVAL, 0);
    }
    if (SUZERAIN_UNLIKELY(i < 0 || i > 2)) {
        SUZERAIN_ERROR_VAL("i < 0 or i > 2", SUZERAIN_EINVAL, 0);
    }
    const underling_extents extents = underling_local_extents(problem, i);
    if (SUZERAIN_UNLIKELY(extents.stride[extents.strideorder[0]] % 2)) {
        SUZERAIN_ERROR_NULL(
                "problem must have an even number of underling_real fields",
                SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(data == NULL)) {
        SUZERAIN_ERROR_NULL("data == NULL", SUZERAIN_EINVAL);
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

    // Prepare the input to fftw_plan_guru_split_dft.  FFTW split interface
    // allows using underling_extents.strides directly.  The tranform is purely
    // in place which sets our output strides equal to our input strides.

    // We transform the long dimension given by extents.strideorder[0]
    const fftw_iodim dims[] = {
        extents.size[extents.strideorder[0]],   // n
        extents.stride[extents.strideorder[0]], // is
        extents.stride[extents.strideorder[0]]  // os
    };
    const int rank = sizeof(dims)/sizeof(dims[0]);

    // We loop over the slowest direction, the second slowest
    // direction, and the individual state fields in row-major
    // order.
    const fftw_iodim howmany_dims[3] = {
        {
            extents.size[extents.strideorder[2]],
            extents.stride[extents.strideorder[2]],
            extents.stride[extents.strideorder[2]]
        },
        {
            extents.size[extents.strideorder[1]],
            extents.stride[extents.strideorder[1]],
            extents.stride[extents.strideorder[1]]
        },
        {
            extents.stride[extents.strideorder[0]] / 2, // howmany/2
            2,                                          // is, interleaved
            2                                           // os, interleaved
        }
    };
    const int howmany_rank = sizeof(howmany_dims)/sizeof(howmany_dims[0]);

    // For interleaved storage, the imaginary data starts one underling_real
    // after data itself.  Per FFTW manual section 4.5.3, FFTW_BACKWARD is
    // FFTW_FORWARD with the real and imaginary parts flipped.
    underling_real * const ri = (fftw_sign == FFTW_FORWARD) ? data : data + 1;
    underling_real * const ii = (fftw_sign == FFTW_FORWARD) ? data + 1 : data;
    underling_real * const ro = ri;
    underling_real * const io = ii;

    fftw_plan fftwp = fftw_plan_guru_split_dft(rank, dims,
                                               howmany_rank, howmany_dims,
                                               ri, ii, ro, io,
                                               fftw_rigor_flags);

    if (SUZERAIN_UNLIKELY(fftwp == NULL)) {
        SUZERAIN_ERROR_NULL("FFTW returned a NULL plan", SUZERAIN_ESANITY);
    }

    // Create and initialize the fftplan workspace
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
