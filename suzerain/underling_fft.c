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
    fftw_plan plan_reorder;
    fftw_plan plan_fft;
};

// ********************************************************************
// INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL
// ********************************************************************

fftw_plan
underling_fftw_plan_nop();

void rotate_left(int *array, int len);

underling_fftplan
underling_fftplan_create_c2c(
        const underling_problem problem,
        int long_ni,
        underling_real * data,
        int fftw_sign,
        unsigned fftw_rigor_flags);

static const unsigned non_rigor_mask =   ~FFTW_ESTIMATE
                                       & ~FFTW_MEASURE
                                       & ~FFTW_PATIENT
                                       & ~FFTW_EXHAUSTIVE
                                       & ~FFTW_WISDOM_ONLY;

// **************************************************************************
// IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION
// **************************************************************************

static
fftw_plan
underling_fftw_plan_nop()
{
    // Create a non-NULL NOP FFTW plan
    fftw_plan nop_plan = fftw_plan_guru_r2r(/*rank*/0,
                                            /*dims*/NULL,
                                            /*howmany_rank*/0,
                                            /*dims*/NULL,
                                            /*in*/NULL,
                                            /*out*/NULL,
                                            /*kind*/NULL,
                                            /*flags*/0);
    assert(nop_plan);
    return nop_plan;
}

void rotate_left(int *array, int len)
{
    if (len > 0) {
        const int tmp = array[0];
        for (int i = 0; i < len - 1; ++i) {
            array[i] = array[i+1];
        }
        array[len - 1] = tmp;
    }
}

underling_fftplan
underling_fftplan_create_c2c_forward(
        const underling_problem problem,
        int long_ni,
        underling_real * data,
        unsigned fftw_rigor_flags)
{
    return underling_fftplan_create_c2c(
            problem, long_ni, data, FFTW_FORWARD, fftw_rigor_flags);
}

underling_fftplan
underling_fftplan_create_c2c_backward(
        const underling_problem problem,
        int long_ni,
        underling_real * data,
        unsigned fftw_rigor_flags)
{
    return underling_fftplan_create_c2c(
            problem, long_ni, data, FFTW_BACKWARD, fftw_rigor_flags);
}

static
underling_fftplan
underling_fftplan_create_c2c(
        const underling_problem problem,
        int long_ni,
        underling_real * data,
        int fftw_sign,
        unsigned fftw_rigor_flags)
{
    // Sanity check input arguments
    if (SUZERAIN_UNLIKELY(problem == NULL)) {
        SUZERAIN_ERROR_VAL("problem == NULL", SUZERAIN_EINVAL, 0);
    }
    if (SUZERAIN_UNLIKELY(long_ni < 0 || long_ni > 2)) {
        SUZERAIN_ERROR_VAL("long_ni < 0 or long_ni > 2", SUZERAIN_EINVAL, 0);
    }
    const underling_extents extents = underling_local_extents(problem, long_ni);
    if (SUZERAIN_UNLIKELY(extents.size[3] % 2)) {
        SUZERAIN_ERROR_NULL(
                "problem must have an even number of underling_real fields",
                SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(extents.order[0] != 3)) {
        SUZERAIN_ERROR_NULL(
                "transformed fields not interleaved: extents.order[0] != 3",
                SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(data == NULL)) {
        SUZERAIN_ERROR_NULL("data == NULL", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(   fftw_sign != FFTW_FORWARD
                          && fftw_sign != FFTW_BACKWARD)) {
        SUZERAIN_ERROR_NULL(
                "fftw_sign not one of FFTW_{FORWARD,BACKWARD}",
                SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(fftw_rigor_flags & non_rigor_mask)) {
        SUZERAIN_ERROR_NULL("FFTW non-rigor bits disallowed", SUZERAIN_EINVAL);
    }

    // Prepare the reordering plan for the input data
    // TODO Fix ESANITY below by reordering for UNDERLING_TRANSPOSED_LONG_N2
    // TODO Fix ESANITY below by reordering for UNDERLING_TRANSPOSED_LONG_N0
    if (SUZERAIN_UNLIKELY(extents.order[1] != long_ni)) {
        SUZERAIN_ERROR_NULL(
                "transformed direction not long: extents.order[1] != long_ni",
                SUZERAIN_ESANITY);
    }
    fftw_plan plan_reorder = underling_fftw_plan_nop();
    if (SUZERAIN_UNLIKELY(plan_reorder == NULL)) {
        SUZERAIN_ERROR_NULL("FFTW returned a NULL NOP plan", SUZERAIN_ESANITY);
    }

    // Prepare the input to fftw_plan_guru_split_dft.  FFTW split interface
    // allows using underling_extents.strides directly.  The tranform is purely
    // in place which sets our output strides equal to our input strides.

    // We transform the long dimension given by extents.order[1]
    const fftw_iodim dims[] = {
        extents.size[extents.order[1]],   // n
        extents.stride[extents.order[1]], // is
        extents.stride[extents.order[1]]  // os
    };
    const int rank = sizeof(dims)/sizeof(dims[0]);

    // We loop over the slowest direction, the second slowest
    // direction, and the individual state fields in row-major
    // order.
    const fftw_iodim howmany_dims[3] = {
        {
            extents.size[extents.order[3]],
            extents.stride[extents.order[3]],
            extents.stride[extents.order[3]]
        },
        {
            extents.size[extents.order[2]],
            extents.stride[extents.order[2]],
            extents.stride[extents.order[2]]
        },
        {
            extents.size[3] / 2, // howmany/2
            2,                   // is, interleaved
            2                    // os, interleaved
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

    fftw_plan plan_fft = fftw_plan_guru_split_dft(rank, dims,
                                                  howmany_rank, howmany_dims,
                                                  ri, ii, ro, io,
                                                  fftw_rigor_flags);

    if (SUZERAIN_UNLIKELY(plan_fft == NULL)) {
        SUZERAIN_ERROR_NULL("FFTW returned a NULL plan", SUZERAIN_ESANITY);
    }

    // Create and initialize the fftplan workspace
    underling_fftplan f = calloc(1, sizeof(struct underling_fftplan_s));
    if (SUZERAIN_UNLIKELY(f == NULL)) {
        SUZERAIN_ERROR_NULL("failed to allocate space for fftplan",
                             SUZERAIN_ENOMEM);
    }
    // Copy the relevant parameters to the fftplan workspace
    f->plan_reorder = plan_reorder;
    f->plan_fft     = plan_fft;

    return f;
}

// FIXME Implement
underling_fftplan
underling_fftplan_create_c2r_backward(
        const underling_problem problem,
        int long_ni,
        underling_real * data,
        unsigned fftw_rigor_flags)
{
    // Sanity check input arguments
    if (SUZERAIN_UNLIKELY(problem == NULL)) {
        SUZERAIN_ERROR_VAL("problem == NULL", SUZERAIN_EINVAL, 0);
    }
    if (SUZERAIN_UNLIKELY(long_ni < 0 || long_ni > 2)) {
        SUZERAIN_ERROR_VAL("long_ni < 0 or long_ni > 2", SUZERAIN_EINVAL, 0);
    }
    const underling_extents extents
        = underling_local_extents(problem, long_ni);
    if (SUZERAIN_UNLIKELY(extents.size[3] % 2)) {
        SUZERAIN_ERROR_NULL(
                "problem must have an even number of underling_real fields",
                SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(extents.order[0] != 3)) {
        SUZERAIN_ERROR_NULL(
                "fields not interleaved: extents.order[0] != 3",
                SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(data == NULL)) {
        SUZERAIN_ERROR_NULL("data == NULL", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(fftw_rigor_flags & non_rigor_mask)) {
        SUZERAIN_ERROR_NULL("FFTW non-rigor bits disallowed", SUZERAIN_EINVAL);
    }

    // Determine the storage ordering necessary for the FFT
    // TODO Fix ESANITY below by reordering for UNDERLING_TRANSPOSED_LONG_N2
    // TODO Fix ESANITY below by reordering for UNDERLING_TRANSPOSED_LONG_N0
    if (SUZERAIN_UNLIKELY(extents.order[1] != long_ni)) {
        SUZERAIN_ERROR_NULL(
                "transformed direction not long: extents.order[1] != long_ni",
                SUZERAIN_ESANITY);
    }

    // Prepare the reordering plan for the input data.  We "rotate" adjacent
    // real and imaginary components so the stride between them is identical.
    fftw_plan plan_reorder = NULL;
    {
        const fftw_iodim howmany_dims[] = {
            { // Loop
                extents.size[extents.order[3]],
                extents.stride[extents.order[3]],
                extents.stride[extents.order[3]]
            },
            { // Loop
                extents.size[extents.order[2]],
                extents.stride[extents.order[2]],
                extents.stride[extents.order[2]]
            },
            { // Loop
                extents.size[extents.order[1]],
                extents.stride[extents.order[1]],
                extents.stride[extents.order[1]]
            },
            { // Transposed
                extents.size[extents.order[0]] / 2,
                2,
                1
            },
            { // Transposed
                2,
                1,
                extents.size[extents.order[0]] / 2
            }
        };
        const int howmany_rank = sizeof(howmany_dims)/sizeof(howmany_dims[0]);

        plan_reorder = fftw_plan_guru_r2r(/*rank*/0, /*dims*/NULL,
                                          howmany_rank, howmany_dims,
                                          data, data,
                                          /*kind*/NULL, fftw_rigor_flags);
    }
    if (SUZERAIN_UNLIKELY(plan_reorder == NULL)) {
        SUZERAIN_ERROR_NULL("FFTW returned a NULL reordering plan",
                SUZERAIN_ESANITY);
    }

    // We transform the long dimension given by extents.order[1] The transform
    // is purely in place which sets our output strides equal to our input
    // strides.
    fftw_plan plan_fft = NULL;
    {
        const fftw_iodim dims[] = {
            {
                2*(extents.size[extents.order[1]] - 1),
                extents.size[extents.order[0]],
                extents.size[extents.order[0]] / 2
            }
        };
        const int rank = sizeof(dims)/sizeof(dims[0]);

        const fftw_iodim howmany_dims[] = {
            {
                extents.size[extents.order[3]],
                extents.stride[extents.order[3]],
                extents.stride[extents.order[3]]
            },
            {
                extents.size[extents.order[2]],
                extents.stride[extents.order[2]],
                extents.stride[extents.order[2]]
            },
            {
                extents.size[extents.order[0]] / 2, // howmany/2
                1,                                  // is, interleaved
                1                                   // os, interleaved
            }
        };
        const int howmany_rank = sizeof(howmany_dims)/sizeof(howmany_dims[0]);

        underling_real * const ri = data;
        underling_real * const ii = data + extents.size[extents.order[0]] / 2;
        underling_real * const out = data;

        plan_fft = fftw_plan_guru_split_dft_c2r(rank, dims,
                                                howmany_rank, howmany_dims,
                                                ri, ii, out,
                                                fftw_rigor_flags);
    }
    if (SUZERAIN_UNLIKELY(plan_fft == NULL)) {
        SUZERAIN_ERROR_NULL("FFTW returned a NULL plan", SUZERAIN_ESANITY);
    }

    // Create and initialize the fftplan workspace
    underling_fftplan f = calloc(1, sizeof(struct underling_fftplan_s));
    if (SUZERAIN_UNLIKELY(f == NULL)) {
        SUZERAIN_ERROR_NULL("failed to allocate space for fftplan",
                             SUZERAIN_ENOMEM);
    }
    // Copy the relevant parameters to the fftplan workspace
    f->plan_reorder = plan_reorder;
    f->plan_fft     = plan_fft;

    return f;

}

int
underling_fftplan_execute(
        const underling_fftplan fftplan)
{
    if (SUZERAIN_UNLIKELY(fftplan == NULL)) {
        SUZERAIN_ERROR("fftplan == NULL", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(fftplan->plan_reorder == NULL)) {
        SUZERAIN_ERROR("fftplan->plan_reorder == NULL", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(fftplan->plan_fft == NULL)) {
        SUZERAIN_ERROR("fftplan->plan_fft == NULL", SUZERAIN_EINVAL);
    }

    fftw_execute(fftplan->plan_reorder);
    fftw_execute(fftplan->plan_fft);

    return SUZERAIN_SUCCESS;
}

void
underling_fftplan_destroy(
        underling_fftplan fftplan)
{
    if (fftplan) {
        if (fftplan->plan_reorder) {
            fftw_destroy_plan(fftplan->plan_reorder);
            fftplan->plan_reorder = NULL;
        }
        if (fftplan->plan_fft) {
            fftw_destroy_plan(fftplan->plan_fft);
            fftplan->plan_fft = NULL;
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
        if (fftplan->plan_reorder) {
            fprintf(output_file, "{plan_reorder:");
            fftw_fprint_plan(fftplan->plan_reorder, output_file);
            fprintf(output_file, "}");
        }
        if (fftplan->plan_fft) {
            fprintf(output_file, "{plan_fft:");
            fftw_fprint_plan(fftplan->plan_fft, output_file);
            fprintf(output_file, "}");
        }
    }
}
