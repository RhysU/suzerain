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

// ********************************************************************
// INTERNAL TYPES INTERNAL TYPES INTERNAL TYPES INTERNAL TYPES INTERNAL
// ********************************************************************

enum transform_type {
    transform_type_unspecified,
    transform_type_c2c_forward,
    transform_type_c2c_backward,
    transform_type_c2r_backward,
    transform_type_r2c_forward
};

struct underling_fft_plan_s {
    int long_ni;                  // Transformed pencil in {0,1,2}
    enum transform_type type;     // Type of transform
    underling_fft_extents input;  // Input data layout
    fftw_plan plan_preorder;      // Executed before the FFT
    fftw_plan plan_fft;           // Performs the FFT
    fftw_plan plan_postorder;     // Executed after the FFT
    underling_fft_extents output; // Output data layout
};

// ********************************************************************
// INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL
// ********************************************************************

fftw_plan
underling_fftw_plan_nop();

fftw_plan
underling_fftw_plan_reorder_complex(
        underling_real * const data,
        const underling_fft_extents * const input,
        const underling_fft_extents * const output);

void
rotate_left(
        int *array,
        int len);

int
adjust_for_fast_stride_in_long_direction(
        underling_extents * const extents,
        const int long_ni);

underling_fft_extents
create_underling_fft_extents_for_complex(
        const underling_extents extents,
        const int long_ni);

underling_fft_extents
create_underling_fft_extents_for_real(
        const underling_extents extents,
        const int long_ni);

underling_fft_plan
underling_fft_plan_create_c2c_internal(
        const int long_ni,
        underling_real * const data,
        const int fftw_sign,
        unsigned fftw_rigor_flags,
        const underling_fft_extents input,
        const underling_fft_extents output);

underling_fft_plan
underling_fft_plan_create_c2r_backward_internal(
        const int long_ni,
        underling_real * data,
        unsigned fftw_rigor_flags,
        const underling_fft_extents input,
        const underling_fft_extents output);

underling_fft_plan
underling_fft_plan_create_r2c_forward_internal(
        const int long_ni,
        underling_real * data,
        unsigned fftw_rigor_flags,
        const underling_fft_extents input,
        const underling_fft_extents output);

void
underling_fft_extents_copy(
        const underling_fft_extents * const e,
        int *start,
        int *size,
        int *stride,
        int *order);

static const unsigned non_rigor_mask =   ~FFTW_ESTIMATE
                                       & ~FFTW_MEASURE
                                       & ~FFTW_PATIENT
                                       & ~FFTW_EXHAUSTIVE
                                       & ~FFTW_WISDOM_ONLY;

// **************************************************************************
// IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION
// **************************************************************************

const underling_fft_extents UNDERLING_FFT_EXTENTS_INVALID = {
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}
};

int
underling_fft_extents_cmp(const underling_fft_extents * const e1,
                          const underling_fft_extents * const e2)
{
    const int start_cmp = memcmp(e1->start, e2->start, sizeof(e1->start));
    if (start_cmp) return start_cmp;

    const int size_cmp = memcmp(e1->size, e2->size, sizeof(e1->size));
    if (size_cmp) return size_cmp;

    const int stride_cmp = memcmp(e1->stride, e2->stride, sizeof(e1->stride));
    if (stride_cmp) return stride_cmp;

    const int order_cmp = memcmp(e1->order, e2->order, sizeof(e1->order));
    if (order_cmp) return order_cmp;

    return 0;
}

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
    if (SUZERAIN_UNLIKELY(nop_plan == NULL)) {
        SUZERAIN_ERROR_NULL("FFTW returned a NULL NOP plan", SUZERAIN_ESANITY);
    }
    return nop_plan;
}

static
fftw_plan
underling_fftw_plan_reorder_complex(
        underling_real * const data,
        const underling_fft_extents * const input,
        const underling_fft_extents * const output)
{
    const int howmany_rank = sizeof(input->size)/sizeof(input->size[0]);
    fftw_iodim howmany_dims[howmany_rank];
    for (int i  = 0; i < howmany_rank; ++i) {
        const int io       = input->order[howmany_rank - 1 - i];
        howmany_dims[i].n  = input->size[io];
        howmany_dims[i].is = input->stride[io];
        howmany_dims[i].os = output->stride[io];
    }
    fftw_plan retval = fftw_plan_guru_r2r(
            0, NULL, howmany_rank, howmany_dims, data, data,
            NULL, FFTW_ESTIMATE);

    if (SUZERAIN_UNLIKELY(retval == NULL)) {
        SUZERAIN_ERROR_NULL("FFTW returned a NULL reorder_complex plan",
                SUZERAIN_ESANITY);
    }

    return retval;
}

static
void
rotate_left(int *array, int len)
{
    if (len) {
        const int tmp = array[0];
        for (int i = 0; i < len - 1; ++i) {
            array[i] = array[i+1];
        }
        array[len - 1] = tmp;
    }
}

static
int
adjust_for_fast_stride_in_long_direction(
        underling_extents * const e,
        int long_ni)
{
    int nrotate = 0;

    while (e->order[1] != long_ni) {
        ++nrotate;
        rotate_left(e->order + 1, 3); // order[0] not touched
    }

    if (nrotate) {
        assert(e->stride[e->order[0]] == 1);
        e->stride[e->order[1]] = e->stride[e->order[0]] * e->size[e->order[0]];
        e->stride[e->order[2]] = e->stride[e->order[1]] * e->size[e->order[1]];
        e->stride[e->order[3]] = e->stride[e->order[2]] * e->size[e->order[2]];
    }

    return nrotate;
}

static
underling_fft_extents
create_underling_fft_extents_for_complex(
        const underling_extents extents,
        const int long_ni)
{
    // Start by copying information from the domain decomposition
    underling_fft_extents retval;
    memcpy(retval.start,  extents.start,  sizeof(extents.start));
    memcpy(retval.size,   extents.size,   sizeof(extents.size));
    memcpy(retval.stride, extents.stride, sizeof(extents.stride));
    memcpy(retval.order,  extents.order,  sizeof(extents.order));

    // Sanity check layout assumptions
    if (SUZERAIN_UNLIKELY(retval.size[3] % 2)) {
        SUZERAIN_ERROR_VAL(
                "problem must have an even number of underling_real fields",
                SUZERAIN_EINVAL,
                UNDERLING_FFT_EXTENTS_INVALID);
    }
    if (SUZERAIN_UNLIKELY(retval.order[0] != 3)) {
        SUZERAIN_ERROR_VAL(
                "transformed fields not interleaved: retval.order[0] != 3",
                SUZERAIN_EINVAL,
                UNDERLING_FFT_EXTENTS_INVALID);
    }
    if (SUZERAIN_UNLIKELY(retval.start[long_ni] != 0)) {
        SUZERAIN_ERROR_VAL(
                "field does not start at zero: retval.start[long_ni] != 0",
                SUZERAIN_EINVAL,
                UNDERLING_FFT_EXTENTS_INVALID);
    }

    // Returned indices 3 and 4 describe interleaved, complex fields
    // built atop extents index 3:
    retval.size[3]   /= 2; // Two adjacent reals make one complex field
    retval.stride[3] *= 2;
    retval.size[4]    = 2; // Each complex value consists of two reals
    retval.stride[4]  = 1; // Real-valued components are adjacent
    retval.start[4]   = 0; // Complex-values are always local
    // Real-valued components are fastest index
    for (int i = 3; i >= 0; --i) {
        retval.order[i+1] = retval.order[i];
    }
    retval.order[0] = 4;

    return retval;
}

static
underling_fft_extents
create_underling_fft_extents_for_real(
        const underling_extents extents,
        const int long_ni)
{
    // Start by copying information from the domain decomposition
    underling_fft_extents retval;
    memcpy(retval.start,  extents.start,  sizeof(extents.start));
    memcpy(retval.size,   extents.size,   sizeof(extents.size));
    memcpy(retval.stride, extents.stride, sizeof(extents.stride));
    memcpy(retval.order,  extents.order,  sizeof(extents.order));

    // Sanity check layout assumptions
    if (SUZERAIN_UNLIKELY(retval.size[3] % 2)) {
        SUZERAIN_ERROR_VAL(
                "problem must have an even number of underling_real fields",
                SUZERAIN_EINVAL,
                UNDERLING_FFT_EXTENTS_INVALID);
    }
    if (SUZERAIN_UNLIKELY(retval.order[0] != 3)) {
        SUZERAIN_ERROR_VAL(
                "transformed fields not interleaved: retval.order[0] != 3",
                SUZERAIN_EINVAL,
                UNDERLING_FFT_EXTENTS_INVALID);
    }
    if (SUZERAIN_UNLIKELY(retval.start[long_ni] != 0)) {
        SUZERAIN_ERROR_VAL(
                "field does not start at zero: retval.start[long_ni] != 0",
                SUZERAIN_EINVAL,
                UNDERLING_FFT_EXTENTS_INVALID);
    }
    if (SUZERAIN_UNLIKELY(retval.start[long_ni] % 2)) {
        SUZERAIN_ERROR_VAL(
                "field does not have even stride",
                SUZERAIN_EINVAL,
                UNDERLING_FFT_EXTENTS_INVALID);
    }

    // The returned layout's index 3 describes real fields built from
    // underling_extents index 3 and the long index.  The long direction does
    // not "cover" the entire underlying memory to allow for real-to-complex
    // padding.  Index 4 reflects that a real scalar value has a single
    // real-valued component.  Definitely a bit goofy, but consistent with
    // complex extents.
    retval.size[long_ni]    = 2*(retval.size[long_ni]-1);
    retval.stride[long_ni] /= 2;
    retval.size[3]         /= 2;
    retval.stride[3]        = 1;
    retval.size[4]          = 1;
    retval.start[4]         = 0;
    retval.stride[4]        = 1;
    // Real-valued components are fastest index
    for (int i = 3; i >= 0; --i) {
        retval.order[i+1] = retval.order[i];
    }
    retval.order[0] = 4;

    return retval;
}

underling_fft_plan
underling_fft_plan_create_c2c_forward(
        const underling_problem problem,
        int long_ni,
        underling_real * data,
        unsigned fftw_rigor_flags)
{
    underling_extents e = underling_local_extents(problem, long_ni);

    const underling_fft_extents input
        = create_underling_fft_extents_for_complex(e, long_ni);

    adjust_for_fast_stride_in_long_direction(&e, long_ni);

    const underling_fft_extents output
        = create_underling_fft_extents_for_complex(e, long_ni);

    return underling_fft_plan_create_c2c_internal(
            long_ni, data, FFTW_FORWARD, fftw_rigor_flags,
            input, output);
}

underling_fft_plan
underling_fft_plan_create_c2c_backward(
        const underling_problem problem,
        int long_ni,
        underling_real * data,
        unsigned fftw_rigor_flags)
{
    underling_extents e = underling_local_extents(problem, long_ni);

    const underling_fft_extents input
        = create_underling_fft_extents_for_complex(e, long_ni);

    adjust_for_fast_stride_in_long_direction(&e, long_ni);

    const underling_fft_extents output
        = create_underling_fft_extents_for_complex(e, long_ni);

    return underling_fft_plan_create_c2c_internal(
            long_ni, data, FFTW_BACKWARD, fftw_rigor_flags,
            input, output);
}

static
underling_fft_plan
underling_fft_plan_create_c2c_internal(
        const int long_ni,
        underling_real * const data,
        const int fftw_sign,
        unsigned fftw_rigor_flags,
        const underling_fft_extents input,
        const underling_fft_extents output)
{
    // Sanity check input arguments
    if (SUZERAIN_UNLIKELY(long_ni < 0 || long_ni > 2)) {
        SUZERAIN_ERROR_VAL("long_ni < 0 or long_ni > 2", SUZERAIN_EINVAL, 0);
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

    // Determine when if/when we reorder relative to the FFT itself
    // Do so by maintaining pointers to the relevant information
    const underling_fft_extents *reorder_in, *reorder_out, *transform;
    const int input_is_long  = (input.order[2]  == long_ni);
    const int output_is_long = (output.order[2] == long_ni);
    if (input_is_long) {
        transform   = &input;  // First, transform the input in-place
        reorder_in  = &input;  // Then reorder the transformed data...
        reorder_out = &output; // to match with the desired output.
    } else if (output_is_long) {
        reorder_in  = &input;  // First, reorder the data ...
        reorder_out = &output; // ...to make it long for the transform.
        transform   = &output; // Then transform it.
    } else {
        SUZERAIN_ERROR_NULL(
                "Neither {input,output}_is_long", SUZERAIN_ESANITY);
    }

    // Cook the reordering plan
    fftw_plan plan_reorder
        = underling_fftw_plan_reorder_complex(data, reorder_in, reorder_out);

    // Prepare the input to fftw_plan_guru_split_dft, which allows using
    // underling_fft_extents.strides directly.  The transform is in place.
    fftw_plan plan_fft = NULL;
    {
        // Transform the long dimension given by transform->order[2]
        assert(transform->order[2] == long_ni);
        const fftw_iodim dims[] = {
            transform->size[transform->order[2]],
            transform->stride[transform->order[2]],
            transform->stride[transform->order[2]]
        };

        // Loop over non-transformed dimensions
        const fftw_iodim howmany_dims[3] = {
            {
                transform->size[transform->order[4]],
                transform->stride[transform->order[4]],
                transform->stride[transform->order[4]]
            },
            {
                transform->size[transform->order[3]],
                transform->stride[transform->order[3]],
                transform->stride[transform->order[3]]
            },
            {
                transform->size[3],
                transform->size[4],
                transform->size[4]
            }
        };

        // FFTW manual section 4.5.3, FFTW_BACKWARD is FFTW_FORWARD with the
        // components flipped.
        underling_real * const ri
            = (fftw_sign == FFTW_FORWARD) ? data : data + transform->stride[4];
        underling_real * const ii
            = (fftw_sign == FFTW_FORWARD) ? data + transform->stride[4] : data;
        underling_real * const ro = ri;
        underling_real * const io = ii;

        plan_fft = fftw_plan_guru_split_dft(
                sizeof(dims)/sizeof(dims[0]), dims,
                sizeof(howmany_dims)/sizeof(howmany_dims[0]), howmany_dims,
                ri, ii, ro, io, fftw_rigor_flags);

        if (SUZERAIN_UNLIKELY(plan_fft == NULL)) {
            SUZERAIN_ERROR_NULL("FFTW returned a NULL FFT plan",
                    SUZERAIN_ESANITY);
        }
    }

    // Create and initialize the plan workspace
    underling_fft_plan f = calloc(1, sizeof(struct underling_fft_plan_s));
    if (SUZERAIN_UNLIKELY(f == NULL)) {
        SUZERAIN_ERROR_NULL("failed to allocate space for plan",
                             SUZERAIN_ENOMEM);
    }

    // Copy the relevant parameters to the plan workspace...
    f->long_ni        = long_ni;
    f->type           = (fftw_sign == FFTW_FORWARD)
                      ? transform_type_c2c_forward
                      : transform_type_c2c_backward;
    f->input          = input;
    f->plan_fft       = plan_fft;
    f->output         = output;
    // ...and fix when the data is reordered relative to the FFT
    if (input_is_long) {
        f->plan_preorder  = underling_fftw_plan_nop();
        f->plan_postorder = plan_reorder;
    } else if (output_is_long) {
        f->plan_preorder  = plan_reorder;
        f->plan_postorder = underling_fftw_plan_nop();
    } else {
        SUZERAIN_ERROR_NULL(
                "Neither {input,output}_is_long", SUZERAIN_ESANITY);
    }

    return f;
}

underling_fft_plan
underling_fft_plan_create_c2r_backward(
        const underling_problem problem,
        int long_ni,
        underling_real * data,
        unsigned fftw_rigor_flags)
{
    underling_extents e = underling_local_extents(problem, long_ni);

    const underling_fft_extents input
        = create_underling_fft_extents_for_complex(e, long_ni);

    adjust_for_fast_stride_in_long_direction(&e, long_ni);

    const underling_fft_extents output
        = create_underling_fft_extents_for_real(e, long_ni);

    return underling_fft_plan_create_c2r_backward_internal(
            long_ni, data, fftw_rigor_flags, input, output);
}

static
underling_fft_plan
underling_fft_plan_create_c2r_backward_internal(
        const int long_ni,
        underling_real * data,
        unsigned fftw_rigor_flags,
        const underling_fft_extents input,
        const underling_fft_extents output)
{
    // Sanity check input arguments
    if (SUZERAIN_UNLIKELY(long_ni < 0 || long_ni > 2)) {
        SUZERAIN_ERROR_VAL("long_ni < 0 or long_ni > 2", SUZERAIN_EINVAL, 0);
    }
    if (SUZERAIN_UNLIKELY(data == NULL)) {
        SUZERAIN_ERROR_NULL("data == NULL", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(fftw_rigor_flags & non_rigor_mask)) {
        SUZERAIN_ERROR_NULL("FFTW non-rigor bits disallowed", SUZERAIN_EINVAL);
    }

    // Prepare the pre-ordering plan.
    fftw_plan plan_preorder = NULL;
    {
        const int howmany_rank = sizeof(input.size)/sizeof(input.size[0]);
        fftw_iodim howmany_dims[howmany_rank];
        // Copy dimensions n{0,1,2}
        for (int i = 0; i < 3; ++i) {
            const int io = input.order[howmany_rank - 1 - i];
            assert(io < 3);
            howmany_dims[i].n  = input.size[io];
            howmany_dims[i].is = input.stride[io];
            howmany_dims[i].os = input.stride[io];
        }
        // Transpose fields to have fixed strides between
        // real and imaginary components.
        howmany_dims[3].n  = input.size[3];
        howmany_dims[3].is = 2;
        howmany_dims[3].os = 1;
        howmany_dims[4].n  = 2;
        howmany_dims[4].is = 1;
        howmany_dims[4].os = input.size[3];

        plan_preorder = fftw_plan_guru_r2r(/*rank*/0, /*dims*/NULL,
                                           howmany_rank, howmany_dims,
                                           data, data,
                                           /*kind*/NULL, fftw_rigor_flags);
        if (SUZERAIN_UNLIKELY(plan_preorder == NULL)) {
            SUZERAIN_ERROR_NULL("FFTW returned a NULL preorder plan",
                    SUZERAIN_ESANITY);
        }
    }

    // Plan in-place transform
    fftw_plan plan_fft = NULL;
    {
        const fftw_iodim dims[] = {       // Transform long_ni
            {
                output.size[long_ni],     // Logical transform size
                input.stride[long_ni],
                input.size[3]             // From preorder plan
            }
        };
        const int rank = sizeof(dims)/sizeof(dims[0]);

        const int howmany_rank = 3;       // Loop over other directions
        fftw_iodim howmany_dims[howmany_rank];
        int j = 0;
        for (const int *oo = output.order+4; oo >= output.order+2; --oo) {
            if (*oo == long_ni) continue; // Skip transformed direction
            assert(*oo < 3);
            assert(j < 2);
            howmany_dims[j].n  = input.size[*oo];
            howmany_dims[j].is = input.stride[*oo];
            howmany_dims[j].os = input.stride[*oo];
            ++j;
        }
        assert(j == 2);
        howmany_dims[2].n  = input.size[3];
        howmany_dims[2].is = 1;
        howmany_dims[2].os = 1;

        underling_real * const ri = data;
        underling_real * const ii = data + input.size[3]; // From preorder plan
        underling_real * const out = data;

        plan_fft = fftw_plan_guru_split_dft_c2r(rank, dims,
                                                howmany_rank, howmany_dims,
                                                ri, ii, out,
                                                fftw_rigor_flags);
        if (SUZERAIN_UNLIKELY(plan_fft == NULL)) {
            SUZERAIN_ERROR_NULL("FFTW returned a NULL FFT plan",
                    SUZERAIN_ESANITY);
        }
    }

    // Prepare output reordering plan
    fftw_plan plan_postorder = NULL;
    {
        const int howmany_rank = sizeof(input.size)/sizeof(input.size[0]);
        fftw_iodim howmany_dims[howmany_rank];
        int j = 0;
        for (const int *oo = output.order+4; oo >= output.order+2; --oo) {
            assert(*oo < 3);
            if (*oo == long_ni) {
                howmany_dims[j].n  = 2*input.size[*oo];   // 2*(n/2+1) pad
                howmany_dims[j].is = input.stride[*oo]/2;
                howmany_dims[j].os = output.stride[*oo];
            } else {
                howmany_dims[j].n  = input.size[*oo];
                howmany_dims[j].is = input.stride[*oo];
                howmany_dims[j].os = output.stride[*oo];
            }
            ++j;
        }
        assert(input.size[3] == output.size[3]);
        howmany_dims[3].n  = output.size[3];
        howmany_dims[3].is = howmany_dims[3].os = 1;
        howmany_dims[4].n  = 1;
        howmany_dims[4].is = howmany_dims[4].os = output.size[3];

        plan_postorder = fftw_plan_guru_r2r(/*rank*/0, /*dims*/NULL,
                                           howmany_rank, howmany_dims,
                                           data, data,
                                           /*kind*/NULL, fftw_rigor_flags);
        if (SUZERAIN_UNLIKELY(plan_postorder == NULL)) {
            SUZERAIN_ERROR_NULL("FFTW returned a NULL postorder plan",
                    SUZERAIN_ESANITY);
        }
    }

    // Create and initialize the plan workspace
    underling_fft_plan f = calloc(1, sizeof(struct underling_fft_plan_s));
    if (SUZERAIN_UNLIKELY(f == NULL)) {
        SUZERAIN_ERROR_NULL("failed to allocate space for plan",
                             SUZERAIN_ENOMEM);
    }
    // Copy the relevant parameters to the plan workspace
    f->long_ni        = long_ni;
    f->type           = transform_type_c2r_backward;
    f->input          = input;
    f->plan_preorder  = plan_preorder;
    f->plan_fft       = plan_fft;
    f->plan_postorder = plan_postorder;
    f->output         = output;

    return f;

}

underling_fft_plan
underling_fft_plan_create_r2c_forward(
        const underling_problem problem,
        int long_ni,
        underling_real * data,
        unsigned fftw_rigor_flags)
{
    underling_extents e = underling_local_extents(problem, long_ni);

    const underling_fft_extents input
        = create_underling_fft_extents_for_real(e, long_ni);

    adjust_for_fast_stride_in_long_direction(&e, long_ni);

    const underling_fft_extents output
        = create_underling_fft_extents_for_complex(e, long_ni);

    return underling_fft_plan_create_r2c_forward_internal(
            long_ni, data, fftw_rigor_flags, input, output);
}

static
underling_fft_plan
underling_fft_plan_create_r2c_forward_internal(
        const int long_ni,
        underling_real * data,
        unsigned fftw_rigor_flags,
        const underling_fft_extents input,
        const underling_fft_extents output)
{
    // Sanity check input arguments
    if (SUZERAIN_UNLIKELY(long_ni < 0 || long_ni > 2)) {
        SUZERAIN_ERROR_VAL("long_ni < 0 or long_ni > 2", SUZERAIN_EINVAL, 0);
    }
    if (SUZERAIN_UNLIKELY(data == NULL)) {
        SUZERAIN_ERROR_NULL("data == NULL", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(fftw_rigor_flags & non_rigor_mask)) {
        SUZERAIN_ERROR_NULL("FFTW non-rigor bits disallowed", SUZERAIN_EINVAL);
    }

    // Prepare the reordering plan for the input data. We must always pay to
    // reorder the output data.  Ignoring the input reordering process avoids
    // touching everything a third separate time.
    // TODO Evaluate performance impact since FFT may be non-stride 1
    const fftw_plan plan_preorder = underling_fftw_plan_nop();

    // The transform is purely in place.
    fftw_plan plan_fft = NULL;
    {
        const fftw_iodim dims[] = {       // Transform long_ni
            {
                input.size[long_ni],
                input.stride[long_ni],
                2*input.stride[long_ni]
            }
        };
        const int rank = sizeof(dims)/sizeof(dims[0]);

        const int howmany_rank = 3;       // Loop over other directions
        fftw_iodim howmany_dims[howmany_rank];
        int j = 0;
        for (const int *io = input.order+4; io >= input.order+1; --io) {
            if (*io == long_ni) continue; // Skip transformed direction
            assert(j < sizeof(howmany_dims)/sizeof(howmany_dims[0]));
            howmany_dims[j].n  = input.size[*io];
            howmany_dims[j].is = input.stride[*io];
            howmany_dims[j].os = input.stride[*io];
            ++j;
        }
        assert(j == sizeof(howmany_dims)/sizeof(howmany_dims[0]));

        underling_real * const in = data;
        underling_real * const ro = data;
        underling_real * const io = data + input.stride[long_ni];

        plan_fft = fftw_plan_guru_split_dft_r2c(rank, dims,
                                                howmany_rank, howmany_dims,
                                                in, ro, io,
                                                fftw_rigor_flags);
        if (SUZERAIN_UNLIKELY(plan_fft == NULL)) {
            SUZERAIN_ERROR_NULL("FFTW returned a NULL FFT plan",
                    SUZERAIN_ESANITY);
        }
    }

    // Prepare the reordering plan for the output data.
    // Use "complex-centric" sizes since they cover the contiguous region.
    fftw_plan plan_postorder = NULL;
    {
        const int howmany_rank = sizeof(output.size)/sizeof(output.size[0]);
        fftw_iodim howmany_dims[howmany_rank];
        int j = 0;
        for (const int *io = input.order+4; io >= input.order+2; --io) {
            assert(j < sizeof(howmany_dims)/sizeof(howmany_dims[0]));
            howmany_dims[j].n  = output.size[*io];
            howmany_dims[j].is = input.stride[*io];
            howmany_dims[j].os = output.stride[*io];
            if (*io == long_ni) {
                howmany_dims[j].is *= 2; // Transform modifies input stride
            }
            ++j;
        }
        assert(j == sizeof(howmany_dims)/sizeof(howmany_dims[0]) - 2);
        howmany_dims[3].n  = output.size[3];
        howmany_dims[3].is = 1;
        howmany_dims[3].os = 2;
        howmany_dims[4].n  = 2;
        howmany_dims[4].is = output.size[3];
        howmany_dims[4].os = 1;

        plan_postorder = fftw_plan_guru_r2r(/*rank*/0, /*dims*/NULL,
                                            howmany_rank, howmany_dims,
                                            data, data,
                                            /*kind*/NULL, fftw_rigor_flags);
        if (SUZERAIN_UNLIKELY(plan_postorder == NULL)) {
            SUZERAIN_ERROR_NULL("FFTW returned a NULL postorder plan",
                    SUZERAIN_ESANITY);
        }
    }

    // Create and initialize the plan workspace
    underling_fft_plan f = calloc(1, sizeof(struct underling_fft_plan_s));
    if (SUZERAIN_UNLIKELY(f == NULL)) {
        SUZERAIN_ERROR_NULL("failed to allocate space for plan",
                             SUZERAIN_ENOMEM);
    }
    // Copy the relevant parameters to the plan workspace
    f->long_ni        = long_ni;
    f->type           = transform_type_r2c_forward;
    f->input          = input;
    f->plan_preorder  = plan_preorder;
    f->plan_fft       = plan_fft;
    f->plan_postorder = plan_postorder;
    f->output         = output;

    return f;
}

underling_fft_plan
underling_fft_plan_create_inverse(
        const underling_fft_plan plan_to_invert,
        underling_real * data,
        unsigned fftw_rigor_flags)
{
    if (SUZERAIN_UNLIKELY(plan_to_invert == NULL)) {
        SUZERAIN_ERROR_NULL("plan_to_invert == NULL", SUZERAIN_EINVAL);
    }

    underling_fft_plan retval = NULL;

    switch (plan_to_invert->type) {
    case transform_type_c2c_forward:
        retval = underling_fft_plan_create_c2c_internal(
                    plan_to_invert->long_ni, data, FFTW_BACKWARD,
                    fftw_rigor_flags,
                    plan_to_invert->output, plan_to_invert->input);
        break;
    case transform_type_c2c_backward:
        retval = underling_fft_plan_create_c2c_internal(
                    plan_to_invert->long_ni, data, FFTW_FORWARD,
                    fftw_rigor_flags,
                    plan_to_invert->output, plan_to_invert->input);
        break;
    case transform_type_c2r_backward:
        retval = underling_fft_plan_create_r2c_forward_internal(
                    plan_to_invert->long_ni, data,
                    fftw_rigor_flags,
                    plan_to_invert->output, plan_to_invert->input);
        break;
    case transform_type_r2c_forward:
        retval = underling_fft_plan_create_c2r_backward_internal(
                    plan_to_invert->long_ni, data,
                    fftw_rigor_flags,
                    plan_to_invert->output, plan_to_invert->input);
        break;
    case transform_type_unspecified:
    default:
        SUZERAIN_ERROR_NULL("Unrecognized plan_to_invert->type",
                SUZERAIN_ESANITY);
    }

    return retval;
}

underling_fft_extents
underling_fft_local_extents_input(
        const underling_fft_plan plan)
{
    if (SUZERAIN_UNLIKELY(plan == NULL)) {
        SUZERAIN_ERROR_VAL("plan == NULL",
                SUZERAIN_EINVAL, UNDERLING_FFT_EXTENTS_INVALID);
    }

    underling_fft_extents retval = plan->input; // Create temporary
    return retval;                              // Return temporary
}

underling_fft_extents
underling_fft_local_extents_output(
        const underling_fft_plan plan)
{
    if (SUZERAIN_UNLIKELY(plan == NULL)) {
        SUZERAIN_ERROR_VAL("plan == NULL",
                SUZERAIN_EINVAL, UNDERLING_FFT_EXTENTS_INVALID);
    }

    underling_fft_extents retval = plan->output; // Create temporary
    return retval;                               // Return temporary
}

static
void
underling_fft_extents_copy(
        const underling_fft_extents * const e,
        int *start,
        int *size,
        int *stride,
        int *order)
{
    if (start) {
        for (int j = 0; j < sizeof(e->start)/sizeof(e->start[0]); ++j)
            start[j] = e->start[j];
    }
    if (size) {
        for (int j = 0; j < sizeof(e->size)/sizeof(e->size[0]); ++j)
            size[j] = e->size[j];
    }
    if (stride) {
        for (int j = 0; j < sizeof(e->stride)/sizeof(e->stride[0]); ++j)
            stride[j] = e->stride[j];
    }
    if (order) {
        for (int j = 0; j < sizeof(e->order)/sizeof(e->order[0]); ++j)
            order[j] = e->order[j];
    }
}

void
underling_fft_local_input(
        const underling_fft_plan plan,
        int *start,
        int *size,
        int *stride,
        int *order)
{
    if (SUZERAIN_UNLIKELY(plan == NULL)) {
        SUZERAIN_ERROR_VOID("plan == NULL", SUZERAIN_EINVAL);
    }

    underling_fft_extents_copy(&plan->input, start, size, stride, order);
}

void
underling_fft_local_output(
        const underling_fft_plan plan,
        int *start,
        int *size,
        int *stride,
        int *order)
{
    if (SUZERAIN_UNLIKELY(plan == NULL)) {
        SUZERAIN_ERROR_VOID("plan == NULL", SUZERAIN_EINVAL);
    }

    underling_fft_extents_copy(&plan->output, start, size, stride, order);
}

int
underling_fft_plan_execute(
        const underling_fft_plan plan)
{
    if (SUZERAIN_UNLIKELY(plan == NULL)) {
        SUZERAIN_ERROR("plan == NULL", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(plan->plan_preorder == NULL)) {
        SUZERAIN_ERROR("plan->plan_preorder == NULL", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(plan->plan_fft == NULL)) {
        SUZERAIN_ERROR("plan->plan_fft == NULL", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(plan->plan_postorder == NULL)) {
        SUZERAIN_ERROR("plan->plan_postorder == NULL", SUZERAIN_EINVAL);
    }

    fftw_execute(plan->plan_preorder);
    fftw_execute(plan->plan_fft);
    fftw_execute(plan->plan_postorder);

    return SUZERAIN_SUCCESS;
}

void
underling_fft_plan_destroy(
        underling_fft_plan plan)
{
    if (plan) {
        if (plan->plan_preorder) {
            fftw_destroy_plan(plan->plan_preorder);
            plan->plan_preorder = NULL;
        }
        if (plan->plan_fft) {
            fftw_destroy_plan(plan->plan_fft);
            plan->plan_fft = NULL;
        }
        if (plan->plan_postorder) {
            fftw_destroy_plan(plan->plan_postorder);
            plan->plan_postorder = NULL;
        }
        free(plan);
    }
}

void
underling_fft_fprint_extents(
        const underling_fft_extents *extents,
        FILE *output_file)
{
    if (!extents) {
        fprintf(output_file, "NULL");
    } else {
        fprintf(output_file, "extents=[%d,%d)x[%d,%d)x[%d,%d)x[%d,%d)x[%d,%d)",
                extents->start[0],
                extents->start[0] + extents->size[0],
                extents->start[1],
                extents->start[1] + extents->size[1],
                extents->start[2],
                extents->start[2] + extents->size[2],
                extents->start[3],
                extents->start[3] + extents->size[3],
                extents->start[4],
                extents->start[4] + extents->size[4]);
        fprintf(output_file, ",strides={%d,%d,%d,%d,%d}",
                extents->stride[0],
                extents->stride[1],
                extents->stride[2],
                extents->stride[3],
                extents->stride[4]);
    }
}

void
underling_fprint_fft_plan(
        const underling_fft_plan plan,
        FILE *output_file)
{
    fprintf(output_file, "{underling_fft_plan:");
    if (!plan) {
        fprintf(output_file, "NULL");
    } else {
        if (plan->plan_preorder) {
            fprintf(output_file, "{plan_preorder:");
            fftw_fprint_plan(plan->plan_preorder, output_file);
            fprintf(output_file, "}");
        }
        if (plan->plan_fft) {
            fprintf(output_file, "{plan_fft:");
            fftw_fprint_plan(plan->plan_fft, output_file);
            fprintf(output_file, "}");
        }
        if (plan->plan_postorder) {
            fprintf(output_file, "{plan_postorder:");
            fftw_fprint_plan(plan->plan_postorder, output_file);
            fprintf(output_file, "}");
        }
    }
}
