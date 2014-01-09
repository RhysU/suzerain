//-----------------------------------------------------------------------bl-
// underling 0.3.1: an FFTW MPI-based library for 3D pencil decompositions
// http://red.ices.utexas.edu/projects/underling
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
//
// This file is part of underling.
//
// underling is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// underling is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with underling.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------el-
// $Id$

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <underling/underling_fftw.h>

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#include <underling/error.h>
#include "common.h"

// TODO Check for memory leaks stemming from planning failures

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

struct underling_fftw_plan_s {
    int long_ni;                   // Transformed pencil in {0,1,2}
    enum transform_type type;      // Type of transform
    underling_fftw_extents input;  // Input data layout
    fftw_plan plan_preorder;       // Executed before the FFT
    fftw_plan plan_fftw;           // Performs the FFT
    fftw_plan plan_postorder;      // Executed after the FFT
    underling_fftw_extents output; // Output data layout
    _Bool in_place;                // Was planning done for in-place transform?
    struct {                       // FFTW new execute interface offsets...
        int ri;                    //   ...real part relative to in buffer
        int ii;                    //   ...imag part relative to in buffer
        int ro;                    //   ...real part relative to out buffer
        int io;                    //   ...imag part relative to out buffer
    } offset;                      // ...not all of which are always used
};

// ********************************************************************
// INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL
// ********************************************************************

static inline void swap(int *a, int *b)
{
    int t = *a; *a = *b; *b = t;
}

static
fftw_plan
underling_fftw_plan_reorder_complex(
        underling_real * const data_in,
        underling_real * const data_out,
        const underling_fftw_extents * const input,
        const underling_fftw_extents * const output,
        const unsigned flags);

static
void
rotate_left(
        int *array,
        int len);

static
int
adjust_for_fast_stride_in_long_direction(
        underling_extents * const extents,
        const int long_ni);

static
void
pack_strides_according_to_order(
        underling_fftw_extents * const e);

static
underling_fftw_extents
create_underling_fftw_extents_for_complex(
        const underling_extents e,
        const int long_ni);

static
underling_fftw_extents
create_underling_fftw_extents_for_real(
        const underling_extents e,
        const int long_ni);

static
underling_fftw_plan
underling_fftw_plan_create_c2c_internal(
        const int long_ni,
        underling_real * const in,
        underling_real * const out,
        const int fftw_sign,
        unsigned fftw_rigor_flags,
        const underling_fftw_extents input,
        const underling_fftw_extents output);

static
underling_fftw_plan
underling_fftw_plan_create_c2r_backward_internal(
        const int long_ni,
        underling_real * const in,
        underling_real * const out,
        unsigned fftw_rigor_flags,
        const underling_fftw_extents input,
        const underling_fftw_extents output);

static
underling_fftw_plan
underling_fftw_plan_create_r2c_forward_internal(
        const int long_ni,
        underling_real * const in,
        underling_real * const out,
        unsigned fftw_rigor_flags,
        const underling_fftw_extents input,
        const underling_fftw_extents output);

static
void
underling_fftw_extents_copy(
        const underling_fftw_extents * const e,
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

const underling_fftw_extents UNDERLING_FFTW_EXTENTS_INVALID = {
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}
};

int
underling_fftw_extents_cmp(const underling_fftw_extents * const e1,
                          const underling_fftw_extents * const e2)
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
underling_fftw_plan_reorder_complex(
        underling_real * const data_in,
        underling_real * const data_out,
        const underling_fftw_extents * const input,
        const underling_fftw_extents * const output,
        const unsigned flags)
{
    const int howmany_rank = sizeof(input->size)/sizeof(input->size[0]);
    fftw_iodim howmany_dims[howmany_rank];
    for (int i = 0; i < howmany_rank; ++i) {
        const int io       = input->order[howmany_rank - 1 - i];
        howmany_dims[i].n  = input->size[io];
        howmany_dims[i].is = input->stride[io];
        howmany_dims[i].os = output->stride[io];
    }
    fftw_plan retval = fftw_plan_guru_r2r(0, NULL, howmany_rank, howmany_dims,
                                          data_in, data_out, NULL, flags);

    if (UNDERLING_UNLIKELY(retval == NULL)) {
        UNDERLING_ERROR_NULL("FFTW returned NULL reorder_complex plan",
                UNDERLING_ESANITY);
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
void
pack_strides_according_to_order(underling_fftw_extents * const e)
{
    // Modify strides to produce packed, contiguous data
    e->stride[e->order[0]] = 1;
    e->stride[e->order[1]] = e->stride[e->order[0]] * e->size[e->order[0]];
    e->stride[e->order[2]] = e->stride[e->order[1]] * e->size[e->order[1]];
    e->stride[e->order[3]] = e->stride[e->order[2]] * e->size[e->order[2]];
    e->stride[e->order[4]] = e->stride[e->order[3]] * e->size[e->order[3]];
}

static
underling_fftw_extents
create_underling_fftw_extents_for_complex(
        const underling_extents e,
        const int long_ni)
{
    // Sanity check layout assumptions
    if (UNDERLING_UNLIKELY(e.size[3] % 2)) {
        UNDERLING_ERROR_VAL(
                "problem must have an even number of underling_real fields",
                UNDERLING_EINVAL,
                UNDERLING_FFTW_EXTENTS_INVALID);
    }
    if (UNDERLING_UNLIKELY(e.order[0] != 3)) {
        UNDERLING_ERROR_VAL("Unknown interleaving: e.order[0] != 3",
                            UNDERLING_ESANITY,
                            UNDERLING_FFTW_EXTENTS_INVALID);
    }
    if (UNDERLING_UNLIKELY(e.start[long_ni] != 0)) {
        UNDERLING_ERROR_VAL(
                "field does not start at zero: e.start[long_ni] != 0",
                UNDERLING_EINVAL,
                UNDERLING_FFTW_EXTENTS_INVALID);
    }

    underling_fftw_extents r;

    // Start by copying information from the domain decomposition
    memcpy(r.start,  e.start,  sizeof(e.start));
    memcpy(r.size,   e.size,   sizeof(e.size));
    memcpy(r.stride, e.stride, sizeof(e.stride));
    memcpy(r.order,  e.order,  sizeof(e.order));

    // Returned indices 3 and 4 describe interleaved, complex fields
    // built atop extents index 3:
    r.size[3]   /= 2;            // Two reals make one complex field
    r.size[4]    = 2;            // Each complex consists of two reals
    r.start[4]   = 0;            // Complex-values are always local
    r.stride[4]  = r.stride[3];  // Real components are adjacent...
    r.stride[3] *= 2;            // while complex values spread out

    // Real-valued components are fastest index
    for (int i = 3; i >= 0; --i) {
        r.order[i+1] = r.order[i];
    }
    r.order[0] = 4;

    return r;
}

static
underling_fftw_extents
create_underling_fftw_extents_for_real(
        const underling_extents e,
        const int long_ni)
{

    // Sanity check layout assumptions
    if (UNDERLING_UNLIKELY(e.size[3] % 2)) {
        UNDERLING_ERROR_VAL(
                "problem must have an even number of underling_real fields",
                UNDERLING_EINVAL,
                UNDERLING_FFTW_EXTENTS_INVALID);
    }
    if (UNDERLING_UNLIKELY(e.order[0] != 3)) {
        UNDERLING_ERROR_VAL("Unknown interleaving: e.order[0] != 3",
                            UNDERLING_ESANITY,
                            UNDERLING_FFTW_EXTENTS_INVALID);
    }
    if (UNDERLING_UNLIKELY(e.start[long_ni] != 0)) {
        UNDERLING_ERROR_VAL(
                "field does not start at zero: e.start[long_ni] != 0",
                UNDERLING_EINVAL,
                UNDERLING_FFTW_EXTENTS_INVALID);
    }
    if (UNDERLING_UNLIKELY(e.start[long_ni] % 2)) {
        UNDERLING_ERROR_VAL(
                "field does not have even stride",
                UNDERLING_EINVAL,
                UNDERLING_FFTW_EXTENTS_INVALID);
    }

    underling_fftw_extents r;

    if (e.order[1] == long_ni) {

        // Start by copying information from the domain decomposition
        memcpy(r.start,  e.start,  sizeof(e.start));
        memcpy(r.size,   e.size,   sizeof(e.size));
        memcpy(r.stride, e.stride, sizeof(e.stride));
        memcpy(r.order,  e.order,  sizeof(e.order));

        // The returned layout's index 3 describes real fields built from
        // underling_extents index 3 and the long index.  The long direction
        // does not "cover" the entire underlying memory to allow for
        // real-to-complex padding.  Index 4 reflects that a real scalar value
        // has a single real-valued component.  Definitely a bit goofy, but
        // consistent with complex extents.
        r.size[long_ni]    = r.size[long_ni] == 1
                             ? /* degenerate complex size == real size == */ 1
                             : 2*(r.size[long_ni]-1);
        r.stride[long_ni] /= 2;
        r.size[3]         /= 2;
        r.stride[3]        = 1;
        r.size[4]          = 1;
        r.start[4]         = 0;
        r.stride[4]        = 1;

        // Real-valued components are fastest index
        for (int i = 3; i >= 0; --i) {
            r.order[i+1] = r.order[i];
        }
        r.order[0] = 4;

    } else {
        UNDERLING_ERROR_VAL("Unknown interleaving: e.order[1] != long_ni",
                            UNDERLING_ESANITY,
                            UNDERLING_FFTW_EXTENTS_INVALID);
    }

    return r;
}

underling_fftw_plan
underling_fftw_plan_create_c2c_forward(
        const underling_problem problem,
        int long_ni,
        underling_real * in,
        underling_real * out,
        unsigned fftw_rigor_flags,
        unsigned packed_flags)
{
    if (!packed_flags) packed_flags = UNDERLING_FFTW_PACKED_NONE; // Default

    underling_extents e = underling_local_extents(problem, long_ni);

    const underling_fftw_extents input
        = create_underling_fftw_extents_for_complex(e, long_ni);

    adjust_for_fast_stride_in_long_direction(&e, long_ni);

    underling_fftw_extents output
        = create_underling_fftw_extents_for_complex(e, long_ni);

    if (    (long_ni == 2 && packed_flags & UNDERLING_FFTW_PACKED_LONG_N2)
         || (long_ni == 0 && packed_flags & UNDERLING_FFTW_PACKED_LONG_N0)) {
        if (UNDERLING_UNLIKELY(in == out)) {
            UNDERLING_ERROR_NULL("invalid packed_flags for in-place transform",
                                 UNDERLING_EINVAL);
        } else {
            pack_strides_according_to_order(&output);
        }
    }

    return underling_fftw_plan_create_c2c_internal(
            long_ni, in, out, FFTW_FORWARD, fftw_rigor_flags,
            input, output);
}

underling_fftw_plan
underling_fftw_plan_create_c2c_backward(
        const underling_problem problem,
        int long_ni,
        underling_real * in,
        underling_real * out,
        unsigned fftw_rigor_flags,
        unsigned packed_flags)
{
    if (!packed_flags) packed_flags = UNDERLING_FFTW_PACKED_NONE; // Default

    underling_extents e = underling_local_extents(problem, long_ni);

    const underling_fftw_extents input
        = create_underling_fftw_extents_for_complex(e, long_ni);

    adjust_for_fast_stride_in_long_direction(&e, long_ni);

    underling_fftw_extents output
        = create_underling_fftw_extents_for_complex(e, long_ni);

    if (    (long_ni == 2 && packed_flags & UNDERLING_FFTW_PACKED_LONG_N2)
         || (long_ni == 0 && packed_flags & UNDERLING_FFTW_PACKED_LONG_N0)) {
        if (UNDERLING_UNLIKELY(in == out)) {
            UNDERLING_ERROR_NULL("invalid packed_flags for in-place transform",
                                 UNDERLING_EINVAL);
        } else {
            pack_strides_according_to_order(&output);
        }
    }

    return underling_fftw_plan_create_c2c_internal(
            long_ni, in, out, FFTW_BACKWARD, fftw_rigor_flags,
            input, output);
}

static
underling_fftw_plan
underling_fftw_plan_create_c2c_internal(
        const int long_ni,
        underling_real * const in,
        underling_real * const out,
        const int fftw_sign,
        unsigned fftw_rigor_flags,
        const underling_fftw_extents input,
        const underling_fftw_extents output)
{
    // Sanity check input arguments
    if (UNDERLING_UNLIKELY(long_ni < 0 || long_ni > 2)) {
        UNDERLING_ERROR_NULL("long_ni < 0 or long_ni > 2", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(in == NULL)) {
        UNDERLING_ERROR_NULL("in == NULL", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(out == NULL)) {
        UNDERLING_ERROR_NULL("out == NULL", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(   fftw_sign != FFTW_FORWARD
                           && fftw_sign != FFTW_BACKWARD)) {
        UNDERLING_ERROR_NULL("fftw_sign not one of FFTW_{FORWARD,BACKWARD}",
                             UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(fftw_rigor_flags & non_rigor_mask)) {
        UNDERLING_ERROR_NULL("FFTW non-rigor bits disallowed",
                             UNDERLING_EINVAL);
    }

    // Check user requested acceptable transform configuration
    const int input_is_long  = (input.order[2]  == long_ni);
    const int output_is_long = (output.order[2] == long_ni);
    if (UNDERLING_UNLIKELY(!input_is_long && !output_is_long)) {
        UNDERLING_ERROR_NULL(
                "Neither {input,output} is long", UNDERLING_EINVAL);
    }

    // Allocate and initialize the plan workspace
    underling_fftw_plan f = malloc(sizeof(struct underling_fftw_plan_s));
    if (UNDERLING_UNLIKELY(f == NULL)) {
        UNDERLING_ERROR_NULL("failed to allocate space for plan",
                             UNDERLING_ENOMEM);
    }
    f->long_ni        = long_ni;
    f->type           = (fftw_sign == FFTW_FORWARD)
                      ? transform_type_c2c_forward
                      : transform_type_c2c_backward;
    f->input          = input;
    f->plan_preorder  = NULL;
    f->plan_fftw      = NULL;
    f->plan_postorder = NULL;
    f->output         = output;
    f->in_place       = (in == out);

    // Determine when if/when we reorder relative to the FFT operation
    const underling_fftw_extents *transform_in, *transform_out;
    if (f->in_place) {
        // In-place requires reordering either before or after the FFT
        if (input_is_long) {
            transform_in = transform_out = &input;
            f->plan_postorder = underling_fftw_plan_reorder_complex(
                    in, out, &input, &output,
                    fftw_rigor_flags | FFTW_DESTROY_INPUT);
            if (UNDERLING_UNLIKELY(!f->plan_postorder)) {
                underling_fftw_plan_destroy(f);
                UNDERLING_ERROR_NULL("!f->plan_postorder", UNDERLING_ESANITY);
            }
        } else if (output_is_long) {
            f->plan_preorder = underling_fftw_plan_reorder_complex(
                    in, out, &input, &output,
                    fftw_rigor_flags | FFTW_DESTROY_INPUT);
            if (UNDERLING_UNLIKELY(!f->plan_preorder)) {
                underling_fftw_plan_destroy(f);
                UNDERLING_ERROR_NULL("!f->plan_preorder", UNDERLING_ESANITY);
            }
            transform_in = transform_out = &output;
        } else {
            UNDERLING_ERROR_NULL("!{input,output}_is_long", UNDERLING_ESANITY);
        }
    } else {
        // Out-of-place does not require additional reordering steps
        transform_in  = &f->input;
        transform_out = &f->output;
    }

    // Prepare the input to fftw_plan_guru_split_dft, which allows using
    // underling_fftw_extents.strides directly from transform_in, transform_out
    {
        const fftw_iodim dims[] = {       // Transform long_ni
            {
                transform_in ->size  [long_ni],
                transform_in ->stride[long_ni],
                transform_out->stride[long_ni]
            }
        };

        const int howmany_rank = 3;       // Loop over other directions
        fftw_iodim howmany_dims[howmany_rank];
        int j = 0;
        for (const int *io = transform_in->order+4;
             io >= transform_in->order+2; --io) {
            if (*io == long_ni) continue; // Skip transformed direction
            assert(*io < 3);
            assert(j < 2);
            howmany_dims[j].n  = transform_in->size[*io];
            howmany_dims[j].is = transform_in->stride[*io];
            howmany_dims[j].os = transform_out->stride[*io];
            ++j;
        }
        assert(j == 2);
        howmany_dims[2].n  = transform_in->size[3];
        howmany_dims[2].is = transform_in->size[4];
        howmany_dims[2].os = transform_out->size[4];

        // Find (in-ri), (in-ii), (out-ro), (out-io) for fftw_execute_split_dft
        // Store offsets to simply using the new array execute interface later
        if (f->in_place) {
            f->offset.ri = f->offset.ro = 0;
            f->offset.ii = f->offset.io = transform_in->stride[4];
        } else {
            f->offset.ri = f->offset.ro = 0;
            f->offset.ii = f->input.stride[4];
            f->offset.io = f->output.stride[4];
        }
        // FFTW_BACKWARD is FFTW_FORWARD with flipped components
        if (fftw_sign == FFTW_BACKWARD) {
            swap(&f->offset.ri, &f->offset.ii);
            swap(&f->offset.ro, &f->offset.io);
        }

        // Finally create the transformation plan
        f->plan_fftw = fftw_plan_guru_split_dft(
                sizeof(dims)/sizeof(dims[0]), dims,
                sizeof(howmany_dims)/sizeof(howmany_dims[0]), howmany_dims,
                in + f->offset.ri, in + f->offset.ii,
                out + f->offset.ro, out + f->offset.io,
                fftw_rigor_flags | FFTW_DESTROY_INPUT);
        if (UNDERLING_UNLIKELY(f->plan_fftw == NULL)) {
            underling_fftw_plan_destroy(f);
            UNDERLING_ERROR_NULL(
                    "FFTW returned NULL FFT plan", UNDERLING_ESANITY);
        }
    }

    return f;
}

underling_fftw_plan
underling_fftw_plan_create_c2r_backward(
        const underling_problem problem,
        int long_ni,
        underling_real * in,
        underling_real * out,
        unsigned fftw_rigor_flags,
        unsigned packed_flags)
{
    if (!packed_flags) packed_flags = UNDERLING_FFTW_PACKED_NONE; // Default

    underling_extents e = underling_local_extents(problem, long_ni);

    const underling_fftw_extents input
        = create_underling_fftw_extents_for_complex(e, long_ni);

    adjust_for_fast_stride_in_long_direction(&e, long_ni);

    underling_fftw_extents output
        = create_underling_fftw_extents_for_real(e, long_ni);

    if (    (long_ni == 2 && packed_flags & UNDERLING_FFTW_PACKED_LONG_N2)
         || (long_ni == 0 && packed_flags & UNDERLING_FFTW_PACKED_LONG_N0)) {
        if (UNDERLING_UNLIKELY(in == out)) {
            UNDERLING_ERROR_NULL("invalid packed_flags for in-place transform",
                                 UNDERLING_EINVAL);
        } else {
            pack_strides_according_to_order(&output);
        }
    }

    if (UNDERLING_UNLIKELY(in == out && input.order[2] != long_ni)) {
        UNDERLING_ERROR_NULL(
                "Creation of in-place c2r_backward plans for"
                " non-stride one directions is currently unimplemented",
                UNDERLING_ESANITY);
    }

    return underling_fftw_plan_create_c2r_backward_internal(
            long_ni, in, out, fftw_rigor_flags, input, output);
}

static
underling_fftw_plan
underling_fftw_plan_create_c2r_backward_internal(
        const int long_ni,
        underling_real * const in,
        underling_real * const out,
        unsigned fftw_rigor_flags,
        const underling_fftw_extents input,
        const underling_fftw_extents output)
{
    // Sanity check input arguments
    if (UNDERLING_UNLIKELY(long_ni < 0 || long_ni > 2)) {
        UNDERLING_ERROR_NULL("long_ni < 0 or long_ni > 2", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(in == NULL)) {
        UNDERLING_ERROR_NULL("in == NULL", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(out == NULL)) {
        UNDERLING_ERROR_NULL("out == NULL", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(fftw_rigor_flags & non_rigor_mask)) {
        UNDERLING_ERROR_NULL("FFTW non-rigor bits disallowed",
                             UNDERLING_EINVAL);
    }

    // Allocate and initialize the plan workspace
    underling_fftw_plan f = malloc(sizeof(struct underling_fftw_plan_s));
    if (UNDERLING_UNLIKELY(f == NULL)) {
        UNDERLING_ERROR_NULL("failed to allocate space for plan",
                             UNDERLING_ENOMEM);
    }
    // Copy the relevant parameters to the plan workspace
    f->long_ni        = long_ni;
    f->type           = transform_type_c2r_backward;
    f->input          = input;
    f->plan_preorder  = NULL;
    f->plan_fftw      = NULL;
    f->plan_postorder = NULL;
    f->output         = output;
    f->in_place       = (in == out);

    // Prepare the pre-ordering plan
    if (f->in_place) {
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

        assert(in == out);
        f->plan_preorder = fftw_plan_guru_r2r(
                0, NULL, howmany_rank, howmany_dims, in, out, NULL,
                fftw_rigor_flags | FFTW_DESTROY_INPUT);
        if (UNDERLING_UNLIKELY(f->plan_preorder == NULL)) {
            underling_fftw_plan_destroy(f);
            UNDERLING_ERROR_NULL(
                    "FFTW returned NULL c2r_backward preorder plan",
                    UNDERLING_ESANITY);
        }
    } else {
        // NOP: No pre-ordering necessary for out-of-place
    }

    // Create in- or out-of-place plan for the FFT
    // Special handling for in-place comes from pre/post-order plans
    {
        const fftw_iodim dims[] = {       // Transform long_ni
            {
                output.size[long_ni],     // Logical transform size
                input.stride[long_ni],
                f->in_place ? input.size[3] : output.stride[long_ni]
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
            howmany_dims[j].os =
                f->in_place ?  input.stride[*oo] : output.stride[*oo];
            ++j;
        }
        assert(j == 2);
        howmany_dims[2].n  = input.size[3];
        howmany_dims[2].is = f->in_place ? 1 : input.stride[3];
        howmany_dims[2].os = f->in_place ? 1 : output.stride[3];

        // Find (in-ri) and (in-ii) offsets for fftw_execute_split_dft_c2r
        // Store offsets to simplify using the new array execute interface
        f->offset.ri = 0;
        f->offset.ii = f->in_place ? input.size[3] : 1;

        f->plan_fftw = fftw_plan_guru_split_dft_c2r(
                rank, dims, howmany_rank, howmany_dims,
                in + f->offset.ri, in + f->offset.ii, out,
                fftw_rigor_flags | FFTW_DESTROY_INPUT);
        if (UNDERLING_UNLIKELY(f->plan_fftw == NULL)) {
            underling_fftw_plan_destroy(f);
            UNDERLING_ERROR_NULL("FFTW returned NULL c2r_backward FFT plan",
                    UNDERLING_ESANITY);
        }
    }

    // Prepare the post-ordering plan for in-place transforms
    if (f->in_place) {
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

        assert(in == out);
        f->plan_postorder = fftw_plan_guru_r2r(
                0, NULL, howmany_rank, howmany_dims, in, out, NULL,
                fftw_rigor_flags | FFTW_DESTROY_INPUT);
        if (UNDERLING_UNLIKELY(f->plan_postorder == NULL)) {
            underling_fftw_plan_destroy(f);
            UNDERLING_ERROR_NULL(
                    "FFTW returned NULL c2r_backward postorder plan",
                    UNDERLING_ESANITY);
        }
    } else {
        // NOP: No post-ordering necessary for out-of-place
    }

    return f;

}

underling_fftw_plan
underling_fftw_plan_create_r2c_forward(
        const underling_problem problem,
        int long_ni,
        underling_real * in,
        underling_real * out,
        unsigned fftw_rigor_flags,
        unsigned packed_flags)
{
    if (!packed_flags) packed_flags = UNDERLING_FFTW_PACKED_NONE; // Default

    underling_extents e = underling_local_extents(problem, long_ni);

    const underling_fftw_extents output
        = create_underling_fftw_extents_for_complex(e, long_ni);

    adjust_for_fast_stride_in_long_direction(&e, long_ni);

    underling_fftw_extents input
        = create_underling_fftw_extents_for_real(e, long_ni);

    if (    (long_ni == 2 && packed_flags & UNDERLING_FFTW_PACKED_LONG_N2)
         || (long_ni == 0 && packed_flags & UNDERLING_FFTW_PACKED_LONG_N0)) {
        if (UNDERLING_UNLIKELY(in == out)) {
            UNDERLING_ERROR_NULL("invalid packed_flags for in-place transform",
                                 UNDERLING_EINVAL);
        } else {
            pack_strides_according_to_order(&input);
        }
    }

    if (UNDERLING_UNLIKELY(in == out && output.order[2] != long_ni)) {
        UNDERLING_ERROR_NULL(
                "Creation of in-place r2c_forward plans for"
                " non-stride one directions is currently unimplemented.",
                UNDERLING_ESANITY);
    }

    return underling_fftw_plan_create_r2c_forward_internal(
            long_ni, in, out, fftw_rigor_flags, input, output);
}

static
underling_fftw_plan
underling_fftw_plan_create_r2c_forward_internal(
        const int long_ni,
        underling_real * const in,
        underling_real * const out,
        unsigned fftw_rigor_flags,
        const underling_fftw_extents input,
        const underling_fftw_extents output)
{
    // Sanity check input arguments
    if (UNDERLING_UNLIKELY(long_ni < 0 || long_ni > 2)) {
        UNDERLING_ERROR_NULL("long_ni < 0 or long_ni > 2", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(in == NULL)) {
        UNDERLING_ERROR_NULL("in == NULL", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(out == NULL)) {
        UNDERLING_ERROR_NULL("out == NULL", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(fftw_rigor_flags & non_rigor_mask)) {
        UNDERLING_ERROR_NULL("FFTW non-rigor bits disallowed",
                             UNDERLING_EINVAL);
    }

    // Create and initialize the plan workspace
    underling_fftw_plan f = malloc(sizeof(struct underling_fftw_plan_s));
    if (UNDERLING_UNLIKELY(f == NULL)) {
        UNDERLING_ERROR_NULL("failed to allocate space for plan",
                             UNDERLING_ENOMEM);
    }
    // Copy the relevant parameters to the plan workspace
    f->long_ni        = long_ni;
    f->type           = transform_type_r2c_forward;
    f->input          = input;
    f->plan_preorder  = NULL;
    f->plan_fftw      = NULL;
    f->plan_postorder = NULL;
    f->output         = output;
    f->in_place       = (in == out);

    // Prepare the pre-ordering plan
    if (f->in_place) {
        // We must always pay to reorder the output data.  Ignoring input
        // reordering process avoids touching all data a separate, third time.
        // TODO Evaluate performance impact since FFT may be non-stride 1
    } else {
        // NOP: No pre-ordering necessary for out-of-place
    }

    // Create in- or out-of-place plan for the FFT
    // Special handling for in-place comes from pre/post-order plans
    {
        const fftw_iodim dims[] = {       // Transform long_ni
            {
                input.size[long_ni],
                input.stride[long_ni],
                f->in_place ? 2*input.stride[long_ni] : output.stride[long_ni]
            }
        };
        const int rank = sizeof(dims)/sizeof(dims[0]);

        const int howmany_rank = 3;       // Loop over other directions
        fftw_iodim howmany_dims[howmany_rank];
        int j = 0;
        for (const int *io = input.order+4; io >= input.order+1; --io) {
            if (*io == long_ni) continue; // Skip transformed direction
            assert((size_t) j < sizeof(howmany_dims)/sizeof(howmany_dims[0]));
            howmany_dims[j].n  = input.size[*io];
            howmany_dims[j].is = input.stride[*io];
            howmany_dims[j].os =
                f->in_place ? input.stride[*io] : output.stride[*io];
            ++j;
        }
        assert((size_t) j == sizeof(howmany_dims)/sizeof(howmany_dims[0]));

        // Find (out-ro) and (out-i0) offsets for fftw_execute_split_dft_r2c
        // Store offsets to simplify using the new array execute interface
        f->offset.ro = 0;
        f->offset.io = f->in_place ? input.stride[long_ni] : 1;

        f->plan_fftw = fftw_plan_guru_split_dft_r2c(
                rank, dims, howmany_rank, howmany_dims,
                in, out + f->offset.ro, out + f->offset.io,
                fftw_rigor_flags | FFTW_DESTROY_INPUT);
        if (UNDERLING_UNLIKELY(f->plan_fftw == NULL)) {
            underling_fftw_plan_destroy(f);
            UNDERLING_ERROR_NULL("FFTW returned NULL r2c_forward FFT plan",
                    UNDERLING_ESANITY);
        }
    }

    // Prepare the post-ordering plan
    if (f->in_place) {
        // Use "complex-centric" sizes since they cover the contiguous region.
        const int howmany_rank = sizeof(output.size)/sizeof(output.size[0]);
        fftw_iodim howmany_dims[howmany_rank];
        int j = 0;
        for (const int *io = input.order+4; io >= input.order+2; --io) {
            assert((size_t) j < sizeof(howmany_dims)/sizeof(howmany_dims[0]));
            howmany_dims[j].n  = output.size[*io];
            howmany_dims[j].is = input.stride[*io];
            howmany_dims[j].os = output.stride[*io];
            if (*io == long_ni) {
                howmany_dims[j].is *= 2; // Transform modifies input stride
            }
            ++j;
        }
        assert((size_t) j == sizeof(howmany_dims)/sizeof(howmany_dims[0]) - 2);
        howmany_dims[3].n  = output.size[3];
        howmany_dims[3].is = 1;
        howmany_dims[3].os = 2;
        howmany_dims[4].n  = 2;
        howmany_dims[4].is = output.size[3];
        howmany_dims[4].os = 1;

        assert(in == out);
        f->plan_postorder = fftw_plan_guru_r2r(
                0, NULL, howmany_rank, howmany_dims, in, out, NULL,
                fftw_rigor_flags | FFTW_DESTROY_INPUT);
        if (UNDERLING_UNLIKELY(f->plan_postorder == NULL)) {
            underling_fftw_plan_destroy(f);
            UNDERLING_ERROR_NULL(
                    "FFTW returned NULL r2c_forward postorder plan",
                    UNDERLING_ESANITY);
        }
    } else {
        // NOP: No post-ordering necessary for out-of-place
    }

    return f;
}

underling_fftw_plan
underling_fftw_plan_create_inverse(
        const underling_fftw_plan plan_to_invert,
        underling_real * in,
        underling_real * out,
        unsigned fftw_rigor_flags)
{
    if (UNDERLING_UNLIKELY(plan_to_invert == NULL)) {
        UNDERLING_ERROR_NULL("plan_to_invert == NULL", UNDERLING_EINVAL);
    }

    underling_fftw_plan retval = NULL;

    switch (plan_to_invert->type) {
    case transform_type_c2c_forward:
        retval = underling_fftw_plan_create_c2c_internal(
                    plan_to_invert->long_ni, in, out, FFTW_BACKWARD,
                    fftw_rigor_flags,
                    plan_to_invert->output, plan_to_invert->input);
        break;
    case transform_type_c2c_backward:
        retval = underling_fftw_plan_create_c2c_internal(
                    plan_to_invert->long_ni, in, out, FFTW_FORWARD,
                    fftw_rigor_flags,
                    plan_to_invert->output, plan_to_invert->input);
        break;
    case transform_type_c2r_backward:
        retval = underling_fftw_plan_create_r2c_forward_internal(
                    plan_to_invert->long_ni, in, out,
                    fftw_rigor_flags,
                    plan_to_invert->output, plan_to_invert->input);
        break;
    case transform_type_r2c_forward:
        retval = underling_fftw_plan_create_c2r_backward_internal(
                    plan_to_invert->long_ni, in, out,
                    fftw_rigor_flags,
                    plan_to_invert->output, plan_to_invert->input);
        break;
    case transform_type_unspecified:
    default:
        UNDERLING_ERROR_NULL("Unrecognized plan_to_invert->type",
                UNDERLING_ESANITY);
    }

    // Check that the inverse indeed inverts storage correctly.  Could have
    // been an assertion, but the runtime cost is low compared to the debugging
    // nightmare if the check ever silently failed in a production build.
    if (UNDERLING_UNLIKELY(underling_fftw_extents_cmp(
                    &plan_to_invert->input,  &retval->output))) {
        UNDERLING_ERROR_REPORT("plan_to_invert->input != retval->output",
                               UNDERLING_ESANITY);
        underling_fftw_plan_destroy(retval);
        return NULL;
    }
    if (UNDERLING_UNLIKELY(underling_fftw_extents_cmp(
                    &plan_to_invert->output, &retval->input))) {
        UNDERLING_ERROR_REPORT("plan_to_invert->output != retval->input",
                               UNDERLING_ESANITY);
        underling_fftw_plan_destroy(retval);
        return NULL;
    }

    return retval;
}

underling_fftw_extents
underling_fftw_local_extents_input(
        const underling_fftw_plan plan)
{
    if (UNDERLING_UNLIKELY(plan == NULL)) {
        UNDERLING_ERROR_VAL("plan == NULL",
                UNDERLING_EINVAL, UNDERLING_FFTW_EXTENTS_INVALID);
    }

    underling_fftw_extents retval = plan->input; // Create temporary
    return retval;                               // Return temporary
}

underling_fftw_extents
underling_fftw_local_extents_output(
        const underling_fftw_plan plan)
{
    if (UNDERLING_UNLIKELY(plan == NULL)) {
        UNDERLING_ERROR_VAL("plan == NULL",
                UNDERLING_EINVAL, UNDERLING_FFTW_EXTENTS_INVALID);
    }

    underling_fftw_extents retval = plan->output; // Create temporary
    return retval;                                // Return temporary
}

static
void
underling_fftw_extents_copy(
        const underling_fftw_extents * const e,
        int *start,
        int *size,
        int *stride,
        int *order)
{
    if (start) {
        for (size_t j = 0; j < sizeof(e->start)/sizeof(e->start[0]); ++j)
            start[j] = e->start[j];
    }
    if (size) {
        for (size_t j = 0; j < sizeof(e->size)/sizeof(e->size[0]); ++j)
            size[j] = e->size[j];
    }
    if (stride) {
        for (size_t j = 0; j < sizeof(e->stride)/sizeof(e->stride[0]); ++j)
            stride[j] = e->stride[j];
    }
    if (order) {
        for (size_t j = 0; j < sizeof(e->order)/sizeof(e->order[0]); ++j)
            order[j] = e->order[j];
    }
}

int
underling_fftw_local_input(
        const underling_fftw_plan plan,
        int *start,
        int *size,
        int *stride,
        int *order)
{
    if (UNDERLING_UNLIKELY(plan == NULL)) {
        UNDERLING_ERROR("plan == NULL", UNDERLING_EINVAL);
    }

    underling_fftw_extents_copy(&plan->input, start, size, stride, order);

    return UNDERLING_SUCCESS;
}

int
underling_fftw_local_output(
        const underling_fftw_plan plan,
        int *start,
        int *size,
        int *stride,
        int *order)
{
    if (UNDERLING_UNLIKELY(plan == NULL)) {
        UNDERLING_ERROR("plan == NULL", UNDERLING_EINVAL);
    }

    underling_fftw_extents_copy(&plan->output, start, size, stride, order);

    return UNDERLING_SUCCESS;
}

int
underling_fftw_plan_execute(
        const underling_fftw_plan plan,
        underling_real * in,
        underling_real * out)
{
    if (UNDERLING_UNLIKELY(plan == NULL)) {
        UNDERLING_ERROR("plan == NULL", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(plan->plan_fftw == NULL)) {
        UNDERLING_ERROR("plan->plan_fftw == NULL", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(in == NULL)) {
        UNDERLING_ERROR("in == NULL", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(out == NULL)) {
        UNDERLING_ERROR("out == NULL", UNDERLING_EINVAL);
    }

    // Check in- vs out-of-place plan applied appropriately
    if (plan->in_place) {
        if (UNDERLING_UNLIKELY(in != out)) {
            UNDERLING_ERROR("in-place plan but in != out", UNDERLING_EINVAL);
        }
    } else if (UNDERLING_UNLIKELY(in == out)) {
        UNDERLING_ERROR("out-of-place plan but in == out", UNDERLING_EINVAL);
    }

    if (plan->plan_preorder) {
        assert(plan->in_place && in == out);
        fftw_execute_r2r(plan->plan_preorder, in, out);
    }
    switch (plan->type) {
        case transform_type_c2c_forward:
        case transform_type_c2c_backward:
            fftw_execute_split_dft(plan->plan_fftw,
                                   in  + plan->offset.ri,
                                   in  + plan->offset.ii,
                                   out + plan->offset.ro,
                                   out + plan->offset.io);
            break;
        case transform_type_c2r_backward:
            fftw_execute_split_dft_c2r(plan->plan_fftw,
                                       in + plan->offset.ri,
                                       in + plan->offset.ii,
                                       out);
            break;
        case transform_type_r2c_forward:
            fftw_execute_split_dft_r2c(plan->plan_fftw,
                                       in,
                                       out + plan->offset.ro,
                                       out + plan->offset.io);
            break;
        default:
            UNDERLING_ERROR("Unknown plan->type", UNDERLING_ESANITY);
    }

    if (plan->plan_postorder) {
        assert(plan->in_place && in == out);
        fftw_execute_r2r(plan->plan_postorder, in, out);
    }

    return UNDERLING_SUCCESS;
}

void
underling_fftw_plan_destroy(
        underling_fftw_plan plan)
{
    if (plan) {
        if (plan->plan_preorder) {
            fftw_destroy_plan(plan->plan_preorder);
            plan->plan_preorder = NULL;
        }
        if (plan->plan_fftw) {
            fftw_destroy_plan(plan->plan_fftw);
            plan->plan_fftw = NULL;
        }
        if (plan->plan_postorder) {
            fftw_destroy_plan(plan->plan_postorder);
            plan->plan_postorder = NULL;
        }
        free(plan);
    }
}

void
underling_fftw_fprint_extents(
        const underling_fftw_extents *extents,
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
underling_fftw_fprint_plan(
        const underling_fftw_plan plan,
        FILE *output_file)
{
    fprintf(output_file, "{underling_fftw_plan:");
    if (!plan) {
        fprintf(output_file, "NULL");
    } else {
        fprintf(output_file,"long_ni=%d,type=%d",plan->long_ni, plan->type);
        fprintf(output_file,",{input=");
        underling_fftw_fprint_extents(&plan->input,output_file);
        fprintf(output_file,"}");
        if (plan->plan_preorder) {
            fprintf(output_file, "{plan_preorder:");
            fftw_fprint_plan(plan->plan_preorder, output_file);
            fprintf(output_file, "}");
        }
        if (plan->plan_fftw) {
            fprintf(output_file, "{plan_fftw:");
            fftw_fprint_plan(plan->plan_fftw, output_file);
            fprintf(output_file, "}");
        }
        if (plan->plan_postorder) {
            fprintf(output_file, "{plan_postorder:");
            fftw_fprint_plan(plan->plan_postorder, output_file);
            fprintf(output_file, "}");
        }
        fprintf(output_file,"{output=");
        underling_fftw_fprint_extents(&plan->output,output_file);
        fprintf(output_file,"},");
        if (plan->in_place) {
            fprintf(output_file,"in-place");
        } else {
            fprintf(output_file,"out-of-place");
        }
        fprintf(output_file,",{offset:ri=%d,ii=%d,ro=%d,io=%d}",
                plan->offset.ri, plan->offset.ii,
                plan->offset.ro, plan->offset.io);
    }
    fprintf(output_file, "}");
}
