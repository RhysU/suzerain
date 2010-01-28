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
 * underling.c: A parallel, three dimensional FFT library atop FFTW3 MPI
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <config.h>
#include <suzerain/common.h>
#pragma hdrstop
#include <suzerain/error.h>
#include <suzerain/underling.h>
#include <fftw3-mpi.h>

// TODO Ensure grid/problem compatibility when both provided to methods!
// TODO Add wisdom broadcasting and FFTW_WISDOM_ONLY handling

// *******************************************************************
// INTERNAL STRUCTS INTERNAL STRUCTS INTERNAL STRUCTS INTERNAL STRUCTS
// *******************************************************************

struct underling_grid_s {
    int      n[3];
    int      pA;
    int      pB;
    MPI_Comm g_comm;                  // Instance owns MPI resource
    MPI_Comm pA_comm;                 // Instance owns MPI resource
    MPI_Comm pB_comm;                 // Instance owns MPI resource
};

typedef struct underling_transpose_s * underling_transpose; // Internal!
struct underling_transpose_s {
    ptrdiff_t d[2];
    ptrdiff_t howmany;
    ptrdiff_t block[2];
    MPI_Comm  comm;                   // Instance owns MPI resource
    unsigned  flags;
    ptrdiff_t local[2];
    ptrdiff_t local_start[2];
    ptrdiff_t local_size;
};

struct underling_problem_s {
    int howmany;                      // # of real values to transpose
    underling_extents long_n[3];      // Layout details for n{0,1,2} long
    underling_transpose backwardA;    // n2 long to n1 long
    underling_transpose backwardB;    // n1 long to n0 long
    underling_transpose forwardB;     // n0 long to n1 long
    underling_transpose forwardA;     // n1 long to n2 long
    size_t local_memory;              // Max of all transpose local sizes
};

struct underling_plan_s {
    fftw_plan plan_backwardA;
    fftw_plan plan_backwardB;
    fftw_plan plan_forwardA;
    fftw_plan plan_forwardB;
    underling_real *data;              // API end-user owns resources
};

// ********************************************************************
// INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL
// ********************************************************************

MPI_Comm
underling_MPI_Comm_dup_with_name(MPI_Comm comm);

fftw_plan
underling_fftw_plan_nop();

underling_transpose
underling_transpose_create(
        ptrdiff_t d0,
        ptrdiff_t d1,
        ptrdiff_t howmany,
        ptrdiff_t block0,
        ptrdiff_t block1,
        MPI_Comm comm,
        unsigned flags);

underling_transpose
underling_transpose_create_inverse(
        const underling_transpose forward);

void
underling_transpose_destroy(
        underling_transpose transpose);

fftw_plan
underling_transpose_fftw_plan(
        const underling_transpose transpose,
        underling_real *in,
        underling_real *out,
        unsigned fftw_flags);

size_t
underling_local_memory_allreduce(
        const underling_grid    grid,
        const underling_problem problem,
        MPI_Op op);

void
underling_fprint_transpose(
        const underling_transpose transpose,
        FILE *output_file);

// **************************************************************************
// IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION
// **************************************************************************

const underling_extents UNDERLING_EXTENTS_INVALID = {
    {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}
};

static
MPI_Comm
underling_MPI_Comm_dup_with_name(MPI_Comm comm)
{
    if (comm == MPI_COMM_NULL) return MPI_COMM_NULL;

    char buffer[MPI_MAX_OBJECT_NAME];
    int resultlen = 0;
    const int get_name_error = MPI_Comm_get_name(comm, buffer, &resultlen);
    if (get_name_error) {
        SUZERAIN_MPICHKR(get_name_error /* MPI_Comm_get_name */);
        return MPI_COMM_NULL;
    }

    MPI_Comm retval = MPI_COMM_NULL;

    const int dup_error = MPI_Comm_dup(comm, &retval);
    if (dup_error) {
        SUZERAIN_MPICHKR(dup_error /* MPI_Comm_dup */);
        return MPI_COMM_NULL;
    }

    if (resultlen > 0) {
        const int set_name_error = MPI_Comm_set_name(retval, buffer);
        if (set_name_error) {
            SUZERAIN_MPICHKR(set_name_error /* MPI_Comm_set_name */);
            SUZERAIN_MPICHKR(MPI_Comm_free(&retval));
            return MPI_COMM_NULL;
        }
    }

    return retval;
}

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

underling_grid
underling_grid_create(
        MPI_Comm comm,
        int n0,
        int n1,
        int n2,
        int pA,
        int pB)
{
    // Sanity check incoming arguments
    if (SUZERAIN_UNLIKELY(comm == MPI_COMM_NULL)) {
        SUZERAIN_ERROR_NULL("comm != MPI_COMM_NULL required", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(n0 < 1)) {
        SUZERAIN_ERROR_NULL("n0 >= 1 required", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(n1 < 1)) {
        SUZERAIN_ERROR_NULL("n1 >= 1 required", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(n2 < 1)) {
        SUZERAIN_ERROR_NULL("n2 >= 1 required", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(pA < 0)) {
        SUZERAIN_ERROR_NULL("pA >= 0 required", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(pB < 0)) {
        SUZERAIN_ERROR_NULL("pB >= 0 required", SUZERAIN_EINVAL);
    }

    // Get number of processors in the communicator
    int nproc;
    SUZERAIN_MPICHKN(MPI_Comm_size(comm, &nproc));

    // Create a balanced processor grid if not specified by pA, pB != 0
    {
        int dims[2] = { pA, pB };
        SUZERAIN_MPICHKN(MPI_Dims_create(nproc, 2, dims));
        // If both directions automatic, ensure dims[0] <= dims[1]
        if (pA == 0 && pB == 0 && dims[0] > dims[1]) {
            const int tmp = dims[0]; dims[0] = dims[1]; dims[1] = tmp;
        }
        pA = dims[0];
        pB = dims[1];
        if (SUZERAIN_UNLIKELY(pA * pB != nproc)) {
            char reason[127];
            snprintf(reason, sizeof(reason)/sizeof(reason[0]),
                    "Invalid processor grid: pA {%d} * pB {%d} != nproc {%d}",
                    pA, pB, nproc);
            SUZERAIN_ERROR_NULL(reason, SUZERAIN_EFAILED);
        }
    }

    // Clone the communicator and create the 2D Cartesian topology
    MPI_Comm g_comm = MPI_COMM_NULL;
    {
        int dims[2]     = { pA, pB };
        int periodic[2] = { 0, 0 };
        SUZERAIN_MPICHKN(MPI_Cart_create(
                comm, 2, dims, periodic, 1/*reordering allowed*/, &g_comm));
    }

    // Create communicator for the PA direction
    MPI_Comm pA_comm = MPI_COMM_NULL;
    {
        int remain_dims[2] = { 1, 0 };
        SUZERAIN_MPICHKN(MPI_Cart_sub(g_comm, remain_dims, &pA_comm));
    }
    // Find the rank of this process within g_comm
    int pA_rank;
    SUZERAIN_MPICHKN(MPI_Comm_rank(pA_comm, &pA_rank));

    // Create communicator for the PB direction
    MPI_Comm pB_comm = MPI_COMM_NULL;
    {
        int remain_dims[2] = { 0, 1 };
        SUZERAIN_MPICHKN(MPI_Cart_sub(g_comm, remain_dims, &pB_comm));
    }
    // Find the rank of this process within g_comm
    int pB_rank;
    SUZERAIN_MPICHKN(MPI_Comm_rank(pB_comm, &pB_rank));

    // Name the three new communicators something mildly descriptive
    char buffer[MPI_MAX_OBJECT_NAME];
    snprintf(buffer, sizeof(buffer)/sizeof(buffer[0]),
            "uGComm%dx%d", pA, pB);
    SUZERAIN_MPICHKN(MPI_Comm_set_name(g_comm, buffer));
    snprintf(buffer, sizeof(buffer)/sizeof(buffer[0]),
            "uPACommXx%d", pB_rank);
    SUZERAIN_MPICHKN(MPI_Comm_set_name(pA_comm, buffer));
    snprintf(buffer, sizeof(buffer)/sizeof(buffer[0]),
            "uPBComm%dxX", pA_rank);
    SUZERAIN_MPICHKN(MPI_Comm_set_name(pB_comm, buffer));

    // Create and initialize the grid workspace
    underling_grid g = calloc(1, sizeof(struct underling_grid_s));
    if (SUZERAIN_UNLIKELY(g == NULL)) {
        SUZERAIN_ERROR_NULL("failed to allocate space for grid",
                             SUZERAIN_ENOMEM);
    }
    // Copy the grid parameters to the grid workspace
    g->n[0]        = n0;
    g->n[1]        = n1;
    g->n[2]        = n2;
    g->pA          = pA;
    g->pB          = pB;
    g->g_comm      = g_comm;
    g->pA_comm     = pA_comm;
    g->pB_comm     = pB_comm;

    return g;
}

void
underling_grid_destroy(underling_grid grid)
{
    if (grid) {
        if (grid->g_comm != MPI_COMM_NULL) {
            SUZERAIN_MPICHKR(MPI_Comm_free(&grid->g_comm));
            grid->g_comm = MPI_COMM_NULL;
        }
        if (grid->pA_comm != MPI_COMM_NULL) {
            SUZERAIN_MPICHKR(MPI_Comm_free(&grid->pA_comm));
            grid->pA_comm = MPI_COMM_NULL;
        }
        if (grid->pB_comm != MPI_COMM_NULL) {
            SUZERAIN_MPICHKR(MPI_Comm_free(&grid->pB_comm));
            grid->pB_comm = MPI_COMM_NULL;
        }
        free(grid);
    }
}

static
underling_transpose
underling_transpose_create(
        ptrdiff_t d0,
        ptrdiff_t d1,
        ptrdiff_t howmany,
        ptrdiff_t block0,
        ptrdiff_t block1,
        MPI_Comm comm,
        unsigned flags)
{
    if (SUZERAIN_UNLIKELY(comm == MPI_COMM_NULL)) {
        SUZERAIN_ERROR_NULL(
                "Unable to create transpose for comm == MPI_COMM_NULL",
                SUZERAIN_EINVAL);
    }

    // Create and initialize the transpose workspace
    underling_transpose t = calloc(1, sizeof(struct underling_transpose_s));
    if (SUZERAIN_UNLIKELY(t == NULL)) {
        SUZERAIN_ERROR_NULL("failed to allocate space for transpose",
                             SUZERAIN_ENOMEM);
    }

    // Fix struct values known from arguments
    t->d[0]     = d0;
    t->d[1]     = d1;
    t->howmany  = howmany;
    t->block[0] = block0;
    t->block[1] = block1;
    t->comm     = underling_MPI_Comm_dup_with_name(comm);
    t->flags    = flags;

    if (t->comm == MPI_COMM_NULL) {
        underling_transpose_destroy(t);
        SUZERAIN_ERROR_NULL("Detected MPI_COMM_NULL in t->comm",
                            SUZERAIN_ESANITY);
    }

    if (SUZERAIN_UNLIKELY(t->d[0] == 0 || t->d[1] == 0)) {
        // Trivial transpose required;
        // fftw_mpi_local_size_many_transpose divides-by-zero on zero size
        t->local[0]       = 0;
        t->local[1]       = 0;
        t->local_start[0] = 0;
        t->local_start[1] = 0;
        t->local_size     = 0;
    } else {
        // Fix transpose details obtained from FFTW MPI call
        t->local_size = fftw_mpi_local_size_many_transposed(
                                                /*rank*/2,
                                                t->d,
                                                t->howmany,
                                                t->block[0],
                                                t->block[1],
                                                t->comm,
                                                &t->local[0],
                                                &t->local_start[0],
                                                &t->local[1],
                                                &t->local_start[1]);
    }

    return t;
}

static
underling_transpose
underling_transpose_create_inverse(
        const underling_transpose forward)
{
    if (SUZERAIN_UNLIKELY(forward == NULL)) {
        SUZERAIN_ERROR_NULL("forward == NULL", SUZERAIN_EINVAL);
    }
    // TODO Handle FFTW_MPI_TRANSPOSED_IN, FFTW_MPI_TRANSPOSED_OUT correctly
    if (SUZERAIN_UNLIKELY(forward->flags)) {
        SUZERAIN_ERROR_NULL("Nontrivial transpose flags not yet implemented!",
                SUZERAIN_ESANITY);
    }

    // Create and initialize the transpose workspace
    underling_transpose backward
        = calloc(1, sizeof(struct underling_transpose_s));
    if (SUZERAIN_UNLIKELY(backward == NULL)) {
        SUZERAIN_ERROR_NULL("failed to allocate space for transpose",
                             SUZERAIN_ENOMEM);
    }

    // Fix struct values known from inverting forward plan
    backward->d[0]     = forward->d[1];
    backward->d[1]     = forward->d[0];
    backward->howmany  = forward->howmany;
    backward->block[0] = forward->block[1];
    backward->block[1] = forward->block[0];
    backward->comm     = underling_MPI_Comm_dup_with_name(forward->comm);
    backward->flags    = forward->flags;

    if (backward->comm == MPI_COMM_NULL) {
        underling_transpose_destroy(backward);
        SUZERAIN_ERROR_NULL("Detected MPI_COMM_NULL in backward->comm",
                            SUZERAIN_ESANITY);
    }

    if (SUZERAIN_UNLIKELY(backward->d[0] == 0 || backward->d[1] == 0)) {
        // Trivial transpose required;
        // fftw_mpi_local_size_many_transpose divides-by-zero on zero size
        backward->local[0]       = 0;
        backward->local[1]       = 0;
        backward->local_start[0] = 0;
        backward->local_start[1] = 0;
        backward->local_size     = 0;
    } else {
        // Fix transpose details obtained from FFTW MPI call
        backward->local_size = fftw_mpi_local_size_many_transposed(
                                                /*rank*/2,
                                                backward->d,
                                                backward->howmany,
                                                backward->block[0],
                                                backward->block[1],
                                                backward->comm,
                                                &backward->local[0],
                                                &backward->local_start[0],
                                                &backward->local[1],
                                                &backward->local_start[1]);
    }

    return backward;
}

static
void
underling_transpose_destroy(
        underling_transpose transpose)
{
    if (transpose) {
        if (transpose->comm != MPI_COMM_NULL) {
            SUZERAIN_MPICHKR(MPI_Comm_free(&transpose->comm));
            transpose->comm = MPI_COMM_NULL;
        }
        free(transpose);
    }
}

static
fftw_plan
underling_transpose_fftw_plan(
        const underling_transpose transpose,
        underling_real *in,
        underling_real *out,
        unsigned fftw_flags)
{
    if (SUZERAIN_UNLIKELY(transpose == NULL)) {
        SUZERAIN_ERROR_NULL("transpose == NULL", SUZERAIN_EINVAL);
    }

    if (SUZERAIN_UNLIKELY(transpose->d[0] == 0 || transpose->d[1] == 0)) {
        // Trivial transpose required;
        // fftw_mpi_plan_many_transpose returns NULL on trivial input
        return underling_fftw_plan_nop();
    } else {
        return fftw_mpi_plan_many_transpose(transpose->d[0],
                                            transpose->d[1],
                                            transpose->howmany,
                                            transpose->block[0],
                                            transpose->block[1],
                                            in,
                                            out,
                                            transpose->comm,
                                            transpose->flags | fftw_flags);
    }
}

underling_problem
underling_problem_create(
        underling_grid grid,
        int howmany)
{
    if (SUZERAIN_UNLIKELY(grid == NULL)) {
        SUZERAIN_ERROR_NULL("grid == NULL", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(howmany < 1)) {
        SUZERAIN_ERROR_NULL("howmany >= 1 required", SUZERAIN_EINVAL);
    }

    // Create and initialize the problem workspace
    underling_problem p = calloc(1, sizeof(struct underling_problem_s));
    if (SUZERAIN_UNLIKELY(p == NULL)) {
        SUZERAIN_ERROR_NULL("failed to allocate space for problem",
                             SUZERAIN_ENOMEM);
    }
    // Copy the problem parameters to the problem workspace
    p->howmany = howmany;

    // Global pencil decomposition details
    // -------------------------------------------------------
    // Long in n2:                        (n0/pB x n1/pA) x n2
    // Long in n1: n2/pA x (n0/pB x n1) = (n2/pA x n0/pB) x n1
    // Long in n0: n1/pB x (n2/pA x n0) = (n1/pB x n2/pA) x n0

    // Fix {n2,n1,n0} dimension details in p->long_{n2,n1,n0}
    p->long_n[2].size[2]  = grid->n[2];
    p->long_n[2].start[2] = 0;
    p->long_n[1].size[1]  = grid->n[1];
    p->long_n[1].start[1] = 0;
    p->long_n[0].size[0]  = grid->n[0];
    p->long_n[0].start[0] = 0;

    // Decompose {n0,n1}/pB and store details in p->long_{(n2,n1),n0}
    {
        ptrdiff_t local_d0, local_d0_start, local_d1, local_d1_start;
        ptrdiff_t dB[2] = {grid->n[0], grid->n[1]}; // Never performed
        fftw_mpi_local_size_many_transposed(2, dB, p->howmany,
                FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, grid->pB_comm,
                &local_d0, &local_d0_start, &local_d1, & local_d1_start);
        assert(local_d0       <= INT_MAX);
        assert(local_d0_start <= INT_MAX);
        assert(local_d1       <= INT_MAX);
        assert(local_d1_start <= INT_MAX);
        p->long_n[2].size[0]  = local_d0;
        p->long_n[2].start[0] = local_d0_start;
        p->long_n[1].size[0]  = local_d0;
        p->long_n[1].start[0] = local_d0_start;
        p->long_n[0].size[1]  = local_d1;
        p->long_n[0].start[1] = local_d1_start;
    }

    // Decompose {n1,n2}/pA and store details in p->long_{n2,(n1,n0)}
    {
        ptrdiff_t local_d0, local_d0_start, local_d1, local_d1_start;
        ptrdiff_t dA[2] = {grid->n[1], grid->n[2]}; // Never performed
        fftw_mpi_local_size_many_transposed(2, dA, p->howmany,
                FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, grid->pA_comm,
                &local_d0, &local_d0_start, &local_d1, &local_d1_start);
        assert(local_d0       <= INT_MAX);
        assert(local_d0_start <= INT_MAX);
        assert(local_d1       <= INT_MAX);
        assert(local_d1_start <= INT_MAX);
        p->long_n[2].size[1]  = local_d0;
        p->long_n[2].start[1] = local_d0_start;
        p->long_n[1].size[2]  = local_d1;
        p->long_n[1].start[2] = local_d1_start;
        p->long_n[0].size[2]  = local_d1;
        p->long_n[0].start[2] = local_d1_start;
    }

    // Determine all necessary strides for row-major storage
    // -----------------------------------------------------
    // Compute strides when long in n2: (n0/pB x n1/pA) x n2
    p->long_n[2].stride[2] = p->howmany;
    p->long_n[2].stride[1] = p->long_n[2].stride[2] * p->long_n[2].size[2];
    p->long_n[2].stride[0] = p->long_n[2].stride[1] * p->long_n[2].size[1];
    p->long_n[2].strideorder[0] = 2; // Fastest
    p->long_n[2].strideorder[1] = 1;
    p->long_n[2].strideorder[2] = 0; // Slowest
    // Compute strides when long in n1: (n2/pA x n0/pB) x n1
    p->long_n[1].stride[1] = p->howmany;
    p->long_n[1].stride[0] = p->long_n[1].stride[1] * p->long_n[1].size[1];
    p->long_n[1].stride[2] = p->long_n[1].stride[0] * p->long_n[1].size[0];
    p->long_n[1].strideorder[0] = 1; // Fastest
    p->long_n[1].strideorder[1] = 0;
    p->long_n[1].strideorder[2] = 2; // Slowest
    // Compute strides when long in n0: (n1/pB x n2/pA) x n0
    p->long_n[0].stride[0] = p->howmany;
    p->long_n[0].stride[2] = p->long_n[0].stride[0] * p->long_n[0].size[0];
    p->long_n[0].stride[1] = p->long_n[0].stride[2] * p->long_n[0].size[2];
    p->long_n[0].strideorder[0] = 0; // Fastest
    p->long_n[0].strideorder[1] = 2;
    p->long_n[0].strideorder[2] = 1; // Slowest

    // Compute extent when long in each direction; redundant but convenient
    for (int i = 0; i < 3; ++i) {
        p->long_n[i].extent =
              p->howmany
            * p->long_n[i].size[0]
            * p->long_n[i].size[1]
            * p->long_n[i].size[2];
    }

    // Transpose pA details: (n0/pB x n1/pA) x n2 to n2/pA x (n0/pB x n1)
    const ptrdiff_t pA_d[2] = { p->long_n[2].size[0] * grid->n[1],
                                grid->n[2] };
    ptrdiff_t pA_block[2]   = { p->long_n[2].size[0] * p->long_n[2].size[1],
                                p->long_n[1].size[2] };
    SUZERAIN_MPICHKN(MPI_Bcast(pA_block, 2, MPI_LONG, 0, grid->pA_comm));

    // Transpose pB details: (n2/pA x n0/pB) x n1 to n1/pB x (n2/pA x n0)
    const ptrdiff_t pB_d[2] = { p->long_n[1].size[2] * grid->n[0],
                                grid->n[1] };
    ptrdiff_t pB_block[2]   = { p->long_n[1].size[2] * p->long_n[1].size[0],
                                p->long_n[0].size[1] };
    SUZERAIN_MPICHKN(MPI_Bcast(pB_block, 2, MPI_LONG, 0, grid->pB_comm));

    // Wave towards physical MPI transpose: long in n2 to long in n1
    p->backwardA = underling_transpose_create(pA_d[0],
                                              pA_d[1],
                                              p->howmany,
                                              pA_block[0],
                                              pA_block[1],
                                              grid->pA_comm,
                                              /*flags*/0);
    if (SUZERAIN_UNLIKELY(p->backwardA == NULL)) {
        underling_problem_destroy(p);
        SUZERAIN_ERROR_NULL("failed creating p->backwardA",
                SUZERAIN_EFAILED);
    }

    // Wave towards physical MPI transpose: long in n1 to long in n0
    p->backwardB = underling_transpose_create(pB_d[0],
                                              pB_d[1],
                                              p->howmany,
                                              pB_block[0],
                                              pB_block[1],
                                              grid->pB_comm,
                                              /*flags*/0);
    if (SUZERAIN_UNLIKELY(p->backwardB == NULL)) {
        underling_problem_destroy(p);
        SUZERAIN_ERROR_NULL("failed creating p->backwardB",
                SUZERAIN_EFAILED);
    }

    // Physical towards wave MPI transpose: long in n0 to long in n1
    p->forwardB = underling_transpose_create_inverse(p->backwardB);
    if (SUZERAIN_UNLIKELY(p->forwardB == NULL)) {
        underling_problem_destroy(p);
        SUZERAIN_ERROR_NULL("failed creating p->forwardB",
                SUZERAIN_EFAILED);
    }

    // Physical towards wave MPI transpose: long in n1 to long in n2
    p->forwardA = underling_transpose_create_inverse(p->backwardA);
    if (SUZERAIN_UNLIKELY(p->forwardA == NULL)) {
        underling_problem_destroy(p);
        SUZERAIN_ERROR_NULL("failed creating p->forwardA",
                SUZERAIN_EFAILED);
    }

    // p->local_memory is overall maximum of all local_size values
    p->local_memory = p->backwardA->local_size;
    if (p->local_memory < p->backwardB->local_size) {
        p->local_memory = p->backwardB->local_size;
    }
    if (p->local_memory < p->forwardB->local_size) {
        p->local_memory = p->forwardB->local_size;
    }
    if (p->local_memory < p->forwardA->local_size) {
        p->local_memory = p->forwardA->local_size;
    }

    return p;
}

size_t
underling_local_memory(
        const underling_problem problem)
{
    if (SUZERAIN_UNLIKELY(problem == NULL)) {
        SUZERAIN_ERROR_VAL("problem == NULL", SUZERAIN_EINVAL, 0);
    }

    return problem->local_memory;
}

static
size_t
underling_local_memory_allreduce(
        const underling_grid    grid,
        const underling_problem problem,
        MPI_Op op)
{
    if (SUZERAIN_UNLIKELY(grid == NULL)) {
        SUZERAIN_ERROR_VAL("grid == NULL", SUZERAIN_EINVAL, 0);
    }
    if (SUZERAIN_UNLIKELY(problem == NULL)) {
        SUZERAIN_ERROR_VAL("problem == NULL", SUZERAIN_EINVAL, 0);
    }

    // Use unsigned long values for safety in heterogeneous environments.
    assert(((size_t)-1) <= ((unsigned long int)-1)); // "Type safety"

    unsigned long int retval = underling_local_memory(problem);
    const int error = MPI_Allreduce(
            MPI_IN_PLACE, &retval, 1, MPI_UNSIGNED_LONG,
            op, grid->g_comm);
    if (SUZERAIN_UNLIKELY(error)) {
        SUZERAIN_MPICHKR(error /* allreduce local_memory */);
        retval = 0;
    }

    return retval;
}

size_t
underling_local_memory_maximum(
        const underling_grid    grid,
        const underling_problem problem)
{
    return underling_local_memory_allreduce(grid, problem, MPI_MAX);
}

size_t
underling_local_memory_minimum(
        const underling_grid    grid,
        const underling_problem problem)
{
    return underling_local_memory_allreduce(grid, problem, MPI_MIN);
}

size_t
underling_global_memory(
        const underling_grid    grid,
        const underling_problem problem)
{
    return underling_local_memory_allreduce(grid, problem, MPI_SUM);
}

size_t
underling_global_memory_optimum(
        const underling_grid    grid,
        const underling_problem problem)
{
    if (SUZERAIN_UNLIKELY(grid == NULL)) {
        SUZERAIN_ERROR_VAL("grid == NULL", SUZERAIN_EINVAL, 0);
    }
    if (SUZERAIN_UNLIKELY(problem == NULL)) {
        SUZERAIN_ERROR_VAL("problem == NULL", SUZERAIN_EINVAL, 0);
    }
    return problem->howmany * grid->n[0] * grid->n[1] * grid->n[2];
}

size_t
underling_local_memory_optimum(
        const underling_grid    grid,
        const underling_problem problem)
{
    const size_t memory      = underling_global_memory_optimum(grid, problem);
    const size_t nprocessors = grid->pA * grid->pB;
    const size_t result      = memory / nprocessors;
    const size_t remainder   = memory % nprocessors;

    int g_rank;
    const int error = MPI_Comm_rank(grid->g_comm, &g_rank);
    if (error) {
        SUZERAIN_MPICHKR(error /* MPI_Comm_rank */);
        return 0;
    }

    // Round up result iff necessary and our 2D rank is low enough
    return result + (!!remainder)*(g_rank <= remainder);
}

size_t
underling_local(
        const underling_problem problem,
        int i,
        int *start,
        int *size,
        int *stride,
        int *strideorder)
{
    if (SUZERAIN_UNLIKELY(i < 0 || i > 2)) {
        SUZERAIN_ERROR_VAL("i < 0 or i > 2", SUZERAIN_EINVAL, 0);
    }
    if (SUZERAIN_UNLIKELY(problem == NULL)) {
        SUZERAIN_ERROR_VAL("problem == NULL", SUZERAIN_EINVAL, 0);
    }

    const underling_extents * const e = &problem->long_n[i];

    if (start) {
        for (int j = 0; j < 3; ++j)
            start[j] = e->start[j];
    }
    if (size) {
        for (int j = 0; j < 3; ++j)
            size[j] = e->size[j];
    }
    if (stride) {
        for (int j = 0; j < 3; ++j)
            stride[j] = e->stride[j];
    }
    if (strideorder) {
        for (int j = 0; j < 3; ++j)
            strideorder[j] = e->strideorder[j];
    }

    return e->extent;
}

underling_extents
underling_local_extents(
        const underling_problem problem,
        int i)
{
    if (SUZERAIN_UNLIKELY(i < 0 || i > 2)) {
        SUZERAIN_ERROR_VAL("i < 0 or i > 2",
                SUZERAIN_EINVAL, UNDERLING_EXTENTS_INVALID);
    }
    if (SUZERAIN_UNLIKELY(problem == NULL)) {
        SUZERAIN_ERROR_VAL("problem == NULL",
                SUZERAIN_EINVAL, UNDERLING_EXTENTS_INVALID);
    }

    underling_extents retval = problem->long_n[i]; // Create temporary
    return retval;                                 // Return temporary
}

void
underling_problem_destroy(
        underling_problem problem)
{
    if (problem) {
        if (problem->backwardA) {
            underling_transpose_destroy(problem->backwardA);
            problem->backwardA = NULL;
        }
        if (problem->backwardB) {
            underling_transpose_destroy(problem->backwardB);
            problem->backwardB = NULL;
        }
        if (problem->forwardB) {
            underling_transpose_destroy(problem->forwardB);
            problem->forwardB = NULL;
        }
        if (problem->forwardA) {
            underling_transpose_destroy(problem->forwardA);
            problem->forwardA = NULL;
        }
        free(problem);
    }
}

underling_plan
underling_plan_create(
        const underling_problem problem,
        underling_real * data,
        unsigned transform_flags,
        unsigned rigor_flags)
{
    if (SUZERAIN_UNLIKELY(problem == NULL)) {
        SUZERAIN_ERROR_NULL("problem == NULL", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(data == NULL)) {
        SUZERAIN_ERROR_NULL("data == NULL", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(transform_flags & ~UNDERLING_TRANSPOSE_ALL)) {
        SUZERAIN_ERROR_NULL(
            "transform_flags contains non-direction bit", SUZERAIN_EINVAL);
    }
    const unsigned non_rigor_mask =   ~FFTW_ESTIMATE
                                    & ~FFTW_MEASURE
                                    & ~FFTW_PATIENT
                                    & ~FFTW_EXHAUSTIVE
                                    & ~FFTW_WISDOM_ONLY;
    if (SUZERAIN_UNLIKELY(rigor_flags & non_rigor_mask)) {
        SUZERAIN_ERROR_NULL("FFTW non-rigor bits disallowed", SUZERAIN_EINVAL);
    }

    // Be ready to execute all transforms if trivial flag provided
    if (transform_flags == 0) {
        transform_flags = UNDERLING_TRANSPOSE_ALL;
    }

    // Create and initialize the plan workspace
    underling_plan p = calloc(1, sizeof(struct underling_plan_s));
    if (SUZERAIN_UNLIKELY(p == NULL)) {
        SUZERAIN_ERROR_NULL("failed to allocate space for plan",
                             SUZERAIN_ENOMEM);
    }
    // Copy the problem parameters to the problem workspace
    p->data = data;

    // Create the requested FFTW MPI plans
    if (transform_flags | UNDERLING_TRANSPOSE_LONG_N2_TO_LONG_N1) {
        p->plan_backwardA = underling_transpose_fftw_plan(
                problem->backwardA, p->data, p->data, rigor_flags);
        if (SUZERAIN_UNLIKELY(p->plan_backwardA == NULL)) {
            underling_plan_destroy(p);
            SUZERAIN_ERROR_NULL(
                    "FFTW MPI returned NULL plan: plan_backwardA",
                    SUZERAIN_EFAILED);
        }
    }

    if (transform_flags | UNDERLING_TRANSPOSE_LONG_N1_TO_LONG_N0) {
        p->plan_backwardB = underling_transpose_fftw_plan(
                problem->backwardB, p->data, p->data, rigor_flags);
        if (SUZERAIN_UNLIKELY(p->plan_backwardB == NULL)) {
            underling_plan_destroy(p);
            SUZERAIN_ERROR_NULL(
                    "FFTW MPI returned NULL plan: plan_backwardB",
                    SUZERAIN_EFAILED);
        }
    }

    if (transform_flags | UNDERLING_TRANSPOSE_LONG_N0_TO_LONG_N1) {
        p->plan_forwardB = underling_transpose_fftw_plan(
                problem->forwardB, p->data, p->data, rigor_flags);
        if (SUZERAIN_UNLIKELY(p->plan_forwardB == NULL)) {
            underling_plan_destroy(p);
            SUZERAIN_ERROR_NULL(
                    "FFTW MPI returned NULL plan: plan_forwardB",
                    SUZERAIN_EFAILED);
        }
    }

    if (transform_flags | UNDERLING_TRANSPOSE_LONG_N1_TO_LONG_N2) {
        p->plan_forwardA = underling_transpose_fftw_plan(
                problem->forwardA, p->data, p->data, rigor_flags);
        if (SUZERAIN_UNLIKELY(p->plan_forwardA == NULL)) {
            underling_plan_destroy(p);
            SUZERAIN_ERROR_NULL(
                    "FFTW MPI returned NULL plan: plan_forwardA",
                    SUZERAIN_EFAILED);
        }
    }

    return p;
}

void
underling_plan_destroy(
        underling_plan plan)
{
    if (plan) {
        if (plan->plan_backwardA) {
            fftw_destroy_plan(plan->plan_backwardA);
            plan->plan_backwardA = NULL;
        }
        if (plan->plan_backwardB) {
            fftw_destroy_plan(plan->plan_backwardB);
            plan->plan_backwardB = NULL;
        }
        if (plan->plan_forwardA) {
            fftw_destroy_plan(plan->plan_forwardA);
            plan->plan_forwardA = NULL;
        }
        if (plan->plan_forwardB) {
            fftw_destroy_plan(plan->plan_forwardB);
            plan->plan_forwardB = NULL;
        }
        plan->data = NULL;
        free(plan);
    }
}

int
underling_execute_long_n2_to_long_n1(
        const underling_plan plan)
{
    if (SUZERAIN_UNLIKELY(plan == NULL)) {
        SUZERAIN_ERROR("plan == NULL", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(plan->plan_backwardA == NULL)) {
        SUZERAIN_ERROR("plan->plan_backwardA == NULL", SUZERAIN_EINVAL);
    }

    fftw_execute(plan->plan_backwardA);

    return SUZERAIN_SUCCESS;
}

int
underling_execute_long_n1_to_long_n0(
        const underling_plan plan)
{
    if (SUZERAIN_UNLIKELY(plan == NULL)) {
        SUZERAIN_ERROR("plan == NULL", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(plan->plan_backwardB == NULL)) {
        SUZERAIN_ERROR("plan->plan_backwardB == NULL", SUZERAIN_EINVAL);
    }

    fftw_execute(plan->plan_backwardB);

    return SUZERAIN_SUCCESS;
}

int
underling_execute_long_n0_to_long_n1(
        const underling_plan plan)
{
    if (SUZERAIN_UNLIKELY(plan == NULL)) {
        SUZERAIN_ERROR("plan == NULL", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(plan->plan_forwardB == NULL)) {
        SUZERAIN_ERROR("plan->plan_forwardB == NULL", SUZERAIN_EINVAL);
    }

    fftw_execute(plan->plan_forwardB);

    return SUZERAIN_SUCCESS;
}

int
underling_execute_long_n1_to_long_n2(
        const underling_plan plan)
{
    if (SUZERAIN_UNLIKELY(plan == NULL)) {
        SUZERAIN_ERROR("plan == NULL", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(plan->plan_forwardA == NULL)) {
        SUZERAIN_ERROR("plan->plan_forwardA == NULL", SUZERAIN_EINVAL);
    }

    fftw_execute(plan->plan_forwardA);

    return SUZERAIN_SUCCESS;
}

void
underling_fprint_grid(
        const underling_grid grid,
        FILE *output_file)
{
    fprintf(output_file, "{underling_grid:");
    if (!grid) {
        fprintf(output_file, "NULL");
    } else {

        fprintf(output_file,
                "{n0=%d,n1=%d,n2=%d}"
                ",{pA=%d,pB=%d}",
                grid->n[0], grid->n[1], grid->n[2],
                grid->pA, grid->pB);

        char buffer[MPI_MAX_OBJECT_NAME];
        int resultlen = 0;
        if (   MPI_Comm_get_name(grid->g_comm, buffer, &resultlen)
            || resultlen == 0) {
            fprintf(output_file, ",{g_comm=%x", grid->g_comm);
        } else {
            fprintf(output_file, ",{g_comm=%s", buffer);
        }

        resultlen = 0;
        if (   MPI_Comm_get_name(grid->pA_comm, buffer, &resultlen)
            || resultlen == 0) {
            fprintf(output_file, ",pA_comm=%x", grid->pA_comm);
        } else {
            fprintf(output_file, ",pA_comm=%s", buffer);
        }

        resultlen = 0;
        if (   MPI_Comm_get_name(grid->pB_comm, buffer, &resultlen)
            || resultlen == 0) {
            fprintf(output_file, ",pB_comm=%x}", grid->pB_comm);
        } else {
            fprintf(output_file, ",pB_comm=%s}", buffer);
        }

    }
    fprintf(output_file, "}");
}

static
void
underling_fprint_transpose(
        const underling_transpose transpose,
        FILE *output_file)
{
    fprintf(output_file, "{underling_transpose:");
    if (!transpose) {
        fprintf(output_file, "NULL");
    } else {
        fprintf(output_file,
                "{d0=%ld,d1=%ld},{block0=%ld,block1=%ld},",
                transpose->d[0], transpose->d[1],
                transpose->block[0], transpose->block[1]);

        char buffer[MPI_MAX_OBJECT_NAME];
        int resultlen = 0;
        if (   MPI_Comm_get_name(transpose->comm, buffer, &resultlen)
            || resultlen == 0) {
            fprintf(output_file, "{comm=%x", transpose->comm);
        } else {
            fprintf(output_file, "{comm=%s", buffer);
        }
        fprintf(output_file,
                ",flags=%u,local_size=%ld},",
                transpose->flags,transpose->local_size);

        fprintf(output_file,
                "{local0=%ld,local_start0=%ld},",
                transpose->local[0],transpose->local_start[0]);
        fprintf(output_file,
                "{local1=%ld,local_start1=%ld}",
                transpose->local[1],transpose->local_start[1]);
    }
    fprintf(output_file, "}");
}

void
underling_fprint_problem(
        const underling_problem problem,
        FILE *output_file)
{
    fprintf(output_file, "{underling_problem:");
    if (!problem) {
        fprintf(output_file, "NULL");
    } else {
        fprintf(output_file,"{howmany=%d,local_memory=%ld}",
                problem->howmany, problem->local_memory);
        for (int i = 2; i >= 0; --i) {
            fprintf(output_file,"{long_n%d:", i);
            underling_fprint_extents(&problem->long_n[i], output_file);
            fprintf(output_file, "}");
        }
        fprintf(output_file,"{backwardA:");
        underling_fprint_transpose(
                problem->backwardA, output_file);
        fprintf(output_file, "}");
        fprintf(output_file,"{backwardB:");
        underling_fprint_transpose(
                problem->backwardB, output_file);
        fprintf(output_file, "}");
        fprintf(output_file,"{forwardB:");
        underling_fprint_transpose(
                problem->forwardB, output_file);
        fprintf(output_file, "}");
        fprintf(output_file,"{forwardA:");
        underling_fprint_transpose(
                problem->forwardA, output_file);
        fprintf(output_file, "}");
    }
    fprintf(output_file, "}");
}

void
underling_fprint_plan(
        const underling_plan plan, FILE *output_file)
{
    fprintf(output_file, "{underling_plan:");
    if (!plan) {
        fprintf(output_file, "NULL");
    } else {

        if (plan->plan_backwardA) {
            fprintf(output_file, "{plan_backwardA:");
            fftw_fprint_plan(plan->plan_backwardA, output_file);
            fprintf(output_file, "}");
        }
        if (plan->plan_backwardB) {
            fprintf(output_file, "{plan_backwardB:");
            fftw_fprint_plan(plan->plan_backwardB, output_file);
            fprintf(output_file, "}");
        }

        if (plan->plan_forwardA) {
            fprintf(output_file, "{plan_forwardB:");
            fftw_fprint_plan(plan->plan_forwardB, output_file);
            fprintf(output_file, "}");
        }
        if (plan->plan_forwardB) {
            fprintf(output_file, "{plan_forwardA:");
            fftw_fprint_plan(plan->plan_forwardA, output_file);
            fprintf(output_file, "}");
        }
    }
}

void
underling_fprint_extents(
        const underling_extents *extents,
        FILE *output_file)
{
    if (!extents) {
        fprintf(output_file, "NULL");
    } else {
        fprintf(output_file, "extents=[%d,%d)x[%d,%d)x[%d,%d)",
                extents->start[0],
                extents->start[0] + extents->size[0],
                extents->start[1],
                extents->start[1] + extents->size[1],
                extents->start[2],
                extents->start[2] + extents->size[2]);
        fprintf(output_file, ",strides={%d,%d,%d}",
                extents->stride[0],
                extents->stride[1],
                extents->stride[2]);
    }
}
