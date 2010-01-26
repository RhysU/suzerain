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

// *******************************************************************
// INTERNAL STRUCTS INTERNAL STRUCTS INTERNAL STRUCTS INTERNAL STRUCTS
// *******************************************************************

struct underling_grid_s {
    int n0;
    int n1;
    int n2;
    int pA;
    int pB;
    MPI_Comm g_comm;
    int g_rank;
    int g_coords[2];
    MPI_Comm pA_comm;
    int pA_rank;
    MPI_Comm pB_comm;
    int pB_rank;
};

typedef struct underling_transpose_s * underling_transpose; // Internal!
struct underling_transpose_s {
    ptrdiff_t d0;
    ptrdiff_t d1;
    ptrdiff_t howmany;
    ptrdiff_t block0;
    ptrdiff_t block1;
    MPI_Comm comm;                    // underling_grid owns resource
    unsigned flags;
    ptrdiff_t local_d0;
    ptrdiff_t local_d0_start;
    ptrdiff_t local_d1;
    ptrdiff_t local_d1_start;
    ptrdiff_t local_size;
};

typedef struct underling_extents { // Internal!
    int local_n0;
    int local_n0_start;
    int local_n1;
    int local_n1_start;
    int local_n2;
    int local_n2_start;
} underling_extents;

struct underling_problem_s {
    underling_grid grid;              // grid owns its resources
    int howmany;                      // # of real values to transpose
    underling_extents long_n2;        // Layout details for n2 long
    underling_extents long_n1;        // Layout details for n1 long
    underling_extents long_n0;        // Layout details for n0 long
    underling_transpose backwardA;    // n2 long to n1 long
    underling_transpose backwardB;    // n1 long to n0 long
    underling_transpose forwardB;     // n0 long to n1 long
    underling_transpose forwardA;     // n1 long to n2 long
    ptrdiff_t local_size;             // Max of all local sizes
};

struct underling_plan_s {
    underling_problem problem;         // problem owns its resources
    fftw_plan plan_backwardA;
    fftw_plan plan_backwardB;
    fftw_plan plan_forwardA;
    fftw_plan plan_forwardB;
    underling_real *data;              // API end-user owns resources
};

// ********************************************************************
// INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL
// ********************************************************************

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
        unsigned flags);

void
underling_fprint_transpose(
        const underling_transpose transpose,
        FILE *output_file);

void
underling_fprint_extents(
        const underling_extents *extents,
        FILE *output_file);

// **************************************************************************
// IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION
// **************************************************************************

underling_grid
underling_grid_create(
        MPI_Comm comm,
        int n0,
        int n1,
        int n2,
        int pA,
        int pB)
{
    // Sanity check incoming, non-MPI arguments
    if (n0 < 1) {
        SUZERAIN_ERROR_NULL("n0 >= 1 required", SUZERAIN_EINVAL);
    }
    if (n1 < 1) {
        SUZERAIN_ERROR_NULL("n1 >= 1 required", SUZERAIN_EINVAL);
    }
    if (n2 < 1) {
        SUZERAIN_ERROR_NULL("n2 >= 1 required", SUZERAIN_EINVAL);
    }
    if (pA < 0) {
        SUZERAIN_ERROR_NULL("pA >= 0 required", SUZERAIN_EINVAL);
    }
    if (pB < 0) {
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
        if (pA * pB != nproc) {
            char reason[127];
            snprintf(reason, sizeof(reason)/sizeof(reason[0]),
                    "Invalid processor grid: pA {%d} * pB {%d} != nproc {%d}",
                    pA, pB, nproc);
            SUZERAIN_ERROR_NULL(reason, SUZERAIN_EFAILED);
        }
    }

    // Clone the communicator and create the 2D Cartesian topology
    MPI_Comm g_comm;
    {
        int dims[2]     = { pA, pB };
        int periodic[2] = { 0, 0 };
        SUZERAIN_MPICHKN(MPI_Cart_create(
                comm, 2, dims, periodic, 1/*reordering allowed*/, &g_comm));
    }
    // Cache the rank and coordinates of this process within g_comm
    int g_rank;
    SUZERAIN_MPICHKN(MPI_Comm_rank(g_comm, &g_rank));
    int g_coords[2];
    SUZERAIN_MPICHKN(MPI_Cart_coords(g_comm, g_rank, 2, g_coords));

    // Create communicator for the PA direction
    MPI_Comm pA_comm;
    {
        int remain_dims[2] = { 1, 0 };
        SUZERAIN_MPICHKN(MPI_Cart_sub(g_comm, remain_dims, &pA_comm));
    }
    // Find the rank of this process within g_comm
    int pA_rank;
    SUZERAIN_MPICHKN(MPI_Comm_rank(pA_comm, &pA_rank));

    // Create communicator for the PB direction
    MPI_Comm pB_comm;
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
    if (g == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate space for grid",
                             SUZERAIN_ENOMEM);
    }
    // Copy the grid parameters to the grid workspace
    g->n0          = n0;
    g->n1          = n1;
    g->n2          = n2;
    g->pA          = pA;
    g->pB          = pB;
    g->g_comm      = g_comm;
    g->g_rank      = g_rank;
    g->g_coords[0] = g_coords[0];
    g->g_coords[1] = g_coords[1];
    g->pA_comm     = pA_comm;
    g->pA_rank     = pA_rank;
    g->pB_comm     = pB_comm;
    g->pB_rank     = pB_rank;

    return g;
}

void
underling_grid_destroy(underling_grid grid)
{
    if (grid) {
        if (grid->g_comm) {
            SUZERAIN_MPICHKV(MPI_Comm_disconnect(&grid->g_comm));
            grid->g_comm = MPI_COMM_NULL;
        }
        if (grid->pA_comm) {
            SUZERAIN_MPICHKV(MPI_Comm_disconnect(&grid->pA_comm));
            grid->pA_comm = MPI_COMM_NULL;
        }
        if (grid->pB_comm) {
            SUZERAIN_MPICHKV(MPI_Comm_disconnect(&grid->pB_comm));
            grid->pB_comm = MPI_COMM_NULL;
        }
        free(grid);
    }
}


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
    // Create and initialize the transpose workspace
    underling_transpose t = calloc(1, sizeof(struct underling_transpose_s));
    if (t == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate space for transpose",
                             SUZERAIN_ENOMEM);
    }

    // Fix struct values known from arguments
    t->d0      = d0;
    t->d1      = d1;
    t->howmany = howmany;
    t->block0  = block0;
    t->block1  = block1;
    t->comm    = comm;
    t->flags   = flags;

    // Fix struct details obtainable via FFTW MPI call
    const ptrdiff_t d[2] = { t->d0, t->d1 };
    t->local_size = fftw_mpi_local_size_many_transposed(
                                           /*rank*/sizeof(d)/sizeof(d[0]),
                                            d,
                                            t->howmany,
                                            t->block0,
                                            t->block1,
                                            t->comm,
                                            &(t->local_d0),
                                            &(t->local_d0_start),
                                            &(t->local_d1),
                                            &(t->local_d1_start));

    return t;
}

underling_transpose
underling_transpose_create_inverse(
        const underling_transpose forward)
{
    // TODO Handle FFTW_MPI_TRANSPOSED_IN, FFTW_MPI_TRANSPOSED_OUT correctly
    if (forward->flags) {
        SUZERAIN_ERROR_NULL("Nontrivial transpose flags not yet implemented!",
                SUZERAIN_ESANITY);
    }

    // Create and initialize the transpose workspace
    underling_transpose backward
        = calloc(1, sizeof(struct underling_transpose_s));
    if (backward == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate space for transpose",
                             SUZERAIN_ENOMEM);
    }

    // Fix struct values known from inverting forward plan
    backward->d0      = forward->d1;
    backward->d1      = forward->d0;
    backward->howmany = forward->howmany;
    backward->block0  = forward->block1;
    backward->block1  = forward->block0;
    backward->comm    = forward->comm;
    backward->flags   = forward->flags;
    // Unsure of how block sizes should be flipped for inversion

    // Fix struct details obtainable via FFTW MPI call
    const ptrdiff_t d[2] = { backward->d0, backward->d1 };
    backward->local_size = fftw_mpi_local_size_many_transposed(
                                           /*rank*/sizeof(d)/sizeof(d[0]),
                                            d,
                                            backward->howmany,
                                            backward->block0,
                                            backward->block1,
                                            backward->comm,
                                            &(backward->local_d0),
                                            &(backward->local_d0_start),
                                            &(backward->local_d1),
                                            &(backward->local_d1_start));

    return backward;
}

void
underling_transpose_destroy(
        underling_transpose transpose)
{
    if (transpose) {
        free(transpose);
    }
}

fftw_plan
underling_transpose_fftw_plan(
        const underling_transpose transpose,
        underling_real *in,
        underling_real *out,
        unsigned flags)
{
    return fftw_mpi_plan_many_transpose(
            transpose->d0,
            transpose->d1,
            transpose->howmany,
            transpose->block0,
            transpose->block1,
            in,
            out,
            transpose->comm,
            transpose->flags | flags);
}

underling_problem
underling_problem_create(
        underling_grid grid,
        int howmany)
{
    if (grid == NULL) {
        SUZERAIN_ERROR_NULL("non-NULL grid required", SUZERAIN_EINVAL);
    }
    if (howmany < 1) {
        SUZERAIN_ERROR_NULL("howmany >= 1 required", SUZERAIN_EINVAL);
    }

    // Create and initialize the problem workspace
    underling_problem p = calloc(1, sizeof(struct underling_problem_s));
    if (p == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate space for problem",
                             SUZERAIN_ENOMEM);
    }
    // Copy the problem parameters to the problem workspace
    p->grid    = grid;
    p->howmany = howmany;

    // Global pencil decomposition details
    // -------------------------------------------------------
    // Long in n2:                        (n0/pB x n1/pA) x n2
    // Long in n1: n2/pA x (n0/pB x n1) = (n2/pA x n0/pB) x n1
    // Long in n0: n1/pB x (n2/pA x n0) = (n1/pB x n2/pA) x n0

    // Fix {n2,n1,n0} dimension details in p->long_{n2,n1,n0}
    p->long_n2.local_n2       = grid->n2;
    p->long_n2.local_n2_start = 0;
    p->long_n1.local_n1       = grid->n1;
    p->long_n1.local_n1_start = 0;
    p->long_n0.local_n0       = grid->n0;
    p->long_n0.local_n0_start = 0;

    // Decompose {n0,n1}/pB and store details in p->long_{(n2,n1),n0}
    {
        ptrdiff_t local_d0, local_d0_start, local_d1, local_d1_start;
        ptrdiff_t dB[2] = {grid->n0, grid->n1}; // Never performed
        fftw_mpi_local_size_many_transposed(2, dB, p->howmany,
                FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, grid->pB_comm,
                &local_d0, &local_d0_start, &local_d1, & local_d1_start);
        assert(local_d0       <= INT_MAX);
        assert(local_d0_start <= INT_MAX);
        assert(local_d1       <= INT_MAX);
        assert(local_d1_start <= INT_MAX);
        p->long_n2.local_n0       = local_d0;
        p->long_n2.local_n0_start = local_d0_start;
        p->long_n1.local_n0       = local_d0;
        p->long_n1.local_n0_start = local_d0_start;
        p->long_n0.local_n1       = local_d1;
        p->long_n0.local_n1_start = local_d1_start;
    }

    // Decompose {n1,n2}/pA and store details in p->long_{n2,(n1,n0)}
    {
        ptrdiff_t local_d0, local_d0_start, local_d1, local_d1_start;
        ptrdiff_t dA[2] = {grid->n1, grid->n2}; // Never performed
        fftw_mpi_local_size_many_transposed(2, dA, p->howmany,
                FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, grid->pA_comm,
                &local_d0, &local_d0_start, &local_d1, &local_d1_start);
        assert(local_d0       <= INT_MAX);
        assert(local_d0_start <= INT_MAX);
        assert(local_d1       <= INT_MAX);
        assert(local_d1_start <= INT_MAX);
        p->long_n2.local_n1       = local_d0;
        p->long_n2.local_n1_start = local_d0_start;
        p->long_n1.local_n2       = local_d1;
        p->long_n1.local_n2_start = local_d1_start;
        p->long_n0.local_n2       = local_d1;
        p->long_n0.local_n2_start = local_d1_start;
    }

    // Transpose pA details: (n0/pB x n1/pA) x n2 to n2/pA x (n0/pB x n1)
    const ptrdiff_t pA_d[2] = { p->long_n2.local_n0 * grid->n1,
                                grid->n2 };
    ptrdiff_t pA_block[2]   = { p->long_n2.local_n0 * p->long_n2.local_n1,
                                p->long_n1.local_n2 };
    SUZERAIN_MPICHKN(MPI_Bcast(pA_block, 2, MPI_LONG, 0, grid->pA_comm));

    // Transpose pB details: (n2/pA x n0/pB) x n1 to n1/pB x (n2/pA x n0)
    const ptrdiff_t pB_d[2] = { p->long_n1.local_n2 * grid->n0,
                                grid->n1 };
    ptrdiff_t pB_block[2]   = { p->long_n1.local_n2 * p->long_n1.local_n0,
                                p->long_n0.local_n1 };
    SUZERAIN_MPICHKN(MPI_Bcast(pB_block, 2, MPI_LONG, 0, grid->pB_comm));

    // Wave towards physical MPI transpose: long in n2 to long in n1
    p->backwardA = underling_transpose_create(pA_d[0],
                                            pA_d[1],
                                            p->howmany,
                                            pA_block[0],
                                            pA_block[1],
                                            grid->pA_comm,
                                            /*flags*/0);
    if (p->backwardA == NULL) {
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
    if (p->backwardB == NULL) {
        underling_problem_destroy(p);
        SUZERAIN_ERROR_NULL("failed creating p->backwardB",
                SUZERAIN_EFAILED);
    }

    // Physical towards wave MPI transpose: long in n0 to long in n1
    p->forwardB = underling_transpose_create_inverse(p->backwardB);
    if (p->forwardB == NULL) {
        underling_problem_destroy(p);
        SUZERAIN_ERROR_NULL("failed creating p->forwardB",
                SUZERAIN_EFAILED);
    }

    // Physical towards wave MPI transpose: long in n1 to long in n2
    p->forwardA = underling_transpose_create_inverse(p->backwardA);
    if (p->forwardA == NULL) {
        underling_problem_destroy(p);
        SUZERAIN_ERROR_NULL("failed creating p->forwardA",
                SUZERAIN_EFAILED);
    }

    // p->local_size is overall maximum of all local_size values
    p->local_size = p->backwardA->local_size;
    if (p->local_size < p->backwardB->local_size) {
        p->local_size = p->backwardB->local_size;
    }
    if (p->local_size < p->forwardB->local_size) {
        p->local_size = p->forwardB->local_size;
    }
    if (p->local_size < p->forwardA->local_size) {
        p->local_size = p->forwardA->local_size;
    }

    return p;
}

size_t
underling_local_size(
        const underling_problem problem)
{
    return problem->local_size;
}

size_t
underling_optimum_local_size(
        const underling_problem problem)
{
    const size_t global_data =   problem->grid->n0
                               * problem->grid->n1
                               * problem->grid->n2
                               * problem->howmany;

    const size_t nprocessors =   problem->grid->pA
                               * problem->grid->pB;

    return global_data/nprocessors;
}

size_t
underling_local_long_n2(
        const underling_problem problem,
        int *start,
        int *size,
        int *stride)
{
    if (problem == NULL) {
        SUZERAIN_ERROR_VAL("problem == NULL", SUZERAIN_EINVAL, 0);
    }

    const underling_extents * const e = &problem->long_n2;
    if (start) {
        start[0] = e->local_n0_start;
        start[1] = e->local_n1_start;
        start[2] = e->local_n2_start;
    }
    if (size) {
        size[0] = e->local_n0;
        size[1] = e->local_n1;
        size[2] = e->local_n2;
    }
    // Long in n2: (n0/pB x  n1/pA) x n2
    // Row-major storage
    if (stride) {
        stride[2] = problem->howmany;
        stride[1] = stride[2] * e->local_n2;
        stride[0] = stride[1] * e->local_n1;
    }

    return problem->howmany * e->local_n0 * e->local_n1 * e->local_n2;
}

size_t
underling_local_long_n1(
        const underling_problem problem,
        int *start,
        int *size,
        int *stride)
{
    if (problem == NULL) {
        SUZERAIN_ERROR_VAL("problem == NULL", SUZERAIN_EINVAL, 0);
    }

    const underling_extents * const e = &problem->long_n1;
    if (start) {
        start[0] = e->local_n0_start;
        start[1] = e->local_n1_start;
        start[2] = e->local_n2_start;
    }
    if (size) {
        size[0] = e->local_n0;
        size[1] = e->local_n1;
        size[2] = e->local_n2;
    }
    // Long in n1: n2/pA x (n0/pB x n1) = (n2/pA x n0/pB) x n1
    // Row-major storage
    if (stride) {
        stride[1] = problem->howmany;
        stride[0] = stride[1] * e->local_n1;
        stride[2] = stride[0] * e->local_n0;
    }

    return problem->howmany * e->local_n0 * e->local_n1 * e->local_n2;
}

size_t
underling_local_long_n0(
        const underling_problem problem,
        int *start,
        int *size,
        int *stride)
{
    if (problem == NULL) {
        SUZERAIN_ERROR_VAL("problem == NULL", SUZERAIN_EINVAL, 0);
    }

    const underling_extents * const e = &problem->long_n0;
    if (start) {
        start[0] = e->local_n0_start;
        start[1] = e->local_n1_start;
        start[2] = e->local_n2_start;
    }
    if (size) {
        size[0] = e->local_n0;
        size[1] = e->local_n1;
        size[2] = e->local_n2;
    }
    // Long in n0: n1/pB x (n2/pA x n0) = ( n1/pB x  n2/pA) x n0
    // Row-major storage
    if (stride) {
        stride[0] = problem->howmany;
        stride[2] = stride[0] * e->local_n0;
        stride[1] = stride[2] * e->local_n2;
    }

    return problem->howmany * e->local_n0 * e->local_n1 * e->local_n2;
}

void
underling_problem_destroy(
        underling_problem problem)
{
    if (problem) {
        problem->grid = NULL;
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
        underling_problem problem,
        underling_real * data,
        unsigned direction_flags,
        unsigned rigor_flags)
{
    if (problem == NULL) {
        SUZERAIN_ERROR_NULL("non-NULL problem required", SUZERAIN_EINVAL);
    }
    if (data == NULL) {
        SUZERAIN_ERROR_NULL("non-NULL data required", SUZERAIN_EINVAL);
    }
    const unsigned non_direction_mask =   ~UNDERLING_DIRECTION_FORWARD
                                        & ~UNDERLING_DIRECTION_BACKWARD;
    if (direction_flags & non_direction_mask) {
        SUZERAIN_ERROR_NULL(
            "direction_flags contains ~UNDERLING_DIRECTION_{FORWARD,BACKWARD}",
            SUZERAIN_EINVAL);
    }
    const unsigned non_rigor_mask =   ~FFTW_ESTIMATE
                                    & ~FFTW_MEASURE
                                    & ~FFTW_PATIENT
                                    & ~FFTW_EXHAUSTIVE
                                    & ~FFTW_WISDOM_ONLY;
    if (rigor_flags & non_rigor_mask) {
        SUZERAIN_ERROR_NULL("FFTW non-rigor bits disallowed", SUZERAIN_EINVAL);
    }

    // Perform both directions if neither specified
    if (direction_flags == 0) {
        direction_flags =   UNDERLING_DIRECTION_FORWARD
                          | UNDERLING_DIRECTION_BACKWARD;
    }

    // Create and initialize the plan workspace
    underling_plan p = calloc(1, sizeof(struct underling_plan_s));
    if (p == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate space for grid",
                             SUZERAIN_ENOMEM);
    }
    p->problem = problem;
    p->data = data;

    // Create FFTW MPI plans to transpose from long in n2 to long in 20
    if (direction_flags | UNDERLING_DIRECTION_BACKWARD) {
        p->plan_backwardA = underling_transpose_fftw_plan(
                p->problem->backwardA, p->data, p->data, rigor_flags);
        if (p->plan_backwardA == NULL) {
            underling_plan_destroy(p);
            SUZERAIN_ERROR_NULL(
                    "FFTW MPI returned NULL plan: plan_backwardA",
                    SUZERAIN_EFAILED);
        }

        p->plan_backwardB = underling_transpose_fftw_plan(
                p->problem->backwardB, p->data, p->data, rigor_flags);
        if (p->plan_backwardB == NULL) {
            underling_plan_destroy(p);
            SUZERAIN_ERROR_NULL(
                    "FFTW MPI returned NULL plan: plan_backwardB",
                    SUZERAIN_EFAILED);
        }

    }

    // Create FFTW MPI plans to transpose from long in n0 to long in n2
    if (direction_flags | UNDERLING_DIRECTION_FORWARD) {
        p->plan_forwardB = underling_transpose_fftw_plan(
                p->problem->forwardB, p->data, p->data, rigor_flags);
        if (p->plan_forwardB == NULL) {
            underling_plan_destroy(p);
            SUZERAIN_ERROR_NULL(
                    "FFTW MPI returned NULL plan: plan_forwardB",
                    SUZERAIN_EFAILED);
        }

        p->plan_forwardA = underling_transpose_fftw_plan(
                p->problem->forwardA, p->data, p->data, rigor_flags);
        if (p->plan_forwardA == NULL) {
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
        plan->problem = NULL;
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
    }
}

int
underling_execute_backward(
        underling_plan plan)
{
    if (plan == NULL) {
        SUZERAIN_ERROR("non-NULL plan required", SUZERAIN_EINVAL);
    }

    const int valid =    plan->plan_backwardA
                      && plan->plan_backwardB;
    if (!valid) {
        SUZERAIN_ERROR("plan has one or more NULL subplans",
                SUZERAIN_EINVAL);
    }

    fftw_execute(plan->plan_backwardA);
    fftw_execute(plan->plan_backwardB);

    return SUZERAIN_SUCCESS;
}

int
underling_execute_forward(
        underling_plan plan)
{
    if (plan == NULL) {
        SUZERAIN_ERROR("non-NULL plan required", SUZERAIN_EINVAL);
    }
    const int valid =    plan->plan_forwardB
                      && plan->plan_forwardA;
    if (!valid) {
        SUZERAIN_ERROR("plan has one or more NULL subplans",
                SUZERAIN_EINVAL);
    }

    fftw_execute(plan->plan_forwardB);
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
                grid->n0, grid->n1, grid->n2,
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
                transpose->d0, transpose->d1,
                transpose->block0, transpose->block1);
        fprintf(output_file,
                "{comm=%x,flags=%u,local_size=%ld},",
                transpose->comm,transpose->flags,transpose->local_size);
        fprintf(output_file,
                "{local_d0=%ld,local_d0_start=%ld},",
                transpose->local_d0,transpose->local_d0_start);
        fprintf(output_file,
                "{local_d1=%ld,local_d1_start=%ld}",
                transpose->local_d1,transpose->local_d1_start);
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
        fprintf(output_file,"{howmany=%d,local_size=%ld}",
                problem->howmany, problem->local_size);
        fprintf(output_file,"{long_n2:");
        underling_fprint_extents(&(problem->long_n2), output_file);
        fprintf(output_file, "}");
        fprintf(output_file,"{long_n1:");
        underling_fprint_extents(&(problem->long_n1), output_file);
        fprintf(output_file, "}");
        fprintf(output_file,"{long_n0:");
        underling_fprint_extents(&(problem->long_n0), output_file);
        fprintf(output_file, "}");
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
        fprintf(output_file, "[%d,%d)x[%d,%d)x[%d,%d)",
                extents->local_n0_start,
                extents->local_n0_start + extents->local_n0,
                extents->local_n1_start,
                extents->local_n1_start + extents->local_n1,
                extents->local_n2_start,
                extents->local_n2_start + extents->local_n2);
    }
}
