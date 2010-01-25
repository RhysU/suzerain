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
    int np0;
    int nw0;
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
    int local_nw0;
    int local_nw0_start;
    int local_n1;
    int local_n1_start;
    int local_n2;
    int local_n2_start;
} underling_extents;

struct underling_problem_s {
    underling_grid grid;              // grid owns its resources
    int nfields;                      // # of complex-valued state...
    int howmany;                      // ...as a # of real values
    underling_extents long_n2;        // Layout details for n2 long
    underling_extents long_n1;        // Layout details for n1 long
    underling_extents long_n0;        // Layout details for n0 long
    underling_transpose tophysA;      // n2 long to n1 long
    underling_transpose tophysB;      // n1 long to n0 long
    underling_transpose towaveB;      // n0 long to n1 long
    underling_transpose towaveA;      // n1 long to n2 long
    ptrdiff_t local_size;             // Max of all local sizes
};

struct underling_plan_s {
    underling_problem problem;         // problem owns its resources
    fftw_plan transpose_tophysA;
    fftw_plan c2c_tophysical_n1;
    fftw_plan transpose_tophysB;
    fftw_plan c2r_tophysical_n0;
    fftw_plan r2c_towave_n0;
    fftw_plan transpose_towaveA;
    fftw_plan c2c_towave_n1;
    fftw_plan transpose_towaveB;
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
        int np0,
        int np1,
        int np2,
        int pA,
        int pB)
{
    // Sanity check incoming, non-MPI arguments
    if (np0 < 1) {
        SUZERAIN_ERROR_NULL("np0 >= 1 required", SUZERAIN_EINVAL);
    }
    if (np1 < 1) {
        SUZERAIN_ERROR_NULL("np1 >= 1 required", SUZERAIN_EINVAL);
    }
    if (np2 < 1) {
        SUZERAIN_ERROR_NULL("np2 >= 1 required", SUZERAIN_EINVAL);
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

    // Determine the wave space dimensions in the n0 where
    // n0 is only about half as long in wave space.
    const int nw0 = np0/2 + 1;
    const int n1  = np1;
    const int n2  = np2;

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
    g->np0         = np0;
    g->nw0         = nw0;
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
        int nfields)
{
    if (grid == NULL) {
        SUZERAIN_ERROR_NULL("non-NULL grid required", SUZERAIN_EINVAL);
    }
    if (nfields < 1) {
        SUZERAIN_ERROR_NULL("nfields >= 1 required", SUZERAIN_EINVAL);
    }

    // Create and initialize the problem workspace
    underling_problem p = calloc(1, sizeof(struct underling_problem_s));
    if (p == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate space for problem",
                             SUZERAIN_ENOMEM);
    }
    // Copy the problem parameters to the problem workspace
    p->grid    = grid;
    p->nfields = nfields;
    p->howmany = nfields * sizeof(underling_complex)/sizeof(underling_real);

    // Global pencil decomposition details
    // ---------------------------------------------------------
    // Long in n2:                         (nw0/pB x  n1/pA) x n2
    // Long in n1: n2/pA x (nw0/pB x n1) = ( n2/pA x nw0/pB) x n1
    // Long in n0: n1/pB x (n2/pA x nw0) = ( n1/pB x  n2/pA) x nw0

    // Fix {n2,n1,n0} dimension details in p->long_{n2,n1,n0}
    p->long_n2.local_n2        = grid->n2;
    p->long_n2.local_n2_start  = 0;
    p->long_n1.local_n1        = grid->n1;
    p->long_n1.local_n1_start  = 0;
    p->long_n0.local_nw0       = grid->nw0;
    p->long_n0.local_nw0_start = 0;

    // Decompose {nw0,n1}/pB and store details in p->long_{(n2,n1),n0}
    {
        ptrdiff_t local_d0, local_d0_start, local_d1, local_d1_start;
        ptrdiff_t dB[2] = {grid->nw0, grid->n1}; // Never performed
        fftw_mpi_local_size_many_transposed(2, dB, p->howmany,
                FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, grid->pB_comm,
                &local_d0, &local_d0_start, &local_d1, & local_d1_start);
        assert(local_d0       <= INT_MAX);
        assert(local_d0_start <= INT_MAX);
        assert(local_d1       <= INT_MAX);
        assert(local_d1_start <= INT_MAX);
        p->long_n2.local_nw0       = local_d0;
        p->long_n2.local_nw0_start = local_d0_start;
        p->long_n1.local_nw0       = local_d0;
        p->long_n1.local_nw0_start = local_d0_start;
        p->long_n0.local_n1        = local_d1;
        p->long_n0.local_n1_start  = local_d1_start;
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

    // Transpose pA details: (nw0/pB x n1/pA) x n2 to n2/pA x (nw0/pB x n1)
    const ptrdiff_t pA_d[2] = { p->long_n2.local_nw0 * grid->n1,
                                grid->n2 };
    ptrdiff_t pA_block[2]   = { p->long_n2.local_nw0 * p->long_n2.local_n1,
                                p->long_n1.local_n2 };
    SUZERAIN_MPICHKN(MPI_Bcast(pA_block, 2, MPI_LONG, 0, grid->pA_comm));

    // Transpose pB details: (n2/pA x nw0/pB) x n1 to n1/pB x (n2/pA x nw0)
    const ptrdiff_t pB_d[2] = { p->long_n1.local_n2 * grid->nw0,
                                grid->n1 };
    ptrdiff_t pB_block[2]   = { p->long_n1.local_n2 * p->long_n1.local_nw0,
                                p->long_n0.local_n1 };
    SUZERAIN_MPICHKN(MPI_Bcast(pB_block, 2, MPI_LONG, 0, grid->pB_comm));

    // Wave towards physical MPI transpose: long in n2 to long in n1
    p->tophysA = underling_transpose_create(pA_d[0],
                                            pA_d[1],
                                            p->howmany,
                                            pA_block[0],
                                            pA_block[1],
                                            grid->pA_comm,
                                            /*flags*/0);
    if (p->tophysA == NULL) {
        underling_problem_destroy(p);
        SUZERAIN_ERROR_NULL("failed creating p->tophysA",
                SUZERAIN_EFAILED);
    }

    // Wave towards physical MPI transpose: long in n1 to long in nw0
    p->tophysB = underling_transpose_create(pB_d[0],
                                            pB_d[1],
                                            p->howmany,
                                            pB_block[0],
                                            pB_block[1],
                                            grid->pB_comm,
                                            /*flags*/0);
    if (p->tophysB == NULL) {
        underling_problem_destroy(p);
        SUZERAIN_ERROR_NULL("failed creating p->tophysB",
                SUZERAIN_EFAILED);
    }

    // Physical towards wave MPI transpose: long in n0 to long in n1
    p->towaveB = underling_transpose_create_inverse(p->tophysB);
    if (p->towaveB == NULL) {
        underling_problem_destroy(p);
        SUZERAIN_ERROR_NULL("failed creating p->towaveB",
                SUZERAIN_EFAILED);
    }

    // Physical towards wave MPI transpose: long in n1 to long in n2
    p->towaveA = underling_transpose_create_inverse(p->tophysA);
    if (p->towaveA == NULL) {
        underling_problem_destroy(p);
        SUZERAIN_ERROR_NULL("failed creating p->towaveA",
                SUZERAIN_EFAILED);
    }

    // p->local_size is overall maximum of all local_size values
    p->local_size = p->tophysA->local_size;
    if (p->local_size < p->tophysB->local_size) {
        p->local_size = p->tophysB->local_size;
    }
    if (p->local_size < p->towaveB->local_size) {
        p->local_size = p->towaveB->local_size;
    }
    if (p->local_size < p->towaveA->local_size) {
        p->local_size = p->towaveA->local_size;
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
    const size_t global_data =   problem->grid->nw0
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
        int *stride_complex)
{
    if (problem == NULL) {
        SUZERAIN_ERROR_VAL("problem == NULL", SUZERAIN_EINVAL, 0);
    }

    const underling_extents * const e = &problem->long_n2;
    if (start) {
        start[0] = e->local_nw0_start;
        start[1] = e->local_n1_start;
        start[2] = e->local_n2_start;
    }
    if (size) {
        size[0] = e->local_nw0;
        size[1] = e->local_n1;
        size[2] = e->local_n2;
    }
    // Long in n2: (nw0/pB x  n1/pA) x n2
    // Row-major storage
    if (stride_complex) {
        stride_complex[2] = problem->nfields;
        stride_complex[1] = stride_complex[2] * e->local_n2;
        stride_complex[0] = stride_complex[1] * e->local_n1;
    }

    return problem->nfields * e->local_nw0 * e->local_n1 * e->local_n2;
}

size_t
underling_local_long_n1(
        const underling_problem problem,
        int *start,
        int *size,
        int *stride_complex)
{
    if (problem == NULL) {
        SUZERAIN_ERROR_VAL("problem == NULL", SUZERAIN_EINVAL, 0);
    }

    const underling_extents * const e = &problem->long_n1;
    if (start) {
        start[0] = e->local_nw0_start;
        start[1] = e->local_n1_start;
        start[2] = e->local_n2_start;
    }
    if (size) {
        size[0] = e->local_nw0;
        size[1] = e->local_n1;
        size[2] = e->local_n2;
    }
    // Long in n1: n2/pA x (nw0/pB x n1) = ( n2/pA x nw0/pB) x n1
    // Row-major storage
    if (stride_complex) {
        stride_complex[1] = problem->nfields;
        stride_complex[0] = stride_complex[1] * e->local_n1;
        stride_complex[2] = stride_complex[0] * e->local_nw0;
    }

    return problem->nfields * e->local_nw0 * e->local_n1 * e->local_n2;
}

size_t
underling_local_long_n0(
        const underling_problem problem,
        int *start,
        int *size,
        int *stride_real)
{
    if (problem == NULL) {
        SUZERAIN_ERROR_VAL("problem == NULL", SUZERAIN_EINVAL, 0);
    }

    const underling_extents * const e = &problem->long_n0;
    if (start) {
        start[0] = e->local_nw0_start;
        start[1] = e->local_n1_start;
        start[2] = e->local_n2_start;
    }
    if (size) {
        size[0] = e->local_nw0;
        size[1] = e->local_n1;
        size[2] = e->local_n2;
    }
    // Long in n0: n1/pB x (n2/pA x nw0) = ( n1/pB x  n2/pA) x nw0
    // Row-major storage
    if (stride_real) {
        stride_real[0] = problem->nfields;
        stride_real[2] = stride_real[0] * e->local_nw0;
        stride_real[1] = stride_real[2] * e->local_n2;
    }

    return problem->nfields * e->local_nw0 * e->local_n1 * e->local_n2;
}

void
underling_problem_destroy(
        underling_problem problem)
{
    if (problem) {
        problem->grid = NULL;
        if (problem->tophysA) {
            underling_transpose_destroy(problem->tophysA);
            problem->tophysA = NULL;
        }
        if (problem->tophysB) {
            underling_transpose_destroy(problem->tophysB);
            problem->tophysB = NULL;
        }
        if (problem->towaveB) {
            underling_transpose_destroy(problem->towaveB);
            problem->towaveB = NULL;
        }
        if (problem->towaveA) {
            underling_transpose_destroy(problem->towaveA);
            problem->towaveA = NULL;
        }
        free(problem);
    }
}

underling_plan
underling_plan_create(
        underling_problem problem,
        underling_real * data,
        int will_perform_c2r,
        int will_perform_r2c,
        unsigned rigor_flags)
{
    if (problem == NULL) {
        SUZERAIN_ERROR_NULL("non-NULL problem required", SUZERAIN_EINVAL);
    }
    if (data == NULL) {
        SUZERAIN_ERROR_NULL("non-NULL data required", SUZERAIN_EINVAL);
    }
    if (!(will_perform_c2r || will_perform_r2c)) {
        SUZERAIN_ERROR_NULL(
                "one or both of will_perform_{c2r,r2c} required",
                SUZERAIN_EINVAL);
    }
    const unsigned non_rigor_flag_mask =   ~FFTW_ESTIMATE
                                         & ~FFTW_MEASURE
                                         & ~FFTW_PATIENT
                                         & ~FFTW_EXHAUSTIVE
                                         & ~FFTW_WISDOM_ONLY;
    if (rigor_flags & non_rigor_flag_mask) {
        SUZERAIN_ERROR_NULL("FFTW non-rigor bits disallowed", SUZERAIN_EINVAL);
    }

    // Create and initialize the plan workspace
    underling_plan p = calloc(1, sizeof(struct underling_plan_s));
    if (p == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate space for grid",
                             SUZERAIN_ENOMEM);
    }
    p->problem = problem;
    p->data = data;

    /* Create FFTW plans to transform from wave to physical space */
    if (will_perform_c2r) {
        p->transpose_tophysA = underling_transpose_fftw_plan(
                p->problem->tophysA, p->data, p->data, rigor_flags);
        if (p->transpose_tophysA == NULL) {
            underling_plan_destroy(p);
            SUZERAIN_ERROR_NULL(
                    "FFTW MPI returned NULL plan: transpose_tophysA",
                    SUZERAIN_EFAILED);
        }

        p->transpose_tophysB = underling_transpose_fftw_plan(
                p->problem->tophysB, p->data, p->data, rigor_flags);
        if (p->transpose_tophysB == NULL) {
            underling_plan_destroy(p);
            SUZERAIN_ERROR_NULL(
                    "FFTW MPI returned NULL plan: transpose_tophysB",
                    SUZERAIN_EFAILED);
        }

        /* TODO Create FFT plan: c2c_tophysical_n1 */
        /* TODO Create FFT plan: c2r_tophysical_n0 */
    }

    /* Create FFTW plans to transform from physical to wave space */
    if (will_perform_r2c) {
        p->transpose_towaveB = underling_transpose_fftw_plan(
                p->problem->towaveB, p->data, p->data, rigor_flags);
        if (p->transpose_towaveB == NULL) {
            underling_plan_destroy(p);
            SUZERAIN_ERROR_NULL(
                    "FFTW MPI returned NULL plan: transpose_towaveB",
                    SUZERAIN_EFAILED);
        }

        p->transpose_towaveA = underling_transpose_fftw_plan(
                p->problem->towaveA, p->data, p->data, rigor_flags);
        if (p->transpose_towaveA == NULL) {
            underling_plan_destroy(p);
            SUZERAIN_ERROR_NULL(
                    "FFTW MPI returned NULL plan: transpose_towaveA",
                    SUZERAIN_EFAILED);
        }

        /* TODO Create FFT plan: r2c_towave_n0 */
        /* TODO Create FFT plan: c2c_towave_n1 */
    }

    return p;
}

void
underling_plan_destroy(
        underling_plan plan)
{
    if (plan) {
        plan->problem = NULL;
        if (plan->transpose_tophysA) {
            fftw_destroy_plan(plan->transpose_tophysA);
            plan->transpose_tophysA = NULL;
        }
        if (plan->c2c_tophysical_n1) {
            fftw_destroy_plan(plan->c2c_tophysical_n1);
            plan->c2c_tophysical_n1 = NULL;
        }
        if (plan->transpose_tophysB) {
            fftw_destroy_plan(plan->transpose_tophysB);
            plan->transpose_tophysB = NULL;
        }
        if (plan->c2r_tophysical_n0) {
            fftw_destroy_plan(plan->c2r_tophysical_n0);
            plan->c2r_tophysical_n0 = NULL;
        }
        if (plan->r2c_towave_n0) {
            fftw_destroy_plan(plan->r2c_towave_n0);
            plan->r2c_towave_n0 = NULL;
        }
        if (plan->transpose_towaveA) {
            fftw_destroy_plan(plan->transpose_towaveA);
            plan->transpose_towaveA = NULL;
        }
        if (plan->c2c_towave_n1) {
            fftw_destroy_plan(plan->c2c_towave_n1);
            plan->c2c_towave_n1 = NULL;
        }
        if (plan->transpose_towaveB) {
            fftw_destroy_plan(plan->transpose_towaveB);
            plan->transpose_towaveB = NULL;
        }
        plan->data = NULL;
    }
}

int
underling_execute_c2r(
        underling_plan plan)
{
    if (plan == NULL) {
        SUZERAIN_ERROR("non-NULL plan required", SUZERAIN_EINVAL);
    }

    /* FIXME check valid c2c_tophysical_n1, c2r_tophysical_n0 as well */
    const int valid =       plan->transpose_tophysA
                         && plan->transpose_tophysB;
    if (!valid) {
        SUZERAIN_ERROR("plan has one or more NULL c2r subplans",
                SUZERAIN_EINVAL);
    }

    fftw_execute(plan->transpose_tophysA);
/* FIXME fftw_execute(plan->c2c_tophysical_n1); */
    fftw_execute(plan->transpose_tophysB);
/* FIXME fftw_execute(plan->c2r_tophysical_n0); */

    return SUZERAIN_SUCCESS;
}

int
underling_execute_r2c(
        underling_plan plan)
{
    if (plan == NULL) {
        SUZERAIN_ERROR("non-NULL plan required", SUZERAIN_EINVAL);
    }
    /* FIXME check valid r2c_towave_n0, c2c_towave_n1 as well */
    const int valid =    plan->transpose_towaveA
                      && plan->transpose_towaveB;
    if (!valid) {
        SUZERAIN_ERROR("plan has one or more NULL r2c subplans",
                SUZERAIN_EINVAL);
    }

/* FIXME fftw_execute(plan->r2c_towave_n0); */
    fftw_execute(plan->transpose_towaveB);
/* FIXME fftw_execute(plan->c2c_towave_n1); */
    fftw_execute(plan->transpose_towaveA);

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
                "{np0=%d,nw0=%d,n1=%d,n2=%d}"
                ",{pA=%d,pB=%d}",
                grid->np0, grid->nw0, grid->n1, grid->n2,
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
        fprintf(output_file,"{nfields=%d,howmany=%d,local_size=%ld}",
                problem->nfields, problem->howmany, problem->local_size);
        fprintf(output_file,"{long_n2:");
        underling_fprint_extents(&(problem->long_n2), output_file);
        fprintf(output_file, "}");
        fprintf(output_file,"{long_n1:");
        underling_fprint_extents(&(problem->long_n1), output_file);
        fprintf(output_file, "}");
        fprintf(output_file,"{long_n0:");
        underling_fprint_extents(&(problem->long_n0), output_file);
        fprintf(output_file, "}");
        fprintf(output_file,"{tophysA:");
        underling_fprint_transpose(
                problem->tophysA, output_file);
        fprintf(output_file, "}");
        fprintf(output_file,"{tophysB:");
        underling_fprint_transpose(
                problem->tophysB, output_file);
        fprintf(output_file, "}");
        fprintf(output_file,"{towaveB:");
        underling_fprint_transpose(
                problem->towaveB, output_file);
        fprintf(output_file, "}");
        fprintf(output_file,"{towaveA:");
        underling_fprint_transpose(
                problem->towaveA, output_file);
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

        fprintf(output_file, "{underling_plan_c2r:");
        if (plan->transpose_tophysA) {
            fprintf(output_file, "{transpose_tophysA:");
            fftw_fprint_plan(plan->transpose_tophysA, output_file);
            fprintf(output_file, "}");
        }
        if (plan->c2c_tophysical_n1) {
            fprintf(output_file, "{c2c_tophysical_n1:");
            fftw_fprint_plan(plan->c2c_tophysical_n1, output_file);
            fprintf(output_file, "}");
        }
        if (plan->transpose_tophysB) {
            fprintf(output_file, "{plan->transpose_tophysB:");
            fftw_fprint_plan(plan->transpose_tophysB, output_file);
            fprintf(output_file, "}");
        }
        if (plan->c2r_tophysical_n0) {
            fprintf(output_file, "{plan->c2r_tophysical_n0:");
            fftw_fprint_plan(plan->c2r_tophysical_n0, output_file);
            fprintf(output_file, "}");
        }
        fprintf(output_file, "}");

        fprintf(output_file, "{underling_plan_r2c:");
        if (plan->r2c_towave_n0) {
            fprintf(output_file, "{r2c_towave_n0:");
            fftw_fprint_plan(plan->r2c_towave_n0, output_file);
            fprintf(output_file, "}");
        }
        if (plan->transpose_towaveA) {
            fprintf(output_file, "{transpose_towaveA:");
            fftw_fprint_plan(plan->transpose_towaveA, output_file);
            fprintf(output_file, "}");
        }
        if (plan->c2c_towave_n1) {
            fprintf(output_file, "{c2c_towave_n1:");
            fftw_fprint_plan(plan->c2c_towave_n1, output_file);
            fprintf(output_file, "}");
        }
        if (plan->transpose_towaveB) {
            fprintf(output_file, "{transpose_towaveB:");
            fftw_fprint_plan(plan->transpose_towaveB, output_file);
            fprintf(output_file, "}");
        }
        fprintf(output_file, "}");
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
                extents->local_nw0_start,
                extents->local_nw0_start + extents->local_nw0,
                extents->local_n1_start,
                extents->local_n1_start + extents->local_n1,
                extents->local_n2_start,
                extents->local_n2_start + extents->local_n2);
    }
}
