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
    int p0;
    int p1;
    MPI_Comm g_comm;
    int g_rank;
    int g_coords[2];
    MPI_Comm p0_comm;
    int p0_rank;
    int p0_coord;
    MPI_Comm p1_comm;
    int p1_rank;
    int p1_coord;
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

struct underling_problem_s {
    underling_grid grid;              // grid owns its resources
    int nfields;                      // # of complex-valued state...
    int howmany;                      // ...as a # of real values
    underling_transpose tophysical_A; // n2 long to n1 long
    underling_transpose tophysical_B; // n1 long to n0 long
    underling_transpose towave_B;     // n0 long to n1 long
    underling_transpose towave_A;     // n1 long to n2 long
    ptrdiff_t local_size;             // Max of all local sizes
};

struct underling_plan_s {
    underling_problem problem;         // problem owns its resources
    fftw_plan transpose_tophysical_A;
    fftw_plan c2c_tophysical_n1;
    fftw_plan transpose_tophysical_B;
    fftw_plan c2r_tophysical_n0;
    fftw_plan r2c_towave_n0;
    fftw_plan transpose_towave_A;
    fftw_plan c2c_towave_n1;
    fftw_plan transpose_towave_B;
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

// **************************************************************************
// IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION
// **************************************************************************

underling_grid
underling_grid_create(
        MPI_Comm comm,
        int np0,
        int np1,
        int np2,
        int p0,
        int p1)
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
    if (p0 < 0) {
        SUZERAIN_ERROR_NULL("p0 >= 0 required", SUZERAIN_EINVAL);
    }
    if (p1 < 0) {
        SUZERAIN_ERROR_NULL("p1 >= 0 required", SUZERAIN_EINVAL);
    }

    // Get number of processors in the communicator
    int nproc;
    SUZERAIN_MPICHKN(MPI_Comm_size(comm, &nproc));

    // Create a balanced processor grid if not specified by p0, p1 != 0
    {
        int dims[2] = { p0, p1 };
        SUZERAIN_MPICHKN(MPI_Dims_create(nproc, 2, dims));
        // If both directions automatic, ensure dims[0] <= dims[1]
        if (p0 == 0 && p1 == 0 && dims[0] > dims[1]) {
            const int tmp = dims[0]; dims[0] = dims[1]; dims[1] = tmp;
        }
        p0 = dims[0];
        p1 = dims[1];
        if (p0 * p1 != nproc) {
            char reason[127];
            snprintf(reason, sizeof(reason)/sizeof(reason[0]),
                    "Invalid processor grid: p0 {%d} * p1 {%d} != nproc {%d}",
                    p0, p1, nproc);
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
        int dims[2]     = { p0, p1 };
        int periodic[2] = { 0, 0 };
        SUZERAIN_MPICHKN(MPI_Cart_create(
                comm, 2, dims, periodic, 1/*reordering allowed*/, &g_comm));
    }
    // Cache the rank and coordinates of this process within g_comm
    int g_rank;
    SUZERAIN_MPICHKN(MPI_Comm_rank(g_comm, &g_rank));
    int g_coords[2];
    SUZERAIN_MPICHKN(MPI_Cart_coords(g_comm, g_rank, 2, g_coords));

    // Create communicator for the P0 direction
    MPI_Comm p0_comm;
    {
        int remain_dims[2] = { 1, 0 };
        SUZERAIN_MPICHKN(MPI_Cart_sub(g_comm, remain_dims, &p0_comm));
    }
    // Cache the rank and coordinates of this process within g_comm
    int p0_rank;
    SUZERAIN_MPICHKN(MPI_Comm_rank(p0_comm, &p0_rank));
    int p0_coord;
    SUZERAIN_MPICHKN(MPI_Cart_coords(p0_comm, p0_rank, 1, &p0_coord));

    // Create communicator for the P1 direction
    MPI_Comm p1_comm;
    {
        int remain_dims[2] = { 0, 1 };
        SUZERAIN_MPICHKN(MPI_Cart_sub(g_comm, remain_dims, &p1_comm));
    }
    // Cache the rank and coordinates of this process within g_comm
    int p1_rank;
    SUZERAIN_MPICHKN(MPI_Comm_rank(p1_comm, &p1_rank));
    int p1_coord;
    SUZERAIN_MPICHKN(MPI_Cart_coords(p1_comm, p1_rank, 1, &p1_coord));

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
    g->p0          = p0;
    g->p1          = p1;
    g->g_comm      = g_comm;
    g->g_rank      = g_rank;
    g->g_coords[0] = g_coords[0];
    g->g_coords[1] = g_coords[1];
    g->p0_comm     = p0_comm;
    g->p0_rank     = p0_rank;
    g->p0_coord    = p0_coord;
    g->p1_comm     = p1_comm;
    g->p1_rank     = p1_rank;
    g->p1_coord    = p1_coord;

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
        if (grid->p0_comm) {
            SUZERAIN_MPICHKV(MPI_Comm_disconnect(&grid->p0_comm));
            grid->p0_comm = MPI_COMM_NULL;
        }
        if (grid->p1_comm) {
            SUZERAIN_MPICHKV(MPI_Comm_disconnect(&grid->p1_comm));
            grid->p1_comm = MPI_COMM_NULL;
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

    // Wave space is initially partitioned across a 2D (p0 x p1) topology
    // according to (nw0 x n1) x n2 with n2 long.  Use FFTW's partitioning to
    // find the portion of the global (nw0 x n1) data that is spread across
    // p1_comm.  We never perform this next transpose, but we do need its
    // local_n0.
    ptrdiff_t p1_specific_nw0n1;
    {
        underling_transpose partition_nw0n1_by_n2_across_p1
            = underling_transpose_create(grid->nw0 * grid->n1,
                                         grid->n2,
                                         p->howmany,
                                         FFTW_MPI_DEFAULT_BLOCK,
                                         FFTW_MPI_DEFAULT_BLOCK,
                                         grid->p1_comm,
                                         /*flags*/0);
        if (partition_nw0n1_by_n2_across_p1 == NULL) {
            underling_problem_destroy(p);
            SUZERAIN_ERROR_NULL(
                    "failed determining partition_nw0n1_by_n2_across_p1",
                    SUZERAIN_EFAILED);
        }
        // Save the partitioning information that we need and destroy the rest
        p1_specific_nw0n1 = partition_nw0n1_by_n2_across_p1->local_d0;
        underling_transpose_destroy(partition_nw0n1_by_n2_across_p1);
    }

    // Wave towards physical MPI transpose: long in n2 to long in n1
    // That is, (nw0 x n1) x n2 becomes n2 x (nw0 x n1)
    p->tophysical_A = underling_transpose_create(
            p1_specific_nw0n1,
            grid->n2,
            p->howmany,
            FFTW_MPI_DEFAULT_BLOCK,
            FFTW_MPI_DEFAULT_BLOCK,
            grid->p0_comm,
            /*flags*/0);
    if (p->tophysical_A == NULL) {
        underling_problem_destroy(p);
        SUZERAIN_ERROR_NULL("failed creating p->tophysical_A",
                SUZERAIN_EFAILED);
    }

    // The global transform size for n2 x (nw0 x n1) = (n2 x nw0) x n1
    // redistribution across p1_comm depends on the rank of the process within
    // p0_comm.  We must accumulate the p0-specific global (n2 x nw0) value.
    ptrdiff_t p0_specific_global_n2nw0;
    const int allreduce_error = MPI_Allreduce(
            &(p->tophysical_A->d0), &p0_specific_global_n2nw0, 1,
            MPI_LONG, MPI_SUM, grid->p1_comm);
    if (allreduce_error) {
        underling_problem_destroy(p);
        SUZERAIN_MPICHKN(allreduce_error);
    }
    assert(p0_specific_global_n2nw0 % grid->n1 == 0);
    p0_specific_global_n2nw0 /= grid->n1;
    p0_specific_global_n2nw0 *= p->tophysical_A->local_d1;

    // Wave towards physical MPI transpose: long in n1 to long in nw0
    // That is, (n2 x nw0) x n1 becomes n1 x (n2 x nw0)
    p->tophysical_B = underling_transpose_create(
            p0_specific_global_n2nw0,
            grid->n1,
            p->howmany,
            FFTW_MPI_DEFAULT_BLOCK,
            FFTW_MPI_DEFAULT_BLOCK,
            grid->p1_comm,
            /*flags*/0);
    if (p->tophysical_B == NULL) {
        underling_problem_destroy(p);
        SUZERAIN_ERROR_NULL("failed creating p->tophysical_B",
                SUZERAIN_EFAILED);
    }

    // Physical towards wave MPI transpose: long in n0 to long in n1
    p->towave_B = underling_transpose_create_inverse(p->tophysical_B);
    if (p->towave_B == NULL) {
        underling_problem_destroy(p);
        SUZERAIN_ERROR_NULL("failed creating p->towave_B",
                SUZERAIN_EFAILED);
    }

    // Physical towards wave MPI transpose: long in n1 to long in n2
    p->towave_A = underling_transpose_create_inverse(p->tophysical_A);
    if (p->towave_A == NULL) {
        underling_problem_destroy(p);
        SUZERAIN_ERROR_NULL("failed creating p->towave_A",
                SUZERAIN_EFAILED);
    }

    // local_size is overall maximum of all transpose local_size values
    p->local_size = p->tophysical_A->local_size;
    if (p->local_size < p->tophysical_B->local_size) {
        p->local_size = p->tophysical_B->local_size;
    }
    if (p->local_size < p->towave_B->local_size) {
        p->local_size = p->towave_B->local_size;
    }
    if (p->local_size < p->towave_A->local_size) {
        p->local_size = p->towave_A->local_size;
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

    const size_t nprocessors =   problem->grid->p0
                               * problem->grid->p1;

    return global_data/nprocessors;
}

void
underling_problem_destroy(
        underling_problem problem)
{
    if (problem) {
        problem->grid = NULL;
        if (problem->tophysical_A) {
            underling_transpose_destroy(problem->tophysical_A);
            problem->tophysical_A = NULL;
        }
        if (problem->tophysical_B) {
            underling_transpose_destroy(problem->tophysical_B);
            problem->tophysical_B = NULL;
        }
        if (problem->towave_B) {
            underling_transpose_destroy(problem->towave_B);
            problem->towave_B = NULL;
        }
        if (problem->towave_A) {
            underling_transpose_destroy(problem->towave_A);
            problem->towave_A = NULL;
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
        p->transpose_tophysical_A = underling_transpose_fftw_plan(
                p->problem->tophysical_A, p->data, p->data, rigor_flags);
        if (p->transpose_tophysical_A == NULL) {
            underling_plan_destroy(p);
            SUZERAIN_ERROR_NULL(
                    "FFTW MPI returned NULL plan: transpose_tophysical_A",
                    SUZERAIN_EFAILED);
        }

        p->transpose_tophysical_B = underling_transpose_fftw_plan(
                p->problem->tophysical_B, p->data, p->data, rigor_flags);
        if (p->transpose_tophysical_B == NULL) {
            underling_plan_destroy(p);
            SUZERAIN_ERROR_NULL(
                    "FFTW MPI returned NULL plan: transpose_tophysical_B",
                    SUZERAIN_EFAILED);
        }

        /* TODO Create FFT plan: c2c_tophysical_n1 */
        /* TODO Create FFT plan: c2r_tophysical_n0 */
    }

    /* Create FFTW plans to transform from physical to wave space */
    if (will_perform_r2c) {
        p->transpose_towave_B = underling_transpose_fftw_plan(
                p->problem->towave_B, p->data, p->data, rigor_flags);
        if (p->transpose_towave_B == NULL) {
            underling_plan_destroy(p);
            SUZERAIN_ERROR_NULL(
                    "FFTW MPI returned NULL plan: transpose_towave_B",
                    SUZERAIN_EFAILED);
        }

        p->transpose_towave_A = underling_transpose_fftw_plan(
                p->problem->towave_A, p->data, p->data, rigor_flags);
        if (p->transpose_towave_A == NULL) {
            underling_plan_destroy(p);
            SUZERAIN_ERROR_NULL(
                    "FFTW MPI returned NULL plan: transpose_towave_A",
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
        if (plan->transpose_tophysical_A) {
            fftw_destroy_plan(plan->transpose_tophysical_A);
            plan->transpose_tophysical_A = NULL;
        }
        if (plan->c2c_tophysical_n1) {
            fftw_destroy_plan(plan->c2c_tophysical_n1);
            plan->c2c_tophysical_n1 = NULL;
        }
        if (plan->transpose_tophysical_B) {
            fftw_destroy_plan(plan->transpose_tophysical_B);
            plan->transpose_tophysical_B = NULL;
        }
        if (plan->c2r_tophysical_n0) {
            fftw_destroy_plan(plan->c2r_tophysical_n0);
            plan->c2r_tophysical_n0 = NULL;
        }
        if (plan->r2c_towave_n0) {
            fftw_destroy_plan(plan->r2c_towave_n0);
            plan->r2c_towave_n0 = NULL;
        }
        if (plan->transpose_towave_A) {
            fftw_destroy_plan(plan->transpose_towave_A);
            plan->transpose_towave_A = NULL;
        }
        if (plan->c2c_towave_n1) {
            fftw_destroy_plan(plan->c2c_towave_n1);
            plan->c2c_towave_n1 = NULL;
        }
        if (plan->transpose_towave_B) {
            fftw_destroy_plan(plan->transpose_towave_B);
            plan->transpose_towave_B = NULL;
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
    const int valid =       plan->transpose_tophysical_A
                         && plan->transpose_tophysical_B;
    if (!valid) {
        SUZERAIN_ERROR("plan has one or more NULL c2r subplans",
                SUZERAIN_EINVAL);
    }

    fftw_execute(plan->transpose_tophysical_A);
/* FIXME fftw_execute(plan->c2c_tophysical_n1); */
    fftw_execute(plan->transpose_tophysical_B);
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
    const int valid =    plan->transpose_towave_A
                      && plan->transpose_towave_B;
    if (!valid) {
        SUZERAIN_ERROR("plan has one or more NULL r2c subplans",
                SUZERAIN_EINVAL);
    }

/* FIXME fftw_execute(plan->r2c_towave_n0); */
    fftw_execute(plan->transpose_towave_B);
/* FIXME fftw_execute(plan->c2c_towave_n1); */
    fftw_execute(plan->transpose_towave_A);

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
                ",{p0=%d,p1=%d}"
                ",{g_comm=%x,p0_comm=%x,p1_comm=%x}",
                grid->np0, grid->nw0, grid->n1, grid->n2,
                grid->p0, grid->p1,
                grid->g_comm, grid->p0_comm, grid->p1_comm);
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
        fprintf(output_file,"{tophysical_A:");
        underling_fprint_transpose(
                problem->tophysical_A, output_file);
        fprintf(output_file, "}");
        fprintf(output_file,"{tophysical_B:");
        underling_fprint_transpose(
                problem->tophysical_B, output_file);
        fprintf(output_file, "}");
        fprintf(output_file,"{towave_B:");
        underling_fprint_transpose(
                problem->towave_B, output_file);
        fprintf(output_file, "}");
        fprintf(output_file,"{towave_A:");
        underling_fprint_transpose(
                problem->towave_A, output_file);
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
        if (plan->transpose_tophysical_A) {
            fprintf(output_file, "{transpose_tophysical_A:");
            fftw_fprint_plan(plan->transpose_tophysical_A, output_file);
            fprintf(output_file, "}");
        }
        if (plan->c2c_tophysical_n1) {
            fprintf(output_file, "{c2c_tophysical_n1:");
            fftw_fprint_plan(plan->c2c_tophysical_n1, output_file);
            fprintf(output_file, "}");
        }
        if (plan->transpose_tophysical_B) {
            fprintf(output_file, "{plan->transpose_tophysical_B:");
            fftw_fprint_plan(plan->transpose_tophysical_B, output_file);
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
        if (plan->transpose_towave_A) {
            fprintf(output_file, "{transpose_towave_A:");
            fftw_fprint_plan(plan->transpose_towave_A, output_file);
            fprintf(output_file, "}");
        }
        if (plan->c2c_towave_n1) {
            fprintf(output_file, "{c2c_towave_n1:");
            fftw_fprint_plan(plan->c2c_towave_n1, output_file);
            fprintf(output_file, "}");
        }
        if (plan->transpose_towave_B) {
            fprintf(output_file, "{transpose_towave_B:");
            fftw_fprint_plan(plan->transpose_towave_B, output_file);
            fprintf(output_file, "}");
        }
        fprintf(output_file, "}");
    }
}
