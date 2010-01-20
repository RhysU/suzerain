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
            int tmp = dims[0]; dims[0] = dims[1]; dims[1] = tmp;
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

    // Transform proceeds like
    // ALL WAVE SPACE
    //      ((nw0/p0 x nw1/p1)) x nw2
    //       \----block_a----/
    //      to
    //      ((nw2/p1 x nw0/p0)) x nw1
    //      ((nw2/p1 x nw0/p0)) x np1
    //       \----block_b----/
    //      to
    //      ((np1/p0 x nw2/p1)) x nw0
    //      ((np1/p0 x nw2/p1)) x np0
    //       \----block_c----/
    // ALL PHYSICAL_SPACE

    // Ensure grid allows the problem data to be distributed evenly.
    // This is a dumb restriction that should probably be lifted.
    if ((nw0*n1) % (p0*p1) != 0) {
        char reason[127];
        snprintf(reason, sizeof(reason)/sizeof(reason[0]),
                "Invalid processor grid: "
                "(nw0 {%d} * n1 {%d}) %% (p0 {%d}* p1 {%d}) != 0",
                nw0, n1, p0, p1);
        SUZERAIN_ERROR_NULL(reason, SUZERAIN_EFAILED);
    }
    const int block_a = (nw0*n1)/(p0*p1);

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
    g->block_a     = block_a;
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

    { /* Wave towards physical MPI transpose: long in n2 to long in n1 */
        underling_transpose_details * const d = &(p->tophysical_A);
        d->n[0]    = grid->p1 * grid->block_a;
        d->n[1]    = grid->n2;
        d->howmany = nfields * sizeof(underling_complex)/sizeof(underling_real);
        d->comm    = grid->p1_comm;
        d->flags   = 0;
        d->local_size = fftw_mpi_local_size_many_transposed(
                    /*rank*/2,
                    d->n,
                    d->howmany,
                    FFTW_MPI_DEFAULT_BLOCK,
                    FFTW_MPI_DEFAULT_BLOCK,
                    d->comm,
                    &(d->local_n0),
                    &(d->local_n0_start),
                    &(d->local_n1),
                    &(d->local_n1_start)
                );
    }
    { /* Wave towards physical MPI transpose: long in n1 to long in n0 */
        underling_transpose_details * const d = &(p->tophysical_B);
        d->n[0]    = grid->p0 * p->tophysical_A.local_n1 * grid->nw0;
        d->n[1]    = grid->n1;
        d->howmany = p->tophysical_A.howmany; /* copy for consistency */
        d->comm    = grid->p0_comm;
        d->flags   = 0;
        d->local_size = fftw_mpi_local_size_many_transposed(
                    /*rank*/2,
                    d->n,
                    d->howmany,
                    FFTW_MPI_DEFAULT_BLOCK,
                    FFTW_MPI_DEFAULT_BLOCK,
                    d->comm,
                    &(d->local_n0),
                    &(d->local_n0_start),
                    &(d->local_n1),
                    &(d->local_n1_start)
                );
    }
    { /* Physical towards wave MPI transpose: long in n0 to long in n1 */
        underling_transpose_details * const d = &(p->towave_B);
        d->n[0]    = p->tophysical_B.n[1];
        d->n[1]    = p->tophysical_B.n[0];
        d->howmany = p->tophysical_A.howmany; /* copy for consistency */
        d->comm    = p->tophysical_B.comm;    /* copy for consistency */
        d->flags   = 0;
        d->local_size = fftw_mpi_local_size_many_transposed(
                    /*rank*/2,
                    d->n,
                    d->howmany,
                    FFTW_MPI_DEFAULT_BLOCK,
                    FFTW_MPI_DEFAULT_BLOCK,
                    d->comm,
                    &(d->local_n0),
                    &(d->local_n0_start),
                    &(d->local_n1),
                    &(d->local_n1_start)
                );
    }
    { /* Physical towards wave MPI transpose: long in n1 to long in n2 */
        underling_transpose_details * const d = &(p->towave_A);
        d->n[0]    = p->tophysical_A.n[1];
        d->n[1]    = p->tophysical_A.n[0];
        d->howmany = p->tophysical_A.howmany; /* copy for consistency */
        d->comm    = p->tophysical_A.comm;    /* copy for consistency */
        d->flags   = 0;
        d->local_size = fftw_mpi_local_size_many_transposed(
                    /*rank*/2,
                    d->n,
                    d->howmany,
                    FFTW_MPI_DEFAULT_BLOCK,
                    FFTW_MPI_DEFAULT_BLOCK,
                    d->comm,
                    &(d->local_n0),
                    &(d->local_n0_start),
                    &(d->local_n1),
                    &(d->local_n1_start)
                );
    }

    /* local_size is overall maximum of all transpose local_size values */
    p->local_size = p->tophysical_A.local_size;
    if (p->local_size < p->tophysical_B.local_size) {
        p->local_size = p->tophysical_B.local_size;
    }
    if (p->local_size < p->towave_B.local_size) {
        p->local_size = p->towave_B.local_size;
    }
    if (p->local_size < p->towave_A.local_size) {
        p->local_size = p->towave_A.local_size;
    }

    return p;
}

size_t
underling_local_size(
        underling_problem problem)
{
    return problem->local_size;
}

void
underling_problem_destroy(
        underling_problem problem)
{
    if (problem) {
        problem->grid = NULL;
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
        {
            underling_transpose_details * const d
                = &(p->problem->tophysical_A);
            p->transpose_tophysical_A = fftw_mpi_plan_many_transpose(
                        d->n[0],
                        d->n[1],
                        d->howmany,
                        FFTW_MPI_DEFAULT_BLOCK,
                        FFTW_MPI_DEFAULT_BLOCK,
                        p->data, /* in-place */
                        p->data, /* in-place */
                        d->comm,
                        d->flags | rigor_flags
                    );
            if (p->transpose_tophysical_A == NULL) {
                underling_plan_destroy(p);
                SUZERAIN_ERROR_NULL(
                        "FFTW MPI returned NULL plan: transpose_tophysical_A",
                        SUZERAIN_EFAILED);
            }
        }
        {
            underling_transpose_details * const d
                = &(p->problem->tophysical_B);
            p->transpose_tophysical_B = fftw_mpi_plan_many_transpose(
                        d->n[0],
                        d->n[1],
                        d->howmany,
                        FFTW_MPI_DEFAULT_BLOCK,
                        FFTW_MPI_DEFAULT_BLOCK,
                        p->data, /* in-place */
                        p->data, /* in-place */
                        d->comm,
                        d->flags | rigor_flags
                    );
            if (p->transpose_tophysical_B == NULL) {
                underling_plan_destroy(p);
                SUZERAIN_ERROR_NULL(
                        "FFTW MPI returned NULL plan: transpose_tophysical_B",
                        SUZERAIN_EFAILED);
            }
        }
        /* TODO Create FFT plan: c2c_tophysical_n1 */
        /* TODO Create FFT plan: c2r_tophysical_n0 */
    }

    /* Create FFTW plans to transform from physical to wave space */
    if (will_perform_r2c) {
        {
            underling_transpose_details * const d
                = &(p->problem->towave_B);
            p->transpose_towave_B = fftw_mpi_plan_many_transpose(
                        d->n[0],
                        d->n[1],
                        d->howmany,
                        FFTW_MPI_DEFAULT_BLOCK,
                        FFTW_MPI_DEFAULT_BLOCK,
                        p->data, /* in-place */
                        p->data, /* in-place */
                        d->comm,
                        d->flags | rigor_flags
                    );
            if (p->transpose_towave_B == NULL) {
                underling_plan_destroy(p);
                SUZERAIN_ERROR_NULL(
                        "FFTW MPI returned NULL plan: transpose_towave_B",
                        SUZERAIN_EFAILED);
            }
        }
        {
            underling_transpose_details * const d
                = &(p->problem->towave_A);
            p->transpose_towave_A = fftw_mpi_plan_many_transpose(
                        d->n[0],
                        d->n[1],
                        d->howmany,
                        FFTW_MPI_DEFAULT_BLOCK,
                        FFTW_MPI_DEFAULT_BLOCK,
                        p->data, /* in-place */
                        p->data, /* in-place */
                        d->comm,
                        d->flags | rigor_flags
                    );
            if (p->transpose_towave_A == NULL) {
                underling_plan_destroy(p);
                SUZERAIN_ERROR_NULL(
                        "FFTW MPI returned NULL plan: transpose_towave_A",
                        SUZERAIN_EFAILED);
            }
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
    fftw_execute(plan->transpose_towave_A);
/* FIXME fftw_execute(plan->c2c_towave_n1); */
    fftw_execute(plan->transpose_towave_B);

    return SUZERAIN_SUCCESS;
}

void
underling_fprint_plan(
        const underling_plan plan, FILE *output_file)
{
    fprintf(output_file, "{underling_plan:");
    if (!plan) {
        fprintf(output_file, "NULL");
    } else {

        fprintf(output_file, "\n{underling_plan_c2r:");
        if (plan->transpose_tophysical_A) {
            fprintf(output_file, "\n{transpose_tophysical_A:");
            fftw_fprint_plan(plan->transpose_tophysical_A, output_file);
            fprintf(output_file, "\n}");
        }
        if (plan->c2c_tophysical_n1) {
            fprintf(output_file, "\n{c2c_tophysical_n1:");
            fftw_fprint_plan(plan->c2c_tophysical_n1, output_file);
            fprintf(output_file, "\n}");
        }
        if (plan->transpose_tophysical_B) {
            fprintf(output_file, "\n{plan->transpose_tophysical_B:");
            fftw_fprint_plan(plan->transpose_tophysical_B, output_file);
            fprintf(output_file, "\n}");
        }
        if (plan->c2r_tophysical_n0) {
            fprintf(output_file, "\n{plan->c2r_tophysical_n0:");
            fftw_fprint_plan(plan->c2r_tophysical_n0, output_file);
            fprintf(output_file, "\n}");
        }
        fprintf(output_file, "\n}");

        fprintf(output_file, "\n{underling_plan_r2c:");
        if (plan->r2c_towave_n0) {
            fprintf(output_file, "\n{r2c_towave_n0:");
            fftw_fprint_plan(plan->r2c_towave_n0, output_file);
            fprintf(output_file, "\n}");
        }
        if (plan->transpose_towave_A) {
            fprintf(output_file, "\n{transpose_towave_A:");
            fftw_fprint_plan(plan->transpose_towave_A, output_file);
            fprintf(output_file, "\n}");
        }
        if (plan->c2c_towave_n1) {
            fprintf(output_file, "\n{c2c_towave_n1:");
            fftw_fprint_plan(plan->c2c_towave_n1, output_file);
            fprintf(output_file, "\n}");
        }
        if (plan->transpose_towave_B) {
            fprintf(output_file, "\n{transpose_towave_B:");
            fftw_fprint_plan(plan->transpose_towave_B, output_file);
            fprintf(output_file, "\n}");
        }
        fprintf(output_file, "\n}");
    }
}

void
underling_print_plan(
        const underling_plan plan)
{
    underling_fprint_plan(plan, stdout);
}
