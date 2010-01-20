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

// Check MPI error code and call SUZERAIN_ERROR if code not MPI_SUCCESS.
// Returns SUZERAIN_EFAILED from the function.
#define MPICHKQ(stmt) \
    do { \
        int _chk_stat = (stmt); \
        if (_chk_stat != MPI_SUCCESS) { \
            char _chk_reason[255]; \
            char *_chk_mpistring = NULL; \
            int _chk_len; \
            const int _chk_string_stat \
                = MPI_Error_string(_chk_stat,_chk_mpistring,&_chk_len); \
            snprintf(_chk_reason, sizeof(_chk_reason)/sizeof(_chk_reason[0]), \
                    "Encountered MPI error code %d: %s", _chk_stat, \
                    (_chk_string_stat == MPI_SUCCESS) \
                    ? _chk_mpistring : "UNKNOWN"); \
            SUZERAIN_ERROR(_chk_reason, SUZERAIN_EFAILED); \
        } \
    } while(0)

// Check MPI error code and call SUZERAIN_ERROR if code not MPI_SUCCESS.
// Returns NULL from the function.
#define MPICHKN(stmt) \
    do { \
        int _chk_stat = (stmt); \
        if (_chk_stat != MPI_SUCCESS) { \
            char _chk_reason[255]; \
            char *_chk_mpistring = NULL; \
            int _chk_len; \
            const int _chk_string_stat \
                = MPI_Error_string(_chk_stat,_chk_mpistring,&_chk_len); \
            snprintf(_chk_reason, sizeof(_chk_reason)/sizeof(_chk_reason[0]), \
                    "Encountered MPI error code %d: %s", _chk_stat, \
                    (_chk_string_stat == MPI_SUCCESS) \
                    ? _chk_mpistring : "UNKNOWN"); \
            SUZERAIN_ERROR_NULL(_chk_reason, SUZERAIN_EFAILED); \
        } \
    } while(0)

// Check MPI error code and call SUZERAIN_ERROR if code not MPI_SUCCESS.
// Returns from a function with a void return value.
#define MPICHKV(stmt) \
    do { \
        int _chk_stat = (stmt); \
        if (_chk_stat != MPI_SUCCESS) { \
            char _chk_reason[255]; \
            char *_chk_mpistring = NULL; \
            int _chk_len; \
            const int _chk_string_stat \
                = MPI_Error_string(_chk_stat,_chk_mpistring,&_chk_len); \
            snprintf(_chk_reason, sizeof(_chk_reason)/sizeof(_chk_reason[0]), \
                    "Encountered MPI error code %d: %s", _chk_stat, \
                    (_chk_string_stat == MPI_SUCCESS) \
                    ? _chk_mpistring : "UNKNOWN"); \
            SUZERAIN_ERROR_VOID(_chk_reason, SUZERAIN_EFAILED); \
        } \
    } while(0)

underling_grid *
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
    MPICHKN(MPI_Comm_size(comm, &nproc));

    // Create a balanced processor grid if not specified by p0, p1 != 0
    {
        int dims[2] = { p0, p1 };
        MPICHKN(MPI_Dims_create(nproc, 2, dims));
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
        MPICHKN(MPI_Cart_create(
                comm, 2, dims, periodic, 1/*reordering allowed*/, &g_comm));
    }
    // Create communicators for the p0, P1 directions
    MPI_Comm p0_comm;
    {
        int remain_dims[2] = { 1, 0 };
        MPICHKN(MPI_Cart_sub(g_comm, remain_dims, &p0_comm));
    }
    MPI_Comm p1_comm;
    {
        int remain_dims[2] = { 0, 1 };
        MPICHKN(MPI_Cart_sub(g_comm, remain_dims, &p1_comm));
    }

    // Create and initialize the grid workspace
    underling_grid * g = malloc(sizeof(underling_grid));
    if (g == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate space for grid",
                             SUZERAIN_ENOMEM);
    }
    // Copy the grid parameters to the grid workspace
    g->np0     = np0;
    g->nw0     = nw0;
    g->n1      = n1;
    g->n2      = n2;
    g->p0      = p0;
    g->p1      = p1;
    g->block_a = block_a;
    g->g_comm  = g_comm;
    g->p0_comm = p0_comm;
    g->p1_comm = p1_comm;

    return g;
}

void
underling_grid_destroy(underling_grid * grid)
{
    if (grid) {
        if (grid->g_comm) {
            MPICHKV(MPI_Comm_disconnect(&grid->g_comm));
            grid->g_comm = MPI_COMM_NULL;
        }
        if (grid->p0_comm) {
            MPICHKV(MPI_Comm_disconnect(&grid->p0_comm));
            grid->p0_comm = MPI_COMM_NULL;
        }
        if (grid->p1_comm) {
            MPICHKV(MPI_Comm_disconnect(&grid->p1_comm));
            grid->p1_comm = MPI_COMM_NULL;
        }
        free(grid);
    }
}


underling_problem *
underling_problem_create(
        underling_grid *grid,
        int nfields)
{
    if (nfields < 1) {
        SUZERAIN_ERROR_NULL("nfields >= 1 required", SUZERAIN_EINVAL);
    }

    // Create and initialize the problem workspace
    underling_problem * p = malloc(sizeof(underling_problem));
    if (p == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate space for problem",
                             SUZERAIN_ENOMEM);
    }
    // Copy the problem parameters to the problem workspace
    p->nfields = nfields;

    {
    }

    return p;
}

void
underling_problem_destroy(
        underling_problem * problem)
{
    if (problem) {
        free(problem);
    }
}
