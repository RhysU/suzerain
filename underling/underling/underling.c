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
#include <underling/underling.h>

#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <fftw3-mpi.h>
#include <underling/error.h>
#include "common.h"

// TODO Ensure grid/problem compatibility when both provided to methods!
// TODO Add wisdom broadcasting and FFTW_WISDOM_ONLY handling
// TODO Check for memory leaks stemming from planning failures

// ********************************************************************
// INTERNAL TYPES INTERNAL TYPES INTERNAL TYPES INTERNAL TYPES INTERNAL
// ********************************************************************

// Quasi-hidden flag used for debugging deadlock conditions
static TLS MPI_Comm debug_comm = MPI_COMM_NULL;

// Employing FFTW_DESTROY_INPUT for transpose planning and/or execution can
// deadlock prior to 3.3-beta1.  See Redmine ticket #1297 for full details.
#ifdef HAVE_FFTW33_BETA1_OR_LATER
static const int destroy_flags = FFTW_DESTROY_INPUT;
#else
static const int destroy_flags = 0;
#endif

// TODO Document internal structures

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
    MPI_Comm  comm;                   // Instance owns MPI resource
    ptrdiff_t d[2];
    ptrdiff_t howmany;
    ptrdiff_t block[2];
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
    _Bool in_place;
    fftw_plan plan_backwardA;
    fftw_plan plan_backwardB;
    fftw_plan plan_forwardA;
    fftw_plan plan_forwardB;
};

// ********************************************************************
// INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL
// ********************************************************************

// TODO Document internal-only methods

static
MPI_Comm
underling_MPI_Comm_dup_with_name(MPI_Comm comm);

static
underling_transpose
underling_transpose_create(
        ptrdiff_t d0,
        ptrdiff_t d1,
        ptrdiff_t howmany,
        ptrdiff_t block0,
        ptrdiff_t block1,
        MPI_Comm comm,
        unsigned flags);

static
underling_transpose
underling_transpose_create_inverse(
        const underling_transpose forward);

static
void
underling_transpose_destroy(
        underling_transpose transpose);

static
fftw_plan
underling_transpose_fftw_plan(
        const underling_transpose transpose,
        underling_real *in,
        underling_real *out,
        unsigned fftw_flags);

static
void
underling_dump_transposes(
        MPI_Comm dump_comm,
        FILE *out,
        const char *prefix,
        const underling_transpose transpose);

static
size_t
underling_local_memory_allreduce(
        const underling_grid    grid,
        const underling_problem problem,
        MPI_Op op);

static
void
underling_fprint_transpose(
        const underling_transpose transpose,
        FILE *output_file);

// **************************************************************************
// IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION
// **************************************************************************

const underling_extents UNDERLING_EXTENTS_INVALID = {
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 0
};

int
underling_extents_cmp(const underling_extents * const e1,
                      const underling_extents * const e2)
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

int
underling_only_init(int *argc, char **argv[])
{
    // Placeholder for future function.  Does nothing currently.

    (void) argc; // Unused
    (void) argv; // Unused

    return UNDERLING_SUCCESS;
}

int
underling_init(int *argc, char **argv[], int nthreads)
{
    int error, flag;

    // Initialize MPI (if necessary)
    if ((error = MPI_Initialized(&flag))) {
        UNDERLING_MPICHKR(error);
        return UNDERLING_ESANITY;
    }
    if (!flag) {
        if ((error = MPI_Init(argc, argv))) {
            UNDERLING_MPICHKR(error);
            return UNDERLING_ESANITY;
        }
    }

    // Initialize FFTW threads (if necessary)
#ifdef HAVE_FFTW3_THREADS
    if (!fftw_init_threads()) {
        UNDERLING_ERROR("fftw_init_threads reported an error",
                        UNDERLING_ESANITY);
    }

    // Possibly lookup nthreads from environment for both OpenMP /and/ pthreads
    if (nthreads == 0) {
        const char * const envstr = getenv("OMP_NUM_THREADS");
        if (envstr) {
            const int envnum = atoi(envstr);
            if (envnum < 1) {
                UNDERLING_ERROR(
                        "Malformed OMP_NUM_THREADS environment variable",
                        UNDERLING_EINVAL);
            }
            nthreads = envnum;
        } else {
            // Value not found in environment so use a sane default
            nthreads = 1;
        }
    }
    fftw_plan_with_nthreads(nthreads);
#endif /* HAVE_FFTW3_THREADS */

    // Initialize FFTW MPI
    fftw_mpi_init();

    // Finally, initialize underling itself
    return underling_only_init(argc, argv);
}

void
underling_only_cleanup()
{
    // Placeholder for future function.  Does nothing currently.
    return;
}

void underling_cleanup()
{
    fftw_mpi_cleanup();
#ifdef HAVE_FFTW3_THREADS
    fftw_cleanup_threads();
#endif
    underling_only_cleanup();

    return;
}

static
MPI_Comm
underling_MPI_Comm_dup_with_name(MPI_Comm comm)
{
    if (comm == MPI_COMM_NULL) return MPI_COMM_NULL;

    char buffer[MPI_MAX_OBJECT_NAME] = ""; // Why does valgrind complain?
    int resultlen = 0;
    const int get_name_error = MPI_Comm_get_name(comm, buffer, &resultlen);
    if (get_name_error) {
        UNDERLING_MPICHKR(get_name_error /* MPI_Comm_get_name */);
        return MPI_COMM_NULL;
    }

    MPI_Comm retval = MPI_COMM_NULL;

    const int dup_error = MPI_Comm_dup(comm, &retval);
    if (dup_error) {
        UNDERLING_MPICHKR(dup_error /* MPI_Comm_dup */);
        return MPI_COMM_NULL;
    }

    if (resultlen > 0) {
        const int set_name_error = MPI_Comm_set_name(retval, buffer);
        if (set_name_error) {
            UNDERLING_MPICHKR(set_name_error /* MPI_Comm_set_name */);
            UNDERLING_MPICHKR(MPI_Comm_free(&retval));
            return MPI_COMM_NULL;
        }
    }

    return retval;
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
    if (UNDERLING_UNLIKELY(comm == MPI_COMM_NULL)) {
        UNDERLING_ERROR_NULL("comm != MPI_COMM_NULL required",
                             UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(n0 < 1)) {
        UNDERLING_ERROR_NULL("n0 >= 1 required", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(n1 < 1)) {
        UNDERLING_ERROR_NULL("n1 >= 1 required", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(n2 < 1)) {
        UNDERLING_ERROR_NULL("n2 >= 1 required", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(pA < 0)) {
        UNDERLING_ERROR_NULL("pA >= 0 required", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(pB < 0)) {
        UNDERLING_ERROR_NULL("pB >= 0 required", UNDERLING_EINVAL);
    }

    // Local scratch variables
    int error;                         // For error reporting
    char buffer[MPI_MAX_OBJECT_NAME];  // For snprintf calls
    int dims[2];                       // For MPI routine dimensions

    // Allocate memory for the grid workspace
    underling_grid g = malloc(sizeof(struct underling_grid_s));
    if (UNDERLING_UNLIKELY(g == NULL)) {
        UNDERLING_ERROR_NULL("failed to allocate space for grid",
                             UNDERLING_ENOMEM);
    }

    // Initialize values that we'll update as we proceed...
    g->n[0]    = n0;
    g->n[1]    = n1;
    g->n[2]    = n2;
    g->pA      = pA;
    g->pB      = pB;
    g->g_comm  = MPI_COMM_NULL;
    g->pA_comm = MPI_COMM_NULL;
    g->pB_comm = MPI_COMM_NULL;

    // Get number of processors in the communicator
    int nproc;
    if ((error = MPI_Comm_size(comm, &nproc))) {
        underling_grid_destroy(g);
        UNDERLING_MPICHKN(error);
    }

    // Create a balanced processor grid if not specified by pA, pB != 0
    dims[0] = g->pA;
    dims[1] = g->pB;
    if ((error = MPI_Dims_create(nproc, 2, dims))) {
        underling_grid_destroy(g);
        UNDERLING_MPICHKN(error);
    }

    // Specifying zero for both pA and pB "aligns" the resulting grid so that
    // the larger of pA and pB decomposes the larger of n0 and n2.
    if (g->pA == 0 && g->pB == 0) {
        if (dims[0] > dims[1]) {           // Ensure dims[0] < dims[1]
            const int tmp = dims[0]; dims[0] = dims[1]; dims[1] = tmp;
        }

        if (g->n[2] < g->n[0]) {
            // NOP: See assignment to g->p{A,B} just below
        } else {
            // Swap: See assignment to g->p{A,B} just below
            const int tmp = dims[0]; dims[0] = dims[1]; dims[1] = tmp;
        }
    }

    // Store new values of pA, pB within workspace
    g->pA = dims[0]; // Most closely associated with n2 decomposition
    g->pB = dims[1]; // Most closely associated with n0 decomposition

    // Sanity check decomposition against processor count
    if (UNDERLING_UNLIKELY(g->pA * g->pB != nproc)) {
        snprintf(buffer, sizeof(buffer),
                "Invalid processor grid: pA {%d} * pB {%d} != nproc {%d}",
                g->pA, g->pB, nproc);
        underling_grid_destroy(g);
        UNDERLING_ERROR_NULL(buffer, UNDERLING_EINVAL);
    }

    // Sanity check decomposition against grid size requirements. See Redmine
    // ticket #1297 for a discussion about handling degenerate transposes.
    // Eliminating this requirement will require updating the Doxygen.
    if (UNDERLING_UNLIKELY(g->n[0] < g->pB)) {
        snprintf(buffer, sizeof(buffer),
                "Decomposition requires n0 {%d} >= pB {%d}",
                g->n[0], g->pB);
        underling_grid_destroy(g);
        UNDERLING_ERROR_NULL(buffer, UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(g->n[1] < g->pA || g->n[1] < g->pB)) {
        snprintf(buffer, sizeof(buffer),
                "Decomposition requires n1 {%d} >= pA {%d}, pB {%d}",
                g->n[1], g->pA, g->pB);
        underling_grid_destroy(g);
        UNDERLING_ERROR_NULL(buffer, UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(g->n[0] < g->pB)) {
        snprintf(buffer, sizeof(buffer),
                "Decomposition requires n2 {%d} >= pA {%d}",
                g->n[2], g->pA);
        underling_grid_destroy(g);
        UNDERLING_ERROR_NULL(buffer, UNDERLING_EINVAL);
    }

    // Clone the communicator and create the 2D Cartesian topology
    {
        dims[0] = g->pA;
        dims[1] = g->pB;
        int periodic[2] = { 0, 0 };
        if ((error = MPI_Cart_create(comm, 2, dims, periodic, 1, &g->g_comm))) {
            underling_grid_destroy(g);
            UNDERLING_MPICHKN(error);
        }
    }

    // Create communicator for the pA direction
    dims[0] = 1;
    dims[1] = 0;
    if ((error = MPI_Cart_sub(g->g_comm, /* remaining */ dims, &g->pA_comm))) {
        underling_grid_destroy(g);
        UNDERLING_MPICHKN(error);
    }
    // Find the rank of this process within pA_comm
    int pA_rank;
    if ((error = MPI_Comm_rank(g->pA_comm, &pA_rank))) {
        underling_grid_destroy(g);
        UNDERLING_MPICHKN(error);
    }

    // Create communicator for the pB direction
    dims[0] = 0;
    dims[1] = 1;
    if ((error = MPI_Cart_sub(g->g_comm, /* remaining */ dims, &g->pB_comm))) {
        underling_grid_destroy(g);
        UNDERLING_MPICHKN(error);
    }
    // Find the rank of this process within pB_comm
    int pB_rank;
    if ((error = MPI_Comm_rank(g->pB_comm, &pB_rank))) {
        underling_grid_destroy(g);
        UNDERLING_MPICHKN(error);
    }

    // Name the three new communicators something mildly descriptive
    // Intra-communicator names can be lexicographically ordered easily
    const int pA_ndigits = snprintf(NULL, 0, "%d", g->pA);
    const int pB_ndigits = snprintf(NULL, 0, "%d", g->pB);
    assert(pA_ndigits > 0);
    assert(pB_ndigits > 0);
    snprintf(buffer, sizeof(buffer),
            "uGComm%0*dx%0*d", pA_ndigits, pA_rank, pB_ndigits, pB_rank);
    if ((error = MPI_Comm_set_name(g->g_comm, buffer))) {
        underling_grid_destroy(g);
        UNDERLING_MPICHKN(error);
    }
    snprintf(buffer, sizeof(buffer),
            "uPACommXx%0*d", pB_ndigits, pB_rank);
    if ((error = MPI_Comm_set_name(g->pA_comm, buffer))) {
        underling_grid_destroy(g);
        UNDERLING_MPICHKN(error);
    }
    snprintf(buffer, sizeof(buffer),
            "uPBComm%0*dxX", pA_ndigits, pA_rank);
    if ((error = MPI_Comm_set_name(g->pB_comm, buffer))) {
        underling_grid_destroy(g);
        UNDERLING_MPICHKN(error);
    }

    return g;
}

int
underling_grid_pA_size(
        const underling_grid grid)
{
    if (UNDERLING_UNLIKELY(grid == NULL)) {
        UNDERLING_ERROR_VAL("grid == NULL", UNDERLING_EINVAL, 0);
    }

    int retval;

    const int error = MPI_Comm_size(grid->pA_comm, &retval);
    if (error) {
        UNDERLING_MPICHKR(error /* MPI_Comm_size */);
        return 0;
    }

    return retval;
}

int
underling_grid_pB_size(
        const underling_grid grid)
{
    if (UNDERLING_UNLIKELY(grid == NULL)) {
        UNDERLING_ERROR_VAL("grid == NULL", UNDERLING_EINVAL, 0);
    }

    int retval;

    const int error = MPI_Comm_size(grid->pB_comm, &retval);
    if (error) {
        UNDERLING_MPICHKR(error /* MPI_Comm_size */);
        return 0;
    }

    return retval;
}

void
underling_grid_destroy(underling_grid grid)
{
    if (grid) {
        if (grid->g_comm != MPI_COMM_NULL) {
            UNDERLING_MPICHKR(MPI_Comm_free(&grid->g_comm));
            grid->g_comm = MPI_COMM_NULL;
        }
        if (grid->pA_comm != MPI_COMM_NULL) {
            UNDERLING_MPICHKR(MPI_Comm_free(&grid->pA_comm));
            grid->pA_comm = MPI_COMM_NULL;
        }
        if (grid->pB_comm != MPI_COMM_NULL) {
            UNDERLING_MPICHKR(MPI_Comm_free(&grid->pB_comm));
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
    if (UNDERLING_UNLIKELY(comm == MPI_COMM_NULL)) {
        UNDERLING_ERROR_NULL(
                "Unable to create transpose for comm == MPI_COMM_NULL",
                UNDERLING_EINVAL);
    }

    // Create and initialize the transpose workspace
    underling_transpose t = calloc(1, sizeof(struct underling_transpose_s));
    if (UNDERLING_UNLIKELY(t == NULL)) {
        UNDERLING_ERROR_NULL("failed to allocate space for transpose",
                             UNDERLING_ENOMEM);
    }

    // Fix struct values known from arguments
    t->comm     = underling_MPI_Comm_dup_with_name(comm);
    t->d[0]     = d0;
    t->d[1]     = d1;
    t->howmany  = howmany;
    t->block[0] = block0;
    t->block[1] = block1;
    t->flags    = flags;

    if (t->comm == MPI_COMM_NULL) {
        underling_transpose_destroy(t);
        UNDERLING_ERROR_NULL("Detected MPI_COMM_NULL in t->comm",
                            UNDERLING_ESANITY);
    }

    if (UNDERLING_UNLIKELY(   t->howmany == 0
                           || t->d[0] == 0
                           || t->d[1] == 0)) {
        // Trivial transpose required;
        // fftw_mpi_local_size_many_transposed divides-by-zero on zero size
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
    if (UNDERLING_UNLIKELY(forward == NULL)) {
        UNDERLING_ERROR_NULL("forward == NULL", UNDERLING_EINVAL);
    }
    const unsigned non_transposed_inout_mask = ~FFTW_MPI_TRANSPOSED_IN
                                             & ~FFTW_MPI_TRANSPOSED_OUT;
    if (UNDERLING_UNLIKELY(forward->flags & non_transposed_inout_mask)) {
        UNDERLING_ERROR_NULL(
                "Flags contains non-FFTW_MPI_TRANSPOSED_{IN,OUT}",
                UNDERLING_ESANITY);
    }

    // Create and initialize the transpose workspace
    underling_transpose backward
        = calloc(1, sizeof(struct underling_transpose_s));
    if (UNDERLING_UNLIKELY(backward == NULL)) {
        UNDERLING_ERROR_NULL("failed to allocate space for transpose",
                             UNDERLING_ENOMEM);
    }

    // Fix struct values known from inverting forward plan
    backward->d[0]     = forward->d[1];
    backward->d[1]     = forward->d[0];
    backward->howmany  = forward->howmany;
    backward->block[0] = forward->block[1];
    backward->block[1] = forward->block[0];
    backward->comm     = underling_MPI_Comm_dup_with_name(forward->comm);
    backward->flags    = 0;

    // Handle transposed input and output flags
    if (forward->flags & FFTW_MPI_TRANSPOSED_IN) {
        backward->flags |= FFTW_MPI_TRANSPOSED_OUT;
    }
    if (forward->flags & FFTW_MPI_TRANSPOSED_OUT) {
        backward->flags |= FFTW_MPI_TRANSPOSED_IN;
    }

    if (backward->comm == MPI_COMM_NULL) {
        underling_transpose_destroy(backward);
        UNDERLING_ERROR_NULL("Detected MPI_COMM_NULL in backward->comm",
                            UNDERLING_ESANITY);
    }

    if (UNDERLING_UNLIKELY(   backward->howmany == 0
                           || backward->d[0]    == 0
                           || backward->d[1]    == 0)) {
        // Trivial transpose required;
        // fftw_mpi_local_size_many_transposed divides-by-zero on zero size
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
            UNDERLING_MPICHKR(MPI_Comm_free(&transpose->comm));
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
    if (UNDERLING_UNLIKELY(transpose == NULL)) {
        UNDERLING_ERROR_NULL("transpose == NULL", UNDERLING_EINVAL);
    }

    if (UNDERLING_UNLIKELY(   transpose->howmany == 0
                           || transpose->d[0]    == 0
                           || transpose->d[1]    == 0)) {
        // fftw_mpi_plan_many_transpose returns NULL on some trivial input.
        // Create a NOP fftw_plan and return it to sidestep the issue.
        return fftw_plan_guru_r2r(0, NULL, 0, NULL, NULL, NULL, NULL, 0);
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

static
void
underling_dump_transposes(
        MPI_Comm dump_comm,
        FILE *out,
        const char *prefix,
        const underling_transpose transpose)
{
    int size, rank;
    UNDERLING_MPICHKV(MPI_Comm_size(dump_comm, &size));
    UNDERLING_MPICHKV(MPI_Comm_rank(dump_comm, &rank));

    UNDERLING_MPICHKV(MPI_Barrier(dump_comm));
    fflush(out);
    for (int i = 0; i < size; ++i) {
        UNDERLING_MPICHKV(MPI_Barrier(dump_comm));
        if (i == rank) {
            fputs(prefix, out);
            fputc(' ', out);
            underling_fprint_transpose(transpose, out);
            fputs("\n", out);
            fflush(out);
        }
    }
    fflush(out);
    UNDERLING_MPICHKV(MPI_Barrier(dump_comm));
}

// FIXME: Leaks memory on pretty much any MPI-based error
underling_problem
underling_problem_create(
        underling_grid grid,
        int howmany,
        unsigned transposed_flags)
{
    if (UNDERLING_UNLIKELY(grid == NULL)) {
        UNDERLING_ERROR_NULL("grid == NULL", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(howmany < 0)) {
        UNDERLING_ERROR_NULL("howmany >= 0 required", UNDERLING_EINVAL);
    }
    const unsigned non_transposed_mask =   ~UNDERLING_TRANSPOSED_LONG_N2
                                         & ~UNDERLING_TRANSPOSED_LONG_N0;
    if (UNDERLING_UNLIKELY(transposed_flags & non_transposed_mask)) {
        UNDERLING_ERROR_NULL(
            "transposed_flags contains non-UNDERLING_TRANSPOSED_LONG_N{0,2}",
            UNDERLING_EINVAL);
    }

    // Create and initialize the problem workspace
    underling_problem p = calloc(1, sizeof(struct underling_problem_s));
    if (UNDERLING_UNLIKELY(p == NULL)) {
        UNDERLING_ERROR_NULL("failed to allocate space for problem",
                             UNDERLING_ENOMEM);
    }
    // Copy the problem parameters to the problem workspace
    p->howmany = howmany;

    // Global pencil decomposition details
    // assuming transposed_flags == 0
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

    // Fix interleaved data field details in p->long_{n2,n1,n0}
    for (int i = 0; i < 3; ++i) {
        p->long_n[i].start[3] = 0;
        p->long_n[i].size[3]  = p->howmany;
    }

    // underling_grid_create() already enforced that n0 >= pB, n1 >= pA, pB,
    // and n2 >= pA.  See Redmine ticket #1297 for a discussion about handling
    // degenerate transposes.

    // Decompose {n0,n1}/pB and store details in p->long_{(n2,n1),n0}
    {
        ptrdiff_t local_d0, local_d0_start, local_d1, local_d1_start;
        ptrdiff_t dB[2] = {grid->n[0], grid->n[1]}; // Never performed
        const int howmany_tmp                       // Allows howmany == 0
            = UNDERLING_UNLIKELY(p->howmany == 0) ? 1 : p->howmany;
        fftw_mpi_local_size_many_transposed(2, dB, howmany_tmp,
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
        const int howmany_tmp                       // Allow howmany == 0
            = UNDERLING_UNLIKELY(p->howmany == 0) ? 1 : p->howmany;
        fftw_mpi_local_size_many_transposed(2, dA, howmany_tmp,
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

    // Determine all necessary stride orders for row-major storage
    // -----------------------------------------------------------
    if (transposed_flags & UNDERLING_TRANSPOSED_LONG_N2) {
        // Set stride order when long in n2:   n2 x (n0/pB x n1/pA) x howmany
        // Confirmed ordering by email to fftw@fftw.org on 30 Nov 2011
        p->long_n[2].order[0] = 3; // Fastest, interleaved data
        p->long_n[2].order[1] = 1;
        p->long_n[2].order[2] = 0;
        p->long_n[2].order[3] = 2; // Slowest, long direction
    } else {
        // Set stride order when long in n2: (n0/pB x n1/pA) x (n2 x howmany)
        p->long_n[2].order[0] = 3; // Fastest, interleaved data
        p->long_n[2].order[1] = 2; // Long direction
        p->long_n[2].order[2] = 1;
        p->long_n[2].order[3] = 0; // Slowest
    }
    // Set stride order when long in n1: (n2/pA x n0/pB) x (n1 x howmany)
    p->long_n[1].order[0] = 3; // Fastest, interleaved data
    p->long_n[1].order[1] = 1; // Long direction
    p->long_n[1].order[2] = 0;
    p->long_n[1].order[3] = 2; // Slowest
    if (transposed_flags & UNDERLING_TRANSPOSED_LONG_N0) {
        // Set stride order when long in n0:   n0 x (n1/pB x n2/pA) x howmany
        // Confirmed ordering by email to fftw@fftw.org on 30 Nov 2011
        p->long_n[0].order[0] = 3; // Fastest, interleaved data
        p->long_n[0].order[1] = 2;
        p->long_n[0].order[2] = 1;
        p->long_n[0].order[3] = 0; // Slowest, long direction
    } else {
        // Set stride order when long in n0: (n1/pB x n2/pA) x (n0 x howmany)
        p->long_n[0].order[0] = 3; // Fastest, interleaved data
        p->long_n[0].order[1] = 0; // Long direction
        p->long_n[0].order[2] = 2;
        p->long_n[0].order[3] = 1; // Slowest
    }

    // Use stride ordering to compute strides in each configuration
    // ------------------------------------------------------------
    for (int i = 0; i < 3; ++i) {
        underling_extents * const e = &p->long_n[i];
        e->stride[e->order[0]] = 1;
        e->stride[e->order[1]] = e->stride[e->order[0]] * e->size[e->order[0]];
        e->stride[e->order[2]] = e->stride[e->order[1]] * e->size[e->order[1]];
        e->stride[e->order[3]] = e->stride[e->order[2]] * e->size[e->order[2]];

        // Compute extent; redundant information but very convenient
        e->extent = e->size[0] * e->size[1] * e->size[2] * e->size[3];
    }

    // Transpose pA details: (n0/pB x n1/pA) x n2 to n2/pA x (n0/pB x n1)
    const ptrdiff_t pA_d[2] = { p->long_n[2].size[0] * grid->n[1],
                                grid->n[2] };
    ptrdiff_t pA_block[2]   = { p->long_n[2].size[0] * p->long_n[2].size[1],
                                p->long_n[1].size[2] };
    UNDERLING_MPICHKN(MPI_Bcast(pA_block, 2, MPI_LONG, 0, grid->pA_comm));

    // Transpose pB details: (n2/pA x n0/pB) x n1 to n1/pB x (n2/pA x n0)
    const ptrdiff_t pB_d[2] = { p->long_n[1].size[2] * grid->n[0],
                                grid->n[1] };
    ptrdiff_t pB_block[2]   = { p->long_n[1].size[2] * p->long_n[1].size[0],
                                p->long_n[0].size[1] };
    UNDERLING_MPICHKN(MPI_Bcast(pB_block, 2, MPI_LONG, 0, grid->pB_comm));

    // Wave towards physical MPI transpose: long in n2 to long in n1
    const unsigned backwardA_flags
        = (transposed_flags & UNDERLING_TRANSPOSED_LONG_N2)
        ? FFTW_MPI_TRANSPOSED_IN
        : 0;
    p->backwardA = underling_transpose_create(pA_d[0],
                                              pA_d[1],
                                              p->howmany,
                                              pA_block[0],
                                              pA_block[1],
                                              grid->pA_comm,
                                              backwardA_flags);
    if (UNDERLING_UNLIKELY(p->backwardA == NULL)) {
        underling_problem_destroy(p);
        UNDERLING_ERROR_NULL("failed creating p->backwardA",
                UNDERLING_EFAILED);
    }

    // Wave towards physical MPI transpose: long in n1 to long in n0
    const unsigned backwardB_flags
        = (transposed_flags & UNDERLING_TRANSPOSED_LONG_N0)
        ? FFTW_MPI_TRANSPOSED_OUT
        : 0;
    p->backwardB = underling_transpose_create(pB_d[0],
                                              pB_d[1],
                                              p->howmany,
                                              pB_block[0],
                                              pB_block[1],
                                              grid->pB_comm,
                                              backwardB_flags);
    if (UNDERLING_UNLIKELY(p->backwardB == NULL)) {
        underling_problem_destroy(p);
        UNDERLING_ERROR_NULL("failed creating p->backwardB",
                UNDERLING_EFAILED);
    }

    // Physical towards wave MPI transpose: long in n0 to long in n1
    p->forwardB = underling_transpose_create_inverse(p->backwardB);
    if (UNDERLING_UNLIKELY(p->forwardB == NULL)) {
        underling_problem_destroy(p);
        UNDERLING_ERROR_NULL("failed creating p->forwardB",
                UNDERLING_EFAILED);
    }

    // Physical towards wave MPI transpose: long in n1 to long in n2
    p->forwardA = underling_transpose_create_inverse(p->backwardA);
    if (UNDERLING_UNLIKELY(p->forwardA == NULL)) {
        underling_problem_destroy(p);
        UNDERLING_ERROR_NULL("failed creating p->forwardA",
                UNDERLING_EFAILED);
    }

    // p->local_memory is overall maximum of all local_size values
    p->local_memory = p->backwardA->local_size;
    if (p->local_memory < (size_t) p->backwardB->local_size) {
        p->local_memory = p->backwardB->local_size;
    }
    if (p->local_memory < (size_t) p->forwardB->local_size) {
        p->local_memory = p->forwardB->local_size;
    }
    if (p->local_memory < (size_t) p->forwardA->local_size) {
        p->local_memory = p->forwardA->local_size;
    }

    return p;
}

size_t
underling_local_memory(
        const underling_problem problem)
{
    if (UNDERLING_UNLIKELY(problem == NULL)) {
        UNDERLING_ERROR_VAL("problem == NULL", UNDERLING_EINVAL, 0);
    }

    return problem->local_memory;
}

size_t
underling_local_memory_optimum(
        const underling_problem problem)
{
    // Optimum local memory is the maximum required be long in any direction

    size_t retval = problem->long_n[0].extent;
    if (retval < problem->long_n[1].extent) {
        retval = problem->long_n[1].extent;
    }
    if (retval < problem->long_n[2].extent) {
        retval = problem->long_n[2].extent;
    }

    return retval;
}

static
size_t
underling_local_memory_allreduce(
        const underling_grid    grid,
        const underling_problem problem,
        MPI_Op op)
{
    if (UNDERLING_UNLIKELY(grid == NULL)) {
        UNDERLING_ERROR_VAL("grid == NULL", UNDERLING_EINVAL, 0);
    }
    if (UNDERLING_UNLIKELY(problem == NULL)) {
        UNDERLING_ERROR_VAL("problem == NULL", UNDERLING_EINVAL, 0);
    }

    // Use unsigned long values for safety in heterogeneous environments.
    assert(((size_t)-1) <= ((unsigned long int)-1)); // "Type safety"

    unsigned long int retval;

    unsigned long int sendbuf = underling_local_memory(problem);
    const int error = MPI_Allreduce(
            &sendbuf, &retval, 1, MPI_UNSIGNED_LONG,
            op, grid->g_comm);
    if (UNDERLING_UNLIKELY(error)) {
        UNDERLING_MPICHKR(error /* allreduce local_memory */);
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
    if (UNDERLING_UNLIKELY(grid == NULL)) {
        UNDERLING_ERROR_VAL("grid == NULL", UNDERLING_EINVAL, 0);
    }
    if (UNDERLING_UNLIKELY(problem == NULL)) {
        UNDERLING_ERROR_VAL("problem == NULL", UNDERLING_EINVAL, 0);
    }
    return problem->howmany * grid->n[0] * grid->n[1] * grid->n[2];
}

size_t
underling_local(
        const underling_problem problem,
        int i,
        int *start,
        int *size,
        int *stride,
        int *order)
{
    if (UNDERLING_UNLIKELY(i < 0 || i > 2)) {
        UNDERLING_ERROR_VAL("i < 0 or i > 2", UNDERLING_EINVAL, 0);
    }
    if (UNDERLING_UNLIKELY(problem == NULL)) {
        UNDERLING_ERROR_VAL("problem == NULL", UNDERLING_EINVAL, 0);
    }

    const underling_extents * const e = &problem->long_n[i];

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

    return e->extent;
}

underling_extents
underling_local_extents(
        const underling_problem problem,
        int i)
{
    if (UNDERLING_UNLIKELY(i < 0 || i > 2)) {
        UNDERLING_ERROR_VAL("i < 0 or i > 2",
                UNDERLING_EINVAL, UNDERLING_EXTENTS_INVALID);
    }
    if (UNDERLING_UNLIKELY(problem == NULL)) {
        UNDERLING_ERROR_VAL("problem == NULL",
                UNDERLING_EINVAL, UNDERLING_EXTENTS_INVALID);
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
        underling_real * in,
        underling_real * out,
        unsigned transpose_flags,
        unsigned rigor_flags)
{
    if (UNDERLING_UNLIKELY(problem == NULL)) {
        UNDERLING_ERROR_NULL("problem == NULL", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(in == NULL)) {
        UNDERLING_ERROR_NULL("in == NULL", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(out == NULL)) {
        UNDERLING_ERROR_NULL("out == NULL", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(transpose_flags
                & ~(UNDERLING_TRANSPOSE_ALL | UNDERLING_TRANSPOSE_NONE))) {
        UNDERLING_ERROR_NULL(
            "transpose_flags contains non-direction bit", UNDERLING_EINVAL);
    }
    const unsigned non_rigor_mask =   ~FFTW_ESTIMATE
                                    & ~FFTW_MEASURE
                                    & ~FFTW_PATIENT
                                    & ~FFTW_EXHAUSTIVE
                                    & ~FFTW_WISDOM_ONLY;
    if (UNDERLING_UNLIKELY(rigor_flags & non_rigor_mask)) {
        UNDERLING_ERROR_NULL(
                "FFTW non-rigor bits disallowed", UNDERLING_EINVAL);
    }

    // Be ready to execute all transforms if trivial flag provided
    if (transpose_flags == 0) {
        transpose_flags = UNDERLING_TRANSPOSE_ALL;
    }
    // Disable all transpose directions whenever transpose none is set
    if (transpose_flags & UNDERLING_TRANSPOSE_NONE) {
        transpose_flags &= ~UNDERLING_TRANSPOSE_ALL;
    }

    // Create and initialize the plan workspace
    underling_plan p = calloc(1, sizeof(struct underling_plan_s));
    if (UNDERLING_UNLIKELY(p == NULL)) {
        UNDERLING_ERROR_NULL("failed to allocate space for plan",
                             UNDERLING_ENOMEM);
    }
    // Copy the problem parameters to the problem workspace
    p->in_place = (in == out);

    // Deadlock debugging logic which can be selectively enabled at runtime
    if (debug_comm != MPI_COMM_NULL) {
        underling_dump_transposes(debug_comm, stdout,
                "UNDERLING_DEBUG backwardA", problem->backwardA);
    }

    // Create the requested FFTW MPI plans
    if (transpose_flags | UNDERLING_TRANSPOSE_LONG_N2_TO_LONG_N1) {
        p->plan_backwardA = underling_transpose_fftw_plan(
                problem->backwardA, in, out, rigor_flags | destroy_flags);
        if (UNDERLING_UNLIKELY(p->plan_backwardA == NULL)) {
            underling_plan_destroy(p);
            UNDERLING_ERROR_NULL(
                    "FFTW MPI returned NULL plan: plan_backwardA",
                    UNDERLING_EFAILED);
        }
    }

    // Deadlock debugging logic which can be selectively enabled at runtime
    if (debug_comm != MPI_COMM_NULL) {
        underling_dump_transposes(debug_comm, stdout,
                "UNDERLING_DEBUG backwardB", problem->backwardB);
    }

    if (transpose_flags | UNDERLING_TRANSPOSE_LONG_N1_TO_LONG_N0) {
        p->plan_backwardB = underling_transpose_fftw_plan(
                problem->backwardB, in, out, rigor_flags | destroy_flags);
        if (UNDERLING_UNLIKELY(p->plan_backwardB == NULL)) {
            underling_plan_destroy(p);
            UNDERLING_ERROR_NULL(
                    "FFTW MPI returned NULL plan: plan_backwardB",
                    UNDERLING_EFAILED);
        }
    }

    // Deadlock debugging logic which can be selectively enabled at runtime
    if (debug_comm != MPI_COMM_NULL) {
        underling_dump_transposes(debug_comm, stdout,
                "UNDERLING_DEBUG forwardB", problem->forwardB);
    }

    if (transpose_flags | UNDERLING_TRANSPOSE_LONG_N0_TO_LONG_N1) {
        p->plan_forwardB = underling_transpose_fftw_plan(
                problem->forwardB, in, out, rigor_flags | destroy_flags);
        if (UNDERLING_UNLIKELY(p->plan_forwardB == NULL)) {
            underling_plan_destroy(p);
            UNDERLING_ERROR_NULL(
                    "FFTW MPI returned NULL plan: plan_forwardB",
                    UNDERLING_EFAILED);
        }
    }

    // Deadlock debugging logic which can be selectively enabled at runtime
    if (debug_comm != MPI_COMM_NULL) {
        underling_dump_transposes(debug_comm, stdout,
                "UNDERLING_DEBUG forwardA", problem->forwardA);
    }

    if (transpose_flags | UNDERLING_TRANSPOSE_LONG_N1_TO_LONG_N2) {
        p->plan_forwardA = underling_transpose_fftw_plan(
                problem->forwardA, in, out, rigor_flags | destroy_flags);
        if (UNDERLING_UNLIKELY(p->plan_forwardA == NULL)) {
            underling_plan_destroy(p);
            UNDERLING_ERROR_NULL(
                    "FFTW MPI returned NULL plan: plan_forwardA",
                    UNDERLING_EFAILED);
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
        free(plan);
    }
}

int
underling_execute_long_n2_to_long_n1(
        const underling_plan plan,
        underling_real * in,
        underling_real * out)
{
    if (UNDERLING_UNLIKELY(plan == NULL)) {
        UNDERLING_ERROR("plan == NULL", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(plan->plan_backwardA == NULL)) {
        UNDERLING_ERROR("plan->plan_backwardA == NULL", UNDERLING_EINVAL);
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

#ifdef HAVE_FFTW33_BETA1_OR_LATER
    fftw_mpi_execute_r2r(plan->plan_backwardA, in, out);
#else
    fftw_execute_r2r(plan->plan_backwardA, in, out);
#endif

    return UNDERLING_SUCCESS;
}

int
underling_execute_long_n1_to_long_n0(
        const underling_plan plan,
        underling_real * in,
        underling_real * out)
{
    if (UNDERLING_UNLIKELY(plan == NULL)) {
        UNDERLING_ERROR("plan == NULL", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(plan->plan_backwardB == NULL)) {
        UNDERLING_ERROR("plan->plan_backwardB == NULL", UNDERLING_EINVAL);
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

#ifdef HAVE_FFTW33_BETA1_OR_LATER
    fftw_mpi_execute_r2r(plan->plan_backwardB, in, out);
#else
    fftw_execute_r2r(plan->plan_backwardB, in, out);
#endif

    return UNDERLING_SUCCESS;
}

int
underling_execute_long_n0_to_long_n1(
        const underling_plan plan,
        underling_real * in,
        underling_real * out)
{
    if (UNDERLING_UNLIKELY(plan == NULL)) {
        UNDERLING_ERROR("plan == NULL", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(plan->plan_forwardB == NULL)) {
        UNDERLING_ERROR("plan->plan_forwardB == NULL", UNDERLING_EINVAL);
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

#ifdef HAVE_FFTW33_BETA1_OR_LATER
    fftw_mpi_execute_r2r(plan->plan_forwardB, in, out);
#else
    fftw_execute_r2r(plan->plan_forwardB, in, out);
#endif

    return UNDERLING_SUCCESS;
}

int
underling_execute_long_n1_to_long_n2(
        const underling_plan plan,
        underling_real * in,
        underling_real * out)
{
    if (UNDERLING_UNLIKELY(plan == NULL)) {
        UNDERLING_ERROR("plan == NULL", UNDERLING_EINVAL);
    }
    if (UNDERLING_UNLIKELY(plan->plan_forwardA == NULL)) {
        UNDERLING_ERROR("plan->plan_forwardA == NULL", UNDERLING_EINVAL);
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

#ifdef HAVE_FFTW33_BETA1_OR_LATER
    fftw_mpi_execute_r2r(plan->plan_forwardA, in, out);
#else
    fftw_execute_r2r(plan->plan_forwardA, in, out);
#endif

    return UNDERLING_SUCCESS;
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
        char buffer[MPI_MAX_OBJECT_NAME];
        int resultlen = 0;
        if (   MPI_Comm_get_name(transpose->comm, buffer, &resultlen)
            || resultlen == 0) {
            fprintf(output_file, "comm=%x", transpose->comm);
        } else {
            fprintf(output_file, "comm=%s", buffer);
        }
        fprintf(output_file,
                ",{d=%td,%td},howmany=%td,{block=%td,%td}",
                transpose->d[0], transpose->d[1], transpose->howmany,
                transpose->block[0], transpose->block[1]);
        fprintf(output_file,
                ",flags=%u,local_size=%td",
                transpose->flags,transpose->local_size);
        fprintf(output_file,
                ",{local=%td,%td},{local_start=%td,%td}",
                transpose->local[0],transpose->local[1],
                transpose->local_start[0],transpose->local_start[1]);
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
        fprintf(output_file,"{howmany=%d,local_memory=%zu}",
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
        fprintf(output_file, "extents=[%d,%d)x[%d,%d)x[%d,%d)x[%d,%d)",
                extents->start[0],
                extents->start[0] + extents->size[0],
                extents->start[1],
                extents->start[1] + extents->size[1],
                extents->start[2],
                extents->start[2] + extents->size[2],
                extents->start[3],
                extents->start[3] + extents->size[3]);
        fprintf(output_file, ",strides={%d,%d,%d,%d}",
                extents->stride[0],
                extents->stride[1],
                extents->stride[2],
                extents->stride[3]);
    }
}

MPI_Comm
underling_debug_transpose(
        MPI_Comm comm)
{
    MPI_Comm last = debug_comm;
    debug_comm = comm;
    return last;
}
