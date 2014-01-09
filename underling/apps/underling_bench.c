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

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sysexits.h>
#include <sys/file.h>
#include <tgmath.h>
#include <unistd.h>

#include "argp.h"
#include "minmax.h"
#include "mpi_argp.h"

#include <mpi.h>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <underling/error.h>
#include <underling/underling.h>
#include <underling/underling_fftw.h>

#ifdef HAVE_GRVY
#include <grvy.h>
#define GRVY_TIMER_BEGIN(id)   grvy_timer_begin(id)
#define GRVY_TIMER_END(id)     grvy_timer_end(id)
#define GRVY_TIMER_FINALIZE()  grvy_timer_finalize()
#define GRVY_TIMER_INIT(id)    grvy_timer_init(id)
#define GRVY_TIMER_RESET()     grvy_timer_reset()
#define GRVY_TIMER_SUMMARIZE() grvy_timer_summarize()
#else
#define GRVY_TIMER_BEGIN(id)
#define GRVY_TIMER_END(id)
#define GRVY_TIMER_FINALIZE()
#define GRVY_TIMER_INIT(id)
#define GRVY_TIMER_RESET()
#define GRVY_TIMER_SUMMARIZE()
#endif

//****************************************************************
// DATA STRUCTURES DATA STRUCTURES DATA STRUCTURES DATA STRUCTURES
//****************************************************************

struct details {
    int world_rank;
    int world_size;
    int verbose;
    int repeat;
    int nfields;
    int howmany;
    int nthreads;
    int n0;
    int n1;
    int n2;
    int pA;
    int pB;
    long bytes;
    unsigned transposed_flags;
    unsigned transform_flags;
    unsigned fftw_rigor_flags;
    unsigned packed_flags;
    char *wisdom_file;
    char fft_n[3];
    int mpi_inplace;
    int fft_inplace;
    int forward;
    int backward;
    double abstol;
};

//*************************************************************
// STATIC PROTOTYPES STATIC PROTOTYPES STATIC PROTOTYPES STATIC
//*************************************************************

static void print_version(FILE *stream, struct argp_state *state);

static void trim(char *a);

static inline int min(int a, int b);

static inline int max(int a, int b);

static void to_human_readable_byte_count(long bytes,
                                         int si,
                                         double *coeff,
                                         const char **units);

static long from_human_readable_byte_count(const char *str);

static underling_real* calloc_field(struct details *d,
                                    const long bytes);

static void fill_field(underling_real *p,
                       struct details *d,
                       const long bytes,
                       unsigned salt);

static underling_real check_field(const underling_real *p,
                                  struct details *d,
                                  const long bytes,
                                  unsigned salt,
                                  underling_fftw_plan plan);

//*******************************************************************
// ARGP DETAILS: http://www.gnu.org/s/libc/manual/html_node/Argp.html
//*******************************************************************

const char *argp_program_version      = "underling_bench " PACKAGE_VERSION;
void (*argp_program_version_hook)(FILE *stream, struct argp_state *state)
                                      = &print_version;
const char *argp_program_bug_address  = PACKAGE_BUGREPORT;
static const char doc[]               =
"Simulate and benchmark underling-based application transformation operations."
"\v"
"Transform parallel, 3D pencil decompositions using underling's data "
"movement capabilities.  Timing information is collected "
"over one or more iterations.  If provided, FFTW wisdom is imported from "
"and accumulated within WISDOM_FILE.\n"
"\n"
"Options taking a 'bytes' parameter can be given common byte-related "
"units.  For example --field-memory=5G indicates that approximately "
"5 gigabytes of memory should be used on each rank to store field data. "
" SI units like 'Ki' or 'MiB' are also accepted.\n"
;

static const char args_doc[] = "\nWISDOM_FILE";

enum {
    KEY_ESTIMATE = 1024, // !isascii
    KEY_MEASURE,
    KEY_PATIENT,
    KEY_EXHAUSTIVE,
    KEY_WISDOM_ONLY,
    KEY_DIR_FORWARD,
    KEY_DIR_BACKWARD,
    KEY_DIR_BOTH,
    KEY_TRANS_N2,
    KEY_TRANS_N0,
    KEY_TRANS_ALL,
    KEY_FFT_IIC,
    KEY_FFT_IIR,
    KEY_FFT_ICC,
    KEY_FFT_ICR,
    KEY_FFT_CCC,
    KEY_FFT_CCR,
    KEY_PACK_N2,
    KEY_PACK_N0,
    KEY_PACK_ALL
};

static struct argp_option options[] = {
    {"check",        'c', "abstol", 0, "absolute round trip tolerance to scaled by transform extents", 0},
    {"verbose",      'v', 0,        0, "produce verbose output",                                       0},
    {"repeat",       'r', "count",  0, "number of repetitions",                                        0},
    {"nthreads",     't', "count",  0, "number of concurrent threads",                                 0},
    {"nfields",      'n', "count",  0, "number of independent fields",                                 0},
    {"howmany",      'h', "count",  0, "howmany components per field",                                 0},
    {"mpi-in-place", 'i', 0,        0, "perform in-place MPI transposes",                              0},
    {"fft-in-place", 'I', 0,        0, "perform in-place FFTs",                                        0},
    {0, 0, 0, 0,
     "Controlling global problem size (specify at most one)",     0 },
    {"field-memory", 'f', "bytes",    0, "per-rank field memory", 0 },
    {"field-global", 'F', "n0xn1xn2", 0, "field global extents",  0 },
    {0, 0, 0, 0,
     "Controlling parallel decomposition per MPI_Dims_create semantics", 0 },
    {"dims",       'P', "pAxpB",    0, "process grid for decomposition", 0 },
    {0, 0, 0, 0,
     "Controlling field transpose and transform directionality", 0},
    {"forward",    KEY_DIR_FORWARD,  0, 0, "go from long_n2 to long_n0", 0},
    {"backward",   KEY_DIR_BACKWARD, 0, 0, "go from long_n0 to long_n2", 0},
    {"both",       KEY_DIR_BOTH,     0, 0, "perform forward then backward (default)", 0},
    {0, 0, 0, 0,
     "Reducing on-node overhead by modifying storage ordering", 0},
    {"trans",    KEY_TRANS_ALL, 0, 0, "use flags TRANSPOSED_LONG_N{0,2}", 0},
    {"trans0",   KEY_TRANS_N0,  0, 0, "use flag  TRANSPOSED_LONG_N0",     0},
    {"trans2",   KEY_TRANS_N2,  0, 0, "use flag  TRANSPOSED_LONG_N2",     0},
    {0, 0, 0, 0,
     "Adding Fourier transformations (--forward shown, --backward inverts)", 0},
    {"ccc", KEY_FFT_CCC, 0, 0, "c2c long_n0, c2c long_n1, c2c long_n2", 0},
    {"ccr", KEY_FFT_CCR, 0, 0, "c2c long_n0, c2c long_n1, r2c long_n2", 0},
    {"icc", KEY_FFT_ICC, 0, 0, "             c2c long_n1, c2c long_n2", 0},
    {"icr", KEY_FFT_ICR, 0, 0, "             c2c long_n1, r2c long_n2", 0},
    {"iic", KEY_FFT_IIC, 0, 0, "                          c2c long_n2", 0},
    {"iir", KEY_FFT_IIR, 0, 0, "                          r2c long_n2", 0},
    {0, 0, 0, 0,
     "Forcing packed, contiguous output from Fourier transforms", 0},
    {"pack",    KEY_PACK_ALL, 0, 0, "use flag FFTW_PACKED_ALL",     0},
    {"pack0",   KEY_PACK_N0,  0, 0, "use flag FFTW_PACKED_LONG_N0", 0},
    {"pack2",   KEY_PACK_N2,  0, 0, "use flag FFTW_PACKED_LONG_N2", 0},
    {0, 0, 0, 0,
     "Changing FFTW planning rigor", 0},
    {"estimate",    KEY_ESTIMATE,    0, 0, "plan with FFTW_ESTIMATE", 0},
    {"measure",     KEY_MEASURE,     0, 0, "plan with FFTW_MEASURE (default)", 0},
    {"patient",     KEY_PATIENT,     0, 0, "plan with FFTW_PATIENT", 0},
    {"exhaustive",  KEY_EXHAUSTIVE,  0, 0, "plan with FFTW_EXHAUSTIVE", 0},
    {"wisdom-only", KEY_WISDOM_ONLY, 0, 0, "plan with FFTW_WISDOM_ONLY", 0},
    {"timelimit",   'T', "seconds", 0, "use fftw_set_timelimit(seconds)", 0},
    { 0, 0, 0, 0, 0, 0 }
};

// Parse a single option following Argp semantics
static error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
    // Get the input argument from argp_parse.
    struct details *d = state->input;

    // Trim any leading/trailing whitespace from arg
    if (arg) trim(arg);

    // Want to ensure we consume the entire argument for many options
    // Many sscanf calls provide an extra sentinel %c dumping into &ignore
    char ignore = '\0';

    switch (key) {
        case ARGP_KEY_ARG:
            if (state->arg_num > 0) {
                argp_usage(state);
            } else {
                d->wisdom_file = strdup(arg);
                assert(d->wisdom_file);
            }
            break;

        case ARGP_KEY_END:
            if (state->arg_num > 1) {
                argp_usage(state);
            }
            break;

        case KEY_ESTIMATE:
            d->fftw_rigor_flags = FFTW_ESTIMATE;
            break;

        case KEY_MEASURE:
            d->fftw_rigor_flags = FFTW_MEASURE;
            break;

        case KEY_PATIENT:
            d->fftw_rigor_flags = FFTW_PATIENT;
            break;

        case KEY_EXHAUSTIVE:
            d->fftw_rigor_flags = FFTW_EXHAUSTIVE;
            break;

        case KEY_WISDOM_ONLY:
            d->fftw_rigor_flags = FFTW_WISDOM_ONLY;
            break;

        case KEY_DIR_FORWARD:
            d->forward  = 1;
            d->backward = 0;
            break;

        case KEY_DIR_BACKWARD:
            d->forward  = 0;
            d->backward = 1;
            break;

        case KEY_DIR_BOTH:
            d->forward  = 1;
            d->backward = 1;
            break;

        case KEY_TRANS_N2:
            d->transposed_flags |= UNDERLING_TRANSPOSED_LONG_N2;
            break;

        case KEY_TRANS_N0:
            d->transposed_flags |= UNDERLING_TRANSPOSED_LONG_N0;
            break;

        case KEY_TRANS_ALL:
            d->transposed_flags |= UNDERLING_TRANSPOSED_LONG_N2;
            d->transposed_flags |= UNDERLING_TRANSPOSED_LONG_N0;
            break;

        case KEY_FFT_IIC:
            d->fft_n[0] = 'i';
            d->fft_n[1] = 'i';
            d->fft_n[2] = 'c';
            break;

        case KEY_FFT_IIR:
            d->fft_n[0] = 'i';
            d->fft_n[1] = 'i';
            d->fft_n[2] = 'r';
            break;

        case KEY_FFT_ICC:
            d->fft_n[0] = 'i';
            d->fft_n[1] = 'c';
            d->fft_n[2] = 'c';
            break;

        case KEY_FFT_ICR:
            d->fft_n[0] = 'i';
            d->fft_n[1] = 'c';
            d->fft_n[2] = 'r';
            break;

        case KEY_FFT_CCC:
            d->fft_n[0] = 'c';
            d->fft_n[1] = 'c';
            d->fft_n[2] = 'c';
            break;

        case KEY_FFT_CCR:
            d->fft_n[0] = 'c';
            d->fft_n[1] = 'c';
            d->fft_n[2] = 'r';
            break;

        case KEY_PACK_N2:
            d->packed_flags |= UNDERLING_FFTW_PACKED_LONG_N2;
            break;

        case KEY_PACK_N0:
            d->packed_flags |= UNDERLING_FFTW_PACKED_LONG_N0;
            break;

        case KEY_PACK_ALL:
            d->packed_flags |= UNDERLING_FFTW_PACKED_ALL;
            break;

        case 'c':
            errno = 0;
            if (1 != sscanf(arg ? arg : "", "%lf %c", &d->abstol, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "check option is malformed: '%s'", arg);
            }
            if (d->abstol < 0) {
                argp_failure(state, EX_USAGE, 0,
                        "abstol value %lf must be positive",
                        d->abstol);
            }
            break;

        case 'T':
            errno = 0;
            {
                double seconds;
                if (1 != sscanf(arg ? arg : "", "%lf %c", &seconds, &ignore)) {
                    argp_failure(state, EX_USAGE, errno,
                            "timelimit option is malformed: '%s'", arg);
                }
                if (seconds < 0) {
                    argp_failure(state, EX_USAGE, 0,
                            "timelimit value %lf must be nonnegative",
                            seconds);
                }
                fftw_set_timelimit(seconds);
            }
            break;

        case 'i':
            d->mpi_inplace = 1;
            break;

        case 'I':
            d->fft_inplace = 1;
            break;

        case 'v':
            ++d->verbose;
            break;

        case 'r':
            errno = 0;
            if (1 != sscanf(arg ? arg : "", "%d %c", &d->repeat, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "repeat option is malformed: '%s'", arg);
            }
            if (d->repeat < 0) {
                argp_failure(state, EX_USAGE, 0,
                        "repeat value %d must be nonnegative",
                        d->repeat);
            }
            break;

        case 't':
            errno = 0;
            if (1 != sscanf(arg ? arg : "", "%d %c", &d->nthreads, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "nthreads is malformed: '%s'", arg);
            }
            if (d->nthreads < 0) {
                argp_failure(state, EX_USAGE, 0,
                        "nthreads value %d must be nonnegative", d->nthreads);
            }
            break;

        case 'n':
            errno = 0;
            if (1 != sscanf(arg ? arg : "", "%d %c", &d->nfields, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "nfields is malformed: '%s'", arg);
            }
            if (d->nfields < 0) {
                argp_failure(state, EX_USAGE, 0,
                        "nfields value %d must be nonnegative", d->nfields);
            }
            break;

        case 'h':
            errno = 0;
            if (1 != sscanf(arg ? arg : "", "%d %c", &d->howmany, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "howmany option is malformed: '%s'", arg);
            }
            if (d->howmany < 1) {
                argp_failure(state, EX_USAGE, 0,
                        "howmany value %d must be strictly positive",
                        d->howmany);
            }
            break;

        case 'f':
            if (d->n0 || d->n1 || d->n2 ) {
                argp_error(state, "only one of --field-{memory,global}"
                           " may be specified");
            }
            if (arg) {
                d->bytes = from_human_readable_byte_count(arg);
                if (d->bytes < 1)
                    argp_failure(state, EX_USAGE, 0,
                            "field-memory option is malformed: '%s'", arg);
            }
            break;

        case 'F':
            if (d->bytes) {
                argp_error(state, "only one of --field-{memory,global}"
                           " may be specified");
            }
            errno = 0;
            if (3 != sscanf(arg ? arg : "", "%d x %d x %d %c",
                            &d->n0, &d->n1, &d->n2, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "field-global option is malformed: '%s'", arg);
            }
            if (d->n0 < 0 || d->n1 < 0 || d->n2 < 0) {
                argp_failure(state, EX_USAGE, 0,
                        "field-global values %dx%dx%d must be nonnegative",
                        d->n0, d->n1, d->n2);
            }
            break;

        case 'P':
            errno = 0;
            if (2 != sscanf(arg ? arg : "", "%d x %d %c",
                            &d->pA, &d->pB, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "dims option not of form pAxpB: '%s'", arg);
            }
            if (d->pA < 0 || d->pB < 0) {
                argp_failure(state, EX_USAGE, 0,
                        "field-global values %dx%d must be nonnegative",
                        d->pA, d->pB);
            }
            break;

        default:
            return ARGP_ERR_UNKNOWN;
    }

    return 0;
}

static struct argp argp = {
    options, parse_opt, args_doc, doc,
    0 /*children*/, 0 /*help_filter*/, 0 /*argp_domain*/
};

// Rank-dependent output streams established in main().
static FILE *rankout, *rankerr;

int main(int argc, char *argv[])
{
    int retval = 0;

    // For miscellaneous timing usage
    double tstart;

    // Initialize default argument storage and default values
    struct details d;
    memset(&d, 0, sizeof(struct details));
    d.repeat   = 1;
    d.nthreads = 0;
    d.nfields  = 1;
    d.howmany  = 2;
    d.forward  = 1;
    d.backward = 1;
    d.fft_n[2] = 'i';
    d.fft_n[1] = 'i';
    d.fft_n[0] = 'i';
    d.abstol   = -1;

    // Initialize/finalize MPI with profiling initially disabled
    MPI_Init(&argc, &argv);
    atexit((void (*) ()) MPI_Finalize);
    MPI_Pcontrol(0);
    MPI_Comm_size(MPI_COMM_WORLD, &d.world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &d.world_rank);

    // Establish rank-dependent output streams
    errno = 0;
    if (d.world_rank == 0) {
        rankout = stdout;
        rankerr = stderr;
    } else {
        rankout = fopen("/dev/null", "w");
        rankerr = rankout;
    }
    if (!rankout) {
        perror("Unable to open rank-dependent output streams"),
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Parse command line arguments using MPI-savvy argp extension
    mpi_argp_parse(d.world_rank, &argp, argc, argv, 0, 0, &d);

    // Initialize/finalize FFTW threads, FFTW MPI, underling
    underling_init(&argc, &argv, d.nthreads);
    atexit(&underling_cleanup);

    // Print program banner that shows version and program arguments
    fprintf(rankout, "%s invoked as\n\t", argp_program_version);
    for (int i = 0; i < argc; ++i) {
        fprintf(rankout, " %s", argv[i]);
    }
    fprintf(rankout, "\n");

    // If available, load wisdom from disk on rank 0 and broadcast it
    // Attempt advisory locking to reduce processes stepping on each other
    if (d.wisdom_file && d.world_rank == 0) {
        FILE *w = fopen(d.wisdom_file, "r");
        if (w) {
            fprintf(rankout, "Loading wisdom from file %s\n", d.wisdom_file);
            if (flock(fileno(w), LOCK_SH)) {
                fprintf(rankout, "WARNING: LOCK_SH failed on wisdom file: %s\n",
                        strerror(errno));
            }
            fftw_import_wisdom_from_file(w);
            if (flock(fileno(w), LOCK_UN)) {
                fprintf(rankout, "WARNING: LOCK_UN failed on wisdom file: %s\n",
                        strerror(errno));
            }
            fclose(w);
        } else {
            fprintf(rankout, "WARNING: Unable to open file %s: %s\n",
                    d.wisdom_file, strerror(errno));
        }
    }
    fftw_mpi_broadcast_wisdom(MPI_COMM_WORLD);

    // If necessary, compute global grid size from per-rank memory constraint
    if (d.bytes) {
        const double nvectors = (d.bytes * d.world_size)
                              / ((double) d.howmany * sizeof(underling_real))
                              / ((double) d.nfields);
        d.n0 = d.n1 = d.n2 = ceil(cbrt(nvectors));

        double coeff;
        const char *units;
        to_human_readable_byte_count(d.bytes, 0, &coeff, &units);
        fprintf(rankout,
            "Per-rank %.2f %s memory requested => %d x %d x %d problem\n",
            coeff, units, d.n0, d.n1, d.n2);
    }

    // Ensure a non-trivial grid was requested
    if (!d.n0 || !d.n1 || !d.n2) {
        fprintf(rankout,
                "You must specify either --field-memory or --field-global!\n");
        MPI_Abort(MPI_COMM_WORLD, EX_USAGE);
    }

    // Initialize underling_grid and print decomposition banner
    fprintf(rankout,
        "\n"
        "Transform and transpose operation to be benchmarked:\n"
        "------------------------------------------------------------------------------\n"
        );
    if (d.forward) {
        fprintf(rankout,
            "Forward:  FFT n2  | n2->n1 | FFT n1  | n1->n0 | FFT n0\n"
            "          %-4s%-3s |   %-3s  | %-4s%-3s |   %-3s  | %-4s%3s\n",
            d.fft_n[2] == 'i' ? "none" :
            d.fft_n[2] == 'r' ? "r2c"  :
            d.fft_n[2] == 'c' ? "c2c"  : "ERR",
            d.fft_n[2] == 'i' ? "" : d.fft_inplace ? "in" : "out",
            d.mpi_inplace ? "in" : "out",
            d.fft_n[1] == 'i' ? "none" :
            d.fft_n[1] == 'r' ? "r2c"  :
            d.fft_n[1] == 'c' ? "c2c"  : "ERR",
            d.fft_n[1] == 'i' ? "" : d.fft_inplace ? "in" : "out",
            d.mpi_inplace ? "in" : "out",
            d.fft_n[0] == 'i' ? "none" :
            d.fft_n[0] == 'r' ? "r2c"  :
            d.fft_n[0] == 'c' ? "c2c"  : "ERR",
            d.fft_n[0] == 'i' ? "" : d.fft_inplace ? "in" : "out");
    }
    if (d.forward && d.backward) fprintf(rankout, "\n");
    if (d.backward) {
        fprintf(rankout,
            "Backward: FFT n0  | n0->n1 | FFT n1  | n1->n2 | FFT n2\n"
            "          %-4s%-3s |   %-3s  | %-4s%-3s |   %-3s  | %-4s%3s\n",
            d.fft_n[0] == 'i' ? "none" :
            d.fft_n[0] == 'r' ? "c2r"  :
            d.fft_n[0] == 'c' ? "c2c"  : "ERR",
            d.fft_n[0] == 'i' ? "" : d.fft_inplace ? "in" : "out",
            d.mpi_inplace ? "in" : "out",
            d.fft_n[1] == 'i' ? "none" :
            d.fft_n[1] == 'r' ? "c2r"  :
            d.fft_n[1] == 'c' ? "c2c"  : "ERR",
            d.fft_n[1] == 'i' ? "" : d.fft_inplace ? "in" : "out",
            d.mpi_inplace ? "in" : "out",
            d.fft_n[2] == 'i' ? "none" :
            d.fft_n[2] == 'r' ? "c2r"  :
            d.fft_n[2] == 'c' ? "c2c"  : "ERR",
            d.fft_n[2] == 'i' ? "" : d.fft_inplace ? "in" : "out");
    }
    if (d.fft_n[0] == 'r' || d.fft_n[1] == 'r' || d.fft_n[2] == 'r') {
        fprintf(rankout,
            "\n"
            "Real-to-complex FFTs imply real-valued grid is %d x %d x %d\n",
            d.fft_n[0] == 'r' ? (2*(d.n0 - 1)+(d.n0 == 1)) : d.n0,
            d.fft_n[1] == 'r' ? (2*(d.n1 - 1)+(d.n1 == 1)) : d.n1,
            d.fft_n[2] == 'r' ? (2*(d.n2 - 1)+(d.n2 == 1)) : d.n2);
    }
    fprintf(rankout,
        "------------------------------------------------------------------------------\n"
        "\n");

    // Initialize underling_grid and print decomposition banner
    underling_grid grid = underling_grid_create(
            MPI_COMM_WORLD, d.n0, d.n1, d.n2, d.pA, d.pB);
    assert(grid);
    d.pA = underling_grid_pA_size(grid); // Obtain any automagic pA value
    d.pB = underling_grid_pB_size(grid); // Obtain any automagic pB value
    fprintf(rankout,
        "Global pencil decomposition details (for transposed_flags == 0)\n"
        "------------------------------------------------------------------------------\n"
        "Long n2:                                     (%1$5d/%5$4d x %2$5d/%4$4d) x %3$5d\n"
        "Long n1: %3$5d/%4$4d x (%1$5d/%5$4d x %2$5d) = (%3$5d/%4$4d x %1$5d/%5$4d) x %2$5d\n"
        "Long n0: %2$5d/%5$4d x (%3$5d/%4$4d x %1$5d) = (%2$5d/%5$4d x %3$5d/%4$4d) x %1$5d\n"
        "------------------------------------------------------------------------------\n"
        "\n", d.n0, d.n1, d.n2, d.pA, d.pB);

    // Initialize underling_problem and find per-field memory requirements
    underling_problem problem = underling_problem_create(
            grid, d.howmany, d.transposed_flags);
    assert(grid);
    const size_t local_memory = underling_local_memory(problem);

    // Display some information about the problem's memory requirements
    {
        double coeff1, coeff2, coeff3;
        const char *units1, *units2, *units3;

        to_human_readable_byte_count(underling_local_memory_optimum(problem)
                * sizeof(underling_real), 0, &coeff1, &units1);
        to_human_readable_byte_count(underling_local_memory_minimum(grid, problem)
                * sizeof(underling_real), 0, &coeff2, &units2);
        to_human_readable_byte_count(underling_local_memory_maximum(grid, problem)
                * sizeof(underling_real), 0, &coeff3, &units3);
        fprintf(rankout, "Optimum per-rank, per-field memory is %.4f %s vs actual %.4f %s -- %.4f %s\n",
                coeff1, units1, coeff2, units2, coeff3, units3);

        to_human_readable_byte_count(d.nfields
                    * underling_local_memory_optimum(problem)
                    * sizeof(underling_real),
                    0, &coeff1, &units1);
        to_human_readable_byte_count(d.nfields
                    * underling_local_memory_maximum(grid, problem)
                    * sizeof(underling_real),
                    0, &coeff2, &units2);
        fprintf(rankout, "Optimum per-rank, total memory is %.4f %s vs worst case %.4f %s\n",
                coeff1, units1, coeff2, units2);

        to_human_readable_byte_count(d.nfields
                    * underling_global_memory_optimum(grid, problem)
                    * sizeof(underling_real),
                    0, &coeff1, &units1);
        to_human_readable_byte_count(d.nfields
                    * underling_global_memory(grid, problem)
                    * sizeof(underling_real),
                    0, &coeff2, &units2);
        fprintf(rankout, "Optimum global, total memory is %.4f %s vs actual %.4f %s\n",
                coeff1, units1, coeff2, units2);
    }

    // Adjust for in- versus out-of-place operation for MPI and FFT substeps.
    // Forward:  FFT n2, n2->n1, FFT n1, n1->n0, FFT n0
    // m[i] is -1 for out-of-place to a lower-indexed buffer, 0 for in-place,
    // and 1 for out-of-place to a higher-indexed buffer.
    int m[10];
    {
        // s tracks the aggregate out-of-place movement that we want to stay in
        // {0, 1}.  To do so, out-of-place operations (!inplace) are negated
        // via a ConditionalNegate bithack before being summed.
        int s;
        s  = (m[0] =  !(d.fft_n[2] == 'i' || d.fft_inplace));
        s += (m[1] = (!d.mpi_inplace ^ -s) + s);
        s += (m[2] = (!(d.fft_n[1] == 'i' || d.fft_inplace) ^ -s) + s);
        s += (m[3] = (!d.mpi_inplace ^ -s) + s);
        s += (m[4] = (!(d.fft_n[0] == 'i' || d.fft_inplace) ^ -s) + s);

        // Backward: FFT n0, n0->n1, FFT n1, n1->n2, FFT n2
        // Backward operations are the inverse of forward operations
        m[5] = -m[4];
        m[6] = -m[3];
        m[7] = -m[2];
        m[8] = -m[1];
        m[9] = -m[0];
    }

    // Do we need an additional field worth of storage for out-of-place ops?
    const int extra = m[0] || m[1] || m[2] || m[3] || m[4] || m[5];

    // Allocate memory for each field plus one optional scratch buffer
    underling_real *f[d.nfields + extra]; // C99
    for (int i = 0; i < d.nfields + extra; ++i) {
        f[i] = calloc_field(&d, local_memory * sizeof(underling_real));
        assert(f[i]);
    }

    // Create the transpose plan
    fprintf(rankout, "\nInvoking underling_plan_create...\n");
    tstart = MPI_Wtime();
    underling_plan t_plan = underling_plan_create(
            problem, f[0], f[!d.mpi_inplace],
            d.transform_flags
                | (d.forward   ? UNDERLING_TRANSPOSE_LONG_N2_TO_LONG_N1
                               | UNDERLING_TRANSPOSE_LONG_N1_TO_LONG_N0 : 0)
                | (d.backward  ? UNDERLING_TRANSPOSE_LONG_N0_TO_LONG_N1
                               | UNDERLING_TRANSPOSE_LONG_N1_TO_LONG_N2 : 0),
            d.fftw_rigor_flags);
    fprintf(rankout, "...underling_plan_create took %lf seconds",
            MPI_Wtime() - tstart);
    if (d.verbose) {
        fputs(" and returned (on rank 0):\n", rankout);
        underling_fprint_plan(t_plan, rankout);
    }
    fputc('\n', rankout);

    // Create any necessary transform plans.
    // Note forward plans are preferred to create corresponding backward plans
    // but we use only FFTW_ESTIMATE if a forward plan is superfluous.
    underling_fftw_plan forward_plan[3]  = { NULL, NULL, NULL };
    underling_fftw_plan backward_plan[3] = { NULL, NULL, NULL };
    for (int i = 3; i --> 0 ;) {
        underling_fftw_plan (*planner)(const underling_problem,
                                       int, underling_real *, underling_real *,
                                       unsigned, unsigned) = NULL;
        const char * planner_name = NULL;
        switch (d.fft_n[i]) {
            case 'c':
                planner = &underling_fftw_plan_create_c2c_forward;
                planner_name = "underling_fftw_plan_create_c2c_forward";
                break;
            case 'r':
                planner = &underling_fftw_plan_create_r2c_forward;
                planner_name = "underling_fftw_plan_create_r2c_forward";
                break;
        }
        if (planner) {
            fprintf(rankout, "\nInvoking %s for long_n%d...\n", planner_name, i);
            tstart = MPI_Wtime();
            forward_plan[i] = planner(problem, i, f[0], f[!d.fft_inplace],
                    d.forward ? d.fftw_rigor_flags : FFTW_ESTIMATE,
                    d.packed_flags);
            fprintf(rankout, "...%s took %lf seconds",
                    planner_name, MPI_Wtime() - tstart);
            if (d.verbose) {
                fputs(" and returned (on rank 0):\n", rankout);
                underling_fftw_fprint_plan(forward_plan[i], rankout);
            }
            fputc('\n', rankout);
        }
        if (d.backward && forward_plan[i]) {
            fprintf(rankout, "\nInvoking underling_fftw_plan_create_inverse for long_n%d...\n", i);
            tstart = MPI_Wtime();
            backward_plan[i] = underling_fftw_plan_create_inverse(
                    forward_plan[i], f[0], f[!d.fft_inplace],
                    d.backward ? d.fftw_rigor_flags : FFTW_ESTIMATE);
            fprintf(rankout, "...underling_fftw_plan_create_inverse took %lf seconds",
                    MPI_Wtime() - tstart);
            if (d.verbose) {
                fputs(" and returned (on rank 0):\n", rankout);
                underling_fftw_fprint_plan(forward_plan[i], rankout);
            }
            fputc('\n', rankout);
        }
    }

    // If requested, gather wisdom and then write to disk on rank 0
    // Attempt advisory locking to reduce processes stepping on each other
    if (d.wisdom_file) {
        fftw_mpi_gather_wisdom(MPI_COMM_WORLD);
    }
    if (d.wisdom_file && d.world_rank == 0) {
        FILE *w = fopen(d.wisdom_file, "w+");
        if (w) {
            fprintf(rankout, "Saving wisdom to file %s\n", d.wisdom_file);
            if (flock(fileno(w), LOCK_EX)) {
                fprintf(rankout, "WARNING: LOCK_EX failed on wisdom file: %s\n",
                        strerror(errno));
            }
            fftw_export_wisdom_to_file(w);
            if (flock(fileno(w), LOCK_UN)) {
                fprintf(rankout, "WARNING: LOCK_UN failed on wisdom file: %s\n",
                        strerror(errno));
            }
            fclose(w);
        } else {
            fprintf(rankout, "WARNING: Unable to open file %s: %s\n",
                    d.wisdom_file, strerror(errno));
        }
    }

    // If sufficiently verbose, output extent information for all stages.
    // Perform one rank at a time to avoid jumbling output.
    if ((d.world_size==1 && d.verbose) || (d.world_size>1 && d.verbose>1)) {
        const int stage[6] = {2, 1, 0, 0, 1, 2};         // Data traversal
        const underling_fftw_plan *dp[6] = {             // across forward and
            forward_plan,  forward_plan,  forward_plan,  // backward directions
            backward_plan, backward_plan, backward_plan  // in a single loop
        };

        for (int r = 0; r < d.world_size; ++r) {         // Walk all ranks
            MPI_Barrier(MPI_COMM_WORLD);                 // Synchronize
            if (r != d.world_rank) continue;             // Take turns in loop

            printf("\nStep-by-step storage info for rank %d:\n", d.world_rank);
            for (size_t i = 0; i < sizeof(stage)/sizeof(stage[0]); ++i) {
                const int l = stage[i];
                if ((dp[i])[l]) {
                    underling_fftw_extents in
                        = underling_fftw_local_extents_input((dp[i])[l]);
                    printf("\tlong_n%d FFT input:  ", l);
                    underling_fftw_fprint_extents(&in, stdout);
                    putchar('\n');

                    underling_fftw_extents out
                        = underling_fftw_local_extents_output((dp[i])[l]);
                    printf("\tlong_n%d FFT output: ", l);
                    underling_fftw_fprint_extents(&out, stdout);
                    putchar('\n');
                } else {
                    underling_extents t = underling_local_extents(problem, l);
                    printf("\tlong_n%d storage:    ", l);
                    underling_fprint_extents(&t, stdout);
                    putchar('\n');
                }
            }
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);  // Let last rank flush
    }

    // Fill state fields with well-defined garbage for check_field purposes
    // We fill all of the field, not just the check_field-relevant parts
    if (d.forward && d.backward) {
        for (int i = 0; i < d.nfields; ++i) {
            fill_field(f[i], &d, local_memory * sizeof(underling_real), i);
        }
    }

    // Precompute the normalization factor required for check_field purposes
    // Normalization factor is the product of all logical Fourier extents
    const underling_real normalization = ((underling_real) 1) / (
            ((long) (
                d.fft_n[0] == 'c' ? d.n0 :
                d.fft_n[0] == 'r' ? (2*(d.n0 - 1)+(d.n0 == 1)) : 1
            )) * ((long) (
                d.fft_n[1] == 'c' ? d.n1 :
                d.fft_n[1] == 'r' ? (2*(d.n1 - 1)+(d.n1 == 1)) : 1
            )) * ((long) (
                d.fft_n[2] == 'c' ? d.n2 :
                d.fft_n[2] == 'r' ? (2*(d.n2 - 1)+(d.n2 == 1)) : 1
            ))
        );

    // During execution field N is pointed to by f[o[N]] for {0, ..., nfields}.
    // Important as m < 0 and m > 0 must run different directions through memory.
    int o[d.nfields]; // C99

    fprintf(rankout, "\nBeginning benchmark main loop...\n");
    GRVY_TIMER_INIT(argp_program_version);
    double tmean = 0;  // Mean, sample variance per Knuth/Welford algorithm from
    double tM2   = 0;  // wikipedia.org/wiki/Algorithms_for_calculating_variance
    for (int i = 0; i < d.repeat; ++i) {
        fprintf(rankout, "\tIteration %d\n", i);
        const double tstart = MPI_Wtime();

        if (d.forward) {

            // Ensure we start from long_n2 (required for forward-only benchmark)
            for (int i = 0; i < d.nfields; ++i)
                o[i] = i;

            // FFT forward transform long_n2
            if (forward_plan[2]) {
                for (int j = 0; j < d.nfields; ++j) {
                    const int k = (m[0] >= 0 ? j : d.nfields - 1 - j);
                    GRVY_TIMER_BEGIN("underling_fftw long_n2 forward");
                    underling_fftw_plan_execute(
                            forward_plan[2], f[o[k]], f[o[k] + m[0]]);
                    GRVY_TIMER_END("underling_fftw long_n2 forward");
                }
                for (int j = 0; j < d.nfields; ++j) o[j] += m[0];
            }

            // MPI transpose long_n2 -> long_n1
            for (int j = 0; j < d.nfields; ++j) {
                const int k = (m[1] >= 0 ? j : d.nfields - 1 - j);
                GRVY_TIMER_BEGIN("underling long_n2_to_long_n1");
                underling_execute_long_n2_to_long_n1(
                        t_plan, f[o[k]], f[o[k] + m[1]]);
                GRVY_TIMER_END("underling long_n2_to_long_n1");
            }
            for (int j = 0; j < d.nfields; ++j) o[j] += m[1];

            // FFT forward transform long_n1
            if (forward_plan[1]) {
                for (int j = 0; j < d.nfields; ++j) {
                    const int k = (m[2] >= 0 ? j : d.nfields - 1 - j);
                    GRVY_TIMER_BEGIN("underling_fftw long_n1 forward");
                    underling_fftw_plan_execute(
                            forward_plan[1], f[o[k]], f[o[k] + m[2]]);
                    GRVY_TIMER_END("underling_fftw long_n1 forward");
                }
            }
            for (int j = 0; j < d.nfields; ++j) o[j] += m[2];

            // MPI transpose long_n1 -> long_n0
            for (int j = 0; j < d.nfields; ++j) {
                const int k = (m[3] >= 0 ? j : d.nfields - 1 - j);
                GRVY_TIMER_BEGIN("underling long_n1_to_long_n0");
                underling_execute_long_n1_to_long_n0(
                        t_plan, f[o[k]], f[o[k] + m[3]]);
                GRVY_TIMER_END("underling long_n1_to_long_n0");
            }
            for (int j = 0; j < d.nfields; ++j) o[j] += m[3];

            // FFT forward transform long_n0
            if (forward_plan[0]) {
                for (int j = 0; j < d.nfields; ++j) {
                    const int k = (m[4] >= 0 ? j : d.nfields - 1 - j);
                    GRVY_TIMER_BEGIN("underling_fftw long_n0 forward");
                    underling_fftw_plan_execute(
                            forward_plan[0], f[o[k]], f[o[k] + m[4]]);
                    GRVY_TIMER_END("underling_fftw long_n0 forward");
                }
                for (int j = 0; j < d.nfields; ++j) o[j] += m[4];
            }

        }

        if (d.backward) {

            // Ensure we start from long_n0 (required for backward-only benchmark)
            for (int i = 0; i < d.nfields; ++i)
                o[i] = m[0] + m[1] + m[2] + m[3] + m[4];

            // FFT backward transform long_n0
            if (backward_plan[0]) {
                for (int j = 0; j < d.nfields; ++j) {
                    const int k = (m[5] >= 0 ? j : d.nfields - 1 - j);
                    GRVY_TIMER_BEGIN("underling_fftw long_n0 backward");
                    underling_fftw_plan_execute(
                            backward_plan[0], f[o[k]], f[o[k] + m[5]]);
                    GRVY_TIMER_END("underling_fftw long_n0 backward");
                }
                for (int j = 0; j < d.nfields; ++j) o[j] += m[5];
            }

            // MPI transpose long_n0 -> long_n1
            for (int j = 0; j < d.nfields; ++j) {
                const int k = (m[6] >= 0 ? j : d.nfields - 1 - j);
                GRVY_TIMER_BEGIN("underling long_n0_to_long_n1");
                underling_execute_long_n0_to_long_n1(
                        t_plan, f[o[k]], f[o[k] + m[6]]);
                GRVY_TIMER_END("underling long_n0_to_long_n1");
            }
            for (int j = 0; j < d.nfields; ++j) o[j] += m[6];

            // FFT backward transform long_n1
            if (backward_plan[1]) {
                for (int j = 0; j < d.nfields; ++j) {
                    const int k = (m[7] >= 0 ? j : d.nfields - 1 - j);
                    GRVY_TIMER_BEGIN("underling_fftw long_n1 backward");
                    underling_fftw_plan_execute(
                            backward_plan[1], f[o[k]], f[o[k] + m[7]]);
                    GRVY_TIMER_END("underling_fftw long_n1 backward");
                }
                for (int j = 0; j < d.nfields; ++j) o[j] += m[7];
            }

            // MPI transpose long_n1 -> long_n2
            for (int j = 0; j < d.nfields; ++j) {
                const int k = (m[8] >= 0 ? j : d.nfields - 1 - j);
                GRVY_TIMER_BEGIN("underling long_n1_to_long_n2");
                underling_execute_long_n1_to_long_n2(
                        t_plan, f[o[k]], f[o[k] + m[8]]);
                GRVY_TIMER_END("underling long_n1_to_long_n2");
            }
            for (int j = 0; j < d.nfields; ++j) o[j] += m[8];

            // FFT backward transform long_n2
            if (backward_plan[2]) {
                for (int j = 0; j < d.nfields; ++j) {
                    const int k = (m[9] >= 0 ? j : d.nfields - 1 - j);
                    GRVY_TIMER_BEGIN("underling_fftw long_n2 backward");
                    underling_fftw_plan_execute(
                            backward_plan[2], f[o[k]], f[o[k] + m[9]]);
                    GRVY_TIMER_END("underling_fftw long_n2 backward");
                }
                for (int j = 0; j < d.nfields; ++j) o[j] += m[9];
            }

        }

        // Update mean and M2 for timing per online Knuth/Welford algorithm
        {
            const double x     = MPI_Wtime() - tstart;
            const double delta = x - tmean;
            tmean += delta/(i + 1);
            tM2   += delta*(x - tmean);
        }

        // Apply forward-and-backward normalization for check_field purposes
        // Ignore for timing as it is not strictly part of the benchmark
        if (d.forward && d.backward) {
            for (int i = 0; i < d.nfields; ++i) {
                for (size_t j = 0; j < local_memory; ++j) {
                    f[i][j] *= normalization;
                }
            }
        }
    }
    GRVY_TIMER_FINALIZE();
    fprintf(rankout, "...ended benchmark main loop\n");
    double tvariance = tM2 / (d.repeat - 1);

    // TODO Get timing information back from multiple ranks
    if (d.world_rank == 0) { GRVY_TIMER_SUMMARIZE(); }

    // Acquire and display worst-observed tmean and tvariance from any rank
    if (d.repeat) {
        struct { double val; int rank; } recvbuf[2], sendbuf[2] = {
            { tmean, d.world_rank }, { tvariance, d.world_rank }
        };
        MPI_Reduce(&sendbuf, &recvbuf, 2, MPI_DOUBLE_INT,
                   MPI_MAXLOC, 0, MPI_COMM_WORLD);

        if (d.world_rank == 0 && d.repeat == 1) {
            printf("Iteration time was %8.6g seconds", recvbuf[0].val);
            if (d.world_size > 1) printf(" on worst rank %d", recvbuf[0].rank);
            putchar('\n');
        } else if (d.world_rank == 0) {
            printf("Iteration mean time was %8.6g seconds across %d iterations",
                   recvbuf[0].val, d.repeat);
            if (d.world_size > 1) printf(" on worst rank %d", recvbuf[0].rank);
            putchar('\n');
            printf("Iteration variance was %8.6g seconds", recvbuf[1].val);
            if (d.world_size > 1) printf(" on worst rank %d", recvbuf[1].rank);
            putchar('\n');
        } else {
            // NOP
        }
    }

    // If we round tripped...
    if (d.repeat && d.forward && d.backward) {

        // ...compute the maximum absolute error observed on any rank...
        // ...(using the same salt as fill_field for consistency)...
        // ...(checking only data and not communication buffers)...
        underling_real maxabserr = 0;
        for (int i = 0; i < d.nfields; ++i) {
            underling_extents e = underling_local_extents(problem, 2);
            assert(e.extent <= local_memory);
            const underling_real abserr = check_field(
                    f[i], &d, e.extent * sizeof(underling_real),
                    i, forward_plan[2]);
            maxabserr = fmax(maxabserr, abserr);
        }

        // ...bring it down to rank zero at known precision and display...
        struct { double val; int rank; } commbuf = { maxabserr, d.world_rank };
        MPI_Allreduce(MPI_IN_PLACE, &commbuf, 1,
                      MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
        fprintf(rankout,
                "Maximum absolute accumulated error %g observed on rank %d\n",
                commbuf.val, commbuf.rank);

        // ...and possibly warn of abstol failure or genuinely fail...
        // ...with failure being uniform across all ranks due to Allreduce.
        if (d.abstol >= 0 && commbuf.val > d.abstol / normalization) {
            fprintf(rankout,
                "Maximum observed error exceeds scaled abstol of %g\n",
                d.abstol / normalization);
            retval = 1;
        }
    }

    // Deallocate memory for each field plus one possible scratch buffer
    for (int i = 0; i < d.nfields + extra; ++i) {
        fftw_free(f[i]);
        f[i] = NULL;
    }

    // Tear down underling_plan, underling_problem, underling_grid
    for (int i = 0; i < 3; ++i) {
        if (backward_plan[i]) underling_fftw_plan_destroy(backward_plan[i]);
    }
    for (int i = 0; i < 3; ++i) {
        if (forward_plan[i]) underling_fftw_plan_destroy(forward_plan[i]);
    }
    underling_plan_destroy(t_plan);
    underling_problem_destroy(problem);
    underling_grid_destroy(grid);

    // Finalizing of MPI, FFTW, underling handled by atexit()

    return retval;
}

void print_version(FILE *stream, struct argp_state *state)
{
    (void) state; // Unused

    fputs(argp_program_version, stream);
    fprintf(stream, " linked against FFTW3 %s", fftw_version);
    int version, subversion;
    if (MPI_SUCCESS == MPI_Get_version(&version, &subversion)) {
        fprintf(stream, " running atop MPI %d.%d", version, subversion);
    }
    fputc('\n', stream);
}


void trim(char *a)
{
    char *b = a;
    while (isspace(*b))   ++b;
    while (*b)            *a++ = *b++;
    *a = '\0';
    while (isspace(*--a)) *a = '\0';
}


static inline int min(int a, int b)
{
    return a < b ? a : b;
}


static inline int max(int a, int b)
{
    return a > b ? a : b;
}


// Adapted from http://stackoverflow.com/questions/3758606/
// how-to-convert-byte-size-into-human-readable-format-in-java
static void to_human_readable_byte_count(long bytes,
                                         int si,
                                         double *coeff,
                                         const char **units)
{
    // Static lookup table of byte-based SI units
    static const char *suffix[][2] = { { "B",  "B"   },
                                       { "kB", "KiB" },
                                       { "MB", "MiB" },
                                       { "GB", "GiB" },
                                       { "TB", "TiB" },
                                       { "EB", "EiB" },
                                       { "ZB", "ZiB" },
                                       { "YB", "YiB" } };
    int unit = si ? 1000 : 1024;
    int exp = 0;
    if (bytes > 0) {
        exp = min( (int) (log(bytes) / log(unit)),
                   (int) sizeof(suffix) / sizeof(suffix[0]) - 1);
    }
    *coeff = bytes / pow(unit, exp);
    *units  = suffix[exp][!!si];
}


// Convert strings like the following into byte counts
//    5MB, 5 MB, 5M, 3.7GB, 123b, 456kBytes
// with some amount of forgiveness baked into the parsing.
static long from_human_readable_byte_count(const char *str)
{
    // Parse leading numeric factor
    char *endptr;
    errno = 0;
    const double coeff = strtod(str, &endptr);
    if (errno) return -1;

    // Skip any intermediate white space
    while (isspace(*endptr)) ++endptr;

    // Read off first character which should be an SI prefix
    int exp  = 0;
    int unit = 1024;
    switch (toupper(*endptr)) {
        case 'B':  exp =  0; break;
        case 'K':  exp =  3; break;
        case 'M':  exp =  6; break;
        case 'G':  exp =  9; break;
        case 'T':  exp = 12; break;
        case 'E':  exp = 15; break;
        case 'Z':  exp = 18; break;
        case 'Y':  exp = 21; break;

        case ' ':
        case '\t':
        case '\0': exp =  0; goto done;

        default:   return -1;
    }
    ++endptr;

    // If an 'i' or 'I' is present use SI factor-of-1000 units
    if (toupper(*endptr) == 'I') {
        ++endptr;
        unit = 1000;
    }

    // Next character must be one of B/empty/whitespace
    switch (toupper(*endptr)) {
        case 'B':
        case ' ':
        case '\t': ++endptr;  break;

        case '\0': goto done;

        default:   return -1;
    }

    // Skip any remaining white space
    while (isspace(*endptr)) ++endptr;

    // Parse error on anything but a null terminator
    if (*endptr) return -1;

done:
    return exp ? coeff * pow(unit, exp / 3) : coeff;
}

static underling_real* calloc_field(struct details *d,
                                    const long bytes)
{
    underling_real *p = fftw_malloc(bytes); // Aligned malloc
    if (!p) {
        double coeff;
        const char *units;
        to_human_readable_byte_count(bytes, 0, &coeff, &units);
        fprintf(stderr, "Unable to malloc %.2f %s on rank %d\n",
                coeff, units, d->world_rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    memset(p, 0, bytes);                    // Make it a calloc

    return p;
}

static void fill_field(underling_real *p,
                       struct details *d,
                       const long bytes,
                       unsigned salt)
{
    assert(p);

    // Store previous random seed and provide a new one based on rank, salt
    char state[64];
    initstate(d->world_rank + salt, state, sizeof(state)/sizeof(state[0]));
    char *previous = setstate(state);

    // Fill using well-defined garbage
    const size_t count = bytes / sizeof(underling_real);
    const underling_real inv_rand_max = ((underling_real) 1) / RAND_MAX;
    for (size_t i = 0; i < count; ++i)
        p[i] = random()*inv_rand_max + salt;

    // Restore previous random state
    setstate(previous);
}

// Check that the data in p matches the values set by fill_field.  One annoying
// complication is that real-value fields may contain padded entries not
// preserved across transforms.  We use data from plan to account for that.
// Another annoyance is that fields may contain communication-only scratch
// space-- we rely on the caller to not specify parameter bytes beyond the
// details we care about.  In practice, this uses underling_extents.extents.
static underling_real check_field(const underling_real *p,
                                  struct details *d,
                                  const long bytes,
                                  unsigned salt,
                                  underling_fftw_plan plan)
{
    assert(p);

    // Store previous random seed and provide a new one based on rank, salt
    char state[64];
    initstate(d->world_rank + salt, state, sizeof(state)/sizeof(state[0]));
    char *previous = setstate(state);

    // Compute maximum observed error using random numbers from fill_field
    underling_real maxabserr = 0;
    size_t count = bytes / sizeof(underling_real);
    const underling_real inv_rand_max = ((underling_real) 1) / RAND_MAX;

    // Check if we need to specially handle for real-valued padding...
    if (plan) {
        underling_fftw_extents e = underling_fftw_local_extents_input(plan);
        if (e.size[4] == 1 /* i.e., real-valued */) {
            // ...Yes we do.
            // Ignore padding between and after any relevant real scalars
            count = MIN(count,(size_t)e.stride[e.order[4]]*e.size[e.order[4]]);
            for (size_t i = 0; i < count; ++i) {
                const underling_real expected = random() * inv_rand_max + salt;
                if (i % e.stride[e.order[3]] < (size_t) e.size[e.order[2]]) {
                    maxabserr = fmax(maxabserr, fabs(p[i] - expected));
                }
            }

            setstate(previous);  // Restore previous random state
            return maxabserr;    // Return to caller
        }
    }

    // ...No we do not.  Just rip through contiguous memory
    for (size_t i = 0; i < count; ++i) {
        const underling_real expected = random() * inv_rand_max + salt;
        maxabserr = fmax(maxabserr, fabs(p[i] - expected));
    }

    setstate(previous);  // Restore previous random state
    return maxabserr;    // Return to caller
}
