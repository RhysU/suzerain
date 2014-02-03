//--------------------------------------------------------------------------
//
// Copyright (C) 2008-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Suzerain is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------

/** @file
 * @copydoc support.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/support/support.hpp>

#include <esio/esio.h>
#include <esio/error.h>
#include <gsl/gsl_errno.h>
#include <sys/file.h>

#include <suzerain/common.hpp>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/coalescing_pool.hpp>
#include <suzerain/countof.h>
#include <suzerain/diffwave.hpp>
#include <suzerain/error.h>
#include <suzerain/exprparse.hpp>
#include <suzerain/htstretch.h>
#include <suzerain/math.hpp>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/rholut.hpp>
#include <suzerain/rngstream.hpp>
#include <suzerain/shared_range.hpp>
#include <suzerain/state.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/validation.hpp>

using boost::numeric_cast;
using std::size_t;

// Helps to identify from whom logging messages are being emitted.
static const std::string who("support");

namespace suzerain {

namespace support {

// Common configuration snippet used in multiple places just below.
// See "Configuration" at http://logging.apache.org/log4cxx/index.html
#define COMMON_CONSOLE_CONFIG                                               \
    "log4j.appender.CONSOLE=org.apache.log4j.ConsoleAppender\n"             \
    "log4j.appender.CONSOLE.layout=org.apache.log4j.PatternLayout\n"        \
    "log4j.appender.CONSOLE.layout.ConversionPattern=%-5p %8r %-11c %m%n\n"

const char log4cxx_config_console[] =
    "## Output INFO or higher messages on CONSOLE\n"
    "log4j.rootLogger=INFO, CONSOLE\n"
    COMMON_CONSOLE_CONFIG;

const char log4cxx_config[] =
    "## Output INFO or higher messages on CONSOLE and in LOG\n"
    "log4j.rootLogger=INFO, CONSOLE, LOG\n"
    COMMON_CONSOLE_CONFIG
    "log4j.appender.LOG=org.apache.log4j.FileAppender\n"
    "log4j.appender.LOG.append=true\n"
    "log4j.appender.LOG.filename=log.dat\n"
    "log4j.appender.LOG.layout=${log4j.appender.CONSOLE.layout}\n"
    "log4j.appender.LOG.layout.ConversionPattern=${log4j.appender.CONSOLE.layout.ConversionPattern}\n"
;

#undef COMMON_CONSOLE_CONFIG

const real_t bsplines_distinct_distance
    = 3*std::numeric_limits<real_t>::epsilon();

void mpi_abort_on_error_handler_gsl(const char * reason,
                                    const char * file,
                                    int line,
                                    int error_code)
{
    return mpi_abort_on_error_handler(reason, file, line,
            error_code, "GSL", gsl_strerror(error_code));
}

void mpi_abort_on_error_handler_suzerain(const char * reason,
                                         const char * file,
                                         int line,
                                         int error_code)
{
    return mpi_abort_on_error_handler(reason, file, line,
            error_code, "Suzerain", suzerain_strerror(error_code));
}

void mpi_abort_on_error_handler_esio(const char * reason,
                                     const char * file,
                                     int line,
                                     int error_code)
{
    return mpi_abort_on_error_handler(reason, file, line,
            error_code, "ESIO", esio_strerror(error_code));
}

#ifdef HAVE_UNDERLING
void mpi_abort_on_error_handler_underling(const char * reason,
                                          const char * file,
                                          int line,
                                          int error_code)
{
    return mpi_abort_on_error_handler(reason, file, line,
            error_code, "underling", underling_strerror(error_code));
}
#endif

void mpi_abort_on_error_handler(const char * reason,
                                const char * file,
                                int line,
                                int error_code,
                                const char * origin,
                                const char * strerror)
{
    FATAL((origin ? origin : "NULLORIGIN")
          << " reports '"
          << (reason ? reason : "NULLREASON")
          << "' as code #"
          << error_code
          << " ('"
          << (strerror ? strerror : "NULLSTRERROR")
          << "') from "
          << (file ? file : "NULLFILE")
          << ':'
          << line);
    MPI_Abort(MPI_COMM_WORLD, errno ? errno : EXIT_FAILURE);
}

bool wisdom_broadcast(const std::string& wisdom_file)
{
    int success = 0;
    if (wisdom_file.empty()) return success; // Short circuit if no path provided

    // If available, load wisdom from disk on rank 0 and broadcast it
    // Attempt advisory locking to avoid multiple jobs stepping on each other
    if (mpi::comm_rank(MPI_COMM_WORLD) == 0) {

        // Import any system-wide wisdom available
        fftw_import_system_wisdom();

        FILE *w = fopen(wisdom_file.c_str(), "r");
        if (w) {
            INFO0(who, "Loading FFTW wisdom from file " << wisdom_file);
            if (flock(fileno(w), LOCK_SH)) {
                WARN0(who, "LOCK_SH failed on wisdom file "
                      << wisdom_file << ": " << strerror(errno));
            }
            success = fftw_import_wisdom_from_file(w);
            if (flock(fileno(w), LOCK_UN)) {
                WARN0(who, "LOCK_UN failed on wisdom file "
                      << wisdom_file << ": " << strerror(errno));
            }
            fclose(w);
        } else {
            WARN0(who, "Unable to open wisdom file "
                  << wisdom_file << ": " << strerror(errno));
        }
    }
    fftw_mpi_broadcast_wisdom(MPI_COMM_WORLD);

    MPI_Bcast(&success, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return success;
}

bool wisdom_gather(const std::string& wisdom_file)
{
    int success = 0;
    if (wisdom_file.empty()) return success; // Short circuit if no path provided

    // If available, gather wisdom and then write to disk on rank 0
    // Attempt advisory locking to reduce processes stepping on each other
    fftw_mpi_gather_wisdom(MPI_COMM_WORLD);
    if (mpi::comm_rank(MPI_COMM_WORLD) == 0) {
        FILE *w = fopen(wisdom_file.c_str(), "w+");
        if (w) {
            INFO0(who, "Saving FFTW wisdom to file " << wisdom_file);
            if (flock(fileno(w), LOCK_EX)) {
                WARN0(who, "LOCK_EX failed on wisdom file "
                      << wisdom_file << ": " << strerror(errno));
            }
            fftw_export_wisdom_to_file(w);
            if (flock(fileno(w), LOCK_UN)) {
                WARN0(who, "LOCK_UN failed on wisdom file "
                      << wisdom_file << ": " << strerror(errno));
            }
            success = !fclose(w);
        } else {
            WARN0(who, "Unable to open wisdom file "
                  << wisdom_file << ": " << strerror(errno));
        }
    }

    MPI_Bcast(&success, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return success;
}

real_t create(const int ndof,
              const int k,
              const double left,
              const double right,
              const double htdelta,
              shared_ptr<bspline>& b,
              shared_ptr<bsplineop>& cop)
{
    INFO0(who, "Creating B-spline basis of order " << k
          << " on [" << left << ", " << right << "] with "
          << ndof << " DOF stretched per htdelta " << htdelta);

////FIXME: Knot vectors are non-increasing for moderate htdelta
/// FIXME: See https://savannah.gnu.org/bugs/index.php?34361
////// Compute collocation point locations using ndof and htdelta
////ArrayXr abscissae(ndof);
////math::linspace(0.0, 1.0, abscissae.size(), abscissae.data());
////for (int i = 0; i < abscissae.size(); ++i) {
////    abscissae[i] = suzerain_htstretch2(htdelta, 1.0, abscissae[i]);
////}
////abscissae = (right - left) * abscissae + left;
////
////// Generate the B-spline workspace based on order and abscissae
////// Maximum non-trivial derivative operators included
////double abserr;
////b = make_shared<bspline>(k, bspline::from_abscissae(),
////                         abscissae.size(), abscissae.data(), &abserr);
////assert(b->n() == ndof);
////cop.reset(new bsplineop(
////            *b, k-2, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE));
////assert(cop->n() == ndof);
////
////INFO0("Created B-spline basis has Greville abscissae abserr of " << abserr);
////
////return abserr;

    // Compute breakpoint point locations using ndof and htdelta
    ArrayXr breakpoints(ndof - k + 2);
    math::linspace(0.0, 1.0, breakpoints.size(), breakpoints.data());
    if (htdelta >= 0) {
        for (int i = 0; i < breakpoints.size(); ++i) {
            breakpoints[i] = suzerain_htstretch2(+htdelta, 1.0, breakpoints[i]);
        }
    } else {
        for (int i = 0; i < breakpoints.size(); ++i) {
            breakpoints[i] = suzerain_htstretch1(-htdelta, 1.0, breakpoints[i]);
        }
    }
    breakpoints = (right - left) * breakpoints + left;

    // Generate the B-spline workspace based on order and breakpoints
    // Maximum non-trivial derivative operators included
    b = make_shared<bspline>(k, bspline::from_breakpoints(),
                             breakpoints.size(), breakpoints.data());
    assert(b->n() == ndof);
    cop.reset(new bsplineop(*b, k-2, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE));
    assert(cop->n() == ndof);

    return 0;
}

void save(const esio_handle h,
          const shared_ptr<bspline>& b,
          const shared_ptr<bsplineop>& cop,
          const shared_ptr<bsplineop>& gop)
{
    // Ensure we were handed the appropriate discrete operators
    SUZERAIN_ENSURE(cop->get()->method == SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE);
    SUZERAIN_ENSURE(gop->get()->method == SUZERAIN_BSPLINEOP_GALERKIN_L2);

    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);

    DEBUG0(who, "Saving B-spline knot details");

    ArrayXr buf(b->nknot());

    for (int i = 0; i < b->nknot(); ++i) buf[i] = b->knot(i);
    esio_line_establish(h, b->nknot(), 0, (procid == 0 ? b->nknot() : 0));
    esio_line_write(h, "knots", buf.data(), 0,
            "Knots used to build B-spline basis");

    for (int i = 0; i < b->nbreak(); ++i) buf[i] = b->breakpoint(i);
    esio_line_establish(h, b->nbreak(), 0, (procid == 0 ? b->nbreak() : 0));
    esio_line_write(h, "breakpoints_y", buf.data(), 0,
            "Breakpoint locations used to build wall-normal B-spline basis");

    for (int i = 0; i < b->n(); ++i) buf[i] = b->collocation_point(i);
    esio_line_establish(h, b->n(), 0, (procid == 0 ? b->n() : 0));
    esio_line_write(h, "collocation_points_y", buf.data(), 0,
            "Collocation points used to build wall-normal discrete operators");

    b->integration_coefficients(0, buf.data());
    esio_line_establish(h, b->n(), 0, (procid == 0 ? b->n() : 0));
    esio_line_write(h, "integration_weights", buf.data(), 0,
            "Integrate by dotting B-spline coefficients against weights");

    char name[8]      = {};
    char comment[127] = {};

    for (int k = 0; k <= cop->nderiv(); ++k) {
        snprintf(name, sizeof(name), "Dy%dT", k);
        snprintf(comment, sizeof(comment),
                "Wall-normal derivative trans(Dy%d(i,j)) = D%dT[j,ku+i-j] for"
                " 0 <= j < n, max(0,j-ku-1) <= i < min(m,j+kl)", k, k);
        const int lda = cop->ku(k) + 1 + cop->kl(k);
        esio_plane_establish(h,
                cop->n(), 0, (procid == 0 ? cop->n() : 0),
                lda,      0, (procid == 0 ? lda          : 0));
        esio_plane_write(h, name, cop->D_T(k), 0, 0, comment);
        esio_attribute_write(h, name, "kl", cop->kl(k));
        esio_attribute_write(h, name, "ku", cop->ku(k));
        esio_attribute_write(h, name, "m",  cop->n());
        esio_attribute_write(h, name, "n",  cop->n());
    }

    DEBUG0(who, "Saving B-spline Galerkin L2 derivative operators");

    for (int k = 0; k <= gop->nderiv(); ++k) {
        snprintf(name, sizeof(name), "Gy%dT", k);
        snprintf(comment, sizeof(comment),
                "Wall-normal Galerkin L2 trans(Gy%d(i,j)) = G%dT[j,ku+i-j] for"
                " 0 <= j < n, max(0,j-ku-1) <= i < min(m,j+kl)", k, k);
        const int lda = gop->ku(k) + 1 + gop->kl(k);
        esio_plane_establish(h,
                gop->n(), 0, (procid == 0 ? gop->n() : 0),
                lda,      0, (procid == 0 ? lda          : 0));
        esio_plane_write(h, name, gop->D_T(k), 0, 0, comment);
        esio_attribute_write(h, name, "kl", gop->kl(k));
        esio_attribute_write(h, name, "ku", gop->ku(k));
        esio_attribute_write(h, name, "m",  gop->n());
        esio_attribute_write(h, name, "n",  gop->n());
    }
}

/**
 * Read an ESIO \c linev of data into line from the first possible named
 * location.  Argument \c first is mutated to return the successful location
 * name.  No suitable location may be detected by checking if <tt>first ==
 * last</tt> on return.
 */
template <typename ForwardIterator>
static void load_linev(
        const esio_handle h,
        ArrayXr &line,
        ForwardIterator& first,
        const ForwardIterator& last)
{
    for ( ; first != last; ++first ) {
        int length;
        int ncomponents;
        if (0 == esio_line_sizev(h, *first, &length, &ncomponents)) {
            line.resize(length, ncomponents);
            esio_line_establish(h, length, 0, length);
            esio_line_readv(h, *first, line.data(), 0);
            return;
        }
    }
}

/**
 * Read an ESIO \c line of data into line from the first possible named
 * location.  Argument \c first is mutated to return the successful location
 * name.  No suitable location may be detected by checking if <tt>first ==
 * last</tt> on return.
 */
template <typename ForwardIterator>
static void load_line(
        const esio_handle h,
        ArrayXr &line,
        ForwardIterator& first,
        const ForwardIterator& last)
{
    for ( ; first != last; ++first ) {
        int length;
        if (0 == esio_line_size(h, *first, &length)) {
            line.resize(length);
            esio_line_establish(h, length, 0, length);
            esio_line_read(h, *first, line.data(), 0);
            return;
        }
    }
}


real_t load(const esio_handle h,
            shared_ptr<bspline>& b,
            shared_ptr<bsplineop>& cop)
{
    using std::abs;
    using std::max;

    real_t abserr = std::numeric_limits<real_t>::quiet_NaN();

    DEBUG0(who, "Loading B-spline workspaces based on restart contents");

    // All ranks load B-spline order
    int k;
    esio_line_establish(h, 1, 0, 1);
    esio_line_read(h, "k", &k, 0);

    // htdelta is ignored

    // knots are ignored

    // All ranks load B-spline breakpoints (as best effort attempt)
    ArrayXr breakpoints;
    array<const char *,2> breakpoints_locs = {{
        "breakpoints_y", "breakpoints"
    }};
    const char **breakpoints_loc = breakpoints_locs.begin();
    load_line(h, breakpoints, breakpoints_loc, breakpoints_locs.end());
    const bool breakpoints_found = (breakpoints_loc != breakpoints_locs.end());

    // All ranks load B-spline collocation points (as best effort attempt)
    ArrayXr colpoints;
    array<const char *,2> colpoints_locs = {{
        "collocation_points_y", "collocation_points"
    }};
    const char **colpoints_loc = colpoints_locs.begin();
    load_line(h, colpoints, colpoints_loc, colpoints_locs.end());
    const bool colpoints_found = (colpoints_loc != colpoints_locs.end());

    // Generally, prefer basis to be formed using collocation points...
    bool abscissae_veto_breakpoints = true;

    // ...unless loaded breakpoints reproduce collocation points very closely.
    // Required because repeated basis calculations at restart not idempotent.
    if (breakpoints_found) {

        b = make_shared<bspline>(k, bspline::from_breakpoints(),
                                 breakpoints.size(), breakpoints.data());

        if (colpoints_found && b->n() == colpoints.size()) {
            double e = 0;
            for (int i = 0; i < colpoints.size(); ++i)
                e = max(e, abs(b->collocation_point(i) - colpoints[i]));

            DEBUG0(who, "Max difference between breakpoint-computed and loaded"
                   " collocation points is " << e);

            if (e < bsplines_distinct_distance)
                abscissae_veto_breakpoints = false;
        }
    }

    if (colpoints_found && abscissae_veto_breakpoints) {
        DEBUG0(who,
               "Collocation points from restart used to build B-spline basis");
        b = make_shared<bspline>(k, bspline::from_abscissae(),
                                 colpoints.size(), colpoints.data(), &abserr);
        DEBUG0(who, "Computed B-spline basis has Greville abscissae abserr of "
               << abserr);
    }

    // Ensure we did get B-spline workspace from the above logic
    if (!b) {
        SUZERAIN_ERROR_VAL("Could not load B-spline workspace from restart",
                SUZERAIN_EFAILED, abserr);
    }

    // Construct B-spline operator workspace from the B-spline workspace
    cop.reset(new bsplineop(
                *b, k-2, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE));

    return abserr;
}

void save_time(const esio_handle h,
               real_t time)
{
    // Root writes details
    int rank;
    esio_handle_comm_rank(h, &rank);
    esio_line_establish(h, 1, 0, (rank == 0) ? 1 : 0);

    esio_line_write(h, "t", &time, 0, "Simulation physical time");

    DEBUG0(who, "Saved simulation time " << time);
}

void load_time(const esio_handle h,
               real_t &time)
{
    // All ranks read details
    esio_line_establish(h, 1, 0, 1);

    esio_line_read(h, "t", &time, 0);

    DEBUG0(who, "Loaded simulation time " << time);
}

// Pooling employed in allocate_padded_state implementations
typedef boost::ptr_map<
        size_t, coalescing_pool<complex_t>
    > padded_state_pools_type;
static padded_state_pools_type padded_state_pools;

template<>
contiguous_state<4,complex_t>* allocate_padded_state(
           const size_t howmany_fields,
           const pencil_grid& dgrid)
{
    typedef coalescing_pool<complex_t> pool_type;

    // Contiguous number of complex_t values necessary to store one field
    // This is sufficient field-to-field padding to allow P3DFFTification
    const size_t blocksize = dgrid.local_wave_storage();

    // Find or create the coalescing_pool matching blocksize
    padded_state_pools.find(blocksize);
    padded_state_pools_type::iterator it = padded_state_pools.find(blocksize);
    if (it == padded_state_pools.end()) {
        std::auto_ptr<pool_type> tmp(new pool_type(blocksize));
        it = padded_state_pools.insert(blocksize, tmp).first;
    }

    // Construct a shared_range for howmany_fields from the pool instance
    // shared_range given boost::bind-based Deleter to invoke release()
    pool_type::blocks blocks = it->second->acquire(howmany_fields);
    shared_range<complex_t> storage(blocks.begin(), blocks.end(),
                                    boost::bind(&pool_type::release,
                                                boost::ref(*(it->second)),
                                                blocks));

    // Create instance using provided storage
    contiguous_state<4,complex_t> * const retval =
        new contiguous_state<4,complex_t>(
            storage,
            to_yxz(howmany_fields, dgrid.local_wave_extent),
            prepend(dgrid.local_wave_storage(),
                    strides_cm(to_yxz(dgrid.local_wave_extent)))
        );

    return retval;
}

} // end namespace suzerain

} // end namespace support
