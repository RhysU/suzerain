//--------------------------------------------------------------------------
//
// Copyright (C) 2008, 2009, 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
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
// support.cpp: Support logic spanning potentially many applications
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <esio/error.h>
#include <gsl/gsl_errno.h>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/coalescing_pool.hpp>
#include <suzerain/countof.h>
#include <suzerain/diffwave.hpp>
#include <suzerain/error.h>
#include <suzerain/exprparse.hpp>
#include <suzerain/htstretch.h>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/rholut.hpp>
#include <suzerain/rngstream.hpp>
#include <suzerain/shared_range.hpp>
#include <suzerain/validation.hpp>
#include <sys/file.h>

#include "logging.hpp"
#include "support.hpp"

using boost::numeric_cast;
using std::size_t;

namespace suzerain {

namespace support {

// Common configuration snippet used in multiple places just below.
// See "Configuration" at http://logging.apache.org/log4cxx/index.html
#define COMMON_CONSOLE_CONFIG                                               \
    "log4j.appender.CONSOLE=org.apache.log4j.ConsoleAppender\n"             \
    "log4j.appender.CONSOLE.layout=org.apache.log4j.PatternLayout\n"        \
    "log4j.appender.CONSOLE.layout.ConversionPattern=%-5p %8r %-10c %m%n\n"

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

void wisdom_broadcast(const std::string& wisdom_file)
{
    if (wisdom_file.empty()) return; // Short circuit if no path provided

    // Only load wisdom from disk if FFTW MPI is available via underling.
    // Otherwise every rank hits the filesystem which is an O(N) bottleneck
    // versus a fixed O(1) planning cost on each rank.
#ifdef HAVE_UNDERLING

    // If available, load wisdom from disk on rank 0 and broadcast it
    // Attempt advisory locking to reduce processes stepping on each other
    if (mpi::comm_rank(MPI_COMM_WORLD) == 0) {

        // Import any system-wide wisdom available
        fftw_import_system_wisdom();

        FILE *w = fopen(wisdom_file.c_str(), "r");
        if (w) {
            INFO0("Loading wisdom from file " << wisdom_file);
            if (flock(fileno(w), LOCK_SH)) {
                WARN0("LOCK_SH failed on wisdom file "
                      << wisdom_file << ": " << strerror(errno));
            }
            fftw_import_wisdom_from_file(w);
            if (flock(fileno(w), LOCK_UN)) {
                WARN0("LOCK_UN failed on wisdom file "
                      << wisdom_file << ": " << strerror(errno));
            }
            fclose(w);
        } else {
            WARN0("Unable to open wisdom file "
                  << wisdom_file << ": " << strerror(errno));
        }
    }
    fftw_mpi_broadcast_wisdom(MPI_COMM_WORLD);

#endif /* HAVE_UNDERLING */
}

void wisdom_gather(const std::string& wisdom_file)
{
    if (wisdom_file.empty()) return; // Short circuit if no path provided

    // Only save wisdom to disk if FFTW MPI is available via underling.
#ifdef HAVE_UNDERLING

    // If available, gather wisdom and then write to disk on rank 0
    // Attempt advisory locking to reduce processes stepping on each other
    fftw_mpi_gather_wisdom(MPI_COMM_WORLD);
    if (mpi::comm_rank(MPI_COMM_WORLD) == 0) {
        FILE *w = fopen(wisdom_file.c_str(), "w+");
        if (w) {
            INFO0("Saving wisdom to file " << wisdom_file);
            if (flock(fileno(w), LOCK_EX)) {
                WARN0("LOCK_EX failed on wisdom file "
                      << wisdom_file << ": " << strerror(errno));
            }
            fftw_export_wisdom_to_file(w);
            if (flock(fileno(w), LOCK_UN)) {
                WARN0("LOCK_UN failed on wisdom file "
                      << wisdom_file << ": " << strerror(errno));
            }
            fclose(w);
        } else {
            WARN0("Unable to open wisdom file "
                  << wisdom_file << ": " << strerror(errno));
        }
    }

#endif /* HAVE_UNDERLING */
}

void store(const esio_handle h,
           const grid_definition& grid)
{
    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);

    DEBUG0("Storing GridDefinition parameters");

    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));

    esio_line_write(h, "Lx", &grid.L.x(), 0,
            grid.options().find("Lx",false).description().c_str());

    esio_line_write(h, "Nx", &grid.N.x(), 0,
            grid.options().find("Nx",false).description().c_str());

    esio_line_write(h, "DAFx", &grid.DAF.x(), 0,
            grid.options().find("DAFx",false).description().c_str());

    esio_line_write(h, "Ly", &grid.L.y(), 0,
            grid.options().find("Ly",false).description().c_str());

    esio_line_write(h, "Ny", &grid.N.y(), 0,
            grid.options().find("Ny",false).description().c_str());

    esio_line_write(h, "k", &grid.k, 0,
            grid.options().find("k",false).description().c_str());

    esio_line_write(h, "htdelta", &grid.htdelta, 0,
            grid.options().find("htdelta",false).description().c_str());

    esio_line_write(h, "Lz", &grid.L.z(), 0,
            grid.options().find("Lz",false).description().c_str());

    esio_line_write(h, "Nz", &grid.N.z(), 0,
            grid.options().find("Nz",false).description().c_str());

    esio_line_write(h, "DAFz", &grid.DAF.z(), 0,
            grid.options().find("DAFz",false).description().c_str());

    DEBUG0("Storing wavenumber vectors for Fourier bases");
    ArrayXc cbuf(std::max(grid.N.x(), grid.N.z()));

    // Obtain wavenumbers via computing 1*(i*kx)/i
    cbuf.fill(complex_t(1,0));
    diffwave::apply(1, 0, complex_t(0,-1), cbuf.data(),
            grid.L.x(), grid.L.z(),
            1, grid.N.x(), grid.N.x(), 0, grid.N.x(), 1, 1, 0, 1);
    esio_line_establish(h, grid.N.x(), 0, (procid == 0 ? grid.N.x() : 0));
    esio_line_write(h, "kx", reinterpret_cast<real_t *>(cbuf.data()),
            2, "Wavenumbers in streamwise X direction"); // Re(cbuf)

    // Obtain wavenumbers via computing 1*(i*kz)/i
    cbuf.fill(complex_t(1,0));
    diffwave::apply(0, 1, complex_t(0,-1), cbuf.data(),
            grid.L.x(), grid.L.z(),
            1, 1, 1, 0, 1, grid.N.z(), grid.N.z(), 0, grid.N.z());
    esio_line_establish(h, grid.N.z(), 0, (procid == 0 ? grid.N.z() : 0));
    esio_line_write(h, "kz", reinterpret_cast<real_t *>(cbuf.data()),
            2, "Wavenumbers in spanwise Z direction"); // Re(cbuf)

    DEBUG0("Storing collocation point vectors for Fourier bases");
    ArrayXr rbuf;

    // Obtain collocation points in x using [-Lx/2, Lx/2]) and dN.x()
    if (grid.dN.x() > 1) {
        rbuf = ArrayXr::LinSpaced(Eigen::Sequential,
                                  grid.dN.x(), 0, grid.dN.x() - 1);
        rbuf *= grid.L.x() / grid.dN.x();
        rbuf -= grid.L.x() / 2;
    } else {
        rbuf = ArrayXr::Constant(grid.dN.x(), 0);
    }
    esio_line_establish(h, rbuf.size(), 0, (procid == 0 ? rbuf.size() : 0));
    esio_line_write(h, "collocation_points_x", rbuf.data(), 0,
            "Collocation points for the dealiased, streamwise X direction");

    // Obtain collocation points in z using [-Lz/2, Lz/2]) and dN.z()
    if (grid.dN.z() > 1) {
        rbuf = ArrayXr::LinSpaced(Eigen::Sequential,
                                  grid.dN.z(), 0, grid.dN.z() - 1);
        rbuf *= grid.L.z() / grid.dN.z();
        rbuf -= grid.L.z() / 2;
    } else {
        rbuf = ArrayXr::Constant(grid.dN.z(), 0);
    }
    esio_line_establish(h, rbuf.size(), 0, (procid == 0 ? rbuf.size() : 0));
    esio_line_write(h, "collocation_points_z", rbuf.data(), 0,
            "Collocation points for the dealiased, spanwise Z direction");
}

void load(const esio_handle h,
          grid_definition& grid)
{
    DEBUG0("Loading GridDefinition parameters");

    esio_line_establish(h, 1, 0, 1); // All ranks load

    if (!(boost::math::isnan)(grid.L.x())) {
        INFO0("Overriding grid using Lx = " << grid.L.x());
    } else {
        esio_line_read(h, "Lx", &grid.L.x(), 0);
    }

    if (grid.N.x()) {
        INFO0("Overriding grid using Nx = " << grid.N.x());
    } else {
        int value;
        esio_line_read(h, "Nx", &value, 0);
        grid.Nx(value);
    }

    if (!((boost::math::isnan)(grid.DAF.x()))) {
        INFO0("Overriding grid using DAFx = " << grid.DAF.x());
    } else {
        double factor;
        esio_line_read(h, "DAFx", &factor, 0);
        grid.DAFx(factor);
    }

    if (!(boost::math::isnan)(grid.L.y())) {
        INFO0("Overriding grid using Ly = " << grid.L.y());
    } else {
        esio_line_read(h, "Ly", &grid.L.y(), 0);
    }

    if (grid.N.y()) {
        INFO0("Overriding grid using Ny = " << grid.N.y());
    } else {
        int value;
        esio_line_read(h, "Ny", &value, 0);
        grid.Ny(value);
    }

    if (grid.k) {
        INFO0("Overriding grid using k = " << grid.k);
    } else {
        esio_line_read(h, "k", &grid.k, 0);
    }

    if (!((boost::math::isnan)(grid.htdelta))) {
        INFO0("Overriding grid using htdelta = " << grid.htdelta);
    } else {
        esio_line_read(h, "htdelta", &grid.htdelta, 0);
    }

    if (!(boost::math::isnan)(grid.L.z())) {
        INFO0("Overriding grid using Lz = " << grid.L.z());
    } else {
        esio_line_read(h, "Lz", &grid.L.z(), 0);
    }

    if (grid.N.z()) {
        INFO0("Overriding grid using Nz = " << grid.N.z());
    } else {
        int value;
        esio_line_read(h, "Nz", &value, 0);
        grid.Nz(value);
    }

    if (!((boost::math::isnan)(grid.DAF.z()))) {
        INFO0("Overriding grid using DAFz = " << grid.DAF.z());
    } else {
        double factor;
        esio_line_read(h, "DAFz", &factor, 0);
        grid.DAFz(factor);
    }
}

void store(const esio_handle h,
           const time_definition& timedef)
{
    DEBUG0("Storing TimeDefinition parameters");

    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);

    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));

    esio_line_write(h, "evmagfactor", &timedef.evmagfactor, 0,
            timedef.options().find("evmagfactor",false).description().c_str());
}

void load(const esio_handle h,
          time_definition& timedef)
{
    DEBUG0("Loading TimeDefinition parameters");

    esio_line_establish(h, 1, 0, 1); // All ranks load

    // evmagfactor not present in pre-revision 24043 restart files
    if (!(boost::math::isnan)(timedef.evmagfactor)) {
        INFO0("Overriding timedef using evmagfactor = " << timedef.evmagfactor);
    } else if (ESIO_NOTFOUND == esio_line_size(h, "evmagfactor", NULL)) {
        timedef.evmagfactor = 0.72;
        INFO0("Employing default evmagfactor = " << timedef.evmagfactor);
    } else {
        esio_line_read(h, "evmagfactor", &timedef.evmagfactor, 0);
    }
}

real_t create(const int ndof,
              const int k,
              const double left,
              const double right,
              const double htdelta,
              boost::shared_ptr<bspline>& b,
              boost::shared_ptr<bsplineop>& bop)
{
    INFO0("Creating B-spline basis of order " << k
          << " on [" << left << ", " << right << "] with "
          << ndof << " DOF stretched per htdelta " << htdelta);

////FIXME: Knot vectors are non-increasing for moderate htdelta
/// FIXME: See https://savannah.gnu.org/bugs/index.php?34361
////// Compute collocation point locations using ndof and htdelta
////Eigen::ArrayXd abscissae(ndof);
////math::linspace(0.0, 1.0, abscissae.size(), abscissae.data());
////for (int i = 0; i < abscissae.size(); ++i) {
////    abscissae[i] = suzerain_htstretch2(htdelta, 1.0, abscissae[i]);
////}
////abscissae = (right - left) * abscissae + left;
////
////// Generate the B-spline workspace based on order and abscissae
////// Maximum non-trivial derivative operators included
////double abserr;
////b = boost::make_shared<bspline>(
////        k, bspline::from_abscissae(),
////        abscissae.size(), abscissae.data(), &abserr);
////assert(b->n() == ndof);
////bop.reset(new bsplineop(
////            *b, k-2, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE));
////assert(bop->n() == ndof);
////
////INFO0("Created B-spline basis has Greville abscissae abserr of " << abserr);
////
////return abserr;

    // Compute breakpoint point locations using ndof and htdelta
    Eigen::ArrayXd breakpoints(ndof - k + 2);
    math::linspace(0.0, 1.0, breakpoints.size(), breakpoints.data());
    for (int i = 0; i < breakpoints.size(); ++i) {
        breakpoints[i] = suzerain_htstretch2(htdelta, 1.0, breakpoints[i]);
    }
    breakpoints = (right - left) * breakpoints + left;

    // Generate the B-spline workspace based on order and breakpoints
    // Maximum non-trivial derivative operators included
    b = boost::make_shared<bspline>(k, bspline::from_breakpoints(),
                                    breakpoints.size(), breakpoints.data());
    assert(b->n() == ndof);
    bop.reset(new bsplineop(*b, k-2, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE));
    assert(bop->n() == ndof);

    return 0;
}

void store(const esio_handle h,
           const boost::shared_ptr<bspline>& b,
           const boost::shared_ptr<bsplineop>& bop,
           const boost::shared_ptr<bsplineop>& gop)
{
    // Ensure we were handed the appropriate discrete operators
    assert(bop->get()->method == SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE);
    assert(gop->get()->method == SUZERAIN_BSPLINEOP_GALERKIN_L2);

    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);

    DEBUG0("Storing B-spline knot details");

    Eigen::ArrayXd buf(b->nknot());

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

    for (int k = 0; k <= bop->nderiv(); ++k) {
        snprintf(name, sizeof(name), "Dy%dT", k);
        snprintf(comment, sizeof(comment),
                "Wall-normal derivative trans(Dy%d(i,j)) = D%dT[j,ku+i-j] for"
                " 0 <= j < n, max(0,j-ku-1) <= i < min(m,j+kl)", k, k);
        const int lda = bop->ku(k) + 1 + bop->kl(k);
        esio_plane_establish(h,
                bop->n(), 0, (procid == 0 ? bop->n() : 0),
                lda,      0, (procid == 0 ? lda          : 0));
        esio_plane_write(h, name, bop->D_T(k), 0, 0, comment);
        esio_attribute_write(h, name, "kl", bop->kl(k));
        esio_attribute_write(h, name, "ku", bop->ku(k));
        esio_attribute_write(h, name, "m",  bop->n());
        esio_attribute_write(h, name, "n",  bop->n());
    }

    DEBUG0("Storing B-spline Galerkin L2 derivative operators");

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

real_t load(const esio_handle h,
            boost::shared_ptr<bspline>& b,
            boost::shared_ptr<bsplineop>& bop)
{
    using std::abs;
    using std::max;

    real_t abserr = std::numeric_limits<real_t>::quiet_NaN();

    DEBUG0("Loading B-spline workspaces based on restart contents");

    // All ranks load B-spline order
    int k;
    esio_line_establish(h, 1, 0, 1);
    esio_line_read(h, "k", &k, 0);

    // htdelta is ignored

    // knots are ignored

    // All ranks load B-spline breakpoints (as best effort attempt)
    ArrayXr breakpoints;
    boost::array<const char *,2> breakpoints_locs = {{
        "breakpoints_y", "breakpoints"
    }};
    const char **breakpoints_loc = breakpoints_locs.begin();
    load_line(h, breakpoints, breakpoints_loc, breakpoints_locs.end());
    const bool breakpoints_found = (breakpoints_loc != breakpoints_locs.end());

    // All ranks load B-spline collocation points (as best effort attempt)
    ArrayXr colpoints;
    boost::array<const char *,2> colpoints_locs = {{
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

        b = boost::make_shared<bspline>(
                k, bspline::from_breakpoints(),
                breakpoints.size(), breakpoints.data());

        if (colpoints_found && b->n() == colpoints.size()) {
            double e = 0;
            for (int i = 0; i < colpoints.size(); ++i)
                e = max(e, abs(b->collocation_point(i) - colpoints[i]));

            DEBUG0("Max difference between breakpoint-computed and loaded"
                   " collocation points is " << e);

            if (e < bsplines_distinct_distance)
                abscissae_veto_breakpoints = false;
        }
    }

    if (colpoints_found && abscissae_veto_breakpoints) {
        DEBUG0("Collocation points from restart used to build B-spline basis");
        b = boost::make_shared<bspline>(
                k, bspline::from_abscissae(),
                colpoints.size(), colpoints.data(), &abserr);
        DEBUG0("Computed B-spline basis has Greville abscissae abserr of "
               << abserr);
    }

    // Ensure we did get B-spline workspace from the above logic
    if (!b) {
        SUZERAIN_ERROR_VAL("Could not load B-spline workspace from restart",
                SUZERAIN_EFAILED, abserr);
    }

    // Construct B-spline operator workspace from the B-spline workspace
    bop.reset(new bsplineop(
                *b, k-2, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE));

    return abserr;
}

void store_time(const esio_handle h,
                real_t time)
{
    // Root writes details
    int rank;
    esio_handle_comm_rank(h, &rank);
    esio_line_establish(h, 1, 0, (rank == 0) ? 1 : 0);

    esio_line_write(h, "t", &time, 0, "Simulation physical time");

    DEBUG0("Stored simulation time " << time);
}

void load_time(const esio_handle h,
               real_t &time)
{
    // All ranks read details
    esio_line_establish(h, 1, 0, 1);

    esio_line_read(h, "t", &time, 0);

    DEBUG0("Loaded simulation time " << time);
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
            boost::bind(&pool_type::release, boost::ref(*(it->second)), blocks));

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

void store_coefficients(
        const esio_handle h,
        const std::vector<field> &fields,
        const contiguous_state<4,complex_t> &swave,
        const grid_definition& grid,
        const pencil_grid& dgrid)
{
    // Ensure swave meets this routine's assumptions
    SUZERAIN_ENSURE(                  swave.shape()[0]  == fields.size()              );
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[1]) == dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[2]) == dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[3]) == dgrid.local_wave_extent.z());

    // Compute wavenumber translation logistics for X direction
    int fxb[2], fxe[2], mxb[2], mxe[2];
    inorder::wavenumber_translate(grid.N.x(),
                                  grid.dN.x(),
                                  dgrid.local_wave_start.x(),
                                  dgrid.local_wave_end.x(),
                                  fxb[0], fxe[0], fxb[1], fxe[1],
                                  mxb[0], mxe[0], mxb[1], mxe[1]);
    // X contains only positive wavenumbers => second range must be empty
    SUZERAIN_ENSURE(fxb[1] == fxe[1]);
    SUZERAIN_ENSURE(mxb[1] == mxe[1]);

    // Compute wavenumber translation logistics for Z direction
    // One or both ranges may be empty
    int fzb[2], fze[2], mzb[2], mze[2];
    inorder::wavenumber_translate(grid.N.z(),
                                  grid.dN.z(),
                                  dgrid.local_wave_start.z(),
                                  dgrid.local_wave_end.z(),
                                  fzb[0], fze[0], fzb[1], fze[1],
                                  mzb[0], mze[0], mzb[1], mze[1]);

    // Save each scalar field in turn...
    for (size_t i = 0; i < fields.size(); ++i) {

        // ...first generate a metadata comment...
        std::string comment = fields[i].description;
        comment += " stored row-major ZXY using a Fourier basis in Z stored"
                   " in-order per /kz; a Fourier basis in X stored using"
                   " Hermitian symmetry per /kx; and a B-spline basis in Y"
                   " defined by /k, /breakpoints_y, and /knots";

        // ...followed by two collective writes per field (once per Z range)
        for (int j = 0; j < 2; ++j) {

            // Collectively establish size of write across all ranks
            esio_field_establish(h,
                                 grid.N.z(),     fzb[j], (fze[j] - fzb[j]),
                                 grid.N.x()/2+1, fxb[0], (fxe[0] - fxb[0]),
                                 grid.N.y(),     0,      (     grid.N.y()));

            // Source of write is NULL for empty WRITE operations
            // Required since MultiArray triggers asserts on invalid indices
            const complex_t * src = NULL;
            if (mxb[0] != mxe[0] && mzb[j] != mze[j]) {
                src = &swave[i][0][mxb[0] - dgrid.local_wave_start.x()]
                                  [mzb[j] - dgrid.local_wave_start.z()];
            }

            // Perform collective write operation
            complex_field_write(h, fields[i].location.c_str(), src,
                                swave.strides()[3],
                                swave.strides()[2],
                                swave.strides()[1],
                                comment.c_str());
        }
    }
}

real_t distance(const bspline& a,
                const bspline& b)
{
    real_t retval = 0;
    if (a.k() != b.k() || a.n() != b.n() || a.nknot() != b.nknot()) {
        retval = std::numeric_limits<real_t>::max();
    } else {
        for (int j = 0; j < b.nknot(); ++j) {
            retval = std::max(retval, std::abs(a.knot(j) - b.knot(j)));
        }
    }
    return retval;
}

void load_coefficients(const esio_handle h,
                       const std::vector<field> &fields,
                       contiguous_state<4,complex_t> &state,
                       const grid_definition& grid,
                       const pencil_grid& dgrid,
                       const bspline& b,
                       const bsplineop& bop)
{
    typedef contiguous_state<4,complex_t> load_type;

    // Ensure local state storage meets this routine's assumptions
    SUZERAIN_ENSURE(                  state.shape()[0]  == fields.size());
    SUZERAIN_ENSURE(numeric_cast<int>(state.shape()[1]) == dgrid.global_wave_extent.y());
    SUZERAIN_ENSURE(numeric_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE(numeric_cast<int>(state.shape()[3]) == dgrid.local_wave_extent.z());

    // Obtain details on the restart field's global sizes
    int Fz, Fx, Fy, ncomponents;
    esio_field_sizev(h, fields[0].location.c_str(), &Fz, &Fx, &Fy, &ncomponents);
    SUZERAIN_ENSURE(ncomponents == 2);

    // Prepare a file-specific B-spline basis
    boost::shared_ptr<bspline> Fb;
    boost::shared_ptr<bsplineop> Fbop;
    load(h, Fb, Fbop);
    SUZERAIN_ENSURE(Fy == Fb->n());

    // Check if the B-spline basis in the file differs from ours.
    const double bsplines_dist = distance(b, *Fb);
    const bool bsplines_same = bsplines_dist < bsplines_distinct_distance;

    // Compute wavenumber translation logistics for X direction.
    // Requires turning a C2R FFT complex-valued coefficient count into a
    // real-valued coefficient count.  Further, need to preserve even- or
    // odd-ness of the coefficient count to handle, for example, Fx == 1.
    int fxb[2], fxe[2], mxb[2], mxe[2];
    inorder::wavenumber_translate(2 * (Fx - 1) + (Fx & 1),
                                  grid.dN.x(),
                                  dgrid.local_wave_start.x(),
                                  dgrid.local_wave_end.x(),
                                  fxb[0], fxe[0], fxb[1], fxe[1],
                                  mxb[0], mxe[0], mxb[1], mxe[1]);
    // X contains only positive wavenumbers => second range must be empty
    SUZERAIN_ENSURE(fxb[1] == fxe[1]);
    SUZERAIN_ENSURE(mxb[1] == mxe[1]);

    // Compute wavenumber translation logistics for Z direction
    // One or both ranges may be empty
    int fzb[2], fze[2], mzb[2], mze[2];
    inorder::wavenumber_translate(Fz,
                                  grid.dN.z(),
                                  dgrid.local_wave_start.z(),
                                  dgrid.local_wave_end.z(),
                                  fzb[0], fze[0], fzb[1], fze[1],
                                  mzb[0], mze[0], mzb[1], mze[1]);

    // Possibly prepare a tmp buffer into which to read each scalar field and a
    // factorization of b's mass matrix.  Used only when !bsplines_same.
    typedef boost::multi_array<
        complex_t, 3, blas::allocator<complex_t>::type
    > tmp_type;
    boost::scoped_ptr<tmp_type> tmp;
    boost::scoped_ptr<bsplineop_luz> mass;
    if (!bsplines_same) {
        INFO0("Differences in B-spline basis require restart projection ("
              << bsplines_dist << " >= " << bsplines_distinct_distance << ")");
        const boost::array<tmp_type::size_type,3> extent = {{
            static_cast<tmp_type::size_type>(Fy),
            state.shape()[2],
            state.shape()[3]
        }};
        tmp.reset(new tmp_type(extent, boost::fortran_storage_order()));
        mass.reset(new bsplineop_luz(bop));
        mass->factor_mass(bop);
    }

    // Load each scalar field in turn
    for (size_t i = 0; i < fields.size(); ++i) {

        // Create a view of the state for just the i-th scalar
        boost::multi_array_types::index_range all;
        boost::array_view_gen<load_type, 3>::type field
                = state[boost::indices[i][all][all][all]];

        // Clear storage prior to load to zero not-loaded coefficients
        if (bsplines_same) {
            multi_array::fill(field, 0);
        } else {
            multi_array::fill(*tmp, 0);
        }

        // Two ESIO read operations per field (once per Z range)
        for (int j = 0; j < 2; ++j) {

            // Collectively establish size of read across all ranks
            esio_field_establish(h, Fz, fzb[j], (fze[j] - fzb[j]),
                                    Fx, fxb[0], (fxe[0] - fxb[0]),
                                    Fy,      0, (            Fy));

            // Destination of read is NULL for empty READ operations.
            // Otherwise find starting point for coefficient loading.
            complex_t * dst = NULL;
            boost::array<load_type::index,3> dst_strides = {{ 0, 0, 0 }};
            if (mxb[0] != mxe[0] && mzb[j] != mze[j]) {
                const boost::array<load_type::index,3> index_list = {{
                        0,
                        mxb[0] - dgrid.local_wave_start.x(),
                        mzb[j] - dgrid.local_wave_start.z()
                }};
                if (bsplines_same) {
                    dst = &field(index_list);
                    std::copy(field.strides(), field.strides() + 3,
                              dst_strides.begin());
                } else {
                    dst = &(*tmp)(index_list);
                    std::copy(tmp->strides(), tmp->strides() + 3,
                              dst_strides.begin());
                }
            }

            // Perform collective read operation into dst
            complex_field_read(h, fields[i].location.c_str(), dst,
                               dst_strides[2], dst_strides[1], dst_strides[0]);
        }


        // If necessary, interpolate between B-spline bases.
        // Relies heavily on both bases being collocation-based.
        // This will change dramatically if we ever go the L_2 route.
        if (!bsplines_same) {

            // Step 0: Obtain collocation points for new basis
            Eigen::ArrayXd points(b.n());
            for (int k = 0; k < b.n(); ++k) points[k] = b.collocation_point(k);

            // Step 1: Affine transformation of points into old basis domain
            // Allows stretching of old range onto new range.
            const double newmin = b.collocation_point(0);
            const double newmax = b.collocation_point(b.n() - 1);
            const double oldmin = Fb->collocation_point(0);
            const double oldmax = Fb->collocation_point(Fb->n() - 1);
            // Only pay for floating point loss if strictly necessary
#pragma warning(push,disable:1572)
            if (oldmin != newmin || oldmax != newmax) {
#pragma warning(pop)
                points = (points - newmin)
                    * ((oldmax - oldmin) / (newmax - newmin)) + oldmin;
            }

            // Step 2: Evaluate old basis + coefficients at points into new
            const int jmax = numeric_cast<int>(field.shape()[1]);
            const int kmax = numeric_cast<int>(field.shape()[2]);
            for (int k = 0; k < kmax; ++k) {
                for (int j = 0; j < jmax; ++j) {
                    Fb->linear_combination(0, &(*tmp)[0][j][k], b.n(),
                                           points.data(), &field[0][j][k], 0);
                }
            }

            // Step 3: Invert for new coefficients using factored mass matrix
            for (int k = 0; k < kmax; ++k) {
                mass->solve(jmax, &field[0][0][k],
                            numeric_cast<int>(field.strides()[0]),
                            numeric_cast<int>(field.strides()[1]));
            }
        }
    }
}

} // end namespace suzerain

} // end namespace support
