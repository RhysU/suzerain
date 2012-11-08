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
// driver_base.hpp: Application driver logic spanning multiple applications
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include "driver_base.hpp"

#include <suzerain/countof.h>
#include <suzerain/error.h>
#include <suzerain/format.hpp>
#include <suzerain/l2.hpp>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/os.h>
#include <suzerain/version.hpp>

#include "logging.hpp"
#include "support.hpp"

namespace suzerain {

namespace support {

driver_base::driver_base(
        const std::string &application_synopsis,
        const std::string &description,
        const std::string &revstr)
    : application_base(application_synopsis,
                       "[RESTART-FILE]",
                       description,
                       revstr)
    , restart(make_shared<restart_definition>(
                /* metadata    */ "metadata.h5.XXXXXX",
                /* uncommitted */ "uncommitted.h5.XXXXXX",
                /* destination */ "restart#.h5",
                /* retain      */ 1,
                /* dt          */ 0,
                /* nt          */ 0))
    , statsdef(make_shared<statistics_definition>(
                /* destination */ "sample#.h5"))
    , timedef(make_shared<time_definition>(
                /* advance_dt  */ 0,
                /* advance_nt  */ 0,
                /* advance_wt  */ 0,
                /* status_dt   */ 0,
                /* status_nt   */ 0,
                /* min_dt      */ 1e-8,
                /* max_dt      */ 1))
    , soft_teardown(false)
    , log_status_L2_show_header(false)
    , log_status_bulk_show_header(false)
    , wtime_load_state(std::numeric_limits<double>::quiet_NaN())
    , wtime_advance_start(std::numeric_limits<double>::quiet_NaN())
    , last_status_nt(std::numeric_limits<std::size_t>::max())
    , last_restart_saved_nt(std::numeric_limits<std::size_t>::max())
{
    std::fill(signal_received.begin(), signal_received.end(), 0);
}

std::string
driver_base::log4cxx_config()
{
    return super::log4cxx_config() + // Append to the default configuration
        "\n"
        "## Collect \"bulk\" messages into bulk.dat mimicking LOG file behavior\n"
        "log4j.logger.bulk=INHERITED, BULK\n"
        "log4j.appender.BULK=${log4j.appender.LOG}\n"
        "log4j.appender.BULK.filename=bulk.dat\n"
        "log4j.appender.BULK.append=${log4j.appender.LOG.append}\n"
        "log4j.appender.BULK.layout=${log4j.appender.LOG.layout}\n"
        "log4j.appender.BULK.layout.ConversionPattern=${log4j.appender.LOG.layout.ConversionPattern}\n"
        "\n"
        "## Collect \"L2.mean\" messages into L2.mean.dat mimicking LOG file behavior\n"
        "log4j.logger.L2.mean=INHERITED, L2MEAN\n"
        "log4j.appender.L2MEAN=${log4j.appender.LOG}\n"
        "log4j.appender.L2MEAN.filename=L2.mean.dat\n"
        "log4j.appender.L2MEAN.append=${log4j.appender.LOG.append}\n"
        "log4j.appender.L2MEAN.layout=${log4j.appender.LOG.layout}\n"
        "log4j.appender.L2MEAN.layout.ConversionPattern=${log4j.appender.LOG.layout.ConversionPattern}\n"
        "\n"
        "## Collect \"rms.fluct\" messages into rms.fluct.dat mimicking LOG file behavior\n"
        "log4j.logger.rms.fluct=INHERITED, RMSFLUCT\n"
        "log4j.appender.RMSFLUCT=${log4j.appender.LOG}\n"
        "log4j.appender.RMSFLUCT.filename=rms.fluct.dat\n"
        "log4j.appender.RMSFLUCT.append=${log4j.appender.LOG.append}\n"
        "log4j.appender.RMSFLUCT.layout=${log4j.appender.LOG.layout}\n"
        "log4j.appender.RMSFLUCT.layout.ConversionPattern=${log4j.appender.LOG.layout.ConversionPattern}\n"
    ;
}

std::vector<std::string>
driver_base::initialize(int argc, char **argv)
{
    // Add problem definitions to options
    options.add_definition(*restart );
    options.add_definition(*statsdef);
    options.add_definition(*timedef );
    options.add_definition( sigdef  );

    // Add additional standalone options
    // TODO

    // Process incoming arguments
    std::vector<std::string> positional = super::initialize(argc, argv);

    return positional;
}

driver_base::~driver_base()
{

    // Remove the metadata file.
    // Preserve restart->uncommitted as it may help post mortem debugging.
    if (mpi::comm_rank(MPI_COMM_WORLD) == 0) {
        if (0 == unlink(restart->metadata.c_str())) {
            DEBUG("Cleaned up temporary file " << restart->metadata);
        } else {
            WARN("Error cleaning up temporary file " << restart->metadata);
        }
    }
}

bool
driver_base::log_status(
        const real_t t,
        const std::size_t nt)
{
    // Notice collective operations are never inside logging macros!

    using std::max;
    using std::floor;
    using std::log10;

    // Defensively avoid multiple invocations with no intervening changes
    if (last_status_nt == nt) {
        DEBUG0("Cowardly refusing to repeatedly show status at nt = " << nt);
        return true;
    }

    SUZERAIN_TIMER_SCOPED("log_status");

    // Build time- and timestep-specific status timeprefix.
    // Precision computations ensure multiple status lines minimally distinct
    std::ostringstream oss;
    real_t np = 0;
    if (timedef->status_dt > 0) {
        np = max(np, -floor(log10(timedef->status_dt)));
    }
    if (timedef->status_nt > 0) {
        np = max(np, -floor(log10(timedef->min_dt * timedef->status_nt)) + 1);
    }
    if (np > 0) {
        oss.setf(std::ios::fixed, std::ios::floatfield);
        const std::streamsize oldprec = oss.precision(np);
        oss << t;
        oss.precision(oldprec);
        oss.unsetf(std::ios::fixed);
    } else {
        oss << t;
    }
    oss << ' ' << std::setw(7) << nt;
    const std::string timeprefix = oss.str();

    // Log information about the various quantities of interest
    log_status_bulk(timeprefix);
    log_status_L2(timeprefix);
    log_status_specific_boundary_state(timeprefix);

    // Permit subclasses to dump arbitrary status information.  E.g. MMS error
    const bool retval = log_status_hook(timeprefix, t, nt);

    last_status_nt = nt; // Maintain last status time step

    return retval;
}

void
driver_base::log_status_L2(
        const std::string& timeprefix,
        const char * const name_L2,
        const char * const name_rms)
{
    // Avoid computational cost when logging is disabled
    logging::logger_type log_L2  = logging::get_logger(name_L2);
    logging::logger_type log_rms = logging::get_logger(name_rms);
    if (!INFO0_ENABLED(log_L2) && !INFO0_ENABLED(log_rms)) return;

    // Show headers only on first invocation
    std::ostringstream msg;
    if (log_status_L2_show_header) {
        msg << timeprefix;
        for (size_t k = 0; k < fields.size(); ++k)
            msg << ' ' << std::setw(fullprec<>::width) << fields[k].identifier;
        INFO0(log_L2, msg.str());
        INFO0(log_rms, msg.str());
        msg.str("");
        log_status_L2_show_header = false;
    }

    // Collective computation of the L_2 norms
    state_nonlinear->assign(*state_linear);
    const std::vector<field_L2> result
        = compute_field_L2(*state_nonlinear, *grid, *dgrid, *gop);

    // Build and log L2 of mean conserved state
    msg << timeprefix;
    for (size_t k = 0; k < result.size(); ++k) {
        msg << ' ' << fullprec<>(result[k].mean());
    }
    INFO0(log_L2, msg.str());

    // Build and log root-mean-squared-fluctuations of conserved state
    // RMS fluctuations are a scaling factor away from L2 fluctuations
    const real_t rms_coeff = 1/std::sqrt(grid->L.x()*grid->L.y()*grid->L.z());
    msg.str("");
    msg << timeprefix;
    for (size_t k = 0; k < result.size(); ++k) {
        msg << ' ' << fullprec<>(rms_coeff*result[k].fluctuating());
    }
    INFO0(log_rms, msg.str());
}

void
driver_base::log_status_bulk(
        const std::string& timeprefix)
{
    // Only continue on the rank housing the zero-zero modes...
    if (!dgrid->has_zero_zero_modes()) return;

    // ...and when logging is enabled.  Notice INFO not INFO0 is used.
    logging::logger_type bulk_state = logging::get_logger("bulk.state");
    if (!INFO_ENABLED(bulk_state)) return;

    // Show headers only on first invocation
    std::ostringstream msg;
    if (log_status_bulk_show_header) {
        msg << timeprefix;
        for (size_t k = 0; k < fields.size(); ++k)
            msg << ' ' << std::setw(fullprec<>::width) << fields[k].identifier;
        INFO0(bulk_state, msg.str());
        msg.str("");
        log_status_bulk_show_header = false;
    }

    // Compute operator for finding bulk quantities from coefficients
    VectorXr bulkcoeff(b->n());
    b->integration_coefficients(0, bulkcoeff.data());
    bulkcoeff /= grid->L.y();

    // Prepare the status message and log it
    msg << timeprefix;
    for (size_t k = 0; k < state_linear->shape()[0]; ++k) {
        Map<VectorXc> mean(
                (*state_linear)[k].origin(), state_linear->shape()[1]);
        msg << ' ' << fullprec<>(bulkcoeff.dot(mean.real()));
    }
    INFO(bulk_state, msg.str());
}

void
driver_base::log_status_specific_boundary_state(
        const std::string& timeprefix)
{
    // Only continue on the rank housing the zero-zero modes.
    if (!dgrid->has_zero_zero_modes()) return;

    logging::logger_type nick[2] = { logging::get_logger("bc.lower"),
                                     logging::get_logger("bc.upper")  };

    // Indices at the lower and upper walls.  Use that bc collocation point
    // values are nothing but the first and last B-spline coefficient values.
    size_t bc[2] = { 0, state_linear->shape()[1] - 1 };

    // Message lists rho, u, v, w, and total energy at walls
    for (size_t l = 0; l < SUZERAIN_COUNTOF(bc); ++l) {

        // Avoid computational cost when logging is disabled
        if (!DEBUG_ENABLED(nick[l])) continue;

        std::ostringstream msg;
        msg << timeprefix;

        const real_t rho = ((*state_linear)[ndx::rho][bc[l]][0][0]).real();
        for (size_t k = 0; k < fields.size(); ++k) {
            real_t val = (k == ndx::rho)
                       ? rho
                       : ((*state_linear)[k][bc[l]][0][0]).real() / rho;
            msg << ' ' << fullprec<>(val);
        }
        DEBUG(nick[l], msg.str());
    }
}

void
driver_base::load_restart(esio_handle esioh, real_t& t)
{
    SUZERAIN_ENSURE(grid);
    SUZERAIN_ENSURE(dgrid);

    // TODO Load everything
    // TODO Permit loading decomposition data too?

    support::load_time(esioh, t);
}

void
driver_base::save_restart(esio_handle esioh, const real_t t)
{
    SUZERAIN_ENSURE(grid);

    // TODO Implement copy of metadata file

    support::store_time(esioh, t);
}

bool
driver_base::log_status_hook(
            const std::string& timeprefix,
            const real_t t,
            const std::size_t nt)
{
    SUZERAIN_UNUSED(timeprefix);
    SUZERAIN_UNUSED(t);
    SUZERAIN_UNUSED(nt);
    return true;
}

// Initialized to zero indicating no signals have been received
driver_base::atomic_signal_received_t atomic_signal_received = {{/*0*/}};

void
driver_base::process_signal(
        const int sig)
{
    // Strictly speaking this handler performs too much work.  The design
    // choice was to have this extra work done on the (rare) signal receipt
    // rather than on the (frequent) polling of signal receipt status.

    std::vector<int>::iterator end;

    // Determine if we should output status due to the signal
    end = sigdef.status.end();
    if (std::find(sigdef.status.begin(), end, sig) != end) {
        atomic_signal_received[0] = sig;
    }

    // Determine if we should write a restart due to the signal
    end = sigdef.restart.end();
    if (std::find(sigdef.restart.begin(), end, sig) != end) {
        atomic_signal_received[1] = sig;
    }

    // Determine if we should tear down the simulation due to the signal
    end = sigdef.teardown.end();
    if (std::find(sigdef.teardown.begin(), end, sig) != end) {
        atomic_signal_received[2] = sig;
    }

    // atomic_signal_received[3] handled outside this routine

    // Determine if we should compute and write statistics due to the signal
    end = sigdef.statistics.end();
    if (std::find(sigdef.statistics.begin(), end, sig) != end) {
        atomic_signal_received[4] = sig;
    }
}

} // end namespace support

} // end namespace suzerain
