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
// driver.hpp: Application driver logic spanning multiple applications
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include "driver.hpp"

#ifdef HAVE_UNDERLING
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <underling/underling.h>
#include <underling/error.h>
#endif

#include <suzerain/mpi_datatype.hpp>
#include <suzerain/mpi.hpp>

#include "logging.hpp"
#include "support.hpp"

namespace suzerain {

namespace support {

Driver::Driver()
    : grid(),
      fftwdef( /* rigor_fft   */ fftw::measure,
               /* rigor_mpi   */ fftw::estimate),
      restart( /* metadata    */ "metadata.h5.XXXXXX",
               /* uncommitted */ "uncommitted.h5.XXXXXX",
               /* destination */ "restart#.h5",
               /* retain      */ 1,
               /* dt          */ 0,
               /* nt          */ 0),
      statsdef(/* destination */ "sample#.h5"),
      timedef( /* advance_dt  */ 0,
               /* advance_nt  */ 0,
               /* advance_wt  */ 0,
               /* status_dt   */ 0,
               /* status_nt   */ 0,
               /* min_dt      */ 1e-8,
               /* max_dt      */ 1),
      b(),
      bop(),
      gop(),
      dgrid(),
      state_linear(),
      state_nonlinear(),
      esioh(NULL),
      soft_teardown(false),
      show_header_information_L2(false),
      show_header_information_bulk(false),
      last_status_nt(std::numeric_limits<std::size_t>::max()),
      last_restart_saved_nt(std::numeric_limits<std::size_t>::max())
{
    std::fill(signal_received.begin(), signal_received.end(), 0);
}

Driver::~Driver()
{

    // Remove the metadata file.
    // Preserve restart.uncommitted as it may help post mortem debugging.
    if (suzerain::mpi::comm_rank(MPI_COMM_WORLD) == 0) {
        if (0 == unlink(restart.metadata.c_str())) {
            DEBUG("Cleaned up temporary file " << restart.metadata);
        } else {
            WARN("Error cleaning up temporary file " << restart.metadata);
        }
    }

    // Finalize any ESIO handle in use
    if (esioh)
        esio_handle_finalize(esioh);

#ifdef HAVE_UNDERLING
    underling_cleanup();
#endif
}

bool Driver::log_status(const real_t t, const std::size_t nt)
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

    // Build time- and timestep-specific status prefix.
    // Precision computations ensure multiple status lines minimally distinct
    std::ostringstream oss;
    real_t np = 0;
    if (timedef.status_dt > 0) {
        np = max(np, -floor(log10(timedef.status_dt)));
    }
    if (timedef.status_nt > 0) {
        np = max(np, -floor(log10(timedef.min_dt * timedef.status_nt)) + 1);
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
    information_bulk(timeprefix);
    information_L2(timeprefix);
    information_boundary_state(timeprefix);

    // Permit subclasses to dump arbitrary status information.  E.g. MMS error
    const bool retval = information_extended(timeprefix, t, nt);

    last_status_nt = nt; // Maintain last status time step

    return retval;
}

// Initialized to zero indicating no signals have been received
Driver::atomic_signal_received_t atomic_signal_received = {{/*0*/}};

void Driver::process_signal(int sig)
{
    // Strictly speaking this handler performs too much work.  The design
    // choice was to have this extra work done on the (rare) signal receipt
    // rather than on the (frequent) polling of signal receipt status.

    std::vector<int>::const_iterator end;

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
