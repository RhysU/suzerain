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
      sigdef(),
      b(),
      bop(),
      gop(),
      dgrid(),
      state_linear(),
      state_nonlinear(),
      esioh(NULL)
{
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

} // end namespace support

} // end namespace suzerain
