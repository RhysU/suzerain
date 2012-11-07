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
// application_base.cpp: building blocks for basic application logic
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include "application_base.hpp"

#ifdef HAVE_UNDERLING
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <underling/underling.h>
#include <underling/error.h>
#endif

#include <esio/error.h>
#include <esio/esio.h>

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

application_base::application_base(
        const std::string &application_synopsis,
        const std::string &description,
        const std::string &revstr)
    : revstr(revstr)
    , options(application_synopsis,
              "FILE",
              description,
              this->revstr)
    , grid(make_shared<grid_definition>())
    , fftwdef(make_shared<fftw_definition>(/* rigor_fft */ fftw::measure,
                                           /* rigor_mpi */ fftw::estimate))
    , b()
    , cop()
    , gop()
    , dgrid()
    , state_linear()
    , state_nonlinear()
    , wtime_mpi_init(std::numeric_limits<double>::quiet_NaN())
    , wtime_fftw_planning(std::numeric_limits<double>::quiet_NaN())
#if defined(SUZERAIN_HAVE_P3DFFT) && defined(SUZERAIN_HAVE_UNDERLING)
    , use_p3dfft(false)
    , use_underling(false)
#endif
{
}

application_base::~application_base()
{
#ifdef HAVE_UNDERLING
    underling_cleanup();
#endif
}

std::vector<std::string>
application_base::initialize(int argc, char **argv)
{
#ifdef SUZERAIN_HAVE_GRVY
    grvy_timer_init(argv[0] ? argv[0] : "NULL");  // Initialize GRVY Timers
#endif
    MPI_Init(&argc, &argv);                          // Initialize MPI...
    wtime_mpi_init = MPI_Wtime();                    // Record MPI_Init time
    atexit((void (*) ()) MPI_Finalize);              // ...finalize at exit
    logging::initialize(MPI_COMM_WORLD,              // Initialize logging
                        this->log4cxx_config().c_str());
#ifdef HAVE_UNDERLING
    underling_init(&argc, &argv, 0);                 // Initialize underling...
    atexit(&underling_cleanup);                      // ...finalize at exit
#endif

    // Hook error handling into logging infrastructure
    gsl_set_error_handler(
            &support::mpi_abort_on_error_handler_gsl);
    suzerain_set_error_handler(
            &support::mpi_abort_on_error_handler_suzerain);
    esio_set_error_handler(
            &support::mpi_abort_on_error_handler_esio);
#ifdef HAVE_UNDERLING
    underling_set_error_handler(
            &support::mpi_abort_on_error_handler_underling);
#endif

    // Add problem definitions to options
    options.add_definition(*grid   );
    options.add_definition(*fftwdef);

    // Add additional standalone options
    options.add_options()
#if defined(SUZERAIN_HAVE_P3DFFT) && defined(SUZERAIN_HAVE_UNDERLING)
        ("p3dfft",    "Use P3DFFT for MPI-parallel FFTs")
        ("underling", "Use underling for MPI-parallel FFTs")
#endif
        ;

    // Process incoming arguments
    std::vector<std::string> positional = options.process(argc, argv);

#if defined(SUZERAIN_HAVE_P3DFFT) && defined(SUZERAIN_HAVE_UNDERLING)
    // Select pencil decomposition and FFT library to use (default p3dfft)
    options.conflicting_options("p3dfft", "underling");
    if (options.variables().count("underling")) {
        use_underling = true;
    } else {
        use_p3dfft = true;
    }
#endif

    // Record build and invocation for posterity and to aid in debugging
    std::ostringstream os;
    std::copy(argv, argv+argc, std::ostream_iterator<const char *>(os," "));
    INFO0("Invocation: " << os.str());
    INFO0("Build:      " << suzerain::version("", revstr));

    switch (options.verbose()) {
        case 0:                   break;
        case 1:  DEBUG0_ENABLE(); break;
        default: TRACE0_ENABLE(); break;
    }
    switch (options.verbose_all()) {
        case 0:                   break;
        case 1:  DEBUG_ENABLE();  break;
        default: TRACE_ENABLE();  break;
    }

    return positional;
}

std::string
application_base::log4cxx_config()
{
    return support::log4cxx_config;
}

void
application_base::load_grid_details(esio_handle esioh, real_t& t)
{
    SUZERAIN_ENSURE(grid);

    // Load the grid parameters from the restart file
    support::load(esioh, *grid);

    // Create the discrete B-spline operators
    support::create(grid->N.y(), grid->k, 0.0,
                    grid->L.y(), grid->htdelta, b, cop);
    gop.reset(new bsplineop(*b, grid->k, SUZERAIN_BSPLINEOP_GALERKIN_L2));

    // Load the simulation time
    support::load_time(esioh, t);
}

void
application_base::store_grid_details(esio_handle esioh, const real_t t)
{
    SUZERAIN_ENSURE(grid);
    SUZERAIN_ENSURE(b);
    SUZERAIN_ENSURE(cop);
    SUZERAIN_ENSURE(gop);

    support::store(esioh, *grid);
    support::store(esioh, b, cop, gop);
    support::store_time(esioh, t);
}

void
application_base::establish_decomposition()
{
    // FIXME Implement
}

} // end namespace support

} // end namespace suzerain
