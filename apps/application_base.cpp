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

#include <suzerain/countof.h>
#include <suzerain/error.h>
#include <suzerain/format.hpp>
#include <suzerain/l2.hpp>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/os.h>
#include <suzerain/pre_gsl.h>
#include <suzerain/state.hpp>
#include <suzerain/version.hpp>

#include "logging.hpp"
#include "support.hpp"

namespace suzerain {

namespace support {

application_base::application_base(
        const std::string &application_synopsis,
        const std::string &argument_synopsis,
        const std::string &description,
        const std::string &revstr)
    : revstr(revstr)
    , options(application_synopsis,
              argument_synopsis,
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

    // Record build and invocation for posterity and to aid in debugging
    std::ostringstream os;
    std::copy(argv, argv+argc, std::ostream_iterator<const char *>(os," "));
    INFO0("Invocation: " << os.str());
    INFO0("Build:      " << suzerain::version("", revstr));

#if defined(SUZERAIN_HAVE_P3DFFT) && defined(SUZERAIN_HAVE_UNDERLING)
    // Select pencil decomposition and FFT library to use (default p3dfft)
    options.conflicting_options("p3dfft", "underling");
    if (options.variables().count("underling")) {
        use_underling = true;
    } else {
        use_p3dfft = true;
    }
#endif

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
application_base::load_grid_and_operators(
        const esio_handle esioh)
{
    SUZERAIN_ENSURE(grid);

    // Possibly load the grid parameters from the restart file
    if (esioh) {
        support::load(esioh, *grid);
    }

    // Create the discrete B-spline operators
    support::create(grid->N.y(), grid->k, 0.0,
                    grid->L.y(), grid->htdelta, b, cop);
    gop.reset(new bsplineop(*b, grid->k, SUZERAIN_BSPLINEOP_GALERKIN_L2));
}

void
application_base::store_grid_and_operators(
        const esio_handle esioh)
{
    SUZERAIN_ENSURE(grid);
    SUZERAIN_ENSURE(b);
    SUZERAIN_ENSURE(cop);
    SUZERAIN_ENSURE(gop);

    support::store(esioh, *grid);
    support::store(esioh, b, cop, gop);
}

void
application_base::establish_decomposition(
        const bool output_size,
        const bool output_plan,
        const bool output_load)
{
    SUZERAIN_ENSURE(grid);

    const std::size_t nranks = suzerain::mpi::comm_size(MPI_COMM_WORLD);

    // Display information on the global degrees of freedom
    if (output_size) {
        INFO0("Number of MPI ranks:               " << nranks);
        INFO0("Grid degrees of freedom    (GDOF): " << grid->N.prod());
        INFO0("GDOF by direction           (XYZ): " << grid->N);
        INFO0("Dealiased GDOF by direction (XYZ): " << grid->dN);
    }

    // Establish the parallel decomposition and output some timing
    if (output_plan) {
        INFO0("Preparing MPI transpose and Fourier transform plans...");
    }
    double begin = MPI_Wtime();
    fftw_set_timelimit(fftwdef->plan_timelimit);
    support::wisdom_broadcast(fftwdef->plan_wisdom);
#if defined(SUZERAIN_HAVE_P3DFFT) && defined(SUZERAIN_HAVE_UNDERLING)
    if (use_p3dfft) {
        dgrid = make_shared<pencil_grid_p3dfft>(
                grid->dN, grid->P, fftwdef->rigor_fft, fftwdef->rigor_mpi);
    } else if (use_underling) {
        dgrid = make_shared<pencil_grid_underling>(
                grid->dN, grid->P, fftwdef->rigor_fft, fftwdef->rigor_mpi);
    } else {
#endif
        dgrid = make_shared<pencil_grid_default>(
                grid->dN, grid->P, fftwdef->rigor_fft, fftwdef->rigor_mpi);
#if defined(SUZERAIN_HAVE_P3DFFT) && defined(SUZERAIN_HAVE_UNDERLING)
    }
#endif
    wtime_fftw_planning = MPI_Wtime() - begin;
    if (output_plan) {
        INFO0("MPI transpose and Fourier transform planning by "
              << dgrid->implementation() << " took "
              << wtime_fftw_planning << " seconds");
        assert((grid->dN == dgrid->global_physical_extent).all());
        INFO0("Rank grid used for decomposition: " << dgrid->processor_grid);
        INFO0("Zero-zero modes located on MPI_COMM_WORLD rank "
              << dgrid->rank_zero_zero_modes);
    }
    begin = MPI_Wtime();
    support::wisdom_gather(fftwdef->plan_wisdom);
    if (output_plan) {
        INFO0("FFTW wisdom gathered and saved in additional "
              << (MPI_Wtime() - begin) << " seconds");
    }

    // Log rank-by-rank decomposition results in debug mode
    DEBUG("Local wave start      (XYZ): " << dgrid->local_wave_start);
    DEBUG("Local wave end        (XYZ): " << dgrid->local_wave_end);
    DEBUG("Local wave extent     (XYZ): " << dgrid->local_wave_extent);
    DEBUG("Local physical start  (XYZ): " << dgrid->local_physical_start);
    DEBUG("Local physical end    (XYZ): " << dgrid->local_physical_end);
    DEBUG("Local physical extent (XYZ): " << dgrid->local_physical_extent);

    // Display normalized workloads metrics relative to zero-zero workload
    if (output_load) {
        real_t sendbuf[4];
        if (dgrid->has_zero_zero_modes()) {
            sendbuf[0] = (real_t) dgrid->local_wave_extent.prod();
            sendbuf[1] = (real_t) dgrid->local_physical_extent.prod();
        }
        SUZERAIN_MPICHKR(MPI_Bcast(sendbuf, 2,
                         suzerain::mpi::datatype<real_t>(), 0,
                         MPI_COMM_WORLD));
        sendbuf[0] = dgrid->local_wave_extent.prod()     / sendbuf[0];
        sendbuf[1] = dgrid->local_physical_extent.prod() / sendbuf[1];
        sendbuf[2] = -sendbuf[0];
        sendbuf[3] = -sendbuf[1];

        real_t recvbuf[4];
        SUZERAIN_MPICHKR(MPI_Reduce(sendbuf, recvbuf, 2,
                         suzerain::mpi::datatype<real_t>(),
                         MPI_SUM, 0, MPI_COMM_WORLD));

        const real_t mean_w = recvbuf[0] / nranks;
        const real_t mean_p = recvbuf[1] / nranks;
        SUZERAIN_MPICHKR(MPI_Reduce(sendbuf, recvbuf, 4,
                         suzerain::mpi::datatype<real_t>(),
                         MPI_MIN, 0, MPI_COMM_WORLD));
        INFO0("Wave space zero-zero normalized workloads     (min/mean/max): "
              << recvbuf[0] << ", " << mean_w << ", " << -recvbuf[2]);
        INFO0("Physical space zero-zero normalized workloads (min/mean/max): "
              << recvbuf[1] << ", " << mean_p << ", " << -recvbuf[3]);
    }
}

void
application_base::establish_state_storage(
        const std::size_t linear_nfields,
        const std::size_t nonlinear_nfields)
{
    SUZERAIN_ENSURE(dgrid);

    // Deallocate any existing unpadded linear state not matching the request
    if (state_linear && state_linear->shape()[0] != linear_nfields) {
        state_linear.reset();
    }

    // Allocate the unpadded linear state to match decomposition
    if (linear_nfields) {
        state_linear = make_shared<
                    state_linear_type
                >(to_yxz(linear_nfields, dgrid->local_wave_extent));
        DEBUG("Linear state shape      (FYXZ): "
              << shape_array(*state_linear));
        DEBUG("Linear state strides    (FYXZ): "
              << strides_array(*state_linear));
    }

    // Deallocate existing transformable nonlinear state not matching request
    if (state_nonlinear && state_nonlinear->shape()[0] != nonlinear_nfields) {
        state_nonlinear.reset();
    }

    // Allocate the transformable nonlinear state to match decomposition
    if (nonlinear_nfields) {
        state_nonlinear.reset(support::allocate_padded_state<
                    state_nonlinear_type
                >(nonlinear_nfields, *dgrid));

        DEBUG("Nonlinear state shape   (FYXZ): "
              << shape_array(*state_nonlinear));
        DEBUG("Nonlinear state strides (FYXZ): "
              << strides_array(*state_nonlinear));
    }
}

void
application_base::establish_ieee_mode()
{
    DEBUG0("Establishing floating point environment from GSL_IEEE_MODE");
    mpi_gsl_ieee_env_setup(suzerain::mpi::comm_rank(MPI_COMM_WORLD));
}

} // end namespace support

} // end namespace suzerain
