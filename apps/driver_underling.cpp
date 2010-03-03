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
 * driver_underling.cpp: An undering library test driver
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop

#include <mpi.h>
#include <fftw3-mpi.h>
#include <log4cxx/logger.h>

// We use HPCT, if available, for some fine-scale timing.
#ifdef HAVE_HPCT
#include <hpct.h>
#endif

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include <suzerain/error.h>
#include <suzerain/fftw.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/profile.hpp>
#include <suzerain/program_options.hpp>
#include <suzerain/types.hpp>
#include <suzerain/underling.hpp>

int main(int argc, char *argv[])
{
    // Initialize external dependencies
    MPI_Init(&argc, &argv);
    MPI_Pcontrol(0);             // Disable MPI profiling on startup
    fftw_mpi_init();

    // Arrange to finalize external dependencies
    atexit((void (*) ()) MPI_Finalize);
    atexit((void (*) ()) fftw_cleanup);
    atexit((void (*) ()) fftw_mpi_cleanup);

    // Initialize logging subsystem
    const int nproc  = suzerain::mpi::comm_size(MPI_COMM_WORLD);
    const int procid = suzerain::mpi::comm_rank(MPI_COMM_WORLD);
    log4cxx::LoggerPtr logger = log4cxx::Logger::getLogger(
            suzerain::mpi::comm_rank_identifier(MPI_COMM_WORLD));

    // Process command-line options
    suzerain::ProgramOptions options;
    suzerain::problem::GridDefinition<> griddef;
    options.add_definition(griddef);
    suzerain::fftw::FFTWDefinition fftwdef;
    options.add_definition(fftwdef);
    int nrep          =  1;
    int howmany       =  1;
    int sleep_barrier = -1;
    namespace po = boost::program_options;
    // TODO Weirdly, add_options() is failing for more than 2 options...
    options.add_options()
        ("howmany", po::value<int>(&howmany)->default_value(1),
        "Number of interleaved real-valued fields to transpose");
    options.add_options()
        ("nrep", po::value<int>(&nrep)->default_value(1),
        "Number of repetitions to perform for timing purposes");
    options.add_options()
        ("sleep_barrier", po::value<int>(&sleep_barrier)->default_value(-1),
        "DEBUG: If > 0, the process rank at which to create a sleep barrier.");
    if (!procid) {
        options.process(argc, argv);
    } else {
        boost::onullstream nullstream;
        options.process(argc, argv,
                        nullstream, nullstream, nullstream, nullstream);
    }
    if (sleep_barrier >= 0) {
        suzerain::mpi::sleep_barrier(MPI_COMM_WORLD, sleep_barrier, std::cerr);
    }
#ifdef HAVE_FFTW3_THREADS
    assert(fftw_init_threads());
    atexit((void (*) ()) fftw_cleanup_threads);
    fftw_plan_with_nthreads(fftwdef.nthreads());
#endif

    // Create a grid and a problem
    namespace underling = suzerain::underling;
    underling::grid grid(MPI_COMM_WORLD,
                         griddef.global_extents().begin(),
                         griddef.processor_grid().begin());
    underling::problem problem(grid, howmany);

    // Allocate automatically freed, aligned storage
    const size_t local_memory = problem.local_memory();
    boost::shared_ptr<underling_real>
        data(static_cast<underling_real *>(
                fftw_malloc(local_memory*sizeof(underling_real))),
             &fftw_free);

    // Obtain some grid-wide memory usage information
    const size_t local_memory_minimum
        = underling::local_memory_minimum(grid, problem);
    const size_t local_memory_maximum
        = underling::local_memory_maximum(grid, problem);

    // Report grid geometry and memory characteristics
    if (!procid) {
        LOG4CXX_INFO(logger,
                     "Number of processors:     " << nproc);
        LOG4CXX_INFO(logger,
                     "Extents:                  " << griddef.global_extents()
                     << " by " << howmany << " scalar field(s)");
        LOG4CXX_INFO(logger,
                     "Requested processor grid: "
                     << griddef.processor_grid());
        LOG4CXX_INFO(logger,
                     "Actual processor grid:    "
                     << grid.pA_size() << " by " << grid.pB_size());
        LOG4CXX_INFO(logger,
                     "Planning rigor:           "
                     << suzerain::fftw::c_str(fftwdef.plan_rigor()));
#ifdef HAVE_FFTW3_THREADS
        LOG4CXX_INFO(logger,
                     "FFTW threads:             " << fftwdef.nthreads());
#else
        LOG4CXX_INFO(logger,
                     "FFTW threads:             Unavailable");
#endif
        LOG4CXX_INFO(logger,
                     "Number of repetitions:    " << nrep);
        LOG4CXX_INFO(logger,
                     "Process storage:          "
                     << (local_memory_minimum*sizeof(underling_real))/1024
                     << " KB min, "
                     << (local_memory_maximum*sizeof(underling_real))/1024
                     << " KB max");
    }

#ifdef HAVE_HPCT
    hpct_timer_init("underling");
    hpct_timer_begin("plan");
#endif
    // Create the plan once
    underling::plan plan(problem,
                         data.get(),
                         underling::transpose::all,
                         fftwdef.plan_rigor());
#ifdef HAVE_HPCT
    hpct_timer_end("plan");
#endif

    // Seed the local memory buffer with processor-dependent, ordered data
    for (int i = 0; i < local_memory; ++i) {
        data.get()[i] = (procid*1e6 + i);
    }

    // Prepare the statistical accumulator for round trip timing
    using namespace boost::accumulators;
    accumulator_set<
            double,
            stats<tag::mean, tag::variance, tag::min, tag::max>
        > roundtrip;

    // Execute the plan the appropriate number of times
    for (int i = 0; i < nrep; ++i) {

        MPI_Barrier(MPI_COMM_WORLD);
        if (!procid) LOG4CXX_DEBUG(logger, "Repetition " << i);

        MPI_Pcontrol(1);                         // Enable MPI profiling
        const double start_trip = MPI_Wtime();   // Start timer
#ifdef HAVE_HPCT
        hpct_timer_begin("long_n2_to_long_n1");
#endif
        plan.execute_long_n2_to_long_n1();
#ifdef HAVE_HPCT
        hpct_timer_end(  "long_n2_to_long_n1");
        hpct_timer_begin("long_n1_to_long_n0");
#endif
        plan.execute_long_n1_to_long_n0();
#ifdef HAVE_HPCT
        hpct_timer_end(  "long_n1_to_long_n0");
        hpct_timer_begin("long_n0_to_long_n1");
#endif
        plan.execute_long_n0_to_long_n1();
#ifdef HAVE_HPCT
        hpct_timer_end(  "long_n0_to_long_n1");
        hpct_timer_begin("long_n1_to_long_n2");
#endif
        plan.execute_long_n1_to_long_n2();
#ifdef HAVE_HPCT
        hpct_timer_end(  "long_n1_to_long_n2");
#endif
        const double end_trip = MPI_Wtime();     // End timer
        MPI_Pcontrol(0);                         // Disable MPI profiling

        // Our round trip time is only as good as the weakest link.
        // In particular, final in-memory transposes may be slower
        // on systems with a poorly balanced problem distribution.
        double sendbuf = end_trip - start_trip;
        double elapsed;
        MPI_Allreduce(&sendbuf, &elapsed, 1,
                      MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        // Accumulate statistics
        roundtrip(elapsed);
    }

#ifdef HAVE_HPCT
    hpct_timer_finalize();
    if (!procid) hpct_timer_summarize();
#endif

    // Dump out the round trip timing statistics
    if (!procid) {
        LOG4CXX_INFO(logger,
                     "Round trip min time:      "
                     << ((min)(roundtrip)) << " seconds");
        LOG4CXX_INFO(logger,
                     "Round trip mean time:     "
                     << ((mean)(roundtrip)) << " seconds");
        LOG4CXX_INFO(logger,
                     "Round trip max time:      "
                     << ((max)(roundtrip)) << " seconds");
        LOG4CXX_INFO(logger,
                     "Round trip stddev:        "
                     << sqrt(variance(roundtrip)) << " seconds");
    }

    // Ensure we recover our seeded processor-dependent, ordered data
    MPI_Barrier(MPI_COMM_WORLD);
    const underling_extents long_n2 = problem.local_extents(2);
    for (int i = 0; i < long_n2.extent; ++i) {
        if (data.get()[i] != (procid*1e6 + i)) {
            LOG4CXX_ERROR(
                    logger,
                    "Did not recover seeded data starting from index " << i);
            break;
        }
    }

    return EXIT_SUCCESS;
}
