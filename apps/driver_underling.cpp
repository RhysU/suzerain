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

#include <suzerain/config.h>
#include <suzerain/common.hpp>
#pragma hdrstop

#include <mpi.h>
#include <fftw3-mpi.h>
#include <log4cxx/logger.h>

#ifdef HAVE_HPCT
#include <hpct.h>
#endif

#include <suzerain/error.h>
#include <suzerain/fftw.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/os.h>
#include <suzerain/problem.hpp>
#include <suzerain/program_options.hpp>
#include <suzerain/types.hpp>
#include <suzerain/underling.hpp>

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);                   // Initialize MPI
    atexit((void (*) ()) MPI_Finalize);       // Finalize MPI at exit
    fftw_mpi_init();                          // Initialize FFTW MPI
    atexit((void (*) ()) fftw_mpi_cleanup);   // Finalize FFTW MPI at exit

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
    int nrep    = 1;
    int howmany = 1;
    namespace po = boost::program_options;
    options.add_options()
        ("howmany", po::value<int>(&nrep)->default_value(1),
        "Number of interleaved real-valued fields to transpose")
        ("rep", po::value<int>(&nrep)->default_value(1),
        "Number of repetitions to perform for timing purposes")
    ;
    if (!procid) {
        options.process(argc, argv);
    } else {
        boost::onullstream nullstream;
        options.process(argc, argv,
                        nullstream, nullstream, nullstream, nullstream);
    }

    // Create a grid and a problem
    namespace underling = suzerain::underling;
    underling::grid grid(MPI_COMM_WORLD,
                         griddef.global_extents().begin(),
                         griddef.processor_grid().begin());
    underling::problem problem(grid, howmany);

    // Allocate automatically free, aligned storage
    boost::shared_ptr<underling_real>
        data(static_cast<underling_real *>(
                fftw_malloc(problem.local_memory()*sizeof(underling_real))),
             &fftw_free);

    // Obtain some grid-wide memory usage information
    const size_t local_memory_minimum
        = underling::local_memory_minimum(grid, problem);
    const size_t local_memory_maximum
        = underling::local_memory_maximum(grid, problem);

    // Report grid geometry and memory characteristics
    if (!procid) {
        LOG4CXX_INFO(logger,
                     "Number of processors:  " << nproc);
        LOG4CXX_INFO(logger,
                     "Extents:               " << griddef.global_extents()
                     << " by " << howmany << " scalar field(s)");
        LOG4CXX_INFO(logger,
                     "Processor grid:        " << griddef.processor_grid());
        LOG4CXX_INFO(logger,
                     "Planning rigor:        "
                     << suzerain::fftw::c_str(fftwdef.plan_rigor()));
        LOG4CXX_INFO(logger,
                     "Number of repetitions: " << nrep);
        LOG4CXX_INFO(logger,
                     "Process storage:       "
                     << (local_memory_maximum*sizeof(underling_real))/1024
                     << " KB min, "
                     << (local_memory_minimum*sizeof(underling_real))/1024
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

    // Execute the plan the appropriate number of times
    for (int i = 0; i < nrep; ++i) {

        MPI_Barrier(MPI_COMM_WORLD);
        if (!procid) LOG4CXX_DEBUG(logger, "Repetition " << i);

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
    }

#ifdef HAVE_HPCT
    hpct_timer_finalize();
    if (!procid) hpct_timer_summarize();
#endif

    return EXIT_SUCCESS;
}
