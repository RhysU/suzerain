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

#include <suzerain/mpi.hpp>
#include <suzerain/types.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/program_options.hpp>
#include <suzerain/underling.h>

#define ONLYPROC0(expr) if (!procid) { expr ; } else

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);                   // Initialize MPI
    atexit((void (*) ()) MPI_Finalize);       // Finalize MPI at exit
    fftw_mpi_init();                          // Initialize FFTW MPI
    atexit((void (*) ()) fftw_mpi_cleanup);   // Finalize FFTW MPI

    const int nproc  = suzerain::mpi::comm_size(MPI_COMM_WORLD);
    const int procid = suzerain::mpi::comm_rank(MPI_COMM_WORLD);

    // Initialize logger with processor number
    std::ostringstream procname;
    procname << "proc"
             << std::setfill('0') << std::setw(ceil(log10(nproc))) << procid;
    log4cxx::LoggerPtr logger = log4cxx::Logger::getLogger(procname.str());

    suzerain::ProgramOptions options;
    suzerain::problem::GridDefinition<> grid;
    options.add_definition(grid);
    namespace po = boost::program_options;
    if (!procid) {
        options.process(argc, argv);
    } else {
        boost::onullstream nullstream;
        options.process(argc, argv,
                        nullstream, nullstream, nullstream, nullstream);
    }

    underling_problem *problem = underling_problem_create(
            MPI_COMM_WORLD, grid.nx(), grid.ny(), grid.nz(), 0, 0);
    assert(problem);

    ONLYPROC0(LOG4CXX_INFO(logger, "Number of processors: " << nproc));

    LOG4CXX_INFO(logger, "problem->np0:     " << problem->np0);
    LOG4CXX_INFO(logger, "problem->nw0:     " << problem->nw0);
    LOG4CXX_INFO(logger, "problem->n1:      " << problem->n1);
    LOG4CXX_INFO(logger, "problem->n2:      " << problem->n2);
    LOG4CXX_INFO(logger, "problem->p0:      " << problem->p0);
    LOG4CXX_INFO(logger, "problem->p1:      " << problem->p1);
    LOG4CXX_INFO(logger, "problem->g_comm:  " << problem->g_comm);
    LOG4CXX_INFO(logger, "problem->p0_comm: " << problem->p0_comm);
    LOG4CXX_INFO(logger, "problem->p1_comm: " << problem->p1_comm);
    LOG4CXX_INFO(logger, "problem->block_a: " << problem->block_a);
    LOG4CXX_INFO(logger, "problem->block_b: " << problem->block_b);
    LOG4CXX_INFO(logger, "problem->block_c: " << problem->block_c);

    underling_problem_destroy(problem);
}
