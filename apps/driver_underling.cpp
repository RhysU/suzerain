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
    suzerain::problem::GridDefinition<> griddef;
    options.add_definition(griddef);
    namespace po = boost::program_options;
    if (!procid) {
        options.process(argc, argv);
    } else {
        boost::onullstream nullstream;
        options.process(argc, argv,
                        nullstream, nullstream, nullstream, nullstream);
    }


    /* Create a grid and a problem*/
    underling_grid grid = underling_grid_create(
            MPI_COMM_WORLD, griddef.nx(), griddef.ny(), griddef.nz(), 0, 0);
    underling_problem problem = underling_problem_create(grid, 1);

    /* Dump some runtime information */
    ONLYPROC0(
        LOG4CXX_INFO(logger, "Number of processors: " << nproc);
        LOG4CXX_INFO(logger, "grid->np0:            " << grid->np0);
        LOG4CXX_INFO(logger, "grid->nw0:            " << grid->nw0);
        LOG4CXX_INFO(logger, "grid->n1:             " << grid->n1);
        LOG4CXX_INFO(logger, "grid->n2:             " << grid->n2);
        LOG4CXX_INFO(logger, "grid->p0:             " << grid->p0);
        LOG4CXX_INFO(logger, "grid->p1:             " << grid->p1);
        LOG4CXX_INFO(logger, "problem->nfields:     " << problem->nfields);
    );
    LOG4CXX_INFO(logger, "grid->p0_comm:        " << grid->p0_comm);
    LOG4CXX_INFO(logger, "grid->p1_comm:        " << grid->p1_comm);
    LOG4CXX_INFO(logger, "problem->local_size:  " << problem->local_size);
    LOG4CXX_INFO(logger, "problem->optimum_size:" << underling_optimum_local_size(problem));

    /* Allocate storage and create a plan */
    const size_t local_size = underling_local_size(problem);
    underling_real * const data
        = (underling_real *) fftw_malloc(local_size*sizeof(underling_real));
    underling_plan plan = underling_plan_create(problem, data, 1, 1, 0);

    MPI_Barrier(MPI_COMM_WORLD);

    /* Initialize test data in wave space */
    for (int i = 0; i < local_size; ++i) {
        data[i] = procid*10000 + i;
        LOG4CXX_DEBUG(logger, "initial data["
                << std::setw(8) << std::setfill('0') << i << "] = "
                << std::setw(8) << std::setfill(' ') << data[i]);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /* Transform to physical space */
    LOG4CXX_DEBUG(logger, "underling_execute_c2r");
    underling_execute_c2r(plan);
    for (int i = 0; i < local_size; ++i) {
        LOG4CXX_DEBUG(logger, "post c2r data["
                << std::setw(8) << std::setfill('0') << i << "] = "
                << std::setw(8) << std::setfill(' ') << data[i]);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /* Transform to wave space */
    LOG4CXX_DEBUG(logger, "underling_execute_r2c");
    underling_execute_r2c(plan);
    for (int i = 0; i < local_size; ++i) {
        LOG4CXX_DEBUG(logger, "post r2c data["
                << std::setw(8) << std::setfill('0') << i << "] = "
                << std::setw(8) << std::setfill(' ') << data[i]);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /* Primitive check for data corruption */
    int corruption = 0;
    for (int i = 0; i < local_size; ++i) {
        if (data[i] != (procid*10000 + i)) {
            LOG4CXX_WARN(logger, "test result discrepancy at index " << i);
            corruption = 1;
            break;
        }
    }
    int retval;
    SUZERAIN_MPICHKQ(MPI_Allreduce(
            &corruption, &retval, 1, MPI_INT, MPI_BOR, MPI_COMM_WORLD));

    /* Clean up after ourselves */
    underling_plan_destroy(plan);
    fftw_free(data);
    underling_problem_destroy(problem);
    underling_grid_destroy(grid);

    return retval;
}
