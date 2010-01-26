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

#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>

#include <suzerain/mpi.hpp>
#include <suzerain/os.h>
#include <suzerain/problem.hpp>
#include <suzerain/program_options.hpp>
#include <suzerain/types.hpp>
#include <suzerain/underling.h>

#define ONLYPROC0(expr) if (!procid) { expr ; } else

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);                   // Initialize MPI
    atexit((void (*) ()) MPI_Finalize);       // Finalize MPI at exit
    fftw_mpi_init();                          // Initialize FFTW MPI
    atexit((void (*) ()) fftw_mpi_cleanup);   // Finalize FFTW MPI at exit

    const int nproc  = suzerain::mpi::comm_size(MPI_COMM_WORLD);
    const int procid = suzerain::mpi::comm_rank(MPI_COMM_WORLD);
    log4cxx::LoggerPtr logger = log4cxx::Logger::getLogger(
            suzerain::mpi::comm_rank_identifier(MPI_COMM_WORLD));

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

    { // Dump grid and problem information
        // TODO: Error checking on these pipe and FILE* operations
        FILE *pipewrite;
        int pipereadfd;
        suzerain_fpipe(&pipereadfd, O_NONBLOCK, &pipewrite, O_NONBLOCK);
        namespace io = boost::iostreams;
        io::stream<io::file_descriptor_source> piperead(pipereadfd);
        ONLYPROC0(
            underling_fprint_grid(grid, pipewrite); fflush(pipewrite);
            LOG4CXX_INFO(logger, piperead.rdbuf());
        );
        underling_fprint_problem(problem, pipewrite); fflush(pipewrite);
        LOG4CXX_DEBUG(logger, piperead.rdbuf());
        fclose(pipewrite);
        close(piperead);
    }

    /* Allocate storage and create a plan */
    const size_t local_size = underling_local_size(problem);
    LOG4CXX_DEBUG(logger, "problem local_size = " << local_size);
    LOG4CXX_DEBUG(logger, "problem optimum_local_size = "
                         << underling_optimum_local_size(problem));
    underling_real * const data
        = (underling_real *) fftw_malloc(local_size*sizeof(underling_real));
    underling_plan plan = underling_plan_create(
            problem, data, UNDERLING_DIRECTION_BOTH, 0);

    MPI_Barrier(MPI_COMM_WORLD);

    /* Initialize test data in wave space */
    for (int i = 0; i < local_size; ++i) {
        data[i] = procid*10000 + i;
        LOG4CXX_TRACE(logger, "initial data["
                << std::setw(8) << std::setfill('0') << i << "] = "
                << std::setw(8) << std::setfill(' ') << data[i]);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /* Transform to physical space */
    LOG4CXX_DEBUG(logger, "underling_execute_backward");
    underling_execute_backward(plan);
    for (int i = 0; i < local_size; ++i) {
        LOG4CXX_TRACE(logger, "post backward data["
                << std::setw(8) << std::setfill('0') << i << "] = "
                << std::setw(8) << std::setfill(' ') << data[i]);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /* Transform to wave space */
    LOG4CXX_DEBUG(logger, "underling_execute_forward");
    underling_execute_forward(plan);
    for (int i = 0; i < local_size; ++i) {
        LOG4CXX_TRACE(logger, "post forward data["
                << std::setw(8) << std::setfill('0') << i << "] = "
                << std::setw(8) << std::setfill(' ') << data[i]);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /* Primitive check for data corruption */
    int corruption = 0;
    const size_t long_n2_data
        = underling_local_long_n2(problem, NULL, NULL, NULL);
    for (int i = 0; i < long_n2_data; ++i) {
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
