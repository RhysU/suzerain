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

    const ptrdiff_t N0 = grid.nx();  // To be rearranged later
    const ptrdiff_t N1 = grid.ny();
    const ptrdiff_t N2 = grid.nz();
    boost::array<ptrdiff_t,2> P = { grid.pg0(), grid.pg1() };
    suzerain::mpi::dims_create(nproc, P.begin(), P.end());
    std::sort(P.begin(), P.end());

    ONLYPROC0(LOG4CXX_INFO(logger, "Number of processors: " << nproc));
    ONLYPROC0(LOG4CXX_INFO(logger, "Computational grid dimensions: "
                           << boost::format("(% 4d, % 4d, % 4d)")
                           % N0 % N1 % N2));
    ONLYPROC0(LOG4CXX_INFO(logger, "Processor grid dimensions: "
                           << boost::format("(%d, %d)") % P[0] % P[1]));

}
