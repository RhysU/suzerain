/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
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
 * driver_p3dfft.cc: A P3DFFT test driver based on work by Dmitry Pekurovsky
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <suzerain/config.h>
#include <suzerain/common.hpp>
#pragma hdrstop

#include <log4cxx/logger.h>
#include <mpi.h>

#include <suzerain/mpi.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/pencil.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/program_options.hpp>

#define ONLYPROC0(expr) if (!procid) { expr ; } else

double real_data(const double x, const double y, const double z) {
    return   1.0* 1.0
           +      2.0*y
           + 2.0* 3.0*sin(x)
           + 2.0* 5.0*sin(z)
           + 2.0* 7.0*sin(x)*y
           + 2.0*11.0*y*sin(z)
           + 4.0*13.0*sin(x)*sin(z)
           + 4.0*17.0*sin(x)*y*sin(z)
           ;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);                   // Initialize MPI on startup
    atexit((void (*) ()) MPI_Finalize);       // Finalize down MPI at exit

    const int nproc  = suzerain::mpi::comm_size(MPI_COMM_WORLD);
    const int procid = suzerain::mpi::comm_rank(MPI_COMM_WORLD);

    // Initialize logger with processor number
    std::ostringstream procname;
    procname << "proc"
             << std::setfill('0') << std::setw(ceil(log10(nproc))) << procid;
    log4cxx::LoggerPtr logger = log4cxx::Logger::getLogger(procname.str());

    // Program-specific option storage
    int nrep;  // Number of times to repeat the test

    suzerain::ProgramOptions options;
    suzerain::problem::GridDefinition<> grid;
    options.add_definition(grid);
    namespace po = boost::program_options;
    options.add_options()
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

    ONLYPROC0(LOG4CXX_INFO(logger, "Number of processors: " << nproc));

    ONLYPROC0(LOG4CXX_INFO(logger, "Physical grid dimensions: "
                           << boost::format("(% 4d, % 4d, % 4d)")
                           % grid.nx() % grid.ny() % grid.nz()));

    ONLYPROC0(LOG4CXX_INFO(logger, "Processor grid dimensions: "
                           << boost::format("(%d, %d)")
                           % grid.pg0() % grid.pg1()));


    // pencil_grid handles P3DFFT setup/clean RAII
    using suzerain::pencil_grid;
    pencil_grid pg(grid.global_extents(), grid.processor_grid());
    // pencil handles memory allocation and storage layout
    using suzerain::pencil;
    pencil<> A(pg);

    // Create a uniform tensor product grid
    std::valarray<double> gridx(A.physical.size_x);
    for (size_t i = 0; i < A.physical.size_x; ++i) {
        gridx[i] = (i+A.physical.start_x) * 2*M_PI/grid.nx();
        LOG4CXX_TRACE(logger, boost::format("gridx[%3d] = % 6g") % i % gridx[i]);
    }

    std::valarray<double> gridy(A.physical.size_y);
    for (size_t j = 0; j < A.physical.size_y; ++j) {
        gridy[j] = (j+A.physical.start_y) * 2*M_PI/grid.ny();
        LOG4CXX_TRACE(logger, boost::format("gridy[%3d] = % 6g") % j % gridy[j]);
    }

    std::valarray<double> gridz(A.physical.size_z);
    for (size_t k = 0; k < A.physical.size_z; ++k) {
        gridz[k] = (k+A.physical.start_z) * 2*M_PI/grid.nz();
        LOG4CXX_TRACE(logger, boost::format("gridz[%3d] = % 6g") % k % gridz[k]);
    }


    LOG4CXX_INFO(logger,
                 "Physical space pencil start and end: "
                 << boost::format("[(%3d, %3d, %3d) ... (%3d, %3d, %3d))")
                 % A.physical.start_x
                 % A.physical.start_y
                 % A.physical.start_z
                 % A.physical.end_x
                 % A.physical.end_y
                 % A.physical.end_z);

    LOG4CXX_INFO(logger,
                 "Wave space pencil start and end:     "
                 << boost::format("[(%3d, %3d, %3d) ... (%3d, %3d, %3d))")
                 % A.wave.start_x
                 % A.wave.start_y
                 % A.wave.start_z
                 % A.wave.end_x
                 % A.wave.end_y
                 % A.wave.end_z);

    MPI_Barrier(MPI_COMM_WORLD);

    for (pencil<>::size_type j = 0; j < A.physical.size_y; ++j) {
        for (pencil<>::size_type k = 0; k < A.physical.size_z; ++k) {
            for (pencil<>::size_type i = 0; i < A.physical.size_x; ++i) {
                const double value = real_data(gridx[i], gridy[j], gridz[k]);
                LOG4CXX_TRACE(logger, boost::format(
                              "Physical space (% 6.4f, % 6.4f, % 6.4f) = % 8.4f")
                              % gridx[i] % gridy[j] % gridz[k] % value);
                A.physical(i,j,k) = value;
            }
        }
    }

    // Scale factor in X and Z directions only
    const double factor = 1.0 / (grid.nx()*grid.nz());

    double rtime1 = 0.0;

    for (int m = 0; m < nrep; m++) {
        MPI_Barrier(MPI_COMM_WORLD);
        rtime1 = rtime1 - MPI_Wtime();
        ONLYPROC0(LOG4CXX_INFO(logger, "Iteration " << m));

        p3dfft_ftran_r2c(A.data(), A.data()); // Physical to wave
        rtime1 = rtime1 + MPI_Wtime();

        std::transform(A.begin(), A.end(), A.begin(),
                       std::bind1st(std::multiplies<pencil<>::real_type>(),factor));

        ONLYPROC0(LOG4CXX_INFO(logger, "Forward transform results "));

        for (pencil<>::complex_iterator it = A.wave.begin();
             it != A.wave.end();
             ++it) {
            if (abs(*it) > 1e-8) {
                pencil<>::index i, j, k;
                A.wave.inverse_global_offset(it - A.wave.begin(), i, j, k);
                LOG4CXX_INFO(logger,
                        boost::format("(%3d, %3d, %3d) = (%12g, %12g) at index %3d")
                        % i % j % k
                        % it->real() % it->imag()
                        % (it-A.wave.begin()) );
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        rtime1 = rtime1 - MPI_Wtime();
        p3dfft_btran_c2r(A.data(), A.data()); // Wave to physical
        rtime1 = rtime1 + MPI_Wtime();

    }

    MPI_Barrier(MPI_COMM_WORLD);

    /* Check results */
    double cdiff = 0.0;
    for (pencil<>::size_type j = 0; j < A.physical.size_y; ++j) {
        for (pencil<>::size_type k = 0; k < A.physical.size_z; ++k) {
            for (pencil<>::size_type i = 0; i < A.physical.size_x; ++i) {
                const double answer = real_data(gridx[i], gridy[j], gridz[k]);
                const double abserr = fabs(A.physical(i,j,k) - answer);
                cdiff               = std::max(cdiff, abserr);
            }
        }
    }

    // Gather error indicator
    double ccdiff = 0.0;
    MPI_Reduce(&cdiff, &ccdiff, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    ONLYPROC0(LOG4CXX_INFO(logger,
                           "Maximum difference: " << std::scientific << ccdiff));

    // Gather timing statistics
    double rtime2 = 0.0;
    MPI_Reduce(&rtime1, &rtime2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    ONLYPROC0(LOG4CXX_INFO(logger, "Time per loop: " << rtime2 / ((double)nrep)));
}
