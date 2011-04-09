//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2008 The PECOS Development Team
//
// Please see http://pecos.ices.utexas.edu for more information.
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
//
// driver_suzerain_p3dfft.cpp: Test driver based on Dmitry Pekurovsky's work
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <suzerain/mpi.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/pencil.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/program_options.hpp>
#include "logger.hpp"

static double real_data(const double x, const double y, const double z)
{
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

    // Establish MPI-savvy, rank-dependent logging names
    name_logger_within_comm_world();

    // Program-specific option storage
    suzerain::ProgramOptions options;
    suzerain::problem::GridDefinition grid(/* default_Nx    */ 16,
                                           /* default_DAFx  */ 3./2.,
                                           /* default_Ny    */ 16,
                                           /* default_k     */ 6,
                                           /* default_Nz    */ 16,
                                           /* default_DAFz  */ 3./2.);
    options.add_definition(grid);

    int nrep = 1;  // Number of times to repeat the test
    namespace po = boost::program_options;
    options.add_options()
        ("nrep", po::value<int>(&nrep)->default_value(nrep),
        "Number of repetitions to perform for timing purposes")
    ;
    options.process(argc, argv);

#pragma warning(push,disable:383)
    const int nproc  = suzerain::mpi::comm_size(MPI_COMM_WORLD);
    INFO0("Number of processors: " << nproc);
    INFO0("Physical grid dimensions: "<< boost::format("(% 4d, % 4d, % 4d)")
          % grid.N.x() % grid.N.y() % grid.N.z());
    INFO0("Processor grid dimensions: " << boost::format("(%d, %d)")
          % grid.P[0] % grid.P[1]);
#pragma warning(pop)

    suzerain::pencil_grid pg(grid.N, grid.P);   // P3DFFT setup/clean RAII
    using suzerain::pencil;
    pencil<> A(pg);                             // Storage management

    // Create a uniform tensor product grid
    std::valarray<double> gridx(A.physical.shape()[0]);
    for (size_t i = 0; i < A.physical.shape()[0]; ++i) {
        gridx[i] = (i+A.physical.index_bases()[0]) * 2*M_PI/grid.N.x();
        TRACE(boost::format("gridx[%3d] = % 6g") % i % gridx[i]);
    }

    std::valarray<double> gridy(A.physical.shape()[1]);
    for (size_t j = 0; j < A.physical.shape()[1]; ++j) {
        gridy[j] = (j+A.physical.index_bases()[1]) * 2*M_PI/grid.N.y();
        TRACE(boost::format("gridy[%3d] = % 6g") % j % gridy[j]);
    }

    std::valarray<double> gridz(A.physical.shape()[2]);
    for (size_t k = 0; k < A.physical.shape()[2]; ++k) {
        gridz[k] = (k+A.physical.index_bases()[2]) * 2*M_PI/grid.N.z();
        TRACE(boost::format("gridz[%3d] = % 6g") % k % gridz[k]);
    }


    INFO("Physical space pencil start and end: "
         << boost::format("[(%3d, %3d, %3d) ... (%3d, %3d, %3d))")
         % A.global_physical.index_bases()[0]
         % A.global_physical.index_bases()[1]
         % A.global_physical.index_bases()[2]
         % (A.global_physical.index_bases()[0] + A.global_physical.shape()[0])
         % (A.global_physical.index_bases()[1] + A.global_physical.shape()[1])
         % (A.global_physical.index_bases()[2] + A.global_physical.shape()[2]));

    INFO("Wave space pencil start and end:     "
         << boost::format("[(%3d, %3d, %3d) ... (%3d, %3d, %3d))")
         % A.global_wave.index_bases()[0]
         % A.global_wave.index_bases()[1]
         % A.global_wave.index_bases()[2]
         % (A.global_wave.index_bases()[0] + A.global_wave.shape()[0])
         % (A.global_wave.index_bases()[1] + A.global_wave.shape()[1])
         % (A.global_wave.index_bases()[2] + A.global_wave.shape()[2]));

    MPI_Barrier(MPI_COMM_WORLD);

    for (pencil<>::size_type j = 0; j < A.physical.shape()[1]; ++j) {
        for (pencil<>::size_type k = 0; k < A.physical.shape()[2]; ++k) {
            for (pencil<>::size_type i = 0; i < A.physical.shape()[0]; ++i) {
                const double value = real_data(gridx[i], gridy[j], gridz[k]);
                TRACE(boost::format(
                      "Physical space (% 6.4f, % 6.4f, % 6.4f) = % 8.4f")
                      % gridx[i] % gridy[j] % gridz[k] % value);
                A.physical[i][j][k] = value;
            }
        }
    }

    // Scale factor in X and Z directions only
    const double factor = 1.0 / (grid.N.x()*grid.N.z());

    double rtime1 = 0.0;

    for (int m = 0; m < nrep; m++) {
        MPI_Barrier(MPI_COMM_WORLD);
        rtime1 = rtime1 - MPI_Wtime();
        INFO0("Iteration " << m);

        p3dfft_ftran_r2c(A.begin(), A.begin()); // Physical to wave
        rtime1 = rtime1 + MPI_Wtime();

        std::transform(A.begin(), A.end(), A.begin(),
                std::bind1st(std::multiplies<pencil<>::real_type>(),factor));

        INFO0("Forward transform results ");

        for (pencil<>::size_type k = A.global_wave.index_bases()[2];
             k < A.global_wave.index_bases()[2] + A.global_wave.shape()[2];
             ++k) {
            for (pencil<>::size_type i = A.global_wave.index_bases()[0];
                i < A.global_wave.index_bases()[0] + A.global_wave.shape()[0];
                ++i) {
                for (pencil<>::size_type j = A.global_wave.index_bases()[1];
                    j < A.global_wave.index_bases()[1] + A.global_wave.shape()[1];
                    ++j) {
                    const pencil<>::complex_type value = A.global_wave[i][j][k];
                    if (abs(value) > 1e-8) {
                        INFO(boost::format("(%3d, %3d, %3d) = (%12g, %12g)")
                             % i % j % k
                             % value.real() % value.imag());
                    }
                }
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        rtime1 = rtime1 - MPI_Wtime();
        p3dfft_btran_c2r(A.begin(), A.begin()); // Wave to physical
        rtime1 = rtime1 + MPI_Wtime();

    }

    MPI_Barrier(MPI_COMM_WORLD);

    /* Check results */
    double cdiff = 0.0;
    for (pencil<>::size_type j = 0; j < A.physical.shape()[1]; ++j) {
        for (pencil<>::size_type k = 0; k < A.physical.shape()[2]; ++k) {
            for (pencil<>::size_type i = 0; i < A.physical.shape()[0]; ++i) {
                const double answer = real_data(gridx[i], gridy[j], gridz[k]);
                const double abserr = fabs(A.physical[i][j][k] - answer);
                cdiff               = std::max(cdiff, abserr);
            }
        }
    }

    // Gather error indicator
    double ccdiff = 0.0;
    MPI_Reduce(&cdiff, &ccdiff, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    INFO0("Maximum difference: " << std::scientific << ccdiff);

    // Gather timing statistics
    double rtime2 = 0.0;
    MPI_Reduce(&rtime1, &rtime2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    INFO0("Time per loop: " << rtime2 / ((double)nrep));
}
