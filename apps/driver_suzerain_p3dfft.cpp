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
// driver_suzerain_p3dfft.cpp: Test driver based on Dmitry Pekurovsky's work
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <fftw3.h>
#include <p3dfft_d.h>
#include <suzerain/fftw_definition.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/pencil.hpp>
#include <suzerain/pre_gsl.h>
#include <suzerain/program_options.hpp>
#include "logging.hpp"

// Provided by driver_suzerain_p3dfft_svnrev.{c,h} to speed recompilation
#pragma warning(push,disable:1419)
extern "C" const char revstr[];
#pragma warning(pop)

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
    logging::initialize(MPI_COMM_WORLD);      // Initialize logging

    DEBUG0("Establishing floating point environment from GSL_IEEE_MODE");
    mpi_gsl_ieee_env_setup(suzerain::mpi::comm_rank(MPI_COMM_WORLD));

    // Process command line options
    suzerain::program_options options(
            "suzerain::pencil_grid_p3dfft performance benchmark",
            "", /* TODO description */ "", revstr);
    suzerain::grid_definition grid(/* Lx UNUSED */ "NaN",
                                  /* Nx        */ 16,
                                  /* DAFx      */ 1.,
                                  /* Ly UNUSED*/ "NaN",
                                  /* Ny        */ 16,
                                  /* k         */ 6,
                                  /* htdelta   */ 0,
                                  /* Lz UNUSED*/ "NaN",
                                  /* Nz        */ 16,
                                  /* DAFz      */ 1.);
    options.add_definition(grid);
    suzerain::fftw_definition fftwdef(
            suzerain::fftw::measure, suzerain::fftw::estimate);
    options.add_definition(fftwdef);

    int  repeat  = 1;
    int  nfields = 1;
    bool inplace = false;
    bool check   = false;
    namespace po = boost::program_options;
    options.add_options()
        ("repeat,r",
         po::value(&repeat)->default_value(repeat),
         "Number of repetitions to perform")
        ("nfields,n",
         po::value(&nfields)->default_value(nfields),
         "Number of independent fields")
        ("in-place,i",
         po::value(&inplace)->default_value(inplace)->zero_tokens(),
         "Perform in-place transposes")
        ("check",
         po::value(&check)->default_value(check)->zero_tokens(),
         "Check results against expected values")
    ;
    options.process(argc, argv);
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
    if (repeat  < 1) throw std::invalid_argument("repeat  < 1");
    if (nfields < 1) throw std::invalid_argument("nfields < 1");

    INFO0("Preparing MPI transpose and Fourier transform execution plans...");
    const double wtime_fftw_planning_start = MPI_Wtime();
    fftw_set_timelimit(fftwdef.plan_timelimit);
    suzerain::pencil_grid_p3dfft pg(
        grid.N, grid.P, fftwdef.rigor_fft, fftwdef.rigor_mpi);
    const double wtime_fftw_planning = MPI_Wtime() - wtime_fftw_planning_start;
    INFO0("MPI transpose and Fourier transform planning took "
          << wtime_fftw_planning << " seconds");

#pragma warning(push,disable:383)
    const int nproc  = suzerain::mpi::comm_size(MPI_COMM_WORLD);

    INFO0("Global physical extents:   " << pg.global_physical_extent.transpose());
    INFO0("Global wave extents:       " << pg.global_wave_extent.transpose());
    INFO0("Processor grid dimensions: " << pg.processor_grid.transpose());
    INFO0("Number of processors:      " << nproc);
    INFO0("Number of fields:          " << nfields);
    INFO0("Performing operations " << (inplace ? "in-place" : "out-of-place"));
#pragma warning(pop)

    // Allocate necessary storage
    using suzerain::pencil;
    const int ploff = inplace ? 0 : 1;
    boost::ptr_vector<pencil<> > pencils(nfields + ploff);
    for (int l = 0; l < nfields + ploff; ++l) {
        pencils.push_back(new pencil<>(pg));
    }

    // Get decomposition and track correctness using first pencil only
    pencil<>& A = pencils[0        ]; // Physical space
    pencil<>& B = pencils[0 + ploff]; // Wave space

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

    DEBUG("Physical space pencil start and end: "
         << boost::format("[(%3d, %3d, %3d) ... (%3d, %3d, %3d))")
         % A.global_physical.index_bases()[0]
         % A.global_physical.index_bases()[1]
         % A.global_physical.index_bases()[2]
         % (A.global_physical.index_bases()[0] + A.global_physical.shape()[0])
         % (A.global_physical.index_bases()[1] + A.global_physical.shape()[1])
         % (A.global_physical.index_bases()[2] + A.global_physical.shape()[2]));

    DEBUG("Wave space pencil start and end:     "
         << boost::format("[(%3d, %3d, %3d) ... (%3d, %3d, %3d))")
         % A.global_wave.index_bases()[0]
         % A.global_wave.index_bases()[1]
         % A.global_wave.index_bases()[2]
         % (A.global_wave.index_bases()[0] + A.global_wave.shape()[0])
         % (A.global_wave.index_bases()[1] + A.global_wave.shape()[1])
         % (A.global_wave.index_bases()[2] + A.global_wave.shape()[2]));

    MPI_Barrier(MPI_COMM_WORLD);

    // Initialize data in first pencil
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
    // ...and copy detail to all other pencils
    for (int l = 1; l < nfields; ++l) pencils[l].physical = A.physical;

    // Scale factor in X and Z directions only
    const double factor = 1.0 / (grid.N.x()*grid.N.z());

    double rtime1 = 0.0;

    for (int m = 0; m < repeat; m++) {
        INFO0("Iteration " << (m + 1));

        // Physical to wave
        MPI_Barrier(MPI_COMM_WORLD);
        rtime1 -= MPI_Wtime();
        for (int l = nfields; l --> 0 ;) {
            p3dfft_ftran_r2c(pencils[l].begin(), pencils[l + ploff].begin());
        }
        rtime1 += MPI_Wtime();

        // Normalization
        for (int l = nfields; l --> 0 ;) {
            std::transform(
                    pencils[l + ploff].begin(), pencils[l + ploff].end(),
                    pencils[l + ploff].begin(),
                    std::bind1st(std::multiplies<pencil<>::real_type>(),factor));
        }

        if (m == repeat - 1 && DEBUG_ENABLED() && check) {
            DEBUG0("Forward transform results ");
            for (pencil<>::size_type k = B.global_wave.index_bases()[2];
                k < B.global_wave.index_bases()[2] + B.global_wave.shape()[2];
                ++k) {
                for (pencil<>::size_type i = B.global_wave.index_bases()[0];
                    i < B.global_wave.index_bases()[0] + B.global_wave.shape()[0];
                    ++i) {
                    for (pencil<>::size_type j = B.global_wave.index_bases()[1];
                        j < B.global_wave.index_bases()[1] + B.global_wave.shape()[1];
                        ++j) {
                        const pencil<>::complex_type value = B.global_wave[i][j][k];
                        if (abs(value) > 1e-8) {
                            DEBUG(boost::format("(%3d, %3d, %3d) = (%12g, %12g)")
                                % i % j % k
                                % value.real() % value.imag());
                        }
                    }
                }
            }
        }

        // Wave to physical
        MPI_Barrier(MPI_COMM_WORLD);
        rtime1 -= MPI_Wtime();
        for (int l = 0; l < nfields; ++l) {
            p3dfft_btran_c2r(pencils[l + ploff].begin(), pencils[l].begin());
        }
        rtime1 += MPI_Wtime();

    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (check) {
        /* Check results */
        INFO0("Checking results against expected values");
        double cdiff = 0.0;
        for (pencil<>::size_type j = 0; j < A.physical.shape()[1]; ++j) {
            for (pencil<>::size_type k = 0; k < A.physical.shape()[2]; ++k) {
                for (pencil<>::size_type i = 0; i < A.physical.shape()[0]; ++i) {
                    const double answer = real_data(gridx[i], gridy[j], gridz[k]);
                    const double abserr = std::abs(A.physical[i][j][k] - answer);
                    cdiff               = std::max(cdiff, abserr);
                }
            }
        }

        // Gather error indicator
        double ccdiff = 0.0;
        MPI_Reduce(&cdiff, &ccdiff, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        INFO0("Maximum difference: " << std::scientific << ccdiff);
    }

    // Gather timing statistics
    double rtime2 = 0.0;
    MPI_Reduce(&rtime1, &rtime2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    INFO0("Seconds per iteration: " << rtime2 / ((double) repeat));
}
