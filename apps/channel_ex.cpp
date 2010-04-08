//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2010 The PECOS Development Team
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
// channel_ex.cpp: A fully explicit channel calculation using Suzerain
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <log4cxx/logger.h>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/bspline_definition.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/math.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/orthonormal.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/program_options.hpp>
#include <suzerain/scenario_definition.hpp>

#pragma warning(disable:383 1572)

// Global scenario parameters, initialized by ProgramOptions in main()
static suzerain::problem::ScenarioDefinition<> def_scenario(100);
static suzerain::problem::ChannelDefinition<>  def_grid;
static suzerain::problem::BsplineDefinition<>  def_bspline;

template<typename FPT, typename SizeType>
static void compute_breakpoints(
        FPT xbegin, FPT xend, SizeType n, FPT alpha, FPT *output)
{
    using suzerain::math::stretchspace;
    const FPT xhalf = (xbegin + xend)/2;

    if (n & 1) {
        stretchspace(xbegin, xhalf, n/2+1, alpha,   output);
        stretchspace(xhalf,  xend,  n/2+1, 1/alpha, output+n/2);
    } else {
        const FPT aitch    = (xend-xbegin)/(n + 1);
        const FPT xhalflen = aitch*(n/2);
        stretchspace(xbegin, xbegin+xhalflen, n/2, alpha,   output);
        stretchspace(xend-xhalflen, xend,     n/2, 1/alpha, output+n/2);
    }

    assert(output[0] == xbegin);
    assert(output[n-1] == xend);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);                   // Initialize MPI on startup
    atexit((void (*) ()) MPI_Finalize);       // Finalize MPI at exit

    // Initialize logger using MPI environment details.
    const int nproc  = suzerain::mpi::comm_size(MPI_COMM_WORLD);
    const int procid = suzerain::mpi::comm_rank(MPI_COMM_WORLD);
    log4cxx::LoggerPtr log = log4cxx::Logger::getLogger(
            suzerain::mpi::comm_rank_identifier(MPI_COMM_WORLD));
    // NB: Log only warnings and above from ranks 1 and higher
    if (procid > 0) log->setLevel(log4cxx::Level::getWarn());

    // Process command line arguments
    suzerain::ProgramOptions options;
    options.add_definition(def_scenario);
    options.add_definition(def_grid);
    options.add_definition(def_bspline);
    if (!procid) {
        options.process(argc, argv);
    } else {
        boost::onullstream nullstream;
        options.process(argc, argv,
                        nullstream, nullstream, nullstream, nullstream);
    }
    assert(def_grid.DAFy() == 1.0);  // Wall normal dealiasing disallowed

    // Dump relevant global scenario parameters
    LOG4CXX_INFO(log, "Number of MPI ranks: " << nproc);
    LOG4CXX_INFO(log, "Global extents:      " << def_grid.global_extents());
    LOG4CXX_INFO(log, "Dealiased extents:   " << def_grid.dealiased_extents());
    LOG4CXX_INFO(log, "B-spline order:      " << def_bspline.k());
    LOG4CXX_INFO(log, "B-spline stretching: " << def_bspline.alpha());
    LOG4CXX_INFO(log, "Reynolds number:     " << def_scenario.Re());
    LOG4CXX_INFO(log, "Prandtl number:      " << def_scenario.Pr());
    LOG4CXX_INFO(log, "gamma = C_p/C_v:     " << def_scenario.gamma());
    LOG4CXX_INFO(log, "beta = ln mu/ln T:   " << def_scenario.beta());

    // Initialize B-spline workspace using [0, Ly] and Ny
    double *breakpoints
        = (double *)suzerain::blas::malloc(def_grid.Ny()*sizeof(double));
    assert(breakpoints);
    compute_breakpoints(0.0, def_grid.Ly(), def_grid.Ny(),
                        def_bspline.alpha(), breakpoints);
    for (int i = 0; i < def_grid.Ny(); ++i) {
        LOG4CXX_TRACE(log,
                      "B-spline breakpoint[" << i << "] = " << breakpoints[i]);
    }
    suzerain::bspline bspw(def_bspline.k(), 2, def_grid.Ny(), breakpoints);
    suzerain::blas::free(breakpoints);

    // Initialize pencil_grid which handles P3DFFT setup/teardown RAII
    suzerain::pencil_grid pg(def_grid.dealiased_extents(),
                             def_grid.processor_grid());
    LOG4CXX_INFO(log, "Processor grid used: " << pg.processor_grid());
}
