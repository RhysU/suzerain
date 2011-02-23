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
 * channel_common.cpp: Channel-related functionality spanning binaries
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <suzerain/diffwave.hpp>
#include "channel_common.hpp"

using boost::numeric_cast;
using boost::scoped_array;

void store(log4cxx::LoggerPtr log,
           const esio_handle esioh,
           const suzerain::problem::ScenarioDefinition<real_t>& scenario,
           MPI_Comm comm)
{
    // Only root writes data
    const int procid = suzerain::mpi::comm_rank(comm);

    LOG4CXX_INFO(log, "Storing ScenarioDefinition parameters");

    esio_line_establish(esioh, 1, 0, (procid == 0 ? 1 : 0));

    esio_line_write(esioh, "Re", &scenario.Re, 0,
            scenario.options().find("Re",false).description().c_str());

    esio_line_write(esioh, "Pr", &scenario.Pr, 0,
            scenario.options().find("Pr",false).description().c_str());

    esio_line_write(esioh, "gamma", &scenario.gamma, 0,
            scenario.options().find("gamma",false).description().c_str());

    esio_line_write(esioh, "beta", &scenario.beta, 0,
            scenario.options().find("beta",false).description().c_str());

    esio_line_write(esioh, "Lx", &scenario.Lx, 0,
            scenario.options().find("Lx",false).description().c_str());

    esio_line_write(esioh, "Ly", &scenario.Ly, 0,
            scenario.options().find("Ly",false).description().c_str());

    esio_line_write(esioh, "Lz", &scenario.Lz, 0,
            scenario.options().find("Lz",false).description().c_str());
}

void load(log4cxx::LoggerPtr log,
          const esio_handle esioh,
          suzerain::problem::ScenarioDefinition<real_t>& scenario)
{
    LOG4CXX_INFO(log, "Loading ScenarioDefinition parameters");

    esio_line_establish(esioh, 1, 0, 1); // All ranks load

    if (scenario.Re) {
        LOG4CXX_INFO(log, "Overriding file using Re = " << scenario.Re);
    } else {
        esio_line_read(esioh, "Re", &scenario.Re, 0);
    }

    if (scenario.Pr) {
        LOG4CXX_INFO(log, "Overriding file using Pr = " << scenario.Pr);
    } else {
        esio_line_read(esioh, "Pr", &scenario.Pr, 0);
    }

    if (scenario.gamma) {
        LOG4CXX_INFO(log, "Overriding file using gamma = " << scenario.gamma);
    } else {
        esio_line_read(esioh, "gamma", &scenario.gamma, 0);
    }

    if (scenario.beta) {
        LOG4CXX_INFO(log, "Overriding file using beta = " << scenario.beta);
    } else {
        esio_line_read(esioh, "beta", &scenario.beta, 0);
    }

    if (scenario.Lx) {
        LOG4CXX_INFO(log, "Overriding file using Lx = " << scenario.Lx);
    } else {
        esio_line_read(esioh, "Lx", &scenario.Lx, 0);
    }

    if (scenario.Ly) {
        LOG4CXX_INFO(log, "Overriding file using Ly = " << scenario.Ly);
    } else {
        esio_line_read(esioh, "Ly", &scenario.Ly, 0);
    }

    if (scenario.Lz) {
        LOG4CXX_INFO(log, "Overriding file using Lz = " << scenario.Lz);
    } else {
        esio_line_read(esioh, "Lz", &scenario.Lz, 0);
    }
}

void store(log4cxx::LoggerPtr log,
           const esio_handle esioh,
           const suzerain::problem::GridDefinition<real_t>& grid,
           MPI_Comm comm,
           const real_t Lx,
           const real_t Lz)
{
    // Only root process writes
    const int procid = suzerain::mpi::comm_rank(comm);

    LOG4CXX_INFO(log, "Storing GridDefinition parameters");

    esio_line_establish(esioh, 1, 0, (procid == 0 ? 1 : 0));

    const int Nx = numeric_cast<int>(grid.Nx);
    esio_line_write(esioh, "Nx", &Nx, 0,
            grid.options().find("Nx",false).description().c_str());

    esio_line_write(esioh, "DAFx", &grid.DAFx, 0,
            grid.options().find("DAFx",false).description().c_str());

    const int Ny = numeric_cast<int>(grid.Ny);
    esio_line_write(esioh, "Ny", &Ny, 0,
            grid.options().find("Ny",false).description().c_str());

    const int k = numeric_cast<int>(grid.k);
    esio_line_write(esioh, "k", &k, 0,
            grid.options().find("k",false).description().c_str());

    const int Nz = numeric_cast<int>(grid.Nz);
    esio_line_write(esioh, "Nz", &Nz, 0,
            grid.options().find("Nz",false).description().c_str());

    esio_line_write(esioh, "DAFz", &grid.DAFz, 0,
            grid.options().find("DAFz",false).description().c_str());

    LOG4CXX_INFO(log, "Storing wavenumber vectors for Fourier bases");

    const int N  = std::max(grid.Nx, grid.Nz);
    scoped_array<complex_t> buf(new complex_t[N]);

    // Obtain wavenumbers via computing 1*(i*kx)/i
    std::fill_n(buf.get(), N, complex_t(1,0));
    suzerain::diffwave::apply(1, 0, complex_t(0,-1), buf.get(),
            Lx, Lz, 1, Nx, Nx, 0, Nx, 1, 1, 0, 1);
    esio_line_establish(esioh, Nx, 0, (procid == 0 ? Nx : 0));
    esio_line_write(esioh, "kx", reinterpret_cast<real_t *>(buf.get()),
            2, "Wavenumbers in streamwise X direction"); // Re(buf)

    // Obtain wavenumbers via computing 1*(i*kz)/i
    std::fill_n(buf.get(), N, complex_t(1,0));
    suzerain::diffwave::apply(1, 0, complex_t(0,-1), buf.get(),
            Lx, Lz, 1, 1, 1, 0, 1, Nz, Nz, 0, Nz);
    esio_line_establish(esioh, Nz, 0, (procid == 0 ? Nz : 0));
    esio_line_write(esioh, "kz", reinterpret_cast<real_t *>(buf.get()),
            2, "Wavenumbers in spanwise Z direction"); // Re(buf)
}

void load(log4cxx::LoggerPtr log,
          const esio_handle esioh,
          suzerain::problem::GridDefinition<real_t>& grid)
{
    LOG4CXX_INFO(log, "Loading GridDefinition parameters");

    esio_line_establish(esioh, 1, 0, 1); // All ranks load

    if (grid.Nx) {
        LOG4CXX_INFO(log, "Overriding file using Nx = " << grid.Nx);
    } else {
        int Nx;
        esio_line_read(esioh, "Nx", &Nx, 0);
        grid.Nx = Nx;
    }

    if (grid.DAFx) {
        LOG4CXX_INFO(log, "Overriding file using DAFx = " << grid.DAFx);
    } else {
        esio_line_read(esioh, "DAFx", &grid.DAFx, 0);
    }

    if (grid.Ny) {
        LOG4CXX_INFO(log, "Overriding file using Ny = " << grid.Ny);
    } else {
        int Ny;
        esio_line_read(esioh, "Ny", &Ny, 0);
        grid.Ny = Ny;
    }

    if (grid.k) {
        LOG4CXX_INFO(log, "Overriding file using k = " << grid.k);
    } else {
        int k;
        esio_line_read(esioh, "k", &k, 0);
        grid.k = k;
    }

    if (grid.Nz) {
        LOG4CXX_INFO(log, "Overriding file using Nz = " << grid.Nz);
    } else {
        int Nz;
        esio_line_read(esioh, "Nz", &Nz, 0);
        grid.Nz = Nz;
    }

    if (grid.DAFz) {
        LOG4CXX_INFO(log, "Overriding file using DAFz = " << grid.DAFz);
    } else {
        esio_line_read(esioh, "DAFz", &grid.DAFz, 0);
    }
}

void store(log4cxx::LoggerPtr log,
           const esio_handle esioh,
           boost::shared_ptr<suzerain::bspline>& bspw, // Yes, a reference
           MPI_Comm comm)
{
    // Only root process writes data
    const int procid = suzerain::mpi::comm_rank(comm);

    LOG4CXX_INFO(log, "Storing B-spline knot details");

    scoped_array<real_t> buf(new real_t[bspw->nknots()]);

    bspw->knots(buf.get(), 1);
    esio_line_establish(esioh, bspw->nknots(),
            0, (procid == 0 ? bspw->nknots() : 0));
    esio_line_write(esioh, "knots", buf.get(), 0,
            "Knots used to build B-spline basis");

    bspw->breakpoints(buf.get(), 1);
    esio_line_establish(esioh, bspw->ndof(),
            0, (procid == 0 ? bspw->ndof() : 0));
    esio_line_write(esioh, "breakpoints", buf.get(), 0,
            "Breakpoint locations used to build B-spline basis");

    bspw->collocation_points(buf.get(), 1);
    esio_line_establish(esioh, bspw->ndof(),
            0, (procid == 0 ? bspw->ndof() : 0));
    esio_line_write(esioh, "collocation_points", buf.get(), 0,
            "Collocation points used to build discrete operators");

    LOG4CXX_INFO(log, "Storing B-spline derivative operators");

    char name[8];
    char comment[127];
    for (int k = 0; k <= bspw->nderivatives(); ++k) {
        snprintf(name, sizeof(name)/sizeof(name[0]), "Dy%d", k);
        snprintf(comment, sizeof(comment)/sizeof(comment[0]),
                "Wall-normal derivative Dy%d(i,j) = D%d[j,ku+i-j] for"
                " 0 <= j < n, max(0,j-ku-1) <= i < min(m,j+kl)", k, k);
        const int lda = bspw->ku(k) + 1 + bspw->kl(k);
        esio_plane_establish(esioh,
                bspw->ndof(), 0, (procid == 0 ? bspw->ndof() : 0),
                lda,          0, (procid == 0 ? lda          : 0));
        esio_plane_write(esioh, name, bspw->D(k), 0, 0, comment);
        esio_attribute_write(esioh, name, "kl", bspw->kl(k));
        esio_attribute_write(esioh, name, "ku", bspw->ku(k));
        esio_attribute_write(esioh, name, "m",  bspw->ndof());
        esio_attribute_write(esioh, name, "n",  bspw->ndof());
    }
}

void load(log4cxx::LoggerPtr log,
          const esio_handle esioh,
          boost::shared_ptr<suzerain::bspline>& bspw, // Yes, a reference
          const suzerain::problem::GridDefinition<real_t>& grid)
{
    LOG4CXX_INFO(log, "Loading B-spline breakpoints");

    const int nbreak = grid.Ny + 2 - grid.k;
    scoped_array<real_t> buf(new double[grid.Ny]);
    esio_line_establish(esioh, nbreak, 0, nbreak); // All ranks load
    esio_line_read(esioh, "breakpoints", buf.get(), 0);
    bspw.reset(new suzerain::bspline(grid.k, grid.k - 2, nbreak, buf.get()));
}
