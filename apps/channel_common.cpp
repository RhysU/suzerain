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
#include "logger.hpp"
#include "channel_common.hpp"

using boost::numeric_cast;
using boost::scoped_array;

const boost::array<const char *,5> field_names = {{
    "rho", "rhou", "rhov", "rhow", "rhoe"
}};

const boost::array<const char *,5> field_descriptions = {{
    "Nondimensional density coefficients stored row-major ZXY",
    "Nondimensional X momentum coefficients stored row-major ZXY",
    "Nondimensional Y momentum coefficients stored row-major ZXY",
    "Nondimensional Z momentum coefficients stored row-major ZXY",
    "Nondimensional total energy coefficients stored row-major ZXY"
}};

void store(const esio_handle esioh,
           const suzerain::problem::ScenarioDefinition<real_t>& scenario)
{
    DEBUG0("Storing ScenarioDefinition parameters");

    // Only root writes data
    int procid;
    esio_handle_comm_rank(esioh, &procid);

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

void load(const esio_handle esioh,
          suzerain::problem::ScenarioDefinition<real_t>& scenario)
{
    DEBUG0("Loading ScenarioDefinition parameters");

    esio_line_establish(esioh, 1, 0, 1); // All ranks load

    if (scenario.Re) {
        INFO0("Overriding scenario using Re = " << scenario.Re);
    } else {
        esio_line_read(esioh, "Re", &scenario.Re, 0);
    }

    if (scenario.Pr) {
        INFO0("Overriding scenario using Pr = " << scenario.Pr);
    } else {
        esio_line_read(esioh, "Pr", &scenario.Pr, 0);
    }

    if (scenario.gamma) {
        INFO0("Overriding scenario using gamma = " << scenario.gamma);
    } else {
        esio_line_read(esioh, "gamma", &scenario.gamma, 0);
    }

    if (scenario.beta) {
        INFO0("Overriding scenario using beta = " << scenario.beta);
    } else {
        esio_line_read(esioh, "beta", &scenario.beta, 0);
    }

    if (scenario.Lx) {
        INFO0("Overriding scenario using Lx = " << scenario.Lx);
    } else {
        esio_line_read(esioh, "Lx", &scenario.Lx, 0);
    }

    if (scenario.Ly) {
        INFO0("Overriding scenario using Ly = " << scenario.Ly);
    } else {
        esio_line_read(esioh, "Ly", &scenario.Ly, 0);
    }

    if (scenario.Lz) {
        INFO0("Overriding scenario using Lz = " << scenario.Lz);
    } else {
        esio_line_read(esioh, "Lz", &scenario.Lz, 0);
    }
}

void store(const esio_handle esioh,
           const suzerain::problem::GridDefinition& grid,
           const real_t Lx,
           const real_t Lz)
{
    // Only root writes data
    int procid;
    esio_handle_comm_rank(esioh, &procid);

    DEBUG0("Storing GridDefinition parameters");

    esio_line_establish(esioh, 1, 0, (procid == 0 ? 1 : 0));

    esio_line_write(esioh, "Nx", &grid.N.x(), 0,
            grid.options().find("Nx",false).description().c_str());

    const double DAFx = grid.DAFx();
    esio_line_write(esioh, "DAFx", &DAFx, 0,
            grid.options().find("DAFx",false).description().c_str());

    esio_line_write(esioh, "Ny", &grid.N.y(), 0,
            grid.options().find("Ny",false).description().c_str());

    esio_line_write(esioh, "k", &grid.k, 0,
            grid.options().find("k",false).description().c_str());

    esio_line_write(esioh, "Nz", &grid.N.z(), 0,
            grid.options().find("Nz",false).description().c_str());

    const double DAFz = grid.DAFz();
    esio_line_write(esioh, "DAFz", &DAFz, 0,
            grid.options().find("DAFz",false).description().c_str());

    DEBUG0("Storing wavenumber vectors for Fourier bases");

    const int N = std::max(grid.N.x(), grid.N.z());
    scoped_array<complex_t> buf(new complex_t[N]);

    // Obtain wavenumbers via computing 1*(i*kx)/i
    std::fill_n(buf.get(), N, complex_t(1,0));
    suzerain::diffwave::apply(1, 0, complex_t(0,-1), buf.get(),
            Lx, Lz, 1, grid.N.x(), grid.N.x(), 0, grid.N.x(), 1, 1, 0, 1);
    esio_line_establish(esioh, grid.N.x(), 0, (procid == 0 ? grid.N.x() : 0));
    esio_line_write(esioh, "kx", reinterpret_cast<real_t *>(buf.get()),
            2, "Wavenumbers in streamwise X direction"); // Re(buf)

    // Obtain wavenumbers via computing 1*(i*kz)/i
    std::fill_n(buf.get(), N, complex_t(1,0));
    suzerain::diffwave::apply(1, 0, complex_t(0,-1), buf.get(),
            Lx, Lz, 1, 1, 1, 0, 1, grid.N.z(), grid.N.z(), 0, grid.N.z());
    esio_line_establish(esioh, grid.N.z(), 0, (procid == 0 ? grid.N.z() : 0));
    esio_line_write(esioh, "kz", reinterpret_cast<real_t *>(buf.get()),
            2, "Wavenumbers in spanwise Z direction"); // Re(buf)
}

void load(const esio_handle esioh,
          suzerain::problem::GridDefinition& grid)
{
    DEBUG0("Loading GridDefinition parameters");

    esio_line_establish(esioh, 1, 0, 1); // All ranks load

    if (grid.N.x()) {
        INFO0("Overriding grid using Nx = " << grid.N.x());
    } else {
        esio_line_read(esioh, "Nx", &grid.N.x(), 0);
    }

    if (grid.DAFx()) {
        INFO0("Overriding grid using DAFx = " << grid.DAFx());
    } else {
        double factor;
        esio_line_read(esioh, "DAFx", &factor, 0);
        grid.DAFx(factor);
    }

    if (grid.N.y()) {
        INFO0("Overriding grid using Ny = " << grid.N.y());
    } else {
        esio_line_read(esioh, "Ny", &grid.N.y(), 0);
    }

    if (grid.k) {
        INFO0("Overriding grid using k = " << grid.k);
    } else {
        esio_line_read(esioh, "k", &grid.k, 0);
    }

    if (grid.N.z()) {
        INFO0("Overriding grid using Nz = " << grid.N.z());
    } else {
        esio_line_read(esioh, "Nz", &grid.N.z(), 0);
    }

    if (grid.DAFz()) {
        INFO0("Overriding grid using DAFz = " << grid.DAFz());
    } else {
        double factor;
        esio_line_read(esioh, "DAFz", &factor, 0);
        grid.DAFz(factor);
    }
}

void store(const esio_handle esioh,
           boost::shared_ptr<suzerain::bspline>& bspw /* Yes, a reference */)
{
    // Only root writes data
    int procid;
    esio_handle_comm_rank(esioh, &procid);

    DEBUG0("Storing B-spline knot details");

    scoped_array<real_t> buf(new real_t[bspw->nknots()]);

    bspw->knots(buf.get(), 1);
    esio_line_establish(esioh, bspw->nknots(),
            0, (procid == 0 ? bspw->nknots() : 0));
    esio_line_write(esioh, "knots", buf.get(), 0,
            "Knots used to build B-spline basis");

    bspw->breakpoints(buf.get(), 1);
    esio_line_establish(esioh, bspw->nbreakpoints(),
            0, (procid == 0 ? bspw->nbreakpoints() : 0));
    esio_line_write(esioh, "breakpoints", buf.get(), 0,
            "Breakpoint locations used to build B-spline basis");

    bspw->collocation_points(buf.get(), 1);
    esio_line_establish(esioh, bspw->ndof(),
            0, (procid == 0 ? bspw->ndof() : 0));
    esio_line_write(esioh, "collocation_points", buf.get(), 0,
            "Collocation points used to build discrete operators");

    DEBUG0("Storing B-spline derivative operators");

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

void load(const esio_handle esioh,
          boost::shared_ptr<suzerain::bspline>& bspw /* Yes, a reference */)
{
    DEBUG0("Loading B-spline breakpoints");

    // All ranks load B-spline order
    int k;
    esio_line_establish(esioh, 1, 0, 1);
    esio_line_read(esioh, "k", &k, 0);

    // All ranks load B-spline breakpoints
    int nbreak;
    esio_line_size(esioh, "breakpoints", &nbreak);
    esio_line_establish(esioh, nbreak, 0, nbreak);
    scoped_array<real_t> buf(new double[nbreak]);
    esio_line_read(esioh, "breakpoints", buf.get(), 0);

    // Construct B-spline workspace
    bspw.reset(new suzerain::bspline(k, k - 2, nbreak, buf.get()));
}

void store_time(const esio_handle esioh,
                real_t time)
{
    // Root writes details
    int rank;
    esio_handle_comm_rank(esioh, &rank);
    esio_line_establish(esioh, 1, 0, (rank == 0) ? 1 : 0);

    esio_line_write(esioh, "t", &time, 0, "Simulation physical time");

    DEBUG0("Stored simulation time " << time);
}

void load_time(const esio_handle esioh,
               real_t &time)
{
    // All ranks read details
    esio_line_establish(esioh, 1, 0, 1);

    esio_line_read(esioh, "t", &time, 0);

    DEBUG0("Loaded simulation time " << time);
}
