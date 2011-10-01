//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2010 The PECOS Development Team
//
// Please see http://pecos.ices.utexas.edu for more information.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------
//
// pencil_grid.cpp: Class to manage data layout concerns for P3DFFT usage
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/mpi.hpp>
#include <suzerain/functional.hpp>
#include <suzerain/pencil_grid.hpp>
#include <p3dfft_d.h>

namespace suzerain {

#pragma warning(push, disable:2022)
void pencil_grid_p3dfft::construct_(int Nx, int Ny, int Nz, int Pa, int Pb)
#pragma warning(pop)
{
    global_physical_extent[0] = Nx;
    global_physical_extent[1] = Ny;
    global_physical_extent[2] = Nz;
    global_wave_extent[0]     = Nx / 2 + 1;
    global_wave_extent[1]     = Ny;
    global_wave_extent[2]     = Nz;
    processor_grid[0]         = Pa;
    processor_grid[1]         = Pb;

    // Ensure MPI is in a sane state
    {
        int flag;
        if (MPI_Initialized(&flag) != MPI_SUCCESS || !flag) {
            throw std::logic_error(
                    std::string("MPI stack not yet initialized"));
        }
        if (MPI_Finalized(&flag) != MPI_SUCCESS || flag) {
            throw std::logic_error(
                    std::string("MPI stack already finalized"));
        }
    }

    // Find the global number of processors
    const int nproc = suzerain::mpi::comm_size(MPI_COMM_WORLD);

    // Find a suitable 2D Cartesian decomposition if necessary
    if (processor_grid[0] == 0 || processor_grid[1] == 0) {
        suzerain::mpi::dims_create(
                nproc,
                processor_grid.data(),
                processor_grid.data() + processor_grid.size());
    }

    // If neither Pa nor Pb was specified, enforce P_2 >= P_1
    if (Pa == 0 && Pb == 0) {
        std::sort(processor_grid.data(),
                  processor_grid.data() + processor_grid.size());
    }

    // Sanity check the 2D decomposition
    if (processor_grid.prod() != nproc) {
        std::ostringstream what;
        what << "Processor grid dimensions " << processor_grid
             << " incompatible with number of processors " << nproc;
        throw std::runtime_error(what.str());
    }

    // Initialize P3DFFT using Y as STRIDE1 direction in wave space
    {
        p3dfft_setup(processor_grid.data(),
                     global_physical_extent[0],
                     global_physical_extent[2],
                     global_physical_extent[1],
                     1); // nuke btrans
        p3dfft_setup_called_ = true;
    }

    // Retrieve information for local input and output pencils
    // P3DFFT uses int types; defensively ensure we do too
    get_dims(local_physical_start.data(),
             local_physical_end.data(),
             local_physical_extent.data(),
             1 /* physical */);
    get_dims(local_wave_start.data(),
             local_wave_end.data(),
             local_wave_extent.data(),
             2 /* wave */);
    // P3DFFT STRIDE1 physical space get_dims returns in (X, Z, Y) ordering
    // Suzerain's pencils require (X, Y, Z) ordering; flip Z and Y data
    std::swap(local_physical_start [1], local_physical_start [2]);
    std::swap(local_physical_end   [1], local_physical_end   [2]);
    std::swap(local_physical_extent[1], local_physical_extent[2]);
    // P3DFFT STRIDE1 wave space get_dims returns in (Y, X, Z) ordering
    // Suzerain's pencils require (X, Y, Z) ordering; flip Y and X data
    std::swap(local_wave_start [0], local_wave_start [1]);
    std::swap(local_wave_end   [0], local_wave_end   [1]);
    std::swap(local_wave_extent[0], local_wave_extent[1]);
    // Transform indices for C conventions; want ranges like [istart, iend)
    local_physical_start -= 1;
    local_wave_start     -= 1;
}

#pragma warning(push,disable:2017)
pencil_grid_p3dfft::~pencil_grid_p3dfft()
{
    if (p3dfft_setup_called_) p3dfft_clean();
}
#pragma warning(pop)

std::size_t pencil_grid_p3dfft::local_physical_storage() const
{
    // Wave space scalars are twice as big as physical space scalars
    return std::max(2*local_wave_extent.prod(), local_physical_extent.prod());
}

std::size_t pencil_grid_p3dfft::local_wave_storage() const
{
    // Physical space scalars are half as big, but we must round up
    const std::size_t prodphys = local_physical_extent.prod();
    const std::size_t prodwave = local_wave_extent.prod();
    return std::max(prodwave, prodphys/2 + (prodphys % 2));
}

} // namespace suzerain
