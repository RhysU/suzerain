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
#include <new>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/functional.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/pencil_grid.hpp>

namespace suzerain {

#ifdef HAVE_P3DFFT ///////////////////////////////////////////////////////////

#pragma warning(push, disable:2022)
void pencil_grid_p3dfft::construct_(int Nx, int Ny, int Nz, int Pa, int Pb,
                                    unsigned rigor_fft, unsigned rigor_mpi)
#pragma warning(pop)
{
    SUZERAIN_UNUSED(rigor_mpi);
    p3dfft_fftw_rigor(rigor_fft);

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

    // Ensure design assumptions hold for linked P3DFFT library
    if (!p3dfft_using_stride1())
        throw std::runtime_error("!p3dfft_using_stride1()");
    if (p3dfft_get_precision() != 2)
        throw std::runtime_error("p3dfft_get_precision() != 2");

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

void pencil_grid_p3dfft::transform_wave_to_physical(double * inout) const
{
    p3dfft_btran_c2r(inout, inout);
}

void pencil_grid_p3dfft::transform_physical_to_wave(double * inout) const
{
    p3dfft_ftran_r2c(inout, inout);
}

#endif // HAVE_P3DFFT ////////////////////////////////////////////////////////


#ifdef HAVE_UNDERLING ////////////////////////////////////////////////////////

pencil_grid_underling::~pencil_grid_underling()
{
    underling_only_cleanup();
}

std::size_t pencil_grid_underling::local_wave_storage() const
{
    // Convert from real to complex sizes, rounding up odd sizes
    const std::size_t n = problem->local_memory();
    return (n / 2) + (n % 2);
}

std::size_t pencil_grid_underling::local_physical_storage() const
{
    // Return twice the computed wave storage (for consistency)
    return 2 * local_wave_storage();
}

void
pencil_grid_underling::transform_wave_to_physical(double * inout) const
{
    transpose->execute_long_n0_to_long_n1(inout, buf.get());
    n1_c2c_backward->execute(buf.get(), inout);
    transpose->execute_long_n1_to_long_n2(inout, buf.get());
    n2_c2r_backward->execute(buf.get(), inout);
}

void
pencil_grid_underling::transform_physical_to_wave(double * inout) const
{
    n2_r2c_forward->execute(inout, buf.get());
    transpose->execute_long_n2_to_long_n1(buf.get(), inout);
    n1_c2c_forward->execute(inout, buf.get());
    transpose->execute_long_n1_to_long_n0(buf.get(), inout);
}

void
pencil_grid_underling::construct_(int Nx, int Ny, int Nz, int Pa, int Pb,
                                  unsigned rigor_fft, unsigned rigor_mpi)
{
    using std::bad_alloc;
    using std::runtime_error;

    // Initialize underling (and all dependencies) as necessary
    underling_init(NULL, NULL, 0);

    // Construct grid and problem
    // Map from pencil_grid to underling's directional expectations
    // long_n0 will be wave space and long_n2 will be physical space
    grid.reset(new underling::grid(MPI_COMM_WORLD, Ny, Nz, Nx/2+1, Pa, Pb));
    if (!grid)
        throw std::runtime_error("underling::grid creation error");

    problem.reset(new underling::problem(*grid, 2,
        underling::transposed::long_n2 & underling::transposed::long_n0));
    if (!problem)
        throw std::runtime_error("underling::problem creation error");

    // Allocate zeroed buffers buf and tmp.  Clearing necessary to avoid NaNs.
    // buf is maintained for execution while tmp is used for planning only.
    namespace blas = suzerain::blas;
    buf.reset(blas::calloc_as<underling::real>(problem->local_memory()),
              blas::free);
    if (!buf) throw bad_alloc();
    boost::shared_array<underling::real> tmp(
            blas::calloc_as<underling::real>(problem->local_memory()),
            blas::free);
    if (!tmp) throw bad_alloc();

    // Prepare execution plans
    transpose.reset(new underling::plan(
            *problem, buf.get(), tmp.get(), 0, rigor_mpi));
    if (!transpose)
        throw runtime_error("underling::plan creation failed");

    n1_c2c_backward.reset(new underling::fftw::plan(
            underling::fftw::plan::c2c_backward(),
            *problem, 1, buf.get(), tmp.get(), rigor_fft,
            underling::fftw::packed::all));
    if (!n1_c2c_backward)
        throw runtime_error("n1_c2c_backward creation failed");

    n2_c2r_backward.reset(new underling::fftw::plan(
            underling::fftw::plan::c2r_backward(),
            *problem, 2, buf.get(), tmp.get(), rigor_fft,
            underling::fftw::packed::all));
    if (!n2_c2r_backward)
        throw runtime_error("n2_c2r_backward creation failed");

    n2_r2c_forward.reset(new underling::fftw::plan(
            *n2_c2r_backward, buf.get(), tmp.get(), rigor_fft));
    if (!n2_r2c_forward)
        throw runtime_error("n2_r2c_forward creation failed");

    n1_c2c_forward.reset(new underling::fftw::plan(
            *n1_c2c_backward, buf.get(), tmp.get(), rigor_fft));
    if (!n1_c2c_forward)
        throw runtime_error("n1_c2c_forward creation failed");

    // Copy underling information into format pencil_grid expects
    global_physical_extent << Nx,           Ny, Nz;
    global_wave_extent     << (Nx / 2 + 1), Ny, Nz;
    processor_grid         << grid->pA_size(), grid->pB_size();

    // Prepare information on wave space storage
    {
        const underling::extents e = problem->local_extents(0);
        local_wave_start[0] = e.start[2];
        local_wave_start[1] = e.start[0];
        local_wave_start[2] = e.start[1];
        local_wave_end[0]   = e.start[2] + e.size[2];
        local_wave_end[1]   = e.start[0] + e.size[0];
        local_wave_end[2]   = e.start[1] + e.size[1];
        local_wave_extent   = local_wave_end - local_wave_start;
    }

    // Prepare information on physical space storage
    {
        const underling::fftw::extents e
                = n2_c2r_backward->local_extents_output();
        local_physical_start[0] = e.start[2];
        local_physical_start[1] = e.start[0];
        local_physical_start[2] = e.start[1];
        local_physical_end[0]   = e.start[2] + e.size[2];
        local_physical_end[1]   = e.start[0] + e.size[0];
        local_physical_end[2]   = e.start[1] + e.size[1];
        local_physical_extent   = local_physical_end - local_physical_start;
    }

    // TODO Ensure wave space storage is column-major, contiguous Y X/2 Z
    // TODO Ensure wave space storage is column-major, contiguous X Z Y
}

#endif // HAVE_UNDERLING /////////////////////////////////////////////////////


} // namespace suzerain
