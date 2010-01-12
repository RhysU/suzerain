/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2009 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * pencil_grid.cpp: Class to manage data layout concerns for P3DFFT usage
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <mpi.h>
#include <p3dfft_d.h>
#include <suzerain/mpi.hpp>
#include <suzerain/pencil_grid.hpp>

namespace suzerain {

pencil_grid::pencil_grid(const pencil_grid::size_type_3d &global_extents,
                         const pencil_grid::size_type_2d &processor_grid)
    : global_extents_(global_extents),
      processor_grid_(processor_grid)
{
    // Ensure MPI is in a sane state at construction time
    {
        int flag;
        if (MPI_Initialized(&flag) != MPI_SUCCESS || !flag) {
            throw std::logic_error("MPI stack not yet initialized");
        }
        if (MPI_Finalized(&flag) != MPI_SUCCESS || flag) {
            throw std::logic_error("MPI stack already finalized");
        }
    }

    // Find the global number of processors
    const int nproc = suzerain::mpi::comm_size(MPI_COMM_WORLD);

    // If processor grid was not fully fixed by the arguments,
    // find a 2D Cartesian decomposition where P_1 <= P_2
    if (processor_grid_[0] == 0 || processor_grid_[1] == 0) {
        suzerain::mpi::dims_create(nproc, processor_grid_);
        std::sort(processor_grid_.begin(), processor_grid_.end());
    }

    // Sanity check the 2D decomposition
    if (processor_grid_[0]*processor_grid_[1] != nproc) {
        std::ostringstream what;
        what << "Processor grid dimensions " << processor_grid_
             << "incompatible with number of processors " << nproc;
        throw std::runtime_error(what.str());
    }

    // Initialize P3DFFT using Y as STRIDE1 direction in wave space
    {
        int pg[2] = {
            boost::numeric_cast<int>(processor_grid_[0]),
            boost::numeric_cast<int>(processor_grid_[1]),
        };
        p3dfft_setup(pg,
                    boost::numeric_cast<int>(global_extents_[0]),
                    boost::numeric_cast<int>(global_extents_[2]),
                    boost::numeric_cast<int>(global_extents_[1]),
                    1 /* nuke btrans */);
        p3dfft_setup_called_ = true;
    }

    // Retrieve information for local input and output pencils
    // P3DFFT uses int types; defensively ensure we do too
    boost::array<int,3> pstart, pextent, pend, wstart, wextent, wend;
    get_dims(pstart.data(), pend.data(), pextent.data(), 1/* physical */);
    get_dims(wstart.data(), wend.data(), wextent.data(), 2/* wave */);
    // P3DFFT STRIDE1 physical space get_dims returns in (X, Z, Y) ordering
    // Suzerain's pencils require (X, Y, Z) ordering; flip Z and Y data
    std::swap(pstart [1], pstart [2]);
    std::swap(pend   [1], pend   [2]);
    std::swap(pextent[1], pextent[2]);
    // P3DFFT STRIDE1 wave space get_dims returns in (Y, X, Z) ordering
    // Suzerain's pencils require (X, Y, Z) ordering; flip Y and X data
    std::swap(wstart [0], wstart [1]);
    std::swap(wend   [0], wend   [1]);
    std::swap(wextent[0], wextent[1]);
    // Transform indices for C conventions; want ranges like [istart, iend)
    std::transform(pstart.begin(), pstart.end(), pstart.begin(),
            std::bind2nd(std::minus<int>(),1));
    std::transform(wstart.begin(), wstart.end(), wstart.begin(),
            std::bind2nd(std::minus<int>(),1));

    // Convert modified P3DFFT dimensions to have appropriate numeric type
    boost::numeric::converter<index,int> converter;
    std::transform(pstart.begin(), pstart.end(), pstart_.begin(), converter);
    std::transform(pend.begin(),   pend.end(),   pend_.begin(),   converter);
    std::transform(pextent.begin(),pextent.end(),pextent_.begin(),converter);
    std::transform(wstart.begin(), wstart.end(), wstart_.begin(), converter);
    std::transform(wend.begin(),   wend.end(),   wend_.begin(),   converter);
    std::transform(wextent.begin(),wextent.end(),wextent_.begin(),converter);
}

pencil_grid::~pencil_grid()
{
    if (p3dfft_setup_called_) p3dfft_clean();
}

} // namespace suzerain
