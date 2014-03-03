//--------------------------------------------------------------------------
//
// Copyright (C) 2008-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
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

/** @file
 * @copydoc extrema.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/extrema.hpp>

#include <suzerain/bspline.hpp>
#include <suzerain/diffwave.hpp>
#include <suzerain/error.h>
#include <suzerain/math.hpp>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/operator_tools.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/physical_view.hpp>
#include <suzerain/specification_grid.hpp>
#include <suzerain/state.hpp>

namespace suzerain {

std::vector<field_extrema_xz>
compute_field_extrema_xz(
        contiguous_state<4,complex_t>& swave,
        const specification_grid& grid,
        const pencil_grid& dgrid,
        const bsplineop& cop)
{
    // State enters as coefficients in X, Y, and Z stored swave[F][Y][X][Z]

    // Contiguous temporary storage for accumulating/broadcasting results
    // Each fast columns stores one wall-normal swaths of data for a scalar
    ArrayXXr buf(swave.shape()[1], 2*swave.shape()[0]);
    ArrayXXr::ColsBlockXpr min = buf.leftCols (swave.shape()[0]);
    ArrayXXr::ColsBlockXpr max = buf.rightCols(swave.shape()[0]);

    // Initialize storage with "uninteresting" min and max values
    // so ranks with "gaps" in portions of the domain do not
    // modify the final Allreduce.
    min.setConstant(std::numeric_limits<ArrayXXr::Scalar>::max());
    max.setConstant(std::numeric_limits<ArrayXXr::Scalar>::min());

    // For each field...
    operator_tools o(grid, dgrid, cop);
    physical_view<> sphys(dgrid, swave);
    for (int f = 0; f < boost::numeric_cast<int>(swave.shape()[0]); ++f) {

        // Transform to physical space (dealiasing along the way)...
        o.zero_dealiasing_modes(swave, f);
        o.bop_apply(0, 1.0, swave, f);
        dgrid.transform_wave_to_physical(&sphys.coeffRef(f,0));

        // Loop over Y...
        for (int offset = 0, j = dgrid.local_physical_start.y();
            j < dgrid.local_physical_end.y();
            ++j) {

            // Loop over X and Z...
            const int last_zxoffset = offset
                                    + dgrid.local_physical_extent.z()
                                    * dgrid.local_physical_extent.x();
            for (; offset < last_zxoffset; ++offset) {

                // Retrieve scalar-valued state at point...
                const real_t val = sphys(f, offset);

                // ...and tracking minimum value at this wall-normal spot
                // being careful to preserve NaNs during local reduction.
                min(j,f) = math::minnan(val, min(j,f));
                max(j,f) = math::maxnan(val, max(j,f));

            }

        }

    }

    // Allreduce the information in a single operation playing a (-max) trick
    max *= -1;
    SUZERAIN_MPICHKR(MPI_Allreduce(MPI_IN_PLACE, buf.data(), buf.size(),
                                   mpi::datatype_of(buf.data()),
                                   MPI_MIN, MPI_COMM_WORLD));
    max *= -1;

    // Copy the wall-normal profiles into the returned vector...
    std::vector<field_extrema_xz> retval(swave.shape()[0]);
    for (size_t f = 0; f < swave.shape()[0]; ++f) {
        retval[f].min = min.col(f);
        retval[f].max = max.col(f);
    }

    // ...and hope the compiler is smart enough to elide the copy.
    return retval;
}

std::vector<field_extrema_xyz>
compute_field_extrema_xyz(
        contiguous_state<4,complex_t>& swave,
        const specification_grid& grid,
        const pencil_grid& dgrid,
        const bsplineop& cop)
{
    // State enters as coefficients in X, Y, and Z stored swave[F][Y][X][Z]

    // TODO Rather than delegating, it may be faster to inline and modify
    std::vector<field_extrema_xz> wallnormal
        = compute_field_extrema_xz(swave, grid, dgrid, cop);

    // Reduce the wall-normal profiles down to a single value
    // TODO Does this honor NaNs correctly?
    std::vector<field_extrema_xyz> retval(wallnormal.size());
    for (size_t f = 0; f < swave.shape()[0]; ++f) {
        retval[f].min = wallnormal[f].min.minCoeff();
        retval[f].max = wallnormal[f].max.maxCoeff();
    }

    return retval;
}

} // end namespace suzerain
