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

struct val_index {
  double val;    /**< Value. */
  int    index;  /**< Global xz index. */
};

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
    int ny = boost::numeric_cast<int>(swave.shape()[1]);
    int nf = boost::numeric_cast<int>(swave.shape()[0]);
    int buf2size = ny * 2*nf;
    suzerain::scoped_array<val_index> buf2(new val_index[buf2size]);

    // Initialize storage with "uninteresting" min and max values
    // so that the final Allreduce ignores contributions from
    // ranks that have no business contributing to certain y_j.
    for (int i = 0; i < buf2size/2; ++i) {
        buf2[i].val   = +std::numeric_limits<real_t>::infinity();
        buf2[i].index = 0;
    }
    for (int i = buf2size/2; i < buf2size; ++i) {
        buf2[i].val   = -std::numeric_limits<real_t>::infinity();
        buf2[i].index = 0;
    }

    // For each field...
    operator_tools o(grid, dgrid, cop);
    physical_view<> sphys(dgrid, swave);
    int dNx = dgrid.global_physical_extent.x();
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
            for (int k = o.dgrid.local_physical_start.z();
                k < o.dgrid.local_physical_end.z();
                ++k) {

                for (int i = o.dgrid.local_physical_start.x();
                    i < o.dgrid.local_physical_end.x();
                    ++i, /* NB */ ++offset) {

                    int xz_global_offset = i + k * dNx;

                    // Retrieve scalar-valued state at point...
                    const real_t val = sphys(f, offset);

                    // ...and tracking minimum value at this wall-normal spot
                    // being careful to preserve NaNs during local reduction.

                    // minimum value
                    if (val >= buf2[j+f*ny].val) {
                        // do nothing
                    } else {
                        buf2[j+f*ny].val   = val;
                        buf2[j+f*ny].index = xz_global_offset;
                    }

                    // maximum value
                    if (val <= buf2[j+(f+nf)*ny].val) {
                        // do nothing
                    } else {
                        buf2[j+(f+nf)*ny].val   = val;
                        buf2[j+(f+nf)*ny].index = xz_global_offset;
                    }

                }

            }

        }

    }

    // Allreduce the information in a single operation playing a (-max) trick
    for (int i = buf2size/2+1; i < buf2size; ++i) {
        buf2[i].val *= -1;
    }
    SUZERAIN_MPICHKR(MPI_Allreduce(MPI_IN_PLACE, &buf2[0], buf2size,
                                   MPI_DOUBLE_INT,
                                   MPI_MINLOC, MPI_COMM_WORLD));
    for (int i = buf2size/2+1; i < buf2size; ++i) {
        buf2[i].val *= -1;
    }

    // Copy the wall-normal profiles into the returned vector...
    std::vector<field_extrema_xz> retval(swave.shape()[0]);
    for (int f = 0; f < nf; ++f) {
        retval[f].min.resize(ny);
        retval[f].max.resize(ny);
        retval[f].imin.resize(ny);
        retval[f].imax.resize(ny); 
        retval[f].kmin.resize(ny);
        retval[f].kmax.resize(ny); 

        for (int j = 0; j < ny; ++j) {
            // value
            retval[f].min(j)  = buf2[    f *ny + j].val;
            retval[f].max(j)  = buf2[(nf+f)*ny + j].val;

            // x-index =  xz_global % dNx
            retval[f].imin(j) = buf2[    f *ny + j].index % dNx;
            retval[f].imax(j) = buf2[(nf+f)*ny + j].index % dNx;

            // z-index = (xz_global - i-index) / dNx
            retval[f].kmin(j) = (buf2[    f *ny + j].index 
                                 - retval[f].imin(j)) / dNx;
            retval[f].kmax(j) = (buf2[(nf+f)*ny + j].index 
                                 - retval[f].imax(j)) / dNx;
        }
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
    std::vector<field_extrema_xyz> retval(wallnormal.size());
    int ny = boost::numeric_cast<int>(swave.shape()[1]);
    for (size_t f = 0; f < swave.shape()[0]; ++f) {
        retval[f].min   = +std::numeric_limits<real_t>::infinity();
        retval[f].max   = -std::numeric_limits<real_t>::infinity();

        for (int j = 0; j < ny; ++j) {

            // minimum value
            if (wallnormal[f].min(j) >= retval[f].min) {
                // do nothing
            } else {
                retval[f].min  = wallnormal[f].min(j);
                retval[f].imin = wallnormal[f].imin(j);
                retval[f].jmin = j;
                retval[f].kmin = wallnormal[f].kmin(j);
            }

            // maximum value
            if (wallnormal[f].max(j) <= retval[f].max) {
                // do nothing
            } else {
                retval[f].max  = wallnormal[f].max(j);
                retval[f].imax = wallnormal[f].imax(j);
                retval[f].jmax = j;
                retval[f].kmax = wallnormal[f].kmax(j);
            }
        }
    }

    return retval;
}

} // end namespace suzerain
