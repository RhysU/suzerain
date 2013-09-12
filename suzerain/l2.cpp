//--------------------------------------------------------------------------
//
// Copyright (C) 2008, 2009, 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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
 * @copydoc l2.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/l2.hpp>

#include <suzerain/common.hpp>
#include <suzerain/error.h>
#include <suzerain/inorder.hpp>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/mpi.hpp>

namespace suzerain {

std::vector<field_L2xyz>
compute_field_L2xyz(
        const contiguous_state<4,complex_t> &state,
        const grid_specification& grid,
        const pencil_grid& dgrid,
        const bsplineop& gop)
{
    // Ensure state storage meets this routine's assumptions
    // Notice state.shape()[0] may be any value
    using boost::numeric_cast;
    assert(numeric_cast<int>(state.shape()[1]) == dgrid.local_wave_extent.y());
    assert(numeric_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
    assert(numeric_cast<int>(state.shape()[3]) == dgrid.local_wave_extent.z());

    // Ensure we were handed Galerkin field_L2xyz operator matrices
    assert(gop.get()->method == SUZERAIN_BSPLINEOP_GALERKIN_L2);

    // Only want non-dealiased X-direction modes to contribute to L2
    // Compute wavenumber translation logistics for X direction
    int fxb[2], fxe[2], mxb[2], mxe[2];
    inorder::wavenumber_translate(grid.N.x(),
                                  grid.dN.x(),
                                  dgrid.local_wave_start.x(),
                                  dgrid.local_wave_end.x(),
                                  fxb[0], fxe[0], fxb[1], fxe[1],
                                  mxb[0], mxe[0], mxb[1], mxe[1]);
    // X contains only positive wavenumbers => second range must be empty
    assert(fxb[1] == fxe[1]);
    assert(mxb[1] == mxe[1]);

    // Only want non-dealiased Z-direction modes to contribute to L2
    // Compute wavenumber translation logistics for Z direction
    // One or both ranges may be empty
    int fzb[2], fze[2], mzb[2], mze[2];
    inorder::wavenumber_translate(grid.N.z(),
                                  grid.dN.z(),
                                  dgrid.local_wave_start.z(),
                                  dgrid.local_wave_end.z(),
                                  fzb[0], fze[0], fzb[1], fze[1],
                                  mzb[0], mze[0], mzb[1], mze[1]);

    // Temporary storage for inner product computations
    VectorXc tmp;
    tmp.setZero(grid.N.y());

    // Contiguous temporary storage for accumulating and broadcasting results
    VectorXc buf(2*state.shape()[0]);
    Map<VectorXc> total2(buf.data() + 0,                state.shape()[0]);
    Map<VectorXc> mean2 (buf.data() + state.shape()[0], state.shape()[0]);

    // Compute the local L2 contribution towards each L^2 norm squared
    // Computation uses partial sums at each loop to reduce swamping which is
    // more-or-less recursive summation using large partitioning factors.
    total2.setZero();
    for (size_t k = 0; k < state.shape()[0]; ++k) {
        complex_t jsum = 0;
        for (int j = 0; j < 2; ++j) {
            complex_t nsum = 0;
            for (int n = mzb[j]; n < mze[j]; ++n) {
                complex_t msum = 0;
                for (int m = mxb[0]; m < mxe[0]; ++m) {
                    const complex_t * u_mn
                        = &state[k][0][m - dgrid.local_wave_start.x()]
                                      [n - dgrid.local_wave_start.z()];
                    gop.accumulate(0, 1.0, u_mn, 1, 0.0, tmp.data(), 1);
                    complex_t dot = blas::dot(grid.N.y(), u_mn,       1,
                                                          tmp.data(), 1);
                    if (m > 0 && m < grid.dN.x()/2) {
                        dot *= 2;
                    }
                    msum += dot;
                }
                nsum += msum;
            }
            jsum += nsum;
        }
        total2[k] += jsum;
    }
    total2 *= grid.L.x() * grid.L.z();

    // Reduce total2 sum onto processor housing the zero-zero mode using
    // mean2 as a scratch buffer to simulate MPI_IN_PLACE
    SUZERAIN_MPICHKR(MPI_Reduce(total2.data(),
            mean2.data(), state.shape()[0]*sizeof(complex_t)/sizeof(real_t),
            mpi::datatype<real_t>(),
            MPI_SUM, dgrid.rank_zero_zero_modes, MPI_COMM_WORLD));
    total2 = mean2;

    // Compute the mean-only L^2 squared for each field using zero-zero modes
    if (dgrid.has_zero_zero_modes()) {
        for (size_t k = 0; k < state.shape()[0]; ++k) {
            const complex_t * u_mn = &state[k][0][0][0];
            gop.accumulate(0, 1.0, u_mn, 1, 0.0, tmp.data(), 1);
            mean2[k] = blas::dot(grid.N.y(), u_mn, 1, tmp.data(), 1);
        }
        mean2 *= grid.L.x() * grid.L.z();
    }

    // Broadcast total2 and mean2 values to all processors
    SUZERAIN_MPICHKR(MPI_Bcast(
            buf.data(), buf.size() * sizeof(complex_t)/sizeof(real_t),
            mpi::datatype<real_t>(),
            dgrid.rank_zero_zero_modes, MPI_COMM_WORLD));

    // Obtain fluctuating2 = total2 - mean2 and pack the return structure
    std::vector<field_L2xyz> retval(state.shape()[0]);
    for (size_t k = 0; k < state.shape()[0]; ++k) {
        retval[k].mean2        = std::abs(mean2[k]);
        retval[k].fluctuating2 = std::abs(total2[k] - mean2[k]);
    }

    return retval;
}

} // end namespace channel
