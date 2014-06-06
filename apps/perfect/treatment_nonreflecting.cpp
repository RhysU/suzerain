//--------------------------------------------------------------------------
//
// Copyright (C) 2013-2014 Rhys Ulerich
// Copyright (C) 2013-2014 The PECOS Development Team
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
 * @copydoc treatment_nonreflecting.hpp
 */

#include "treatment_nonreflecting.hpp"

#include <suzerain/bspline.hpp>
#include <suzerain/inorder.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/specification_grid.hpp>
#include <suzerain/timers.h>

#include "definition_scenario.hpp"

#pragma float_control(precise, on)
#pragma fenv_access(on)
#pragma float_control(except, on)
#pragma fp_contract(off)
static inline
suzerain::real_t twopiover(const suzerain::real_t L)
{
    return 2*M_PI/L;
}
#pragma float_control(except, off)
#pragma fenv_access(off)
#pragma float_control(precise, off)
#pragma fp_contract(on)

namespace suzerain {

namespace perfect {

treatment_nonreflecting::treatment_nonreflecting(
        const linearize::type& linearization,
        const definition_scenario &scenario,
        const specification_isothermal &isothermal,
        const specification_grid &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        bspline &b)
    : operator_base(grid, dgrid, cop, b)
    , scenario(scenario)
    , isothermal(isothermal)
    , linearization(linearization)
    , who("treatment_nonreflecting")
{
    // Defensively NaN working storage to reduce misuse possibility
    VL_S_RY      .setConstant(std::numeric_limits<real_t>::quiet_NaN());
    PG_BG_VL_S_RY.setConstant(std::numeric_limits<real_t>::quiet_NaN());
    PG_CG_VL_S_RY.setConstant(std::numeric_limits<real_t>::quiet_NaN());
    PG_VL_S_RY   .setConstant(std::numeric_limits<real_t>::quiet_NaN());
    inv_VL_S_RY  .setConstant(std::numeric_limits<real_t>::quiet_NaN());
}

std::vector<real_t>
treatment_nonreflecting::apply_operator(
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const lowstorage::method_interface<complex_t> &method,
            const std::size_t substep_index) const
{
    // Simply delegate all processing unless the grid is one-sided.
    if (!grid.one_sided()) {
        return N->apply_operator(time, swave, method, substep_index);
    }
    // On a one-sided grid, the work is much more interesting...
    SUZERAIN_TIMER_SCOPED("treatment_nonreflecting::apply_operator");

    // Implementation approach:
    //   1) Preserve the wavenumber-dependent state at the upper boundary
    //   2) Invoke the wrapped nonlinear operator in the usual fashion
    //   3) Build the various matrices we need from the reference state.
    //   4) Prepare pre-computable products of the various matrices
    //   5) Modify the right hand side in a wave-number dependent fashion.
    //
    // Implementation currently works for only the 5-equation case
    // (though, given correct Giles matrices, nothing prevents extension).
    SUZERAIN_ENSURE(swave.shape()[0] == 5U);

    // State enters method as coefficients in X, Y, and Z directions

    // If the k_x and k_z first-order NRBC terms are to be applied explicitly,
    // make a local, stride-1 copy of -I times upper boundary state values.
    // Notice that boundary coefficients are 1-1 with boundary point values.
    // That is, applying the mass matrix is an ignorable NOP at the boundary.
    Matrix5Xc negI_stash;
    if (linearization != linearize::rhome_xyz) {
        negI_stash.resize(NoChange, swave.shape()[2] * swave.shape()[3]);
        const complex_t negI_unit(0, -1);
        const int ku = boost::numeric_cast<int>(swave.shape()[0]);
        const int l  = boost::numeric_cast<int>(swave.shape()[1]) - 1; // Upper
        const int mu = boost::numeric_cast<int>(swave.shape()[2]);
        const int nu = boost::numeric_cast<int>(swave.shape()[3]);
        for (int k = 0; k < ku; ++k) {
            for (int n = 0; n < nu; ++n) {
                for (int m = 0; m < mu; ++m) {
                    negI_stash(k, m + mu*n) = negI_unit * swave[k][l][m][n];
                }
            }
        }
    }

    // Invoke the wrapped nonlinear operator (which may compute references!)
    const std::vector<real_t> retval = N->apply_operator(
            time, swave, method, substep_index);

    // State is now coefficients in X and Z directions
    // State is now collocation point values in Y direction

    // Prepare the matrices required to implement the boundary condition.
    // The hideous const_cast is required due to timestepping API.
    if (substep_index == 0) {
        const_cast<treatment_nonreflecting*>(this)
                ->compute_giles_matrices_upper();
    }

    // Wavenumber traversal modeled after those found in suzerain/diffwave.c
    const int Ny   = dgrid.global_wave_extent.y();
    const int dNx  = grid.dN.x();
    const int dkbx = dgrid.local_wave_start.x();
    const int dkex = dgrid.local_wave_end.x();
    const int dNz  = grid.dN.z();
    const int dkbz = dgrid.local_wave_start.z();
    const int dkez = dgrid.local_wave_end.z();
    const real_t twopioverLx_by_chi = twopiover(grid.L.x()) / dgrid.chi();
    const real_t twopioverLz_by_chi = twopiover(grid.L.z()) / dgrid.chi();

    // Traverse wavenumbers updating the RHS with the Giles boundary condition
    // TODO Could short circuit some of this traversal using dealiasing
    const Matrix5r ImPG_VL_S_RY = VL_S_RY - PG_VL_S_RY;
    if (negI_stash.size()) {

        // "Implementation primarily within the nonlinear explicit operator"
        // which is also appropriate for wall-normal implicit treatment
        assert(   linearization == linearize::none
               || linearization == linearize::rhome_y);
        using suzerain::inorder::wavenumber;
        const int mu = boost::numeric_cast<int>(swave.shape()[2]);
        for (int n = dkbz; n < dkez; ++n) {
            const int wn = wavenumber(dNz, n);
            const real_t kn_by_chi = twopioverLz_by_chi*wn;

            for (int m = dkbx; m < dkex; ++m) {
                const int wm = wavenumber(dNx, m);
                const real_t km_by_chi = twopioverLx_by_chi*wm;

                // Accumulates kn/km NRBC terms before possible ImPG swamping.
                // The -I scaling has been accommodated within negI_stash.
                Vector5c tmp   = kn_by_chi
                               * PG_BG_VL_S_RY.cast<complex_t>()
                               * negI_stash.col((m - dkbx) + mu*(n - dkbz));
                tmp.noalias() += km_by_chi
                               * PG_CG_VL_S_RY.cast<complex_t>()
                               * negI_stash.col((m - dkbx) + mu*(n - dkbz));

                // Pack upper boundary nonlinear RHS into contiguous buffer
                Vector5c N;
                N[0] = swave[ndx::e  ][Ny - 1][m - dkbx][n - dkbz];
                N[1] = swave[ndx::mx ][Ny - 1][m - dkbx][n - dkbz];
                N[2] = swave[ndx::my ][Ny - 1][m - dkbx][n - dkbz];
                N[3] = swave[ndx::mz ][Ny - 1][m - dkbx][n - dkbz];
                N[4] = swave[ndx::rho][Ny - 1][m - dkbx][n - dkbz];

                // Project out unwanted characteristics from RHS and accumulate
                // followed by projecting result back to conserved state
                tmp.noalias() += ImPG_VL_S_RY.cast<complex_t>() * N;
                N.noalias()    = inv_VL_S_RY.cast<complex_t>()  * tmp;

                // Unpack new upper boundary RHS from the contiguous buffer
                swave[ndx::e  ][Ny - 1][m - dkbx][n - dkbz] = N[0];
                swave[ndx::mx ][Ny - 1][m - dkbx][n - dkbz] = N[1];
                swave[ndx::my ][Ny - 1][m - dkbx][n - dkbz] = N[2];
                swave[ndx::mz ][Ny - 1][m - dkbx][n - dkbz] = N[3];
                swave[ndx::rho][Ny - 1][m - dkbx][n - dkbz] = N[4];

            }
        }

    } else {

        // "Implementation primarily within the linear implicit operator"
        // which is appropriate for linear implicit work in three directions.
        assert(linearization == linearize::rhome_xyz);
        const Matrix5r upper_nrbc_n = inv_VL_S_RY * ImPG_VL_S_RY;
        for (int n = dkbz; n < dkez; ++n) {
            for (int m = dkbx; m < dkex; ++m) {

                // Pack upper boundary RHS into a contiguous buffer
                Vector5c N;
                N[0] = swave[ndx::e  ][Ny - 1][m - dkbx][n - dkbz];
                N[1] = swave[ndx::mx ][Ny - 1][m - dkbx][n - dkbz];
                N[2] = swave[ndx::my ][Ny - 1][m - dkbx][n - dkbz];
                N[3] = swave[ndx::mz ][Ny - 1][m - dkbx][n - dkbz];
                N[4] = swave[ndx::rho][Ny - 1][m - dkbx][n - dkbz];

                // Modify the packed RHS to remove unwanted characteristics
                N.applyOnTheLeft(upper_nrbc_n.cast<complex_t>());

                // Unpack new upper boundary RHS from the contiguous buffer
                swave[ndx::e  ][Ny - 1][m - dkbx][n - dkbz] = N[0];
                swave[ndx::mx ][Ny - 1][m - dkbx][n - dkbz] = N[1];
                swave[ndx::my ][Ny - 1][m - dkbx][n - dkbz] = N[2];
                swave[ndx::mz ][Ny - 1][m - dkbx][n - dkbz] = N[3];
                swave[ndx::rho][Ny - 1][m - dkbx][n - dkbz] = N[4];

            }
        }

    }

    return retval;
}

void
treatment_nonreflecting::compute_giles_matrices_upper()
{
    SUZERAIN_TIMER_SCOPED("compute_giles_matrices_upper");

    giles_matrices_upper(scenario.Ma,
                         scenario.gamma,
                         isothermal.upper_rho,
                         isothermal.upper_u,
                         isothermal.upper_v,
                         isothermal.upper_w,
                         std::sqrt(isothermal.upper_T),
                         VL_S_RY,
                         PG_BG_VL_S_RY,
                         PG_CG_VL_S_RY,
                         PG_VL_S_RY,
                         inv_VL_S_RY);
}

} // namespace perfect

} // namespace suzerain
