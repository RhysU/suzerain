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
        const definition_scenario &scenario,
        const specification_isothermal &isothermal,
        const specification_grid &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        bspline &b)
    : operator_base(grid, dgrid, cop, b)
    , scenario(scenario)
    , isothermal(isothermal)
    , who("treatment_nonreflecting")
{
#ifndef NDEBUG
    // Defensively NaN working storage to reduce misuse possibility
    VL_S_RY             .setConstant(std::numeric_limits<real_t>::quiet_NaN());
    PG_BG_VL_S_RY_by_chi.setConstant(std::numeric_limits<real_t>::quiet_NaN());
    PG_CG_VL_S_RY_by_chi.setConstant(std::numeric_limits<real_t>::quiet_NaN());
    ImPG_VL_S_RY        .setConstant(std::numeric_limits<real_t>::quiet_NaN());
    inv_VL_S_RY         .setConstant(std::numeric_limits<real_t>::quiet_NaN());
#endif
}

std::vector<real_t>
treatment_nonreflecting::apply_operator(
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const lowstorage::method_interface<complex_t> &method,
            const std::size_t substep_index) const
{
    // Implementation approach:
    //   1) Preserve the wavenumber-dependent state at the upper boundary
    //   2) Invoke the wrapped nonlinear operator in the usual fashion
    //   3) Build the various matrices we need from the reference state.
    //   4) Prepare pre-computable products of the various matrices
    //   5) Modify the right hand side in a wave-number dependent fashion.

    // State enters method as coefficients in X, Y, and Z directions

    // Make a local, stride-1 copy of -I times upper boundary state values.
    // Notice that boundary coefficients are 1-1 with boundary point values.
    // That is, applying the mass matrix is an ignorable NOP at the boundary.
    Matrix5Xc negI_stash(5, swave.shape()[2] * swave.shape()[3]);
    {
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
    const real_t twopioverLx = twopiover(grid.L.x());  // Weird looking...
    const real_t twopioverLz = twopiover(grid.L.z());  // ...for FP control

    // Traverse wavenumbers updating the RHS with the Giles boundary condition
    // TODO Could short circuit some of this traversal using dealiasing
    using suzerain::inorder::wavenumber;
    const int mu = boost::numeric_cast<int>(swave.shape()[2]);
    for (int n = dkbz; n < dkez; ++n) {
        const int wn = wavenumber(dNz, n);
        const real_t kn = twopioverLz*wn;

        for (int m = dkbx; m < dkex; ++m) {
            const int wm = wavenumber(dNx, m);
            const real_t km = twopioverLx*wm;

            // Unpack upper boundary RHS into a contiguous buffer
            Vector5c N;
            for (int f = 0; f < N.size(); ++f) {
                N(f) = swave[f][Ny - 1][m - dkbx][n - dkbz];
            }

            // Modify the packed RHS per
            // "Implementation primarily within the nonlinear explicit operator"
            // Accumulates kn/km terms before any possible ImPG swamping.
            // The -I factors have already been included within negI_stash.
            Vector5c tmp   = kn
                           * PG_BG_VL_S_RY_by_chi.cast<complex_t>()
                           * negI_stash.col((m - dkbx) + mu*(n - dkbz));
            tmp.noalias() += km
                           * PG_CG_VL_S_RY_by_chi.cast<complex_t>()
                           * negI_stash.col((m - dkbx) + mu*(n - dkbz));
            tmp.noalias() += ImPG_VL_S_RY.cast<complex_t>() * N;
            N.noalias()    = inv_VL_S_RY.cast<complex_t>() * tmp;

            // Pack new upper boundary RHS from the contiguous buffer
            for (int f = 0; f < N.size(); ++f) {
                swave[f][Ny - 1][m - dkbx][n - dkbz] = N(f);
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
                         PG_BG_VL_S_RY_by_chi, // Rescaled below
                         PG_CG_VL_S_RY_by_chi, // Rescaled below
                         ImPG_VL_S_RY,         // Adjusted below
                         inv_VL_S_RY);
    PG_BG_VL_S_RY_by_chi /= dgrid.chi();
    PG_CG_VL_S_RY_by_chi /= dgrid.chi();
    ImPG_VL_S_RY = VL_S_RY - ImPG_VL_S_RY;
}

} // namespace perfect

} // namespace suzerain
