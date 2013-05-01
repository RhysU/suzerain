//--------------------------------------------------------------------------
//
// Copyright (C) 2013 Rhys Ulerich
// Copyright (C) 2013 The PECOS Development Team
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
 * @copydoc nonreflecting_treatment.hpp
 */

#include "nonreflecting_treatment.hpp"

#include <suzerain/bspline.hpp>
#include <suzerain/grid_specification.hpp>
#include <suzerain/inorder.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/pencil_grid.hpp>

#include "common_block.hpp"
#include "scenario_definition.hpp"

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

nonreflecting_treatment::nonreflecting_treatment(
        const scenario_definition &scenario,
        const grid_specification &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        bspline &b,
        operator_common_block &common)
    : operator_base(grid, dgrid, cop, b)
    , scenario(scenario)
    , common(common)
    , who("nonreflecting_treatment")
{
    // NOP
}

std::vector<real_t> nonreflecting_treatment::apply_operator(
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag,
            const std::size_t substep_index) const
{
    // Implementation approach:
    //   1) Preserve the wavenumber-dependent state at the upper boundary
    //   2) Invoke the wrapped nonlinear operator in the usual fashion
    //   3) Build the various matrices we need from the reference state.
    //   4) Prepare pre-computable products of the various matrices
    //   5) Modify the right hand side in a wave-number dependent fashion.

    // State enters method as coefficients in X, Y, and Z directions

    // Make a local, stride-1 copy of I times upper boundary state point values.
    // Notice that boundary coefficients are 1-1 with boundary point values.
    // That is, applying the mass matrix to the boundary is an ignorable NOP.
    Matrix5Xc i_stash(5, swave.shape()[2] * swave.shape()[3]);
    {
        const complex_t imag_unit(0, 1);
        const int ku = boost::numeric_cast<int>(swave.shape()[0]);
        const int l  = boost::numeric_cast<int>(swave.shape()[1]) - 1; // Upper
        const int mu = boost::numeric_cast<int>(swave.shape()[2]);
        const int nu = boost::numeric_cast<int>(swave.shape()[3]);
        for (int k = 0; k < ku; ++k) {
            for (int n = 0; n < nu; ++n) {
                for (int m = 0; m < mu; ++m) {
                    i_stash(k, m + n * mu) = imag_unit * swave[k][l][m][n];
                }
            }
        }
    }

    // Invoke the wrapped nonlinear operator (which may compute references!)
    const std::vector<real_t> retval = N->apply_operator(
            time, swave, evmaxmag_real, evmaxmag_imag, substep_index);

    // State is now coefficients in X and Z directions
    // State is now collocation point values in Y direction

    // Prepare the rotation and its inverse that reorders from
    // ndx::{e, mx, my, mz, rho} to {rho = 0, my = 1, mz = 2, mx = 3, e = 4}.
    // All remaining matrices/logic within the routine uses the latter order!
    Matrix5r RY(Matrix5r::Zero());
    {
        RY(0, ndx::rho) = 1;
        RY(1, ndx::my ) = 1;
        RY(2, ndx::mz ) = 1;
        RY(3, ndx::mx ) = 1;
        RY(4, ndx::e  ) = 1;
    }
    Matrix5r inv_RY(Matrix5r::Zero());
    {
        inv_RY(ndx::e  , 4) = 1;
        inv_RY(ndx::mx , 3) = 1;
        inv_RY(ndx::my , 1) = 1;
        inv_RY(ndx::mz , 2) = 1;
        inv_RY(ndx::rho, 0) = 1;
    }

    // Retrieve upper boundary reference state from the common block.
    // These assignments are odd looking to permit others as documented.
    const real_t rho = common.ref_rho().tail<1>()[0];
    const real_t u   = common.ref_uy ().tail<1>()[0]; // Ref u = v' per RY
    const real_t v   = common.ref_uz ().tail<1>()[0]; // Ref v = w' per RY
    const real_t w   = common.ref_ux ().tail<1>()[0]; // Ref w = u' per RY
          real_t a   = common.ref_a  ().tail<1>()[0];

    // Prepare oft-used quantities derived from the reference state
    const real_t inv_rho = 1 / rho;
    const real_t u2      = u * u;
    const real_t v2      = v * v;
    const real_t w2      = w * w;
          real_t inv_a   = 1 / a;
          real_t a2      = a * a;
          real_t inv_a2  = 1 / a2;

    // Prepare oft-used scenario-related constants
    const real_t chi        = dgrid.chi();
    const real_t half       = static_cast<real_t>(1) / 2;
    const real_t gamma      = scenario.gamma;
    const real_t inv_gamma  = 1 / gamma;
    const real_t gamma1     = gamma - 1;
    const real_t inv_gamma1 = 1 / gamma1;
    const real_t Ma         = scenario.Ma;
    const real_t inv_Ma     = 1 / Ma;
    const real_t Ma2        = Ma * Ma;
    const real_t inv_Ma2    = 1 / Ma2;

    // Build the conserved -> primitive transformation and its inverse
    Matrix5r S(Matrix5r::Zero());
    {
        S(0, 0) =   1;
        S(1, 0) = - u * inv_rho;
        S(2, 0) = - v * inv_rho;
        S(3, 0) = - w * inv_rho;
        S(4, 0) =   a2 * inv_gamma * inv_Ma2;
        S(1, 1) =   inv_rho;
        S(4, 1) = - gamma1 * u;
        S(2, 2) =   inv_rho;
        S(4, 2) = - gamma1 * v;
        S(3, 3) =   inv_rho;
        S(4, 3) = - gamma1 * w;
        S(4, 4) =   gamma1 * inv_Ma2;
    }
    Matrix5r inv_S(Matrix5r::Zero());
    {
        inv_S(0, 0) =   1;
        inv_S(1, 0) =   u;
        inv_S(2, 0) =   v;
        inv_S(3, 0) =   w;
        inv_S(4, 0) =   Ma2 * (u2 + v2 + w2) + a2 * inv_gamma * inv_gamma1;
        inv_S(1, 1) =   rho;
        inv_S(4, 1) =   Ma2 * rho * u;
        inv_S(2, 2) =   rho;
        inv_S(4, 2) =   Ma2 * rho * v;
        inv_S(3, 3) =   rho;
        inv_S(4, 3) =   Ma2 * rho * w;
        inv_S(4, 4) =   Ma2 * inv_gamma1;
    }

    // Per the model document
    //   When reusing [V^L, B^G, and C^G] derived for U ... for V^\ast
    //   every sound speed must be scaled by 1 / Ma because
    //   \bar{a} / u_0 = a_0 \bar{a}^\ast / u_0 = \bar{a}^\ast / Ma.
    // so we now once-and-for-all adjust a and a2 appropriately.
    a      *= inv_Ma;
    inv_a  *= Ma;
    a2     *= inv_Ma2;
    inv_a2 *= Ma2;

    // Build the primitive -> characteristic transformation and its inverse
    Matrix5r VL(Matrix5r::Zero());
    {
        VL(0, 0) = - a2;
        VL(3, 1) =   rho * a;
        VL(4, 1) = - rho * a;
        VL(1, 2) =   rho * a;
        VL(2, 3) =   rho * a;
        VL(0, 4) =   1;
        VL(3, 4) =   1;
        VL(4, 4) =   1;
    }
    Matrix5r inv_VL(Matrix5r::Zero());
    {
        inv_VL(0, 0) = - inv_a2;
        inv_VL(2, 1) =   inv_rho * inv_a;
        inv_VL(3, 2) =   inv_rho * inv_a;
        inv_VL(0, 3) =   half * inv_a2;
        inv_VL(1, 3) =   half * inv_rho * inv_a;
        inv_VL(4, 3) =   half;
        inv_VL(0, 4) =   half * inv_a2;
        inv_VL(1, 4) = - half * inv_rho * inv_a;
        inv_VL(4, 4) =   half;
    }

    // The upper boundary is an inflow when the reference velocity is negative.
    // In the event of u == 0, treat the boundary like an outflow.
    const bool inflow = u < 0;

    // Build the in-vs-outflow characteristic-preserving projection
    Matrix5r PG(Matrix5r::Zero());
    if (inflow) {
        PG(0, 0) = 1;
        PG(1, 1) = 1;
        PG(2, 2) = 1;
        PG(3, 3) = 1;
    } else {
        PG(4, 4) = 1;
    }

    Matrix5r BG(Matrix5r::Zero());  // Medida's B^G_1
    if (inflow) {
        BG(1, 1) = v;
        BG(3, 1) = half * (a - u);
        BG(2, 2) = v;
        BG(1, 3) = half * (a + u);
        BG(3, 3) = v;
        BG(1, 4) = half * (a - u);
    } else {
        BG(4, 1) = u;
        BG(4, 4) = v;
    }

    Matrix5r CG(Matrix5r::Zero());  // Medida's C^G_1
    if (inflow) {
        CG(1, 1) = w;
        CG(2, 2) = w;
        CG(3, 2) = half * (a - u);
        CG(2, 3) = half * (a + u);
        CG(3, 3) = w;
        CG(2, 4) = half * (a - u);
    } else {
        CG(4, 2) = u;
        CG(4, 4) = w;
    }

    // Prepare all necessary real-valued products of the above matrices
    const Matrix5r VL_S_RY           = VL * S * RY;
    const Matrix5r BG_VL_S_RY_by_chi = BG * VL_S_RY / chi;
    const Matrix5r CG_VL_S_RY_by_chi = CG * VL_S_RY / chi;
    const Matrix5r ImPG_VL_S_RY      = (Matrix5r::Identity() - PG) * VL_S_RY;
    const Matrix5r inv_VL_S_RY       = inv_RY * inv_S * inv_VL;

    // Wavenumber traversal modeled after those found in suzerain/diffwave.c
    const int Ny   = dgrid.global_wave_extent.y();
    const int Nx   = grid.N.x();
    const int dNx  = grid.dN.x();
    const int dkbx = dgrid.local_wave_start.x();
    const int dkex = dgrid.local_wave_end.x();
    const int Nz   = grid.N.z();
    const int dNz  = grid.dN.z();
    const int dkbz = dgrid.local_wave_start.z();
    const int dkez = dgrid.local_wave_end.z();
    const real_t twopioverLx = twopiover(grid.L.x());  // Weird looking...
    const real_t twopioverLz = twopiover(grid.L.z());  // ...for FP control

    // Traverse wavenumbers updating the RHS with the Giles boundary condition
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
            // The imaginary unit has already been included within i_stash.
            Vector5c tmp;
            tmp.noalias()  = ImPG_VL_S_RY.cast<complex_t>() * N;
            tmp.noalias() += kn
                           * BG_VL_S_RY_by_chi.cast<complex_t>()
                           * i_stash.col((m - dkbx) + mu*(n - dkbz));
            tmp.noalias() += km
                           * CG_VL_S_RY_by_chi.cast<complex_t>()
                           * i_stash.col((m - dkbx) + mu*(n - dkbz));
            N.noalias()    = inv_VL_S_RY.cast<complex_t>() * tmp;

            // Pack new upper boundary RHS from the contiguous buffer
            for (int f = 0; f < N.size(); ++f) {
                swave[f][Ny - 1][m - dkbx][n - dkbz] = N(f);
            }

        }
    }

    return retval;
}

} // namespace perfect

} // namespace suzerain