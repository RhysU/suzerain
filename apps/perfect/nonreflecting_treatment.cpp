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
#include <suzerain/ndx.hpp>
#include <suzerain/pencil_grid.hpp>

#include "common_block.hpp"
#include "scenario_definition.hpp"

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
    // State enters method as coefficients in X, Y, and Z directions

    // Prepare the rotation and its inverse that reorders from
    // ndx::{e, mx, my, mz, rho} to {rho = 0, my = 1, mz = 2, mx = 3, e = 4}.
    // All remaining matrices/logic within the routine uses the latter order!
    Matrix5r RY = Matrix5r::Zero();
    {
        RY(0, ndx::rho) = 1;
        RY(1, ndx::my ) = 1;
        RY(2, ndx::mz ) = 1;
        RY(3, ndx::mx ) = 1;
        RY(4, ndx::e  ) = 1;
    }
    Matrix5r inv_RY = Matrix5r::Zero();
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

    // Prepare oft-used derived quantities
    const real_t inv_rho = 1 / rho;
    const real_t u2      = u * u;
    const real_t v2      = v * v;
    const real_t w2      = w * w;
          real_t inv_a   = 1 / a;
          real_t a2      = a * a;
          real_t inv_a2  = 1 / a2;

    // Prepare oft-used scenario-related constants
    const real_t half       = real_t(1) / 2;
    const real_t gamma      = scenario.gamma;
    const real_t inv_gamma  = 1 / gamma;
    const real_t gamma1     = gamma - 1;
    const real_t inv_gamma1 = 1 / gamma1;
    const real_t Ma         = scenario.Ma;
    const real_t inv_Ma     = 1 / Ma;
    const real_t Ma2        = Ma * Ma;
    const real_t inv_Ma2    = 1 / Ma2;

    // Build the variable transformation matrices
    Matrix5r S = Matrix5r::Zero();
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

    Matrix5r inv_S  = Matrix5r::Zero();
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

    Matrix5r VL = Matrix5r::Zero();
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

    Matrix5r inv_VL = Matrix5r::Zero();
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
    // In the event of v == 0, also treat the boundary like an inflow.
    const bool inflow = u <= 0;

    Matrix5r PG = Matrix5r::Zero();  // Characteristic-preserving projection
    if (inflow) {
        // TODO
    } else {
        // TODO
    }

    Matrix5r BG = Matrix5r::Zero();  // Medida's B^G_1
    if (inflow) {
        // TODO
    } else {
        // TODO
    }

    Matrix5r CG = Matrix5r::Zero();  // Medida's C^G_1
    if (inflow) {
        // TODO
    } else {
        // TODO
    }

    // TODO Implement
    return N->apply_operator(
            time, swave, evmaxmag_real, evmaxmag_imag, substep_index);

    // State leaves method as coefficients in X and Z directions
    // State leaves method as collocation point values in Y direction
}

} // namespace perfect

} // namespace suzerain
