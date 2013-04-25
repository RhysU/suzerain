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

    // Retrieve upper boundary reference state from the common block
    const real_t rho = common.ref_rho().tail<1>()[0];
    const real_t u   = common.ref_ux ().tail<1>()[0];
    const real_t v   = common.ref_uy ().tail<1>()[0];
    const real_t w   = common.ref_uz ().tail<1>()[0];
    const real_t a   = common.ref_a  ().tail<1>()[0];

    // Prepare oft-used derived quantities
    const real_t inv_rho = 1 / rho;
    const real_t u2      = u * u;
    const real_t v2      = v * v;
    const real_t w2      = w * w;
    const real_t a2      = a * a;

    // Prepare oft-used scenario-related constants
    const real_t gamma      = scenario.gamma;
    const real_t inv_gamma  = 1 / gamma;
    const real_t gamma1     = gamma - 1;
    const real_t inv_gamma1 = 1 / gamma1;
    const real_t Ma         = scenario.Ma;
    const real_t inv_Ma     = 1 / Ma;
    const real_t Ma2        = Ma * Ma;
    const real_t inv_Ma2    = 1 / Ma2;

    // Build the appropriate 5x5 projection matrices:
    //   1) The documentation lists these as rho, rho_u, rho_v, rho_w, rho_E.
    //   2) The code uses a different order, so logically access entries.
    //   3) Form piece parts for readability and then pre-compute products.
    Matrix5r RY = Matrix5r::Zero();
    {
        RY(ndx::e  , ndx::e  ) = 1;
        RY(ndx::mx , ndx::my ) = 1;
        RY(ndx::my , ndx::mz ) = 1;
        RY(ndx::mz , ndx::mx ) = 1;
        RY(ndx::rho, ndx::rho) = 1;
    }

    Matrix5r inv_RY = Matrix5r::Zero();
    {
        inv_RY(ndx::e  , ndx::e  ) = 1;
        inv_RY(ndx::mx , ndx::mz ) = 1;
        inv_RY(ndx::my , ndx::mx ) = 1;
        inv_RY(ndx::mz , ndx::my ) = 1;
        inv_RY(ndx::rho, ndx::rho) = 1;
    }

    Matrix5r S = Matrix5r::Zero();
    {
        S(ndx::e  , ndx::e )  =   gamma1 * inv_Ma2;
        S(ndx::e  , ndx::mx)  = - gamma1 * u;
        S(ndx::e  , ndx::my)  = - gamma1 * v;
        S(ndx::e  , ndx::mz)  = - gamma1 * w;
        S(ndx::e  , ndx::rho) =   a2 * inv_gamma * inv_Ma2;
        S(ndx::mx , ndx::mx)  =   inv_rho;
        S(ndx::mx , ndx::rho) = - u * inv_rho;
        S(ndx::my , ndx::my)  =   inv_rho;
        S(ndx::my , ndx::rho) = - v * inv_rho;
        S(ndx::mz , ndx::mz)  =   inv_rho;
        S(ndx::mz , ndx::rho) = - w * inv_rho;
        S(ndx::rho, ndx::rho) =   1;
    }

    Matrix5r inv_S  = Matrix5r::Zero();
    {
        S(ndx::e  , ndx::e )  =   Ma2 * inv_gamma1;
        S(ndx::e  , ndx::mx)  =   Ma2 * rho * u;
        S(ndx::e  , ndx::my)  =   Ma2 * rho * v;
        S(ndx::e  , ndx::mz)  =   Ma2 * rho * w;
        S(ndx::e  , ndx::rho) =   Ma2 * (u2 + v2 + w2)
                              +   a2 * inv_gamma * inv_gamma1;
        S(ndx::mx , ndx::mx)  =   rho;
        S(ndx::mx , ndx::rho) =   u;
        S(ndx::my , ndx::my)  =   rho;
        S(ndx::my , ndx::rho) =   v;
        S(ndx::mz , ndx::mz)  =   rho;
        S(ndx::mz , ndx::rho) =   w;
        S(ndx::rho, ndx::rho) =   1;
    }

    Matrix5r VL = Matrix5r::Zero(), inv_VL = Matrix5r::Zero();
    Matrix5r BG = Matrix5r::Zero();
    Matrix5r CG = Matrix5r::Zero();
    Matrix5r PG = Matrix5r::Zero();

    // TODO Implement
    return N->apply_operator(
            time, swave, evmaxmag_real, evmaxmag_imag, substep_index);

    // State leaves method as coefficients in X and Z directions
    // State leaves method as collocation point values in Y direction
}

} // namespace perfect

} // namespace suzerain
