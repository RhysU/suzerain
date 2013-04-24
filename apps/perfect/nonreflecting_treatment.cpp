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

    // TODO Implement
    return N->apply_operator(
            time, swave, evmaxmag_real, evmaxmag_imag, substep_index);

    // State leaves method as coefficients in X and Z directions
    // State leaves method as collocation point values in Y direction
}

} // namespace perfect

} // namespace suzerain
