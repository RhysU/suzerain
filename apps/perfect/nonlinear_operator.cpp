//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
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
 * @copydoc nonlinear_operator.hpp
 */

#include "nonlinear_operator.hpp"

#include <suzerain/common.hpp>
#include <suzerain/grid_specification.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/state.hpp>

#include "common_block.hpp"
#include "manufactured_solution.hpp"
#include "navier_stokes.hpp"
#include "scenario_definition.hpp"

#pragma warning(disable:383 1572)

namespace suzerain {

namespace perfect {

nonlinear_operator::nonlinear_operator(
        const scenario_definition &scenario,
        const grid_specification &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        bspline &b,
        operator_common_block &common,
        const shared_ptr<const manufactured_solution>& msoln)
    : operator_base(grid, dgrid, cop, b)
    , scenario(scenario)
    , common(common)
    , msoln(msoln)
    , who("operator.N")
{
    // NOP
}

std::vector<real_t> nonlinear_operator::apply_operator(
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const lowstorage::method_interface<complex_t> &method,
            const std::size_t substep_index) const
{
#define ARGUMENTS this->scenario.alpha,                     \
                  this->scenario.beta,                      \
                  this->scenario.gamma,                     \
                  this->scenario.Ma,                        \
                  this->scenario.Pr,                        \
                  this->scenario.Re,                        \
                  *this, common, msoln, time, swave, method

    // Dispatch to an optimized right hand side per substep/linearization:
    switch (common.linearization) {
    case linearize::rhome_xyz:
        return (substep_index == 0)
             ? apply_navier_stokes_spatial_operator<true,
                    linearize::rhome_xyz>(ARGUMENTS)
             : apply_navier_stokes_spatial_operator<false,
                    linearize::rhome_xyz>(ARGUMENTS);

    case linearize::rhome_y:
        return (substep_index == 0)
             ? apply_navier_stokes_spatial_operator<true,
                    linearize::rhome_y>(ARGUMENTS)
             : apply_navier_stokes_spatial_operator<false,
                    linearize::rhome_y>(ARGUMENTS);

    case linearize::none:
        return (substep_index == 0)
             ? apply_navier_stokes_spatial_operator<true,
                    linearize::none>(ARGUMENTS)
             : apply_navier_stokes_spatial_operator<false,
                    linearize::none>(ARGUMENTS);

    default:
        SUZERAIN_ERROR_VAL_UNIMPLEMENTED(std::vector<real_t>());
        break;
    }

#undef ARGUMENTS
}

} // namespace perfect

} // namespace suzerain
