//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
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
 * @copydoc nonlinear_operator.hpp
 */

#include "nonlinear_operator.hpp"

#include <suzerain/common.hpp>
#include <suzerain/grid_specification.hpp>
#include <suzerain/largo_specification.hpp>
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
        largo_specification& sg,
        const shared_ptr<const manufactured_solution>& msoln)
    : operator_base(grid, dgrid, cop, b)
    , scenario(scenario)
    , common(common)
    , sg(sg)
    , msoln(msoln)
    , who("operator.N")
{
}

std::vector<real_t> nonlinear_operator::apply_operator(
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const lowstorage::method_interface<complex_t> &method,
            const std::size_t substep_index) const
{
#define ARGUMENTS this->scenario.alpha,                         \
                  this->scenario.beta,                          \
                  this->scenario.gamma,                         \
                  this->scenario.Ma,                            \
                  this->scenario.Pr,                            \
                  this->scenario.Re,                            \
                  *this, common, sg, msoln, time, swave,        \
                  method, substep_index

    // Dispatch to an optimized implementation depending on runtime settings.
    // Done since the compiler can theoretically hammer out much savings when
    // many small, localized jumps can be hoisted up to this level.  In
    // practice, with ZerothSubstep, this has saved 1-2% of wall time.
    switch (common.slow_treatment) {

    case slowgrowth::none:
        switch (common.linearization) {
        case linearize::rhome_xyz:
            return (substep_index == 0)
                 ? apply_navier_stokes_spatial_operator<true,
                        linearize::rhome_xyz, slowgrowth::none>(ARGUMENTS)
                 : apply_navier_stokes_spatial_operator<false,
                        linearize::rhome_xyz, slowgrowth::none>(ARGUMENTS);

        case linearize::rhome_y:
            return (substep_index == 0)
                 ? apply_navier_stokes_spatial_operator<true,
                        linearize::rhome_y, slowgrowth::none>(ARGUMENTS)
                 : apply_navier_stokes_spatial_operator<false,
                        linearize::rhome_y, slowgrowth::none>(ARGUMENTS);

        case linearize::none:
            return (substep_index == 0)
                 ? apply_navier_stokes_spatial_operator<true,
                        linearize::none, slowgrowth::none>(ARGUMENTS)
                 : apply_navier_stokes_spatial_operator<false,
                        linearize::none, slowgrowth::none>(ARGUMENTS);

        default:
            SUZERAIN_ERROR_VAL_UNIMPLEMENTED(std::vector<real_t>());
            break;
        }
        break;

    case slowgrowth::largo:
        switch (common.linearization) {
        case linearize::rhome_xyz:
            return (substep_index == 0)
                 ? apply_navier_stokes_spatial_operator<true,
                        linearize::rhome_xyz, slowgrowth::largo>(ARGUMENTS)
                 : apply_navier_stokes_spatial_operator<false,
                        linearize::rhome_xyz, slowgrowth::largo>(ARGUMENTS);

        case linearize::rhome_y:
            return (substep_index == 0)
                 ? apply_navier_stokes_spatial_operator<true,
                        linearize::rhome_y, slowgrowth::largo>(ARGUMENTS)
                 : apply_navier_stokes_spatial_operator<false,
                        linearize::rhome_y, slowgrowth::largo>(ARGUMENTS);

        case linearize::none:
            return (substep_index == 0)
                 ? apply_navier_stokes_spatial_operator<true,
                        linearize::none, slowgrowth::largo>(ARGUMENTS)
                 : apply_navier_stokes_spatial_operator<false,
                        linearize::none, slowgrowth::largo>(ARGUMENTS);

        default:
            SUZERAIN_ERROR_VAL_UNIMPLEMENTED(std::vector<real_t>());
            break;
        }
        break;

    default:
        SUZERAIN_ERROR_VAL_UNIMPLEMENTED(std::vector<real_t>());
        break;

    }

#undef ARGUMENTS
}

} // namespace perfect

} // namespace suzerain
