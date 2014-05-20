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
 * @copydoc operator_nonlinear.hpp
 */

#include "operator_nonlinear.hpp"

#include <suzerain/common.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/specification_grid.hpp>
#include <suzerain/specification_largo.hpp>
#include <suzerain/state.hpp>

#include "definition_scenario.hpp"
#include "manufactured_solution.hpp"
#include "navier_stokes.hpp"
#include "operator_common_block.hpp"

#pragma warning(disable:383 1572)

namespace suzerain {

namespace perfect {

operator_nonlinear::operator_nonlinear(
        const definition_scenario &scenario,
        const specification_grid &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        bspline &b,
        operator_common_block &common,
        specification_largo& sg,
        const shared_ptr<const manufactured_solution>& msoln)
    : operator_base(grid, dgrid, cop, b)
    , scenario(scenario)
    , common(common)
    , sg(sg)
    , msoln(msoln)
    , who("operator.N")
{
}

std::vector<real_t> operator_nonlinear::apply_operator(
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const lowstorage::method_interface<complex_t> &method,
            const std::size_t substep_index) const
{
#define ARGUMENTS scenario, *this, common, sg, msoln, time, swave, \
                  method, substep_index

    // Dispatch to an optimized implementation depending on runtime settings.
    // Done since the compiler can theoretically hammer out much savings when
    // many small, localized jumps can be hoisted up to this level.  In
    // practice, with ZerothSubstep, explicitly this saved 1-2% of wall time.
    if (sg.formulation.enabled()) {

        switch (common.linearization) {
        case linearize::rhome_xyz:
            return (substep_index == 0)
                 ? apply_navier_stokes_spatial_operator<
                        true,  linearize::rhome_xyz, true>(ARGUMENTS)
                 : apply_navier_stokes_spatial_operator<
                        false, linearize::rhome_xyz, true>(ARGUMENTS);

        case linearize::rhome_y:
            return (substep_index == 0)
                 ? apply_navier_stokes_spatial_operator<
                        true,  linearize::rhome_y, true>(ARGUMENTS)
                 : apply_navier_stokes_spatial_operator<
                        false, linearize::rhome_y, true>(ARGUMENTS);

        case linearize::none:
            return (substep_index == 0)
                 ? apply_navier_stokes_spatial_operator<
                        true,  linearize::none, true>(ARGUMENTS)
                 : apply_navier_stokes_spatial_operator<
                        false, linearize::none, true>(ARGUMENTS);

        default:
            SUZERAIN_ERROR_VAL_UNIMPLEMENTED(std::vector<real_t>());
            break;
        }

    } else {

        switch (common.linearization) {
        case linearize::rhome_xyz:
            return (substep_index == 0)
                 ? apply_navier_stokes_spatial_operator<
                        true,  linearize::rhome_xyz, false>(ARGUMENTS)
                 : apply_navier_stokes_spatial_operator<
                        false, linearize::rhome_xyz, false>(ARGUMENTS);

        case linearize::rhome_y:
            return (substep_index == 0)
                 ? apply_navier_stokes_spatial_operator<
                        true,  linearize::rhome_y, false>(ARGUMENTS)
                 : apply_navier_stokes_spatial_operator<
                        false, linearize::rhome_y, false>(ARGUMENTS);

        case linearize::none:
            return (substep_index == 0)
                 ? apply_navier_stokes_spatial_operator<
                        true,  linearize::none, false>(ARGUMENTS)
                 : apply_navier_stokes_spatial_operator<
                        false, linearize::none, false>(ARGUMENTS);

        default:
            SUZERAIN_ERROR_VAL_UNIMPLEMENTED(std::vector<real_t>());
            break;
        }

    }

#undef ARGUMENTS
}

} // namespace perfect

} // namespace suzerain
