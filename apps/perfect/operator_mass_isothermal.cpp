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
 * @copydoc operator_mass_isothermal.hpp
 */

#include "operator_mass_isothermal.hpp"

#include <suzerain/common.hpp>

#include "definition_scenario.hpp"
#include "operator_common_block.hpp"

#pragma warning(disable:383 1572)

namespace suzerain {

namespace perfect {

operator_mass_isothermal::operator_mass_isothermal(
        const definition_scenario &scenario,
        const specification_isothermal &spec,
        const specification_grid &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        bspline &b,
        operator_common_block &common)
    : suzerain::operator_mass_isothermal(spec, grid, dgrid, cop, b)
    , scenario(scenario)
    , common(common)
    , who("operator.L")
{
}

void operator_mass_isothermal::invert_mass_plus_scaled_operator(
        const complex_t &phi,
        multi_array::ref<complex_t,4> &state,
        const lowstorage::method_interface<complex_t> &method,
        const component delta_t,
        const std::size_t substep_index,
        multi_array::ref<complex_t,4> *ic0) const
{
    SUZERAIN_ENSURE(common.linearization == linearize::none);
    return suzerain::operator_mass_isothermal::invert_mass_plus_scaled_operator(
            phi, state, method, delta_t, substep_index, ic0);
}

real_t operator_mass_isothermal::lower_E(
        const real_t lower_T,
        const real_t lower_u,
        const real_t lower_v,
        const real_t lower_w,
        const std::vector<real_t> lower_cs) const
{
    // FIXME Gather and use instantaneous averages when isnan(lower_[uvw]).
    SUZERAIN_ENSURE(lower_cs.size() == 1U && lower_cs[0] == 1.0);
    const real_t E_internal = lower_T / (scenario.gamma*(scenario.gamma - 1));
    const real_t E_kinetic  = (scenario.Ma*scenario.Ma)
                            * (  lower_u*lower_u
                               + lower_v*lower_v
                               + lower_w*lower_w ) / 2;
    return E_internal + E_kinetic;
}

real_t operator_mass_isothermal::upper_E(
        const real_t upper_T,
        const real_t upper_u,
        const real_t upper_v,
        const real_t upper_w,
        const std::vector<real_t> upper_cs) const
{
    // FIXME Gather and use instantaneous averages when isnan(lower_[uvw]).
    return lower_E(upper_T, upper_u, upper_v, upper_w, upper_cs);
}

} // namespace perfect

} // namespace suzerain
