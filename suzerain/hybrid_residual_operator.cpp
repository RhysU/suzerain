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
 * @copydoc hybrid_residual_operator.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/hybrid_residual_operator.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/state.hpp>

namespace suzerain {

hybrid_residual_operator::hybrid_residual_operator(
        const real_t chi)
    : chi(chi)
{
}

std::vector<real_t> hybrid_residual_operator::apply_operator(
            const real_t time,
            contiguous_state<4,complex_t> &state,
            const lowstorage::method_interface<element>& method,
            const std::size_t substep_index) const
{
    // Allocate (potentially large) extra working storage
    state_linear_type extra(shape_array(state));

    // The following steps are taken
    //     (1) extra <- state
    //     (2) state <- R(state)
    //     (3) state <- (M +      0  * L) extra - 1 * state
    //     (4) state <- (M - (1/chi) * L) extra - 1 * state
    // and they ultimately result in
    //         state <- R(u) - (1/chi) Lu = N(u)
    // There are lots of ways to make this process less computationally
    // expensive if it turns out to be problematic, but this procedure
    // works in the context of the current APIs.

    extra.assign_from(state);
    std::vector<real_t> retval = R->apply_operator(
            time, state, method, substep_index);
    L->accumulate_mass_plus_scaled_operator(
            0,      extra, -1, state, substep_index);
    L->accumulate_mass_plus_scaled_operator(
            -1/chi, extra, -1, state, substep_index);

    return retval;
}

} // namespace suzerain
