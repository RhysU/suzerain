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
 * @copydoc operator_mass.hpp
 */

#include <suzerain/operator_mass.hpp>

#include <suzerain/pencil_grid.hpp>
#include <suzerain/specification_grid.hpp>
#include <suzerain/state.hpp>

namespace suzerain {

operator_mass::operator_mass(
        const specification_grid &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        bspline &b)
    : operator_base(grid, dgrid, cop, b)
{
    SUZERAIN_UNUSED(grid);
    SUZERAIN_UNUSED(dgrid);
    SUZERAIN_UNUSED(b);

    // Ensure cached mass matrix factorized prior to operator application
    // Strictly speaking unnecessary, but reduces timing variability
    this->massluz();
}

operator_mass::~operator_mass()
{
}

void operator_mass::apply_mass_plus_scaled_operator(
        const complex_t &phi,
        multi_array::ref<complex_t,4> &state,
        const std::size_t substep_index) const
{
    SUZERAIN_UNUSED(phi);
    SUZERAIN_UNUSED(substep_index);

    // Verify required assumptions
    SUZERAIN_ENSURE(state.strides()[1] == 1);
    SUZERAIN_ENSURE(state.shape()[1]   == static_cast<unsigned>(cop.n()));
    SUZERAIN_ENSURE(multi_array::is_contiguous(state));

    // Those assumptions holding, apply operator to each wall-normal pencil.
    // Speedup may be possible by not applying to dealiased modes
    const int nrhs = state.shape()[0]*state.shape()[2]*state.shape()[3];
    cop.apply(0, nrhs, 1, state.data(), 1, state.shape()[1]);
}

void operator_mass::accumulate_mass_plus_scaled_operator(
        const complex_t &phi,
        const multi_array::ref<complex_t,4> &input,
        const complex_t &beta,
        contiguous_state<4,complex_t> &output,
        const std::size_t substep_index) const
{
    SUZERAIN_UNUSED(phi);
    SUZERAIN_UNUSED(substep_index);

    SUZERAIN_ENSURE(output.is_isomorphic(input));

    const multi_array::ref<complex_t,4> &x = input;  // Shorthand
    contiguous_state<4,complex_t>       &y = output; // Shorthand
    const complex_t c_one = 1;

    // Sidesteps assertions triggered by dereferencing trivial input and output
    if (SUZERAIN_UNLIKELY(0U == x.shape()[1] * x.shape()[2])) return;

    // Loops go from slower to faster indices for contiguous_state<4,complex_t>
    // Speedup may be possible by not accumulating dealiased modes
    typedef contiguous_state<4,complex_t>::index index;
    for (index ix = x.index_bases()[0], iy = y.index_bases()[0];
        ix < static_cast<index>(x.index_bases()[0] + x.shape()[0]);
        ++ix, ++iy) {

        for (index lx = x.index_bases()[3], ly = y.index_bases()[3];
            lx < static_cast<index>(x.index_bases()[3] + x.shape()[3]);
            ++lx, ++ly) {

            cop.accumulate(0, x.shape()[2],
                    c_one,
                    &x[ix][x.index_bases()[1]][x.index_bases()[2]][lx],
                    x.strides()[1], x.strides()[2],
                    beta,
                    &y[iy][y.index_bases()[1]][y.index_bases()[2]][ly],
                    y.strides()[1], y.strides()[2]);
        }
    }
}

void operator_mass::invert_mass_plus_scaled_operator(
        const complex_t &phi,
        multi_array::ref<complex_t,4> &state,
        const lowstorage::method_interface<complex_t> &method,
        const component delta_t,
        const std::size_t substep_index,
        multi_array::ref<complex_t,4> *ic0) const
{
    SUZERAIN_UNUSED(phi);
    SUZERAIN_UNUSED(method);
    SUZERAIN_UNUSED(delta_t);
    SUZERAIN_UNUSED(substep_index);

    // Verify required assumptions
    SUZERAIN_ENSURE(state.strides()[1] == 1);
    SUZERAIN_ENSURE(state.shape()[1]   == static_cast<unsigned>(cop.n()));
    SUZERAIN_ENSURE(multi_array::is_contiguous(state));

    // Those assumptions holding, invert operator on each wall-normal pencil
    massluz()->solve(state.shape()[0]*state.shape()[2]*state.shape()[3],
                     state.data(), 1, state.shape()[1]);

    if (ic0) {

        // Likewise, verify required assumptions
        SUZERAIN_ENSURE(ic0->strides()[1] == 1);
        SUZERAIN_ENSURE(ic0->shape()[1]   == static_cast<unsigned>(cop.n()));
        SUZERAIN_ENSURE(multi_array::is_contiguous(*ic0));

        // Likewise, invert operator on each wall-normal pencil
        massluz()->solve(ic0->shape()[0]*ic0->shape()[2]*ic0->shape()[3],
                         ic0->data(), 1, ic0->shape()[1]);

    }
}

} // namespace suzerain
