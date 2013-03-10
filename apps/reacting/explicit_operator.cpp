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
 * @copydoc explicit_operator.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include "explicit_operator.hpp"

#include <suzerain/common.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/state.hpp>

#include "nonlinear_operator.hpp"

#pragma warning(disable:383 1572)

namespace suzerain {

namespace reacting {

// A helper class for implementing isothermal, no-slip boundary conditions
class IsothermalNoSlipFunctor
{
private:
    const ptrdiff_t field_stride;
    const real_t    Cv;
    const real_t    T_wall;

public:
    IsothermalNoSlipFunctor(ptrdiff_t field_stride, real_t Cv, real_t T_wall)
        : field_stride(field_stride)
	, Cv(Cv)
	, T_wall(T_wall)
    {}

    void operator()(complex_t &rho) const
    {
        (&rho)[(ndx::mx - ndx::rho)*field_stride] = 0;
        (&rho)[(ndx::my - ndx::rho)*field_stride] = 0;
        (&rho)[(ndx::mz - ndx::rho)*field_stride] = 0;
        (&rho)[(ndx::e  - ndx::rho)*field_stride] = rho*Cv*T_wall; //rho*inv_gamma_gamma1;
    }
};

isothermal_mass_operator::isothermal_mass_operator(
        const single_ideal_gas_constitutive& cmods,
        const channel_definition &chdef,
        const grid_specification &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        bspline &b,
        operator_common_block &common)
    : mass_operator(grid, dgrid, cop, b)
    , cmods(cmods)
    , chdef(chdef)
    , common(common)
    , who("operator.L")
{
    // NOP
}

void isothermal_mass_operator::invert_mass_plus_scaled_operator(
        const complex_t &phi,
        multi_array::ref<complex_t,4> &state,
        const timestepper::lowstorage::method_interface<complex_t> &method,
        const component delta_t,
        const std::size_t substep_index,
        multi_array::ref<complex_t,4> *ic0) const
{
    // State enters method as coefficients in X and Z directions
    // State enters method as collocation point values in Y direction

    // Shorthand
    using boost::indices;
    typedef boost::multi_array_types::index_range range;

    // Indexes only the first and last collocation point
    const std::size_t Ny         = state.shape()[1];
    const std::size_t wall_lower = 0;
    const std::size_t wall_upper = Ny - 1;
    range walls(wall_lower, wall_upper + 1, wall_upper - wall_lower);

    // Prepare a state view of density locations at lower and upper walls
    multi_array::ref<complex_t,4>::array_view<3>::type state_view
            = state[indices[ndx::rho][walls][range()][range()]];

    // Prepare functor setting pointwise BCs given density locations
    const IsothermalNoSlipFunctor bc_functor(
	  state.strides()[0], cmods.Cv, chdef.T_wall);

    // Apply the functor to all wall-only density locations
    multi_array::for_each(state_view, bc_functor);

    // Apply boundary conditions to any requested constraint problems
    if (ic0) {
        multi_array::ref<complex_t,4>::array_view<3>::type ic0_view
                = (*ic0)[indices[ndx::rho][walls][range()][range()]];
        SUZERAIN_ENSURE(state.strides()[0] == ic0->strides()[0]); // NB!
        multi_array::for_each(ic0_view, bc_functor);
    }

    // channel_treatment step (3) performs the usual operator solve
    base::invert_mass_plus_scaled_operator(
            phi, state, method, delta_t, substep_index, ic0);

    // State leaves method as coefficients in X, Y, and Z directions
}

explicit_nonlinear_operator::explicit_nonlinear_operator(
        const single_ideal_gas_constitutive& cmods,
        const grid_specification &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        bspline &b,
        operator_common_block &common,
        const filter_definition &fsdef,
        const shared_ptr<const manufactured_solution>& msoln)
    : operator_base(grid, dgrid, cop, b)
    , cmods(cmods)
    , common(common)
    , fsdef(fsdef)
    , msoln(msoln)
    , who("operator.N")
{
    // NOP
}

std::vector<real_t> explicit_nonlinear_operator::apply_operator(
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag,
            const std::size_t substep_index) const
{
    // Dispatch to implementation paying nothing for substep-related ifs
    if (substep_index == 0) {
        return apply_navier_stokes_spatial_operator<true,  linearize::none>
            (*this, common, fsdef, msoln, cmods, time, swave, evmaxmag_real, evmaxmag_imag);
    } else {
        return apply_navier_stokes_spatial_operator<false, linearize::none>
            (*this, common, fsdef, msoln, cmods, time, swave, evmaxmag_real, evmaxmag_imag);
    }
}

} // namespace reacting

} // namespace suzerain
