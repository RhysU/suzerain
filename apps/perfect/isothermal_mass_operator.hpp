//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012, 2013 Rhys Ulerich
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

#ifndef SUZERAIN_PERFECT_ISOTHERMAL_MASS_OPERATOR_HPP
#define SUZERAIN_PERFECT_ISOTHERMAL_MASS_OPERATOR_HPP

/** @file
 * Fully explicit, linearization-ready Navier--Stokes operators.
 */

#include <suzerain/mass_operator.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/state_fwd.hpp>

#pragma warning(disable:383 1572)

namespace suzerain {

// Forward declarations
class grid_specification;
class pencil_grid;

namespace perfect {

// Forward declarations
class operator_common_block;
class scenario_definition;

/**
 * A mass operator that provides no slip, isothermal walls.  It requires
 * interoperation with nonlinear_operator via operator_common_block.
 */
class isothermal_mass_operator : public mass_operator
{

    typedef mass_operator base;

public:

    isothermal_mass_operator(
            const scenario_definition &scenario,
            const grid_specification &grid,
            const pencil_grid &dgrid,
            const bsplineop &cop,
            bspline &b,
            operator_common_block &common);

    /**
     * Performs the following steps:
     * <ul>
     * <li>
     *     channel_treatment step (3) performs the operator solve which for
     *     the implicit treatment must be combined with boundary conditions
     * </li><li>
     *     channel_treatment step (8) sets no-slip conditions
     *     on wall collocation points.
     * </li><li>
     * </li>
     *     channel_treatment step (9) sets isothermal conditions at walls
     *     using rho_wall = e_wall * gamma * (gamma - 1).
     * </ul>
     */
    virtual void invert_mass_plus_scaled_operator(
            const complex_t &phi,
            multi_array::ref<complex_t,4> &state,
            const timestepper::lowstorage::method_interface<complex_t> &method,
            const component delta_t,
            const std::size_t substep_index,
            multi_array::ref<complex_t,4> *ic0 = NULL) const;

protected:

    /** The scenario in which the operator is used */
    const scenario_definition &scenario;

    /** Houses data required for \ref invert_mass_plus_scaled_operator */
    operator_common_block &common;

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;

};

} // namespace perfect

} // namespace suzerain

#endif  /* SUZERAIN_PERFECT_ISOTHERMAL_MASS_OPERATOR_HPP */
