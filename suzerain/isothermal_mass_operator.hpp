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

#ifndef SUZERAIN_ISOTHERMAL_MASS_OPERATOR_HPP
#define SUZERAIN_ISOTHERMAL_MASS_OPERATOR_HPP

/** @file
 * Provides \ref isothermal_mass_operator.
 */

#include <suzerain/mass_operator.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/state_fwd.hpp>

namespace suzerain {

// Forward declarations
class grid_specification;
class isothermal_specification;
class pencil_grid;

/**
 * An abstract mass operator that provides isothermal boundary conditions.
 * Subclasses must provide a way to compute specific total energy given local
 * state information.  In this way, different equations of state may be
 * accommodated.
 */
class isothermal_mass_operator : public mass_operator
{

    typedef mass_operator base;

public:

    /** Construct an instance using the given specification. */
    isothermal_mass_operator(
            const isothermal_specification &spec,
            const grid_specification &grid,
            const pencil_grid &dgrid,
            const bsplineop &cop,
            bspline &b);

////FIXME Draft
////virtual void invert_mass_plus_scaled_operator(
////        const complex_t &phi,
////        multi_array::ref<complex_t,4> &state,
////        const timestepper::lowstorage::method_interface<complex_t> &method,
////        const component delta_t,
////        const std::size_t substep_index,
////        multi_array::ref<complex_t,4> *ic0 = NULL) const;

protected:

    /** The lower and upper isothermal boundary specification. */
    const isothermal_specification &spec;

};

} // namespace suzerain

#endif  /* SUZERAIN_ISOTHERMAL_MASS_OPERATOR_HPP */
