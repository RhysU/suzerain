//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
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

#ifndef SUZERAIN_PERFECT_OPERATOR_MASS_ISOTHERMAL_HPP
#define SUZERAIN_PERFECT_OPERATOR_MASS_ISOTHERMAL_HPP

/** @file
 * Provides \ref operator_mass_isothermal
 */

#include <suzerain/operator_mass_isothermal.hpp>

#pragma warning(disable:383 1572)

namespace suzerain {

// Forward declarations
class specification_grid;
class specification_isothermal;
class pencil_grid;

namespace perfect {

// Forward declarations
class operator_common_block;
class definition_scenario;

/**
 * A mass operator that provides various isothermal boundaries.  It requires
 * interoperation with operator_nonlinear via operator_common_block.
 */
class operator_mass_isothermal : public suzerain::operator_mass_isothermal
{

public:

    operator_mass_isothermal(
            const definition_scenario &scenario,
            const specification_isothermal &spec,
            const specification_grid &grid,
            const pencil_grid &dgrid,
            const bsplineop &cop,
            bspline &b,
            operator_common_block &common);

    virtual void invert_mass_plus_scaled_operator(
            const complex_t &phi,
            multi_array::ref<complex_t,4> &state,
            const lowstorage::method_interface<complex_t> &method,
            const component delta_t,
            const std::size_t substep_index,
            multi_array::ref<complex_t,4> *ic0 = NULL) const;

    virtual real_t lower_E(const real_t lower_T,
                           const real_t lower_u,
                           const real_t lower_v,
                           const real_t lower_w,
                           const std::vector<real_t> lower_cs) const;

    virtual real_t upper_E(const real_t upper_T,
                           const real_t upper_u,
                           const real_t upper_v,
                           const real_t upper_w,
                           const std::vector<real_t> upper_cs) const;
protected:

    /** The scenario in which the operator is used */
    const definition_scenario &scenario;

    /** Houses data required for \ref invert_mass_plus_scaled_operator */
    operator_common_block &common;

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;

};

} // namespace perfect

} // namespace suzerain

#endif  /* SUZERAIN_PERFECT_OPERATOR_MASS_ISOTHERMAL_HPP */
