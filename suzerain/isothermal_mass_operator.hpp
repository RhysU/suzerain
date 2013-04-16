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

#include <suzerain/common.hpp>
#include <suzerain/isothermal_specification.hpp>
#include <suzerain/mass_operator.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/state_fwd.hpp>

namespace suzerain {

// Forward declarations
class grid_specification;
class pencil_grid;

/**
 * An abstract mass operator that provides isothermal boundary conditions for
 * equations ordered per \ref suzerain::ndx.  Subclasses must provide a way to
 * compute specific total energy given local state information.  In this way,
 * different equations of state may be accommodated.
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

    /**
     * Compute the specific total energy \f$E\f$ on the lower boundary.  The
     * specified condition is provided but subclasses should account for any
     * relevant auxiliary information.  For example, when enforcing \c lower_v
     * is disabled using \c NaN per \ref isothermal_specification#lower_v for a
     * non-reflecting boundary condition, the instantaneous mean velocity
     * should instead be used when computing \f$E\f$.
     *
     * @param lower_T  Specified per \ref isothermal_specification#lower_T.
     * @param lower_u  Specified per \ref isothermal_specification#lower_u.
     * @param lower_v  Specified per \ref isothermal_specification#lower_v.
     * @param lower_w  Specified per \ref isothermal_specification#lower_w.
     * @param lower_cs Specified per \ref isothermal_specification#lower_cs.
     *
     * @return the value of \f$E\f$ to use when adjusting the lower boundary.
     */
    virtual real_t lower_E(const real_t lower_T,
                           const real_t lower_u,
                           const real_t lower_v,
                           const real_t lower_w,
                           const std::vector<real_t> lower_cs) const = 0;

    /**
     * Compute the specific total energy \f$E\f$ on the upper boundary.
     * See \ref lower_E for additional comments.
     *
     * @param upper_T  Specified per \ref isothermal_specification#upper_T.
     * @param upper_u  Specified per \ref isothermal_specification#upper_u.
     * @param upper_v  Specified per \ref isothermal_specification#upper_v.
     * @param upper_w  Specified per \ref isothermal_specification#upper_w.
     * @param upper_cs Specified per \ref isothermal_specification#upper_cs.
     *
     * @return the value of \f$E\f$ to use when adjusting the upper boundary.
     */
    virtual real_t upper_E(const real_t upper_T,
                           const real_t upper_u,
                           const real_t upper_v,
                           const real_t upper_w,
                           const std::vector<real_t> upper_cs) const = 0;

    virtual void invert_mass_plus_scaled_operator(
            const complex_t &phi,
            multi_array::ref<complex_t,4> &state,
            const timestepper::lowstorage::method_interface<complex_t> &method,
            const component delta_t,
            const std::size_t substep_index,
            multi_array::ref<complex_t,4> *ic0 = NULL) const;

private:

    /** The lower and upper isothermal boundary specification. */
    const isothermal_specification spec;

};

} // namespace suzerain

#endif  /* SUZERAIN_ISOTHERMAL_MASS_OPERATOR_HPP */
