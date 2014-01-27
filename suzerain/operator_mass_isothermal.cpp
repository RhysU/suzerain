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

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/operator_mass_isothermal.hpp>

#include <suzerain/inorder.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/specification_grid.hpp>
#include <suzerain/state.hpp>

namespace suzerain {

operator_mass_isothermal::operator_mass_isothermal(
        const specification_isothermal &spec,
        const specification_grid &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        bspline &b)
    : operator_mass(grid, dgrid, cop, b)
    , spec(spec)
{
}

/**
 * A helper functor for implementing isothermal boundary conditions
 * per \ref specification_isothermal.
 */
class isothermal_enforcer
    : public boost::noncopyable // To defend against expense in for_each
{

    const std::size_t               Ny;       ///< # of points across y
    const std::ptrdiff_t            incf;     ///< Stride between fields
    const specification_isothermal& spec;     ///< Boundary specification
    const real_t                    lower_E;  ///< Specific total energy
    const real_t                    upper_E;  ///< Specific total energy

    /**
     * Bit fields marking exactly which BCs should be enforced.
     * @{*/
    const int enforce_lower_e  : 1;
    const int enforce_upper_e  : 1;
    const int enforce_lower_u  : 1;
    const int enforce_upper_u  : 1;
    const int enforce_lower_v  : 1;
    const int enforce_upper_v  : 1;
    const int enforce_lower_w  : 1;
    const int enforce_upper_w  : 1;
    const int enforce_lower_cs : 1;
    const int enforce_upper_cs : 1;

    /**@}*/

public:

    /** Prepare an instance for the given strides and parameters. */
    isothermal_enforcer(const specification_grid&       grid,
                        const std::ptrdiff_t            incf,
                        const specification_isothermal& spec,
                        const real_t                    lower_E,
                        const real_t                    upper_E)
        : Ny(grid.N.y())
        , incf(incf)
        , spec(spec)
        , lower_E(lower_E)
        , upper_E(upper_E)
        , enforce_lower_e (!(boost::math::isnan)(lower_E     ))
        , enforce_upper_e (!(boost::math::isnan)(upper_E     )
                           && grid.two_sided())
        , enforce_lower_u (!(boost::math::isnan)(spec.lower_u))
        , enforce_upper_u (!(boost::math::isnan)(spec.upper_u)
                           && grid.two_sided())
        , enforce_lower_v (!(boost::math::isnan)(spec.lower_v))
        , enforce_upper_v (!(boost::math::isnan)(spec.upper_v)
                           && grid.two_sided())
        , enforce_lower_w (!(boost::math::isnan)(spec.lower_w))
        , enforce_upper_w (!(boost::math::isnan)(spec.upper_w)
                           && grid.two_sided())
        , enforce_lower_cs(true)
        , enforce_upper_cs(grid.two_sided())
    {}

    /**
     * Apply the conditions to wave-space mode right hand sides.
     * Application order follows the expected state sequencing in memory.
     */
    void operator()(complex_t& lower_rho) const
    {
        // Locate the upper density relative to the provided lower density
        complex_t& upper_rho = (&lower_rho)[Ny - 1];

#define LOWER(var) ((&lower_rho)[((var) - ndx::rho)*incf])
#define UPPER(var) ((&upper_rho)[((var) - ndx::rho)*incf])

        // Set total energy to be density times specific total energy
        if (enforce_lower_e) LOWER(ndx::e) = lower_rho * lower_E;
        if (enforce_upper_e) UPPER(ndx::e) = upper_rho * upper_E;

        // Set streamwise momentum to be density times prescribed velocity
        if (enforce_lower_u) LOWER(ndx::mx) = lower_rho * spec.lower_u;
        if (enforce_upper_u) UPPER(ndx::mx) = upper_rho * spec.upper_u;

        // Set wall-normal momentum to be density times prescribed velocity
        if (enforce_lower_v) LOWER(ndx::my) = lower_rho * spec.lower_v;
        if (enforce_upper_v) UPPER(ndx::my) = upper_rho * spec.upper_v;

        // Set spanwise momentum to be density times prescribed velocity
        if (enforce_lower_w) LOWER(ndx::mz) = lower_rho * spec.lower_w;
        if (enforce_upper_w) UPPER(ndx::mz) = upper_rho * spec.upper_w;

        // Do nothing to the density equation

        // Set species partial densities to be density times mixture fraction.
        if (enforce_lower_cs & enforce_upper_cs) {
            assert(spec.lower_cs.size() == spec.upper_cs.size());
            const std::size_t num_species = spec.lower_cs.size();
            for (std::size_t s = 1; s < num_species; ++s) {
                LOWER(ndx::rho_(s)) = lower_rho * spec.lower_cs[s];
                UPPER(ndx::rho_(s)) = upper_rho * spec.upper_cs[s];
            }
        } else if (enforce_lower_cs) {
            const std::size_t num_species = spec.lower_cs.size();
            for (std::size_t s = 1; s < num_species; ++s) {
                LOWER(ndx::rho_(s)) = lower_rho * spec.lower_cs[s];
            }
        } else if (enforce_upper_cs) {
            const std::size_t num_species = spec.upper_cs.size();
            for (std::size_t s = 1; s < num_species; ++s) {
                UPPER(ndx::rho_(s)) = upper_rho * spec.upper_cs[s];
            }
        } else {
            // NOP
        }

#undef LOWER
#undef UPPER

    }

};

void operator_mass_isothermal::invert_mass_plus_scaled_operator(
        const complex_t &phi,
        multi_array::ref<complex_t,4> &state,
        const lowstorage::method_interface<complex_t> &method,
        const component delta_t,
        const std::size_t substep_index,
        multi_array::ref<complex_t,4> *ic0) const
{
    // State enters method as coefficients in X and Z directions
    // State enters method as collocation point values in Y direction

    // Prepare functor setting pointwise BCs given lower density locations.
    // Applies these to BOTH lower and upper boundaries given only lower one!
    SUZERAIN_ENSURE(state.shape()[1] == (unsigned) grid.N.y());
    const isothermal_enforcer enforcer(this->grid,
                                       state.strides()[0], // field_stride
                                       this->spec,
                                       this->lower_E(spec.lower_T,
                                                     spec.lower_u,
                                                     spec.lower_v,
                                                     spec.lower_w,
                                                     spec.lower_cs),
                                       this->upper_E(spec.upper_T,
                                                     spec.upper_u,
                                                     spec.upper_v,
                                                     spec.upper_w,
                                                     spec.upper_cs));

    // Apply the enforcer to the lower boundary densities.
    typedef boost::multi_array_types::index_range range;
    multi_array::ref<complex_t,4>::array_view<2>::type state_view
            = state[boost::indices[ndx::rho][0][range()][range()]];
    multi_array::for_each(state_view,
                          boost::bind<void>(boost::cref(enforcer), _1));

    // Apply zero-zero mode boundary conditions to any requested constraints.
    if (ic0) {
        SUZERAIN_ENSURE(ic0->shape()[0] == state.shape()[0]);
        SUZERAIN_ENSURE(ic0->shape()[1] == state.shape()[1]);
        multi_array::ref<complex_t,4>::array_view<2>::type ic0_view
                = (*ic0)[boost::indices[ndx::rho][0][range()][range()]];
        multi_array::for_each(ic0_view,
                              boost::bind<void>(boost::cref(enforcer), _1));
    }

    // Perform the usual operator_mass solve across all equations
    base::invert_mass_plus_scaled_operator(
            phi, state, method, delta_t, substep_index, ic0);

    // State leaves method as coefficients in X, Y, and Z directions
}

} // namespace suzerain
