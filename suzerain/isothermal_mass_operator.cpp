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
 * @copydoc isothermal_mass_operator.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/isothermal_mass_operator.hpp>

#include <suzerain/grid_specification.hpp>
#include <suzerain/inorder.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/state.hpp>

namespace suzerain {

isothermal_mass_operator::isothermal_mass_operator(
        const isothermal_specification &spec,
        const grid_specification &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        bspline &b)
    : mass_operator(grid, dgrid, cop, b)
    , spec(spec)
{
    // NOP
}

/**
 * A helper functor for implementing isothermal boundary conditions
 * per \ref isothermal_specification.
 */
class isothermal_enforcer
    : public boost::noncopyable // To defend against expense in for_each
{

    /** Bitmasks used to precompute and cache which conditions to #flags. */
    enum {
        ENFORCE_LOWER_E  = 1 << 0,
        ENFORCE_LOWER_U  = 1 << 1,
        ENFORCE_LOWER_V  = 1 << 2,
        ENFORCE_LOWER_W  = 1 << 3,
        ENFORCE_LOWER_CS = 1 << 4,
        ENFORCE_UPPER_E  = 1 << 5,
        ENFORCE_UPPER_U  = 1 << 6,
        ENFORCE_UPPER_V  = 1 << 7,
        ENFORCE_UPPER_W  = 1 << 8,
        ENFORCE_UPPER_CS = 1 << 9,
    };

    const std::size_t               Ny;       ///< # of points across y
    const std::ptrdiff_t            incf;     ///< Stride between fields
    const isothermal_specification& spec;     ///< Boundary specification
    const real_t                    lower_E;  ///< Specific total energy
    const real_t                    upper_E;  ///< Specific total energy
    const unsigned int              flags;    ///< Which BCs are enforced?

public:

    /** Prepare an instance for the given strides and parameters. */
    isothermal_enforcer(const grid_specification&       grid,
                        const std::ptrdiff_t            incf,
                        const isothermal_specification& spec,
                        const real_t                    lower_E,
                        const real_t                    upper_E)
        : Ny(grid.N.y())
        , incf(incf)
        , spec(spec)
        , lower_E(lower_E)
        , upper_E(upper_E)
        , flags(
            ENFORCE_LOWER_E  * !(boost::math::isnan)(lower_E     )
          | ENFORCE_LOWER_U  * !(boost::math::isnan)(spec.lower_u)
          | ENFORCE_LOWER_V  * !(boost::math::isnan)(spec.lower_v)
          | ENFORCE_LOWER_W  * !(boost::math::isnan)(spec.lower_w)
          | ENFORCE_LOWER_CS
          | ENFORCE_UPPER_E  * !(boost::math::isnan)(upper_E     )
          | ENFORCE_UPPER_U  * !(boost::math::isnan)(spec.upper_u)
          | ENFORCE_UPPER_V  * !(boost::math::isnan)(spec.upper_v)
          | ENFORCE_UPPER_W  * !(boost::math::isnan)(spec.upper_w)
          | ENFORCE_UPPER_CS
          )
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
        if (flags & ENFORCE_LOWER_E) LOWER(ndx::e) = lower_rho * lower_E;
        if (flags & ENFORCE_UPPER_E) UPPER(ndx::e) = upper_rho * upper_E;

        // Set streamwise momentum to be density times prescribed velocity
        if (flags & ENFORCE_LOWER_U) LOWER(ndx::mx) = lower_rho * spec.lower_u;
        if (flags & ENFORCE_UPPER_U) UPPER(ndx::mx) = upper_rho * spec.upper_u;

        // Set wall-normal momentum to be density times prescribed velocity
        if (flags & ENFORCE_LOWER_V) LOWER(ndx::my) = lower_rho * spec.lower_v;
        if (flags & ENFORCE_UPPER_V) UPPER(ndx::my) = upper_rho * spec.upper_v;

        // Set spanwise momentum to be density times prescribed velocity
        if (flags & ENFORCE_LOWER_W) LOWER(ndx::mz) = lower_rho * spec.lower_w;
        if (flags & ENFORCE_UPPER_W) UPPER(ndx::mz) = upper_rho * spec.upper_w;

        // Do nothing to the density equation

        // Set species partial densities to be density times mixture fraction.
        if (flags & (ENFORCE_LOWER_CS | ENFORCE_UPPER_CS)) {
            assert(spec.lower_cs.size() == spec.upper_cs.size());
            const std::size_t num_species = spec.lower_cs.size();
            for (std::size_t s = 1; s < num_species; ++s) {
                LOWER(ndx::rho_(s)) = lower_rho * spec.lower_cs[s];
                UPPER(ndx::rho_(s)) = upper_rho * spec.upper_cs[s];
            }
        } else if (flags & ENFORCE_LOWER_CS) {
            const std::size_t num_species = spec.lower_cs.size();
            for (std::size_t s = 1; s < num_species; ++s) {
                LOWER(ndx::rho_(s)) = lower_rho * spec.lower_cs[s];
            }
        } else if (flags & ENFORCE_UPPER_CS) {
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

    // Perform the usual mass_operator solve across all equations
    base::invert_mass_plus_scaled_operator(
            phi, state, method, delta_t, substep_index, ic0);

    // State leaves method as coefficients in X, Y, and Z directions
}

} // namespace suzerain
