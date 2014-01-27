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

#ifndef SUZERAIN_REACTING_EXPLICIT_OPERATOR_HPP
#define SUZERAIN_REACTING_EXPLICIT_OPERATOR_HPP

/** @file
 * Fully explicit, linearization-ready Navier--Stokes operators.
 */

#include "nonlinear_operator_fwd.hpp"

#include <suzerain/isothermal_mass_operator.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/specification_grid.hpp>
#include <suzerain/state_fwd.hpp>

#include "reacting.hpp"
#include "manufactured_solution.hpp"
#include "definition_channel.hpp"
#include "antioch_constitutive.hpp"
#include "definition_filter.hpp"

#pragma warning(disable:383 1572)

namespace suzerain {

// Forward declarations
class specification_largo;

namespace reacting {

/**
 * A boundary-condition agnostic, fully explicit Navier&ndash;Stokes operator.
 *
 * @see apply_navier_stokes_spatial_operator for the guts of the implementation.
 */
class explicit_nonlinear_operator
    : public operator_base,
      public lowstorage::nonlinear_operator< contiguous_state<4,complex_t> >
{
public:

    explicit_nonlinear_operator(
            const antioch_constitutive& cmods,
            const specification_grid &grid,
            const pencil_grid &dgrid,
            const bsplineop &cop,
            bspline &b,
            operator_common_block &common,
            const definition_filter &fsdef,
            const specification_largo &sgdef,
            const shared_ptr<const manufactured_solution>& msoln);

    virtual std::vector<real_t> apply_operator(
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const lowstorage::method_interface<complex_t> &method,
            const std::size_t substep_index) const;

protected:

    /** Contains constitutive models and related parameters */
    const antioch_constitutive& cmods;

    /** Houses data additionally required for some linear operators */
    operator_common_block &common;

    /** Holds optional manufactured solution forcing details */
    const shared_ptr<const manufactured_solution> msoln;

    /** The filter source definition */
    const definition_filter &fsdef;

    /** The slow growth definition */
    const specification_largo &sgdef;

private:

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;

    // boost::noncopyable trips Intel non-virtual base destructor warnings.
    explicit_nonlinear_operator(const explicit_nonlinear_operator&);
    explicit_nonlinear_operator& operator=(const explicit_nonlinear_operator&);

};

/**
 * A mass operator that provides no slip, isothermal walls.  It requires
 * interoperation with explicit_nonlinear_operator via operator_common_block.
 */
class isothermal_mass_operator : public suzerain::isothermal_mass_operator
{

public:

    isothermal_mass_operator(
            const antioch_constitutive& cmods,
            const specification_isothermal &isospec,
            const definition_channel &chdef,
            const specification_grid &grid,
            const pencil_grid &dgrid,
            const bsplineop &cop,
            bspline &b,
            operator_common_block &common);

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

    /** Contains constitutive models and related parameters */
    const antioch_constitutive& cmods;

    /** The channel flow case for which the operator is used */
    const definition_channel &chdef;

    /** Houses data required for \ref invert_mass_plus_scaled_operator */
    operator_common_block &common;

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;
};

} // namespace reacting

} // namespace suzerain

#endif  /* SUZERAIN_REACTING_EXPLICIT_OPERATOR_HPP */
