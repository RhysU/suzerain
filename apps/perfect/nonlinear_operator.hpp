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

#ifndef SUZERAIN_PERFECT_NONLINEAR_OPERATOR_HPP
#define SUZERAIN_PERFECT_NONLINEAR_OPERATOR_HPP

/** @file
 * A nonlinear Navier--Stokes spatial operator.
 */

#include <suzerain/mass_operator.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/state_fwd.hpp>

namespace suzerain {

// Forward declarations
class grid_specification;
class largo_specification;
class pencil_grid;

namespace perfect {

// Forward declarations
class operator_common_block;
class manufactured_solution;
class scenario_definition;

/**
 * A boundary-condition agnostic Navier&ndash;Stokes operator.
 * The operator can be used in fully explicit or hybrid explicit/implicit mode.
 *
 * @see apply_navier_stokes_spatial_operator for the guts of the implementation.
 */
class nonlinear_operator
    : public operator_base,
      public lowstorage::nonlinear_operator< contiguous_state<4,complex_t> >
{
public:

    nonlinear_operator(
            const scenario_definition &scenario,
            const grid_specification &grid,
            const pencil_grid &dgrid,
            const bsplineop &cop,
            bspline &b,
            operator_common_block &common,
            largo_specification& sg,
            const shared_ptr<const manufactured_solution>& msoln);

    virtual std::vector<real_t> apply_operator(
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const lowstorage::method_interface<complex_t> &method,
            const std::size_t substep_index) const;

protected:

    /** The scenario in which the operator is used */
    const scenario_definition &scenario;

    /** Houses data additionally required for some linear operators */
    operator_common_block &common;

    /** Houses data around slow growth forcing computed via Largo */
    largo_specification& sg;

    /** Holds optional manufactured solution forcing details */
    const shared_ptr<const manufactured_solution> msoln;

private:

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;

    // boost::noncopyable trips Intel non-virtual base destructor warnings.
    nonlinear_operator(const nonlinear_operator&);
    nonlinear_operator& operator=(const nonlinear_operator&);

};

} // namespace perfect

} // namespace suzerain

#endif  /* SUZERAIN_PERFECT_NONLINEAR_OPERATOR_HPP */
