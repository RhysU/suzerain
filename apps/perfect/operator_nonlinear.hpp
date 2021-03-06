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

#ifndef SUZERAIN_PERFECT_OPERATOR_NONLINEAR_HPP
#define SUZERAIN_PERFECT_OPERATOR_NONLINEAR_HPP

/** @file
 * A nonlinear Navier--Stokes spatial operator.
 */

#include <suzerain/operator_mass.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/state_fwd.hpp>

namespace suzerain {

// Forward declarations
class specification_grid;
class specification_largo;
class pencil_grid;

namespace perfect {

// Forward declarations
class operator_common_block;
class manufactured_solution;
class definition_scenario;

/**
 * A boundary-condition agnostic Navier&ndash;Stokes operator.
 * The operator can be used in fully explicit or hybrid explicit/implicit mode.
 *
 * @see apply_navier_stokes_spatial_operator for the guts of the implementation.
 */
class operator_nonlinear
    : public operator_base,
      public lowstorage::operator_nonlinear< contiguous_state<4,complex_t> >
{
public:

    operator_nonlinear(
            const definition_scenario &scenario,
            const specification_grid &grid,
            const pencil_grid &dgrid,
            const bsplineop &cop,
            bspline &b,
            operator_common_block &common,
            specification_largo& sg,
            const shared_ptr<const manufactured_solution>& msoln);

    virtual std::vector<real_t> apply_operator(
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const lowstorage::method_interface<complex_t> &method,
            const std::size_t substep_index) const;

protected:

    /** The scenario in which the operator is used */
    const definition_scenario &scenario;

    /** Houses data additionally required for some linear operators */
    operator_common_block &common;

    /** Houses data around slow growth forcing computed via Largo */
    specification_largo& sg;

    /** Holds optional manufactured solution forcing details */
    const shared_ptr<const manufactured_solution> msoln;

private:

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;

    // boost::noncopyable trips Intel non-virtual base destructor warnings.
    operator_nonlinear(const operator_nonlinear&);
    operator_nonlinear& operator=(const operator_nonlinear&);

};

} // namespace perfect

} // namespace suzerain

#endif  /* SUZERAIN_PERFECT_OPERATOR_NONLINEAR_HPP */
