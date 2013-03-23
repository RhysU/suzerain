//--------------------------------------------------------------------------
//
// Copyright (C) 2012, 2013 Rhys Ulerich
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

#ifndef SUZERAIN_PERFECT_HYBRID_OPERATOR_HPP
#define SUZERAIN_PERFECT_HYBRID_OPERATOR_HPP

/** @file
 * Hybrid implicit/explicit Navier--Stokes operators.
 */

#include <suzerain/bsmbsm_solver.hpp>
#include <suzerain/grid_specification.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/state_fwd.hpp>
#include <suzerain/timestepper.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/zgbsv_specification.hpp>

#include "common_block.hpp"
#include "manufactured_solution.hpp"
#include "navier_stokes.hpp"

#pragma warning(disable:383 1572)

namespace suzerain {

namespace perfect {

/**
 * A hybrid implicit operator that provides no slip, isothermal walls.  It
 * requires interoperation with nonlinear_operator via operator_common_block.
 */
class isothermal_hybrid_linear_operator
  : public operator_base,
    public timestepper::lowstorage::linear_operator<
        multi_array::ref<complex_t,4>,
        contiguous_state<4,complex_t>
    >
{
public:

    isothermal_hybrid_linear_operator(
            const zgbsv_specification& spec,
            const scenario_definition &scenario,
            const grid_specification &grid,
            const pencil_grid &dgrid,
            const bsplineop &cop,
            bspline &b,
            operator_common_block &common);

    ~isothermal_hybrid_linear_operator();

    virtual void apply_mass_plus_scaled_operator(
             const complex_t &phi,
             multi_array::ref<complex_t,4> &state,
             const timestepper::lowstorage::method_interface<complex_t> &method,
             const component delta_t,
             const std::size_t substep_index) const;

    virtual void accumulate_mass_plus_scaled_operator(
            const complex_t &phi,
            const multi_array::ref<complex_t,4> &input,
            const complex_t &beta,
            contiguous_state<4,complex_t> &output,
            const timestepper::lowstorage::method_interface<complex_t> &method,
            const component delta_t,
            const std::size_t substep_index) const;

    /**
     * Sets no-slip, isothermal boundary conditions at first and last
     * wall-normal collocation point using that
     * rho_wall = e_wall * gamma * (gamma - 1).
     */
    virtual void invert_mass_plus_scaled_operator(
            const complex_t &phi,
            multi_array::ref<complex_t,4> &state,
            const timestepper::lowstorage::method_interface<complex_t> &method,
            const component delta_t,
            const std::size_t substep_index,
            multi_array::ref<complex_t,4> *ic0 = NULL) const;

protected:

    /** Controls the solves performed during invert_mass_plus_scaled_operator */
    shared_ptr<bsmbsm_solver> solver;

    /** The scenario in which the operator is used */
    const scenario_definition &scenario;

    /** Houses data required for operator application and inversion */
    operator_common_block &common;

private:

    /** Pack scenario parameters for rholut_imexop.h usage */
    suzerain_rholut_imexop_scenario imexop_s() const;

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;

};

} // namespace perfect

} // namespace suzerain

#endif  /* SUZERAIN_PERFECT_HYBRID_OPERATOR_HPP */
