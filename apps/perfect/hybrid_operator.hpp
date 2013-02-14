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

#include "nonlinear_operator_fwd.hpp"

#include <suzerain/grid_specification.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/state_fwd.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/zgbsv_specification.hpp>

#include "manufactured_solution.hpp"

#pragma warning(disable:383 1572)

namespace suzerain { namespace perfect {

/**
 * A hybrid implicit operator that provides no slip, isothermal walls.  It
 * requires interoperation with hybrid_nonlinear_operator via
 * operator_common_block.
 */
class isothermal_hybrid_linear_operator
  : public operator_base,
    public timestepper::lowstorage::linear_operator<
        multi_array::ref<complex_t,4>,
        contiguous_state<4,complex_t>
    >
{
public:

    /** Determines the extent of the implicit treatment */
    static const linearize::type linearization = linearize::rhome;

    isothermal_hybrid_linear_operator(
            const zgbsv_specification& spec,
            const scenario_definition &scenario,
            const grid_specification &grid,
            const pencil_grid &dgrid,
            const bsplineop &cop,
            bspline &b,
            operator_common_block &common)
        : operator_base(grid, dgrid, cop, b),
          spec(spec),
          scenario(scenario),
          common(common)
    {
        INFO0("Linear isothermal_hybrid_linear_operator using "
              << static_cast<std::string>(spec));
    }

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
    const zgbsv_specification spec;

    /** The scenario in which the operator is used */
    const scenario_definition &scenario;

    /** Houses data required for operator application and inversion */
    operator_common_block &common;

private:

    /** Pack scenario parameters for rholut_imexop.h usage */
    suzerain_rholut_imexop_scenario imexop_s() const
    {
        suzerain_rholut_imexop_scenario retval;
        retval.Re    = scenario.Re;
        retval.Pr    = scenario.Pr;
        retval.Ma    = scenario.Ma;
        retval.alpha = scenario.alpha;
        retval.gamma = scenario.gamma;
        return retval;
    }

};

/**
 * A boundary-condition agnostic, hybrid explicit Navier&ndash;Stokes operator.
 *
 * @see apply_navier_stokes_spatial_operator for the guts of the implementation.
 */
class hybrid_nonlinear_operator
    : public operator_base,
      public timestepper::nonlinear_operator<
            contiguous_state<4,complex_t>
      >
{
public:

    /** Determines the implicit treatment of the paired linear_operator. */
    static const linearize::type linearization
            = isothermal_hybrid_linear_operator::linearization;

    hybrid_nonlinear_operator(
            const scenario_definition &scenario,
            const grid_specification &grid,
            const pencil_grid &dgrid,
            const bsplineop &cop,
            bspline &b,
            operator_common_block &common,
            const shared_ptr<const manufactured_solution>& msoln)
        : operator_base(grid, dgrid, cop, b),
          scenario(scenario),
          common(common),
          msoln(msoln)
    {}

    virtual std::vector<real_t> apply_operator(
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag,
            const std::size_t substep_index) const;

protected:

    /** The scenario in which the operator is used */
    const scenario_definition &scenario;

    /** Houses data additionally required for some linear operators */
    operator_common_block &common;

    /** Holds optional manufactured solution forcing details */
    const shared_ptr<const manufactured_solution> msoln;

private:

    // Using boost::noncopyable trips Intel non-virtual base destructor warnings.
    hybrid_nonlinear_operator(const hybrid_nonlinear_operator&);
    hybrid_nonlinear_operator& operator=(const hybrid_nonlinear_operator&);

};

} /* namespace perfect */ } /* namespace suzerain */

#endif  /* SUZERAIN_PERFECT_HYBRID_OPERATOR_HPP */
