//--------------------------------------------------------------------------
//
// Copyright (C) 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
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
// hybrid_operator.hpp: Hybrid implicit/explicit Navier--Stokes operators
// $Id$

#ifndef HYBRID_OPERATOR_HPP
#define HYBRID_OPERATOR_HPP

#include <suzerain/grid_definition.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/spec_zgbsv.hpp>
#include <suzerain/state_fwd.hpp>

#include "../logging.hpp"
#include "nonlinear_operator_fwd.hpp"

#pragma warning(disable:383 1572)

namespace suzerain { namespace perfect {

/**
 * A hybrid implicit operator that provides no slip, isothermal walls.  It
 * requires interoperation with HybridNonlinearOperator via
 * operator_common_block.
 */
class HybridIsothermalLinearOperator
  : public operator_base,
    public timestepper::lowstorage::linear_operator<
        multi_array::ref<complex_t,4>,
        contiguous_state<4,complex_t>
    >
{
public:

    /** Determines the extent of the implicit treatment */
    static const linearize::type linearization = linearize::rhome;

    HybridIsothermalLinearOperator(
            const spec_zgbsv& spec,
            const scenario_definition &scenario,
            const grid_definition &grid,
            const pencil_grid &dgrid,
            bspline &b,
            const bsplineop &bop,
            operator_common_block &common)
        : operator_base(grid, dgrid, b, bop),
          spec(spec),
          scenario(scenario),
          common(common)
    {
        INFO0("HybridIsothermalLinearOperator solving using "
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
    const spec_zgbsv spec;

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
class HybridNonlinearOperator
    : public operator_base,
      public timestepper::nonlinear_operator<
            contiguous_state<4,complex_t>
      >
{
public:

    /** Determines the implicit treatment of the paired linear_operator. */
    static const linearize::type linearization
            = HybridIsothermalLinearOperator::linearization;

    HybridNonlinearOperator(
            const scenario_definition &scenario,
            const grid_definition &grid,
            const pencil_grid &dgrid,
            bspline &b,
            const bsplineop &bop,
            operator_common_block &common,
            const boost::shared_ptr<const manufactured_solution>& msoln)
        : operator_base(grid, dgrid, b, bop),
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
    const boost::shared_ptr<const manufactured_solution> msoln;

private:

    // Using boost::noncopyable trips Intel non-virtual base destructor warnings.
    HybridNonlinearOperator(const HybridNonlinearOperator&);
    HybridNonlinearOperator& operator=(const HybridNonlinearOperator&);

};

} /* namespace perfect */ } /* namespace suzerain */

#endif  /* HYBRID_OPERATOR_HPP */
