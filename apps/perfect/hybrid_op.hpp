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
// hybrid_op.hpp: Operators for hybrid implicit/explicit time advance
// $Id$

#ifndef HYBRID_OP_HPP
#define HYBRID_OP_HPP

#include <suzerain/grid_definition.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/spec_zgbsv.hpp>
#include <suzerain/state_fwd.hpp>

#include "../logging.hpp"
#include "nonlinear_fwd.hpp"

#pragma warning(disable:383 1572)

namespace suzerain { namespace perfect {

/**
 * A hybrid implicit operator that provides no slip, isothermal walls.  It
 * requires interoperation with HybridNonlinearOperator via
 * OperatorCommonBlock.
 */
class HybridIsothermalLinearOperator
  : public operator_base,
    public timestepper::lowstorage::ILinearOperator<
        multi_array::ref<complex_t,4>,
        contiguous_state<4,complex_t>
    >
{
public:

    /** Determines the extent of the implicit treatment */
    static const linearize::type linearization = linearize::rhome;

    HybridIsothermalLinearOperator(
            const spec_zgbsv& spec,
            const ScenarioDefinition &scenario,
            const problem::GridDefinition &grid,
            const pencil_grid &dgrid,
            bspline &b,
            const bsplineop &bop,
            OperatorCommonBlock &common)
        : operator_base(grid, dgrid, b, bop),
          spec(spec),
          scenario(scenario),
          common(common)
    {
        INFO0("HybridIsothermalLinearOperator solving using "
              << static_cast<std::string>(spec));
    }

    virtual void applyMassPlusScaledOperator(
             const complex_t &phi,
             multi_array::ref<complex_t,4> &state,
             const timestepper::lowstorage::IMethod<complex_t> &method,
             const component delta_t,
             const std::size_t substep_index) const;

    virtual void accumulateMassPlusScaledOperator(
            const complex_t &phi,
            const multi_array::ref<complex_t,4> &input,
            const complex_t &beta,
            contiguous_state<4,complex_t> &output,
            const timestepper::lowstorage::IMethod<complex_t> &method,
            const component delta_t,
            const std::size_t substep_index) const;

    /**
     * Performs the following steps:
     * <ul>
     * <li>
     *     channel_treatment step (3) performs the operator solve which for
     *     the implicit treatment must be combined with boundary conditions
     * </li><li>
     *     channel_treatment step (8) sets no-slip conditions
     *     on wall collocation points.
     * </li><li>
     * </li>
     *     channel_treatment step (9) sets isothermal conditions at walls
     *     using rho_wall = e_wall * gamma * (gamma - 1).
     * </ul>
     */
    virtual void invertMassPlusScaledOperator(
            const complex_t &phi,
            multi_array::ref<complex_t,4> &state,
            const timestepper::lowstorage::IMethod<complex_t> &method,
            const component delta_t,
            const std::size_t substep_index,
            multi_array::ref<complex_t,4> *ic0 = NULL) const;

protected:

    /** Controls the solves performed during invertMassPlusScaledOperator */
    const spec_zgbsv spec;

    /** The scenario in which the operator is used */
    const ScenarioDefinition &scenario;

    /** Houses data required for operator application and inversion */
    OperatorCommonBlock &common;

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
 * @see perfect::applyNonlinearOperator for the guts of the implementation.
 */
class HybridNonlinearOperator
    : public operator_base,
      public timestepper::nonlinear_operator<
            contiguous_state<4,complex_t>
      >
{
public:

    /** Determines the implicit treatment of the paired ILinearOperator. */
    static const linearize::type linearization
            = HybridIsothermalLinearOperator::linearization;

    HybridNonlinearOperator(
            const ScenarioDefinition &scenario,
            const problem::GridDefinition &grid,
            const pencil_grid &dgrid,
            bspline &b,
            const bsplineop &bop,
            OperatorCommonBlock &common,
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
    const ScenarioDefinition &scenario;

    /** Houses data additionally required for some linear operators */
    OperatorCommonBlock &common;

    /** Holds optional manufactured solution forcing details */
    const boost::shared_ptr<const manufactured_solution> msoln;

private:

    // Using boost::noncopyable trips Intel non-virtual base destructor warnings.
    HybridNonlinearOperator(const HybridNonlinearOperator&);
    HybridNonlinearOperator& operator=(const HybridNonlinearOperator&);

};

} /* namespace perfect */ } /* namespace suzerain */

#endif  /* HYBRID_OP_HPP */
