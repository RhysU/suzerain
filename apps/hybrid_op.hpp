//--------------------------------------------------------------------------
//
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
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
#include <suzerain/scenario_definition.hpp>
#include <suzerain/state_fwd.hpp>

#include "precision.hpp"
#include "channel.hpp"
#include "nonlinear_fwd.hpp"

#pragma warning(disable:383 1572)

namespace channel {

/**
 * A hybrid implicit operator that forces bulk momentum and provides no slip,
 * isothermal walls.  It requires interoperation with HybridNonlinearOperator
 * via OperatorCommonBlock.
 *
 * Implicit treatment follows "Numerical considerations" within
 * <tt>writeups/derivation.tex</tt>.  During \ref invertMassPlusScaledOperator
 * implicit momentum forcing is applied following the section of
 * <tt>writeups/channel_treatment.tex</tt> titled "Enforcing a target bulk
 * momentum via the linear operator" and using information from
 * OperatorCommonBlock::u().
 *
 * Also during \ref invertMassPlusScaledOperator, OperatorCommonBlock::f(),
 * OperatorCommonBlock::f_dot_u(), and OperatorCommonBlock::qb() are
 * accumulated using ILowStorageMethod::iota().
 */
class HybridIsothermalLinearOperator
  : public suzerain::OperatorBase<real_t>,
    public suzerain::timestepper::lowstorage::ILinearOperator<
        suzerain::multi_array::ref<complex_t,4>,
        suzerain::ContiguousState<4,complex_t>
    >
{
public:

    /** Determines the extent of the implicit treatment */
    static const linearize::type linearization = linearize::rhome;

    HybridIsothermalLinearOperator(
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop);

    virtual void applyMassPlusScaledOperator(
             const complex_t &phi,
             suzerain::multi_array::ref<complex_t,4> &state,
             const component delta_t,
             const std::size_t substep_index) const;

     virtual void accumulateMassPlusScaledOperator(
             const complex_t &phi,
             const suzerain::multi_array::ref<complex_t,4> &input,
             const complex_t &beta,
             suzerain::ContiguousState<4,complex_t> &output,
             const component delta_t,
             const std::size_t substep_index) const;

     virtual void invertMassPlusScaledOperator(
             const complex_t &phi,
             suzerain::multi_array::ref<complex_t,4> &state,
             const component delta_t,
             const std::size_t substep_index,
             const real_t iota) const;

private:

    /** Precomputed integration coefficients */
    Eigen::VectorXr bulkcoeff;

    /** Houses data required for \ref invertMassPlusScaledOperator */
    OperatorCommonBlock &common;

};

/**
 * A boundary-condition agnostic, hybrid explicit Navier&ndash;Stokes operator.
 *
 * @see channel::applyNonlinearOperator for the guts of the implementation.
 */
class HybridNonlinearOperator
    : public suzerain::OperatorBase<real_t>,
      public suzerain::timestepper::INonlinearOperator<
            suzerain::ContiguousState<4,complex_t>
      >
{
public:

    /** Determines the implicit treatment of the paired ILinearOperator. */
    static const linearize::type linearization
            = HybridIsothermalLinearOperator::linearization;

    HybridNonlinearOperator(
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop,
            OperatorCommonBlock &common,
            const boost::shared_ptr<
                  const channel::manufactured_solution>& msoln)
        : suzerain::OperatorBase<real_t>(scenario, grid, dgrid, b, bop),
          common(common),
          msoln(msoln)
    {}

    virtual std::vector<real_t> applyOperator(
            const real_t time,
            suzerain::ContiguousState<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag,
            const std::size_t substep_index) const;

protected:

    /** Houses data additionally required for some linear operators */
    OperatorCommonBlock &common;

    /** Holds optional manufactured solution forcing details */
    const boost::shared_ptr<const channel::manufactured_solution> msoln;

private:

    // Using boost::noncopyable trips Intel non-virtual base destructor warnings.
    HybridNonlinearOperator(const HybridNonlinearOperator&);
    HybridNonlinearOperator& operator=(const HybridNonlinearOperator&);

};

} // namespace channel

#endif  /* HYBRID_OP_HPP */
