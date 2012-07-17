//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// explicit_op.hpp: Fully explicit Navier--Stokes operators
// $Id$

#ifndef EXPLICIT_OP_HPP
#define EXPLICIT_OP_HPP

#include <suzerain/grid_definition.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/scenario_definition.hpp>
#include <suzerain/state_fwd.hpp>

#include "../precision.hpp"
#include "../support.hpp"

#include "nonlinear_fwd.hpp"

#pragma warning(disable:383 1572)

namespace channel {

/**
 * A boundary-condition agnostic, fully explicit Navier&ndash;Stokes operator.
 *
 * @see channel::applyNonlinearOperator for the guts of the implementation.
 */
class NonlinearOperator
    : public suzerain::OperatorBase<real_t>,
      public suzerain::timestepper::INonlinearOperator<
            suzerain::ContiguousState<4,complex_t>
      >
{
public:

    NonlinearOperator(
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
    NonlinearOperator(const NonlinearOperator&);
    NonlinearOperator& operator=(const NonlinearOperator&);

};

/** An operator which applies or inverts a B-spline mass matrix */
class BsplineMassOperator
  : public suzerain::OperatorBase<real_t>,
    public suzerain::timestepper::lowstorage::ILinearOperator<
        suzerain::multi_array::ref<complex_t,4>,
        suzerain::ContiguousState<4,complex_t>
    >
{
public:

    BsplineMassOperator(
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

     /** Precomputed mass matrix factorization */
    suzerain::bsplineop_luz massluz;

};

/**
 * A mass operator that provides no slip, isothermal walls.  It requires
 * interoperation with NonlinearOperator via OperatorCommonBlock.
 */
class BsplineMassOperatorIsothermal
    : public BsplineMassOperator
{

    typedef BsplineMassOperator base;

public:

    BsplineMassOperatorIsothermal(
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop,
            OperatorCommonBlock &common)
        : BsplineMassOperator(scenario, grid, dgrid, b, bop),
          common(common)
    {}

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
            suzerain::multi_array::ref<complex_t,4> &state,
            const component delta_t,
            const std::size_t substep_index,
            const real_t iota) const;

protected:

    /** Houses data required for \ref invertMassPlusScaledOperator */
    OperatorCommonBlock &common;

};

} // namespace channel

#endif  /* EXPLICIT_OP_HPP */
