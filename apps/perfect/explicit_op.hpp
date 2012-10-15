//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012 Rhys Ulerich
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
// explicit_op.hpp: Fully explicit Navier--Stokes operators
// $Id$

#ifndef EXPLICIT_OP_HPP
#define EXPLICIT_OP_HPP

#include "nonlinear_fwd.hpp"

#include <suzerain/grid_definition.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/state_fwd.hpp>

#include "perfect.hpp"

#pragma warning(disable:383 1572)

namespace suzerain { namespace perfect {

/**
 * A boundary-condition agnostic, fully explicit Navier&ndash;Stokes operator.
 *
 * @see support::applyNonlinearOperator for the guts of the implementation.
 */
class NonlinearOperator
    : public suzerain::OperatorBase,
      public suzerain::timestepper::INonlinearOperator<
            suzerain::ContiguousState<4,complex_t>
      >
{
public:

    NonlinearOperator(
            const ScenarioDefinition &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop,
            OperatorCommonBlock &common,
            const boost::shared_ptr<
                  const perfect::manufactured_solution>& msoln)
        : suzerain::OperatorBase(grid, dgrid, b, bop),
          scenario(scenario),
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

    /** The scenario in which the operator is used */
    const ScenarioDefinition &scenario;

    /** Houses data additionally required for some linear operators */
    OperatorCommonBlock &common;

    /** Holds optional manufactured solution forcing details */
    const boost::shared_ptr<const perfect::manufactured_solution> msoln;

private:

    // Using boost::noncopyable trips Intel non-virtual base destructor warnings.
    NonlinearOperator(const NonlinearOperator&);
    NonlinearOperator& operator=(const NonlinearOperator&);

};

/** An operator which applies or inverts a B-spline mass matrix */
class BsplineMassOperator
  : public suzerain::OperatorBase,
    public suzerain::timestepper::lowstorage::ILinearOperator<
        suzerain::multi_array::ref<complex_t,4>,
        suzerain::ContiguousState<4,complex_t>
    >
{
public:

    BsplineMassOperator(
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop);

    virtual void applyMassPlusScaledOperator(
             const complex_t &phi,
             suzerain::multi_array::ref<complex_t,4> &state,
             const suzerain::timestepper::lowstorage::IMethod<complex_t> &method,
             const component delta_t,
             const std::size_t substep_index) const;

     virtual void accumulateMassPlusScaledOperator(
             const complex_t &phi,
             const suzerain::multi_array::ref<complex_t,4> &input,
             const complex_t &beta,
             suzerain::ContiguousState<4,complex_t> &output,
             const suzerain::timestepper::lowstorage::IMethod<complex_t> &method,
             const component delta_t,
             const std::size_t substep_index) const;

     virtual void invertMassPlusScaledOperator(
             const complex_t &phi,
             suzerain::multi_array::ref<complex_t,4> &state,
             const suzerain::timestepper::lowstorage::IMethod<complex_t> &method,
             const component delta_t,
             const std::size_t substep_index) const;

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
            const ScenarioDefinition &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop,
            OperatorCommonBlock &common)
        : BsplineMassOperator(grid, dgrid, b, bop),
          scenario(scenario),
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
            const suzerain::timestepper::lowstorage::IMethod<complex_t> &method,
            const component delta_t,
            const std::size_t substep_index) const;

protected:

    /** The scenario in which the operator is used */
    const ScenarioDefinition &scenario;

    /** Houses data required for \ref invertMassPlusScaledOperator */
    OperatorCommonBlock &common;

};

} /* namespace perfect */ } /* namespace suzerain */

#endif  /* EXPLICIT_OP_HPP */
