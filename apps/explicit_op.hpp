//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2011 The PECOS Development Team
//
// Please see http://pecos.ices.utexas.edu for more information.
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
//
// explicit_op.hpp: Operators for channel_explicit
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef EXPLICIT_OP_HPP
#define EXPLICIT_OP_HPP

#include <suzerain/grid_definition.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/scenario_definition.hpp>
#include <suzerain/state_fwd.hpp>

#include "precision.hpp"
#include "channel.hpp"

#pragma warning(disable:383 1572)

namespace channel {

/**
 * Storage for holding quantities computed during nonlinear operator
 * application which are required for linear operator application.
 */
class OperatorCommonBlock
{
public:
    // See http://eigen.tuxfamily.org/dox-devel/TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    OperatorCommonBlock() {}

    Eigen::ArrayXXr data;

    // Prepare logical indices using a struct for scoping (e.g. mean::u).
    struct mean { enum {
          u,
          count // Sentry
    }; };

    /** Zero storage within common block and reset its size as specified. */
    void reset(const std::size_t Ny) { data.setZero(mean::count, Ny); }

private:

    // Using boost::noncopyable trips Intel non-virtual base destructor warnings.
    OperatorCommonBlock(const OperatorCommonBlock&);
    OperatorCommonBlock& operator=(const OperatorCommonBlock&);
};

/** An operator which merely applies or inverts a B-spline mass matrix */
class BsplineMassOperator
  : public suzerain::OperatorBase<real_t>,
    public suzerain::timestepper::lowstorage::ILinearOperator<
        suzerain::ContiguousState<4,complex_t>
    >
{
public:

    typedef suzerain::ContiguousState<4,complex_t> state_type;

    BsplineMassOperator(
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop);

    virtual void applyMassPlusScaledOperator(
             const complex_t &phi,
             state_type &state,
             const std::size_t substep_index) const;

     virtual void accumulateMassPlusScaledOperator(
             const complex_t &phi,
             const state_type &input,
             const complex_t &beta,
             state_type &output,
             const std::size_t substep_index) const;

     virtual void invertMassPlusScaledOperator(
             const complex_t &phi,
             state_type &state,
             const std::size_t substep_index) const;

private:

    suzerain::bsplineop_luz massluz;
};

/**
 * An mass operator that forces bulk momentum and provides no slip, isothermal
 * walls.
 */
class BsplineMassOperatorIsothermal
  : public BsplineMassOperator
{
public:

    typedef BsplineMassOperator base;

    BsplineMassOperatorIsothermal(
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop,
            OperatorCommonBlock &common);

    virtual void invertMassPlusScaledOperator(
            const complex_t &phi,
            state_type &state,
            const std::size_t substep_index) const;

protected:

    Eigen::VectorXr bulkcoeff;

    OperatorCommonBlock &common;

};

/**
 * A boundary-condition agnostic, fully explicit Navier&ndash;Stokes operator.
 */
class NonlinearOperator
    : public suzerain::OperatorBase<real_t>,
      public suzerain::timestepper::INonlinearOperator<
            suzerain::ContiguousState<4,complex_t>
      >
{
public:

    typedef suzerain::ContiguousState<4,complex_t> state_type;

    NonlinearOperator(
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop,
            OperatorCommonBlock &common,
            const boost::shared_ptr<
                  const channel::manufactured_solution>& msoln);

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

/**
 * A fully explicit Navier&ndash;Stokes operator for isothermal boundaries.
 */
class NonlinearOperatorIsothermal
    : public NonlinearOperator
{
public:

    typedef NonlinearOperator base;

    NonlinearOperatorIsothermal(
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop,
            OperatorCommonBlock &common,
            const boost::shared_ptr<
                      const channel::manufactured_solution>& msoln);

    virtual std::vector<real_t> applyOperator(
            const real_t time,
            suzerain::ContiguousState<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag,
            const std::size_t substep_index) const;

protected:

    Eigen::VectorXr bulkcoeff;

};


} // namespace channel

#endif  /* EXPLICIT_OP_HPP */
