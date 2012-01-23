//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012 The PECOS Development Team
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
#include <suzerain/multi_array.hpp>
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
 * application which either are required for linear operator application or for
 * statistics sampling purposes.
 *
 * The quantities are as follows:
 * \li \c u  The nonlinear operator computes the instantaneous spatial (x, z)
 *     mean streamwise velocity profile as collocation point values.
 *     The linear operator then uses the information to compute
 *     the implicit \f$f\cdot{}u\f$ term in the total energy equation.
 * \li \c f The linear operator accumulates the time-step-specific
 *     temporal mean streamwise (x) component of the implicit \f$f\f$
 *     term in the momentum equation stored as coefficients.
 * \li \c f_dot_u The linear operator accumulates the time-step-specific
 *     temporal mean the implicit \f$f\cdot{}u\f$ term in the energy equation
 *     stored as coefficients.
 * \li \c qb The linear operator accumulates the time-step-specific
 *     temporal mean the implicit \f$q_b\f$ term in the energy equation
 *     stored as coefficients.
 *
 * "Time-step-specific temporal means" are time averages taken across a single
 * time step of quantities which vary on each substep.  As the substeps are all
 * of equal length, a simple running mean reset on substep zero and then
 * accumulated.
 */
class OperatorCommonBlock
{
public:
    // See http://eigen.tuxfamily.org/dox-devel/TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    OperatorCommonBlock() {}

private:

    /** Type of the contiguous storage used to house all scalars */
    typedef Eigen::Array<real_t, Eigen::Dynamic, 4> storage_type;

public:

    /** Contiguous storage used to house all scalars */
    storage_type storage;

    // Declare a named, mutable "view" into \c storage for each quantity
    storage_type::ColXpr u()       { return storage.col(0); }
    storage_type::ColXpr f()       { return storage.col(1); }
    storage_type::ColXpr f_dot_u() { return storage.col(2); }
    storage_type::ColXpr qb()      { return storage.col(3); }

    // Declare a named, mutable "view" into \c storage for each quantity
    storage_type::ConstColXpr u()       const { return storage.col(0); }
    storage_type::ConstColXpr f()       const { return storage.col(1); }
    storage_type::ConstColXpr f_dot_u() const { return storage.col(2); }
    storage_type::ConstColXpr qb()      const { return storage.col(3); }

private:

    // Using boost::noncopyable trips Intel non-virtual base destructor warnings.
    OperatorCommonBlock(const OperatorCommonBlock&);
    OperatorCommonBlock& operator=(const OperatorCommonBlock&);
};

/**
 * A boundary-condition agnostic, fully explicit Navier&ndash;Stokes operator.
 *
 * During \ref applyOperator the instantaneous wall-normal velocity is averaged
 * across the streamwise and spanwise directions and stored into
 * OperatorCommonBlock::u() using an instance provided at construction time.
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
            const std::size_t substep_index) const
    {
        // Dispatch to implementation paying nothing for substep-related ifs
        if (substep_index == 0) {
            return applyOperator<true >(time, swave, evmaxmag_real, evmaxmag_imag);
        } else {
            return applyOperator<false>(time, swave, evmaxmag_real, evmaxmag_imag);
        }
    }

protected:

    /** Houses data additionally required for some linear operators */
    OperatorCommonBlock &common;

    /** Holds optional manufactured solution forcing details */
    const boost::shared_ptr<const channel::manufactured_solution> msoln;

private:

    // Internal implementation of INonlinearOperator::applyOperator logic
    template< bool zeroth_substep >
    std::vector<real_t> applyOperator(
            const real_t time,
            suzerain::ContiguousState<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag) const;

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
 * A mass operator that forces bulk momentum and provides no slip, isothermal
 * walls.  It requires interoperation with NonlinearOperator via
 * OperatorCommonBlock.
 *
 * During \ref invertMassPlusScaledOperator implicit momentum forcing is
 * applied following the section of <tt>writeups/channel_treatment.tex</tt>
 * titled "Enforcing a target bulk momentum via the linear operator" and using
 * information from OperatorCommonBlock::u().
 *
 * Also during \ref invertMassPlusScaledOperator, OperatorCommonBlock::f(),
 * OperatorCommonBlock::f_dot_u(), and OperatorCommonBlock::qb() are
 * accumulated using ILowStorageMethod::iota().
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
            suzerain::multi_array::ref<complex_t,4> &state,
            const component delta_t,
            const std::size_t substep_index,
            const real_t iota) const;

protected:

    /** Precomputed integration coefficients */
    Eigen::VectorXr bulkcoeff;

    /** Houses data required for \ref invertMassPlusScaledOperator */
    OperatorCommonBlock &common;

};

} // namespace channel

#endif  /* EXPLICIT_OP_HPP */
