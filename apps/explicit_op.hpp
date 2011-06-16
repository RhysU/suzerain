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

/** An operator which merely applies or inverts a B-spline mass matrix */
class BsplineMassOperator
  : public suzerain::OperatorBase<real_t>,
    public suzerain::timestepper::lowstorage::ILinearOperator<
        suzerain::NoninterleavedState<4,complex_t>
    >
{
public:

    typedef suzerain::NoninterleavedState<4,complex_t> state_type;

    BsplineMassOperator(
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop);

    virtual void applyMassPlusScaledOperator(
             const complex_t &phi,
             state_type &state) const;

     virtual void accumulateMassPlusScaledOperator(
             const complex_t &phi,
             const state_type &input,
             const complex_t &beta,
             state_type &output) const;

     virtual void invertMassPlusScaledOperator(
             const complex_t &phi,
             state_type &state) const;

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
            const suzerain::bsplineop &bop);

    virtual void applyMassPlusScaledOperator(
            const complex_t &phi,
            state_type &state) const;

    virtual void accumulateMassPlusScaledOperator(
            const complex_t &phi,
            const state_type &input,
            const complex_t &beta,
            state_type &output) const;

    virtual void invertMassPlusScaledOperator(
            const complex_t &phi,
            state_type &state) const;

protected:

    Eigen::VectorXr bulkcoeff;

private:

    void save_mean_state_at_collocation_points(const state_type &state) const;
    mutable Eigen::VectorXr saved_mean_rho;
    mutable Eigen::VectorXr saved_mean_rhou;

};


/**
 * A boundary-condition agnostic, fully explicit Navier&ndash;Stokes operator.
 */
class NonlinearOperator
    : public boost::noncopyable,                        // Nontrivial storage
      public suzerain::OperatorBase<real_t>,
      public suzerain::timestepper::INonlinearOperator<
            suzerain::NoninterleavedState<4,complex_t>
      >
{

public:
    typedef suzerain::NoninterleavedState<4,complex_t> state_type;

    NonlinearOperator(
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop,
            const boost::shared_ptr<
                const nsctpl_rholut::manufactured_solution<real_t> >& msoln);

    virtual real_t applyOperator(
            const real_t time,
            suzerain::NoninterleavedState<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag,
            const bool delta_t_requested = false) const;

protected:

    /** Holds optional manufactured solution forcing details */
    const boost::shared_ptr<
                const nsctpl_rholut::manufactured_solution<real_t>
          > msoln;

    /** Auxiliary scalar-field storage used within applyOperator */
    mutable state_type auxw;

    // Inner struct purely for name-scoping purposes
    struct aux {

        // Logical indexes into auxiliary scalar-field storage
        // Seemingly goofy ordering matches usage in applyOperator
        enum {
            rho_y, rho_yy,
            rho_x, rho_xx, rho_xz, rho_z, rho_zz, rho_xy, rho_yz,
            mx_y, mx_yy,
            mx_x, mx_xx, mx_xz, mx_z, mx_zz, mx_xy, mx_yz,
            my_y, my_yy,
            my_x, my_xx, my_xz, my_z, my_zz, my_xy, my_yz,
            mz_y, mz_yy,
            mz_x, mz_xx, mz_xz, mz_z, mz_zz, mz_xy, mz_yz,
            e_y, div_grad_e, e_x, e_z,

            count // Sentry
        };

        private: aux();
    };

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
            const boost::shared_ptr<
                const nsctpl_rholut::manufactured_solution<real_t> >& msoln);

    virtual real_t applyOperator(
            const real_t time,
            suzerain::NoninterleavedState<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag,
            const bool delta_t_requested = false) const;

protected:

    Eigen::VectorXr bulkcoeff;
    mutable Eigen::VectorXr rho_fm;
    mutable Eigen::VectorXr fm_dot_m;

};


} // namespace channel

#endif  /* EXPLICIT_OP_HPP */
