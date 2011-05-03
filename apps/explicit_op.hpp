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
#include <suzerain/pencil_grid.hpp>
#include <suzerain/scenario_definition.hpp>
#include <suzerain/state_fwd.hpp>

#include "precision.hpp"
#include "channel.hpp"

#pragma warning(disable:383 1572)

namespace channel {

/** Helper infrastructure for operator implementations. */
class OperatorBase
{
public:

    OperatorBase(
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop);

protected:

    const suzerain::problem::ScenarioDefinition<real_t> &scenario;
    const suzerain::problem::GridDefinition &grid;
    const suzerain::pencil_grid &dgrid;
    const suzerain::bsplineop &bop;

    /** Shorthand for scaled operator accumulation */
    template<typename AlphaType, typename MultiArrayX,
             typename BetaType,  typename MultiArrayY>
    int bop_accumulate(
            int nderiv,
            const AlphaType& alpha, const MultiArrayX &x, int ndx_x,
            const BetaType& beta,         MultiArrayY &y, int ndx_y) const
    {
        assert(x.shape()[1] == (unsigned) bop.n());
        assert((unsigned) x.strides()[3] == x.shape()[2] * x.strides()[2] );
        assert((unsigned) y.strides()[3] == y.shape()[2] * y.strides()[2] );
        assert(std::equal(x.shape() + 1, x.shape() + 4, y.shape() + 1));

        return bop.accumulate(
                nderiv, x.shape()[2] * x.shape()[3],
                alpha,  x[ndx_x].origin(), x.strides()[1], x.strides()[2],
                beta,   y[ndx_y].origin(), y.strides()[1], y.strides()[2]);
    }

    /** Shorthand for scaled operator application */
    template<typename AlphaType, typename MultiArray>
    int bop_apply(
            int nderiv, const AlphaType& alpha, MultiArray &x, int ndx) const
    {
        assert(x.shape()[1] == (unsigned) bop.n());
        assert((unsigned) x.strides()[3] == x.shape()[2] * x.strides()[2]);

        return bop.apply(
                nderiv, x.shape()[2] * x.shape()[3],
                alpha,  x[ndx].origin(), x.strides()[1], x.strides()[2]);
    }

    /** Shorthand for wave space-based differentiation accumulation */
    template<typename MultiArrayX, typename MultiArrayY>
    void diffwave_accumulate(int dxcnt,
                             int dzcnt,
                             const typename MultiArrayX::element& alpha,
                             const MultiArrayX &x,
                             int ndx_x,
                             const typename MultiArrayY::element& beta,
                             MultiArrayY &y,
                             int ndx_y) const
    {
        assert(std::equal(x.shape()   + 1, x.shape()   + 4, y.shape()   + 1));
        assert(std::equal(x.strides() + 1, x.strides() + 4, y.strides() + 1));

        return suzerain::diffwave::accumulate(
                dxcnt, dzcnt,
                alpha, x[ndx_x].origin(),
                beta,  y[ndx_y].origin(),
                scenario.Lx, scenario.Lz,
                dgrid.global_wave_extent.y(),
                grid.N.x(),
                grid.dN.x(),
                dgrid.local_wave_start.x(),
                dgrid.local_wave_end.x(),
                grid.N.z(),
                grid.dN.z(),
                dgrid.local_wave_start.z(),
                dgrid.local_wave_end.z());
    }

    /** Shorthand for wave space-based differentiation application */
    template<typename MultiArray>
    void diffwave_apply(int dxcnt,
                        int dzcnt,
                        const typename MultiArray::element& alpha,
                        const MultiArray &x,
                        int ndx_x) const
    {
        return suzerain::diffwave::accumulate(
                dxcnt, dzcnt,
                alpha, x[ndx_x].origin(),
                scenario.Lx, scenario.Lz,
                dgrid.global_wave_extent.y(),
                grid.N.x(),
                grid.dN.x(),
                dgrid.local_wave_start.x(),
                dgrid.local_wave_end.x(),
                grid.N.z(),
                grid.dN.z(),
                dgrid.local_wave_start.z(),
                dgrid.local_wave_end.z());
    }

    // Only valid for j \in dgrid.local_physical_{start,end}.y()
    real_t x(std::size_t i) const {
        return i * scenario.Lx / grid.dN.x() - scenario.Lx / 2;
    }

    real_t y(std::size_t j) const {
        return y_[j];
    }

    real_t z(std::size_t k) const {
        return k * scenario.Lz / grid.dN.z() - scenario.Lz / 2;
    }

    // Only valid for j \in dgrid.local_physical_{start,end}.y()
    real_t one_over_delta_y(std::size_t j) const {
        return one_over_delta_y_[j];
    }

    const bool has_zero_zero_mode;
    const real_t one_over_delta_x;
    const real_t one_over_delta_z;

private:
    boost::multi_array<real_t,1> y_;
    boost::multi_array<real_t,1> one_over_delta_y_;
};

/** An operator which merely applies or inverts a B-spline mass matrix */
class BsplineMassOperator
  : public OperatorBase,
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
    Eigen::VectorXr massinv_elower;
    Eigen::VectorXr massinv_eupper;

private:

    void save_mean_state_at_collocation_points(const state_type &state) const;
    mutable Eigen::VectorXr mean_rho;
    mutable Eigen::VectorXr mean_rhou;

};


/**
 * A boundary-condition agnostic, fully explicit Navier&ndash;Stokes operator.
 */
class NonlinearOperator
    : public OperatorBase,
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
            const suzerain::bsplineop &bop);

    virtual real_t applyOperator(
            suzerain::NoninterleavedState<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag,
            const bool delta_t_requested = false) const;

protected:

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
            const suzerain::bsplineop &bop);

    virtual real_t applyOperator(
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
