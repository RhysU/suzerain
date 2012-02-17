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
// channel_treatment.hpp: mixin providing channel problem forcing
// $Id$

#ifndef CHANNEL_TREATMENT_HPP
#define CHANNEL_TREATMENT_HPP

#include "nonlinear.hpp"

#include <suzerain/bspline.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/scenario_definition.hpp>

#include "precision.hpp"

namespace channel {

/**
 * A mixin providing channel problem treatment atop any
 * <tt>
 *    public suzerain::timestepper::INonlinearOperator<
 *          suzerain::ContiguousState<4,complex_t>
 *    >
 * </tt>.  During \ref invertMassPlusScaledOperator implicit momentum forcing
 * is applied following the section of <tt>writeups/channel_treatment.tex</tt>
 * titled "Enforcing a target bulk momentum via the linear operator" and using
 * information from OperatorCommonBlock::u() via instance provided at
 * construction time.
 *
 * Means of the implicit momentum and energy forcing coefficients are also
 * maintained across each individual time step for sampling the statistics
 * /bar_f, /bar_f_dot_u, and /bar_qb using OperatorCommonBlock.
 */
template< typename BaseClass >
class ChannelTreatment : public BaseClass
{
public:

    /**
     * Constructor delegating to BaseClass.
     *
     * BaseClass must make its constructor arguments available as member
     * variables under the same name as those found in this constructor.
     */
    ChannelTreatment(
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop,
            OperatorCommonBlock &common);

    /**
     * Force the channel problem delegating to BaseClass when appropriate.
     * The BaseClass must perform the following steps.
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
     * The remainder of the channel treatment is done by this mixin.
     */
    virtual void invertMassPlusScaledOperator(
            const complex_t &phi,
            suzerain::multi_array::ref<complex_t,4> &state,
            const real_t delta_t,
            const std::size_t substep_index,
            const real_t iota) const;

private:

    /** Precomputed integration coefficients */
    Eigen::VectorXr bulkcoeff;

};

template< typename BaseClass >
ChannelTreatment<BaseClass>::ChannelTreatment(
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop,
            OperatorCommonBlock &common)
    : BaseClass(scenario, grid, dgrid, b, bop, common),
      bulkcoeff(b.n())
{
    // Precompute operator for finding bulk quantities from coefficients
    b.integration_coefficients(0, bulkcoeff.data());
    bulkcoeff /= scenario.Ly;
}

template< typename BaseClass >
void ChannelTreatment<BaseClass>::invertMassPlusScaledOperator(
        const complex_t &phi,
        suzerain::multi_array::ref<complex_t,4> &state,
        const real_t delta_t,
        const std::size_t substep_index,
        const real_t iota) const
{
    // State enters method as coefficients in X and Z directions
    // State enters method as collocation point values in Y direction

    // Shorthand
    using Eigen::ArrayXc;
    using Eigen::ArrayXr;
    using Eigen::InnerStride;
    using Eigen::Map;
    using suzerain::inorder::wavenumber;
    namespace field = channel::field;
    namespace ndx   = field::ndx;

    // Bring some BaseClass-dependent superclass member names into scope
    OperatorCommonBlock &common
            = this->common;
    const suzerain::problem::ScenarioDefinition<real_t> &scenario
            = this->scenario;
    const suzerain::pencil_grid &dgrid
            = this->dgrid;

    // More shorthand
    const int Ny = dgrid.global_wave_extent.y();
    const std::size_t wall_lower    = 0;
    const std::size_t wall_upper    = Ny - 1;
    const bool has_zero_zero_modes  = dgrid.has_zero_zero_modes();

    // Sidesteps assertions when local rank contains no wavespace information
    if (SUZERAIN_UNLIKELY(0U == state.shape()[1])) return;

    // Incoming state has wall-normal pencils of interleaved state scalars?
    SUZERAIN_ENSURE(state.shape()  [1] == (unsigned) Ny);
    SUZERAIN_ENSURE(state.strides()[1] ==             1);
    SUZERAIN_ENSURE(state.strides()[0] == (unsigned) Ny);
    SUZERAIN_ENSURE(state.shape()  [0] ==  field::count);

    // Means of the implicit momentum and energy forcing coefficients(!) are
    // maintained across each individual time step for sampling the statistics
    // /bar_f, /bar_f_dot_u, and /bar_qb using OperatorCommonBlock via
    // ILowStorageMethod::iota a la mean += iota*(sample - mean).  Note that
    // the sample must be divided by delta_t to account for step sizes.

    // See channel_treatment writeup for information on the steps below.
    // Steps may appear out of order relative to the writeup.
    // For phi real, the implicit operator is purely real for km == kn == 0.
    // Success of the zero-zero mode treatment relies on real-valued phi!
    SUZERAIN_ENSURE(phi.imag() == 0);

    // channel_treatment step (1) done during nonlinear operator application
    // via shared OperatorCommonBlock storage space.

    // Integral constraints enabled only when parameters are non-inf, non-NaN.
    // Allow disabling these to meet manufactured solution verification needs.
    const bool constrain_bulk_rhou
            = (boost::math::isnormal)(scenario.bulk_rhou);
    const bool constrain_bulk_rho
            = (boost::math::isnormal)(scenario.bulk_rho);

    // channel_treatment step (2) loads ones and local mean streamwise velocity
    // into the imaginary part of the zero-zero mode non-wall collocation
    // points of the streamwise momentum and total energy, respectively.  These
    // wall collocation points are set to zero to avoid forcing at the walls
    // when enforcing integral constraints.
    //
    // The remaining equations have their zero-zero mode imaginary portions
    // zeroed to reduce cross-talk between equations during implicit solves.
    // These must be re-zeroed after the solve.
    if (constrain_bulk_rhou && has_zero_zero_modes) {

        const std::size_t start = wall_lower + 1;
        const std::size_t size  = wall_upper - start;

        Map<ArrayXc>(state[ndx::rho ].origin(), Ny).imag().setZero();

        Map<ArrayXc> mean_rhou(state[ndx::rhou].origin(), Ny);
        mean_rhou.imag()[wall_lower] = 0;
        mean_rhou.imag().segment(start, size).setOnes();
        mean_rhou.imag()[wall_upper] = 0;

        Map<ArrayXc>(state[ndx::rhov].origin(), Ny).imag().setZero();

        Map<ArrayXc>(state[ndx::rhow].origin(), Ny).imag().setZero();

        Map<ArrayXc> mean_rhoe(state[ndx::rhoe].origin(), Ny);
        mean_rhoe.imag()[wall_lower]          = 0;
        mean_rhoe.imag().segment(start, size) = common.u().segment(start, size);
        mean_rhoe.imag()[wall_upper]          = 0;
    }

    // The BaseClass is responsible for all of the following steps:
    //
    // channel_treatment step (3) performs the operator solve which for the
    // implicit treatment must be combined with boundary conditions
    //
    // channel_treatment step (8) sets no-slip conditions
    // on wall collocation points.
    //
    // channel_treatment step (9) sets isothermal conditions at walls
    // using rho_wall = e_wall * gamma * (gamma - 1).
    BaseClass::invertMassPlusScaledOperator(
            phi, state, delta_t, substep_index, iota);

    if (constrain_bulk_rhou && has_zero_zero_modes) {

        Map<ArrayXc>(state[ndx::rho ].origin(), Ny).imag().setZero();

        // channel_treatment steps (4), (5), and (6) determine and apply the
        // appropriate bulk momentum forcing to achieve a target value.
        Map<ArrayXc> mean_rhou(state[ndx::rhou].origin(), Ny);
        const complex_t bulk = bulkcoeff.cast<complex_t>().dot(mean_rhou.matrix());
        const real_t varphi = (scenario.bulk_rhou - bulk.real()) / bulk.imag();
        mean_rhou.real() += varphi * mean_rhou.imag();
        common.f() += iota*((varphi/delta_t)*mean_rhou.imag() - common.f());
        mean_rhou.imag().setZero();

        Map<ArrayXc>(state[ndx::rhov].origin(), Ny).imag().setZero();

        Map<ArrayXc>(state[ndx::rhow].origin(), Ny).imag().setZero();

        // channel_treatment step (7) accounts for the momentum forcing
        // within the total energy equation including the Mach squared
        // factor arising from the nondimensionalization choices.
        Map<ArrayXc> mean_rhoe(state[ndx::rhoe].origin(), Ny);
        mean_rhoe.real() += (varphi*scenario.Ma*scenario.Ma) * mean_rhoe.imag();
        common.f_dot_u() += iota*(   (varphi/delta_t)*mean_rhoe.imag()
                                   - common.f_dot_u());
        mean_rhoe.imag().setZero();

        // channel_treatment steps (8) and (9) already performed above
    }

    // B-spline numerics are not conservative when the mean density field has
    // any asymmetry.  Combat the bulk density's very, very small tendency to
    // drift by (usually) dribbling in mass so that we maintain a target bulk
    // density.
    //
    // When nonzero, this mass forcing does perturb the thermodynamic state.  A
    // constant value is added to all density coefficients to avoid introducing
    // any additional density gradient (using that the basis is a partition of
    // unity).  However, examining the equation of state shows that this
    // approach does muck with pressure and temperature gradients (albeit also
    // uniformly).
    //
    // We neglect the associated contributions to the momentum and energy
    // equations as this mass dribbling is small and vanishes when the mean
    // field is symmetric in the wall-normal direction.
    //
    // Enforcing using an integral constraint (a la the bulk momentum target)
    // would be cleaner.  However, it would require additional implicit solves
    // to implement in the general case.  The particulars would be sensitive to
    // boundary conditions in the presence of equation cross-talk.
    if (constrain_bulk_rho && has_zero_zero_modes) {

        // Statistics for the bulk density forcing are not tracked.
        Map<ArrayXr,0,InnerStride<2> > mean_rho(
                (real_t *)state[ndx::rho].origin(), Ny);
        mean_rho *= (scenario.bulk_rho / bulkcoeff.dot(mean_rho.matrix()));
    }

    // No volumetric energy forcing in performed current substep
    common.qb() += iota*(/* (zero/delta_t) */ - common.qb());

    // State leaves method as coefficients in X, Y, and Z directions
}

} // namespace channel

#endif /* CHANNEL_TREATMENT_HPP */
