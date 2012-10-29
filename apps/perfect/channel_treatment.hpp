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
#include <suzerain/spec_zgbsv.hpp>
#include <suzerain/storage.hpp>

#include "perfect.hpp"

namespace suzerain { namespace perfect {

/**
 * A mixin providing channel problem treatment atop any
 * <tt>
 *    public timestepper::INonlinearOperator< ContiguousState<4,complex_t> >
 * </tt>.  During \ref invertMassPlusScaledOperator implicit momentum forcing
 * is applied following the section of <tt>writeups/channel_treatment.tex</tt>
 * titled "Enforcing a target bulk momentum via the linear operator" and using
 * information from OperatorCommonBlock::u() via instance provided at
 * construction time.
 *
 * Means of the implicit momentum and energy forcing coefficients are also
 * maintained across each individual time step for sampling the statistics \c
 * /bar_f, \c /bar_f_dot_u, \c and /bar_qb using OperatorCommonBlock.  Integral
 * constraint means are also tracked for sampling \c /bar_Crho, \c /bar_Crhou,
 * \c /bar_Crhov, \c /bar_Crhow, \c /bar_CrhoE, and \c /bar_Crhou_dot_u.
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
            const ScenarioDefinition &scenario,
            const problem::GridDefinition &grid,
            const pencil_grid &dgrid,
            bspline &b,
            const bsplineop &bop,
            OperatorCommonBlock &common);

    /**
     * Constructor delegating to BaseClass.
     *
     * BaseClass must make its constructor arguments available as member
     * variables under the same name as those found in this constructor.
     */
    ChannelTreatment(
            const spec_zgbsv& spec,
            const ScenarioDefinition &scenario,
            const problem::GridDefinition &grid,
            const pencil_grid &dgrid,
            bspline &b,
            const bsplineop &bop,
            OperatorCommonBlock &common);

    /**
     * Force the channel problem delegating to BaseClass when appropriate.
     * The BaseClass is responsible for enforcing all boundary conditions.
     */
    virtual void invertMassPlusScaledOperator(
            const complex_t &phi,
            multi_array::ref<complex_t,4> &state,
            const timestepper::lowstorage::IMethod<complex_t> &method,
            const real_t delta_t,
            const std::size_t substep_index,
            multi_array::ref<complex_t,4> *ic0 = NULL) const;

private:

    /** Common initialization code to be called at end of constructor */
    void finish_construction(
            const problem::GridDefinition &grid,
            const pencil_grid &dgrid,
            bspline &b,
            const bsplineop &bop);

    /** Should bulk streamwise momentum constraint be enforced? */
    bool constrain_bulk_rho_u() const
    {
        // Yes when non-inf, non-NaN.
        return (boost::math::isnormal)(this->scenario.bulk_rho_u);
    }

    /** Should bulk density constraint be enforced? */
    bool constrain_bulk_rho() const
    {
        // Yes when non-inf, non-NaN.
        return (boost::math::isnormal)(this->scenario.bulk_rho);
    }

    /** Precomputed integration coefficients */
    VectorXr bulkcoeff;
};

template< typename BaseClass >
ChannelTreatment<BaseClass>::ChannelTreatment(
            const ScenarioDefinition &scenario,
            const problem::GridDefinition &grid,
            const pencil_grid &dgrid,
            bspline &b,
            const bsplineop &bop,
            OperatorCommonBlock &common)
    : BaseClass(scenario, grid, dgrid, b, bop, common)
{
    this->finish_construction(grid, dgrid, b, bop);
}

template< typename BaseClass >
ChannelTreatment<BaseClass>::ChannelTreatment(
            const spec_zgbsv& spec,
            const ScenarioDefinition &scenario,
            const problem::GridDefinition &grid,
            const pencil_grid &dgrid,
            bspline &b,
            const bsplineop &bop,
            OperatorCommonBlock &common)
    : BaseClass(spec, scenario, grid, dgrid, b, bop, common)
{
    this->finish_construction(grid, dgrid, b, bop);
}

template< typename BaseClass >
void ChannelTreatment<BaseClass>::finish_construction(
            const problem::GridDefinition &grid,
            const pencil_grid &dgrid,
            bspline &b,
            const bsplineop &bop)
{
    SUZERAIN_UNUSED(bop);

    // Precomputed results only necessary on rank with zero-zero modes
    if (!dgrid.has_zero_zero_modes()) return;

    // Precompute operator for finding bulk quantities from coefficients
    bulkcoeff.resize(b.n());
    b.integration_coefficients(0, bulkcoeff.data());
    bulkcoeff /= grid.L.y();
}

template< typename BaseClass >
void ChannelTreatment<BaseClass>::invertMassPlusScaledOperator(
        const complex_t &phi,
        multi_array::ref<complex_t,4> &state,
        const timestepper::lowstorage::IMethod<complex_t> &method,
        const real_t delta_t,
        const std::size_t substep_index,
        multi_array::ref<complex_t,4> *ic0) const
{
    // State enters method as coefficients in X and Z directions
    // State enters method as collocation point values in Y direction

    SUZERAIN_TIMER_SCOPED("ChannelTreatment");

    // Shorthand
    namespace ndx = support::field::ndx;
    OperatorCommonBlock &common = this->common;
    const ScenarioDefinition &scenario = this->scenario;
    const int Ny = this->dgrid.global_wave_extent.y();

    // Sidesteps assertions when local rank contains no wavespace information
    if (SUZERAIN_UNLIKELY(0U == state.shape()[1])) return;

    // Incoming state has wall-normal pencils of interleaved state scalars?
    SUZERAIN_ENSURE(state.shape()  [1] == (unsigned) Ny);
    SUZERAIN_ENSURE(state.strides()[1] ==             1);
    SUZERAIN_ENSURE(state.strides()[0] == (unsigned) Ny);
    SUZERAIN_ENSURE(state.shape()  [0] == support::field::count);

    // See channel_treatment writeup (redux) for information on the steps below

    // Have a tantrum if caller expects us to compute any constraints
    SUZERAIN_ENSURE(ic0 == NULL);

    // Constraint data wrapped is kept within cdata:
    //   cdata.col(0) is for the bulk density constraint
    //   cdata.col(1) is for the bulk momentum constraint
    //
    // On zero-zero rank, re-use ic0 to wrap data for BaseClass invocation
    MatrixX2c cdata;
    if (this->dgrid.has_zero_zero_modes()) {

        // Prepare data for bulk density and bulk momentum constraints
        // Nondimensionalization requires scaling forcing work by Mach^2
        cdata.setZero(state.shape()[0]*Ny, cdata.cols());
        cdata.col(0).segment(ndx::rho * Ny, Ny).setOnes();
        cdata.col(1).segment(ndx::e   * Ny, Ny).real() = scenario.Ma
                                                       * scenario.Ma
                                                       * common.u();
        cdata.col(1).segment(ndx::mx  * Ny, Ny).setOnes();

        // Wrap data into appropriately digestible format
        const boost::array<std::size_t,4> sizes
            = {{ state.shape()[0], Ny, cdata.cols(), 1 }};
        ic0 = new multi_array::ref<complex_t,4>(
                cdata.data(), sizes, storage::interleaved<4>());
    }

    // Delegate to BaseClass for the solution procedure
    BaseClass::invertMassPlusScaledOperator(
            phi, state, method, delta_t, substep_index, ic0);

    // Clean up any integral constraint data wrapper we may have employed
    delete ic0;

    // Only the rank with zero-zero modes proceeds!
    if (!this->dgrid.has_zero_zero_modes()) return;

    // Get an Eigen-friendly map of the zero-zero mode coefficients
    Map<VectorXc> mean(state.origin(), state.shape()[0]*Ny);

    // Assemble the matrix problem for simultaneous integral constraints
    // Notice Matrix2r::operator<< depicts the matrix in row-major fashion
    Vector2r crhs;
    crhs <<    scenario.bulk_rho
             - bulkcoeff.dot(mean.segment(ndx::rho * Ny, Ny).real())
         ,     scenario.bulk_rho_u
             - bulkcoeff.dot(mean.segment(ndx::mx  * Ny, Ny).real())
         ;
    Matrix2r cmat;
    cmat << bulkcoeff.dot(cdata.col(0).segment(ndx::rho * Ny, Ny).real())
         ,  bulkcoeff.dot(cdata.col(1).segment(ndx::rho * Ny, Ny).real())
         ,  bulkcoeff.dot(cdata.col(0).segment(ndx::mx  * Ny, Ny).real())
         ,  bulkcoeff.dot(cdata.col(1).segment(ndx::mx  * Ny, Ny).real())
         ;
    Vector2r cphi;

    // Solve the requested, possibly simultaneous constraint problem
    if (this->constrain_bulk_rho() && this->constrain_bulk_rho_u()) {
        enum {options = Eigen::ComputeFullU | Eigen::ComputeFullV};
        cphi = cmat.jacobiSvd(options).solve(crhs);
    } else if (this->constrain_bulk_rho()  ) {
        cphi(0) = crhs(0) / cmat(0,0);        // Solve 1x1 system
        cphi(1) = 0;                          // bulk_rho_u not applied
    } else if (this->constrain_bulk_rho_u()) {
        cphi(0) = 0;                          // bulk_rho not applied
        cphi(1) = crhs(1) / cmat(1,1);        // Solve 1x1 system
    } else {
        cphi.setZero();                       // Neither constraint applied
    }

    // Mutate constraint data by scaling and then add to the mean state
    cdata.col(0) *= cphi(0);
    mean += cdata.col(0);
    cdata.col(1) *= cphi(1);
    mean += cdata.col(1);

    // The implicitly applied integral constraints, as coefficients, must be
    // averaged across each substep to permit accounting for their impact on
    // the Reynolds averaged equations using IMethod::iota_alpha similar to
    //
    //    mean += iota_alpha*(sample/(alpha*delta_t) - mean).
    //
    // The alpha*delta_t factor accounts for time step sizes.
    const real_t iota_alpha   = method.iota_alpha(substep_index);
    const real_t inv_alpha_dt = 1 / (method.alpha(substep_index)*delta_t);

    // Fully-coupled implicit solves can cause either constraint, when coupled
    // with boundary conditions, to impact every equation and so we must
    // accommodate that unpleasant possibility.
    //
    // Separately tracking the bulk_rho_u constraint as a
    // pressure-gradient-like forcing f doing work f_dot_u is desirable.
    // All other impacts are lumped into Crho{,u,v,w,E,u_dot_u}.  That is,
    //
    //     -----   --------Constraint--------
    //      Eqn    bulk_rho  bulk_rho_u  zero
    //     -----   --------------------------
    //     rho_E   CrhoE     f_dot_u
    //     rho_u   Crhou     f
    //     rho_v   Crhov     Crhov
    //     rho_w   Crhow     Crhow
    //     rho     Crho      Crho        qb
    //     -----   --------------------------
    //
    // The logic is a bit convoluted to avoid introducing temporaries.  Sorry.

    // First, track physically-oriented forcing
    common.f()       += iota_alpha*(
                            inv_alpha_dt*cdata.col(1).segment(ndx::mx * Ny, Ny)
                                                     .real().array()
                          - common.f()
                        );
    cdata.col(1).segment(ndx::mx * Ny, Ny).setZero();  // Mark done
    common.f_dot_u() += iota_alpha*(
                            inv_alpha_dt*cdata.col(1).segment(ndx::e  * Ny, Ny)
                                                     .real().array()
                          - common.f_dot_u()
                        );
    cdata.col(1).segment(ndx::e  * Ny, Ny).setZero();  // Mark done
    common.qb()      += iota_alpha*(
                            /* inv_alpha_dt * zero */  // No bulk heating
                          - common.qb()
                        );

    // Second, merge and track any remaining numerically-oriented quantities
    cdata.col(0) = cdata.rowwise().sum();
    common.CrhoE() += iota_alpha*(
                          inv_alpha_dt*cdata.col(0).segment(ndx::e  * Ny, Ny)
                                                   .real().array()
                        - common.CrhoE()
                      );
    common.Crhou() += iota_alpha*(
                          inv_alpha_dt*cdata.col(0).segment(ndx::mx * Ny, Ny)
                                                   .real().array()
                        - common.Crhou()
                      );
    common.Crhov() += iota_alpha*(
                          inv_alpha_dt*cdata.col(0).segment(ndx::my * Ny, Ny)
                                                   .real().array()
                        - common.Crhov()
                      );
    common.Crhow() += iota_alpha*(
                          inv_alpha_dt*cdata.col(0).segment(ndx::mz * Ny, Ny)
                                                   .real().array()
                        - common.Crhow()
                      );
    common.Crho()  += iota_alpha*(
                          inv_alpha_dt*cdata.col(0).segment(ndx::rho* Ny, Ny)
                                                   .real().array()
                        - common.Crho()
                      );

    // Last, track the impact of Crho{u,v,w} on the TKE equation
    // Nondimensionalization requires scaling constraint work by Mach^2
    common.Crhou_dot_u() += iota_alpha*(
                              (inv_alpha_dt*scenario.Ma*scenario.Ma)
                                *cdata.col(0).segment(ndx::mx*Ny, Ny)
                                             .real().array()
                                *common.u()
                            + (inv_alpha_dt*scenario.Ma*scenario.Ma)
                                *cdata.col(0).segment(ndx::my*Ny, Ny)
                                             .real().array()
                                *common.v()
                            + (inv_alpha_dt*scenario.Ma*scenario.Ma)
                                *cdata.col(0).segment(ndx::mz*Ny, Ny)
                                             .real().array()
                                *common.w()
                            - common.Crho()
                          );

    // State leaves method as coefficients in X, Y, and Z directions
}

} /* namespace perfect */ } /* namespace suzerain */

#endif /* CHANNEL_TREATMENT_HPP */
