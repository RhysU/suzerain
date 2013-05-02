//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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

#ifndef SUZERAIN_REACTING_CHANNEL_TREATMENT_HPP
#define SUZERAIN_REACTING_CHANNEL_TREATMENT_HPP

/** @file
 * A mixin providing implicit channel forcing.
 */

#include "nonlinear_operator_fwd.hpp"

#include <suzerain/bspline.hpp>
#include <suzerain/grid_specification.hpp>
#include <suzerain/isothermal_specification.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/storage.hpp>
#include <suzerain/zgbsv_specification.hpp>

#include "reacting.hpp"
#include "antioch_constitutive.hpp"

namespace suzerain {

namespace reacting {

/**
 * A mixin providing channel problem treatment atop any
 * <tt>
 *    public timestepper::nonlinear_operator< contiguous_state<4,complex_t> >
 * </tt>.  During \ref invert_mass_plus_scaled_operator implicit momentum forcing
 * is applied following the section of <tt>writeups/channel_treatment.tex</tt>
 * titled "Enforcing a target bulk momentum via the linear operator" and using
 * information from operator_common_block::u() via instance provided at
 * construction time.
 *
 * Means of the implicit momentum and energy forcing coefficients are also
 * maintained across each individual time step for sampling the statistics \c
 * /bar_f, \c /bar_f_dot_u, \c and /bar_qb using operator_common_block.  Integral
 * constraint means are also tracked for sampling \c /bar_Crho, \c /bar_Crhou,
 * \c /bar_Crhov, \c /bar_Crhow, \c /bar_CrhoE, and \c /bar_Crhou_dot_u.
 */
template< typename BaseClass >
class channel_treatment : public BaseClass
{
public:

    /**
     * Constructor delegating to BaseClass.
     *
     * BaseClass must make its constructor arguments available as member
     * variables under the same name as those found in this constructor.
     */
    channel_treatment(
            const antioch_constitutive& cmods,
            const isothermal_specification &isospec,
            const channel_definition &chdef,
            const grid_specification &grid,
            const pencil_grid &dgrid,
            const bsplineop &cop,
            bspline &b,
            operator_common_block &common);

    /**
     * Constructor delegating to BaseClass.
     *
     * BaseClass must make its constructor arguments available as member
     * variables under the same name as those found in this constructor.
     */
    channel_treatment(
            const zgbsv_specification& spec,
            const antioch_constitutive& cmods,
            const channel_definition &chdef,
            const grid_specification &grid,
            const pencil_grid &dgrid,
            const bsplineop &cop,
            bspline &b,
            operator_common_block &common);

    /** Virtual destructor as the class has virtual methods. */
    virtual ~channel_treatment() {}

    /**
     * Force the channel problem delegating to BaseClass when appropriate.
     * The BaseClass is responsible for enforcing all boundary conditions.
     */
    virtual void invert_mass_plus_scaled_operator(
            const complex_t &phi,
            multi_array::ref<complex_t,4> &state,
            const timestepper::lowstorage::method_interface<complex_t> &method,
            const real_t delta_t,
            const std::size_t substep_index,
            multi_array::ref<complex_t,4> *ic0 = NULL) const;

private:

    /** Common initialization code to be called at end of constructor */
    void finish_construction(
            const grid_specification &grid,
            const pencil_grid &dgrid,
            const bsplineop &cop,
            bspline &b);

    /** Should bulk streamwise momentum constraint be enforced? */
    bool constrain_bulk_rho_u() const
    {
        // Yes when non-inf, non-NaN.
        return (boost::math::isnormal)(this->chdef.bulk_rho_u);
    }

    /** Should bulk density constraint be enforced? */
    bool constrain_bulk_rho() const
    {
        // Yes when non-inf, non-NaN.
        return (boost::math::isnormal)(this->chdef.bulk_rho);
    }

    /** Precomputed integration coefficients */
    VectorXr bulkcoeff;

    /**
     * Constraint data passed to BaseClass.
     *
     * \li cdata.col(0) is for the bulk density constraint.
     * \li cdata.col(1) is for the bulk momentum constraint.
     *
     * Mutable member avoids repeated allocation/deallocation.
     */
    mutable MatrixX2c cdata;

    /**
     * Least squares constraint solver.
     * Mutable member avoids repeated allocation/deallocation.
     */
    mutable Eigen::JacobiSVD<Matrix2r,Eigen::NoQRPreconditioner> jacobiSvd;
};

template< typename BaseClass >
channel_treatment<BaseClass>::channel_treatment(
            const antioch_constitutive& cmods,
            const isothermal_specification &isospec,
            const channel_definition &chdef,
            const grid_specification &grid,
            const pencil_grid &dgrid,
            const bsplineop &cop,
            bspline &b,
            operator_common_block &common)
    : BaseClass(cmods, isospec, chdef, grid, dgrid, cop, b, common),
      jacobiSvd(2, 2, Eigen::ComputeFullU | Eigen::ComputeFullV)
{
    this->finish_construction(grid, dgrid, cop, b);
}

template< typename BaseClass >
channel_treatment<BaseClass>::channel_treatment(
            const zgbsv_specification& spec,
            const antioch_constitutive& cmods,
            const channel_definition &chdef,
            const grid_specification &grid,
            const pencil_grid &dgrid,
            const bsplineop &cop,
            bspline &b,
            operator_common_block &common)
    : BaseClass(spec, cmods, chdef, grid, dgrid, cop, b, common),
      jacobiSvd(2, 2, Eigen::ComputeFullU | Eigen::ComputeFullV)
{
    this->finish_construction(grid, dgrid, cop, b);
}

template< typename BaseClass >
void channel_treatment<BaseClass>::finish_construction(
            const grid_specification &grid,
            const pencil_grid &dgrid,
            const bsplineop &cop,
            bspline &b)
{
    SUZERAIN_UNUSED(cop);

    // Precomputed results only necessary on rank with zero-zero modes
    if (!dgrid.has_zero_zero_modes()) return;

    // Precompute operator for finding bulk quantities from coefficients
    bulkcoeff.resize(b.n());
    b.integration_coefficients(0, bulkcoeff.data());
    bulkcoeff /= grid.L.y();
}

template< typename BaseClass >
void channel_treatment<BaseClass>::invert_mass_plus_scaled_operator(
        const complex_t &phi,
        multi_array::ref<complex_t,4> &state,
        const timestepper::lowstorage::method_interface<complex_t> &method,
        const real_t delta_t,
        const std::size_t substep_index,
        multi_array::ref<complex_t,4> *ic0) const
{
    // State enters method as coefficients in X and Z directions
    // State enters method as collocation point values in Y direction

    SUZERAIN_TIMER_SCOPED("ChannelTreatment");

    // Shorthand
    using std::size_t;
    operator_common_block &common = this->common;
    const channel_definition &chdef = this->chdef;
    const int Ny = this->dgrid.global_wave_extent.y();

    // Sidesteps assertions when local rank contains no wavespace information
    if (SUZERAIN_UNLIKELY(0U == state.shape()[1])) return;

    // Incoming state has wall-normal pencils of interleaved state scalars?
    // Any amount of incoming state is valid so long as there's enough there
    SUZERAIN_ENSURE(state.shape()  [1] == (unsigned) Ny      );
    SUZERAIN_ENSURE(state.strides()[1] ==             1      );
    SUZERAIN_ENSURE(state.strides()[0] == (unsigned) Ny      );
    SUZERAIN_ENSURE(state.shape()  [0] >  (unsigned) ndx::e  );
    SUZERAIN_ENSURE(state.shape()  [0] >  (unsigned) ndx::mx );
    SUZERAIN_ENSURE(state.shape()  [0] >  (unsigned) ndx::mx );
    SUZERAIN_ENSURE(state.shape()  [0] >  (unsigned) ndx::mx );
    SUZERAIN_ENSURE(state.shape()  [0] >  (unsigned) ndx::rho);

    // See channel_treatment writeup (redux) for information on the steps below

    // Have a tantrum if caller expects us to compute any constraints
    SUZERAIN_ENSURE(ic0 == NULL);

    // On zero-zero rank, re-use ic0 to wrap cdata for BaseClass invocation
    if (this->dgrid.has_zero_zero_modes()) {

        // Prepare data for bulk density and bulk momentum constraints
        //
        // Notice scaling by Mach^2 to cause bulk_rho_u-related forcing
        //     work to have nondimensional energy "units" because we will
        //     directly add the result to the total energy equation.
        cdata.setZero(state.shape()[0]*Ny, cdata.cols());
        cdata.col(0).segment(ndx::rho * Ny, Ny).setOnes();
        cdata.col(1).segment(ndx::e   * Ny, Ny).real() = common.u();
        cdata.col(1).segment(ndx::mx  * Ny, Ny).setOnes();

        // Wrap data into appropriately digestible format
        const array<size_t,4> sizes = {{
                state.shape()[0], (size_t) Ny, (size_t) cdata.cols(), 1
        }};
        ic0 = new multi_array::ref<complex_t,4>(
                cdata.data(), sizes, storage::interleaved<4>());
    }

    // Delegate to BaseClass for the solution procedure
    BaseClass::invert_mass_plus_scaled_operator(
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
    crhs <<    chdef.bulk_rho
             - bulkcoeff.dot(mean.segment(ndx::rho * Ny, Ny).real())
         ,     chdef.bulk_rho_u
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
        jacobiSvd.compute(cmat);              // Fancy decomposition for
        cphi = jacobiSvd.solve(crhs);         // a robust 2x2 constraint solve
    } else if (this->constrain_bulk_rho()  ) {
        cphi(0) = crhs(0) / cmat(0,0);        // Solve 1x1 system
        cphi(1) = 0;                          // bulk_rho_u not applied
    } else if (this->constrain_bulk_rho_u()) {
        cphi(0) = 0;                          // bulk_rho not applied
        cphi(1) = crhs(1) / cmat(1,1);        // Solve 1x1 system
    } else {
        cphi.setZero();                       // Neither constraint applied
    }

    // Add scaled constraints to the mean state
    mean += cphi(0)*cdata.col(0) + cphi(1)*cdata.col(1);

    // The implicitly applied integral constraints, as coefficients, must be
    // averaged across each substep to permit accounting for their impact on
    // the Reynolds averaged equations using method_interface::iota similarly to
    //
    //    mean += iota * ((sample / delta_t) - mean).
    //
    // The delta_t accounts for step sizes already implicitly included in cphi.
    //
    // Notice bulk_rho_u-related forcing is NOT scaled by Mach^2 when tracked
    // because our post-processing routines will account for Mach^2 factor.
    const real_t iota     = method.iota(substep_index);
    common.f()           += iota * (
                                ArrayX1r::Constant(Ny, cphi(1) / delta_t)
                              - common.f()
                            );
    common.f_dot_u()     += iota * (
                                (cphi(1) / delta_t) * common.u()
                              - common.f_dot_u()
                            );
    common.qb()          += iota * (/* zero */ - common.qb());
    common.CrhoE()       += iota * (/* zero */ - common.CrhoE());
    common.Crhou()       += iota * (/* zero */ - common.Crhou());
    common.Crhov()       += iota * (/* zero */ - common.Crhov());
    common.Crhow()       += iota * (/* zero */ - common.Crhow());
    common.Crho()        += iota * (
                                ArrayX1r::Constant(Ny, cphi(0) / delta_t)
                              - common.Crho()
                            );
    common.Crhou_dot_u() += iota * (/* zero */ - common.Crhou_dot_u());

    // State leaves method as coefficients in X, Y, and Z directions
}

} // namespace reacting

} // namespace suzerain

#endif /* SUZERAIN_REACTING_CHANNEL_TREATMENT_HPP */
