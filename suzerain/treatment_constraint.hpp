//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
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

#ifndef SUZERAIN_TREATMENT_CONSTRAINT_HPP
#define SUZERAIN_TREATMENT_CONSTRAINT_HPP

/** @file
 * Provides \ref treatment_constraint.
 */

#include <suzerain/bspline.h>
#include <suzerain/common.hpp>
#include <suzerain/constraint.hpp>
#include <suzerain/lowstorage.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/state_fwd.hpp>

namespace suzerain {

namespace constraint {

//// TODO? Use Eigen Maps to decouple operator_common_block if sensible

/**
 * A wrapper applying integral constraint treatment atop any linear operator.
 * During \ref invert_mass_plus_scaled_operator implicit momentum forcing is
 * applied following the section of <tt>writeups/treatment_channel.tex</tt>
 * titled "Enforcing a target bulk momentum via the linear operator" and using
 * information from \ref operator_common_block via an instance provided at
 * construction time.
 */
template<typename CommonBlock>
class treatment
    : public lowstorage::linear_operator<
          multi_array::ref<complex_t,4>,
          contiguous_state<4,complex_t>
      >
    , public boost::noncopyable
{

public:

    /**
    * Abstract interface for tracking instantaneous collocation point
    * profiles taken to be inputs to constraint treatment logic.
    */
    class inputs
    {
    public:
        virtual ~inputs() {} ///< Virtual for interface
        virtual ArrayXXr::ConstColXpr u()  const = 0; ///< Streamwise velocity
        virtual ArrayXXr::ConstColXpr v()  const = 0; ///< Wall-normal velocity
        virtual ArrayXXr::ConstColXpr w()  const = 0; ///< Spanwise velocity
        virtual ArrayXXr::ConstColXpr uu() const = 0; ///< Reynolds stress component uu
        virtual ArrayXXr::ConstColXpr uv() const = 0; ///< Reynolds stress component uv
        virtual ArrayXXr::ConstColXpr uw() const = 0; ///< Reynolds stress component uw
        virtual ArrayXXr::ConstColXpr vv() const = 0; ///< Reynolds stress component vv
        virtual ArrayXXr::ConstColXpr vw() const = 0; ///< Reynolds stress component vw
        virtual ArrayXXr::ConstColXpr ww() const = 0; ///< Reynolds stress component ww
    };

    /**
    * Abstract interface for tracking forcing collocation point profile
    * outputs maintained by constraint treatment logic.
    */
    class outputs
    {
    public:
        virtual ~outputs() {}                        ///< Virtual for interface
        virtual ArrayXXr::ColXpr fx()           = 0; ///< Streamwise momentum
        virtual ArrayXXr::ColXpr fy()           = 0; ///< Wall-normal momentum
        virtual ArrayXXr::ColXpr fz()           = 0; ///< Spanwise momentum
        virtual ArrayXXr::ColXpr f_dot_u()      = 0; ///< Forcing work
        virtual ArrayXXr::ColXpr qb()           = 0; ///< Heating
        virtual ArrayXXr::ColXpr CrhoE()        = 0; ///< Constrained total energy
        virtual ArrayXXr::ColXpr Crhou()        = 0; ///< Constrained streamwise momentum
        virtual ArrayXXr::ColXpr Crhov()        = 0; ///< Constrained wall-normal momentum
        virtual ArrayXXr::ColXpr Crhow()        = 0; ///< Constrained spanwise momentum
        virtual ArrayXXr::ColXpr Crhou_dot_u()  = 0; ///< Constrained forcing work
        virtual ArrayXXr::ColXpr Crho()         = 0; ///< Constrained density
        virtual ArrayXXr::ColXpr C2rhoE()       = 0; ///< Squared \ref CrhoE
        virtual ArrayXXr::ColXpr C2rhou()       = 0; ///< Squared \ref Crhou
        virtual ArrayXXr::ColXpr C2rhov()       = 0; ///< Squared \ref Crhov
        virtual ArrayXXr::ColXpr C2rhow()       = 0; ///< Squared \ref Crhow
        virtual ArrayXXr::ColXpr C2rhou_dot_u() = 0; ///< Squared \ref Crhou_dot_u
        virtual ArrayXXr::ColXpr C2rho()        = 0; ///< Squared \ref Crho
    };

    /**
     * Constructor.  After construction, #L must be provided.  Generally, one
     * or more constraints will be supplied via #physical or #numerical before
     * usage.
     *
     * @param Ma     Nondimensional Mach number for kinetic energy computations.
     *               In a dimensional setting, this should be one.  Notice that
     *               it is a \e reference permitting the value to track another
     *               setting.  For example, one in a \ref definition_scenario.
     * @param dgrid  Parallel decomposition details used to determine which
     *               rank houses the "zero-zero" Fourier modes.
     * @param common Storage from which mean velocity profiles will be read and
     *               to which implicit forcing averages will be accumulated.
     */
    treatment(
            const real_t& Ma,
            const pencil_grid& dgrid,
            bspline& b,
            CommonBlock& common);

    /** Delegates invocation to #L */
    virtual void apply_mass_plus_scaled_operator(
             const complex_t &phi,
             multi_array::ref<complex_t,4> &state,
             const std::size_t substep_index) const;

    /** Delegates invocation to #L */
    virtual void accumulate_mass_plus_scaled_operator(
             const complex_t &phi,
             const multi_array::ref<complex_t,4> &input,
             const complex_t &beta,
             contiguous_state<4,complex_t> &output,
             const std::size_t substep_index) const;

    /**
     * Force the channel problem delegating to #L when appropriate.
     * #L is responsible for enforcing all boundary conditions.
     */
    virtual void invert_mass_plus_scaled_operator(
            const complex_t& phi,
            multi_array::ref<complex_t,4>& state,
            const lowstorage::method_interface<complex_t>& method,
            const real_t delta_t,
            const std::size_t substep_index,
            multi_array::ref<complex_t,4> *ic0 = NULL) const;

    /** The operator whose behavior is modified by this instance. */
    shared_ptr<lowstorage::linear_operator<
                multi_array::ref<complex_t,4>,
                contiguous_state<4,complex_t>
            > > L;

    /**
     * Catalog of physically-oriented constraints to be indexed by equation
     * number ndx::e, ndx::mx, ndx::my, ndx::mz, or ndx::rho.  These are
     * constraints which interpretable as physics-related source terms which
     * may include contributions to multiple equations.  For example, an
     * ndx::mx momentum constraint will also contribute forcing work terms to
     * the ndx::e total energy.
     *
     * \warning The ndx::rho entry is not treated "physically" in the sense
     * that density added by the constraint updates neither the momentum nor
     * total energy equations.  Therefore, when the physical[ndx::rho]
     * constraint modifies the density it also adjusts the velocity and
     * specific total energy.
     *
     * Means of the implicit ndx::mx, ndx::my, ndx::mz, and ndx::e forcing are
     * maintained across each individual time step for sampling the statistics
     * \c /bar_f, \c /bar_f_dot_u, and \c /bar_qb using #com.  The mean
     * implicit numerical[ndx::rho] forcing is combined with
     * physical[ndx::rho]'s effect and stored within \c /bar_Crho.
     */
    array<shared_ptr<constraint::base>, 5> physical;

    /**
     * Catalog of numerically-oriented constraints to be indexed by equation
     * number ndx::e, ndx::mx, ndx::my, ndx::mz, or ndx::rho.  These are
     * numerical implementation artifacts to control the behavior of single
     * equations rather than physics-related sources.  For example, an
     * ndx::mx momentum constraint will not not contribute forcing work terms
     * to the ndx::e total energy.
     *
     * Means of the implicit ndx::mx, ndx::my, ndx::mz, and ndx::e forcing are
     * maintained across each individual time step for sampling the statistics
     * \c /bar_Crhou, \c /bar_Crhou_dot_u, and \c /bar_CrhoE using #com.
     * The mean implicit numerical[ndx::rho] forcing is combined with
     * physical[ndx::rho]'s effect and stored within \c /bar_Crho.
     */
    array<shared_ptr<constraint::base>, 5> numerical;

    /** An appropriately-size, do-nothing constraint usable by callers. */
    const shared_ptr<constraint::disabled> none;

protected:

    /**
     * What Mach number should be used to scale kinetic
     * energy contributions to the total energy equation?
     */
    const real_t& Ma;

    /** Does this rank possess the "zero-zero" Fourier modes? */
    const bool rank_has_zero_zero_modes;

    /**
     * Used to obtain mean primitive state profiles to compute
     * implicit forcing work contributions to the total energy equation.
     */
    CommonBlock& com;

private:

    /**
     * Constraint data passed to #L indexed by ndx::type.  Entries <tt>[ndx::e,
     * ndx::rho]</tt> pertain to #physical.  Entries <tt>[ndx::rho_(1) +
     * ndx::e, ndx::rho_(1) + ndx::rho]</tt> pertain to #numerical.  Mutable
     * member avoids repeated allocation/deallocation.
     */
    mutable MatrixXAc cdata;

    /**
     * Least squares constraint solver.
     * Mutable member avoids repeated allocation/deallocation.
     */
    mutable Eigen::JacobiSVD<
            MatrixXXr, Eigen::FullPivHouseholderQRPreconditioner
        > jacobiSvd;

};

template <typename CommonBlock>
treatment<CommonBlock>::treatment(
            const real_t& Ma,
            const pencil_grid& dgrid,
            bspline& b,
            CommonBlock& common)
    : none(new constraint::disabled(b))
    , Ma(Ma)
    , rank_has_zero_zero_modes(dgrid.has_zero_zero_modes())
    , com(common)
    , jacobiSvd(0, 0, Eigen::ComputeFullU | Eigen::ComputeFullV)
{
    using std::fill;
    fill(physical.begin(),  physical.end(),  none);
    fill(numerical.begin(), numerical.end(), none);
}

template <typename CommonBlock>
void treatment<CommonBlock>::apply_mass_plus_scaled_operator(
        const complex_t &phi,
        multi_array::ref<complex_t,4> &state,
        const std::size_t substep_index) const
{
    return L->apply_mass_plus_scaled_operator(
            phi, state, substep_index);
}

template <typename CommonBlock>
void
treatment<CommonBlock>::accumulate_mass_plus_scaled_operator(
        const complex_t &phi,
        const multi_array::ref<complex_t,4> &input,
        const complex_t &beta,
        contiguous_state<4,complex_t> &output,
        const std::size_t substep_index) const
{
    return L->accumulate_mass_plus_scaled_operator(
            phi, input, beta, output, substep_index);
}

template <typename CommonBlock>
void
treatment<CommonBlock>::invert_mass_plus_scaled_operator(
        const complex_t& phi,
        multi_array::ref<complex_t,4>& state,
        const lowstorage::method_interface<complex_t>& method,
        const real_t delta_t,
        const std::size_t substep_index,
        multi_array::ref<complex_t,4> *ic0) const
{
    // State enters method as coefficients in X and Z directions
    // State enters method as collocation point values in Y direction

    SUZERAIN_TIMER_SCOPED("constraint::treatment");

    // Shorthand
    using ndx::e;
    using ndx::mx;
    using ndx::my;
    using ndx::mz;
    using ndx::rho;
    const int rho_1 = ndx::rho_(1);

    // Sidesteps assertions when local rank contains no wavespace information
    const int N = (int) state.shape()[1];
    if (SUZERAIN_UNLIKELY(0 == N)) return;

    // Incoming state has wall-normal pencils of interleaved state scalars?
    // Any amount of incoming state is valid so long as there's enough there
    SUZERAIN_ENSURE(state.strides()[1] ==            1  );
    SUZERAIN_ENSURE(state.strides()[0] == (unsigned) N  );
    SUZERAIN_ENSURE(state.shape()  [0] >  (unsigned) e  );
    SUZERAIN_ENSURE(state.shape()  [0] >  (unsigned) mx );
    SUZERAIN_ENSURE(state.shape()  [0] >  (unsigned) mx );
    SUZERAIN_ENSURE(state.shape()  [0] >  (unsigned) mx );
    SUZERAIN_ENSURE(state.shape()  [0] >  (unsigned) rho);

    // See treatment_channel writeup (redux) for information on steps below.

    // Have a tantrum if caller expects us to compute any constraints
    SUZERAIN_ENSURE(ic0 == NULL);

    // On zero-zero rank, re-use ic0 to wrap cdata for BaseClass invocation
    if (rank_has_zero_zero_modes) {

        // Prepare RHS data for density, momentum, and total energy constraints
        //
        // Notice scaling by Mach^2 to cause mx-, my-, and mz-related forcing
        //     work to have nondimensional energy "units" because we will
        //     directly add the result to the total energy equation.
        //
        // Notice ndx::rho constraint ignores other equations BY DESIGN.
        // This choice is reflected in the later update of com.Crho.
        cdata.setZero(state.shape()[0]*N, cdata.cols());
        cdata.col(e  ).segment(N*e  , N).real() = physical[e  ]->shape;
        cdata.col(mx ).segment(N*e  , N).real() = physical[mx ]->shape
                                                * com.u() * (Ma * Ma);
        cdata.col(mx ).segment(N*mx , N).real() = physical[mx ]->shape;
        cdata.col(my ).segment(N*e  , N).real() = physical[my ]->shape
                                                * com.v() * (Ma * Ma);
        cdata.col(my ).segment(N*my , N).real() = physical[my ]->shape;
        cdata.col(mz ).segment(N*e  , N).real() = physical[mz ]->shape
                                                * com.w() * (Ma * Ma);
        cdata.col(mz ).segment(N*mz , N).real() = physical[mz ]->shape;
        cdata.col(rho).segment(N*rho, N).real() = physical[rho]->shape;

        // Prepare RHS for one decoupled numerical constraint per equation
        cdata.col(rho_1+e  ).segment(N*e  , N).real() = numerical[e  ]->shape;
        cdata.col(rho_1+mx ).segment(N*mx , N).real() = numerical[mx ]->shape;
        cdata.col(rho_1+my ).segment(N*my , N).real() = numerical[my ]->shape;
        cdata.col(rho_1+mz ).segment(N*mz , N).real() = numerical[mz ]->shape;
        cdata.col(rho_1+rho).segment(N*rho, N).real() = numerical[rho]->shape;

        // Wrap data into appropriately digestible format
        using std::size_t;
        const array<size_t,4> sizes = {{
                state.shape()[0], (size_t) N, (size_t) cdata.cols(), 1
        }};
        ic0 = new multi_array::ref<complex_t,4>(
                cdata.data(), sizes, storage::interleaved<4>());
    }

    // Delegate to wrapped operator for the solution procedure
    L->invert_mass_plus_scaled_operator(
            phi, state, method, delta_t, substep_index, ic0);

    // Clean up any integral constraint data wrapper we may have employed
    delete ic0;

    // Only the rank with zero-zero modes proceeds!
    if (!rank_has_zero_zero_modes) return;

    // Get an Eigen-friendly map of the zero-zero mode coefficients
    Map<VectorXc> mean(state.origin(), state.shape()[0]*N);

    // Solve the requested, possibly simultaneous constraint problem.  A fancy
    // decomposition is used for a simple 10x10 solve or to permit one or more
    // inactive constraints via least squares.  Least squares also adds
    // robustness whenever the constraints are somehow incompatible.
    //
    // 1) Assemble the matrix problem for simultaneous integral constraints.
    //    Nonzero only when both the i-th and j-th constraints are enabled.
    MatrixAr cmat(MatrixAr::Zero());
    for (int j = 0; j < cdata.cols(); ++j) {
        if (j<rho_1 ? physical[j]->enabled() : numerical[j-rho_1]->enabled()) {
            for (int i = 0; i < physical.static_size;  ++i) {
                if (physical[i]->enabled()) {
                    cmat(i,j) = physical[i]->coeff.dot(
                            cdata.col(j).segment(N*i,N).real());
                }
            }
            for (int i = 0; i < numerical.static_size; ++i) {
                if (numerical[i]->enabled()) {
                    cmat(rho_1+i,j) = numerical[i]->coeff.dot(
                            cdata.col(j).segment(N*i,N).real());
                }
            }
        }
    }
    // 2) Prepare initial values relative to desired targets when active.
    //    Otherwise, zero the associated row and column in matrix.
    VectorAr crhs(VectorAr::Zero());
    for (int i = 0; i < physical.static_size;  ++i) {
        if (physical[i]->enabled()) {
            crhs[i] = physical[i]->target()
                    - physical[i]->coeff.dot(
                            mean.segment(N*i,N).real());
        }
    }
    for (int i = 0; i < numerical.static_size;  ++i) {
        if (numerical[i]->enabled()) {
            crhs[rho_1+i] = numerical[i]->target()
                          - numerical[i]->coeff.dot(
                                  mean.segment(N*i,N).real());
        }
    }
    // 3) Prepare matrix decomposition and solve using least squares
    VectorAr cphi = jacobiSvd.compute(cmat).solve(crhs);
    // 4) Add correctly scaled results to the mean state to satisfy constraints
    mean += (cdata * cphi.asDiagonal()).rowwise().sum();

    // The implicitly applied integral constraints, as point values, must be
    // averaged across each substep to permit accounting for their impact on
    // the Reynolds averaged equations using method_interface::iota as in
    //
    //    mean += iota * ((sample / delta_t) - mean).
    //
    // The delta_t accounts for step sizes already implicitly included in cphi.
    //
    // Notice mx-related forcing is NOT scaled by Mach^2 when tracked
    // because our post-processing routines will account for Mach^2 factor.
    //
    // Notice physical[ndx::rho] constraint lumped into Crho BY DESIGN.
    cphi               /= delta_t; // Henceforth includes 1/delta_t scaling!
    const real_t iota   = method.iota(substep_index);
    com.fx()           += iota*(
                              cphi[mx]*physical[mx]->shape
                            - com.fx()
                          );
    com.fy()           += iota*(
                              cphi[my]*physical[my]->shape
                            - com.fy()
                          );
    com.fz()           += iota*(
                              cphi[mz]*physical[mz]->shape
                            - com.fz()
                          );
    com.f_dot_u()      += iota*(
                              cphi[mx]*physical[mx]->shape*com.u()
                            + cphi[my]*physical[my]->shape*com.v()
                            + cphi[mz]*physical[mz]->shape*com.w()
                            - com.f_dot_u()
                          );
    com.qb()           += iota*(
                              cphi[e]*physical[e]->shape
                            - com.qb()
                          );
    com.CrhoE()        += iota*(
                              cphi[rho_1+e  ]*numerical[e  ]->shape
                            - com.CrhoE()
                          );
    com.C2rhoE()       += iota*(
                              (cphi[rho_1+e  ]*numerical[e  ]->shape).square()
                            - com.C2rhoE()
                          );
    com.Crhou()        += iota*(
                              cphi[rho_1+mx ]*numerical[mx ]->shape
                            - com.Crhou()
                          );
    com.C2rhou()       += iota*(
                              (cphi[rho_1+mx ]*numerical[mx ]->shape).square()
                            - com.C2rhou()
                          );
    com.Crhov()        += iota*(
                              cphi[rho_1+my ]*numerical[my ]->shape
                            - com.Crhov()
                          );
    com.C2rhov()       += iota*(
                              (cphi[rho_1+my ]*numerical[my ]->shape).square()
                            - com.C2rhov()
                          );
    com.Crhow()        += iota*(
                              cphi[rho_1+mz ]*numerical[mz ]->shape
                            - com.Crhow()
                          );
    com.C2rhow()       += iota*(
                              (cphi[rho_1+mz ]*numerical[mz ]->shape).square()
                            - com.C2rhow()
                          );
    com.Crho()         += iota*(
                              cphi[rho]      *physical [rho]->shape
                            + cphi[rho_1+rho]*numerical[rho]->shape
                            - com.Crho()
                          );
    com.C2rho()        += iota*( // Terms like aa + 2ab + bb
                                (cphi[rho]      *physical [rho]->shape).square()
                            + 2*(cphi[rho]      *physical [rho]->shape)
                               *(cphi[rho_1+rho]*numerical[rho]->shape)
                            +   (cphi[rho_1+rho]*numerical[rho]->shape).square()
                            - com.C2rho()
                          );
    com.Crhou_dot_u()  += iota*(
                              cphi[rho_1+mx]*numerical[mx]->shape*com.u()
                            + cphi[rho_1+my]*numerical[my]->shape*com.v()
                            + cphi[rho_1+mz]*numerical[mz]->shape*com.w()
                            - com.Crhou_dot_u()
                          );
    com.C2rhou_dot_u() += iota*( // Terms like uu + 2uv + 2uw + vv + 2vw + ww
                                (cphi[rho_1+mx]*numerical[mx]->shape).square()
                               *com.uu()
                            + 2*(cphi[rho_1+mx]*cphi[rho_1+my])
                               *numerical[mx]->shape
                               *numerical[my]->shape
                               *com.uv()
                            + 2*(cphi[rho_1+mx]*cphi[rho_1+mz])
                               *numerical[mx]->shape
                               *numerical[mz]->shape
                               *com.uw()
                            +   (cphi[rho_1+my]*numerical[my]->shape).square()
                               *com.vv()
                            + 2*(cphi[rho_1+my]*cphi[rho_1+mz])
                               *numerical[my]->shape
                               *numerical[mz]->shape
                               *com.vw()
                            +   (cphi[rho_1+mz]*numerical[mz]->shape).square()
                               *com.ww()
                            - com.C2rhou_dot_u()
                          );

    // State leaves method as coefficients in X, Y, and Z directions
}

} // namespace constraint

} // namespace suzerain

#endif /* SUZERAIN_TREATMENT_CONSTRAINT_HPP */
