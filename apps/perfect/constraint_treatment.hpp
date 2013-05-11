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

#ifndef SUZERAIN_PERFECT_CONSTRAINT_TREATMENT_HPP
#define SUZERAIN_PERFECT_CONSTRAINT_TREATMENT_HPP

/** @file
 * Provides \ref constraint_treatment.
 */

#include <suzerain/common.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/state_fwd.hpp>
#include <suzerain/timestepper.hpp>

namespace suzerain {

// Forward declarations
class bspline;
class bsplineop;
class grid_specification;
class pencil_grid;

namespace perfect {

// Forward declarations
class operator_common_block;
class scenario_definition;

/**
 * A wrapper applying integral constraints treatment atop any linear operator.
 * This class began as a way to provide integral constraints for driving a
 * Coleman-like channel.
 *
 * During \ref invert_mass_plus_scaled_operator implicit momentum forcing is
 * applied following the section of <tt>writeups/channel_treatment.tex</tt>
 * titled "Enforcing a target bulk momentum via the linear operator" and using
 * information from operator_common_block::u() via an instance provided at
 * construction time.
 *
 * Means of the implicit momentum and energy forcing coefficients are also
 * maintained across each individual time step for sampling the statistics \c
 * /bar_f, \c /bar_f_dot_u, \c and /bar_qb using operator_common_block.
 * Integral constraint means are also tracked for sampling \c /bar_Crho, \c
 * /bar_Crhou, \c /bar_Crhov, \c /bar_Crhow, \c /bar_CrhoE, and \c
 * /bar_Crhou_dot_u.
 */
class constraint_treatment
    : public operator_base
    , public timestepper::lowstorage::linear_operator<
          multi_array::ref<complex_t,4>,
          contiguous_state<4,complex_t>
      >
{
public:

    /**
     * Constructor.
     * After construction, #L should be provided.
     */
    constraint_treatment(
            const scenario_definition& scenario,
            const grid_specification& grid,
            const pencil_grid& dgrid,
            const bsplineop& cop,
            bspline& b,
            operator_common_block& common);

    /** Will bulk total energy be constrained? */
    bool constrain_bulk_rho_E() const;

    /** Will bulk streamwise momentum be constrained? */
    bool constrain_bulk_rho_u() const;

    /** Will bulk density be constrained enforced? */
    bool constrain_bulk_rho() const;

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
            const timestepper::lowstorage::method_interface<complex_t>& method,
            const real_t delta_t,
            const std::size_t substep_index,
            multi_array::ref<complex_t,4> *ic0 = NULL) const;

    /** The operator whose behavior is modified by this instance. */
    shared_ptr<timestepper::lowstorage::linear_operator<
                multi_array::ref<complex_t,4>,
                contiguous_state<4,complex_t>
            > > L;

protected:

    /** The scenario in which the operator is used */
    const scenario_definition &scenario;

    /** Houses data additionally required for some linear operators */
    operator_common_block &common;

    /** Precomputed integration coefficients */
    VectorXr bulkcoeff;

    /**
     * Constraint data passed to #L.
     *
     * \li cdata.col(0) is for the bulk density constraint.
     * \li cdata.col(1) is for the bulk momentum constraint.
     * \li cdata.col(2) is for the bulk total energy constraint.
     *
     * Mutable member avoids repeated allocation/deallocation.
     */
    mutable MatrixX3c cdata;

    /**
     * Least squares constraint solver.
     * Mutable member avoids repeated allocation/deallocation.
     */
    mutable Eigen::JacobiSVD<
            MatrixXXr, Eigen::FullPivHouseholderQRPreconditioner
        > jacobiSvd;
};

} // namespace perfect

} // namespace suzerain

#endif /* SUZERAIN_PERFECT_CONSTRAINT_TREATMENT_HPP */
