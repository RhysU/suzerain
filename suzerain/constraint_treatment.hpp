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
#include <suzerain/lowstorage.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/state_fwd.hpp>

namespace suzerain {

// Forward declarations
class bspline;
class pencil_grid;
namespace constraint { class base; }

namespace perfect {

// Forward declarations
class operator_common_block;

// TODO Use Eigen Maps to decouple operator_common_block if sensible

/**
 * A wrapper applying integral constraint treatment atop any linear operator.
 * During \ref invert_mass_plus_scaled_operator implicit momentum forcing is
 * applied following the section of <tt>writeups/channel_treatment.tex</tt>
 * titled "Enforcing a target bulk momentum via the linear operator" and using
 * information from \ref operator_common_block via an instance provided at
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
    : public lowstorage::linear_operator<
          multi_array::ref<complex_t,4>,
          contiguous_state<4,complex_t>
      >
    , public  boost::noncopyable
    , private array<shared_ptr<constraint::base>, 5>
{

    /** Hide the composition-related superclass from client code */
    typedef array<shared_ptr<constraint::base>, 5> implementation_defined;

public:

    /**
     * Constructor.  After construction, #L must be provided.  One or more
     * <tt>operator[](ndx::type)</tt> calls should be made to establish
     * constraints when constraints are desired.
     *
     * @param Ma     Nondimensional Mach number for kinetic energy computations.
     *               In a dimensional setting, this should be one.  Notice that
     *               it is a \e reference permitting the value to track another
     *               setting.  For example, one in a \ref scenario_definition.
     * @param dgrid  Parallel decomposition details used to determine which
     *               rank houses the "zero-zero" Fourier modes.
     * @param common Storage from which mean velocity profiles will be read and
     *               to which implicit forcing averages will be accumulated.
     */
    constraint_treatment(
            const real_t& Ma,
            const pencil_grid& dgrid,
            operator_common_block& common);

    // Permit subscripting/iterating to access equation-specific constraints
    using implementation_defined::operator[];
    using implementation_defined::size;

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

protected:

    /**
     * What Mach number should be used to scale kinetic
     * energy contributions to the total energy equation?
     */
    const real_t& Ma;

    /**
     * Used to obtain mean primitive state profiles to compute
     * implicit forcing work contributions to the total energy equation.
     */
    operator_common_block& common;

private:

    /** Does this rank possess the "zero-zero" Fourier modes? */
    const bool rank_has_zero_zero_modes;

    /**
     * Constraint data passed to #L indexed by ndx::type.
     * Mutable member avoids repeated allocation/deallocation.
     */
    mutable MatrixX5c cdata;

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
