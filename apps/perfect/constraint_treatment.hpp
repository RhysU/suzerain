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

/**
 * Encapsulates the ways \ref constraint_treatment can constrain a value.
 */
class constraint
{
public:

    /** A particular value can be constrained in what ways? */
    enum what_type {
          nothing = 0 ///< Enforce nothing.  That is, no constraint.
        , value_lower ///< Enforce collocation value at \f$y=0\f$
        , value_upper ///< Enforce collocation value at \f$y=L_y\f$
        , value_bulk  ///< Enforce bulk value across \f$y=\left[0,L_y\right]\f$
    } what;

    /** What numeric value should be targeted? */
    real_t target;

    /** Default constructor creates an unenforced constraint. */
    constraint()
        : what(nothing)
        , target(std::numeric_limits<real_t>::quiet_NaN())
    {}

    /** Construct an instance enforcing \c target in \c what manner. */
    constraint(const what_type what, const real_t target)
        : what(what)
        , target(target)
    {}

    /** Is this constraint enforceable as specified? */
    bool enabled() const
    { return what != nothing && !(boost::math::isnan)(target); }

};

/**
 * A wrapper applying integral constraint treatment atop any linear operator.
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
     * Constructor.  After construction, #L must be provided.
     * One or more of #rho, #rho_u, or #rho_E should be called when
     * constraints are desired.
     *
     * @param Ma Nondimensional Mach number for kinetic energy computations.
     *           In a dimensional setting, this should be one.  Notice that it
     *           is a \e reference permitting the value to track other
     *           settings.  For example, those in a \ref scenario_definition.
     */
    constraint_treatment(
            real_t& Ma,
            const grid_specification& grid,
            const pencil_grid& dgrid,
            const bsplineop& cop,
            bspline& b,
            operator_common_block& common);

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

    /** Specify the desired constraint treatment for density \f$\rho\f$. */
    constraint_treatment& specify_rho(const constraint& c);

    /** Specify the desired constraint treatment for momentum \f$\rho{}u\f$. */
    constraint_treatment& specify_rho_u(const constraint& c);

    /** Specify the desired constraint treatment for energy \f$\rho{}E\f$. */
    constraint_treatment& specify_rho_E(const constraint& c);

protected:

    /** Helper used to implement \c specify_XXX methods */
    constraint_treatment& specify(const constraint& src,
                                        constraint& dst,
                                        VectorXr&   coeff);

    /**
     * What Mach number should be used to scale kinetic
     * energy contributions to the total energy equation?
     */
    real_t &Ma;

    /** Houses data additionally required for some linear operators */
    operator_common_block &common;

    /** In what manner should \f$\rho\f$ be constrained? */
    constraint rho;

    /** In what manner should \f$\rho{}u\f$ be constrained? */
    constraint mx;

    /** In what manner should \f$\rho{}E\f$ be constrained? */
    constraint e;

    /** Precomputed integration coefficients */
    VectorXr coeff_bulk;

    /** Precomputed integration-like coefficients for \c rho constraint */
    VectorXr coeff_rho;

    /** Precomputed integration-like coefficients for \c rho_u constraint */
    VectorXr coeff_mx;

    /** Precomputed integration-like coefficients for \c rho_E constraint */
    VectorXr coeff_e;

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
