//--------------------------------------------------------------------------
//
// Copyright (C) 2013 Rhys Ulerich
// Copyright (C) 2013 The PECOS Development Team
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

#ifndef SUZERAIN_SUPPORT_HYBRID_RESIDUAL_OPERATOR_HPP
#define SUZERAIN_SUPPORT_HYBRID_RESIDUAL_OPERATOR_HPP

/** @file
 * Provides \ref hybrid_residual_operator.
 */

#include <suzerain/common.hpp>
#include <suzerain/state_fwd.hpp>
#include <suzerain/timestepper.hpp>

namespace suzerain {

namespace support {

/**
 * Builds a hybrid implicit/explicit-ready \ref nonlinear_operator from an
 * implicit-agnostic \ref nonlinear_operator and a compatible \ref
 * linear_operator.
 *
 * More concretely, suppose one wants to use the hybrid implicit/explicit
 * timestepping schemes defined in \ref timestepper::lowstorage which the state
 * vector \f$ u(t) \f$ to \f$u(t+\Delta{}t)\f$ according to
 * \f[
 *   M u_{t} = Lu + N(u)
 * \f]
 * where \f$M\f$, \f$L\f$, and \f$N\f$ are linear, linear, and nonlinear
 * operators, respectively.  Instead of \f$N\F$ you have a fully explicit
 * operator implementation $R$ satisfying
 * \f[
 *   M u_{t} = R(u).
 * \f]
 * and an implementation of \f$M + \phi L\f$ per \ref
 * timestepper::lowstorage::linear_operator.  To use the low-storage
 * implicit-explicit interface you need \f$N(u) = R(u) - L(u)\f$.
 *
 * This class programmatically implements such an \f$N\f$ given implementations
 * of \f$R\f$ and \f$M + \phi L\f$.  It permits having your \f$N\f$
 * automatically adjust to any \f$L\f$ but it comes at the cost of additional
 * memory overhead and runtime cost.
 */
class hybrid_residual_operator
    : public timestepper::nonlinear_operator< contiguous_state<4,complex_t> >
{
public:

    /**
     * Default constructor.  After construction, #L and #R must be specified
     * prior to #apply_operator invocation.
     * */
    hybrid_residual_operator();

    /** State type associated with the linear operator #L. */
    typedef interleaved_state<4, complex_t> state_linear_type;

    /** State type associated with the nonlinear operator #R. */
    typedef contiguous_state<4, complex_t> state_nonlinear_type;

    /**
     * An ancestor common to \ref state_linear_type and \ref
     * state_nonlinear_type.  A common ancestor is necessary for
     * interoperability of #L and #N.
     */
    typedef multi_array::ref<complex_t, 4> state_common_type;

    /**
     * The hybrid implicit/explicit linear operator which is able to
     * interoperate between #state_linear and #state_nonlinear.
     */
    shared_ptr<timestepper::lowstorage::linear_operator<
                state_common_type, state_nonlinear_type
            > > L;

    /**
     * The fully-explicit nonlinear operator which operators on
     * #state_nonlinear.
     */
    shared_ptr<timestepper::nonlinear_operator<
                state_nonlinear_type
            > > R;

    // Inherits documentation
    virtual std::vector<real_t> apply_operator(
            const real_t time,
            state_nonlinear_type &state,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag,
            const std::size_t substep_index) const;

private:

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;

    // Using boost::noncopyable trips Intel non-virtual base destructor warnings.
    hybrid_residual_operator(const hybrid_residual_operator&);
    hybrid_residual_operator& operator=(const hybrid_residual_operator&);

};

} // namespace support

} // namespace suzerain

#endif  /* SUZERAIN_SUPPORT_HYBRID_RESIDUAL_OPERATOR_HPP */
