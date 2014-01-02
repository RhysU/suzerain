//--------------------------------------------------------------------------
//
// Copyright (C) 2013-2014 Rhys Ulerich
// Copyright (C) 2013-2014 The PECOS Development Team
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

#ifndef SUZERAIN_HYBRID_RESIDUAL_OPERATOR_HPP
#define SUZERAIN_HYBRID_RESIDUAL_OPERATOR_HPP

/** @file
 * Provides \ref hybrid_residual_operator.
 */

#include <suzerain/common.hpp>
#include <suzerain/lowstorage.hpp>
#include <suzerain/state_fwd.hpp>

namespace suzerain {

/**
 * Builds a hybrid implicit/explicit-ready \ref nonlinear_operator from an
 * implicit-agnostic \ref nonlinear_operator and a compatible \ref
 * linear_operator.
 *
 * More concretely, suppose one wants to use the hybrid implicit/explicit
 * timestepping schemes defined in \ref lowstorage which the state
 * vector \f$ u(t) \f$ is advanced to \f$u(t+\Delta{}t)\f$ according to
 * \f[
 *   M u_{t} = Lu + \chi N(u)
 * \f]
 * where \f$M\f$, \f$L\f$, and \f$N\f$ are linear, linear, and nonlinear
 * operators, respectively.  Scalar \f$\chi\f$ is a fixed multiplicative factor
 * often used to hide FFT normalization costs.  Instead of \f$N\f$ you have a
 * fully explicit operator implementation $R$ per
 * \f[
 *   M u_{t} = \chi R(u).
 * \f]
 * and an implementation of \f$M + \phi L\f$ per \ref
 * lowstorage::linear_operator.  To use the low-storage
 * implicit-explicit interface you need
 * \f[
 *   N(u) = R(u) - \frac{1}{\chi} Lu.
 * \f].
 *
 * This class programmatically implements such an \f$N\f$ given implementations
 * of \f$R\f$ and \f$M + \phi L\f$.  It permits having your \f$N\f$
 * automatically adjust to any sane \f$L\f$ implementation but it comes at the
 * cost of nontrivial additional memory and moderate runtime overhead.
 *
 * When a good \f$L\f$ and \f$R\f$ combination is found, it is well worthwhile
 * to encode the corresponding \f$N\f$ directly rather than to continue using
 * this adaptor.
 */
class hybrid_residual_operator
    : public lowstorage::nonlinear_operator< contiguous_state<4,complex_t> >
{
public:

    /**
     * Default constructor.
     *
     * After construction, #L and #R must be specified
     * prior to #apply_operator invocation.
     *
     * @param chi   Scaling factor \f$\chi\f used to form
     *              \f$ N(u) = R(u) - \frac{1}{\chi} L(u) f$.
     */
    explicit hybrid_residual_operator(const real_t chi);

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
    shared_ptr<lowstorage::linear_operator<
                state_common_type, state_nonlinear_type
            > > L;

    /**
     * The fully-explicit nonlinear operator which operates on
     * #state_nonlinear.
     */
    shared_ptr<lowstorage::nonlinear_operator<
                state_nonlinear_type
            > > R;

    // Inherits documentation
    virtual std::vector<real_t> apply_operator(
            const real_t time,
            state_nonlinear_type &state,
            const lowstorage::method_interface<element>& method,
            const std::size_t substep_index) const;

private:

    /**
     * Scaling factor \f$\chi\f used to form
     * \f$ N(u) = R(u) - \frac{1}{\chi} L(u) f$.
     */
    const real_t chi;

    // boost::noncopyable trips Intel non-virtual base destructor warnings.
    hybrid_residual_operator(const hybrid_residual_operator&);
    hybrid_residual_operator& operator=(const hybrid_residual_operator&);

};

} // namespace suzerain

#endif  /* SUZERAIN_HYBRID_RESIDUAL_OPERATOR_HPP */
