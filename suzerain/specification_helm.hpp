//--------------------------------------------------------------------------
//
// Copyright (C) 2014 Rhys Ulerich
// Copyright (C) 2014 The PECOS Development Team
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

#ifndef SUZERAIN_SPECIFICATION_HELM_HPP
#define SUZERAIN_SPECIFICATION_HELM_HPP

/** @file
 * Provides \ref specification_helm.
 */

#include <suzerain/helm.h>

namespace suzerain {

/**
 * Encapsulates PID controller settings and operation using \ref helm.h.
 * Hides aspects of the incremental control from \ref helm.h to ease use.
 * Const-ness in this class refers to logical const-ness as
 * operating the control mutates internal state.  The controller
 * tracks whether or not it has been #enabled.
 */
class specification_helm
{

public:

    /**
     * Construct an instance with standard-form control settings.
     * That is, with coefficients from the equation
     * \f{align}{
     *     \frac{\mathrm{d}}{\mathrm{d}t} v(t) &= k_p \left[
     *               - \frac{\mathrm{d}}{\mathrm{d}t} y(t)
     *               + \frac{r(t) - y(t)}{T_i}
     *               + \frac{u(t) - v(t)}{T_t}
     *               + \frac{T_d}{T_f}\left(
     *                   \frac{\mathrm{d}}{\mathrm{d}t} f(t)
     *                 - \frac{\mathrm{d}}{\mathrm{d}t} y(t)
     *                 \right)
     *             \right]
     *             .
     * \f}
     * To be enabled, argument \c kp must be nonzero.
     * After construction, \ref approach must be called once before \ref steady.
     *
     * \param kp \copydoc helm_state#kp
     * \param Td \copydoc helm_state#Td
     * \param Tf \copydoc helm_state#Tf
     * \param Ti \copydoc helm_state#Ti
     * \param Tt \copydoc helm_state#Tp
     */
    explicit
    specification_helm(const double kp = 1.00,
                       const double Td = 1.00,
                       const double Tf = 0.01,
                       const double Ti = 1.00,
                       const double Tt = 1.00);

    /** Is the controller turned on? */
    bool enabled;

    /**
     * \brief Reset all tuning parameters, but \e not transient state.
     * \see \ref helm_reset
     */
    specification_helm *
    reset();

    /**
     * \brief Reset any transient state, but \e not tuning parameters.
     * Has the side-effect of setting <code>enabled == true</code>.
     *
     * \param t Absolute time at which control is being established.
     * \param v Absolute actuator position to establish.
     * \param r Desired reference value, often called the setpoint.
     * \see \ref helm_approach
     */
    specification_helm *
    approach(const double t,
             const double v,
             const double r);

    /**
     * \brief Find the control signal necessary to steady unsteady process y(t).
     * \param t Current absolute simulation time.
     * \param u Currently observed actuator signal.
     * \param y Currently observed process output.
     * \return Incremental suggested change to control signal \c u.
     * \see \ref helm_steady
     */
    double
    steady(const double t,
           const double u,
           const double y) const;

protected:

    /** The state and settings for underlying incremental controller. */
    mutable helm_state h;

    /** The last simulation time \c t at which steady was called. */
    mutable double t;

    /** Tracks the position-based control signal \c v. */
    mutable double v;

    /** Desired reference signal, often called the setpoint. */
    mutable double r;

};

} // namespace suzerain

#endif // SUZERAIN_SPECIFICATION_HELM_HPP
