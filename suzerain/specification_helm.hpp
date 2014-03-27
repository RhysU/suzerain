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
 */
class specification_helm : public helm_state
{

public:

    /**
     * Construct an instance with default values.
     * To enable, argument \c kp must be nonzero.
     */
    explicit
    specification_helm(const double kp,
                       const double Td = 1,
                       const double Tf = 0.01,
                       const double Ti = 1,
                       const double Tt = 1);

    /**
     * \brief Is the controller turned on mean is #kp nonzero?
     * \return True if #kp is nonzero.  False otherwise.
     */
    bool
    enabled();

    /**
     * \brief Reset all tuning parameters, but \e not transient state.
     * \see \ref helm_reset
     */
    specification_helm *
    reset();

    /**
     * \brief Reset any transient state, but \e not tuning parameters.
     * \param r Desired reference value, often called the setpoint.
     * \param v Absolute actuator position to establish.
     * \see \ref helm_approach
     */
    specification_helm *
    approach(const double r,
             const double v);

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
           const double y);

    /** Desired reference signal, often called the setpoint. */
    double r;

    /** The last simulation time \c t at which steady was called. */
    double t;

    /** Tracks the position-based control signal \c v. */
    double v;

};

} // namespace suzerain

#endif // SUZERAIN_SPECIFICATION_HELM_HPP
