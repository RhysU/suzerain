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

#ifndef SUZERAIN_SUPPORT_DEFINITION_HELM_HPP
#define SUZERAIN_SUPPORT_DEFINITION_HELM_HPP

/** @file
 * Provides \ref definition_helm.
 */

#include <suzerain/specification_helm.hpp>
#include <suzerain/support/definition_base.hpp>
#include <suzerain/support/savable.hpp>

namespace suzerain {

namespace support {

/**
 * Upgrades a \ref specification_helm with \ref definition_base behavior.
 * This permits using the instance with \ref program_options.
 */
class definition_helm
    : public specification_helm
    , public virtual definition_base
    , public virtual savable
{
public:

    /**
     * \copydoc specification_helm()
     * \param r Reference value, often called the setpoint.
     */
    explicit
    definition_helm(const double kp,
                    const double Td = 1.00,
                    const double Tf = 0.01,
                    const double Ti = 1.00,
                    const double Tt = 1.00,
                    const double r  = 1.00);

    /** @copydoc definition_base::options_description() */
    virtual boost::program_options::options_description options_description();

    // Make specification_helm::approach(t, v, r) visible
    using specification_helm::approach;

    /**
     * \brief Reset any transient state, but \e not tuning parameters.
     * Reference from construction, possibly modified by options parsing, used.
     * \param t Absolute time at which control is being established.
     * \param v Absolute actuator position to establish.
     * \see \ref specification_helm::approach.
     */
    definition_helm *
    approach(const double t,
             const double v)
    {
        specification_helm::approach(t, v, this->r);
        return this;
    }

    /** @copydoc savable::save */
    virtual void save(
            const esio_handle h) const;

};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_DEFINITION_HELM_HPP
