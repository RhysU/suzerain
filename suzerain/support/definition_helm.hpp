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

namespace suzerain {

namespace support {

/**
 * Upgrades a \ref specification_helm with \ref definition_base behavior.
 * This permits using the instance with \ref program_options.
 */
class definition_helm
    : public specification_helm
    , public virtual definition_base
{
public:

    /** \copydoc specification_helm() */
    explicit
    definition_helm(const double kp,
                    const double Td = 1.00,
                    const double Tf = 0.01,
                    const double Ti = 1.00,
                    const double Tt = 1.00);

    /** @copydoc definition_base::options_description() */
    virtual boost::program_options::options_description options_description();

};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_DEFINITION_HELM_HPP
