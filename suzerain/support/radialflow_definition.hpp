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

#ifndef SUZERAIN_SUPPORT_RADIALFLOW_DEFINITION_HPP
#define SUZERAIN_SUPPORT_RADIALFLOW_DEFINITION_HPP

/** @file
 * Routines for adding noise/perturbations to state fields
 */

#include <suzerain/common.hpp>
#include <suzerain/radialflow_specification.hpp>
#include <suzerain/support/definition_base.hpp>
#include <suzerain/support/esio_fwd.hpp>
#include <suzerain/support/loadable.hpp>
#include <suzerain/support/overridable.hpp>
#include <suzerain/support/populatable.hpp>
#include <suzerain/support/savable.hpp>

namespace suzerain {

namespace support {

/**
 * Upgrades a \ref radialflow_specification with \ref definition_base behavior.
 * This permits using the instance with \ref program_options.
 */
class radialflow_definition
    : public virtual definition_base
    , public virtual loadable
    , public virtual overridable<radialflow_specification>
    , public virtual populatable<radialflow_specification>
    , public virtual savable
    , public radialflow_specification
{
public:

    /** Construct an instance with the given default values */
    explicit radialflow_definition(
            double deltae = std::numeric_limits<double>::quiet_NaN(),
            double gam0   = std::numeric_limits<double>::quiet_NaN(),
            double Mae    = std::numeric_limits<double>::quiet_NaN(),
            double pexi   = std::numeric_limits<double>::quiet_NaN(),
            double Te     = std::numeric_limits<double>::quiet_NaN());

    /** @copydoc populatable::populate */
    virtual void populate(
        const radialflow_specification& that,
        const bool verbose = false);

    /** @copydoc overridable::override */
    virtual void override(
        const radialflow_specification& that,
        const bool verbose = false);

    /** @copydoc savable::save */
    virtual void save(
        const esio_handle h) const;

    /** @copydoc loadable::load */
    virtual void load(
        const esio_handle h,
        const bool verbose = true);

    /** @copydoc definition_base::options_description() */
    virtual boost::program_options::options_description options_description();

};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_RADIALFLOW_DEFINITION_HPP
