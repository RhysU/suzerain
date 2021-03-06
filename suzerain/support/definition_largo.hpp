//--------------------------------------------------------------------------
//
// Copyright (C) 2012-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
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

#ifndef SUZERAIN_SUPPORT_DEFINITION_LARGO_HPP
#define SUZERAIN_SUPPORT_DEFINITION_LARGO_HPP

/** @file
 * Provides \ref definition_largo.
 */

#include <suzerain/common.hpp>
#include <suzerain/specification_largo.hpp>
#include <suzerain/support/definition_base.hpp>
#include <suzerain/support/esio_fwd.hpp>
#include <suzerain/support/loadable.hpp>
#include <suzerain/support/overridable.hpp>
#include <suzerain/support/populatable.hpp>
#include <suzerain/support/savable.hpp>

namespace suzerain
{

namespace support
{

/**
 * Holds parameters defining largo boundary cases.
 */
class definition_largo
    : public virtual definition_base
    , public virtual loadable
    , public virtual overridable<definition_largo> // Not specification_largo!
    , public virtual populatable<definition_largo> // Not specification_largo!
    , public virtual savable
    , public specification_largo
{
public:

    /**
     * Construct an instance with all parameters set to NaN.
     * Clients can use NaN as a not-yet-specified or use-the-default value.
     */
    definition_largo();

    /** @copydoc populatable::populate */
    virtual void populate(
        const definition_largo& that,
        const bool verbose = false);

    /** @copydoc overridable::override */
    virtual void override(
        const definition_largo& that,
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

#endif // SUZERAIN_SUPPORT_DEFINITION_LARGO_HPP
