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

#ifndef SUZERAIN_SUPPORT_DEFINITION_ARSEL_HPP
#define SUZERAIN_SUPPORT_DEFINITION_ARSEL_HPP

/** @file
 * Provides \ref definition_arsel.
 */

#include <suzerain/specification_arsel.hpp>
#include <suzerain/support/definition_base.hpp>
#include <suzerain/support/loadable.hpp>
#include <suzerain/support/savable.hpp>

namespace suzerain {

namespace support {

/**
 * Upgrades a \ref specification_arsel with \ref definition_base behavior.
 * This permits using the instance with \ref program_options.
 */
class definition_arsel
    : public specification_arsel
    , public virtual definition_base
////, public virtual loadable
    , public virtual savable
{
public:

    /** Construct an instance with the recommended default values */
    definition_arsel();

    /** @copydoc definition_base::options_description() */
    virtual boost::program_options::options_description options_description();

    /** @copydoc savable::save */
    virtual void save(
            const esio_handle h) const;

/////** @copydoc loadable::load */
////virtual void load(
////        const esio_handle h,
////        const bool verbose = true);

};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_DEFINITION_ARSEL_HPP
