//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
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

#ifndef SUZERAIN_SUPPORT_LOADABLE_HPP
#define SUZERAIN_SUPPORT_LOADABLE_HPP

/** @file
 * Provides \ref loadable.
 */

#include <suzerain/support/esio_fwd.hpp>

namespace suzerain {

namespace support {

/** Abstract interface indicating details may be loaded from an ESIO handle. */
class loadable
{
public:

    /**
     * Load details from an ESIO-based file.
     *
     * Descendants should override this method adding any desired functionality
     * either before or after invoking the superclass version.
     *
     * @param h       Open, readable handle from which details will be loaded.
     * @param verbose Should logging be emitted when a value is retained?
     */
    virtual void load(
            const esio_handle h,
            const bool verbose = true) = 0;

    /** Virtual destructor to permit deleting through subclass. */
    virtual ~loadable() {}

};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_LOADABLE_HPP
