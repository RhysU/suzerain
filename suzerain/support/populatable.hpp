//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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

#ifndef SUZERAIN_SUPPORT_POPULATABLE_HPP
#define SUZERAIN_SUPPORT_POPULATABLE_HPP

/** @file
 * Provides \ref populatable.
 */

namespace suzerain {

namespace support {

/**
 * Abstract interface indicating details may be populated from another instance.
 */
template <class Derived>
class populatable
{
public:

    /**
     * Populate any NaN members in \c this with values from \c that.
     *
     * Descendants should override this method adding any desired functionality
     * either before or after invoking the superclass version.
     *
     * @param that    Instance from which information is taken.
     * @param verbose Should logging be emitted when a value is retained?
     */
    virtual void populate(
            const Derived& that,
            const bool verbose = false) = 0;

};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_POPULATABLE_HPP
