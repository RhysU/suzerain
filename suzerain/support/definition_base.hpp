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

#ifndef SUZERAIN_SUPPORT_DEFINITION_BASE_HPP
#define SUZERAIN_SUPPORT_DEFINITION_BASE_HPP

/** @file
 * Provides an abstract base class for handling problem definitions.
 */

#include <suzerain/common.hpp>

// https://svn.boost.org/trac/boost/ticket/7568
SUZERAIN_GCC_DIAG_OFF(unused-parameter);
#include <boost/program_options.hpp>
SUZERAIN_GCC_DIAG_ON(unused-parameter);

namespace suzerain {

namespace support {

/** An abstract base class for defining a collection of related options. */
class definition_base
{
public:

    /**
     * Build a Boost \c options_description encompassing all information in
     * this definition instance.
     *
     * Values copied by Boost Program Options, for example default values, are
     * copied at the time of this invocation.  This method may be expensive to
     * invoke and users are encouraged to cache its return value rather than
     * call it repeatedly.
     *
     * @return An instance suitable for <tt>add</tt>-ing to a
     *         <tt>boost::program_options::options_description</tt> instance.
     *
     * @see <a href="http://www.boost.org/doc/libs/release/libs/program_options">
     *      Boost.Program_options</a> for more information.
     */
    virtual boost::program_options::options_description options_description()
        = 0;

    /** Virtual destructor as appropriate for a base class. */
    virtual ~definition_base() {}

};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_DEFINITION_BASE_HPP
