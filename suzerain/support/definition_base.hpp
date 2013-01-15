//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#ifndef SUZERAIN_SUPPORT_DEFINITION_BASE_HPP
#define SUZERAIN_SUPPORT_DEFINITION_BASE_HPP

/** @file
 * Provides an abstract base class for handling problem definitions.
 */

#include <boost/program_options.hpp>

#include <suzerain/common.hpp>

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
