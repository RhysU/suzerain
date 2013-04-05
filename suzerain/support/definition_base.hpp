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

protected:

    /**@name Helpers for subclass implementations */
    /**@{*/

    /**
     * If \c destination is NaN, populate it with the value from \c source.
     * When \c verbose, log an informative message using \c name and \c
     * description.
     *
     * @param name        Name to use for any logging.
     * @param description A short description of \c name to use for logging.
     *                    If \c NULL, no description will be logged.
     * @param destination Location to possibly assign from \c source.
     * @param source      Source data to possibly use.
     * @param verbose     Should a human-readable message be logged?
     *
     * @return \c True whenever \c destination was assigned from \c source.
     */
    static bool maybe_populate(const char*   name,
                               const char*   description,
                                     real_t& destination,
                               const real_t& source,
                               const bool    verbose);

    /**
     * If \c source is non-NaN, override \c destination with its value.
     * @copydetails maybe_populate(const char*,const char*,real_t&,const real_t&,const bool)
     */
    static bool maybe_override(const char*   name,
                               const char*   description,
                                     real_t& destination,
                               const real_t& source,
                               const bool    verbose);

    /**
     * If \c destination is zero, populate it with the value from \c source.
     * @copydetails maybe_populate(const char*,const char*,real_t&,const real_t&,const bool)
     */
    static bool maybe_populate(const char* name,
                               const char* description,
                                     int&  destination,
                               const int&  source,
                               const bool  verbose);

    /**
     * If \c source is non-zero, override \c destination with its value.
     * @copydetails maybe_populate(const char*,const char*,real_t&,const real_t&,const bool)
     */
    static bool maybe_override(const char* name,
                               const char* description,
                                     int&  destination,
                               const int&  source,
                               const bool  verbose);

    /**@}*/
};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_DEFINITION_BASE_HPP
