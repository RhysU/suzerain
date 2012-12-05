//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// definition_base.hpp: base class for handling problem definitions
// $Id$

#ifndef SUZERAIN_SUPPORT_DEFINITION_BASE_HPP
#define SUZERAIN_SUPPORT_DEFINITION_BASE_HPP

#include <boost/program_options.hpp>

#include <suzerain/common.hpp>

/** @file
 * Provides \ref definition_base
 */

namespace suzerain {

namespace support {

/**
 * An abstract base class for problem-related definitions.  These are
 * collections of runtime arguments or scenario parameters used to perform
 * calculations.
 */
class definition_base
{
public:
    /**
     * Constructor setting no overall options description.
     */
    definition_base() : options_() {};

    /**
     * Constructor setting an overall options description.
     *
     * @param caption Caption to use when describing options
     **/
    definition_base(const std::string &caption) : options_(caption) {};

    /** Virtual destructor */
    virtual ~definition_base() {};

    /**
     * Obtain a Boost options_description encompassing all information
     * in the definition.
     *
     * @return A reference suitable for <tt>add</tt>-ing to a
     *         <tt>boost::program_options::options_description</tt> instance.
     *
     * @see <a href="http://www.boost.org/doc/libs/release/libs/program_options">
     *      Boost.Program_options</a> for more information.
     */
    const boost::program_options::options_description& options() const
    {
        return options_;
    }

    /**
     * Obtain a Boost options_description encompassing all information
     * in the definition.
     *
     * @return A reference suitable for <tt>add</tt>-ing to a
     *         <tt>boost::program_options::options_description</tt> instance.
     *
     * @see <a href="http://www.boost.org/doc/libs/release/libs/program_options">
     *      Boost.Program_options</a> for more information.
     */
    boost::program_options::options_description& options()
    {
        return options_;
    }

    /**
     * Return a Boost program_options_easy_init suitable for adding options
     * to this definition.
     *
     * @return An
     *         <tt>boost::program_options::options_description_easy_init</tt>
     *         instance.
     *
     * @see <a href="http://www.boost.org/doc/libs/release/libs/program_options">
     *      Boost.Program_options</a> for more information.
     */
    boost::program_options::options_description_easy_init add_options()
    {
        return options_.add_options();
    }

protected:

    /** Stores the program options processing information */
    boost::program_options::options_description options_;
};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_DEFINITION_BASE_HPP
