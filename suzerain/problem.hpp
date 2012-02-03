//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// problem.hpp: classes handling problem definitions
// $Id$

#ifndef __SUZERAIN_PROBLEM_HPP
#define __SUZERAIN_PROBLEM_HPP

#include <suzerain/common.hpp>

/** @file
 * Provides classes handling problem definitions, which are runtime
 * arguments or scenario parameters used to perform calculations.
 */

namespace suzerain {

/**
 * Provides classes handling problem definitions, which are runtime
 * arguments or scenario parameters used to perform calculations.
 */
namespace problem {

/**
 * An interface describing a problem definition.
 */
class IDefinition
{
public:
    /**
     * Constructor setting no overall options description.
     */
    IDefinition() : options_() {};

    /**
     * Constructor setting an overall options description.
     *
     * @param caption Caption to use when describing options
     **/
    IDefinition(const std::string &caption) : options_(caption) {};

    /** Virtual destructor */
    virtual ~IDefinition() {};

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
    const boost::program_options::options_description& options() const {
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
    boost::program_options::options_description& options() {
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
    boost::program_options::options_description_easy_init add_options() {
        return options_.add_options();
    }

protected:

    /** Stores the program options processing information */
    boost::program_options::options_description options_;
};

} // namespace problem

} // namespace suzerain

#endif // __SUZERAIN_PROBLEM_HPP
