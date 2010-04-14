/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Suzerain is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * problem.hpp: classes handling problem definitions
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
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
    /** Constructor appropriate for an abstract base class */
    IDefinition() {};

    /** Virtual destructor appropriate for an abstract base class */
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
    virtual const boost::program_options::options_description& options() = 0;
};

} // namespace problem

} // namespace suzerain

#endif // __SUZERAIN_PROBLEM_HPP
