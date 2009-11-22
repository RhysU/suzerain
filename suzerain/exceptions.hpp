/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
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
 * exceptions.hpp: Exceptions used within Suzerain
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_EXCEPTIONS_H
#define __SUZERAIN_EXCEPTIONS_H

#include <suzerain/common.hpp>

/** \file
 * Provides exception types used throughout Suzerain.
 */

namespace suzerain
{

/** Reports arguments to functions that are outside the valid input range.
 *
 * \internal Intended to have the same semantics as \c std::domain_error
 * but with the benefit of boost::exception as a base class.
 */
class domain_error: public std::domain_error, public boost::exception {
public:

    /*!
     * Construct an instance reporting \c what_arg through its what member
     * function.
     *
     * \param what_arg details describing the cause of the exception.
     */
    domain_error(const std::string &what_arg)
        : std::domain_error(what_arg) {};
};

/** Reports an invalid argument was provided to a function.
 *
 * \internal Intended to have the same semantics as \c std::invalid_argument
 * but with the benefit of boost::exception as a base class.
 */
class invalid_argument: public std::invalid_argument, public boost::exception {
public:

    /*! \copydoc suzerain::domain_error::domain_error(const std::string&) */
    invalid_argument(const std::string &what_arg)
        : std::invalid_argument(what_arg) {};
};

/** Reports function arguments indicative of a programming error.
 *
 * \internal Intended to have the same semantics as \c std::logic_error
 * but with the benefit of boost::exception as a base class.
 */
class logic_error: public std::logic_error, public boost::exception {
public:

    /*! \copydoc domain_error::domain_error(const std::string&) */
    logic_error(const std::string &what_arg)
        : std::logic_error(what_arg) {};
};

} // namespace suzerain

#endif // __SUZERAIN_EXCEPTIONS_H
