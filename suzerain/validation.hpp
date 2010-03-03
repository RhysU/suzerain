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
 * validation.hpp: helpers that ensure values meet some validation criteria
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_VALIDATION_HPP
#define __SUZERAIN_VALIDATION_HPP

#include <suzerain/common.hpp>

/** @file
 * Provides helpers to ensure values meet simple validation criteria.
 */

namespace suzerain {

/** Provides helpers to ensure values meet simple validation criteria. */
namespace validation {

/**
 * Throws an \c invalid_argument exception if a value is not positive.
 *
 * @param t    Value to check.
 * @param name Optional name to use in the exception message
 */
template< typename T >
void ensure_positive(T t, const char * name = NULL)
throw(std::invalid_argument)
{
    // NOTE: t has non-const, non-reference type to avoid std::bind2nd issues

    BOOST_STATIC_ASSERT(boost::is_arithmetic<T>::value);
    if (!(t > T(0))) {
        std::ostringstream msg;
        if (name) {
            msg << "Value of '" << name << "' (" << t << ") must be positive";
        } else {
            msg << "Value " << t << " must be positive";
        }
        throw std::invalid_argument(msg.str());
    }
}

/**
 * Throws an \c invalid_argument exception if a value is not nonnegative.
 *
 * @param t    Value to check.
 * @param name Optional name to use in the exception message
 */
template< typename T >
void ensure_nonnegative(T t, const char * name = NULL)
throw(std::invalid_argument)
{
    // NOTE: t has non-const, non-reference type to avoid std::bind2nd issues

    BOOST_STATIC_ASSERT(boost::is_arithmetic<T>::value);
    const bool applies =    boost::is_floating_point<T>::value
                         || boost::is_signed<T>::value;
    if (applies && !(t >= T(0))) {
        std::ostringstream msg;
        if (name) {
            msg << "Value of '" << name << "' (" << t << ") must be nonnegative";
        } else {
            msg << "Value " << t << " must be nonnegative";
        }
        throw std::invalid_argument(msg.str());
    }
}

} // namespace validation

} // namespace suzerain

#endif // __SUZERAIN_VALIDATION_HPP
