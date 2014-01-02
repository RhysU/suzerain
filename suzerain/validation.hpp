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

#ifndef SUZERAIN_VALIDATION_HPP
#define SUZERAIN_VALIDATION_HPP

/** @file
 * Helpers to ensure values meet simple validation criteria.
 */

#include <suzerain/common.hpp>

namespace suzerain {

/** Helpers to ensure values meet simple validation criteria. */
namespace validation {

/**
 * Throws an \c invalid_argument exception if a value is neither positive
 * nor NaN.
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
    if (!(t > T(0) || (boost::math::isnan)(t))) {
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
 * Throws an \c invalid_argument exception if a value is neither nonnegative
 * nor NaN.
 *
 * @param t    Value to check.
 * @param name Optional name to use in the exception message
 */
template< typename T >
typename boost::disable_if<boost::is_unsigned<T> >::type
    ensure_nonnegative(T t, const char * name = NULL)
throw(std::invalid_argument)
{
    // NOTE: t has non-const, non-reference type to avoid std::bind2nd issues

    BOOST_STATIC_ASSERT(boost::is_arithmetic<T>::value);
    if (!(t >= T(0) || (boost::math::isnan)(t))) {
        std::ostringstream msg;
        if (name) {
            msg << "Value of '" << name << "' (" << t << ") must be nonnegative";
        } else {
            msg << "Value " << t << " must be nonnegative";
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
typename boost::enable_if<boost::is_unsigned<T> >::type
ensure_nonnegative(T t, const char * name = NULL)
throw(std::invalid_argument)
{
    // NOTE: t has non-const, non-reference type to avoid std::bind2nd issues

    BOOST_STATIC_ASSERT(boost::is_arithmetic<T>::value);

    // NOP: Unsigned types are trivially nonnegative
    SUZERAIN_UNUSED(t);
    SUZERAIN_UNUSED(name);
}

/**
 * Throws an \c invalid_argument exception if a value is outside bounds
 * specified and non NaN.
 *
 * @param t    Value to check.
 * @param name Optional name to use in the exception message
 */
template< typename T >
void ensure_bounded(T t, T lower, T upper, bool lower_inclusive = true,
    bool upper_inclusive = true, const char * name = NULL)
throw(std::invalid_argument)
{
    // NOTE: t has non-const, non-reference type to avoid std::bind2nd issues

    BOOST_STATIC_ASSERT(boost::is_arithmetic<T>::value);
    if (   (   (lower_inclusive ? t >= lower : t > lower)
            && (upper_inclusive ? t <= upper : t < upper))
        || (boost::math::isnan)(t))
      return;

    // Otherwise...
    std::ostringstream msg;
    if (name) {
      msg << "Value of " << name << " (" << t << ") not inside "
          << (lower_inclusive ? '[' : '(') << lower << ", " << upper
          << (upper_inclusive ? ']' : ')');
    } else {
      msg << "Value " << t << " not inside "
          << (lower_inclusive ? '[' : '(') << lower << ", " << upper
          << (upper_inclusive ? ']' : ')');
    }
    throw std::invalid_argument(msg.str());
}

} // namespace validation

} // namespace suzerain

#endif // SUZERAIN_VALIDATION_HPP
