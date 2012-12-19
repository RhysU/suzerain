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

} // namespace validation

} // namespace suzerain

#endif // SUZERAIN_VALIDATION_HPP
