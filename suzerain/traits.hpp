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

#ifndef SUZERAIN_TRAITS_HPP
#define SUZERAIN_TRAITS_HPP

/** @file
 * Traits classes for manipulating various types
 */

#include <suzerain/common.hpp>
#include <suzerain/complex.hpp>

namespace suzerain {

/**
 * Provides type traits suitable for manipulating real- and complex-valued
 * types.
 *
 * @see <a href="http://www.boost.org/doc/libs/release/libs/type_traits/">
 * Boost.TypeTraits</a> for more details.
 */
namespace traits {

/**
 * A type trait that, given a real- or complex-valued type,
 * returns the underlying scalar component type as <tt>type</tt>.
 */
template<class T, class Enable = void> struct component {};

/** A specialization to handle recognized real-valued arithmetic types. */
template<class T>
struct component<T, typename boost::enable_if<
    boost::is_arithmetic<T> >::type
> {
    typedef T type;
};

/**
 * A specialization to handle all recognized complex types.
 * @see suzerain::complex::traits::is_complex for details.
 */
template<class T>
struct component<T, typename boost::enable_if<
    suzerain::complex::traits::is_complex<T> >::type
> {
    typedef typename suzerain::complex::traits::real<T>::type type;
};

} // namespace traits

} // namespace suzerain

#endif // SUZERAIN_TRAITS_HPP
