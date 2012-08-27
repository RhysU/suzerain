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
// traits.hpp: utilities for manipulating various types
// $Id$

#ifndef SUZERAIN_TRAITS_HPP
#define SUZERAIN_TRAITS_HPP

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
 * @see ::suzerain::complex::traits::is_complex for details.
 */
template<class T>
struct component<T, typename boost::enable_if<
    ::suzerain::complex::traits::is_complex<T> >::type
> {
    typedef typename ::suzerain::complex::traits::real<T>::type type;
};

} // namespace traits

} // namespace suzerain

#endif // SUZERAIN_TRAITS_HPP
