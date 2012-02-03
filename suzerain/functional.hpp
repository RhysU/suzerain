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
// functional.hpp: functional tools
// $Id$

#ifndef __SUZERAIN_FUNCTIONAL_HPP
#define __SUZERAIN_FUNCTIONAL_HPP

#include <suzerain/common.hpp>
#include <suzerain/utility.hpp>
// Pull in libstdc++'s STL algorithm extensions.
// We can reimplement the necessary if this becomes a portability issue.
#include <ext/algorithm>

namespace suzerain {

/**
 * Provides general functional programming tools.
 */
namespace functional {

/**
 * Import SGI STL's copy_n implementation.
 * @see The SGI STL <a href="http://www.sgi.com/tech/stl/copy_n.html">copy_n</a>
 *      or GNU
 *      <a href="http://gcc.gnu.org/onlinedocs/libstdc++/manual/ext_algorithms.html">
 *      libstdc++'s algorithm extensions</a> for more information.
 */
using __gnu_cxx::copy_n;

/**
 * Compute the product of the values in the range <tt>[first,last)</tt>. This
 * combines <tt>std::accumulate</tt> with <tt>std::multiplies<typename
 * std::iterator_traits<InputIterator>::value_type ></tt>.
 *
 * @param first Initial position in a sequence.
 * @param last Final position in a sequence.
 *
 * @return The result of taking the product of all values
 *         in the range <tt>[first,last)</tt>.
 */
template<class InputIterator>
typename std::iterator_traits<InputIterator>::value_type
product(InputIterator first, InputIterator last) {
    typedef typename
        std::iterator_traits<InputIterator>::value_type value_type;
    return std::accumulate(
            first, last, value_type(1), std::multiplies<value_type>());
}

/**
 * Compute the product of the values in the range <tt>[first,last)</tt> with \c
 * init.  This combines <tt>std::accumulate</tt> with
 * <tt>std::multiplies<T></tt>.  Note that \c init controls the multiplication
 * rules in use.
 *
 * @param first Initial position in a sequence.
 * @param last Final position in a sequence.
 * @param init Initial value for the product.
 *
 * @return The result of taking the product of \c init with all values
 *         in the range <tt>[first,last)</tt>.
 */
template<class InputIterator, class T>
T product(InputIterator first, InputIterator last, T init) {
    return std::accumulate(first, last, init, std::multiplies<T>());
}

/**
 * Compute the product of values in the range <tt>[c.begin(),c.end())</tt>.
 *
 * @param c Collection whose values should be computed.
 *
 * @return The result of taking the product of all values in the range
 *         <tt>[c.begin(),c.end())</tt>.
 * @see The Boost <a
 * href="http://www.boost.org/doc/libs/release/libs/utility/Collection.html">
 * Collection concept documentation</a> for more information.
 */
template<class Collection>
typename Collection::value_type product(const Collection& c) {
    return product(c.begin(), c.end());
}

/** A binary operator that computes the minimum of two values. */
template<class T>
struct minimum : std::binary_function <T,T,T> {
    T operator() (const T& x, const T& y) const {
        return std::min(x,y);
    }
};

/** A binary operator that computes the maximum of two values. */
template<class T>
struct maximum : std::binary_function <T,T,T> {
    T operator() (const T& x, const T& y) const {
        return std::max(x,y);
    }
};

/** An functor computing suzerain::any */
struct any {
    typedef bool result_type;

    template<typename InputIterator>
    bool operator()(InputIterator first, InputIterator last) const
    {
        return suzerain::any(first, last);
    }
};

/** An functor computing suzerain::all */
struct all {
    typedef bool result_type;

    template<typename InputIterator>
    bool operator()(InputIterator first, InputIterator last) const
    {
        return suzerain::all(first, last);
    }
};

/**
 * A functor that performs assignment to type \c Target from type \c Source.
 * \c Target must be assignable from \c Source.
 *
 * The \c Enable template parameter allows using <tt>boost::enable_if</tt> to
 * specialize or partially specialize the functor per <tt>enable_if</tt>'s
 * documentation section <a
 * href="http://www.boost.org/doc/libs/1_41_0/libs/utility/enable_if.html">
 * 3.1: Enabling template class specializations</a>.
 */
template<class Target, class Source, class Enable = void>
struct assign {

    /**
     * Create an instance which assigns \c s when applied.
     *
     * @param s source of assignment operation occurring via
     *          <tt>operator()</tt>.
     */
    assign(const Source &s) : s_(s) {}

    /**
     * Assign the value provided at construction to \c t.
     *
     * @param t to be assigned.
     */
    void operator()(Target& t) const { t = s_; }

private:
    const Source &s_; /**< Source for assignment operations */
};

} // namespace functional

} // namespace suzerain

#endif // __SUZERAIN_FUNCTIONAL_HPP
