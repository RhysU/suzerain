/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2009 The PECOS Development Team
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
 * multi_array.hpp: Utilities for Boost.MultiArray
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_MULTIARRAY_H
#define __SUZERAIN_MULTIARRAY_H

#include <suzerain/common.hpp>

namespace suzerain
{

/**
 * Provides utilities atop the Boost.MultiArray concept
 * and its concrete implementations.
 *
 * @see <a href="http://www.boost.org/doc/libs/release/libs/multi_array">
 *      Boost.MultiArray</a> for more information on the MultiArray concept.
 */
namespace multi_array
{

namespace { // anonymous

template<std::size_t D>
struct for_each_functor {

    BOOST_STATIC_ASSERT(D != 0); // Nonsensical behavior for zero dimensions
    BOOST_STATIC_ASSERT(D > 1);  // Ensure not instantiated for specialized values

    // Cannot use Array reference in signature below, why?  See
    // http://groups.google.com/group/boost-list/browse_thread/thread/e16f32c4411dea08
    // for any potential changes to the reference-related equations below
    template<class Array,
             class UnaryFunction >
    void operator()(Array /*&*/x, UnaryFunction &f) const {
        for_each_functor<D-1> functor;
        for (typename Array::iterator i = x.begin(); i != x.end(); ++i) {
            functor(*i,f);
        }
    }
};

template<> struct for_each_functor<1> {
    // Cannot use Array reference in signature below, why?  See
    // http://groups.google.com/group/boost-list/browse_thread/thread/e16f32c4411dea08
    // for any potential changes to the reference-related equations below
    template<class Array,
             class UnaryFunction >
    void operator()(Array /*&*/x, UnaryFunction &f) const {
        std::for_each(x.begin(), x.end(), f);
    }
};

} // namespace anonymous

template<class Array,
         class UnaryFunction>
inline
UnaryFunction for_each(Array &x,
                       UnaryFunction f) {
    for_each_functor<Array::dimensionality>()(x,f);
    return f;
}

template<class ValueType, std::size_t NumDims, class Allocator,
         class UnaryFunction>
inline
UnaryFunction for_each(boost::multi_array<ValueType,NumDims,Allocator> &x,
                       UnaryFunction f) {
    return std::for_each(x.data(), x.data() + x.num_elements(), f);
}

template<class ValueType, std::size_t NumDims,
         class UnaryFunction>
inline
UnaryFunction for_each(boost::multi_array_ref<ValueType,NumDims> &x,
                       UnaryFunction f)
{
    return std::for_each(x.data(), x.data() + x.num_elements(), f);
}

namespace { // anonymous

template<class S>
struct assign_value_functor {
    assign_value_functor(const S &s) : s_(s) {}

    template<class T> void operator()(T& t) const { t = s_; }

private:
    const S &s_;
};

} // namespace anonymous

template<class Array, class V>
void fill(Array &x, const V &v) {
    for_each(x, assign_value_functor<V>(v));
}

} // namespace multi_array

} // namespace suzerain

#endif // __SUZERAIN_MULTIARRAY_H
