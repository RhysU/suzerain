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
#include <suzerain/complex.hpp>

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

    // Cannot use MultiArray reference in signature below, why?  See
    // http://groups.google.com/group/boost-list/browse_thread/thread/e16f32c4411dea08
    // for any potential changes to the reference-related equations below
    template<class MultiArray,
             class UnaryFunction >
    void operator()(MultiArray /*&*/x, UnaryFunction &f) const {
        for_each_functor<D-1> functor;
        for (typename MultiArray::iterator i = x.begin(); i != x.end(); ++i) {
            functor(*i,f);
        }
    }
};

template<> struct for_each_functor<1> {
    // Cannot use MultiArray reference in signature below, why?  See
    // http://groups.google.com/group/boost-list/browse_thread/thread/e16f32c4411dea08
    // for any potential changes to the reference-related equations below
    template<class MultiArray,
             class UnaryFunction >
    void operator()(MultiArray /*&*/x, UnaryFunction &f) const {
        std::for_each(x.begin(), x.end(), f);
    }
};

} // namespace anonymous

/**
 * Invoke \c f for each element in MultiArray \c x.  The
 * order in which the invocations take place is undefined.
 *
 * @param x MultiArray over which to iterate.
 * @param f UnaryFunction to invoke on each element of \c x.
 *
 * @return \c f.
 *
 * @see SGI's <a href="http://www.sgi.com/tech/stl/UnaryFunction.html">
 *      UnaryFunction</a> concept for more information.
 */
template<class MultiArray,
         class UnaryFunction>
inline
UnaryFunction for_each(MultiArray &x,
                       UnaryFunction f) {
    for_each_functor<MultiArray::dimensionality>()(x,f);
    return f;
}

/**
 * Invoke \c f for each element in <tt>boost::multi_array</tt> \c x.  The order
 * in which the invocations take place is undefined.  This specialization
 * takes advantage of <tt>boost::multi_array</tt>'s contiguous storage.
 *
 * @param x <tt>boost::multi_array</tt> over which to iterate.
 * @param f UnaryFunction to invoke on each element of \c x.
 *
 * @return \c f.
 *
 * @see SGI's <a href="http://www.sgi.com/tech/stl/UnaryFunction.html">
 *      UnaryFunction</a> concept for more information.
 */
template<class ValueType, std::size_t NumDims, class Allocator,
         class UnaryFunction>
inline
UnaryFunction for_each(boost::multi_array<ValueType,NumDims,Allocator> &x,
                       UnaryFunction f) {
    return std::for_each(x.data(), x.data() + x.num_elements(), f);
}

/**
 * Invoke \c f for each element in <tt>boost::multi_array_ref</tt> \c x.  The
 * order in which the invocations take place is undefined.  This specialization
 * takes advantage of <tt>boost::multi_array_ref</tt>'s contiguous storage.
 *
 * @param x <tt>boost::multi_array_ref</tt> over which to iterate.
 * @param f UnaryFunction to invoke on each element of \c x.
 *
 * @return \c f.
 *
 * @see SGI's <a href="http://www.sgi.com/tech/stl/UnaryFunction.html">
 *      UnaryFunction</a> concept for more information.
 */
template<class ValueType, std::size_t NumDims,
         class UnaryFunction>
inline
UnaryFunction for_each(boost::multi_array_ref<ValueType,NumDims> &x,
                       UnaryFunction f)
{
    return std::for_each(x.data(), x.data() + x.num_elements(), f);
}

namespace { // anonymous

// General assignment functor
template<class Source, class Enable = void>
struct assign_functor {
    assign_functor(const Source &s) : s_(s) {}

    template<class Target> void operator()(Target& t) const { t = s_; }

private:
    const Source &s_;
};

// Specialized assignment functor only accepting recognized complex types
template<class Source>
struct assign_functor<
    Source,
    typename boost::enable_if<
        ::suzerain::complex::traits::is_complex<Source>
    >::type >
{
    assign_functor(const Source &s) : s_(s) {};

    template<class Target>
    void operator()(Target& t) const {
        BOOST_STATIC_ASSERT(::suzerain::complex::traits::is_complex<Target>::value);
        ::suzerain::complex::assign_complex(t, s_);
    }

private:
    const Source &s_;
};

} // namespace anonymous

/**
 * Fill MultiArray \c x with the value \c v.  MultiArray <tt>x</tt>'s elements
 * must be assignable from \c v.
 *
 * @param x MultiArray to fill.
 * @param v Value with which to fill \c x.
 */
template<class MultiArray, class V>
void fill(MultiArray &x, const V &v) {
    for_each(x, assign_functor<V>(v));
}

/**
 * Fill MultiArray \c x with real NaN.
 * MultiArray::element must be a floating point value.
 *
 * @param x MultiArray to fill.
 */
template<class MultiArray>
typename boost::disable_if<
    ::suzerain::complex::traits::is_complex<typename MultiArray::element>,
    void
>::type fill_with_NaN(MultiArray &x) {
    typedef typename MultiArray::element real_type;
    BOOST_STATIC_ASSERT(std::numeric_limits<real_type>::has_quiet_NaN);
    fill(x, std::numeric_limits<real_type>::quiet_NaN());
}

/**
 * Fill MultiArray \c x with complex NaN.  MultiArray::element
 * must be a complex type built from two floating point values.
 *
 * @param x MultiArray to fill.
 */
template<class MultiArray>
typename boost::enable_if<
    ::suzerain::complex::traits::is_complex<typename MultiArray::element>,
    void
>::type fill_with_NaN(MultiArray &x) {
    fill(x, ::suzerain::complex::NaN<typename MultiArray::element>());
}

} // namespace multi_array

} // namespace suzerain

#endif // __SUZERAIN_MULTIARRAY_H
