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
#include <suzerain/functional.hpp>

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

/**
 * Obtain a <tt>boost::array</tt> containing the shape information from a
 * MultiArray implementation.
 *
 * @param a MultiArray to interrogate.
 *
 * @return A copy of the MultiArray's shape information with a length
 *         equal to the MultiArray's dimensionality.
 */
template< typename MultiArray >
boost::array<typename MultiArray::size_type,MultiArray::dimensionality>
    shape_array(const MultiArray &a)
{
    typedef typename MultiArray::size_type size_type;
    boost::array<size_type,MultiArray::dimensionality> retval;
    const size_type *shapeBegin = a.shape();
    std::copy(shapeBegin, shapeBegin + MultiArray::dimensionality,
              retval.begin());
    return retval;
}

namespace { // anonymous

template<std::size_t D>
struct for_each_functor {

    BOOST_STATIC_ASSERT(D != 0); // Nonsensical behavior for zero dimensions
    BOOST_STATIC_ASSERT(D > 1);  // Ensure not instantiated for special values

    // See http://groups.google.com/group/boost-list/browse_thread/thread/e16f32c4411dea08
    // for details about why MultiArray::iterator::reference is used below.
    template<class MultiArray,
             class UnaryFunction >
    void operator()(MultiArray &x, UnaryFunction &f) const {
        for_each_functor<D-1> functor;
        for (typename MultiArray::iterator i = x.begin(); i != x.end(); ++i) {
            typename MultiArray::iterator::reference ri = *i;
            functor(ri,f);
        }
    }
};

template<> struct for_each_functor<1> {
    template<class MultiArray,
             class UnaryFunction >
    void operator()(MultiArray &x, UnaryFunction &f) const {
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

/**
 * Fill MultiArray \c x with the value \c v.  MultiArray <tt>x</tt>'s elements
 * must be assignable from \c v.  The underlying assignment uses
 * ::suzerain::functional::assign to allow specializations of that functor to
 * be found.
 *
 * @param x MultiArray to fill.
 * @param v Value with which to fill \c x.
 */
template<class MultiArray, class V>
void fill(MultiArray &x, const V &v) {
    using suzerain::functional::assign;
    for_each(x, assign<typename MultiArray::element,V>(v));
}

} // namespace multi_array

} // namespace suzerain

#endif // __SUZERAIN_MULTIARRAY_H
