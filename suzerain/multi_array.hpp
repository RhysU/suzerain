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
 * multi_array.hpp: Utilities for Boost.MultiArray
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_MULTIARRAY_H
#define __SUZERAIN_MULTIARRAY_H

#include <suzerain/common.hpp>
#include <suzerain/functional.hpp>
#include <suzerain/storage.hpp>

namespace suzerain
{

/**
 * Provides utilities atop the Boost.MultiArray concept and its concrete
 * implementations.
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

/**
 * Obtain a <tt>boost::array</tt> containing the stride information from a
 * MultiArray implementation.
 *
 * @param a MultiArray to interrogate.
 *
 * @return A copy of the MultiArray's stride information with a length
 *         equal to the MultiArray's dimensionality.
 */
template< typename MultiArray >
boost::array<typename MultiArray::index,MultiArray::dimensionality>
    strides_array(const MultiArray &a)
{
    typedef typename MultiArray::index index;
    boost::array<index,MultiArray::dimensionality> retval;
    const index *stridesBegin = a.strides();
    std::copy(stridesBegin, stridesBegin + MultiArray::dimensionality,
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

/**
 * A <tt>boost::multi_array_ref</tt>-like class for strided indexing into a
 * contiguous region of memory provided to the constructor.
 *
 * @internal Nearly all of the implementation is inherited, but we choose to
 * avoid publicly being a <tt>multi_array_ref</tt> subclass since the
 * underlying elements are not strictly contiguous.
 *
 * @see <a href="http://www.boost.org/doc/libs/release/libs/multi_array">
 *      Boost.MultiArray</a> for more information on the MultiArray concept.
 */
template<typename ValueType, std::size_t NumDims>
class ref
    : private boost::multi_array_ref<ValueType, NumDims>
{
private:
    typedef typename boost::multi_array_ref<ValueType, NumDims> base;

public:

    // MultiArray Table 2. Associated Types
    typedef typename base::value_type             value_type;
    typedef typename base::reference              reference;
    typedef typename base::const_reference        const_reference;
    typedef typename base::size_type              size_type;
    typedef typename base::difference_type        difference_type;
    typedef typename base::iterator               iterator;
    typedef typename base::const_iterator         const_iterator;
    typedef typename base::reverse_iterator       reverse_iterator;
    typedef typename base::const_reverse_iterator const_reverse_iterator;
    typedef typename base::element                element;
    typedef typename base::index                  index;
    typedef typename base::index_gen              index_gen;
    typedef typename base::index_range            index_range;
    using base::subarray;
    using base::const_subarray;
    using base::array_view;
    using base::const_array_view;

    // MultiArray Table 3. Valid Expressions
    using base::dimensionality;
    using base::shape;
    using base::strides;
    using base::index_bases;
    using base::origin;
    using base::num_dimensions;
    using base::num_elements;
    using base::size;
    using base::operator();
    using base::begin;
    using base::end;
    using base::rbegin;
    using base::rend;
    using base::operator[];
    using base::operator==;
    using base::operator<;
    using base::operator<=;
    using base::operator>;
    using base::operator>=;

    // Useful miscellany not strictly required by MultiArray
    typedef typename base::extent_gen   extent_gen;
    typedef typename base::extent_range extent_range;
    using base::data;
    using base::storage_order;

    template<typename ExtentList,
             typename StorageOrderSequence>
    explicit ref(
            element* data,
            const ExtentList& sizes,
            const suzerain::storage::general<StorageOrderSequence>& storage)
        : base(data, sizes, storage.storage_order())
    {
        // stride_list_ protected in boost::const_multi_array_ref ancestor
        storage.compute_strides(this->shape(), this->stride_list_.begin());
    }

    template<typename StorageOrderSequence>
    explicit ref(
            element* data,
            typename base::extent_gen ranges,
            const suzerain::storage::general<StorageOrderSequence>& storage)
        : base(data, ranges, storage.storage_order())
    {
        // stride_list_ protected in boost::const_multi_array_ref ancestor
        storage.compute_strides(this->shape(), this->stride_list_.begin());
    }

    template<typename ExtentList,
             typename MinStrideList,
             typename StorageOrderSequence>
    explicit ref(
            element* data,
            const ExtentList& sizes,
            const MinStrideList& minstrides,
            const suzerain::storage::general<StorageOrderSequence>& storage)
        : base(data, sizes, storage.storage_order())
    {
        // stride_list_ protected in boost::const_multi_array_ref ancestor
        storage.compute_strides(this->shape(),
                                minstrides.begin(),
                                this->stride_list_.begin());
    }

    template<typename MinStrideList,
             typename StorageOrderSequence>
    explicit ref(
            element* data,
            typename base::extent_gen ranges,
            const MinStrideList& minstrides,
            const suzerain::storage::general<StorageOrderSequence>& storage)
        : base(data, ranges, storage.storage_order())
    {
        // stride_list_ protected in boost::const_multi_array_ref ancestor
        storage.compute_strides(this->shape(),
                                minstrides.begin(),
                                this->stride_list_.begin());
    }

    explicit ref(const ref&  other) : base(other /* shallow */) {}

    explicit ref(const base& other) : base(other /* shallow */) {}

    // Assignment operator provided by 'using base::operator=' above
};

} // namespace multi_array

} // namespace suzerain

#endif // __SUZERAIN_MULTIARRAY_H
