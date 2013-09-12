//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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

#ifndef SUZERAIN_MULTIARRAY_HPP
#define SUZERAIN_MULTIARRAY_HPP

/** @file
 * Utilities for Boost.MultiArray.
 */

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
 * Obtain an array containing the shape information from a MultiArray
 * implementation.
 *
 * @param a MultiArray to interrogate.
 *
 * @return A copy of the MultiArray's shape information with a length
 *         equal to the MultiArray's dimensionality.
 */
template< typename MultiArray >
array<typename MultiArray::size_type,MultiArray::dimensionality>
    shape_array(const MultiArray &a)
{
    typedef typename MultiArray::size_type size_type;
    array<size_type,MultiArray::dimensionality> retval;
    const size_type *shapeBegin = a.shape();
    std::copy(shapeBegin, shapeBegin + MultiArray::dimensionality,
              retval.begin());
    return retval;
}

/**
 * Obtain an array containing the stride information from a MultiArray
 * implementation.
 *
 * @param a MultiArray to interrogate.
 *
 * @return A copy of the MultiArray's stride information with a length
 *         equal to the MultiArray's dimensionality.
 */
template< typename MultiArray >
array<typename MultiArray::index,MultiArray::dimensionality>
    strides_array(const MultiArray &a)
{
    typedef typename MultiArray::index index;
    array<index,MultiArray::dimensionality> retval;
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

// Detects storage order to walk memory linearly for 2D MultiArrays
template<> struct for_each_functor<2> {

    template<class MultiArray,
             class UnaryFunction >
    void operator()(MultiArray &x, UnaryFunction &f) const {

        typedef typename MultiArray::index     index;
        typedef typename MultiArray::size_type size_type;

        // Compute lower and upper index bounds
        const index     * const base   = x.index_bases();
        const size_type * const shape  = x.shape();
        const index il = base[0], iu = il + (index) shape[0];
        const index jl = base[1], ju = jl + (index) shape[1];

        // Obtain absolute strides for storage order determination
        const index * const stride = x.strides();
        const index si = std::abs(stride[0]);
        const index sj = std::abs(stride[1]);

        // Traverse according to runtime-determined storage order
        if (si >= sj) {
            for (index i = il; i < iu; ++i)
            for (index j = jl; j < ju; ++j)
                f(x[i][j]);
        } else {
            for (index j = jl; j < ju; ++j)
            for (index i = il; i < iu; ++i)
                f(x[i][j]);
        }

    }

};

// Detects storage order to walk memory linearly for 3D MultiArrays
template<> struct for_each_functor<3> {

    template<class MultiArray,
             class UnaryFunction >
    void operator()(MultiArray &x, UnaryFunction &f) const {

        typedef typename MultiArray::index     index;
        typedef typename MultiArray::size_type size_type;

        // Compute lower and upper index bounds
        const index     * const base   = x.index_bases();
        const size_type * const shape  = x.shape();
        const index il = base[0], iu = il + (index) shape[0];
        const index jl = base[1], ju = jl + (index) shape[1];
        const index kl = base[2], ku = kl + (index) shape[2];

        // Obtain absolute strides for storage order determination
        const index * const stride = x.strides();
        const index si = std::abs(stride[0]);
        const index sj = std::abs(stride[1]);
        const index sk = std::abs(stride[2]);

        // Perform all possible stride comparisons
        const bool oij = si >= sj, oji = !oij;
        const bool oik = si >= sk, oki = !oik;
        const bool ojk = sj >= sk, okj = !ojk;

#ifndef SUZERAIN_PARSED_BY_DOXYGEN
#define FOR(ndx) for (index ndx = ndx##l; ndx < ndx##u; ++ndx)
#endif

        // Traverse according to runtime-determined storage order
        if      (oij && ojk) // boost::c_storage_order
            FOR(i) FOR(j) FOR(k) f(x[i][j][k]);
        else if (okj && oji) // boost::fortran_storage_order
            FOR(k) FOR(j) FOR(i) f(x[i][j][k]);
        else if (oki && oij) // suzerain::storage::interleaved<3>
            FOR(k) FOR(i) FOR(j) f(x[i][j][k]);
        else if (oik && okj) // suzerain::storage::contiguous<3>
            FOR(i) FOR(k) FOR(j) f(x[i][j][k]);
        else if (oji && oik) // lexicographic from here onward...
            FOR(j) FOR(i) FOR(k) f(x[i][j][k]);
        else // (ojk && oki)
            FOR(j) FOR(k) FOR(i) f(x[i][j][k]);

#ifndef SUZERAIN_PARSED_BY_DOXYGEN
#undef FOR
#endif

    }

};

// Detects storage order to walk memory linearly for 4D MultiArrays
template<> struct for_each_functor<4> {

    template<class MultiArray,
             class UnaryFunction >
    void operator()(MultiArray &x, UnaryFunction &f) const {

        typedef typename MultiArray::index     index;
        typedef typename MultiArray::size_type size_type;

        // Compute lower and upper index bounds
        const index     * const base   = x.index_bases();
        const size_type * const shape  = x.shape();
        const index il = base[0], iu = il + (index) shape[0];
        const index jl = base[1], ju = jl + (index) shape[1];
        const index kl = base[2], ku = kl + (index) shape[2];
        const index ll = base[3], lu = ll + (index) shape[3];

        // Obtain absolute strides for storage order determination
        const index * const stride = x.strides();
        const index si = std::abs(stride[0]);
        const index sj = std::abs(stride[1]);
        const index sk = std::abs(stride[2]);
        const index sl = std::abs(stride[3]);

        // Perform all possible stride comparisons
        const bool oij = si >= sj, oji = !oij;
        const bool oik = si >= sk, oki = !oik;
        const bool oil = si >= sl, oli = !oil;
        const bool ojk = sj >= sk, okj = !ojk;
        const bool ojl = sj >= sl, olj = !ojl;
        const bool okl = sk >= sl, olk = !okl;

#ifndef SUZERAIN_PARSED_BY_DOXYGEN
#define FOR(ndx) for (index ndx = ndx##l; ndx < ndx##u; ++ndx)
#endif

        // Traverse according to runtime-determined storage order
        if      (oij && ojk && okl) // boost::c_storage_order
            FOR(i) FOR(j) FOR(k) FOR(l) f(x[i][j][k][l]);
        else if (olk && okj && oji) // boost::fortran_storage_order
            FOR(l) FOR(k) FOR(j) FOR(i) f(x[i][j][k][l]);
        else if (olk && oki && oij) // suzerain::storage::interleaved<4>
            FOR(l) FOR(k) FOR(i) FOR(j) f(x[i][j][k][l]);
        else if (oil && olk && okj) // suzerain::storage::contiguous<4>
            FOR(i) FOR(l) FOR(k) FOR(j) f(x[i][j][k][l]);
        else if (oij && ojl && olk) // lexicographic from here onward...
            FOR(i) FOR(j) FOR(l) FOR(k) f(x[i][j][k][l]);
        else if (oik && okj && ojl)
            FOR(i) FOR(k) FOR(j) FOR(l) f(x[i][j][k][l]);
        else if (oik && okl && olj)
            FOR(i) FOR(k) FOR(l) FOR(j) f(x[i][j][k][l]);
        else if (oil && olj && ojk)
            FOR(i) FOR(l) FOR(j) FOR(k) f(x[i][j][k][l]);
        else if (oji && oik && okl)
            FOR(j) FOR(i) FOR(k) FOR(l) f(x[i][j][k][l]);
        else if (oji && oil && olk)
            FOR(j) FOR(i) FOR(l) FOR(k) f(x[i][j][k][l]);
        else if (ojk && oki && oil)
            FOR(j) FOR(k) FOR(i) FOR(l) f(x[i][j][k][l]);
        else if (ojk && okl && oli)
            FOR(j) FOR(k) FOR(l) FOR(i) f(x[i][j][k][l]);
        else if (ojl && oli && oik)
            FOR(j) FOR(l) FOR(i) FOR(k) f(x[i][j][k][l]);
        else if (ojl && olk && oki)
            FOR(j) FOR(l) FOR(k) FOR(i) f(x[i][j][k][l]);
        else if (oki && oij && ojl)
            FOR(k) FOR(i) FOR(j) FOR(l) f(x[i][j][k][l]);
        else if (oki && oil && olj)
            FOR(k) FOR(i) FOR(l) FOR(j) f(x[i][j][k][l]);
        else if (okj && oji && oil)
            FOR(k) FOR(j) FOR(i) FOR(l) f(x[i][j][k][l]);
        else if (okj && ojl && oli)
            FOR(k) FOR(j) FOR(l) FOR(i) f(x[i][j][k][l]);
        else if (okl && oli && oij)
            FOR(k) FOR(l) FOR(i) FOR(j) f(x[i][j][k][l]);
        else if (okl && olj && oji)
            FOR(k) FOR(l) FOR(j) FOR(i) f(x[i][j][k][l]);
        else if (oli && oij && ojk)
            FOR(l) FOR(i) FOR(j) FOR(k) f(x[i][j][k][l]);
        else if (oli && oik && okj)
            FOR(l) FOR(i) FOR(k) FOR(j) f(x[i][j][k][l]);
        else if (olj && oji && oik)
            FOR(l) FOR(j) FOR(i) FOR(k) f(x[i][j][k][l]);
        else // (olj && ojk && oki)
            FOR(l) FOR(j) FOR(k) FOR(i) f(x[i][j][k][l]);

#ifndef SUZERAIN_PARSED_BY_DOXYGEN
#undef FOR
#endif

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
                       UnaryFunction f)
{
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
                       UnaryFunction f)
{
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
 * suzerain::functional::assign to allow specializations of that functor to
 * be found.
 *
 * @param x MultiArray to fill.
 * @param v Value with which to fill \c x.
 */
template<class MultiArray, class V>
void fill(MultiArray &x, const V &v)
{
    using suzerain::functional::assign;
    for_each(x, assign<typename MultiArray::element,V>(v));
}

/**
 * Is the given MultiArray contiguous in memory?
 *
 * A MultiArray is contiguous when the <tt>product(shape(), shape() +
 * dimensionality)</tt> equals the distance between the first and last element
 * according to <tt>index_bases()</tt>, <tt>strides()</tt> and
 * <tt>shape()</tt>.
 *
 *
 * @param x MultiArray to check for contiguity.
 *
 * @return True if \c is contiguous.  False otherwise.
 *
 * @see SGI's <a href="http://www.sgi.com/tech/stl/UnaryFunction.html">
 *      UnaryFunction</a> concept for more information.
 */
template<class MultiArray>
bool is_contiguous(const MultiArray &x)
{
    // Using ptrdiff_t avoids unsigned/signed warnings when dimensionality == 0.
    static const std::ptrdiff_t dimensionality = MultiArray::dimensionality;

    typename MultiArray::size_type min_extent = (dimensionality == 0) ? 0 : 1;
    typename MultiArray::size_type abs_extent = 1;

    // Obtain minimum possible extent and actual extent using abs(strides).
    // (shape()[i] == 0) overflows abs_extent but min_extent becomes zero.
    for (std::ptrdiff_t i = 0; i < dimensionality; ++i) {
        const typename MultiArray::size_type shape = x.shape()[i];
        min_extent *= shape;
        abs_extent += (shape - 1)*std::abs(x.strides()[i]);
    }

    return (min_extent == 0) ? true : (min_extent == abs_extent);
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

    /**@name MultiArray Associated Types
     * Following Table 2 in the MultiArray concept specification, the
     * following types are available:
     *@{
     */
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
    /*@}*/

    /**@name MultiArray Valid Expressions
     * Following Table 3 in the MultiArray concept specification, the
     * following expressions are valid:
     *@{
     */
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
    /*@}*/

    // Useful miscellany not strictly required by MultiArray
    /**@name MultiArray Extensions
     * Though not strictly required by the MultiArray concept, the following
     * queries are useful and are also provided by <tt>boost::multi_array</tt>
     * and <tt>boost::multi_array_ref</tt>:
     *@{
     */
    typedef typename base::extent_gen         extent_gen;
    typedef typename base::extent_range       extent_range;
    typedef typename base::storage_order_type storage_order_type;
    using base::data;
    using base::storage_order;
    /*@}*/

    /**
     * Construct an instance wrapping \c data with dimensions \c sizes
     * stored according to \c storage.
     *
     * @param data    Raw elements to wrap
     * @param sizes   A length ref::dimensionality collection of the
     *                valid extents
     * @param storage A marker type specifying the storage ordering
     */
    template<typename ExtentList,
             typename StorageOrderSequence>
    explicit ref(
            element* data,
            const ExtentList& sizes,
            const storage::general<StorageOrderSequence>& storage)
        : base(data, sizes, storage.storage_order())
    {
        // stride_list_ protected in boost::const_multi_array_ref ancestor
        storage.compute_strides(this->shape(), this->stride_list_.begin());
    }

    /**
     * Construct an instance wrapping \c data with sizes and indexing
     * requirements based on \c ranges stored according to \c storage.
     *
     * @param data    Raw elements to wrap
     * @param ranges  Index range and size information for \c data
     * @param storage A marker type specifying the storage ordering
     */
    template<typename StorageOrderSequence>
    explicit ref(
            element* data,
            typename extent_gen::template gen_type<NumDims>::type &ranges,
            const storage::general<StorageOrderSequence>& storage)
        : base(data, ranges, storage.storage_order())
    {
        // stride_list_ protected in boost::const_multi_array_ref ancestor
        storage.compute_strides(this->shape(), this->stride_list_.begin());
    }

    /**
     * Construct an instance wrapping \c data with dimensions \c sizes
     * stored obeying \c minstrides according to \c storage.
     *
     * @param data       Raw elements to wrap
     * @param sizes      A length ref::dimensionality collection
     *                   of the valid extents
     * @param minstrides A length ref::dimensionality collection
     *                   of the minimum acceptable strides
     * @param storage    A marker type specifying the storage ordering
     */
    template<typename ExtentList,
             typename MinStrideList,
             typename StorageOrderSequence>
    explicit ref(
            element* data,
            const ExtentList& sizes,
            const MinStrideList& minstrides,
            const storage::general<StorageOrderSequence>& storage)
        : base(data, sizes, storage.storage_order())
    {
        // stride_list_ protected in boost::const_multi_array_ref ancestor
        storage.compute_strides(this->shape(),
                                minstrides.begin(),
                                this->stride_list_.begin());
    }

    /**
     * Construct an instance wrapping \c data with sizes and indexing
     * requirements based on \c ranges stored obeying \c minstrides according
     * to \c storage.
     *
     * @param data       Raw elements to wrap
     * @param ranges     Index range and size information for \c data
     * @param minstrides A length ref::dimensionality collection
     *                   of the minimum acceptable strides
     * @param storage    A marker type specifying the storage ordering
     */
    template<typename MinStrideList,
             typename StorageOrderSequence>
    explicit ref(
            element* data,
            typename extent_gen::template gen_type<NumDims>::type &ranges,
            const MinStrideList& minstrides,
            const storage::general<StorageOrderSequence>& storage)
        : base(data, ranges, storage.storage_order())
    {
        // stride_list_ protected in boost::const_multi_array_ref ancestor
        storage.compute_strides(this->shape(),
                                minstrides.begin(),
                                this->stride_list_.begin());
    }

    /**
     * Construct a shallow copy of \c other.
     *
     * @param other Instance to mimic.
     */
    explicit ref(const ref& other)
        : base(other /* shallow */)
    {
        // stride_list_ protected in boost::const_multi_array_ref ancestor
        std::copy(other.strides(), other.strides() + other.dimensionality,
                  this->stride_list_.begin());
    }

    /**
     * Construct a shallow copy of \c other.
     *
     * @param other Instance to mimic.
     */
    explicit ref(const base& other)
        : base(other /* shallow */)
    {
        // stride_list_ protected in boost::const_multi_array_ref ancestor
        std::copy(other.strides(), other.strides() + other.dimensionality,
                  this->stride_list_.begin());
    }

    // Assignment operator provided by 'using base::operator=' above
    // base::operator= seems to appropriately copy stride information
};

} // namespace multi_array

} // namespace suzerain

#endif // SUZERAIN_MULTIARRAY_HPP
