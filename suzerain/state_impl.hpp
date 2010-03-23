//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2010 The PECOS Development Team
//
// Please see http://pecos.ices.utexas.edu for more information.
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
//
// state_impl.hpp: Implementations of the IState interface
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef __SUZERAIN_STATE_IMPL_HPP
#define __SUZERAIN_STATE_IMPL_HPP

#include <suzerain/common.hpp>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/functional.hpp>
#include <suzerain/mpl.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/storage.hpp>
#include <suzerain/state.hpp>

/** @file
 * Provides implementations of the IState interface.
 */

namespace suzerain
{

template<
    typename Element,
    typename Allocator = std::allocator<Element>
>
class RawMemory
{
public:
    template< typename I >
    RawMemory(
        I count,
        typename boost::enable_if<boost::is_integral<I> >::type* dummy = 0)
        : count_(boost::numeric_cast<std::size_t>(count)),
          a_(Allocator()),
          p_(a_.allocate(count)) {}

    RawMemory(const RawMemory &other)
        : count_(other.count_), a_(other.a_), p_(a_.allocate(count_))
    {
        std::memcpy(p_, other.p_, count_*sizeof(Element));
    }

    ~RawMemory() { a_.deallocate(p_, count_); }

    typename Allocator::pointer raw_memory() const {
        return typename Allocator::pointer(p_); // defensive copy
    }

    std::size_t raw_memory_count() const { return count_; }

private:
    const RawMemory& operator=( const RawMemory& ); // Disable

    const std::size_t count_;
    Allocator a_;
    typename Allocator::pointer p_;
};

// Forward declarations
template< std::size_t NumDims, typename Element, typename Allocator >
    class InterleavedState;
template< std::size_t NumDims, typename Element, typename Allocator >
    class NoninterleavedState;

template<
    std::size_t NumDims,
    typename Element,
    typename Allocator = typename suzerain::blas::allocator<Element>::type
>
class InterleavedState
    : public IState<NumDims,Element,suzerain::storage::interleaved<NumDims> >,
      public RawMemory<Element,Allocator>,
      public boost::multi_array_ref<Element, NumDims>
{
public:
    typedef typename boost::multi_array_ref<Element, NumDims>
            multi_array_type;
    typedef typename suzerain::storage::interleaved<NumDims>
            storage_interleaved;

    template< typename ExtentList >
    explicit InterleavedState(
            const ExtentList& sizes,
            typename Allocator::size_type min_total_contiguous_count = 0);

    InterleavedState(const InterleavedState& other);

    virtual ~InterleavedState() {}

    virtual void scale(const Element& factor);

    virtual bool isConformant(
        const IState<NumDims,Element,storage_interleaved>& other) const
        throw(std::bad_cast);

    virtual void addScaled(
            const Element& factor,
            const IState<NumDims,Element,storage_interleaved>& other)
            throw(std::bad_cast, std::logic_error);

    virtual void assign(
            const IState<NumDims,Element,storage_interleaved>& other)
            throw(std::bad_cast, std::logic_error);

    virtual void exchange(
            IState<NumDims,Element,storage_interleaved>& other)
            throw(std::bad_cast, std::logic_error);

private:
    // Disable assignment operators
    const multi_array_type& operator=( const multi_array_type& );
    const InterleavedState& operator=( const InterleavedState& );
};

template< std::size_t NumDims, typename Element, typename Allocator >
template< typename ExtentList >
InterleavedState<NumDims,Element,Allocator>::InterleavedState(
        const ExtentList& sizes,
        typename Allocator::size_type min_total_contiguous_count)
    : IStateBase<NumDims,Element>(),
      IState<NumDims,Element,storage_interleaved>(),
      RawMemory<Element,Allocator>(std::max<typename Allocator::size_type>(
                    suzerain::functional::product(sizes.begin(), sizes.end()),
                    min_total_contiguous_count)),
      multi_array_type(RawMemory<Element,Allocator>::raw_memory(),
                       sizes,
                       storage_interleaved::storage_order())
{
    // NOP
}

template< std::size_t NumDims, typename Element, typename Allocator >
InterleavedState<NumDims,Element,Allocator>::InterleavedState(
        const InterleavedState &other)
    : IStateBase<NumDims,Element>(other),
      IState<NumDims,Element,storage_interleaved>(other),
      RawMemory<Element,Allocator>(other),
      multi_array_type(RawMemory<Element,Allocator>::raw_memory(),
                       suzerain::multi_array::shape_array(other),
                       storage_interleaved::storage_order())
{
    // Data copied by RawMemory's copy constructor
}

template< std::size_t NumDims, typename Element, typename Allocator >
void InterleavedState<NumDims,Element,Allocator>::scale(
        const Element& factor)
{
    suzerain::blas::scal(this->num_elements(), factor, this->data(), 1);
}

template< std::size_t NumDims, typename Element, typename Allocator >
bool InterleavedState<NumDims,Element,Allocator>::isConformant(
            const IState<NumDims,Element,storage_interleaved>& other) const
throw(std::bad_cast)
{
    const InterleavedState& o = dynamic_cast<const InterleavedState&>(other);

    return std::equal(this->shape(), this->shape() + NumDims, o.shape());
}

template< std::size_t NumDims, typename Element, typename Allocator >
void InterleavedState<NumDims,Element,Allocator>::addScaled(
            const Element& factor,
            const IState<NumDims,Element,storage_interleaved>& other)
throw(std::bad_cast, std::logic_error)
{
    if (SUZERAIN_UNLIKELY(this == boost::addressof(other)))
        throw std::logic_error("Unable to handle this->addScaled(...,this)");

    if (SUZERAIN_UNLIKELY(!isConformant(other))) throw std::logic_error(
            std::string("Nonconformant other in ") + __PRETTY_FUNCTION__);

    const InterleavedState& o = dynamic_cast<const InterleavedState&>(other);
    assert(o.num_elements() == this->num_elements());

    suzerain::blas::axpy(
            this->num_elements(), factor, o.data(), 1, this->data(), 1);
}

template< std::size_t NumDims, typename Element, typename Allocator >
void InterleavedState<NumDims,Element,Allocator>::assign(
            const IState<NumDims,Element,storage_interleaved>& other)
throw(std::bad_cast, std::logic_error)
{
    if (SUZERAIN_UNLIKELY(this == boost::addressof(other))) return; // Self?

    if (SUZERAIN_UNLIKELY(!isConformant(other))) throw std::logic_error(
            std::string("Nonconformant other in ") + __PRETTY_FUNCTION__);

    const InterleavedState& o = dynamic_cast<const InterleavedState&>(other);
    assert(o.num_elements() == this->num_elements());

    suzerain::blas::copy(
            this->num_elements(), o.data(), 1, this->data(), 1);
}

template< std::size_t NumDims, typename Element, typename Allocator >
void InterleavedState<NumDims,Element,Allocator>::exchange(
            IState<NumDims,Element,storage_interleaved>& other)
throw(std::bad_cast, std::logic_error)
{
    if (SUZERAIN_UNLIKELY(this == boost::addressof(other))) return; // Self?

    if (SUZERAIN_UNLIKELY(!isConformant(other))) throw std::logic_error(
            std::string("Nonconformant other in ") + __PRETTY_FUNCTION__);

    InterleavedState& o = dynamic_cast<InterleavedState&>(other);
    assert(o.num_elements() == this->num_elements());

    suzerain::blas::swap(
            this->num_elements(), o.data(), 1, this->data(), 1);
}

template<
    std::size_t NumDims,
    typename Element,
    typename Allocator = typename suzerain::blas::allocator<Element>::type
>
class NoninterleavedState
    : public IState<NumDims,
                    Element,
                    suzerain::storage::noninterleaved<NumDims> >,
      public RawMemory<Element,Allocator>,
      public boost::detail::multi_array::multi_array_view<Element, NumDims>
{
public:
    typedef typename
        boost::detail::multi_array::multi_array_view<Element, NumDims>
        multi_array_type;
    typedef typename
        suzerain::storage::noninterleaved<NumDims>
        storage_noninterleaved;

    template< typename ExtentList >
    explicit NoninterleavedState(
            const ExtentList& sizes,
            typename Allocator::size_type min_variable_stride = 0);

    NoninterleavedState(const NoninterleavedState& other);

    virtual ~NoninterleavedState() {}

    virtual void scale(const Element& factor);

    virtual bool isConformant(
        const IState<NumDims,Element,storage_noninterleaved>& other) const
        throw(std::bad_cast);

    virtual void addScaled(
            const Element& factor,
            const IState<NumDims,Element,storage_noninterleaved>& other)
            throw(std::bad_cast, std::logic_error);

    virtual void assign(
            const IState<NumDims,Element,storage_noninterleaved>& other)
            throw(std::bad_cast, std::logic_error);

    virtual void exchange(
            IState<NumDims,Element,storage_noninterleaved>& other)
            throw(std::bad_cast, std::logic_error);

private:
    // Disable assignment operators
    const multi_array_type& operator=( const multi_array_type& );
    const NoninterleavedState& operator=( const NoninterleavedState& );

    // Compute stride information, including variable padding requirements
    template< typename ExtentList >
    static boost::array<typename multi_array_type::index,NumDims>
        computeStrides(
            const ExtentList& extents_list,
            typename Allocator::size_type min_variable_stride = 0);

    template< typename ExtentList >
    static std::size_t computeRawMemoryCount(
            const ExtentList& extents_list,
            typename Allocator::size_type min_variable_stride = 0);
};

template< std::size_t NumDims, typename Element, typename Allocator >
template< typename ExtentList >
boost::array< typename NoninterleavedState<
        NumDims,Element,Allocator
    >::multi_array_type::index,NumDims>
NoninterleavedState<NumDims,Element,Allocator>::computeStrides(
        const ExtentList &extents,
        typename Allocator::size_type min_variable_stride)
{
    // Slow, but we need random access to extents_list's contents
    boost::array<typename ExtentList::value_type,NumDims> extents;
    std::copy(extents.begin(), extents.end(), extents_list.begin());

    // Obtain storage ordering type information in runtime array
    const suzerain::mpl::sequence_array<storage_noninterleaved> storage;
    assert(storage.size() == stride.size());
    assert(storage.size > 0);
    assert(storage[NumDims-1] == 0);

    // Compute the usual stride information by the usual means
    boost::array<index,NumDims> stride;
    stride[storage[0]] = 1;
    for (index i = 1; i < NumDims; ++i) {
        stride[storage[i]] = stride[storage[i-1]]*extents[storage[i-1]];
    }
    // Enforce the minimum stride requirement on the final stride
    stride[storage[NumDims-1]] = std::max(
            boost::numeric_cast<index>(min_variable_stride),
            stride[storage[NumDims-1]]);
}

template< std::size_t NumDims, typename Element, typename Allocator >
template< typename ExtentList >
std::size_t
NoninterleavedState<NumDims,Element,Allocator>::computeRawMemoryCount(
        const ExtentList &extents,
        typename Allocator::size_type min_variable_stride)
{
    assert(extents_list.begin() != extents_list.end())

    using suzerain::functional::product;
    const std::size_t nonpadded_count
        = product(extents_list.begin(), extents_list.end());
    const std::size_t padded_count
        = *(extents_list.begin()) * min_variable_stride;

    return std::max(padded, nonpadded);
}

} // namespace suzerain

#endif // __SUZERAIN_STATE_IMPL_HPP
