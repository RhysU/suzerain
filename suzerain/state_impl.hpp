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
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::size_type size_type;

    template< typename I >
    RawMemory(
        I count,
        typename boost::enable_if<boost::is_integral<I> >::type* dummy = 0)
        : count_(boost::numeric_cast<size_type>(count)),
          a_(Allocator()),
          p_(a_.allocate(count)) {}

    RawMemory(const RawMemory &other)
        : count_(other.count_), a_(other.a_), p_(a_.allocate(count_))
    {
        std::memcpy(p_, other.p_, count_*sizeof(Element));
    }

    ~RawMemory() { a_.deallocate(p_, count_); }

    pointer raw_memory() const { return pointer(p_); /* defensive copy */ }

    size_type raw_memory_count() const { return count_; }

private:
    const RawMemory& operator=( const RawMemory& ); // Disable

    const size_type count_;
    Allocator a_;
    pointer p_;
};

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
            typename Allocator::size_type min_contiguous_count = 0);

    InterleavedState(const InterleavedState& other);

    virtual ~InterleavedState() {}

    virtual void scale(const Element &factor);

    virtual bool isConformant(
        const IState<NumDims,Element,storage_interleaved>& other) const
        throw(std::bad_cast);

    virtual void addScaled(
            const Element &factor,
            const IState<NumDims,Element,storage_interleaved>& other)
            throw(std::bad_cast, std::logic_error);

    virtual void assign(
            const IState<NumDims,Element,storage_interleaved>& other)
            throw(std::bad_cast, std::logic_error);

    virtual void exchange(
            IState<NumDims,Element,storage_interleaved>& other)
            throw(std::bad_cast, std::logic_error);

protected:
    virtual boost::array<std::size_t,NumDims> shapeContainer() const {
        boost::array<std::size_t,NumDims> a;
        std::copy(this->shape(), this->shape() + NumDims, a.begin());
        return a;
    }

private:
    const multi_array_type& operator=( const multi_array_type& ); // Disable
    const InterleavedState& operator=( const InterleavedState& ); // Disable
};

template< std::size_t NumDims, typename Element, typename Allocator >
template< typename ExtentList >
InterleavedState<NumDims,Element,Allocator>::InterleavedState(
        const ExtentList& sizes,
        typename Allocator::size_type min_contiguous_count)
    : IStateBase<NumDims,Element>(),
      IState<NumDims,Element,storage_interleaved>(),
      RawMemory<Element,Allocator>(std::max<typename Allocator::size_type>(
                    suzerain::functional::product(sizes.begin(), sizes.end()),
                    min_contiguous_count)),
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
                       other.shapeContainer(),
                       storage_interleaved::storage_order())
{
    // Data copied by RawMemory's copy constructor
}

template< std::size_t NumDims, typename Element, typename Allocator >
void InterleavedState<NumDims,Element,Allocator>::scale(
        const Element &factor)
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
            const Element &factor,
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

} // namespace suzerain

#endif // __SUZERAIN_STATE_IMPL_HPP
