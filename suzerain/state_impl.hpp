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
class ContiguousMemory
    : private Allocator // permits empty base class optimization
{
public:
    template< typename I >
    ContiguousMemory(
        I count,
        typename boost::enable_if<boost::is_integral<I> >::type* dummy = 0)
        : Allocator(),
          pbegin_(Allocator::allocate(boost::numeric_cast<
                      typename Allocator::size_type>(count))),
          pend_(pbegin_ + count) {
        SUZERAIN_UNUSED(dummy);
    }

    ContiguousMemory(const ContiguousMemory &other)
        : Allocator(other),
          pbegin_(Allocator::allocate(boost::numeric_cast<
                      typename Allocator::size_type>(
                      std::distance(other.pbegin_,other.pend_)))),
          pend_(pbegin_ + std::distance(other.pbegin_,other.pend_))
    {
        std::memcpy(pbegin_, other.pbegin_,
            std::distance(other.pbegin_, other.pend_) * sizeof(Element));
    }

    ~ContiguousMemory() {
        Allocator::deallocate(pbegin_, std::distance(pbegin_,pend_));
    }

    typename Allocator::pointer const memory_begin() const {
        return pbegin_;
    }

    typename Allocator::pointer const memory_end() const {
        return pend_;
    }

private:
    const ContiguousMemory& operator=( const ContiguousMemory& ); // Disable

    typename Allocator::pointer const pbegin_;
    typename Allocator::pointer const pend_;
};

// Forward declarations
template< std::size_t NumDims, typename Element, typename Allocator >
    class InterleavedState;
template< std::size_t NumDims, typename Element, typename Allocator >
    class NoninterleavedState;

// FIXME Add inheritance from IState<...,interleaved,noninterleaved>
template<
    std::size_t NumDims,
    typename Element,
    typename Allocator = typename suzerain::blas::allocator<Element>::type
>
class InterleavedState
    : public IState<NumDims,Element,suzerain::storage::interleaved<NumDims>,
                                    suzerain::storage::interleaved<NumDims> >,
      public ContiguousMemory<Element,Allocator>,
      public boost::multi_array_ref<Element, NumDims>
{
public:
    typedef typename suzerain::storage::interleaved<NumDims>
            storage_interleaved;
    typedef typename suzerain::storage::noninterleaved<NumDims>
            storage_noninterleaved;

    typedef typename boost::multi_array_ref<Element, NumDims> multi_array_type;
    typedef typename multi_array_type::value_type value_type;
    typedef typename multi_array_type::reference reference;
    typedef typename multi_array_type::const_reference const_reference;
    typedef typename multi_array_type::size_type size_type;
    typedef typename multi_array_type::difference_type difference_type;
    typedef typename multi_array_type::index index;

    template< typename ExtentList >
    explicit InterleavedState(
            const ExtentList& sizes,
            size_type min_total_contiguous_count = 0);

    InterleavedState(const InterleavedState& other);

    virtual ~InterleavedState() {}

    virtual void scale(const Element& factor);

    virtual bool isConformant(
        const IState<NumDims,Element,storage_interleaved,
                                     storage_interleaved>& other) const
        throw(std::bad_cast);

// FIXME Signature and implementation
//  virtual bool isConformant(
//      const IState<NumDims,Element,storage_noninterleaved,
//                                   storage_interleaved>& other) const
//      throw(std::bad_cast);

    virtual void addScaled(
            const Element& factor,
            const IState<NumDims,Element,storage_interleaved,
                                         storage_interleaved>& other)
            throw(std::bad_cast, std::logic_error);

    virtual void assign(
            const IState<NumDims,Element,storage_interleaved,
                                         storage_interleaved>& other)
            throw(std::bad_cast, std::logic_error);

    virtual void exchange(
            IState<NumDims,Element,storage_interleaved,
                                   storage_interleaved>& other)
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
        size_type min_total_contiguous_count)
    : IStateBase<NumDims,Element>(),
      IState<NumDims,Element,storage_interleaved,
                             storage_interleaved>(),
      ContiguousMemory<Element,Allocator>(std::max<size_type>(
                    suzerain::functional::product(sizes.begin(), sizes.end()),
                    min_total_contiguous_count)),
      multi_array_type(ContiguousMemory<Element,Allocator>::memory_begin(),
                       sizes,
                       storage_interleaved::storage_order())
{
    assert(NumDims == std::distance(sizes.begin(),sizes.end()));
    // NOP
}

template< std::size_t NumDims, typename Element, typename Allocator >
InterleavedState<NumDims,Element,Allocator>::InterleavedState(
        const InterleavedState &other)
    : IStateBase<NumDims,Element>(other),
      IState<NumDims,Element,storage_interleaved,
                             storage_interleaved>(other),
      ContiguousMemory<Element,Allocator>(other),
      multi_array_type(ContiguousMemory<Element,Allocator>::memory_begin(),
                       suzerain::multi_array::shape_array(other),
                       storage_interleaved::storage_order())
{
    // Data copied by ContiguousMemory's copy constructor
}

template< std::size_t NumDims, typename Element, typename Allocator >
void InterleavedState<NumDims,Element,Allocator>::scale(
        const Element& factor)
{
    suzerain::blas::scal(this->num_elements(), factor, this->data(), 1);
}

template< std::size_t NumDims, typename Element, typename Allocator >
bool InterleavedState<NumDims,Element,Allocator>::isConformant(
            const IState<NumDims,Element,storage_interleaved,
                                         storage_interleaved>& other) const
throw(std::bad_cast)
{
    const InterleavedState& o = dynamic_cast<const InterleavedState&>(other);

    return std::equal(this->shape(), this->shape() + NumDims, o.shape());
}

template< std::size_t NumDims, typename Element, typename Allocator >
void InterleavedState<NumDims,Element,Allocator>::addScaled(
            const Element& factor,
            const IState<NumDims,Element,storage_interleaved,
                                         storage_interleaved>& other)
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
            const IState<NumDims,Element,storage_interleaved,
                                         storage_interleaved>& other)
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
            IState<NumDims,Element,storage_interleaved,
                                   storage_interleaved>& other)
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
                    suzerain::storage::noninterleaved<NumDims>,
                    suzerain::storage::noninterleaved<NumDims> >,
      public ContiguousMemory<Element,Allocator>,
      public boost::multi_array_ref<Element, NumDims>
{
public:
    typedef typename suzerain::storage::interleaved<NumDims>
            storage_interleaved;
    typedef typename suzerain::storage::noninterleaved<NumDims>
            storage_noninterleaved;

    typedef typename boost::multi_array_ref<Element, NumDims> multi_array_type;
    typedef typename multi_array_type::value_type value_type;
    typedef typename multi_array_type::reference reference;
    typedef typename multi_array_type::const_reference const_reference;
    typedef typename multi_array_type::size_type size_type;
    typedef typename multi_array_type::difference_type difference_type;
    typedef typename multi_array_type::index index;

    template< typename ExtentList >
    explicit NoninterleavedState(
            const ExtentList& sizes);

    template< typename ExtentList, typename MinStrideList >
    explicit NoninterleavedState(
            const ExtentList& sizes,
            const MinStrideList& minstrides);

    NoninterleavedState(const NoninterleavedState& other);

    virtual ~NoninterleavedState() {}

    virtual void scale(const Element& factor);

    virtual bool isConformant(
        const IState<NumDims,Element,storage_noninterleaved,
                                     storage_noninterleaved>& other) const
        throw(std::bad_cast);

    virtual void addScaled(
            const Element& factor,
            const IState<NumDims,Element,storage_noninterleaved,
                                         storage_noninterleaved>& other)
            throw(std::bad_cast, std::logic_error);

    virtual void assign(
            const IState<NumDims,Element,storage_noninterleaved,
                                         storage_noninterleaved>& other)
            throw(std::bad_cast, std::logic_error);

    virtual void exchange(
            IState<NumDims,Element,storage_noninterleaved,
                                   storage_noninterleaved>& other)
            throw(std::bad_cast, std::logic_error);

private:
    // Disable assignment operators
    const multi_array_type& operator=( const multi_array_type& );
    const NoninterleavedState& operator=( const NoninterleavedState& );
};

template< std::size_t NumDims, typename Element, typename Allocator >
template< typename ExtentList >
NoninterleavedState<NumDims,Element,Allocator>::NoninterleavedState(
        const ExtentList& sizes)
    : IStateBase<NumDims,Element>(),
      IState<NumDims,Element,storage_noninterleaved,
                             storage_noninterleaved>(),
      ContiguousMemory<Element,Allocator>(
              storage_noninterleaved::compute_storage(sizes.begin())),
      multi_array_type(ContiguousMemory<Element,Allocator>::memory_begin(),
                       sizes,
                       storage_noninterleaved::storage_order())
{
    assert(NumDims == std::distance(sizes.begin(),sizes.end()));

    // Abuse encapsulation; stride_list_ is protected from our
    // boost::const_multi_array_ref ancestor.
    storage_noninterleaved::compute_strides(
            sizes.begin(), this->stride_list_.begin());
}

template< std::size_t NumDims, typename Element, typename Allocator >
template< typename ExtentList, typename MinStrideList >
NoninterleavedState<NumDims,Element,Allocator>::NoninterleavedState(
        const ExtentList& sizes,
        const MinStrideList& minstrides)
    : IStateBase<NumDims,Element>(),
      IState<NumDims,Element,storage_noninterleaved,
                             storage_noninterleaved>(),
      ContiguousMemory<Element,Allocator>(
              storage_noninterleaved::compute_storage(
                  sizes.begin(), minstrides.begin())),
      multi_array_type(ContiguousMemory<Element,Allocator>::memory_begin(),
                       sizes,
                       storage_noninterleaved::storage_order())
{
    assert(NumDims == std::distance(sizes.begin(),sizes.end()));
    assert(NumDims == std::distance(minstrides.begin(),minstrides.end()));

    // Abuse encapsulation; stride_list_ is protected from our
    // boost::const_multi_array_ref ancestor.
    storage_noninterleaved::compute_strides(
            sizes.begin(), minstrides.begin(), this->stride_list_.begin());
}

template< std::size_t NumDims, typename Element, typename Allocator >
NoninterleavedState<NumDims,Element,Allocator>::NoninterleavedState(
        const NoninterleavedState &other)
    : IStateBase<NumDims,Element>(other),
      IState<NumDims,Element,storage_noninterleaved,
                             storage_noninterleaved>(other),
      ContiguousMemory<Element,Allocator>(other),
      multi_array_type(this->memory_begin(),
                       suzerain::multi_array::shape_array(other),
                       storage_noninterleaved::storage_order())
{
    // Data copied by ContiguousMemory's copy constructor

    // Beware that boost::multi_array_ref has a shallow copy constructor; we
    // must use non-copy constructor and now explicitly copy strides.  Abuse
    // encapsulation; stride_list_ is protected from our
    // boost::const_multi_array_ref ancestor.
    std::copy(other.stride_list_.begin(), other.stride_list_.end(),
              this->stride_list_.begin());
}

template< std::size_t NumDims, typename Element, typename Allocator >
void NoninterleavedState<NumDims,Element,Allocator>::scale(
        const Element& factor)
{
    suzerain::blas::scal(
            std::distance(this->memory_begin(),this->memory_end()),
            factor, this->memory_begin(), 1);
}

template< std::size_t NumDims, typename Element, typename Allocator >
bool NoninterleavedState<NumDims,Element,Allocator>::isConformant(
            const IState<NumDims,Element,storage_noninterleaved,
                                         storage_noninterleaved>& other) const
throw(std::bad_cast)
{
    const NoninterleavedState& o
        = dynamic_cast<const NoninterleavedState&>(other);

    return std::equal(this->shape(), this->shape() + NumDims, o.shape());
}

namespace detail {

template< typename BLASFunctor, typename Element, typename Allocator >
void apply(BLASFunctor functor,
           NoninterleavedState<1,Element,Allocator>& x,
           NoninterleavedState<1,Element,Allocator>& y)
{
    assert(std::equal(x.shape(), x.shape() + 1, y.shape()));

    functor(x.shape()[0],
            x.memory_begin(), x.strides()[0],
            y.memory_begin(), y.strides()[0]);
}

template< typename BLASFunctor, typename Element, typename Allocator >
void apply(BLASFunctor functor,
           NoninterleavedState<2,Element,Allocator>& x,
           NoninterleavedState<2,Element,Allocator>& y)
{
    typedef typename NoninterleavedState<2,Element,Allocator>::index index;
    assert(std::equal(x.shape(), x.shape() + 2, y.shape()));
    using boost::numeric_cast;
    const index iu = numeric_cast<index>(x.index_bases()[0] + x.shape()[0]);

    for (index ix = x.index_bases()[0], iy = y.index_bases()[0];
         ix < iu;
         ++ix, ++iy) {

        functor(x.shape()[1],
                &(x[ix][x.index_bases()[1]]), x.strides()[1],
                &(y[iy][y.index_bases()[1]]), y.strides()[1]);
    }
}

template< typename BLASFunctor, typename Element, typename Allocator >
void apply(BLASFunctor functor,
           NoninterleavedState<3,Element,Allocator>& x,
           NoninterleavedState<3,Element,Allocator>& y)
{
    typedef typename NoninterleavedState<3,Element,Allocator>::index index;
    assert(std::equal(x.shape(), x.shape() + 3, y.shape()));
    using boost::numeric_cast;
    const index iu = numeric_cast<index>(x.index_bases()[0] + x.shape()[0]);
    const index ku = numeric_cast<index>(x.index_bases()[2] + x.shape()[2]);

    for (index ix = x.index_bases()[0], iy = y.index_bases()[0];
        ix < iu;
        ++ix, ++iy) {

        for (index kx = x.index_bases()[2], ky = y.index_bases()[2];
             kx < ku;
             ++kx, ++ky) {

            functor(x.shape()[1],
                    &(x[ix][x.index_bases()[1]][kx]), x.strides()[1],
                    &(y[iy][y.index_bases()[1]][ky]), y.strides()[1]);
        }
    }
}

template< typename BLASFunctor, typename Element, typename Allocator >
void apply(BLASFunctor functor,
           NoninterleavedState<4,Element,Allocator>& x,
           NoninterleavedState<4,Element,Allocator>& y)
{
    typedef typename NoninterleavedState<4,Element,Allocator>::index index;
    assert(std::equal(x.shape(), x.shape() + 4, y.shape()));
    using boost::numeric_cast;
    const index iu = numeric_cast<index>(x.index_bases()[0] + x.shape()[0]);
    const index ku = numeric_cast<index>(x.index_bases()[2] + x.shape()[2]);
    const index lu = numeric_cast<index>(x.index_bases()[3] + x.shape()[3]);

    for (index ix = x.index_bases()[0], iy = y.index_bases()[0];
        ix < iu;
        ++ix, ++iy) {

        for (index lx = x.index_bases()[3], ly = y.index_bases()[3];
            lx < lu;
            ++lx, ++ly) {

            for (index kx = x.index_bases()[2], ky = y.index_bases()[2];
                kx < ku;
                ++kx, ++ky) {

                functor(x.shape()[1],
                        &(x[ix][x.index_bases()[1]][kx][lx]), x.strides()[1],
                        &(y[iy][y.index_bases()[1]][ky][ly]), y.strides()[1]);
            }
        }
    }
}

template< typename BLASFunctor, typename Element, typename Allocator >
void apply(BLASFunctor functor,
           NoninterleavedState<5,Element,Allocator>& x,
           NoninterleavedState<5,Element,Allocator>& y)
{
    typedef typename NoninterleavedState<5,Element,Allocator>::index index;
    assert(std::equal(x.shape(), x.shape() + 5, y.shape()));
    using boost::numeric_cast;
    const index iu = numeric_cast<index>(x.index_bases()[0] + x.shape()[0]);
    const index ku = numeric_cast<index>(x.index_bases()[2] + x.shape()[2]);
    const index lu = numeric_cast<index>(x.index_bases()[3] + x.shape()[3]);
    const index mu = numeric_cast<index>(x.index_bases()[4] + x.shape()[4]);

    for (index ix = x.index_bases()[0], iy = y.index_bases()[0];
        ix < iu;
        ++ix, ++iy) {

        for (index mx = x.index_bases()[4], my = y.index_bases()[4];
            mx < mu;
            ++mx, ++my) {

            for (index lx = x.index_bases()[3], ly = y.index_bases()[3];
                lx < lu;
                ++lx, ++ly) {

                for (index kx = x.index_bases()[2], ky = y.index_bases()[2];
                    kx < ku;
                    ++kx, ++ky) {

                    functor(x.shape()[1],
                            &(x[ix][x.index_bases()[1]][kx][lx][mx]),
                            x.strides()[1],
                            &(y[iy][y.index_bases()[1]][ky][ly][my]),
                            y.strides()[1]);
                }
            }
        }
    }
}

} // namespace detail

template< std::size_t NumDims, typename Element, typename Allocator >
void NoninterleavedState<NumDims,Element,Allocator>::addScaled(
            const Element& factor,
            const IState<NumDims,Element,storage_noninterleaved,
                                         storage_noninterleaved>& other)
throw(std::bad_cast, std::logic_error)
{
    if (SUZERAIN_UNLIKELY(this == boost::addressof(other)))
        throw std::logic_error("Unable to handle this->addScaled(...,this)");

    if (SUZERAIN_UNLIKELY(!isConformant(other))) throw std::logic_error(
            std::string("Nonconformant other in ") + __PRETTY_FUNCTION__);

    const NoninterleavedState& o
        = dynamic_cast<const NoninterleavedState&>(other);

    if (SUZERAIN_UNLIKELY(
        std::equal(o.strides(), o.strides() + NumDims, this->strides()))) {
        // Identical padding between variables: use a single BLAS call
        suzerain::blas::axpy(
                std::distance(this->memory_begin(),this->memory_end()),
                factor, o.memory_begin(), 1, this->memory_begin(), 1);
    } else {
        // Different padding between variables: loop over BLAS calls
        detail::apply(::suzerain::blas::functor::axpy<Element>(factor),
                      const_cast<NoninterleavedState&>(o), *this);
    }
}

template< std::size_t NumDims, typename Element, typename Allocator >
void NoninterleavedState<NumDims,Element,Allocator>::assign(
            const IState<NumDims,Element,storage_noninterleaved,
                                         storage_noninterleaved>& other)
throw(std::bad_cast, std::logic_error)
{
    if (SUZERAIN_UNLIKELY(this == boost::addressof(other))) return; // Self?

    if (SUZERAIN_UNLIKELY(!isConformant(other))) throw std::logic_error(
            std::string("Nonconformant other in ") + __PRETTY_FUNCTION__);

    const NoninterleavedState& o
        = dynamic_cast<const NoninterleavedState&>(other);

    if (SUZERAIN_UNLIKELY(
        std::equal(o.strides(), o.strides() + NumDims, this->strides()))) {
        // Identical padding between variables: use a single BLAS call
        suzerain::blas::copy(
                std::distance(this->memory_begin(),this->memory_end()),
                o.memory_begin(), 1, this->memory_begin(), 1);
    } else {
        // Different padding between variables: loop over BLAS calls
        detail::apply(::suzerain::blas::functor::copy(),
                      const_cast<NoninterleavedState&>(o), *this);
    }
}

template< std::size_t NumDims, typename Element, typename Allocator >
void NoninterleavedState<NumDims,Element,Allocator>::exchange(
            IState<NumDims,Element,storage_noninterleaved,
                                   storage_noninterleaved>& other)
throw(std::bad_cast, std::logic_error)
{
    if (SUZERAIN_UNLIKELY(this == boost::addressof(other))) return; // Self?

    if (SUZERAIN_UNLIKELY(!isConformant(other))) throw std::logic_error(
            std::string("Nonconformant other in ") + __PRETTY_FUNCTION__);

    NoninterleavedState& o = dynamic_cast<NoninterleavedState&>(other);

    if (SUZERAIN_UNLIKELY(
        std::equal(o.strides(), o.strides() + NumDims, this->strides()))) {
        // Identical padding between variables: use a single BLAS call
        suzerain::blas::swap(
                std::distance(this->memory_begin(),this->memory_end()),
                o.memory_begin(), 1, this->memory_begin(), 1);
    } else {
        // Different padding between variables: loop over BLAS calls
        detail::apply(::suzerain::blas::functor::swap(), o, *this);
    }
}

} // namespace suzerain

#endif // __SUZERAIN_STATE_IMPL_HPP
