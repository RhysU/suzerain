//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2011 The PECOS Development Team
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
// state.hpp: Implementations of State classes
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef __SUZERAIN_STATE_HPP
#define __SUZERAIN_STATE_HPP

#include <suzerain/common.hpp>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/functional.hpp>
#include <suzerain/mpl.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/state_fwd.hpp>
#include <suzerain/storage.hpp>

/** @file
 * Implementations of StateBase subclasses.
 */

namespace suzerain
{

template< std::size_t NumDims, typename Element, typename Allocator >
template< typename ExtentList >
NoninterleavedState<NumDims,Element,Allocator>::NoninterleavedState(
        const ExtentList& sizes)
    : ContiguousMemory<Element,Allocator>(
              storage_type::compute_storage(sizes.begin())),
      multi_array_type(
              ContiguousMemory<Element,Allocator>::memory_begin(),
              sizes, storage_type())
{
    // NOP
}

template< std::size_t NumDims, typename Element, typename Allocator >
template< typename ExtentList, typename MinStrideList >
NoninterleavedState<NumDims,Element,Allocator>::NoninterleavedState(
        const ExtentList& sizes,
        const MinStrideList& minstrides)
    : ContiguousMemory<Element,Allocator>(
            storage_type::compute_storage(sizes.begin(), minstrides.begin())),
      multi_array_type(
            ContiguousMemory<Element,Allocator>::memory_begin(),
            sizes, minstrides, storage_type())
{
    // NOP
}

template< std::size_t NumDims, typename Element, typename Allocator >
NoninterleavedState<NumDims,Element,Allocator>::NoninterleavedState(
        const NoninterleavedState &other)
    : ContiguousMemory<Element,Allocator>(other),
      multi_array_type(this->memory_begin(),
                       suzerain::multi_array::shape_array(other),
                       suzerain::multi_array::strides_array(other),
                       storage_type())
{
    // Data copied by ContiguousMemory's copy constructor
    // Strides copied by chosen multi_array_type constructor
}

template< std::size_t NumDims, typename Element, typename Allocator >
void NoninterleavedState<NumDims,Element,Allocator>::scale(
        const Element& factor)
{
    suzerain::blas::scal(
            std::distance(this->memory_begin(),this->memory_end()),
            factor, this->memory_begin(), 1);
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
                &x[ix][x.index_bases()[1]], x.strides()[1],
                &y[iy][y.index_bases()[1]], y.strides()[1]);
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
                    &x[ix][x.index_bases()[1]][kx], x.strides()[1],
                    &y[iy][y.index_bases()[1]][ky], y.strides()[1]);
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
                        &x[ix][x.index_bases()[1]][kx][lx], x.strides()[1],
                        &y[iy][y.index_bases()[1]][ky][ly], y.strides()[1]);
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
                            &x[ix][x.index_bases()[1]][kx][lx][mx],
                            x.strides()[1],
                            &y[iy][y.index_bases()[1]][ky][ly][my],
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
            const NoninterleavedState& other)
{
    if (SUZERAIN_UNLIKELY(this == boost::addressof(other)))
        throw std::logic_error("Unable to handle this->addScaled(...,this)");

    if (SUZERAIN_UNLIKELY(!isIsomorphic(other))) throw std::logic_error(
            std::string("Non-isomorphic other in ") + __PRETTY_FUNCTION__);

    if (SUZERAIN_UNLIKELY(std::equal(other.strides(),
                    other.strides() + NumDims, this->strides()))) {
        // Identical strides between elements: use a single BLAS call
        suzerain::blas::axpy(
                std::distance(this->memory_begin(),this->memory_end()),
                factor, other.memory_begin(), 1, this->memory_begin(), 1);
    } else {
        // Different strides between elements: loop over BLAS calls
        detail::apply(::suzerain::blas::functor::axpy<Element>(factor),
                      const_cast<NoninterleavedState&>(other), *this);
    }
}

template< std::size_t NumDims, typename Element, typename Allocator >
void NoninterleavedState<NumDims,Element,Allocator>::assign(
            const NoninterleavedState& other)
{
    if (SUZERAIN_UNLIKELY(this == boost::addressof(other))) return; // Self?

    if (SUZERAIN_UNLIKELY(!isIsomorphic(other))) throw std::logic_error(
            std::string("Non-isomorphic other in ") + __PRETTY_FUNCTION__);

    if (SUZERAIN_UNLIKELY(std::equal(other.strides(),
                    other.strides() + NumDims, this->strides()))) {
        // Identical strides between elements: use a single BLAS call
        suzerain::blas::copy(
                std::distance(this->memory_begin(),this->memory_end()),
                other.memory_begin(), 1, this->memory_begin(), 1);
    } else {
        // Different strides between elements: loop over BLAS calls
        detail::apply(::suzerain::blas::functor::copy(),
                      const_cast<NoninterleavedState&>(other), *this);
    }
}

template< std::size_t NumDims, typename Element, typename Allocator >
void NoninterleavedState<NumDims,Element,Allocator>::exchange(
            NoninterleavedState& other)
{
    if (SUZERAIN_UNLIKELY(this == boost::addressof(other))) return; // Self?

    if (SUZERAIN_UNLIKELY(!isIsomorphic(other))) throw std::logic_error(
            std::string("Non-isomorphic other in ") + __PRETTY_FUNCTION__);

    if (SUZERAIN_UNLIKELY(std::equal(other.strides(),
                    other.strides() + NumDims, this->strides()))) {
        // Identical strides between elements: use a single BLAS call
        suzerain::blas::swap(
                std::distance(this->memory_begin(),this->memory_end()),
                other.memory_begin(), 1, this->memory_begin(), 1);
    } else {
        // Different strides between elements: loop over BLAS calls
        detail::apply(::suzerain::blas::functor::swap(), other, *this);
    }
}


template< std::size_t NumDims, typename Element, typename Allocator >
template< typename ExtentList >
InterleavedState<NumDims,Element,Allocator>::InterleavedState(
        const ExtentList& sizes,
        size_type min_total_contiguous_count)
    : ContiguousMemory<Element,Allocator>(std::max<size_type>(
                    suzerain::functional::product(sizes.begin(), sizes.end()),
                    min_total_contiguous_count)),
      multi_array_type(ContiguousMemory<Element,Allocator>::memory_begin(),
                       sizes, storage_type())
{
    // NOP
}

template< std::size_t NumDims, typename Element, typename Allocator >
InterleavedState<NumDims,Element,Allocator>::InterleavedState(
        const InterleavedState &other)
    : ContiguousMemory<Element,Allocator>(other),
      multi_array_type(ContiguousMemory<Element,Allocator>::memory_begin(),
                       suzerain::multi_array::shape_array(other),
                       storage_type())
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
void InterleavedState<NumDims,Element,Allocator>::addScaled(
            const Element& factor,
            const InterleavedState& other)
{
    if (SUZERAIN_UNLIKELY(this == boost::addressof(other)))
        throw std::logic_error("Unable to handle this->addScaled(...,this)");

    if (SUZERAIN_UNLIKELY(!isIsomorphic(other))) throw std::logic_error(
            std::string("Non-isomorphic other in ") + __PRETTY_FUNCTION__);

    assert(other.num_elements() == this->num_elements());
    suzerain::blas::axpy(
            this->num_elements(), factor, other.data(), 1, this->data(), 1);
}

template< std::size_t NumDims, typename Element, typename Allocator >
void InterleavedState<NumDims,Element,Allocator>::assign(
            const InterleavedState& other)
{
    if (SUZERAIN_UNLIKELY(this == boost::addressof(other))) return; // Self?

    if (SUZERAIN_UNLIKELY(!isIsomorphic(other))) throw std::logic_error(
            std::string("Non-isomorphic other in ") + __PRETTY_FUNCTION__);

    assert(other.num_elements() == this->num_elements());
    suzerain::blas::copy(
            this->num_elements(), other.data(), 1, this->data(), 1);
}

template< std::size_t NumDims, typename Element, typename Allocator >
void InterleavedState<NumDims,Element,Allocator>::exchange(
            InterleavedState& other)
{
    if (SUZERAIN_UNLIKELY(this == boost::addressof(other))) return; // Self?

    if (SUZERAIN_UNLIKELY(!isIsomorphic(other))) throw std::logic_error(
            std::string("Non-isomorphic other in ") + __PRETTY_FUNCTION__);

    assert(other.num_elements() == this->num_elements());
    suzerain::blas::swap(
            this->num_elements(), other.data(), 1, this->data(), 1);
}

} // namespace suzerain

#endif // __SUZERAIN_STATE_HPP
