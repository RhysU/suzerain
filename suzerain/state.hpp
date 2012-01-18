//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012 The PECOS Development Team
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
#include <suzerain/shared_range.hpp>
#include <suzerain/state_fwd.hpp>
#include <suzerain/storage.hpp>

/** @file
 * Implementations of StateBase subclasses.
 */

namespace suzerain
{

template< std::size_t NumDims, typename Element >
template< typename ExtentList >
ContiguousState<NumDims,Element>::ContiguousState(
        const ExtentList& sizes)
    : shared_range_type(suzerain::allocate_shared_range(
                typename suzerain::blas::allocator<Element>::type(),
                storage_order_type::compute_storage(sizes.begin()))),
      multi_array_type(shared_range_type::begin(), sizes, storage_order_type())
{
    // NOP
}

template< std::size_t NumDims, typename Element >
template< typename ExtentList >
ContiguousState<NumDims,Element>::ContiguousState(
        const shared_range_type& storage,
        const ExtentList& sizes)
    : shared_range_type(storage),
      multi_array_type(shared_range_type::begin(), sizes, storage_order_type())
{
    const std::ptrdiff_t minimum_size
            = storage_order_type::compute_storage(sizes.begin());
    if (shared_range_type::size() < minimum_size) {
        std::ostringstream oss;
        oss << __PRETTY_FUNCTION__ << " requires storage.size() >= "
            << minimum_size << " but only " << shared_range_type::size()
            << " was provided.";
        throw std::invalid_argument(oss.str());
    }
}

template< std::size_t NumDims, typename Element >
template< typename ExtentList, typename MinStrideList >
ContiguousState<NumDims,Element>::ContiguousState(
        const ExtentList& sizes,
        const MinStrideList& minstrides)
    : shared_range_type(suzerain::allocate_shared_range(
                typename suzerain::blas::allocator<Element>::type(),
                storage_order_type::compute_storage(
                        sizes.begin(), minstrides.begin()))),
      multi_array_type(shared_range_type::begin(),
                       sizes, minstrides, storage_order_type())
{
    // NOP
}

template< std::size_t NumDims, typename Element >
template< typename ExtentList, typename MinStrideList >
ContiguousState<NumDims,Element>::ContiguousState(
        const shared_range_type& storage,
        const ExtentList& sizes,
        const MinStrideList& minstrides)
    : shared_range_type(storage),
      multi_array_type(shared_range_type::begin(),
                       sizes, minstrides, storage_order_type())
{
    const std::ptrdiff_t minimum_size = storage_order_type::compute_storage(
            sizes.begin(), minstrides.begin());
    if (shared_range_type::size() < minimum_size) {
        std::ostringstream oss;
        oss << __PRETTY_FUNCTION__ << " requires storage.size() >= "
            << minimum_size << " but only " << shared_range_type::size()
            << " was provided.";
        throw std::invalid_argument(oss.str());
    }
}

template< std::size_t NumDims, typename Element >
ContiguousState<NumDims,Element>::ContiguousState(
        const ContiguousState &other)
    : shared_range_type(suzerain::clone_shared_range(
                typename suzerain::blas::allocator<Element>::type(),
                other.range())),
      multi_array_type(shared_range_type::begin(),
                       suzerain::multi_array::shape_array(other),
                       suzerain::multi_array::strides_array(other),
                       storage_order_type())
{
    // Data copied by suzerain::clone_shared_range
    // Strides copied by chosen multi_array_type constructor
}

template< std::size_t NumDims, typename Element >
void ContiguousState<NumDims,Element>::scale(
        const Element& factor)
{
    suzerain::blas::scal(
            shared_range_type::size(), factor, shared_range_type::begin(), 1);
}

namespace detail {

template< typename BLASFunctor, typename Element >
void apply(BLASFunctor functor,
           multi_array::ref<Element,1>& x,
           ContiguousState<1,Element>&  y)
{
    assert(std::equal(x.shape(), x.shape() + 1, y.shape()));

    functor(x.shape()[0],
            x.range().begin(), x.strides()[0],
            y.range().begin(), y.strides()[0]);
}

template< typename BLASFunctor, typename Element >
void apply(BLASFunctor functor,
           multi_array::ref<Element,2>& x,
           ContiguousState<2,Element>&  y)
{
    typedef typename ContiguousState<2,Element>::index index;
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

template< typename BLASFunctor, typename Element >
void apply(BLASFunctor functor,
           multi_array::ref<Element,3>& x,
           ContiguousState<3,Element>&  y)
{
    typedef typename ContiguousState<3,Element>::index index;
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

template< typename BLASFunctor, typename Element >
void apply(BLASFunctor functor,
           multi_array::ref<Element,4>& x,
           ContiguousState<4,Element>&  y)
{
    typedef typename ContiguousState<4,Element>::index index;
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

template< typename BLASFunctor, typename Element >
void apply(BLASFunctor functor,
           multi_array::ref<Element,5>& x,
           ContiguousState<5,Element>&  y)
{
    typedef typename ContiguousState<5,Element>::index index;
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

template< typename BLASFunctor, typename Element >
void apply(BLASFunctor functor,
           multi_array::ref<Element,1>& x,
           InterleavedState<1,Element>& y)
{
    assert(std::equal(x.shape(), x.shape() + 1, y.shape()));

    functor(x.shape()[0],
            x.range().begin(), x.strides()[0],
            y.range().begin(), y.strides()[0]);
}

template< typename BLASFunctor, typename Element >
void apply(BLASFunctor functor,
           multi_array::ref<Element,2>& x,
           InterleavedState<2,Element>& y)
{
    typedef typename InterleavedState<2,Element>::index index;
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

template< typename BLASFunctor, typename Element >
void apply(BLASFunctor functor,
           multi_array::ref<Element,3>& x,
           InterleavedState<3,Element>& y)
{
    typedef typename InterleavedState<3,Element>::index index;
    assert(std::equal(x.shape(), x.shape() + 3, y.shape()));
    using boost::numeric_cast;
    const index iu = numeric_cast<index>(x.index_bases()[0] + x.shape()[0]);
    const index ku = numeric_cast<index>(x.index_bases()[2] + x.shape()[2]);

    for (index kx = x.index_bases()[2], ky = y.index_bases()[2];
         kx < ku;
         ++kx, ++ky) {

        for (index ix = x.index_bases()[0], iy = y.index_bases()[0];
             ix < iu;
             ++ix, ++iy) {

            functor(x.shape()[1],
                    &x[ix][x.index_bases()[1]][kx], x.strides()[1],
                    &y[iy][y.index_bases()[1]][ky], y.strides()[1]);
        }
    }
}

template< typename BLASFunctor, typename Element >
void apply(BLASFunctor functor,
           multi_array::ref<Element,4>& x,
           InterleavedState<4,Element>& y)
{
    typedef typename InterleavedState<4,Element>::index index;
    assert(std::equal(x.shape(), x.shape() + 4, y.shape()));
    using boost::numeric_cast;
    const index iu = numeric_cast<index>(x.index_bases()[0] + x.shape()[0]);
    const index ku = numeric_cast<index>(x.index_bases()[2] + x.shape()[2]);
    const index lu = numeric_cast<index>(x.index_bases()[3] + x.shape()[3]);

    for (index lx = x.index_bases()[3], ly = y.index_bases()[3];
         lx < lu;
         ++lx, ++ly) {

        for (index kx = x.index_bases()[2], ky = y.index_bases()[2];
             kx < ku;
             ++kx, ++ky) {

            for (index ix = x.index_bases()[0], iy = y.index_bases()[0];
                 ix < iu;
                 ++ix, ++iy) {

                functor(x.shape()[1],
                        &x[ix][x.index_bases()[1]][kx][lx], x.strides()[1],
                        &y[iy][y.index_bases()[1]][ky][ly], y.strides()[1]);
            }
        }
    }
}

template< typename BLASFunctor, typename Element >
void apply(BLASFunctor functor,
           multi_array::ref<Element,5>& x,
           InterleavedState<5,Element>& y)
{
    typedef typename InterleavedState<5,Element>::index index;
    assert(std::equal(x.shape(), x.shape() + 5, y.shape()));
    using boost::numeric_cast;
    const index iu = numeric_cast<index>(x.index_bases()[0] + x.shape()[0]);
    const index ku = numeric_cast<index>(x.index_bases()[2] + x.shape()[2]);
    const index lu = numeric_cast<index>(x.index_bases()[3] + x.shape()[3]);
    const index mu = numeric_cast<index>(x.index_bases()[4] + x.shape()[4]);

    for (index mx = x.index_bases()[4], my = y.index_bases()[4];
         mx < mu;
         ++mx, ++my) {

        for (index lx = x.index_bases()[3], ly = y.index_bases()[3];
             lx < lu;
             ++lx, ++ly) {

            for (index kx = x.index_bases()[2], ky = y.index_bases()[2];
                 kx < ku;
                 ++kx, ++ky) {

                for (index ix = x.index_bases()[0], iy = y.index_bases()[0];
                     ix < iu;
                     ++ix, ++iy) {

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

template< std::size_t NumDims, typename Element >
void ContiguousState<NumDims,Element>::addScaled(
            const Element& factor,
            const ContiguousState& other)
{
    if (SUZERAIN_UNLIKELY(this == boost::addressof(other)))
        throw std::logic_error("Unable to handle this->addScaled(...,this)");

    if (SUZERAIN_UNLIKELY(!isIsomorphic(other))) throw std::logic_error(
            std::string("Non-isomorphic other in ") + __PRETTY_FUNCTION__);

    if (SUZERAIN_UNLIKELY(std::equal(other.strides(),
                    other.strides() + NumDims, this->strides()))) {
        // Identical strides between elements: use a single BLAS call
        suzerain::blas::axpy(shared_range_type::size(), factor,
                             other.shared_range_type::begin(), 1,
                             shared_range_type::begin(), 1);
    } else {
        // Different strides between elements: loop over BLAS calls
        detail::apply(::suzerain::blas::functor::axpy<Element>(factor),
                      const_cast<ContiguousState&>(other), *this);
    }
}

template< std::size_t NumDims, typename Element >
void ContiguousState<NumDims,Element>::assign(
            const ContiguousState& other)
{
    if (SUZERAIN_UNLIKELY(this == boost::addressof(other))) return; // Self?

    if (SUZERAIN_UNLIKELY(!isIsomorphic(other))) throw std::logic_error(
            std::string("Non-isomorphic other in ") + __PRETTY_FUNCTION__);

    if (SUZERAIN_UNLIKELY(std::equal(other.strides(),
                    other.strides() + NumDims, this->strides()))) {
        // Identical strides between elements: use a single BLAS call
        suzerain::blas::copy(shared_range_type::size(),
                other.shared_range_type::begin(), 1,
                shared_range_type::begin(), 1);
    } else {
        // Different strides between elements: loop over BLAS calls
        detail::apply(::suzerain::blas::functor::copy(),
                      const_cast<ContiguousState&>(other), *this);
    }
}

template< std::size_t NumDims, typename Element >
void ContiguousState<NumDims,Element>::exchange(
            ContiguousState& other)
{
    if (SUZERAIN_UNLIKELY(this == boost::addressof(other))) return; // Self?

    if (SUZERAIN_UNLIKELY(!isIsomorphic(other))) throw std::logic_error(
            std::string("Non-isomorphic other in ") + __PRETTY_FUNCTION__);

    if (SUZERAIN_UNLIKELY(std::equal(other.strides(),
                    other.strides() + NumDims, this->strides()))) {
        // Identical strides between elements: use a single BLAS call
        suzerain::blas::swap(shared_range_type::size(),
                other.shared_range_type::begin(), 1,
                shared_range_type::begin(), 1);
    } else {
        // Different strides between elements: loop over BLAS calls
        detail::apply(::suzerain::blas::functor::swap(), other, *this);
    }
}


template< std::size_t NumDims, typename Element >
template< typename ExtentList >
InterleavedState<NumDims,Element>::InterleavedState(
        const ExtentList& sizes,
        size_type min_total_contiguous_count)
    : shared_range_type(suzerain::allocate_shared_range(
                typename suzerain::blas::allocator<Element>::type(),
                std::max<size_type>(
                    suzerain::functional::product(sizes.begin(), sizes.end()),
                    min_total_contiguous_count))),
      multi_array_type(shared_range_type::begin(), sizes, storage_order_type())
{
    // NOP
}

template< std::size_t NumDims, typename Element >
InterleavedState<NumDims,Element>::InterleavedState(
        const InterleavedState &other)
    : shared_range_type(suzerain::clone_shared_range(
                typename suzerain::blas::allocator<Element>::type(),
                other.range())),
      multi_array_type(shared_range_type::begin(),
                       suzerain::multi_array::shape_array(other),
                       storage_order_type())
{
    // Data copied by suzerain::clone_shared_range
}

template< std::size_t NumDims, typename Element >
void InterleavedState<NumDims,Element>::scale(
        const Element& factor)
{
    // Data guaranteed to be contiguous in first num_elements.
    // Any padding from min_total_contiguous_count is unmodified.
    suzerain::blas::scal(this->num_elements(), factor, this->data(), 1);
}

template< std::size_t NumDims, typename Element >
void InterleavedState<NumDims,Element>::addScaled(
            const Element& factor,
            const InterleavedState& other)
{
    if (SUZERAIN_UNLIKELY(this == boost::addressof(other)))
        throw std::logic_error("Unable to handle this->addScaled(...,this)");

    if (SUZERAIN_UNLIKELY(!isIsomorphic(other))) throw std::logic_error(
            std::string("Non-isomorphic other in ") + __PRETTY_FUNCTION__);

    // Data in this and other guaranteed to be contiguous in num_elements.
    // Any padding from min_total_contiguous_count is unmodified.
    suzerain::blas::axpy(
            this->num_elements(), factor, other.data(), 1, this->data(), 1);
}

template< std::size_t NumDims, typename Element >
void InterleavedState<NumDims,Element>::assign(
            const InterleavedState& other)
{
    if (SUZERAIN_UNLIKELY(this == boost::addressof(other))) return; // Self?

    if (SUZERAIN_UNLIKELY(!isIsomorphic(other))) throw std::logic_error(
            std::string("Non-isomorphic other in ") + __PRETTY_FUNCTION__);

    // Data in this and other guaranteed to be contiguous in num_elements.
    // Any padding from min_total_contiguous_count is unmodified.
    suzerain::blas::copy(
            this->num_elements(), other.data(), 1, this->data(), 1);
}

template< std::size_t NumDims, typename Element >
void InterleavedState<NumDims,Element>::exchange(
            InterleavedState& other)
{
    if (SUZERAIN_UNLIKELY(this == boost::addressof(other))) return; // Self?

    if (SUZERAIN_UNLIKELY(!isIsomorphic(other))) throw std::logic_error(
            std::string("Non-isomorphic other in ") + __PRETTY_FUNCTION__);

    // Data in this and other guaranteed to be contiguous in num_elements.
    // Any padding from min_total_contiguous_count is unmodified.
    suzerain::blas::swap(
            this->num_elements(), other.data(), 1, this->data(), 1);
}

} // namespace suzerain

#endif // __SUZERAIN_STATE_HPP
