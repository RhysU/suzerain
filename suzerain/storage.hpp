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
// storage.hpp: Marker classes indicating storage ordering
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef __SUZERAIN_STORAGE_HPP
#define __SUZERAIN_STORAGE_HPP

#include <suzerain/common.hpp>
#include <suzerain/iterator.hpp>
#include <suzerain/mpl.hpp>

/** @file
 * Provides marker types indicating state element type and element storage
 * configurations.
 **/

namespace suzerain
{

/**
 * Provides marker types indicating state element type and element storage
 * configurations.
 **/
namespace storage
{

// FIXME Document
template< typename StorageOrderSequence >
class general {
private:
    typedef typename suzerain::mpl::sequence_array<StorageOrderSequence>
            sequence_array_type;

public:
    typedef StorageOrderSequence storage_order_sequence;
    static const std::size_t dimensionality = sequence_array_type::static_size;
    typedef boost::general_storage_order<dimensionality> storage_order_type;

    static storage_order_type storage_order() {
        storage_order_type result(
                sequence_array_type().begin(),
                suzerain::iterator::make_infinite_constant(true));
        return result;
    }

    template< typename RAIterator1,
              typename RAIterator2,
              typename OutputIterator >
    static std::size_t compute_strides(RAIterator1 sizes,
                                       RAIterator2 minstrides,
                                       OutputIterator strides);

    template< typename RAIterator,
              typename OutputIterator >
    static std::size_t compute_strides(RAIterator sizes,
                                       OutputIterator strides);

    template< typename RAIterator1,
              typename RAIterator2 >
    static std::size_t compute_storage(RAIterator1 sizes,
                                       RAIterator2 minstrides);

    template< typename RAIterator >
    static std::size_t compute_storage(RAIterator sizes);
};

template< typename StorageOrderSequence >
template< typename RAIterator1, typename RAIterator2, typename OutputIterator >
std::size_t general<StorageOrderSequence>::compute_strides(
        RAIterator1 sizes,
        RAIterator2 minstrides,
        OutputIterator strides)
{
    typedef typename std::iterator_traits<OutputIterator>::value_type
        output_value_type;
    boost::array<output_value_type, dimensionality> stride_array;

    // Compute stride details using known storage order
    const sequence_array_type seq;
    stride_array[seq[0]] = std::max<output_value_type>(1, minstrides[seq[0]]);
    for (std::size_t i = 1; i < dimensionality; ++i) {
        stride_array[seq[i]] = std::max<output_value_type>(
                sizes[seq[i-1]] * stride_array[seq[i-1]],
                minstrides[seq[i]]);
    }

    // Copy stride information to output parameter
    std::copy(stride_array.begin(), stride_array.end(), strides);

    // Return storage requirements from the slowest stride and size
    return sizes[seq[dimensionality-1]] * stride_array[seq[dimensionality-1]];
}

template< typename StorageOrderSequence >
template< typename RAIterator, typename OutputIterator >
std::size_t general<StorageOrderSequence>::compute_strides(
        RAIterator sizes,
        OutputIterator strides)
{
    boost::array<
        typename std::iterator_traits<OutputIterator>::value_type,
        dimensionality
    > scratch;
    std::fill(scratch.begin(), scratch.end(), 1);
    return compute_strides(sizes, scratch.begin(), strides);
}

template< typename StorageOrderSequence >
template< typename RAIterator1, typename RAIterator2 >
std::size_t general<StorageOrderSequence>::compute_storage(
        RAIterator1 sizes,
        RAIterator2 minstrides)
{
    boost::array<std::size_t,dimensionality> scratch;
    return compute_strides(sizes, minstrides, scratch.begin());
}

template< typename StorageOrderSequence >
template< typename RAIterator >
std::size_t general<StorageOrderSequence>::compute_storage(
        RAIterator sizes)
{
    boost::array<std::size_t,dimensionality> scratch;
    return compute_strides(sizes, scratch.begin());
}

// FIXME Document
template< std::size_t NumDims > class interleaved {
    BOOST_STATIC_ASSERT( NumDims != NumDims ); // Never to be instantiated
};

template<> class interleaved<1>
    : public general< boost::mpl::vector_c<std::size_t,0> > {};
template<> class interleaved<2>
    : public general< boost::mpl::vector_c<std::size_t,0,1> > {};
template<> class interleaved<3>
    : public general< boost::mpl::vector_c<std::size_t,0,1,2> > {};
template<> class interleaved<4>
    : public general< boost::mpl::vector_c<std::size_t,0,1,2,3> > {};
template<> class interleaved<5>
    : public general< boost::mpl::vector_c<std::size_t,0,1,2,3,4> > {};
template<> class interleaved<6>
    : public general< boost::mpl::vector_c<std::size_t,0,1,2,3,4,5> > {};

// FIXME Document
template< std::size_t NumDims > class noninterleaved {
    BOOST_STATIC_ASSERT( NumDims != NumDims ); // Never to be instantiated
};

template<> class noninterleaved<1>
    : public general< boost::mpl::vector_c<std::size_t,0> > {};
template<> class noninterleaved<2>
    : public general< boost::mpl::vector_c<std::size_t,1,0> > {};
template<> class noninterleaved<3>
    : public general< boost::mpl::vector_c<std::size_t,1,2,0> > {};
template<> class noninterleaved<4>
    : public general< boost::mpl::vector_c<std::size_t,1,2,3,0> > {};
template<> class noninterleaved<5>
    : public general< boost::mpl::vector_c<std::size_t,1,2,3,4,0> > {};
template<> class noninterleaved<6>
    : public general< boost::mpl::vector_c<std::size_t,1,2,3,4,5,0> > {};

} // namespace storage

} // namespace suzerain

#endif // __SUZERAIN_STORAGE_HPP
