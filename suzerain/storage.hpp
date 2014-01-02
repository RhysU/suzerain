//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
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

#ifndef SUZERAIN_STORAGE_HPP
#define SUZERAIN_STORAGE_HPP

/** @file
 * Provides marker types indicating element storage configurations.
 **/

#include <suzerain/common.hpp>
#include <suzerain/iterator.hpp>
#include <suzerain/mpl.hpp>

namespace suzerain {

/**
 * Provides marker types indicating element storage configurations.
 **/
namespace storage {

/**
 * Provides storage and stride computations given a storage ordering.
 *
 * @tparam StorageOrderSequence Storage ordering specified using a type
 * adhering to the Boost MPL's <a
 * href="http://www.boost.org/doc/libs/release/libs/mpl/doc/refmanual/forward-sequence.html">
 * Forward Sequence</a> concept.
 */
template< typename StorageOrderSequence >
class general {
private:
    typedef typename suzerain::mpl::sequence_array<StorageOrderSequence>
            sequence_array_type;

public:


    /** The Boost MPL Forward Sequence used to specify the storage order */
    typedef StorageOrderSequence storage_order_sequence;

    /** The dimensionality of the storage order */
    static const std::size_t dimensionality = sequence_array_type::static_size;

    /**
     * A Boost.MultiArray-friendly type for specifying the storage order.
     * @see <a href="http://www.boost.org/doc/libs/release/libs/multi_array">
     * Boost.MultiArray</a> for more information on the MultiArray concept.
     */
    typedef boost::general_storage_order<dimensionality> storage_order_type;

    /**
     * Retrieve a Boost.MultiArray-friendly instance for specifying the storage
     * order.
     * @see <a
     * href="http://www.boost.org/doc/libs/release/libs/multi_array">
     * Boost.MultiArray</a> for more information on the MultiArray concept.
     */
    static storage_order_type storage_order() {
        storage_order_type result(
                sequence_array_type().begin(),
                suzerain::iterator::make_infinite_constant(true));
        return result;
    }

    /**
     * Compute the dimension-by-dimension \c strides given extents with \c
     * sizes where \c minstrides must be obeyed.  All arguments must be
     * iterators with at least general::dimensionality slots available.
     *
     * @param sizes      The extents to contain in each dimension.
     * @param minstrides The minimum stride in each dimension.
     * @param strides    Output strides meeting the constraints set by
     *                   \c sizes and \c minstrides.
     *
     * @return Number of elements of storage necessary to contain
     *         \c sizes per output \c strides.
     */
    template< typename RandomAccessIterator1,
              typename RandomAccessIterator2,
              typename OutputIterator >
    static std::size_t compute_strides(RandomAccessIterator1 sizes,
                                       RandomAccessIterator2 minstrides,
                                       OutputIterator strides);

    /**
     * Compute the dimension-by-dimension \c strides given extents with \c
     * sizes.  All arguments must be iterators with at least
     * general::dimensionality slots available.
     *
     * @param sizes      The extents to contain in each dimension.
     * @param strides    Output strides matching \c sizes.
     *
     * @return Number of elements of storage necessary to contain \c sizes.
     */
    template< typename RandomAccessIterator,
              typename OutputIterator >
    static std::size_t compute_strides(RandomAccessIterator sizes,
                                       OutputIterator strides);

    /**
     * Compute the amount of storage necessary to contain \c sizes elements
     * where \c minstrides must be obeyed.  All arguments must be iterators
     * with at least general::dimensionality slots available.
     *
     * @param sizes      The extents to contain in each dimension.
     * @param minstrides The minimum stride in each dimension.
     *
     * @return Number of elements of storage necessary to contain
     *         \c sizes per output \c strides.
     */
    template< typename RandomAccessIterator1,
              typename RandomAccessIterator2 >
    static std::size_t compute_storage(RandomAccessIterator1 sizes,
                                       RandomAccessIterator2 minstrides);

    /**
     * Compute the amount of storage necessary to contain \c sizes elements.
     * All arguments must be iterators with at least general::dimensionality
     * slots available.
     *
     * @param sizes      The extents to contain in each dimension.
     *
     * @return Number of elements of storage necessary to contain \c sizes.
     */
    template< typename RandomAccessIterator >
    static std::size_t compute_storage(RandomAccessIterator sizes);
};

template< typename StorageOrderSequence >
template< typename RandomAccessIterator1,
          typename RandomAccessIterator2,
          typename OutputIterator >
std::size_t general<StorageOrderSequence>::compute_strides(
        RandomAccessIterator1 sizes,
        RandomAccessIterator2 minstrides,
        OutputIterator strides)
{
    typedef typename std::iterator_traits<OutputIterator>::value_type
        output_value_type;
    array<output_value_type, dimensionality> stride_array;

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
template< typename RandomAccessIterator, typename OutputIterator >
std::size_t general<StorageOrderSequence>::compute_strides(
        RandomAccessIterator sizes,
        OutputIterator strides)
{
    array<
        typename std::iterator_traits<OutputIterator>::value_type,
        dimensionality
    > scratch;
    std::fill(scratch.begin(), scratch.end(), 1);
    return compute_strides(sizes, scratch.begin(), strides);
}

template< typename StorageOrderSequence >
template< typename RandomAccessIterator1, typename RandomAccessIterator2 >
std::size_t general<StorageOrderSequence>::compute_storage(
        RandomAccessIterator1 sizes,
        RandomAccessIterator2 minstrides)
{
    array<std::size_t,dimensionality> scratch;
    return compute_strides(sizes, minstrides, scratch.begin());
}

template< typename StorageOrderSequence >
template< typename RandomAccessIterator >
std::size_t general<StorageOrderSequence>::compute_storage(
        RandomAccessIterator sizes)
{
    array<std::size_t,dimensionality> scratch;
    return compute_strides(sizes, scratch.begin());
}

/**
 * A marker type for specifying interleaved state storage.  The storage is
 * Fortran's row-major ordering \em except that the first two indices are
 * exchanged.  The second index is fastest, the first index next fastest, and
 * indices three and higher follow regular row-major ordering.  Regions with
 * the first index held constant are therefore "interleaved" with the second
 * index and hence the name.
 */
template< std::size_t NumDims > class interleaved;

/** Interleaved storage specification for one dimension */
template<> class interleaved<1>
    : public general< boost::mpl::vector_c<std::size_t,0> > {};

/** Interleaved storage specification for two dimensions */
template<> class interleaved<2>
    : public general< boost::mpl::vector_c<std::size_t,1,0> > {};

/** Interleaved storage specification for three dimensions */
template<> class interleaved<3>
    : public general< boost::mpl::vector_c<std::size_t,1,0,2> > {};

/** Interleaved storage specification for four dimensions */
template<> class interleaved<4>
    : public general< boost::mpl::vector_c<std::size_t,1,0,2,3> > {};

/** Interleaved storage specification for five dimensions */
template<> class interleaved<5>
    : public general< boost::mpl::vector_c<std::size_t,1,0,2,3,4> > {};

/** Interleaved storage specification for six dimensions */
template<> class interleaved<6>
    : public general< boost::mpl::vector_c<std::size_t,1,0,2,3,4,5> > {};

/**
 * A marker type for specifying non-interleaved state storage.  This is storage
 * for which the first index is slowest and the remaining indices follow
 * Fortran row-major storage ordering.  Regions with the first index held
 * constant are therefore "contiguous" and hence the name.
 */
template< std::size_t NumDims > class contiguous;

/** Contiguous storage specification for one dimension */
template<> class contiguous<1>
    : public general< boost::mpl::vector_c<std::size_t,0> > {};

/** Contiguous storage specification for two dimensions */
template<> class contiguous<2>
    : public general< boost::mpl::vector_c<std::size_t,1,0> > {};

/** Contiguous storage specification for three dimensions */
template<> class contiguous<3>
    : public general< boost::mpl::vector_c<std::size_t,1,2,0> > {};

/** Contiguous storage specification for four dimensions */
template<> class contiguous<4>
    : public general< boost::mpl::vector_c<std::size_t,1,2,3,0> > {};

/** Contiguous storage specification for five dimensions */
template<> class contiguous<5>
    : public general< boost::mpl::vector_c<std::size_t,1,2,3,4,0> > {};

/** Contiguous storage specification for six dimensions */
template<> class contiguous<6>
    : public general< boost::mpl::vector_c<std::size_t,1,2,3,4,5,0> > {};

} // namespace storage

} // namespace suzerain

#endif // SUZERAIN_STORAGE_HPP
