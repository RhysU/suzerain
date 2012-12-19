//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#ifndef SUZERAIN_UTILITY_HPP
#define SUZERAIN_UTILITY_HPP

/** @file
 * Miscellaneous utilities that seemingly fit nowhere else.
 */

#include <suzerain/common.hpp>

namespace suzerain {

/**
 * Convert a 3 element array from XYZ to YXZ ordering.
 *
 * @param xyz to convert.
 *
 * @return A copy of <tt>xyz</tt> with the elements reordered.
 */
template< typename RandomAccessContainer >
array<
    typename RandomAccessContainer::value_type, 3
> to_yxz(const RandomAccessContainer& xyz)
{
    array<typename RandomAccessContainer::value_type,3> retval
            = {{ xyz[1], xyz[0], xyz[2] }};
    return retval;
}

/**
 * Convert a 3 element array from XYZ to YXZ ordering.
 *
 * @param xyz to convert.
 *
 * @return A copy of <tt>xyz</tt> with the elements reordered.
 */
template< typename EigenDenseType >
array<
    typename EigenDenseType::Scalar, 3
> to_yxz(const EigenDenseType& xyz)
{
    array<typename EigenDenseType::Scalar,3> retval
            = {{ xyz[1], xyz[0], xyz[2] }};
    return retval;
}

/**
 * Convert a 3 element array from XYZ to YXZ ordering with an
 * additional element prepended.
 *
 * @param prepend to prepend.
 * @param xyz to convert.
 *
 * @return A copy of <tt>xyz</tt> with the elements reordered
 *         and <tt>prepend</tt> prepended.
 */
template< typename RandomAccessContainer >
array<
    typename RandomAccessContainer::value_type, 4
> to_yxz(const typename RandomAccessContainer::value_type& prepend,
         const RandomAccessContainer& xyz)
{
    array<typename RandomAccessContainer::value_type,4> retval
            = {{ prepend, xyz[1], xyz[0], xyz[2] }};
    return retval;
}

/**
 * Convert a 3 element array from XYZ to YXZ ordering with an
 * additional element prepended.
 *
 * @param prepend to prepend.
 * @param xyz to convert.
 *
 * @return A copy of <tt>xyz</tt> with the elements reordered
 *         and <tt>prepend</tt> prepended.
 */
template< typename EigenDenseType >
array<
    typename EigenDenseType::Scalar, 4
> to_yxz(const typename EigenDenseType::Scalar& prepend,
         const EigenDenseType& xyz)
{
    array<typename EigenDenseType::Scalar,4> retval
            = {{ prepend, xyz[1], xyz[0], xyz[2] }};
    return retval;
}

/**
 * Convert a 3 element array from XYZ to XZY ordering.
 *
 * @param xyz to convert.
 *
 * @return A copy of <tt>xyz</tt> with the elements reordered.
 */
template< typename RandomAccessContainer >
array<
    typename RandomAccessContainer::value_type, 3
> to_xzy(const RandomAccessContainer& xyz)
{
    array<typename RandomAccessContainer::value_type,3> retval
            = {{ xyz[0], xyz[2], xyz[1] }};
    return retval;
}

/**
 * Convert a 3 element array from XYZ to XZY ordering.
 *
 * @param xyz to convert.
 *
 * @return A copy of <tt>xyz</tt> with the elements reordered.
 */
template< typename EigenDenseType >
array<
    typename EigenDenseType::Scalar, 3
> to_xzy(const EigenDenseType& xyz)
{
    array<typename EigenDenseType::Scalar,3> retval
            = {{ xyz[0], xyz[2], xyz[1] }};
    return retval;
}

/**
 * Convert a 3 element array from XYZ to XZY ordering with an additional
 * element prepended.
 *
 * @param prepend to prepend.
 * @param xyz to convert.
 *
 * @return A copy of <tt>xyz</tt> with the elements reordered
 *         and <tt>prepend</tt> prepended.
 */
template< typename RandomAccessContainer, typename U >
array<typename RandomAccessContainer::value_type, 4> to_xzy(
        const U& prepend,
        const RandomAccessContainer& xyz)
{
    array<typename RandomAccessContainer::value_type,4> retval
        = {{ prepend, xyz[0], xyz[2], xyz[1] }};
    return retval;
}

/**
 * Convert a 3 element array from XYZ to XZY ordering with an additional
 * element prepended.
 *
 * @param prepend to prepend.
 * @param xyz to convert.
 *
 * @return A copy of <tt>xyz</tt> with the elements reordered
 *         and <tt>prepend</tt> prepended.
 */
template< typename EigenDenseType, typename U >
array<typename EigenDenseType::Scalar, 4> to_xzy(
        const U& prepend,
        const EigenDenseType& xyz)
{
    array<typename EigenDenseType::Scalar,4> retval
        = {{ prepend, xyz[0], xyz[2], xyz[1] }};
    return retval;
}

/**
 * Return a copy of an array with an element prepended.
 *
 * @param prepend to prepend.
 * @param a to copy.
 *
 * @return An array with one more entry than \c array.
 */
template<typename T, typename U, std::size_t N>
array<T, N+1> prepend(
        const U& prepend,
        const array<T,N>& a)
{
    array<T,N+1> retval;
    retval[0] = prepend;
    std::copy(a.begin(), a.end(), retval.begin() + 1);
    return retval;
}

/**
 * Compute column-major strides given an array of extents.
 *
 * @param extents to use.
 *
 * @return Column-major strides for <tt>extents</tt> elements
 *         stored contiguously in memory.
 */
template< typename T, std::size_t N >
array<T,N> strides_cm(const array<T,N>& extents)
{
    array<T,N> retval = {{ 1 }};
    for (typename array<T,N>::size_type i = 1; i < N; ++i) {
        retval[i] = retval[i-1] * extents[i-1];
    }
    return retval;
}

/** A specialization to handle degenerate zero element input. */
template< typename T >
array<T,0> strides_cm(const array<T,0>& extents)
{
    array<T,0> retval;
    return retval;
}

/**
 * Compute row-major strides given an array of extents.
 *
 * @param extents to use.
 *
 * @return Row-major strides for <tt>extents</tt> elements
 *         stored contiguously in memory.
 */
template< typename T, std::size_t N >
array<T,N> strides_rm(const array<T,N>& extents)
{
    array<T,N> retval;
    retval[N-1] = 1;
    for (typename array<T,N>::size_type i = N-1; i > 0; --i) {
        retval[i-1] = retval[i] * extents[i];
    }
    return retval;
}

/** A specialization to handle degenerate zero element input. */
template< typename T, typename U >
array<T,0> strides_rm(const array<T,0>& extents)
{
    array<T,0> retval;
    return retval;
}

/**
 * Return the minimum of three arguments.
 *
 * @return <tt>std::min(a,std::min(b,c))</tt>.
 */
template< typename T >
const T& min(const T& a,const T& b, const T& c) {
    return ::std::min(a, ::std::min(b,c));
}

/**
* Check if the argument is nonnegative.
*
* @return True if <tt>t >= 0</tt> and false otherwise.
*/
template< typename T >
typename boost::disable_if<boost::is_unsigned<T>, bool>::type
is_nonnegative(T t) {
    return t >= 0;
}

/**
* Check if the argument is nonnegative.
*
* @return True if <tt>t >= 0</tt> and false otherwise.
*/
template< typename T >
typename boost::enable_if<boost::is_unsigned<T>, bool>::type
is_nonnegative(T t) {
    SUZERAIN_UNUSED(t);
    return true;
}

/**
 * Check if any of the values in the sequence is true.
 *
 * @return True if any value in the sequence is true.  False otherwise.
 */
template<class InputIterator>
bool any(InputIterator first, InputIterator last)
{
    using std::iterator_traits;
    typedef typename iterator_traits<InputIterator>::value_type value_type;
    return std::accumulate(first,
                           last,
                           false,
                           std::logical_or<value_type>());
}

/**
 * Check if all of the values in the sequence are true.
 *
 * @return True if all values in the sequence are true.  False otherwise.
 */
template<class InputIterator>
bool all(InputIterator first, InputIterator last)
{
    using std::iterator_traits;
    typedef typename iterator_traits<InputIterator>::value_type value_type;
    return std::accumulate(first,
                           last,
                           true,
                           std::logical_and<value_type>());
}

} // namespace suzerain

#endif // SUZERAIN_UTILITY_HPP
