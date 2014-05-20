//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
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


/**
 * Check if all accumulators with an \ref array of boost::accumulators
 * saw an identical number of samples.  Often useful when it is known
 * <em>a priori</em> that this should be true permitting asserting
 * that was <em>a posteriori</em> the case as a way of looking for
 * logic errors.
 *
 * @tparam FPT   Floating point type for the accumulators.
 * @tparam Stats The Accumulators framework details for the accumulators
 *               which \e must specify <tt>boost::accumulators::tag::count</tt>
 *               for this method to compile.
 * @tparam N     The number of accumulators in the array.
 * @param  acc   The array of accumulators to check.
 *
 * @return Zero if all accumulators saw the same number of invocations.
 *         Nonzero indicating the first accumulator possessing a count
 *         not matching that of <code>acc[i]</code>.
 */
template <typename FPT, typename Stats, std::size_t N>
std::size_t inconsistent_accumulation_count(
    const array<boost::accumulators::accumulator_set<FPT, Stats>, N> acc)
{
    using boost::accumulators::count;
    const std::size_t expected = count(acc.front());
    for (std::size_t i = 1; i < N; ++i) {
        const std::size_t observed = count(acc[i]);
        if (expected != observed) {
            return i;
        }
    }
    return 0;
}

// Specialization of inconsistent_accumulation_count for degenerate length-0
template <typename FPT, typename Stats>
std::size_t inconsistent_accumulation_count(
    const array<boost::accumulators::accumulator_set<FPT, Stats>, 0> acc)
{
    return 0;
}

// Specialization of inconsistent_accumulation_count for degenerate length-1
template <typename FPT, typename Stats>
std::size_t inconsistent_accumulation_count(
    const array<boost::accumulators::accumulator_set<FPT, Stats>, 1> acc)
{
    return 0;
}

} // namespace suzerain

#endif // SUZERAIN_UTILITY_HPP
