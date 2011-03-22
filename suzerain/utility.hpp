/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * utility.hpp: miscellaneous utilities that seemingly fit nowhere else
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_UTILITY_HPP
#define __SUZERAIN_UTILITY_HPP

#include <suzerain/common.hpp>

namespace suzerain {

namespace { // anonymous

// Look up scalar type from a Container, e.g. a boost::array
template< typename T >
struct scalar {
    typedef typename T::value_type type;
};

// Look up scalar type from an Eigen base class
template< typename Derived >
struct scalar< Eigen::DenseBase<Derived> > {
    typedef typename Eigen::internal::traits<Derived>::Scalar type;
};

} // end namespace anonymous

/**
 * Convert a 3 element array from XYZ to YXZ ordering.
 *
 * @param xyz to convert.
 *
 * @return A copy of <tt>xyz</tt> with the elements reordered.
 */
template< typename RandomAccessContainer >
boost::array<
    typename scalar<RandomAccessContainer>::type, 3
> to_yxz(const RandomAccessContainer& xyz)
{
    boost::array<typename scalar<RandomAccessContainer>::type,3> retval
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
template< typename RandomAccessContainer, typename U >
boost::array<
    typename scalar<RandomAccessContainer>::type, 4
> to_yxz(const U& prepend, const RandomAccessContainer& xyz)
{
    boost::array<typename scalar<RandomAccessContainer>::type,4> retval
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
boost::array<
    typename scalar<RandomAccessContainer>::type, 3
> to_xzy(const RandomAccessContainer& xyz)
{
    boost::array<typename scalar<RandomAccessContainer>::type,3> retval
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
boost::array<
    typename scalar<RandomAccessContainer>::type, 4
> to_xzy(const U& prepend, const RandomAccessContainer& xyz)
{
    boost::array<typename scalar<RandomAccessContainer>::type,4> retval
        = {{ prepend, xyz[0], xyz[2], xyz[1] }};
    return retval;
}

/**
 * Convert a 3 element array from containing wave space XYZ extents to
 * containing physical space XYZ extents assuming the X-direction is
 * transformed via a complex-to-real Fourier transform.  The X extent's even or
 * odd nature is preserved.  That is, @{Nx, Ny, Nz}@ becomes @{2*(Nx-1), Ny,
 * Nz}@ when @Nx@ is even or it becomes @{2*(Nx-1)+1, Ny, Nz}@ when @Nx is odd.
 * Preserving oddness is necessary for @Nx == 1@ to return non-trivial extents.
 *
 * @param xyz to transform.
 *
 * @return A copy of <tt>xyz</tt> with the elements transformed.
 */
template< typename RandomAccessContainer >
boost::array<
    typename scalar<RandomAccessContainer>::type, 3
> to_physical_xc2r(const RandomAccessContainer& xyz)
{
    boost::array<typename scalar<RandomAccessContainer>::type,3> retval
        = {{ 2*(xyz[0]-1) + (xyz[0] & 1), xyz[1], xyz[2] }};
    assert((retval[0]/2)+1 == xyz[0]);
    return retval;
}

/**
 * Return a copy of a boost::array with an element prepended.
 *
 * @param prepend to prepend.
 * @param array to copy.
 *
 * @return A boost::array with one more entry than \c array.
 */
template<typename T, typename U, std::size_t N>
boost::array<T, N+1> prepend(const U& prepend, const boost::array<T,N>& array)
{
    boost::array<T,N+1> retval;
    retval[0] = prepend;
    std::copy(array.begin(), array.end(), retval.begin() + 1);
    return retval;
}

/**
 * Compute column-major strides given a boost::array of extents.
 *
 * @param extents to use.
 *
 * @return Column-major strides for <tt>extents</tt> elements
 *         stored contiguously in memory.
 */
template< typename T, std::size_t N >
boost::array<T,N> strides_cm(const boost::array<T,N>& extents)
{
    boost::array<T,N> retval = {{ 1 }};
    for (typename boost::array<T,N>::size_type i = 1; i < N; ++i) {
        retval[i] = retval[i-1] * extents[i-1];
    }
    return retval;
}

/** A specialization to handle degenerate zero element input. */
template< typename T >
boost::array<T,0> strides_cm(const boost::array<T,0>& extents)
{
    boost::array<T,0> retval;
    return retval;
}

/**
 * Compute row-major strides given a boost::array of extents.
 *
 * @param extents to use.
 *
 * @return Row-major strides for <tt>extents</tt> elements
 *         stored contiguously in memory.
 */
template< typename T, std::size_t N >
boost::array<T,N> strides_rm(const boost::array<T,N>& extents)
{
    boost::array<T,N> retval;
    retval[N-1] = 1;
    for (typename boost::array<T,N>::size_type i = N-1; i > 0; --i) {
        retval[i-1] = retval[i] * extents[i];
    }
    return retval;
}

/** A specialization to handle degenerate zero element input. */
template< typename T, typename U >
boost::array<T,0> strides_rm(const boost::array<T,0>& extents)
{
    boost::array<T,0> retval;
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

#endif // __SUZERAIN_UTILITY_HPP
