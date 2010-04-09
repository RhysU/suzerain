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

/**
 * Convert a boost::array from XYZ to YXZ ordering.
 *
 * @param xyz to convert.
 *
 * @return A copy of <tt>xyz</tt> with the elements reordered.
 */
template< typename T >
boost::array<T,3> to_yxz(const boost::array<T,3>& xyz)
{
    boost::array<T,3> retval = { xyz[1], xyz[0], xyz[2] };
    return retval;
}

/**
 * Convert a 3D boost::array from XYZ to YXZ ordering with an
 * additional element prepended.
 *
 * @param prepend to prepend.
 * @param xyz to convert.
 *
 * @return A copy of <tt>xyz</tt> with the elements reordered
 *         and <tt>prepend</tt> prepended.
 */
template< typename T, typename U >
boost::array<T,4> to_yxz(const U& prepend, const boost::array<T,3>& xyz)
{
    boost::array<T,4> retval = { prepend, xyz[1], xyz[0], xyz[2] };
    return retval;
}

/**
 * Convert a boost::array from XYZ to XZY ordering.
 *
 * @param xyz to convert.
 *
 * @return A copy of <tt>xyz</tt> with the elements reordered.
 */
template< typename T >
boost::array<T,3> to_xzy(const boost::array<T,3>& xyz)
{
    boost::array<T,3> retval = { xyz[0], xyz[2], xyz[1] };
    return retval;
}

/**
 * Convert a 3D boost::array from XYZ to XZY ordering with an
 * additional element prepended.
 *
 * @param prepend to prepend.
 * @param xyz to convert.
 *
 * @return A copy of <tt>xyz</tt> with the elements reordered
 *         and <tt>prepend</tt> prepended.
 */
template< typename T, typename U >
boost::array<T,4> to_xzy(const U& prepend, const boost::array<T,3>& xyz)
{
    boost::array<T,4> retval = { prepend, xyz[0], xyz[2], xyz[1] };
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
    boost::array<T,N> retval = { 1 };
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

} // namespace suzerain

#endif // __SUZERAIN_UTILITY_HPP
