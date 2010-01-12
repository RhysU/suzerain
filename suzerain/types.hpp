/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2009 The PECOS Development Team
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
 * pencil_types.hpp: Class providing basic type traits for P3DFFT pencils
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_TYPES_H
#define __SUZERAIN_TYPES_H

#include <suzerain/common.hpp>

namespace suzerain
{

/**
 * An empty base class providing basic integral types used pervasively in
 * Suzerain's pencil decomposition routines.
 */
class integral_types {
public:

    /**
     * @name Index and index offset types
     * Type choices mimic those found in <tt>boost::multi_array_types</tt>.
     * @see <a href="http://www.boost.org/doc/libs/release/libs/multi_array">
     *      Boost.MultiArray</a> for more information on
     *      <tt>boost::multi_array_types</tt>.
     * @{ */
    /**
     * This is a signed integral type used for indexing into pencils.
     * It is also used to represent strides and index bases.
     */
    typedef boost::multi_array_types::index index;

    /**
     * This is an unsigned integral type.
     * It is primarily used to specify pencil shapes.
     */
    typedef boost::multi_array_types::size_type size_type;

    /**
     * This is a signed integral type used to represent the distance
     * between two iterators.  It is the same type as
     * <tt>std::iterator_traits<iterator>::difference_type</tt>.
     **/
    typedef boost::multi_array_types::difference_type difference_type;
    /**  @} */


    /**
     * Multidimensional array types to hold multiple indices and/or offsets.
     * @name Multidimensional index and index offset types
     * @see <a href="http://www.boost.org/doc/libs/release/libs/array">
     *      Boost.Array</a> for more information on <tt>boost::array</tt>.
     * @{ */
    /**
     * An array of three index values representing X, Y, and Z information,
     * respectively.
     **/
    typedef boost::array<index,3> index_3d;

    /** An array holding two index values. */
    typedef boost::array<index,2> index_2d;

    /**
     * An array of three size_type values representing X, Y, and Z information,
     * respectively.
     **/
    typedef boost::array<size_type,3> size_type_3d;

    /** An array holding two size_type values. */
    typedef boost::array<size_type,2> size_type_2d;
    /**  @} */

protected:
    /** Prevent direct instantiation of this traits-like class */
    integral_types() {}
};

} // namespace suzerain

#endif // __SUZERAIN_TYPES_H
