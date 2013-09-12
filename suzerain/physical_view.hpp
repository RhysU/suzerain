//--------------------------------------------------------------------------
//
// Copyright (C) 2008, 2009, 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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

#ifndef SUZERAIN_PHYSICAL_VIEW_HPP
#define SUZERAIN_PHYSICAL_VIEW_HPP

/** @file
 * Simplifies efficiently accessing contiguous_state in physical space
 */

#include <suzerain/common.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/state.hpp>

namespace suzerain {

/**
 * Constant providing Eigen \c MapOptions for physical_view templates.
 * FIXME Unaligned is defensive but inefficient and likely unnecessary
 */
enum { physical_view_map_options = Unaligned };

/**
 * An \ref ArrayXXr-like class providing a physical space view of wave space
 * \ref contiguous_state given pencil decomposition information.  In physical
 * space, the view reshapes the 4D row-major (NFields, Y, Z, X) with contiguous
 * (Y, Z, X) into a 2D (NFields, Y*Z*X) layout.  Reducing the dimensionality
 * encourages linear access and eases indexing overhead.
 *
 * When NFields is known <em>a priori</em>, it may be supplied as the first
 * template parameter.  Doing so should provide a slight performance gain.
 */
template <int NFields = Dynamic>
struct physical_view;

/**
 * Logic for the general case when <tt>NFields == Dynamic</tt>.
 * The Map superclass has a dynamic number of rows.
 */
template <>
struct physical_view<Dynamic> : Map<Array<real_t, Dynamic, Dynamic, RowMajor>,
                                    physical_view_map_options,
                                    OuterStride<Dynamic> >
{
    /**
     * Create a view given state and the parallel decomposition.
     *
     * @param dgrid Parallel pencil composition to use for determining local
     *              portion of the global state field.
     * @param state Padded state storage used to house data when
     *              viewed from wave space perspective.
     */
    physical_view(const pencil_grid &dgrid,
                  contiguous_state<4,complex_t> &state)
        : Map<Array<real_t, Dynamic, Dynamic, RowMajor>,
              physical_view_map_options,
              OuterStride<Dynamic> >
          (reinterpret_cast<real_t *>(state.origin()),
           state.shape()[0],                   // NFields
           dgrid.local_physical_extent.prod(), // Y*Z*X
           OuterStride<>(  sizeof(complex_t)/sizeof(real_t)
                         * state.strides()[0]))
    {
        SUZERAIN_ENSURE(static_cast<int>(state.shape()[1]) == dgrid.local_wave_extent.y());
        SUZERAIN_ENSURE(static_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
        SUZERAIN_ENSURE(static_cast<int>(state.shape()[3]) == dgrid.local_wave_extent.z());
    }
};

/**
 * Logic for the case when \c NFields is strictly positive.
 * The Map superclass has a constant number of rows.
 */
template <int NFields>
struct physical_view : Map<Array<real_t, Dynamic, Dynamic, RowMajor>,
                           physical_view_map_options,
                           OuterStride<Dynamic> >
{
#ifdef SUZERAIN_PARSED_BY_DOXYGEN
    BOOST_STATIC_ASSERT(NFields >= 0);
#endif

    /**
     * Create a view given state and the parallel decomposition.
     *
     * @param dgrid Parallel pencil composition to use for determining local
     *              portion of the global state field.
     * @param state Padded state storage used to house data when
     *              viewed from wave space perspective.
     *
     * @throw std::logic_error if <tt>state.shape()[0] != NFields</tt>.
     */
    physical_view(const pencil_grid &dgrid,
                  contiguous_state<4,complex_t> &state)
        : Map<Array<real_t, Dynamic, Dynamic, RowMajor>,
              physical_view_map_options,
              OuterStride<Dynamic> >
          (reinterpret_cast<real_t *>(state.origin()),
           NFields,                            // NFields
           dgrid.local_physical_extent.prod(), // Y*Z*X
           OuterStride<>(  sizeof(complex_t)/sizeof(real_t)
                         * state.strides()[0]))
    {
        if (SUZERAIN_UNLIKELY(static_cast<int>(state.shape()[0]) != NFields)) {
            std::ostringstream oss;
            oss << "state.shape()[0] = " << state.shape()[0]
                << " mismatch with NFields = " << NFields
                << " in " << __PRETTY_FUNCTION__;
            throw std::logic_error(oss.str());
        }

        SUZERAIN_ENSURE(static_cast<int>(state.shape()[1]) == dgrid.local_wave_extent.y());
        SUZERAIN_ENSURE(static_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
        SUZERAIN_ENSURE(static_cast<int>(state.shape()[3]) == dgrid.local_wave_extent.z());
    }
};

} // end namespace suzerain

#endif // SUZERAIN_PHYSICAL_VIEW_HPP
