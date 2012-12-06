//--------------------------------------------------------------------------
//
// Copyright (C) 2008, 2009, 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
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
// physical_view.hpp: simplifies viewing contiguous_state in physical space
// $Id$

#ifndef SUZERAIN_PHYSICAL_VIEW_HPP
#define SUZERAIN_PHYSICAL_VIEW_HPP

#include <suzerain/common.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/state.hpp>

namespace suzerain {

// FIXME Unaligned is defensive but inefficient and likely unnecessary
namespace { enum {
    physical_view_map_options = Unaligned
}; }

/**
 * A template typedef for how to view multiple \ref contiguous_state fields in
 * physical space, including a convenient method for constructing such an
 * instance.  The optional first template parameter may be specified to provide
 * a number of fields known at compile time.
 *
 * In physical space, the view reshapes the 4D row-major (F, Y, Z, X) with
 * contiguous (Y, Z, X) into a 2D (F, Y*Z*X) layout where we know F a priori.
 * Reducing the dimensionality encourages linear access and eases indexing
 * overhead.
 */
template <int NFields = Dynamic>
struct physical_view;

/**
 * Specialization for the case when <tt>NFields == Dynamic</tt>.
 * This specialization returns a view with a dynamic number of fields.
 */
template <>
struct physical_view<Dynamic> {

    /** Type of the view returned by \ref create. */
    typedef Map<
                Array<real_t, Dynamic, Dynamic, RowMajor>,
                physical_view_map_options,
                OuterStride<Dynamic>
            > type;

    /** Create a view given state and the parallel decomposition. */
    static inline type create(
            const pencil_grid &dgrid,
            contiguous_state<4,complex_t> &state)
    {
        SUZERAIN_ENSURE(static_cast<int>(state.shape()[1]) == dgrid.local_wave_extent.y());
        SUZERAIN_ENSURE(static_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
        SUZERAIN_ENSURE(static_cast<int>(state.shape()[3]) == dgrid.local_wave_extent.z());

        type retval(reinterpret_cast<real_t *>(state.origin()),
                    state.shape()[0],                   // F
                    dgrid.local_physical_extent.prod(), // Y*Z*X
                    OuterStride<>(  state.strides()[0]
                                  * sizeof(complex_t)/sizeof(real_t)));

        return retval;
    }
};

/**
 * General case for when NFields is known at compile time.  This case returns a
 * view with compile-time-fixed number of fields.  It is more efficient and
 * should be preferred when possible.
 */
template <int NFields>
struct physical_view {

#ifdef DOXYGEN_SHOULD_SKIP_THIS
    BOOST_STATIC_ASSERT(NFields >= 0);
#endif

    /** Type of the view returned by \ref create. */
    typedef Map<
                Array<real_t, NFields, Dynamic, RowMajor>,
                physical_view_map_options,
                OuterStride<Dynamic>
            > type;

    /** Create a view given state and the parallel decomposition. */
    static inline type create(
            const pencil_grid &dgrid,
            contiguous_state<4,complex_t> &state)
    {
        SUZERAIN_ENSURE(static_cast<int>(state.shape()[1]) == dgrid.local_wave_extent.y());
        SUZERAIN_ENSURE(static_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
        SUZERAIN_ENSURE(static_cast<int>(state.shape()[3]) == dgrid.local_wave_extent.z());

        if (static_cast<int>(state.shape()[0]) != NFields) {
            std::ostringstream oss;
            oss << "state.shape()[0] = " << state.shape()[0]
                << " mismatch with NFields = " << NFields
                << " in " << __PRETTY_FUNCTION__;
            throw std::logic_error(oss.str());
        }

        type retval(reinterpret_cast<real_t *>(state.origin()),
                    NFields,                            // F
                    dgrid.local_physical_extent.prod(), // Y*Z*X
                    OuterStride<>(  state.strides()[0]
                                  * sizeof(complex_t)/sizeof(real_t)));

        return retval;
    }
};

} // end namespace suzerain

#endif // SUZERAIN_PHYSICAL_VIEW_HPP
