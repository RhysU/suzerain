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

/**
 * A template typedef for how to view multiple \ref contiguous_state fields in
 * physical space, Including a convenient method for constructing such an
 * instance.  The optional first template parameter may be specified to provide
 * a number of fields known at compile time.
 */
template <int NFields = Dynamic>
struct physical_view {

    BOOST_STATIC_ASSERT(NFields == Dynamic || NFields >= 0);

    /**
     * In physical space, we'll employ a view to reshape the 4D row-major (F,
     * Y, Z, X) with contiguous (Y, Z, X) into a 2D (F, Y*Z*X) layout where we
     * know F a priori.  Reducing the dimensionality encourages linear access
     * and eases indexing overhead.
     */
    typedef Map<
                Array<real_t, NFields, Dynamic, RowMajor>,
                Unaligned, // FIXME Defensive but likely unnecessary
                OuterStride<Dynamic>
            > type;

    /**
     * Create a view instance given state storage and sufficient information
     * about the parallel decomposition.  The default value of \c nfields may
     * only be used when the template parameter \c NFields was not Dynamic.
     */
    static inline type create(
            const pencil_grid &dgrid,
            contiguous_state<4,complex_t> &state,
            const int nfields = NFields)
    {
        if (NFields == Dynamic || NFields == nfields) {
            type retval(reinterpret_cast<real_t *>(state.origin()),
                        nfields,                            // F
                        dgrid.local_physical_extent.prod(), // Y*Z*X
                        OuterStride<>(  state.strides()[0]
                                     * sizeof(complex_t)/sizeof(real_t)));

            return retval;
        }

        throw std::invalid_argument(
                "NFields, nfields mismatch in physical_view::create");
    }

};

} // end namespace suzerain

#endif // SUZERAIN_PHYSICAL_VIEW_HPP
