//--------------------------------------------------------------------------
//
// Copyright (C) 2008-2014 Rhys Ulerich
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

#ifndef SUZERAIN_EXTREMA_HPP
#define SUZERAIN_EXTREMA_HPP

/** @file
 * Compute pointwise extrema of quantities on distributed data.
 */

#include <suzerain/common.hpp>

#include <suzerain/state_fwd.hpp>

namespace suzerain {

// Forward declarations
class bsplineop;
class specification_grid;
class pencil_grid;

/**
 * Holds information on the pointwise extrema within a scalar field as a
 * function of wall-normal $y$ collocation point.
 */
struct field_extrema_xz {

    ArrayXr min;  /**< Pointwise minimum versus collocation point index. */
    ArrayXr max;  /**< Pointwise maximum versus collocation point index. */

};

/**
 * Compute pointwise extrema of all given scalar fields from \c swave
 * destroying \c swave in the process.
 *
 * Incoming fields are represented as Fourier coefficients in $x$ and $z$ \e and
 * B-spline coefficients in $y$ at the B-spline collocation points defining the
 * mass operator in \c cop.
 */
std::vector<field_extrema_xz>
compute_field_extrema_xz(
        contiguous_state<4,complex_t>& swave,
        const specification_grid& grid,
        const pencil_grid& dgrid,
        const bsplineop& cop);

/**
 * Holds information on global extrema within a scalar field.
 */
struct field_extrema_xyz {

    real_t min;  /**< Global pointwise minimum. */
    real_t max;  /**< Global pointwise maximum. */

};

/**
 * Compute pointwise extrema for all given scalar fields from \c swave
 * destroying \c swave in the process.
 *
 * Incoming fields are represented as Fourier coefficients in $x$ and $z$ \e and
 * B-spline coefficients in $y$ at the B-spline collocation points defining the
 * mass operator in \c cop.
 */
std::vector<field_extrema_xyz>
compute_field_extrema_xyz(
        contiguous_state<4,complex_t>& swave,
        const specification_grid& grid,
        const pencil_grid& dgrid,
        const bsplineop& cop);

} // end namespace suzerain

#endif // SUZERAIN_EXTREMA_HPP
