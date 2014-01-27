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

#ifndef SUZERAIN_L2_HPP
#define SUZERAIN_L2_HPP

/** @file
 * Compute \f$L^2\f$ norms and related quantities on distributed data.
 */

#include <suzerain/common.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/specification_grid.hpp>
#include <suzerain/state.hpp>

namespace suzerain {

/**
 * Holds information on the \f$L^2_{xyz}\f$ norm of a scalar field
 * or inner product of two scalar fields.
 *
 * To compute root-mean-square (RMS) values instead of \f$L^2\f$ quantities,
 * divide #mean2, #fluctuating2, and #total2 by \f$L_x L_y L_z\f$.  To obtain
 * RMS values from #mean, #fluctuating, and #total, divide each of them by
 * \f$\sqrt{L_x L_y L_z}\f$.
 */
struct field_L2xyz {

    real_t mean;         /**< \f$L^2\f$ of mean state.          */
    real_t fluctuating;  /**< \f$L^2\f$ of state less the mean. */

    real_t mean2       () const { return mean       *mean;         }
    real_t fluctuating2() const { return fluctuating*fluctuating;  }
    real_t total2      () const { return mean2() + fluctuating2(); }
    real_t total       () const { return std::sqrt(total2());      }

};

/**
 * Compute the \f$L^2_{xyz}\f$ norm of all given scalar fields.
 * See writeup/L2.tex for full details.
 */
std::vector<field_L2xyz>
compute_field_L2xyz(
        const contiguous_state<4,complex_t> &state,
        const grid_specification& grid,
        const pencil_grid& dgrid,
        const bsplineop& gop);

/**
 * Holds information on the \f$L^2_{xz}\f$ norm of a scalar field or inner
 * product of two scalar fields.  Provides both coefficient-wise and vectorized
 * access to quantities of interest.
 *
 * To compute root-mean-square (RMS) values instead of \f$L^2\f$ quantities,
 * divide #mean2, #fluctuating2, and #total2 by \f$L_x L_z\f$.  To obtain RMS
 * values from #mean, #fluctuating, and #total, divide each of them by
 * \f$\sqrt{L_x L_z}\f$.
 */
struct field_L2xz {

    ArrayXr mean;         /**< \f$L^2\f$ of mean state.          */
    ArrayXr fluctuating;  /**< \f$L^2\f$ of state less the mean. */

    real_t mean2       (int i) const { return mean       [i]*mean       [i]; }
    real_t fluctuating2(int i) const { return fluctuating[i]*fluctuating[i]; }
    real_t total2      (int i) const { return mean2(i) + fluctuating2(i);    }
    real_t total       (int i) const { return std::sqrt(total2(i));          }

};

/**
 * Compute the \f$L^2_{xz}\f$ norm of all given scalar fields from \c state
 * represented as Fourier coefficients in $x$ and $z$ \e and B-spline
 * coefficients in $y$ at the B-spline collocation points defining the mass
 * operator in \c cop.  See writeup/L2.tex for full details.
 */
std::vector<field_L2xz>
compute_field_L2xz(
        const contiguous_state<4,complex_t> &state,
        const grid_specification& grid,
        const pencil_grid& dgrid,
        const bsplineop& cop);

/**
 * Compute the \f$L^2_{xz}\f$ norm of all given scalar fields from \c state
 * represented as Fourier coefficients in $x$ and $z$ \e but collocation point
 * values in $y$.  See writeup/L2.tex for full details.
 */
std::vector<field_L2xz>
compute_field_L2xz(
        const contiguous_state<4,complex_t> &state,
        const grid_specification& grid,
        const pencil_grid& dgrid);

} // end namespace suzerain

#endif // SUZERAIN_L2_HPP
