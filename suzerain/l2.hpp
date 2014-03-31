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
#include <suzerain/state_fwd.hpp>

namespace suzerain {

// Forward declarations
class bsplineop;
class pencil_grid;
class specification_grid;

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
        const specification_grid& grid,
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
 * represented as Fourier coefficients in \f$x\f$ and \f$z\f$ \e and B-spline
 * coefficients in \f$y\f$ at the B-spline collocation points defining the mass
 * operator in \c cop.  See writeup/L2.tex for full details.
 */
std::vector<field_L2xz>
compute_field_L2xz(
        const contiguous_state<4,complex_t> &state,
        const specification_grid& grid,
        const pencil_grid& dgrid,
        const bsplineop& cop);

/**
 * Compute the \f$L^2_{xz}\f$ norm of all given scalar fields from \c state
 * represented as Fourier coefficients in \f$x\f$ and \f$z\f$ \e but collocation
 * point values in \f$y\f$.  See writeup/L2.tex for full details.
 */
std::vector<field_L2xz>
compute_field_L2xz(
        const contiguous_state<4,complex_t> &state,
        const specification_grid& grid,
        const pencil_grid& dgrid);

/**
 * Compute the <em>rank-local</em> Fourier transform of two-point correlation
 * versus separation in the Hermitian-symmetric \f$x\f$ direction.  See
 * writeup/twopoint.tex for full details.  Obtaining a globally correct answer
 * requires using <code>MPI_Reduce</code> (or similar) with <code>MPI_SUM</code>
 * to sum the resulting buffer across all ranks.  Dealiasing modes are \e not
 * included in computation.
 *
 * @param state[in  ] Scalar fields represented as Fourier coefficients in
 *                    \f$x\f$ and \f$z\f$ \e but point values in \f$y\f$.
 * @param si   [in  ] Scalar index of interest within \c state.
 * @param sj   [in  ] Scalar index of interest within \c state.
 * @param grid [in  ] Domain specification.
 * @param dgrid[in  ] Fourier-based domain specification.
 * @param out  [out ] Output to be stored as a row-major, contiguous
 *                    matrix of size <code>grid.N.y()</code> by
 *                    <code>grid.N.x()</code> indexed by \f$(y_j, k_x)\f$.
 */
void
compute_twopoint_xlocal(
        const contiguous_state<4,complex_t> &state,
        const contiguous_state<4,complex_t>::index &si,
        const contiguous_state<4,complex_t>::index &sj,
        const specification_grid& grid,
        const pencil_grid& dgrid,
        complex_t * const out);

/**
 * Compute the <em>rank-local</em> Fourier transform of two-point correlation
 * versus separation in the \f$z\f$ direction.  See writeup/twopoint.tex for
 * full details.  Obtaining a globally correct answer requires using
 * <code>MPI_Reduce</code> (or similar) with <code>MPI_SUM</code> to sum the
 * resulting buffer across all ranks.  Dealiasing modes are \e not included in
 * computation.
 *
 * @param state[in  ] Scalar fields represented as Fourier coefficients in
 *                    \f$x\f$ and \f$z\f$ \e but point values in \f$y\f$.
 * @param si   [in  ] Scalar index of interest within \c state.
 * @param sj   [in  ] Scalar index of interest within \c state.
 * @param grid [in  ] Domain specification.
 * @param dgrid[in  ] Fourier-based domain specification.
 * @param out  [out ] Output to be stored as a row-major, contiguous
 *                    matrix of size <code>grid.N.y()</code> by
 *                    <code>grid.N.z()</code> indexed by \f$(y_j, k_z)\f$.
 */
void
compute_twopoint_zlocal(
        const contiguous_state<4,complex_t> &state,
        const contiguous_state<4,complex_t>::index &si,
        const contiguous_state<4,complex_t>::index &sj,
        const specification_grid& grid,
        const pencil_grid& dgrid,
        complex_t * const out);

} // end namespace suzerain

#endif // SUZERAIN_L2_HPP
