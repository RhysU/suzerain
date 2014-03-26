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

#ifndef SUZERAIN_SPECIFICATION_GRID_HPP
#define SUZERAIN_SPECIFICATION_GRID_HPP

// TODO Distinguish between two- versus one-sided stretching in specification_grid

/** @file
 * Provides classes handling three-dimensional, dealiased grid specifications.
 */

#include <suzerain/common.hpp>

namespace suzerain {

/**
 * Holds basic three dimensional computational grid details for a distributed,
 * mixed Fourier/B-spline method.  The B-spline representation is used in the
 * wall-normal Y direction.  Logical grid sizes should be specified in terms of
 * physical space coefficient counts.
 */
class specification_grid
{
public:
    // See http://eigen.tuxfamily.org/dox/TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /**
     * Construct an instance containing default values intended to be
     * overwritten.  Integer values will be zeros and floating point
     * values will be NaNs.
     */
    specification_grid();

    /**
     * Construct an instance with the given default values.
     *
     * @param Lx      Physical domain extent in the X direction.
     * @param Nx      Logical grid size in the X direction.
     * @param DAFx    Dealiasing factor in the X direction.
     * @param Ly      Physical domain extent in the Y direction.
     * @param Ny      Logical grid size in the Y direction.
     *                This is the number of B-spline basis functions
     *                (equivalently, wall-normal degrees of freedom)
     *                to use.
     * @param k       Uniform B-spline basis order plus one.
     *                Piecewise cubics correspond to
     *                <tt>default_k == 4</tt>.
     * @param htdelta Hyperbolic tangent stretching parameter
     *                to use when computing breakpoint locations.
     *                Positive values indicate two-sided stretching.
     *                Negative values indicate one-sided stretching.
     * @param Lz      Physical domain extent in the Z direction.
     * @param Nz      Logical grid size in the Z direction.
     * @param DAFz    Dealiasing factor in the Z direction.
     */
    specification_grid(real_t Lx,
                       int    Nx,
                       real_t DAFx,
                       real_t Ly,
                       int    Ny,
                       int    k,
                       real_t htdelta,
                       real_t Lz,
                       int    Nz,
                       real_t DAFz);

    /**@{*/

    /** Physical domain extents in the X, Y, and Z directions. */
    Array3r L;

    /**
     * Global logical extents in the X, Y, and Z directions.
     * Mutate via \ref Nx(), \ref Ny(), and \ref Nz() members.
     */
    const Array3i N;

    /**
     * Set the logical extents in the X direction.
     *
     * @param value New, nonnegative value to set
     * @return <tt>*this</tt>
     */
    specification_grid& Nx(int value);

    /**
     * Set the logical extents in the Y direction.
     *
     * @param value New, nonnegative value to set
     * @return <tt>*this</tt>
     */
    specification_grid& Ny(int value);

    /**
     * Set the logical extents in the Z direction.
     *
     * @param value New, nonnegative value to set
     * @return <tt>*this</tt>
     */
    specification_grid& Nz(int value);

    /**@}*/

    /**@{*/

    /**
     * Dealiasing factors in the X, Y, and Z directions
     * Mutate via DAFx(), DAFy(), and DAFz() members.
     */
    const Array3r DAF;

    /** Global dealiased logical extents in the X, Y, and Z directions. */
    const Array3i dN;

    /**
     * Set the dealiasing factor in the X direction.
     *
     * @param factor New, nonnegative factor to set.
     * @return <tt>*this</tt>
     */
    specification_grid& DAFx(real_t factor);

    /**
     * Set the dealiasing factor in the Z direction.
     *
     * @param factor New, nonnegative factor to set.
     * @return <tt>*this</tt>
     */
    specification_grid& DAFz(real_t factor);

    /**@}*/

    /**
     * The two dimensional processor grid.
     *
     * In physical space, \f$ P[0] \f$ is the grid extent in the Z direction
     * and \f$ P[1] \f$ is the grid extent in the Y direction.  In wave space,
     * \f$ P[0] \f$ is the grid extent in the X direction and \f$ P[1] \f$  is
     * the grid extent in the Z direction.
     */
    Vector2i P;

    /**
     * The B-spline basis order plus one.  For example, piecewise cubics have
     * <tt>k() == 4</tt>.
     */
    int k;

    /**
     * The hyperbolic tangent stretching parameter to use when computing
     * breakpoint locations.  Positive values indicate two-sided stretching.
     * Negative values indicate one-sided stretching.
     *
     * @see suzerain_htstretch1() and suzerain_htstretch2() for examples
     *      of how this parameter is used.
     */
    real_t htdelta;

    /** Is two-sided grid stretching via suzerain_htstretch2() in effect? */
    bool two_sided() const { return htdelta >= 0; }

    /** Is one-sided grid stretching via suzerain_htstretch1() in effect? */
    bool one_sided() const { return !two_sided(); }

    /** @copydoc Nx(int) */
    specification_grid& Nx(const std::string& value);

    /** @copydoc Ny(int) */
    specification_grid& Ny(const std::string& value);

    /** @copydoc Nz(int) */
    specification_grid& Nz(const std::string& value);

    /** @copydoc DAFx(real_t) */
    specification_grid& DAFx(const std::string& factor);

    /** @copydoc DAFz(real_t) */
    specification_grid& DAFz(const std::string& factor);

    /**
     * Return the <tt>i</tt>th \c globally-indexed x collocation point, taking
     * into account the chosen dealiasing factors.
     *
     * On a degenerate grid with only a single point, the collocation point is
     * taken to be at zero.
     */
    real_t x(std::size_t i) const
    {
        return dN.x() > 1 ? ((i * L.x()) / dN.x() - L.x() / 2) : 0;
    }

    /**
     * Return the <tt>k</tt>th \c globally-indexed z collocation point, taking
     * into account the chosen dealiasing factors.
     *
     * On a degenerate grid with only a single point, the collocation point is
     * taken to be at zero.
     */
    real_t z(std::size_t k) const
    {
        return dN.z() > 1 ? ((k * L.z()) / dN.z() - L.z() / 2) : 0;
    }

};

} // namespace suzerain

#endif // SUZERAIN_SPECIFICATION_GRID_HPP
