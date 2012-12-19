//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#ifndef SUZERAIN_GRID_SPECIFICATION_HPP
#define SUZERAIN_GRID_SPECIFICATION_HPP

// TODO Distinguish between two- versus one-sided stretching in grid_specification

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
class grid_specification
{
public:
    // See http://eigen.tuxfamily.org/dox/TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /**
     * Construct an instance containing default values intended to be
     * overwritten.  Integer values will be zeros and floating point
     * values will be NaNs.
     */
    grid_specification();

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
     * @param Lz      Physical domain extent in the Z direction.
     * @param Nz      Logical grid size in the Z direction.
     * @param DAFz    Dealiasing factor in the Z direction.
     */
    grid_specification(real_t Lx,
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
    grid_specification& Nx(int value);

    /**
     * Set the logical extents in the Y direction.
     *
     * @param value New, nonnegative value to set
     * @return <tt>*this</tt>
     */
    grid_specification& Ny(int value);

    /**
     * Set the logical extents in the Z direction.
     *
     * @param value New, nonnegative value to set
     * @return <tt>*this</tt>
     */
    grid_specification& Nz(int value);

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
    grid_specification& DAFx(real_t factor);

    /**
     * Set the dealiasing factor in the Z direction.
     *
     * @param factor New, nonnegative factor to set.
     * @return <tt>*this</tt>
     */
    grid_specification& DAFz(real_t factor);

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
     * breakpoint locations.
     *
     * @see suzerain_htstretch1() and suzerain_htstretch2() for examples
     *      of how this parameter is used.
     */
    real_t htdelta;

    /** @copydoc Nx(int) */
    grid_specification& Nx(const std::string& value);

    /** @copydoc Ny(int) */
    grid_specification& Ny(const std::string& value);

    /** @copydoc Nz(int) */
    grid_specification& Nz(const std::string& value);

    /** @copydoc DAFx(real_t) */
    grid_specification& DAFx(const std::string& factor);

    /** @copydoc DAFz(real_t) */
    grid_specification& DAFz(const std::string& factor);
};

} // namespace suzerain

#endif // SUZERAIN_GRID_SPECIFICATION_HPP
