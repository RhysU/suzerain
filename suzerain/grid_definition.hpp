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
// grid_definition.hpp: classes handling grid definitions
// $Id$

#ifndef SUZERAIN_GRID_DEFINITION_HPP
#define SUZERAIN_GRID_DEFINITION_HPP

#include <suzerain/common.hpp>
#include <suzerain/problem.hpp>

// TODO Distinguish between two- versus one-sided stretching in GridDefinition

/** @file
 * Provides classes handling problem grid definitions.
 */

namespace suzerain {

namespace problem {

/**
 * Holds basic three dimensional computational grid details for a distributed,
 * mixed Fourier/B-spline method.  The B-spline representation is used in the
 * wall-normal Y direction.  Logical grid sizes should be specified in terms
 * of physical space coefficient counts.
 */
class GridDefinition : public IDefinition
{
public:
    // See http://eigen.tuxfamily.org/dox/TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /**
     * Construct an instance containing default values intended to be
     * overwritten.  Integer values will be zeros and floating point
     * values will be NaNs.
     */
    GridDefinition();

    /**
     * Construct an instance with the given default values.
     *
     * @param Lx      Physical domain extent in the X direction
     *                which is evaluated using exprparse().
     * @param Nx      Logical grid size in the X direction.
     * @param DAFx    Dealiasing factor in the X direction.
     * @param Ly      Physical domain extent in the Y direction
     *                which is evaluated using exprparse().
     * @param Ny      Logical grid size in the Y direction.
     *                This is the number of B-spline basis functions
     *                (equivalently, wall-normal degrees of freedom)
     *                to use.
     * @param k       Uniform B-spline basis order plus one.
     *                Piecewise cubics correspond to
     *                <tt>default_k == 4</tt>.
     * @param htdelta Hyperbolic tangent stretching parameter
     *                to use when computing breakpoint locations.
     * @param Lz      Physical domain extent in the Z direction
     *                which is evaluated using exprparse().
     * @param Nz      Logical grid size in the Z direction.
     * @param DAFz    Dealiasing factor in the Z direction.
     */
    GridDefinition(const char * Lx,
                   int          Nx,
                   real_t       DAFx,
                   const char * Ly,
                   int          Ny,
                   int          k,
                   real_t       htdelta,
                   const char * Lz,
                   int          Nz,
                   real_t       DAFz);

    /**@{*/

    /** Physical domain extents in the X, Y, and Z directions. */
    Array3r L;

    /** Global logical extents in the X, Y, and Z directions. */
    Array3i N;

    /**
     * Set the logical extents in the X direction.
     *
     * @param value New, nonnegative value to set
     * @return <tt>*this</tt>
     */
    GridDefinition& Nx(int value);

    /**
     * Set the logical extents in the Y direction.
     *
     * @param value New, nonnegative value to set
     * @return <tt>*this</tt>
     */
    GridDefinition& Ny(int value);

    /**
     * Set the logical extents in the Z direction.
     *
     * @param value New, nonnegative value to set
     * @return <tt>*this</tt>
     */
    GridDefinition& Nz(int value);

    /**@}*/

    /**@{*/

    /** Dealiasing factors in the X, Y, and Z directions */
    Array3r DAF;

    /** Global dealiased logical extents in the X, Y, and Z directions. */
    Array3i dN;

    /**
     * Set the dealiasing factor in the X direction.
     *
     * @param factor New, nonnegative factor to set.
     * @return <tt>*this</tt>
     */
    GridDefinition& DAFx(real_t factor);

    /**
     * Set the dealiasing factor in the Z direction.
     *
     * @param factor New, nonnegative factor to set.
     * @return <tt>*this</tt>
     */
    GridDefinition& DAFz(real_t factor);

    /**@}*/

    /**
     * The two dimensional processor grid.
     *
     * In physical space, \f$ P[0] \f$ is the grid extent in the Z direction
     * and \f$ P[1] \f$ is the grid extent in the Y direction.  In wave space,
     * \f$ P[0] \f$ is the grid extent in the X direction and \f$ P[1] \f$  is
     * the grid extent in the Z direction.
     */
    Eigen::Vector2i P;

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

private:
    /** Options initialization common to all constructors */
    void initialize_options(const char * default_Lx,
                            const char * default_Ly,
                            const char * default_Lz);

    /** @copydoc Nx(int) */
    GridDefinition& Nx(const std::string& value);

    /** @copydoc Ny(int) */
    GridDefinition& Ny(const std::string& value);

    /** @copydoc Nz(int) */
    GridDefinition& Nz(const std::string& value);

    /** @copydoc DAFx(real_t) */
    GridDefinition& DAFx(const std::string& value);

    /** @copydoc DAFz(real_t) */
    GridDefinition& DAFz(const std::string& value);
};

} // namespace problem

} // namespace suzerain

#endif // SUZERAIN_GRID_DEFINITION_HPP
