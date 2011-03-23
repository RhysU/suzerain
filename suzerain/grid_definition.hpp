/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Suzerain is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * grid_definition.hpp: classes handling grid definitions
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_GRID_DEFINITION_HPP
#define __SUZERAIN_GRID_DEFINITION_HPP

#include <suzerain/common.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/validation.hpp>

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
     * Construct an instance with the given default values.  Setting default
     * values of zero changes option semantics so the parameters become
     * optional.
     *
     * @param default_Nx      Default logical grid size in the X direction.
     * @param default_DAFx    Default dealiasing factor in the X direction.
     * @param default_Ny      Default logical grid size in the Y direction.
     *                        This is the number of B-spline basis functions
     *                        (equivalently, wall-normal degrees of freedom)
     *                        to use.
     * @param default_k       Default uniform B-spline basis order plus one.
     *                        Piecewise cubics correspond to
     *                        <tt>default_k == 4</tt>.
     * @param default htdelta Default hyperbolic tangent stretching parameter
     *                        to use when computing breakpoint locations.
     * @param default_Nz      Default logical grid size in the Z direction.
     * @param default_DAFz    Default dealiasing factor in the Z direction.
     *
     * @note Since <tt>htdelta == 0</tt> is a potentially useful value,
     *       <tt>default_htdelta == -0.0</tt> allows detecting when
     *       <tt>+0.0</tt> is explicitly provided.
     */
    explicit GridDefinition(int    default_Nx      =  0,
                            double default_DAFx    =  0.0,
                            int    default_Ny      =  0,
                            int    default_k       =  0,
                            double default_htdelta = -0.0,
                            int    default_Nz      =  0  ,
                            double default_DAFz    =  0.0);

    /**@{*/

    /** Global logical extents in the X, Y, and Z directions. */
    const Eigen::Array3i N;

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
    const Eigen::Array3d DAF;

    /** Global dealiased logical extents in the X, Y, and Z directions. */
    const Eigen::Array3i dN;

    /**
     * Set the dealiasing factor in the X direction.
     *
     * @param factor New, nonnegative factor to set.
     * @return <tt>*this</tt>
     */
    GridDefinition& DAFx(double factor);

    /**
     * Set the dealiasing factor in the Z direction.
     *
     * @param factor New, nonnegative factor to set.
     * @return <tt>*this</tt>
     */
    GridDefinition& DAFz(double factor);

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
    double htdelta;
};

} // namespace problem

} // namespace suzerain

#endif // __SUZERAIN_GRID_DEFINITION_HPP
