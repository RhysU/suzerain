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
 * wall-normal Y direction.
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
     * @param default_Nx   Default logical grid size in the X direction.
     * @param default_DAFx Default dealiasing factor in the X direction.
     * @param default_Ny   Default logical grid size in the Y direction.
     *                     This is the number of B-spline basis functions
     *                     (equivalently, wall-normal degrees of freedom)
     *                     to use.
     * @param default_k    Default uniform B-spline basis order plus one.
     *                     Piecewise cubics correspond to
     *                     <tt>default_k == 4</tt>.
     * @param default_Nz   Default logical grid size in the Z direction.
     * @param default_DAFz Default dealiasing factor in the Z direction.
     */
    explicit GridDefinition(int    default_Nx   = 0,
                            double default_DAFx = 0,
                            int    default_Ny   = 0,
                            int    default_k    = 0,
                            int    default_Nz   = 0,
                            double default_DAFz = 0);

    /** Global logical extents in the X, Y, and Z directions. */
    Eigen::Array3i N;

    /** Global dealiased logical extents in the X, Y, and Z directions.  */
    Eigen::Array3i dN;

    /**
     * The B-spline basis order plus one.  For example, piecewise cubics have
     * <tt>k() == 4</tt>.
     */
    int k;

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
     * Obtain the dealiasing factor in the X direction.
     *
     * @return <tt>dN.x()/N.x()</tt> performed using floating point division.
     */
    double DAFx() const { return static_cast<double>(dN.x()) / N.x(); }

    /**
     * Set the dealiasing factor in the X direction.
     *
     * @param factor New, nonnegative factor to set.
     * @return <tt>*this</tt>
     */
    GridDefinition& DAFx(double factor);

    /**
     * Obtain the dealiasing factor in the Z direction.
     *
     * @return <tt>dN.z()/N.z()</tt> performed using floating point division.
     */
    double DAFz() const { return static_cast<double>(dN.z()) / N.z(); }

    /**
     * Set the dealiasing factor in the Z direction.
     *
     * @param factor New, nonnegative factor to set.
     * @return <tt>*this</tt>
     */
    GridDefinition& DAFz(double factor);
};

} // namespace problem

} // namespace suzerain

#endif // __SUZERAIN_GRID_DEFINITION_HPP
