//--------------------------------------------------------------------------
//
// Copyright (C) 2013 Rhys Ulerich
// Copyright (C) 2013 The PECOS Development Team
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

#ifndef SUZERAIN_ISOTHERMAL_SPECIFICATION_HPP
#define SUZERAIN_ISOTHERMAL_SPECIFICATION_HPP

/** @file
 * Provides \ref isothermal_specification.
 */

#include <suzerain/common.hpp>

namespace suzerain {

/**
 * Holds parameters for specifying simple isothermal boundary conditions on the
 * \f$y=0\f$ and \f$y=L_y\f$ planes.  These conditions include no-slip walls,
 * transpiring walls, and constant freestream boundaries.
 */
class isothermal_specification
{
public:

    /**
     * Construct an instance with all parameters set to NaN and zero-length
     * #lower_cs and #upper_cs.  Clients can use NaN as a not-yet-specified or
     * use-the-default value.
     */
    isothermal_specification();

    /**
     * Specify a no-slip wall with given \f$T\f$ at both boundaries.
     *
     * @param wall_T \f$T\f$ for both walls.
     */
    isothermal_specification(real_t wall_T);

    /**
     * Specify a transpiring wall with given inflow velocity \f$v\f$ and
     * temperature \f$T\f$ at both boundaries.
     *
     * @param wall_T          \f$T\f$ for both walls.
     * @param inflow_velocity \f$v\f$ for both walls with positive values
     *                        oriented blowing into the domain interior.
     */
    isothermal_specification(real_t wall_T,
                             real_t inflow_velocity);

    /**
     * Specify two different temperatures and wall blowing velocities.
     *
     * @param lower_T \f$T\f$ at \f$y=0\f$
     * @param lower_v \f$v\f$ at \f$y=0\f$
     * @param upper_T \f$T\f$ at \f$y=L_y\f$
     * @param upper_v \f$v\f$ at \f$y=L_y\f$
     */
    isothermal_specification(real_t lower_T,
                             real_t lower_v,
                             real_t upper_T,
                             real_t upper_v);

    /**
     * Specify two different temperatures, wall blowing velocities,
     * and associated streamwise velocities.
     *
     * @param lower_T \f$T\f$ at \f$y=0\f$
     * @param lower_u \f$u\f$ at \f$y=0\f$
     * @param lower_v \f$v\f$ at \f$y=0\f$
     * @param upper_T \f$T\f$ at \f$y=L_y\f$
     * @param upper_u \f$u\f$ at \f$y=L_y\f$
     * @param upper_v \f$v\f$ at \f$y=L_y\f$
     */
    isothermal_specification(real_t lower_T,
                             real_t lower_u,
                             real_t lower_v,
                             real_t upper_T,
                             real_t upper_u,
                             real_t upper_v);

    /**
     * Conditions on the \f$y=0\f$ boundary.
     * @{
     */

    /** Temperature \f$T\f$ at \f$y=0\f$. */
    real_t lower_T;

    /** Streamwise velocity \f$u\f$ at \f$y=0\f$. */
    real_t lower_u;

    /** Wall-normal velocity \f$v\f$ at \f$y=0\f$. */
    real_t lower_v;

    /** Spanwise velocity \f$w\f$ at \f$y=0\f$. */
    real_t lower_w;

    /** Species mass fractions \f$c_s\f$ at \f$y=0\f$. */
    std::vector<real_t> lower_cs;

    /**@}*/

    /**
     * Conditions on the \f$y=L_y\f$ boundary.
     * @{
     */

    /** Temperature \f$T\f$ at \f$y=L_y\f$. */
    real_t upper_T;

    /** Streamwise velocity \f$u\f$ at \f$y=L_y\f$. */
    real_t upper_u;

    /** Wall-normal velocity \f$v\f$ at \f$y=L_y\f$. */
    real_t upper_v;

    /** Spanwise velocity \f$w\f$ at \f$y=L_y\f$. */
    real_t upper_w;

    /** Species mass fractions \f$c_s\f$ at \f$y=L_y\f$. */
    std::vector<real_t> upper_cs;

    /**@}*/

};

} // namespace suzerain

#endif // SUZERAIN_ISOTHERMAL_SPECIFICATION_HPP
