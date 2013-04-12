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
     * #lower_mass_fractions and #upper_mass_fractions.  Clients can use NaN as
     * a not-yet-specified or use-the-default value.
     */
    isothermal_specification();


    /** Virtual destructor to permit use as a base class */
    virtual ~isothermal_specification();

    /**
     * Conditions on the \f$y=0\f$ boundary.
     * @{
     */

    /** Streamwise velocity \f$u\f$ at \f$y=0\f$. */
    real_t lower_u;

    /** Wall-normal velocity \f$v\f$ at \f$y=0\f$. */
    real_t lower_v;

    /** Spanwise velocity \f$w\f$ at \f$y=0\f$. */
    real_t lower_w;

    /** Temperature \f$T\f$ at \f$y=0\f$. */
    real_t lower_T;

    /** Species mass fractions \f$c_s\f$ at \f$y=0\f$. */
    std::vector<real_t> lower_mass_fractions;

    /**@}*/

    /**
     * Conditions on the \f$y=L_y\f$ boundary.
     * @{
     */

    /** Streamwise velocity \f$u\f$ at \f$y=L_y\f$. */
    real_t upper_u;

    /** Wall-normal velocity \f$v\f$ at \f$y=L_y\f$. */
    real_t upper_v;

    /** Spanwise velocity \f$w\f$ at \f$y=L_y\f$. */
    real_t upper_w;

    /** Temperature \f$T\f$ at \f$y=L_y\f$. */
    real_t upper_T;

    /** Species mass fractions \f$c_s\f$ at \f$y=L_y\f$. */
    std::vector<real_t> upper_mass_fractions;

    /**@}*/

};

} // namespace suzerain

#endif // SUZERAIN_ISOTHERMAL_SPECIFICATION_HPP
