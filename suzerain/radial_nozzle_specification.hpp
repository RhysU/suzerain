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

#ifndef SUZERAIN_RADIAL_NOZZLE_SPECIFICATION_HPP
#define SUZERAIN_RADIAL_NOZZLE_SPECIFICATION_HPP

/** @file
 * Classes for handling radial nozzle problems.
 */

#include <suzerain/common.hpp>

namespace suzerain {

/**
 * Encapsulates the scenario parameters necessary when invoking \ref
 * suzerain_radial_nozzle_solver().  See \ref radial_nozzle.h for how to
 * interpret these settings and baseflow.tex for a writeup on the general
 * problem class.
 */
class radial_nozzle_specification
{

public:

    /** Construct an instance with the given default values. */
    explicit radial_nozzle_specification(
            double Ma0  = std::numeric_limits<double>::quiet_NaN(),
            double gam0 = std::numeric_limits<double>::quiet_NaN(),
            double rho1 = std::numeric_limits<double>::quiet_NaN(),
            double u1   = std::numeric_limits<double>::quiet_NaN(),
            double p1   = std::numeric_limits<double>::quiet_NaN(),
            double R1   = std::numeric_limits<double>::quiet_NaN());

    double Ma0;   //!< Reference Mach number         \f$\mbox{Ma}_0\f$
    double gam0;  //!< Reference specific heat ratio \f$\gamma_0   \f$
    double rho1;  //!< Inner density                 \f$\rho\left(R_1\right)\f$
    double u1;    //!< Inner radial velocity         \f$u   \left(R_1\right)\f$
    double p1;    //!< Inner pressure                \f$p   \left(R_1\right)\f$
    double R1;    //!< Inner radius of interest      \f$R_1\f$

};

} // namespace suzerain

#endif // SUZERAIN_RADIAL_NOZZLE_SPECIFICATION_HPP
