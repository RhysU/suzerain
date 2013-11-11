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

    /** Construct an instance with all NaN initial values. */
    radial_nozzle_specification();

    /** Construct an instance with the given default values. */
    explicit radial_nozzle_specification(
            double Ma0,
            double gam0,
            double rho1,
            double u1,
            double p1,
            double R1);

    double Ma0;   //!< Reference Mach number         \f$\mbox{Ma}_0\f$
    double gam0;  //!< Reference specific heat ratio \f$\gamma_0   \f$
    double rho1;  //!< Initial density               \f$\rho\left(r_1\right)\f$
    double u1;    //!< Initial radial velocity       \f$u   \left(r_1\right)\f$
    double p1;    //!< Initial pressure              \f$p   \left(r_1\right)\f$
    double R1;    //!< Inner radius of interest      \f$R_1\f$

};

} // namespace suzerain

#endif // SUZERAIN_RADIAL_NOZZLE_SPECIFICATION_HPP
