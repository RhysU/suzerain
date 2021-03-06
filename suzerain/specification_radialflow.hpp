//--------------------------------------------------------------------------
//
// Copyright (C) 2013-2014 Rhys Ulerich
// Copyright (C) 2013-2014 The PECOS Development Team
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

#ifndef SUZERAIN_SPECIFICATION_RADIALFLOW_HPP
#define SUZERAIN_SPECIFICATION_RADIALFLOW_HPP

/** @file
 * Classes for handling radial nozzle problems.
 */

#include <suzerain/common.hpp>

namespace suzerain {

/**
 * Encapsulates the non-scenario parameters necessary when invoking \ref
 * suzerain_radialflow_qoi_match(). See \ref radialflow.h for how to interpret
 * these settings and radialflow.tex for a writeup on the general problem class.
 */
class specification_radialflow
{

public:

    /** Construct an instance with the given default values. */
    explicit specification_radialflow(
            double deltae = std::numeric_limits<double>::quiet_NaN(),
            double gamma  = std::numeric_limits<double>::quiet_NaN(),
            double Mae    = std::numeric_limits<double>::quiet_NaN(),
            double pexi   = std::numeric_limits<double>::quiet_NaN());

    /**
     * Edge distance \f$\delta_e\f$ above the \f$x\f$-axis.
     */
    double deltae;

    /**
     * Reference specific heat ratio \f$\gamma\f$.
     *
     * If either 0 or \c NaN, some scenario value should be used.
     */
    double gamma;

    /**
     * Edge Mach number \f$\mbox{Ma}_e\f$ defined per \ref
     * suzerain_radialflow_qoi_Mae().
     *
     * If either 0 or \c NaN, some scenario value should be used.
     */
    double Mae;

    /**
     * Pressure gradient parameter \f$p^\ast_{e,\xi}\f$ defined per \ref
     * suzerain_radialflow_qoi_pexi().
     *
     * If either 0 or \c NaN, the radial flow is considered \ref trivial().
     */
    double pexi;

    /**
     * Is the instance devoid of any interesting data?
     */
    bool trivial() const
    {
        return (boost::math::isnan)(pexi) || (pexi == 0.0);
    }
};

} // namespace suzerain

#endif // SUZERAIN_SPECIFICATION_RADIALFLOW_HPP
