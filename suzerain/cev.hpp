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

#ifndef SUZERAIN_CEV_HPP
#define SUZERAIN_CEV_HPP

/** @file
 * Provides tabulated data for Crew Exploration Vehicle-related simulations.
 */

namespace suzerain {

/** Provides tabulated data for Crew Exploration Vehicle-related simulations. */
namespace cev {

/**
 * Scenario information from fully-laminar CEV simulations at International
 * Space Station return peak heating conditions.  Simulations by Paul Bauman
 * using the FIN-S hypersonic flow code.  Postprocessing performed by Rhys
 * Ulerich per <a
 * href="https://svn.ices.utexas.edu/repos/pecos/turbulence/heatshield_bl/">
 * https://svn.ices.utexas.edu/repos/pecos/turbulence/heatshield_bl/</a>.
 * Reduced data is tracked within <tt>notebooks/cev_iss_lam.in</tt>.  It has
 * been reordered relative to that data file to permit interpolation.
 *
 * @param dstag   Leeward distance from the stagnation point in meters.
 * @param gammae  Ratio of specific heats at the boundary layer edge.
 * @param gammaw  Ratio of specific heats at the wall.
 * @param Mae     Mach number at the boundary layer edge.
 * @param Pre     Prandtl number at the boundary layer edge.
 * @param Prw     Prandtl number at the wall.
 * @param pexi    The inviscid-friendly pressure parameter
 *                \f$p_{e,x}^\ast = \frac{\delta}{\rho_e u_e^2}
 *                \frac{\partial p}{\partial x}\f$.
 * @param T_ratio The edge-to-wall temperature ratio, \f$T_e / T_w\f$.
 *
 * \return ::SUZERAIN_SUCCESS on success.  On error throws an exception.
 */
int iss_laminar(
    const double dstag,
    double& gammae,
    double& gammaw,
    double& Mae,
    double& Pre,
    double& Prw,
    double& pexi,
    double& T_ratio);

} // namespace cev

} // namespace suzerain

#endif /* SUZERAIN_CEV_HPP */
