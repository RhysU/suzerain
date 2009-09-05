/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2009 The PECOS Development Team
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
 * classical.hpp: Obtains classical state from conservative variables
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef PECOS_SUZERAIN_CLASSICAL_H
#define PECOS_SUZERAIN_CLASSICAL_H

#include <cmath>
#include <Eigen/Core>

/** @file
 * Provides routines that compute classical state quantities (based on
 * nondimensional \f$p\f$, \f$T\f$, \f$\vec{u}\f$) from conservative state
 * quantities (based on nondimensional \f$\rho\f$, \f$\vec{m}\f$, \f$e\f$)
 * using the ideal gas equations of state nondimensionalized by a reference
 * length, temperature, and density:
 * \f{align*}
 *     p &= \left(\gamma-1\right) \left(
 *       e - \frac{m\cdot{}m}{2\rho}
 *     \right)
 *     \\
 *     T &= \gamma{} \frac{p}{\rho}
 *     \\
 *     \mu &= T^{\beta}
 *     \\
 *     \lambda &= - \frac{2}{3} \mu
 * \f}
 *
 * The following nondimensional variable definitions are used:
 *   - \f$p\f$ or \c p is pressure.
 *   - \f$T\f$ or \c T is temperature.
 *   - \f$\vec{u}\f$ or \c u is the velocity vector.
 *   - \f$\rho\f$ or \c pho is density.
 *   - \f$\vec{m}\f$ or \c m is the momentum vector \f$\rho\vec{u}\f$.
 *   - \f$e\f$ or \c e is the total energy \f$\rho\tilde{e}\f$ where
 *         \f$\tilde{e}\f$ is the sum of the internal and kinetic
 *         energy densities.
 *   - \f$\mu\f$ or \c mu is the first viscosity.
 *   - \f$\lambda\f$ or \c lambda is the second viscosity.
 *   - \f$\gamma\f$ or \c gamma is the constant ratio of specific heats
 *         \f$C_p/C_v\f$.
 */

namespace pecos
{

namespace suzerain
{

namespace cartesian
{

namespace rhome
{

template < typename Scalar >
void p_T_mu_lambda(const Scalar beta,
        const Scalar gamma,
        const Scalar rho,
        const Eigen::Matrix<Scalar,3,1> &m,
        const Scalar e,
        Scalar &p,
        Scalar &T,
        Scalar &mu,
        Scalar &lambda)
{
    const Scalar rho_inverse = 1.0/rho;

    p      = (gamma-1)*(e - (1.0/2.0)*rho_inverse*m.squaredNorm());
    T      = gamma * p * rho_inverse;
    mu     = pow(T, beta);
    lambda = -2.0/3.0*mu;
}

} // namespace rhome

} // namespace cartesian

} // namespace suzerain

} // namespace pecos


#endif // PECOS_SUZERAIN_CLASSICAL_H
