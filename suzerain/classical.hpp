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
 * nondimensional \f$p\f$, \f$T\f$, \f$\vec{u}\f$) from other state quantities.
 */

namespace pecos
{

namespace suzerain
{

/**
 * Provides routines that compute classical state quantities (based on
 * nondimensional \f$p\f$, \f$T\f$, \f$\vec{u}\f$) from other state quantities
 * under the assumption of an orthonormal coordinate system with an identity
 * metric tensor.
 */
namespace orthonormal
{

/**
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
 *   - \f$\beta\f$ or \c beta is the constant coefficient in the power
 *          viscosity law.
 *   - \f$\gamma\f$ or \c gamma is the constant ratio of specific heats
 *         \f$C_p/C_v\f$.
 */
namespace rhome
{

/**
 * Compute \f$p\f$, \f$T\f$, \f$\mu\f$, and \f$\lambda\f$ and their gradients
 * using the equation of state.  The gradients are computed using these
 * expansions:
 * \f{align*}
 *      \vec{\nabla}p &= (\gamma-1)\left[
 *            \vec{\nabla}e
 *          + \frac{1}{2}\rho^{-2}\left(\vec{m}\cdot\vec{m}\right) \vec{\nabla}\rho
 *          - \rho^{-1} \left(\vec{\nabla}\vec{m}\right)^{\mathrm{T}} \vec{m}
 *      \right]
 *      \\
 *      \vec{\nabla}T &= \gamma\rho^{-1}\vec{\nabla}p
 *                     - \gamma\rho^{-2} p \vec{\nabla}\rho
 *      \\
 *      \vec{\nabla}\mu &= \beta{}T^{\beta-1}\vec{\nabla}T
 *      \\
 *      \vec{\nabla}\lambda &= -\frac{2}{3}\vec{\nabla}\mu
 * \f}
 *
 * @param[in]  beta \f$\beta\f$
 * @param[in]  gamma \f$\gamma\f$
 * @param[in]  rho \f$rho\f$
 * @param[in]  grad_rho \f$\vec{\nabla}rho\f$
 * @param[in]  m \f$\vec{m}\f$
 * @param[in]  grad_m \f$\vec{\nabla}\vec{m}\f$
 * @param[in]  e \f$e\f$
 * @param[in]  grad_e \f$\vec{\nabla}e\f$
 * @param[out] p \f$p\f$
 * @param[out] grad_p \f$\vec{\nabla}p\f$
 * @param[out] T \f$T\f$
 * @param[out] grad_T \f$\vec{\nabla}T\f$
 * @param[out] mu \f$\mu\f$
 * @param[out] grad_mu \f$\vec{\nabla}\mu\f$
 * @param[out] lambda \f$\lambda\f$
 * @param[out] grad_lambda \f$\vec{\nabla}\lambda\f$
 */
template<typename Scalar>
void p_T_mu_lambda(
        const Scalar &beta,
        const Scalar &gamma,
        const Scalar &rho,
        const Eigen::Matrix<Scalar,3,1> &grad_rho,
        const Eigen::Matrix<Scalar,3,1> &m,
        const Eigen::Matrix<Scalar,3,3> &grad_m,
        const Scalar &e,
        const Eigen::Matrix<Scalar,3,1> &grad_e,
        Scalar &p,
        Eigen::Matrix<Scalar,3,1> &grad_p,
        Scalar &T,
        Eigen::Matrix<Scalar,3,1> &grad_T,
        Scalar &mu,
        Eigen::Matrix<Scalar,3,1> &grad_mu,
        Scalar &lambda,
        Eigen::Matrix<Scalar,3,1> &grad_lambda)
{
    const Scalar rho_inverse     = 1.0/rho;

    // Compute scalar quantities
    p      = (gamma-1.0)*(e - (1.0/2.0)*rho_inverse*m.squaredNorm());
    T      = gamma * p * rho_inverse;
    mu     = pow(T, beta);
    lambda = -2.0/3.0*mu;

    // Compute vector quantities
    grad_p = (gamma-1.0)*(
                grad_e + rho_inverse*(
                      0.5*rho_inverse*m.squaredNorm()*grad_rho
                    - grad_m.transpose()*m
                )
             );
    grad_T      = gamma*rho_inverse*(grad_p - rho_inverse*p*grad_rho);
    grad_mu     = beta*pow(T,beta-1.0)*grad_T;
    grad_lambda = -2.0/3.0*grad_mu;
}

/**
 * Compute \f$\vec{\nabla}\vec{u}\f$.  Uses the expansion
 * \f[
 *      \vec{\nabla}\vec{u} =
 *          \rho^{-1} \vec{\nabla}\vec{m} - \rho^{-2}
 *          \vec{m}\otimes\vec{\nabla}\rho.
 * \f]
 *
 * @param[in]  rho \f$\rho\f$
 * @param[in]  grad_rho \f$\vec{\nabla}\rho\f$.
 * @param[in]  m \f$\vec{m}\f$
 * @param[in]  grad_m \f$\vec{\nabla}\vec{m}\f$
 *
 * @return The gradient of the velocity based on the provided fields.
 */
template<typename Scalar>
Eigen::Matrix<Scalar,3,3> grad_u(
        const Scalar &rho,
        const Eigen::Matrix<Scalar,3,1> &grad_rho,
        const Eigen::Matrix<Scalar,3,1> &m,
        const Eigen::Matrix<Scalar,3,3> &grad_m)
{
    const Scalar rho_inverse = 1.0/rho;

    return rho_inverse*(grad_m - rho_inverse*m*grad_rho.transpose());
}

/**
 * Compute \f$\vec{\nabla}\cdot\vec{u}\f$.  Uses the expansion
 * \f[
 *      \vec{\nabla}\cdot\vec{u} =
 *            \rho^{-1} \vec{\nabla}\cdot\vec{m}
 *          - \rho^{-2} \vec{\nabla}\rho\cdot\vec{m}.
 * \f]
 *
 * @param[in]  rho \f$\rho\f$
 * @param[in]  grad_rho \f$\vec{\nabla}\rho\f$.
 * @param[in]  m \f$\vec{m}\f$
 * @param[in]  div_m \f$\vec{\nabla}\cdot\vec{m}\f$
 *
 * @return The divergence of the velocity based on the provided fields.
 */
template<typename Scalar>
Scalar div_u(
        const Scalar &rho,
        const Eigen::Matrix<Scalar,3,1> &grad_rho,
        const Eigen::Matrix<Scalar,3,1> &m,
        const Scalar &div_m)
{
    const Scalar rho_inverse = 1.0/rho;

    return rho_inverse*(div_m - rho_inverse*grad_rho.dot(m));
}

/**
 * Compute \f$\vec{\nabla}\vec{\nabla}\cdot{}\vec{u}\f$.  Uses the expansion
 * \f[
 * \vec{\nabla}\vec{\nabla}\cdot{}\vec{u} =
 *      2\rho^{-3}\left(\vec{\nabla}\rho\cdot\vec{m}\right)\vec{\nabla}\rho
 *     - \rho^{-2}\left(\vec{\nabla}\vec{\nabla}\rho\right)\vec{m}
 *     - \rho^{-2}\left(\vec{\nabla}\vec{m}\right)^{\mathrm{T}}\vec{\nabla}\rho
 *     - \rho^{-2}\left(\vec{\nabla}\cdot\vec{m}\right)\vec{\nabla}\rho
 *     + \rho^{-1}\vec{\nabla}\vec{\nabla}\cdot\vec{m}.
 * \f]
 *
 * @param[in]  rho \f$\rho\f$
 * @param[in]  grad_rho \f$\vec{\nabla}\rho\f$.
 * @param[in]  grad_grad_rho \f$\vec{\nabla}\vec{\nabla}\rho\f$.
 * @param[in]  m \f$\vec{m}\f$
 * @param[in]  div_m \f$\vec{\nabla}\cdot\vec{m}\f$
 * @param[in]  grad_m \f$\vec{\nabla}\vec{m}\f$
 * @param[in]  grad_div_m \f$\vec{\nabla}\vec{\nabla}\cdot\vec{m}\f$
 *
 * @return The gradient of the divergence of the velocity based on the provided
 *         fields.
 */
template<typename Scalar>
Eigen::Matrix<Scalar,3,1> grad_div_u(
        const Scalar &rho,
        const Eigen::Matrix<Scalar,3,1> &grad_rho,
        const Eigen::Matrix<Scalar,3,3> &grad_grad_rho,
        const Eigen::Matrix<Scalar,3,1> &m,
        const Scalar &div_m,
        const Eigen::Matrix<Scalar,3,3> &grad_m,
        const Eigen::Matrix<Scalar,3,1> &grad_div_m)
{
    const Scalar rho_inverse = 1.0/rho;

    return rho_inverse * (
                grad_div_m + rho_inverse*(
                      (2*rho_inverse*grad_rho.dot(m) - div_m)*grad_rho
                    - grad_grad_rho*m - grad_m.transpose()*grad_rho
                )
            );
}

/**
 * Compute \f$\vec{\nabla}\cdot{}\vec{\nabla}\vec{u}\f$.  Uses the expansion
 * \f[
 * \vec{\nabla}\cdot{}\vec{\nabla}\vec{u} =
 *      2\rho^{-3}\left(\vec{\nabla}\rho\cdot\vec{\nabla}\rho\right)\vec{m}
 *     -2\rho^{-2}\left(\vec{\nabla}\vec{m}\right)\vec{\nabla}\rho
 *     - \rho^{-2}\left(\vec{\nabla}\cdot\vec{\nabla}\rho\right)\vec{m}
 *     + \rho^{-1}\vec{\nabla}\cdot\vec{\nabla}\vec{m}.
 * \f]
 *
 * @param[in]  rho \f$\rho\f$
 * @param[in]  grad_rho \f$\vec{\nabla}\rho\f$.
 * @param[in]  div_grad_rho \f$\vec{\nabla}\cdot\vec{\nabla}\rho\f$.
 * @param[in]  m \f$\vec{m}\f$
 * @param[in]  grad_m \f$\vec{\nabla}\vec{m}\f$
 * @param[in]  div_grad_m \f$\vec{\nabla}\cdot\vec{\nabla}\vec{m}\f$
 *
 * @return The divergence of the gradient of the velocity based on the provided
 *         fields.
 */
template<typename Scalar>
Eigen::Matrix<Scalar,3,1> div_grad_u(
        const Scalar &rho,
        const Eigen::Matrix<Scalar,3,1> &grad_rho,
        const Scalar &div_grad_rho,
        const Eigen::Matrix<Scalar,3,1> &m,
        const Eigen::Matrix<Scalar,3,3> &grad_m,
        const Eigen::Matrix<Scalar,3,1> &div_grad_m)
{
    const Scalar rho_inverse = 1.0/rho;

    return rho_inverse*(
                div_grad_m + rho_inverse*(
                    (2.0*rho_inverse*grad_rho.squaredNorm() - div_grad_rho)*m
                   - 2.0*grad_m*grad_rho
                )
            );
}

// FIXME tau is symmetric, which should be enforced through Eigen return type
/**
 * Compute \f$\accentset{\leftrightarrow}{\tau} =
 *    \mu\left(
 *        \vec{\nabla}\vec{u}
 *      + \left(\vec{\nabla}\vec{u}\right)^{\mathrm{T}}
 *    \right)
 *  + \lambda \left( \vec{\nabla}\cdot\vec{u} \right)
 *            \accentset{\leftrightarrow}{I}\f$.
 *
 * @param[in]  mu \f$\mu\f$
 * @param[in]  lambda \f$\lambda\f$
 * @param[in]  div_u \f$\vec{\nabla}\cdot\vec{u}\f$
 * @param[in]  grad_u \f$\vec{\nabla}\vec{u}\f$
 *
 * @return The viscous stress tensor based on the provided fields.
 */
template<typename Scalar>
Eigen::Matrix<Scalar,3,3> tau(
        const Scalar &mu,
        const Scalar &lambda,
        const Scalar &div_u,
        const Eigen::Matrix<Scalar,3,3> &grad_u)
{
    return mu*(grad_u + grad_u.transpose())
        + lambda*div_u*Eigen::Matrix<Scalar,3,3>::Identity();
}

} // namespace rhome

} // namespace orthonormal

} // namespace suzerain

} // namespace pecos


#endif // PECOS_SUZERAIN_CLASSICAL_H
