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
 * orthonormal.hpp: Computes classical quantities in an orthonormal frame
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef PECOS_SUZERAIN_ORTHONORMAL_H
#define PECOS_SUZERAIN_ORTHONORMAL_H

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
 * nondimensional \f$p\f$, \f$T\f$, \f$\vec{u}\f$, \f$mu\f$, \f$\lambda\f$)
 * from other conserved and classical state quantities under the assumption of
 * an orthonormal coordinate system with an identity metric tensor.
 *
 * The following nondimensional variable definitions are used:
 *   - \f$p\f$ or \c p is pressure.
 *   - \f$T\f$ or \c T is temperature.
 *   - \f$\vec{u}\f$ or \c u is the velocity vector.
 *   - \f$e\f$ or \c e is the total energy \f$\rho\tilde{e}\f$ where
 *         \f$\tilde{e}\f$ is the sum of the internal and kinetic
 *         energy densities.
 *   - \f$\mu\f$ or \c mu is the first viscosity.
 *   - \f$\lambda\f$ or \c lambda is the second viscosity.
 *   - \f$\accentset{\leftrightarrow}{\tau}\f$ or \c tau is the
 *          viscous stress tensor.
 */
namespace orthonormal
{

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
template<typename Scalar,
         typename Tensor = Eigen::Matrix<Scalar,3,3> >
Tensor tau(
        const Scalar &mu,
        const Scalar &lambda,
        const Scalar &div_u,
        const Tensor &grad_u)
{
// FIXME tau's symmetry should be used/enforced in the computation
    return mu*(grad_u + grad_u.transpose())
        + lambda*div_u*Tensor::Identity();
}

/**
 * Compute \f$\vec{\nabla}\cdot\accentset{\leftrightarrow}{\tau}\f$.
 * Uses the expansion
 * \f[
 *      \vec{\nabla}\cdot\accentset{\leftrightarrow}{\tau} =
 *        \left[
 *           \vec{\nabla}\vec{u}
 *         + \left(\vec{\nabla}\vec{u}\right)^{\mathrm{T}}
 *        \right] \vec{\nabla}\mu
 *      + \mu \vec{\nabla}\cdot\vec{\nabla}\vec{u}
 *      + \left(\mu+\lambda\right)\vec{\nabla}\vec{\nabla}\cdot\vec{u}
 *      + \left(\vec{\nabla}\cdot\vec{u}\right)\vec{\nabla}\lambda
 * \f]
 *
 * @param[in]  mu \f$\mu\f$
 * @param[in]  grad_mu \f$\vec{\nabla}\mu\f$
 * @param[in]  lambda \f$\lambda\f$
 * @param[in]  grad_lambda \f$\vec{\nabla}\lambda\f$
 * @param[in]  div_u \f$\vec{\nabla}\cdot\vec{u}\f$
 * @param[in]  grad_u \f$\vec{\nabla}\vec{u}\f$
 * @param[in]  div_grad_u \f$\vec{\nabla}\cdot\vec{\nabla}\vec{u}\f$
 * @param[in]  grad_div_u \f$\vec{\nabla}\vec{\nabla}\cdot\vec{u}\f$
 *
 * @return The divergence of the viscous stress tensor based on the provided
 *         fields.
 */
template<typename Scalar,
         typename Vector = Eigen::Matrix<Scalar,3,1>,
         typename Tensor = Eigen::Matrix<Scalar,3,3> >
Vector div_tau(
        const Scalar &mu,
        const Vector &grad_mu,
        const Scalar &lambda,
        const Vector &grad_lambda,
        const Scalar &div_u,
        const Tensor &grad_u,
        const Vector &div_grad_u,
        const Vector &grad_div_u)
{
// FIXME symmetry of grad_u + grad_u^T should be used in the computation
    return (grad_u + grad_u.transpose())*grad_mu
        + mu*div_grad_u
        + (mu+lambda)*grad_div_u
        + div_u*grad_lambda;
}

/**
 * Compute \f$\vec{\nabla}\cdot{}e\vec{u}\f$.
 * Uses the expansion
 * \f[
 *      \vec{\nabla}\cdot{}e\vec{u} =
 *          e\vec{\nabla}\cdot\vec{u} + \vec{u}\cdot\vec{\nabla}e
 * \f]
 *
 * @param[in] e \f$e\f$
 * @param[in] grad_e \f$\vec{\nabla}e\f$
 * @param[in] u \f$\vec{u}\f$
 * @param[in] div_u \f$\vec{\nabla}\cdot\vec{u}\f$
 *
 * @return The divergence of the product of total energy and velocity.
 */
template<typename Scalar,
         typename Vector = Eigen::Matrix<Scalar,3,1> >
Scalar div_e_u(
        const Scalar &e,
        const Vector &grad_e,
        const Vector &u,
        const Scalar &div_u)
{
    return e*div_u + u.dot(grad_e);
}

/**
 * Compute \f$\vec{\nabla}\cdot{}\mu\vec{\nabla}T\f$.
 * Uses the expansion
 * \f[
 *      \vec{\nabla}\cdot{}\mu\vec{\nabla}T =
 *          \vec{\nabla}\mu\cdot\vec{\nabla}T
 *          + \mu \vec{\nabla}\cdot\vec{\nabla}T
 * \f]
 *
 * @param[in] grad_T \f$\vec{\nabla}T\f$
 * @param[in] div_grad_T \f$\vec{\nabla}\cdot\vec{\nabla}T\f$
 * @param[in] mu \f$\mu\f$
 * @param[in] grad_mu \f$\vec{\nabla}\mu\f$
 *
 * @return The divergence of the product of viscosity and
 *      the temperature gradient.
 */
template<typename Scalar,
         typename Vector = Eigen::Matrix<Scalar,3,1> >
Scalar div_mu_grad_T(
        const Vector &grad_T,
        const Scalar &div_grad_T,
        const Scalar &mu,
        const Vector &grad_mu)
{
    return grad_mu.dot(grad_T) + mu*div_grad_T;
}

/**
 * Compute \f$\vec{\nabla}\cdot{}p\vec{u}\f$.
 * Uses the expansion
 * \f[
 *      \vec{\nabla}\cdot{}p\vec{u} =
 *            p\vec{\nabla}\cdot\vec{u}
 *          + \vec{u}\cdot\vec{\nabla}p
 * \f]
 *
 * @param[in] p \f$p\f$
 * @param[in] grad_p \f$\vec{\nabla}p\f$
 * @param[in] u \f$\vec{u}\f$
 * @param[in] div_u \f$\vec{\nabla}\cdot\vec{u}\f$
 *
 * @return The divergence of the product of pressure and velocity.
 */
template<typename Scalar,
         typename Vector = Eigen::Matrix<Scalar,3,1> >
Scalar div_p_u(
        const Scalar &p,
        const Vector &grad_p,
        const Vector &u,
        const Scalar &div_u)
{
    return p*div_u + u.dot(grad_p);
}

/**
 * Compute \f$\vec{\nabla}\cdot\accentset{\leftrightarrow}{\tau}\vec{u}\f$.
 * Uses the expansion
 * \f[
 *      \vec{\nabla}\cdot\accentset{\leftrightarrow}{\tau}\vec{u} =
 *      \vec{u}\cdot\left(
 *          \vec{\nabla}\cdot\accentset{\leftrightarrow}{\tau}
 *      \right)
 *      +
 *      \mathrm{trace}\!\left(
 *          \accentset{\leftrightarrow}{\tau}\,\vec{\nabla}\vec{u}
 *      \right)
 * \f]
 *
 * @param[in] u \f$\vec{u}\f$
 * @param[in] grad_u \f$\vec{\nabla}\vec{u}\f$
 * @param[in] tau \f$\accentset{\leftrightarrow}{\tau}\f$
 * @param[in] div_tau \f$\vec{\nabla}\cdot\accentset{\leftrightarrow}{\tau}\f$
 *
 * @return The divergence of the viscous stress tensor applied to the
 *      velocity.
 */
template<typename Scalar,
         typename Vector = Eigen::Matrix<Scalar,3,1>,
         typename Tensor = Eigen::Matrix<Scalar,3,3> >
Scalar div_tau_u(
        const Vector &u,
        const Tensor &grad_u,
        const Tensor &tau,
        const Vector &div_tau)
{
    return u.dot(div_tau) + (tau*grad_u).lazy().trace();
}

/**
 * Provides routines that compute classical state quantities (based on
 * nondimensional \f$p\f$, \f$T\f$, \f$\vec{u}\f$, \f$\mu\f$, \f$\lambda\f$)
 * from conservative state quantities (based on nondimensional \f$\rho\f$,
 * \f$\vec{m}\f$, \f$e\f$) using the ideal gas equations of state
 * nondimensionalized by a reference length, temperature, and density:
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
 *   - \f$\rho\f$ or \c rho is density.
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
 *          + \frac{1}{2}\rho^{-2}\left(\vec{m}\cdot\vec{m}\right)
 *            \vec{\nabla}\rho
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
template<typename Scalar,
         typename Vector = Eigen::Matrix<Scalar,3,1>,
         typename Tensor = Eigen::Matrix<Scalar,3,3> >
void p_T_mu_lambda(
        const Scalar &beta,
        const Scalar &gamma,
        const Scalar &rho,
        const Vector &grad_rho,
        const Vector &m,
        const Tensor &grad_m,
        const Scalar &e,
        const Vector &grad_e,
        Scalar &p,
        Vector &grad_p,
        Scalar &T,
        Vector &grad_T,
        Scalar &mu,
        Vector &grad_mu,
        Scalar &lambda,
        Vector &grad_lambda)
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
 * Compute \f$\vec{u}=\rho^{-1}\vec{m}\f$
 *
 * @param rho \f$\rho\f$
 * @param m \f$\vec{m}\f$
 *
 * @return The velocity vector.
 */
template<typename Scalar,
         typename Vector = Eigen::Matrix<Scalar,3,1> >
Vector u(const Scalar &rho,
         const Vector &m)
{
    const Scalar rho_inverse = 1.0/rho;

    return rho_inverse*m;
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
template<typename Scalar,
         typename Vector = Eigen::Matrix<Scalar,3,1>,
         typename Tensor = Eigen::Matrix<Scalar,3,3> >
Tensor grad_u(
        const Scalar &rho,
        const Vector &grad_rho,
        const Vector &m,
        const Tensor &grad_m)
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
template<typename Scalar,
         typename Vector = Eigen::Matrix<Scalar,3,1> >
Scalar div_u(
        const Scalar &rho,
        const Vector &grad_rho,
        const Vector &m,
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
template<typename Scalar,
         typename Vector = Eigen::Matrix<Scalar,3,1>,
         typename Tensor = Eigen::Matrix<Scalar,3,3> >
Vector grad_div_u(
        const Scalar &rho,
        const Vector &grad_rho,
        const Tensor &grad_grad_rho,
        const Vector &m,
        const Scalar &div_m,
        const Tensor &grad_m,
        const Vector &grad_div_m)
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
template<typename Scalar,
         typename Vector = Eigen::Matrix<Scalar,3,1>,
         typename Tensor = Eigen::Matrix<Scalar,3,3> >
Vector div_grad_u(
        const Scalar &rho,
        const Vector &grad_rho,
        const Scalar &div_grad_rho,
        const Vector &m,
        const Tensor &grad_m,
        const Vector &div_grad_m)
{
    const Scalar rho_inverse = 1.0/rho;

    return rho_inverse*(
                div_grad_m + rho_inverse*(
                    (2.0*rho_inverse*grad_rho.squaredNorm() - div_grad_rho)*m
                   - 2.0*grad_m*grad_rho
                )
            );
}

} // namespace rhome

} // namespace orthonormal

} // namespace suzerain

} // namespace pecos


#endif // PECOS_SUZERAIN_ORTHONORMAL_H
