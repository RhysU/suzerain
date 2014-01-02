//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
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

#ifndef SUZERAIN_RHOLT_HPP
#define SUZERAIN_RHOLT_HPP

/** @file
 * Provides compute kernels for formulations using a reference density
 * \f$\rho_0\f$, length \f$l_0\f$, and temperature \f$T_0\f$.
 */

#include <suzerain/common.hpp>

namespace suzerain {

/**
 * Provides compute kernels for formulations using a reference density
 * \f$\rho_0\f$, length \f$l_0\f$, and temperature \f$T_0\f$.  The ideal gas
 * equation of state is employed:
 * \f{align*}
 *     p &= \left(\gamma-1\right) \left(
 *       e - \frac{m\cdot{}m}{2\rho}
 *     \right)
 *     \\
 *     T &= \gamma{} \frac{p}{\rho}
 *     \\
 *     \mu &= T^{\beta}
 *     \\
 *     \lambda &= \left(\alpha - \frac{2}{3}\right) \mu
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
 *   - \f$\alpha\f$ or \c alpha is the scaling factor relating \f$\lambda\f$
 *          and \f$ \mu \f$.
 *   - \f$\beta\f$ or \c beta is the constant coefficient in the power
 *          viscosity law.
 *   - \f$\gamma\f$ or \c gamma is the constant ratio of specific heats
 *         \f$C_p/C_v\f$.
 *
 * Templated \c Scalar, \c Vector, \c and Tensor types
 * are used.  These functions are intended to be used with <a
 * href="http://eigen.tuxfamily.org/">Eigen</a>'s vector and matrix
 * types.  For example, when <tt>Scalar == double</tt> \c Vector
 * and \c Tensor may be <tt>Eigen::Matrix<Scalar,3,1></tt> and
 * <tt>Eigen::Matrix<Scalar,3,3></tt>, respectively.
 */
namespace rholt {

/**
 * Compute \f$\tau =
 *    \mu\left(
 *        \vec{\nabla}\vec{u}
 *      + \left(\vec{\nabla}\vec{u}\right)^{\mathsf{T}}
 *    \right)
 *  + \lambda \left( \vec{\nabla}\cdot\vec{u} \right) I\f$.
 *
 * @param[in] mu \f$\mu\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] lambda \f$\lambda\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] div_u \f$\vec{\nabla}\cdot\vec{u}\f$
 *            computed from, for example, div_u()
 * @param[in] grad_u \f$\vec{\nabla}\vec{u}\f$
 *            computed from, for example, grad_u()
 *
 * @return The viscous stress tensor based on the provided fields.
 */
template<typename Scalar,
         typename Tensor >
inline
Tensor tau(
        const Scalar &mu,
        const Scalar &lambda,
        const Scalar &div_u,
        const Tensor &grad_u)
{
    return mu*(grad_u + grad_u.transpose()) + lambda*div_u*Tensor::Identity();
}

/**
 * Compute
 * \f$\vec{\nabla}\cdot\left(\frac{1}{\rho}\vec{m}\otimes\vec{m}\right)\f$.
 * Uses the expansion
 * \f[
 *      \vec{\nabla}\cdot\left(\frac{1}{\rho}\vec{m}\otimes\vec{m}\right)
 *      = \rho^{-1} \left[
 *            \left(\vec{\nabla}\vec{m}\right)\vec{m}
 *          + \left(
 *              \vec{\nabla}\cdot\vec{m}
 *              - \rho^{-1} \vec{m}\cdot\vec{\nabla}\rho
 *            \right)\vec{m}
 *        \right]
 * \f]
 *
 * @param[in] rho \f$\rho\f$
 * @param[in] grad_rho \f$\vec{\nabla}\rho\f$
 * @param[in] m \f$\vec{m}\f$
 * @param[in] div_m \f$\vec{\nabla}\cdot\vec{m}\f$
 * @param[in] grad_m \f$\vec{\nabla}\vec{m}\f$
 *
 * @return The convective derivative computed from the divergence form.
 */
template<typename Scalar,
         typename Vector,
         typename Tensor >
inline
Vector div_rho_inverse_m_outer_m(
        const Scalar &rho,
        const Vector &grad_rho,
        const Vector &m,
        const Scalar &div_m,
        const Tensor &grad_m)
{
    const Scalar rho_inverse = 1/rho;

    return rho_inverse*(
                  grad_m*m
                + (div_m - rho_inverse*m.dot(grad_rho))*m
            );
}

/**
 * Compute
 * \f$\vec{\nabla}\cdot\left(\vec{u}\otimes\vec{m}\right)\f$.
 * Uses the expansion
 * \f[
 *      \vec{\nabla}\cdot\left(\vec{u}\otimes\vec{m}\right)
 *      =   \left(\vec{\nabla}\vec{m}\right)\vec{u}
 *        + \left(\vec{\nabla}\cdot\vec{u}\right)\vec{m}
 * \f]
 *
 * @param[in] m \f$\vec{m}\f$
 * @param[in] grad_m \f$\vec{\nabla}\vec{m}\f$
 * @param[in] u \f$\vec{u}\f$
 *            computed from, for example, u()
 * @param[in] div_u \f$\vec{\nabla}\cdot\vec{u}\f$
 *            computed from, for example, div_u()
 *
 * @return The convective derivative computed from the divergence form.
 */
template<typename Scalar,
         typename Vector,
         typename Tensor >
inline
Vector div_u_outer_m(
        const Vector &m,
        const Tensor &grad_m,
        const Vector &u,
        const Scalar &div_u)
{
    return grad_m*u + div_u*m;
}

/**
 * Compute \f$\vec{\nabla}\cdot\tau\f$.
 * Uses the expansion
 * \f[
 *      \vec{\nabla}\cdot\tau =
 *        \left[
 *           \vec{\nabla}\vec{u}
 *         + \left(\vec{\nabla}\vec{u}\right)^{\mathsf{T}}
 *        \right] \vec{\nabla}\mu
 *      + \mu \vec{\nabla}\cdot\vec{\nabla}\vec{u}
 *      + \left(\mu+\lambda\right)\vec{\nabla}\vec{\nabla}\cdot\vec{u}
 *      + \left(\vec{\nabla}\cdot\vec{u}\right)\vec{\nabla}\lambda
 * \f]
 *
 * @param[in] mu \f$\mu\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] grad_mu \f$\vec{\nabla}\mu\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] lambda \f$\lambda\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] grad_lambda \f$\vec{\nabla}\lambda\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] div_u \f$\vec{\nabla}\cdot\vec{u}\f$
 *            computed from, for example, div_u()
 * @param[in] grad_u \f$\vec{\nabla}\vec{u}\f$
 *            computed from, for example, grad_u()
 * @param[in] div_grad_u \f$\vec{\nabla}\cdot\vec{\nabla}\vec{u}\f$
 *            computed from, for example, div_grad_u()
 * @param[in] grad_div_u \f$\vec{\nabla}\vec{\nabla}\cdot\vec{u}\f$
 *            computed from, for example, grad_div_u()
 *
 * @return The divergence of the viscous stress tensor based on the provided
 *         fields.
 */
template<typename Scalar,
         typename Vector,
         typename Tensor >
inline
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
 *            computed from, for example, u()
 * @param[in] div_u \f$\vec{\nabla}\cdot\vec{u}\f$
 *            computed from, for example, div_u()
 *
 * @return The divergence of the product of total energy and velocity.
 */
template<typename Scalar,
         typename Vector >
inline
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
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] div_grad_T \f$\vec{\nabla}\cdot\vec{\nabla}T\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] mu \f$\mu\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] grad_mu \f$\vec{\nabla}\mu\f$
 *            computed from, for example, p_T_mu_lambda()
 *
 * @return The divergence of the product of viscosity and
 *      the temperature gradient.
 */
template<typename Scalar,
         typename Vector >
inline
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
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] grad_p \f$\vec{\nabla}p\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] u \f$\vec{u}\f$
 *            computed from, for example, u()
 * @param[in] div_u \f$\vec{\nabla}\cdot\vec{u}\f$
 *            computed from, for example, div_u()
 *
 * @return The divergence of the product of pressure and velocity.
 */
template<typename Scalar,
         typename Vector >
inline
Scalar div_p_u(
        const Scalar &p,
        const Vector &grad_p,
        const Vector &u,
        const Scalar &div_u)
{
    return p*div_u + u.dot(grad_p);
}

/**
 * Compute \f$\vec{\nabla}\cdot\tau\vec{u}\f$.
 * Uses the expansion
 * \f[
 *      \vec{\nabla}\cdot\tau\vec{u} =
 *      \vec{u}\cdot\left(
 *          \vec{\nabla}\cdot\tau
 *      \right)
 *      +
 *      \operatorname{trace}\left(
 *          \tau^{\mathsf{T}}\,\vec{\nabla}\vec{u}
 *      \right)
 * \f]
 *
 * @param[in] u \f$\vec{u}\f$
 *            computed from, for example, u()
 * @param[in] grad_u \f$\vec{\nabla}\vec{u}\f$
 *            computed from, for example, grad_u()
 * @param[in] tau \f$\tau\f$
 *            computed from, for example, tau()
 * @param[in] div_tau \f$\vec{\nabla}\cdot\tau\f$
 *            computed from, for example, div_tau()
 *
 * @note Compilers may have trouble automatically deducing the Scalar type
 * because no Scalar appears in the function signature.  Explicitly providing
 * it via <tt>div_tau_u<Scalar></tt> should resolve any 'no instance of
 * function template' errors encountered.
 *
 * @return The divergence of the viscous stress tensor applied to the
 *      velocity.
 */
template<typename Scalar,
         typename Vector,
         typename Tensor >
inline
Scalar div_tau_u(
        const Vector &u,
        const Tensor &grad_u,
        const Tensor &tau,
        const Vector &div_tau)
{
    return u.dot(div_tau) + tau.cwiseProduct(grad_u).sum();
}

/**
 * Compute the explicit part of
 * \f$\vec{\nabla}\cdot\left(\frac{1}{\rho}\vec{m}\otimes\vec{m}\right)\f$.
 * Uses the expansion
 * \f[
 *      \vec{\nabla}\cdot\left(\frac{1}{\rho}\vec{m}\otimes\vec{m}\right)
 *      = \vec{\nabla}\vec{m}\left(
 *            \vec{u} - \left\{\vec{u}\right\}_0
 *        \right)
 *      + \left(\vec{\nabla}\cdot\vec{m}\right)\left(
 *            \vec{u} - \left\{\vec{u}\right\}_0
 *        \right)
 *      - \left(
 *            \vec{u}\otimes\vec{u} - \left\{\vec{u}\otimes\vec{u}\right\}_0
 *        \right) \vec{\nabla}\rho
 * \f]
 * where \f$\left\{\rho^{-1}\vec{m}\right\}_0\f$ and
 * \f$\left\{\rho^{-1}\vec{m}\otimes\rho^{-1}\vec{m}\right\}_0\f$ are based on
 * the same reference velocity.  The remaining linear portion of
 * \f$\vec{\nabla}\cdot\left(\frac{1}{\rho}\vec{m}\otimes\vec{m}\right)\f$ is
 * \f[
 *        \vec{\nabla}\vec{m} \left\{\vec{u}\right\}_0
 *      + \left(\vec{\nabla}\cdot\vec{m}\right) \left\{\vec{u}\right\}_0
 *      - \left\{\vec{u}\otimes\vec{u}\right\}_0 \vec{\nabla}\rho
 * \f]
 *
 * @param[in] grad_rho           \f$\vec{\nabla}\rho\f$
 * @param[in] div_m              \f$\vec{\nabla}\cdot\vec{m}\f$
 * @param[in] grad_m             \f$\vec{\nabla}\vec{m}\f$
 * @param[in] u                  \f$\vec{u}\f$ computed from, for example, u()
 * @param[in] refcoeff_u         \f$\left\{\vec{u}\right\}_0\f$
 * @param[in] refcoeff_u_outer_u \f$\left\{\vec{u}\otimes\vec{u}\right\}_0\f$
 *
 * @return The convective derivative computed from the divergence form.
 */
template<typename Scalar,
         typename Vector,
         typename Tensor,
         typename VectorCoefficient,
         typename TensorCoefficient >
inline
Vector explicit_div_rho_inverse_m_outer_m(
        const Vector &grad_rho,
        const Scalar &div_m,
        const Tensor &grad_m,
        const Vector &u,
        const VectorCoefficient &refcoeff_u,
        const TensorCoefficient &refcoeff_u_outer_u)
{
// FIXME symmetry of refcoeff_u_outer_u should be used in the computation
    const Vector coeff_u = u - refcoeff_u;

    return grad_m*coeff_u
         + div_m*coeff_u
         - (u*u.transpose() - refcoeff_u_outer_u)*grad_rho;
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}\cdot\vec{m}\f$ in the explicit portion
 * of \f$\vec{\nabla}\cdot{}p\vec{u}\f$.
 *
 * @param[in] rho \f$\rho\f$
 * @param[in] p \f$p\f$
 *
 * @return \f$\rho^{-1}p\f$
 * @see explicit_div_p_u for more details on the explicit
 *      operator.
 */
template<typename Scalar>
inline
Scalar explicit_div_p_u_refcoeff_div_m(
        const Scalar &rho,
        const Scalar &p)
{
    return p/rho;
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}\rho\f$ in the explicit portion
 * of \f$\vec{\nabla}\cdot{}p\vec{u}\f$.
 *
 * @param[in] rho \f$\rho\f$
 * @param[in] m \f$\vec{m}\f$
 * @param[in] p \f$p\f$
 *
 * @return \f$\rho^{-2}p\vec{m}\f$
 * @see explicit_div_p_u for more details on the explicit
 *      operator.
 */
template<typename Scalar,
         typename Vector >
inline
Vector explicit_div_p_u_refcoeff_grad_rho(
        const Scalar &rho,
        const Vector &m,
        const Scalar &p)
{
    const Scalar rho_inverse = 1/rho;
    return rho_inverse*rho_inverse*p*m;
}

/**
 * Compute the explicit portion of
 * \f$\vec{\nabla}\cdot{}p\vec{u}\f$.
 * Uses the expansion
 * \f{align*}
 * \vec{\nabla}\cdot{}p\vec{u} &=
 *     \rho^{-1}\vec{m}\cdot\vec{\nabla}p
 *   + \left(
 *           \rho^{-1}p
 *         - \left\{\rho^{-1}p\right\}_0
 *     \right)\vec{\nabla}\cdot\vec{m}
 *   - \left(
 *           \rho^{-2}p\vec{m}
 *         - \left\{\rho^{-2}p\vec{m}\right\}_0
 *     \right)\cdot\vec{\nabla}\rho
 * \f}
 * where \f$\left\{\rho^{-1}p\right\}_0\f$
 * and \f$\left\{\rho^{-2}p\vec{m}\right\}_0\f$ are fixed by
 * \c refcoeff_div_m and \c refcoeff_grad_rho, respectively.
 * The remaining linear portion of
 * \f$\vec{\nabla}\cdot{}p\vec{u}\f$ is
 * \f[
 *      \left\{\rho^{-1}p\right\}_0 \vec{\nabla}\vec{m}
 *    - \left\{\rho^{-2}p\vec{m}\right\}_0\cdot\vec{\nabla}\rho
 * \f]
 *
 * @param rho \f$\rho\f$
 * @param grad_rho \f$\vec{\nabla}\rho\f$
 * @param m \f$\vec{m}\f$
 * @param div_m \f$\vec{\nabla}\cdot\vec{m}\f$
 * @param p \f$p\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param grad_p \f$\vec{\nabla}p\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param refcoeff_div_m the reference coefficient
 *        on \f$\vec{\nabla}\cdot\vec{m}\f$ which may
 *        be computed using explicit_div_p_u_refcoeff_div_m()
 * @param refcoeff_grad_rho the reference coefficient
 *        on \f$\vec{\nabla}\rho\f$ which may
 *        be computed using explicit_div_p_u_refcoeff_grad_rho()
 *
 * @return The explicit portion of divergence of the pressure
 *         times velocity.
 */
template<typename Scalar,
         typename Vector,
         typename ScalarCoefficient,
         typename VectorCoefficient >
inline
Scalar explicit_div_p_u(
        const Scalar            &rho,
        const Vector            &grad_rho,
        const Vector            &m,
        const Scalar            &div_m,
        const Scalar            &p,
        const Vector            &grad_p,
        const ScalarCoefficient &refcoeff_div_m,
        const VectorCoefficient &refcoeff_grad_rho)
{
    const Scalar coeff_div_m(
            explicit_div_p_u_refcoeff_div_m(rho, p) - refcoeff_div_m);
    const Vector coeff_grad_rho(
            explicit_div_p_u_refcoeff_grad_rho(rho, m, p) - refcoeff_grad_rho);

    return   m.dot(grad_p)/rho
           + coeff_div_m*div_m
           - coeff_grad_rho.dot(grad_rho);
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}\rho\f$ in the explicit portion of \f$\vec{\nabla}p\f$.
 *
 * @param[in] rho \f$\rho\f$
 * @param[in] m \f$\vec{m}\f$
 *
 * @return \f$\frac{\vec{m}\cdot\vec{m}}{\rho^{2}}\f$
 * @see explicit_grad_p for more details on the explicit operator.
 */
template<typename Scalar,
         typename Vector >
inline
Scalar explicit_grad_p_refcoeff_grad_rho(
        const Scalar &rho,
        const Vector &m)
{
    const Scalar rho_inverse = 1/rho;

    return m.squaredNorm()*rho_inverse*rho_inverse;
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}\vec{m}\f$ in the explicit portion of \f$\vec{\nabla}p\f$.
 *
 * @param[in] rho \f$\rho\f$
 * @param[in] m \f$\vec{m}\f$
 *
 * @return \f$\frac{\vec{m}}{\rho}\f$
 * @see explicit_grad_p for more details on the explicit operator.
 */
template<typename Scalar,
         typename Vector >
inline
Vector explicit_grad_p_refcoeff_grad_m(
        const Scalar &rho,
        const Vector &m)
{
    return m / rho;
}

/**
 * Compute the explicit portion of \f$\vec{\nabla}p\f$.
 * Uses the expansion
 * \f{align*}
 *   \nabla{}p &= \left(\gamma-1\right)\left(
 *       \frac{1}{2}\left(
 *           \frac{\vec{m}\cdot\vec{m}}{\rho^{2}}
 *         - \left\{\frac{\vec{m}\cdot\vec{m}}{\rho^{2}}\right\}_0
 *       \right)\vec{\nabla}\rho
 *     - \left(\vec{\nabla}\vec{m}\right)^{\mathsf{T}}\left(
 *           \frac{\vec{m}}{\rho}
 *         - \left\{\frac{\vec{m}}{\rho}\right\}_0
 *       \right)
 *   \right)
 * \f}
 * where \f$
 * \left\{\frac{\vec{m}\cdot\vec{m}}{\rho^{2}}\right\}_{0} \f$ and \f$
 * \left\{\frac{m}{\rho}\right\}_{0} \f$ are fixed by \c refcoeff_grad_rho and
 * \c refcoeff_grad_m, respectively.
 * The remaining linear portion of \f$\vec{\nabla}p\f$ is
 * \f[
 *     \left(\gamma-1\right)\vec{\nabla}e
 *   + \frac{\gamma-1}{2}
 *     \left\{\frac{\vec{m}\cdot\vec{m}}{\rho^{2}}\right\}_{0} \vec{\nabla}\rho
 *   - \left(\gamma-1\right)
 *     \left(\vec{\nabla}\vec{m}\right)^{\mathsf{T}}
 *     \left\{\frac{m}{\rho}\right\}_{0}
 * \f]
 * Unlike many other methods beginning with <tt>explicit_</tt>, when zero
 * reference coefficients are used the full \f$\vec{\nabla}p\f$ term is
 * <i>not</i> recovered.
 *
 * @param gamma \f$\gamma\f$
 * @param rho \f$\rho\f$
 * @param grad_rho \f$\vec{\nabla}\rho\f$
 * @param m \f$\vec{m}\f$
 * @param grad_m \f$\vec{\nabla}\vec{m}\f$
 * @param refcoeff_grad_rho the reference coefficient on \f$\vec{\nabla}\rho\f$
 *        which may be computed using explicit_grad_p_refcoeff_grad_rho.
 * @param refcoeff_grad_m the reference coefficient on \f$\vec{\nabla}\vec{m}\f$
 *        which may be computed using explicit_grad_p_refcoeff_grad_m.
 *
 * @return The explicit portion of the gradient of the pressure.
 */
template<typename Scalar,
         typename Vector,
         typename Tensor,
         typename ScalarCoefficient,
         typename VectorCoefficient >
inline
Vector explicit_grad_p(
        const Scalar            &gamma,
        const Scalar            &rho,
        const Vector            &grad_rho,
        const Vector            &m,
        const Tensor            &grad_m,
        const ScalarCoefficient &refcoeff_grad_rho,
        const VectorCoefficient &refcoeff_grad_m)
{
    const Scalar coeff_grad_rho(
            explicit_grad_p_refcoeff_grad_rho(rho, m) - refcoeff_grad_rho);
    const Vector coeff_grad_m(
            explicit_grad_p_refcoeff_grad_m(rho, m) - refcoeff_grad_m);

    return (gamma-1)*(
                  (coeff_grad_rho/2)*grad_rho
                - grad_m.transpose()*coeff_grad_m
            );
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}\cdot\vec{m}\f$ in the explicit portion
 * of \f$\vec{\nabla}\cdot\left(e+p\right)\vec{u}\f$.
 *
 * @param[in] rho \f$\rho\f$
 * @param[in] e   \f$e\f$
 * @param[in] p   \f$p\f$
 *
 * @return \f$\rho^{-1}\left(e+p\right)\f$
 * @see explicit_div_e_plus_p_u for more details on the explicit
 *      operator.
 */
template<typename Scalar>
inline
Scalar explicit_div_e_plus_p_u_refcoeff_div_m(
        const Scalar &rho,
        const Scalar &e,
        const Scalar &p)
{
    return (e + p) / rho;
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}\rho\f$ in the explicit portion
 * of \f$\vec{\nabla}\cdot\left(e+p\right)\vec{u}\f$.
 *
 * @param[in] gamma \f$\gamma\f$
 * @param[in] rho   \f$\rho\f$
 * @param[in] m     \f$\vec{m}\f$
 * @param[in] e     \f$e\f$
 * @param[in] p     \f$p\f$
 *
 * @return \f$\rho^{-2}\vec{m}\left((\gamma-2)e-2p\right)\f$
 * @see explicit_div_e_plus_p_u for more details on the explicit
 *      operator.
 */
template<typename Scalar,
         typename Vector >
inline
Vector explicit_div_e_plus_p_u_refcoeff_grad_rho(
        const Scalar &gamma,
        const Scalar &rho,
        const Vector &m,
        const Scalar &e,
        const Scalar &p)
{
    return (((gamma-2)*e-2*p)/(rho*rho))*m;
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}e\f$ in the explicit portion
 * of \f$\vec{\nabla}\cdot\left(e+p\right)\vec{u}\f$.
 *
 * @param[in] rho \f$\rho\f$
 * @param[in] m \f$\vec{m}\f$
 *
 * @return \f$\rho^{-1}\vec{m}\f$
 * @see explicit_div_e_plus_p_u for more details on the explicit
 *      operator.
 */
template<typename Scalar,
         typename Vector >
inline
Vector explicit_div_e_plus_p_u_refcoeff_grad_e(
        const Scalar &rho,
        const Vector &m)
{
    return m/rho;
}

/**
 * Compute the explicit portion of
 * \f$\vec{\nabla}\cdot\left(e+p\right)\vec{u}\f$.
 * Uses the expansion
 * \f{align*}
 * \vec{\nabla}\cdot\left(e+p\right)\vec{u} =
 *  &- \left(\gamma-1\right)\rho^{-2}\vec{m}\cdot
 *     \left(\vec{\nabla}\vec{m}\right)^{\mathsf{T}}\vec{m}
 * \\
 *  &+ \left(
 *                 \rho^{-1}\left(e+p\right)
 *        - \left\{\rho^{-1}\left(e+p\right)\right\}_0
 *     \right)\vec{\nabla}\cdot\vec{m}
 * \\
 *  &+ \left(
 *                 \rho^{-2}\vec{m}\left((\gamma-2)e-2p\right)
 *        - \left\{\rho^{-2}\vec{m}\left((\gamma-2)e-2p\right)\right\}_0
 *     \right)\cdot\vec{\nabla}\rho
 * \\
 *  &+ \gamma\left(
 *          \rho^{-1}\vec{m}
 *        - \left\{\rho^{-1}\vec{m}\right\}_0
 *     \right)\cdot\vec{\nabla}e
 * \f}
 * where
 * \f$\left\{\rho^{-1}\left(e+p\right)\right\}_0\f$
 *  \f$\left\{\rho^{-2}\vec{m}\left((\gamma-2)e-2p\right)\right\}_0\f$, and
 *  \f$\left\{\rho^{-1}\vec{m}\right\}_0\f$ are fixed by \c
 *  refcoeff_div_m, \c refcoeff_grad_rho, and
 * \c refcoeff_grad_e, respectively.
 * The remaining linear portion of
 * \f$\vec{\nabla}\cdot\left(e+p\right)\vec{u}\f$ is
 * \f[
 *      \left\{\rho^{-1}\left(e+p\right)\right\}_0
 *        \vec{\nabla}\cdot\vec{m}
 *    + \left\{\rho^{-2}\vec{m}\left((\gamma-2)e-2p\right)\right\}_0
 *        \cdot\vec{\nabla}\rho
 *    + \gamma\left\{\rho^{-1}\vec{m}\right\}_0
 *        \cdot\vec{\nabla}e
 * \f]
 *
 * @param gamma \f$\gamma\f$
 * @param rho \f$\rho\f$
 * @param grad_rho \f$\vec{\nabla}\rho\f$
 * @param m \f$\vec{m}\f$
 * @param div_m \f$\vec{\nabla}\cdot\vec{m}\f$
 * @param grad_m \f$\vec{\nabla}\vec{m}\f$
 * @param e \f$e\f$
 * @param grad_e \f$\vec{\nabla}e\f$
 * @param p \f$p\f$
 * @param refcoeff_div_m the reference coefficient
 *        on \f$\vec{\nabla}\cdot\vec{m}\f$ which may
 *        be computed using explicit_div_e_plus_p_u_refcoeff_div_m()
 * @param refcoeff_grad_rho the reference coefficient
 *        on \f$\vec{\nabla}\rho\f$ which may
 *        be computed using explicit_div_e_plus_p_u_refcoeff_grad_rho()
 * @param refcoeff_grad_e the reference coefficient
 *        on \f$\vec{\nabla}e\f$ which may
 *        be computed using explicit_div_e_plus_p_u_refcoeff_grad_e()
 *
 * @return The explicit portion of divergence of the energy plus
 *         the pressure times velocity.
 */
template<typename Scalar,
         typename Vector,
         typename Tensor,
         typename ScalarCoefficient,
         typename VectorCoefficient1,
         typename VectorCoefficient2 >
inline
Scalar explicit_div_e_plus_p_u(
        const Scalar             &gamma,
        const Scalar             &rho,
        const Vector             &grad_rho,
        const Vector             &m,
        const Scalar             &div_m,
        const Tensor             &grad_m,
        const Scalar             &e,
        const Vector             &grad_e,
        const Scalar             &p,
        const ScalarCoefficient  &refcoeff_div_m,
        const VectorCoefficient1 &refcoeff_grad_rho,
        const VectorCoefficient2 &refcoeff_grad_e)
{
    const Scalar coeff_div_m(
            explicit_div_e_plus_p_u_refcoeff_div_m(rho, e, p)
          - refcoeff_div_m);
    const Vector coeff_grad_rho(
            explicit_div_e_plus_p_u_refcoeff_grad_rho(gamma, rho, m, e, p)
          - refcoeff_grad_rho);
    const Vector coeff_grad_e(
            explicit_div_e_plus_p_u_refcoeff_grad_e(rho, m)
          - refcoeff_grad_e);

    return   (1-gamma)/(rho*rho)*m.dot(grad_m.transpose()*m)
           + coeff_div_m*div_m
           + coeff_grad_rho.dot(grad_rho)
           + gamma*coeff_grad_e.dot(grad_e);
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}\cdot\vec{\nabla}\vec{m}\f$ in the explicit portion
 * of \f$\mu\vec{\nabla}\cdot\vec{\nabla}\vec{u}\f$.
 *
 * @param[in] mu \f$\mu\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] rho \f$\rho\f$
 *
 * @return \f$\mu\rho^{-1}\f$
 * @see explicit_mu_div_grad_u for more details on the explicit
 *      operator.
 */
template<typename Scalar>
inline
Scalar explicit_mu_div_grad_u_refcoeff_div_grad_m(
        const Scalar &mu,
        const Scalar &rho)
{
    return mu/rho;
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}\cdot\vec{\nabla}\rho\f$ in the explicit portion
 * of \f$\mu\vec{\nabla}\cdot\vec{\nabla}\vec{u}\f$.
 *
 * @param[in] mu \f$\mu\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] rho \f$\rho\f$
 * @param[in] m \f$\vec{m}\f$
 *
 * @return \f$\mu\rho^{-2}\vec{m}\f$
 * @see explicit_mu_div_grad_u for more details on the explicit
 *      operator.
 */
template<typename Scalar,
         typename Vector >
inline
Vector explicit_mu_div_grad_u_refcoeff_div_grad_rho(
        const Scalar &mu,
        const Scalar &rho,
        const Vector &m)
{
    const Scalar rho_inverse = 1/rho;
    return mu*rho_inverse*rho_inverse*m;
}

/**
 * Compute the explicit portion of
 * \f$\mu\vec{\nabla}\cdot\vec{\nabla}\vec{u}\f$.
 * Uses the expansion
 * \f{align*}
 * \mu\vec{\nabla}\cdot\vec{\nabla}\vec{u} = &\phantom{{}+}
 *     2\mu\rho^{-2}\left[
 *           \rho^{-1}\vec{m}\left({\vec\nabla}\rho\right)^{2}
 *         - \left(\vec{\nabla}\vec{m}\right)\vec{\nabla}\rho
 *     \right]
 * \\
 *   &{}+ \left(
 *              \mu\rho^{-1}
 *            - \left\{\mu\rho^{-1}\right\}_0
 *        \right)\vec{\nabla}\cdot\vec{\nabla}\vec{m}
 *      - \left(
 *              \mu\rho^{-2}\vec{m}
 *            - \left\{\mu\rho^{-2}\vec{m}\right\}_0
 *        \right)\vec{\nabla}\cdot\vec{\nabla}\rho
 * \f}
 * where \f$\left\{\mu\rho^{-1}\right\}_0\f$
 * and \f$\left\{\mu\rho^{-2}\vec{m}\right\}_0\f$ are fixed by
 * \c refcoeff_div_grad_m and \c refcoeff_div_grad_rho, respectively.
 * The remaining linear portion of
 * \f$\mu\vec{\nabla}\cdot\vec{\nabla}\vec{u}\f$ is
 * \f[
 *      \left\{\mu\rho^{-1}\right\}_0 \vec{\nabla}\cdot\vec{\nabla}\vec{m}
 *    - \left\{\mu\rho^{-2}\vec{m}\right\}_0 \vec{\nabla}\cdot\vec{\nabla}\rho
 * \f]
 *
 * @param mu \f$\mu\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param rho \f$\rho\f$
 * @param grad_rho \f$\vec{\nabla}\rho\f$
 * @param div_grad_rho \f$\vec{\nabla}\cdot\vec{\nabla}\rho\f$
 * @param m \f$\vec{m}\f$
 * @param grad_m \f$\vec{\nabla}\vec{m}\f$
 * @param div_grad_m \f$\vec{\nabla}\cdot\vec{\nabla}\vec{m}\f$
 * @param refcoeff_div_grad_m the reference coefficient
 *        on \f$\vec{\nabla}\cdot\vec{\nabla}\vec{m}\f$ which may
 *        be computed using explicit_mu_div_grad_u_refcoeff_div_grad_m()
 * @param refcoeff_div_grad_rho the reference coefficient
 *        on \f$\vec{\nabla}\cdot\vec{\nabla}\rho\f$ which may
 *        be computed using explicit_mu_div_grad_u_refcoeff_div_grad_rho()
 *
 * @return The explicit portion of the viscosity times the Laplacian
 *         of velocity.
 */
template<typename Scalar,
         typename Vector,
         typename Tensor,
         typename ScalarCoefficient,
         typename VectorCoefficient >
inline
Vector explicit_mu_div_grad_u(
        const Scalar            &mu,
        const Scalar            &rho,
        const Vector            &grad_rho,
        const Scalar            &div_grad_rho,
        const Vector            &m,
        const Tensor            &grad_m,
        const Vector            &div_grad_m,
        const ScalarCoefficient &refcoeff_div_grad_m,
        const VectorCoefficient &refcoeff_div_grad_rho)
{
    const Scalar rho_inverse  = 1/rho;
    const Scalar rho_inverse2 = rho_inverse*rho_inverse;
    const Scalar coeff_div_grad_m(
            explicit_mu_div_grad_u_refcoeff_div_grad_m(mu, rho)
          - refcoeff_div_grad_m);
    const Vector coeff_div_grad_rho(
            explicit_mu_div_grad_u_refcoeff_div_grad_rho(mu, rho, m)
          - refcoeff_div_grad_rho);

    return   2*mu*rho_inverse2*(
                    rho_inverse*grad_rho.squaredNorm()*m
                  - grad_m*grad_rho
             )
           + coeff_div_grad_m  *div_grad_m
           - coeff_div_grad_rho*div_grad_rho;
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}\cdot\vec{\nabla}\vec{m}\f$ in the explicit portion
 * of \f$\vec{u}\cdot\mu\vec{\nabla}\cdot\vec{\nabla}\vec{u}\f$.
 *
 * @param[in] mu \f$\mu\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] rho \f$\rho\f$
 * @param[in] m  \f$m\f$
 *
 * @return \f$\mu\rho^{-2}\vec{m}\f$
 * @see explicit_u_dot_mu_div_grad_u for more details on the explicit
 *      operator.
 */
template<typename Scalar,
         typename Vector>
inline
Vector explicit_u_dot_mu_div_grad_u_refcoeff_div_grad_m(
        const Scalar &mu,
        const Scalar &rho,
        const Vector &m)
{
    return mu/(rho*rho)*m;
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}\cdot\vec{\nabla}\rho\f$ in the explicit portion
 * of \f$\vec{u}\cdot\mu\vec{\nabla}\cdot\vec{\nabla}\vec{u}\f$.
 *
 * @param[in] mu \f$\mu\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] rho \f$\rho\f$
 * @param[in] m \f$\vec{m}\f$
 *
 * @return \f$\mu\rho^{-3}\vec{m}^2\f$
 * @see explicit_u_dot_mu_div_grad_u for more details on the explicit
 *      operator.
 */
template<typename Scalar,
         typename Vector >
inline
Scalar explicit_u_dot_mu_div_grad_u_refcoeff_div_grad_rho(
        const Scalar &mu,
        const Scalar &rho,
        const Vector &m)
{
    return mu*m.squaredNorm()/(rho*rho*rho);
}

/**
 * Compute the explicit portion of
 * \f$\vec{u}\cdot\mu\vec{\nabla}\cdot\vec{\nabla}\vec{u}\f$.
 * Uses the expansion
 * \f{align*}
 * \vec{u}\cdot\mu\vec{\nabla}\cdot\vec{\nabla}\vec{u} = &\phantom{{}+}
 *     2\mu\rho^{-3}\vec{m}\cdot\left[
 *           \rho^{-1}\vec{m}\left({\vec\nabla}\rho\right)^{2}
 *         - \left(\vec{\nabla}\vec{m}\right)\vec{\nabla}\rho
 *     \right]
 * \\
 *   &{}+ \left(
 *              \mu\rho^{-2}\vec{m}
 *            - \left\{\mu\rho^{-2}m\right\}_0
 *        \right)\cdot\vec{\nabla}\cdot\vec{\nabla}\vec{m}
 *      - \left(
 *              \mu\rho^{-3}\vec{m}^2
 *            - \left\{\mu\rho^{-3}\vec{m}^2\right\}_0
 *        \right)\vec{\nabla}\cdot\vec{\nabla}\rho
 * \f}
 * where \f$\left\{\mu\rho^{-2}\vec{m}\right\}_0\f$
 * and \f$\left\{\mu\rho^{-3}\vec{m}^2\right\}_0\f$ are fixed by
 * \c refcoeff_div_grad_m and \c refcoeff_div_grad_rho, respectively.
 * The remaining linear portion of
 * \f$\vec{u}\cdot\mu\vec{\nabla}\cdot\vec{\nabla}\vec{u}\f$ is
 * \f[
 *      \left\{\mu\rho^{-1}\vec{m}\right\}_0
 *      \cdot\vec{\nabla}\cdot\vec{\nabla}\vec{m}
 *    - \left\{\mu\rho^{-3}\vec{m}^2\right\}_0
 *      \vec{\nabla}\cdot\vec{\nabla}\rho
 * \f]
 *
 * @param mu \f$\mu\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param rho \f$\rho\f$
 * @param grad_rho \f$\vec{\nabla}\rho\f$
 * @param div_grad_rho \f$\vec{\nabla}\cdot\vec{\nabla}\rho\f$
 * @param m \f$\vec{m}\f$
 * @param grad_m \f$\vec{\nabla}\vec{m}\f$
 * @param div_grad_m \f$\vec{\nabla}\cdot\vec{\nabla}\vec{m}\f$
 * @param refcoeff_div_grad_m the reference coefficient
 *        on \f$\vec{\nabla}\cdot\vec{\nabla}\vec{m}\f$ which may
 *        be computed using explicit_u_dot_mu_div_grad_u_refcoeff_div_grad_m()
 * @param refcoeff_div_grad_rho the reference coefficient
 *        on \f$\vec{\nabla}\cdot\vec{\nabla}\rho\f$ which may
 *        be computed using explicit_u_dot_mu_div_grad_u_refcoeff_div_grad_rho()
 *
 * @return The explicit portion of the velocity times the viscosity dotted
 *         against the Laplacian of velocity.
 */
template<typename Scalar,
         typename Vector,
         typename Tensor,
         typename ScalarCoefficient,
         typename VectorCoefficient >
inline
Scalar explicit_u_dot_mu_div_grad_u(
        const Scalar            &mu,
        const Scalar            &rho,
        const Vector            &grad_rho,
        const Scalar            &div_grad_rho,
        const Vector            &m,
        const Tensor            &grad_m,
        const Vector            &div_grad_m,
        const VectorCoefficient &refcoeff_div_grad_m,
        const ScalarCoefficient &refcoeff_div_grad_rho)
{
    const Scalar rho_inverse  = 1/rho;
    const Scalar rho_inverse3 = rho_inverse*rho_inverse*rho_inverse;
    const Vector coeff_div_grad_m(
            explicit_u_dot_mu_div_grad_u_refcoeff_div_grad_m(mu, rho, m)
          - refcoeff_div_grad_m);
    const Scalar coeff_div_grad_rho(
            explicit_u_dot_mu_div_grad_u_refcoeff_div_grad_rho(mu, rho, m)
          - refcoeff_div_grad_rho);

    return   2*mu*rho_inverse3*m.dot(
                    rho_inverse*grad_rho.squaredNorm()*m
                  - grad_m*grad_rho
             )
           + coeff_div_grad_m.dot(div_grad_m)
           - coeff_div_grad_rho*div_grad_rho;
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}\vec{\nabla}\cdot\vec{m}\f$ in the explicit portion
 * of \f$\left(\mu+\lambda\right)\vec{\nabla}\vec{\nabla}\cdot\vec{u}\f$.
 *
 * @param[in] mu \f$\mu\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] lambda \f$\lambda\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] rho \f$\rho\f$
 *
 * @return \f$\left(\mu+\lambda\right)\rho^{-1}\f$
 * @see explicit_mu_plus_lambda_grad_div_u() for more details on the explicit
 *      operator.
 */
template<typename Scalar>
inline
Scalar explicit_mu_plus_lambda_grad_div_u_refcoeff_grad_div_m(
        const Scalar &mu,
        const Scalar &lambda,
        const Scalar &rho)
{
    return (mu+lambda)/rho;
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}\vec{\nabla}\rho\f$ in the explicit portion
 * of \f$\left(\mu+\lambda\right)\vec{\nabla}\vec{\nabla}\cdot\vec{u}\f$.
 *
 * @param[in] mu \f$\mu\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] lambda \f$\lambda\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] rho \f$\rho\f$
 * @param[in] m \f$\vec{m}\f$
 *
 * @return \f$\left(\mu+\lambda\right)\rho^{-2}\vec{m}\f$
 * @see explicit_mu_plus_lambda_grad_div_u() for more details on the explicit
 *      operator.
 */
template<typename Scalar,
         typename Vector >
inline
Vector explicit_mu_plus_lambda_grad_div_u_refcoeff_grad_grad_rho(
        const Scalar &mu,
        const Scalar &lambda,
        const Scalar &rho,
        const Vector &m)
{
    const Scalar rho_inverse = 1/rho;
    return (mu+lambda)*rho_inverse*rho_inverse*m;
}

/**
 * Compute the explicit portion of
 * \f$\left(\mu+\lambda\right)\vec{\nabla}\vec{\nabla}\cdot\vec{u}\f$.
 * Uses the expansion
 * \f{align*}
 * \left(\mu+\lambda\right)\vec{\nabla}\vec{\nabla}\cdot\vec{u} =
 *   &\phantom{{}+}
 *    \left(\mu+\lambda\right)\rho^{-2}\left[
 *        \left(
 *              2\rho^{-1}\vec{\nabla}\rho\cdot\vec{m}
 *            - \vec{\nabla}\cdot\vec{m}
 *        \right)\vec{\nabla}\rho
 *      - {\vec{\nabla}\vec{m}}^{\mathsf{T}}\vec{\nabla}\rho
 *    \right]
 * \\
 *   &{}+ \left(
 *              \left(\mu+\lambda\right)\rho^{-1}
 *            - \left\{\left(\mu+\lambda\right)\rho^{-1}\right\}_0
 *        \right) \vec{\nabla}\vec{\nabla}\cdot\vec{m}
 *      - \left(
 *              \left(\mu+\lambda\right)\rho^{-2}\vec{m}
 *            - \left\{\left(\mu+\lambda\right)\rho^{-2}\vec{m}\right\}_0
 *        \right) \vec{\nabla}\vec{\nabla}\rho
 * \\
 * \f}
 * where \f$\left\{\left(\mu+\lambda\right)\rho^{-1}\right\}_0\f$
 * and \f$\left\{\left(\mu+\lambda\right)\rho^{-2}\vec{m}\right\}_0\f$
 * are fixed by \c refcoeff_grad_div_m and \c refcoeff_grad_grad_rho,
 * respectively.
 * The remaining linear portion of
 * \f$\left(\mu+\lambda\right)\vec{\nabla}\vec{\nabla}\cdot\vec{u}\f$ is
 * \f[
 *        \left\{\left(\mu+\lambda\right)\rho^{-1}\right\}_0
 *        \vec{\nabla}\vec{\nabla}\cdot\vec{m}
 *      - \left\{\left(\mu+\lambda\right)\rho^{-2}\vec{m}\right\}_0
 *        \vec{\nabla}\vec{\nabla}\rho
 * \f]
 *
 * @param mu \f$\mu\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param lambda \f$\lambda\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param rho \f$\rho\f$
 * @param grad_rho \f$\vec{\nabla}\rho\f$
 * @param grad_grad_rho \f$\vec{\nabla}\vec{\nabla}\rho\f$
 * @param m \f$\vec{m}\f$
 * @param div_m \f$\vec{\nabla}\cdot\vec{m}\f$
 * @param grad_m \f$\vec{\nabla}\vec{m}\f$
 * @param grad_div_m \f$\vec{\nabla}\vec{\nabla}\cdot\vec{m}\f$
 * @param refcoeff_grad_div_m the reference coefficient
 *        on \f$\vec{\nabla}\vec{\nabla}\cdot\vec{m}\f$ which may
 *        be computed using
 *        explicit_mu_plus_lambda_grad_div_u_refcoeff_grad_div_m()
 * @param refcoeff_grad_grad_rho the reference coefficient
 *        on \f$\vec{\nabla}\vec{\nabla}\rho\f$ which may
 *        be computed using
 *        explicit_mu_plus_lambda_grad_div_u_refcoeff_grad_grad_rho()
 *
 * @return The explicit portion of the sum of the viscosities times
 *         the gradient of the divergence of velocity.
 */
template<typename Scalar,
         typename Vector,
         typename Tensor,
         typename ScalarCoefficient,
         typename VectorCoefficient >
inline
Vector explicit_mu_plus_lambda_grad_div_u(
        const Scalar            &mu,
        const Scalar            &lambda,
        const Scalar            &rho,
        const Vector            &grad_rho,
        const Tensor            &grad_grad_rho,
        const Vector            &m,
        const Scalar            &div_m,
        const Tensor            &grad_m,
        const Vector            &grad_div_m,
        const ScalarCoefficient &refcoeff_grad_div_m,
        const VectorCoefficient &refcoeff_grad_grad_rho)
{
    const Scalar rho_inverse  = 1/rho;
    const Scalar rho_inverse2 = rho_inverse*rho_inverse;
    const Scalar coeff_grad_div_m(
            explicit_mu_plus_lambda_grad_div_u_refcoeff_grad_div_m(
              mu, lambda, rho)
          - refcoeff_grad_div_m);
    const Vector coeff_grad_grad_rho(
            explicit_mu_plus_lambda_grad_div_u_refcoeff_grad_grad_rho(
              mu, lambda, rho, m)
          - refcoeff_grad_grad_rho);

    // TODO Use symmetry of grad_grad_rho
    return   (mu+lambda)*rho_inverse2*(
                  (2*rho_inverse*grad_rho.dot(m) - div_m)*grad_rho
                - grad_m.transpose()*grad_rho
             )
           + coeff_grad_div_m*grad_div_m
           - grad_grad_rho*coeff_grad_grad_rho;
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}\vec{\nabla}\cdot\vec{m}\f$ in the explicit portion of
 * \f$\vec{u}\cdot\left(\mu+\lambda\right)\vec{\nabla}\vec{\nabla}\cdot\vec{u}\f$.
 *
 * @param[in] mu \f$\mu\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] lambda \f$\lambda\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] rho \f$\rho\f$
 * @param[in] m \f$\vec{m}\f$
 *
 * @return \f$\left(\mu+\lambda\right)\rho^{-2}\vec{m}\f$
 * @see explicit_u_dot_mu_plus_lambda_grad_div_u() for more details on the
 *      explicit operator.
 */
template<typename Scalar,
         typename Vector>
inline
Vector explicit_u_dot_mu_plus_lambda_grad_div_u_refcoeff_grad_div_m(
        const Scalar &mu,
        const Scalar &lambda,
        const Scalar &rho,
        const Vector &m)
{
    return (mu+lambda)/(rho*rho)*m;
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}\vec{\nabla}\rho\f$ in the explicit portion of \f$\vec{u}
 * \cdot \left(\mu+\lambda\right) \vec{\nabla}\vec{\nabla}\cdot\vec{u}\f$.
 *
 * @param[in] mu \f$\mu\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] lambda \f$\lambda\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] rho \f$\rho\f$
 * @param[in] m \f$\vec{m}\f$
 *
 * @return \f$\left(\mu+\lambda\right)\rho^{-3}\vec{m}\otimes\vec{m}\f$
 * @see explicit_u_dot_mu_plus_lambda_grad_div_u() for more details on the
 *      explicit operator.
 */
template<typename Scalar,
         typename Vector,
         typename Tensor>
inline
Tensor explicit_u_dot_mu_plus_lambda_grad_div_u_refcoeff_grad_grad_rho(
        const Scalar &mu,
        const Scalar &lambda,
        const Scalar &rho,
        const Vector &m)
{
    return ((mu+lambda)/(rho*rho*rho)*m)*m.transpose();
}

/**
 * Compute the explicit portion of \f$\vec{u} \cdot \left(\mu+\lambda\right)
 * \vec{\nabla}\vec{\nabla}\cdot\vec{u}\f$.  Uses the expansion
 * \f{align*}
 * \vec{u}\cdot\left(\mu+\lambda\right)\vec{\nabla}\vec{\nabla}\cdot\vec{u} =
 *   &\phantom{{}+}
 *    \left(\mu+\lambda\right)\rho^{-3}\vec{m}\cdot\left[
 *        \left(
 *              2\rho^{-1}\vec{\nabla}\rho\cdot\vec{m}
 *            - \vec{\nabla}\cdot\vec{m}
 *        \right)\vec{\nabla}\rho
 *      - {\vec{\nabla}\vec{m}}^{\mathsf{T}}\vec{\nabla}\rho
 *    \right]
 * \\
 *   &{}+ \left(
 *              \left(\mu+\lambda\right)\rho^{-2}\vec{m}
 *            - \left\{\left(\mu+\lambda\right)\rho^{-2}\vec{m}\right\}_0
 *        \right)\cdot\vec{\nabla}\vec{\nabla}\cdot\vec{m}
 * \\
 *   &{}- \operatorname{trace}\left[
 *          {\vec{\nabla}\vec{\nabla}\rho}^{\mathsf{T}}
 *          \left(
 *                \left(\mu+\lambda\right)\rho^{-3}
 *                \vec{m}\otimes\vec{m}
 *              - \left\{\left(\mu+\lambda\right)
 *                \rho^{-3}\vec{m}\otimes\vec{m}\right\}_0
 *          \right)
 *        \right]
 * \\
 * \f}
 * where \f$\left\{\left(\mu+\lambda\right)\rho^{-2}\vec{m}\right\}_0\f$ and
 * \f$\left\{\left(\mu+\lambda\right)
 * \rho^{-3}\vec{m}\otimes\vec{m}\right\}_0\f$ are fixed by \c
 * refcoeff_grad_div_m and \c refcoeff_grad_grad_rho, respectively.
 * The remaining linear portion of \f$\vec{u}\cdot\left(\mu+\lambda\right)
 * \vec{\nabla}\vec{\nabla}\cdot\vec{u}\f$ is
 * \f[
 *        \left\{\left(\mu+\lambda\right)\rho^{-2}\vec{m}\right\}_0
 *        \cdot\vec{\nabla}\vec{\nabla}\cdot\vec{m}
 *      - \operatorname{trace}\left[
 *          {\vec{\nabla}\vec{\nabla}\rho}^{\mathsf{T}}
 *          \left\{
 *            \left(\mu+\lambda\right)\rho^{-3}\vec{m}\otimes\vec{m}
 *          \right\}_0
 *        \right]
 * \f]
 *
 * @param mu \f$\mu\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param lambda \f$\lambda\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param rho \f$\rho\f$
 * @param grad_rho \f$\vec{\nabla}\rho\f$
 * @param grad_grad_rho \f$\vec{\nabla}\vec{\nabla}\rho\f$
 * @param m \f$\vec{m}\f$
 * @param div_m \f$\vec{\nabla}\cdot\vec{m}\f$
 * @param grad_m \f$\vec{\nabla}\vec{m}\f$
 * @param grad_div_m \f$\vec{\nabla}\vec{\nabla}\cdot\vec{m}\f$
 * @param refcoeff_grad_div_m the reference coefficient
 *        on \f$\vec{\nabla}\vec{\nabla}\cdot\vec{m}\f$ which may
 *        be computed using
 *        explicit_u_dot_mu_plus_lambda_grad_div_u_refcoeff_grad_div_m()
 * @param refcoeff_grad_grad_rho the reference coefficient
 *        on \f$\vec{\nabla}\vec{\nabla}\rho\f$ which may
 *        be computed using
 *        explicit_u_dot_mu_plus_lambda_grad_div_u_refcoeff_grad_grad_rho()
 *
 * @return The explicit portion of the sum of the viscosities times the
 *         velocity dotted against the gradient of the divergence of velocity.
 */
template<typename Scalar,
         typename Vector,
         typename Tensor,
         typename VectorCoefficient,
         typename TensorCoefficient >
inline
Scalar explicit_u_dot_mu_plus_lambda_grad_div_u(
        const Scalar            &mu,
        const Scalar            &lambda,
        const Scalar            &rho,
        const Vector            &grad_rho,
        const Tensor            &grad_grad_rho,
        const Vector            &m,
        const Scalar            &div_m,
        const Tensor            &grad_m,
        const Vector            &grad_div_m,
        const VectorCoefficient &refcoeff_grad_div_m,
        const TensorCoefficient &refcoeff_grad_grad_rho)
{
    const Vector coeff_grad_div_m(
            explicit_u_dot_mu_plus_lambda_grad_div_u_refcoeff_grad_div_m(
              mu, lambda, rho, m)
          - refcoeff_grad_div_m);
    const Tensor coeff_grad_grad_rho(
            explicit_u_dot_mu_plus_lambda_grad_div_u_refcoeff_grad_grad_rho<
                Scalar, Vector, Tensor
            >(mu, lambda, rho, m)
          - refcoeff_grad_grad_rho);

    return   ((mu+lambda)/(rho*rho*rho))*m.dot(
                  (2*grad_rho.dot(m)/rho - div_m)*grad_rho
                - grad_m.transpose()*grad_rho
             )
           + coeff_grad_div_m.dot(grad_div_m)
           - grad_grad_rho.cwiseProduct(coeff_grad_grad_rho).sum();
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}\cdot\vec{\nabla}e\f$ in the explicit portion
 * of \f$\mu\vec{\nabla}\cdot\vec{\nabla}T\f$.
 *
 * @param[in] mu \f$\mu\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] rho \f$\rho\f$
 *
 * @return \f$\mu\rho^{-1}\f$
 * @see explicit_mu_div_grad_T() for more details on the explicit
 *      operator.
 */
template<typename Scalar>
inline
Scalar explicit_mu_div_grad_T_refcoeff_div_grad_e(
        const Scalar &mu,
        const Scalar &rho)
{
    return mu/rho;
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}\cdot\vec{\nabla}\vec{m}\f$ in the explicit portion
 * of \f$\mu\vec{\nabla}\cdot\vec{\nabla}T\f$.
 *
 * @param[in] mu \f$\mu\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] rho \f$\rho\f$
 * @param[in] m \f$\vec{m}\f$
 *
 * @return \f$\mu\rho^{-2}\vec{m}\f$
 * @see explicit_mu_div_grad_T() for more details on the explicit
 *      operator.
 */
template<typename Scalar,
         typename Vector >
inline
Vector explicit_mu_div_grad_T_refcoeff_div_grad_m(
        const Scalar &mu,
        const Scalar &rho,
        const Vector &m)
{
    return (mu/(rho*rho))*m;
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}\cdot\vec{\nabla}\rho\f$ in the explicit portion
 * of \f$\mu\vec{\nabla}\cdot\vec{\nabla}T\f$.
 *
 * @param[in] gamma \f$\gamma\f$
 * @param[in] mu \f$\mu\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param[in] rho \f$\rho\f$
 * @param[in] e \f$e\f$
 * @param[in] p \f$p\f$
 *            computed from, for example, p_T_mu_lambda()
 *
 * @return \f$\mu\rho^{-2}\left((\gamma-1)e-2p\right)\f$
 * @see explicit_mu_div_grad_T() for more details on the explicit
 *      operator.
 */
template<typename Scalar>
inline
Scalar explicit_mu_div_grad_T_refcoeff_div_grad_rho(
        const Scalar &gamma,
        const Scalar &mu,
        const Scalar &rho,
        const Scalar &e,
        const Scalar &p)
{
    return mu/(rho*rho)*((gamma-1)*e-2*p);
}

/**
 * Compute the explicit portion of \f$\mu\vec{\nabla}\cdot\vec{\nabla}T\f$.
 * Uses the expansion
 * \f{align*}
 * \mu\vec{\nabla}\cdot\vec{\nabla}T =
 *   &{}- 2\gamma\mu\rho^{-2}\vec{\nabla}\rho\cdot
 *        \left(\vec{\nabla}p-\rho^{-1}p\vec{\nabla}\rho\right)
 * \\
 *   &{}- \gamma\left(\gamma-1\right)\mu\rho^{-2}\left[
 *              \operatorname{trace}\left(
 *                  {\vec{\nabla}\vec{m}}^{\mathsf{T}}\vec{\nabla}\vec{m}
 *              \right)
 *            - \rho^{-1}\left[
 *                2{\vec{\nabla}\vec{m}}^{\mathsf{T}}\vec{m}
 *                      \cdot\vec{\nabla}\rho
 *              - \rho^{-1} \vec{m}^2 \left(\vec{\nabla}\rho\right)^{2}
 *            \right]
 *        \right]
 * \\
 *   &{}+ \gamma\left(\gamma-1\right)\left(
 *              \mu\rho^{-1}
 *            - \left\{\mu\rho^{-1}\right\}_0
 *        \right)\vec{\nabla}\cdot\vec{\nabla}e
 *      - \gamma\left(\gamma-1\right)\left(
 *               \mu\rho^{-2}m
 *             - \left\{\mu\rho^{-2}\vec{m}\right\}_0
 *        \right)\cdot\vec{\nabla}\cdot\vec{\nabla}\vec{m}
 * \\
 *   &{}+ \gamma\left(
 *                   \mu\rho^{-2}\left((\gamma-1)e-2p\right)
 *          - \left\{\mu\rho^{-2}\left((\gamma-1)e-2p\right)\right\}_0
 *        \right)\vec{\nabla}\cdot\vec{\nabla}\rho
 * \f}
 * where \f$\left\{\mu\rho^{-1}\right\}_0\f$,
 * \f$\left\{\mu\rho^{-2}\vec{m}\right\}_0\f$, and
 * \f$\left\{\mu\rho^{-2}\left((\gamma-1)e-2p\right)\right\}_0\f$ are fixed by
 * \c refcoeff_div_grad_e, \c refcoeff_div_grad_m, and \c
 * refcoeff_div_grad_rho, respectively.  The remaining linear portion of
 * \f$\mu\vec{\nabla}\cdot\vec{\nabla}T\f$ is
 * \f[
 *      \gamma\left(\gamma-1\right)
 *      \left\{\mu\rho^{-1}\right\}_0
 *      \vec{\nabla}\cdot\vec{\nabla}e
 *    - \gamma\left(\gamma-1\right)
 *      \left\{\mu\rho^{-2}\vec{m}\right\}_0
 *      \cdot\vec{\nabla}\cdot\vec{\nabla}\vec{m}
 *    + \gamma
 *      \left\{\mu\rho^{-2}\left((\gamma-1)e-2p\right)\right\}_0
 *      \vec{\nabla}\cdot\vec{\nabla}\rho
 * \f]
 *
 * @param gamma \f$\gamma\f$
 * @param mu \f$\mu\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param rho \f$\rho\f$
 * @param grad_rho \f$\vec{\nabla}\rho\f$
 * @param div_grad_rho \f$\vec{\nabla}\cdot\vec{\nabla}\rho\f$
 * @param m \f$\vec{m}\f$
 * @param grad_m \f$\vec{\nabla}\vec{m}\f$
 * @param div_grad_m \f$\vec{\nabla}\cdot\vec{\nabla}\vec{m}\f$
 * @param e \f$e\f$
 * @param div_grad_e \f$\vec{\nabla}\cdot\vec{\nabla}e\f$
 * @param p \f$p\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param grad_p \f$\vec{\nabla}p\f$
 *            computed from, for example, p_T_mu_lambda()
 * @param refcoeff_div_grad_e the reference coefficient
 *        on \f$\vec{\nabla}\cdot\vec{\nabla}e\f$ which may
 *        be computed using
 *        explicit_mu_div_grad_T_refcoeff_div_grad_e()
 * @param refcoeff_div_grad_m the reference coefficient
 *        on \f$\vec{\nabla}\cdot\vec{\nabla}\vec{m}\f$ which may
 *        be computed using
 *        explicit_mu_div_grad_T_refcoeff_div_grad_m()
 * @param refcoeff_div_grad_rho the reference coefficient
 *        on \f$\vec{\nabla}\cdot\vec{\nabla}\rho\f$ which may
 *        be computed using
 *        explicit_mu_div_grad_T_refcoeff_div_grad_rho()
 *
 * @return The explicit portion of the viscosity times
 *         the Laplacian of temperature.
 */
template<typename Scalar,
         typename Vector,
         typename Tensor,
         typename ScalarCoefficient1,
         typename ScalarCoefficient2,
         typename VectorCoefficient >
inline
Scalar explicit_mu_div_grad_T(
        const Scalar             &gamma,
        const Scalar             &mu,
        const Scalar             &rho,
        const Vector             &grad_rho,
        const Scalar             &div_grad_rho,
        const Vector             &m,
        const Tensor             &grad_m,
        const Vector             &div_grad_m,
        const Scalar             &e,
        const Scalar             &div_grad_e,
        const Scalar             &p,
        const Vector             &grad_p,
        const ScalarCoefficient1 &refcoeff_div_grad_e,
        const VectorCoefficient  &refcoeff_div_grad_m,
        const ScalarCoefficient2 &refcoeff_div_grad_rho)
{
    const Scalar rho_inverse  = 1/rho;
    const Scalar rho_inverse2 = rho_inverse*rho_inverse;
    const Scalar coeff_div_grad_e(
            explicit_mu_div_grad_T_refcoeff_div_grad_e(mu, rho)
          - refcoeff_div_grad_e);
    const Vector coeff_div_grad_m(
            explicit_mu_div_grad_T_refcoeff_div_grad_m(mu, rho, m)
          - refcoeff_div_grad_m);
    const Scalar coeff_div_grad_rho(
            explicit_mu_div_grad_T_refcoeff_div_grad_rho(gamma, mu, rho, e, p)
          - refcoeff_div_grad_rho);

    return gamma*(
                coeff_div_grad_rho*div_grad_rho
              - 2*mu*rho_inverse2*grad_rho.dot(
                    grad_p - rho_inverse*p*grad_rho
                )
              + (gamma-1)*(
                  - mu*rho_inverse2*(
                        grad_m.squaredNorm() // Frobenius norm
                      - rho_inverse*(
                            2*grad_rho.dot(grad_m.transpose()*m)
                          - rho_inverse*m.squaredNorm()*grad_rho.squaredNorm()
                        )
                    )
                  + coeff_div_grad_e*div_grad_e
                  - div_grad_m.dot(coeff_div_grad_m)
                )
           );
}

/**
 * Compute \f$p\f$ using the equation of state.
 *
 * @param[in]  alpha \f$\alpha\f$
 * @param[in]  beta \f$\beta\f$
 * @param[in]  gamma \f$\gamma\f$
 * @param[in]  rho \f$\rho\f$
 * @param[in]  m \f$\vec{m}\f$
 * @param[in]  e \f$e\f$
 * @param[out] p \f$p\f$
 */
template<typename Scalar,
         typename Vector  >
inline
void p(const Scalar &alpha,
       const Scalar &beta,
       const Scalar &gamma,
       const Scalar &rho,
       const Vector &m,
       const Scalar &e,
       Scalar &p)
{
    SUZERAIN_UNUSED(alpha);
    SUZERAIN_UNUSED(beta);

    // Compute scalar quantities
    p = (gamma - 1)*(e - m.squaredNorm()/2/rho);
}

/**
 * Compute \f$p\f$ using the equation of state.
 *
 * @param[in]  gamma \f$\gamma\f$
 * @param[in]  rho   \f$\rho\f$
 * @param[in]  T     \f$T\f$
 * @param[out] p     \f$p\f$
 */
template<typename Scalar>
inline
void p(const Scalar &gamma,
       const Scalar &rho,
       const Scalar &T,
       Scalar &p)
{
    p = rho*T/gamma;
}

/**
 * Compute \f$p\f$ and \f$T\f$ using the equation of state.
 *
 * @param[in]  alpha \f$\alpha\f$
 * @param[in]  beta \f$\beta\f$
 * @param[in]  gamma \f$\gamma\f$
 * @param[in]  rho \f$\rho\f$
 * @param[in]  m \f$\vec{m}\f$
 * @param[in]  e \f$e\f$
 * @param[out] p \f$p\f$
 * @param[out] T \f$T\f$
 */
template<typename Scalar,
         typename Vector  >
inline
void p_T(
        const Scalar &alpha,
        const Scalar &beta,
        const Scalar &gamma,
        const Scalar &rho,
        const Vector &m,
        const Scalar &e,
        Scalar &p,
        Scalar &T)
{
    SUZERAIN_UNUSED(alpha);
    SUZERAIN_UNUSED(beta);

    const Scalar rho_inverse = 1/rho;

    // Compute scalar quantities
    p = (gamma - 1)*(e - rho_inverse*m.squaredNorm()/2);
    T = gamma * p * rho_inverse;
}

/**
 * Compute \f$p\f$, \f$T\f$, and their gradients using the equation of state.
 * The gradients are computed using these expansions:
 * \f{align*}
 *      \vec{\nabla}p &= (\gamma-1)\left[
 *            \vec{\nabla}e
 *          + \frac{1}{2}\rho^{-2}\left(\vec{m}\cdot\vec{m}\right)
 *            \vec{\nabla}\rho
 *          - \rho^{-1} \left(\vec{\nabla}\vec{m}\right)^{\mathsf{T}}\vec{m}
 *      \right]
 *      \\
 *      \vec{\nabla}T &= \gamma\rho^{-1}\vec{\nabla}p
 *                     - \gamma\rho^{-2} p \vec{\nabla}\rho
 * \f}
 *
 * @param[in]  alpha \f$\alpha\f$
 * @param[in]  beta \f$\beta\f$
 * @param[in]  gamma \f$\gamma\f$
 * @param[in]  rho \f$\rho\f$
 * @param[in]  grad_rho \f$\vec{\nabla}\rho\f$
 * @param[in]  m \f$\vec{m}\f$
 * @param[in]  grad_m \f$\vec{\nabla}\vec{m}\f$
 * @param[in]  e \f$e\f$
 * @param[in]  grad_e \f$\vec{\nabla}e\f$
 * @param[out] p \f$p\f$
 * @param[out] grad_p \f$\vec{\nabla}p\f$
 * @param[out] T \f$T\f$
 * @param[out] grad_T \f$\vec{\nabla}T\f$
 */
template<typename Scalar,
         typename Vector,
         typename Tensor >
inline
void p_T(
        const Scalar &alpha,
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
        Vector &grad_T)
{
    SUZERAIN_UNUSED(alpha);  // Present for API consistency
    SUZERAIN_UNUSED(beta);   // Present for API consistency

    const Scalar rho_inverse      = 1/rho;
    const Scalar half_rho_inverse = rho_inverse/2;
    const Scalar gamma1           = gamma - 1;

    // Compute scalar quantities
    p = gamma1*(e - half_rho_inverse*m.squaredNorm());
    T = gamma * p * rho_inverse;

    // Compute vector quantities
    grad_p = gamma1*(
                grad_e + rho_inverse*(
                      (half_rho_inverse*m.squaredNorm())*grad_rho
                    - grad_m.transpose()*m
                )
             );
    grad_T = gamma*rho_inverse*(grad_p - rho_inverse*p*grad_rho);
}

/**
 * Compute \f$p\f$, \f$T\f$, \f$\mu\f$, and \f$\lambda\f$
 * using the equation of state.
 *
 * @param[in]  alpha \f$\alpha\f$
 * @param[in]  beta \f$\beta\f$
 * @param[in]  gamma \f$\gamma\f$
 * @param[in]  rho \f$\rho\f$
 * @param[in]  m \f$\vec{m}\f$
 * @param[in]  e \f$e\f$
 * @param[out] p \f$p\f$
 * @param[out] T \f$T\f$
 * @param[out] mu \f$\mu\f$
 * @param[out] lambda \f$\lambda\f$
 */
template<typename Scalar,
         typename Vector  >
inline
void p_T_mu_lambda(
        const Scalar &alpha,
        const Scalar &beta,
        const Scalar &gamma,
        const Scalar &rho,
        const Vector &m,
        const Scalar &e,
        Scalar &p,
        Scalar &T,
        Scalar &mu,
        Scalar &lambda)
{
    using std::pow;
    const Scalar rho_inverse = 1/rho;

    // Compute scalar quantities
    p      = (gamma - 1)*(e - rho_inverse*m.squaredNorm()/2);
    T      = gamma * p * rho_inverse;
    mu     = pow(T, beta);
    lambda = (alpha - Scalar(2)/3)*mu;
}

/**
 * Compute \f$p\f$, \f$T\f$, \f$\mu\f$, and \f$\lambda\f$ and their gradients
 * using the equation of state.  The gradients are computed using these
 * expansions:
 * \f{align*}
 *      \vec{\nabla}p &= (\gamma-1)\left[
 *            \vec{\nabla}e
 *          + \frac{1}{2}\rho^{-2}\left(\vec{m}\cdot\vec{m}\right)
 *            \vec{\nabla}\rho
 *          - \rho^{-1} \left(\vec{\nabla}\vec{m}\right)^{\mathsf{T}}\vec{m}
 *      \right]
 *      \\
 *      \vec{\nabla}T &= \gamma\rho^{-1}\vec{\nabla}p
 *                     - \gamma\rho^{-2} p \vec{\nabla}\rho
 *      \\
 *      \vec{\nabla}\mu &= \beta{}T^{\beta-1}\vec{\nabla}T
 *      \\
 *      \vec{\nabla}\lambda &= \left(\alpha-\frac{2}{3}\right)\vec{\nabla}\mu
 * \f}
 *
 * @param[in]  alpha \f$\alpha\f$
 * @param[in]  beta \f$\beta\f$
 * @param[in]  gamma \f$\gamma\f$
 * @param[in]  rho \f$\rho\f$
 * @param[in]  grad_rho \f$\vec{\nabla}\rho\f$
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
         typename Vector,
         typename Tensor >
inline
void p_T_mu_lambda(
        const Scalar &alpha,
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
    using std::pow;
    const Scalar rho_inverse      = 1/rho;
    const Scalar half_rho_inverse = rho_inverse/2;
    const Scalar gamma1           = gamma - 1;
    const Scalar alpha23          = alpha - Scalar(2)/3;

    // Compute scalar quantities
    p      = gamma1*(e - half_rho_inverse*m.squaredNorm());
    T      = gamma * p * rho_inverse;

    const Scalar T_to_beta1 = pow(T, beta - 1); // Avoid div for grad_mu

    mu     = T_to_beta1 * T;
    lambda = alpha23*mu;

    // Compute vector quantities
    grad_p = gamma1*(
                grad_e + rho_inverse*(
                      (half_rho_inverse*m.squaredNorm())*grad_rho
                    - grad_m.transpose()*m
                )
             );
    grad_T      = gamma*rho_inverse*(grad_p - rho_inverse*p*grad_rho);
    grad_mu     = (beta * T_to_beta1)*grad_T;
    grad_lambda = alpha23*grad_mu;
}

/**
 * Compute the internal energy by subtracting the kinetic
 * energy from the total energy.
 *
 * @param[in] rho \f$\rho\f$
 * @param[in] m   \f$\vec{m}\f$
 * @param[in] e   \f$e\f$
 */
template<typename Scalar,
         typename Vector  >
inline
Scalar energy_internal(
        const Scalar &rho,
        const Vector &m,
        const Scalar &e)
{
    return e - m.squaredNorm()/rho/2;
}

/**
 * Compute the internal energy using the pressure.
 *
 * @param[in] gamma \f$\gamma\f$
 * @param[in] p     \f$p\f$
 */
template<typename Scalar>
inline
Scalar energy_internal(
        const Scalar &gamma,
        const Scalar &p)
{
    return p / (gamma - 1);
}

/**
 * Compute the gradient of the internal energy from the pressure gradient.
 *
 * @param[in] gamma  \f$\gamma\f$
 * @param[in] grad_p \f$\vec{\nabla}p\f$
 */
template<typename Scalar,
         typename Vector>
inline
Vector energy_internal_gradient(
        const Scalar &gamma,
        const Vector &grad_p)
{
    return grad_p / (gamma - 1);
}

/**
 * Compute the kinetic energy.
 *
 * @param[in] rho \f$\rho\f$
 * @param[in] m   \f$\vec{m}\f$
 */
template<typename Scalar,
         typename Vector  >
inline
Scalar energy_kinetic(
        const Scalar &rho,
        const Vector &m)
{
    return m.squaredNorm()/rho/2;
}

/**
 * Compute the gradient of the kinetic energy.
 *
 * @param[in] rho      \f$\rho\f$
 * @param[in] grad_rho \f$\vec{\nabla}\rho\f$
 * @param[in] m        \f$\vec{m}\f$
 * @param[in] grad_m   \f$\vec{\nabla}\vec{m}\f$
 */
template<typename Scalar,
         typename Vector,
         typename Tensor  >
inline
Vector energy_kinetic_gradient(
        const Scalar &rho,
        const Vector &grad_rho,
        const Vector &m,
        const Tensor &grad_m)
{
    const Scalar rho_inverse  = 1/rho;
    return rho_inverse*(
               grad_m.transpose()*m
             - rho_inverse/2*m.squaredNorm()*grad_rho
           );
}

/**
 * Compute \f$\vec{\nabla}\cdot\vec\nabla{}p\f$
 * using the equation of state.  Uses the expansion
 * \f{align*}
 *      \vec{\nabla}\cdot\vec{\nabla}p =
 *      \left(\gamma-1\right)\left[
 *          \vec{\nabla}\cdot\vec{\nabla}e
 *          - \rho^{-1}\left[
 *                \operatorname{trace}\left(
 *                    \vec{\nabla}\vec{m}^{\mathsf{T}}\vec{\nabla}\vec{m}
 *                \right)
 *              + \left(\vec{\nabla}\cdot\vec{\nabla}\vec{m}\right)\cdot\vec{m}
 *              - \rho^{-1}\left[
 *                    2\vec{\nabla}\vec{m}^{\mathsf{T}}\vec{m}
 *                        \cdot\vec{\nabla}\rho
 *                  + \frac{1}{2}\left(\vec{m}\cdot\vec{m}\right)
 *                        \vec{\nabla}\cdot\vec{\nabla}\rho
 *                  - \rho^{-1}\left(\vec{m}\cdot\vec{m}\right)
 *                        \left(\vec{\nabla}\rho\cdot\vec{\nabla}\rho\right)
 *              \right]
 *          \right]
 *      \right]
 *      .
 * \f}
 * @param[in]  gamma \f$\gamma\f$
 * @param[in]  rho \f$\rho\f$
 * @param[in]  grad_rho \f$\vec{\nabla}\rho\f$
 * @param[in]  div_grad_rho \f$\vec{\nabla}\cdot\vec{\nabla}\rho\f$
 * @param[in]  m \f$\vec{m}\f$
 * @param[in]  grad_m \f$\vec{\nabla}\vec{m}\f$
 * @param[in]  div_grad_m \f$\vec{\nabla}\cdot\vec{\nabla}\vec{m}\f$
 * @param[in]  e \f$e\f$
 * @param[in]  grad_e \f$\vec{\nabla}e\f$
 * @param[in]  div_grad_e \f$\vec{\nabla}\cdot\vec{\nabla}e\f$
 *
 * @return The Laplacian of pressure.
 */
template<typename Scalar,
         typename Vector,
         typename Tensor >
inline
Scalar div_grad_p(
        const Scalar &gamma,
        const Scalar &rho,
        const Vector &grad_rho,
        const Scalar &div_grad_rho,
        const Vector &m,
        const Tensor &grad_m,
        const Vector &div_grad_m,
        const Scalar &e,
        const Vector &grad_e,
        const Scalar &div_grad_e)
{
    SUZERAIN_UNUSED(e);
    SUZERAIN_UNUSED(grad_e);

    const Scalar rho_inverse = 1/rho;
    const Scalar m_dot_m     = m.squaredNorm();

    return (gamma - 1)*(
                  div_grad_e
                - rho_inverse*(
                          grad_m.squaredNorm() // Frobenius norm
                        + div_grad_m.dot(m)
                        - rho_inverse*(
                                2*grad_rho.dot(grad_m.transpose()*m)
                              + m_dot_m*div_grad_rho/2
                              - rho_inverse*m_dot_m*grad_rho.squaredNorm()
                            )
                    )
            );
}

/**
 * Compute \f$\vec{\nabla}\cdot\vec\nabla{}T\f$
 * using the equation of state.  Uses the expansion
 * \f[
 *      \vec{\nabla}\cdot\vec\nabla{}T =
 *          \gamma\rho^{-1}\left[
 *                \vec{\nabla}\cdot\vec{\nabla}p
 *              - \rho^{-1}\left[
 *                    p\vec{\nabla}\cdot\vec{\nabla}\rho
 *                  + 2\vec{\nabla}\rho\cdot\left(
 *                      \vec{\nabla}p - \rho^{-1}p\vec{\nabla}\rho
 *                  \right)
 *              \right]
 *          \right]
 *      .
 * \f]
 * @param[in]  gamma \f$\gamma\f$
 * @param[in]  rho \f$\rho\f$
 * @param[in]  grad_rho \f$\vec{\nabla}\rho\f$
 * @param[in]  div_grad_rho \f$\vec{\nabla}\cdot\vec{\nabla}\rho\f$
 * @param[in]  p \f$p\f$
 *             computed from, for example, p_T_mu_lambda()
 * @param[in]  grad_p \f$\vec{\nabla}p\f$
 *             computed from, for example, p_T_mu_lambda()
 * @param[in]  div_grad_p \f$\vec{\nabla}\cdot\vec{\nabla}p\f$
 *             computed from, for example, div_grad_p()
 *
 * @return The Laplacian of temperature.
 */
template<typename Scalar,
         typename Vector >
inline
Scalar div_grad_T(
        const Scalar &gamma,
        const Scalar &rho,
        const Vector &grad_rho,
        const Scalar &div_grad_rho,
        const Scalar &p,
        const Vector &grad_p,
        const Scalar &div_grad_p)
{
    const Scalar rho_inverse = 1/rho;

    return gamma*rho_inverse*(
                  div_grad_p
                - rho_inverse*(
                          p*div_grad_rho
                        + 2*grad_rho.dot(grad_p - rho_inverse*p*grad_rho)
                    )
            );
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
         typename Vector >
inline
Vector u(const Scalar &rho,
         const Vector &m)
{
    const Scalar rho_inverse = 1/rho;

    return rho_inverse*m;
}

/**
 * Compute \f$\vec{\nabla}\vec{u}\f$.  Uses the expansion
 * \f[
 *      \vec{\nabla}\vec{u} =
 *      \rho^{-1} \left[
 *              \vec{\nabla}\vec{m}
 *            - \rho^{-1} \vec{m}\otimes\vec{\nabla}\rho
 *      \right]
 *      .
 * \f]
 *
 * @param[in]  rho \f$\rho\f$
 * @param[in]  grad_rho \f$\vec{\nabla}\rho\f$
 * @param[in]  m \f$\vec{m}\f$
 * @param[in]  grad_m \f$\vec{\nabla}\vec{m}\f$
 *
 * @return The gradient of the velocity based on the provided fields.
 */
template<typename Scalar,
         typename Vector,
         typename Tensor >
inline
Tensor grad_u(
        const Scalar &rho,
        const Vector &grad_rho,
        const Vector &m,
        const Tensor &grad_m)
{
    // Expression mangled to avoid GCC 4.3 per Redmine issue #1562
    const Scalar rho_inverse = 1/rho;
    Tensor retval(-rho_inverse*m*grad_rho.transpose());
    retval += grad_m;
    retval *= rho_inverse;
    return retval;
}

/**
 * Compute \f$\vec{\nabla}\cdot\vec{u}\f$.  Uses the expansion
 * \f[
 *      \vec{\nabla}\cdot\vec{u} =
 *      \rho^{-1}\left[
 *          \vec{\nabla}\cdot\vec{m} - \rho^{-1}\vec{m}\cdot\vec{\nabla}\rho
 *      \right]
 *      .
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
         typename Vector >
inline
Scalar div_u(
        const Scalar &rho,
        const Vector &grad_rho,
        const Vector &m,
        const Scalar &div_m)
{
    const Scalar rho_inverse = 1/rho;

    return rho_inverse*(div_m - rho_inverse*grad_rho.dot(m));
}

/**
 * Compute \f$\vec{\nabla}\vec{\nabla}\cdot{}\vec{u}\f$.  Uses the expansion
 * \f[
 * \vec{\nabla}\vec{\nabla}\cdot{}\vec{u} =
 *      \rho^{-1}\left[
 *            \vec{\nabla}\vec{\nabla}\cdot\vec{m}
 *          + \rho^{-1}\left[
 *                \left(
 *                     2 \rho^{-1} \vec{\nabla}\rho\cdot\vec{m}
 *                   - \vec{\nabla}\cdot\vec{m}
 *                \right)\vec{\nabla}\rho
 *              - \left(\vec{\nabla}\vec{\nabla}\rho\right)\vec{m}
 *              - \vec{\nabla}\vec{m}^{\mathsf{T}}\vec{\nabla}\rho
 *            \right]
 *      \right]
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
         typename Vector,
         typename Tensor >
inline
Vector grad_div_u(
        const Scalar &rho,
        const Vector &grad_rho,
        const Tensor &grad_grad_rho,
        const Vector &m,
        const Scalar &div_m,
        const Tensor &grad_m,
        const Vector &grad_div_m)
{
    const Scalar rho_inverse = 1/rho;

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
 *     \rho^{-1}\left[
 *        \vec{\nabla}\cdot\vec{\nabla}\vec{m}
 *        +
 *        \rho^{-1}\left[
 *            \left(
 *                2\rho^{-1}\left(\vec{\nabla}\rho\cdot\vec{\nabla}\rho\right)
 *              - \vec{\nabla}\cdot\vec{\nabla}\rho
 *            \right) \vec{m}
 *          - 2 \left(\vec{\nabla}\vec{m}\right)\vec{\nabla}\rho
 *        \right]
 *     \right]
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
         typename Vector,
         typename Tensor >
inline
Vector div_grad_u(
        const Scalar &rho,
        const Vector &grad_rho,
        const Scalar &div_grad_rho,
        const Vector &m,
        const Tensor &grad_m,
        const Vector &div_grad_m)
{
    const Scalar rho_inverse = 1/rho;

    return rho_inverse*(
                div_grad_m + rho_inverse*(
                    (2*rho_inverse*grad_rho.squaredNorm() - div_grad_rho)*m
                   - 2*grad_m*grad_rho
                )
            );
}

/**
 * Compute \f$\vec{\nabla}\times\vec{\nabla}\times\vec{u}\f$.  Uses the
 * expansion
 * \f[
 * \vec{\nabla}\times\vec{\nabla}\times\vec{u} =
 *      \rho^{-1}\left[
 *            \vec{\nabla}\times\vec{\nabla}\times\vec{m}
 *          + \rho^{-1} \left[
 *                \left(
 *                      2\vec{\nabla}\vec{m}
 *                    - \vec{\nabla}{\vec{m}}^{\mathsf{T}}
 *                \right) \vec{\nabla}\rho
 *              + \left(
 *                      \vec{\nabla}\cdot\vec{\nabla}\rho
 *                    - 2 \rho^{-1} \left(\vec{\nabla}\rho\right)^2
 *                \right) \vec{m}
 *              - \left( \vec{\nabla}\vec{\nabla}\rho \right) \vec{m}
 *              + \left(
 *                      2 \rho^{-1} \vec{\nabla}\rho \cdot \vec{m}
 *                    - \vec{\nabla}\cdot\vec{m}
 *                \right) \vec{\nabla}\rho
 *          \right]
 *      \right]
 * \f]
 *
 * @param[in]  rho \f$\rho\f$
 * @param[in]  grad_rho \f$\vec{\nabla}\rho\f$.
 * @param[in]  grad_grad_rho \f$\vec{\nabla}\vec{\nabla}\rho\f$.
 * @param[in]  m \f$\vec{m}\f$
 * @param[in]  div_m \f$\vec{\nabla}\cdot\vec{m}\f$
 * @param[in]  grad_m \f$\vec{\nabla}\vec{m}\f$
 * @param[in]  curl_curl_m \f$\vec{\nabla}\times\vec{\nabla}\times\vec{m}\f$
 *
 * @return The curl of the curl of the velocity based on the provided fields.
 */
template<typename Scalar,
         typename Vector,
         typename Tensor >
inline
Vector curl_curl_u(
        const Scalar &rho,
        const Vector &grad_rho,
        const Tensor &grad_grad_rho,
        const Vector &m,
        const Scalar &div_m,
        const Tensor &grad_m,
        const Vector &curl_curl_m)
{
    const Scalar rho_inverse = 1/rho;

    return rho_inverse*(
                 curl_curl_m
               + rho_inverse*(
                     (   2*grad_m
                       - grad_m.transpose()
                     ) * grad_rho
                   + (   grad_grad_rho.trace()
                       - 2*rho_inverse*grad_rho.squaredNorm()
                     ) * m
                   - grad_grad_rho * m
                   + (   2*rho_inverse*grad_rho.dot(m)
                       - div_m
                     ) * grad_rho
               )
           );
}

} // namespace rholt

} // namespace suzerain

#endif // SUZERAIN_RHOLT_HPP
