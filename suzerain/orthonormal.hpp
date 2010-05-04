/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
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
#ifndef __SUZERAIN_ORTHONORMAL_H
#define __SUZERAIN_ORTHONORMAL_H

#include <suzerain/common.hpp>

/** @file
 * Provides routines that compute classical state quantities (based on
 * nondimensional \f$p\f$, \f$T\f$, \f$\vec{u}\f$) from other state quantities.
 */

namespace suzerain
{

/**
 * Provides routines that compute classical quantities (based on
 * nondimensional \f$p\f$, \f$T\f$, \f$\vec{u}\f$, \f$\mu\f$, \f$\lambda\f$)
 * from other conserved and classical quantities under the assumption of
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
 *      + \left(\vec{\nabla}\vec{u}\right)^{\mathsf{T}}
 *    \right)
 *  + \lambda \left( \vec{\nabla}\cdot\vec{u} \right)
 *            \accentset{\leftrightarrow}{I}\f$.
 *
 * @param[in] mu \f$\mu\f$
 *            computed from, for example, rhome::p_T_mu_lambda()
 * @param[in] lambda \f$\lambda\f$
 *            computed from, for example, rhome::p_T_mu_lambda()
 * @param[in] div_u \f$\vec{\nabla}\cdot\vec{u}\f$
 *            computed from, for example, rhome::div_u()
 * @param[in] grad_u \f$\vec{\nabla}\vec{u}\f$
 *            computed from, for example, rhome::grad_u()
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
         typename Vector = Eigen::Matrix<Scalar,3,1>,
         typename Tensor = Eigen::Matrix<Scalar,3,3> >
Vector div_rho_inverse_m_outer_m(
        const Scalar &rho,
        const Vector &grad_rho,
        const Vector &m,
        const Scalar &div_m,
        const Tensor &grad_m)
{
    const Scalar rho_inverse = 1.0/rho;

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
 *            computed from, for example, rhome::u()
 * @param[in] div_u \f$\vec{\nabla}\cdot\vec{u}\f$
 *            computed from, for example, rhome::div_u()
 *
 * @return The convective derivative computed from the divergence form.
 */
template<typename Scalar,
         typename Vector = Eigen::Matrix<Scalar,3,1>,
         typename Tensor = Eigen::Matrix<Scalar,3,3> >
Vector div_u_outer_m(
        const Vector &m,
        const Tensor &grad_m,
        const Vector &u,
        const Scalar &div_u)
{
    return grad_m*u + div_u*m;
}

/**
 * Compute \f$\vec{\nabla}\cdot\accentset{\leftrightarrow}{\tau}\f$.
 * Uses the expansion
 * \f[
 *      \vec{\nabla}\cdot\accentset{\leftrightarrow}{\tau} =
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
 *            computed from, for example, rhome::p_T_mu_lambda()
 * @param[in] grad_mu \f$\vec{\nabla}\mu\f$
 *            computed from, for example, rhome::p_T_mu_lambda()
 * @param[in] lambda \f$\lambda\f$
 *            computed from, for example, rhome::p_T_mu_lambda()
 * @param[in] grad_lambda \f$\vec{\nabla}\lambda\f$
 *            computed from, for example, rhome::p_T_mu_lambda()
 * @param[in] div_u \f$\vec{\nabla}\cdot\vec{u}\f$
 *            computed from, for example, rhome::div_u()
 * @param[in] grad_u \f$\vec{\nabla}\vec{u}\f$
 *            computed from, for example, rhome::grad_u()
 * @param[in] div_grad_u \f$\vec{\nabla}\cdot\vec{\nabla}\vec{u}\f$
 *            computed from, for example, rhome::div_grad_u()
 * @param[in] grad_div_u \f$\vec{\nabla}\vec{\nabla}\cdot\vec{u}\f$
 *            computed from, for example, rhome::grad_div_u()
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
 *            computed from, for example, rhome::u()
 * @param[in] div_u \f$\vec{\nabla}\cdot\vec{u}\f$
 *            computed from, for example, rhome::div_u()
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
 *            computed from, for example, rhome::p_T_mu_lambda()
 * @param[in] div_grad_T \f$\vec{\nabla}\cdot\vec{\nabla}T\f$
 *            computed from, for example, rhome::p_T_mu_lambda()
 * @param[in] mu \f$\mu\f$
 *            computed from, for example, rhome::p_T_mu_lambda()
 * @param[in] grad_mu \f$\vec{\nabla}\mu\f$
 *            computed from, for example, rhome::p_T_mu_lambda()
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
 *            computed from, for example, rhome::p_T_mu_lambda()
 * @param[in] grad_p \f$\vec{\nabla}p\f$
 *            computed from, for example, rhome::p_T_mu_lambda()
 * @param[in] u \f$\vec{u}\f$
 *            computed from, for example, rhome::u()
 * @param[in] div_u \f$\vec{\nabla}\cdot\vec{u}\f$
 *            computed from, for example, rhome::div_u()
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
 *      \operatorname{trace}\left(
 *          \accentset{\leftrightarrow}{\tau}\,\vec{\nabla}\vec{u}
 *      \right)
 * \f]
 *
 * @param[in] u \f$\vec{u}\f$
 *            computed from, for example, rhome::u()
 * @param[in] grad_u \f$\vec{\nabla}\vec{u}\f$
 *            computed from, for example, rhome::grad_u()
 * @param[in] tau \f$\accentset{\leftrightarrow}{\tau}\f$
 *            computed from, for example, tau()
 * @param[in] div_tau \f$\vec{\nabla}\cdot\accentset{\leftrightarrow}{\tau}\f$
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
         typename Vector = Eigen::Matrix<Scalar,3,1> >
inline
Vector explicit_div_p_u_refcoeff_grad_rho(
        const Scalar &rho,
        const Vector &m,
        const Scalar &p)
{
    const Scalar rho_inverse = 1.0/rho;
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
 *            computed from, for example, rhome::p_T_mu_lambda()
 * @param grad_p \f$\vec{\nabla}p\f$
 *            computed from, for example, rhome::p_T_mu_lambda()
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
         typename Vector            = Eigen::Matrix<Scalar,3,1>,
         typename ScalarCoefficient = Scalar,
         typename VectorCoefficient = Vector >
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
 * \f$\vec{\nabla}\cdot\vec{m}\f$ in the explicit portion
 * of \f$\vec{\nabla}\cdot\left(e+p\right)\vec{u}\f$.
 *
 * @param[in] gamma \f$\gamma\f$
 * @param[in] rho \f$\rho\f$
 * @param[in] m \f$\vec{m}\f$
 * @param[in] e \f$e\f$
 *
 * @return \f$\rho^{-1}\left( \gamma{}e
 *              - \frac{\gamma-1}{2}\rho^{-1}\vec{m}^2\right)\f$
 * @see explicit_div_e_plus_p_u for more details on the explicit
 *      operator.
 */
template<typename Scalar,
         typename Vector = Eigen::Matrix<Scalar,3,1> >
inline
Scalar explicit_div_e_plus_p_u_refcoeff_div_m(
        const Scalar &gamma,
        const Scalar &rho,
        const Vector &m,
        const Scalar &e)
{
    const Scalar rho_inverse = 1.0/rho;
    return rho_inverse*(gamma*e - (gamma-1)/2*rho_inverse*m.squaredNorm());
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}\rho\f$ in the explicit portion
 * of \f$\vec{\nabla}\cdot\left(e+p\right)\vec{u}\f$.
 *
 * @param[in] gamma \f$\gamma\f$
 * @param[in] rho \f$\rho\f$
 * @param[in] m \f$\vec{m}\f$
 * @param[in] e \f$e\f$
 *
 * @return \f$\rho^{-2}\left(
 *              \left(\gamma-1\right) \rho^{-1}\vec{m}^2
 *            - \gamma{}e
 *         \right)\vec{m}\f$
 * @see explicit_div_e_plus_p_u for more details on the explicit
 *      operator.
 */
template<typename Scalar,
         typename Vector = Eigen::Matrix<Scalar,3,1> >
inline
Vector explicit_div_e_plus_p_u_refcoeff_grad_rho(
        const Scalar &gamma,
        const Scalar &rho,
        const Vector &m,
        const Scalar &e)
{
    const Scalar rho_inverse = 1.0/rho;
    return rho_inverse*rho_inverse*(
              (gamma-1)*rho_inverse*m.squaredNorm()
            - gamma*e
           )*m;
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}e\f$ in the explicit portion
 * of \f$\vec{\nabla}\cdot\left(e+p\right)\vec{u}\f$.
 *
 * @param[in] gamma \f$\gamma\f$
 * @param[in] rho \f$\rho\f$
 * @param[in] m \f$\vec{m}\f$
 *
 * @return \f$\gamma\rho^{-1}\vec{m}\f$
 * @see explicit_div_e_plus_p_u for more details on the explicit
 *      operator.
 */
template<typename Scalar,
         typename Vector = Eigen::Matrix<Scalar,3,1> >
inline
Vector explicit_div_e_plus_p_u_refcoeff_grad_e(
        const Scalar &gamma,
        const Scalar &rho,
        const Vector &m)
{
    return gamma/rho*m;
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
 *          \rho^{-1}\left(\gamma{}e
 *             -\frac{\gamma-1}{2}\rho^{-1}\vec{m}^2\right)
 *        - \left\{\rho^{-1}\left(\gamma{}e
 *             -\frac{\gamma-1}{2}\rho^{-1}\vec{m}^2\right)\right\}_0
 *     \right)\vec{\nabla}\cdot\vec{m}
 * \\
 *  &+ \left(
 *          \rho^{-2}\left(\left(\gamma-1\right)\rho^{-1}\vec{m}^2
 *                  -\gamma{}e\right)\vec{m}
 *        - \left\{\rho^{-2}\left(\left(\gamma-1\right)\rho^{-1}\vec{m}^2
 *                  -\gamma{}e\right)\vec{m}\right\}_0
 *     \right)\cdot\vec{\nabla}\rho
 * \\
 *  &+ \left(
 *          \gamma\rho^{-1}\vec{m}
 *        - \left\{\gamma\rho^{-1}\vec{m}\right\}_0
 *     \right)\cdot\vec{\nabla}e
 * \f}
 * where
 * \f$\left\{\rho^{-1}\left(\gamma{}e
 * -\frac{\gamma-1}{2}\rho^{-1}\vec{m}^2\right)\right\}_0\f$
 *  \f$\left\{\rho^{-2}\left(\left(\gamma-1\right)\rho^{-1}\vec{m}^2
 *  -\gamma{}e\right)\vec{m}\right\}_0\f$, and
 *  \f$\left\{\gamma\rho^{-1}\vec{m}\right\}_0\f$ are fixed by \c
 *  refcoeff_div_m, \c refcoeff_grad_rho, and
 * \c refcoeff_grad_e, respectively.
 * The remaining linear portion of
 * \f$\vec{\nabla}\cdot\left(e+p\right)\vec{u}\f$ is
 * \f[
 *      \left\{\rho^{-1}\left(\gamma{}e
 *        -\frac{\gamma-1}{2}\rho^{-1}\vec{m}^2\right)\right\}_0
 *        \vec{\nabla}\cdot\vec{m}
 *    + \left\{\rho^{-2}\left(\left(\gamma-1\right)\rho^{-1}\vec{m}^2
 *        -\gamma{}e\right)\vec{m}\right\}_0
 *        \cdot\vec{\nabla}\rho
 *    + \left\{\gamma\rho^{-1}\vec{m}\right\}_0
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
         typename Vector             = Eigen::Matrix<Scalar,3,1>,
         typename Tensor             = Eigen::Matrix<Scalar,3,3>,
         typename ScalarCoefficient  = Scalar,
         typename VectorCoefficient1 = Vector,
         typename VectorCoefficient2 = Vector >
Scalar explicit_div_e_plus_p_u(
        const Scalar             &gamma,
        const Scalar             &rho,
        const Vector             &grad_rho,
        const Vector             &m,
        const Scalar             &div_m,
        const Tensor             &grad_m,
        const Scalar             &e,
        const Vector             &grad_e,
        const ScalarCoefficient  &refcoeff_div_m,
        const VectorCoefficient1 &refcoeff_grad_rho,
        const VectorCoefficient2 &refcoeff_grad_e)
{
    const Scalar rho_inverse = 1.0/rho;
    const Scalar coeff_div_m(
            explicit_div_e_plus_p_u_refcoeff_div_m(gamma, rho, m, e)
          - refcoeff_div_m);
    const Vector coeff_grad_rho(
            explicit_div_e_plus_p_u_refcoeff_grad_rho(gamma, rho, m, e)
          - refcoeff_grad_rho);
    const Vector coeff_grad_e(
            explicit_div_e_plus_p_u_refcoeff_grad_e(gamma, rho, m)
          - refcoeff_grad_e);
    // FIXME What about linearizing the grad_m term?

    return -(gamma-1)*rho_inverse*rho_inverse*m.dot(grad_m.transpose()*m)
           + coeff_div_m*div_m
           + coeff_grad_rho.dot(grad_rho)
           + coeff_grad_e.dot(grad_e);
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}\cdot\vec{\nabla}\vec{m}\f$ in the explicit portion
 * of \f$\mu\vec{\nabla}\cdot\vec{\nabla}\vec{u}\f$.
 *
 * @param[in] mu \f$\mu\f$
 *            computed from, for example, rhome::p_T_mu_lambda()
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
 *            computed from, for example, rhome::p_T_mu_lambda()
 * @param[in] rho \f$\rho\f$
 * @param[in] m \f$\vec{m}\f$
 *
 * @return \f$\mu\rho^{-2}\vec{m}\f$
 * @see explicit_mu_div_grad_u for more details on the explicit
 *      operator.
 */
template<typename Scalar,
         typename Vector = Eigen::Matrix<Scalar,3,1> >
inline
Vector explicit_mu_div_grad_u_refcoeff_div_grad_rho(
        const Scalar &mu,
        const Scalar &rho,
        const Vector &m)
{
    const Scalar rho_inverse = 1.0/rho;
    return mu*rho_inverse*rho_inverse*m;
}

/**
 * Compute the explicit portion of
 * \f$\mu\vec{\nabla}\cdot\vec{\nabla}\vec{u}\f$.
 * Uses the expansion
 * \f{align*}
 * \mu\vec{\nabla}\cdot\vec{\nabla}\vec{u} = &\phantom{{}+}
 *     \mu\rho^{-2}\left[
 *           2\rho^{-1}\vec{m}\left({\vec\nabla}\rho\right)^{2}
 *         - 2 \left(\vec{\nabla}\vec{m}\right)\vec{\nabla}\rho
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
 * and \f$\left\{\mu\rho^{-2}m\right\}_0\f$ are fixed by
 * \c refcoeff_div_grad_m and \c refcoeff_div_grad_rho, respectively.
 * The remaining linear portion of
 * \f$\mu\vec{\nabla}\cdot\vec{\nabla}\vec{u}\f$ is
 * \f[
 *      \left\{\mu\rho^{-1}\right\}_0 \vec{\nabla}\cdot\vec{\nabla}\vec{m}
 *    - \left\{\mu\rho^{-2}\vec{m}\right\}_0 \vec{\nabla}\cdot\vec{\nabla}\rho
 * \f]
 *
 * @param mu \f$\mu\f$
 *            computed from, for example, rhome::p_T_mu_lambda()
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
         typename Vector            = Eigen::Matrix<Scalar,3,1>,
         typename Tensor            = Eigen::Matrix<Scalar,3,3>,
         typename ScalarCoefficient = Scalar,
         typename VectorCoefficient = Vector >
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
    const Scalar rho_inverse  = 1.0/rho;
    const Scalar rho_inverse2 = rho_inverse*rho_inverse;
    const Scalar coeff_div_grad_m(
            explicit_mu_div_grad_u_refcoeff_div_grad_m(mu, rho)
          - refcoeff_div_grad_m);
    const Vector coeff_div_grad_rho(
            explicit_mu_div_grad_u_refcoeff_div_grad_rho(mu, rho, m)
          - refcoeff_div_grad_rho);

    return   mu*rho_inverse2*(
                    2*rho_inverse*grad_rho.squaredNorm()*m
                  - 2*grad_m*grad_rho
             )
           + coeff_div_grad_m  *div_grad_m
           - coeff_div_grad_rho*div_grad_rho;
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}\vec{\nabla}\cdot\vec{m}\f$ in the explicit portion
 * of \f$\left(\mu+\lambda\right)\vec{\nabla}\vec{\nabla}\cdot\vec{u}\f$.
 *
 * @param[in] mu \f$\mu\f$
 *            computed from, for example, rhome::p_T_mu_lambda()
 * @param[in] lambda \f$\lambda\f$
 *            computed from, for example, rhome::p_T_mu_lambda()
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
 *            computed from, for example, rhome::p_T_mu_lambda()
 * @param[in] lambda \f$\lambda\f$
 *            computed from, for example, rhome::p_T_mu_lambda()
 * @param[in] rho \f$\rho\f$
 * @param[in] m \f$\vec{m}\f$
 *
 * @return \f$\left(\mu+\lambda\right)\rho^{-2}\vec{m}\f$
 * @see explicit_mu_plus_lambda_grad_div_u() for more details on the explicit
 *      operator.
 */
template<typename Scalar,
         typename Vector = Eigen::Matrix<Scalar,3,1> >
inline
Vector explicit_mu_plus_lambda_grad_div_u_refcoeff_grad_grad_rho(
        const Scalar &mu,
        const Scalar &lambda,
        const Scalar &rho,
        const Vector &m)
{
    const Scalar rho_inverse = 1.0/rho;
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
 *      - \vec{\nabla}\vec{\nabla}\rho \left(
 *              \left(\mu+\lambda\right)\rho^{-2}\vec{m}
 *            - \left\{\left(\mu+\lambda\right)\rho^{-2}\vec{m}\right\}_0
 *        \right)
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
 *      - \vec{\nabla}\vec{\nabla}\rho
 *        \left\{\left(\mu+\lambda\right)\rho^{-2}\vec{m}\right\}_0
 * \f]
 *
 * @param mu \f$\mu\f$
 *            computed from, for example, rhome::p_T_mu_lambda()
 * @param lambda \f$\lambda\f$
 *            computed from, for example, rhome::p_T_mu_lambda()
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
         typename Vector            = Eigen::Matrix<Scalar,3,1>,
         typename Tensor            = Eigen::Matrix<Scalar,3,3>,
         typename ScalarCoefficient = Scalar,
         typename VectorCoefficient = Vector >
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
    const Scalar rho_inverse  = 1.0/rho;
    const Scalar rho_inverse2 = rho_inverse*rho_inverse;
    const Scalar coeff_grad_div_m(
            explicit_mu_plus_lambda_grad_div_u_refcoeff_grad_div_m(
              mu, lambda, rho)
          - refcoeff_grad_div_m);
    const Vector coeff_grad_grad_rho(
            explicit_mu_plus_lambda_grad_div_u_refcoeff_grad_grad_rho(
              mu, lambda, rho, m)
          - refcoeff_grad_grad_rho);

    return   (mu+lambda)*rho_inverse2*(
                  (2*rho_inverse*grad_rho.dot(m) - div_m)*grad_rho
                - grad_m.transpose()*grad_rho
             )
           + coeff_grad_div_m*grad_div_m
           - grad_grad_rho*coeff_grad_grad_rho;
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}\cdot\vec{\nabla}e\f$ in the explicit portion
 * of \f$\mu\vec{\nabla}\cdot\vec{\nabla}T\f$.
 *
 * @param[in] mu \f$\mu\f$
 *            computed from, for example, rhome::p_T_mu_lambda()
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
 *            computed from, for example, rhome::p_T_mu_lambda()
 * @param[in] rho \f$\rho\f$
 * @param[in] m \f$\vec{m}\f$
 *
 * @return \f$\mu\rho^{-2}\vec{m}\f$
 * @see explicit_mu_div_grad_T() for more details on the explicit
 *      operator.
 */
template<typename Scalar,
         typename Vector = Eigen::Matrix<Scalar,3,1> >
inline
Vector explicit_mu_div_grad_T_refcoeff_div_grad_m(
        const Scalar &mu,
        const Scalar &rho,
        const Vector &m)
{
    const Scalar rho_inverse = 1.0/rho;
    return mu*rho_inverse*rho_inverse*m;
}

/**
 * Compute the reference coefficient for the term containing
 * \f$\vec{\nabla}\cdot\vec{\nabla}\rho\f$ in the explicit portion
 * of \f$\mu\vec{\nabla}\cdot\vec{\nabla}T\f$.
 *
 * @param[in] gamma \f$\gamma\f$
 * @param[in] mu \f$\mu\f$
 *            computed from, for example, rhome::p_T_mu_lambda()
 * @param[in] rho \f$\rho\f$
 * @param[in] m \f$\vec{m}\f$
 * @param[in] p \f$p\f$
 *            computed from, for example, rhome::p_T_mu_lambda()
 *
 * @return \f$\mu\rho^{-2}\left(
 *            \frac{\gamma-1}{2}\rho^{-1}\vec{m}^2-p\right)\f$
 * @see explicit_mu_div_grad_T() for more details on the explicit
 *      operator.
 */
template<typename Scalar,
         typename Vector = Eigen::Matrix<Scalar,3,1> >
inline
Scalar explicit_mu_div_grad_T_refcoeff_div_grad_rho(
        const Scalar &gamma,
        const Scalar &mu,
        const Scalar &rho,
        const Vector &m,
        const Scalar &p)
{
    const Scalar rho_inverse = 1.0/rho;
    return mu*rho_inverse*rho_inverse*(
             0.5*(gamma-1.0)*rho_inverse*m.squaredNorm() - p
           );
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
 *            \mu\rho^{-2}\left(\frac{\gamma-1}{2}\rho^{-1}\vec{m}^2-p\right)
 *          - \left\{
 *              \mu\rho^{-2}\left(\frac{\gamma-1}{2}\rho^{-1}\vec{m}^2-p\right)
 *            \right\}_0
 *        \right)\vec{\nabla}\cdot\vec{\nabla}\rho
 * \f}
 * where \f$\left\{\mu\rho^{-1}\right\}_0\f$,
 *  \f$\left\{\mu\rho^{-2}\vec{m}\right\}_0\f$,
 * and \f$
 *   \left\{
 *     \mu\rho^{-2}\left(\frac{\gamma-1}{2}\rho^{-1}\vec{m}^2-p\right)
 *   \right\}_0
 * \f$
 * are fixed by \c refcoeff_div_grad_e, \c refcoeff_div_grad_m,
 * and \c refcoeff_div_grad_rho, respectively.
 * The remaining linear portion of
 * \f$\mu\vec{\nabla}\cdot\vec{\nabla}T\f$ is
 * \f[
 *      \gamma\left(\gamma-1\right)
 *      \left\{\mu\rho^{-1}\right\}_0
 *      \vec{\nabla}\cdot\vec{\nabla}e
 *    - \gamma\left(\gamma-1\right)
 *      \left\{\mu\rho^{-2}\vec{m}\right\}_0
 *      \cdot\vec{\nabla}\cdot\vec{\nabla}\vec{m}
 *    + \gamma
 *      \left\{
 *        \mu\rho^{-2}\left(\frac{\gamma-1}{2}\rho^{-1}\vec{m}^2-p\right)
 *      \right\}_0
 *      \vec{\nabla}\cdot\vec{\nabla}\rho
 * \f]
 *
 * @param gamma \f$\gamma\f$
 * @param mu \f$\mu\f$
 *            computed from, for example, rhome::p_T_mu_lambda()
 * @param rho \f$\rho\f$
 * @param grad_rho \f$\vec{\nabla}\rho\f$
 * @param div_grad_rho \f$\vec{\nabla}\cdot\vec{\nabla}\rho\f$
 * @param m \f$\vec{m}\f$
 * @param grad_m \f$\vec{\nabla}\vec{m}\f$
 * @param div_grad_m \f$\vec{\nabla}\cdot\vec{\nabla}\vec{m}\f$
 * @param div_grad_e \f$\vec{\nabla}\cdot\vec{\nabla}e\f$
 * @param p \f$p\f$
 *            computed from, for example, rhome::p_T_mu_lambda()
 * @param grad_p \f$\vec{\nabla}p\f$
 *            computed from, for example, rhome::p_T_mu_lambda()
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
         typename Vector             = Eigen::Matrix<Scalar,3,1>,
         typename Tensor             = Eigen::Matrix<Scalar,3,3>,
         typename ScalarCoefficient1 = Scalar,
         typename ScalarCoefficient2 = Scalar,
         typename VectorCoefficient  = Vector >
Scalar explicit_mu_div_grad_T(
        const Scalar             &gamma,
        const Scalar             &mu,
        const Scalar             &rho,
        const Vector             &grad_rho,
        const Scalar             &div_grad_rho,
        const Vector             &m,
        const Tensor             &grad_m,
        const Vector             &div_grad_m,
        const Scalar             &div_grad_e,
        const Scalar             &p,
        const Vector             &grad_p,
        const ScalarCoefficient1 &refcoeff_div_grad_e,
        const VectorCoefficient  &refcoeff_div_grad_m,
        const ScalarCoefficient2 &refcoeff_div_grad_rho)
{
    const Scalar rho_inverse  = 1.0/rho;
    const Scalar rho_inverse2 = rho_inverse*rho_inverse;
    const Scalar coeff_div_grad_e(
            explicit_mu_div_grad_T_refcoeff_div_grad_e(mu, rho)
          - refcoeff_div_grad_e);
    const Vector coeff_div_grad_m(
            explicit_mu_div_grad_T_refcoeff_div_grad_m(mu, rho, m)
          - refcoeff_div_grad_m);
    const Scalar coeff_div_grad_rho(
            explicit_mu_div_grad_T_refcoeff_div_grad_rho(gamma, mu, rho, m, p)
          - refcoeff_div_grad_rho);

    return gamma*(
                coeff_div_grad_rho*div_grad_rho
              - 2.0*mu*rho_inverse2*grad_rho.dot(
                    grad_p - rho_inverse*p*grad_rho
                )
              + (gamma-1.0)*(
                  - mu*rho_inverse2*(
                        (grad_m.transpose()*grad_m).trace()
                      - rho_inverse*(
                            2.0*grad_rho.dot(grad_m.transpose()*m)
                          - rho_inverse*m.squaredNorm()*grad_rho.squaredNorm()
                        )
                    )
                  + coeff_div_grad_e*div_grad_e
                  - div_grad_m.dot(coeff_div_grad_m)
                )
           );
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
 *          - \rho^{-1} \left(\vec{\nabla}\vec{m}\right)^{\mathsf{T}}\vec{m}
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
    const Scalar rho_inverse = 1.0/rho;

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
         typename Vector = Eigen::Matrix<Scalar,3,1>,
         typename Tensor = Eigen::Matrix<Scalar,3,3> >
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
    const Scalar rho_inverse = 1.0/rho;
    const Scalar m_dot_m     = m.squaredNorm();

    return (gamma - 1.0)*(
                  div_grad_e
                - rho_inverse*(
                          (grad_m.transpose()*grad_m).trace()
                        + div_grad_m.dot(m)
                        - rho_inverse*(
                                2.0*grad_rho.dot(grad_m.transpose()*m)
                              + 0.5*m_dot_m*div_grad_rho
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
         typename Vector = Eigen::Matrix<Scalar,3,1> >
Scalar div_grad_T(
        const Scalar &gamma,
        const Scalar &rho,
        const Vector &grad_rho,
        const Scalar &div_grad_rho,
        const Scalar &p,
        const Vector &grad_p,
        const Scalar &div_grad_p)
{
    const Scalar rho_inverse = 1.0/rho;

    return gamma*rho_inverse*(
                  div_grad_p
                - rho_inverse*(
                          p*div_grad_rho
                        + 2.0*grad_rho.dot(grad_p - rho_inverse*p*grad_rho)
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

#endif // __SUZERAIN_ORTHONORMAL_H
