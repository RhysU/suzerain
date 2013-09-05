//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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

#ifndef SUZERAIN_RHOLUT_HPP
#define SUZERAIN_RHOLUT_HPP

/** @file
 * Provides compute kernels for formulations using a reference density
 * \f$\rho_0\f$, length \f$l_0\f$, velocity \f$u_0\f$, and temperature
 * \f$T_0\f$.
 */

#include <suzerain/common.hpp>
#include <suzerain/rholt.hpp>

namespace suzerain {

/**
 * Provides compute kernels for formulations using a reference density
 * \f$\rho_0\f$, length \f$l_0\f$, velocity \f$u_0\f$, and temperature
 * \f$T_0\f$.  The ideal gas equation of state is employed:
 * \f{align*}
 *     p &= \left(\gamma-1\right) \left(
 *       e - \mbox{Ma}^{2} \frac{m\cdot{}m}{2\rho}
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
 *   - \f$\mbox{Ma}\f$ or \c Ma is the Mach number \f$u_0/a_0\f$.
 *
 * Templated \c Scalar, \c Vector, \c and Tensor types
 * are used.  These functions are intended to be used with <a
 * href="http://eigen.tuxfamily.org/">Eigen</a>'s vector and matrix
 * types.  For example, when <tt>Scalar == double</tt> \c Vector
 * and \c Tensor may be <tt>Eigen::Matrix<Scalar,3,1></tt> and
 * <tt>Eigen::Matrix<Scalar,3,3></tt>, respectively.
 *
 * Many compute kernels in this nondimensionalization are identical to ones
 * used in suzerain::rholt.  Those kernels have been brought into this
 * namespace with using declarations.
 */
namespace rholut {


/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::tau
 */
using suzerain::rholt::tau;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::div_rho_inverse_m_outer_m
 */
using suzerain::rholt::div_rho_inverse_m_outer_m;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::div_u_outer_m
 */
using suzerain::rholt::div_u_outer_m;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::div_tau
 */
using suzerain::rholt::div_tau;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::div_e_u
 */
using suzerain::rholt::div_e_u;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::div_mu_grad_T
 */
using suzerain::rholt::div_mu_grad_T;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::div_p_u
 */
using suzerain::rholt::div_p_u;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::div_tau_u
 */
using suzerain::rholt::div_tau_u;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_div_rho_inverse_m_outer_m
 */
using suzerain::rholt::explicit_div_rho_inverse_m_outer_m;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_grad_p_refcoeff_grad_rho
 */
using suzerain::rholt::explicit_grad_p_refcoeff_grad_rho;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_grad_p_refcoeff_grad_m
 */
using suzerain::rholt::explicit_grad_p_refcoeff_grad_m;

/**
 * Compute the explicit portion of \f$\vec{\nabla}p\f$.
 * Uses the expansion
 * \f{align*}
 *   \nabla{}p &= \left(\gamma-1\right)\mbox{Ma}^{2}\left(
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
 *   + \frac{\gamma-1}{2}\mbox{Ma}^{2}
 *     \left\{\frac{\vec{m}\cdot\vec{m}}{\rho^{2}}\right\}_{0} \vec{\nabla}\rho
 *   - \left(\gamma-1\right)\mbox{Ma}^{2}
 *     \left(\vec{\nabla}\vec{m}\right)^{\mathsf{T}}
 *     \left\{\frac{m}{\rho}\right\}_{0}
 * \f]
 * Unlike many other methods beginning with <tt>explicit_</tt>, when zero
 * reference coefficients are used the full \f$\vec{\nabla}p\f$ term is
 * <i>not</i> recovered.
 *
 * @param gamma \f$\gamma\f$
 * @param Ma \f$\mbox{Ma}\f$
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
        const Scalar            &Ma,
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

    return (gamma-1)*Ma*Ma*(
                  (coeff_grad_rho/2)*grad_rho
                - grad_m.transpose()*coeff_grad_m
            );
}

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_div_p_u_refcoeff_div_m
 */
using suzerain::rholt::explicit_div_p_u_refcoeff_div_m;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_div_p_u_refcoeff_grad_rho
 */
using suzerain::rholt::explicit_div_p_u_refcoeff_grad_rho;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_div_p_u
 */
using suzerain::rholt::explicit_div_p_u;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_div_e_plus_p_u_refcoeff_div_m
 */
using suzerain::rholt::explicit_div_e_plus_p_u_refcoeff_div_m;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_div_e_plus_p_u_refcoeff_grad_rho
 */
using suzerain::rholt::explicit_div_e_plus_p_u_refcoeff_grad_rho;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_div_e_plus_p_u_refcoeff_grad_e
 */
using suzerain::rholt::explicit_div_e_plus_p_u_refcoeff_grad_e;

/**
 * Compute the explicit portion of
 * \f$\vec{\nabla}\cdot\left(e+p\right)\vec{u}\f$.
 * Uses the expansion
 * \f{align*}
 * \vec{\nabla}\cdot\left(e+p\right)\vec{u} =
 *  &- \left(\gamma-1\right)\mbox{Ma}^{2}\rho^{-2}\vec{m}\cdot
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
 * where \f$\rho^{-1}\left(e+p\right)\f$,
 * \f$\rho^{-2}\vec{m}\left((\gamma-2)e-2p\right)\f$, and
 * \f$\left\{\rho^{-1}\vec{m}\right\}_0\f$ are fixed by \c refcoeff_div_m, \c
 * refcoeff_grad_rho, and \c refcoeff_grad_e, respectively.  The remaining
 * linear portion of
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
 * @param Ma \f$\mbox{Ma}\f$
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
        const Scalar             &Ma,
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

    return   ((1-gamma)*Ma*Ma)/(rho*rho)*m.dot(grad_m.transpose()*m)
           + coeff_div_m*div_m
           + coeff_grad_rho.dot(grad_rho)
           + gamma*coeff_grad_e.dot(grad_e);
}

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_mu_div_grad_u_refcoeff_div_grad_m
 */
using suzerain::rholt::explicit_mu_div_grad_u_refcoeff_div_grad_m;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_mu_div_grad_u_refcoeff_div_grad_rho
 */
using suzerain::rholt::explicit_mu_div_grad_u_refcoeff_div_grad_rho;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_mu_div_grad_u
 */
using suzerain::rholt::explicit_mu_div_grad_u;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_u_dot_mu_div_grad_u_refcoeff_div_grad_m
 */
using suzerain::rholt::explicit_u_dot_mu_div_grad_u_refcoeff_div_grad_m;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_u_dot_mu_div_grad_u_refcoeff_div_grad_rho
 */
using suzerain::rholt::explicit_u_dot_mu_div_grad_u_refcoeff_div_grad_rho;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_u_dot_mu_div_grad_u
 */
using suzerain::rholt::explicit_u_dot_mu_div_grad_u;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_mu_plus_lambda_grad_div_u_refcoeff_grad_div_m
 */
using suzerain::rholt::explicit_mu_plus_lambda_grad_div_u_refcoeff_grad_div_m;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_mu_plus_lambda_grad_div_u_refcoeff_grad_grad_rho
 */
using suzerain::rholt::explicit_mu_plus_lambda_grad_div_u_refcoeff_grad_grad_rho;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_mu_plus_lambda_grad_div_u
 */
using suzerain::rholt::explicit_mu_plus_lambda_grad_div_u;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_u_dot_mu_plus_lambda_grad_div_u_refcoeff_grad_div_m
 */
using suzerain::rholt::explicit_u_dot_mu_plus_lambda_grad_div_u_refcoeff_grad_div_m;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_u_dot_mu_plus_lambda_grad_div_u_refcoeff_grad_grad_rho
 */
using suzerain::rholt::explicit_u_dot_mu_plus_lambda_grad_div_u_refcoeff_grad_grad_rho;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_u_dot_mu_plus_lambda_grad_div_u
 */
using suzerain::rholt::explicit_u_dot_mu_plus_lambda_grad_div_u;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_mu_div_grad_T_refcoeff_div_grad_e
 */
using suzerain::rholt::explicit_mu_div_grad_T_refcoeff_div_grad_e;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_mu_div_grad_T_refcoeff_div_grad_m
 */
using suzerain::rholt::explicit_mu_div_grad_T_refcoeff_div_grad_m;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::explicit_mu_div_grad_T_refcoeff_div_grad_rho
 */
using suzerain::rholt::explicit_mu_div_grad_T_refcoeff_div_grad_rho;

/**
 * Compute the explicit portion of \f$\mu\vec{\nabla}\cdot\vec{\nabla}T\f$.
 * Uses the expansion
 * \f{align*}
 * \mu\vec{\nabla}\cdot\vec{\nabla}T =
 *   &{}- 2\gamma\mu\rho^{-2}\vec{\nabla}\rho\cdot
 *        \left(\vec{\nabla}p-\rho^{-1}p\vec{\nabla}\rho\right)
 * \\
 *   &{}- \gamma\left(\gamma-1\right)\mbox{Ma}^{2}\mu\rho^{-2}\left[
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
 *      - \gamma\left(\gamma-1\right)\mbox{Ma}^{2}\left(
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
 *    - \gamma\left(\gamma-1\right)\mbox{Ma}^{2}
 *      \left\{\mu\rho^{-2}\vec{m}\right\}_0
 *      \cdot\vec{\nabla}\cdot\vec{\nabla}\vec{m}
 *    + \gamma
 *      \left\{\mu\rho^{-2}\left((\gamma-1)e-2p\right)\right\}_0
 *      \vec{\nabla}\cdot\vec{\nabla}\rho
 * \f]
 *
 * @param gamma \f$\gamma\f$
 * @param Ma \f$\mbox{Ma}\f$
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
        const Scalar             &Ma,
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
                    coeff_div_grad_e*div_grad_e
                  - Ma * Ma * (
                      mu*rho_inverse2*(
                          grad_m.squaredNorm() // Frobenius norm
                        - rho_inverse*(
                              2*grad_rho.dot(grad_m.transpose()*m)
                            - rho_inverse*m.squaredNorm()*grad_rho.squaredNorm()
                          )
                      )
                    + div_grad_m.dot(coeff_div_grad_m)
                  )
                )
           );
}

/**
 * Compute \f$p\f$, \f$T\f$, \f$\mu\f$, and \f$\lambda\f$
 * using the equation of state.
 *
 * @param[in]  alpha \f$\alpha\f$
 * @param[in]  beta \f$\beta\f$
 * @param[in]  gamma \f$\gamma\f$
 * @param[in]  Ma \f$\mbox{Ma}\f$
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
        const Scalar &Ma,
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
    p      = (gamma - 1)*(e - Ma*Ma*rho_inverse*m.squaredNorm()/2);
    T      = gamma * p * rho_inverse;
    mu     = pow(T, beta);
    lambda = (alpha - Scalar(2)/3)*mu;
}

/**
 * Compute \f$p\f$ and \f$T\f$ using the equation of state.
 *
 * @param[in]  alpha \f$\alpha\f$
 * @param[in]  beta \f$\beta\f$
 * @param[in]  gamma \f$\gamma\f$
 * @param[in]  Ma \f$\mbox{Ma}\f$
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
        const Scalar &Ma,
        const Scalar &rho,
        const Vector &m,
        const Scalar &e,
        Scalar &p,
        Scalar &T)
{
    SUZERAIN_UNUSED(alpha);  // Present for API consistency
    SUZERAIN_UNUSED(beta);   // Present for API consistency
    const Scalar rho_inverse = 1/rho;

    // Compute scalar quantities
    p = (gamma - 1)*(e - Ma*Ma*rho_inverse*m.squaredNorm()/2);
    T = gamma * p * rho_inverse;
}

/**
 * Compute \f$p\f$ using the equation of state.
 *
 * @param[in]  alpha \f$\alpha\f$
 * @param[in]  beta \f$\beta\f$
 * @param[in]  gamma \f$\gamma\f$
 * @param[in]  Ma \f$\mbox{Ma}\f$
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
       const Scalar &Ma,
       const Scalar &rho,
       const Vector &m,
       const Scalar &e,
       Scalar &p)
{
    SUZERAIN_UNUSED(alpha);  // Present for API consistency
    SUZERAIN_UNUSED(beta);   // Present for API consistency

    p = (gamma - 1)*(e - Ma*Ma*m.squaredNorm()/2/rho);
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
 * Compute the internal energy by subtracting the kinetic
 * energy from the total energy.
 *
 * @param[in] Ma \f$\mbox{Ma}\f$
 * @param[in] rho \f$\rho\f$
 * @param[in] m \f$\vec{m}\f$
 * @param[in] e \f$e\f$
 */
template<typename Scalar,
         typename Vector  >
inline
Scalar energy_internal(
        const Scalar &Ma,
        const Scalar &rho,
        const Vector &m,
        const Scalar &e)
{
    return e - Ma*Ma*m.squaredNorm()/rho/2;
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
 * @param[in] Ma \f$\mbox{Ma}\f$
 * @param[in] rho \f$\rho\f$
 * @param[in] m \f$\vec{m}\f$
 */
template<typename Scalar,
         typename Vector  >
inline
Scalar energy_kinetic(
        const Scalar &Ma,
        const Scalar &rho,
        const Vector &m)
{
    return Ma*Ma*m.squaredNorm()/rho/2;
}

/**
 * Compute the gradient of the kinetic energy.
 *
 * @param[in] Ma       \f$\mbox{Ma}\f$
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
        const Scalar &Ma,
        const Scalar &rho,
        const Vector &grad_rho,
        const Vector &m,
        const Tensor &grad_m)
{
    const Scalar rho_inverse  = 1/rho;
    return rho_inverse*Ma*Ma*(
               grad_m.transpose()*m
             - rho_inverse/2*m.squaredNorm()*grad_rho
           );
}

/**
 * Compute \f$p\f$, \f$T\f$, \f$\mu\f$, and \f$\lambda\f$ and their gradients
 * using the equation of state.  The gradients are computed using these
 * expansions:
 * \f{align*}
 *      \vec{\nabla}p &= (\gamma-1)\left[
 *            \vec{\nabla}e
 *          + \frac{\mbox{Ma}^{2}}{\rho} \left[
 *               \frac{\vec{m}\cdot\vec{m}}{2\rho} \vec{\nabla}{\rho}
 *             - \left(\vec{\nabla}\vec{m}\right)^{\mathsf{T}}\vec{m}
 *            \right]
 *      \right]
 *      \\
 *      \vec{\nabla}T &= \frac{\gamma}{\rho}\left(
 *                          \vec{\nabla}p - \frac{p}{\rho}\vec{\nabla}\rho
 *                       \right)
 *                     = \rho^{-1} \left(
 *                          \gamma\vec{\nabla}p - T\vec{\nabla}\rho
 *                       \right)
 *      \\
 *      \vec{\nabla}\mu &= \beta{}T^{\beta-1}\vec{\nabla}T
 *      \\
 *      \vec{\nabla}\lambda &= \left(\alpha-\frac{2}{3}\right)\vec{\nabla}\mu
 * \f}
 *
 * @param[in]  alpha \f$\alpha\f$
 * @param[in]  beta \f$\beta\f$
 * @param[in]  gamma \f$\gamma\f$
 * @param[in]  Ma \f$\mbox{Ma}\f$
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
        const Scalar &Ma,
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
    const Scalar squared_Ma       = Ma * Ma;

    // Compute scalar quantities
    p      = gamma1*(e - squared_Ma*half_rho_inverse*m.squaredNorm());
    T      = gamma * p * rho_inverse;

    const Scalar T_to_beta1 = pow(T, beta - 1); // Avoid div for grad_mu

    mu     = T_to_beta1 * T;
    lambda = alpha23*mu;

    // Compute vector quantities
    grad_p = gamma1*(
                grad_e + squared_Ma*rho_inverse*(
                      (half_rho_inverse*m.squaredNorm())*grad_rho
                    - grad_m.transpose()*m
                )
             );
    grad_T      = rho_inverse*(gamma*grad_p - T*grad_rho);
    grad_mu     = (beta * T_to_beta1)*grad_T;
    grad_lambda = alpha23*grad_mu;
}

/**
 * Compute \f$\vec{\nabla}\cdot\vec\nabla{}p\f$
 * using the equation of state.  Uses the expansion
 * \f{align*}
 *      \vec{\nabla}\cdot\vec{\nabla}p =
 *      \left(\gamma-1\right)\left[
 *          \vec{\nabla}\cdot\vec{\nabla}e
 *          - \frac{\mbox{Ma}^2}{\rho}\left[
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
 * @param[in]  Ma \f$\mbox{Ma}\f$
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
        const Scalar &Ma,
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
                - Ma*Ma*rho_inverse*(
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
 * \namespace suzerain::rholut
 * \see suzerain::rholt::div_grad_T
 */
using suzerain::rholt::div_grad_T;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::u
 */
using suzerain::rholt::u;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::grad_u
 */
using suzerain::rholt::grad_u;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::div_u
 */
using suzerain::rholt::div_u;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::grad_div_u
 */
using suzerain::rholt::grad_div_u;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::div_grad_u
 */
using suzerain::rholt::div_grad_u;

/**
 * \namespace suzerain::rholut
 * \see suzerain::rholt::curl_curl_u
 */
using suzerain::rholt::curl_curl_u;

} // namespace rholut

} // namespace suzerain

#endif // SUZERAIN_RHOLUT_HPP
