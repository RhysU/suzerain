/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2013 Rhys Ulerich
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
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
 */

#ifndef SUZERAIN_RADIAL_NOZZLE_H
#define SUZERAIN_RADIAL_NOZZLE_H

#include <stddef.h>

/** @file
 * A GNU Scientific Library-based axisymmetric radial nozzle solver.  Use \ref
 * suzerain_radial_nozzle_solver to produce a \ref
 * suzerain_radial_nozzle_solution.
 *
 * @see Model documentation in <tt>writeups/baseflow.tex</tt> for details.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Stores pointwise radial nozzle extended state at some point \f$R\f$.
 */
typedef struct suzerain_radial_nozzle_state {
    double R;       /**< Radius                     \f$R              \f$ */
    double u;       /**< Radial velocity            \f$u              \f$ */
    double rho;     /**< Density                    \f$\rho           \f$ */
    double p;       /**< Pressure                   \f$p              \f$ */
    double a2;      /**< Sound speed squared        \f$a^2            \f$ */
    double up;      /**< Radial velocity derivative \f$\partial_R u   \f$ */
    double rhop;    /**< Density derivative         \f$\partial_R \rho\f$ */
    double pp;      /**< Pressure derivative        \f$\partial_R p   \f$ */
} suzerain_radial_nozzle_state;

/**
 * A radial nozzle solution produced by \ref suzerain_radial_nozzle_solver.
 * The struct hack is used to provide solution-wide parameters, e.g. \c Ma0,
 * followed by \c size local details in \c state.
 */
typedef struct suzerain_radial_nozzle_solution {
    double Ma0;     /**< Reference Mach number         \f$\mbox{Ma}_0\f$ */
    double gam0;    /**< Reference specific heat ratio \f$\gamma_0   \f$ */
    size_t size;    /**< Number of entries in flexible member \c state */

    /** Pointwise solution details indexable per the struct hack. */
    suzerain_radial_nozzle_state state[0];
} suzerain_radial_nozzle_solution;

/**
 * Solve the requested radial nozzle problem.  That is, advance the system
 * \f{align}{
 *     u^\prime
 *   &=
 *     -\frac{u}{r}
 *     \cdot
 *     \frac{
 *        2
 *      + \mbox{Ma}_{0}^2\left(\gamma_0-1\right)
 *      - \mbox{Ma}_{0}^2\left(\gamma_0-1\right) u^2
 *     }{
 *        2
 *      + \mbox{Ma}_{0}^2\left(\gamma_0-1\right)
 *      - \mbox{Ma}_{0}^2\left(\gamma_0+1\right) u^2
 *     }
 *   \\
 *     \left(\log\rho\right)^\prime
 *   &=
 *     -\frac{\mbox{Ma}_{0}^2}{2}\frac{\left(u^2\right)^\prime}{a^2}
 *   \\
 *     p^\prime
 *   &=
 *     - \frac{1}{2}\mbox{Ma}_{0}^2 \, \rho\left(u^2\right)^\prime
 * \f}
 * employing the equation of state
 * \f{align}{
 *     a^2
 *   &=
 *       \frac{\gamma-1}{\gamma_0-1}
 *     + \mbox{Ma}_{0}^2\frac{\gamma-1}{2}\left(1-u^2\right)
 * \f}
 * to solve for the quantities in \ref suzerain_radial_nozzle_state as a
 * function of \f$R\f$ given initial conditions \f$\rho\left(r_1\right)\f$,
 * \f$u\left(r_1\right)\f$, and \f$p\left(r_1\right)\f$ satisfying the
 * realizability constraint
 * \f{align}{
 *     u^2
 *   &<
 *     \frac{2}{\mbox{Ma}_{0}^2\left(\gamma_0-1\right)} + 1
 * .
 * \f}
 *
 * @param Ma0   Reference Mach number         \f$\mbox{Ma}_0\f$
 * @param gam0  Reference specific heat ratio \f$\gamma_0   \f$
 * @param rho1  Initial density               \f$\rho\left(r_1\right)\f$
 * @param u1    Initial radial velocity       \f$u   \left(r_1\right)\f$
 * @param p1    Initial pressure              \f$p   \left(r_1\right)\f$
 * @param R     Radii of interest with \f$R_1\f$ taken from \c R[0].
 *              Must be a contiguous array of length \c size.
 *              Entries must be strictly increasing.
 * @param size  Number of radii within \c R.
 *
 * @return On success returns a \ref suzerain_radial_nozzle_solution
 *         encapsulating the solution details.  The caller is responsible for
 *         subsequently <tt>free</tt>ing the returned pointer.
 *         On failure calls \ref suzerain_error() and returns \c NULL.
 *
 * @see Model documentation in <tt>writeups/baseflow.tex</tt> for details.
 */
suzerain_radial_nozzle_solution *
suzerain_radial_nozzle_solver(
    const double         Ma0,
    const double         gam0,
    const double         rho1,
    const double         u1,
    const double         p1,
    const double * const R,
    const size_t         size);

/**
 * Compute distance \f$\delta_i\f$ to the \f$x\f$-axis under Cartesian base
 * flow geometry assumptions.  Takes locations <tt>s->state[j].R</tt> to be
 * along an \f$x\f$-oriented line segment beginning at Cartesian position
 * \f$\left(R_0,0\right)\f$ and extending until \f$\left(R_0,\delta\right)\f$
 * where \f$\delta=\sqrt{R_i^2 - R_0^2}\f$.
 *
 * @param s A solution generated by \ref suzerain_radial_nozzle_solver.
 * @param i An indexed location within \c s in range <tt>[0,s->size)</tt>.
 *
 * @return\f$\delta=\sqrt{R_i^2 - R_0^2}\f$.
 *
 * @see Model documentation in <tt>writeups/baseflow.tex</tt> for details.
 */
double
suzerain_radial_nozzle_delta(
    const suzerain_radial_nozzle_solution * s,
    const size_t i);

/**
 * Compute \f$\mbox{Ma}_e\f$, the edge Mach number, at location index \c i
 * within solution \c s.  That is,
 * \f{align}{
 *   \mbox{Ma}_{e}
 *   &\equiv
 *   \left. \frac{u_0 u_\xi}{a_0 a} \right|_{\left(R_0,\delta\right)}
 *   =
 *   \left.
 *     \frac{\mbox{Ma}_{0} R_0}{R}
 *     \frac{\left|u\left(R\right)\right|}
 *          {      a\left(R\right)       }
 *   \right|_{R = \sqrt{R_0^2 + \delta^2}}
 * \f}
 * where \f$\xi\f$ is denotes the \f$x\f$ direction possibly reflected so that
 * streamwise velocity always has positive sign. \f$\delta\f$ is computed from
 * \ref suzerain_radial_nozzle_delta using \c s and \c i.  Given this geometric
 * assumption, only the \f$\xi\f$-component of the radial velocity
 * \f$u\left(R\right)\f$ contributes to \f$\mbox{Ma}_e\f$.
 *
 * @param s A solution generated by \ref suzerain_radial_nozzle_solver.
 * @param i An indexed location within \c s in range <tt>[0,s->size)</tt>.
 *
 * @return \f$\mbox{Ma}_e\f$ computed using described geometric assumptions.
 *
 * @see Model documentation in <tt>writeups/baseflow.tex</tt> for details.
 */
double
suzerain_radial_nozzle_qoi_Mae(
    const suzerain_radial_nozzle_solution * s,
    const size_t i);

/**
 * Compute \f$p^\ast_{e,\xi}\f$, a nondimensional pressure gradient parameter,
 * at location index \c i within solution \c s.  That is,
 * \f{align}{
 *   p^\ast_{e,\xi}
 *   &\equiv
 *   \left.
 *   \frac{l_0 \delta}{\rho_0 \rho \, u_0^2 u^2}
 *     \frac{\partial\left(p_0 p\right)}{\partial\left(l_0 \xi\right)}
 *   \right|_{\left(R_0,\delta\right)}
 *   =
 *   \left.
 *     \frac{\operatorname{sgn}(u) \, \delta}{\mbox{Ma}_{0}^2 \rho u^2}
 *       \frac{\partial{}p}{\partial{}x}
 *   \right|_{\left(R_0,\delta\right)}
 *   =
 *   \left.
 *     \frac{\operatorname{sgn}(u) R \, \delta \, p^\prime\left(R\right)}
 *          {\mbox{Ma}_{0}^2 R_0 \, \rho\left(R\right) u^2\left(R\right)}
 *   \right|_{R=\sqrt{R_0^2 + \delta^2}}
 * \f}
 * where \f$\xi\f$ is denotes the \f$x\f$ direction possibly reflected so that
 * streamwise velocity always has positive sign. \f$\delta\f$ is computed from
 * \ref suzerain_radial_nozzle_delta using \c s and \c i.  Given this geometric
 * assumption, only the \f$\xi\f$-component of the pressure gradient
 * \f$p^\prime\left(R\right)\f$ contributes to \f$p^\ast_{e\,\xi}\f$.
 *
 * @param s A solution generated by \ref suzerain_radial_nozzle_solver.
 * @param i An indexed location within \c s in range <tt>[0,s->size)</tt>.
 *
 * @return \f$p^\ast_{e,\xi}\f$ computed using described geometric assumptions.
 *
 * @see Model documentation in <tt>writeups/baseflow.tex</tt> for details.
 *
 */
double
suzerain_radial_nozzle_qoi_pexi(
    const suzerain_radial_nozzle_solution * s,
    const size_t i);

/**
 * Compute Cartesian base flow primitive state at \f$\left(R_0,
 * \delta_i\right)\f$, including streamwise derivatives, given a radial nozzle
 * solution.  Downstream is oriented to be in the positive \f$x\f$ direction
 * regardless of the sub- versus supersonic nature of the radial solution.
 *
 * More concretely, compute
 * \f{align}{
 *       \rho  &= \rho\left(R\right)
 *   \\  u_\xi &= \left|u\left(R\right)\right| \frac{x}{R}
 *   \\  u_y   &=       u\left(R\right)        \frac{y}{R}
 *   \\  p     &= \frac{\mbox{Ma}_0^2}{\mbox{Ma}^2} \, p\left(R\right)
 * \f}
 * and
 * \f{align}{
 *     \frac{\partial}{\partial\xi} \rho
 *   &=
 *     \frac{x\operatorname{sgn}(u)}{R} \rho^\prime\left(R\right)
 *   \\
 *     \frac{\partial}{\partial\xi} u_\xi
 *   &=
 *     \operatorname{sgn}(u) \left(
 *         \frac{x^2 u^\prime\left(R\right)}{R^2}
 *       + \frac{y^2 u       \left(R\right)}{R^3}
 *     \right)
 *   \\
 *     \frac{\partial}{\partial\xi} u_y
 *   &=
 *     x y \operatorname{sgn}(u) \left(
 *         \frac{u^\prime\left(R\right)}{R^2}
 *       - \frac{u       \left(R\right)}{R^3}
 *     \right)
 *   \\
 *     \frac{\partial}{\partial\xi} p
 *   &=
 *     \frac{x\operatorname{sgn}(u)}{R}
 *     \frac{\mbox{Ma}_0^2}{\mbox{Ma}^2}
 *     \,
 *     p^\prime\left(R\right)
 * \f}
 * where \f$x=R_0\f$ and \f$y=\delta\f$.  Direction \f$\xi\f$ is nothing but
 * \f$x\f$ possibly reflected so that streamwise velocity always has positive
 * sign.  Coordinate \f$\delta\f$ is computed from \ref
 * suzerain_radial_nozzle_delta using \c s and \c i.  Parameter \f$\mbox{Ma}\f$
 * permits translating the nondimensional results into a setting where
 * \f$\mbox{Ma} \ne \mbox{Ma}_0\f$.
 *
 * @param[in ] s A solution generated by \ref suzerain_radial_nozzle_solver.
 * @param[in ] i An indexed location within \c s in range <tt>[0,s->size)</tt>.
 * @param[in ] Ma   Target Mach number \f$Ma\f$ for returned pressure values.
 * @param[out] rho  Result \f$\rho \f$
 * @param[out] u    Result \f$u_\xi\f$
 * @param[out] v    Result \f$u_y  \f$
 * @param[out] p    Result \f$p    \f$
 * @param[out] rhop Result \f$\frac{\partial}{\partial\xi} \rho \f$
 * @param[out] up   Result \f$\frac{\partial}{\partial\xi} u_\xi\f$
 * @param[out] vp   Result \f$\frac{\partial}{\partial\xi} u_y  \f$
 * @param[out] pp   Result \f$\frac{\partial}{\partial\xi} p    \f$
 *
 * @see Model documentation in <tt>writeups/baseflow.tex</tt> for details.
 */
void
suzerain_radial_nozzle_cartesian_primitive(
    const suzerain_radial_nozzle_solution * s,
    const size_t i,
    const double Ma,
    double *rho,
    double *u,
    double *v,
    double *p,
    double *rhop,
    double *up,
    double *vp,
    double *pp);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_RADIAL_NOZZLE_H */
