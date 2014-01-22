/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2013-2014 Rhys Ulerich
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

#ifndef SUZERAIN_RADIALFLOW_H
#define SUZERAIN_RADIALFLOW_H

#include <stddef.h>

/** @file
 * A GNU Scientific Library-based axisymmetric radial flow solver. Solves
 * isenthalpic flows in nondimensional primitive variables. Use \ref
 * suzerain_radialflow_solver to produce a \ref suzerain_radialflow_solution.
 *
 * @see Model documentation in <tt>writeups/radialflow.tex</tt> for details.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Stores pointwise radial nozzle extended state at some point \f$R\f$.
 */
typedef struct suzerain_radialflow_state {
    double R;       /**< Radius                     \f$R              \f$ */
    double u;       /**< Radial velocity            \f$u              \f$ */
    double rho;     /**< Density                    \f$\rho           \f$ */
    double p;       /**< Pressure                   \f$p              \f$ */
    double a2;      /**< Sound speed squared        \f$a^2            \f$ */
    double up;      /**< Radial velocity derivative \f$\partial_R u   \f$ */
    double rhop;    /**< Density derivative         \f$\partial_R \rho\f$ */
    double pp;      /**< Pressure derivative        \f$\partial_R p   \f$ */
} suzerain_radialflow_state;

/**
 * A radial nozzle solution produced by \ref suzerain_radialflow_solver.
 * The struct hack is used to provide solution-wide parameters, e.g. \c Ma0,
 * followed by \c size local details in \c state.  The entries in the flexible
 * member \c state \e must be in either increasing or decreasing order by \c R.
 */
typedef struct suzerain_radialflow_solution {
    double Ma0;     /**< Reference Mach number         \f$\mbox{Ma}_0\f$ */
    double gam0;    /**< Reference specific heat ratio \f$\gamma_0   \f$ */
    size_t size;    /**< Number of entries in flexible member \c state */

    /** Pointwise solution details indexable per the struct hack. */
    suzerain_radialflow_state state[0];
} suzerain_radialflow_solution;

/**
 * Solve the requested radial nozzle problem.  That is, advance the system
 * \f{align}{
 *     u^\prime
 *  &=
 *     \frac{u}{r}
 *     \,
 *     \frac{
 *       \left[
 *          2 \mbox{Ma}_0^{-2} + \left(\gamma_0-1\right) \left(1 - u^2\right)
 *       \right]
 *     }{
 *         2 u^2
 *       - \left[
 *            2 \mbox{Ma}0^{-2} + \left(\gamma_0-1\right) \left(1 - u^2\right)
 *         \right]
 *     }
 *   \\
 *     \rho^\prime
 *   &=
 *     \frac{1}{a^2} p^\prime
 *   \\
 *     p^\prime
 *   &=
 *     -\frac{\mbox{Ma}_0^2}{2} \rho u u^\prime
 * \f}
 * employing the equation of state
 * \f{align}{
 *     a^2
 *   &=
 *     1 + \mbox{Ma}_{0}^2\frac{\gamma_0-1}{2}\left(1-u^2\right)
 * \f}
 * to solve for the quantities in \ref suzerain_radialflow_state as a
 * function of \f$R\f$ given initial conditions \f$\rho\left(r_1\right)\f$ and
 * \f$u\left(r_1\right)\f$ satisfying the realizability constraint
 * \f{align}{
 *     u^2
 *   &<
 *     \frac{2}{\mbox{Ma}_{0}^2\left(\gamma_0-1\right)} + 1
 * .
 * \f}
 *
 * @param Ma0   Reference Mach number         \f$\mbox{Ma}_0\f$.
 * @param gam0  Reference specific heat ratio \f$\gamma_0   \f$.
 * @param u1    Initial radial velocity       \f$u   \left(R_1\right)\f$.
 * @param rho1  Initial density               \f$\rho\left(R_1\right)\f$.
 * @param p1    Initial pressure              \f$p   \left(R_1\right)\f$.
 * @param R     Radii of interest with \f$R_1\f$ taken from \c R[0].
 *              Must be a contiguous array of length \c size
 *              with either increasing or decreasing values.
 * @param size  Number of radii within \c R.
 *
 * @return On success returns a \ref suzerain_radialflow_solution
 *         encapsulating the solution details.  The caller is responsible for
 *         subsequently <tt>free</tt>ing the returned pointer.
 *         On failure calls \ref suzerain_error() and returns \c NULL.
 *
 * @see Model documentation in <tt>writeups/radialflow.tex</tt> for details.
 */
suzerain_radialflow_solution *
suzerain_radialflow_solver(
    const double         Ma0,
    const double         gam0,
    const double         u1,
    const double         rho1,
    const double         p1,
    const double * const R,
    const size_t         size);

/**
 * Compute distance \f$\delta_i\f$ to the \f$x\f$-axis under Cartesian base
 * flow geometry assumptions.  Takes locations <tt>s->state[j].R</tt> to be
 * along a positive-\f$y\f$-oriented line segment beginning at Cartesian
 * position \f$\left(R_0,0\right)\f$ and extending until
 * \f$\left(R_0,\delta\right)\f$ where \f$\delta=\sqrt{R_i^2 - R_0^2}\f$ and
 * \f$R_0\f$ is taken from the \e minimum <tt>s->state[j].R</tt>.
 *
 * @param s A solution generated by \ref suzerain_radialflow_solver.
 * @param i An indexed location within \c s in range <tt>[0,s->size)</tt>.
 *
 * @return\f$\delta=\sqrt{R_i^2 - R_0^2}\f$.
 *
 * @see Model documentation in <tt>writeups/radialflow.tex</tt> for details.
 */
double
suzerain_radialflow_delta(
    const suzerain_radialflow_solution * s,
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
 * \ref suzerain_radialflow_delta using \c s and \c i.  Given these geometric
 * assumptions, only the \f$\xi\f$-component of the radial velocity
 * \f$u\left(R\right)\f$ contributes to \f$\mbox{Ma}_e\f$.
 *
 * @param s A solution generated by \ref suzerain_radialflow_solver.
 * @param i An indexed location within \c s in range <tt>[0,s->size)</tt>.
 *
 * @return \f$\mbox{Ma}_e\f$ computed using described geometric assumptions.
 *
 * @see Model documentation in <tt>writeups/radialflow.tex</tt> for details.
 */
double
suzerain_radialflow_qoi_Mae(
    const suzerain_radialflow_solution * s,
    const size_t i);

/**
 * Compute \f$p^\ast_{e,\xi}\f$, a nondimensional pressure gradient parameter,
 * at location index \c i within solution \c s.  That is,
 * \f{align}{
 *    p^\ast_{e,\xi}
 *    &\equiv
 *    \left.
 *    \frac{l_0 \delta}{\rho_0 \rho \, u_0^2 u^2}
 *      \frac{\partial\left(p_0 p\right)}{\partial\left(l_0 \xi\right)}
 *    \right|_{\left(R_0,\delta\right)}
 * \\&=
 *    \left.
 *      \frac{\operatorname{sgn}(u) \, \delta}{\mbox{Ma}_{0}^2 \rho u^2}
 *        \frac{\partial{}p}{\partial{}x}
 *    \right|_{\left(R_0,\delta\right)}
 * \\&=
 *    \left.
 *      \frac{\operatorname{sgn}(u) R \, \delta \, p^\prime\left(R\right)}
 *           {\mbox{Ma}_{0}^2 R_0 \, \rho\left(R\right) u^2\left(R\right)}
 *    \right|_{R=\sqrt{R_0^2 + \delta^2}}
 * \\&=
 *    - \left.
 *        \frac{R \, \delta \, u'\left(R\right)}
 *             {R_0 \, \left|u\left(R\right)\right|}
 *    \right|_{R=\sqrt{R_0^2 + \delta^2}}
 * \f}
 * where \f$\xi\f$ is denotes the \f$x\f$ direction possibly reflected so that
 * streamwise velocity always has positive sign. \f$\delta\f$ is computed from
 * \ref suzerain_radialflow_delta using \c s and \c i.  Given these geometric
 * assumptions, only the \f$\xi\f$-component of the pressure gradient
 * \f$p^\prime\left(R\right)\f$ contributes to \f$p^\ast_{e\,\xi}\f$.
 *
 * @param s A solution generated by \ref suzerain_radialflow_solver.
 * @param i An indexed location within \c s in range <tt>[0,s->size)</tt>.
 *
 * @return \f$p^\ast_{e,\xi}\f$ computed using described geometric assumptions.
 *
 * @see Model documentation in <tt>writeups/radialflow.tex</tt> for details.
 *
 */
double
suzerain_radialflow_qoi_pexi(
    const suzerain_radialflow_solution * s,
    const size_t i);

/**
 * Compute \f$\mbox{T}_e\f$, the edge temperature under ideal gas assumptions,
 * at location index \c i within solution \c s.  That is,
 * \f{align}{
 *   \mbox{T}_{e}
 *   &\equiv
 *   \gamma_0 \frac{p  \left(x,y;\mbox{Ma}\right)}
 *                 {rho\left(x,y          \right)}
 *   =
 *   \left.
 *   \gamma_0 \frac{\mbox{Ma}^2   p  \left(R\right)}
 *                 {\mbox{Ma}^2_0 rho\left(R\right)}
 *   \right|_{R = \sqrt{R_0^2 + \delta^2}}
 *   .
 * \f}
 * Here, \f$\mbox{Ma}\f$ is the Mach number for some target
 * nondimensionalization sharing references \f$u_0\f$, \f$\rho_0\f$, and
 * relationship \f$p_0=\rho_0 a_0^2\f$ but not necessarily \f$a_0\f$.
 *
 * @param s  A solution generated by \ref suzerain_radialflow_solver.
 * @param i  An indexed location within \c s in range <tt>[0,s->size)</tt>.
 * @param Ma Target Mach number \f$Ma\f$ for returned temperature value.
 *
 * @return \f$\mbox{T}_e\f$ computed using described assumptions.
 *
 * @see Model documentation in <tt>writeups/radialflow.tex</tt> for details.
 */
double
suzerain_radialflow_qoi_Te(
    const suzerain_radialflow_solution * s,
    const size_t i,
    const double Ma);

/**
 * Find flow state \f$u\left(R\right)$, \f$\rho\left(R\right)$, and
 * \f$p\left(R\right)$ and radial extents \f$\left[R_0, R\right]\f$ reproducing
 * the specified boundary layer edge parameters.
 *
 * @param[in ] delta Edge distance \f$\delta\f$ above the \f$x\f$-axis.
 * @param[in ] gam0  Reference specific heat ratio \f$\gamma_0\f$.
 * @param[in ] Ma_e  Edge Mach number \f$\mbox{Ma}_e\f$ defined
 *                   per \ref suzerain_radialflow_qoi_Mae.
 * @param[in ] p_exi Pressure gradient parameter \f$p^\ast_{e,\xi}\f$ defined
 *                   per \ref suzerain_radialflow_qoi_pexi.
 * @param[in ] T_e   Edge temperature \f$T_e\f$ defined
 *                   per \ref suzerain_radialflow_qoi_Te.
 * @param[out] Ma0   Reference Mach number \f$\mbox{Ma}_0\f$.
 * @param[out] R0    Radius causing edge radius \f$R\f$ to be at
 *                   Cartesian coordinate \f$\left(R_0,\delta\right)\f$.
 * @param[out] R     Edge radius \f$R\f$.
 * @param[out] uR    Radial velocity at edge, \f$u\left(R\right)\f$
 * @param[out] rhoR  Density at edge \f$\rho\left(R\right)\f$
 * @param[out] pR    Pressure at edge \f$p\left(R\right)\f$
 *
 * @return SUZERAIN_SUCCESS (zero) on success and nonzero on failure.
 *         On failure, all output values are additionally set to NaN.
 *
 * @see Model documentation in <tt>writeups/radialflow.tex</tt> for details.
 */
int
suzerain_radialflow_qoi_match(
    const double delta,
    const double gam0,
    const double Ma_e,
    const double p_exi,
    const double T_e,
    double *Ma0,
    double *R0,
    double *R,
    double *uR,
    double *rhoR,
    double *pR);

/**
 * Compute Cartesian base flow primitive state at \f$\left(R_0,
 * \delta_i\right)\f$, including streamwise derivatives, given a radial nozzle
 * solution.  Downstream \f$\xi\f$ is oriented to be in the positive \f$x\f$
 * direction regardless of the sub- versus supersonic nature of the radial
 * solution.  Coordinate \f$\delta\f$ is computed from \ref
 * suzerain_radialflow_delta using \c s and \c i.  Parameter \f$\mbox{Ma}\f$
 * permits translating the nondimensional results into a setting where
 * \f$\mbox{Ma} \ne \mbox{Ma}_0\f$.
 *
 * @param[in ] s A solution generated by \ref suzerain_radialflow_solver.
 * @param[in ] i An indexed location within \c s in range <tt>[0,s->size)</tt>.
 * @param[in ] Ma     Target Mach number \f$Ma\f$ for returned pressure values.
 * @param[out] rho    Result \f$\rho \f$
 * @param[out] u      Result \f$u_\xi\f$
 * @param[out] v      Result \f$u_y  \f$
 * @param[out] p      Result \f$p    \f$
 * @param[out] p      Result \f$a    \f$
 * @param[out] rho_xi Result \f$\frac{\partial}{\partial\xi} \rho \f$
 * @param[out] u_xi   Result \f$\frac{\partial}{\partial\xi} u_\xi\f$
 * @param[out] v_xi   Result \f$\frac{\partial}{\partial\xi} u_y  \f$
 * @param[out] p_xi   Result \f$\frac{\partial}{\partial\xi} p    \f$
 * @param[out] a_xi   Result \f$\frac{\partial}{\partial\xi} a    \f$
 * @param[out] rho_y  Result \f$\frac{\partial}{\partial{}y} \rho \f$
 * @param[out] u_y    Result \f$\frac{\partial}{\partial{}y} u_\xi\f$
 * @param[out] v_y    Result \f$\frac{\partial}{\partial{}y} u_y  \f$
 * @param[out] p_y    Result \f$\frac{\partial}{\partial{}y} p    \f$
 * @param[out] a_y    Result \f$\frac{\partial}{\partial{}y} a    \f$
 *
 * @see Model documentation in <tt>writeups/radialflow.tex</tt> for details.
 */
void
suzerain_radialflow_cartesian_primitive(
    const suzerain_radialflow_solution * s,
    const size_t i,
    const double Ma,
    double *rho,
    double *u,
    double *v,
    double *p,
    double *a,
    double *rho_xi,
    double *u_xi,
    double *v_xi,
    double *p_xi,
    double *a_xi,
    double *rho_y,
    double *u_y,
    double *v_y,
    double *p_y,
    double *a_y);

/**
 * Compute Cartesian base flow conserved state and pressure at \f$\left(R_0,
 * \delta_i\right)\f$, including streamwise derivatives, given a radial nozzle
 * solution.  This collection of information happens to be what Largo requires
 * for a slow growth baseflow definition per \ref largo_state.
 *
 * @param[in ] s A solution generated by \ref suzerain_radialflow_solver.
 * @param[in ] i An indexed location within \c s in range <tt>[0,s->size)</tt>.
 * @param[in ] Ma    Target Mach number \f$Ma\f$ for returned total energy.
 * @param[out] r     Result \f$\rho      \f$
 * @param[out] ru    Result \f$\rho u_\xi\f$
 * @param[out] rv    Result \f$\rho u_y  \f$
 * @param[out] rE    Result \f$\rho E    \f$
 * @param[out] p     Result \f$p         \f$
 * @param[out] r_xi  Result \f$\frac{\partial}{\partial\xi} \rho      \f$
 * @param[out] ru_xi Result \f$\frac{\partial}{\partial\xi} \rho u_\xi\f$
 * @param[out] rv_xi Result \f$\frac{\partial}{\partial\xi} \rho u_y  \f$
 * @param[out] rE_xi Result \f$\frac{\partial}{\partial\xi} \rho E    \f$
 * @param[out] p_xi  Result \f$\frac{\partial}{\partial\xi} p         \f$
 * @param[out] r_y   Result \f$\frac{\partial}{\partial{}y} \rho      \f$
 * @param[out] ru_y  Result \f$\frac{\partial}{\partial{}y} \rho u_\xi\f$
 * @param[out] rv_y  Result \f$\frac{\partial}{\partial{}y} \rho u_y  \f$
 * @param[out] rE_y  Result \f$\frac{\partial}{\partial{}y} \rho E    \f$
 * @param[out] p_y   Result \f$\frac{\partial}{\partial{}y} p         \f$
 *
 * @see Method \ref suzerain_radialflow_cartesian_primitive for
 *      exactly the primitive state to which this results correspond.
 */
void
suzerain_radialflow_cartesian_conserved(
    const suzerain_radialflow_solution * s,
    const size_t i,
    const double Ma,
    double *r,
    double *ru,
    double *rv,
    double *rE,
    double *p,
    double *r_xi,
    double *ru_xi,
    double *rv_xi,
    double *rE_xi,
    double *p_xi,
    double *r_y,
    double *ru_y,
    double *rv_y,
    double *rE_y,
    double *p_y);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_RADIALFLOW_H */
