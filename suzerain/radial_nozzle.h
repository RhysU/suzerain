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
    double rho;     /**< Density                    \f$\rho           \f$ */
    double u;       /**< Radial velocity            \f$u              \f$ */
    double p;       /**< Pressure                   \f$p              \f$ */
    double a2;      /**< Sound speed squared        \f$a^2            \f$ */
    double rhop;    /**< Density derivative         \f$\partial_R \rho\f$ */
    double up;      /**< Radial velocity derivative \f$\partial_R u   \f$ */
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
    suzerain_radial_nozzle_state state[0]; /**< Pointwise solution details */
} suzerain_radial_nozzle_solution;

/**
 * Solve the requested radial nozzle problem.
 *
 * @param Ma0   Reference Mach number         \f$\mbox{Ma}_0\f$
 * @param gam0  Reference specific heat ratio \f$\gamma_0   \f$
 * @param rho1  Initial density               \f$\rho\!\left(R_1\right)\f$
 * @param u1    Initial radial velocity       \f$u   \!\left(R_1\right)\f$
 * @param p1    Initial pressure              \f$p   \!\left(R_1\right)\f$
 * @param R     Radii of interest with $R_1$ taken from \c R[0].
 *              Must be a contiguous array of length \c size.
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

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_RADIAL_NOZZLE_H */
