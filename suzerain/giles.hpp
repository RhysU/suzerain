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

#ifndef SUZERAIN_GILES_HPP
#define SUZERAIN_GILES_HPP

/** @file
 * Forms Medida's Giles-like nonreflecting boundary condition matrices.
 * Background and notation set within <tt>writeups/perfect_gas.tex</tt> section
 * entitled "Nonreflecting freestream boundary conditions".
 */

#include <suzerain/common.hpp>

namespace suzerain {

/**
 * Compute Medida's Giles-like boundary condition matrices ordered per \ref
 * suzerain::ndx for a nonreflecting \f$y\f$ boundary given reference state.
 * Inflow vs outflow and subsonic vs supersonic conditions are determined using
 * \c ref_v and \c normal_sign.  All matrices are independent of wavenumber in a
 * wavespace context.
 *
 * @param[in]  Ma      Reference Mach number \f$\mbox{Ma} = u_0 / a_0\f$
 *                     if matrices are to be used in a nondimensional context.
 *                     To use them dimensionally, provide the value one.
 * @param[in]  gamma   Reference ratio of specific heats.
 * @param[in]  ref_rho Reference density \f$rho\f$.
 * @param[in]  ref_u   Reference streamwise velocity \f$u\f$.
 * @param[in]  ref_v   Reference velocity \f$v\f$.
 * @param[in]  ref_w   Reference spanwise velocity \f$w\f$.
 * @param[in]  ref_a   Reference sound speed \f$a\f$.
 * @param[out] VL_S_RY       \f$\left[V^L S\right] R^Y\f
 * @param[out] PG_BG_VL_S_RY \f$\left[B^G\right] \left[V^L S\right] R^Y\f$
 * @param[out] PG_CG_VL_S_RY \f$\left[C^G\right] \left[V^L S\right] R^Y\f$
 * @param[out] ImPG_VL_S_RY  \f$\left(I - P^G\right) \left[V^L S\right] R^Y\f$
 * @param[out] inv_VL_S_RY   \f${R^Y}^{-1} \left[V^L S\right]^{-1}\f$
 * @param[in]  normal_sign Sign of a boundary-normal vector.
 *                    For the boundary at \f$y=0\f$ this should be negative.
 *                    For the boundary at \f$y=L_y\f$, it must be positive.
 */
void
giles_matrices(
        const real_t Ma,
        const real_t gamma,
        const real_t ref_rho,
        const real_t ref_u,
        const real_t ref_v,
        const real_t ref_w,
        const real_t ref_a,
        Matrix5r& VL_S_RY,
        Matrix5r& PG_BG_VL_S_RY,
        Matrix5r& PG_CG_VL_S_RY,
        Matrix5r& ImPG_VL_S_RY,
        Matrix5r& inv_VL_S_RY,
        const real_t normal_sign);

/**
 * Invokes giles_matrices() for lower \f$y=L_y\f$ boundary.
 */
inline
void
giles_matrices_lower(
        const real_t Ma,
        const real_t gamma,
        const real_t ref_rho,
        const real_t ref_u,
        const real_t ref_v,
        const real_t ref_w,
        const real_t ref_a,
        Matrix5r& VL_S_RY,
        Matrix5r& PG_BG_VL_S_RY,
        Matrix5r& PG_CG_VL_S_RY,
        Matrix5r& ImPG_VL_S_RY,
        Matrix5r& inv_VL_S_RY)
{
    return giles_matrices(Ma, gamma, ref_rho, ref_u, ref_v, ref_w, ref_a,
            VL_S_RY, PG_BG_VL_S_RY, PG_CG_VL_S_RY, ImPG_VL_S_RY, inv_VL_S_RY,
            -1);
}

/**
 * Invokes giles_matrices() for upper \f$y=0\f$ boundary.
 */
inline
void
giles_matrices_upper(
        const real_t Ma,
        const real_t gamma,
        const real_t ref_rho,
        const real_t ref_u,
        const real_t ref_v,
        const real_t ref_w,
        const real_t ref_a,
        Matrix5r& VL_S_RY,
        Matrix5r& PG_BG_VL_S_RY,
        Matrix5r& PG_CG_VL_S_RY,
        Matrix5r& ImPG_VL_S_RY,
        Matrix5r& inv_VL_S_RY)
{
    return giles_matrices(Ma, gamma, ref_rho, ref_u, ref_v, ref_w, ref_a,
            VL_S_RY, PG_BG_VL_S_RY, PG_CG_VL_S_RY, ImPG_VL_S_RY, inv_VL_S_RY,
            +1);
}


} // namespace suzerain

#endif  /* SUZERAIN_GILES_HPP */
