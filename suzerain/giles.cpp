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

/** @file
 * @copydoc giles.hpp
 */

#include <suzerain/giles.hpp>

#include <suzerain/ndx.hpp>

namespace suzerain {

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
        Matrix5r& PG_VL_S_RY,
        Matrix5r& inv_VL_S_RY,
        const real_t normal_sign)
{
    // Prepare a rotation (differing from the model document), reordering,
    // ndx::{e, mx, my, mz, rho} to {rho = 0, my = 1, mz = 2, mx = 3, e = 4}.
    // All remaining matrices/logic within the routine uses the latter order!
    // The inverse rotation RY^{-1} is accomplished using RY.transpose().
    Matrix5r RY(Matrix5r::Zero());
    {
        RY(0, ndx::rho) = 1;
        RY(1, ndx::my ) = 1;
        RY(2, ndx::mz ) = 1;
        RY(3, ndx::mx ) = 1;
        RY(4, ndx::e  ) = 1;
    }

    // Prepare upper boundary reference state per rotated frame
    const real_t rho = ref_rho;
    const real_t u   = ref_v  ; // Ref u = Ref v' per RY
    const real_t v   = ref_w  ; // Ref v = Ref w' per RY
    const real_t w   = ref_u  ; // Ref w = Ref u' per RY
          real_t a   = ref_a  ;

    // Prepare oft-used quantities derived from the reference state
    const real_t inv_rho = 1 / rho;
    const real_t u2      = u * u;
    const real_t v2      = v * v;
    const real_t w2      = w * w;
          real_t inv_a   = 1 / a;
          real_t a2      = a * a;
          real_t inv_a2  = 1 / a2;

    // Prepare oft-used scenario-related constants
    const real_t half    = static_cast<real_t>(1) / 2;
    const real_t gamma1  = gamma - 1;
    const real_t inv_Ma  = 1 / Ma;
    const real_t Ma2     = Ma * Ma;
    const real_t inv_Ma2 = 1 / Ma2;

    // Build the conserved -> primitive transformation and its inverse
    Matrix5r S(Matrix5r::Zero());
    {
        S(0, 0) =   1;
        S(1, 0) = - u * inv_rho;
        S(2, 0) = - v * inv_rho;
        S(3, 0) = - w * inv_rho;
        S(4, 0) =   gamma1 * (u2 + v2 + w2) / 2;
        S(1, 1) =   inv_rho;
        S(4, 1) = - gamma1 * u;
        S(2, 2) =   inv_rho;
        S(4, 2) = - gamma1 * v;
        S(3, 3) =   inv_rho;
        S(4, 3) = - gamma1 * w;
        S(4, 4) =   gamma1 * inv_Ma2;
    }
    Matrix5r inv_S(Matrix5r::Zero());
    {
        inv_S(0, 0) =   1;
        inv_S(1, 0) =   u;
        inv_S(2, 0) =   v;
        inv_S(3, 0) =   w;
        inv_S(4, 0) =   Ma2 * (u2 + v2 + w2) / 2;
        inv_S(1, 1) =   rho;
        inv_S(4, 1) =   Ma2 * rho * u;
        inv_S(2, 2) =   rho;
        inv_S(4, 2) =   Ma2 * rho * v;
        inv_S(3, 3) =   rho;
        inv_S(4, 3) =   Ma2 * rho * w;
        inv_S(4, 4) =   Ma2 / gamma1;
    }

    // Per the model document
    //   When reusing [V^L, B^G, and C^G] derived for U ... for V^\ast
    //   every sound speed must be scaled by 1 / Ma because
    //   \bar{a} / u_0 = a_0 \bar{a}^\ast / u_0 = \bar{a}^\ast / Ma.
    // so we now once-and-for-all adjust a and a2 appropriately.
    a      *= inv_Ma;
    inv_a  *= Ma;
    a2     *= inv_Ma2;
    inv_a2 *= Ma2;

    // Build the primitive -> characteristic transformation and its inverse
    Matrix5r VL(Matrix5r::Zero());
    {
        VL(0, 0) = - a2;
        VL(3, 1) =   rho * a;
        VL(4, 1) = - rho * a;
        VL(1, 2) =   rho * a;
        VL(2, 3) =   rho * a;
        VL(0, 4) =   1;
        VL(3, 4) =   1;
        VL(4, 4) =   1;
    }
    Matrix5r inv_VL(Matrix5r::Zero());
    {
        inv_VL(0, 0) = - inv_a2;
        inv_VL(2, 1) =   inv_rho * inv_a;
        inv_VL(3, 2) =   inv_rho * inv_a;
        inv_VL(0, 3) =   half * inv_a2;
        inv_VL(1, 3) =   half * inv_rho * inv_a;
        inv_VL(4, 3) =   half;
        inv_VL(0, 4) =   half * inv_a2;
        inv_VL(1, 4) = - half * inv_rho * inv_a;
        inv_VL(4, 4) =   half;
    }

    // Build the in-vs-outflow characteristic-preserving projection
    // automatically accounting for sub- versus supersonic boundaries.
    Vector5r PG;
    PG(0) = static_cast<real_t>(normal_sign * (u    ) < 0);  // Entropy
    PG(1) = static_cast<real_t>(normal_sign * (u    ) < 0);  // Vorticity
    PG(2) = static_cast<real_t>(normal_sign * (u    ) < 0);  // Vorticity
    PG(3) = static_cast<real_t>(normal_sign * (u + a) < 0);  // Pressure
    PG(4) = static_cast<real_t>(normal_sign * (u - a) < 0);  // Pressure

    Matrix5r BG(Matrix5r::Zero());          // Medida equations 5.82 and 5.83
    BG(1, 1) = v;
    BG(3, 1) = PG(0) ? half * (a - u) : u;  // Inflow modifications to...
    BG(4, 1) = PG(0) ? half * (a + u) : u;  // ...produce well-posed result.
    BG(2, 2) = v;
    BG(1, 3) = half * (a + u);
    BG(3, 3) = v;
    BG(1, 4) = half * (a - u);
    BG(4, 4) = v;

    Matrix5r CG(Matrix5r::Zero());          // Medida equations 5.82 and 5.83
    CG(1, 1) = w;
    CG(2, 2) = w;
    CG(3, 2) = PG(0) ? half * (a - u) : u;  // Inflow modifications to...
    CG(4, 2) = PG(0) ? half * (a + u) : u;  // ...produce well-posed result.
    CG(2, 3) = half * (a + u);
    CG(3, 3) = w;
    CG(2, 4) = half * (a - u);
    CG(4, 4) = w;

    // Prepare all necessary real-valued products of the above matrices
    VL_S_RY       = VL * S * RY;
    PG_BG_VL_S_RY = PG.asDiagonal() * BG * VL_S_RY;
    PG_CG_VL_S_RY = PG.asDiagonal() * CG * VL_S_RY;
    PG_VL_S_RY    = PG.asDiagonal()      * VL_S_RY;
    inv_VL_S_RY   = RY.transpose() * inv_S * inv_VL;  // Symbolic inverse
////inv_VL_S_RY   = VL_S_RY.inverse();                // Numeric inverse
}

} // namespace suzerain
