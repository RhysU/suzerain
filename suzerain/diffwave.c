/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2014 Rhys Ulerich
 * Copyright (C) 2012-2014 The PECOS Development Team
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

/** @file
 * @copydoc diffwave.h
 */

#include <suzerain/diffwave.h>

#include <gsl/gsl_sf_pow_int.h>

#include <suzerain/common.h>
#include <suzerain/blas_et_al.h>
#include <suzerain/inorder.h>

static inline
void scale_by_imaginary_power(const complex_double in,
                              complex_double * const out, int p)
{
    // Modulo-four-like operation for 2's complement p
    switch (p & 3) {
    case 0: *out =  in           ; break; // I^0 = 1
    case 1: *out =  in*_Complex_I; break; // I^1 = I = I^-3
    case 2: *out = -in           ; break; // I^2 = -1 = I^-2
    case 3: *out = -in*_Complex_I; break; // I^3 = -I = I^-1
    }
}

#pragma float_control(precise, on)
#pragma fenv_access(on)
#pragma float_control(except, on)
#pragma fp_contract(off)
static inline
double twopiover(const double L)
{
    return 2
         * 3.1415926535897932384626433832795028841971693993751058209
         / L;
}
#pragma float_control(except, off)
#pragma fenv_access(off)
#pragma float_control(precise, off)
#pragma fp_contract(on)

void suzerain_diffwave_apply(
    const int dxcnt,
    const int dzcnt,
    const complex_double alpha, complex_double * const x,
    const double Lx,
    const double Lz,
    const int Ny,
    const int Nx, const int dNx, const int dkbx, const int dkex,
    const int Nz, const int dNz, const int dkbz, const int dkez)
{
    assert(dxcnt >= 0);             // Only differentiation supported
    assert(dzcnt >= 0);
    assert(Ny >= 0);                // Sanity check Y direction
    assert(Nx >= 0);                // Sanity check X direction
    assert(dNx >= Nx);
    assert(dkex <= dNx);
    assert(dkbx <= dkex);
    assert(Nz >= 0);                // Sanity check Z direction
    assert(dNz >= Nz);
    assert(dkez <= dNz);
    assert(dkbz <= dkez);

    // Compute X, Z strides based on locally contiguous storage assumption
    const int sx = Ny, sz = (dkex - dkbx)*sx;

    // Precompute loop independent constants
    const int absmin_wm = suzerain_inorder_wavenumber_absmin(Nx);
    const int absmin_wn = suzerain_inorder_wavenumber_absmin(Nz);
    const double twopioverLx = twopiover(Lx);  // Weird looking for FP control
    const double twopioverLz = twopiover(Lz);  // Weird looking for FP control
    complex_double alpha_ipow;
    scale_by_imaginary_power(alpha, &alpha_ipow, dxcnt + dzcnt);

    // Overwrite x with alpha*D*x for storage, dealiasing details
    // Generalized loops are hideous but permit one linear pass through memory
    int nfreqidx = INT_MAX;  // Hoisted out of loops to avoid uninitialized
    int mfreqidx = INT_MAX;  // usage warnings.  Convince yourself it's OK.
    for (int n = dkbz; n < dkez; ++n) {
        const int noff = sz*(n - dkbz);
        const int nkeeper  = (dzcnt > 0)
            ? (nfreqidx = suzerain_inorder_wavenumber_diff(Nz, dNz, n))
            : suzerain_inorder_wavenumber_abs(dNz, n) <= absmin_wn;
        if (nkeeper) {
            // Relies on gsl_sf_pow_int(0.0, 0) == 1.0
            const double nscale = gsl_sf_pow_int(twopioverLz*nfreqidx, dzcnt);
            for (int m = dkbx; m < dkex; ++m) {
                const int moff = noff + sx*(m - dkbx);
                const int mkeeper = (dxcnt > 0)
                    ? (mfreqidx = suzerain_inorder_wavenumber_diff(Nx, dNx, m))
                    : suzerain_inorder_wavenumber_abs(dNx, m) <= absmin_wm;
                if (mkeeper) {
                    // Relies on gsl_sf_pow_int(0.0, 0) == 1.0
                    const double mscale
                        = nscale*gsl_sf_pow_int(twopioverLx*mfreqidx, dxcnt);
                    const complex_double malpha = mscale*alpha_ipow;
                    suzerain_blas_zscal(Ny, malpha, x+moff, 1);
                } else {
                    memset(x+moff, 0, Ny*sizeof(x[0])); // Scale by zero
                }
            }
        } else {
            memset(x+noff, 0, sz*sizeof(x[0]));  // Scale by zero
        }
    }
}

void suzerain_diffwave_accumulate(
    const int dxcnt,
    const int dzcnt,
    const complex_double alpha, const complex_double * const x,
    const complex_double beta,        complex_double * const y,
    const double Lx,
    const double Lz,
    const int Ny,
    const int Nx, const int dNx, const int dkbx, const int dkex,
    const int Nz, const int dNz, const int dkbz, const int dkez)
{
    assert(dkex == dkbx || dkez == dkbz || x != y); // Trivial or not aliased?
    assert(dxcnt >= 0);                             // No integration!
    assert(dzcnt >= 0);
    assert(Ny >= 0);                                // Sanity check Y direction
    assert(Nx >= 0);                                // Sanity check X direction
    assert(dNx >= Nx);
    assert(dkex <= dNx);
    assert(dkbx <= dkex);
    assert(Nz >= 0);                                // Sanity check Z direction
    assert(dNz >= Nz);
    assert(dkez <= dNz);
    assert(dkbz <= dkez);

    // Compute X, Z strides based on locally contiguous storage assumption
    const int sx = Ny, sz = (dkex - dkbx)*sx;

    // Precompute loop independent constants
    const int absmin_wm = suzerain_inorder_wavenumber_absmin(Nx);
    const int absmin_wn = suzerain_inorder_wavenumber_absmin(Nz);
    const double twopioverLx = twopiover(Lx);  // Weird looking for FP control
    const double twopioverLz = twopiover(Lz);  // Weird looking for FP control
    complex_double alpha_ipow;
    scale_by_imaginary_power(alpha, &alpha_ipow, dxcnt + dzcnt);

    // Accumulate y <- alpha*D*x + beta*y for storage, dealiasing details
    // Generalized loops are hideous but permit one linear pass through memory
    int nfreqidx = INT_MAX;  // Hoisted out of loops to avoid uninitialized
    int mfreqidx = INT_MAX;  // usage warnings.  Convince yourself it's OK.
    for (int n = dkbz; n < dkez; ++n) {
        const int noff = sz*(n - dkbz);
        const int nkeeper = (dzcnt > 0)
            ? (nfreqidx = suzerain_inorder_wavenumber_diff(Nz, dNz, n))
            : suzerain_inorder_wavenumber_abs(dNz, n) <= absmin_wn;
        if (nkeeper) {
            // Relies on gsl_sf_pow_int(0.0, 0) == 1.0
            const double nscale = gsl_sf_pow_int(twopioverLz*nfreqidx, dzcnt);
            for (int m = dkbx; m < dkex; ++m) {
                const int moff = noff + sx*(m - dkbx);
                const int mkeeper  = (dxcnt > 0)
                    ? (mfreqidx = suzerain_inorder_wavenumber_diff(Nx, dNx, m))
                    : suzerain_inorder_wavenumber_abs(dNx, m) <= absmin_wm;
                if (mkeeper) {
                    // Relies on gsl_sf_pow_int(0.0, 0) == 1.0
                    const double mscale
                        = nscale*gsl_sf_pow_int(twopioverLx*mfreqidx, dxcnt);
                    const complex_double malpha = mscale*alpha_ipow;
                    suzerain_blas_zaxpby(Ny, malpha, x+moff, 1, beta, y+moff, 1);
                } else {
                    suzerain_blas_zscal(Ny, beta, y+moff, 1);
                }
            }
        } else {
            suzerain_blas_zscal(sz,beta,y+noff,1);
        }
    }
}
