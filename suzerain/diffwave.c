/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * diffwave.c: Computational kernels for differentiating in wave space
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.h>
#pragma hdrstop
#include <gsl/gsl_sf_pow_int.h>
#include <suzerain/inorder.h>
#include <suzerain/diffwave.h>

static inline
double _Complex scale_by_imaginary_power(const double in[2], int p)
{
    // Modulo-four-like operation for 2's complement p
    switch (p & 3) {
        default: // unreachable default silences warnings
        case 0:  return  in[0] + in[1]*_Complex_I; // I^0 =  1
        case 1:  return -in[1] + in[0]*_Complex_I; // I^1 =  I = I^-3
        case 2:  return -in[0] - in[1]*_Complex_I; // I^2 = -1 = I^-2
        case 3:  return  in[1] - in[0]*_Complex_I; // I^3 = -I = I^-1
    }
}

#pragma float_control(precise, on)
#pragma fenv_access(on)
#pragma float_control(except, on)
#pragma fp_contract(off)
static inline
double twopiover(const double L)
{
    return 2*M_PI/L;
}
#pragma float_control(except, off)
#pragma fenv_access(off)
#pragma float_control(precise, off)
#pragma fp_contract(on)

// Specialization of BLAS zscal for exactly our needs
static inline
void internal_zscal(const int n,
                    const double _Complex alpha,
                    double _Complex * const restrict x)
{
    for (int i = 0; i < n; ++i)
        x[i] *= alpha;
}

// Specialization of BLAS zaxpby for exactly our needs
static inline
void internal_zaxpby(const int n,
                     const double _Complex alpha,
                     const double _Complex * const restrict x,
                     const double _Complex beta,
                     double _Complex * const restrict y)
{
    for (int i = 0; i < n; ++i)
        y[i] = alpha*x[i] + beta*y[i];
}

// Special case for diffwave_apply when dxcnt = dzcnt = 0 and alpha = 1
static void zero_wavenumbers_used_only_for_dealiasing(
    double (* const x)[2],
    const int Ny,
    const int Nx, const int dNx, const int dkbx, const int dkex,
    const int Nz, const int dNz, const int dkbz, const int dkez);

void suzerain_diffwave_apply(
    const int dxcnt,
    const int dzcnt,
    const double alpha[2], double (* const x)[2],
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
    assert(dkex <= dNx);            // Often dkex <= (dNx/2+1)
    assert(dkbx <= dkex);
    assert(Nz >= 0);                // Sanity check Z direction
    assert(dNz >= Nz);
    assert(dkez <= dNz);
    assert(dkbz <= dkez);

    // Special logic for common-but-degenerate case
#pragma warning(push,disable:1572)
    if (dxcnt == 0 && dzcnt == 0 && alpha[0] == 1 && alpha[1] == 0) {
#pragma warning(pop)
        return zero_wavenumbers_used_only_for_dealiasing(
                x, Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
    }

    // Compute loop independent constants
    const double twopioverLx = twopiover(Lx);  // Weird looking for FP control
    const double twopioverLz = twopiover(Lz);  // Weird looking for FP control
    const double _Complex alpha_ipow
            = scale_by_imaginary_power(alpha, dxcnt + dzcnt);

    // Compute X, Z strides
    const int sx = Ny, sz = (dkex - dkbx)*sx;

    // Overwrite x with alpha*D*x for storage, dealiasing assumptions
    for (int n = dkbz; n < dkez; ++n) {
        const int noff     = sz*(n - dkbz);
        const int nfreqidx = suzerain_inorder_wavenumber_diff(Nz, dNz, n);
        const int nkeeper  = (dzcnt > 0)
                ? nfreqidx
                : suzerain_inorder_wavenumber_abs(dNz, n) <= (Nz-1)/2;
        if (nkeeper) {
            // Relies on gsl_sf_pow_int(0.0, 0) == 1.0
            const double nscale = gsl_sf_pow_int(twopioverLz*nfreqidx, dzcnt);
            for (int m = dkbx; m < dkex; ++m) {
                const int moff     = noff + sx*(m - dkbx);
                const int mfreqidx
                        = suzerain_inorder_wavenumber_diff(Nx, dNx, m);
                const int mkeeper  = (dxcnt > 0)
                        ? mfreqidx
                        : suzerain_inorder_wavenumber_abs(dNx, m) <= (Nx-1)/2;
                if (mkeeper) {
                    // Relies on gsl_sf_pow_int(0.0, 0) == 1.0
                    const double mscale
                        = nscale*gsl_sf_pow_int(twopioverLx*mfreqidx, dxcnt);
                    internal_zscal(Ny, alpha_ipow*mscale, (void*)&x[moff][0]);
                } else {
                    memset(x+moff, 0, Ny*sizeof(x[0])); // Scale by zero
                }
            }
        } else {
            memset(x+noff, 0, sz*sizeof(x[0]));  // Scale by zero
        }
    }
}

static void zero_wavenumbers_used_only_for_dealiasing(
    double (* const x)[2],
    const int Ny,
    const int Nx, const int dNx, const int dkbx, const int dkex,
    const int Nz, const int dNz, const int dkbz, const int dkez)
{
    assert(Ny >= 0);                // Sanity check Y direction
    assert(Nx >= 0);                // Sanity check X direction
    assert(dNx >= Nx);
    assert(dkex <= dNx);            // Often dkex <= (dNx/2+1)
    assert(dkbx <= dkex);
    assert(Nz >= 0);                // Sanity check Z direction
    assert(dNz >= Nz);
    assert(dkez <= dNz);
    assert(dkbz <= dkez);

    // Compute X, Z strides
    const int sx = Ny, sz = (dkex - dkbx)*sx;

    // Zero only non-dealiased wavenumbers in x using storage assumptions
    // Examine suzerain_diffwave_apply() to help unravel this loop logic
    for (int n = dkbz; n < dkez; ++n) {
        const int noff = sz*(n - dkbz);
        if (suzerain_inorder_wavenumber_abs(dNz, n) <= (Nz-1)/2) {
            for (int m = dkbx; m < dkex; ++m) {
                if (suzerain_inorder_wavenumber_abs(dNx, m) > (Nx-1)/2) {
                    const int moff = noff + sx*(m - dkbx);
                    memset(x + moff, 0, Ny*sizeof(x[0]));
                }
            }
        } else {
            memset(x+noff, 0, sz*sizeof(x[0]));
        }
    }
}

void suzerain_diffwave_accumulate(
    const int dxcnt,
    const int dzcnt,
    const double alpha[2], const double (* const x)[2],
    const double beta[2],        double (* const y)[2],
    const double Lx,
    const double Lz,
    const int Ny,
    const int Nx, const int dNx, const int dkbx, const int dkex,
    const int Nz, const int dNz, const int dkbz, const int dkez)
{
    assert(dxcnt >= 0);             // Only differentiation supported
    assert(dzcnt >= 0);
    assert((void*) x != (void *)y); // Sanity check aliasing requirements
    assert(Ny >= 0);                // Sanity check Y direction
    assert(Nx >= 0);                // Sanity check X direction
    assert(dNx >= Nx);
    assert(dkex <= dNx);            // Often dkex <= (dNx/2+1)
    assert(dkbx <= dkex);
    assert(Nz >= 0);                // Sanity check Z direction
    assert(dNz >= Nz);
    assert(dkez <= dNz);
    assert(dkbz <= dkez);

    // Compute loop independent constants
    const double twopioverLx = twopiover(Lx);  // Weird looking for FP control
    const double twopioverLz = twopiover(Lz);  // Weird looking for FP control
    const double _Complex alpha_ipow = scale_by_imaginary_power(alpha, dxcnt + dzcnt);
    double _Complex beta_c;
    memcpy(&beta_c, beta, sizeof(beta));

    // Compute X, Z strides
    const int sx = Ny, sz = (dkex - dkbx)*sx;

    // Accumulate y <- alpha*D*x + beta*y for storage, dealiasing assumptions
    for (int n = dkbz; n < dkez; ++n) {
        const int noff     = sz*(n - dkbz);
        const int nfreqidx = suzerain_inorder_wavenumber_diff(Nz, dNz, n);
        const int nkeeper  = (dzcnt > 0)
                ? nfreqidx
                : suzerain_inorder_wavenumber_abs(dNz, n) <= (Nz-1)/2;
        if (nkeeper) {
            // Relies on gsl_sf_pow_int(0.0, 0) == 1.0
            const double nscale = gsl_sf_pow_int(twopioverLz*nfreqidx, dzcnt);
            for (int m = dkbx; m < dkex; ++m) {
                const int moff = noff + sx*(m - dkbx);
                const int mfreqidx
                    = suzerain_inorder_wavenumber_diff(Nx, dNx, m);
                const int mkeeper  = (dxcnt > 0)
                        ? mfreqidx
                        : suzerain_inorder_wavenumber_abs(dNx, m) <= (Nx-1)/2;
                if (mkeeper) {
                    // Relies on gsl_sf_pow_int(0.0, 0) == 1.0
                    const double mscale
                        = nscale*gsl_sf_pow_int(twopioverLx*mfreqidx, dxcnt);
                    internal_zaxpby(Ny, mscale*alpha_ipow, (void*)&x[moff][0],
                                        beta_c,            (void*)&y[moff][0]);
                } else {
                    internal_zscal(Ny, beta_c, (void*)&y[moff][0]);
                }
            }
        } else {
            internal_zscal(sz, beta_c, (void*)&y[noff][0]);
        }
    }
}
