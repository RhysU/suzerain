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
#include <suzerain/blas_et_al.h>
#include <suzerain/diffwave.h>

inline
void scale_by_imaginary_power(const double in[2], double out[2], int p)
{
    // Modulo-four-like operation for 2's complement p
    switch (p & 3) {
        case 0: // I^0 = 1
            out[0] = in[0];
            out[1] = in[1];
            break;
        case 1: // I^1 = I = I^-3
            out[0] = -in[1];
            out[1] =  in[0];
            break;
        case 2: // I^2 = -1 = I^-2
            out[0] = -in[0];
            out[1] = -in[1];
            break;
        case 3: // I^3 = -I = I^-1
            out[0] =  in[1];
            out[1] = -in[0];
            break;
    }
}

static
void suzerain_diffwave_accumulate_y0x0z0(
    const double alpha[2], const double (* const x)[2],
    const double beta[2],        double (* const y)[2],
    const int Ny,
    const int Nx, const int dNx, const int dkbx, const int dkex,
    const int Nz, const int dNz, const int dkbz, const int dkez)
{
    // {n,m}keeper complexity because we must include zero and Nyquist modes
    const int sx = Ny, sz = (dkex - dkbx)*sx; // Compute X, Z strides
    for (int n = dkbz; n < dkez; ++n) {
        const int noff = sz*(n - dkbz);
        const int nkeeper =    suzerain_diffwave_freqindex(Nz, dNz, n)
                            || (n == 0) || (!(Nz & 1) && n == Nz/2);
        if (nkeeper) {
            for (int m = dkbx; m < dkex; ++m) {
                const int moff = noff + sx*(m - dkbx);
                const int mkeeper =    suzerain_diffwave_freqindex(Nx, dNx, m)
                                    || (m == 0) || (!(Nx & 1) && m == Nx/2);
                if (mkeeper) {
                    suzerain_blas_zaxpby(Ny,alpha,x+moff,1,beta,y+moff,1);
                } else {
                    suzerain_blas_zscal(Ny,beta,y+moff,1);
                }
            }
        } else {
            suzerain_blas_zscal(sz,beta,y+noff,1);
        }
    }
}

#pragma float_control(precise, on)
#pragma fenv_access(on)
#pragma float_control(except, on)
#pragma fp_contract(off)
inline double twopiover(const double L)
{
    return 2*M_PI/L;
}
#pragma float_control(except, off)
#pragma fenv_access(off)
#pragma float_control(precise, off)
#pragma fp_contract(on)

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

    // No derivatives case is sufficiently weird to warrant special handling
    if (dxcnt == 0 && dzcnt == 0) {
        return suzerain_diffwave_accumulate_y0x0z0(alpha, x, beta,  y,
                                                   Ny,
                                                   Nx, dNx, dkbx, dkex,
                                                   Nz, dNz, dkbz, dkez);
    }

    // Compute loop independent constants
    const double twopioverLx = twopiover(Lx);  // Weird looking for FP control
    const double twopioverLz = twopiover(Lz);  // Weird looking for FP control
    double alpha_ipow[2];
    scale_by_imaginary_power(alpha, alpha_ipow, dxcnt + dzcnt);

    const int sx = Ny, sz = (dkex - dkbx)*sx; // Compute X, Z strides
    for (int n = dkbz; n < dkez; ++n) {
        const int noff = sz*(n - dkbz);
        const int nfreqidx = suzerain_diffwave_freqindex(Nz, dNz, n);
        if (nfreqidx) {
            const double nscale = gsl_sf_pow_int(twopioverLz*nfreqidx, dzcnt);
            for (int m = dkbx; m < dkex; ++m) {
                const int moff = noff + sx*(m - dkbx);
                const int mfreqidx = suzerain_diffwave_freqindex(Nx, dNx, m);
                if (mfreqidx) {
                    const double mscale
                        = nscale*gsl_sf_pow_int(twopioverLx*mfreqidx, dxcnt);
                    const double malpha[2] = { mscale*alpha_ipow[0],
                                               mscale*alpha_ipow[1] };
                    suzerain_blas_zaxpby(Ny,malpha,x+moff,1,beta,y+moff,1);
                } else {
                    suzerain_blas_zscal(Ny,beta,y+moff,1);
                }
            }
        } else {
            suzerain_blas_zscal(sz,beta,y+noff,1);
        }
    }
}
