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
#include <suzerain/blas_et_al.h>
#include <suzerain/diffwave.h>

#define Y_ASSERT(N) \
    do {assert((N) >= 0);} while (0)

#define XZ_ASSERT(N,dN,dkb,dke) \
    do { \
        assert((N) >= 0); \
        assert((dN) >= (N)); \
        assert((dke) <= (dN));  \
        assert((dkb) <= (dke)); \
    } while (0)

void suzerain_diffwave_accumulate_y0x0z0(
    const double alpha[2], const double (* const x)[2],
    const double beta[2],        double (* const y)[2],
    const double Lx,
    const double Lz,
    const int Ny,
    const int Nx, const int dNx, const int dkbx, const int dkex,
    const int Nz, const int dNz, const int dkbz, const int dkez)
{
    SUZERAIN_UNUSED(Lx);
    SUZERAIN_UNUSED(Lz);

    assert((void*) x != (void *)y);
    Y_ASSERT(Ny);
    XZ_ASSERT(Nx, dNx, dkbx, dkex);
    XZ_ASSERT(Nz, dNz, dkbz, dkez);

    // {n,m}keeper complexity because we must not nuke zero and Nyquist modes
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
