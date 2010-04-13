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

#include <suzerain/common.h>
#include <suzerain/diffwave.h>

inline
void y_assert(const int Ny)
{
    assert(Ny >= 0);
}

inline
void xz_assert(const int N, const int dN, const int dkb, const int dke)
{
    assert(N >= 0);
    assert(dN >= N);
    assert(dke <= dN);  // Usually dkex <= dNx/2+1 for X, but not required.
    assert(dkb <= dke);
}

void suzerain_diffwave_accumulate_y0x0z0(
    const double alpha, const double * SUZERAIN_RESTRICT x,
    const double beta,  double * SUZERAIN_RESTRICT y,
    const double Lx,
    const double Lz,
    const int Ny,
    const int Nx, const int dNx, const int dkbx, const int dkex,
    const int Nz, const int dNz, const int dkbz, const int dkez)
{
    SUZERAIN_UNUSED(Lx);
    SUZERAIN_UNUSED(Lz);

    assert((void*) x != (void *)y);
    y_assert(Ny);
    xz_assert(Nx, dNx, dkbx, dkex);
    xz_assert(Nz, dNz, dkbz, dkez);

    // {n,m}keeper complexity because we must not nuke zero and Nyquist modes
    for (int n = dkbz; n < dkez; ++n) {
        const int nkeeper =     (n == 0)
                             || (!(Nz & 1) && n == Nz/2)
                             || suzerain_diffwave_freqindex(Nz, dNz, n);
        if (nkeeper) {
            for (int m = dkbx; m < dkex; ++m) {
                const int mkeeper =    (m == 0)
                                    || (!(Nx & 1) && m == Nx/2)
                                    || suzerain_diffwave_freqindex(Nx, dNx, m);
                if (mkeeper) {
                    for (int l = 0; l < Ny; ++l) {
                        *y++ = beta*(*y) + alpha*(*x++); // Real
                        *y++ = beta*(*y) + alpha*(*x++); // Imag
                    }
                } else {
                    const size_t fillcount = Ny*2;
                    memset(y, 0, fillcount*sizeof(y[0]));
                    y += fillcount;
                    x += fillcount;
                }
            }
        } else {
            const size_t fillcount = (dkex-dkbx)*Ny*2;
            memset(y, 0, fillcount*sizeof(y[0]));
            y += fillcount;
            x += fillcount;
        }
    }
}
