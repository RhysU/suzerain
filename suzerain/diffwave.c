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

#include <assert.h>
#include <suzerain/diffwave.h>

inline
void y_assert(const int Ny)
{
    assert(Ny >= 0);
}

inline
void x_assert(const int N, const int dN, const int dkb, const int dke)
{
    assert(N >= 0);
    assert(dN >= N);
    assert(dke <= dN/2);
    assert(dkb <= dke);
}

inline
void z_assert(const int N, const int dN, const int dkb, const int dke)
{
    assert(N >= 0);
    assert(dN >= N);
    assert(dke <= dN);
    assert(dkb <= dke);
}

int suzerain_diffwave_accumulate_y0x0z0(
    const double alpha, const double * restrict x,
    const double beta,  double * restrict y,
    const double Lx,
    const double Ly,
    const int Ny,
    const int Nx, const int dNx, const int dkbx, const int dkex,
    const int Nz, const int dNz, const int dkbz, const int dkez)
{
    assert((void*) x != (void *)y);
    y_assert(Ny);
    x_assert(Nx, dNx, dkbx, dkex);
    z_assert(Nz, dNz, dkbz, dkez);

}
