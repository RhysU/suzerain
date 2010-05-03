/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
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
 *
 * diffwave.h: Computational kernels for differentiating in wave space
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_DIFFWAVE_H
#define __SUZERAIN_DIFFWAVE_H

/** @file
 * Provides computational kernels for differentiating a three dimensional
 * wavespace field stored (Y, X, Z) in column-major order.  Only
 * differentiation of the X and Z directions is supported.  Additional padding
 * (e.g. due to dealiasing) may be present and will be zeroed during
 * processing.
 */

#ifdef __cplusplus
extern "C" {
#endif

inline
int suzerain_diffwave_freqindex(const int N, const int i)
{
    assert(0 <= i && i < N);
    return (i < N/2+1) ? i : -N + i;
}

inline
int suzerain_diffwave_absfreqindex(const int N, const int i)
{
    assert(0 <= i && i < N);
    return (i < N/2+1) ? i : N - i;
}

inline
int suzerain_diffwave_freqdiffindex(const int N, const int dN, const int i)
{
    assert(0 <= i && i < dN && N <= dN);
    if (i < (N+1)/2) {
        return i;
    } else if (i >= dN - (N-1)/2) {
        return -dN+i;
    } else {
        return 0;
    }
}

void suzerain_diffwave_apply(
    const int dxcnt,
    const int dzcnt,
    const double alpha[2], double (*x)[2],
    const double Lx,
    const double Lz,
    const int Ny,
    const int Nx, const int dNx, const int dkbx, const int dkex,
    const int Nz, const int dNz, const int dkbz, const int dkez);

void suzerain_diffwave_accumulate(
    const int dxcnt,
    const int dzcnt,
    const double alpha[2], const double (*x)[2],
    const double beta[2],        double (*y)[2],
    const double Lx,
    const double Lz,
    const int Ny,
    const int Nx, const int dNx, const int dkbx, const int dkex,
    const int Nz, const int dNz, const int dkbz, const int dkez);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // __SUZERAIN_DIFFWAVE_H
