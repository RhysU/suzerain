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

#ifndef SUZERAIN_DIFFWAVE_H
#define SUZERAIN_DIFFWAVE_H

/** @file
 * Computational kernels for differentiating in wave space.
 *
 * Provides computational kernels for differentiating a three dimensional
 * wavespace field stored (Y, X, Z) in column-major order.  Only
 * differentiation of the X and Z directions is supported.  Additional padding
 * (e.g. due to dealiasing) may be present and will be zeroed during
 * processing.
 */

#include <suzerain/complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Given a complex-valued, wave-space field \f$x\f$, compute \f$ x \leftarrow{}
 * \alpha \partial{}x^{\mbox{dxcnt}} \partial{}z^{\mbox{dzcnt}} x\f$.  The
 * implementation accounts for the field potentially being dealiased,
 * distributed across multiple machines, and representing a domain of arbitrary
 * length in the X and Z directions.  This routine may also be used to
 * efficiently zero higher wavenumbers used only for dealiasing when <tt>dxcnt
 * == 0</tt>, <tt>dzcnt == 0</tt> (where usually \f$\alpha = 1\f$).  All
 * Nyquist modes \e are zeroed when <tt>dxcnt == 0</tt>, <tt>dzcnt == 0</tt>.
 *
 * The input and output data \c x is stored column-major over the Y direction
 * (index range <tt>0</tt> to <tt>Ny-1</tt>), X direction (index range
 * <tt>dkbx</tt> to <tt>dkex</tt>), and Z direction (index range <tt>dkbz</tt>
 * to <tt>dkez</tt>).  This layout is equivalent to Dmitry Pekurovsky's P3DFFT
 * storage order when Y is specified to be STRIDE1 in wave space.  Complex
 * values are stored as C arrays of length two with the real part preceding
 * the imaginary part.
 *
 * @param[in]     dxcnt Partial derivative order to compute in the X direction
 * @param[in]     dzcnt Partial derivative order to compute in the Z direction
 * @param[in]     alpha Complex-valued scaling factor \f$\alpha\f$
 * @param[in,out] x     Input and output wave-space field
 * @param[in]     Lx    Length of the domain in the X direction
 * @param[in]     Lz    Length of the domain in the Y direction
 * @param[in]     Ny    Number of points in the Y direction
 * @param[in]     Nx    Number of points in the X direction, which
 *                      determines the maximum wavenumbers which are
 *                      retained when differentiating.
 * @param[in]     dNx   Number of dealiased points in the X direction,
 *                      which determines of offsets are translated into
 *                      frequencies.
 * @param[in]     dkbx  The first (inclusive) in-order frequency contained
 *                      in field \c x in the X direction.
 * @param[in]     dkex  The last (exclusive) in-order frequency contained
 *                      in field \c x in the X direction.
 * @param[in]     Nz    Number of points in the Z direction, which
 *                      determines the maximum wavenumbers which are
 *                      retained when differentiating.
 * @param[in]     dNz   Number of dealiased points in the Z direction,
 *                      which determines of offsets are translated into
 *                      frequencies.
 * @param[in]     dkbz  The first (inclusive) in-order frequency contained
 *                      in field \c z in the Z direction.
 * @param[in]     dkez  The last (exclusive) in-order frequency contained
 *                      in field \c z in the Z direction.
 */
void suzerain_diffwave_apply(
    const int dxcnt,
    const int dzcnt,
    const complex_double alpha, complex_double *x,
    const double Lx,
    const double Lz,
    const int Ny,
    const int Nx, const int dNx, const int dkbx, const int dkex,
    const int Nz, const int dNz, const int dkbz, const int dkez);

/**
 * Given two complex-valued, wave-space fields \f$x\f$ and \f$y\f$, compute \f$
 * y \leftarrow{} \alpha \partial{}x^{\mbox{dxcnt}} \partial{}z^{\mbox{dzcnt}}
 * x + \beta{}y\f$.  The implementation accounts for the field potentially
 * being dealiased, distributed across multiple machines, and representing a
 * domain of arbitrary length in the X and Z directions.
 *
 * @param[in]     dxcnt Partial derivative order to compute in the X direction
 * @param[in]     dzcnt Partial derivative order to compute in the Z direction
 * @param[in]     alpha Complex-valued scaling factor \f$\alpha\f$
 * @param[in]     x     Input wave-space field to be differentiated
 * @param[in]     beta  Complex-valued scaling factor \f$\beta\f$
 * @param[in,out] y     Input and output wave-space field where
 *                      accumulation takes place.
 * @param[in]     Lx    Length of the domain in the X direction
 * @param[in]     Lz    Length of the domain in the Y direction
 * @param[in]     Ny    Number of points in the Y direction
 * @param[in]     Nx    Number of points in the X direction, which
 *                      determines the maximum wavenumbers which are
 *                      retained when differentiating.
 * @param[in]     dNx   Number of dealiased points in the X direction,
 *                      which determines of offsets are translated into
 *                      frequencies.
 * @param[in]     dkbx  The first (inclusive) in-order frequency contained
 *                      in field \c x in the X direction.
 * @param[in]     dkex  The last (exclusive) in-order frequency contained
 *                      in field \c x in the X direction.
 * @param[in]     Nz    Number of points in the Z direction, which
 *                      determines the maximum wavenumbers which are
 *                      retained when differentiating.
 * @param[in]     dNz   Number of dealiased points in the Z direction,
 *                      which determines of offsets are translated into
 *                      frequencies.
 * @param[in]     dkbz  The first (inclusive) in-order frequency contained
 *                      in field \c z in the Z direction.
 * @param[in]     dkez  The last (exclusive) in-order frequency contained
 *                      in field \c z in the Z direction.
 *
 * @see suzerain_diffwave_apply() for more information on the storage layout
 *      for field \c x.  Field \c y has storage requirements identical to
 *      those of \c x.
 */
void suzerain_diffwave_accumulate(
    const int dxcnt,
    const int dzcnt,
    const complex_double alpha, const complex_double *x,
    const complex_double beta,        complex_double *y,
    const double Lx,
    const double Lz,
    const int Ny,
    const int Nx, const int dNx, const int dkbx, const int dkex,
    const int Nz, const int dNz, const int dkbz, const int dkez);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // SUZERAIN_DIFFWAVE_H
