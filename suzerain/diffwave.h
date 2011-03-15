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

/**
 * For a one dimensional DFT of length \c N, compute the frequency index
 * associated with "in-order" entry \c i.  That is, the integer \f$i\f$ such
 * that \f$k_i = \frac{2\pi{}i}{L}\f$ for \f$i\in\left\{0,1,\dots,N\right\}\f$
 * ranges over the frequencies supported on a domain of length \f$L\f$.
 *
 * For example, for <tt>N=6</tt> the values <tt>{ 0, 1, 2, 3, -2, -1 }</tt>
 * will be returned for <tt>i = 0, 1, 2, 3, 4, 5</tt>.
 *
 * @param N The length of the DFT.
 * @param i The entry of interest where <tt>0 <= i < N</tt>.
 *
 * @return The frequency index associated with in-order entry <tt>i</tt>.
 * @see FFTW's discussion of <a href="http://www.fftw.org/fftw3.3alpha_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html">The 1d Discrete Fourier Transform</a>
 *      for an extended description of in-order storage.
 */
inline
int suzerain_diffwave_freqindex(const int N, const int i)
{
    assert(N > 0);
    assert(0 <= i && i < N);
    return (i < N/2+1) ? i : -N + i;
}

/**
 * For a one dimensional DFT of length \c N, determine if the given frequency
 * index is supportable by the grid.  For example, for <tt>N=6</tt> the values
 * <tt>{ 0, 1, 2, 3, -2, -1 }</tt> are supportable.  Values less than -2 or
 * greater than 3 are not.
 *
 * @param N The length of the DFT.
 * @param i The frequency index of interest.
 *
 * @return True if the frequency index is supportable.  False otherwise.
 * @see FFTW's discussion of <a href="http://www.fftw.org/fftw3.3alpha_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html">The 1d Discrete Fourier Transform</a>
 *      for an extended description of in-order storage.
 */
inline
int suzerain_diffwave_freqindexsupported(const int N, const int i)
{
    assert(N > 0);
    return (-N+1)/2 <= i && i <= N/2;
}

/**
 * For a one dimensional DFT of length \c N, compute the "in-order" entry
 * associated with the given frequency index \c i.  This is the inverse
 * function of suzerain_diffwave_freqindex().
 *
 * For example, for <tt>N=6</tt> the values <tt>i = 0, 1, 2, 3, 4, 5</tt>
 * will be returned for <tt>{ 0, 1, 2, 3, -2, -1 }</tt>
 *
 * @param N The length of the DFT.
 * @param i The frequency index of interest.
 *
 * @return The in-order entry index associated with the frequency <tt>i</tt>
 * @see FFTW's discussion of <a href="http://www.fftw.org/fftw3.3alpha_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html">The 1d Discrete Fourier Transform</a>
 *      for an extended description of in-order storage.
 */
inline
int suzerain_diffwave_indexfreq(const int N, const int i)
{
    assert(N > 0);
    assert(suzerain_diffwave_freqindexsupported(N,i));
    return (i >= 0) ? i : N + i;
}

/**
 * Compute the absolute value of suzerain_diffwave_freqindex().  This can be
 * computed at reduced cost relative to taking the absolute value of
 * suzerain_diffwave_freqindex().
 *
 * For example, for <tt>N=6</tt> the values <tt>{ 0, 1,  2,  3, 2, 1 }</tt>
 * will be returned for <tt>i = 0, 1, 2, 3, 4, 5</tt>.
 *
 * @param N The length of the DFT.
 * @param i The entry of interest where <tt>0 <= i < N</tt>.
 *
 * @return The absolute value of the frequency index associated with in-order
 *         entry <tt>i</tt>.
 */
inline
int suzerain_diffwave_absfreqindex(const int N, const int i)
{
    assert(N > 0);
    assert(0 <= i && i < N);
    return (i < N/2+1) ? i : N - i;
}

/**
 * Compute the frequency index of the scaling factor used to differentiate
 * a wave-space signal on a dealiased grid of length <tt>dN</tt>.  This
 * method returns zero for modes which cannot be supported on a grid of
 * length <tt>N</tt> or which are removed by differentiation.
 *
 * For example, for <tt>dN=9</tt> and <tt>N=6</tt> the values <tt> {0, 1, 2, 0,
 * 0, 0, 0, -2, -1}</tt> will be returned for <tt>i=0, 1, 2, 3, 4, 5, 6, 7,
 * 8</tt>.
 *
 * @param N The length of the DFT.
 * @param dN The dealiased length of the DFT.
 * @param i The entry of interest where <tt>0 <= i < dN</tt>.
 *
 * @return The frequency index of the scaling factor
 *         associated with in-order entry <tt>i</tt>.
 */
inline
int suzerain_diffwave_freqdiffindex(const int N, const int dN, const int i)
{
    assert(N  > 0);
    assert(dN > 0);
    assert(0 <= i && i < dN && N <= dN);
    if (i < (N+1)/2) {
        return i;
    } else if (i >= dN - (N-1)/2) {
        return -dN+i;
    } else {
        return 0;
    }
}

/**
 * Given a complex-valued, wave-space field \f$x\f$, compute \f$ x \leftarrow{}
 * \alpha \partial{}x^{\mbox{dxcnt}} \partial{}z^{\mbox{dzcnt}} x\f$.  The
 * implementation accounts for the field potentially being dealiased,
 * distributed across multiple machines, and representing a domain of arbitrary
 * length in the X and Z directions.
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
    const double alpha[2], double (*x)[2],
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
