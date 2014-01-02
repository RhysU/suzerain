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

#ifndef SUZERAIN_INORDER_H
#define SUZERAIN_INORDER_H

#include <assert.h>

/** @file
 * Provides utilities for manipulating one dimensional discrete Fourier
 * transform coefficients stored "in-order".  From
 * href="http://www.fftw.org/fftw3_doc/One_002dDimensional-DFTs-of-Real-Data.html">FFTW's
 * documentation</a>:
 *     "...the <tt>k</tt>-th output corresponds to the frequency
 *     <tt>k/n</tt> (or <tt>k/T</tt>, where <tt>T</tt> is your total
 *     sampling period). For those who like to think in terms of
 *     positive and negative frequencies, this means that the positive
 *     frequencies are stored in the first half of the output and the
 *     negative frequencies are stored in backwards order in the second
 *     half of the output. (The frequency <tt>-k/n</tt> is the same as
 *     the frequency <tt>(n-k)/n</tt>.)"
 *
 * The following concepts are used throughout these routine descriptions.  They
 * are defined assuming one is working with a one-dimensional discrete Fourier
 * transform (DFT) data of length <tt>N</tt>:
 * <dl>
 *   <dt>index</dt>
 *   <dd>
 *      The indices are the nonnegative integers <tt>0</tt>, <tt>1</tt>, ...,
 *      <tt>N-1</tt>.  These are the usual C indices used to access an array of
 *      length <tt>N</tt>.
 *   </dd>
 *   <dt>wavenumber</dt>
 *   <dd>
 *      The wavenumbers are the integers <tt>(-N+1)/2</tt>, ..., <tt>-1</tt>,
 *      <tt>0</tt>, <tt>1</tt>, ..., <tt>N/2</tt> such that \f$k_i =
 *      \frac{2\pi{}i}{L}\f$ ranges over the frequencies supported on a domain
 *      of length \f$L\f$.
 *   </dd>
 *   <dt>range</dt>
 *   <dd>
 *      A range <tt>[a,b)</tt> contains the integers <tt>a</tt>, <tt>a+1</tt>,
 *      ..., <tt>b-2</tt>, <tt>b-1</tt> assuming <tt>b > a</tt>.  For <tt>b ==
 *      a</tt> the range contains nothing and is said to be empty.
 *   </dd>
 * </dl>
 * To make these concepts concrete, for <tt>N = 8</tt> the indices run from
 * <tt>0</tt> to <tt>7</tt>, inclusive.  This is index range <tt>[0,8)</tt>.
 * The wavenumbers corresponding to these indices are <tt>0</tt>, <tt>1</tt>,
 * <tt>2</tt>, <tt>3</tt>, <tt>4</tt>, <tt>-3</tt>, <tt>-2</tt>, <tt>-1</tt>.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * For a length <tt>N</tt> DFT, find the wavenumber associated with index
 * <tt>i</tt>.
 *
 * For example, for <tt>N=6</tt> the values <tt>{ 0, 1, 2, 3, -2, -1 }</tt>
 * will be returned for <tt>i = 0, 1, 2, 3, 4, 5</tt>.
 *
 * @param N The length of the DFT.
 * @param i The index of interest.
 *
 * @return The wavenumber associated with index <tt>i</tt>.
 * @see The documentation for inorder.h for terminology details.
 */
inline
int suzerain_inorder_wavenumber(const int N, const int i)
{
    assert(0 <= i && i < N);
    return (i < N/2+1) ? i : -N + i;
}

/**
 * For a length <tt>N</tt> DFT, find the absolute value of the wavenumber
 * associated with index <tt>i</tt>.  This can be computed at reduced cost
 * relative to taking the absolute value of suzerain_inorder_wavenumber().
 *
 * For example, for <tt>N=6</tt> the values <tt>{ 0, 1, 2, 3, 2, 1 }</tt>
 * will be returned for <tt>i = 0, 1, 2, 3, 4, 5</tt>.
 *
 * @param N The length of the DFT.
 * @param i The index of interest.
 *
 * @return The absolute value of the wavenumber for index <tt>i</tt>.
 * @see The documentation for inorder.h for terminology details.
 */
inline
int suzerain_inorder_wavenumber_abs(const int N, const int i)
{
    assert(0 <= i && i < N);
    return (i < N/2+1) ? i : N - i;
}

/**
 * Find the minimum wavenumber contained in a length <tt>N</tt> DFT.
 *
 * @param N The length of the DFT.
 *
 * @return the minimum wavenumber <tt>(1-N)/2</tt>.
 * @see The documentation for inorder.h for terminology details.
 */
inline
int suzerain_inorder_wavenumber_min(const int N)
{
    assert(0 < N);
    return (1-N)/2;
}

/**
 * Find the absolute value of the minimum wavenumber contained in a length
 * <tt>N</tt> DFT.  For even or odd <tt>N</tt>, this is the highest wavenumber
 * for which both a positive and negative mode is supported (e.g. 3 for N = 7
 * and N = 8).  For even <tt>N</tt>, this is one less than the Nyquist mode.
 *
 * @param N The length of the DFT.
 *
 * @return the absolute of the minimum wavenumber <tt>(N-1)/2</tt>.
 * @see The documentation for inorder.h for terminology details.
 */
inline
int suzerain_inorder_wavenumber_absmin(const int N)
{
    assert(0 < N);
    return (N-1)/2;
}

/**
 * Find the maximum wavenumber contained in a length <tt>N</tt> DFT.
 * For even <tt>N</tt> this is the Nyquist mode.
 *
 * @param N The length of the DFT.
 *
 * @return the maximum wavenumber <tt>N/2</tt>.
 * @see The documentation for inorder.h for terminology details.
 */
inline
int suzerain_inorder_wavenumber_max(const int N)
{
    assert(0 < N);
    return N/2;
}

/**
 * For a \e logically even length <tt>N</tt> DFT of a real-valued function, is
 * the wavenumber <tt>w</tt> the Nyquist mode?
 *
 * @param N The \e logical length of the DFT.
 * @param w The wavenumber of interest.
 *
 * @return One if the wavenumber is the Nyquist mode for <tt>N</tt> even.
 *         Zero otherwise.
 * @see The documentation for inorder.h for terminology details.
 */
inline
int suzerain_inorder_wavenumber_nyquist(const int N, const int w)
{
    // Deliberately no assertions on out-of-range wavenumbers!
    return (w == N/2)*(1 - (N % 2));
}

/**
 * For a \e logically even length <tt>N</tt> DFT of a real-valued function, is
 * the wavenumber associated with index <tt>i</tt> the Nyquist mode?
 *
 * For example, for <tt>N=6</tt> the values <tt>{ 0, 0, 0, 1, 0, 0 }</tt> will
 * be returned for <tt>i = 0, 1, 2, 3, 4, 5</tt>.  For <tt>N=5</tt> the values
 * <tt>{ 0, 0, 0, 0, 0 }</tt> will be returned for <tt>i = 0, 1, 2, 3, 4</tt>.
 *
 * @param N The \e logical length of the DFT.
 * @param i The index of interest.
 *
 * @return One if the indexed wavenumber is the Nyquist mode for <tt>N</tt>
 *         even.  Zero otherwise.
 * @see The documentation for inorder.h for terminology details.
 */
inline
int suzerain_inorder_wavenumber_nyquist_index(const int N, const int i)
{
    assert(0 <= i && i < N);
    return (i == N/2)*(1 - (N % 2));
}

/**
 * For a \e logically length <tt>N</tt> DFT of a real-valued function, are the
 * imaginary components of the wavenumber <tt>w</tt> \e necessarily zero due to
 * conjugate symmetry?
 *
 * @param N The \e logical length of the DFT.
 * @param w The wavenumber of interest.
 *
 * @return One if the wavenumber must have zero imaginary components.
 *         Zero otherwise.
 * @see The documentation for inorder.h for terminology details.
 */
inline
int suzerain_inorder_wavenumber_imagzero(const int N, const int w)
{
    return (w == 0) + suzerain_inorder_wavenumber_nyquist(N, w);
}

/**
 * For a \e logically length <tt>N</tt> DFT of a real-valued function, are the
 * imaginary components of the wavenumber associated with index <tt>i</tt> \e
 * necessarily zero due to conjugate symmetry?
 *
 * For example, for <tt>N=6</tt> the values <tt>{ 1, 0, 0, 1, 0, 0 }</tt> will
 * be returned for <tt>i = 0, 1, 2, 3, 4, 5</tt>.  For <tt>N=5</tt> the values
 * <tt>{ 1, 0, 0, 0, 0 }</tt> will be returned for <tt>i = 0, 1, 2, 3, 4</tt>.
 *
 * @param N The \e logical length of the DFT.
 * @param i The index of interest.
 *
 * @return One if the indexed wavenumber must have zero imaginary components.
 *         Zero otherwise.
 * @see The documentation for inorder.h for terminology details.
 */
inline
int suzerain_inorder_wavenumber_imagzero_index(const int N, const int i)
{
    return (i == 0) + suzerain_inorder_wavenumber_nyquist_index(N, i);
}

/**
 * Determine if the provided wavenumber is valid for a length <tt>N</tt> DFT.
 *
 * @param N The length of the DFT.
 * @param w The wavenumber of interest.
 *
 * @return True if the wavenumber is contained in such a DFT.  False otherwise.
 * @see The documentation for inorder.h for terminology details.
 */
inline
int suzerain_inorder_wavenumber_valid(const int N, const int w)
{
    return    suzerain_inorder_wavenumber_min(N) <= w
           && w <= suzerain_inorder_wavenumber_max(N);
}

/**
 * Compute the real-valued scaling factor used to differentiate a wave-space
 * signal on a dealiased DFT of length <tt>dN</tt>.  This method returns zero
 * for modes which cannot be supported on a DFT of length <tt>N</tt> or which
 * are removed by differentiation.
 *
 * For example, for <tt>dN=9</tt> and <tt>N=6</tt> the values <tt> {0, 1, 2, 0,
 * 0, 0, 0, -2, -1}</tt> will be returned for <tt>i=0, 1, 2, 3, 4, 5, 6, 7,
 * 8</tt>.
 *
 * @param N  The length of the DFT used to compute what wavenumbers are valid.
 * @param dN The dealiased length of the DFT.
 * @param i  The index of interest where <tt>0 <= i && i < dN</tt>.
 *
 * @return The wavenumber-like scaling factor associated with index <tt>i</tt>.
 * @see The documentation for inorder.h for terminology details.
 */
inline
int suzerain_inorder_wavenumber_diff(const int N, const int dN, const int i)
{
    assert(N <= dN);
    assert(0 <= i && i < dN);
    if (i < (N+1)/2) {
        return i;
    } else if (i >= dN - (N-1)/2) {
        return -dN+i;
    } else {
        return 0;
    }
}

/**
 * Compute an indicator for which indices from a length <tt>dN</tt> DFT are
 * translatable (as wavenumbers) to a smaller length <tt>N</tt> DFT.
 *
 * For example, for <tt>dN=9</tt> and <tt>N=6</tt> the values <tt> {1, 1, 1, 0,
 * 0, 0, 0, 1, 1}</tt> will be returned for <tt>i=0, 1, 2, 3, 4, 5, 6, 7,
 * 8</tt>.
 *
 * @param N  The length of the smaller DFT.
 * @param dN The length of the larger DFT.
 * @param i  The index of interest where <tt>0 <= i < dN</tt>.
 *
 * @return True if the wavenumber corresponding to index <tt>i</tt> on a
 *         length <tt>dN</tt> DFT is also present on a length <tt>N</tt> DFT.
 * @see The documentation for inorder.h for terminology details.
 */
inline
int suzerain_inorder_wavenumber_translatable(
        const int N, const int dN, const int i)
{
    assert(N <= dN);
    return suzerain_inorder_wavenumber_valid(
            N, suzerain_inorder_wavenumber(dN, i));
}

/**
 * Translate the index range <tt>[tb, te)</tt> from a length <tt>T</tt> DFT to
 * a length <tt>S</tt> DFT.  In-order storage for <tt>T != S</tt> implies the
 * index range may be broken into at most two contiguous ranges.
 *
 * More specifically, when ignoring wavenumbers not found in a length
 * <tt>S</tt> DFT, the data kept in <tt>[tb,te)</tt> is equivalent to
 * <tt>[tb1,te1)</tt> and <tt>[tb2,te2)</tt>.  Further, the range
 * <tt>[tb1,te1)</tt> contains data equivalent to <tt>[sb1,se1)</tt> and
 * <tt>[tb2,te2)</tt> contains data equivalent to <tt>[sb2,se2)</tt>.
 *
 * Such information could be used to copy data prepared for DFTs of disparate
 * lengths.  Two steps would be necessary.  First, a copy from
 * <tt>[sb1,se1)</tt> to <tt>[tb1,te1)</tt>.  Second, a copy from
 * <tt>[sb2,se2)</tt> to <tt>[tb2,te2)</tt>.  Either copy may be of zero
 * length.
 *
 * @param S    Length of source DFT.
 * @param T    Length of target DFT.
 * @param tb   Beginning index into target data.
 * @param te   Ending index into target data.
 * @param sb1  Beginning index of first range into source data.
 * @param se1  Ending index of first range into source data.
 * @param sb2  Beginning index of second range into source data.
 * @param se2  Ending index of second range into source data.
 * @param tb1  Target index corresponding to \c sb1
 * @param te1  Target index corresponding to \c se1
 * @param tb2  Target index corresponding to \c sb2
 * @param te2  Target index corresponding to \c se2
 */
void suzerain_inorder_wavenumber_translate(
        const int S, const int T, const int tb, const int te,
        int* sb1, int* se1, int* sb2, int* se2,
        int* tb1, int* te1, int* tb2, int* te2);

/**
 * For a length <tt>N</tt> DFT, find the index associated with wavenumber
 * <tt>w</tt>.
 *
 * For example, for <tt>N=6</tt> the values <tt>i = 0, 1, 2, 3, 4, 5</tt>
 * will be returned for <tt>{ 0, 1, 2, 3, -2, -1 }</tt>
 *
 * @param N The length of the DFT.
 * @param w The wavenumber of interest.
 *
 * @return The index associated with wavenumber <tt>w</tt>.
 * @see The documentation for inorder.h for terminology details.
 */
inline
int suzerain_inorder_index(const int N, const int w)
{
    assert(suzerain_inorder_wavenumber_valid(N, w));
    return (w >= 0) ? w : N + w;
}


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // SUZERAIN_INORDER_H
