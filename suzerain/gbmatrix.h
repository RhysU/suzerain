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

#ifndef SUZERAIN_GBMATRIX_H
#define SUZERAIN_GBMATRIX_H

/** @file
 * General band storage matrix helpers
 */

#include <assert.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Compute the BLAS-compatible offset to <tt>a(i,j)</tt> for general
 * banded matrices <tt>a(i,j) -> storage(ku+i-j,j)</tt> where storage
 * is column-major with leading dimension ld.  When comparing with BLAS
 * documentation, note missing constant one compared to Fortran due to
 * C's 0-indexing.
 *
 * @param ld Leading dimension for column-major ordering
 * @param kl Number of superdiagonals
 * @param ku Number of subdiagonals
 * @param i  Row index desired
 * @param j  Column index desired
 *
 * @return Column-major offset where entry <tt>(i,j)</tt> is stored.
 */
static inline
int suzerain_gbmatrix_offset(int ld, int kl, int ku, int i, int j) {
    /* Unused parameters present to be consistent with gbmatrix_in_band */
    (void) kl;
    return j*ld+(ku+i-j);
}

/**
 * Determine if indices fall within the band of a general banded matrix.
 *
 * @param kl Number of superdiagonals
 * @param ku Number of subdiagonals
 * @param i  Row index desired
 * @param j  Column index desired
 *
 * @return true if entry <tt>(i,j)</tt> falls within the matrix's band.
 */
static inline
int suzerain_gbmatrix_in_band(int kl, int ku, int i, int j) {
    return j-ku <= i && i <= j+kl;
}

/**
 * Compute the starting point and extents of the <tt>j</tt>-th column within a
 * general banded matrix.
 *
 * Assuming <tt>s == sizeof(double)</tt> on entry and the returned value was \c
 * p, then <tt>((double *)p)[il]</tt> is the first nonzero entry in column \c j
 * and <tt>((double *)p)[iu]</tt> is the one past the last nonzero entry of the
 * same column.
 *
 * @param[in]  m  Number of rows.
 * @param[in]  n  Number of columns.
 * @param[in]  kl Number of subdiagonals.
 * @param[in]  ku Number of superdiagonals.
 * @param[in]  a  BLAS-compatible general banded storage.
 * @param[in]  ld Leading dimension for column-major ordering.
 * @param[in]  s  Size of each matrix entry in bytes.
 *                For example, <tt>sizeof(double)</tt>.
 * @param[in]  j  Desired column in the half-open range <tt>[0,n)</tt>.
 * @param[out] il First nonzero entry's row index \c i in column \c j.
 * @param[out] iu Last plus one nonzero entry's row index \c i in column \c j.
 *
 * @return A pointer which, after being cast to match \c s, may be
 *         dereferenced to access entries in column \c j of \c a.
 */
static inline
void *suzerain_gbmatrix_col(int m, int n, int kl, int ku,
                            void *a, int ld, size_t s, int j,
                            int *il, int *iu)
{
    (void) n;
    assert(m  >= 0);
    assert(n  >= 0);
    assert(kl >= 0);
    assert(ku >= 0);
    assert(ld >= kl + ku);

    *il = j - ku;     if (*il < 0) *il = 0;
    *iu = j + kl + 1; if (*iu > m) *iu = m;
    return ((char*)a) + s*(ku + j*(ld - 1));
}

/**
 * Compute the starting point, stride, and extents of the <tt>i</tt>-th
 * row within a general banded matrix.
 *
 * Assuming <tt>s == sizeof(double)</tt> on entry and the returned value was \c
 * p, then <tt>((double *)p)[jl*inc]</tt> is the first nonzero entry in row \c
 * i and <tt>((double *)p)[ju*inc]</tt> is the one past the last nonzero entry
 * of the same row.
 *
 * @param[in]  m   Number of rows.
 * @param[in]  n   Number of columns.
 * @param[in]  kl  Number of subdiagonals.
 * @param[in]  ku  Number of superdiagonals.
 * @param[in]  a   BLAS-compatible general banded storage.
 * @param[in]  ld  Leading dimension for column-major ordering.
 * @param[in]  s   Size of each matrix entry in bytes.
 *                 For example, <tt>sizeof(double)</tt>.
 * @param[in]  i   Desired row in the half-open range <tt>[0,m)</tt>.
 * @param[out] jl  First nonzero entry's column index \c j in row \c i.
 * @param[out] ju  Last plus one nonzero entry's column index \c j in row \c i.
 * @param[out] inc Increment between entries in row \c i.
 *
 * @return A pointer which, after being cast to match \c s, may be
 *         dereferenced to access entries in row \c i of \c a.
 */
static inline
void *suzerain_gbmatrix_row(int m, int n, int kl, int ku,
                            void *a, int ld, size_t s, int i,
                            int *jl, int *ju, int *inc)
{
    (void) m;
    assert(m  >= 0);
    assert(n  >= 0);
    assert(kl >= 0);
    assert(ku >= 0);
    assert(ld >= kl + ku);

    // Column-major banded matrix layout: a[(ku + i)*inc + j*(ld - inc)]
    //
    // Logically, here's what needs to happen...
    //
    // Transpose the matrix storage information to traverse by rows
    // *inc = ld - 1;                          // Start from column-major...
    // a = ((char*)a) + s*(ku - kl*(*inc));    // ...and traverse a by rows
    // kl ^= ku; ku ^= kl; kl ^= ku;           // Swap kl and ku for A^T
    // m  ^= n;  n  ^= m;  m  ^= n;            // Swap m and n for A^T
    //
    // Now our row problem looks just like indexing the i-th column of A
    // *jl = i - ku;     if (*jl < 0) *jl = 0;
    // *ju = i + kl + 1; if (*ju > m) *ju = m;
    // return ((char*)a) + s*(ku*(*inc) + i*(ld - *inc));
    //
    // In practice, here's the shortest path there...

    *inc = ld - 1;
    *jl  = i - kl;     if (*jl < 0) *jl = 0;
    *ju  = i + ku + 1; if (*ju > n) *ju = n;
    return ((char*)a) + s*(ku + i);
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // SUZERAIN_GBMATRIX_H
