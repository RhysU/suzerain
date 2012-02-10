/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2010, 2011, 2012 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *--------------------------------------------------------------------------
 * gbmatrix.h: General band storage matrix helpers
 * $Id$
 */

#ifndef __SUZERAIN_GBMATRIX_H
#define __SUZERAIN_GBMATRIX_H

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
 * @param ld Leading dimension for column-major ordering
 * @param kl Number of superdiagonals
 * @param ku Number of subdiagonals
 * @param i  Row index desired
 * @param j  Column index desired
 *
 * @return true if entry <tt>(i,j)</tt> falls within the matrix's band.
 */
static inline
int suzerain_gbmatrix_in_band(int ld, int kl, int ku, int i, int j) {
    /* Unused parameters present to be consistent with gbmatrix_offset */
    (void) ld;
    return j-ku <= i && i <= j+kl;
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
 * @param[in]  ld Leading dimension for column-major ordering.
 * @param[in]  s   Size of each matrix entry in bytes.
 *                 For example, <tt>sizeof(double)</tt>.
 * @param[in]  i   Desired row in the half-open range <tt>[0,m)</tt>.
 * @param[out] jl  First nonzero entry's column index \c j in row \c i.
 * @param[out] ju  Last plus one nonzero entry's column index \c in row \c i.
 * @param[out] inc Increment between entries in row \c i.
 *
 * @return A pointer which, after been cast to match \c s, may be
 *         dereferenced to access entries in rho \c i of \c a.
 */
static inline
void *suzerain_gbmatrix_row(int m, int n, int kl, int ku,
                            void *a, int ld, size_t s, int i,
                            int *jl, int *ju, int *inc)
{
    // Column-major banded matrix layout: a[(ku + i)*inc + j*(ld - inc)]

    // Transpose the matrix storage information to traverse by rows
    *inc = ld - 1;                          // Start from column-major...
    a = ((char*)a) + s*(ku - kl*(*inc));    // ...and traverse a by rows
    kl ^= ku; ku ^= kl; kl ^= ku;           // Swap kl and ku for A^T
    m  ^= n;  n  ^= m;  m  ^= n;            // Swap m and n for A^T

    // Now our row problem looks just like indexing the i-th column of A
    *jl = i - ku;     if (*jl < 0) *jl = 0;
    *ju = i + kl + 1; if (*ju > m) *ju = m;
    return ((char*)a) + s*(ku*(*inc) + i*(ld - *inc));
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // __SUZERAIN_GBMATRIX_H
