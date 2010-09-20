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
 * gbmatrix.h: General band storage matrix helpers
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_GBMATRIX_H
#define __SUZERAIN_GBMATRIX_H

#include <suzerain/common.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Compute the BLAS-compatible offset to <tt>a(i,j)</tt> for general
 * banded matrices <tt>a(i,j) -> storage(ku+i-j,j)</tt> where storage
 * is column-major with LDA lda.  When comparing with BLAS documentation, note
 * missing constant one compared to Fortran due to C's 0-indexing.
 *
 * @param lda Leading dimension for column-major ordering
 * @param kl  Number of superdiagonals
 * @param ku  Number of subdiagonals
 * @param i   Row index desired
 * @param j   Column index desired
 *
 * @return Column-major offset where entry <tt>(i,j)</tt> is stored.
 */
inline
int suzerain_gbmatrix_offset(int lda, int kl, int ku, int i, int j) {
    /* Unused parameters present to be consistent with gbmatrix_in_band */
    SUZERAIN_UNUSED(kl);
    return j*lda+(ku+i-j);
}

/** Determine if indices fall within the band of a general banded matrix.
 *
 * @param lda Leading dimension for column-major ordering
 * @param kl  Number of superdiagonals
 * @param ku  Number of subdiagonals
 * @param i   Row index desired
 * @param j   Column index desired
 *
 * @return true if entry <tt>(i,j)</tt> falls within the matrix's band.
 */
inline
int suzerain_gbmatrix_in_band(int lda, int kl, int ku, int i, int j) {
    /* Unused parameters present to be consistent with gbmatrix_offset */
    SUZERAIN_UNUSED(lda);
    return j-ku <= i && i <= j+kl;
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // __SUZERAIN_GBMATRIX_H
