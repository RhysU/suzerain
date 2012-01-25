/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2012 The PECOS Development Team
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
 * bsmbsm.h: routines for blocked square matrices with banded submatrices
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_BSMBSM_H__
#define __SUZERAIN_BSMBSM_H__

/** @file
 * Utilities for working with blocked square matrices with banded submatrices
 * (BSMBSM).
 *
 * Matrix \f$A\f$ is a blocked square matrix with banded submatrices (BSMBSM)
 * when
 * \f[A = \begin{pmatrix}
 *     B^{0\,0}  & \cdots & B^{0\,S-1}   \\
 *     \vdots    & \ddots & \vdots       \\
 *     B^{S-1,0} & \cdots & B^{S-1\,S-1}
 * \end{pmatrix}\f]
 * where every \f$B^{i\,j}\f$ is an <tt>n</tt>x<tt>n</tt> banded submatrix
 * containing \c kl subdiagonals and \c ku superdiagonals.  We set the
 * convention that lowercase identifiers indicate submatrix details while
 * uppercase identifiers indicate global matrix details.  The structure of a
 * BSMBSM is defined completely by the parameters <tt>S</tt>, <tt>n</tt>,
 * <tt>kl</tt>, and <tt>ku</tt>.  The number of rows and columns is <tt>N :=
 * S*n</tt>.
 *
 * Building \f$A\f$ from individually contiguous, banded submatrices
 * \f$B_{i\,j}\f$ within \f$A\f$ is both convenient and efficient.  For
 * example, banded matrix accumulation operations and boundary condition
 * imposition are simple in such a storage format.  However, using individually
 * contiguous, banded submatrices is grossly inefficient for solving linear
 * equations.
 *
 * With appropriate renumbering of \f$A\f$, solving linear equations can be
 * done efficiently.  The zero-indexed permutation vector \f[q(i) \coloneqq
 * \left(i\bmod{}S\right)n + \lfloor{}i/S\rfloor{}\f] may always be used to
 * convert a BSMBSM into a globally banded <tt>N</tt>x<tt>N</tt> matrix with
 * minimum bandwidth.  More concretely, the permutation matrix \f$P\f$ uniquely
 * defined by vector \f$q\f$ causes \f$P A P^{\mbox{T}}\f$ to have <tt>KL :=
 * S*(kl+1)-1</tt> subdiagonals and <tt>KU := S*(ku+1)-1</tt> superdiagonals
 * summing to overall bandwidth <tt>KL + 1 + KU = S*(kl + ku + 2)-1</tt>.  The
 * reverse permutation vector has a simple closed form \f[q^{-1}(i) \coloneqq
 * \left(i\bmod{}n\right)S + \lfloor{}i/n\rfloor{}.\f] With \f$A_{i\,j}\f$ in
 * hand, the banded renumbering may formed using the relationships
 * \f{align*}{
 *        \left.A\right|_{i\,j}
 *     &= \left.P A P^{\mbox{T}}\right|_{q(i)\,q(j)}
 *     &
 *        \left.P A P^{\mbox{T}}\right|_{i\,j}
 *     &= \left.A\right|_{q^{-1}(i)\,q^{-1}(j)}.
 * \f}
 * This renumbering can be LU factorized in order <tt>N*(KL + 1 + KU)^2 =
 * S*n*(S*(kl + ku + 2)-1)^2</tt> floating point operations to find \f$LU = P A
 * P^{\mbox{T}}\f$.  The linear equation \f$AX=B\f$, which is equivalent to
 * \f$LUPX=PB\f$, has the solution \f[X = A^{-1}B =
 * P^{\mbox{T}}\left(LU\right)^{-1}PB\f] where inversion has been used as a
 * notational convenience representing triangular backsubstitution.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Storage details for a banded square matrix with banded submatrices (BSMBSM).
 * No data is stored within this POD type.  However, the struct contains
 * redundant data which should be computed consistently using
 * suzerain_bsmbsm_construct().
 *
 * @see bsmbsm.h for details on BSMBSMs.
 */
typedef struct suzerain_bsmbsm {

    int S;  /**< Number of rows and columns of banded submatrices */
    int n;  /**< Number or rows and columns in each banded submatrix */
    int kl; /**< Number of subdiagonals in each banded submatrix */
    int ku; /**< Number of superdiagonals in each banded submatrix */
    int ld; /**< Minimum leading dimension for each submatrix */
    int N;  /**< Number of rows and colums in global matrix */
    int KL; /**< Number of subdiagonals in matrix after renumbering */
    int KU; /**< Number of superdiagonals in matrix after renumbering */
    int LD; /**< Minimum leading dimension for matrix after renumbering */

} suzerain_bsmbsm;

/**
 * Populate a suzerain_bsmbsm instance with all BSMBSM storage details.
 *
 * @param S   Number of rows and columns of banded submatrices
 * @param n   Number or rows and columns in each banded submatrix
 * @param kl  Number of subdiagonals in each banded submatrix
 * @param ku  Number of superdiagonals in each banded submatrix
 *
 * @return a fully populated BSMBSM storage description.
 *
 * @see bsmbsm.h for details on BSMBSMs.
 */
static inline
suzerain_bsmbsm suzerain_bsmbsm_construct(int S, int n, int kl, int ku)
{
    suzerain_bsmbsm tmp;
    tmp.S  = S;
    tmp.n  = n;
    tmp.kl = kl;
    tmp.ku = ku;
    tmp.ld = tmp.kl + 1 + tmp.ku;
    tmp.N  = tmp.S*tmp.n;
    tmp.KL = tmp.S*(tmp.kl + 1) - 1;
    tmp.KU = tmp.S*(tmp.ku + 1) - 1;
    tmp.LD = tmp.KL + 1 + tmp.KU;
    return tmp;
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __SUZERAIN_BSMBSM_H__ */
