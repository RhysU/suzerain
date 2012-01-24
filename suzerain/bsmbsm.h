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
 * (BSMBSM).  In particular, data structures to aid building and efficiently
 * solving linear equations with BSMBSM operators are provided.
 *
 * Matrix \f$A\f$ is a BSMBSM when
 * \f[A = \begin{pmatrix}
 *     B_{0\,0}  & \cdots & B_{0\,S-1}   \\
 *     \vdots    & \ddots & \vdots       \\
 *     B_{S-1,0} & \cdots & B_{S-1\,S-1}
 * \end{pmatrix}\f]
 * where every \f$B_{i\,j}\f$ is an <tt>n</tt>x<tt>n</tt> banded submatrix
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
 * The zero-indexed permutation vector \f$q(i) \coloneqq
 * \left(i\bmod{}S\right)n + \lfloor{}i/S\rfloor{}\f$ may always be used to
 * convert a BSMBSM into a globally banded <tt>N</tt>x<tt>N</tt> matrix with
 * minimum bandwidth.  More concretely, the permutation matrix \f$P\f$ uniquely
 * defined by vector \f$q\f$ causes \f$P A P^{\mbox{T}}\f$ to have <tt>KL :=
 * S*(kl+1)-1</tt> subdiagonals and <tt>KU := S*(ku+1)-1</tt> superdiagonals
 * summing to overall bandwidth <tt>KL + 1 + KU = S*(kl + ku + 2)-1</tt>.  Note
 * that the reverse permutation vector has a simple closed form \f$q^{-1}(i)
 * \coloneqq \left(i\bmod{}n\right)S + \lfloor{}i/n\rfloor{}\f$.
 *
 * With \f$A\f$ in hand, the banded renumbering \f$P A P^{\mbox{T}}\f$ may be
 * formed using \f$\left.P A P^{\mbox{T}}\right|_{i\,j} =
 * \left.A\right|_{q^{-1}(i)\,q^{-1}(j)}\f$.  This renumbering may be LU
 * factorized in order <tt>S*n*(S*(kl + ku + 2)-1)^2</tt> floating point
 * operations to find \f$LU = P A P^{\mbox{T}}\f$.  The linear equation
 * \f$AX=B\f$, which is equivalent to \f$LUPX=PB\f$, has the solution \f$X =
 * A^{-1}B = P^{\mbox{T}}\left(LU\right)^{-1}PB\f$.  Here inversion has been
 * used as a notational convenience for solve operations.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Encapsulates storage of a complex-valued BSMBSM matrix.
 * @see bsmbsm.h for details on BSMBSMs.
 */
typedef struct suzerain_zbsmbsm {

    int n;

    int kl;

    int ku;

    int ld;

    int S;

    int N;

    int KL;

    int KU;

    int LD;

} suzerain_zbsmbsm;

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __SUZERAIN_BSMBSM_H__ */
