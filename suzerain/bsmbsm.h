/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2012-2014 Rhys Ulerich
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

#ifndef SUZERAIN_BSMBSM_H
#define SUZERAIN_BSMBSM_H

/** @file
 * Routines for blocked square matrices with banded submatrices.
 */

#include <assert.h>

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute.h>

#include <suzerain/complex.h>

/** @file
 * Utilities for working with blocked square matrices with banded submatrices
 * (BSMBSM).
 *
 * Matrix \f$A\f$ is a blocked square matrix with banded submatrices (BSMBSM)
 * when
 * \f[A = \begin{pmatrix}
 *     B^{0,0}  & \cdots & B^{0,S-1}   \\
 *     \vdots    & \ddots & \vdots       \\
 *     B^{S-1,0} & \cdots & B^{S-1,S-1}
 * \end{pmatrix}\f]
 * where every \f$B^{i,j}\f$ is an <tt>n</tt>x<tt>n</tt> banded submatrix
 * containing \c kl subdiagonals and \c ku superdiagonals.  We set the
 * convention that lowercase identifiers indicate submatrix details while
 * uppercase identifiers indicate global matrix details.  The structure of a
 * BSMBSM is defined completely by the parameters <tt>S</tt>, <tt>n</tt>,
 * <tt>kl</tt>, and <tt>ku</tt>.  The number of rows and columns is <tt>N :=
 * S*n</tt>.
 *
 * Building \f$A\f$ from individually contiguous, banded submatrices
 * \f$B_{i,j}\f$ within \f$A\f$ is both convenient and efficient.  For
 * example, banded matrix accumulation operations and boundary condition
 * imposition are simple in such a storage format.  However, using individually
 * contiguous, banded submatrices is grossly inefficient for solving linear
 * equations.
 *
 * With appropriate renumbering of \f$A\f$, solving linear equations can be
 * done efficiently.  The zero-indexed permutation vector \f[q(i) =
 * \left(i\bmod{}S\right)n + \lfloor{}i/S\rfloor{}\f] may always be used to
 * convert a BSMBSM into a globally banded <tt>N</tt>x<tt>N</tt> matrix with
 * minimum bandwidth.  More concretely, the permutation matrix \f$P\f$ uniquely
 * defined by vector \f$q\f$ causes \f$P A P^{\mbox{T}}\f$ to have <tt>KL :=
 * S*(kl+1)-1</tt> subdiagonals and <tt>KU := S*(ku+1)-1</tt> superdiagonals
 * summing to overall bandwidth <tt>KL + 1 + KU = S*(kl + ku + 2)-1</tt>.  The
 * reverse permutation vector has a simple closed form \f[q^{-1}(i) =
 * \left(i\bmod{}n\right)S + \lfloor{}i/n\rfloor{}.\f] With \f$A_{i,j}\f$ in
 * hand, the banded renumbering may formed using the relationships
 * \f{align*}{
 *        \left.A\right|_{i,j}
 *     &= \left.P A P^{\mbox{T}}\right|_{q^{-1}(i),q^{-1}(j)}
 *     &
 *        \left.P A P^{\mbox{T}}\right|_{i,j}
 *     &= \left.A\right|_{q(i),q(j)}.
 * \f}
 * This renumbering can be LU factorized in order <tt>N*(KL + 1 + KU)^2 =
 * S*n*(S*(kl + ku + 2)-1)^2</tt> floating point operations to find \f$LU = P A
 * P^{\mbox{T}}\f$.  The linear equation \f$AX=B\f$, which is equivalent to
 * \f$LUPX=PB\f$, has the solution \f[X = A^{-1}B =
 * P^{\mbox{T}}\left(LU\right)^{-1}PB\f] where inversion has been used as a
 * notational convenience representing triangular back substitution.
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
 * @see bsmbsm.h for full details on BSMBSMs.
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
 * @param[in] S  Number of rows and columns of banded submatrices
 * @param[in] n  Number or rows and columns in each banded submatrix
 * @param[in] kl Number of subdiagonals in each banded submatrix
 * @param[in] ku Number of superdiagonals in each banded submatrix
 *
 * @return a fully populated BSMBSM storage description.
 *
 * @see bsmbsm.h for full details on BSMBSMs.
 */
static inline
suzerain_bsmbsm suzerain_bsmbsm_construct(int S, int n, int kl, int ku)
{
    assert(S  >  0);
    assert(n  >  0);
    assert(kl >= 0);
    assert(ku >= 0);

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

/**
 * Compute the <tt>i</tt>th entry in the permutation vector \f$q\f$.
 *
 * @param[in] S Number of rows and columns of banded submatrices
 * @param[in] n Number or rows and columns in each banded submatrix
 * @param[in] i Index \f$i\f$ in the range <tt>[0, S*n)</tt>.
 *
 * @return the value of \f$q(i)\f$.
 *
 * @see bsmbsm.h for full details on the permutation vector \f$q\f$.
 */
static inline
int suzerain_bsmbsm_q(int S, int n, int i)
{
    assert(0 <= i && i < S*n);

    div_t t = div(i, S);
    return t.rem*n + t.quot;
}

/**
 * Compute the <tt>i</tt>th entry in the permutation vector \f$q^{-1}\f$.
 *
 * @param[in] S Number of rows and columns of banded submatrices
 * @param[in] n Number or rows and columns in each banded submatrix
 * @param[in] i Index \f$i\f$ in the range <tt>[0, S*n)</tt>.
 *
 * @return the value of \f$q^{-1}(i)\f$.
 *
 * @see bsmbsm.h for full details on the permutation vector \f$q^{-1}\f$.
 */
static inline
int suzerain_bsmbsm_qinv(int S, int n, int i)
{
    return suzerain_bsmbsm_q(n, S, i);
}

/**
 * @brief Compute \f$ y \leftarrow{} \alpha{} P x + \beta{}y \f$ or \f$ y
 * \leftarrow{} \alpha{} P^{\mbox{T}} x + \beta{}y \f$ where \f$P\f$ is defined
 * by permutation vector \f$q\f$.  Storage \c y must not alias storage \c x.
 * Negative strides may be used and are interpreted as in the BLAS.
 *
 * @param[in]  trans Either 'N' for \f$P x\f$ or 'T' for \f$P^{\mbox{T}}\f$.
 * @param[in]  S Number of rows and columns of banded submatrices
 * @param[in]  n Number or rows and columns in each banded submatrix
 * @param[in]  alpha Multiplicative scalar \f$ \alpha \f$
 * @param[in]  x First source vector of length <tt>S*n</tt>.
 * @param[in]  incx First source vector stride.
 * @param[in]  beta Multiplicative scalar \f$ \beta \f$
 * @param[out] y Second source vector and target vector of length <tt>S*n</tt>.
 * @param[in]  incy Second source vector and target vector stride.
 *
 * @return Zero on success.
 *         Otherwise calls suzerain_blas_xerbla() and returns nonzero.
 * @see bsmbsm.h for full details on the permutation vector \f$q\f$.
 */
int
suzerain_bsmbsm_saPxpby(
    char trans,
    int S,
    int n,
    const float alpha,
    const float *x,
    int incx,
    const float beta,
    float *y,
    int incy);

/** @copydoc suzerain_bsmbsm_saPxpby */
int
suzerain_bsmbsm_daPxpby(
    char trans,
    int S,
    int n,
    const double alpha,
    const double *x,
    int incx,
    const double beta,
    double *y,
    int incy);

/** @copydoc suzerain_bsmbsm_saPxpby */
int
suzerain_bsmbsm_caPxpby(
    char trans,
    int S,
    int n,
    const complex_float alpha,
    const complex_float *x,
    int incx,
    const complex_float beta,
    complex_float *y,
    int incy);

/** @copydoc suzerain_bsmbsm_saPxpby */
int
suzerain_bsmbsm_zaPxpby(
    char trans,
    int S,
    int n,
    const complex_double alpha,
    const complex_double *x,
    int incx,
    const complex_double beta,
    complex_double *y,
    int incy);

/**
 * Create a <a
 * href="http://www.gnu.org/software/gsl/">GNU Scientific Library</a> (GSL) <a
 * href="http://www.gnu.org/software/gsl/manual/html_node/Permutations.html">
 * gsl_permutation</a> representing permutation vector \f$q\f$.
 * The returned <tt>gsl_permutation*</tt> must be deallocated using
 * <tt>gsl_permutation_free</tt>.
 *
 * After creation, the GSL methods
 * <ul>
 *     <li><code>gsl_permute_float</code>,
 *         <code>gsl_permute_inverse</code></li>
 *     <li><code>gsl_permute</code>,
 *         <code>gsl_permute_float_inverse</code></li>
 *     <li><code>gsl_permute_complex_float</code>,
 *         <code>gsl_permute_complex_float_inverse</code></li>
 *     <li><code>gsl_permute_complex</code>,
 *         <code>gsl_permute_complex_inverse</code></li>
 * </ul>
 * may be used to apply \f$q\f$ or \f$q^{-1}\f$ <em>in-place</em> on
 * <tt>float</tt>, <tt>double</tt>, <tt>float[2]</tt>, or <tt>double[2]</tt>
 * data (with the complex cases requiring casting).
 *
 * Wherever possible other methods like suzerain_bsmbsm_dPcopy(),
 * suzerain_bsmbsm_dPcopy(), suzerain_bsmbsm_daPxpby(), or
 * suzerain_bsmbsm_daPxpby() should be used to perform out-of-place
 * permutation is is much, much faster than the in-place methods available
 * in the GSL.
 *
 * @param[in] S Number of rows and columns of banded submatrices
 * @param[in] n Number or rows and columns in each banded submatrix
 *
 * @return A <tt>gsl_permutation*</tt> representing \f$q\f$ on success.
 *         NULL otherwise (e.g. on memory allocation failure).
 */
gsl_permutation * suzerain_bsmbsm_permutation(int S, int n);

/**
 * Pack scaled, banded submatrix \f$\alpha{}B^{\hat{\imath},\hat{\jmath}}\f$
 * into the corresponding locations within banded storage of \f$P A
 * P^{\mbox{T}}\f$.  More specifically,
 * \f[
 *   \left.\alpha{}B^{\hat{\imath},\hat{\jmath}}\right|_{i,j}
 *   =
 *   \left.A\right|_{\hat{\imath}n + i, \hat{\jmath}n + j}
 *   =
 *   \left.PAP^{\mbox{T}}\right|_{
 *       q^{-1}\left(\hat{\imath}n + i\right)
 *       ,
 *       q^{-1}\left(\hat{\jmath}n + j\right)
 *   }
 * \f]
 * relates the source and target values through the permutation vector \f$q\f$.
 * Storage \c b must not alias storage \c papt.  All calls packing data into
 * one \f$P A P^{\mbox{T}}\f$ must supply the same parameters for \c KL and \c
 * KU.
 *
 * @param[in]  S      Number of rows and columns of banded submatrices.
 * @param[in]  n      Number or rows and columns in each banded submatrix.
 * @param[in]  ihat   Submatrix row offset \f$\hat{\imath}\f$.
 * @param[in]  jhat   Submatrix column offset \f$\hat{\jmath}\f$.
 * @param[in]  kl     Number of subdiagonals in each banded submatrix.
 * @param[in]  ku     Number of superdiagonals in each banded submatrix.
 * @param[in]  alpha  Scaling factor \f$\alpha\f$ to apply to submatrix.
 * @param[in]  b      Band storage of submatrix
 *                    \f$B^{\hat{\imath},\hat{\jmath}}\f$.
 * @param[in]  ldb    Leading dimension of storage \c b.
 * @param[in]  KL     Number of subdiagonals in the renumbered matrix
 *                    which must be at least <tt>S*(kl + 1) - 1</tt>.
 * @param[in]  KU     Number of superdiagonals in the renumbered matrix
 *                    which must be at least <tt>S*(ku + 1) - 1</tt>.
 * @param[out] papt   Band storage of renumbered matrix \f$PAP^{\mbox{T}}\f$
 * @param[in]  ldpapt Leading dimension of storage \c papt.
 *                    which must be <tt>KU + 1 + KL</tt>.
 *
 * @return Zero on success.
 *         Otherwise calls suzerain_blas_xerbla() and returns nonzero.
 * @see bsmbsm.h for full details on the permutation vector \f$q\f$.
 */
int
suzerain_bsmbsm_spack(
        int S,  int n, int ihat, int jhat,
        int kl, int ku, const float alpha, const float *b,    int ldb,
        int KL, int KU,                          float *papt, int ldpapt);

/** @copydoc suzerain_bsmbsm_spack */
int
suzerain_bsmbsm_dpack(
        int S,  int n, int ihat, int jhat,
        int kl, int ku, const double alpha, const double *b,    int ldb,
        int KL, int KU,                           double *papt, int ldpapt);

/** @copydoc suzerain_bsmbsm_spack */
int
suzerain_bsmbsm_cpack(
        int S,  int n, int ihat, int jhat,
        int kl, int ku, const complex_float alpha,
                        const complex_float *b,    int ldb,
        int KL, int KU,       complex_float *papt, int ldpapt);

/** @copydoc suzerain_bsmbsm_spack */
int
suzerain_bsmbsm_zpack(
        int S,  int n, int ihat, int jhat,
        int kl, int ku, const complex_double alpha,
                        const complex_double *b,    int ldb,
        int KL, int KU,       complex_double *papt, int ldpapt);

/** @copydoc suzerain_bsmbsm_spack */
int
suzerain_bsmbsm_cspack(
        int S,  int n, int ihat, int jhat,
        int kl, int ku, const complex_float alpha,
                        const         float *b,    int ldb,
        int KL, int KU,       complex_float *papt, int ldpapt);

/** @copydoc suzerain_bsmbsm_spack */
int
suzerain_bsmbsm_zdpack(
        int S,  int n, int ihat, int jhat,
        int kl, int ku, const complex_double alpha,
                        const         double *b,    int ldb,
        int KL, int KU,       complex_double *papt, int ldpapt);

/**
 * Pack contiguous, scaled, banded submatrix
 * \f$\alpha{}B^{\hat{\imath},\hat{\jmath}}\f$ into the corresponding
 * locations within contiguous storage of \f$P A P^{\mbox{T}}\f$.  This is a
 * convenience method to simplify preparing storage for use by the BLAS'
 * <tt>gbmv</tt> or LAPACK's <tt>gbsvx</tt>.
 *
 * @param[in]  A     Storage details for the BSMBSM matrix.
 * @param[in]  ihat  Submatrix row offset \f$\hat{\imath}\f$.
 * @param[in]  jhat  Submatrix column offset \f$\hat{\jmath}\f$.
 * @param[in]  alpha Scaling factor \f$\alpha\f$ to apply to submatrix.
 * @param[in]  b     Band storage of submatrix
 *                   \f$B^{\hat{\imath},\hat{\jmath}}\f$
 *                   which <em>must</em> have <tt>A->kl</tt> and
 *                   <tt>A->ku</tt> diagonals and leading dimension
 *                   <tt>A->ld = A->ku + 1 + A->kl</tt>.
 * @param[out] papt  Band storage of renumbered matrix \f$PAP^{\mbox{T}}\f$
 *                   which <em>must</em> have <tt>A->KL</tt> and
 *                   <tt>A->KU</tt> diagonals and leading dimension
 *                   <tt>A->LD</tt>.
 *
 * @return Zero on success.
 *         Otherwise calls suzerain_blas_xerbla() and returns nonzero.
 */
inline static int
suzerain_bsmbsm_spackc(const suzerain_bsmbsm *A, int ihat, int jhat,
                       const float alpha, const float *b, float *papt)
{
    return suzerain_bsmbsm_spack(A->S, A->n, ihat, jhat,
                                 A->kl, A->ku, alpha, b,    A->ld,
                                 A->KL, A->KU,        papt, A->LD);
}

/** @copydoc suzerain_bsmbsm_spackc */
inline static int
suzerain_bsmbsm_dpackc(const suzerain_bsmbsm *A, int ihat, int jhat,
                       const double alpha, const double *b, double *papt)
{
    return suzerain_bsmbsm_dpack(A->S, A->n, ihat, jhat,
                                 A->kl, A->ku, alpha, b,    A->ld,
                                 A->KL, A->KU,        papt, A->LD);
}

/** @copydoc suzerain_bsmbsm_spackc */
inline static int
suzerain_bsmbsm_cpackc(const suzerain_bsmbsm *A, int ihat, int jhat,
                       const complex_float alpha, const complex_float *b,
                                                        complex_float *papt)
{
    return suzerain_bsmbsm_cpack(A->S, A->n, ihat, jhat,
                                 A->kl, A->ku, alpha, b,    A->ld,
                                 A->KL, A->KU,        papt, A->LD);
}

/** @copydoc suzerain_bsmbsm_spackc */
inline static int
suzerain_bsmbsm_zpackc(const suzerain_bsmbsm *A, int ihat, int jhat,
                       const complex_double alpha, const complex_double *b,
                                                         complex_double *papt)
{
    return suzerain_bsmbsm_zpack(A->S, A->n, ihat, jhat,
                                 A->kl, A->ku, alpha, b,    A->ld,
                                 A->KL, A->KU,        papt, A->LD);
}

/** @copydoc suzerain_bsmbsm_spackc */
inline static int
suzerain_bsmbsm_cspackc(const suzerain_bsmbsm *A, int ihat, int jhat,
                        const complex_float alpha, const float *b,
                        complex_float *papt)
{
    return suzerain_bsmbsm_cspack(A->S, A->n, ihat, jhat,
                                  A->kl, A->ku, alpha, b,    A->ld,
                                  A->KL, A->KU,        papt, A->LD);
}

/** @copydoc suzerain_bsmbsm_spackc */
inline static int
suzerain_bsmbsm_zdpackc(const suzerain_bsmbsm *A, int ihat, int jhat,
                        const complex_double alpha, const double *b,
                        complex_double *papt)
{
    return suzerain_bsmbsm_zdpack(A->S, A->n, ihat, jhat,
                                  A->kl, A->ku, alpha, b,    A->ld,
                                  A->KL, A->KU,        papt, A->LD);
}

/**
 * Pack contiguous, scaled, banded submatrix
 * \f$\alpha{}B^{\hat{\imath},\hat{\jmath}}\f$ into the corresponding
 * locations within contiguous, LU factorization-ready storage of \f$P A
 * P^{\mbox{T}}\f$.  This is a convenience method to simplify preparing storage
 * for use by LAPACK's <tt>gbtrf</tt> or <tt>gbsv</tt>.
 *
 * @param[in]  A     Storage details for the BSMBSM matrix.
 * @param[in]  ihat  Submatrix row offset \f$\hat{\imath}\f$.
 * @param[in]  jhat  Submatrix column offset \f$\hat{\jmath}\f$.
 * @param[in]  alpha Scaling factor \f$\alpha\f$ to apply to submatrix.
 * @param[in]  b     Band storage of submatrix
 *                   \f$B^{\hat{\imath},\hat{\jmath}}\f$
 *                   which <em>must</em> have <tt>A->kl</tt> and
 *                   <tt>A->ku</tt> diagonals and leading dimension
 *                   <tt>A->ld = A->ku + 1 + A->kl</tt>.
 * @param[out] papt  Band storage of renumbered matrix \f$PAP^{\mbox{T}}\f$
 *                   which <em>must</em> have <tt>A->KL</tt> and <tt>A->KU</tt>
 *                   diagonals and leading dimension
 *                   <tt>A->LD + A->KL = 2*A->KL + A->KU + 1</tt>.
 *
 * @return Zero on success.
 *         Otherwise calls suzerain_blas_xerbla() and returns nonzero.
 */
inline static int
suzerain_bsmbsm_spackf(const suzerain_bsmbsm *A, int ihat, int jhat,
                       const float alpha, const float *b, float *papt)
{
    return suzerain_bsmbsm_spack(
            A->S, A->n, ihat, jhat,
            A->kl, A->ku, alpha, b,            A->ld,
            A->KL, A->KU,        papt + A->KL, A->LD + A->KL);
}

/** @copydoc suzerain_bsmbsm_spackf */
inline static int
suzerain_bsmbsm_dpackf(const suzerain_bsmbsm *A, int ihat, int jhat,
                       const double alpha, const double *b, double *papt)
{
    return suzerain_bsmbsm_dpack(
            A->S, A->n, ihat, jhat,
            A->kl, A->ku, alpha, b,            A->ld,
            A->KL, A->KU,        papt + A->KL, A->LD + A->KL);
}

/** @copydoc suzerain_bsmbsm_spackf */
inline static int
suzerain_bsmbsm_cpackf(const suzerain_bsmbsm *A, int ihat, int jhat,
                       const complex_float alpha, const complex_float *b,
                                                        complex_float *papt)
{
    return suzerain_bsmbsm_cpack(
            A->S, A->n, ihat, jhat,
            A->kl, A->ku, alpha, b,            A->ld,
            A->KL, A->KU,        papt + A->KL, A->LD + A->KL);
}

/** @copydoc suzerain_bsmbsm_spackf */
inline static int
suzerain_bsmbsm_zpackf(const suzerain_bsmbsm *A, int ihat, int jhat,
                       const complex_double alpha, const complex_double *b,
                                                         complex_double *papt)
{
    return suzerain_bsmbsm_zpack(
            A->S, A->n, ihat, jhat,
            A->kl, A->ku, alpha, b,            A->ld,
            A->KL, A->KU,        papt + A->KL, A->LD + A->KL);
}

/** @copydoc suzerain_bsmbsm_spackf */
inline static int
suzerain_bsmbsm_cspackf(const suzerain_bsmbsm *A, int ihat, int jhat,
                        const complex_float alpha, const float *b,
                        complex_float *papt)
{
    return suzerain_bsmbsm_cspack(
            A->S, A->n, ihat, jhat,
            A->kl, A->ku, alpha, b,            A->ld,
            A->KL, A->KU,        papt + A->KL, A->LD + A->KL);
}

/** @copydoc suzerain_bsmbsm_spackf */
inline static int
suzerain_bsmbsm_zdpackf(const suzerain_bsmbsm *A, int ihat, int jhat,
                        const complex_double alpha, const double *b,
                        complex_double *papt)
{
    return suzerain_bsmbsm_zdpack(
            A->S, A->n, ihat, jhat,
            A->kl, A->ku, alpha, b,            A->ld,
            A->KL, A->KU,        papt + A->KL, A->LD + A->KL);
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_BSMBSM_H */
