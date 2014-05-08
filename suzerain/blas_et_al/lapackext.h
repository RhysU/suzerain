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

#ifndef SUZERAIN_BLAS_ET_AL_LAPACKEXT_H
#define SUZERAIN_BLAS_ET_AL_LAPACKEXT_H

/*!\file
 * LAPACK-like extensions built atop LAPACK, the BLAS, and/or on custom coded
 * loops.  Some of these extensions resemble newer LAPACK features which
 * have not trickled down to banded matrices from general matrices.
 */

#include <suzerain/common.h>
#include <suzerain/complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief Compute the solution to a banded system of linear equations \f$ AX=B
 * \f$ using an in-place LU factorization.  This is implemented as a
 * possible call to GBTRF followed by a call to GBTRS.
 *
 * Transposes of \f$ A \f$ can be taken using the \c trans parameter.
 *
 * \param[in,out] fact  One of 'N' or 'F'.  When 'N' on entry, \f$ A \f$ is
 *                      factorized and \c fact is set to 'F' on exit.
 *                      Otherwise, \f$ A \f$ must have been factorized in-place
 *                      by a previous call to this method employing 'N'.
 * \param[in]     trans One of 'N', 'T', or 'C' for no transpose,
 *                      a transpose, or a conjugate transpose, respectively.
 * \param[in]     n
 * \param[in]     kl
 * \param[in]     ku
 * \param[in]     nrhs
 * \param[in,out] ab    Dimension (<tt>ldab</tt>,<tt>n</tt>)
 * \param[in]     ldab  Minimum <tt>kl + ku + 1</tt>
 * \param[in,out] ipiv  Dimension \c n
 * \param[in,out] b     Dimension (<tt>ldb</tt>,<tt>nrhs</tt>)
 * \param[in]     ldb   Minimum \c n
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 */
int
suzerain_lapackext_sgbsv(
        char * const fact,
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        float *ab,
const int ldab,
        int *ipiv,
        float *b,
        const int ldb);

/*! \copydoc suzerain_lapackext_sgbsv */
int
suzerain_lapackext_dgbsv(
        char * const fact,
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        double *ab,
        const int ldab,
        int *ipiv,
        double *b,
        const int ldb);

/*! \copydoc suzerain_lapackext_sgbsv */
int
suzerain_lapackext_cgbsv(
        char * const fact,
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        complex_float *ab,
        const int ldab,
        int *ipiv,
        complex_float *b,
        const int ldb);

/*! \copydoc suzerain_lapackext_sgbsv */
int
suzerain_lapackext_zgbsv(
        char * const fact,
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        complex_double *ab,
        const int ldab,
        int *ipiv,
        complex_double *b,
        const int ldb);

/*!
 * \brief Compute the solution to a banded system of linear equations \f$A x =
 * b \f$ using single precision LU factorization followed by double precision
 * iterative refinement.
 *
 * The overall functionality is an extension of the capability of LAPACK's <a
 * href="http://www.netlib.org/lapack/double/dsgesv.f">DSGESV</a> but written
 * for general banded matrices.  The approach follows DSGESV which is based on
 * the DGESIRSV algorithm presented in June 2006 as <a
 * href="http://www.netlib.org/lapack/lawnspdf/lawn175.pdf">LAWN 175:
 * Exploiting the Performance of 32 bit Floating Point Arithmetic in Obtaining
 * 64 bit Accuracy (Revisiting Iterative Refinement for Linear Systems)</a> by
 * Langou et al.  Unlike their algorithm, this routine
 * <ol>
 *   <li>permits the user to provide a preexisting single
 *       or double precision factorization,</li>
 *   <li>returns the single or double precision factorization
 *       for subsequent reuse,</li>
 *   <li>permits the caller to provide the Frobenius norm of \f$A\f$,</li>
 *   <li>assumes only one right hand side is of interest
 *       (above changes, however, permit efficient repeated invocation),</li>
 *   <li>requires that matrix storage is contiguous (i.e.
 *       <tt>ldab == kl + ku + 1</tt>, <tt>ldafb == 2*kl + ku + 1</tt>),</li>
 *       for reasons of cache-friendliness,</li>
 *   <li>performs precision demotion and promotion in-place
 *       for reasons of cache-friendliness,</li>
 *   <li>performs iterative refinement after the double-precision
 *       fallback LU factorization if the direct solution fails
 *       to meet the solution tolerance,</li>
 *   <li>permits the user to scale the required tolerance by some
 *       amount to accommodate situations where lower precision
 *       results are acceptable,</li>
 *   <li>requires the user to specify the maximum number of refinements
 *       to attempt (Langou et al. suggested thirty, <tt>?gbrfsx</tt> suggests
 *       up to one hundred in aggressive circumstances like providing
 *       approximate factorizations),</li>
 *   <li>makes accessible the number iterative refinements attempted,</li>
 *   <li>and makes accessible the final residual vector and residual.</li>
 * </ol>
 * All in all, these changes are designed to permit combining single precision
 * factorization along with possibly using approximate factorizations and
 * refinement as a means to Newton iteration.  The consistent error bound
 * testing in either single precision factorization or double precision
 * factorization is intended to permit the results to be \e indistinguishable
 * regardless of solution method.
 *
 * \param[in,out] fact  On input, 'N' if factorization must be performed;
 *                      'S' if a single precision factorization is
 *                      provided in \c afb and \c ipiv; 'D' if a double
 *                      precision factorization is provided in \c afb
 *                      and \c ipiv.  On output, either 'S' or 'D'
 *                      indicating which factorization precision has
 *                      been returned in \c afb and \c ipiv.
 * \param[in,out] apprx On input, if \c fact is 'S' or 'D' and \c apprx
 *                      is nonzero, assume the factorization provided in
 *                      \c afb and \c ipiv is \f$ LUP = A + \Delta{}A\f$
 *                      for some perturbation \f$ \Delta{}A \f$.
 *                      If the provided factorization fails to deliver
 *                      adequate convergence as described under \c aiter,
 *                      \c ab will be factorized and \c apprx set to
 *                      zero on return.  If \c fact is 'N' on input,
 *                      \c apprx is ignored and set to zero on return.
 * \param[in]     aiter At most \c aiter refinement steps will be attempted
 *                      before it can be confirmed that the residual
 *                      is decaying by a factor of two at each step.
 *                      Parameter \c aiter is used during mixed-precision
 *                      refinements and during double-precision refinements
 *                      when \c apprx is nonzero.
 * \param[in]     trans If 'N' solve \f$A      x = b\f$.
 *                      If 'T' solve \f$A^\top x = b\f$.
 *                      If 'C' solve \f$A^H    x = b\f$.
 * \param[in]     n     Number of rows and columns in \f$ A \f$.
 * \param[in]     kl    Number of subdiagonals in \f$ A \f$.
 * \param[in]     ku    Number of superdiagonals in \f$ A \f$.
 * \param[in]     ab    Double precision matrix in banded storage of
 *                      dimension (<tt>ldab == kl + ku + 1</tt>, <tt>n</tt>).
 * \param[in,out] afrob The Frobenius norm of \f$ A \f$.
 *                      If negative on entry and \c tolsc is strictly
 *                      positive, the norm is computed and returned to
 *                      permit caching the result.  Otherwise, \c afrob
 *                      is not modified.
 * \param[in,out] afb   Single or double precision LU factorization of \f$A\f$
 *                      of dimension (<tt>ldafb == 2*kl + ku + 1</tt>,
 *                      <tt>n</tt>).  Whether a single or double
 *                      precision result is returned is communicated by
 *                      \c fact on return.
 * \param[in,out] ipiv  LU pivoting information for the factorization
 *                      in \c afb of dimension \c n.
 * \param[in]     b     Right hand side \f$ b \f$ of length \c n.
 * \param[out]    x     Computed solution \f$ x \f$ of length \c n.
 * \param[in,out] siter On input, the maximum number of single precision
 *                      iterative refinements to attempt.  Zero specifies
 *                      no single precision refinements are to occur.
 *                      If negative, no single precision factorization
 *                      or solve will be formed and double precision
 *                      processing occurs.  On output, the number of
 *                      such refinements performed.
 * \param[in,out] diter On input, the maximum number of double precision
 *                      iterative refinements to attempt.  Zero specifies
 *                      no double precision refinements are to occur.
 *                      If negative, no double precision factorization or
 *                      solve will be performed.  On output, the number
 *                      of such refinements performed.
 * \param[in,out] tolsc On input, a nonnegative multiplicative factor
 *                      used to scale the Langou et al. tolerance \f$
 *                      \sqrt{n} \text{eps} \|A\|_\text{Fro} \|x\|_2
 *                      \f$ against which the residual is compared as a
 *                      stopping criterion.  The recommended value from
 *                      Langou et al is \c 1.0 to regain full accuracy
 *                      as measured per backward stability.  The special
 *                      value \c 0.0 can be provided to specify that
 *                      double precision machine epsilon should be used
 *                      which, in conjunction with <tt>siter < 0</tt>,
 *                      makes the refinement process act like LAPACK's
 *                      <tt>?gbrfs</tt>.  On output, the fraction of the
 *                      tolerance represented by the returned solution.
 *                      Values greater than one indicate the desired
 *                      tolerance could not be met.
 * \param[out]    r     Solution residual \f$ b - A x \f$.
 * \param[out]    res   2-norm of the residual.
 *                      That is, \f$\|b - A x\|_2\f$.
 *
 * One or both of \c siter or \c diter must be nonnegative on entry.
 *
 * If \c tolsc is identically zero, \c siter should be strictly negative.
 * This requirement is not enforced, but it is madness to think that
 * mixed-precision iterative refinement will provide residuals better
 * than those returned by a full-precision procedure.
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 *         Errors related to the <tt>i</tt>-th argument are indicated
 *         by a return value of <tt>-i</tt>.
 */
int
suzerain_lapackext_dsgbsvx(
        char * const fact,
        int * const apprx,
        const int aiter,
        char trans,
        const int n,
        const int kl,
        const int ku,
        double * const ab,
        double * const afrob,
        double * const afb,
        int * const ipiv,
        double * const b,
        double * const x,
        int * const siter,
        int * const diter,
        double * const tolsc,
        double * const r,
        double * const res);

/*! \copydoc suzerain_lapackext_dsgbsvx */
int
suzerain_lapackext_zcgbsvx(
        char * const fact,
        int * const apprx,
        const int aiter,
        char trans,
        const int n,
        const int kl,
        const int ku,
        complex_double * const ab,
        double * const afrob,
        complex_double * const afb,
        int * const ipiv,
        complex_double * const b,
        complex_double * const x,
        int * const siter,
        int * const diter,
        double * const tolsc,
        complex_double * const r,
        double * const res);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_BLAS_ET_AL_LAPACKEXT_H */
