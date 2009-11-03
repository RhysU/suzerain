/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2009 The PECOS Development Team
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
 * timestepper.h: SMR91 timestepping implementation
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __PECOS_SUZERAIN_TIMESTEPPER_H
#define __PECOS_SUZERAIN_TIMESTEPPER_H

#include <suzerain/common.h>

/** @file
 * Provides a low-storage, hybrid explicit/implicit Runge-Kutta time integrator
 * suitable for use from both C and Fortran.  The integrator advances the state
 * vector \f$ u(t) \f$ to \f$u(t+\Delta{}t)\f$ according to \f$ u_t = Lu + N(u)
 * \f$ where \f$L\f$ and \f$N\f$ are linear and nonlinear operators,
 * respectively.  Neither operator may depend on time.
 *
 * Different order schemes may be computed using suzerain_lsrk_substep() given
 * appropriate constants in an instance of ::suzerain_lsrk_method.  Currently
 * only ::suzerain_lsrk_smr91 is provided.
 *
 * @see ::suzerain_lsrk_method for the equation used to advance between
 * substeps.
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif
__BEGIN_DECLS
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 * Encapsulates a low-storage, hybrid explicit/implicit Runge-Kutta scheme's
 * coefficients.  Schemes advance between substeps following the equation
 * \f[
 *   u^{i+1} = u^i + \Delta{}t \left[
 *         \alpha_{i} L u^i
 *       + \beta_{i}  L u^{i+1}
 *       + \gamma_{i} N\left( u^{i} \right)
 *       + \zeta_{i} N\left( u^{i} \right)
 *   \right]
 * \f]
 * The coefficients must further obey the relationship \f$ \alpha_{i} +
 * \beta_{i} = \gamma_{i} + \zeta_{i}\f$ so that the implicit and explicit
 * operator substeps are equal.  Index \f$i\f$ runs from 0 to
 * <tt>substeps-1</tt> inclusive.
 */
typedef struct suzerain_lsrk_method {

    /** Human-readable name of the method suitable for output */
    const char * const name;

    /**
     * Number of substeps required to advance \f$u(t)\f$ to
     * \f$u(t+\Delta{}t)\f$.
     */
    const int substeps;

    /** Coefficients \f$\alpha_{i}\f$ in the substep equation. */
    const double * const alpha;

    /** Coefficients \f$\beta{i}\f$ in the substep equation. */
    const double * const beta;

    /** Coefficients \f$\gamma{i}\f$ in the substep equation. */
    const double * const gamma;

    /** Coefficients \f$\zeta{i}\f$ in the substep equation. */
    const double * const zeta;

} suzerain_lsrk_method;

/**
 * Encapsulates the three stage, third order scheme from Appendix A of Spalart,
 * Moser, and Rogers' 1991 ``Spectral Methods for the Navier-Stokes Equations
 * with One Infinite and Two Periodic Directions'' published in the
 * <em>Journal of Computational Physics</em> volume 96 pages 297-324.
 */
extern const suzerain_lsrk_method * const suzerain_lsrk_smr91;

/**
 * Provides a low-storage, hybrid explicit/implicit Runge-Kutta time substepper
 * where the linear operator has the form \f$L = M^{-1} \xi_j D_j\f$ for
 * \f$j\in[0, \dots, \mbox{nD}]\f$ and must be provided using BLAS-compatible
 * general band storage matrices.  This linear operator form facilitates using
 * B-spline derivative operators obtained using bspline.h.  The nonlinear
 * operator must be computed by the caller between Runge-Kutta substeps.
 *
 * More specifically, storage \c a is advanced from \f$u^{i}\f$ to
 * \f$u^{i+1}\f$.  Storage \c b must contain \f$N(u^{i})\f$ on entry; it is
 * preserved for use in a subsequent substep.  Storage \c c must contain
 * \f$N(u^{i-1})\f$ on entry; it is overwritten during the invocation.  Code
 * that uses the substepper should look something like
 * \code
 *      for (int i = 0; i < nsteps; ++i) {
 *          for (int substep = 0;
 *               substep < suzerain_lsrk_smr91->substeps;
 *               ++substep) {
 *              b[0] = nonlinear_operator(a[0]);
 *              if (suzerain_lsrk_substep(
 *                          suzerain_lsrk_smr91,
 *                          n, kl, ku,
 *                          M, ldM,
 *                          nD, xi, D, ldD,
 *                          delta_t, nrhs,
 *                          a, inca, lda,
 *                          b, incb, ldb,
 *                          c, incc, ldc,
 *                          substep)) {
 *                  report_error_somehow();
 *              }
 *              c[0] = b[0];
 *          }
 *      }
 * \endcode
 * It may be more performant to copy \c a into \c b prior to the substep and
 * then compute the nonlinear operator in place.  The caller may also wish to
 * alternate invocations using <tt>a, inca, lda, b, incb, ldb, c, incc,
 * ldc</tt> with ones using <tt>a, inca, lda, c, incc, ldc, b, incb, ldb</tt>
 * rather than copying \c b into \c c after each substep.
 *
 * Multiple, noncontiguous state vectors may be updated within a single
 * substep.  They are parameterized by the number of state vectors (\c nrhs),
 * the leading dimension between state vectors (\c lda, \c ldb, \c ldc), and
 * the increment between elements in each state vector (\c inca, \c incb, \c
 * incc).  For real-valued \f$M\f$, \f$\xi_j\f$, and \f$D_j\f$, these
 * parameters allow computing the nonlinear operator once for a complex field
 * and then separately performing the substep on its real and imaginary
 * subfields using two invocations with different \f$\xi_j\f$ coefficients and
 * imaginary offset storage locations (<tt>a+1</tt>, <tt>b+1</tt>, <tt>c+1</tt>).
 *
 * The implementation employs column-major <a
 * href="http://www.intel.com/software/products/mkl/docs/mklqref/matrixst.htm">
 * general banded matrix</a> BLAS and LAPACK functionality where possible.  The
 * algorithm uses a rearrangement of low-storage substep equation given in the
 * documentation for ::suzerain_lsrk_method under the substitution \f$L =
 * M^{-1} \xi_j D_j\f$:
 * \f[
 *     \underbrace{
 *       \left(M-\sum_{j}\Delta{}t\beta_{i}\xi_{j}D_{j}\right)
 *     }_{\hat{M}}
 *     \left(\frac{u^{i+1} - u^{i}}{\Delta{}t}\right)
 *     =
 *       \underbrace{
 *         \left(\sum_{j}\left(\gamma_i+\zeta_{i-1}\right)\xi_{j}D_{j}\right)
 *       }_{\hat{D}} u^i
 *     + M \left[
 *           \gamma_{i} N\left( u^{i} \right)
 *         + \zeta_{i} N\left( u^{i-1} \right)
 *       \right]
 * \f]
 * The number of state vectors is assumed to be large enough that it is
 * worthwhile to form \f$\hat{D}\f$ prior to computing \f$\hat{D}u^{i}\f$.  The
 * computational algorithm is
 * -# Preconditions:
 *   Storage \f$a = u^i\f$;
 *   storage \f$b = N\left(u^{i}\right)\f$;
 *   storage \f$c = N\left(u^{i-1}\right)\f$
 * -# Allocate \f$d\f$ to be the length of one state vector
 * -# GB_ACC: Form operator \f$\hat{D}\f$ using \f$\gamma_i\f$,
 *           \f$\zeta_{i}\f$, \f$\xi_j\f$, and \f$D_j\f$
 * -# GB_ACC: Form operator \f$\hat{M}\f$ using \f$M\f$,
 *           \f$\Delta{}t\f$, \f$\beta_i\f$, \f$\xi_j\f$, and \f$D_j\f$
 * -# GBTRF:  \f$\hat{M}\leftarrow{}\mbox{lu}\left( \hat{M} \right)\f$
 * -# For each right hand side \f$a_j\f$, \f$b_j\f$, and \f$c_j\f$
 *      in storage \c a, \c b, and \c c, respectively:
 *   -# AXPBY: \f$c_{j}\leftarrow{}\gamma_{i}b_{j}+\zeta_{i}c_{j}\f$
 *   -# GBMV: \f$d\leftarrow{}1 M c_{j} + 0 d\f$
 *   -# GBMV: \f$d\leftarrow{}1 \hat{D} a_{j} + 1 d\f$
 *   -# GBTRS: \f$d\leftarrow{}\hat{M}^{-1} d\f$
 *   -# AXPY: \f$a_{j}\leftarrow{}\Delta{}t\,d + a_{j}\f$
 * -# Deallocate \f$d\f$
 * -# Postconditions:
 *   Storage \f$a = u^{i+1}\f$;
 *   storage \f$b = N\left(u^{i}\right)\f$
 * .
 * Note that the routine chooses to assemble the right hand side of the
 * equation one state vector at a time because
 * - GBMV can only operate on a single vector per invocation, and because
 * - GBTRS does not allow an arbitrary increment between elements in a vector.
 *
 * @param[in] method Low-storage method to use.
 * @param[in] n Number of rows and columns in \f$M\f$ and \f$D_j\f$.
 * @param[in] kl Number of subdiagonals in \f$M\f$ and \f$D_j\f$.
 * @param[in] ku Number of superdiagonals in \f$M\f$ and \f$D_j\f$.
 * @param[in] M Column-major band storage of \f$M\f$ using \c kl, \c ku,
 *      and \c ldM.
 * @param[in] ldM Leading dimension between columns of \c M.
 * @param[in] nD Number of operators stored in \c D.
 * @param[in] xi Leading coefficients \f$\xi_j\f$ for operators \f$D_j\f$
 *      indexed from \c 0 to <tt>nD-1</tt> inclusive.
 * @param[in] D Column-major band storage of operators \f$D_j\f$ using
 *      \c kl, \c ku, and \c ldD.  Operator \f$D_0\f$ is stored in
 *      <tt>D[0]</tt>, etc.  Indices run from \c 0 to <tt>nD-1</tt>
 *      inclusive.
 * @param[in] ldD Leading dimension between columns of each operator
 *      <tt>D[j]</tt>.
 * @param[in] delta_t Size of the Runge-Kutta step to take.  The same
 *      value must be provided for all substeps when advancing from
 *      \f$u(t)\f$ to \f$u(t+\Delta{}t)\f$.
 * @param[in] nrhs The number of right hand sides or state vectors stored
 *      in \c a, \c b, and \c c.
 * @param[in,out] a Contains \f$u^{i}\f$ on entry and \f$u^{i+1}\f$ on
 *      exit where \f$i\f$ is the provided \c substep.
 * @param[in] inca Increment between two consecutive elements in a
 *      column vector in \c a.
 * @param[in] lda Leading dimension between two state vectors in \c a.
 *      That is, \f$a_j\f$ starts at <tt>a+j*inca</tt>.
 * @param[in] b Contains \f$N(u^{i})\f$ on entry where \f$i\f$ is
 *      the provided \c substep.
 * @param[in] incb Increment between two consecutive elements in a
 *      column vector in \c b.
 * @param[in] ldb Leading dimension between two state vectors in \c b.
 *      That is, \f$b_j\f$ starts at <tt>b+j*incb</tt>.
 * @param[in,out] c Contains \f$N(u^{i+1})\f$ on entry and contains
 *      undefined values on exit where \f$i\f$ is the provided \c substep.
 * @param[in] incc Increment between two consecutive elements in a
 *      column vector in \c c.
 * @param[in] ldc Leading dimension between two state vectors in \c c.
 *      That is, \f$c_j\f$ starts at <tt>b+j*incc</tt>.
 * @param[in] substep The substep to take which must be in the inclusive
 *      range 0 to <tt>method.substeps</tt>.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_lsrk_substep(
        const suzerain_lsrk_method * const method,
        const int n,
        const int kl,
        const int ku,
        const double * const M,
        const int ldM,
        const int nD,
        const double * const xi,
        const double * const * const D,
        const int ldD,
        double delta_t,
        const int nrhs,
        double       * const a, const int inca, const int lda,
        const double * const b, const int incb, const int ldb,
        double       * const c, const int incc, const int ldc,
        const int substep);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
__END_DECLS
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif // __PECOS_SUZERAIN_TIMESTEPPER_H
