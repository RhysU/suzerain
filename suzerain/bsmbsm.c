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

/** @file
 * @copydoc bsmbsm.h
 */

#include <suzerain/bsmbsm.h>

#include <suzerain/common.h>
#include <suzerain/blas_et_al.h>

// Shorthand
#define UNLIKELY(expr) SUZERAIN_UNLIKELY(expr)
static inline int imin(int a, int b) { return a < b ? a : b; }
static inline int imax(int a, int b) { return a > b ? a : b; }

// C99 extern declaration for static inline function
extern
suzerain_bsmbsm suzerain_bsmbsm_construct(int S, int n, int kl, int ku);

// C99 extern declaration for static inline function
extern
int suzerain_bsmbsm_q(int S, int n, int i);

// C99 extern declaration for static inline function
extern
int suzerain_bsmbsm_qinv(int S, int n, int i);

// suzerain_bsmbsm_saPxpby
#define COMPONENT float
#define AFFIXPREC(pre,post) pre ## s ## post
#include "bsmbsm_aPxpby_real.def"

// suzerain_bsmbsm_daPxpby
#define COMPONENT double
#define AFFIXPREC(pre,post) pre ## d ## post
#include "bsmbsm_aPxpby_real.def"

// suzerain_bsmbsm_caPxpby
#define COMPONENT float
#define COMPLEX   complex_float
#define AFFIXPREC(pre,post) pre ## c ## post
#include "bsmbsm_aPxpby_complex.def"

// suzerain_bsmbsm_zaPxpby
#define COMPONENT double
#define COMPLEX   complex_double
#define AFFIXPREC(pre,post) pre ## z ## post
#include "bsmbsm_aPxpby_complex.def"

// suzerain_bsmbsm_spack
#define SOURCE float
#define TARGET float
#define AFFIXPREC(pre,post) pre ## s ## post
#include "bsmbsm_pack.def"

// suzerain_bsmbsm_dpack
#define SOURCE double
#define TARGET double
#define AFFIXPREC(pre,post) pre ## d ## post
#include "bsmbsm_pack.def"

// suzerain_bsmbsm_cpack
#define SOURCE complex_float
#define TARGET complex_float
#define AFFIXPREC(pre,post) pre ## c ## post
#include "bsmbsm_pack.def"

// suzerain_bsmbsm_zpack
#define SOURCE complex_double
#define TARGET complex_double
#define AFFIXPREC(pre,post) pre ## z ## post
#include "bsmbsm_pack.def"

// suzerain_bsmbsm_cspack
#define SOURCE float
#define TARGET complex_float
#define AFFIXPREC(pre,post) pre ## cs ## post
#include "bsmbsm_pack.def"

// suzerain_bsmbsm_zdpack
#define SOURCE double
#define TARGET complex_double
#define AFFIXPREC(pre,post) pre ## zd ## post
#include "bsmbsm_pack.def"

gsl_permutation *
suzerain_bsmbsm_permutation(int S, int n)
{
    assert(S > 0);
    assert(n > 0);

    const ptrdiff_t N = S*n;
    gsl_permutation * const p = gsl_permutation_alloc(N);
    if (p) {
        size_t * const q = p->data;
        for (ptrdiff_t i = 0; i < N; ++i) {
            q[i] = suzerain_bsmbsm_q(S, n, i);
        }
    }
    return p;
}

// C99 extern declaration for static inline function
extern int
suzerain_bsmbsm_spackc(const suzerain_bsmbsm *A, int ihat, int jhat,
                       const float alpha, const float *b, float *papt);

// C99 extern declaration for static inline function
extern int
suzerain_bsmbsm_dpackc(const suzerain_bsmbsm *A, int ihat, int jhat,
                       const double alpha, const double *b, double *papt);

// C99 extern declaration for static inline function
extern int
suzerain_bsmbsm_cpackc(const suzerain_bsmbsm *A, int ihat, int jhat,
                       const complex_float alpha, const complex_float *b,
                                                        complex_float *papt);

// C99 extern declaration for static inline function
extern int
suzerain_bsmbsm_zpackc(const suzerain_bsmbsm *A, int ihat, int jhat,
                       const complex_double alpha, const complex_double *b,
                                                         complex_double *papt);

// C99 extern declaration for static inline function
extern int
suzerain_bsmbsm_cspackc(const suzerain_bsmbsm *A, int ihat, int jhat,
                        const complex_float alpha, const float *b,
                        complex_float *papt);

// C99 extern declaration for static inline function
extern int
suzerain_bsmbsm_zdpackc(const suzerain_bsmbsm *A, int ihat, int jhat,
                        const complex_double alpha, const double *b,
                        complex_double *papt);

// C99 extern declaration for static inline function
extern float*
suzerain_bsmbsm_spackc_ptr(const suzerain_bsmbsm *A, int qinv_i, int qinv_j,
                           float *papt);

// C99 extern declaration for static inline function
extern double*
suzerain_bsmbsm_dpackc_ptr(const suzerain_bsmbsm *A, int qinv_i, int qinv_j,
                           double *papt);

// C99 extern declaration for static inline function
extern complex_float*
suzerain_bsmbsm_cpackc_ptr(const suzerain_bsmbsm *A, int qinv_i, int qinv_j,
                           complex_float *papt);

// C99 extern declaration for static inline function
extern complex_double*
suzerain_bsmbsm_zpackc_ptr(const suzerain_bsmbsm *A, int qinv_i, int qinv_j,
                           complex_double *papt);

// C99 extern declaration for static inline function
extern int
suzerain_bsmbsm_spackf(const suzerain_bsmbsm *A, int ihat, int jhat,
                       const float alpha, const float *b, float *papt);

// C99 extern declaration for static inline function
extern int
suzerain_bsmbsm_dpackf(const suzerain_bsmbsm *A, int ihat, int jhat,
                       const double alpha, const double *b, double *papt);

// C99 extern declaration for static inline function
extern int
suzerain_bsmbsm_cpackf(const suzerain_bsmbsm *A, int ihat, int jhat,
                       const complex_float alpha, const complex_float *b,
                                                        complex_float *papt);

// C99 extern declaration for static inline function
extern int
suzerain_bsmbsm_zpackf(const suzerain_bsmbsm *A, int ihat, int jhat,
                       const complex_double alpha, const complex_double *b,
                                                         complex_double *papt);

// C99 extern declaration for static inline function
extern int
suzerain_bsmbsm_cspackf(const suzerain_bsmbsm *A, int ihat, int jhat,
                        const complex_float alpha, const float *b,
                        complex_float *papt);

// C99 extern declaration for static inline function
extern int
suzerain_bsmbsm_zdpackf(const suzerain_bsmbsm *A, int ihat, int jhat,
                        const complex_double alpha, const double *b,
                        complex_double *papt);

// C99 extern declaration for static inline function
extern float*
suzerain_bsmbsm_spackf_ptr(const suzerain_bsmbsm *A, int qinv_i, int qinv_j,
                           float *papt);

// C99 extern declaration for static inline function
extern double*
suzerain_bsmbsm_dpackf_ptr(const suzerain_bsmbsm *A, int qinv_i, int qinv_j,
                           double *papt);

// C99 extern declaration for static inline function
extern complex_float*
suzerain_bsmbsm_cpackf_ptr(const suzerain_bsmbsm *A, int qinv_i, int qinv_j,
                           complex_float *papt);

// C99 extern declaration for static inline function
extern complex_double*
suzerain_bsmbsm_zpackf_ptr(const suzerain_bsmbsm *A, int qinv_i, int qinv_j,
                           complex_double *papt);
