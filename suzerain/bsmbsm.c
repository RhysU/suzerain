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

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/bsmbsm.h>

#include <suzerain/common.h>
#include <suzerain/blas_et_al.h>

// Shorthand
#define UNLIKELY(expr) SUZERAIN_UNLIKELY(expr)
static inline int imin(int a, int b) { return a < b ? a : b; }
static inline int imax(int a, int b) { return a > b ? a : b; }

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
