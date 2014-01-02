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
 * Optimized zero-zero mode implementations from rholut_imexop.h.
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/rholut_imexop.h>

#include <suzerain/common.h>
#include <suzerain/blas_et_al.h>
#include <suzerain/bsmbsm.h>

// The special case suzerain_rholut_imexop_accumulate00()
// is distilled from suzerain_rholut_imexop_accumulate()
void
suzerain_rholut_imexop_accumulate00(
        const complex_double phi,
        const suzerain_rholut_imexop_scenario * const s,
        const suzerain_rholut_imexop_ref      * const r,
        const suzerain_rholut_imexop_refld    * const ld,
        const suzerain_bsplineop_workspace    * const w,
        const complex_double *in_rho_E,
        const complex_double *in_rho_u,
        const complex_double *in_rho_v,
        const complex_double *in_rho_w,
        const complex_double *in_rho ,
        const complex_double beta,
        complex_double *out_rho_E,
        complex_double *out_rho_u,
        complex_double *out_rho_v,
        complex_double *out_rho_w,
        complex_double *out_rho )
{
    // The special case suzerain_rholut_imexop_accumulate00()
    // is distilled from suzerain_rholut_imexop_accumulate().
    // See maintenance-related comments at the top of that routine.

    // Sanity checks
    assert(w->nderiv >= 2);            // Adequate workspace?
    assert(!(!in_rho_E ^ !out_rho_E)); // Both non-NULL or NULL?
    assert(!(!in_rho_u ^ !out_rho_u)); // ditto
    assert(!(!in_rho_v ^ !out_rho_v)); // ditto
    assert(!(!in_rho_w ^ !out_rho_w)); // ditto
    assert(!(!in_rho   ^ !out_rho  )); // ditto

    // Prepare shorthand for some useful derived values
    const double gm1         = s->gamma - 1;
    const double gm3         = s->gamma - 3;
    const double ap43        = s->alpha + 4.0/3.0;
    const double ap13        = s->alpha + 1.0/3.0;
    const double Ma2         = s->Ma*s->Ma;
    const double invRe       = 1 / s->Re;
    const double invMa2      = 1 / Ma2;
    const double ginvPr      = s->gamma / s->Pr;
    const double ginvRePr    = s->gamma / (s->Re * s->Pr);

    // Accumulate the requested portions of the M + \varphi L operator.  Scale
    // output by beta, accumulate non-mass contributions, and finally
    // accumulate the mass contributions.  Mass contributions come last as they
    // are expected to have magnitudes much larger than phi-- this helps to
    // ensure better rounding of \varphi L contributions.

    // Readable shorthand for the common code patterns appearing below
    //
    // Notice we need to account for suzerain_bsplineop_workspace storing the
    // transpose of the operators when we invoke suzerain_blaseext_* routines.
    enum { M = 0, D1 = 1, D2 = 2 };
#   define LIKELY(expr)    SUZERAIN_LIKELY(expr)
#   define IN(quantity)    in_##quantity,  1
#   define OUT(quantity)   out_##quantity, 1
#   define REF(quantity)   r->quantity, ld->quantity
#   define PREAMBLE_N(op)  'T', w->n, w->kl[op], w->ku[op]
#   define PREAMBLE_NN(op) 'T', w->n, w->n, w->kl[op], w->ku[op]

    if (LIKELY(in_rho_E)) {  // Accumulate total energy terms into out_rho_E

        suzerain_blas_zscal(w->n, beta, OUT(rho_E));

        /* in_rho_E */ {
            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D1),
                -phi*s->gamma,               REF(uy),
                w->D_T[D1], w->ld, IN(rho_E), 1.0, OUT(rho_E));

            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D2),
                phi*ginvRePr,                REF(nu),
                w->D_T[D2], w->ld, IN(rho_E), 1.0, OUT(rho_E));
        }

        if (LIKELY(in_rho_u)) {
            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D2),
                phi*Ma2*invRe*(1-ginvPr),    REF(nuux),
                w->D_T[D2], w->ld, IN(rho_u), 1.0, OUT(rho_E));
        }

        if (LIKELY(in_rho_v)) {
            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D1),
                -phi,                        REF(e_divm),
                w->D_T[D1], w->ld, IN(rho_v), 1.0, OUT(rho_E));

            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D2),
                phi*Ma2*invRe*(ap43-ginvPr), REF(nuuy),
                w->D_T[D2], w->ld, IN(rho_v), 1.0, OUT(rho_E));
        }

        if (LIKELY(in_rho_w)) {
            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D2),
                phi*Ma2*invRe*(1-ginvPr),    REF(nuuz),
                w->D_T[D2], w->ld, IN(rho_w), 1.0, OUT(rho_E));
        }

        if (LIKELY(in_rho)) {

            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D1),
                -phi,                        REF(ey_gradrho),
                w->D_T[D1], w->ld, IN(rho), 1.0, OUT(rho_E));

            suzerain_blasext_zgbdddmv_d_z(PREAMBLE_N(D2),
                -phi*Ma2*invRe,              REF(nuu2),
                -phi*Ma2*invRe*ap13,         REF(nuuyuy),
                phi*ginvRePr/gm1,            REF(e_deltarho),
                w->D_T[D2], w->ld, IN(rho), 1.0, OUT(rho_E));
        }

        suzerain_blas_zgbmv_d_z(PREAMBLE_NN(M),
            1.0, w->D_T[M], w->ld, IN(rho_E), 1.0, OUT(rho_E));
    }

    if (LIKELY(in_rho_u)) {  // Accumulate X momentum terms into out_rho_u

        suzerain_blas_zscal(w->n, beta, OUT(rho_u));

        if (LIKELY(in_rho_E)) { /* NOP */ }

        /* in_rho_u */ {
            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D1),
                -phi,                      REF(uy),
                w->D_T[D1], w->ld, IN(rho_u), 1.0, OUT(rho_u));

            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D2),
                phi*invRe,                 REF(nu),
                w->D_T[D2], w->ld, IN(rho_u), 1.0, OUT(rho_u));
        }

        if (LIKELY(in_rho_v)) {
            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D1),
                -phi,                      REF(ux),
                w->D_T[D1], w->ld, IN(rho_v), 1.0, OUT(rho_u));
        }

        if (LIKELY(in_rho_w)) { /* NOP */ }

        if (LIKELY(in_rho)) {
            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D1),
                phi,                       REF(uxuy),
                w->D_T[D1], w->ld, IN(rho), 1.0, OUT(rho_u));

            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D2),
                -phi*invRe,                REF(nuux),
                w->D_T[D2], w->ld, IN(rho), 1.0, OUT(rho_u));
        }

        suzerain_blas_zgbmv_d_z(PREAMBLE_NN(M),
            1.0, w->D_T[M], w->ld, IN(rho_u), 1.0, OUT(rho_u));
    }

    if (LIKELY(in_rho_v)) {  // Accumulate Y momentum terms into out_rho_v

        suzerain_blas_zscal(w->n, beta, OUT(rho_v));

        if (LIKELY(in_rho_E)) suzerain_blas_zgbmv_d_z(PREAMBLE_NN(D1),
                -phi*gm1*invMa2, w->D_T[D1], w->ld, IN(rho_E),
                1.0, OUT(rho_v));

        if (LIKELY(in_rho_u)) {
            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D1),
                phi*gm1,              REF(ux),
                w->D_T[D1],  w->ld, IN(rho_u), 1.0, OUT(rho_v));
        }

        /* in_rho_v */ {
            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D1),
                phi*gm3,              REF(uy),
                w->D_T[D1], w->ld, IN(rho_v), 1.0, OUT(rho_v));

            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D2),
                phi*ap43*invRe,       REF(nu),
                w->D_T[D2], w->ld, IN(rho_v), 1.0, OUT(rho_v));
        }

        if (LIKELY(in_rho_w)) {
            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D1),
                phi*gm1,              REF(uz),
                w->D_T[D1],  w->ld, IN(rho_w), 1.0, OUT(rho_v));
        }

        if (LIKELY(in_rho)) {
            suzerain_blasext_zgbddmv_d_z(PREAMBLE_N(D1),
                -phi*0.5*gm1,         REF(u2),
                 phi,                 REF(uyuy),
                w->D_T[D1], w->ld, IN(rho), 1.0, OUT(rho_v));

            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D2),
                -phi*ap43*invRe,      REF(nuuy),
                w->D_T[D2], w->ld, IN(rho), 1.0, OUT(rho_v));
        }

        suzerain_blas_zgbmv_d_z(PREAMBLE_NN(M),
            1.0, w->D_T[M], w->ld, IN(rho_v), 1.0, OUT(rho_v));
    }

    if (LIKELY(in_rho_w)) {  // Accumulate Z momentum terms into out_rho_w

        suzerain_blas_zscal(w->n, beta, OUT(rho_w));

        if (LIKELY(in_rho_E)) { /* NOP */ }

        if (LIKELY(in_rho_u)) { /* NOP */ }

        if (LIKELY(in_rho_v)) {
            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D1),
                -phi,                      REF(uz),
                w->D_T[D1], w->ld, IN(rho_v), 1.0, OUT(rho_w));
        }

        /* in_rho_w */ {
            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D1),
                -phi,                      REF(uy),
                w->D_T[D1], w->ld, IN(rho_w), 1.0, OUT(rho_w));

            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D2),
                phi*invRe,                 REF(nu),
                w->D_T[D2], w->ld, IN(rho_w), 1.0, OUT(rho_w));
        }

        if (LIKELY(in_rho)) {
            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D1),
                 phi,                      REF(uyuz),
                w->D_T[D1], w->ld, IN(rho), 1.0, OUT(rho_w));

            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D2),
                -phi*invRe,                REF(nuuz),
                w->D_T[D2], w->ld, IN(rho), 1.0, OUT(rho_w));
        }

        suzerain_blas_zgbmv_d_z(PREAMBLE_NN(M),
            1.0, w->D_T[M], w->ld, IN(rho_w), 1.0, OUT(rho_w));
    }

    if (LIKELY(in_rho )) {  // Accumulate density terms into out_rho

        suzerain_blas_zscal(w->n, beta, OUT(rho));

        if (LIKELY(in_rho_E)) {/* NOP */}

        if (LIKELY(in_rho_u)) {/* NOP */}

        if (LIKELY(in_rho_v)) suzerain_blas_zgbmv_d_z(PREAMBLE_NN(D1),
                -phi,     w->D_T[D1], w->ld, IN(rho_v), 1.0, OUT(rho));

        if (LIKELY(in_rho_w)) {/* NOP */}

        suzerain_blas_zgbmv_d_z(PREAMBLE_NN(M),
            1.0, w->D_T[M], w->ld, IN(rho), 1.0, OUT(rho));
    }

#   undef PREAMBLE_NN
#   undef PREAMBLE_N
#   undef LIKELY
#   undef REF
#   undef IN
#   undef OUT

}

// suzerain_rholut_imexop_pack{c,f}00 differ trivially
// use preprocessor to generate both from the same source template

// suzerain_rholut_imexop_packc00
#define FUNCNAME00() suzerain_rholut_imexop_packc00
#define PACK(x)      suzerain_bsmbsm_ ## x ## packc
#include "rholut_imexop00.def"

// suzerain_rholut_imexop_packf00
#define FUNCNAME00() suzerain_rholut_imexop_packf00
#define PACK(x)      suzerain_bsmbsm_ ## x ## packf
#include "rholut_imexop00.def"
