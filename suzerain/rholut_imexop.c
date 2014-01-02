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
 * @copydoc rholut_imexop.h
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/rholut_imexop.h>

#include <suzerain/common.h>
#include <suzerain/blas_et_al.h>
#include <suzerain/bsmbsm.h>

void
suzerain_rholut_imexop_accumulate(
        const complex_double phi,
        const double km,
        const double kn,
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
    // When you modify this routine, you must also modify rholut_imexop.def so
    // that operator accumulation-without-assembly and assembly match.  The
    // test cases in tests/test_rholut_imexop.cpp are invaluable in checking
    // the coherence of these two pieces of logic.
    //
    // When you modify this routine, you may also need to modify
    // rholut_imexop00.c and rholut_imexop00.def so the zero-zero
    // mode case is handled consistently.

    // Sanity checks
    assert(w->nderiv >= 2);            // Adequate workspace?
    assert(!(!in_rho_E ^ !out_rho_E)); // Both non-NULL or NULL?
    assert(!(!in_rho_u ^ !out_rho_u)); // ditto
    assert(!(!in_rho_v ^ !out_rho_v)); // ditto
    assert(!(!in_rho_w ^ !out_rho_w)); // ditto
    assert(!(!in_rho   ^ !out_rho  )); // ditto

    // Prepare shorthand for some useful derived values
    const complex_double ikm = _Complex_I*km;
    const complex_double ikn = _Complex_I*kn;
    const double km2         = km*km;
    const double kn2         = kn*kn;
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
            suzerain_blasext_zgbdddmv_d_z(PREAMBLE_N(M),
                -phi*s->gamma*ikm,           REF(ux),
                -phi*s->gamma*ikn,           REF(uz),
                -phi*ginvRePr*(km2+kn2),     REF(nu),
                w->D_T[M],  w->ld, IN(rho_E), 1.0, OUT(rho_E));

            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D1),
                -phi*s->gamma,               REF(uy),
                w->D_T[D1], w->ld, IN(rho_E), 1.0, OUT(rho_E));

            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D2),
                phi*ginvRePr,                REF(nu),
                w->D_T[D2], w->ld, IN(rho_E), 1.0, OUT(rho_E));
        }

        if (LIKELY(in_rho_u)) {
            const double coeff_nuux
                = Ma2*invRe*((ginvPr-ap43)*km2 + (ginvPr-1)*kn2);
            suzerain_blasext_zgbdddmv_d_z(PREAMBLE_N(M),
                phi*coeff_nuux,              REF(nuux),
                -phi*Ma2*invRe*ap13*km*kn,   REF(nuuz),
                -phi*ikm,                    REF(e_divm),
                w->D_T[M],  w->ld, IN(rho_u), 1.0, OUT(rho_E));

            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D1),
                phi*Ma2*invRe*ap13*ikm,      REF(nuuy),
                w->D_T[D1], w->ld, IN(rho_u), 1.0, OUT(rho_E));

            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D2),
                phi*Ma2*invRe*(1-ginvPr),    REF(nuux),
                w->D_T[D2], w->ld, IN(rho_u), 1.0, OUT(rho_E));
        }

        if (LIKELY(in_rho_v)) {
            const double coeff_nuuy
                = Ma2*invRe*(ginvPr-1)*(km2+kn2);
            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(M),
                phi*coeff_nuuy,              REF(nuuy),
                w->D_T[M],  w->ld, IN(rho_v), 1.0, OUT(rho_E));

            suzerain_blasext_zgbdddmv_d_z(PREAMBLE_N(D1),
                phi*Ma2*invRe*ap13*ikm,      REF(nuux),
                phi*Ma2*invRe*ap13*ikn,      REF(nuuz),
                -phi,                        REF(e_divm),
                w->D_T[D1], w->ld, IN(rho_v), 1.0, OUT(rho_E));

            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D2),
                phi*Ma2*invRe*(ap43-ginvPr), REF(nuuy),
                w->D_T[D2], w->ld, IN(rho_v), 1.0, OUT(rho_E));
        }

        if (LIKELY(in_rho_w)) {
            const double coeff_nuuz
                = Ma2*invRe*((ginvPr-1)*km2 + (ginvPr-ap43)*kn2);
            suzerain_blasext_zgbdddmv_d_z(PREAMBLE_N(M),
                -phi*Ma2*invRe*ap13*km*kn,   REF(nuux),
                phi*coeff_nuuz,              REF(nuuz),
                -phi*ikn,                    REF(e_divm),
                w->D_T[M],  w->ld, IN(rho_w), 1.0, OUT(rho_E));

            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D1),
                phi*Ma2*invRe*ap13*ikn,      REF(nuuy),
                w->D_T[D1], w->ld, IN(rho_w), 1.0, OUT(rho_E));

            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D2),
                phi*Ma2*invRe*(1-ginvPr),    REF(nuuz),
                w->D_T[D2], w->ld, IN(rho_w), 1.0, OUT(rho_E));
        }

        if (LIKELY(in_rho)) {

            // Mass terms done in 2 passes to avoid zgbdddddddmv_d.
            // Writing such a beast may provide a tiny speedup.
            suzerain_blasext_zgbddddmv_d_z(PREAMBLE_N(M),
                phi*Ma2*invRe*(km2+kn2),     REF(nuu2),
                phi*Ma2*invRe*ap13*km2,      REF(nuuxux),
                phi*Ma2*invRe*ap13*2*km*kn,  REF(nuuxuz),
                phi*Ma2*invRe*ap13*kn2,      REF(nuuzuz),
                w->D_T[M],  w->ld, IN(rho), 1.0, OUT(rho_E));
            suzerain_blasext_zgbdddmv_d_z(PREAMBLE_N(M),
                -phi*ikm,                    REF(ex_gradrho),
                -phi*ikn,                    REF(ez_gradrho),
                -phi*ginvRePr/gm1*(km2+kn2), REF(e_deltarho),
                w->D_T[M],  w->ld, IN(rho), 1.0, OUT(rho_E));

            suzerain_blasext_zgbdddmv_d_z(PREAMBLE_N(D1),
                -phi*Ma2*invRe*ap13*2*ikm,   REF(nuuxuy),
                -phi*Ma2*invRe*ap13*2*ikn,   REF(nuuyuz),
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

        if (LIKELY(in_rho_E)) suzerain_blas_zgbmv_d_z(PREAMBLE_NN(M),
                -phi*gm1*invMa2*ikm, w->D_T[M], w->ld, IN(rho_E),
                1.0, OUT(rho_u));

        /* in_rho_u */ {
            suzerain_blasext_zgbdddmv_d_z(PREAMBLE_N(M),
                phi*gm3*ikm,               REF(ux),
                -phi*ikn,                  REF(uz),
                -phi*invRe*(ap43*km2+kn2), REF(nu),
                w->D_T[M],  w->ld, IN(rho_u), 1.0, OUT(rho_u));

            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D1),
                -phi,                      REF(uy),
                w->D_T[D1], w->ld, IN(rho_u), 1.0, OUT(rho_u));

            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D2),
                phi*invRe,                 REF(nu),
                w->D_T[D2], w->ld, IN(rho_u), 1.0, OUT(rho_u));
        }

        if (LIKELY(in_rho_v)) {
            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(M),
                phi*gm1*ikm,               REF(uy),
                w->D_T[M],  w->ld, IN(rho_v), 1.0, OUT(rho_u));

            suzerain_blasext_zgbddmv_d_z(PREAMBLE_N(D1),
                -phi,                      REF(ux),
                phi*ap13*invRe*ikm,        REF(nu),
                w->D_T[D1], w->ld, IN(rho_v), 1.0, OUT(rho_u));
        }

        if (LIKELY(in_rho_w)) {
            suzerain_blasext_zgbdddmv_d_z(PREAMBLE_N(M),
                -phi*ikn,                  REF(ux),
                 phi*gm1*ikm,              REF(uz),
                -phi*ap13*invRe*km*kn,     REF(nu),
                w->D_T[M],  w->ld, IN(rho_w), 1.0, OUT(rho_u));
        }

        if (LIKELY(in_rho)) {
            suzerain_blasext_zgbdddddmv_d_z(PREAMBLE_N(M),
                -phi*0.5*gm1*ikm,          REF(u2),
                phi*ikm,                   REF(uxux),
                phi*ikn,                   REF(uxuz),
                phi*invRe*(ap43*km2+kn2),  REF(nuux),
                phi*ap13*invRe*km*kn,      REF(nuuz),
                w->D_T[M],  w->ld, IN(rho), 1.0, OUT(rho_u));

            suzerain_blasext_zgbddmv_d_z(PREAMBLE_N(D1),
                phi,                       REF(uxuy),
                -phi*ap13*invRe*ikm,       REF(nuuy),
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
            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(M),
                -phi*ikm,             REF(uy),
                w->D_T[M],   w->ld, IN(rho_u), 1.0, OUT(rho_v));

            suzerain_blasext_zgbddmv_d_z(PREAMBLE_N(D1),
                phi*gm1,              REF(ux),
                phi*ap13*invRe*ikm,   REF(nu),
                w->D_T[D1],  w->ld, IN(rho_u), 1.0, OUT(rho_v));
        }

        /* in_rho_v */ {
            suzerain_blasext_zgbdddmv_d_z(PREAMBLE_N(M),
                -phi*ikm,             REF(ux),
                -phi*ikn,             REF(uz),
                -phi*invRe*(km2+kn2), REF(nu),
                w->D_T[M],  w->ld, IN(rho_v), 1.0, OUT(rho_v));

            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D1),
                phi*gm3,              REF(uy),
                w->D_T[D1], w->ld, IN(rho_v), 1.0, OUT(rho_v));

            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D2),
                phi*ap43*invRe,       REF(nu),
                w->D_T[D2], w->ld, IN(rho_v), 1.0, OUT(rho_v));
        }

        if (LIKELY(in_rho_w)) {
            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(M),
                -phi*ikn,             REF(uy),
                w->D_T[M],   w->ld, IN(rho_w), 1.0, OUT(rho_v));

            suzerain_blasext_zgbddmv_d_z(PREAMBLE_N(D1),
                phi*gm1,              REF(uz),
                phi*ap13*invRe*ikn,   REF(nu),
                w->D_T[D1],  w->ld, IN(rho_w), 1.0, OUT(rho_v));
        }

        if (LIKELY(in_rho)) {
            suzerain_blasext_zgbdddmv_d_z(PREAMBLE_N(M),
                phi*ikm,              REF(uxuy),
                phi*ikn,              REF(uyuz),
                phi*invRe*(km2+kn2),  REF(nuuy),
                w->D_T[M],  w->ld, IN(rho), 1.0, OUT(rho_v));

            suzerain_blasext_zgbddddmv_d_z(PREAMBLE_N(D1),
                -phi*0.5*gm1,         REF(u2),
                 phi,                 REF(uyuy),
                -phi*ap13*invRe*ikm,  REF(nuux),
                -phi*ap13*invRe*ikn,  REF(nuuz),
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

        if (LIKELY(in_rho_E)) suzerain_blas_zgbmv_d_z(PREAMBLE_NN(M),
                -phi*gm1*invMa2*ikn, w->D_T[M], w->ld, IN(rho_E),
                1.0, OUT(rho_w));

        if (LIKELY(in_rho_u)) {
            suzerain_blasext_zgbdddmv_d_z(PREAMBLE_N(M),
                 phi*gm1*ikn,              REF(ux),
                -phi*ikm,                  REF(uz),
                -phi*ap13*invRe*km*kn,     REF(nu),
                w->D_T[M],  w->ld, IN(rho_u), 1.0, OUT(rho_w));
        }

        if (LIKELY(in_rho_v)) {
            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(M),
                phi*gm1*ikn,               REF(uy),
                w->D_T[M],  w->ld, IN(rho_v), 1.0, OUT(rho_w));

            suzerain_blasext_zgbddmv_d_z(PREAMBLE_N(D1),
                -phi,                      REF(uz),
                phi*ap13*invRe*ikn,        REF(nu),
                w->D_T[D1], w->ld, IN(rho_v), 1.0, OUT(rho_w));
        }

        /* in_rho_w */ {
            suzerain_blasext_zgbdddmv_d_z(PREAMBLE_N(M),
                -phi*ikm,                  REF(ux),
                phi*gm3*ikn,               REF(uz),
                -phi*invRe*(km2+ap43*kn2), REF(nu),
                w->D_T[M],  w->ld, IN(rho_w), 1.0, OUT(rho_w));

            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D1),
                -phi,                      REF(uy),
                w->D_T[D1], w->ld, IN(rho_w), 1.0, OUT(rho_w));

            suzerain_blasext_zgbdmv_d_z(PREAMBLE_N(D2),
                phi*invRe,                 REF(nu),
                w->D_T[D2], w->ld, IN(rho_w), 1.0, OUT(rho_w));
        }

        if (LIKELY(in_rho)) {
            suzerain_blasext_zgbdddddmv_d_z(PREAMBLE_N(M),
                -phi*0.5*gm1*ikn,          REF(u2),
                phi*ikm,                   REF(uxuz),
                phi*ikn,                   REF(uzuz),
                phi*invRe*(km2+ap43*kn2),  REF(nuuz),
                phi*ap13*invRe*km*kn,      REF(nuux),
                w->D_T[M],  w->ld, IN(rho), 1.0, OUT(rho_w));

            suzerain_blasext_zgbddmv_d_z(PREAMBLE_N(D1),
                 phi,                      REF(uyuz),
                -phi*ap13*invRe*ikn,       REF(nuuy),
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

        if (LIKELY(in_rho_u)) suzerain_blas_zgbmv_d_z(PREAMBLE_NN(M),
                -phi*ikm, w->D_T[M], w->ld, IN(rho_u),  1.0, OUT(rho));

        if (LIKELY(in_rho_v)) suzerain_blas_zgbmv_d_z(PREAMBLE_NN(D1),
                -phi,     w->D_T[D1], w->ld, IN(rho_v), 1.0, OUT(rho));

        if (LIKELY(in_rho_w)) suzerain_blas_zgbmv_d_z(PREAMBLE_NN(M),
                -phi*ikn, w->D_T[M], w->ld, IN(rho_w),  1.0, OUT(rho));

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

// suzerain_rholut_imexop_pack{c,f} differ trivially
// use preprocessor to generate both from the same source template

// suzerain_rholut_imexop_packc
#define FUNCNAME()   suzerain_rholut_imexop_packc
#define PACK(x)      suzerain_bsmbsm_ ## x ## packc
#include "rholut_imexop.def"

// suzerain_rholut_imexop_packf
#define FUNCNAME()   suzerain_rholut_imexop_packf
#define PACK(x)      suzerain_bsmbsm_ ## x ## packf
#include "rholut_imexop.def"
