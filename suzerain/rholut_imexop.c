/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2012 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *--------------------------------------------------------------------------
 * rholut_imexop.c: hybrid implicit/explicit operator apply and solve
 * $Id$
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/rholut_imexop.h>
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
        const complex_double *in_rho ,
        const complex_double *in_rhou,
        const complex_double *in_rhov,
        const complex_double *in_rhow,
        const complex_double *in_rhoe,
        const complex_double beta,
        complex_double *out_rho ,
        complex_double *out_rhou,
        complex_double *out_rhov,
        complex_double *out_rhow,
        complex_double *out_rhoe)
{
    // When you modify this routine, you must also modify rholut_imexop.def so
    // that operator accumulation-without-assembly and assembly match.  The
    // test cases in tests/test_rholut_imexop.cpp are invaluable in checking
    // the coherence of these two pieces of logic.

    // Sanity checks
    assert(w->nderiv >= 2);          // Adequate workspace?
    assert(!(!in_rho  ^ !out_rho )); // Both non-NULL or NULL?
    assert(!(!in_rhou ^ !out_rhou)); // ditto
    assert(!(!in_rhov ^ !out_rhov)); // ditto
    assert(!(!in_rhow ^ !out_rhow)); // ditto
    assert(!(!in_rhoe ^ !out_rhoe)); // ditto

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

    // We need to account for suzerain_bsplineop_workspace storing the
    // transpose of the operators when we invoke suzerain_blaseext_* routines.
    static const char trans = 'T';

    // Prepare several oft-used constants to aid readability
    static const int M       = 0;
    static const int D1      = 1;
    static const int D2      = 2;
    const int        n       = w->n;

    // Shorthand for the common pattern of providing a "r->foo, ld->foo" pair.
#   define REF(quantity) r->quantity, ld->quantity

    // Shorthand for the common pattern of providing "in_foo, inc" pairs.
#   define IN(quantity)  in_##quantity,  1
#   define OUT(quantity) out_##quantity, 1

    if (in_rho ) {  // Accumulate density terms into out_rho

        suzerain_blas_zscal(n, beta, OUT(rho));

        if (in_rhou) suzerain_blas_zgbmv_d(trans, n, n, w->kl[M],  w->ku[M],
                -phi*ikm, w->D_T[M], w->ld, IN(rhou),  1.0, OUT(rho));

        if (in_rhov) suzerain_blas_zgbmv_d(trans, n, n, w->kl[D1], w->ku[D1],
                -phi,     w->D_T[D1], w->ld, IN(rhov), 1.0, OUT(rho));

        if (in_rhow) suzerain_blas_zgbmv_d(trans, n, n, w->kl[M],  w->ku[M],
                -phi*ikn, w->D_T[M], w->ld, IN(rhow),  1.0, OUT(rho));

        if (in_rhoe) {/* NOP */};

        suzerain_blas_zgbmv_d(trans, n, n, w->kl[M], w->ku[M],
            1.0, w->D_T[M], w->ld, IN(rho), 1.0, OUT(rho));
    }

    if (in_rhou) {  // Accumulate X momentum terms into out_rhou

        suzerain_blas_zscal(n, beta, OUT(rhou));

        if (in_rho) {
            suzerain_blasext_zgbdddddmv_d(trans, n, w->kl[M], w->ku[M],
                -phi*0.5*gm1*ikm,          REF(u2),
                phi*ikm,                   REF(uxux),
                phi*ikn,                   REF(uxuz),
                phi*invRe*(ap43*km2+kn2),  REF(nuux),
                phi*ap13*invRe*km*kn,      REF(nuuz),
                w->D_T[M],  w->ld, IN(rho), 1.0, OUT(rhou));

            suzerain_blasext_zgbddmv_d(trans, n, w->kl[D1], w->ku[D1],
                phi,                       REF(uxuy),
                -phi*ap13*invRe*ikm,       REF(nuuy),
                w->D_T[D1], w->ld, IN(rho), 1.0, OUT(rhou));

            suzerain_blasext_zgbdmv_d(trans, n, w->kl[D2], w->ku[D2],
                -phi*invRe,                REF(nuux),
                w->D_T[D2], w->ld, IN(rho), 1.0, OUT(rhou));
        }

        /* in_rhou */ {
            suzerain_blasext_zgbdddmv_d(trans, n, w->kl[M], w->ku[M],
                phi*gm3*ikm,               REF(ux),
                -phi*ikn,                  REF(uz),
                -phi*invRe*(ap43*km2+kn2), REF(nu),
                w->D_T[M],  w->ld, IN(rhou), 1.0, OUT(rhou));

            suzerain_blasext_zgbdmv_d(trans, n, w->kl[D1], w->ku[D1],
                -phi,                      REF(uy),
                w->D_T[D1], w->ld, IN(rhou), 1.0, OUT(rhou));

            suzerain_blasext_zgbdmv_d(trans, n, w->kl[D2], w->ku[D2],
                phi*invRe,                 REF(nu),
                w->D_T[D2], w->ld, IN(rhou), 1.0, OUT(rhou));
        }

        if (in_rhov) {
            suzerain_blasext_zgbdmv_d(trans, n, w->kl[M], w->ku[M],
                phi*gm1*ikm,               REF(uy),
                w->D_T[M],  w->ld, IN(rhov), 1.0, OUT(rhou));

            suzerain_blasext_zgbddmv_d(trans, n, w->kl[D1], w->ku[D1],
                -phi,                      REF(ux),
                phi*ap13*invRe*ikm,        REF(nu),
                w->D_T[D1], w->ld, IN(rhov), 1.0, OUT(rhou));
        }

        if (in_rhow) {
            suzerain_blasext_zgbdddmv_d(trans, n, w->kl[M], w->ku[M],
                -phi*ikn,                  REF(ux),
                 phi*gm1*ikm,              REF(uz),
                -phi*ap13*invRe*km*kn,     REF(nu),
                w->D_T[M],  w->ld, IN(rhow), 1.0, OUT(rhou));
        }

        if (in_rhoe) suzerain_blas_zgbmv_d(trans, n, n, w->kl[M], w->ku[M],
                -phi*gm1*invMa2*ikm, w->D_T[M], w->ld, IN(rhoe),
                1.0, OUT(rhou));

        suzerain_blas_zgbmv_d(trans, n, n, w->kl[M], w->ku[M],
            1.0, w->D_T[M], w->ld, IN(rhou), 1.0, OUT(rhou));
    }

    if (in_rhov) {  // Accumulate Y momentum terms into out_rhov

        suzerain_blas_zscal(n, beta, OUT(rhov));

        if (in_rho) {
            suzerain_blasext_zgbdddmv_d(trans, n, w->kl[M], w->ku[M],
                phi*ikm,              REF(uxuy),
                phi*ikn,              REF(uyuz),
                phi*invRe*(km2+kn2),  REF(nuuy),
                w->D_T[M],  w->ld, IN(rho), 1.0, OUT(rhov));

            suzerain_blasext_zgbddddmv_d(trans, n, w->kl[D1], w->ku[D1],
                -phi*0.5*gm1,         REF(u2),
                 phi,                 REF(uyuy),
                -phi*ap13*invRe*ikm,  REF(nuux),
                -phi*ap13*invRe*ikn,  REF(nuuz),
                w->D_T[D1], w->ld, IN(rho), 1.0, OUT(rhov));

            suzerain_blasext_zgbdmv_d(trans, n, w->kl[D2], w->ku[D2],
                -phi*ap43*invRe,      REF(nuuy),
                w->D_T[D2], w->ld, IN(rho), 1.0, OUT(rhov));
        }

        if (in_rhou) {
            suzerain_blasext_zgbdmv_d(trans, n, w->kl[M], w->ku[M],
                -phi*ikm,             REF(uy),
                w->D_T[M],   w->ld, IN(rhou), 1.0, OUT(rhov));

            suzerain_blasext_zgbddmv_d(trans, n, w->kl[D1], w->ku[D1],
                phi*gm1,              REF(ux),
                phi*ap13*invRe*ikm,   REF(nu),
                w->D_T[D1],  w->ld, IN(rhou), 1.0, OUT(rhov));
        }

        /* in_rhov */ {
            suzerain_blasext_zgbdddmv_d(trans, n, w->kl[M], w->ku[M],
                -phi*ikm,             REF(ux),
                -phi*ikn,             REF(uz),
                -phi*invRe*(km2+kn2), REF(nu),
                w->D_T[M],  w->ld, IN(rhov), 1.0, OUT(rhov));

            suzerain_blasext_zgbdmv_d(trans, n, w->kl[D1], w->ku[D1],
                phi*gm3,              REF(uy),
                w->D_T[D1], w->ld, IN(rhov), 1.0, OUT(rhov));

            suzerain_blasext_zgbdmv_d(trans, n, w->kl[D2], w->ku[D2],
                phi*ap43*invRe,       REF(nu),
                w->D_T[D2], w->ld, IN(rhov), 1.0, OUT(rhov));
        }

        if (in_rhow) {
            suzerain_blasext_zgbdmv_d(trans, n, w->kl[M], w->ku[M],
                -phi*ikn,             REF(uy),
                w->D_T[M],   w->ld, IN(rhow), 1.0, OUT(rhov));

            suzerain_blasext_zgbddmv_d(trans, n, w->kl[D1], w->ku[D1],
                phi*gm1,              REF(uz),
                phi*ap13*invRe*ikn,   REF(nu),
                w->D_T[D1],  w->ld, IN(rhow), 1.0, OUT(rhov));
        }

        if (in_rhoe) suzerain_blas_zgbmv_d(trans, n, n, w->kl[D1], w->ku[D1],
                -phi*gm1*invMa2, w->D_T[D1], w->ld, IN(rhoe),
                1.0, OUT(rhov));

        suzerain_blas_zgbmv_d(trans, n, n, w->kl[M], w->ku[M],
            1.0, w->D_T[M], w->ld, IN(rhov), 1.0, OUT(rhov));
    }

    if (in_rhow) {  // Accumulate Z momentum terms into out_rhow

        suzerain_blas_zscal(n, beta, OUT(rhow));

        if (in_rho) {
            suzerain_blasext_zgbdddddmv_d(trans, n, w->kl[M], w->ku[M],
                -phi*0.5*gm1*ikn,          REF(u2),
                phi*ikm,                   REF(uxuz),
                phi*ikn,                   REF(uzuz),
                phi*invRe*(km2+ap43*kn2),  REF(nuuz),
                phi*ap13*invRe*km*kn,      REF(nuux),
                w->D_T[M],  w->ld, IN(rho), 1.0, OUT(rhow));

            suzerain_blasext_zgbddmv_d(trans, n, w->kl[D1], w->ku[D1],
                 phi,                      REF(uyuz),
                -phi*ap13*invRe*ikn,       REF(nuuy),
                w->D_T[D1], w->ld, IN(rho), 1.0, OUT(rhow));

            suzerain_blasext_zgbdmv_d(trans, n, w->kl[D2], w->ku[D2],
                -phi*invRe,                REF(nuuz),
                w->D_T[D2], w->ld, IN(rho), 1.0, OUT(rhow));
        }

        if (in_rhou) {
            suzerain_blasext_zgbdddmv_d(trans, n, w->kl[M], w->ku[M],
                 phi*gm1*ikn,              REF(ux),
                -phi*ikm,                  REF(uz),
                -phi*ap13*invRe*km*kn,     REF(nu),
                w->D_T[M],  w->ld, IN(rhou), 1.0, OUT(rhow));
        }

        if (in_rhov) {
            suzerain_blasext_zgbdmv_d(trans, n, w->kl[M], w->ku[M],
                phi*gm1*ikn,               REF(uy),
                w->D_T[M],  w->ld, IN(rhov), 1.0, OUT(rhow));

            suzerain_blasext_zgbddmv_d(trans, n, w->kl[D1], w->ku[D1],
                -phi,                      REF(uz),
                phi*ap13*invRe*ikn,        REF(nu),
                w->D_T[D1], w->ld, IN(rhov), 1.0, OUT(rhow));
        }

        /* in_rhow */ {
            suzerain_blasext_zgbdddmv_d(trans, n, w->kl[M], w->ku[M],
                -phi*ikm,                  REF(ux),
                phi*gm3*ikn,               REF(uz),
                -phi*invRe*(km2+ap43*kn2), REF(nu),
                w->D_T[M],  w->ld, IN(rhow), 1.0, OUT(rhow));

            suzerain_blasext_zgbdmv_d(trans, n, w->kl[D1], w->ku[D1],
                -phi,                      REF(uy),
                w->D_T[D1], w->ld, IN(rhow), 1.0, OUT(rhow));

            suzerain_blasext_zgbdmv_d(trans, n, w->kl[D2], w->ku[D2],
                phi*invRe,                 REF(nu),
                w->D_T[D2], w->ld, IN(rhow), 1.0, OUT(rhow));
        }

        if (in_rhoe) suzerain_blas_zgbmv_d(trans, n, n, w->kl[M], w->ku[M],
                -phi*gm1*invMa2*ikn, w->D_T[M], w->ld, IN(rhoe),
                1.0, OUT(rhow));

        suzerain_blas_zgbmv_d(trans, n, n, w->kl[M], w->ku[M],
            1.0, w->D_T[M], w->ld, IN(rhow), 1.0, OUT(rhow));
    }

    if (in_rhoe) {  // Accumulate total energy terms into out_rhoe

        suzerain_blas_zscal(n, beta, OUT(rhoe));

        if (in_rho) {

            // Mass terms done in 2 passes to avoid zgbdddddddmv_d.
            // Writing such a beast may provide a tiny speedup.
            suzerain_blasext_zgbddddmv_d(trans, n, w->kl[M], w->ku[M],
                phi*Ma2*invRe*(km2+kn2),     REF(nuu2),
                phi*Ma2*invRe*ap13*km2,      REF(nuuxux),
                phi*Ma2*invRe*ap13*2*km*kn,  REF(nuuxuz),
                phi*Ma2*invRe*ap13*kn2,      REF(nuuzuz),
                w->D_T[M],  w->ld, IN(rho), 1.0, OUT(rhoe));
            suzerain_blasext_zgbdddmv_d( trans, n, w->kl[M], w->ku[M],
                -phi*ikm,                    REF(ex_gradrho),
                -phi*ikn,                    REF(ez_gradrho),
                -phi*ginvRePr/gm1*(km2+kn2), REF(e_deltarho),
                w->D_T[M],  w->ld, IN(rho), 1.0, OUT(rhoe));

            suzerain_blasext_zgbdddmv_d(trans, n, w->kl[D1], w->ku[D1],
                -phi*Ma2*invRe*ap13*2*ikm,   REF(nuuxuy),
                -phi*Ma2*invRe*ap13*2*ikn,   REF(nuuyuz),
                -phi,                        REF(ey_gradrho),
                w->D_T[D1], w->ld, IN(rho), 1.0, OUT(rhoe));

            suzerain_blasext_zgbdddmv_d(trans, n, w->kl[D2], w->ku[D2],
                -phi*Ma2*invRe,              REF(nuu2),
                -phi*Ma2*invRe*ap13,         REF(nuuyuy),
                phi*ginvRePr/gm1,            REF(e_deltarho),
                w->D_T[D2], w->ld, IN(rho), 1.0, OUT(rhoe));
        }

        if (in_rhou) {
            const double coeff_nuux
                = Ma2*invRe*((ginvPr-ap43)*km2 + (ginvPr-1)*kn2);
            suzerain_blasext_zgbdddmv_d(trans, n, w->kl[M], w->ku[M],
                phi*coeff_nuux,              REF(nuux),
                -phi*Ma2*invRe*ap13*km*kn,   REF(nuuz),
                -phi*ikm,                    REF(e_divm),
                w->D_T[M],  w->ld, IN(rhou), 1.0, OUT(rhoe));

            suzerain_blasext_zgbdmv_d(trans, n, w->kl[D1], w->ku[D1],
                phi*Ma2*invRe*ap13*ikm,      REF(nuuy),
                w->D_T[D1], w->ld, IN(rhou), 1.0, OUT(rhoe));

            suzerain_blasext_zgbdmv_d(trans, n, w->kl[D2], w->ku[D2],
                phi*Ma2*invRe*(1-ginvPr),    REF(nuux),
                w->D_T[D2], w->ld, IN(rhou), 1.0, OUT(rhoe));
        }

        if (in_rhov) {
            const double coeff_nuuy
                = Ma2*invRe*(ginvPr-1)*(km2+kn2);
            suzerain_blasext_zgbdmv_d(trans, n, w->kl[M], w->ku[M],
                phi*coeff_nuuy,              REF(nuuy),
                w->D_T[M],  w->ld, IN(rhov), 1.0, OUT(rhoe));

            suzerain_blasext_zgbdddmv_d(trans, n, w->kl[D1], w->ku[D1],
                phi*Ma2*invRe*ap13*ikm,      REF(nuux),
                phi*Ma2*invRe*ap13*ikn,      REF(nuuz),
                -phi,                        REF(e_divm),
                w->D_T[D1], w->ld, IN(rhov), 1.0, OUT(rhoe));

            suzerain_blasext_zgbdmv_d(trans, n, w->kl[D2], w->ku[D2],
                phi*Ma2*invRe*(ap43-ginvPr), REF(nuuy),
                w->D_T[D2], w->ld, IN(rhov), 1.0, OUT(rhoe));
        }

        if (in_rhow) {
            const double coeff_nuuz
                = Ma2*invRe*((ginvPr-1)*km2 + (ginvPr-ap43)*kn2);
            suzerain_blasext_zgbdddmv_d(trans, n, w->kl[M], w->ku[M],
                -phi*Ma2*invRe*ap13*km*kn,   REF(nuux),
                phi*coeff_nuuz,              REF(nuuz),
                -phi*ikn,                    REF(e_divm),
                w->D_T[M],  w->ld, IN(rhow), 1.0, OUT(rhoe));

            suzerain_blasext_zgbdmv_d(trans, n, w->kl[D1], w->ku[D1],
                phi*Ma2*invRe*ap13*ikn,      REF(nuuy),
                w->D_T[D1], w->ld, IN(rhow), 1.0, OUT(rhoe));

            suzerain_blasext_zgbdmv_d(trans, n, w->kl[D2], w->ku[D2],
                phi*Ma2*invRe*(1-ginvPr),    REF(nuuz),
                w->D_T[D2], w->ld, IN(rhow), 1.0, OUT(rhoe));
        }

        /* in_rhoe */ {
            suzerain_blasext_zgbdddmv_d(trans, n, w->kl[M], w->ku[M],
                -phi*s->gamma*ikm,           REF(ux),
                -phi*s->gamma*ikn,           REF(uz),
                -phi*ginvRePr*(km2+kn2),     REF(nu),
                w->D_T[M],  w->ld, IN(rhoe), 1.0, OUT(rhoe));

            suzerain_blasext_zgbdmv_d(trans, n, w->kl[D1], w->ku[D1],
                -phi*s->gamma,               REF(uy),
                w->D_T[D1], w->ld, IN(rhoe), 1.0, OUT(rhoe));

            suzerain_blasext_zgbdmv_d(trans, n, w->kl[D2], w->ku[D2],
                phi*ginvRePr,                REF(nu),
                w->D_T[D2], w->ld, IN(rhoe), 1.0, OUT(rhoe));
        }

        suzerain_blas_zgbmv_d(trans, n, n, w->kl[M], w->ku[M],
            1.0, w->D_T[M], w->ld, IN(rhoe), 1.0, OUT(rhoe));
    }

#   undef REF
#   undef IN
#   undef OUT

}

// suzerain_rholut_imexop_pack{c,f} differ trivially
// use preprocessor to generate both from the same source template

// suzerain_rholut_imexop_packc
#define FUNCNAME()  suzerain_rholut_imexop_packc
#define PACK(x)     suzerain_bsmbsm_ ## x ## packc
#include "rholut_imexop.def"

// suzerain_rholut_imexop_packf
#define FUNCNAME()  suzerain_rholut_imexop_packf
#define PACK(x)     suzerain_bsmbsm_ ## x ## packf
#include "rholut_imexop.def"
