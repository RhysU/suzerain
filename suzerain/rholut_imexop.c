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

    // Prepare several oft-used constants to aid readability
    static const int inc        = 1;
    static const int nrhs       = 1;
    static const int M          = 0;
    static const int D1         = 1;
    static const int D2         = 2;
    const int    n              = w->n;

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
    const double ginvRePr    = s->gamma / (s->Re * s->Pr);

    // Accumulate the requested portions of the M + \varphi L operator.  Scale
    // output by beta, accumulate non-mass contributions, and finally
    // accumulate the mass contributions.  Mass contributions come last as they
    // are expected to have magnitudes much larger than phi-- this helps to
    // ensure better rounding of \varphi L contributions.

    // Shorthand for the common pattern of providing a r->foo, ld->foo pair.
#   define REF(quantity) r->quantity, ld->quantity

    if (in_rho ) {  // Accumulate density terms into out_rho

        suzerain_blas_zscal(n, beta, out_rho, inc);

        if (in_rhou) suzerain_bsplineop_accumulate_complex(M, nrhs,
                -phi*ikm, in_rhou, inc, n, 1.0, out_rho, inc, n, w);

        if (in_rhov) suzerain_bsplineop_accumulate_complex(D1,nrhs,
                -phi*1.0, in_rhov, inc, n, 1.0, out_rho, inc, n, w);

        if (in_rhow) suzerain_bsplineop_accumulate_complex(M, nrhs,
                -phi*ikn, in_rhow, inc, n, 1.0, out_rho, inc, n, w);

        if (in_rhoe) {/* NOP */};

        suzerain_bsplineop_accumulate_complex(
            M, nrhs, 1.0, in_rho, inc, n, 1.0, out_rho, inc, n, w);
    }

    if (in_rhou) {  // Accumulate X momentum terms into out_rhou

        suzerain_blas_zscal(n, beta, out_rhou, inc);

        if (in_rho) {
            suzerain_blasext_zgbdddddmv_d('N', n, w->kl[M], w->ku[M],
                phi*ikm,                   REF(uxux),
                phi*ikn,                   REF(uxuz),
                phi*invRe*(ap43*km2+kn2),  REF(nuux),
                phi*ap13*invRe*km*kn,      REF(nuuz),
                -phi*0.5*gm1*ikm,          REF(m_gradrho),
                w->D[M],  w->ld, in_rho, inc, 1.0, out_rhou, inc);

            suzerain_blasext_zgbddmv_d('N', n, w->kl[D1], w->ku[D1],
                phi,                       REF(uxuy),
                -phi*ap13*invRe*ikm,       REF(nuuy),
                w->D[D1], w->ld, in_rho, inc, 1.0, out_rhou, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D2], w->ku[D2],
                -phi*invRe,                REF(nuux),
                w->D[D2], w->ld, in_rho, inc, 1.0, out_rhou, inc);
        }

        /* in_rhou */ {
            suzerain_blasext_zgbdddmv_d('N', n, w->kl[M], w->ku[M],
                phi*gm3*ikm,               REF(ux),
                -phi*ikn,                  REF(uz),
                -phi*invRe*(ap43*km2+kn2), REF(nu),
                w->D[M],  w->ld, in_rhou, inc, 1.0, out_rhou, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D1], w->ku[D1],
                -phi,                      REF(uy),
                w->D[D1], w->ld, in_rhou, inc, 1.0, out_rhou, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D2], w->ku[D2],
                phi*invRe,                 REF(nu),
                w->D[D2], w->ld, in_rhou, inc, 1.0, out_rhou, inc);
        }

        if (in_rhov) {
            suzerain_blasext_zgbdmv_d('N', n, w->kl[M], w->ku[M],
                phi*gm1*ikm,               REF(uy),
                w->D[M],  w->ld, in_rhov, inc, 1.0, out_rhou, inc);

            suzerain_blasext_zgbddmv_d('N', n, w->kl[D1], w->ku[D1],
                -phi,                      REF(ux),
                phi*ap13*invRe*ikm,        REF(nu),
                w->D[D1], w->ld, in_rhov, inc, 1.0, out_rhou, inc);
        }

        if (in_rhow) {
            suzerain_blasext_zgbdddmv_d('N', n, w->kl[M], w->ku[M],
                -phi*ikn,                  REF(ux),
                 phi*gm1*ikm,              REF(uz),
                -phi*ap13*invRe*km*kn,     REF(nu),
                w->D[M],  w->ld, in_rhow, inc, 1.0, out_rhou, inc);
        }

        if (in_rhoe) suzerain_bsplineop_accumulate_complex(M, nrhs,
                -phi*gm1*invMa2*ikm, in_rhoe, inc, n,
                1.0, out_rhou, inc, n, w);

        suzerain_bsplineop_accumulate_complex(
            M, nrhs, 1.0, in_rhou, inc, n, 1.0, out_rhou, inc, n, w);
    }

    if (in_rhov) {  // Accumulate Y momentum terms into out_rhov

        suzerain_blas_zscal(n, beta, out_rhov, inc);

        if (in_rho) {
            suzerain_blasext_zgbdddmv_d('N', n, w->kl[M], w->ku[M],
                phi*ikm,              REF(uxuy),
                phi*ikn,              REF(uyuz),
                phi*invRe*(km2+kn2),  REF(nuuy),
                w->D[M],  w->ld, in_rho, inc, 1.0, out_rhov, inc);

            suzerain_blasext_zgbddddmv_d('N', n, w->kl[D1], w->ku[D1],
                 phi,                 REF(uyuy),
                -phi*ap13*invRe*ikm,  REF(nuux),
                -phi*ap13*invRe*ikn,  REF(nuuz),
                -phi*0.5*gm1,         REF(m_gradrho),
                w->D[D1], w->ld, in_rho, inc, 1.0, out_rhov, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D2], w->ku[D2],
                -phi*ap43*invRe,      REF(nuuy),
                w->D[D2], w->ld, in_rho, inc, 1.0, out_rhov, inc);
        }

        if (in_rhou) {
            suzerain_blasext_zgbdmv_d('N', n, w->kl[M], w->ku[M],
                -phi*ikm,             REF(uy),
                w->D[M],   w->ld, in_rhou, inc, 1.0, out_rhov, inc);

            suzerain_blasext_zgbddmv_d('N', n, w->kl[D1], w->ku[D1],
                phi*gm1,              REF(ux),
                phi*ap13*invRe*ikm,   REF(nu),
                w->D[D1],  w->ld, in_rhou, inc, 1.0, out_rhov, inc);
        }

        /* in_rhov */ {
            suzerain_blasext_zgbdddmv_d('N', n, w->kl[M], w->ku[M],
                -phi*ikm,             REF(ux),
                -phi*ikn,             REF(uz),
                -phi*invRe*(km2+kn2), REF(nu),
                w->D[M],  w->ld, in_rhov, inc, 1.0, out_rhov, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D1], w->ku[D1],
                phi*gm3,              REF(uy),
                w->D[D1], w->ld, in_rhov, inc, 1.0, out_rhov, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D2], w->ku[D2],
                phi*ap43*invRe,       REF(nu),
                w->D[D2], w->ld, in_rhov, inc, 1.0, out_rhov, inc);
        }

        if (in_rhow) {
            suzerain_blasext_zgbdmv_d('N', n, w->kl[M], w->ku[M],
                -phi*ikn,             REF(uy),
                w->D[M],   w->ld, in_rhow, inc, 1.0, out_rhov, inc);

            suzerain_blasext_zgbddmv_d('N', n, w->kl[D1], w->ku[D1],
                phi*gm1,              REF(uz),
                phi*ap13*invRe*ikn,   REF(nu),
                w->D[D1],  w->ld, in_rhow, inc, 1.0, out_rhov, inc);
        }

        if (in_rhoe) suzerain_bsplineop_accumulate_complex(D1, nrhs,
                -phi*gm1*invMa2, in_rhoe, inc, n,
                1.0, out_rhov, inc, n, w);

        suzerain_bsplineop_accumulate_complex(
            M, nrhs, 1.0, in_rhov, inc, n, 1.0, out_rhov, inc, n, w);
    }

    if (in_rhow) {  // Accumulate Z momentum terms into out_rhow

        suzerain_blas_zscal(n, beta, out_rhow, inc);

        if (in_rho) {
            suzerain_blasext_zgbdddddmv_d('N', n, w->kl[M], w->ku[M],
                phi*ikm,                   REF(uxuz),
                phi*ikn,                   REF(uzuz),
                phi*invRe*(km2+ap43*kn2),  REF(nuuz),
                phi*ap13*invRe*km*kn,      REF(nuux),
                -phi*0.5*gm1*ikn,          REF(m_gradrho),
                w->D[M],  w->ld, in_rho, inc, 1.0, out_rhow, inc);

            suzerain_blasext_zgbddmv_d('N', n, w->kl[D1], w->ku[D1],
                 phi,                      REF(uyuz),
                -phi*ap13*invRe*ikn,       REF(nuuy),
                w->D[D1], w->ld, in_rho, inc, 1.0, out_rhow, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D2], w->ku[D2],
                -phi*invRe,                REF(nuuz),
                w->D[D2], w->ld, in_rho, inc, 1.0, out_rhow, inc);
        }

        if (in_rhou) {
            suzerain_blasext_zgbdddmv_d('N', n, w->kl[M], w->ku[M],
                 phi*gm1*ikn,              REF(ux),
                -phi*ikm,                  REF(uz),
                -phi*ap13*invRe*km*kn,     REF(nu),
                w->D[M],  w->ld, in_rhou, inc, 1.0, out_rhow, inc);
        }

        if (in_rhov) {
            suzerain_blasext_zgbdmv_d('N', n, w->kl[M], w->ku[M],
                phi*gm1*ikn,               REF(uy),
                w->D[M],  w->ld, in_rhov, inc, 1.0, out_rhow, inc);

            suzerain_blasext_zgbddmv_d('N', n, w->kl[D1], w->ku[D1],
                -phi,                      REF(uz),
                phi*ap13*invRe*ikn,        REF(nu),
                w->D[D1], w->ld, in_rhov, inc, 1.0, out_rhow, inc);
        }

        /* in_rhow */ {
            suzerain_blasext_zgbdddmv_d('N', n, w->kl[M], w->ku[M],
                -phi*ikm,                  REF(ux),
                phi*gm3*ikn,               REF(uz),
                -phi*invRe*(km2+ap43*kn2), REF(nu),
                w->D[M],  w->ld, in_rhow, inc, 1.0, out_rhow, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D1], w->ku[D1],
                -phi,                      REF(uy),
                w->D[D1], w->ld, in_rhow, inc, 1.0, out_rhow, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D2], w->ku[D2],
                phi*invRe,                 REF(nu),
                w->D[D2], w->ld, in_rhow, inc, 1.0, out_rhow, inc);
        }

        if (in_rhoe) suzerain_bsplineop_accumulate_complex(M, nrhs,
                -phi*gm1*invMa2*ikn, in_rhoe, inc, n,
                1.0, out_rhow, inc, n, w);

        suzerain_bsplineop_accumulate_complex(
            M, nrhs, 1.0, in_rhow, inc, n, 1.0, out_rhow, inc, n, w);
    }

    if (in_rhoe) {  // Accumulate total energy terms into out_rhoe

        suzerain_blas_zscal(n, beta, out_rhoe, inc);

        if (in_rho) {
            suzerain_blasext_zgbdddmv_d('N', n, w->kl[M], w->ku[M],
                -phi*ikm,                    REF(ex_gradrho),
                -phi*ikn,                    REF(ez_gradrho),
                -phi*ginvRePr/gm1*(km2+kn2), REF(e_deltarho),
                w->D[M],  w->ld, in_rho, inc, 1.0, out_rhoe, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D1], w->ku[D1],
                -phi,                        REF(ey_gradrho),
                w->D[D1], w->ld, in_rho, inc, 1.0, out_rhoe, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D2], w->ku[D2],
                phi*ginvRePr/gm1,            REF(e_deltarho),
                w->D[D2], w->ld, in_rho, inc, 1.0, out_rhoe, inc);
        }

        if (in_rhou) {
            suzerain_blasext_zgbddmv_d('N', n, w->kl[M], w->ku[M],
                phi*ginvRePr*Ma2*(km2+kn2),  REF(nuux),
                -phi*ikm,                    REF(e_divm),
                w->D[M],  w->ld, in_rhou, inc, 1.0, out_rhoe, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D2], w->ku[D2],
                -phi*ginvRePr*Ma2,           REF(nuux),
                w->D[D2], w->ld, in_rhou, inc, 1.0, out_rhoe, inc);
        }

        if (in_rhov) {
            suzerain_blasext_zgbdmv_d('N', n, w->kl[M], w->ku[M],
                phi*ginvRePr*Ma2*(km2+kn2),  REF(nuuy),
                w->D[M],  w->ld, in_rhov, inc, 1.0, out_rhoe, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D1], w->ku[D1],
                -phi,                        REF(e_divm),
                w->D[D1], w->ld, in_rhov, inc, 1.0, out_rhoe, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D2], w->ku[D2],
                -phi*ginvRePr*Ma2,           REF(nuuy),
                w->D[D2], w->ld, in_rhov, inc, 1.0, out_rhoe, inc);
        }

        if (in_rhow) {
            suzerain_blasext_zgbddmv_d('N', n, w->kl[M], w->ku[M],
                phi*ginvRePr*Ma2*(km2+kn2),  REF(nuuz),
                -phi*ikn,                    REF(e_divm),
                w->D[M],  w->ld, in_rhow, inc, 1.0, out_rhoe, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D2], w->ku[D2],
                -phi*ginvRePr*Ma2,           REF(nuuz),
                w->D[D2], w->ld, in_rhow, inc, 1.0, out_rhoe, inc);
        }

        /* in_rhoe */ {
            suzerain_blasext_zgbdddmv_d('N', n, w->kl[M], w->ku[M],
                -phi*s->gamma*ikm,           REF(ux),
                -phi*s->gamma*ikn,           REF(uz),
                -phi*ginvRePr*(km2+kn2),     REF(nu),
                w->D[M],  w->ld, in_rhoe, inc, 1.0, out_rhoe, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D1], w->ku[D1],
                -phi*s->gamma,               REF(uy),
                w->D[D1], w->ld, in_rhoe, inc, 1.0, out_rhoe, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D2], w->ku[D2],
                phi*ginvRePr,                REF(nu),
                w->D[D2], w->ld, in_rhoe, inc, 1.0, out_rhoe, inc);
        }

        suzerain_bsplineop_accumulate_complex(
            M, nrhs, 1.0, in_rhoe, inc, n, 1.0, out_rhoe, inc, n, w);
    }

#   undef REF

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
