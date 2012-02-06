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
suzerain_rholut_imexop_apply(
        const double phi,
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
    const double ap43        = s->alpha + 4.0/3.0;
    const double ap13        = s->alpha + 1.0/3.0;
    const double Ma2         = s->Ma * s->Ma;
    const double invRe       = 1 / s->Re;
    const double invMa2      = 1 / Ma2;
    const double ginvRePr    = s->gamma / (s->Re * s->Pr);

    // Accumulate the requested portions of the M + \varphi L operator.
    // The zscal beta/phi and M gbmv phi coefficients scale output by beta.
    // Notice that the phi in (M + phi*L) is achieved using final M gbmv.

    // Shorthand for the common pattern of providing a r->foo, ld->foo pair.
#   define REF(quantity) r->quantity, ld->quantity

    if (in_rho ) {  // Accumulate density terms into out_rho
        suzerain_blas_zscal(n, beta/phi, out_rho, inc);

        if (in_rhou) suzerain_bsplineop_accumulate_complex(M, nrhs,
                -ikm, in_rhou, inc, n, 1.0, out_rho, inc, n, w);

        if (in_rhov) suzerain_bsplineop_accumulate_complex(D1,nrhs,
                -1.0, in_rhov, inc, n, 1.0, out_rho, inc, n, w);

        if (in_rhow) suzerain_bsplineop_accumulate_complex(M, nrhs,
                -ikn, in_rhow, inc, n, 1.0, out_rho, inc, n, w);

        if (in_rhoe) {/* NOP */};

        suzerain_bsplineop_accumulate_complex(
                M, nrhs, 1.0, in_rho,  inc, n, phi, out_rho, inc, n, w);
    }

    if (in_rhou) {  // Accumulate X momentum terms into out_rhou
        suzerain_blas_zscal(n, beta/phi, out_rhou, inc);

        if (in_rho) {
            suzerain_blasext_zgbdddmv_d('N', n, w->kl[M], w->ku[M],
                -0.5*gm1*ikm,          REF(m_gradrho),
                invRe*(ap43*km2+kn2),  REF(nuux),
                -ap13*invRe*(ikm+ikn), REF(nuuz),
                w->D[M],  w->ld, in_rho, inc, 1.0, out_rhou, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D1], w->ku[D1],
                -ap13*invRe*ikm,       REF(nuuy),
                w->D[D1], w->ld, in_rho, inc, 1.0, out_rhou, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D2], w->ku[D2],
                -invRe,                REF(nuux),
                w->D[D2], w->ld, in_rho, inc, 1.0, out_rhou, inc);
        }

        /* in_rhou */ {
            suzerain_blasext_zgbddmv_d('N', n, w->kl[M], w->ku[M],
                -gm1*ikm,              REF(ux),
                -invRe*(ap43*km2+kn2), REF(nu),
                w->D[M],  w->ld, in_rhou, inc, 1.0, out_rhou, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D2], w->ku[D2],
                invRe,                 REF(nu),
                w->D[D2], w->ld, in_rhou, inc, 1.0, out_rhou, inc);
        }

        if (in_rhov) {
            suzerain_blasext_zgbdmv_d('N', n, w->kl[M], w->ku[M],
                -gm1*ikm,       REF(uy),
                w->D[M],  w->ld, in_rhov, inc, 1.0, out_rhou, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D1], w->ku[D1],
                ap13*invRe*ikm, REF(nu),
                w->D[D1], w->ld, in_rhov, inc, 1.0, out_rhou, inc);
        }

        if (in_rhow) {
            suzerain_blasext_zgbddmv_d('N', n, w->kl[M], w->ku[M],
                -gm1*ikm,             REF(uz),
                ap13*invRe*(ikm+ikn), REF(nu),
                w->D[M],  w->ld, in_rhow, inc, 1.0, out_rhou, inc);
        }

        if (in_rhoe) suzerain_bsplineop_accumulate_complex(M, nrhs,
                -gm1*invMa2*ikm, in_rhoe, inc, n, 1.0, out_rhou, inc, n, w);

        suzerain_bsplineop_accumulate_complex(
                M, nrhs, 1.0, in_rhou, inc, n, phi, out_rhou, inc, n, w);
    }

    if (in_rhov) {  // Accumulate Y momentum terms into out_rhov
        suzerain_blas_zscal(n, beta/phi, out_rhov, inc);

        if (in_rho) {
            suzerain_blasext_zgbdmv_d('N', n, w->kl[M], w->ku[M],
                invRe*(km2+kn2), REF(nuuy),
                w->D[M],  w->ld, in_rho, inc, 1.0, out_rhov, inc);

            suzerain_blasext_zgbdddmv_d('N', n, w->kl[D1], w->ku[D1],
                -0.5*gm1,        REF(m_gradrho),
                -ap13*invRe*ikm, REF(nuux),
                -ap13*invRe*ikn, REF(nuuz),
                w->D[D1], w->ld, in_rho, inc, 1.0, out_rhov, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D2], w->ku[D2],
                -ap43*invRe,     REF(nuuy),
                w->D[D2], w->ld, in_rho, inc, 1.0, out_rhov, inc);
        }

        if (in_rhou) {
            suzerain_blasext_zgbddmv_d('N', n, w->kl[D1], w->ku[D1],
                gm1,            REF(ux),
                ap13*invRe*ikm, REF(nu),
                w->D[D1],  w->ld, in_rhou, inc, 1.0, out_rhov, inc);
        }

        /* in_rhov */ {
            suzerain_blasext_zgbdmv_d('N', n, w->kl[M], w->ku[M],
                -invRe*(km2+kn2), REF(nu),
                w->D[M],  w->ld, in_rhov, inc, 1.0, out_rhov, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D1], w->ku[D1],
                gm1,             REF(uy),
                w->D[D1], w->ld, in_rhov, inc, 1.0, out_rhov, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D2], w->ku[D2],
                ap43*invRe,      REF(nu),
                w->D[D2], w->ld, in_rhov, inc, 1.0, out_rhov, inc);
        }

        if (in_rhow) {
            suzerain_blasext_zgbddmv_d('N', n, w->kl[D1], w->ku[D1],
                gm1,            REF(uz),
                ap13*invRe*ikn, REF(nu),
                w->D[D1],  w->ld, in_rhow, inc, 1.0, out_rhov, inc);
        }

        if (in_rhoe) suzerain_bsplineop_accumulate_complex(D1, nrhs,
                -gm1*invMa2, in_rhoe, inc, n, 1.0, out_rhov, inc, n, w);

        suzerain_bsplineop_accumulate_complex(
                M, nrhs, 1.0, in_rhov, inc, n, phi, out_rhov, inc, n, w);
    }

    if (in_rhow) {  // Accumulate Z momentum terms into out_rhow
        suzerain_blas_zscal(n, beta/phi, out_rhow, inc);

        if (in_rho) {
            suzerain_blasext_zgbdddmv_d('N', n, w->kl[M], w->ku[M],
                -0.5*gm1*ikn,          REF(m_gradrho),
                invRe*(km2+ap43*kn2),  REF(nuuz),
                -ap13*invRe*(ikm+ikn), REF(nuux),
                w->D[M],  w->ld, in_rho, inc, 1.0, out_rhow, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D1], w->ku[D1],
                -ap13*invRe*ikn,       REF(nuuy),
                w->D[D1], w->ld, in_rho, inc, 1.0, out_rhow, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D2], w->ku[D2],
                -invRe,                REF(nuuz),
                w->D[D2], w->ld, in_rho, inc, 1.0, out_rhow, inc);
        }

        if (in_rhou) {
            suzerain_blasext_zgbddmv_d('N', n, w->kl[M], w->ku[M],
                -gm1*ikn,             REF(ux),
                ap13*invRe*(ikm+ikn), REF(nu),
                w->D[M],  w->ld, in_rhou, inc, 1.0, out_rhow, inc);
        }

        if (in_rhov) {
            suzerain_blasext_zgbdmv_d('N', n, w->kl[M], w->ku[M],
                -gm1*ikn,       REF(uy),
                w->D[M],  w->ld, in_rhov, inc, 1.0, out_rhow, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D1], w->ku[D1],
                ap13*invRe*ikn, REF(nu),
                w->D[D1], w->ld, in_rhov, inc, 1.0, out_rhow, inc);
        }

        /* in_rhow */ {
            suzerain_blasext_zgbddmv_d('N', n, w->kl[M], w->ku[M],
                -gm1*ikn,              REF(uz),
                -invRe*(km2+ap43*kn2), REF(nu),
                w->D[M],  w->ld, in_rhow, inc, 1.0, out_rhow, inc);

            suzerain_blasext_zgbdmv_d('N', n, w->kl[D2], w->ku[D2],
                invRe,                 REF(nu),
                w->D[D2], w->ld, in_rhow, inc, 1.0, out_rhow, inc);
        }

        if (in_rhoe) suzerain_bsplineop_accumulate_complex(M, nrhs,
                -gm1*invMa2*ikn, in_rhoe, inc, n, 1.0, out_rhow, inc, n, w);

        suzerain_bsplineop_accumulate_complex(
                M, nrhs, 1.0, in_rhow, inc, n, phi, out_rhow, inc, n, w);
    }

    if (in_rhoe) {  // Accumulate total energy terms into out_rhoe
        suzerain_blas_zscal(n, beta/phi, out_rhoe, inc);

        suzerain_bsplineop_accumulate_complex(
                M, nrhs, 1.0, in_rhoe, inc, n, phi, out_rhoe, inc, n, w);
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
