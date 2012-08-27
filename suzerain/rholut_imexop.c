/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2012 The PECOS Development Team
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
        const int imagzero,
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

    // Two use cases are accommodated.  The first (and vastly more common) case
    // is applying the operator to in_rho{,u,v,w,e} as provided.  The second
    // case applies the operator for Im(in_rho{,u,v,w,e}) == 0.  Imaginary
    // parts are entirely ignored.  This is useful to handle the 0th and
    // Nyquist modes in an even-length DFT without resorting to auxiliary
    // buffers, different memory access patterns, or duplicated code.
    //
    // For brevity, believe it or not, typedef some gbmv-like function
    // signatures operations but without a (const void *) x vector type:
    typedef int gbmv_t(const char, const int, const int, const int, const int,
                            const complex_double,   const double *, const int,
                                             /* NB */ const void *, const int,
                            const complex_double, complex_double *, const int);
    typedef int gbdmv_t    (const char, const int, const int, const int,
                            const complex_double,   const double *, const int,
                                                    const double *, const int,
                                             /* NB */ const void *, const int,
                            const complex_double, complex_double *, const int);
    typedef int gbddmv_t   (const char, const int, const int, const int,
                            const complex_double,   const double *, const int,
                            const complex_double,   const double *, const int,
                                                    const double *, const int,
                                             /* NB */ const void *, const int,
                            const complex_double, complex_double *, const int);
    typedef int gbdddmv_t  (const char, const int, const int, const int,
                            const complex_double,   const double *, const int,
                            const complex_double,   const double *, const int,
                            const complex_double,   const double *, const int,
                                                    const double *, const int,
                                             /* NB */ const void *, const int,
                            const complex_double, complex_double *, const int);
    typedef int gbddddmv_t (const char, const int, const int, const int,
                            const complex_double,   const double *, const int,
                            const complex_double,   const double *, const int,
                            const complex_double,   const double *, const int,
                            const complex_double,   const double *, const int,
                                                    const double *, const int,
                                             /* NB */ const void *, const int,
                            const complex_double, complex_double *, const int);
    typedef int gbdddddmv_t(const char, const int, const int, const int,
                            const complex_double,   const double *, const int,
                            const complex_double,   const double *, const int,
                            const complex_double,   const double *, const int,
                            const complex_double,   const double *, const int,
                            const complex_double,   const double *, const int,
                                                    const double *, const int,
                                             /* NB */ const void *, const int,
                            const complex_double, complex_double *, const int);

    // In the first case, treat x arguments as complex with stride one...
    int inc_in = 1;
    gbmv_t      *p_gbmv      = (gbmv_t *)      &        suzerain_blas_zgbmv_d_z;
    gbdmv_t     *p_gbdmv     = (gbdmv_t *)     &    suzerain_blasext_zgbdmv_d_z;
    gbddmv_t    *p_gbddmv    = (gbddmv_t *)    &   suzerain_blasext_zgbddmv_d_z;
    gbdddmv_t   *p_gbdddmv   = (gbdddmv_t *)   &  suzerain_blasext_zgbdddmv_d_z;
    gbddddmv_t  *p_gbddddmv  = (gbddddmv_t *)  & suzerain_blasext_zgbddddmv_d_z;
    gbdddddmv_t *p_gbdddddmv = (gbdddddmv_t *) &suzerain_blasext_zgbdddddmv_d_z;

    // ...but in the second case, treat x arguments as real with stride two.
    if (SUZERAIN_UNLIKELY(imagzero)) {
        inc_in = 2;
        p_gbmv      = (gbmv_t *)      &        suzerain_blas_zgbmv_d_d;
        p_gbdmv     = (gbdmv_t *)     &    suzerain_blasext_zgbdmv_d_d;
        p_gbddmv    = (gbddmv_t *)    &   suzerain_blasext_zgbddmv_d_d;
        p_gbdddmv   = (gbdddmv_t *)   &  suzerain_blasext_zgbdddmv_d_d;
        p_gbddddmv  = (gbddddmv_t *)  & suzerain_blasext_zgbddddmv_d_d;
        p_gbdddddmv = (gbdddddmv_t *) &suzerain_blasext_zgbdddddmv_d_d;
    }

    // Shorthand for the common pattern of providing "in_foo, inc" pairs.
#   define IN(quantity)  in_##quantity, inc_in
#   define OUT(quantity) out_##quantity, 1

    // Shorthand for the common pattern of providing a "r->foo, ld->foo" pair.
#   define REF(quantity) r->quantity, ld->quantity

    // Just plain shorthand
#   define LIKELY(expr) SUZERAIN_LIKELY(expr)

    // We need to account for suzerain_bsplineop_workspace storing the
    // transpose of the operators when we invoke suzerain_blaseext_* routines.
    static const char trans = 'T';

    // Prepare several oft-used constants to aid readability
    static const int M       = 0;
    static const int D1      = 1;
    static const int D2      = 2;
    const int        n       = w->n;

    if (LIKELY(in_rho )) {  // Accumulate density terms into out_rho

        suzerain_blas_zscal(n, beta, OUT(rho));

        if (LIKELY(in_rhou)) (*p_gbmv)(trans, n, n, w->kl[M],  w->ku[M],
                -phi*ikm, w->D_T[M], w->ld, IN(rhou),  1.0, OUT(rho));

        if (LIKELY(in_rhov)) (*p_gbmv)(trans, n, n, w->kl[D1], w->ku[D1],
                -phi,     w->D_T[D1], w->ld, IN(rhov), 1.0, OUT(rho));

        if (LIKELY(in_rhow)) (*p_gbmv)(trans, n, n, w->kl[M],  w->ku[M],
                -phi*ikn, w->D_T[M], w->ld, IN(rhow),  1.0, OUT(rho));

        if (LIKELY(in_rhoe)) {/* NOP */};

        (*p_gbmv)(trans, n, n, w->kl[M], w->ku[M],
            1.0, w->D_T[M], w->ld, IN(rho), 1.0, OUT(rho));
    }

    if (LIKELY(in_rhou)) {  // Accumulate X momentum terms into out_rhou

        suzerain_blas_zscal(n, beta, OUT(rhou));

        if (LIKELY(in_rho)) {
            (*p_gbdddddmv)(trans, n, w->kl[M], w->ku[M],
                -phi*0.5*gm1*ikm,          REF(u2),
                phi*ikm,                   REF(uxux),
                phi*ikn,                   REF(uxuz),
                phi*invRe*(ap43*km2+kn2),  REF(nuux),
                phi*ap13*invRe*km*kn,      REF(nuuz),
                w->D_T[M],  w->ld, IN(rho), 1.0, OUT(rhou));

            (*p_gbddmv)(trans, n, w->kl[D1], w->ku[D1],
                phi,                       REF(uxuy),
                -phi*ap13*invRe*ikm,       REF(nuuy),
                w->D_T[D1], w->ld, IN(rho), 1.0, OUT(rhou));

            (*p_gbdmv)(trans, n, w->kl[D2], w->ku[D2],
                -phi*invRe,                REF(nuux),
                w->D_T[D2], w->ld, IN(rho), 1.0, OUT(rhou));
        }

        /* in_rhou */ {
            (*p_gbdddmv)(trans, n, w->kl[M], w->ku[M],
                phi*gm3*ikm,               REF(ux),
                -phi*ikn,                  REF(uz),
                -phi*invRe*(ap43*km2+kn2), REF(nu),
                w->D_T[M],  w->ld, IN(rhou), 1.0, OUT(rhou));

            (*p_gbdmv)(trans, n, w->kl[D1], w->ku[D1],
                -phi,                      REF(uy),
                w->D_T[D1], w->ld, IN(rhou), 1.0, OUT(rhou));

            (*p_gbdmv)(trans, n, w->kl[D2], w->ku[D2],
                phi*invRe,                 REF(nu),
                w->D_T[D2], w->ld, IN(rhou), 1.0, OUT(rhou));
        }

        if (LIKELY(in_rhov)) {
            (*p_gbdmv)(trans, n, w->kl[M], w->ku[M],
                phi*gm1*ikm,               REF(uy),
                w->D_T[M],  w->ld, IN(rhov), 1.0, OUT(rhou));

            (*p_gbddmv)(trans, n, w->kl[D1], w->ku[D1],
                -phi,                      REF(ux),
                phi*ap13*invRe*ikm,        REF(nu),
                w->D_T[D1], w->ld, IN(rhov), 1.0, OUT(rhou));
        }

        if (LIKELY(in_rhow)) {
            (*p_gbdddmv)(trans, n, w->kl[M], w->ku[M],
                -phi*ikn,                  REF(ux),
                 phi*gm1*ikm,              REF(uz),
                -phi*ap13*invRe*km*kn,     REF(nu),
                w->D_T[M],  w->ld, IN(rhow), 1.0, OUT(rhou));
        }

        if (LIKELY(in_rhoe)) (*p_gbmv)(trans, n, n, w->kl[M], w->ku[M],
                -phi*gm1*invMa2*ikm, w->D_T[M], w->ld, IN(rhoe),
                1.0, OUT(rhou));

        (*p_gbmv)(trans, n, n, w->kl[M], w->ku[M],
            1.0, w->D_T[M], w->ld, IN(rhou), 1.0, OUT(rhou));
    }

    if (LIKELY(in_rhov)) {  // Accumulate Y momentum terms into out_rhov

        suzerain_blas_zscal(n, beta, OUT(rhov));

        if (LIKELY(in_rho)) {
            (*p_gbdddmv)(trans, n, w->kl[M], w->ku[M],
                phi*ikm,              REF(uxuy),
                phi*ikn,              REF(uyuz),
                phi*invRe*(km2+kn2),  REF(nuuy),
                w->D_T[M],  w->ld, IN(rho), 1.0, OUT(rhov));

            (*p_gbddddmv)(trans, n, w->kl[D1], w->ku[D1],
                -phi*0.5*gm1,         REF(u2),
                 phi,                 REF(uyuy),
                -phi*ap13*invRe*ikm,  REF(nuux),
                -phi*ap13*invRe*ikn,  REF(nuuz),
                w->D_T[D1], w->ld, IN(rho), 1.0, OUT(rhov));

            (*p_gbdmv)(trans, n, w->kl[D2], w->ku[D2],
                -phi*ap43*invRe,      REF(nuuy),
                w->D_T[D2], w->ld, IN(rho), 1.0, OUT(rhov));
        }

        if (LIKELY(in_rhou)) {
            (*p_gbdmv)(trans, n, w->kl[M], w->ku[M],
                -phi*ikm,             REF(uy),
                w->D_T[M],   w->ld, IN(rhou), 1.0, OUT(rhov));

            (*p_gbddmv)(trans, n, w->kl[D1], w->ku[D1],
                phi*gm1,              REF(ux),
                phi*ap13*invRe*ikm,   REF(nu),
                w->D_T[D1],  w->ld, IN(rhou), 1.0, OUT(rhov));
        }

        /* in_rhov */ {
            (*p_gbdddmv)(trans, n, w->kl[M], w->ku[M],
                -phi*ikm,             REF(ux),
                -phi*ikn,             REF(uz),
                -phi*invRe*(km2+kn2), REF(nu),
                w->D_T[M],  w->ld, IN(rhov), 1.0, OUT(rhov));

            (*p_gbdmv)(trans, n, w->kl[D1], w->ku[D1],
                phi*gm3,              REF(uy),
                w->D_T[D1], w->ld, IN(rhov), 1.0, OUT(rhov));

            (*p_gbdmv)(trans, n, w->kl[D2], w->ku[D2],
                phi*ap43*invRe,       REF(nu),
                w->D_T[D2], w->ld, IN(rhov), 1.0, OUT(rhov));
        }

        if (LIKELY(in_rhow)) {
            (*p_gbdmv)(trans, n, w->kl[M], w->ku[M],
                -phi*ikn,             REF(uy),
                w->D_T[M],   w->ld, IN(rhow), 1.0, OUT(rhov));

            (*p_gbddmv)(trans, n, w->kl[D1], w->ku[D1],
                phi*gm1,              REF(uz),
                phi*ap13*invRe*ikn,   REF(nu),
                w->D_T[D1],  w->ld, IN(rhow), 1.0, OUT(rhov));
        }

        if (LIKELY(in_rhoe)) (*p_gbmv)(trans, n, n, w->kl[D1], w->ku[D1],
                -phi*gm1*invMa2, w->D_T[D1], w->ld, IN(rhoe),
                1.0, OUT(rhov));

        (*p_gbmv)(trans, n, n, w->kl[M], w->ku[M],
            1.0, w->D_T[M], w->ld, IN(rhov), 1.0, OUT(rhov));
    }

    if (LIKELY(in_rhow)) {  // Accumulate Z momentum terms into out_rhow

        suzerain_blas_zscal(n, beta, OUT(rhow));

        if (LIKELY(in_rho)) {
            (*p_gbdddddmv)(trans, n, w->kl[M], w->ku[M],
                -phi*0.5*gm1*ikn,          REF(u2),
                phi*ikm,                   REF(uxuz),
                phi*ikn,                   REF(uzuz),
                phi*invRe*(km2+ap43*kn2),  REF(nuuz),
                phi*ap13*invRe*km*kn,      REF(nuux),
                w->D_T[M],  w->ld, IN(rho), 1.0, OUT(rhow));

            (*p_gbddmv)(trans, n, w->kl[D1], w->ku[D1],
                 phi,                      REF(uyuz),
                -phi*ap13*invRe*ikn,       REF(nuuy),
                w->D_T[D1], w->ld, IN(rho), 1.0, OUT(rhow));

            (*p_gbdmv)(trans, n, w->kl[D2], w->ku[D2],
                -phi*invRe,                REF(nuuz),
                w->D_T[D2], w->ld, IN(rho), 1.0, OUT(rhow));
        }

        if (LIKELY(in_rhou)) {
            (*p_gbdddmv)(trans, n, w->kl[M], w->ku[M],
                 phi*gm1*ikn,              REF(ux),
                -phi*ikm,                  REF(uz),
                -phi*ap13*invRe*km*kn,     REF(nu),
                w->D_T[M],  w->ld, IN(rhou), 1.0, OUT(rhow));
        }

        if (LIKELY(in_rhov)) {
            (*p_gbdmv)(trans, n, w->kl[M], w->ku[M],
                phi*gm1*ikn,               REF(uy),
                w->D_T[M],  w->ld, IN(rhov), 1.0, OUT(rhow));

            (*p_gbddmv)(trans, n, w->kl[D1], w->ku[D1],
                -phi,                      REF(uz),
                phi*ap13*invRe*ikn,        REF(nu),
                w->D_T[D1], w->ld, IN(rhov), 1.0, OUT(rhow));
        }

        /* in_rhow */ {
            (*p_gbdddmv)(trans, n, w->kl[M], w->ku[M],
                -phi*ikm,                  REF(ux),
                phi*gm3*ikn,               REF(uz),
                -phi*invRe*(km2+ap43*kn2), REF(nu),
                w->D_T[M],  w->ld, IN(rhow), 1.0, OUT(rhow));

            (*p_gbdmv)(trans, n, w->kl[D1], w->ku[D1],
                -phi,                      REF(uy),
                w->D_T[D1], w->ld, IN(rhow), 1.0, OUT(rhow));

            (*p_gbdmv)(trans, n, w->kl[D2], w->ku[D2],
                phi*invRe,                 REF(nu),
                w->D_T[D2], w->ld, IN(rhow), 1.0, OUT(rhow));
        }

        if (LIKELY(in_rhoe)) (*p_gbmv)(trans, n, n, w->kl[M], w->ku[M],
                -phi*gm1*invMa2*ikn, w->D_T[M], w->ld, IN(rhoe),
                1.0, OUT(rhow));

        (*p_gbmv)(trans, n, n, w->kl[M], w->ku[M],
            1.0, w->D_T[M], w->ld, IN(rhow), 1.0, OUT(rhow));
    }

    if (LIKELY(in_rhoe)) {  // Accumulate total energy terms into out_rhoe

        suzerain_blas_zscal(n, beta, OUT(rhoe));

        if (LIKELY(in_rho)) {

            // Mass terms done in 2 passes to avoid zgbdddddddmv_d.
            // Writing such a beast may provide a tiny speedup.
            (*p_gbddddmv)(trans, n, w->kl[M], w->ku[M],
                phi*Ma2*invRe*(km2+kn2),     REF(nuu2),
                phi*Ma2*invRe*ap13*km2,      REF(nuuxux),
                phi*Ma2*invRe*ap13*2*km*kn,  REF(nuuxuz),
                phi*Ma2*invRe*ap13*kn2,      REF(nuuzuz),
                w->D_T[M],  w->ld, IN(rho), 1.0, OUT(rhoe));
            (*p_gbdddmv)( trans, n, w->kl[M], w->ku[M],
                -phi*ikm,                    REF(ex_gradrho),
                -phi*ikn,                    REF(ez_gradrho),
                -phi*ginvRePr/gm1*(km2+kn2), REF(e_deltarho),
                w->D_T[M],  w->ld, IN(rho), 1.0, OUT(rhoe));

            (*p_gbdddmv)(trans, n, w->kl[D1], w->ku[D1],
                -phi*Ma2*invRe*ap13*2*ikm,   REF(nuuxuy),
                -phi*Ma2*invRe*ap13*2*ikn,   REF(nuuyuz),
                -phi,                        REF(ey_gradrho),
                w->D_T[D1], w->ld, IN(rho), 1.0, OUT(rhoe));

            (*p_gbdddmv)(trans, n, w->kl[D2], w->ku[D2],
                -phi*Ma2*invRe,              REF(nuu2),
                -phi*Ma2*invRe*ap13,         REF(nuuyuy),
                phi*ginvRePr/gm1,            REF(e_deltarho),
                w->D_T[D2], w->ld, IN(rho), 1.0, OUT(rhoe));
        }

        if (LIKELY(in_rhou)) {
            const double coeff_nuux
                = Ma2*invRe*((ginvPr-ap43)*km2 + (ginvPr-1)*kn2);
            (*p_gbdddmv)(trans, n, w->kl[M], w->ku[M],
                phi*coeff_nuux,              REF(nuux),
                -phi*Ma2*invRe*ap13*km*kn,   REF(nuuz),
                -phi*ikm,                    REF(e_divm),
                w->D_T[M],  w->ld, IN(rhou), 1.0, OUT(rhoe));

            (*p_gbdmv)(trans, n, w->kl[D1], w->ku[D1],
                phi*Ma2*invRe*ap13*ikm,      REF(nuuy),
                w->D_T[D1], w->ld, IN(rhou), 1.0, OUT(rhoe));

            (*p_gbdmv)(trans, n, w->kl[D2], w->ku[D2],
                phi*Ma2*invRe*(1-ginvPr),    REF(nuux),
                w->D_T[D2], w->ld, IN(rhou), 1.0, OUT(rhoe));
        }

        if (LIKELY(in_rhov)) {
            const double coeff_nuuy
                = Ma2*invRe*(ginvPr-1)*(km2+kn2);
            (*p_gbdmv)(trans, n, w->kl[M], w->ku[M],
                phi*coeff_nuuy,              REF(nuuy),
                w->D_T[M],  w->ld, IN(rhov), 1.0, OUT(rhoe));

            (*p_gbdddmv)(trans, n, w->kl[D1], w->ku[D1],
                phi*Ma2*invRe*ap13*ikm,      REF(nuux),
                phi*Ma2*invRe*ap13*ikn,      REF(nuuz),
                -phi,                        REF(e_divm),
                w->D_T[D1], w->ld, IN(rhov), 1.0, OUT(rhoe));

            (*p_gbdmv)(trans, n, w->kl[D2], w->ku[D2],
                phi*Ma2*invRe*(ap43-ginvPr), REF(nuuy),
                w->D_T[D2], w->ld, IN(rhov), 1.0, OUT(rhoe));
        }

        if (LIKELY(in_rhow)) {
            const double coeff_nuuz
                = Ma2*invRe*((ginvPr-1)*km2 + (ginvPr-ap43)*kn2);
            (*p_gbdddmv)(trans, n, w->kl[M], w->ku[M],
                -phi*Ma2*invRe*ap13*km*kn,   REF(nuux),
                phi*coeff_nuuz,              REF(nuuz),
                -phi*ikn,                    REF(e_divm),
                w->D_T[M],  w->ld, IN(rhow), 1.0, OUT(rhoe));

            (*p_gbdmv)(trans, n, w->kl[D1], w->ku[D1],
                phi*Ma2*invRe*ap13*ikn,      REF(nuuy),
                w->D_T[D1], w->ld, IN(rhow), 1.0, OUT(rhoe));

            (*p_gbdmv)(trans, n, w->kl[D2], w->ku[D2],
                phi*Ma2*invRe*(1-ginvPr),    REF(nuuz),
                w->D_T[D2], w->ld, IN(rhow), 1.0, OUT(rhoe));
        }

        /* in_rhoe */ {
            (*p_gbdddmv)(trans, n, w->kl[M], w->ku[M],
                -phi*s->gamma*ikm,           REF(ux),
                -phi*s->gamma*ikn,           REF(uz),
                -phi*ginvRePr*(km2+kn2),     REF(nu),
                w->D_T[M],  w->ld, IN(rhoe), 1.0, OUT(rhoe));

            (*p_gbdmv)(trans, n, w->kl[D1], w->ku[D1],
                -phi*s->gamma,               REF(uy),
                w->D_T[D1], w->ld, IN(rhoe), 1.0, OUT(rhoe));

            (*p_gbdmv)(trans, n, w->kl[D2], w->ku[D2],
                phi*ginvRePr,                REF(nu),
                w->D_T[D2], w->ld, IN(rhoe), 1.0, OUT(rhoe));
        }

        (*p_gbmv)(trans, n, n, w->kl[M], w->ku[M],
            1.0, w->D_T[M], w->ld, IN(rhoe), 1.0, OUT(rhoe));
    }

#   undef LIKELY
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
