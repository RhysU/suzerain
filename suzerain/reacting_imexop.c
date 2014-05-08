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
 * @copydoc reacting_imexop.h
 */

#include <suzerain/reacting_imexop.h>

#include <suzerain/common.h>
#include <suzerain/blas_et_al/blas_et_al.h>
#include <suzerain/bsmbsm.h>

void
suzerain_reacting_flow_imexop_accumulate(
        const complex_double phi,
        const suzerain_reacting_imexop_scenario * const s,
        const suzerain_reacting_imexop_ref      * const r,
        const suzerain_reacting_imexop_refld    * const ld,
        const suzerain_bsplineop_workspace    * const w,
        const int imagzero,
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
    // When you modify this routine, you must also modify reacting_imexop.def so
    // that operator accumulation-without-assembly and assembly match.  The
    // test cases in tests/test_reacting_imexop.cpp are invaluable in checking
    // the coherence of these two pieces of logic.

    // Sanity checks
    assert(w->nderiv >= 2);            // Adequate workspace?
    assert(!(!in_rho_E ^ !out_rho_E)); // Both non-NULL or NULL?
    assert(!(!in_rho_u ^ !out_rho_u)); // ditto
    assert(!(!in_rho_v ^ !out_rho_v)); // ditto
    assert(!(!in_rho_w ^ !out_rho_w)); // ditto
    assert(!(!in_rho   ^ !out_rho  )); // ditto

    // Prepare shorthand for some useful derived values
    const double ap43        = s->alpha + 4.0/3.0;

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

    // In the first case, treat x arguments as complex with stride one...
    int inc_in = 1;
    gbmv_t      *p_gbmv      = (gbmv_t *)      &        suzerain_blas_zgbmv_d_z;
    gbdmv_t     *p_gbdmv     = (gbdmv_t *)     &    suzerain_blasext_zgbdmv_d_z;
    gbddmv_t    *p_gbddmv    = (gbddmv_t *)    &   suzerain_blasext_zgbddmv_d_z;

    // ...but in the second case, treat x arguments as real with stride two.
    if (SUZERAIN_UNLIKELY(imagzero)) {
        inc_in = 2;
        p_gbmv      = (gbmv_t *)      &        suzerain_blas_zgbmv_d_d;
        p_gbdmv     = (gbdmv_t *)     &    suzerain_blasext_zgbdmv_d_d;
        p_gbddmv    = (gbddmv_t *)    &   suzerain_blasext_zgbddmv_d_d;
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


    if (LIKELY(in_rho_E)) {  // Accumulate total energy terms into out_rho_E

        suzerain_blas_zscal(n, beta, OUT(rho_E));

        /* in_rho_E */ {
            (*p_gbdmv)(trans, n, w->kl[D1], w->ku[D1],
                -phi,               REF(vp_rE),
                w->D_T[D1], w->ld, IN(rho_E), 1.0, OUT(rho_E));

            (*p_gbdmv)(trans, n, w->kl[D2], w->ku[D2],
                phi,                         REF(korCv),
                w->D_T[D2], w->ld, IN(rho_E), 1.0, OUT(rho_E));
        }

        if (LIKELY(in_rho_u)) {

            (*p_gbdmv)(trans, n, w->kl[D1], w->ku[D1],
                -phi,              REF(vp_ru),
                w->D_T[D1],  w->ld, IN(rho_u), 1.0, OUT(rho_E));

        }

        if (LIKELY(in_rho_v)) {

            (*p_gbdmv)(trans, n, w->kl[D1], w->ku[D1],
                phi,                        REF(Ce_rv),
                w->D_T[D1], w->ld, IN(rho_v), 1.0, OUT(rho_E));

        }

        if (LIKELY(in_rho_w)) {

            (*p_gbdmv)(trans, n, w->kl[D1], w->ku[D1],
                -phi,              REF(vp_rw),
                w->D_T[D1],  w->ld, IN(rho_w), 1.0, OUT(rho_E));

        }

        if (LIKELY(in_rho)) {

            (*p_gbdmv)(trans, n, w->kl[D1], w->ku[D1],
                phi,                        REF(Ce_rho),
                w->D_T[D1], w->ld, IN(rho), 1.0, OUT(rho_E));

        }

        (*p_gbmv)(trans, n, n, w->kl[M], w->ku[M],
            1.0, w->D_T[M], w->ld, IN(rho_E), 1.0, OUT(rho_E));
    }

    if (LIKELY(in_rho_u)) {  // Accumulate X momentum terms into out_rho_u

        suzerain_blas_zscal(n, beta, OUT(rho_u));

        if (LIKELY(in_rho_E)) {/* NOP */};

        /* in_rho_u */ {

            (*p_gbdmv)(trans, n, w->kl[D2], w->ku[D2],
                phi,                      REF(nu),
                w->D_T[D2], w->ld, IN(rho_u), 1.0, OUT(rho_u));
        }

        if (LIKELY(in_rho_v)) {

            (*p_gbdmv)(trans, n, w->kl[D1], w->ku[D1],
                -phi,                      REF(ux),
                w->D_T[D1], w->ld, IN(rho_v), 1.0, OUT(rho_u));
        }

        if (LIKELY(in_rho_w)) {/* NOP */};

        if (LIKELY(in_rho)) {

            (*p_gbdmv)(trans, n, w->kl[D1], w->ku[D1],
                phi,                       REF(uxuy),
                w->D_T[D1], w->ld, IN(rho), 1.0, OUT(rho_u));

        }


        (*p_gbmv)(trans, n, n, w->kl[M], w->ku[M],
            1.0, w->D_T[M], w->ld, IN(rho_u), 1.0, OUT(rho_u));
    }

    if (LIKELY(in_rho_v)) {  // Accumulate Y momentum terms into out_rho_v

        suzerain_blas_zscal(n, beta, OUT(rho_v));

        if (LIKELY(in_rho_E)) {

            (*p_gbdmv)(trans, n, w->kl[D1], w->ku[D1],
                -phi,              REF(p_rE),
                w->D_T[D1],  w->ld, IN(rho_E), 1.0, OUT(rho_v));
        }


        if (LIKELY(in_rho_u)) {

            (*p_gbdmv)(trans, n, w->kl[D1], w->ku[D1],
                -phi,              REF(p_ru),
                w->D_T[D1],  w->ld, IN(rho_u), 1.0, OUT(rho_v));
        }

        /* in_rho_v */ {

            (*p_gbddmv)(trans, n, w->kl[D1], w->ku[D1],
                -phi,              REF(uy),
                 phi,              REF(vp_rE), // FIXME: correct but misleading (since vp_rE = v*(gam-1))
                w->D_T[D1], w->ld, IN(rho_v), 1.0, OUT(rho_v));

            (*p_gbdmv)(trans, n, w->kl[D2], w->ku[D2],
                phi*ap43,          REF(nu),
                w->D_T[D2], w->ld, IN(rho_v), 1.0, OUT(rho_v));
        }

        if (LIKELY(in_rho_w)) {

            (*p_gbdmv)(trans, n, w->kl[D1], w->ku[D1],
                -phi,                 REF(p_rw),
                w->D_T[D1],  w->ld, IN(rho_w), 1.0, OUT(rho_v));
        }

        if (LIKELY(in_rho)) {

            (*p_gbdmv)(trans, n, w->kl[D1], w->ku[D1],
                 phi,                 REF(Cmy_rho),
                w->D_T[D1], w->ld, IN(rho), 1.0, OUT(rho_v));

        }

        (*p_gbmv)(trans, n, n, w->kl[M], w->ku[M],
            1.0, w->D_T[M], w->ld, IN(rho_v), 1.0, OUT(rho_v));
    }

    if (LIKELY(in_rho_w)) {  // Accumulate Z momentum terms into out_rho_w

        suzerain_blas_zscal(n, beta, OUT(rho_w));

        if (LIKELY(in_rho_E)) {/* NOP */};

        if (LIKELY(in_rho_u)) {/* NOP */};

        if (LIKELY(in_rho_v)) {

            (*p_gbdmv)(trans, n, w->kl[D1], w->ku[D1],
                -phi,                      REF(uz),
                w->D_T[D1], w->ld, IN(rho_v), 1.0, OUT(rho_w));
        }

        /* in_rho_w */ {

            (*p_gbdmv)(trans, n, w->kl[D2], w->ku[D2],
                phi,                      REF(nu),
                w->D_T[D2], w->ld, IN(rho_w), 1.0, OUT(rho_w));
        }

        if (LIKELY(in_rho)) {

            (*p_gbdmv)(trans, n, w->kl[D1], w->ku[D1],
                 phi,                      REF(uyuz),
                w->D_T[D1], w->ld, IN(rho), 1.0, OUT(rho_w));

        }

        (*p_gbmv)(trans, n, n, w->kl[M], w->ku[M],
            1.0, w->D_T[M], w->ld, IN(rho_w), 1.0, OUT(rho_w));
    }

    if (LIKELY(in_rho )) {  // Accumulate density terms into out_rho

        suzerain_blas_zscal(n, beta, OUT(rho));

        if (LIKELY(in_rho_E)) {/* NOP */};

        if (LIKELY(in_rho_u)) {/* NOP */};

        if (LIKELY(in_rho_v)) {/* NOP */};

        if (LIKELY(in_rho_w)) {/* NOP */};

        (*p_gbmv)(trans, n, n, w->kl[M], w->ku[M],
            1.0, w->D_T[M], w->ld, IN(rho), 1.0, OUT(rho));
    }

#   undef LIKELY
#   undef REF
#   undef IN
#   undef OUT

}


void
suzerain_reacting_species_imexop_accumulate(
        const complex_double phi,
        const suzerain_reacting_imexop_scenario * const s,
        const suzerain_reacting_imexop_ref      * const r,
        const suzerain_reacting_imexop_refld    * const ld,
        const suzerain_bsplineop_workspace    * const w,
        const int imagzero,
        const complex_double *in_rho_s,
        const complex_double beta,
        complex_double *out_rho_s )
{
    SUZERAIN_UNUSED(s);

    // When you modify this routine, you must also modify reacting_imexop.def so
    // that operator accumulation-without-assembly and assembly match.  The
    // test cases in tests/test_reacting_imexop.cpp are invaluable in checking
    // the coherence of these two pieces of logic.

    // Sanity checks
    assert(w->nderiv >= 2);  // Adequate workspace?
    assert( in_rho_s!=NULL); // Require both input and
    assert(out_rho_s!=NULL); // output rho_s

    // Accumulate the requested portions of the M + \varphi L operator.  Scale
    // output by beta, accumulate non-mass contributions, and finally
    // accumulate the mass contributions.  Mass contributions come last as they
    // are expected to have magnitudes much larger than phi-- this helps to
    // ensure better rounding of \varphi L contributions.

    // Two use cases are accommodated.  The first (and vastly more common) case
    // is applying the operator to in_rho_s as provided.  The second
    // case applies the operator for Im(in_rho_s) == 0.  Imaginary
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

    // In the first case, treat x arguments as complex with stride one...
    int inc_in = 1;
    gbmv_t      *p_gbmv      = (gbmv_t *)      &        suzerain_blas_zgbmv_d_z;
    gbdmv_t     *p_gbdmv     = (gbdmv_t *)     &    suzerain_blasext_zgbdmv_d_z;

    // ...but in the second case, treat x arguments as real with stride two.
    if (SUZERAIN_UNLIKELY(imagzero)) {
        inc_in = 2;
        p_gbmv      = (gbmv_t *)      &        suzerain_blas_zgbmv_d_d;
        p_gbdmv     = (gbdmv_t *)     &    suzerain_blasext_zgbdmv_d_d;
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
    SUZERAIN_UNUSED(D1);


    suzerain_blas_zscal(n, beta, OUT(rho_s));

    (*p_gbdmv)(trans, n, w->kl[D2], w->ku[D2],
               phi,                         REF(Ds),
               w->D_T[D2], w->ld, IN(rho_s), 1.0, OUT(rho_s));

    (*p_gbmv)(trans, n, n, w->kl[M], w->ku[M],
              1.0, w->D_T[M], w->ld, IN(rho_s), 1.0, OUT(rho_s));

#   undef LIKELY
#   undef REF
#   undef IN
#   undef OUT

}


// suzerain_reacting_imexop_pack{c,f} differ trivially
// use preprocessor to generate both from the same source template

// suzerain_reacting_imexop_packc
#define FLOW_FUNCNAME()    suzerain_reacting_flow_imexop_packc
#define SPECIES_FUNCNAME() suzerain_reacting_species_imexop_packc
#define PACK(x)            suzerain_bsmbsm_ ## x ## packc
#include "reacting_imexop.def"

// suzerain_reacting_imexop_pack
#define FLOW_FUNCNAME()    suzerain_reacting_flow_imexop_packf
#define SPECIES_FUNCNAME() suzerain_reacting_species_imexop_packf
#define PACK(x)            suzerain_bsmbsm_ ## x ## packf
#include "reacting_imexop.def"
