/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2012 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * rholut_imexop.c: hybrid implicit/explicit operator apply and solve
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

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

    // Accumulate the requested portions of the M + \varphi L operator.
    // The zscal beta/phi and M gbmv phi coefficients scale output by beta.

    if (out_rho ) {  // Accumulate density terms into out_rho
        suzerain_blas_zscal(n, beta/phi, out_rho, inc);

        suzerain_bsplineop_accumulate_complex(
                M, nrhs, 1.0, in_rho, inc, n, phi, out_rho, inc, n, w);
    }

    if (out_rhou) {  // Accumulate X momentum terms into out_rhou
        suzerain_blas_zscal(n, beta/phi, out_rhou, inc);

        suzerain_bsplineop_accumulate_complex(
                M, nrhs, 1.0, in_rhou, inc, n, phi, out_rhou, inc, n, w);
    }

    if (out_rhov) {  // Accumulate Y momentum terms into out_rhov
        suzerain_blas_zscal(n, beta/phi, out_rhov, inc);

        suzerain_bsplineop_accumulate_complex(
                M, nrhs, 1.0, in_rhov, inc, n, phi, out_rhov, inc, n, w);
    }

    if (out_rhow) {  // Accumulate Z momentum terms into out_rhow
        suzerain_blas_zscal(n, beta/phi, out_rhow, inc);

        suzerain_bsplineop_accumulate_complex(
                M, nrhs, 1.0, in_rhow, inc, n, phi, out_rhow, inc, n, w);
    }

    if (out_rhoe) {  // Accumulate total energy terms into out_rhoe
        suzerain_blas_zscal(n, beta/phi, out_rhoe, inc);

        suzerain_bsplineop_accumulate_complex(
                M, nrhs, 1.0, in_rhoe, inc, n, phi, out_rhoe, inc, n, w);
    }
}

// suzerain_rholut_imexop_pack{c,f} differ trivially
// use preprocessor to generate both from the same source template

// suzerain_rholut_imexop_packc
#define FUNCNAME()  suzerain_rholut_imexop_packc
#define PACKFUNC(x) suzerain_bsmbsm_ ## x ## packc
#include "rholut_imexop.def"

// suzerain_rholut_imexop_packf
#define FUNCNAME()  suzerain_rholut_imexop_packf
#define PACKFUNC(x) suzerain_bsmbsm_ ## x ## packf
#include "rholut_imexop.def"
