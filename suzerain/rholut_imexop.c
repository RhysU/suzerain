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
        const double (*in_rho )[2],
        const double (*in_rhou)[2],
        const double (*in_rhov)[2],
        const double (*in_rhow)[2],
        const double (*in_rhoe)[2],
        const double beta,
        double (*out_rho )[2],
        double (*out_rhou)[2],
        double (*out_rhov)[2],
        double (*out_rhow)[2],
        double (*out_rhoe)[2])
{
    // Prepare several oft-used constants to aid readability
    static const int inc        = 1;
    static const int nrhs       = 1;
    static const int M          = 0;
    static const int D1         = 1;
    static const int D2         = 2;
    const int    n              = w->n;
    const double c_one[2]       = {          1, 0 };
    const double c_phi[2]       = {        phi, 0 };
    const double beta_by_phi[2] = { beta / phi, 0 };

    // Incrementally build out_rho
    suzerain_blas_zscal( // TODO incorporate beta_by_rho into first GBMV
            n, beta_by_phi, out_rho, inc);
    suzerain_bsplineop_accumulate_complex(
            M, nrhs, c_one, in_rho, inc, n, c_phi, out_rho, inc, n, w);

    // Incrementally build out_rhou
    suzerain_blas_zscal( // TODO incorporate beta_by_rho into first GBMV
            n, beta_by_phi, out_rhou, inc);
    suzerain_bsplineop_accumulate_complex(
            M, nrhs, c_one, in_rhou, inc, n, c_phi, out_rhou, inc, n, w);

    // Incrementally build out_rhov
    suzerain_blas_zscal( // TODO incorporate beta_by_rho into first GBMV
            n, beta_by_phi, out_rhov, inc);
    suzerain_bsplineop_accumulate_complex(
            M, nrhs, c_one, in_rhov, inc, n, c_phi, out_rhov, inc, n, w);

    // Incrementally build out_rhow
    suzerain_blas_zscal( // TODO incorporate beta_by_rho into first GBMV
            n, beta_by_phi, out_rhow, inc);
    suzerain_bsplineop_accumulate_complex(
            M, nrhs, c_one, in_rhow, inc, n, c_phi, out_rhow, inc, n, w);

    // Incrementally build out_rhoe
    suzerain_blas_zscal( // TODO incorporate beta_by_rho into first GBMV
            n, beta_by_phi, out_rhoe, inc);
    suzerain_bsplineop_accumulate_complex(
            M, nrhs, c_one, in_rhoe, inc, n, c_phi, out_rhoe, inc, n, w);
}
