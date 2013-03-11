/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2012, 2013 The PECOS Development Team
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
 * @copydoc filterop.h
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/filterop.h>

#include <suzerain/common.h>
#include <suzerain/blas_et_al.h>
#include <suzerain/error.h>
#include <suzerain/gbmatrix.h>
#include <suzerain/inorder.h>

// Double precision as suzerain_filterop_workspace
#define SCALAR         double           /* Scalar type        */
#define BLAS(pre,post) pre ## d ## post /* Call BLAS routines */
#define WORK(pre,post) pre ##      post /* Workspace          */
#include "filterop.def"

// Complex double as suzerain_filteropz_workspace
#define SCALAR         complex_double   /* Scalar type        */
#define BLAS(pre,post) pre ## z ## post /* Call BLAS routines */
#define WORK(pre,post) pre ## z ## post /* Workspace          */
#include "filterop.def"


void suzerain_filteropz_source_apply(
    const complex_double alpha, complex_double * const x,
    const int Ny,
    const int Nx, const int dNx, const int dkbx, const int dkex,
    const int Nz, const int dNz, const int dkbz, const int dkez)
{
    assert(Ny >= 0);                // Sanity check Y direction
    assert(Nx >= 0);                // Sanity check X direction
    assert(dNx >= Nx);
    assert(dkex <= dNx);
    assert(dkbx <= dkex);
    assert(Nz >= 0);                // Sanity check Z direction
    assert(dNz >= Nz);
    assert(dkez <= dNz);
    assert(dkbz <= dkez);

    // Compute X, Z strides based on locally contiguous storage assumption
    const int sx = Ny, sz = (dkex - dkbx)*sx;

    // Precompute loop independent constants
    const int absmin_wm = suzerain_inorder_wavenumber_absmin(Nx);
    const int absmin_wn = suzerain_inorder_wavenumber_absmin(Nz);

    // Overwrite x with alpha*D*x for storage, dealiasing details
    // Generalized loops are hideous but permit one linear pass through memory
    int nfreqidx = INT_MAX;  // Hoisted out of loops to avoid uninitialized
    int mfreqidx = INT_MAX;  // usage warnings.  Convince yourself it's OK.
    for (int n = dkbz; n < dkez; ++n) {
        const int noff = sz*(n - dkbz);
        const int nkeeper  = suzerain_inorder_wavenumber_abs(dNz, n) <= absmin_wn;
        if (nkeeper) {
            // Relies on gsl_sf_pow_int(0.0, 0) == 1.0
//             const double nscale = gsl_sf_pow_int(twopioverLz*nfreqidx, dzcnt);
            for (int m = dkbx; m < dkex; ++m) {
                const int moff = noff + sx*(m - dkbx);
                const int mkeeper = suzerain_inorder_wavenumber_abs(dNx, m) <= absmin_wm;
                if (mkeeper) {
                    // FIXME: Call method to compute filter/filter source

                    // Relies on gsl_sf_pow_int(0.0, 0) == 1.0
//                     const double mscale
//                         = nscale*gsl_sf_pow_int(twopioverLx*mfreqidx, dxcnt);
//                     const complex_double malpha = mscale*alpha_ipow;
//                     suzerain_blas_zscal(Ny, malpha, x+moff, 1);
                } else {
                    memset(x+moff, 0, Ny*sizeof(x[0])); // Scale by zero
                }
            }
        } else {
            memset(x+noff, 0, sz*sizeof(x[0]));  // Scale by zero
        }
    }
}


void suzerain_filteropz_source_accumulate(
    const complex_double alpha, const complex_double * const x,
    const complex_double beta,        complex_double * const y,
    const int Ny,
    const int Nx, const int dNx, const int dkbx, const int dkex,
    const int Nz, const int dNz, const int dkbz, const int dkez)
{
    assert(dkex == dkbx || dkez == dkbz || x != y); // Trivial or not aliased?
    assert(Ny >= 0);                                // Sanity check Y direction
    assert(Nx >= 0);                                // Sanity check X direction
    assert(dNx >= Nx);
    assert(dkex <= dNx);
    assert(dkbx <= dkex);
    assert(Nz >= 0);                                // Sanity check Z direction
    assert(dNz >= Nz);
    assert(dkez <= dNz);
    assert(dkbz <= dkez);

    // Compute X, Z strides based on locally contiguous storage assumption
    const int sx = Ny, sz = (dkex - dkbx)*sx;

    // Precompute loop independent constants
    const int absmin_wm = suzerain_inorder_wavenumber_absmin(Nx);
    const int absmin_wn = suzerain_inorder_wavenumber_absmin(Nz);

    // Accumulate y <- alpha*D*x + beta*y for storage, dealiasing details
    // Generalized loops are hideous but permit one linear pass through memory
    int nfreqidx = INT_MAX;  // Hoisted out of loops to avoid uninitialized
    int mfreqidx = INT_MAX;  // usage warnings.  Convince yourself it's OK.
    for (int n = dkbz; n < dkez; ++n) {
        const int noff = sz*(n - dkbz);
        const int nkeeper = suzerain_inorder_wavenumber_abs(dNz, n) <= absmin_wn;
        if (nkeeper) {
            // Relies on gsl_sf_pow_int(0.0, 0) == 1.0
//             const double nscale = gsl_sf_pow_int(twopioverLz*nfreqidx, dzcnt);
            for (int m = dkbx; m < dkex; ++m) {
                const int moff = noff + sx*(m - dkbx);
                const int mkeeper = suzerain_inorder_wavenumber_abs(dNx, m) <= absmin_wm;
                if (mkeeper) {
                    // FIXME: Call method to compute filter/filter source

                    // Relies on gsl_sf_pow_int(0.0, 0) == 1.0
//                     const double mscale
//                         = nscale*gsl_sf_pow_int(twopioverLx*mfreqidx, dxcnt);
//                     const complex_double malpha = mscale*alpha_ipow;
//                     suzerain_blas_zaxpby(Ny, malpha, x+moff, 1, beta, y+moff, 1);
                } else {
                    suzerain_blas_zscal(Ny, beta, y+moff, 1);
                }
            }
        } else {
            suzerain_blas_zscal(sz,beta,y+noff,1);
        }
    }
}

