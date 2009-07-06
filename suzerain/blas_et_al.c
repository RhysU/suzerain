/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2009 The PECOS Development Team
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
 * blas_et_al.c: wraps external implementations of BLAS, LAPACK, et al.
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_BLAS_ET_AL_H__
#define __SUZERAIN_BLAS_ET_AL_H__

#include "config.h"

#ifdef HAVE_MKL
#include <mkl_types.h>
#include <mkl_blas.h>
#include <mkl_lapack.h>
#else
#error "No suitable BLAS and/or LAPACK library found during configuration"
#endif

#include <suzerain/blas_et_al.h>

int suzerain_lapack_dgbtrf(
        int m,
        int n,
        int kl,
        int ku,
        double *ab,
        int ldab,
        int *ipiv)
{
#ifdef HAVE_MKL
    MKL_INT _m    = m;
    MKL_INT _n    = n;
    MKL_INT _kl   = kl;
    MKL_INT _ku   = ku;
    MKL_INT _ldab = ldab;
    MKL_INT _info = 0;
#else
#error "Sanity failure"
#endif

    dgbtrf(&_m, &_n, &_kl, &_ku, ab, &_ldab, ipiv, &_info);

    return _info;
}

int suzerain_lapack_dgbtrs(
        char trans,
        int n,
        int kl,
        int ku,
        int nrhs,
        double *ab,
        int ldab,
        int *ipiv,
        double *b,
        int ldb)
{
#ifdef HAVE_MKL
    char _trans   = (char) trans;
    MKL_INT _n    = n;
    MKL_INT _kl   = kl;
    MKL_INT _ku   = ku;
    MKL_INT _nrhs = nrhs;
    MKL_INT _ldab = ldab;
    MKL_INT _ldb  = ldb;
    MKL_INT _info = 0;
#else
#error "Sanity failure"
#endif

    dgbtrs(&_trans, &_n, &_kl, &_ku, &_nrhs,
           ab, &_ldab, ipiv, b, &_ldb, &_info);

    return _info;
}

#endif /* __SUZERAIN_BLAS_ET_AL_H__ */
