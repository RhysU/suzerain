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
 * blas_et_al.h: wraps external implementations of BLAS, LAPACK, et al.
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_BLAS_ET_AL_H__
#define __SUZERAIN_BLAS_ET_AL_H__

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

enum suzerain_lapack_transpose {
    SUZERAIN_LAPACK_TRANSPOSE_NONE      = 'N',
    SUZERAIN_LAPACK_TRANSPOSE_NORMAL    = 'T',
    SUZERAIN_LAPACK_TRANSPOSE_CONJUGATE = 'C'
};

int suzerain_lapack_dgbtrf(
        int m,
        int n,
        int kl,
        int ku,
        double *ab,
        int ldab,
        int *ipiv);

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
        int ldb);

__END_DECLS

#endif /* __SUZERAIN_BLAS_ET_AL_H__ */
