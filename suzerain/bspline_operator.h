/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
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
 *
 * bspline_operator.h: operator creation routines for a bspline basis
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef PECOS_SUZERAIN_BSPLINE_OPERATOR_H
#define PECOS_SUZERAIN_BSPLINE_OPERATOR_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS   /* empty */
#endif

__BEGIN_DECLS

enum suzerain_bspline_operator_method {
    SUZERAIN_BSPLINE_OPERATOR_COLLOCATION_GREVILLE = 1
};

typedef struct {
    int order;
    int nbreak;
    int nderiv;
    int M_kl, M_ku, M_lda, M_storagesize;
    int D_kl, D_ku, D_lda, D_storagesize;
    double *M;
    double **D;
} suzerain_bspline_operator_workspace;

suzerain_bspline_operator_workspace *
suzerain_bspline_operator_alloc(int order,
                                int nbreak,
                                int nderiv,
                                enum suzerain_bspline_operator_method method);

void
suzerain_bspline_operator_free(suzerain_bspline_operator_workspace *w);

__END_DECLS

#endif // PECOS_SUZERAIN_BSPLINE_OPERATOR_H
