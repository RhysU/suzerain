/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2014 Rhys Ulerich
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
 * @copydoc gbmatrix.h
 */

#include <suzerain/gbmatrix.h>

// C99 extern declaration for inlined function in gbmatrix.h
extern
int suzerain_gbmatrix_offset(int lda, int kl, int ku, int i, int j);

// C99 extern declaration for inlined function in gbmatrix.h
extern
int suzerain_gbmatrix_in_band(int kl, int ku, int i, int j);

// C99 extern declaration for inlined function in gbmatrix.h
extern
void *suzerain_gbmatrix_col(int m, int n, int kl, int ku,
                            void *a, int ld, size_t s, int j,
                            int *il, int *iu);

// C99 extern declaration for inlined function in gbmatrix.h
extern
void *suzerain_gbmatrix_row(int m, int n, int kl, int ku,
                            void *a, int ld, size_t s, int i,
                            int *jl, int *ju, int *inc);
