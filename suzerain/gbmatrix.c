/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2010, 2011, 2012 Rhys Ulerich
 * Copyright (C) 2012 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *--------------------------------------------------------------------------
 * gbmatrix.c: General band storage matrix helpers
 * $Id$
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
