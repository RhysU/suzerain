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
 * filterop.c: Discrete filtering operator routines
 * $Id$
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.h>
#pragma hdrstop
#include <suzerain/filterop.h>
#include <suzerain/blas_et_al.h>
#include <suzerain/error.h>
#include <suzerain/gbmatrix.h>

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
