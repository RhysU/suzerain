/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2009 The PECOS Development Team
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
 * richardson.c: generalized Richardson extrapolation routines
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include "config.h"

#include <gsl/gsl_math.h>
#include <suzerain/error.h>
#include <suzerain/richardson.h>

int
suzerain_richardson_step(
    const gsl_vector * Aih,
    const gsl_vector * Aiht,
    const int ki,
    const double t,
    gsl_vector * Aip1h)
{
    const double tki       = gsl_pow_int(t, ki);
    const double inv_tkim1 = 1.0/(tki-1.0);
    int error;

    error = gsl_blas_dcopy(Aiht, Aip1h);
    if (error) return error;

    gsl_blas_dscal(tki*inv_tkim1, Aip1h);

    error = gsl_blas_daxpy(-inv_tkim1, Aih, Aip1h);
    if (error) return error;

    return SUZERAIN_SUCCESS;
}
