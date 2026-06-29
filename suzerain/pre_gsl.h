/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2011-2014 Rhys Ulerich
 * Copyright (C) 2012-2014 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This file is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This file is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this file.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 */

#ifndef SUZERAIN_PRE_GSL_H
#define SUZERAIN_PRE_GSL_H

/** @file
 * An MPI-friendly wrapper around GSL's IEEE environment setup.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Call <a
 * href="http://www.gnu.org/software/gsl/manual/html_node/Setting-up-your-IEEE-environment.html">
 * gsl_ieee_env_setup()</a> in an MPI-friendly way.  Processes with nonzero
 * rank will have their \c stdout and \c stderr redirected to
 * <tt>/dev/null</tt> during \c gsl_ieee_env_setup.
 *
 * @param rank MPI rank of this process.  Output from \c gsl_ieee_env_setup
 *             will only be observable from rank zero.
 */
void
mpi_gsl_ieee_env_setup(const int rank);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // SUZERAIN_PRE_GSL_H
