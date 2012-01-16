/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010, 2011, 2012 The PECOS Development Team
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
 * mpi.h: Provides an incompatibility-hiding wrapper over the usual <mpi.h>
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifdef __GNUC__
/* Some older implementations trigger many warnings.  Suppress them. */
#pragma GCC system_header
#endif
#include <mpi.h>

#ifndef __SUZERAIN_MPI_H__
#define __SUZERAIN_MPI_H__

/** @file
 * Provides a small wrapper over the MPI vendor's mpi.h that hides
 * incompatibilities as much as possible.
 */

#ifdef __cplusplus
extern "C" {
#endif

/* The MPI standard uses MPI_MAX_OBJECT_NAME for MPI_Comm_get_name, etc.  */
#ifndef MPI_MAX_OBJECT_NAME
/* mpvapich only provides MPI_MAX_NAME_STRING */
#ifdef MPI_MAX_NAME_STRING
#define MPI_MAX_OBJECT_NAME MPI_MAX_NAME_STRING
#endif /* MPI_MAX_NAME_STRING */
#endif /* MPI_MAX_OBJECT_NAME */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __SUZERAIN_MPI_H__ */
