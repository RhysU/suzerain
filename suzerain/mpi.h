/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2010, 2011, 2012 The PECOS Development Team
 *
 * This file is part of Suzerain.
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *--------------------------------------------------------------------------
 * mpi.h: Provides an incompatibility-hiding wrapper over the usual <mpi.h>
 * $Id$
 */

#ifndef __SUZERAIN_MPI_H
#define __SUZERAIN_MPI_H

#ifdef __GNUC__
/* Some older implementations trigger many warnings.  Suppress them. */
#pragma GCC system_header
#endif
#include <mpi.h>

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

#endif /* __SUZERAIN_MPI_H */
