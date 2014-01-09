/***************************************************************************
//-----------------------------------------------------------------------bl-
// underling 0.3.1: an FFTW MPI-based library for 3D pencil decompositions
// http://red.ices.utexas.edu/projects/underling
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
//
// This file is part of underling.
//
// underling is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// underling is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with underling.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------el-
// $Id$
***************************************************************************/

#ifndef UNDERLING_COMMON_H
#define UNDERLING_COMMON_H

#if defined(__GNUC__) && __GNUC__ >= 3
/**
Provides hint to the compiler to optimize for the expression being false.
@param expr boolean expression.
@returns value of expression.
*/
#define UNDERLING_UNLIKELY(expr) __builtin_expect(expr, 0)
#else
/**
Provides hint to the compiler to optimize for the expression being false.
@param expr boolean expression.
@returns value of expression.
**/
#define UNDERLING_UNLIKELY(expr) expr
#endif

/**
 * Explicitly marks a variable as unused to avoid compiler warnings.
 *
 * @param variable Variable to be marked.
 */
#define UNDERLING_UNUSED(variable) do {(void)(variable);} while (0)

/* The MPI standard uses MPI_MAX_OBJECT_NAME for MPI_Comm_get_name, etc.  */
#ifndef MPI_MAX_OBJECT_NAME
/* mpvapich only provides MPI_MAX_NAME_STRING */
#ifdef MPI_MAX_NAME_STRING
#define MPI_MAX_OBJECT_NAME MPI_MAX_NAME_STRING
#endif /* MPI_MAX_NAME_STRING */
#endif /* MPI_MAX_OBJECT_NAME */

#endif // UNDERLING_COMMON_H
