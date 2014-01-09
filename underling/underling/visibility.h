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

#ifndef UNDERLING_VISIBILITY_H
#define UNDERLING_VISIBILITY_H

/** \cond INTERNAL */

// Content and ideas modified from http://gcc.gnu.org/wiki/Visibility

// Generic helper definitions for shared library support
#if __GNUC__ >= 4
#  define UNDERLING_HELPER_SHARED_IMPORT __attribute__ ((visibility("default")))
#  define UNDERLING_HELPER_SHARED_EXPORT __attribute__ ((visibility("default")))
#  define UNDERLING_HELPER_SHARED_LOCAL  __attribute__ ((visibility("hidden")))
#else
#  define UNDERLING_HELPER_SHARED_IMPORT
#  define UNDERLING_HELPER_SHARED_EXPORT
#  define UNDERLING_HELPER_SHARED_LOCAL
#endif

// Now we use the generic helper definitions above to define UNDERLING_API and
// UNDERLING_LOCAL.  UNDERLING_API is used for the public API symbols.  It
// either imports or exports the relevant symbol.  UNDERLING_LOCAL is used for
// non-API symbols.

#ifdef UNDERLING_SHARED_EXPORTS // defined if building UNDERLING versus just using it
#  define UNDERLING_API UNDERLING_HELPER_SHARED_EXPORT
#else
#  define UNDERLING_API UNDERLING_HELPER_SHARED_IMPORT
#endif // UNDERLING_SHARED_EXPORTS
#define UNDERLING_LOCAL UNDERLING_HELPER_SHARED_LOCAL

/** \endcond */

#endif /* UNDERLING_VISIBILITY_H__ */
