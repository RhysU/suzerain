//--------------------------------------------------------------------------
//
// Copyright (C) 2013-2014 Rhys Ulerich
// Copyright (C) 2013-2014 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Suzerain is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------

#ifndef SUZERAIN_PROFILE_H
#define SUZERAIN_PROFILE_H

/** @file
 * <a
 * href="www.boost.org/doc/libs/release/libs/preprocessor/">Boost.Preprocessor</a>
 * definitions driving \ref suzerain::profile. Isolated to ease use, to simplify
 * Isolated to ease use, to simplify re-use, and to facilitate debugging.
 */

/**
 * A Boost.Preprocessor sequence of tuples of quantities generally sampled in
 * wave space.  These are part of the data found in \ref suzerain::profile.
 */
#define SUZERAIN_PROFILE_WAVE \
    ((rho,   1)) /* scalar */ \
    ((rho_u, 2)) /* vector */

/**
 * A Boost.Preprocessor sequence of tuples of quantities generally sampled  in
 * physical space.  These are part of the data found in \ref suzerain::profile.
 */
#define SUZERAIN_PROFILE_PHYSICAL \
    ((a,     1))  /* scalar */    \
    ((H0,    1))  /* scalar */    \
    ((ke,    1))  /* scalar */    \
    ((mu,    1))  /* scalar */    \
    ((T,     1))  /* scalar */    \
    ((u,     2))  /* vector */

/**
 * A Boost.Preprocessor sequence of tuples of all profile quantities.
 * These are nothing but the data found in \ref suzerain::profile.
 */
#define SUZERAIN_PROFILE      \
    SUZERAIN_PROFILE_WAVE     \
    SUZERAIN_PROFILE_PHYSICAL

#endif // SUZERAIN_PROFILE_H
