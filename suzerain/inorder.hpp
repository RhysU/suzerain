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
 * inorder.hpp: C++ wrappers for the C-based suzerain_inorder_* API
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_INORDER_HPP
#define __SUZERAIN_INORDER_HPP

#include <suzerain/inorder.h>

/** @file
 * Provides C++ wrappers for the C-based API in inorder.h.
 */

namespace suzerain {

/**
 * Provides C++ wrappers for C-based API in inorder.h.
 */
namespace inorder {

/** @see suzerain_inorder_wavenumber() */
inline
int wavenumber(const int N, const int i)
{
    return suzerain_inorder_wavenumber(N, i);
}

/** @see suzerain_inorder_wavenumber_abs() */
inline
int wavenumber_abs(const int N, const int i)
{
    return suzerain_inorder_wavenumber_abs(N, i);
}

/** @see suzerain_inorder_wavenumber_min() */
inline
int wavenumber_min(const int N)
{
    return suzerain_inorder_wavenumber_min(N);
}

/** @see suzerain_inorder_wavenumber_max() */
inline
int wavenumber_max(const int N)
{
    return suzerain_inorder_wavenumber_max(N);
}

/** see suzerain_inorder_wavenumber_valid() */
inline
bool wavenumber_valid(const int N, const int w)
{
    return suzerain_inorder_wavenumber_valid(N, w);
}

/** @see suzerain_inorder_wavenumber_diff() */
inline
int wavenumber_diff(const int N, const int dN, const int i)
{
    return suzerain_inorder_wavenumber_diff(N, dN, i);
}

/** @see suzerain_inorder_wavenumber_translatable() */
inline
int wavenumber_translatable(const int N, const int dN, const int i)
{
    return suzerain_inorder_wavenumber_translatable(N, dN, i);
}

/** @see suzerain_inorder_wavenumber_translate() */
inline
void wavenumber_translate(const int S, const int T, const int tb, const int te,
                          int& sb1, int& se1, int& sb2, int& se2,
                          int& tb1, int& te1, int& tb2, int& te2)
{
    return suzerain_inorder_wavenumber_translate(
        S, T, tb, te, &sb1, &se1, &sb2, &se2, &tb1, &te1, &tb2, &te2);
}

/** @see suzerain_inorder_index() */
inline
int index(const int N, const int w)
{
    return suzerain_inorder_index(N, w);
}

} // namespace inorder

} // namespace suzerain

#endif // __SUZERAIN_INORDER_HPP
