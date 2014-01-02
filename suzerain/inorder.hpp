//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
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

#ifndef SUZERAIN_INORDER_HPP
#define SUZERAIN_INORDER_HPP

/** @file
 * Provides C++ wrappers for the C-based API in inorder.h.
 */

#include <suzerain/inorder.h>

namespace suzerain {

/**
 * Provides C++ wrappers for C-based API in inorder.h.
 */
namespace inorder {

/** @copydoc suzerain_inorder_wavenumber() */
inline
int wavenumber(const int N, const int i)
{
    return suzerain_inorder_wavenumber(N, i);
}

/** @copydoc suzerain_inorder_wavenumber_abs() */
inline
int wavenumber_abs(const int N, const int i)
{
    return suzerain_inorder_wavenumber_abs(N, i);
}

/** @copydoc suzerain_inorder_wavenumber_min() */
inline
int wavenumber_min(const int N)
{
    return suzerain_inorder_wavenumber_min(N);
}

/** @copydoc suzerain_inorder_wavenumber_absmin() */
inline
int wavenumber_absmin(const int N)
{
    return suzerain_inorder_wavenumber_absmin(N);
}

/** @copydoc suzerain_inorder_wavenumber_max() */
inline
int wavenumber_max(const int N)
{
    return suzerain_inorder_wavenumber_max(N);
}

/** @copydoc suzerain_inorder_wavenumber_nyquist() */
inline
int wavenumber_nyquist(const int N, const int w)
{
    return suzerain_inorder_wavenumber_nyquist(N, w);
}

/** @copydoc suzerain_inorder_wavenumber_nyquist_index() */
inline
int wavenumber_nyquist_index(const int N, const int i)
{
    return suzerain_inorder_wavenumber_nyquist_index(N, i);
}

/** @copydoc suzerain_inorder_wavenumber_imagzero() */
inline
int wavenumber_imagzero(const int N, const int w)
{
    return suzerain_inorder_wavenumber_imagzero(N, w);
}

/** @copydoc suzerain_inorder_wavenumber_imagzero_index() */
inline
int wavenumber_imagzero_index(const int N, const int i)
{
    return suzerain_inorder_wavenumber_imagzero_index(N, i);
}

/** @copydoc suzerain_inorder_wavenumber_valid() */
inline
bool wavenumber_valid(const int N, const int w)
{
    return suzerain_inorder_wavenumber_valid(N, w);
}

/** @copydoc suzerain_inorder_wavenumber_diff() */
inline
int wavenumber_diff(const int N, const int dN, const int i)
{
    return suzerain_inorder_wavenumber_diff(N, dN, i);
}

/** @copydoc suzerain_inorder_wavenumber_translatable() */
inline
int wavenumber_translatable(const int N, const int dN, const int i)
{
    return suzerain_inorder_wavenumber_translatable(N, dN, i);
}

/** @copydoc suzerain_inorder_wavenumber_translate() */
inline
void wavenumber_translate(const int S, const int T, const int tb, const int te,
                          int& sb1, int& se1, int& sb2, int& se2,
                          int& tb1, int& te1, int& tb2, int& te2)
{
    return suzerain_inorder_wavenumber_translate(
        S, T, tb, te, &sb1, &se1, &sb2, &se2, &tb1, &te1, &tb2, &te2);
}

/** @copydoc suzerain_inorder_index() */
inline
int index(const int N, const int w)
{
    return suzerain_inorder_index(N, w);
}

} // namespace inorder

} // namespace suzerain

#endif // SUZERAIN_INORDER_HPP
