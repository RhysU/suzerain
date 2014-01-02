/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2014 Rhys Ulerich
 * Copyright (C) 2012-2014 The PECOS Development Team
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
 */

/** @file
 * @copydoc inorder.h
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/inorder.h>

#include <suzerain/common.h>

// C99 extern declarations for inline functions from inorder.h
// See "C99 model" in http://www.greenend.org.uk/rjk/2003/03/inline.html
extern
int suzerain_inorder_wavenumber(const int N, const int i);

extern
int suzerain_inorder_wavenumber_abs(const int N, const int i);

extern
int suzerain_inorder_wavenumber_min(const int N);

extern
int suzerain_inorder_wavenumber_absmin(const int N);

extern
int suzerain_inorder_wavenumber_max(const int N);

extern
int suzerain_inorder_wavenumber_nyquist(const int N, const int w);

extern
int suzerain_inorder_wavenumber_nyquist_index(const int N, const int i);

extern
int suzerain_inorder_wavenumber_imagzero(const int N, const int w);

extern
int suzerain_inorder_wavenumber_imagzero_index(const int N, const int i);

extern
int suzerain_inorder_wavenumber_valid(const int N, const int w);

extern
int suzerain_inorder_wavenumber_diff(const int N, const int dN, const int i);

extern
int suzerain_inorder_wavenumber_translatable(
        const int N, const int dN, const int i);

extern
int suzerain_inorder_index(const int N, const int w);

void suzerain_inorder_wavenumber_translate(
        const int S, const int T, const int tb, const int te,
        int* sb1, int* se1, int* sb2, int* se2,
        int* tb1, int* te1, int* tb2, int* te2)
{
    // This could not be much uglier, could it?  The inline routines used to
    // cobble together this logic are a minefield of assertions.  Beware.

    // Sanity check incoming arguments
    assert(0 < S);
    assert(0 < T);
    assert(0 <= tb && tb <= te && te <= T);

    // If necessary, transform problem so S <= T and solve it
    if (S > T) {
        int sb = S;
        if (tb < T) {
            sb = suzerain_inorder_index(
                    S, suzerain_inorder_wavenumber(T, tb));
        }
        int se = S;
        if (te < T) {
            se = suzerain_inorder_index(
                    S, suzerain_inorder_wavenumber(T, te));
        }
        return suzerain_inorder_wavenumber_translate(
                T, S, sb, se, tb1, te1, tb2, te2, sb1, se1, sb2, se2);
    }

    // Sanity check preconditions hold for following algorithm
    assert(S <= T);

    // Find <tt>[tb1,te1)</tt> and <tt>[tb2,te2)</tt>
    *tb1 = tb;
    while (*tb1 < te && !suzerain_inorder_wavenumber_translatable(S, T, *tb1)) {
        ++*tb1;
    }
    *te1 = *tb1;
    while (*te1 < te &&  suzerain_inorder_wavenumber_translatable(S, T, *te1)) {
        ++*te1;
    }
    *tb2 = *te1;
    while (*tb2 < te && !suzerain_inorder_wavenumber_translatable(S, T, *tb2)) {
        ++*tb2;
    }
    *te2 = *tb2;
    while (*te2 < te &&  suzerain_inorder_wavenumber_translatable(S, T, *te2)) {
        ++*te2;
    }

    // Sanity check intermediate results
    assert(0 <= *tb1 && *tb1 <= T);
    assert(0 <= *te1 && *te1 <= T);
    assert(0 <= *tb2 && *tb2 <= T);
    assert(0 <= *te2 && *te2 <= T);
    assert(*tb1 <= *te1 && *te1 <= *tb2 && *tb2 <= *te2);

    // Find <tt>[sb1,se1)</tt>
    *sb1 = S;
    if (*tb1 < T) {
        const int f = suzerain_inorder_wavenumber(T, *tb1);
        if (suzerain_inorder_wavenumber_valid(S, f)) {
            *sb1 = suzerain_inorder_index(S, f);
        }
    }
    *se1 = *sb1 + (*te1 - *tb1);

    // Find <tt>[sb2,se2)</tt>
    *sb2 = *se1;
    if (*tb2 < T) {
        const int f = suzerain_inorder_wavenumber(T, *tb2);
        if (suzerain_inorder_wavenumber_valid(S, f)) {
            *sb2 = suzerain_inorder_index(S, f);
        }
    }
    *se2 = *sb2 + (*te2 - *tb2);
}
