/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
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
 * diffwave.hpp: C++ wrappers for the C-based suzerain_diffwave_* API
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_DIFFWAVE_HPP
#define __SUZERAIN_DIFFWAVE_HPP

#include <suzerain/diffwave.h>
#include <suzerain/complex.hpp>

/** @file
 * Provides C++ wrappers for the C-based API in diffwave.h.
 */

namespace suzerain {

/**
 * Provides C++ wrappers for C-based API in diffwave.h.
 * @see diffwave.h
 */
namespace diffwave {

/** @see suzerain_diffwave_apply */
template< typename Complex1,
          typename Complex2 >
typename boost::enable_if<boost::mpl::and_<
    suzerain::complex::traits::is_complex_double<Complex1>,
    suzerain::complex::traits::is_complex_double<Complex2>
>, void>::type apply(
    const int dxcnt,
    const int dzcnt,
    const Complex1 &alpha, Complex2 *x,
    const double Lx,
    const double Lz,
    const int Ny,
    const int Nx, const int dNx, const int dkbx, const int dkex,
    const int Nz, const int dNz, const int dkbz, const int dkez)
{
    return suzerain_diffwave_apply(
            dxcnt, dzcnt,
            reinterpret_cast<const double *>(&alpha),
            reinterpret_cast<double (*)[2]>(x),
            Lx,
            Lz,
            Ny,
            Nx, dNx, dkbx, dkex,
            Nz, dNz, dkbz, dkez);
}

/** @see suzerain_diffwave_accumulate */
template< typename Complex1,
          typename Complex2,
          typename Complex3,
          typename Complex4 >
typename boost::enable_if<boost::mpl::and_<
    suzerain::complex::traits::is_complex_double<Complex1>,
    suzerain::complex::traits::is_complex_double<Complex2>,
    suzerain::complex::traits::is_complex_double<Complex3>,
    suzerain::complex::traits::is_complex_double<Complex4>
>, void>::type accumulate(
    const int dxcnt,
    const int dzcnt,
    const Complex1 &alpha, const Complex2 *x,
    const Complex3 &beta, Complex4 *y,
    const double Lx,
    const double Lz,
    const int Ny,
    const int Nx, const int dNx, const int dkbx, const int dkex,
    const int Nz, const int dNz, const int dkbz, const int dkez)
{
    return suzerain_diffwave_accumulate(
            dxcnt, dzcnt,
            reinterpret_cast<const double *>(&alpha),
            reinterpret_cast<const double (*)[2]>(x),
            reinterpret_cast<const double *>(&beta),
            reinterpret_cast<double (*)[2]>(y),
            Lx,
            Lz,
            Ny,
            Nx, dNx, dkbx, dkex,
            Nz, dNz, dkbz, dkez);
}

} // namespace diffwave

} // namespace suzerain

#endif // __SUZERAIN_DIFFWAVE_HPP
