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

#ifndef SUZERAIN_DIFFWAVE_HPP
#define SUZERAIN_DIFFWAVE_HPP

/** @file
 * Provides C++ wrappers for the C-based API in diffwave.h.
 */

#include <suzerain/complex.hpp>
#include <suzerain/diffwave.h>

namespace suzerain {

/**
 * Provides C++ wrappers for C-based API in diffwave.h.
 * @see diffwave.h
 */
namespace diffwave {

/** @see suzerain_diffwave_apply */
template< typename AlphaType, typename Complex >
typename boost::enable_if<
    suzerain::complex::traits::is_complex_t<Complex>
>::type apply(
    const int dxcnt,
    const int dzcnt,
    const AlphaType &alpha, Complex *x,
    const real_t Lx,
    const real_t Lz,
    const int Ny,
    const int Nx, const int dNx, const int dkbx, const int dkex,
    const int Nz, const int dNz, const int dkbz, const int dkez)
{
    complex_t alpha_complex;
    suzerain::complex::assign_complex(alpha_complex, alpha);
    return suzerain_diffwave_apply(
            dxcnt, dzcnt,
            alpha_complex,
            reinterpret_cast<complex_t *>(x),
            Lx,
            Lz,
            Ny,
            Nx, dNx, dkbx, dkex,
            Nz, dNz, dkbz, dkez);
}

/** @see suzerain_diffwave_accumulate */
template< typename AlphaType, typename BetaType,
          typename Complex1,  typename Complex2 >
typename boost::enable_if<boost::mpl::and_<
    suzerain::complex::traits::is_complex_t<Complex1>,
    suzerain::complex::traits::is_complex_t<Complex2>
>, void>::type accumulate(
    const int dxcnt,
    const int dzcnt,
    const AlphaType &alpha, const Complex1 *x,
    const BetaType  &beta,        Complex2 *y,
    const real_t Lx,
    const real_t Lz,
    const int Ny,
    const int Nx, const int dNx, const int dkbx, const int dkex,
    const int Nz, const int dNz, const int dkbz, const int dkez)
{
    complex_t alpha_complex;
    suzerain::complex::assign_complex(alpha_complex, alpha);
    complex_t beta_complex;
    suzerain::complex::assign_complex(beta_complex, beta);
    return suzerain_diffwave_accumulate(
            dxcnt, dzcnt,
            alpha_complex,
            reinterpret_cast<const complex_t *>(x),
            beta_complex,
            reinterpret_cast<complex_t *>(y),
            Lx,
            Lz,
            Ny,
            Nx, dNx, dkbx, dkex,
            Nz, dNz, dkbz, dkez);
}

} // namespace diffwave

} // namespace suzerain

#endif // SUZERAIN_DIFFWAVE_HPP
