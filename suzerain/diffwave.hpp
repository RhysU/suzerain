//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// diffwave.hpp: C++ wrappers for the C-based suzerain_diffwave_* API
// $Id$

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
template< typename AlphaType, typename Complex >
typename boost::enable_if<
    suzerain::complex::traits::is_complex_double<Complex>
>::type apply(
    const int dxcnt,
    const int dzcnt,
    const AlphaType &alpha, Complex *x,
    const double Lx,
    const double Lz,
    const int Ny,
    const int Nx, const int dNx, const int dkbx, const int dkex,
    const int Nz, const int dNz, const int dkbz, const int dkez)
{
    std::complex<double> alpha_complex;
    suzerain::complex::assign_complex(alpha_complex, alpha);
    return suzerain_diffwave_apply(
            dxcnt, dzcnt,
            alpha_complex,
            reinterpret_cast<std::complex<double> *>(x),
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
    suzerain::complex::traits::is_complex_double<Complex1>,
    suzerain::complex::traits::is_complex_double<Complex2>
>, void>::type accumulate(
    const int dxcnt,
    const int dzcnt,
    const AlphaType &alpha, const Complex1 *x,
    const BetaType  &beta,        Complex2 *y,
    const double Lx,
    const double Lz,
    const int Ny,
    const int Nx, const int dNx, const int dkbx, const int dkex,
    const int Nz, const int dNz, const int dkbz, const int dkez)
{
    std::complex<double> alpha_complex;
    suzerain::complex::assign_complex(alpha_complex, alpha);
    std::complex<double> beta_complex;
    suzerain::complex::assign_complex(beta_complex, beta);
    return suzerain_diffwave_accumulate(
            dxcnt, dzcnt,
            alpha_complex,
            reinterpret_cast<const std::complex<double> *>(x),
            beta_complex,
            reinterpret_cast<std::complex<double> *>(y),
            Lx,
            Lz,
            Ny,
            Nx, dNx, dkbx, dkex,
            Nz, dNz, dkbz, dkez);
}

} // namespace diffwave

} // namespace suzerain

#endif // __SUZERAIN_DIFFWAVE_HPP
