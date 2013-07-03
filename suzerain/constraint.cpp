//--------------------------------------------------------------------------
//
// Copyright (C) 2013 Rhys Ulerich
// Copyright (C) 2013 The PECOS Development Team
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

/** @file
 * @copydoc constraint.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/constraint.hpp>

#include <suzerain/bspline.hpp>
#include <suzerain/gbmatrix.h>

namespace suzerain {

namespace constraint {

base::~base()
{
    // NOP
}

uniform::uniform(bspline &b)
{
    shape.setOnes(b.n());
}

uniform::uniform(const bsplineop &bop)
{
    shape.setOnes(bop.n());
}

coefficient::coefficient(const bsplineop &bop, const int i)
{
    // Apply the mass matrix to express coefficient i via collocation points
    coeff.setZero(bop.n());
    coeff[i] = 1;
    bop.apply(0, 1.0, coeff.data(), coeff.innerStride());
}

lower::lower(const bsplineop &bop, const int nderiv)
{
    // Reconstructing the first row of D^{(nderiv)} in coeff
    // requires working with the first column of D_T
    coeff.setZero(bop.n());
    int il, iu;
    real_t * const col = (real_t *) suzerain_gbmatrix_col(
            bop.n(), bop.n(), bop.kl(nderiv), bop.ku(nderiv),
            (void *) bop.D_T(nderiv), bop.ld(), sizeof(real_t),
            /*first*/0, &il, &iu);
    for (int i = il; i < iu; ++i) {
        coeff[i] = col[i];
    }
}

upper::upper(const bsplineop &bop, const int nderiv)
{
    // Reconstructing the last row of D^{(nderiv)} in coeff
    // requires working with the last column of D_T
    coeff.setZero(bop.n());
    int il, iu;
    real_t * const col = (real_t *) suzerain_gbmatrix_col(
            bop.n(), bop.n(), bop.kl(nderiv), bop.ku(nderiv),
            (void *) bop.D_T(nderiv), bop.ld(), sizeof(real_t),
            /*last*/bop.n()-1, &il, &iu);
    for (int i = il; i < iu; ++i) {
        coeff[i] = col[i];
    }
}

bulk::bulk(bspline &b)
{
    coeff.resize(b.n());
    b.integration_coefficients(0, coeff.data());
    real_t L = b.breakpoint(b.nbreak()-1) - b.breakpoint(0);
    coeff /= L;
}

constant::constant(const real_t target)
    : t(target)
{
    // NOP
}

real_t
constant::target() const
{
    return t;
}

bool
constant::enabled() const
{
    return !(boost::math::isnan)(t);
}

reference::reference(const real_t& target)
    : t(&target)
{
    // NOP
}

real_t
reference::target() const
{
    return *t;
}

bool
reference::enabled() const
{
    return !(boost::math::isnan)(*t);
}

disabled::disabled(bspline &b)
    : uniform(b)
{
    // NOP
}

real_t
disabled::target() const
{
    return std::numeric_limits<real_t>::quiet_NaN();
}

bool
disabled::enabled() const
{
    return false;
}

constant_lower::constant_lower(
        const real_t target,
        const bsplineop &bop,
        const int nderiv)
    : constant(target)
    , lower(bop, nderiv)
    , uniform(bop)
{
    // NOP
}

constant_upper::constant_upper(
        const real_t target,
        const bsplineop &bop,
        const int nderiv)
    : constant(target)
    , upper(bop, nderiv)
    , uniform(bop)
{
    // NOP
}

constant_bulk::constant_bulk(const real_t target, bspline &b)
    : constant(target)
    , bulk(b)
    , uniform(b)
{
    // NOP
}

reference_lower::reference_lower(
        const real_t target,
        const bsplineop &bop,
        const int nderiv)
    : reference(target)
    , lower(bop, nderiv)
    , uniform(bop)
{
    // NOP
}

reference_upper::reference_upper(
        const real_t target,
        const bsplineop &bop,
        const int nderiv)
    : reference(target)
    , upper(bop, nderiv)
    , uniform(bop)
{
    // NOP
}

reference_bulk::reference_bulk(const real_t& target, bspline &b)
    : reference(target)
    , bulk(b)
    , uniform(b)
{
    // NOP
}

constant_upper_derivative::constant_upper_derivative(
        const real_t target,
        const bsplineop &bop,
        const int nderiv)
    : constant(target)
    , upper(bop, nderiv)
    , coefficient(bop, nderiv)
{
    // NOP
}

} // namespace constraint

} // namespace suzerain
