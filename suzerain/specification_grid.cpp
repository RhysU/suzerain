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

/** @file
 * @copydoc specification_grid.hpp
 */

#include <suzerain/specification_grid.hpp>

#include <suzerain/exprparse.hpp>
#include <suzerain/validation.hpp>

namespace suzerain {

specification_grid::specification_grid()
    : L(std::numeric_limits<real_t>::quiet_NaN(),
        std::numeric_limits<real_t>::quiet_NaN(),
        std::numeric_limits<real_t>::quiet_NaN()),
      N(0, 0, 0),
      DAF(std::numeric_limits<real_t>::quiet_NaN(),
          1 /* Never dealiased */,
          std::numeric_limits<real_t>::quiet_NaN()),
      dN(0, 0, 0),
      P(0, 0),
      k(0),
      htdelta(std::numeric_limits<real_t>::quiet_NaN())
{
}

specification_grid::specification_grid(real_t Lx,
                                       int    Nx,
                                       real_t DAFx,
                                       real_t Ly,
                                       int    Ny,
                                       int    k,
                                       real_t htdelta,
                                       real_t Lz,
                                       int    Nz,
                                       real_t DAFz)
    : L(Lx, Ly, Lz)
    , N(Nx, Ny, Nz)
    , DAF(DAFx, 1 /* Never dealiased */, DAFz)
    , dN(Nx * DAFx, Ny, Nz * DAFz)
    , P(0, 0)
    , k(k)
    , htdelta(htdelta)
{
}

specification_grid& specification_grid::Nx(int value)
{
    if (N.x())
    {
        validation::ensure_positive(value, "Nx");
    }
    else
    {
        validation::ensure_nonnegative(value, "Nx");
    }
    const_cast<int&>(N.x())  = value;
#pragma warning(push,disable:2259)
    const_cast<int&>(dN.x()) = N.x() * DAF.x();
#pragma warning(pop)
    return *this;
}

specification_grid& specification_grid::Ny(int value)
{
    if (N.y())
    {
        validation::ensure_positive(value, "Ny");
    }
    else
    {
        validation::ensure_nonnegative(value, "Ny");
    }
    const_cast<int&>(N.y())  = value;
#pragma warning(push,disable:2259)
    const_cast<int&>(dN.y()) = N.y() * DAF.y();
#pragma warning(pop)
    return *this;
}

specification_grid& specification_grid::Nz(int value)
{
    if (N.z())
    {
        validation::ensure_positive(value, "Nz");
    }
    else
    {
        validation::ensure_nonnegative(value, "Nz");
    }
    const_cast<int&>(N.z())  = value;
#pragma warning(push,disable:2259)
    const_cast<int&>(dN.z()) = N.z() * DAF.z();
#pragma warning(pop)
    return *this;
}

specification_grid& specification_grid::DAFx(real_t factor)
{
#pragma warning(push,disable:1572)
    if (DAF.x() != 0)
    {
#pragma warning(pop)
        validation::ensure_positive(factor, "DAFx");
    }
    else
    {
        validation::ensure_nonnegative(factor, "DAFx");
    }
    const_cast<real_t&>(DAF.x()) = factor;
#pragma warning(push,disable:2259)
    const_cast<int&   >(dN.x())  = N.x() * DAF.x();
#pragma warning(pop)
    return *this;
}

specification_grid& specification_grid::DAFz(real_t factor)
{
#pragma warning(push,disable:1572)
    if (DAF.z() != 0)
    {
#pragma warning(pop)
        validation::ensure_positive(factor, "DAFz");
    }
    else
    {
        validation::ensure_nonnegative(factor, "DAFz");
    }
    const_cast<real_t&>(DAF.z()) = factor;
#pragma warning(push,disable:2259)
    const_cast<int&   >(dN.z())  = N.z() * DAF.z();
#pragma warning(pop)
    return *this;
}

specification_grid& specification_grid::Nx(const std::string& value)
{
#pragma warning(push,disable:2259)
    return Nx(static_cast<int>(exprparse<real_t>(value, "Nx")));
#pragma warning(pop)
}

specification_grid& specification_grid::Ny(const std::string& value)
{
#pragma warning(push,disable:2259)
    return Ny(static_cast<int>(exprparse<real_t>(value, "Ny")));
#pragma warning(pop)
}

specification_grid& specification_grid::Nz(const std::string& value)
{
#pragma warning(push,disable:2259)
    return Nz(static_cast<int>(exprparse<real_t>(value, "Nz")));
#pragma warning(pop)
}

specification_grid& specification_grid::DAFx(const std::string& factor)
{
    return DAFx(exprparse<real_t>(factor, "DAFx"));
}

specification_grid& specification_grid::DAFz(const std::string& factor)
{
    return DAFz(exprparse<real_t>(factor, "DAFz"));
}

} // end namespace suzerain
