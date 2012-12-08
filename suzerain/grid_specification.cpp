//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// grid_specification.cpp: classes handling 3D, dealiased grid specifications
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/grid_specification.hpp>

#include <suzerain/exprparse.hpp>
#include <suzerain/validation.hpp>

namespace suzerain {

grid_specification::grid_specification()
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

grid_specification::grid_specification(real_t Lx,
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

grid_specification& grid_specification::Nx(int value)
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

grid_specification& grid_specification::Ny(int value)
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

grid_specification& grid_specification::Nz(int value)
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

grid_specification& grid_specification::DAFx(real_t factor)
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

grid_specification& grid_specification::DAFz(real_t factor)
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

grid_specification& grid_specification::Nx(const std::string& value)
{
#pragma warning(push,disable:2259)
    return Nx(static_cast<int>(exprparse<real_t>(value, "Nx")));
#pragma warning(pop)
}

grid_specification& grid_specification::Ny(const std::string& value)
{
#pragma warning(push,disable:2259)
    return Ny(static_cast<int>(exprparse<real_t>(value, "Ny")));
#pragma warning(pop)
}

grid_specification& grid_specification::Nz(const std::string& value)
{
#pragma warning(push,disable:2259)
    return Nz(static_cast<int>(exprparse<real_t>(value, "Nz")));
#pragma warning(pop)
}

grid_specification& grid_specification::DAFx(const std::string& value)
{
    return DAFx(exprparse<real_t>(value, "DAFx"));
}

grid_specification& grid_specification::DAFz(const std::string& value)
{
    return DAFz(exprparse<real_t>(value, "DAFz"));
}

} // end namespace suzerain