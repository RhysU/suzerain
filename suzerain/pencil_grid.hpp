/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
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
 *
 * pencil_grid.hpp: Class to manage data layout concerns for P3DFFT usage
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef PECOS_SUZERAIN_PENCILGRID
#define PECOS_SUZERAIN_PENCILGRID

#include <suzerain/exceptions.hpp>

namespace pecos
{

namespace suzerain
{

template < typename I = int >
class pencil_grid
{

public:
    typedef I dim_type;

    pencil_grid(const I proc_dims[2], const I nx, const I ny, const I nz)
    throw(domain_error);

    const dim_type pg1;
    const dim_type pg2;

    const dim_type nx;
    const dim_type ny;
    const dim_type nz;

};

template<typename I>
pencil_grid<I>::pencil_grid(
    const I proc_dims[2],
    const I nx, const I ny, const I nz)
throw(domain_error)
        : pg1(proc_dims[0]), pg2(proc_dims[1]),
        nx(nx), ny(ny), nz(nz)
{
    if (pg1 < 0) throw domain_error();

    if (pg2 < 0) throw domain_error();

    if (ny  < 0) throw domain_error();

    if (ny  < 0) throw domain_error();

    if (ny  < 0) throw domain_error();
}

}

}

#endif // PECOS_SUZERAIN_PENCILGRID
