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
 * pencil_grid.h: Class to manage data layout concerns for P3DFFT usage
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef PECOS_SUZERAIN_PENCILGRID
#define PECOS_SUZERAIN_PENCILGRID

#include <algorithm>
#include <stdexcept>

#include "utility.h"

namespace pecos { namespace suzerain {

template<typename T = int>
class pencil_grid {

  public:
    typedef T dim_type;

    pencil_grid(
        const T proc_dims[2],
        const T nx, const T ny, const T nz,
        const T pstart[3], const T pend[3], const T psize[3],
        const T wstart[3], const T wend[3], const T wsize[3]);

    const dim_type pg1;
    const dim_type pg2;

    const dim_type nx;
    const dim_type ny;
    const dim_type nz;

    const dim_type pstart_x;
    const dim_type pstart_y;
    const dim_type pstart_z;
    const dim_type pend_x;
    const dim_type pend_y;
    const dim_type pend_z;
    const dim_type psize_x;
    const dim_type psize_y;
    const dim_type psize_z;

    const dim_type wstart_x;
    const dim_type wstart_y;
    const dim_type wstart_z;
    const dim_type wend_x;
    const dim_type wend_y;
    const dim_type wend_z;
    const dim_type wsize_x;
    const dim_type wsize_y;
    const dim_type wsize_z;

};

template<typename T>
pencil_grid<T>::pencil_grid(
    const T proc_dims[2],
    const T nx, const T ny, const T nz,
    const T pstart[3], const T pend[3], const T psize[3],
    const T wstart[3], const T wend[3], const T wsize[3]) 
  : pg1(proc_dims[0]), pg2(proc_dims[1]), 
    nx(nx), ny(ny), nz(nz),
    pstart_x(pstart[0]), pstart_y(pstart[1]), pstart_z(pstart[2]),
    pend_x(pend[0]), pend_y(pend[1]), pend_z(pend[2]),
    psize_x(psize[0]), psize_y(psize[1]), psize_z(psize[2]),
    wstart_x(wstart[0]), wstart_y(wstart[1]), wstart_z(wstart[2]),
    wend_x(wend[0]), wend_y(wend[1]), wend_z(wend[2]),
    wsize_x(wsize[0]), wsize_y(wsize[1]), wsize_z(wsize[2])
{
  if (pg1 < 0) throw domain_error();
  if (pg2 < 0) throw domain_error();

  if (ny  < 0) throw domain_error();
  if (ny  < 0) throw domain_error();
  if (ny  < 0) throw domain_error();

  if (pstart_x < 0) throw domain_error();
  if (pstart_y < 0) throw domain_error();
  if (pstart_z < 0) throw domain_error();
  if (pend_x   < 0) throw domain_error();
  if (pend_y   < 0) throw domain_error();
  if (pend_z   < 0) throw domain_error();
  if (psize_x  < 0) throw domain_error();
  if (psize_y  < 0) throw domain_error();
  if (psize_z  < 0) throw domain_error();

  if (wstart_x < 0) throw domain_error();
  if (wstart_y < 0) throw domain_error();
  if (wstart_z < 0) throw domain_error();
  if (wend_x   < 0) throw domain_error();
  if (wend_y   < 0) throw domain_error();
  if (wend_z   < 0) throw domain_error();
  if (wsize_x  < 0) throw domain_error();
  if (wsize_y  < 0) throw domain_error();
  if (wsize_z  < 0) throw domain_error();
}

} } 

#endif // PECOS_SUZERAIN_PENCILGRID
