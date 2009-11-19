/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2009 The PECOS Development Team
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
 * state.hpp: Class to manage mutable state vectors
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_STATE_H
#define __SUZERAIN_STATE_H

#include <suzerain/common.hpp>
#include <suzerain/blas_et_al.h>

namespace suzerain
{

template < typename Type >
class state : boost::noncopyable
{

public:
    typedef boost::multi_array_types::index index;
    typedef boost::multi_array_types::size_type size_type;
    typedef boost::multi_array_types::difference_type difference_type;
    typedef typename boost::multi_array_ref<Type, 3> data_type;

    explicit state(size_type variable_count,
                   size_type vector_length,
                   size_type vector_count)
    throw(std::bad_alloc);

private:

    boost::shared_array<Type> raw;  /**< Raw, aligned storage */

public:

    data_type data;                 /**< MultiArray of state data */

    const size_type variable_count; /**< Number of state variables */
    const size_type vector_length;  /**< Length of each state vector */
    const size_type vector_count;   /**< Number of state vectors */
};

template< typename Type >
state<Type>::state(size_type variable_count,
                   size_type vector_length,
                   size_type vector_count)
throw(std::bad_alloc)
    : raw(reinterpret_cast<Type *>(suzerain_blas_malloc(
                sizeof(Type)*variable_count*vector_length*vector_count)),
          std::ptr_fun(free)),
      data(raw.get(),
           boost::extents[variable_count][vector_length][vector_count],
           boost::fortran_storage_order()),
      variable_count(variable_count),
      vector_length(vector_length),
      vector_count(vector_count)
{
    if (!raw) throw std::bad_alloc();
}

} // namespace suzerain

#endif // __SUZERAIN_STATE_H
