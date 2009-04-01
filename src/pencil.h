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
 * pencil.h: Class to support P3DFFT's physical and wave space data layout
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef PECOS_SUZERAIN_PENCIL
#define PECOS_SUZERAIN_PENCIL

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <complex>

namespace pecos { namespace suzerain {

template<typename T = double>
class pencil {

  private:
    typedef pencil<T> self_type;

    typedef typename boost::numeric::ublas::shallow_array_adaptor<T> adaptor_pspace;
    typedef typename boost::numeric::ublas::vector<T, adaptor_pspace> vector_pspace;

    typedef typename boost::numeric::ublas::shallow_array_adaptor<std::complex<T> > adaptor_wspace;
    typedef typename boost::numeric::ublas::vector<std::complex<T>, adaptor_wspace> vector_wspace;

  public:
    typedef typename vector_pspace::size_type pspace_size_type;
    typedef typename vector_pspace::difference_type pspace_difference_type;
    typedef typename vector_pspace::value_type pspace_value_type;
    typedef typename vector_pspace::const_reference pspace_const_reference;
    typedef typename vector_pspace::reference pspace_reference;
    typedef typename vector_pspace::const_pointer pspace_const_pointer;
    typedef typename vector_pspace::pointer pspace_pointer;

    typedef typename vector_wspace::size_type wspace_size_type;
    typedef typename vector_wspace::difference_type wspace_difference_type;
    typedef typename vector_wspace::value_type wspace_value_type;
    typedef typename vector_wspace::const_reference wspace_const_reference;
    typedef typename vector_wspace::reference wspace_reference;
    typedef typename vector_wspace::const_pointer wspace_const_pointer;
    typedef typename vector_wspace::pointer wspace_pointer;

};

} }

#endif // PECOS_SUZERAIN_PENCIL
