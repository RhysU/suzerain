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
 * pencil.hpp: Class to support P3DFFT's physical and wave space data layout
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef PECOS_SUZERAIN_PENCIL
#define PECOS_SUZERAIN_PENCIL

#include <algorithm>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/shared_array.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <complex>

#include "pencil_grid.hpp"

namespace pecos
  {
  namespace suzerain
    {

    template<typename T = double, typename G = pencil_grid<> >
    class pencil
      {
      public:
        typedef typename G::dim_type dim_type;

      private:
        typedef typename boost::numeric::ublas::shallow_array_adaptor<T> adaptor_pspace_type;
        typedef typename boost::numeric::ublas::vector<T, adaptor_pspace_type> vector_pspace_type;

        typedef typename boost::numeric::ublas::shallow_array_adaptor<std::complex<T> > adaptor_wspace_type;
        typedef typename boost::numeric::ublas::vector<std::complex<T>, adaptor_wspace_type> vector_wspace_type;

      public:
        typedef typename vector_pspace_type::value_type pspace_value_type;
        typedef typename vector_pspace_type::const_reference pspace_const_reference;
        typedef typename vector_pspace_type::reference pspace_reference;
        typedef typename vector_pspace_type::const_pointer pspace_const_pointer;
        typedef typename vector_pspace_type::pointer pspace_pointer;

        typedef typename vector_wspace_type::value_type wspace_value_type;
        typedef typename vector_wspace_type::const_reference wspace_const_reference;
        typedef typename vector_wspace_type::reference wspace_reference;
        typedef typename vector_wspace_type::const_pointer wspace_const_pointer;
        typedef typename vector_wspace_type::pointer wspace_pointer;

      private:
        // Ensure design assumptions valid when instantiated
        BOOST_STATIC_ASSERT(
            2*sizeof(pspace_value_type) == sizeof(wspace_value_type));

        BOOST_STATIC_ASSERT((boost::is_same<
              typename vector_pspace_type::size_type,
              typename vector_wspace_type::size_type>::value));

        BOOST_STATIC_ASSERT((boost::is_same<
              typename vector_pspace_type::difference_type,
              typename vector_wspace_type::difference_type>::value));

      public:
        // By static assertion, same as also typename vector_wspace_type::size_type
        typedef typename vector_pspace_type::size_type size_type;

        // By static assertion, same as typename vector_wspace_type::difference_type
        typedef typename vector_pspace_type::difference_type difference_type;

      public:
        pencil(const dim_type pstart[3], const dim_type psize[3],
               const dim_type wstart[3], const dim_type wsize[3])
        throw(domain_error);

        const dim_type pstart_x;
        const dim_type pstart_y;
        const dim_type pstart_z;
        const dim_type psize_x;
        const dim_type psize_y;
        const dim_type psize_z;

        const dim_type wstart_x;
        const dim_type wstart_y;
        const dim_type wstart_z;
        const dim_type wsize_x;
        const dim_type wsize_y;
        const dim_type wsize_z;

      private:
        size_type pspace_nelem_;
        size_type wspace_nelem_;
        size_type array_nelem_;  
        boost::shared_array<pspace_value_type> array_;
        adaptor_pspace_type adaptor_pspace_;
        vector_pspace_type vector_pspace_;
        adaptor_wspace_type adaptor_wspace_;
        vector_wspace_type vector_wspace_;
      };

    template<typename T, typename G>
    pencil<T,G>::pencil(
      const dim_type pstart[3], const dim_type psize[3],
      const dim_type wstart[3], const dim_type wsize[3])
    throw(domain_error)
        : 
        pstart_x(pstart[0]), pstart_y(pstart[1]), pstart_z(pstart[2]),
        psize_x(psize[0]), psize_y(psize[1]), psize_z(psize[2]),
        wstart_x(wstart[0]), wstart_y(wstart[1]), wstart_z(wstart[2]),
        wsize_x(wsize[0]), wsize_y(wsize[1]), wsize_z(wsize[2]),
        pspace_nelem_(psize_x*psize_y*psize_z),
        wspace_nelem_(wsize_x*wsize_y*wsize_z),
        array_nelem_(std::max(pspace_nelem_, 2*wspace_nelem_)),
        array_(new pspace_value_type[array_nelem_]),
        adaptor_pspace_(pspace_nelem_, array_.get()),
        vector_pspace_(pspace_nelem_, adaptor_pspace_),
        adaptor_wspace_(wspace_nelem_, reinterpret_cast<wspace_pointer>(array_.get())),
        vector_wspace_(wspace_nelem_, adaptor_wspace_)
    {
      if (pstart_x < 0) throw domain_error();
      if (pstart_y < 0) throw domain_error();
      if (pstart_z < 0) throw domain_error();
      if (psize_x  < 0) throw domain_error();
      if (psize_y  < 0) throw domain_error();
      if (psize_z  < 0) throw domain_error();

      if (wstart_x < 0) throw domain_error();
      if (wstart_y < 0) throw domain_error();
      if (wstart_z < 0) throw domain_error();
      if (wsize_x  < 0) throw domain_error();
      if (wsize_y  < 0) throw domain_error();
      if (wsize_z  < 0) throw domain_error();
    }

  }
}

#endif // PECOS_SUZERAIN_PENCIL
