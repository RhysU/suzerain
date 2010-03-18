//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2010 The PECOS Development Team
//
// Please see http://pecos.ices.utexas.edu for more information.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------
//
// mpl.hpp: tools related to Boost.MPL
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef __SUZERAIN_MPL_HPP
#define __SUZERAIN_MPL_HPP

#include <suzerain/common.hpp>

namespace suzerain {

/**
 * Provides tools related to Boost.MPL.
 * @see <a href="http://www.boost.org/doc/libs/release/libs/mpl">Boost.MPL</a>
 * for more information.
 */
namespace mpl {

/**
 * A boost::array subclass that contains the contents of Sequence after
 * construction.  The template parameter Sequence must model MPL's <a
 * href="http://www.boost.org/doc/libs/release/libs/mpl/doc/refmanual/forward-sequence.html">
 * Forward Sequence</a>.
 *
 * @tparam Sequence The sequence to use, which determines the
 *                  type and size of the boost::array ancestor as well as
 *                  an instance's default constructed content.
 *
 * @see <a href="http://www.boost.org/doc/libs/release/libs/array/">boost::array</a>
 *      for more information on the ancestor class.
 */
template< typename Sequence >
class sequence_array
    : public boost::array<
            typename Sequence::value_type,
            boost::mpl::size<Sequence>::type::value
      >
{
    typedef typename boost::array<
                typename Sequence::value_type,
                boost::mpl::size<Sequence>::type::value
            >::iterator iterator;

    struct copier_ {
        copier_(iterator it) : it_(it) {}
        template<typename U> void operator()(U u) { *(it_++) = u; }
        iterator it_;
    };

public:
    sequence_array() {
        boost::mpl::for_each<Sequence>(copier_(this->begin()));
    }
};

} // namespace mpl

} // namespace suzerain

#endif // __SUZERAIN_MPL_HPP
