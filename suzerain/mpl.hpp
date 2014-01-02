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

#ifndef SUZERAIN_MPL_HPP
#define SUZERAIN_MPL_HPP

/** @file
 * Tools related to Boost.MPL.
 */

#include <suzerain/common.hpp>

namespace suzerain {

/**
 * Provides tools related to Boost.MPL.
 * @see <a href="http://www.boost.org/doc/libs/release/libs/mpl">Boost.MPL</a>
 * for more information.
 */
namespace mpl {

/**
 * An array subclass that contains the contents of Sequence after construction.
 * The template parameter Sequence must model MPL's <a
 * href="http://www.boost.org/doc/libs/release/libs/mpl/doc/refmanual/forward-sequence.html">
 * Forward Sequence</a>.
 *
 * @tparam Sequence The sequence to use, which determines the
 *                  type and size of the \ref array ancestor as well as
 *                  an instance's default constructed content.
 */
template< typename Sequence >
class sequence_array
    : public array<
            typename Sequence::value_type,
            boost::mpl::size<Sequence>::type::value
      >
{
    typedef typename array<
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

#endif // SUZERAIN_MPL_HPP
