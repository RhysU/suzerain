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
// mpl.hpp: tools related to Boost.MPL
// $Id$

#ifndef SUZERAIN_MPL_HPP
#define SUZERAIN_MPL_HPP

#include <suzerain/common.hpp>

namespace suzerain {

/**
 * Provides tools related to Boost.MPL.
 * @see <a href="http://www.boost.org/doc/libs/release/libs/mpl">Boost.MPL</a>
 * for more information.
 */
namespace mpl {

/**
 * An \ref array subclass that contains the contents of Sequence after
 * construction.  The template parameter Sequence must model MPL's <a
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
