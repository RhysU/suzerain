/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
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
 * common.hpp: C++ common definitions, utility macros, and inline functions
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_COMMON_HPP
#define __SUZERAIN_COMMON_HPP

// Include all of the C common material
#include <suzerain/common.h>

// Required standard C++ functionality used throughout Suzerain
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <ios>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <locale>
#include <memory>
#include <new>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <utility>
#include <valarray>
#include <vector>

// Include Boost functionality used throughout Suzerain
#ifdef SUZERAIN_HAVE_BOOST
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/array.hpp>
#include <boost/bind.hpp>
#include <boost/concept/assert.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/function.hpp>
#include <boost/integer_traits.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/make_shared.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/list_c.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/multi_array.hpp>
#include <boost/noncopyable.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/numeric/conversion/converter.hpp>
#include <boost/preprocessor.hpp>
#include <boost/program_options.hpp>
#include <boost/ptr_container/ptr_list.hpp>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ref.hpp>
#include <boost/scoped_array.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/static_assert.hpp>
#include <boost/swap.hpp>
#include <boost/test/utils/nullstream.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility.hpp>

// Provide an operator<<(basic_ostream, boost::array) template in ::boost
namespace boost {
template< typename charT, typename traits, typename T, ::std::size_t N >
::std::basic_ostream<charT,traits>& operator<<(
        ::std::basic_ostream<charT,traits> &os,
        const ::boost::array<T,N> &array)
{
    os << '[' << N << "]{ ";
    ::std::copy(array.begin(),
                array.end(),
                ::std::ostream_iterator<T,charT,traits>(os, " "));
    os << '}';
    return os;
}
} // namespace boost

// SHIFTED_SUM taken from http://lists.boost.org/boost-users/2009/10/53245.php

/**
 * Helper macro used within SUZERAIN_SHIFTED_SUM.
 * @see SUZERAIN_SHIFTED_SUM for more details.
 */
#define SUZERAIN_SHIFTED_SUM_OP(s, state, seq) \
    (BOOST_PP_SEQ_PUSH_BACK(                   \
        BOOST_PP_TUPLE_ELEM(2, 0, state),      \
        (BOOST_PP_TUPLE_ELEM(2, 0, seq),       \
         BOOST_PP_TUPLE_ELEM(2, 1, state))),   \
     BOOST_PP_ADD(                             \
        BOOST_PP_TUPLE_ELEM(2, 1, state),      \
        BOOST_PP_TUPLE_ELEM(2, 1, seq)))

/**
 * A Boost.Preprocessor shifted sum operation written by Steven Watanabe which
 * operates on the second element of each tuple in a sequence
 * (http://lists.boost.org/boost-users/2009/10/53245.php).  It, for example,
 * takes the sequence of tuples <tt> ((A, 1)) ((B, 1)) ((C, 1)) ((D, 2)) ((E,
 * 4)) </tt> to the sequence of tuples <tt> ((A, 0)) ((B, 1)) ((C, 2)) ((D, 3))
 * ((E, 5)) </tt> and is handy for computing absolute offsets given a sequence
 * containing names and sizes.
 */
#define SUZERAIN_SHIFTED_SUM(seq)    \
    BOOST_PP_TUPLE_ELEM(2, 0,        \
        BOOST_PP_SEQ_FOLD_LEFT(      \
            SUZERAIN_SHIFTED_SUM_OP, \
            (BOOST_PP_SEQ_NIL, 0),   \
            seq))

#endif // SUZERAIN_HAVE_BOOST

// Include other functionality used throughout Suzerain
#ifdef SUZERAIN_HAVE_EIGEN
#include <Eigen/Core>
#endif // SUZERAIN_HAVE_EIGEN

#endif // __SUZERAIN_COMMON_HPP
