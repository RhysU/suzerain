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
// common.hpp: C++ common definitions, utility macros, and inline functions
// $Id$

#ifndef SUZERAIN_COMMON_HPP
#define SUZERAIN_COMMON_HPP

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
// Boost.Preprocessor was included in common.h
#ifdef SUZERAIN_HAVE_BOOST
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/array.hpp>
#include <boost/bind.hpp>
#include <boost/concept/assert.hpp>
#include <boost/current_function.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/function.hpp>
#include <boost/integer_traits.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/make_shared.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/list_c.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/multi_array.hpp>
#include <boost/noncopyable.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/numeric/conversion/converter.hpp>
#include <boost/program_options.hpp>
#include <boost/ptr_container/ptr_list.hpp>
// https://svn.boost.org/trac/boost/ticket/4276
SUZERAIN_GCC_DIAG_ON(ignored-qualifiers);
#include <boost/ptr_container/ptr_map.hpp>
SUZERAIN_GCC_DIAG_OFF(ignored-qualifiers);
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

/**
 * Ensure that \c expr evaluates to boolean \c true at runtime.  If \c expr
 * evaluates to boolean \c false, then an exception \c except is thrown with
 * message \c msg.
 *
 * This macro is intended for <tt>assert</tt>-like checks which should always
 * be performed regardless of whether or not \c NDEBUG is <tt>#define</tt>d.
 */
#define SUZERAIN_ENSURE_MSGEXCEPT(expr, msg, except) \
    if (SUZERAIN_UNLIKELY(!(expr)))                  \
        throw except(::std::string(msg " in ") + BOOST_CURRENT_FUNCTION)

/**
 * Ensure that \c expr evaluates to boolean \c true at runtime.  If \c expr
 * evaluates to boolean \c false, then a <tt>std::logic_error</tt> is thrown
 * with message \c msg.
 *
 * This macro is intended for <tt>assert</tt>-like checks which should always
 * be performed regardless of whether or not \c NDEBUG is <tt>#define</tt>d.
 */
#define SUZERAIN_ENSURE_MSG(expr, msg) \
    SUZERAIN_ENSURE_MSGEXCEPT(expr, msg, ::std::logic_error)

/**
 * Ensure that \c expr evaluates to boolean \c true at runtime.  If \c expr
 * evaluates to boolean \c false, then a <tt>std::logic_error</tt> is thrown.
 *
 * This macro is intended for <tt>assert</tt>-like checks which should always
 * be performed regardless of whether or not \c NDEBUG is <tt>#define</tt>d.
 */
#define SUZERAIN_ENSURE(expr) \
    SUZERAIN_ENSURE_MSG(expr, BOOST_PP_STRINGIZE(expr) " false")

/**
 * Ensure that \c expr evaluates to boolean \c true at runtime.  If \c expr
 * evaluates to boolean \c false, then an exception \c except is thrown.
 *
 * This macro is intended for <tt>assert</tt>-like checks which should always
 * be performed regardless of whether or not \c NDEBUG is <tt>#define</tt>d.
 */
#define SUZERAIN_ENSURE_EXCEPT(expr, except) \
    SUZERAIN_ENSURE_MSGEXCEPT(expr, BOOST_PP_STRINGIZE(expr) " false", except)


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

#endif // SUZERAIN_COMMON_HPP
