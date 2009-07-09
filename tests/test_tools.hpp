/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2009 The PECOS Development Team
 * Based heavily on the Boost Test predicate utilities
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * test_tools.hpp: one-off test predicates and checks
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef PECOS_SUZERAIN_TEST_TOOLS_H
#define PECOS_SUZERAIN_TEST_TOOLS_H

#include <sstream>
#include <boost/test/test_tools.hpp>
#include <boost/test/floating_point_comparison.hpp>

// Like BOOST_CHECK_EQUAL_COLLECTION but uses BOOST_CHECK_CLOSE, which allows
// for a floating point tolerance on each comparison.  Have submitted a request
// to Boost Test maintainer to add an official BOOST_<level>_CLOSE_COLLECTION
// test predicate.
template<typename FPT>
bool check_close_collections(const FPT *left_begin, const FPT *left_end,
                             const FPT *right_begin, const FPT *right_end,
                             FPT percent_tolerance)
{
    const ::boost::test_tools::close_at_tolerance<FPT> is_close
        = ::boost::test_tools::close_at_tolerance<FPT>(
                ::boost::test_tools::percent_tolerance(percent_tolerance));

    int pos = 0;
    bool res = true;
    std::ostringstream msg;

    for( ;
         left_begin != left_end && right_begin != right_end;
         ++left_begin, ++right_begin, ++pos ) {

        if( !is_close(*left_begin,*right_begin) ) {
            msg << "\nMismatch to tolerance "
                << percent_tolerance << "%% at position " << pos << ": "
                << *left_begin << " != " << *right_begin;
            res = false;
        }
    }

    if( left_begin != left_end ) {
        std::size_t r_size = pos;
        while( left_begin != left_end ) {
            ++pos;
            ++left_begin;
        }

        msg << "\nCollections size mismatch: " << pos << " != " << r_size;
        res = false;
    }

    if( right_begin != right_end ) {
        std::size_t l_size = pos;
        while( right_begin != right_end ) {
            ++pos;
            ++right_begin;
        }

        msg << "\nCollections size mismatch: " << l_size << " != " << pos;
        res = false;
    }

    BOOST_CHECK_MESSAGE(res, msg.str());

    return res;
}

#endif // PECOS_SUZERAIN_TEST_TOOLS_H
