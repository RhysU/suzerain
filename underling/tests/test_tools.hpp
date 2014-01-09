//-----------------------------------------------------------------------bl-
// underling 0.3.1: an FFTW MPI-based library for 3D pencil decompositions
// http://red.ices.utexas.edu/projects/underling
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
//
// This file is part of underling.
//
// underling is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// underling is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with underling.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------el-
// $Id$

#ifndef UNDERLING_TEST_TOOLS_HPP
#define UNDERLING_TEST_TOOLS_HPP

#include <complex>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#pragma warning(push,disable:1418)

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
    std::ostringstream errormsg;

    for( ;
         left_begin != left_end && right_begin != right_end;
         ++left_begin, ++right_begin, ++pos ) {

        bool pos_okay = true;

        if( *left_begin == 0 && *right_begin == 0 ) {
            // NOP on both identically zero
        } else if ( *left_begin == 0 ) {
            if (std::abs(*right_begin) > percent_tolerance) pos_okay = false;
        } else if ( *right_begin == 0 ) {
            if (std::abs(*left_begin) > percent_tolerance) pos_okay = false;
        } else if( !is_close(*left_begin,*right_begin) ) {
            pos_okay = false;
        }

        if (!pos_okay) {
            errormsg << "\nMismatch to tolerance "
                << percent_tolerance << "% at position " << pos << ": ";
            const ::std::ios_base::fmtflags flags = errormsg.flags();
            const ::std::streamsize         prec  = errormsg.precision();
            errormsg.flags(::std::ios::scientific | ::std::ios::showpos);
            errormsg.precision(::std::numeric_limits<FPT>::digits10 + 1);
            errormsg << *left_begin << " != " << *right_begin;
            errormsg.flags(flags);
            errormsg.precision(prec);
            res = false;
        }
    }

    if( left_begin != left_end ) {
        std::size_t r_size = pos;
        while( left_begin != left_end ) {
            ++pos;
            ++left_begin;
        }

        errormsg << "\nCollections size mismatch: " << pos << " != " << r_size;
        res = false;
    }

    if( right_begin != right_end ) {
        std::size_t l_size = pos;
        while( right_begin != right_end ) {
            ++pos;
            ++right_begin;
        }

        errormsg << "\nCollections size mismatch: " << l_size << " != " << pos;
        res = false;
    }

    BOOST_CHECK_MESSAGE(res, errormsg.str());

    return res;
}

// Like BOOST_CHECK_EQUAL_COLLECTION but uses BOOST_CHECK_SMALL and
// compares the absolute value of a complex difference.
template<typename FPT>
bool check_close_complex_collections(
        const FPT (*left_begin )[2], const FPT (*left_end )[2],
        const FPT (*right_begin)[2], const FPT (*right_end)[2],
        FPT abs_tolerance)
{
    int pos = 0;
    bool res = true;
    std::ostringstream errormsg;

    for( ;
         left_begin != left_end && right_begin != right_end;
         ++left_begin, ++right_begin, ++pos ) {

        const FPT diff_re  = (*left_begin)[0] - (*right_begin)[0];
        const FPT diff_im  = (*left_begin)[1] - (*right_begin)[1];
        const FPT diff_abs = std::sqrt(diff_re*diff_re + diff_im*diff_im);

        const bool pos_okay = diff_abs <= abs_tolerance;

        if (!pos_okay) {
            errormsg << "\nMismatch to abs tolerance "
                << abs_tolerance << " at position " << pos << ": ";
            const ::std::ios_base::fmtflags flags = errormsg.flags();
            const ::std::streamsize         prec  = errormsg.precision();
            errormsg.flags(::std::ios::scientific | ::std::ios::showpos);
            errormsg.precision(::std::numeric_limits<FPT>::digits10 + 1);
            errormsg <<'{'<<(*left_begin)[0] <<','<<(*left_begin)[1] <<'}'
                     << " != "
                     <<'{'<<(*right_begin)[0]<<','<<(*right_begin)[1]<<'}';
            errormsg.flags(flags);
            errormsg.precision(prec);
            res = false;
        }
    }

    if( left_begin != left_end ) {
        std::size_t r_size = pos;
        while( left_begin != left_end ) {
            ++pos;
            ++left_begin;
        }

        errormsg << "\nCollections size mismatch: " << pos << " != " << r_size;
        res = false;
    }

    if( right_begin != right_end ) {
        std::size_t l_size = pos;
        while( right_begin != right_end ) {
            ++pos;
            ++right_begin;
        }

        errormsg << "\nCollections size mismatch: " << l_size << " != " << pos;
        res = false;
    }

    BOOST_CHECK_MESSAGE(res, errormsg.str());

    return res;
}

#pragma warning(pop)
#endif // UNDERLING_TEST_TOOLS_HPP
