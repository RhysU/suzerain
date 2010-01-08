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

#ifndef PECOS_SUZERAIN_TEST_TOOLS_HPP
#define PECOS_SUZERAIN_TEST_TOOLS_HPP

#include <suzerain/common.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#define CHECK_GBMATRIX_CLOSE(                                        \
            e_m, e_n, e_kl, e_ku, e, e_ld,                           \
            r_m, r_n, r_kl, r_ku, r, r_ld,                           \
            percent_tolerance)                                       \
        {                                                            \
            ::std::string errormsg(_suzerain_check_gbmatrix_close(   \
                    e_m, e_n, e_kl, e_ku, e, e_ld,                   \
                    r_m, r_n, r_kl, r_ku, r, r_ld,                   \
                    percent_tolerance));                             \
            BOOST_CHECK_MESSAGE(!errormsg.length(), errormsg);       \
        }

template<typename FPT>
std::string
_suzerain_check_gbmatrix_close(
    int e_m, int e_n, int e_kl, int e_ku, const FPT * const e, int e_ld,
    int r_m, int r_n, int r_kl, int r_ku, const FPT * const r, int r_ld,
    FPT percent_tolerance)
{
    bool checkequality = true;
    std::ostringstream errors;

    // Test sanity checks for expected values
    if (e_m < 0) {
        errors << "\nParameter e_m = " << e_m << " < 0";
        checkequality = false;
    }
    if (e_n < 0) {
        errors << "\nParameter e_n = " << e_n << " < 0";
        checkequality = false;
    }
    if (e_kl < 0) {
        errors << "\nParameter e_kl = " << e_kl << " < 0";
        checkequality = false;
    }
    if (e_ku < 0) {
        errors << "\nParameter e_ku = " << e_ku << " < 0";
        checkequality = false;
    }
    if (e_ld < 0) {
        errors << "\nParameter e_ld = " << e_ld << " < 0";
        checkequality = false;
    }
    if (e == NULL) {
        errors << "\nParameter e == NULL";
        checkequality = false;
    }

    // Test sanity checks for result values
    if (r_m < 0) {
        errors << "\nParameter r_m = " << r_m << " < 0";
        checkequality = false;
    }
    if (r_n < 0) {
        errors << "\nParameter r_n = " << r_n << " < 0";
        checkequality = false;
    }
    if (r_kl < 0) {
        errors << "\nParameter r_kl = " << r_kl << " < 0";
        checkequality = false;
    }
    if (r_ku < 0) {
        errors << "\nParameter r_ku = " << r_ku << " < 0";
        checkequality = false;
    }
    if (r_ld < 0) {
        errors << "\nParameter r_ld = " << r_ld << " < 0";
        checkequality = false;
    }
    if (r == NULL) {
        errors << "\nParameter r == NULL";
        checkequality = false;
    }

    // Check that expected and results have compatible shapes
    if (e_m != r_m) {
        errors << "\nMismatch in number of rows: " << e_m << " vs " << r_m;
        checkequality = false;
    }
    if (e_n != r_n) {
        errors << "\nMismatch in number of columns: " << e_n << " vs " << r_n;
        checkequality = false;
    }

    // Any further error messages are useless if the above tests fail
    // so short circuit the remainder of the test if any did.
    if (checkequality) {
        const boost::test_tools::close_at_tolerance<FPT> is_close
            = boost::test_tools::close_at_tolerance<FPT>(
                boost::test_tools::percent_tolerance(percent_tolerance));

        for (int j = 0; j < e_n; ++j) {
            for (int i = 0; i < e_m; ++i) {

                const bool e_in_band
                    = suzerain_gbmatrix_in_band(e_ld, e_kl, e_ku, i, j);
                const bool r_in_band
                    = suzerain_gbmatrix_in_band(r_ld, r_kl, r_ku, i, j);
                const int e_offset
                    = suzerain_gbmatrix_offset(e_ld, e_kl, e_ku, i, j);
                const int r_offset
                    = suzerain_gbmatrix_offset(r_ld, r_kl, r_ku, i, j);

                 if (e_in_band && r_in_band) {
                    const FPT r_value = r[r_offset];
                    const FPT e_value = e[e_offset];
                    if (!is_close(e_value, r_value)) {
                        errors << "\nMismatch to "
                            << percent_tolerance << "% at index ("
                            << std::setw(2) << i << ","
                            << std::setw(2) << j << "): ";
                        const std::ios_base::fmtflags flags = errors.flags();
                        const std::streamsize prec = errors.precision();
                        errors.flags(std::ios::scientific | std::ios::showpos);
                        errors.precision(std::numeric_limits<FPT>::digits10);
                        errors << e_value << " != " << r_value;
                        errors.flags(flags);
                        errors.precision(prec);
                    }
                } else if (!e_in_band && !r_in_band) {
                    // NOP
                } else if (!e_in_band) {
                    const FPT r_value = r[r_offset];
                    if (r_value != 0.0) {
                        errors << "\nNonzero outside expected's band "
                            << "at result index (" << std::setw(2) << i
                            << "," << std::setw(2) << j << "): ";
                        const std::ios_base::fmtflags flags = errors.flags();
                        const std::streamsize prec = errors.precision();
                        errors.flags(std::ios::scientific | std::ios::showpos);
                        errors.precision(std::numeric_limits<FPT>::digits10);
                        errors << r_value;
                        errors.flags(flags);
                        errors.precision(prec);
                    }
                } else if (!r_in_band) {
                    const FPT e_value = e[e_offset];
                    if (e_value != 0.0) {
                        errors << "\nNonzero outside result's band "
                            << "at expected index (" << std::setw(2) << i
                            << "," << std::setw(2) << j << "): ";
                        const std::ios_base::fmtflags flags = errors.flags();
                        const std::streamsize prec = errors.precision();
                        errors.flags(std::ios::scientific | std::ios::showpos);
                        errors.precision(std::numeric_limits<FPT>::digits10);
                        errors << e_value;
                        errors.flags(flags);
                        errors.precision(prec);
                    }
                }
            }
        }
    }

    return errors.str();
}

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
            if (abs(*right_begin) > percent_tolerance) pos_okay = false;
        } else if ( *right_begin == 0 ) {
            if (abs(*left_begin) > percent_tolerance) pos_okay = false;
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


// Provide a general operator<<(basic_ostream, boost::array) template in ::boost
// Required for many of the Boost.Test predicates to compile correctly
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
}

#endif // PECOS_SUZERAIN_TEST_TOOLS_HPP
