//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
// Functionality is based heavily on the Boost Test predicate utilities.
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
//--------------------------------------------------------------------------

#ifndef PECOS_SUZERAIN_TEST_TOOLS_HPP
#define PECOS_SUZERAIN_TEST_TOOLS_HPP

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include <suzerain/common.hpp>
#include <suzerain/complex.hpp>
#include <suzerain/error.h>
#include <suzerain/gbmatrix.h>
#include <suzerain/multi_array.hpp>

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#ifdef HAVE_MKL
#include <mkl.h>
#endif

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

#define CHECK_GBMATRIX_SYMMETRIC(m, n, kl, ku, r, ld)                  \
        {                                                              \
            ::std::string errormsg(_suzerain_check_gbmatrix_symmetric( \
                    m, n, kl, ku, r, ld));                             \
            BOOST_CHECK_MESSAGE(!errormsg.length(), errormsg);         \
        }

#pragma warning(push,disable:1418)

// TODO Two versions _suzerain_check_gbmatrix_close have mucho copy'n'paste

template<typename FPT>
std::string
_suzerain_check_gbmatrix_close(
    int e_m, int e_n, int e_kl, int e_ku, const FPT * const e, int e_ld,
    int r_m, int r_n, int r_kl, int r_ku, const FPT * const r, int r_ld,
    FPT percent_tolerance,
    FPT small_tolerance = -1)
{
    if (small_tolerance == -1) small_tolerance = percent_tolerance;

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
        using boost::test_tools::check_is_small;

        for (int j = 0; j < e_n; ++j) {
            for (int i = 0; i < e_m; ++i) {

                const bool e_in_band
                    = suzerain_gbmatrix_in_band(e_kl, e_ku, i, j);
                const bool r_in_band
                    = suzerain_gbmatrix_in_band(r_kl, r_ku, i, j);
                const int e_offset
                    = suzerain_gbmatrix_offset(e_ld, e_kl, e_ku, i, j);
                const int r_offset
                    = suzerain_gbmatrix_offset(r_ld, r_kl, r_ku, i, j);

                 if (e_in_band && r_in_band) {
                    const FPT r_value = r[r_offset];
                    const FPT e_value = e[e_offset];
                    if (e_value == FPT(0)) {
                        if (!check_is_small(r_value, small_tolerance)) {
                            errors << "\nMismatch of expected zero to "
                                << small_tolerance << " at index ("
                                << std::setw(2) << i << ","
                                << std::setw(2) << j << "): ";
                            const std::ios_base::fmtflags flags = errors.flags();
                            const std::streamsize prec = errors.precision();
                            errors.flags(std::ios::scientific | std::ios::showpos);
                            errors.precision(std::numeric_limits<FPT>::digits10);
                            errors << r_value;
                            errors.flags(flags);
                            errors.precision(prec);
                        }
                    } else if (!is_close(e_value, r_value)) {
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

template<typename FPT>
std::string
_suzerain_check_gbmatrix_close(
    int e_m, int e_n, int e_kl, int e_ku, const FPT (* const e)[2], int e_ld,
    int r_m, int r_n, int r_kl, int r_ku, const FPT (* const r)[2], int r_ld,
    FPT abs_tolerance)
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

        for (int j = 0; j < e_n; ++j) {
            for (int i = 0; i < e_m; ++i) {

                const bool e_in_band
                    = suzerain_gbmatrix_in_band(e_kl, e_ku, i, j);
                const bool r_in_band
                    = suzerain_gbmatrix_in_band(r_kl, r_ku, i, j);
                const int e_offset
                    = suzerain_gbmatrix_offset(e_ld, e_kl, e_ku, i, j);
                const int r_offset
                    = suzerain_gbmatrix_offset(r_ld, r_kl, r_ku, i, j);

                 if (e_in_band && r_in_band) {
                    const FPT r_value_re = r[r_offset][0];
                    const FPT r_value_im = r[r_offset][1];
                    const FPT e_value_re = e[e_offset][0];
                    const FPT e_value_im = e[e_offset][1];

                    const FPT diff_re  = r_value_re - e_value_re;
                    const FPT diff_im  = r_value_im - e_value_im;
                    const FPT diff_abs
                        = std::sqrt(diff_re*diff_re+ diff_im*diff_im);
                    const bool is_close = diff_abs <= abs_tolerance;

                    if (!is_close) {
                        errors << "\nMismatch to "
                            << abs_tolerance << "% at index ("
                            << std::setw(2) << i << ","
                            << std::setw(2) << j << "): ";
                        const std::ios_base::fmtflags flags = errors.flags();
                        const std::streamsize prec = errors.precision();
                        errors.flags(std::ios::scientific | std::ios::showpos);
                        errors.precision(std::numeric_limits<FPT>::digits10);
                        errors <<'{'<< e_value_re <<','<< e_value_im <<'}'
                               << " != "
                               <<'{'<< r_value_re <<','<< r_value_im <<'}';
                        errors.flags(flags);
                        errors.precision(prec);
                    }
                } else if (!e_in_band && !r_in_band) {
                    // NOP
                } else if (!e_in_band) {
                    const FPT r_value_re = r[r_offset][0];
                    const FPT r_value_im = r[r_offset][1];
                    if (r_value_re != 0.0 || r_value_im != 0.0) {
                        errors << "\nNonzero outside expected's band "
                            << "at result index (" << std::setw(2) << i
                            << "," << std::setw(2) << j << "): ";
                        const std::ios_base::fmtflags flags = errors.flags();
                        const std::streamsize prec = errors.precision();
                        errors.flags(std::ios::scientific | std::ios::showpos);
                        errors.precision(std::numeric_limits<FPT>::digits10);
                        errors <<'{'<< r_value_re <<','<< r_value_im <<'}';
                        errors.flags(flags);
                        errors.precision(prec);
                    }
                } else if (!r_in_band) {
                    const FPT e_value_re = e[e_offset][0];
                    const FPT e_value_im = e[e_offset][1];
                    if (e_value_re != 0.0 || e_value_im != 0.0) {
                        errors << "\nNonzero outside result's band "
                            << "at expected index (" << std::setw(2) << i
                            << "," << std::setw(2) << j << "): ";
                        const std::ios_base::fmtflags flags = errors.flags();
                        const std::streamsize prec = errors.precision();
                        errors.flags(std::ios::scientific | std::ios::showpos);
                        errors.precision(std::numeric_limits<FPT>::digits10);
                        errors <<'{'<< e_value_re <<','<< e_value_im <<'}';
                        errors.flags(flags);
                        errors.precision(prec);
                    }
                }
            }
        }
    }

    return errors.str();
}

template<typename FPT>
std::string
_suzerain_check_gbmatrix_close(
    int e_m, int e_n, int e_kl, int e_ku,
    const ::std::complex<FPT> * const e, int e_ld,
    int r_m, int r_n, int r_kl, int r_ku,
    const ::std::complex<FPT> * const r, int r_ld,
    FPT abs_tolerance)
{
    return _suzerain_check_gbmatrix_close(
            e_m, e_n, e_kl, e_ku, (const FPT (*)[2])e, e_ld,
            r_m, r_n, r_kl, r_ku, (const FPT (*)[2])r, r_ld,
            abs_tolerance);
}

template<typename FPT>
std::string
_suzerain_check_gbmatrix_symmetric(
        int m, int n, int kl, int ku, const FPT *r, int ld)
{
    bool checkequality = true;
    std::ostringstream errors;

    // Test sanity checks
    if (m < 0) {
        errors << "\nParameter m = " << m << " < 0";
        checkequality = false;
    }
    if (n < 0) {
        errors << "\nParameter n = " << n << " < 0";
        checkequality = false;
    }
    if (m != n) {
        errors << "\nParameter m != n";
        checkequality = false;
    }
    if (kl < 0) {
        errors << "\nParameter kl = " << kl << " < 0";
        checkequality = false;
    }
    if (ku < 0) {
        errors << "\nParameter ku = " << ku << " < 0";
        checkequality = false;
    }
    if (ld < 0) {
        errors << "\nParameter ld = " << ld << " < 0";
        checkequality = false;
    }
    if (r == NULL) {
        errors << "\nParameter r == NULL";
        checkequality = false;
    }

    // Any further error messages are useless if the above tests fail
    // so short circuit the remainder of the test if any did.
    if (checkequality) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < i; ++j) {
                if (   suzerain_gbmatrix_in_band(kl, ku, i, j)
                    && suzerain_gbmatrix_in_band(kl, ku, j, i)) {
                    int o_ij = suzerain_gbmatrix_offset(ld, kl, ku, i, j);
                    int o_ji = suzerain_gbmatrix_offset(ld, kl, ku, j, i);
                    double v_ij = r[o_ij];
                    double v_ji = r[o_ji];

                    if (v_ij != v_ji) {
                        errors << "\nAsymmetry in matrix at index ("
                            << std::setw(2) << i << ","
                            << std::setw(2) << j << "): ";
                        const std::ios_base::fmtflags flags = errors.flags();
                        const std::streamsize prec = errors.precision();
                        errors.flags(std::ios::scientific | std::ios::showpos);
                        errors.precision(std::numeric_limits<FPT>::digits10);
                        errors <<'{'<< v_ij <<'}'
                               << " != "
                               <<'{'<< v_ji <<'}';
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
bool check_close_collections(
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

template<typename FPT>
bool check_close_collections(const std::complex<FPT> *left_begin,
                             const std::complex<FPT> *left_end,
                             const std::complex<FPT> *right_begin,
                             const std::complex<FPT> *right_end,
                             FPT abs_tolerance)
{
    return check_close_collections((const FPT (*)[2]) left_begin,
                                   (const FPT (*)[2]) left_end,
                                   (const FPT (*)[2]) right_begin,
                                   (const FPT (*)[2]) right_end,
                                   abs_tolerance);
}

/**
 * Fill a floating point MultiArray \c x with real NaN.
 *
 * @param x MultiArray to fill.
 */
template<class MultiArray>
typename boost::disable_if<
    ::suzerain::complex::traits::is_complex<typename MultiArray::element>,
    void
>::type fill_with_NaN(MultiArray &x) {
    typedef typename MultiArray::element real_type;
    BOOST_STATIC_ASSERT(std::numeric_limits<real_type>::has_quiet_NaN);
    using namespace ::suzerain::multi_array;
    fill(x, std::numeric_limits<real_type>::quiet_NaN());
}

/**
 * Fill a complex-valued MultiArray \c x with complex NaN.
 *
 * @param x MultiArray to fill.
 */
template<class MultiArray>
typename boost::enable_if<
    ::suzerain::complex::traits::is_complex<typename MultiArray::element>,
    void
>::type fill_with_NaN(MultiArray &x) {
    using namespace ::suzerain::multi_array;
    fill(x, ::suzerain::complex::NaN<typename MultiArray::element>());
}



/** Provides a periodic function useful for testing FFT behavior */
template<typename FPT = double, typename Integer = int>
class periodic_function {
public:

    /**
     * Produce a periodic real signal with known frequency content on domain of
     * supplied length.
     *
     * @param N Number of points in the physical domain.
     * @param max_mode_exclusive Exclusive upper bound on the signal's
     *        frequency content.  Providing a value less than one is
     *        shorthand for the maximum mode supportable given \c N.
     * @param shift Phase shift in the signal
     * @param length Domain size over which the signal is periodic
     * @param constant The constant content of the signal and an amplitude
     *                 factor used to scale all modes.
     */
     explicit periodic_function(const Integer N,
                                const Integer max_mode_exclusive = 0,
                                const FPT shift = M_PI/3,
                                const FPT length = 2*M_PI,
                                const FPT constant = 17)
        : N(N),
          max_mode_exclusive(
                    max_mode_exclusive > 0 ? max_mode_exclusive : N/2+1),
          shift(shift),
          length(length),
          constant(constant)
        {
            assert(max_mode_exclusive <= (N/2+1));
        }

    /**
     * Retrieve the real-valued signal amplitude at the given location.
     *
     * @param x The desired grid location.
     *          Should be within <tt>[0,length)</tt>.
     * @param derivative Desired derivative of the signal with zero
     *        indicating the signal itself.
     *
     * @param Returns the requested value.
     */
    FPT physical_evaluate(const FPT x, const Integer derivative = 0) const;

    /**
     * Retrieve the real-valued signal amplitude at the given gridpoint.
     *
     * @param i Zero-indexed value of the desired grid location.
     *        Must be within <tt>[0,N)</tt>
     * @param derivative Desired derivative of the signal with zero
     *        indicating the signal itself.
     *
     * @param Returns the requested value.
     */
    FPT physical(const Integer i, const Integer derivative = 0) const;

    /**
     * Retrieve the requested complex-valued signal modes in wave space.
     *
     * @param i Zero-indexed mode number, where zero is the constant
     *        mode.
     * @param derivative Desired derivative of the signal with zero
     *        indicating the signal itself.
     *
     * @param Returns the requested value.
     */
    typename std::complex<FPT> wave(const Integer i,
                                    const Integer derivative = 0) const;

    const Integer N;
    const Integer max_mode_exclusive;
    const FPT shift;
    const FPT length;
    const FPT constant;
};

template<typename FPT, typename Integer>
FPT periodic_function<FPT,Integer>::physical_evaluate(
        const FPT x,
        const Integer derivative) const
{
    assert(0 <= (derivative % 4) && (derivative % 4) <= 3);

    FPT retval = (max_mode_exclusive > 0 && derivative == 0 )
        ? constant : 0;
    for (Integer j = 1; j < max_mode_exclusive; ++j) {
        switch (derivative % 4) {
        case 0:
            retval +=   j * pow(j*(2.0*M_PI/length), derivative)
                      * constant
                      * sin(j*(2.0*M_PI/length)*x + shift);
            break;
        case 1:
            retval +=   j * pow(j*(2.0*M_PI/length), derivative)
                      * constant
                      * cos(j*(2.0*M_PI/length)*x + shift);
            break;
        case 2:
            retval -=   j * pow(j*(2.0*M_PI/length), derivative)
                      * constant
                      * sin(j*(2.0*M_PI/length)*x + shift);
            break;
        case 3:
            retval -=   j * pow(j*(2.0*M_PI/length), derivative)
                      * constant
                      * cos(j*(2.0*M_PI/length)*x + shift);
            break;
        }
    }
    return retval;
}

template<typename FPT, typename Integer>
FPT periodic_function<FPT,Integer>::physical(
        const Integer i,
        const Integer derivative) const
{
    assert(0 <= i);
    assert(i <  N);

    const FPT xi = i*length/N;
    return physical_evaluate(xi, derivative);
}

template<typename FPT, typename Integer>
typename std::complex<FPT> periodic_function<FPT,Integer>::wave(
        const Integer i,
        const Integer derivative) const
{
    assert(0 <= i && i < N);

    typedef typename std::complex<FPT> complex_type;
    complex_type retval;
    if (i == 0) {
        if (i < max_mode_exclusive) {
            retval = complex_type(1, 0);
        }
    } else if (i < (N/2)) {
        if (i < max_mode_exclusive) {
            retval = complex_type(   (i)*sin(shift)/2, -  (i)*cos(shift)/2 );
        }
    } else if (i == (N/2) &&  (N%2)) {
        if (i < max_mode_exclusive) {
            retval = complex_type(   (i)*sin(shift)/2, -  (i)*cos(shift)/2 );
        }
    } else if (i == (N/2) && !(N%2)) { // Highest half mode on even grid
        if (i < max_mode_exclusive) {
            retval = complex_type(   (i)*sin(shift), 0    );
        }
    } else {
        if ((N-i) < max_mode_exclusive) {
            retval = complex_type( (N-i)*sin(shift)/2, +(N-i)*cos(shift)/2 );
        }
    }
    retval *= constant;

    for (int j = 0; j < derivative; ++j) {
        retval *= complex_type(0, 2*M_PI/length);
    }

    return retval;
}

// Relative error routine pulled from <boost/math/tools/test.hpp>.
// Could not use original directly without adding a link-time dependency.
// Slightly simplified from the source version.
template <class T>
T relative_error(T a, T b)
{
#pragma warning(push,disable:1572)
    T min_val = std::numeric_limits<T>::min();
    T max_val = std::numeric_limits<T>::max();

    if((a != 0) && (b != 0)) {
        if(std::abs(b) >= max_val) {
            if(std::abs(a) >= max_val)
                return 0;  // one infinity is as good as another!
        }
        // If the result is denormalised, treat all denorms as equivalent:
        if((a < min_val) && (a > 0))
            a = min_val;
        else if((a > -min_val) && (a < 0))
            a = -min_val;
        if((b < min_val) && (b > 0))
            b = min_val;
        else if((b > -min_val) && (b < 0))
            b = -min_val;
        return (std::max)(std::abs((a-b)/a), std::abs((a-b)/b));
    }

    // Handle special case where one or both are zero:
    if(min_val == 0)
        return std::abs(a-b);
    if(std::abs(a) < min_val)
        a = min_val;
    if(std::abs(b) < min_val)
        b = min_val;
    return (std::max)(std::abs((a-b)/a), std::abs((a-b)/b));
#pragma warning(pop)
}

/** A fixture for the Boost.Test that replaces suzerain_error */
#pragma warning(push,disable:2017 2021)
class BoostFailErrorHandlerFixture {
public:
    /** A suzerain_error_handler_t that invokes BOOST_FAIL */
    static void boost_fail_error_handler(
            const char *reason, const char *file, int line, int suzerain_errno)
    {
        std::ostringstream oss;
        oss << "Encountered '"
            << suzerain_strerror(suzerain_errno)
            << "' at "
            << file
            << ':'
            << line
            << " with reason '"
            << reason
            << "'";
        BOOST_FAIL(oss.str());
    }

    BoostFailErrorHandlerFixture()
        : previous_(suzerain_set_error_handler(&boost_fail_error_handler)) {}

    ~BoostFailErrorHandlerFixture() {
        suzerain_set_error_handler(previous_);
    }
private:
    suzerain_error_handler_t * previous_;
};
#pragma warning(pop)

#ifdef SUZERAIN_HAVE_MKL
#include <mkl_service.h>
#endif
class BlasCleanupFixture {
public:
    BlasCleanupFixture() {
#ifdef SUZERAIN_HAVE_MKL
#if INTEL_MKL_VERSION < 110002
        MKL_FreeBuffers();
#else
        mkl_free_buffers();
#endif
#endif
    }
};

#pragma warning(pop)
#endif // PECOS_SUZERAIN_TEST_TOOLS_HPP
