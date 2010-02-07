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

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>
#include <suzerain/common.hpp>
#include <suzerain/complex.hpp>
#include <suzerain/multi_array.hpp>

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
     * @param N Number of points in the physical domain
     * @param max_mode_exclusive Exclusive upper bound on the signal's
     *        frequency content.
     * @param shift Phase shift in the signal
     * @param length Domain size over which the signal is periodic
     * @param constant The constant content of the signal
     */
     periodic_function(const Integer N,
                       const Integer max_mode_exclusive = -1,
                       const FPT shift = M_PI/3,
                       const FPT length = 2*M_PI,
                       const FPT constant = 17)
        : N(N),
          max_mode_exclusive(
                    max_mode_exclusive >= 0 ? max_mode_exclusive : N/2+1),
          shift(shift),
          length(length),
          constant(constant)
        {
            assert(max_mode_exclusive <= (N/2+1));
        }

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
    const Integer constant;
};


template<typename FPT, typename Integer>
FPT periodic_function<FPT,Integer>::physical(
        const Integer i,
        const Integer derivative) const
{
    assert(0 <= i);
    assert(i < N);
    assert(0 <= (derivative % 4) && (derivative % 4) <= 3);

    const FPT xi = i*length/N;
    FPT retval = (max_mode_exclusive > 0 && derivative == 0 )
        ? constant : 0;
    for (Integer i = 1; i < max_mode_exclusive; ++i) {
        switch (derivative % 4) {
            case 0:
                retval +=   i * pow(i*(2.0*M_PI/length), derivative)
                          * sin(i*(2.0*M_PI/length)*xi + shift);
                break;
            case 1:
                retval +=   i * pow(i*(2.0*M_PI/length), derivative)
                          * cos(i*(2.0*M_PI/length)*xi + shift);
                break;
            case 2:
                retval -=   i * pow(i*(2.0*M_PI/length), derivative)
                          * sin(i*(2.0*M_PI/length)*xi + shift);
                break;
            case 3:
                retval -=   i * pow(i*(2.0*M_PI/length), derivative)
                          * cos(i*(2.0*M_PI/length)*xi + shift);
                break;
        }
    }
    return retval;
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
            retval = complex_type(constant, 0);
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

    for (int j = 0; j < derivative; ++j) {
        retval *= complex_type(0, 2*M_PI/length);
    }

    return retval;
}

#endif // PECOS_SUZERAIN_TEST_TOOLS_HPP
