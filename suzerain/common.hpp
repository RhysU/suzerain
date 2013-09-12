//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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

#ifndef SUZERAIN_COMMON_HPP
#define SUZERAIN_COMMON_HPP

/** @file
 * Common includes, definitions, utility macros, and inline functions for C++.
 */

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

// If possible, begin suppressing warnings from Eigen and Boost headers
SUZERAIN_GCC_DIAG_OFF(ignored-qualifiers);
SUZERAIN_GCC_DIAG_OFF(suggest-attribute=const);
SUZERAIN_GCC_DIAG_OFF(suggest-attribute=noreturn);
SUZERAIN_GCC_DIAG_OFF(suggest-attribute=pure);
SUZERAIN_GCC_DIAG_OFF(unused-parameter);
SUZERAIN_GCC_DIAG_OFF(unused-variable);

// Include Eigen functionality used throughout Suzerain
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>

// Include Boost functionality used throughout Suzerain
// Boost.Preprocessor was included in common.h
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
#include <boost/math/special_functions.hpp>
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
#include <boost/typeof/typeof.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility.hpp>

// If possible, stop suppressing warnings from Eigen and Boost headers
SUZERAIN_GCC_DIAG_ON(unused-variable);
SUZERAIN_GCC_DIAG_ON(unused-parameter);
SUZERAIN_GCC_DIAG_ON(suggest-attribute=pure);
SUZERAIN_GCC_DIAG_ON(suggest-attribute=noreturn);
SUZERAIN_GCC_DIAG_ON(suggest-attribute=const);
SUZERAIN_GCC_DIAG_ON(ignored-qualifiers);

/**
 * Ensure that \c expr evaluates to boolean \c true at runtime.  If \c expr
 * evaluates to boolean \c false, then an exception \c except is thrown with
 * message \c msg.
 *
 * This macro is intended for <tt>assert</tt>-like checks which should always
 * be performed regardless of whether or not \c NDEBUG is <tt>define</tt>d.
 */
#define SUZERAIN_ENSURE_MSGEXCEPT(expr, msg, except) \
    if (SUZERAIN_UNLIKELY(!(expr)))                  \
        throw except(::std::string(msg " in ") + BOOST_CURRENT_FUNCTION)

/**
 * Ensure that \c expr evaluates to boolean \c true at runtime.  If \c expr
 * evaluates to boolean \c false, then a <tt>std::runtime_error</tt> is thrown
 * with message \c msg.
 *
 * This macro is intended for <tt>assert</tt>-like checks which should always
 * be performed regardless of whether or not \c NDEBUG is <tt>define</tt>d.
 */
#define SUZERAIN_ENSURE_MSG(expr, msg) \
    SUZERAIN_ENSURE_MSGEXCEPT(expr, msg, ::std::runtime_error)

/**
 * Ensure that \c expr evaluates to boolean \c true at runtime.  If \c expr
 * evaluates to boolean \c false, then a <tt>std::runtime_error</tt> is thrown.
 *
 * This macro is intended for <tt>assert</tt>-like checks which should always
 * be performed regardless of whether or not \c NDEBUG is <tt>define</tt>d.
 */
#define SUZERAIN_ENSURE(expr) \
    SUZERAIN_ENSURE_MSG(expr, BOOST_PP_STRINGIZE(expr) " false")

/**
 * Ensure that \c expr evaluates to boolean \c true at runtime.  If \c expr
 * evaluates to boolean \c false, then an exception \c except is thrown.
 *
 * This macro is intended for <tt>assert</tt>-like checks which should always
 * be performed regardless of whether or not \c NDEBUG is <tt>define</tt>d.
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

// Provide an operator<<(basic_ostream, boost::array) template in namespace boost
// http://agentzlerich.blogspot.com/2009/12/small-gotcha-when-combining-boostarray.html
namespace boost {
template< typename CharT, typename Traits, typename T, ::std::size_t N >
::std::basic_ostream<CharT,Traits>& operator<<(
        ::std::basic_ostream<CharT,Traits> &os,
        const ::boost::array<T,N> &array)
{
    os << '[' << N << "]{ ";
    ::std::copy(array.begin(), array.end(),
                ::std::ostream_iterator<T,CharT,Traits>(os, " "));
    os << '}';
    return os;
}
}

//////////////////////////////////////////////////////////
// Miscellaneous functionality used throughout Suzerain //
//////////////////////////////////////////////////////////

/** Topmost namespace for Suzerain: a spectral, compressible DNS framework. */
namespace suzerain {

/**
 * Typedefs for real and complex-valued scalars.
 * @{
 */

/**
 * The default real-valued scalar type within Suzerain.  Currently only
 * <tt>real_t == double</tt> is supported by many, many components.
 */
typedef double real_t;

/**
 * The default complex-valued scalar type within Suzerain based on \ref real_t.
 */
typedef std::complex<real_t> complex_t;

/** @} */

// Bring ubiquitous smart pointer definitions into this namespace
// Done because these may be part of public APIs, as well as for convenience.
using boost::make_shared;   /**< \namespace suzerain */
using boost::scoped_array;  /**< \namespace suzerain */
using boost::scoped_ptr;    /**< \namespace suzerain */
using boost::shared_array;  /**< \namespace suzerain */
using boost::shared_ptr;    /**< \namespace suzerain */

// Likewise, bring boost::array into this namespace
using boost::array;         /**< \namespace suzerain */

/** \weakgroup EigenTypedefs Typedefs and declarations to simplify Eigen usage.
 *
 * Defines typedefs for <a href="http://eigen.tuxfamily.org/">Eigen</a>-related
 * classes within the \ref suzerain namespace.  Typedefs like Matrix3r and
 * Vector3c permit simplifying using \ref suzerain::real_t and \ref
 * suzerain::complex_t with Eigen's Matrix and Array templates.  They
 * also permit codifying Suzerain's expectations regarding column- versus
 * row-major storage regardless of the value of \c EIGEN_DEFAULT_TO_ROW_MAJOR.
 *
 * @{
 */

// Bring common Eigen names into this namespace
using Eigen::Aligned;        /**< \namespace suzerain */
using Eigen::Array;          /**< \namespace suzerain */
using Eigen::ColMajor;       /**< \namespace suzerain */
using Eigen::DenseBase;      /**< \namespace suzerain */
using Eigen::DiagonalMatrix; /**< \namespace suzerain */
using Eigen::Dynamic;        /**< \namespace suzerain */
using Eigen::InnerStride;    /**< \namespace suzerain */
using Eigen::Map;            /**< \namespace suzerain */
using Eigen::Matrix;         /**< \namespace suzerain */
using Eigen::NoChange;       /**< \namespace suzerain */
using Eigen::OuterStride;    /**< \namespace suzerain */
using Eigen::RowMajor;       /**< \namespace suzerain */
using Eigen::Sequential;     /**< \namespace suzerain */
using Eigen::Stride;         /**< \namespace suzerain */
using Eigen::Unaligned;      /**< \namespace suzerain */

/**
 * Typedefs for \ref real_t, column-major <a
 * href="http://eigen.tuxfamily.org/">Eigen</a> matrices.  The names have been
 * extended to handle fixed sizes up to 15 using hexadecimal notation.
 * @{
 */
typedef Matrix<    real_t,       1,       1, ColMajor> Matrix1r ;
typedef Matrix<    real_t,       2,       2, ColMajor> Matrix2r ;
typedef Matrix<    real_t,       3,       3, ColMajor> Matrix3r ;
typedef Matrix<    real_t,       4,       4, ColMajor> Matrix4r ;
typedef Matrix<    real_t,       5,       5, ColMajor> Matrix5r ;
typedef Matrix<    real_t,       6,       6, ColMajor> Matrix6r ;
typedef Matrix<    real_t,       7,       7, ColMajor> Matrix7r ;
typedef Matrix<    real_t,       8,       8, ColMajor> Matrix8r ;
typedef Matrix<    real_t,       9,       9, ColMajor> Matrix9r ;
typedef Matrix<    real_t,      10,      10, ColMajor> MatrixAr ;
typedef Matrix<    real_t,      11,      11, ColMajor> MatrixBr ;
typedef Matrix<    real_t,      12,      12, ColMajor> MatrixCr ;
typedef Matrix<    real_t,      13,      13, ColMajor> MatrixDr ;
typedef Matrix<    real_t,      14,      14, ColMajor> MatrixEr ;
typedef Matrix<    real_t,      15,      15, ColMajor> MatrixFr ;
typedef Matrix<    real_t,       1, Dynamic, ColMajor> Matrix1Xr;
typedef Matrix<    real_t,       2, Dynamic, ColMajor> Matrix2Xr;
typedef Matrix<    real_t,       3, Dynamic, ColMajor> Matrix3Xr;
typedef Matrix<    real_t,       4, Dynamic, ColMajor> Matrix4Xr;
typedef Matrix<    real_t,       5, Dynamic, ColMajor> Matrix5Xr;
typedef Matrix<    real_t,       6, Dynamic, ColMajor> Matrix6Xr;
typedef Matrix<    real_t,       7, Dynamic, ColMajor> Matrix7Xr;
typedef Matrix<    real_t,       8, Dynamic, ColMajor> Matrix8Xr;
typedef Matrix<    real_t,       9, Dynamic, ColMajor> Matrix9Xr;
typedef Matrix<    real_t,      10, Dynamic, ColMajor> MatrixAXr;
typedef Matrix<    real_t,      11, Dynamic, ColMajor> MatrixBXr;
typedef Matrix<    real_t,      12, Dynamic, ColMajor> MatrixCXr;
typedef Matrix<    real_t,      13, Dynamic, ColMajor> MatrixDXr;
typedef Matrix<    real_t,      14, Dynamic, ColMajor> MatrixEXr;
typedef Matrix<    real_t,      15, Dynamic, ColMajor> MatrixFXr;
typedef Matrix<    real_t, Dynamic,       1, ColMajor> MatrixX1r;
typedef Matrix<    real_t, Dynamic,       2, ColMajor> MatrixX2r;
typedef Matrix<    real_t, Dynamic,       3, ColMajor> MatrixX3r;
typedef Matrix<    real_t, Dynamic,       4, ColMajor> MatrixX4r;
typedef Matrix<    real_t, Dynamic,       5, ColMajor> MatrixX5r;
typedef Matrix<    real_t, Dynamic,       6, ColMajor> MatrixX6r;
typedef Matrix<    real_t, Dynamic,       7, ColMajor> MatrixX7r;
typedef Matrix<    real_t, Dynamic,       8, ColMajor> MatrixX8r;
typedef Matrix<    real_t, Dynamic,       9, ColMajor> MatrixX9r;
typedef Matrix<    real_t, Dynamic,      10, ColMajor> MatrixXAr;
typedef Matrix<    real_t, Dynamic,      11, ColMajor> MatrixXBr;
typedef Matrix<    real_t, Dynamic,      12, ColMajor> MatrixXCr;
typedef Matrix<    real_t, Dynamic,      13, ColMajor> MatrixXDr;
typedef Matrix<    real_t, Dynamic,      14, ColMajor> MatrixXEr;
typedef Matrix<    real_t, Dynamic,      15, ColMajor> MatrixXFr;
typedef Matrix<    real_t, Dynamic, Dynamic, ColMajor> MatrixXXr;
/** @} */

/**
 * Typedefs for \ref complex_t, column-major <a
 * href="http://eigen.tuxfamily.org/">Eigen</a> matrices.  The names have been
 * extended to handle fixed sizes up to 15 using hexadecimal notation.
 * @{
 */
typedef Matrix< complex_t,       1,       1, ColMajor> Matrix1c ;
typedef Matrix< complex_t,       2,       2, ColMajor> Matrix2c ;
typedef Matrix< complex_t,       3,       3, ColMajor> Matrix3c ;
typedef Matrix< complex_t,       4,       4, ColMajor> Matrix4c ;
typedef Matrix< complex_t,       5,       5, ColMajor> Matrix5c ;
typedef Matrix< complex_t,       6,       6, ColMajor> Matrix6c ;
typedef Matrix< complex_t,       7,       7, ColMajor> Matrix7c ;
typedef Matrix< complex_t,       8,       8, ColMajor> Matrix8c ;
typedef Matrix< complex_t,       9,       9, ColMajor> Matrix9c ;
typedef Matrix< complex_t,      10,      10, ColMajor> MatrixAc ;
typedef Matrix< complex_t,      11,      11, ColMajor> MatrixBc ;
typedef Matrix< complex_t,      12,      12, ColMajor> MatrixCc ;
typedef Matrix< complex_t,      13,      13, ColMajor> MatrixDc ;
typedef Matrix< complex_t,      14,      14, ColMajor> MatrixEc ;
typedef Matrix< complex_t,      15,      15, ColMajor> MatrixFc ;
typedef Matrix< complex_t,       1, Dynamic, ColMajor> Matrix1Xc;
typedef Matrix< complex_t,       2, Dynamic, ColMajor> Matrix2Xc;
typedef Matrix< complex_t,       3, Dynamic, ColMajor> Matrix3Xc;
typedef Matrix< complex_t,       4, Dynamic, ColMajor> Matrix4Xc;
typedef Matrix< complex_t,       5, Dynamic, ColMajor> Matrix5Xc;
typedef Matrix< complex_t,       6, Dynamic, ColMajor> Matrix6Xc;
typedef Matrix< complex_t,       7, Dynamic, ColMajor> Matrix7Xc;
typedef Matrix< complex_t,       8, Dynamic, ColMajor> Matrix8Xc;
typedef Matrix< complex_t,       9, Dynamic, ColMajor> Matrix9Xc;
typedef Matrix< complex_t,      10, Dynamic, ColMajor> MatrixAXc;
typedef Matrix< complex_t,      11, Dynamic, ColMajor> MatrixBXc;
typedef Matrix< complex_t,      12, Dynamic, ColMajor> MatrixCXc;
typedef Matrix< complex_t,      13, Dynamic, ColMajor> MatrixDXc;
typedef Matrix< complex_t,      14, Dynamic, ColMajor> MatrixEXc;
typedef Matrix< complex_t,      15, Dynamic, ColMajor> MatrixFXc;
typedef Matrix< complex_t, Dynamic,       1, ColMajor> MatrixX1c;
typedef Matrix< complex_t, Dynamic,       2, ColMajor> MatrixX2c;
typedef Matrix< complex_t, Dynamic,       3, ColMajor> MatrixX3c;
typedef Matrix< complex_t, Dynamic,       4, ColMajor> MatrixX4c;
typedef Matrix< complex_t, Dynamic,       5, ColMajor> MatrixX5c;
typedef Matrix< complex_t, Dynamic,       6, ColMajor> MatrixX6c;
typedef Matrix< complex_t, Dynamic,       7, ColMajor> MatrixX7c;
typedef Matrix< complex_t, Dynamic,       8, ColMajor> MatrixX8c;
typedef Matrix< complex_t, Dynamic,       9, ColMajor> MatrixX9c;
typedef Matrix< complex_t, Dynamic,      10, ColMajor> MatrixXAc;
typedef Matrix< complex_t, Dynamic,      11, ColMajor> MatrixXBc;
typedef Matrix< complex_t, Dynamic,      12, ColMajor> MatrixXCc;
typedef Matrix< complex_t, Dynamic,      13, ColMajor> MatrixXDc;
typedef Matrix< complex_t, Dynamic,      14, ColMajor> MatrixXEc;
typedef Matrix< complex_t, Dynamic,      15, ColMajor> MatrixXFc;
typedef Matrix< complex_t, Dynamic, Dynamic, ColMajor> MatrixXXc;
/** @} */

/**
 * Typedefs for \ref real_t, column-major <a
 * href="http://eigen.tuxfamily.org/">Eigen</a> matrices.  The names have been
 * extended to handle fixed sizes up to 15 using hexadecimal notation.
 * @{
 */
typedef Matrix<       int,       1,       1, ColMajor> Matrix1i ;
typedef Matrix<       int,       2,       2, ColMajor> Matrix2i ;
typedef Matrix<       int,       3,       3, ColMajor> Matrix3i ;
typedef Matrix<       int,       4,       4, ColMajor> Matrix4i ;
typedef Matrix<       int,       5,       5, ColMajor> Matrix5i ;
typedef Matrix<       int,       6,       6, ColMajor> Matrix6i ;
typedef Matrix<       int,       7,       7, ColMajor> Matrix7i ;
typedef Matrix<       int,       8,       8, ColMajor> Matrix8i ;
typedef Matrix<       int,       9,       9, ColMajor> Matrix9i ;
typedef Matrix<       int,      10,      10, ColMajor> MatrixAi ;
typedef Matrix<       int,      11,      11, ColMajor> MatrixBi ;
typedef Matrix<       int,      12,      12, ColMajor> MatrixCi ;
typedef Matrix<       int,      13,      13, ColMajor> MatrixDi ;
typedef Matrix<       int,      14,      14, ColMajor> MatrixEi ;
typedef Matrix<       int,      15,      15, ColMajor> MatrixFi ;
typedef Matrix<       int,       1, Dynamic, ColMajor> Matrix1Xi;
typedef Matrix<       int,       2, Dynamic, ColMajor> Matrix2Xi;
typedef Matrix<       int,       3, Dynamic, ColMajor> Matrix3Xi;
typedef Matrix<       int,       4, Dynamic, ColMajor> Matrix4Xi;
typedef Matrix<       int,       5, Dynamic, ColMajor> Matrix5Xi;
typedef Matrix<       int,       6, Dynamic, ColMajor> Matrix6Xi;
typedef Matrix<       int,       7, Dynamic, ColMajor> Matrix7Xi;
typedef Matrix<       int,       8, Dynamic, ColMajor> Matrix8Xi;
typedef Matrix<       int,       9, Dynamic, ColMajor> Matrix9Xi;
typedef Matrix<       int,      10, Dynamic, ColMajor> MatrixAXi;
typedef Matrix<       int,      11, Dynamic, ColMajor> MatrixBXi;
typedef Matrix<       int,      12, Dynamic, ColMajor> MatrixCXi;
typedef Matrix<       int,      13, Dynamic, ColMajor> MatrixDXi;
typedef Matrix<       int,      14, Dynamic, ColMajor> MatrixEXi;
typedef Matrix<       int,      15, Dynamic, ColMajor> MatrixFXi;
typedef Matrix<       int, Dynamic,       1, ColMajor> MatrixX1i;
typedef Matrix<       int, Dynamic,       2, ColMajor> MatrixX2i;
typedef Matrix<       int, Dynamic,       3, ColMajor> MatrixX3i;
typedef Matrix<       int, Dynamic,       4, ColMajor> MatrixX4i;
typedef Matrix<       int, Dynamic,       5, ColMajor> MatrixX5i;
typedef Matrix<       int, Dynamic,       6, ColMajor> MatrixX6i;
typedef Matrix<       int, Dynamic,       7, ColMajor> MatrixX7i;
typedef Matrix<       int, Dynamic,       8, ColMajor> MatrixX8i;
typedef Matrix<       int, Dynamic,       9, ColMajor> MatrixX9i;
typedef Matrix<       int, Dynamic,      10, ColMajor> MatrixXAi;
typedef Matrix<       int, Dynamic,      11, ColMajor> MatrixXBi;
typedef Matrix<       int, Dynamic,      12, ColMajor> MatrixXCi;
typedef Matrix<       int, Dynamic,      13, ColMajor> MatrixXDi;
typedef Matrix<       int, Dynamic,      14, ColMajor> MatrixXEi;
typedef Matrix<       int, Dynamic,      15, ColMajor> MatrixXFi;
typedef Matrix<       int, Dynamic, Dynamic, ColMajor> MatrixXXi;
/** @} */

/**
 * Typedefs for \ref real_t, column-major <a
 * href="http://eigen.tuxfamily.org/">Eigen</a> arrays.  The names have been
 * extended to handle fixed sizes up to 15 using hexadecimal notation.
 * @{
 */
typedef Array<     real_t,       1,       1, ColMajor> Array1r ;
typedef Array<     real_t,       2,       1, ColMajor> Array2r ;
typedef Array<     real_t,       3,       1, ColMajor> Array3r ;
typedef Array<     real_t,       4,       1, ColMajor> Array4r ;
typedef Array<     real_t,       5,       1, ColMajor> Array5r ;
typedef Array<     real_t,       6,       1, ColMajor> Array6r ;
typedef Array<     real_t,       7,       1, ColMajor> Array7r ;
typedef Array<     real_t,       8,       1, ColMajor> Array8r ;
typedef Array<     real_t,       9,       1, ColMajor> Array9r ;
typedef Array<     real_t,      10,       1, ColMajor> ArrayAr ;
typedef Array<     real_t,      11,       1, ColMajor> ArrayBr ;
typedef Array<     real_t,      12,       1, ColMajor> ArrayCr ;
typedef Array<     real_t,      13,       1, ColMajor> ArrayDr ;
typedef Array<     real_t,      14,       1, ColMajor> ArrayEr ;
typedef Array<     real_t,      15,       1, ColMajor> ArrayFr ;
typedef Array<     real_t, Dynamic,       1, ColMajor> ArrayXr ;
typedef Array<     real_t,       1, Dynamic, ColMajor> Array1Xr;
typedef Array<     real_t,       2, Dynamic, ColMajor> Array2Xr;
typedef Array<     real_t,       3, Dynamic, ColMajor> Array3Xr;
typedef Array<     real_t,       4, Dynamic, ColMajor> Array4Xr;
typedef Array<     real_t,       5, Dynamic, ColMajor> Array5Xr;
typedef Array<     real_t,       6, Dynamic, ColMajor> Array6Xr;
typedef Array<     real_t,       7, Dynamic, ColMajor> Array7Xr;
typedef Array<     real_t,       8, Dynamic, ColMajor> Array8Xr;
typedef Array<     real_t,       9, Dynamic, ColMajor> Array9Xr;
typedef Array<     real_t,      10, Dynamic, ColMajor> ArrayAXr;
typedef Array<     real_t,      11, Dynamic, ColMajor> ArrayBXr;
typedef Array<     real_t,      12, Dynamic, ColMajor> ArrayCXr;
typedef Array<     real_t,      13, Dynamic, ColMajor> ArrayDXr;
typedef Array<     real_t,      14, Dynamic, ColMajor> ArrayEXr;
typedef Array<     real_t,      15, Dynamic, ColMajor> ArrayFXr;
typedef Array<     real_t, Dynamic,       1, ColMajor> ArrayX1r;
typedef Array<     real_t, Dynamic,       2, ColMajor> ArrayX2r;
typedef Array<     real_t, Dynamic,       3, ColMajor> ArrayX3r;
typedef Array<     real_t, Dynamic,       4, ColMajor> ArrayX4r;
typedef Array<     real_t, Dynamic,       5, ColMajor> ArrayX5r;
typedef Array<     real_t, Dynamic,       6, ColMajor> ArrayX6r;
typedef Array<     real_t, Dynamic,       7, ColMajor> ArrayX7r;
typedef Array<     real_t, Dynamic,       8, ColMajor> ArrayX8r;
typedef Array<     real_t, Dynamic,       9, ColMajor> ArrayX9r;
typedef Array<     real_t, Dynamic,      10, ColMajor> ArrayXAr;
typedef Array<     real_t, Dynamic,      11, ColMajor> ArrayXBr;
typedef Array<     real_t, Dynamic,      12, ColMajor> ArrayXCr;
typedef Array<     real_t, Dynamic,      13, ColMajor> ArrayXDr;
typedef Array<     real_t, Dynamic,      14, ColMajor> ArrayXEr;
typedef Array<     real_t, Dynamic,      15, ColMajor> ArrayXFr;
typedef Array<     real_t, Dynamic, Dynamic, ColMajor> ArrayXXr;
/** @} */

/**
 * Typedefs for \ref complex_t, column-major <a
 * href="http://eigen.tuxfamily.org/">Eigen</a> arrays.  The names have been
 * extended to handle fixed sizes up to 15 using hexadecimal notation.
 * @{
 */
typedef Array<  complex_t,       1,       1, ColMajor> Array1c ;
typedef Array<  complex_t,       2,       1, ColMajor> Array2c ;
typedef Array<  complex_t,       3,       1, ColMajor> Array3c ;
typedef Array<  complex_t,       4,       1, ColMajor> Array4c ;
typedef Array<  complex_t,       5,       1, ColMajor> Array5c ;
typedef Array<  complex_t,       6,       1, ColMajor> Array6c ;
typedef Array<  complex_t,       7,       1, ColMajor> Array7c ;
typedef Array<  complex_t,       8,       1, ColMajor> Array8c ;
typedef Array<  complex_t,       9,       1, ColMajor> Array9c ;
typedef Array<  complex_t,      10,       1, ColMajor> ArrayAc ;
typedef Array<  complex_t,      11,       1, ColMajor> ArrayBc ;
typedef Array<  complex_t,      12,       1, ColMajor> ArrayCc ;
typedef Array<  complex_t,      13,       1, ColMajor> ArrayDc ;
typedef Array<  complex_t,      14,       1, ColMajor> ArrayEc ;
typedef Array<  complex_t,      15,       1, ColMajor> ArrayFc ;
typedef Array<  complex_t, Dynamic,       1, ColMajor> ArrayXc ;
typedef Array<  complex_t,       1, Dynamic, ColMajor> Array1Xc;
typedef Array<  complex_t,       2, Dynamic, ColMajor> Array2Xc;
typedef Array<  complex_t,       3, Dynamic, ColMajor> Array3Xc;
typedef Array<  complex_t,       4, Dynamic, ColMajor> Array4Xc;
typedef Array<  complex_t,       5, Dynamic, ColMajor> Array5Xc;
typedef Array<  complex_t,       6, Dynamic, ColMajor> Array6Xc;
typedef Array<  complex_t,       7, Dynamic, ColMajor> Array7Xc;
typedef Array<  complex_t,       8, Dynamic, ColMajor> Array8Xc;
typedef Array<  complex_t,       9, Dynamic, ColMajor> Array9Xc;
typedef Array<  complex_t,      10, Dynamic, ColMajor> ArrayAXc;
typedef Array<  complex_t,      11, Dynamic, ColMajor> ArrayBXc;
typedef Array<  complex_t,      12, Dynamic, ColMajor> ArrayCXc;
typedef Array<  complex_t,      13, Dynamic, ColMajor> ArrayDXc;
typedef Array<  complex_t,      14, Dynamic, ColMajor> ArrayEXc;
typedef Array<  complex_t,      15, Dynamic, ColMajor> ArrayFXc;
typedef Array<  complex_t, Dynamic,       1, ColMajor> ArrayX1c;
typedef Array<  complex_t, Dynamic,       2, ColMajor> ArrayX2c;
typedef Array<  complex_t, Dynamic,       3, ColMajor> ArrayX3c;
typedef Array<  complex_t, Dynamic,       4, ColMajor> ArrayX4c;
typedef Array<  complex_t, Dynamic,       5, ColMajor> ArrayX5c;
typedef Array<  complex_t, Dynamic,       6, ColMajor> ArrayX6c;
typedef Array<  complex_t, Dynamic,       7, ColMajor> ArrayX7c;
typedef Array<  complex_t, Dynamic,       8, ColMajor> ArrayX8c;
typedef Array<  complex_t, Dynamic,       9, ColMajor> ArrayX9c;
typedef Array<  complex_t, Dynamic,      10, ColMajor> ArrayXAc;
typedef Array<  complex_t, Dynamic,      11, ColMajor> ArrayXBc;
typedef Array<  complex_t, Dynamic,      12, ColMajor> ArrayXCc;
typedef Array<  complex_t, Dynamic,      13, ColMajor> ArrayXDc;
typedef Array<  complex_t, Dynamic,      14, ColMajor> ArrayXEc;
typedef Array<  complex_t, Dynamic,      15, ColMajor> ArrayXFc;
typedef Array<  complex_t, Dynamic, Dynamic, ColMajor> ArrayXXc;
/** @} */

/**
 * Typedefs for \c int, column-major <a
 * href="http://eigen.tuxfamily.org/">Eigen</a> arrays.  The names have been
 * extended to handle fixed sizes up to 15 using hexadecimal notation.
 * @{
 */
typedef Array<        int,       1,       1, ColMajor> Array1i ;
typedef Array<        int,       2,       1, ColMajor> Array2i ;
typedef Array<        int,       3,       1, ColMajor> Array3i ;
typedef Array<        int,       4,       1, ColMajor> Array4i ;
typedef Array<        int,       5,       1, ColMajor> Array5i ;
typedef Array<        int,       6,       1, ColMajor> Array6i ;
typedef Array<        int,       7,       1, ColMajor> Array7i ;
typedef Array<        int,       8,       1, ColMajor> Array8i ;
typedef Array<        int,       9,       1, ColMajor> Array9i ;
typedef Array<        int,      10,       1, ColMajor> ArrayAi ;
typedef Array<        int,      11,       1, ColMajor> ArrayBi ;
typedef Array<        int,      12,       1, ColMajor> ArrayCi ;
typedef Array<        int,      13,       1, ColMajor> ArrayDi ;
typedef Array<        int,      14,       1, ColMajor> ArrayEi ;
typedef Array<        int,      15,       1, ColMajor> ArrayFi ;
typedef Array<        int, Dynamic,       1, ColMajor> ArrayXi ;
typedef Array<        int,       1, Dynamic, ColMajor> Array1Xi;
typedef Array<        int,       2, Dynamic, ColMajor> Array2Xi;
typedef Array<        int,       3, Dynamic, ColMajor> Array3Xi;
typedef Array<        int,       4, Dynamic, ColMajor> Array4Xi;
typedef Array<        int,       5, Dynamic, ColMajor> Array5Xi;
typedef Array<        int,       6, Dynamic, ColMajor> Array6Xi;
typedef Array<        int,       7, Dynamic, ColMajor> Array7Xi;
typedef Array<        int,       8, Dynamic, ColMajor> Array8Xi;
typedef Array<        int,       9, Dynamic, ColMajor> Array9Xi;
typedef Array<        int,      10, Dynamic, ColMajor> ArrayAXi;
typedef Array<        int,      11, Dynamic, ColMajor> ArrayBXi;
typedef Array<        int,      12, Dynamic, ColMajor> ArrayCXi;
typedef Array<        int,      13, Dynamic, ColMajor> ArrayDXi;
typedef Array<        int,      14, Dynamic, ColMajor> ArrayEXi;
typedef Array<        int,      15, Dynamic, ColMajor> ArrayFXi;
typedef Array<        int, Dynamic,       1, ColMajor> ArrayX1i;
typedef Array<        int, Dynamic,       2, ColMajor> ArrayX2i;
typedef Array<        int, Dynamic,       3, ColMajor> ArrayX3i;
typedef Array<        int, Dynamic,       4, ColMajor> ArrayX4i;
typedef Array<        int, Dynamic,       5, ColMajor> ArrayX5i;
typedef Array<        int, Dynamic,       6, ColMajor> ArrayX6i;
typedef Array<        int, Dynamic,       7, ColMajor> ArrayX7i;
typedef Array<        int, Dynamic,       8, ColMajor> ArrayX8i;
typedef Array<        int, Dynamic,       9, ColMajor> ArrayX9i;
typedef Array<        int, Dynamic,      10, ColMajor> ArrayXAi;
typedef Array<        int, Dynamic,      11, ColMajor> ArrayXBi;
typedef Array<        int, Dynamic,      12, ColMajor> ArrayXCi;
typedef Array<        int, Dynamic,      13, ColMajor> ArrayXDi;
typedef Array<        int, Dynamic,      14, ColMajor> ArrayXEi;
typedef Array<        int, Dynamic,      15, ColMajor> ArrayXFi;
typedef Array<        int, Dynamic, Dynamic, ColMajor> ArrayXXi;
/** @} */

/**
 * Typedefs for \ref real_t <a href="http://eigen.tuxfamily.org/">Eigen</a>
 * column vectors.  The names have been extended to handle fixed sizes up to 15
 * using hexadecimal notation.
 * @{
 */
typedef Matrix<    real_t,       1, 1, ColMajor> Vector1r;
typedef Matrix<    real_t,       2, 1, ColMajor> Vector2r;
typedef Matrix<    real_t,       3, 1, ColMajor> Vector3r;
typedef Matrix<    real_t,       4, 1, ColMajor> Vector4r;
typedef Matrix<    real_t,       5, 1, ColMajor> Vector5r;
typedef Matrix<    real_t,       6, 1, ColMajor> Vector6r;
typedef Matrix<    real_t,       7, 1, ColMajor> Vector7r;
typedef Matrix<    real_t,       8, 1, ColMajor> Vector8r;
typedef Matrix<    real_t,       9, 1, ColMajor> Vector9r;
typedef Matrix<    real_t,      10, 1, ColMajor> VectorAr;
typedef Matrix<    real_t,      11, 1, ColMajor> VectorBr;
typedef Matrix<    real_t,      12, 1, ColMajor> VectorCr;
typedef Matrix<    real_t,      13, 1, ColMajor> VectorDr;
typedef Matrix<    real_t,      14, 1, ColMajor> VectorEr;
typedef Matrix<    real_t,      15, 1, ColMajor> VectorFr;
typedef Matrix<    real_t, Dynamic, 1, ColMajor> VectorXr;
/** @} */

/**
 * Typedefs for \ref complex_t <a href="http://eigen.tuxfamily.org/">Eigen</a>
 * column vectors.  The names have been extended to handle fixed sizes up to 15
 * using hexadecimal notation.
 * @{
 */
typedef Matrix< complex_t,       1, 1, ColMajor> Vector1c;
typedef Matrix< complex_t,       2, 1, ColMajor> Vector2c;
typedef Matrix< complex_t,       3, 1, ColMajor> Vector3c;
typedef Matrix< complex_t,       4, 1, ColMajor> Vector4c;
typedef Matrix< complex_t,       5, 1, ColMajor> Vector5c;
typedef Matrix< complex_t,       6, 1, ColMajor> Vector6c;
typedef Matrix< complex_t,       7, 1, ColMajor> Vector7c;
typedef Matrix< complex_t,       8, 1, ColMajor> Vector8c;
typedef Matrix< complex_t,       9, 1, ColMajor> Vector9c;
typedef Matrix< complex_t,      10, 1, ColMajor> VectorAc;
typedef Matrix< complex_t,      11, 1, ColMajor> VectorBc;
typedef Matrix< complex_t,      12, 1, ColMajor> VectorCc;
typedef Matrix< complex_t,      13, 1, ColMajor> VectorDc;
typedef Matrix< complex_t,      14, 1, ColMajor> VectorEc;
typedef Matrix< complex_t,      15, 1, ColMajor> VectorFc;
typedef Matrix< complex_t, Dynamic, 1, ColMajor> VectorXc;
/** @} */

/**
 * Typedefs for \c int <a href="http://eigen.tuxfamily.org/">Eigen</a> column
 * vectors.  The names have been extended to handle fixed sizes up to 15 using
 * hexadecimal notation.
 * @{
 */
typedef Matrix<       int,       1, 1, ColMajor> Vector1i;
typedef Matrix<       int,       2, 1, ColMajor> Vector2i;
typedef Matrix<       int,       3, 1, ColMajor> Vector3i;
typedef Matrix<       int,       4, 1, ColMajor> Vector4i;
typedef Matrix<       int,       5, 1, ColMajor> Vector5i;
typedef Matrix<       int,       6, 1, ColMajor> Vector6i;
typedef Matrix<       int,       7, 1, ColMajor> Vector7i;
typedef Matrix<       int,       8, 1, ColMajor> Vector8i;
typedef Matrix<       int,       9, 1, ColMajor> Vector9i;
typedef Matrix<       int,      10, 1, ColMajor> VectorAi;
typedef Matrix<       int,      11, 1, ColMajor> VectorBi;
typedef Matrix<       int,      12, 1, ColMajor> VectorCi;
typedef Matrix<       int,      13, 1, ColMajor> VectorDi;
typedef Matrix<       int,      14, 1, ColMajor> VectorEi;
typedef Matrix<       int,      15, 1, ColMajor> VectorFi;
typedef Matrix<       int, Dynamic, 1, ColMajor> VectorXi;
/** @} */

/** @} */

} // end namespace suzerain

#endif // SUZERAIN_COMMON_HPP
