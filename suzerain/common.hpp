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

// Include Eigen functionality used through Suzerain
#include <Eigen/Core>
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
#include <boost/typeof/typeof.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility.hpp>

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
 * evaluates to boolean \c false, then a <tt>std::logic_error</tt> is thrown
 * with message \c msg.
 *
 * This macro is intended for <tt>assert</tt>-like checks which should always
 * be performed regardless of whether or not \c NDEBUG is <tt>define</tt>d.
 */
#define SUZERAIN_ENSURE_MSG(expr, msg) \
    SUZERAIN_ENSURE_MSGEXCEPT(expr, msg, ::std::logic_error)

/**
 * Ensure that \c expr evaluates to boolean \c true at runtime.  If \c expr
 * evaluates to boolean \c false, then a <tt>std::logic_error</tt> is thrown.
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
using Eigen::Array;        /**< \namespace suzerain */
using Eigen::ColMajor;     /**< \namespace suzerain */
using Eigen::Dynamic;      /**< \namespace suzerain */
using Eigen::InnerStride;  /**< \namespace suzerain */
using Eigen::Map;          /**< \namespace suzerain */
using Eigen::Matrix;       /**< \namespace suzerain */
using Eigen::NoChange;     /**< \namespace suzerain */
using Eigen::OuterStride;  /**< \namespace suzerain */
using Eigen::RowMajor;     /**< \namespace suzerain */
using Eigen::Sequential;   /**< \namespace suzerain */
using Eigen::Unaligned;    /**< \namespace suzerain */

/**
 * Typedefs for \ref real_t, column-major <a
 * href="http://eigen.tuxfamily.org/">Eigen</a> matrices.
 * @{
 */
typedef Matrix<    real_t,       1,       1, ColMajor> Matrix1r ;
typedef Matrix<    real_t,       2,       2, ColMajor> Matrix2r ;
typedef Matrix<    real_t,       3,       3, ColMajor> Matrix3r ;
typedef Matrix<    real_t,       4,       4, ColMajor> Matrix4r ;
typedef Matrix<    real_t,       5,       5, ColMajor> Matrix5r ;
typedef Matrix<    real_t,       1, Dynamic, ColMajor> Matrix1Xr;
typedef Matrix<    real_t,       2, Dynamic, ColMajor> Matrix2Xr;
typedef Matrix<    real_t,       3, Dynamic, ColMajor> Matrix3Xr;
typedef Matrix<    real_t,       4, Dynamic, ColMajor> Matrix4Xr;
typedef Matrix<    real_t,       5, Dynamic, ColMajor> Matrix5Xr;
typedef Matrix<    real_t, Dynamic,       1, ColMajor> MatrixX1r;
typedef Matrix<    real_t, Dynamic,       2, ColMajor> MatrixX2r;
typedef Matrix<    real_t, Dynamic,       3, ColMajor> MatrixX3r;
typedef Matrix<    real_t, Dynamic,       4, ColMajor> MatrixX4r;
typedef Matrix<    real_t, Dynamic,       5, ColMajor> MatrixX5r;
typedef Matrix<    real_t, Dynamic, Dynamic, ColMajor> MatrixXXr;
/** @} */

/**
 * Typedefs for \ref complex_t, column-major <a
 * href="http://eigen.tuxfamily.org/">Eigen</a> matrices.
 * @{
 */
typedef Matrix< complex_t,       1,       1, ColMajor> Matrix1c ;
typedef Matrix< complex_t,       2,       2, ColMajor> Matrix2c ;
typedef Matrix< complex_t,       3,       3, ColMajor> Matrix3c ;
typedef Matrix< complex_t,       4,       4, ColMajor> Matrix4c ;
typedef Matrix< complex_t,       5,       5, ColMajor> Matrix5c ;
typedef Matrix< complex_t,       1, Dynamic, ColMajor> Matrix1Xc;
typedef Matrix< complex_t,       2, Dynamic, ColMajor> Matrix2Xc;
typedef Matrix< complex_t,       3, Dynamic, ColMajor> Matrix3Xc;
typedef Matrix< complex_t,       4, Dynamic, ColMajor> Matrix4Xc;
typedef Matrix< complex_t,       5, Dynamic, ColMajor> Matrix5Xc;
typedef Matrix< complex_t, Dynamic,       1, ColMajor> MatrixX1c;
typedef Matrix< complex_t, Dynamic,       2, ColMajor> MatrixX2c;
typedef Matrix< complex_t, Dynamic,       3, ColMajor> MatrixX3c;
typedef Matrix< complex_t, Dynamic,       4, ColMajor> MatrixX4c;
typedef Matrix< complex_t, Dynamic,       5, ColMajor> MatrixX5c;
typedef Matrix< complex_t, Dynamic, Dynamic, ColMajor> MatrixXXc;
/** @} */

/**
 * Typedefs for \ref real_t, column-major <a
 * href="http://eigen.tuxfamily.org/">Eigen</a> matrices.
 * @{
 */
typedef Matrix<       int,       1,       1, ColMajor> Matrix1i ;
typedef Matrix<       int,       2,       2, ColMajor> Matrix2i ;
typedef Matrix<       int,       3,       3, ColMajor> Matrix3i ;
typedef Matrix<       int,       4,       4, ColMajor> Matrix4i ;
typedef Matrix<       int,       5,       5, ColMajor> Matrix5i ;
typedef Matrix<       int,       1, Dynamic, ColMajor> Matrix1Xi;
typedef Matrix<       int,       2, Dynamic, ColMajor> Matrix2Xi;
typedef Matrix<       int,       3, Dynamic, ColMajor> Matrix3Xi;
typedef Matrix<       int,       4, Dynamic, ColMajor> Matrix4Xi;
typedef Matrix<       int,       5, Dynamic, ColMajor> Matrix5Xi;
typedef Matrix<       int, Dynamic,       1, ColMajor> MatrixX1i;
typedef Matrix<       int, Dynamic,       2, ColMajor> MatrixX2i;
typedef Matrix<       int, Dynamic,       3, ColMajor> MatrixX3i;
typedef Matrix<       int, Dynamic,       4, ColMajor> MatrixX4i;
typedef Matrix<       int, Dynamic,       5, ColMajor> MatrixX5i;
typedef Matrix<       int, Dynamic, Dynamic, ColMajor> MatrixXXi;
/** @} */

/**
 * Typedefs for \ref real_t, column-major <a
 * href="http://eigen.tuxfamily.org/">Eigen</a> arrays.
 * @{
 */
typedef Array<     real_t,       1,       1, ColMajor> Array1r ;
typedef Array<     real_t,       2,       1, ColMajor> Array2r ;
typedef Array<     real_t,       3,       1, ColMajor> Array3r ;
typedef Array<     real_t,       4,       1, ColMajor> Array4r ;
typedef Array<     real_t,       5,       1, ColMajor> Array5r ;
typedef Array<     real_t, Dynamic,       1, ColMajor> ArrayXr ;
typedef Array<     real_t,       1, Dynamic, ColMajor> Array1Xr;
typedef Array<     real_t,       2, Dynamic, ColMajor> Array2Xr;
typedef Array<     real_t,       3, Dynamic, ColMajor> Array3Xr;
typedef Array<     real_t,       4, Dynamic, ColMajor> Array4Xr;
typedef Array<     real_t,       5, Dynamic, ColMajor> Array5Xr;
typedef Array<     real_t, Dynamic,       1, ColMajor> ArrayX1r;
typedef Array<     real_t, Dynamic,       2, ColMajor> ArrayX2r;
typedef Array<     real_t, Dynamic,       3, ColMajor> ArrayX3r;
typedef Array<     real_t, Dynamic,       4, ColMajor> ArrayX4r;
typedef Array<     real_t, Dynamic,       5, ColMajor> ArrayX5r;
typedef Array<     real_t, Dynamic, Dynamic, ColMajor> ArrayXXr;
/** @} */

/**
 * Typedefs for \ref complex_t, column-major <a
 * href="http://eigen.tuxfamily.org/">Eigen</a> arrays.
 * @{
 */
typedef Array<  complex_t,       1,       1, ColMajor> Array1c ;
typedef Array<  complex_t,       2,       2, ColMajor> Array2c ;
typedef Array<  complex_t,       3,       3, ColMajor> Array3c ;
typedef Array<  complex_t,       4,       4, ColMajor> Array4c ;
typedef Array<  complex_t,       5,       5, ColMajor> Array5c ;
typedef Array<  complex_t, Dynamic,       1, ColMajor> ArrayXc ;
typedef Array<  complex_t,       1, Dynamic, ColMajor> Array1Xc;
typedef Array<  complex_t,       2, Dynamic, ColMajor> Array2Xc;
typedef Array<  complex_t,       3, Dynamic, ColMajor> Array3Xc;
typedef Array<  complex_t,       4, Dynamic, ColMajor> Array4Xc;
typedef Array<  complex_t,       5, Dynamic, ColMajor> Array5Xc;
typedef Array<  complex_t, Dynamic,       1, ColMajor> ArrayX1c;
typedef Array<  complex_t, Dynamic,       2, ColMajor> ArrayX2c;
typedef Array<  complex_t, Dynamic,       3, ColMajor> ArrayX3c;
typedef Array<  complex_t, Dynamic,       4, ColMajor> ArrayX4c;
typedef Array<  complex_t, Dynamic,       5, ColMajor> ArrayX5c;
typedef Array<  complex_t, Dynamic, Dynamic, ColMajor> ArrayXXc;
/** @} */

/**
 * Typedefs for \c int, column-major <a
 * href="http://eigen.tuxfamily.org/">Eigen</a> arrays.
 * @{
 */
typedef Array<        int,       1,       1, ColMajor> Array1i ;
typedef Array<        int,       2,       1, ColMajor> Array2i ;
typedef Array<        int,       3,       1, ColMajor> Array3i ;
typedef Array<        int,       4,       1, ColMajor> Array4i ;
typedef Array<        int,       5,       1, ColMajor> Array5i ;
typedef Array<        int, Dynamic,       1, ColMajor> ArrayXi ;
typedef Array<        int,       1, Dynamic, ColMajor> Array1Xi;
typedef Array<        int,       2, Dynamic, ColMajor> Array2Xi;
typedef Array<        int,       3, Dynamic, ColMajor> Array3Xi;
typedef Array<        int,       4, Dynamic, ColMajor> Array4Xi;
typedef Array<        int,       5, Dynamic, ColMajor> Array5Xi;
typedef Array<        int, Dynamic,       1, ColMajor> ArrayX1i;
typedef Array<        int, Dynamic,       2, ColMajor> ArrayX2i;
typedef Array<        int, Dynamic,       3, ColMajor> ArrayX3i;
typedef Array<        int, Dynamic,       4, ColMajor> ArrayX4i;
typedef Array<        int, Dynamic,       5, ColMajor> ArrayX5i;
typedef Array<        int, Dynamic, Dynamic, ColMajor> ArrayXXi;
/** @} */

/**
 * Typedefs for \ref real_t <a href="http://eigen.tuxfamily.org/">Eigen</a>
 * column vectors.
 * @{
 */
typedef Matrix<    real_t,       1, 1, ColMajor> Vector1r;
typedef Matrix<    real_t,       2, 1, ColMajor> Vector2r;
typedef Matrix<    real_t,       3, 1, ColMajor> Vector3r;
typedef Matrix<    real_t,       4, 1, ColMajor> Vector4r;
typedef Matrix<    real_t,       5, 1, ColMajor> Vector5r;
typedef Matrix<    real_t, Dynamic, 1, ColMajor> VectorXr;
/** @} */

/**
 * Typedefs for \ref complex_t <a href="http://eigen.tuxfamily.org/">Eigen</a>
 * column vectors.
 * @{
 */
typedef Matrix< complex_t,       1, 1, ColMajor> Vector1c;
typedef Matrix< complex_t,       2, 1, ColMajor> Vector2c;
typedef Matrix< complex_t,       3, 1, ColMajor> Vector3c;
typedef Matrix< complex_t,       4, 1, ColMajor> Vector4c;
typedef Matrix< complex_t,       5, 1, ColMajor> Vector5c;
typedef Matrix< complex_t, Dynamic, 1, ColMajor> VectorXc;
/** @} */

/**
 * Typedefs for \c int <a href="http://eigen.tuxfamily.org/">Eigen</a>
 * column vectors.
 * @{
 */
typedef Matrix<       int,       1, 1, ColMajor> Vector1i;
typedef Matrix<       int,       2, 1, ColMajor> Vector2i;
typedef Matrix<       int,       3, 1, ColMajor> Vector3i;
typedef Matrix<       int,       4, 1, ColMajor> Vector4i;
typedef Matrix<       int,       5, 1, ColMajor> Vector5i;
typedef Matrix<       int, Dynamic, 1, ColMajor> VectorXi;
/** @} */

/** @} */

} // end namespace suzerain

#endif // SUZERAIN_COMMON_HPP
