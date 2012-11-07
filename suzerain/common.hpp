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

/** \namespace suzerain */
using boost::shared_ptr;

/** \namespace suzerain */
using boost::scoped_ptr;

/** \namespace suzerain */
using boost::shared_array;

/** \namespace suzerain */
using boost::scoped_array;

/** \weakgroup EigenTypedefs Typedefs and declarations to simplify Eigen usage.
 *
 * Defines typedefs for <a href="http://eigen.tuxfamily.org/">Eigen</a>-related
 * classes within the \ref suzerain namespace.  Typedefs like \ref Matrix3r and
 * \ref Vector3c permit simplifying using \ref suzerain::real_t and \ref
 * suzerain::complex_t with Eigen's #Matrix and #Array templates.  They also
 * permit codifying Suzerain's expectations regarding column- versus row-major
 * storage regardless of the value of \c EIGEN_DEFAULT_TO_ROW_MAJOR.
 *
 * @{
 */

/**
 * Typedefs for \ref real_t, column-major <a
 * href="http://eigen.tuxfamily.org/">Eigen</a> matrices.
 * @{
 */
typedef Eigen::Matrix<    real_t,              1,              1, Eigen::ColMajor> Matrix1r ;
typedef Eigen::Matrix<    real_t,              2,              2, Eigen::ColMajor> Matrix2r ;
typedef Eigen::Matrix<    real_t,              3,              3, Eigen::ColMajor> Matrix3r ;
typedef Eigen::Matrix<    real_t,              4,              4, Eigen::ColMajor> Matrix4r ;
typedef Eigen::Matrix<    real_t,              5,              5, Eigen::ColMajor> Matrix5r ;
typedef Eigen::Matrix<    real_t,              1, Eigen::Dynamic, Eigen::ColMajor> Matrix1Xr;
typedef Eigen::Matrix<    real_t,              2, Eigen::Dynamic, Eigen::ColMajor> Matrix2Xr;
typedef Eigen::Matrix<    real_t,              3, Eigen::Dynamic, Eigen::ColMajor> Matrix3Xr;
typedef Eigen::Matrix<    real_t,              4, Eigen::Dynamic, Eigen::ColMajor> Matrix4Xr;
typedef Eigen::Matrix<    real_t,              5, Eigen::Dynamic, Eigen::ColMajor> Matrix5Xr;
typedef Eigen::Matrix<    real_t, Eigen::Dynamic,              1, Eigen::ColMajor> MatrixX1r;
typedef Eigen::Matrix<    real_t, Eigen::Dynamic,              2, Eigen::ColMajor> MatrixX2r;
typedef Eigen::Matrix<    real_t, Eigen::Dynamic,              3, Eigen::ColMajor> MatrixX3r;
typedef Eigen::Matrix<    real_t, Eigen::Dynamic,              4, Eigen::ColMajor> MatrixX4r;
typedef Eigen::Matrix<    real_t, Eigen::Dynamic,              5, Eigen::ColMajor> MatrixX5r;
typedef Eigen::Matrix<    real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> MatrixXXr;
/** @} */

/**
 * Typedefs for \ref complex_t, column-major <a
 * href="http://eigen.tuxfamily.org/">Eigen</a> matrices.
 * @{
 */
typedef Eigen::Matrix< complex_t,              1,              1, Eigen::ColMajor> Matrix1c ;
typedef Eigen::Matrix< complex_t,              2,              2, Eigen::ColMajor> Matrix2c ;
typedef Eigen::Matrix< complex_t,              3,              3, Eigen::ColMajor> Matrix3c ;
typedef Eigen::Matrix< complex_t,              4,              4, Eigen::ColMajor> Matrix4c ;
typedef Eigen::Matrix< complex_t,              5,              5, Eigen::ColMajor> Matrix5c ;
typedef Eigen::Matrix< complex_t,              1, Eigen::Dynamic, Eigen::ColMajor> Matrix1Xc;
typedef Eigen::Matrix< complex_t,              2, Eigen::Dynamic, Eigen::ColMajor> Matrix2Xc;
typedef Eigen::Matrix< complex_t,              3, Eigen::Dynamic, Eigen::ColMajor> Matrix3Xc;
typedef Eigen::Matrix< complex_t,              4, Eigen::Dynamic, Eigen::ColMajor> Matrix4Xc;
typedef Eigen::Matrix< complex_t,              5, Eigen::Dynamic, Eigen::ColMajor> Matrix5Xc;
typedef Eigen::Matrix< complex_t, Eigen::Dynamic,              1, Eigen::ColMajor> MatrixX1c;
typedef Eigen::Matrix< complex_t, Eigen::Dynamic,              2, Eigen::ColMajor> MatrixX2c;
typedef Eigen::Matrix< complex_t, Eigen::Dynamic,              3, Eigen::ColMajor> MatrixX3c;
typedef Eigen::Matrix< complex_t, Eigen::Dynamic,              4, Eigen::ColMajor> MatrixX4c;
typedef Eigen::Matrix< complex_t, Eigen::Dynamic,              5, Eigen::ColMajor> MatrixX5c;
typedef Eigen::Matrix< complex_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> MatrixXXc;
/** @} */

/**
 * Typedefs for \ref real_t, column-major <a
 * href="http://eigen.tuxfamily.org/">Eigen</a> matrices.
 * @{
 */
typedef Eigen::Matrix<       int,              1,              1, Eigen::ColMajor> Matrix1i ;
typedef Eigen::Matrix<       int,              2,              2, Eigen::ColMajor> Matrix2i ;
typedef Eigen::Matrix<       int,              3,              3, Eigen::ColMajor> Matrix3i ;
typedef Eigen::Matrix<       int,              4,              4, Eigen::ColMajor> Matrix4i ;
typedef Eigen::Matrix<       int,              5,              5, Eigen::ColMajor> Matrix5i ;
typedef Eigen::Matrix<       int,              1, Eigen::Dynamic, Eigen::ColMajor> Matrix1Xi;
typedef Eigen::Matrix<       int,              2, Eigen::Dynamic, Eigen::ColMajor> Matrix2Xi;
typedef Eigen::Matrix<       int,              3, Eigen::Dynamic, Eigen::ColMajor> Matrix3Xi;
typedef Eigen::Matrix<       int,              4, Eigen::Dynamic, Eigen::ColMajor> Matrix4Xi;
typedef Eigen::Matrix<       int,              5, Eigen::Dynamic, Eigen::ColMajor> Matrix5Xi;
typedef Eigen::Matrix<       int, Eigen::Dynamic,              1, Eigen::ColMajor> MatrixX1i;
typedef Eigen::Matrix<       int, Eigen::Dynamic,              2, Eigen::ColMajor> MatrixX2i;
typedef Eigen::Matrix<       int, Eigen::Dynamic,              3, Eigen::ColMajor> MatrixX3i;
typedef Eigen::Matrix<       int, Eigen::Dynamic,              4, Eigen::ColMajor> MatrixX4i;
typedef Eigen::Matrix<       int, Eigen::Dynamic,              5, Eigen::ColMajor> MatrixX5i;
typedef Eigen::Matrix<       int, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> MatrixXXi;
/** @} */

/**
 * Typedefs for \ref real_t, column-major <a
 * href="http://eigen.tuxfamily.org/">Eigen</a> arrays.
 * @{
 */
typedef Eigen::Array<     real_t,              1,              1, Eigen::ColMajor> Array1r ;
typedef Eigen::Array<     real_t,              2,              1, Eigen::ColMajor> Array2r ;
typedef Eigen::Array<     real_t,              3,              1, Eigen::ColMajor> Array3r ;
typedef Eigen::Array<     real_t,              4,              1, Eigen::ColMajor> Array4r ;
typedef Eigen::Array<     real_t,              5,              1, Eigen::ColMajor> Array5r ;
typedef Eigen::Array<     real_t, Eigen::Dynamic,              1, Eigen::ColMajor> ArrayXr ;
typedef Eigen::Array<     real_t,              1, Eigen::Dynamic, Eigen::ColMajor> Array1Xr;
typedef Eigen::Array<     real_t,              2, Eigen::Dynamic, Eigen::ColMajor> Array2Xr;
typedef Eigen::Array<     real_t,              3, Eigen::Dynamic, Eigen::ColMajor> Array3Xr;
typedef Eigen::Array<     real_t,              4, Eigen::Dynamic, Eigen::ColMajor> Array4Xr;
typedef Eigen::Array<     real_t,              5, Eigen::Dynamic, Eigen::ColMajor> Array5Xr;
typedef Eigen::Array<     real_t, Eigen::Dynamic,              1, Eigen::ColMajor> ArrayX1r;
typedef Eigen::Array<     real_t, Eigen::Dynamic,              2, Eigen::ColMajor> ArrayX2r;
typedef Eigen::Array<     real_t, Eigen::Dynamic,              3, Eigen::ColMajor> ArrayX3r;
typedef Eigen::Array<     real_t, Eigen::Dynamic,              4, Eigen::ColMajor> ArrayX4r;
typedef Eigen::Array<     real_t, Eigen::Dynamic,              5, Eigen::ColMajor> ArrayX5r;
typedef Eigen::Array<     real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ArrayXXr;
/** @} */

/**
 * Typedefs for \ref complex_t, column-major <a
 * href="http://eigen.tuxfamily.org/">Eigen</a> arrays.
 * @{
 */
typedef Eigen::Array<  complex_t,              1,              1, Eigen::ColMajor> Array1c ;
typedef Eigen::Array<  complex_t,              2,              2, Eigen::ColMajor> Array2c ;
typedef Eigen::Array<  complex_t,              3,              3, Eigen::ColMajor> Array3c ;
typedef Eigen::Array<  complex_t,              4,              4, Eigen::ColMajor> Array4c ;
typedef Eigen::Array<  complex_t,              5,              5, Eigen::ColMajor> Array5c ;
typedef Eigen::Array<  complex_t, Eigen::Dynamic,              1, Eigen::ColMajor> ArrayXc ;
typedef Eigen::Array<  complex_t,              1, Eigen::Dynamic, Eigen::ColMajor> Array1Xc;
typedef Eigen::Array<  complex_t,              2, Eigen::Dynamic, Eigen::ColMajor> Array2Xc;
typedef Eigen::Array<  complex_t,              3, Eigen::Dynamic, Eigen::ColMajor> Array3Xc;
typedef Eigen::Array<  complex_t,              4, Eigen::Dynamic, Eigen::ColMajor> Array4Xc;
typedef Eigen::Array<  complex_t,              5, Eigen::Dynamic, Eigen::ColMajor> Array5Xc;
typedef Eigen::Array<  complex_t, Eigen::Dynamic,              1, Eigen::ColMajor> ArrayX1c;
typedef Eigen::Array<  complex_t, Eigen::Dynamic,              2, Eigen::ColMajor> ArrayX2c;
typedef Eigen::Array<  complex_t, Eigen::Dynamic,              3, Eigen::ColMajor> ArrayX3c;
typedef Eigen::Array<  complex_t, Eigen::Dynamic,              4, Eigen::ColMajor> ArrayX4c;
typedef Eigen::Array<  complex_t, Eigen::Dynamic,              5, Eigen::ColMajor> ArrayX5c;
typedef Eigen::Array<  complex_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ArrayXXc;
/** @} */

/**
 * Typedefs for \c int, column-major <a
 * href="http://eigen.tuxfamily.org/">Eigen</a> arrays.
 * @{
 */
typedef Eigen::Array<        int,              1,              1, Eigen::ColMajor> Array1i ;
typedef Eigen::Array<        int,              2,              1, Eigen::ColMajor> Array2i ;
typedef Eigen::Array<        int,              3,              1, Eigen::ColMajor> Array3i ;
typedef Eigen::Array<        int,              4,              1, Eigen::ColMajor> Array4i ;
typedef Eigen::Array<        int,              5,              1, Eigen::ColMajor> Array5i ;
typedef Eigen::Array<        int, Eigen::Dynamic,              1, Eigen::ColMajor> ArrayXi ;
typedef Eigen::Array<        int,              1, Eigen::Dynamic, Eigen::ColMajor> Array1Xi;
typedef Eigen::Array<        int,              2, Eigen::Dynamic, Eigen::ColMajor> Array2Xi;
typedef Eigen::Array<        int,              3, Eigen::Dynamic, Eigen::ColMajor> Array3Xi;
typedef Eigen::Array<        int,              4, Eigen::Dynamic, Eigen::ColMajor> Array4Xi;
typedef Eigen::Array<        int,              5, Eigen::Dynamic, Eigen::ColMajor> Array5Xi;
typedef Eigen::Array<        int, Eigen::Dynamic,              1, Eigen::ColMajor> ArrayX1i;
typedef Eigen::Array<        int, Eigen::Dynamic,              2, Eigen::ColMajor> ArrayX2i;
typedef Eigen::Array<        int, Eigen::Dynamic,              3, Eigen::ColMajor> ArrayX3i;
typedef Eigen::Array<        int, Eigen::Dynamic,              4, Eigen::ColMajor> ArrayX4i;
typedef Eigen::Array<        int, Eigen::Dynamic,              5, Eigen::ColMajor> ArrayX5i;
typedef Eigen::Array<        int, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ArrayXXi;
/** @} */

/**
 * Typedefs for \ref real_t <a href="http://eigen.tuxfamily.org/">Eigen</a>
 * column vectors.
 * @{
 */
typedef Eigen::Matrix<    real_t,              1, 1, Eigen::ColMajor> Vector1r;
typedef Eigen::Matrix<    real_t,              2, 1, Eigen::ColMajor> Vector2r;
typedef Eigen::Matrix<    real_t,              3, 1, Eigen::ColMajor> Vector3r;
typedef Eigen::Matrix<    real_t,              4, 1, Eigen::ColMajor> Vector4r;
typedef Eigen::Matrix<    real_t,              5, 1, Eigen::ColMajor> Vector5r;
typedef Eigen::Matrix<    real_t, Eigen::Dynamic, 1, Eigen::ColMajor> VectorXr;
/** @} */

/**
 * Typedefs for \ref complex_t <a href="http://eigen.tuxfamily.org/">Eigen</a>
 * column vectors.
 * @{
 */
typedef Eigen::Matrix< complex_t,              1, 1, Eigen::ColMajor> Vector1c;
typedef Eigen::Matrix< complex_t,              2, 1, Eigen::ColMajor> Vector2c;
typedef Eigen::Matrix< complex_t,              3, 1, Eigen::ColMajor> Vector3c;
typedef Eigen::Matrix< complex_t,              4, 1, Eigen::ColMajor> Vector4c;
typedef Eigen::Matrix< complex_t,              5, 1, Eigen::ColMajor> Vector5c;
typedef Eigen::Matrix< complex_t, Eigen::Dynamic, 1, Eigen::ColMajor> VectorXc;
/** @} */

/**
 * Typedefs for \c int <a href="http://eigen.tuxfamily.org/">Eigen</a>
 * column vectors.
 * @{
 */
typedef Eigen::Matrix<       int,              1, 1, Eigen::ColMajor> Vector1i;
typedef Eigen::Matrix<       int,              2, 1, Eigen::ColMajor> Vector2i;
typedef Eigen::Matrix<       int,              3, 1, Eigen::ColMajor> Vector3i;
typedef Eigen::Matrix<       int,              4, 1, Eigen::ColMajor> Vector4i;
typedef Eigen::Matrix<       int,              5, 1, Eigen::ColMajor> Vector5i;
typedef Eigen::Matrix<       int, Eigen::Dynamic, 1, Eigen::ColMajor> VectorXi;
/** @} */

// Make <code>suzerain::Map</code> behave like <code>Eigen::Map</code>.
/** \namespace suzerain */
using Eigen::Map;

/** @} */

} // end namespace suzerain

#endif // SUZERAIN_COMMON_HPP
