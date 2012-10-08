//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// precision.hpp: Typedefs related to the working precision
// $Id$

#ifndef SUZERAIN_PRECISION_HPP
#define SUZERAIN_PRECISION_HPP

#include <Eigen/Core>

// TODO real_t and complex_t are not, strictly speaking, POSIX compliant names
// But they are, strictly speaking, shorter than real_type and complex_type.

/** @file
 * Typedefs related to the working precision.  Many of these typedefs are for
 * classes within Eigen.
 *
 * @see <a href="http://eigen.tuxfamily.org/">Eigen</a> for more information.
 */

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

// Make <code>suzerain::Map</code> behave like <code>Eigen::Map</code>.
using Eigen::Map;

}  // end namespace suzerain

#endif // SUZERAIN_PRECISION_HPP
