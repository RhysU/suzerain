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
// precision.hpp: Typedefs selecting the working precision
// $Id$

#ifndef PRECISION_HPP
#define PRECISION_HPP

#include <Eigen/Core>

// Scalar- and complex-valued typedefs
// Currently only real_t == double is supported by many, many components
typedef double                real_t;
typedef std::complex<real_t>  complex_t;

// Augment Eigen with some typedefs based on real_t and complex_t
// Allows us to write precision agnostic Eigen-base code
namespace Eigen {

    // Typedefs for real-valued matrices
    typedef Matrix<    real_t,       1,       1> Matrix1r ;
    typedef Matrix<    real_t,       2,       2> Matrix2r ;
    typedef Matrix<    real_t,       3,       3> Matrix3r ;
    typedef Matrix<    real_t,       4,       4> Matrix4r ;
    typedef Matrix<    real_t,       1, Dynamic> Matrix1Xr;
    typedef Matrix<    real_t,       2, Dynamic> Matrix2Xr;
    typedef Matrix<    real_t,       3, Dynamic> Matrix3Xr;
    typedef Matrix<    real_t,       4, Dynamic> Matrix4Xr;
    typedef Matrix<    real_t, Dynamic,       1> MatrixX1r;
    typedef Matrix<    real_t, Dynamic,       2> MatrixX2r;
    typedef Matrix<    real_t, Dynamic,       3> MatrixX3r;
    typedef Matrix<    real_t, Dynamic,       4> MatrixX4r;
    typedef Matrix<    real_t, Dynamic, Dynamic> MatrixXXr;

    // Typedefs for complex-valued matrices
    typedef Matrix< complex_t,       1,       1> Matrix1c ;
    typedef Matrix< complex_t,       2,       2> Matrix2c ;
    typedef Matrix< complex_t,       3,       3> Matrix3c ;
    typedef Matrix< complex_t,       4,       4> Matrix4c ;
    typedef Matrix< complex_t,       1, Dynamic> Matrix1Xc;
    typedef Matrix< complex_t,       2, Dynamic> Matrix2Xc;
    typedef Matrix< complex_t,       3, Dynamic> Matrix3Xc;
    typedef Matrix< complex_t,       4, Dynamic> Matrix4Xc;
    typedef Matrix< complex_t, Dynamic,       1> MatrixX1c;
    typedef Matrix< complex_t, Dynamic,       2> MatrixX2c;
    typedef Matrix< complex_t, Dynamic,       3> MatrixX3c;
    typedef Matrix< complex_t, Dynamic,       4> MatrixX4c;
    typedef Matrix< complex_t, Dynamic, Dynamic> MatrixXXc;

    // Typedefs for real-valued arrays
    typedef Array<    real_t,       1,       1> Array1r ;
    typedef Array<    real_t,       2,       1> Array2r ;
    typedef Array<    real_t,       3,       1> Array3r ;
    typedef Array<    real_t,       4,       1> Array4r ;
    typedef Array<    real_t, Dynamic,       1> ArrayXr ;
    typedef Array<    real_t,       1, Dynamic> Array1Xr;
    typedef Array<    real_t,       2, Dynamic> Array2Xr;
    typedef Array<    real_t,       3, Dynamic> Array3Xr;
    typedef Array<    real_t,       4, Dynamic> Array4Xr;
    typedef Array<    real_t, Dynamic,       1> ArrayX1r;
    typedef Array<    real_t, Dynamic,       2> ArrayX2r;
    typedef Array<    real_t, Dynamic,       3> ArrayX3r;
    typedef Array<    real_t, Dynamic,       4> ArrayX4r;
    typedef Array<    real_t, Dynamic, Dynamic> ArrayXXr;

    // Typedefs for complex-valued arrays
    typedef Array< complex_t,       1,       1> Array1c ;
    typedef Array< complex_t,       2,       2> Array2c ;
    typedef Array< complex_t,       3,       3> Array3c ;
    typedef Array< complex_t,       4,       4> Array4c ;
    typedef Array< complex_t, Dynamic,       1> ArrayXc ;
    typedef Array< complex_t,       1, Dynamic> Array1Xc;
    typedef Array< complex_t,       2, Dynamic> Array2Xc;
    typedef Array< complex_t,       3, Dynamic> Array3Xc;
    typedef Array< complex_t,       4, Dynamic> Array4Xc;
    typedef Array< complex_t, Dynamic,       1> ArrayX1c;
    typedef Array< complex_t, Dynamic,       2> ArrayX2c;
    typedef Array< complex_t, Dynamic,       3> ArrayX3c;
    typedef Array< complex_t, Dynamic,       4> ArrayX4c;
    typedef Array< complex_t, Dynamic, Dynamic> ArrayXXc;

    // Typedefs for real-valued vectors
    typedef Matrix<    real_t,       1, 1> Vector1r;
    typedef Matrix<    real_t,       2, 1> Vector2r;
    typedef Matrix<    real_t,       3, 1> Vector3r;
    typedef Matrix<    real_t,       4, 1> Vector4r;
    typedef Matrix<    real_t, Dynamic, 1> VectorXr;

    // Typedefs for complex-valued vectors
    typedef Matrix< complex_t,       1, 1> Vector1c;
    typedef Matrix< complex_t,       2, 1> Vector2c;
    typedef Matrix< complex_t,       3, 1> Vector3c;
    typedef Matrix< complex_t,       4, 1> Vector4c;
    typedef Matrix< complex_t, Dynamic, 1> VectorXc;
}

#endif // PRECISION_HPP
