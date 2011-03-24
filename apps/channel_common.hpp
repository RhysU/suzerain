/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Suzerain is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * channel_common.hpp: Channel-related functionality spanning binaries
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef CHANNEL_COMMON_HPP
#define CHANNEL_COMMON_HPP

#include <Eigen/Core>
#include <esio/esio.h>
#include <suzerain/bspline.hpp>
#include <suzerain/diffwave.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/inorder.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/scenario_definition.hpp>
#include <suzerain/state.hpp>
#include <suzerain/state_impl.hpp>
#include <suzerain/storage.hpp>
#include <suzerain/timestepper.hpp>

// Scalar- and complex-valued typedefs
// Currently only real_t == double is supported by many, many components
typedef double                real_t;
typedef std::complex<real_t>  complex_t;

// Augment Eigen with some typedefs based on real_t and complex_t
// Allows us to write precision agnostic code
namespace Eigen {

    // Typedefs for real-valued matrices
    typedef Matrix<    real_t,       2,       2> Matrix2r ;
    typedef Matrix<    real_t,       3,       3> Matrix3r ;
    typedef Matrix<    real_t,       4,       4> Matrix4r ;
    typedef Matrix<    real_t,       2, Dynamic> Matrix2Xr;
    typedef Matrix<    real_t,       3, Dynamic> Matrix3Xr;
    typedef Matrix<    real_t,       4, Dynamic> Matrix4Xr;
    typedef Matrix<    real_t, Dynamic,       2> MatrixX2r;
    typedef Matrix<    real_t, Dynamic,       3> MatrixX3r;
    typedef Matrix<    real_t, Dynamic,       4> MatrixX4r;
    typedef Matrix<    real_t, Dynamic, Dynamic> MatrixXXr;

    // Typedefs for complex-valued matrices
    typedef Matrix< complex_t,       2,       2> Matrix2c ;
    typedef Matrix< complex_t,       3,       3> Matrix3c ;
    typedef Matrix< complex_t,       4,       4> Matrix4c ;
    typedef Matrix< complex_t,       2, Dynamic> Matrix2Xc;
    typedef Matrix< complex_t,       3, Dynamic> Matrix3Xc;
    typedef Matrix< complex_t,       4, Dynamic> Matrix4Xc;
    typedef Matrix< complex_t, Dynamic,       2> MatrixX2c;
    typedef Matrix< complex_t, Dynamic,       3> MatrixX3c;
    typedef Matrix< complex_t, Dynamic,       4> MatrixX4c;
    typedef Matrix< complex_t, Dynamic, Dynamic> MatrixXXc;

    // Typedefs for real-valued arrays
    typedef Array<    real_t,       2,       1> Array2r ;
    typedef Array<    real_t,       3,       1> Array3r ;
    typedef Array<    real_t,       4,       1> Array4r ;
    typedef Array<    real_t, Dynamic,       1> ArrayXr ;
    typedef Array<    real_t,       2, Dynamic> Array2Xr;
    typedef Array<    real_t,       3, Dynamic> Array3Xr;
    typedef Array<    real_t,       4, Dynamic> Array4Xr;
    typedef Array<    real_t, Dynamic,       2> ArrayX2r;
    typedef Array<    real_t, Dynamic,       3> ArrayX3r;
    typedef Array<    real_t, Dynamic,       4> ArrayX4r;
    typedef Array<    real_t, Dynamic, Dynamic> ArrayXXr;

    // Typedefs for complex-valued arrays
    typedef Array< complex_t,       2,       2> Array2c ;
    typedef Array< complex_t,       3,       3> Array3c ;
    typedef Array< complex_t,       4,       4> Array4c ;
    typedef Array< complex_t, Dynamic,       1> ArrayXc ;
    typedef Array< complex_t,       2, Dynamic> Array2Xc;
    typedef Array< complex_t,       3, Dynamic> Array3Xc;
    typedef Array< complex_t,       4, Dynamic> Array4Xc;
    typedef Array< complex_t, Dynamic,       2> ArrayX2c;
    typedef Array< complex_t, Dynamic,       3> ArrayX3c;
    typedef Array< complex_t, Dynamic,       4> ArrayX4c;
    typedef Array< complex_t, Dynamic, Dynamic> ArrayXXc;

    // Typedefs for real-valued vectors
    typedef Matrix<    real_t,       2, 1> Vector2r;
    typedef Matrix<    real_t,       3, 1> Vector3r;
    typedef Matrix<    real_t,       4, 1> Vector4r;
    typedef Matrix<    real_t, Dynamic, 1> VectorXr;

    // Typedefs for complex-valued vectors
    typedef Matrix< complex_t,       2, 1> Vector2c;
    typedef Matrix< complex_t,       3, 1> Vector3c;
    typedef Matrix< complex_t,       4, 1> Vector4c;
    typedef Matrix< complex_t, Dynamic, 1> VectorXc;
}

/** Field names as stored in restart files */
extern const boost::array<const char *,5> field_names;

/** Field descriptions for use in restart file comments */
extern const boost::array<const char *,5> field_descriptions;

/** Store a ScenarioDefinition in a restart file */
void store(const esio_handle esioh,
           const suzerain::problem::ScenarioDefinition<real_t>& scenario);

/** Load a ScenarioDefinition from a restart file */
void load(const esio_handle esioh,
          suzerain::problem::ScenarioDefinition<real_t>& scenario);

/** Store a GridDefinition in a restart file */
void store(const esio_handle esioh,
           const suzerain::problem::GridDefinition& grid,
           const real_t Lx,
           const real_t Lz);

/** Load a GridDefinition from a restart file */
void load(const esio_handle esioh,
          suzerain::problem::GridDefinition& grid);

/** Create a B-spline workspace on [a,b] per ndof, k, and htdelta */
void create(const int ndof,
            const int k,
            const double a,
            const double b,
            const double htdelta,
            boost::shared_ptr<const suzerain::bspline>& bspw);

/** Store a suzerain::bspline workspace in a restart file */
void store(const esio_handle esioh,
           const boost::shared_ptr<const suzerain::bspline>& bspw);

/** Load a suzerain::bspline workspace from a restart file */
void load(const esio_handle esioh,
          boost::shared_ptr<const suzerain::bspline>& bspw);

/** Store the current simulation time information */
void store_time(const esio_handle esioh,
                real_t time);

/** Load the current simulation time information */
void load_time(const esio_handle esioh,
               real_t &time);

/** Load the current simulation state from an open restart file in h */
void load(const esio_handle esioh,
          suzerain::NoninterleavedState<4,complex_t> &state,
          const suzerain::problem::GridDefinition& grid,
          const suzerain::pencil_grid& dgrid,
          const suzerain::bspline& bspw);

/** Read a complex-valued field via ESIO */
template< typename I >
inline
void complex_field_read(esio_handle h, const char *name, complex_t *field,
                        I cstride = 0, I bstride = 0, I astride = 0)
{
    using boost::numeric_cast;

    esio_field_readv(h, name, reinterpret_cast<real_t *>(field),
                     2*numeric_cast<int>(cstride),
                     2*numeric_cast<int>(bstride),
                     2*numeric_cast<int>(astride),
                     2);
}

/** Read a complex-valued field via ESIO */
inline
void complex_field_read(esio_handle h, const char *name, complex_t *field)
{
    // When no strides are provided, we must specify the stride type.
    return complex_field_read<int>(h, name, field);
}

/** Write a complex-valued field via ESIO */
template< typename I >
inline
void complex_field_write(esio_handle h,
                         const char *name, const complex_t *field,
                         I cstride = 0, I bstride = 0, I astride = 0,
                         const char * comment = 0)
{
    using boost::numeric_cast;

    esio_field_writev(h, name, reinterpret_cast<const real_t *>(field),
                      2*numeric_cast<int>(cstride),
                      2*numeric_cast<int>(bstride),
                      2*numeric_cast<int>(astride),
                      2, comment);
}

/** Write a complex-valued field via ESIO */
inline
void complex_field_write(esio_handle h,
                         const char *name, const complex_t *field)
{
    // When no strides are provided, we must specify the stride type.
    return complex_field_write<int>(h, name, field);
}

#endif // CHANNEL_COMMON_HPP
