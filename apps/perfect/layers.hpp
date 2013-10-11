//--------------------------------------------------------------------------
//
// Copyright (C) 2013 Rhys Ulerich
// Copyright (C) 2013 The PECOS Development Team
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

#ifndef SUZERAIN_PERFECT_LAYERS_HPP
#define SUZERAIN_PERFECT_LAYERS_HPP

/** @file
 * Sampling logistics for boundary layer profiles.
 */

#include <suzerain/common.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/state_fwd.hpp>

namespace suzerain {

// Forward declarations
class grid_specification;
class pencil_grid;

namespace perfect {

// Forward declarations
class scenario_definition;

/**
 * Encapsulates boundary layer profiles necessary for using \ref bl.h.
 *
 * Samples of each quantity are made available through two-dimensional
 * column-major arrays.  The row index iterates over wall-normal B-spline
 * coefficients and the column index iterates over tensor indices.  Scalars
 * have only a single tensor index.  Vector layers have two indices
 * corresponding to the streamwise x and wall-normal y directions.  The
 * homogeneous spanwise direction z is not reported.
 */
class layers
{
public:

/**
 * A Boost.Preprocessor sequence of tuples of layers computed in wave
 * space.
 */
#define SUZERAIN_PERFECT_LAYERS_WAVE             \
    ((rho,                      1)) /* scalar */ \
    ((rho_u,                    2)) /* vector */

/**
 * A Boost.Preprocessor sequence of tuples of layers computed in physical
 * space.
 */
#define SUZERAIN_PERFECT_LAYERS_PHYSICAL   \
    ((a,                 1))  /* scalar */ \
    ((H0,                1))  /* scalar */ \
    ((mu,                1))  /* scalar */ \
    ((T,                 1))  /* scalar */ \
    ((u,                 2))  /* vector */

/** A Boost.Preprocessor sequence of tuples of all sampled layers. */
#define SUZERAIN_PERFECT_LAYERS      \
    SUZERAIN_PERFECT_LAYERS_WAVE     \
    SUZERAIN_PERFECT_LAYERS_PHYSICAL

    /* Compile-time totals of the number of scalars sampled at each point */
    struct nscalars { enum {
#define EXTRACT(r, data, tuple) BOOST_PP_TUPLE_ELEM(2, 1, tuple)
#define SUM(s, state, x) BOOST_PP_ADD(state, x)

        wave = BOOST_PP_SEQ_FOLD_LEFT(SUM, 0, BOOST_PP_SEQ_TRANSFORM(
                    EXTRACT,,SUZERAIN_PERFECT_LAYERS_WAVE)),

        physical = BOOST_PP_SEQ_FOLD_LEFT(SUM, 0, BOOST_PP_SEQ_TRANSFORM(
                    EXTRACT,,SUZERAIN_PERFECT_LAYERS_PHYSICAL)),

        total = BOOST_PP_SEQ_FOLD_LEFT(SUM, 0, BOOST_PP_SEQ_TRANSFORM(
                    EXTRACT,,SUZERAIN_PERFECT_LAYERS))

#undef EXTRACT
#undef SUM
    }; };

    /** Simulation time when layer profiles were obtained */
    real_t t;

    /** Type of the contiguous storage used to house all scalars */
    typedef Array<real_t, Dynamic, nscalars::total> storage_type;

    /** Contiguous storage used to house all means */
    storage_type storage;

    /**
     * Constructor setting <tt>this->t = NaN</tt>.
     * Caller will need to resize <tt>this->storage</tt> prior to use.
     */
    layers();

    /**
     * Constructor setting <tt>this->t = t</tt>.
     * Caller will need to resize <tt>this->storage</tt> prior to use.
     */
    explicit layers(real_t t);

    /**
     * Constructor setting <tt>this->t = t</tt> and preparing a zero-filled \c
     * storage containing \c Ny rows.
     */
    layers(real_t t, storage_type::Index Ny);

#define OP(r, data, tuple)                                              \
    BOOST_PP_TUPLE_ELEM(2, 0, tuple) = BOOST_PP_TUPLE_ELEM(2, 1, tuple)

    /** Compile-time offsets for each quantity within \c storage */
    struct start { enum {
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(
                OP,,SUZERAIN_SHIFTED_SUM(SUZERAIN_PERFECT_LAYERS)))
    }; };

    /** Compile-time sizes for each quantity within \c storage */
    struct size { enum {
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP,,
                    SUZERAIN_PERFECT_LAYERS))
    }; };

#undef OP

    // Declare a named, mutable "view" into storage for each quantity
#define DECLARE(r, data, tuple)                                               \
    storage_type::NColsBlockXpr<size::BOOST_PP_TUPLE_ELEM(2, 0, tuple)>::Type \
    BOOST_PP_TUPLE_ELEM(2, 0, tuple)()                                        \
    {                                                                         \
        return storage.middleCols<size::BOOST_PP_TUPLE_ELEM(2, 0, tuple)>(    \
                start::BOOST_PP_TUPLE_ELEM(2, 0, tuple));                     \
    }
    BOOST_PP_SEQ_FOR_EACH(DECLARE,,SUZERAIN_PERFECT_LAYERS)
#undef DECLARE

    // Declare a named, immutable "view" into storage for each quantity
#define DECLARE(r, data, tuple)                                                    \
    storage_type::ConstNColsBlockXpr<size::BOOST_PP_TUPLE_ELEM(2, 0, tuple)>::Type \
    BOOST_PP_TUPLE_ELEM(2, 0, tuple)() const                                       \
    {                                                                              \
        return storage.middleCols<size::BOOST_PP_TUPLE_ELEM(2, 0, tuple)>(         \
                start::BOOST_PP_TUPLE_ELEM(2, 0, tuple));                          \
    }
    BOOST_PP_SEQ_FOR_EACH(DECLARE,,SUZERAIN_PERFECT_LAYERS)
#undef DECLARE

private:

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;
};

/**
 * Using the provided state, sample the mean boundary layer profiles declared
 * in \ref layers.  This is an expensive, collective method producing valid
 * results <em>only on rank zero</em>.
 *
 * @param[in]     scenario Scenario parameters.
 * @param[in]     grid     Grid parameters.
 * @param[in]     dgrid    Pencil decomposition parameters.
 * @param[in]     cop      B-spline operator workspace.
 * @param[in,out] swave    Destroyed in the computation
 * @param[in]     t        Current simulation time
 *
 * @return Mean boundary layer layers as B-spline coefficients.
 */
layers sample_layers(
        const scenario_definition &scenario,
        const grid_specification &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        contiguous_state<4,complex_t> &swave,
        const real_t t);

} // end namespace perfect

} // end namespace suzerain

#endif // SUZERAIN_PERFECT_LAYERS_HPP
