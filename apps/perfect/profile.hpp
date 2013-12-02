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

#ifndef SUZERAIN_PERFECT_PROFILE_HPP
#define SUZERAIN_PERFECT_PROFILE_HPP

/** @file
 * Sampling logistics for solution profiles.
 */

#include <suzerain/common.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/state_fwd.hpp>

// Forward declarations
struct suzerain_bl_local;
struct suzerain_bl_pg;
struct suzerain_bl_qoi;
struct suzerain_bl_reynolds;
struct suzerain_bl_thicknesses;
struct suzerain_bl_viscous;
struct suzerain_channel_local;
struct suzerain_channel_qoi;
struct suzerain_channel_viscous;

namespace suzerain {

// Forward declarations
class grid_specification;
class largo_specification;
class pencil_grid;

namespace perfect {

// Forward declarations
class scenario_definition;
class quantities;

/**
 * Encapsulates solution profiles necessary for using \ref bl.h
 * or \ref channel.h.
 *
 * Samples of each quantity are made available through two-dimensional
 * column-major arrays.  The row index iterates over wall-normal B-spline
 * coefficients and the column index iterates over tensor indices.  Scalars
 * have only a single tensor index.  Vector quantities have two indices
 * corresponding to the streamwise x and wall-normal y directions.  The
 * homogeneous spanwise direction z is not reported.
 */
class profile
{
public:

/**
 * A Boost.Preprocessor sequence of tuples of quantities computed in wave
 * space.
 */
#define SUZERAIN_PERFECT_PROFILE_WAVE             \
    ((rho,                      1)) /* scalar */ \
    ((rho_u,                    2)) /* vector */

/**
 * A Boost.Preprocessor sequence of tuples of quantities computed in physical
 * space.
 */
#define SUZERAIN_PERFECT_PROFILE_PHYSICAL   \
    ((a,                 1))  /* scalar */ \
    ((H0,                1))  /* scalar */ \
    ((mu,                1))  /* scalar */ \
    ((T,                 1))  /* scalar */ \
    ((u,                 2))  /* vector */

/** A Boost.Preprocessor sequence of tuples of all sampled quantity profiles. */
#define SUZERAIN_PERFECT_PROFILE      \
    SUZERAIN_PERFECT_PROFILE_WAVE     \
    SUZERAIN_PERFECT_PROFILE_PHYSICAL

    /* Compile-time totals of the number of scalars sampled at each point */
    struct nscalars { enum {
#define EXTRACT(r, data, tuple) BOOST_PP_TUPLE_ELEM(2, 1, tuple)
#define SUM(s, state, x) BOOST_PP_ADD(state, x)

        wave = BOOST_PP_SEQ_FOLD_LEFT(SUM, 0, BOOST_PP_SEQ_TRANSFORM(
                    EXTRACT,,SUZERAIN_PERFECT_PROFILE_WAVE)),

        physical = BOOST_PP_SEQ_FOLD_LEFT(SUM, 0, BOOST_PP_SEQ_TRANSFORM(
                    EXTRACT,,SUZERAIN_PERFECT_PROFILE_PHYSICAL)),

        total = BOOST_PP_SEQ_FOLD_LEFT(SUM, 0, BOOST_PP_SEQ_TRANSFORM(
                    EXTRACT,,SUZERAIN_PERFECT_PROFILE))

#undef EXTRACT
#undef SUM
    }; };

    /** Type of the contiguous storage used to house all scalars */
    typedef Array<real_t, Dynamic, nscalars::total> storage_type;

    /** Contiguous storage used to house all quantity profiles */
    storage_type storage;

    /** Default constructor. Resize <tt>this->storage</tt> prior to use. */
    profile();

    /** Constructor prepares zero-filled \c storage containing \c Ny rows. */
    explicit profile(storage_type::Index Ny);

#define OP(r, data, tuple)                                              \
    BOOST_PP_TUPLE_ELEM(2, 0, tuple) = BOOST_PP_TUPLE_ELEM(2, 1, tuple)

    /** Compile-time offsets for each quantity within \c storage */
    struct start { enum {
        wave     = 0,                      // Start of wave block
        physical = wave + nscalars::wave,  // Start of physical block
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(
                OP,,SUZERAIN_SHIFTED_SUM(SUZERAIN_PERFECT_PROFILE)))
    }; };

    /** Compile-time sizes for each quantity within \c storage */
    struct size { enum {
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP,,
                    SUZERAIN_PERFECT_PROFILE))
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
    BOOST_PP_SEQ_FOR_EACH(DECLARE,,SUZERAIN_PERFECT_PROFILE)
#undef DECLARE

    // Declare a named, immutable "view" into storage for each quantity
#define DECLARE(r, data, tuple)                                                    \
    storage_type::ConstNColsBlockXpr<size::BOOST_PP_TUPLE_ELEM(2, 0, tuple)>::Type \
    BOOST_PP_TUPLE_ELEM(2, 0, tuple)() const                                       \
    {                                                                              \
        return storage.middleCols<size::BOOST_PP_TUPLE_ELEM(2, 0, tuple)>(         \
                start::BOOST_PP_TUPLE_ELEM(2, 0, tuple));                          \
    }
    BOOST_PP_SEQ_FOR_EACH(DECLARE,,SUZERAIN_PERFECT_PROFILE)
#undef DECLARE

    /**
     * Copy information from a \ref quantities instance.
     * May permit avoiding \ref sample_profile calls in some circumstances.
     */
    profile& operator=(const quantities &q);

};

/**
 * Using the provided state, sample the mean solution profiles declared
 * in \ref profile.  This is a mildly expensive, collective method.
 *
 * @param[in]     scenario Scenario parameters.
 * @param[in]     grid     Grid parameters.
 * @param[in]     dgrid    Pencil decomposition parameters.
 * @param[in]     cop      B-spline operator workspace.
 * @param[in,out] swave    Destroyed in the computation
 *
 * @return Mean quantity profiles as B-spline coefficients.
 */
profile sample_profile(
        const scenario_definition &scenario,
        const grid_specification &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        contiguous_state<4,complex_t> &swave);

/**
 * Use the boundary layer information in \c prof and possibly base flow
 * information in \c sg to compute many quantities of interest.  This is
 * a purely local computation requiring no communication.
 *
 * @param[in]  prof     Profile information from \ref sample_profile().
 * @param[in]  scenario Scenario of interest.
 * @param[in]  sg       Slow growth definition optionally in use
 *                      which provides base flow details for
 *                      streamwise pressure and velocity gradients.
 * @param[out] wall     Populated on return.
 * @param[out] viscous  Populated on return.
 * @param[out] edge     Populated on return.
 * @param[out] thick    Populated on return.
 * @param[out] reynolds Populated on return.
 * @param[out] qoi      Populated on return.
 * @param[out] pg       Populated on return.
 */
void summarize_boundary_layer_nature(
        const profile &prof,
        const scenario_definition &scenario,
        const shared_ptr<largo_specification> &sg,
        bspline &b,
        suzerain_bl_local       &wall,
        suzerain_bl_viscous     &viscous,
        suzerain_bl_local       &edge,
        suzerain_bl_thicknesses &thick,
        suzerain_bl_reynolds    &reynolds,
        suzerain_bl_qoi         &qoi,
        suzerain_bl_pg          &pg);

/**
 * Use the boundary layer information in \c prof to compute many quantities of
 * interest.  This is a purely local computation requiring no communication.
 *
 * @param[in]  prof     Profile information from \ref sample_profile().
 * @param[in]  scenario Scenario of interest.
 * @param[out] wall     Populated on return.
 * @param[out] viscous  Populated on return.
 * @param[out] center   Populated on return.
 * @param[out] qoi      Populated on return.
 */
void summarize_channel_nature(
        const profile &prof,
        const scenario_definition &scenario,
        bspline &b,
        suzerain_channel_local   &wall,
        suzerain_channel_viscous &viscous,
        suzerain_channel_local   &center,
        suzerain_channel_qoi     &qoi);

} // end namespace perfect

} // end namespace suzerain

#endif // SUZERAIN_PERFECT_PROFILE_HPP
