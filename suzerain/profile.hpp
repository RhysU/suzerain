//--------------------------------------------------------------------------
//
// Copyright (C) 2013-2014 Rhys Ulerich
// Copyright (C) 2013-2014 The PECOS Development Team
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

#ifndef SUZERAIN_PROFILE_HPP
#define SUZERAIN_PROFILE_HPP

/** @file
 * Sampling logistics for boundary layer and channel solution profiles.
 */

#include <suzerain/common.hpp>
#include <suzerain/profile.h>

namespace suzerain {

// Forward declarations
class samples;
class summary;

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

    /* Compile-time totals of the number of scalars sampled at each point */
    struct nscalars { enum {
#define EXTRACT(r, data, tuple) BOOST_PP_TUPLE_ELEM(2, 1, tuple)
#define SUM(s, state, x) BOOST_PP_ADD(state, x)

        wave = BOOST_PP_SEQ_FOLD_LEFT(SUM, 0, BOOST_PP_SEQ_TRANSFORM(
                    EXTRACT,,SUZERAIN_PROFILE_WAVE)),

        physical = BOOST_PP_SEQ_FOLD_LEFT(SUM, 0, BOOST_PP_SEQ_TRANSFORM(
                    EXTRACT,,SUZERAIN_PROFILE_PHYSICAL)),

        total = BOOST_PP_SEQ_FOLD_LEFT(SUM, 0, BOOST_PP_SEQ_TRANSFORM(
                    EXTRACT,,SUZERAIN_PROFILE))

#undef EXTRACT
#undef SUM
    }; };

    /** Type of the contiguous storage used to house all scalars */
    typedef Array<real_t, Dynamic, nscalars::total, ColMajor> storage_type;

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
                OP,,SUZERAIN_SHIFTED_SUM(SUZERAIN_PROFILE)))
    }; };

    /**
     * Provide access to contiguous subregions within storage.
     * @{
     */
    storage_type::NColsBlockXpr <nscalars::wave>::Type wave()
    { return storage.middleCols<nscalars::wave>(start::wave); }

    storage_type::NColsBlockXpr <nscalars::physical>::Type physical()
    { return storage.middleCols<nscalars::physical>(start::physical); }

    storage_type::ConstNColsBlockXpr<nscalars::wave>::Type wave() const
    { return storage.middleCols<nscalars::wave>(start::wave); }

    storage_type::ConstNColsBlockXpr<nscalars::physical>::Type physical() const
    { return storage.middleCols<nscalars::physical>(start::physical); }
    /** @} */

    /** Compile-time sizes for each quantity within \c storage */
    struct size { enum {
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP,, SUZERAIN_PROFILE))
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
    BOOST_PP_SEQ_FOR_EACH(DECLARE,,SUZERAIN_PROFILE)
#undef DECLARE

    // Declare a named, immutable "view" into storage for each quantity
#define DECLARE(r, data, tuple)                                                    \
    storage_type::ConstNColsBlockXpr<size::BOOST_PP_TUPLE_ELEM(2, 0, tuple)>::Type \
    BOOST_PP_TUPLE_ELEM(2, 0, tuple)() const                                       \
    {                                                                              \
        return storage.middleCols<size::BOOST_PP_TUPLE_ELEM(2, 0, tuple)>(         \
                start::BOOST_PP_TUPLE_ELEM(2, 0, tuple));                          \
    }
    BOOST_PP_SEQ_FOR_EACH(DECLARE,,SUZERAIN_PROFILE)
#undef DECLARE

    /**
     * Copy information from a \ref samples instance.
     * May permit avoiding \ref sample_profile calls in some circumstances.
     */
    profile& operator=(const samples &q);

    /**
     * Copy information from a \ref summary instance.
     * May permit avoiding \ref sample_profile calls in some circumstances.
     */
    profile& operator=(const summary &q);

};

} // namespace suzerain

#endif // SUZERAIN_PROFILE_HPP
