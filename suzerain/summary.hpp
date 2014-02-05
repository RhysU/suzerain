//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
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

#ifndef SUZERAIN_SUMMARY_HPP
#define SUZERAIN_SUMMARY_HPP

/** @file
 * Named storage for summarizing sampled quantity components and their
 * derivatives of interest for 3D Navier--Stokes.
 */

#include <suzerain/common.hpp>
#include <suzerain/summary.h>

namespace suzerain {

/**
 * Summarizes quantity components from \ref suzerain::sample and their
 * derivatives on collocation points.  Instances are meant as a human-readable,
 * B-spline-agnostic data interchange and post-processing contiguous storage
 * format.  Accordingly, easy access to names and descriptions are provided.
 *
 * \internal Many memory overhead and implementation consistency issues have
 * been traded for the headache of reading Boost.Preprocessor-based logic.
 */
class summary
{
public:

    /* Compile-time scalar quantity component counts  */
    struct nscalars { enum {

        grid        = BOOST_PP_SEQ_SIZE(SUZERAIN_SUMMARY_GRID),
        sampled     = BOOST_PP_SEQ_SIZE(SUZERAIN_SUMMARY_SAMPLED),
        sampled__y  = BOOST_PP_SEQ_SIZE(SUZERAIN_SUMMARY_SAMPLED__Y),
        sampled__yy = BOOST_PP_SEQ_SIZE(SUZERAIN_SUMMARY_SAMPLED__YY),
        total       = SUZERAIN_SUMMARY_COUNT

    }; };

    /** Type of the underlying storage mechanism. */
    typedef Array<real_t, Dynamic, nscalars::total> storage_type;

    /**
     * Contiguous storage used to house all data.
     * Notice each scalar comprises a contiguous column.
     */
    storage_type storage;

    /** Compile-time offsets for each component within \c storage */
    struct offset { enum {

        grid        = 0,
        sampled     = grid       + nscalars::grid,
        sampled__y  = sampled    + nscalars::sampled,
        sampled__yy = sampled__y + nscalars::sampled__y  // Comma below

#define DECLARE(data, name, description, offset) , name = offset
        SUZERAIN_SUMMARY_FOR_EACH(DECLARE,)
#undef DECLARE

    }; };

    /** Provides human-readable descriptions indexed on \ref offset */
    static const char * description[nscalars::total];

    /** Caller will need to resize <tt>this->storage</tt> prior to use. */
    summary();

    /** Constructor preparing zero-filled storage containing \c Ny rows. */
    explicit summary(storage_type::Index Ny);

    /**
     * Provide access to contiguous subregions within storage.
     * @{
     */
    storage_type::NColsBlockXpr<nscalars::grid>::Type grid()
    { return storage.middleCols<nscalars::grid>(offset::grid); }

    storage_type::NColsBlockXpr<nscalars::sampled>::Type sampled()
    { return storage.middleCols<nscalars::sampled>(offset::sampled); }

    storage_type::NColsBlockXpr<nscalars::sampled__y>::Type sampled__y()
    { return storage.middleCols<nscalars::sampled__y>(offset::sampled__y); }

    storage_type::NColsBlockXpr<nscalars::sampled__yy>::Type sampled__yy()
    { return storage.middleCols<nscalars::sampled__yy>(offset::sampled__yy); }

    storage_type::ConstNColsBlockXpr<nscalars::grid>::Type grid() const
    { return storage.middleCols<nscalars::grid>(offset::grid); }

    storage_type::ConstNColsBlockXpr<nscalars::sampled>::Type sampled() const
    { return storage.middleCols<nscalars::sampled>(offset::sampled); }

    storage_type::ConstNColsBlockXpr<nscalars::sampled__y>::Type sampled__y() const
    { return storage.middleCols<nscalars::sampled__y>(offset::sampled__y); }

    storage_type::ConstNColsBlockXpr<nscalars::sampled__yy>::Type sampled__yy() const
    { return storage.middleCols<nscalars::sampled__yy>(offset::sampled__yy); }
    /** @} */

    // Declare mutable column views into storage for each component
#define DECLARE(data, name, description, offset) \
    storage_type::ColXpr name() { return storage.col(offset); }
    SUZERAIN_SUMMARY_FOR_EACH(DECLARE,)
#undef DECLARE

    // Declare immutable column views into storage for each component
#define DECLARE(data, name, description, offset) \
    storage_type::ConstColXpr name() const { return storage.col(offset); }
    SUZERAIN_SUMMARY_FOR_EACH(DECLARE,)
#undef DECLARE

};

} // end namespace suzerain

#endif // SUZERAIN_SUMMARY_HPP
