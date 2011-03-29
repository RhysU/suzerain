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
#include <suzerain/timestepper.hpp>

#include "precision.hpp"

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

/**
 * Store the current simulation state into an open restart file
 * Only non-dealiased state content is saved.
 */
void store(const esio_handle esioh,
           const suzerain::NoninterleavedState<4,complex_t> &state,
           const suzerain::problem::GridDefinition& grid,
           const suzerain::pencil_grid& dgrid);

/**
 * Load the current simulation state from an open restart file.
 * Handles the very non-trivial task of adjusting the restart
 * to match the provided \c grid, \c dgrid, and \c bspw.
 */
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
