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

#include <esio/esio.h>
#include <log4cxx/logger.h>
#include <suzerain/bspline.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/scenario_definition.hpp>
#include <suzerain/state.hpp>
#include <suzerain/state_impl.hpp>
#include <suzerain/storage.hpp>
#include <suzerain/timestepper.hpp>

// Introduce scalar- and complex-valued typedefs
// Currently only real_t == double is supported by many, many components
typedef double               real_t;
typedef std::complex<real_t> complex_t;

/** Field names as stored in restart files */
extern const boost::array<const char *,5> field_names;

/** Field descriptions for use in restart file comments */
extern const boost::array<const char *,5> field_descriptions;

/** Store a ScenarioDefinition in a restart file */
void store(log4cxx::LoggerPtr log,
           const esio_handle esioh,
           const suzerain::problem::ScenarioDefinition<real_t>& scenario);

/** Load a ScenarioDefinition from a restart file */
void load(log4cxx::LoggerPtr log,
          const esio_handle esioh,
          suzerain::problem::ScenarioDefinition<real_t>& scenario);

/** Store a GridDefinition in a restart file */
void store(log4cxx::LoggerPtr log,
           const esio_handle esioh,
           const suzerain::problem::GridDefinition<real_t>& grid,
           const real_t Lx,
           const real_t Lz);

/** Load a GridDefinition from a restart file */
void load(log4cxx::LoggerPtr log,
          const esio_handle esioh,
          suzerain::problem::GridDefinition<real_t>& grid);

/** Store a suzerain::bspline workspace in a restart file */
void store(log4cxx::LoggerPtr log,
           const esio_handle esioh,
           boost::shared_ptr<suzerain::bspline>& bspw /* Yes, a reference */);

/** Load a suzerain::bspline workspace from a restart file */
void load(log4cxx::LoggerPtr log,
          const esio_handle esioh,
          boost::shared_ptr<suzerain::bspline>& bspw, // Yes, a reference
          const suzerain::problem::GridDefinition<real_t>& grid);

/** Store the current simulation time information */
void store_time(log4cxx::LoggerPtr log,
                const esio_handle esioh,
                real_t time);

/** Load the current simulation time information */
void load_time(log4cxx::LoggerPtr log,
               const esio_handle esioh,
               real_t &time);

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

/** Write a complex-valued field via ESIO */
template< typename I >
inline
void complex_field_write(esio_handle h, const char *name, complex_t *field,
                         I cstride = 0, I bstride = 0, I astride = 0,
                         const char * comment = 0)
{
    using boost::numeric_cast;

    esio_field_writev(h, name, reinterpret_cast<real_t *>(field),
                      2*numeric_cast<int>(cstride),
                      2*numeric_cast<int>(bstride),
                      2*numeric_cast<int>(astride),
                      2, comment);
}

#endif // CHANNEL_COMMON_HPP
