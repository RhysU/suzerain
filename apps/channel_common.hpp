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

#include <log4cxx/logger.h>
#include <esio/esio.h>
#include <suzerain/bspline.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/scenario_definition.hpp>

// Introduce scalar- and complex-valued typedefs
// Currently only real_t == double is supported by many, many components
typedef double               real_t;
typedef std::complex<real_t> complex_t;

void store(log4cxx::LoggerPtr log,
           const esio_handle esioh,
           const suzerain::problem::ScenarioDefinition<real_t>& scenario);

void load(log4cxx::LoggerPtr log,
          const esio_handle esioh,
          suzerain::problem::ScenarioDefinition<real_t>& scenario);

void store(log4cxx::LoggerPtr log,
           const esio_handle esioh,
           const suzerain::problem::GridDefinition<real_t>& grid,
           const real_t Lx,
           const real_t Lz);

void load(log4cxx::LoggerPtr log,
          const esio_handle esioh,
          suzerain::problem::GridDefinition<real_t>& grid);

void store(log4cxx::LoggerPtr log,
           const esio_handle esioh,
           boost::shared_ptr<suzerain::bspline>& bspw /* Yes, a reference */);

void load(log4cxx::LoggerPtr log,
          const esio_handle esioh,
          boost::shared_ptr<suzerain::bspline>& bspw, // Yes, a reference
          const suzerain::problem::GridDefinition<real_t>& grid);

#endif // CHANNEL_COMMON_HPP
