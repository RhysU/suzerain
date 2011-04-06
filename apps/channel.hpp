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
 * channel.hpp: Channel-related functionality spanning binaries
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef CHANNEL_HPP
#define CHANNEL_HPP

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

/**
 * Contains cross-cutting functionality used within the channel binaries.
 */
namespace channel {

/** Contains basic details about the scalar state fields employed */
namespace field {

/** Contains state variable indices within state storage */
namespace ndx {

// Anonymous enum to declare our state variable storage indices.
// Update count just below if you modify this enum!
enum {
    rho,  /**< Nondimensional density */
    rhou, /**< Nondimensional streamwise momentum */
    rhov, /**< Nondimensional wall-normal momentum */
    rhow, /**< Nondimensional spanwise momentum */
    rhoe  /**< Nondimensional total energy */
};

} // end namespace ndx;

/** Contains the number of distinct state variables we track */
const std::size_t count = static_cast<std::size_t>(ndx::rhoe) + 1;

/** Field names as stored in restart files */
extern const boost::array<const char *, count> name;

} // end namespace field


/** Store a ScenarioDefinition in a restart file */
void store(const esio_handle h,
           const suzerain::problem::ScenarioDefinition<real_t>& scenario);

/** Load a ScenarioDefinition from a restart file */
void load(const esio_handle h,
          suzerain::problem::ScenarioDefinition<real_t>& scenario);

/** Store a GridDefinition in a restart file */
void store(const esio_handle h,
           const suzerain::problem::GridDefinition& grid,
           const real_t Lx,
           const real_t Lz);

/** Load a GridDefinition from a restart file */
void load(const esio_handle h,
          suzerain::problem::GridDefinition& grid);

/** Create a B-spline workspace on [left,right] per ndof, k, and htdelta */
void create(const int ndof,
            const int k,
            const double left,
            const double right,
            const double htdelta,
            boost::shared_ptr<suzerain::bspline>& b,
            boost::shared_ptr<suzerain::bsplineop>& bop);

/** Store a suzerain::bspline workspace in a restart file */
void store(const esio_handle h,
           const boost::shared_ptr<suzerain::bspline>& b,
           const boost::shared_ptr<suzerain::bsplineop>& bop);

/** Load a suzerain::bspline workspace from a restart file */
void load(const esio_handle h,
          boost::shared_ptr<suzerain::bspline>& b,
          boost::shared_ptr<suzerain::bsplineop>& bop);

/** Store the current simulation time information */
void store_time(const esio_handle h,
                real_t time);

/** Load the current simulation time information */
void load_time(const esio_handle h,
               real_t &time);

/**
 * Store the current simulation state into an open restart file
 * Only non-dealiased state content is saved.
 */
void store(const esio_handle h,
           const suzerain::NoninterleavedState<4,complex_t> &state,
           const suzerain::problem::GridDefinition& grid,
           const suzerain::pencil_grid& dgrid);

/**
 * Load the current simulation state from an open restart file.
 * Handles the very non-trivial task of adjusting the restart
 * to match the provided \c grid, \c dgrid, \c b, and \c bop
 */
void load(const esio_handle h,
          suzerain::NoninterleavedState<4,complex_t> &state,
          const suzerain::problem::GridDefinition& grid,
          const suzerain::pencil_grid& dgrid,
          const suzerain::bspline& b,
          const suzerain::bsplineop& bop);

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

/** Holds information on the \f$L^2\f$ norm of a scalar field */
struct L2 {
    real_t mean2;
    real_t fluctuating2;
    real_t total2()        { return mean2 + fluctuating2;    };
    real_t total()         { return std::sqrt(total2());     };
    real_t mean()          { return std::sqrt(mean2);        };
    real_t fluctuating()   { return std::sqrt(fluctuating2); };
};

/**
 * Compute information about the \f$L^2\f$ norm of all scalar fields.
 * See writeup/L2.tex for full details.
 */
boost::array<L2,field::count>
field_L2(const suzerain::NoninterleavedState<4,complex_t> &state,
         const suzerain::problem::ScenarioDefinition<real_t>& scenario,
         const suzerain::problem::GridDefinition& grid,
         const suzerain::pencil_grid& dgrid,
         suzerain::bspline& b);

} // end namespace channel

#endif // CHANNEL_HPP
