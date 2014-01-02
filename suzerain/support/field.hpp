//--------------------------------------------------------------------------
//
// Copyright (C) 2008-2014 Rhys Ulerich
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

#ifndef SUZERAIN_SUPPORT_FIELD_HPP
#define SUZERAIN_SUPPORT_FIELD_HPP

/** @file
 * Logic for saving and loading distributed state fields.
 */

#include <suzerain/common.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/exprparse.hpp>
#include <suzerain/grid_specification.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/state_fwd.hpp>
#include <suzerain/support/esio_fwd.hpp>

namespace suzerain {

namespace support {

/**
 * Collects details about scalar-valued fields.
 * For example, spanwise momentum.
 */
class field
{
public:

    /**
     * A brief string used in status displays.
     * For example, "rho_w".
     */
    std::string identifier;

    /**
     * Human-readable text used to generate descriptions.
     * For example, "spanwise momentum".
     */
    std::string description;

    /**
     * The ESIO location (that is, HDF5 dataset) in which the field is stored
     * in files.  Often this will just be \c identifier.
     */
    std::string location;

};

/**
 * Save the current simulation conserved state as expansion coefficients into
 * an open restart file.   Only non-dealiased, conserved state is saved as
 * "wave space" coefficients.  This is the most efficient and flexible way to
 * save state to disk.
 */
void save_coefficients(
        const esio_handle h,
        const std::vector<field> &fields,
        const contiguous_state<4,complex_t> &swave,
        const grid_specification& grid,
        const pencil_grid& dgrid);

/**
 * Load the current simulation state from an open coefficient-based restart
 * file.  Handles the non-trivial task of adjusting the restart to match the
 * provided \c grid, \c dgrid, \c b, and \c cop.
 */
void load_coefficients(
        const esio_handle h,
        const std::vector<field> &fields,
        contiguous_state<4,complex_t> &state,
        const grid_specification& grid,
        const pencil_grid& dgrid,
        const bsplineop& cop,
        const bspline& b);

/**
 * Save the current simulation state as collocation point values into an open
 * restart file.  Note that <tt>state</tt>'s contents are destroyed.
 * Collocation point values required only for dealiasing purposes <i>are</i>
 * stored but only those points informed by non-dealiased state.  This method
 * is less efficient and the restart data less flexible than that produced by
 * save_coefficients().
 */
void save_collocation_values(
        const esio_handle h,
        const std::vector<field> &fields,
        contiguous_state<4,complex_t>& swave,
        const grid_specification& grid,
        const pencil_grid& dgrid,
        const bsplineop& cop);

/**
 * Load the current simulation state from an open collocation point value
 * restart file.  Cannot handle interpolating onto a different grid.
 */
void load_collocation_values(
        const esio_handle h,
        const std::vector<field> &fields,
        contiguous_state<4,complex_t>& state,
        const grid_specification& grid,
        const pencil_grid& dgrid,
        const bsplineop& cop,
        const bspline& b);

} // end namespace support

} // end namespace suzerain

#endif // SUZERAIN_SUPPORT_FIELD_HPP
