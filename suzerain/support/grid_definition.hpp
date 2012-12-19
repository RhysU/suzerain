//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#ifndef SUZERAIN_SUPPORT_GRID_DEFINITION_HPP
#define SUZERAIN_SUPPORT_GRID_DEFINITION_HPP

/** @file
 * Provides classes handling problem grid definitions.
 */

#include <esio/esio.h>

#include <suzerain/common.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/grid_specification.hpp>
#include <suzerain/support/definition_base.hpp>

namespace suzerain {

namespace support {

/**
 * Upgrades a \ref grid_specification with \ref definition_base behavior.
 * This permits using the instance with \ref program_options.
 */
class grid_definition : public grid_specification, public definition_base
{
public:

    /**
     * Construct an instance containing default values intended to be
     * overwritten.  Integer values will be zeros and floating point
     * values will be NaNs.
     */
    grid_definition();

    /**
     * Construct an instance with the given default values.
     *
     * @param Lx      Physical domain extent in the X direction
     *                which is evaluated using exprparse().
     * @param Nx      Logical grid size in the X direction.
     * @param DAFx    Dealiasing factor in the X direction.
     * @param Ly      Physical domain extent in the Y direction
     *                which is evaluated using exprparse().
     * @param Ny      Logical grid size in the Y direction.
     *                This is the number of B-spline basis functions
     *                (equivalently, wall-normal degrees of freedom)
     *                to use.
     * @param k       Uniform B-spline basis order plus one.
     *                Piecewise cubics correspond to
     *                <tt>default_k == 4</tt>.
     * @param htdelta Hyperbolic tangent stretching parameter
     *                to use when computing breakpoint locations.
     * @param Lz      Physical domain extent in the Z direction
     *                which is evaluated using exprparse().
     * @param Nz      Logical grid size in the Z direction.
     * @param DAFz    Dealiasing factor in the Z direction.
     */
    grid_definition(const real_t Lx,
                    const int    Nx,
                    const real_t DAFx,
                    const real_t Ly,
                    const int    Ny,
                    const int    k,
                    const real_t htdelta,
                    const real_t Lz,
                    const int    Nz,
                    const real_t DAFz);

    /** @copydoc support::definition_base::options_description() */
    boost::program_options::options_description options_description();

};

/**
 * Save a grid_definition in an ESIO-based file.
 *
 * @param h    Open, writable handle in which details will be saved.
 * @param grid Grid to be saved.
 */
void save(const esio_handle h,
          const grid_definition& grid);

/**
 * Load a grid_definition from an ESIO-based file.
 *
 * @param h    Open, readable handle from which details will be loaded.
 * @param grid Grid to be saved.
 */
void load(const esio_handle h,
          grid_definition& grid);

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_GRID_DEFINITION_HPP
