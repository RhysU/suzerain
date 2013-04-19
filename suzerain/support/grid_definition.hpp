//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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
class grid_definition
    : public virtual definition_base
    , public grid_specification
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

    /** Virtual destructor to permit use as a base class */
    virtual ~grid_definition();

    /**
     * Populate any NaN members in \c this with values from \c that.
     * Subclasses should override this method adding any desired functionality
     * either before or after invoking the superclass version.
     *
     * @param that    Instance from which information is taken.
     * @param verbose Should logging be emitted when a value is retained?
     */
    virtual void populate(
            const grid_definition& that,
            const bool verbose = false);

    /**
     * Override members in \c this with non-NaN values from \c that.
     * Subclasses should override this method adding any desired functionality
     * either before or after invoking the superclass version.
     *
     * @param that    Instance from which information is taken.
     * @param verbose Should logging be emitted when an override occurs?
     */
    virtual void override(
            const grid_definition& that,
            const bool verbose = false);

    /**
     * Save grid details into an ESIO-based file.
     * Subclasses should override this method adding any desired functionality
     * either before or after invoking the superclass version.
     *
     * @param h Open, writable handle in which details will be saved.
     */
    virtual void save(
            const esio_handle h) const;

    /**
     * Populate grid details from an ESIO-based file.
     * Subclasses should override this method adding any desired functionality
     * either before or after invoking the superclass version.
     *
     * @param h       Open, readable handle from which details will be loaded.
     * @param verbose Should logging be emitted when a value is retained?
     */
    virtual void load(
            const esio_handle h,
            const bool verbose = true);

    /** @copydoc support::definition_base::options_description() */
    virtual boost::program_options::options_description options_description();

};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_GRID_DEFINITION_HPP
