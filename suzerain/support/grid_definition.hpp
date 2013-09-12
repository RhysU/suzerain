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

#include <suzerain/common.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/grid_specification.hpp>
#include <suzerain/support/definition_base.hpp>
#include <suzerain/support/esio_fwd.hpp>
#include <suzerain/support/loadable.hpp>
#include <suzerain/support/overridable.hpp>
#include <suzerain/support/populatable.hpp>
#include <suzerain/support/savable.hpp>

namespace suzerain {

namespace support {

/**
 * Upgrades a \ref grid_specification with \ref definition_base behavior.
 * This permits using the instance with \ref program_options.
 */
class grid_definition
    : public virtual definition_base
    , public virtual loadable
    , public virtual overridable<grid_specification>
    , public virtual populatable<grid_specification>
    , public virtual savable
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

    /** @copydoc populatable::populate */
    virtual void populate(
            const grid_specification& that,
            const bool verbose = false);

    /** @copydoc overridable::override */
    virtual void override(
            const grid_specification& that,
            const bool verbose = false);

    /** @copydoc savable::save */
    virtual void save(
            const esio_handle h) const;

    /** @copydoc loadable::load */
    virtual void load(
            const esio_handle h,
            const bool verbose = true);

    /** @copydoc definition_base::options_description() */
    virtual boost::program_options::options_description options_description();

};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_GRID_DEFINITION_HPP
