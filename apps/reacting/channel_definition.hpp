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

#ifndef SUZERAIN_CHANNEL_DEFINITION_HPP
#define SUZERAIN_CHANNEL_DEFINITION_HPP

/** @file
 * Provides classes handling problem definition for channel
 * flow, e.g., bulk density and momentum.
 */

#include <esio/esio.h>

#include <suzerain/common.hpp>
#include <suzerain/support/definition_base.hpp>
#include <suzerain/support/loadable.hpp>
#include <suzerain/support/overridable.hpp>
#include <suzerain/support/populatable.hpp>
#include <suzerain/support/savable.hpp>

namespace suzerain {

namespace reacting {

/**
 * Holds parameters defining channel flow case.
 */
class channel_definition
    : public virtual support::definition_base
    , public virtual support::loadable
    , public virtual support::overridable<channel_definition>
    , public virtual support::populatable<channel_definition>
    , public virtual support::savable
{
public:

    /**
     * Construct an instance with all parameters set to NaN.
     * Clients can use NaN as a not-yet-specified or use-the-default value.
     */
    channel_definition();

    /**
     * Construct an instance with the given parameter values.
     *
     * @param bulk_rho   Bulk density target.
     * @param bulk_rho_u Bulk streamwise momentum target.
     */
    channel_definition(const real_t bulk_rho,
                       const real_t bulk_rho_u,
                       const real_t T_wall,
                       const std::vector<real_t> wall_mass_fractions);


    /** @copydoc support::populatable::populate */
    virtual void populate(
            const channel_definition& that,
            const bool verbose = false);

    /** @copydoc support::overridable::override */
    virtual void override(
            const channel_definition& that,
            const bool verbose = false);

    /** @copydoc support::savable::save */
    virtual void save(
            const esio_handle h) const;

    /** @copydoc support::loadable::load */
    virtual void load(
            const esio_handle h,
            const bool verbose = true);

    /** @copydoc support::definition_base::options_description() */
    virtual boost::program_options::options_description options_description();

    /**
     * The bulk density used as a target for integral constraints.
     */
    real_t bulk_rho;

    /**
     * The bulk streamwise momentum used as a target for integral constraints.
     */
    real_t bulk_rho_u;

    /**
     * The desired wall temperature.
     */
    real_t T_wall;

    /**
     * The desired wall mass fractions
     */
    std::vector<real_t> wall_mass_fractions;

};

} // namespace reacting

} // namespace suzerain

#endif // SUZERAIN_CHANNEL_DEFINITION_HPP
