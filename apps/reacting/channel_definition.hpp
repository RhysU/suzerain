//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#ifndef SUZERAIN_CHANNEL_DEFINITION_HPP
#define SUZERAIN_CHANNEL_DEFINITION_HPP

/** @file
 * Classes handling reacting channel flow problem scenario parameters.
 */

#include <esio/esio.h>

#include <suzerain/common.hpp>
#include <suzerain/support/definition_base.hpp>

/** @file 
 * Provides classes handling problem defintion for channel
 * flow, e.g., bulk density and momentum.
 */

namespace suzerain {

namespace reacting {

/**
 * Holds parameters defining channel flow case.
 */
class channel_definition : public support::definition_base
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
		       const real_t T_wall);


    /** Virtual destructor to permit use as a base class */
    virtual ~channel_definition();

    /**
     * Populate any NaN members in \c this with values from \c that.
     * Subclasses should override this method adding any desired functionality
     * either before or after invoking the superclass version.
     *
     * @param that    Instance from which information is taken.
     * @param verbose Should logging be emitted when a value is retained?
     */
    virtual void populate(
            const channel_definition& that,
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
            const channel_definition& that,
            const bool verbose = false);

    /**
     * Save scenario into an ESIO-based file.
     * Subclasses should override this method adding any desired functionality
     * either before or after invoking the superclass version.
     *
     * @param h Open, writable handle in which details will be saved.
     */
    virtual void save(
            const esio_handle h) const;

    /**
     * Populate scenario from an ESIO-based file.
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
    boost::program_options::options_description options_description();

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
