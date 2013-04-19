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

#ifndef SUZERAIN_SUPPORT_ISOTHERMAL_DEFINITION_HPP
#define SUZERAIN_SUPPORT_ISOTHERMAL_DEFINITION_HPP

/** @file
 * Provides \ref isothermal_definition.
 */

#include <esio/esio.h>

#include <suzerain/common.hpp>
#include <suzerain/isothermal_specification.hpp>
#include <suzerain/support/definition_base.hpp>

namespace suzerain {

namespace support {

// FIXME: Declare other constructors per isothermal_specification.

/**
 * Holds parameters defining isothermal boundary cases.
 */
class isothermal_definition
    : public virtual support::definition_base
    , public isothermal_specification
{
public:

    /**
     * Construct an instance with all parameters set to NaN.
     * Clients can use NaN as a not-yet-specified or use-the-default value.
     */
    isothermal_definition();

    /** Virtual destructor to permit use as a base class */
    virtual ~isothermal_definition();

    /**
     * Populate any NaN members in \c this with values from \c that.
     * Subclasses should override this method adding any desired functionality
     * either before or after invoking the superclass version.
     *
     * @param that    Instance from which information is taken.
     * @param verbose Should logging be emitted when a value is retained?
     */
    virtual void populate(
            const isothermal_definition& that,
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
            const isothermal_definition& that,
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
    virtual boost::program_options::options_description options_description();

};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_ISOTHERMAL_DEFINITION_HPP
