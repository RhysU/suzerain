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

#ifndef SUZERAIN_FILTER_DEFINITION_HPP
#define SUZERAIN_FILTER_DEFINITION_HPP

/** @file
 * Classes handling filter source parameters.
 */

#include <esio/esio.h>

#include <suzerain/common.hpp>
#include <suzerain/filterop.h>
#include <suzerain/timers.h>
#include <suzerain/support/definition_base.hpp>

/** @file
 * Provides classes handling problem definition for filter
 * source term, e.g., strength coefficient.
 */

namespace suzerain {

namespace reacting {

/**
 * Holds parameters defining filter source.
 */
class filter_definition : public support::definition_base
{
public:

    /**
     * Construct an instance with all parameters set to NaN.
     * Clients can use NaN as a not-yet-specified or use-the-default value.
     */
    filter_definition();

    /**
     * Construct an instance with the given parameter values.
     *
     * @param filter_phi Filter source strength.
     */
    explicit filter_definition(const real_t filter_phi);


    /** Virtual destructor to permit use as a base class */
    virtual ~filter_definition();

    /**
     * Populate any NaN members in \c this with values from \c that.
     * Subclasses should override this method adding any desired functionality
     * either before or after invoking the superclass version.
     *
     * @param that    Instance from which information is taken.
     * @param verbose Should logging be emitted when a value is retained?
     */
    virtual void populate(
            const filter_definition& that,
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
            const filter_definition& that,
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
     * The filter source strength coefficient.
     */
    real_t filter_phi;

    /**
     * @copydoc suzerain_filterop_filter
     * @param n Vector length for \c x and \c y
     */
    int filter(const int     n,
               const double  alpha,
               const double* x,
               const int     incx,
                     double* y) const
    {
        if (SUZERAIN_UNLIKELY(!r || r->n != n)) {
            SUZERAIN_TIMER_SCOPED("filter_definition::prepare_real");
            const int err = prepare_real(n);
            if (SUZERAIN_UNLIKELY(err)) return err;
        }
        return suzerain_filterop_filter(alpha, x, incx, y, r.get());
    }

    /**
     * @copydoc suzerain_filterop_filter
     * @param n Vector length for \c x and \c y
     */
    int filter(const int     n,
               const complex_double  alpha,
               const complex_double* x,
               const int     incx,
                     complex_double* y) const
    {
        if (SUZERAIN_UNLIKELY(!z || z->n != n)) {
            SUZERAIN_TIMER_SCOPED("filter_definition::prepare_complex");
            const int err = prepare_complex(n);
            if (SUZERAIN_UNLIKELY(err)) return err;
        }
        return suzerain_filteropz_filter(alpha, x, incx, y, z.get());
    }

    /** Reset any previously initialized filtering parameters. */
    void reset();

protected:

    /** Prepare a real-valued filter workspace of length \c n. */
    int prepare_real(const int n) const;

    /** Prepare a complex-valued filter workspace of length \c n. */
    int prepare_complex(const int n) const;

    /** A filtering workspace for real-valued state. */
    mutable shared_ptr<suzerain_filterop_workspace> r;

    /** A filtering workspace for complex-valued state. */
    mutable shared_ptr<suzerain_filteropz_workspace> z;

};

} // namespace reacting

} // namespace suzerain

#endif // SUZERAIN_FILTER_DEFINITION_HPP
