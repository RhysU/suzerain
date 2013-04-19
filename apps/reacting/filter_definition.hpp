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

#ifndef SUZERAIN_FILTER_DEFINITION_HPP
#define SUZERAIN_FILTER_DEFINITION_HPP

/** @file
 * Provides classes handling problem definition for filter
 * source term, e.g., strength coefficient.
 */

#include <esio/esio.h>

#include <suzerain/common.hpp>
#include <suzerain/filterop.h>
#include <suzerain/timers.h>
#include <suzerain/support/definition_base.hpp>
#include <suzerain/grid_specification.hpp>
#include <suzerain/pencil_grid.hpp>

namespace suzerain {

namespace reacting {

/**
 * Holds parameters defining filter source.
 */
class filter_definition
    : public virtual support::definition_base
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
        SUZERAIN_TIMER_SCOPED("filter_definition::filter(real)");

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
        SUZERAIN_TIMER_SCOPED("filter_definition::filter(complex)");

        if (SUZERAIN_UNLIKELY(!z || z->n != n)) {
            SUZERAIN_TIMER_SCOPED("filter_definition::prepare_complex");
            const int err = prepare_complex(n);
            if (SUZERAIN_UNLIKELY(err)) return err;
        }

        return suzerain_filteropz_filter(alpha, x, incx, y, z.get());
    }

    /**
     * @copydoc suzerain_filteropz_source_apply
     */
    template<typename MultiArrayX, typename MultiArrayY>
    int source_apply(const grid_specification &grid,
                     const pencil_grid &dgrid,
                     const typename MultiArrayX::element& alpha,
                     const MultiArrayX &x,
                     int ndx_x) const
    {
        SUZERAIN_TIMER_SCOPED("filter_definition::source_apply");

        if (SUZERAIN_UNLIKELY(!z || z->n != dgrid.global_wave_extent.y())) {
            SUZERAIN_TIMER_SCOPED("filter_definition::prepare_complex");
            const int err = prepare_complex(dgrid.global_wave_extent.y());
            if (SUZERAIN_UNLIKELY(err)) return err;
        }

        return suzerain_filteropz_source_apply(
                alpha, x[ndx_x].origin(),
                dgrid.global_wave_extent.y(),
                grid.N.x(),
                grid.dN.x(),
                dgrid.local_wave_start.x(),
                dgrid.local_wave_end.x(),
                grid.N.z(),
                grid.dN.z(),
                dgrid.local_wave_start.z(),
                dgrid.local_wave_end.z(),
                z.get());
    }

    /**
     * @copydoc suzerain_filteropz_source_accumulate
     */
    template<typename MultiArrayX, typename MultiArrayY>
    int source_accumulate(const grid_specification &grid,
                          const pencil_grid &dgrid,
                          const typename MultiArrayX::element& alpha,
                          const MultiArrayX &x,
                          int ndx_x,
                          const typename MultiArrayY::element& beta,
                          MultiArrayY &y,
                          int ndx_y) const
    {
        SUZERAIN_TIMER_SCOPED("filter_definition::source_accumulate");

        assert(std::equal(x.shape()   + 1, x.shape()   + 4, y.shape()   + 1));
        assert(std::equal(x.strides() + 1, x.strides() + 4, y.strides() + 1));

        if (SUZERAIN_UNLIKELY(!z || z->n != dgrid.global_wave_extent.y())) {
            SUZERAIN_TIMER_SCOPED("filter_definition::prepare_complex");
            const int err = prepare_complex(dgrid.global_wave_extent.y());
            if (SUZERAIN_UNLIKELY(err)) return err;
        }

        return suzerain_filteropz_source_accumulate(
                alpha, x[ndx_x].origin(),
                beta,  y[ndx_y].origin(),
                dgrid.global_wave_extent.y(),
                grid.N.x(),
                grid.dN.x(),
                dgrid.local_wave_start.x(),
                dgrid.local_wave_end.x(),
                grid.N.z(),
                grid.dN.z(),
                dgrid.local_wave_start.z(),
                dgrid.local_wave_end.z(),
                z.get());
    }

    /** Reset any previously initialized filtering parameters. */
    void reset();

    /** Save filter operators. */
    void save_filteropz(const esio_handle h,
                        const int n,
                        const char *location = "filter_operators");

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
