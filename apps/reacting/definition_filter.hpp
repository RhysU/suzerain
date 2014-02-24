//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
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

#ifndef SUZERAIN_REACTING_DEFINITION_FILTER_HPP
#define SUZERAIN_REACTING_DEFINITION_FILTER_HPP

/** @file
 * Provides classes handling problem definition for filter
 * source term, e.g., strength coefficient.
 */

#include <suzerain/common.hpp>
#include <suzerain/filterop.h>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/specification_grid.hpp>
#include <suzerain/support/definition_base.hpp>
#include <suzerain/support/esio_fwd.hpp>
#include <suzerain/support/loadable.hpp>
#include <suzerain/support/overridable.hpp>
#include <suzerain/support/populatable.hpp>
#include <suzerain/support/savable.hpp>
#include <suzerain/timers.h>

namespace suzerain {

namespace reacting {

/**
 * Holds parameters defining filter source.
 */
class definition_filter
    : public virtual support::definition_base
    , public virtual support::loadable
    , public virtual support::overridable<definition_filter>
    , public virtual support::populatable<definition_filter>
    , public virtual support::savable
{
public:

    /**
     * Construct an instance with all parameters set to NaN.
     * Clients can use NaN as a not-yet-specified or use-the-default value.
     */
    definition_filter();

    /**
     * Construct an instance with the given parameter values.
     *
     * @param filter_phi Filter source strength.
     */
    explicit definition_filter(const real_t filter_phi);

    /** @copydoc support::populatable::populate */
    virtual void populate(
            const definition_filter& that,
            const bool verbose = false);

    /** @copydoc support::overridable::override */
    virtual void override(
            const definition_filter& that,
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
        SUZERAIN_TIMER_SCOPED("definition_filter::filter(real)");

        if (SUZERAIN_UNLIKELY(!r || r->n != n)) {
            SUZERAIN_TIMER_SCOPED("definition_filter::prepare_real");
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
        SUZERAIN_TIMER_SCOPED("definition_filter::filter(complex)");

        if (SUZERAIN_UNLIKELY(!z || z->n != n)) {
            SUZERAIN_TIMER_SCOPED("definition_filter::prepare_complex");
            const int err = prepare_complex(n);
            if (SUZERAIN_UNLIKELY(err)) return err;
        }

        return suzerain_filteropz_filter(alpha, x, incx, y, z.get());
    }

    /**
     * @copydoc suzerain_filteropz_source_apply
     */
    template<typename MultiArrayX, typename MultiArrayY>
    int source_apply(const specification_grid &grid,
                     const pencil_grid &dgrid,
                     const typename MultiArrayX::element& alpha,
                     const MultiArrayX &x,
                     int ndx_x) const
    {
        SUZERAIN_TIMER_SCOPED("definition_filter::source_apply");

        if (SUZERAIN_UNLIKELY(!z || z->n != dgrid.global_wave_extent.y())) {
            SUZERAIN_TIMER_SCOPED("definition_filter::prepare_complex");
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
    int source_accumulate(const specification_grid &grid,
                          const pencil_grid &dgrid,
                          const typename MultiArrayX::element& alpha,
                          const MultiArrayX &x,
                          int ndx_x,
                          const typename MultiArrayY::element& beta,
                          MultiArrayY &y,
                          int ndx_y) const
    {
        SUZERAIN_TIMER_SCOPED("definition_filter::source_accumulate");

        assert(std::equal(x.shape()   + 1, x.shape()   + 4, y.shape()   + 1));
        assert(std::equal(x.strides() + 1, x.strides() + 4, y.strides() + 1));

        if (SUZERAIN_UNLIKELY(!z || z->n != dgrid.global_wave_extent.y())) {
            SUZERAIN_TIMER_SCOPED("definition_filter::prepare_complex");
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
                        const int n);

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

#endif // SUZERAIN_REACTING_DEFINITION_FILTER_HPP
