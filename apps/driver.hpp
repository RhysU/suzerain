//--------------------------------------------------------------------------
//
// Copyright (C) 2008, 2009, 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
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
// driver.hpp: Application driver logic spanning multiple applications
// $Id$

#ifndef SUZERAIN_SUPPORT_DRIVER_HPP
#define SUZERAIN_SUPPORT_DRIVER_HPP

#include <esio/esio.h>

#include <suzerain/common.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/fftw.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/restart_definition.hpp>
#include <suzerain/signal_definition.hpp>
#include <suzerain/state.hpp>
#include <suzerain/statistics_definition.hpp>
#include <suzerain/time_definition.hpp>
#include <suzerain/timestepper.hpp>

namespace suzerain {

namespace support {

class field;

/**
 * An extensible driver class for managing a Suzerain application.
 * Intended for time-varying, three-dimensional problems.
 * Instantiate from within \c main().
 *
 * Signal handling capabilities may misbehave if multiple instances
 * are executing within the same process.
 */
class Driver
{
public:

    typedef InterleavedState<4,complex_t> linear_state_type;

    typedef ContiguousState<4,complex_t> nonlinear_state_type;

    Driver();

    virtual ~Driver();

    std::vector<support::field> fields;

    problem::GridDefinition grid;

    fftw::FFTWDefinition fftwdef;

    problem::RestartDefinition restart;

    problem::StatisticsDefinition statsdef;

    problem::TimeDefinition timedef;

    boost::shared_ptr<bspline> b;

    boost::shared_ptr<bsplineop> bop; // Collocation operators

    boost::shared_ptr<bsplineop> gop; // Galerkin L2 operators

    boost::shared_ptr<pencil_grid> dgrid;

    boost::shared_ptr<linear_state_type> state_linear;

    boost::shared_ptr<nonlinear_state_type> state_nonlinear;

    esio_handle esioh;

    /** Controls which signals trigger which processing. */
    static const problem::SignalDefinition sigdef;

    /**
     * Type of atomic locations used to track local receipt of the following
     * signal-based actions:
     *
     * \li \c 0 Output a status message
     * \li \c 1 Write a restart file
     * \li \c 2 Tear down the simulation (reactively  due to an incoming signal)
     * \li \c 3 Tear down the simulation (proactively due to --advance_wt limit)
     * \li \c 4 Compute and write a statistics file
     */
    typedef boost::array<
            volatile sig_atomic_t, 5
        > atomic_signal_received_t;

    /** Atomic locations used to track local signal receipt. */
    static atomic_signal_received_t atomic_signal_received;

    /**
     * When \c true, any time advance should be stopped as soon
     * as reasonably possible.
     */
    bool soft_teardown;

    /**
     * Routine to output status, generally called via the TimeController.
     *
     * Invokes \ref log_status_bulk, \ref log_status_L2, \ref
     * log_status_boundary_state, and \ref log_status_hook.
     */
    virtual bool log_status(
            const real_t t,
            const std::size_t nt);

    /** Log messages containing mean L2 and RMS fluctuation information. */
    virtual void log_status_L2(
            const std::string& timeprefix,
            const char * const name_L2  = "L2.mean",
            const char * const name_rms = "rms.fluct");

    /** Log messages containing bulk quantities. */
    virtual void log_status_bulk(
            const std::string& timeprefix);

    /**
     * Log messages containing specific state quantities at the upper and lower
     * boundaries.  Density is reported as-is.  All other scalars are divided
     * by density.
     */
    virtual void log_status_specific_boundary_state(
            const std::string& timeprefix);

    /**
     * Routine to save a restart file, generally called via the TimeController.
     *
     * The restart saves the data in \ref state_linear.
     * The data in \ref state_nonlinear is destroyed by this call.
     */
    virtual bool save_restart(
            real_t t,
            size_t nt);

protected:

    /**
     * Hook permitting subclasses to output additional status information.
     * Invoked at the end of \ref log_status.  Returning \c false causes the
     * TimeController to halt.
     */
    virtual bool log_status_hook(
            const std::string& timeprefix,
            const real_t t,
            const std::size_t nt);

private:

    /**
     * Flag used to control whether \ref log_status_L2 shows headers.
     * The default implementation disables headers after the first invocation.
     */
    bool log_status_L2_show_header;

    /**
     * Flag used to control whether \ref log_status_bulk shows headers.
     * The default implementation disables headers after the first invocation.
     */
    bool log_status_bulk_show_header;

    /** Signal handler which mutates \c atomic_signal_received. */
    static void process_signal(const int sig);

    /** Tracks last time a status line was output */
    std::size_t last_status_nt;

    /** Tracks last time a restart file was written successfully */
    std::size_t last_restart_saved_nt;

    /**
     * Type of non-atomic locations used to track global receipt of the
     * same actions as \ref atomic_signal_received_t.
     */
    typedef boost::array<
            int, atomic_signal_received_t::static_size
        > signal_received_t;

    /** Non-atomic locations used to track global signal receipt. */
    signal_received_t signal_received;

};

} // end namespace support

} // end namespace suzerain

#endif // SUZERAIN_SUPPORT_DRIVER_HPP
