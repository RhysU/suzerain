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

    /** Log messages containing mean L2 and RMS fluctuation information. */
    virtual void information_L2(
            const std::string& prefix,
            const char * const name_L2  = "L2.mean",
            const char * const name_rms = "rms.fluct")
        = 0;

    /** Log messages containing bulk quantities. */
    virtual void information_bulk(
            const std::string& prefix)
        = 0;

    /** Log messages containing state quantities at the boundaries. */
    virtual void information_boundary_state(
            const std::string& prefix)
        = 0;

    /**
     * Hook permitting subclasses to output additional status information.
     * Returning \c false causes the TimeController to halt.
     */
    virtual bool information_extended(
            const std::string& prefix,
            const real_t simulation_time,
            const std::size_t nt)
        = 0;

    /**
     * Routine to output status, generally via the TimeController.
     *
     * Invokes \ref information_bulk, \ref information_L2, \ref
     * information_boundary_state, and \ref information_extended.
     */
    bool log_status(
            const real_t t,
            const std::size_t nt);

protected:

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

    /** Atomic locations used to track local signal receipt. */
    static atomic_signal_received_t atomic_signal_received;

    /** Flag used to indicate early time advancement stopping is legit. */
    bool soft_teardown;

    /** Flag used to control whether information_L2 shows headers. */
    bool show_header_information_L2;

    /** Flag used to control whether information_bulk shows headers. */
    bool show_header_information_bulk;

private:

    /** Signal handler which mutates \c atomic_signal_received. */
    static void process_signal(int sig);

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
