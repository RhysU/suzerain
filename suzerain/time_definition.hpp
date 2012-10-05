//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// time_definition.hpp: classes handling advance definitions
// $Id$

#ifndef SUZERAIN_TIME_DEFINITION_HPP
#define SUZERAIN_TIME_DEFINITION_HPP

#include <suzerain/common.hpp>
#include <suzerain/problem.hpp>

/** @file
 * Provides classes handling time advancement settings.
 */

namespace suzerain {

namespace problem {

/**
 * Encapsulates flags related to time advancement.  Includes details on how
 * far the simulation should be advanced as long as how frequently status
 * updates should occur.
 */
class TimeDefinition : public IDefinition
{
public:
    /**
     * Construct an instance with the given default values.  All of these can
     * be overridden by command line options.  Values not appearing here
     * <i>will</i> appear on the command line with integers and floats
     * defaulting to zero and NaN, respectively.
     *
     * @param advance_dt  Maximum amount of physical time to advance
     *                    the simulation.
     * @param advance_nt  Maximum number of discrete time steps to
     *                    advance the simulation.
     * @param advance_wt  Maximum amount of wall time to advance
     *                    the simulation.
     *                    advance the simulation.
     * @param status_dt   Maximum physical time between status updates.
     * @param status_nt   Maximum number of discrete time steps between
     *                    status updates.
     * @param min_dt      Minimum allowable physically-driven time step.
     * @param max_dt      Maximum allowable physically-driven time step.
     */
    TimeDefinition(double advance_dt,
                   int    advance_nt,
                   double advance_wt,
                   double status_dt,
                   int    status_nt,
                   double min_dt,
                   double max_dt);

    /**
     * Construct an instance with the given default values.  All of these can
     * be overridden by command line options.  Values not appearing here
     * <i>will not</i> appear on the command line with integers and floats
     * defaulting to zero and NaN, respectively.
     *
     * @param evmagfactor Safety factor used to pre-multiply a scheme's
     *                    maximum pure real and imaginary eigenvalue
     *                    magnitudes.  Usually in <tt>(0,1]</tt>, this
     *                    increases the conservativeness of the time stepping.
     */
    TimeDefinition(const char * evmagfactor);

    /** Maximum amount of physical time to advance the simulation. */
    double advance_dt;

    /** Maximum number of discrete time steps to advance the simulation. */
    int advance_nt;

    /** Maximum amount of wall time to advance the simulation. */
    double advance_wt;

    /** Maximum physical time between status updates. */
    double status_dt;

    /** Maximum number of discrete time steps between status updates. */
    int status_nt;

    /** Minimum allowable physically-driven time step */
    double min_dt;

    /** Maximum allowable physically-driven time step */
    double max_dt;

    /**
     * Factor used to pre-multiply a scheme's pure real and imaginary
     * eigenvalue magnitudes.
     *
     * @see suzerain::timestepper::lowstorage::Method for details.
     */
    double evmagfactor;

private:

    /** Prepare one-time options for time advancement */
    void initialize_advancement(double default_advance_dt,
                                int    default_advance_nt,
                                double default_advance_wt,
                                double default_status_dt,
                                int    default_status_nt,
                                double default_min_dt,
                                double default_max_dt);

    /** Prepare repeatedly-used options similar to scenario parameters */
    void initialize_scenario(const char * default_evmagfactor);
};

} // namespace problem

} // namespace suzerain

#endif // SUZERAIN_TIME_DEFINITION_HPP
