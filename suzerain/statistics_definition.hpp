//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// statistics_definition.hpp: classes handling statistics output definitions
// $Id$

#ifndef SUZERAIN_STATISTICS_DEFINITION_HPP
#define SUZERAIN_STATISTICS_DEFINITION_HPP

#include <suzerain/common.hpp>
#include <suzerain/problem.hpp>

/** @file
 * Provides classes handling arguments used to control statistics file output.
 */

namespace suzerain {

namespace problem {

/**
 * Encapsulates flags related to statistics output behavior for the simulation.
 * Includes the sampling rate and particulars around writing and archiving
 * statistics files.
 */
class StatisticsDefinition : public IDefinition
{
public:
    /**
     * Construct an instance with the given default values.
     * All of these can be overridden by command line options.
     *
     * @param destination  Archiving pattern to use when
     *                     committing statistics output files.
     * @param retain       Maximum number of committed statistics files
     *                     to retain.
     * @param nt           Number of simulation steps between sampling.
     * @param dt           Amount of simulation time between sampling.
     *
     * @see ESIO's esio_file_close_restart() for the semantics of
     *      \c desttemplate and \c retain.
     */
    explicit StatisticsDefinition(
            const std::string& destination,
            int retain = (1<<15),
            real_t dt  = 0,
            int nt     = 0);

    /**
     * The archiving pattern to use when committing statistics output files.
     *
     * @see ESIO's esio_file_close_restart() for the semantics of
     *      destination.
     */
    std::string destination;

    /**
     * The maximum number of committed statistics files to retain.
     */
    int retain;

    /**
     * The maximum amount of simulation time between sampling statistics.
     */
    real_t dt;

    /**
     * The maximum number of simulation steps to take between sampling
     * statistics.
     */
    int nt;
};

} // namespace problem

} // namespace suzerain

#endif // SUZERAIN_STATISTICS_DEFINITION_HPP
