//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// restart_definition.hpp: classes handling restart definitions
// $Id$

#ifndef SUZERAIN_RESTART_DEFINITION_HPP
#define SUZERAIN_RESTART_DEFINITION_HPP

#include <suzerain/common.hpp>
#include <suzerain/problem.hpp>

/** @file
 * Provides classes handling restart definitions, which are runtime
 * arguments used to control simulation restart behavior.
 */

namespace suzerain {

namespace problem {

/**
 * Encapsulates flags related to restart file behavior for the simulation.
 * Includes the particulars around writing and archiving restart files.
 */
class RestartDefinition : public IDefinition
{
public:
    /**
     * Construct an instance with the given default values.
     * All of these can be overridden by command line options.
     *
     * @param metadata     Restart file path to use when saving
     *                     common restart file metadata.
     * @param uncommitted  Restart file path to use when saving
     *                     uncommitted restart data.
     * @param destination  Restart archiving pattern to use when
     *                     committing restart files.
     * @param retain       Maximum number of committed restart files
     *                     to retain.
     * @param nt           Number of simulation steps between restart writes.
     * @param dt           Amount of simulation time between restart writes.
     *
     * @see ESIO's esio_file_close_restart() for the semantics of
     *      \c destination and \c retain.
     */
    RestartDefinition(const std::string& metadata,
                      const std::string& uncommitted,
                      const std::string& destination,
                      int retain = 1,
                      double dt  = 0,
                      int nt     = 0);

    /**
     * The restart file path to use when saving common restart file
     * metadata.
     */
    std::string metadata;

    /**
     * The file path to use when saving uncommitted restart data.
     */
    std::string uncommitted;

    /**
     * The archiving pattern to use when committing restart files.
     *
     * @see ESIO's esio_file_close_restart() for the semantics of
     *      destination.
     */
    std::string destination;

    /**
     * The maximum number of committed restart files to retain.
     */
    int retain;

    /**
     * The maximum amount of simulation time between writing restart files.
     */
    double dt;

    /**
     * The maximum number of simulation steps to take between writing restart
     * files.
     */
    int nt;

    /**
     * Save restart fields as collocation point values in physical space?
     */
    bool physical;
};

} // namespace problem

} // namespace suzerain

#endif // SUZERAIN_RESTART_DEFINITION_HPP
