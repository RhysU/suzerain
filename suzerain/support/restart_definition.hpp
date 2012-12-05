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
// restart_definition.hpp: classes handling restart definitions
// $Id$

#ifndef SUZERAIN_SUPPORT_RESTART_DEFINITION_HPP
#define SUZERAIN_SUPPORT_RESTART_DEFINITION_HPP

#include <suzerain/common.hpp>
#include <suzerain/support/definition_base.hpp>

/** @file
 * Provides classes handling restart definitions, which are runtime
 * arguments used to control simulation restart behavior.
 */

namespace suzerain {

namespace support {

/**
 * Encapsulates flags related to restart file behavior for the simulation.
 * Includes the particulars around writing and archiving restart files.
 */
class restart_definition : public definition_base
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
    restart_definition(const std::string& metadata,
                       const std::string& uncommitted,
                       const std::string& destination,
                       std::size_t retain = 1,
                       real_t dt          = 0,
                       std::size_t nt     = 0);

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
    std::size_t retain;

    /**
     * The maximum amount of simulation time between writing restart files.
     */
    real_t dt;

    /**
     * The maximum number of simulation steps to take between writing restart
     * files.
     */
    std::size_t nt;

    /**
     * Save restart fields as collocation point values in physical space?
     */
    bool physical;
};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_RESTART_DEFINITION_HPP
